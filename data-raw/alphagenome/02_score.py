"""AlphaGenome pilot benchmark, step 2: score the targets of 01_targets.R.

For every gene (one 1 Mb window) and every condition it requests splice-site
probabilities, splice-site usage, splice junctions and RNA-seq coverage, and
keeps only the values at the annotated targets:

  scores/sites/<gene>.parquet      splice-site tracks at each annotated site,
                                   at offsets -2..2 around it (the offset that
                                   holds the signal is settled in step 3)
  scores/regions/<gene>.parquet    mean RNA-seq coverage over each exon/intron
  scores/junctions/<gene>.parquet  every predicted junction inside the gene
  scores/tracks.parquet            metadata of the tracks, per output type

Coordinates in the outputs are those of the inputs (1-based, inclusive) for
sites and regions, and AlphaGenome's own (0-based start, end) for junctions.
Nothing is aggregated over tracks here.

The API key is read from the environment variable ALPHAGENOME_API_KEY.

Usage (from data-raw/alphagenome/):
  pixi run python 02_score.py --work-dir <work_dir> --list-conditions
  pixi run python 02_score.py --work-dir <work_dir> --smoke
  pixi run python 02_score.py --work-dir <work_dir> \
      --conditions UBERON:0002107 EFO:0002067
"""

import argparse
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

OUTPUTS = {
    "splice_sites": dna_client.OutputType.SPLICE_SITES,
    "splice_site_usage": dna_client.OutputType.SPLICE_SITE_USAGE,
    "splice_junctions": dna_client.OutputType.SPLICE_JUNCTIONS,
    "rna_seq": dna_client.OutputType.RNA_SEQ,
}
SITE_OFFSETS = (-2, -1, 0, 1, 2)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--work-dir", required=True, type=pathlib.Path)
    parser.add_argument(
        "--conditions",
        nargs="*",
        default=None,
        help="Condition curies (see --list-conditions); default: no condition.",
    )
    parser.add_argument("--list-conditions", action="store_true")
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Score the first gene only and describe the output.",
    )
    parser.add_argument("--max-genes", type=int, default=None)
    parser.add_argument("--retries", type=int, default=4)
    return parser.parse_args()


def metadata_table(track_data, output, condition):
    meta = track_data.metadata.reset_index(drop=True).copy()
    meta.insert(0, "track", np.arange(len(meta)))
    meta.insert(0, "condition", condition)
    meta.insert(0, "output", output)
    return meta.astype({c: str for c in meta.columns if meta[c].dtype == object})


def site_scores(track_data, output, sites, window_start):
    """Track values at each site and offset; long table."""
    if track_data.resolution != 1:
        raise ValueError(f"{output}: expected 1 bp resolution")
    values = track_data.values
    idx0 = sites["pos"].to_numpy() - 1 - window_start
    frames = []
    for offset in SITE_OFFSETS:
        idx = idx0 + offset
        ok = (idx >= 0) & (idx < values.shape[0])
        block = values[idx[ok], :]
        n_sites, n_tracks = block.shape
        frames.append(
            pd.DataFrame({
                "site": np.repeat(sites["site"].to_numpy()[ok], n_tracks),
                "pos": np.repeat(sites["pos"].to_numpy()[ok], n_tracks),
                "offset": offset,
                "output": output,
                "track": np.tile(np.arange(n_tracks), n_sites),
                "value": block.reshape(-1).astype("float32"),
            })
        )
    return pd.concat(frames, ignore_index=True)


def region_scores(track_data, regions, window_start):
    """Mean coverage of each region per track, from a cumulative sum."""
    res = track_data.resolution
    values = track_data.values.astype("float64")
    csum = np.vstack([np.zeros((1, values.shape[1])), np.cumsum(values, axis=0)])
    start = (regions["start_pos"].to_numpy() - 1 - window_start) // res
    end = -(-(regions["end_pos"].to_numpy() - window_start) // res)
    start = np.clip(start, 0, values.shape[0])
    end = np.clip(end, 0, values.shape[0])
    width = np.maximum(end - start, 1)
    means = (csum[end, :] - csum[start, :]) / width[:, None]
    n_regions, n_tracks = means.shape
    return pd.DataFrame({
        "region": np.repeat(regions["region"].to_numpy(), n_tracks),
        "start_pos": np.repeat(regions["start_pos"].to_numpy(), n_tracks),
        "end_pos": np.repeat(regions["end_pos"].to_numpy(), n_tracks),
        "output": "rna_seq",
        "track": np.tile(np.arange(n_tracks), n_regions),
        "value": means.reshape(-1).astype("float32"),
    })


def junction_scores(junction_data, gene):
    """Every predicted junction on the gene strand inside the gene."""
    junctions = junction_data.junctions
    keep = [
        i
        for i, j in enumerate(junctions)
        if j.strand == gene.strand
        and j.start >= gene.start_pos - 2
        and j.end <= gene.end_pos + 2
    ]
    values = junction_data.values[keep, :]
    n_junctions, n_tracks = values.shape
    return pd.DataFrame({
        "junction_start": np.repeat([junctions[i].start for i in keep], n_tracks),
        "junction_end": np.repeat([junctions[i].end for i in keep], n_tracks),
        "strand": gene.strand,
        "track": np.tile(np.arange(n_tracks), n_junctions),
        "value": values.reshape(-1).astype("float32"),
    }).astype({"junction_start": "int64", "junction_end": "int64"})


def predict(model, interval, condition, retries):
    for attempt in range(retries + 1):
        try:
            return model.predict_interval(
                interval,
                requested_outputs=list(OUTPUTS.values()),
                condition=condition,
            )
        except Exception as error:  # gRPC errors: quota, deadline, unavailable
            if attempt == retries:
                raise
            wait = 15 * 2**attempt
            print(f"  retry in {wait}s: {error}", file=sys.stderr)
            time.sleep(wait)


def main():
    args = parse_args()
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        sys.exit("Set ALPHAGENOME_API_KEY.")
    model = dna_client.create(api_key)

    if args.list_conditions:
        conditions = model.supported_conditions()
        out = args.work_dir / "supported_conditions.parquet"
        conditions.astype(
            {c: str for c in conditions.columns if conditions[c].dtype == object}
        ).to_parquet(out)
        with pd.option_context("display.max_rows", 60, "display.width", 200):
            print(conditions)
        print(f"{len(conditions)} conditions written to {out}")
        return

    targets = args.work_dir / "targets"
    genes = pd.read_parquet(targets / "genes.parquet")
    sites = pd.read_parquet(targets / "sites.parquet")
    regions = pd.read_parquet(targets / "regions.parquet")
    conditions = args.conditions or [None]
    if args.smoke:
        genes = genes.head(1)
    elif args.max_genes:
        genes = genes.head(args.max_genes)

    scores = args.work_dir / "scores"
    for sub in ("sites", "regions", "junctions"):
        (scores / sub).mkdir(parents=True, exist_ok=True)
    tracks = {}
    tracks_file = scores / "tracks.parquet"

    for n, gene in enumerate(genes.itertuples(index=False), start=1):
        done = scores / "junctions" / f"{gene.gene_id}.parquet"
        if done.exists() and tracks_file.exists() and not args.smoke:
            continue
        interval = genome.Interval(gene.seqid, gene.window_start, gene.window_end)
        gene_sites = sites[sites["gene_id"] == gene.gene_id]
        gene_regions = regions[regions["gene_id"] == gene.gene_id]
        site_frames, region_frames, junction_frames = [], [], []
        for condition in conditions:
            label = str(condition) if condition is not None else "none"
            output = predict(model, interval, condition, args.retries)
            if args.smoke:
                describe(output)
            for name in ("splice_sites", "splice_site_usage"):
                data = getattr(output, name)
                tracks.setdefault((name, label), metadata_table(data, name, label))
                frame = site_scores(data, name, gene_sites, gene.window_start)
                site_frames.append(frame.assign(condition=label))
            tracks.setdefault(
                ("rna_seq", label), metadata_table(output.rna_seq, "rna_seq", label)
            )
            region_frames.append(
                region_scores(output.rna_seq, gene_regions, gene.window_start)
                .assign(condition=label)
            )
            tracks.setdefault(
                ("splice_junctions", label),
                metadata_table(output.splice_junctions, "splice_junctions", label),
            )
            junction_frames.append(
                junction_scores(output.splice_junctions, gene).assign(condition=label)
            )
        for sub, frames in (
            ("sites", site_frames),
            ("regions", region_frames),
            ("junctions", junction_frames),
        ):
            # Junctions are written last: their file marks the gene as done
            frame = pd.concat(frames, ignore_index=True).assign(gene_id=gene.gene_id)
            frame.to_parquet(scores / sub / f"{gene.gene_id}.parquet", index=False)
        if not tracks_file.exists():
            pd.concat(tracks.values(), ignore_index=True).to_parquet(
                tracks_file, index=False
            )
        print(f"[{n}/{len(genes)}] {gene.gene_name} {gene.gene_id}", flush=True)


def describe(output):
    """Print the shape and track metadata of each output (smoke mode)."""
    with pd.option_context("display.max_rows", 12, "display.width", 200):
        for name in ("splice_sites", "splice_site_usage", "rna_seq"):
            data = getattr(output, name)
            print(f"\n== {name}: values {data.values.shape}, resolution {data.resolution}")
            print(data.metadata)
        junctions = output.splice_junctions
        print(f"\n== splice_junctions: {len(junctions)} junctions, values {junctions.values.shape}")
        print(junctions.metadata)
        print([str(j) for j in junctions.junctions[:5]])


if __name__ == "__main__":
    main()
