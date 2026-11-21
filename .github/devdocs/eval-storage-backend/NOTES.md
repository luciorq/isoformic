# Storage backend: keep DuckDB or replace it?

Evaluation of 2026-11-18, at HEAD 4a6e983 plus the working tree. Every number
below was measured or checked on that date; where a measurement is missing the
text says so. The benchmark figures already include the corrections of an
independent review (section 5.6).

## 1. Question

Is it worth dropping DuckDB as the storage and query engine of isoformic? The
alternatives considered: a Rust crate (extendr or savvy with arrow-rs
`parquet` or polars), the `arrow` R package, `nanoarrow`, `nanoparquet`,
`polars` (r-polars), `data.table`. Related questions: what to do about the
`duckdb (>= 1.5.4)` pin while CRAN has 2.0.0, and the PLAN.md backlog item
"DuckDB limits" (threads and memory defaults under Slurm, `ulimit -v`).

## 2. How isoformic uses DuckDB

All access goes through `R/utils-duckdb.R`: `duckdb_connect()` (a connection
that closes itself, with the extension, secret and spill directories under
`tempdir()`, and no `threads` or `memory_limit` setting), `duckdb_run()`,
and the exported `read_parquet()`, `write_parquet()` (zstd) and
`read_tabular()`, plus `parquet_nrow()` and `read_delimited()`. Outside these
helpers there are 23 `duckdb::duckdb_register()` calls. No package code opens
a file-backed database.

Call sites in `R/`, by kind:

| Kind | Sites | Notes |
|---|---|---|
| A. Write a whole data frame to Parquet | 32 | experiment IO, merged annotations, frozen data, ORF cache, sequences, splice sites |
| B. Read a whole Parquet file | 16 | `read_isoformic`, `read_assay_dir`, caches, frozen data |
| C. Read with column projection | 14 | 13 of 29 `read_annotation_table()` callers, plus `read_treeterminus` |
| Row-filter (WHERE) pushdown on Parquet | 0 | `read_parquet()` has no filter argument |
| D. SQL on registered R data frames | 10 | range joins, window functions, overlap logic, ordered `string_agg` (reverse complement) |
| E. Text parsing with `read_csv` + SQL | 7 direct, 9 via helpers | GFF3/GTF annotation build (about 317 lines of SQL, data never enters R), FASTA headers, multi-file `quant.sf` with `PIVOT`, Matrix Market, BED/GTF lines, the GTEx GCT |

Outside `R/`: 69 matching lines in 11 test files, 44 `read_parquet(` calls in
the articles, 141 matching lines in `data-raw/`. In `R/` the SQL has 17
`glue_sql()` calls in 11 files, 42 JOINs, 11 window functions, 28 `COPY`
statements.

Data on disk (103 Parquet files surveyed: GENCODE 50 annotation directory,
an experiment directory, `inst/extdata`, `vignettes/articles/data`, the frozen
cache) uses only UTF8 strings, DOUBLE, INT32, BOOLEAN and FLOAT (the latter in
one bundled file). No nested, int64, date or dictionary columns are stored.

Hot paths: the GFF3/GTF annotation build (the only heavy DuckDB job, the only
one that keeps its intermediate data out of R and can spill to disk),
experiment IO (plain reads and writes), `read_salmon` (multi-file read and
PIVOT), and the analysis functions that use the D queries
(`transcript_features()`, `splicing_events()`, polyadenylation, productivity
evidence, diversity). Draws are read in base R; DuckDB is not on that path.

Dependency weight: in the pixi environment duckdb 1.5.5 is 56 MB installed
(55 MB in `libs/`), about 55% of the 101 MB `Imports` closure. A CRAN build of
duckdb 2.0.0 in a separate library was 24 MB; the difference between the two
builds was not investigated. Compile time from source was not measured.

A replacement would have to cover plain Parquet IO (every candidate does),
the 16 text-parsing sites, the out-of-memory annotation build, and the 10
SQL queries on in-memory frames.

## 3. Candidates

| Candidate | CRAN | conda-forge (linux / osx / win-64) | Installed | Build | Licence | Fit |
|---|---|---|---|---|---|---|
| duckdb 2.0.0 | 2.0.0 (2026-10-23); install ERROR on macOS x86_64 release and macOS oldrel | 1.5.5 on linux-64, osx-64, osx-arm64; win-64 stops at 1.1.2; feedstock PRs for 1.5.6, 2.0.0 and "Re-enable Windows builds" open | 24 MB (CRAN build), 56 MB (conda-forge 1.5.5) | vendored C++ engine, C++17 | MIT | covers everything in use |
| nanoparquet 0.5.2 | yes, 13/13 OK | 0.5.2 on all, including win-64 | 8.2 MB | self-contained C/C++, no Imports | MIT | plain read/write and projection; no filter pushdown, no SQL, no CSV, fails on STRUCT columns |
| arrow 26.0.1 | yes, 1 check FAIL (macOS x86_64) | 25.0.1 on all, including win-64 | 25 MB | C++20 plus libarrow, 11 Imports | Apache-2.0 | Parquet, datasets with filter/projection, CSV; SQL-like work only through dplyr/acero |
| nanoarrow 0.9.0 | yes | 0.8.0 on all | 0.9 MB | C | Apache-2.0 | no Parquet at all; not a candidate on its own |
| polars 2.1.0 (r-polars) | **no** (0.7.0 archived 2023 for policy violations); R-multiverse and r-universe only | none | 77 MB | Rust; prebuilt download, source build took 563 s wall in an earlier attempt | MIT | Parquet, lazy scans, CSV, a SQL context |
| data.table 1.18.6.1 | yes | yes, including win-64 | 6.0 MB tarball | C, OpenMP | MPL-2.0 | no Parquet; could only replace in-memory transforms and CSV parsing |
| Custom Rust crate (extendr/savvy + arrow-rs `parquet` 60.1.0) | would ship inside isoformic | possible (other Rust R packages exist on conda-forge) | probe: 5.4 MB `vendor.tar.xz`, 81 crates | rustc >= 1.88; offline build at `-j 2`: 1 min 54 s, 746 MB peak RSS | Apache-2.0 / MIT | all IO and transform code written and maintained in-house; CRAN Rust policy (vendoring, version checks, `-j 2`, author list) |
| fst | last release 2022 | - | - | - | AGPL-3 | not Parquet; ruled out |

DuckDB 2.0.0 compared with 1.5.x: 2.0 is an LTS release (support until at
least October 2027; 1.5 extended to March 2027). The single-arrow lambda
syntax is disabled (not used in `R/`). Database files written by 2.0 are not
readable by older versions (isoformic writes none). macOS x86_64 is being
phased out and 2.1 drops it. `PIVOT` without `IN (...)` warns (already fixed
in `read_salmon`). 2.0.0 prints a storage-directory notice unless the
storage settings are given; `duckdb_connect()` already sets them.

## 4. Resource limits (backlog "DuckDB limits")

Tested with duckdb 1.5.5 and 2.0.0, which behave the same:

| Limit applied | DuckDB `memory_limit` | DuckDB `threads` |
|---|---|---|
| cgroup `MemoryMax=4G` plus `taskset -c 0-3` | 3.1 GiB (80% of the cgroup cap) | 16 (all cores; CPU affinity ignored) |
| cgroup 12G plus `ulimit -v` 4 GB | 9.5 GiB (from the cgroup; ulimit ignored) | not recorded |
| `ulimit -v` 2 GB | R itself fails to start (OpenBLAS) before DuckDB | - |

So the memory default already follows a cgroup limit (which Slurm sets when
it enforces memory with cgroups), but the thread default ignores CPU affinity
and therefore Slurm's CPU allocation. The earlier report "the default
configuration failed to connect under `ulimit -v`" did not reproduce at 4 GB.

## 5. Benchmarks

### 5.1 Method

- Real data: the GENCODE 50 annotation directory (exons 5.23M rows,
  sequences 646k rows with 1.45 GB of text, ORFs and transcripts 670k rows),
  a draws table (295k x 301 doubles), and a 9-file experiment directory
  (1.96M rows). GFF3: the GENCODE 50 `chr_patch_hapl_scaff` file.
- Each measurement is a fresh `Rscript` (R 4.5.3), 3 repeats, median
  reported. Op time measured inside R; op peak is VmHWM after a reset, so it
  includes data already in memory (agrees with `/usr/bin/time -v` within 1%).
  Every read result is touched so lazy (ALTREP) vectors are materialised.
- Threads: duckdb `threads`, arrow `set_cpu_count` (arrow kept at least 2 IO
  threads even at "1 thread"), `POLARS_MAX_THREADS`. nanoparquet is
  single-threaded and ran only at 1 thread, so it is compared at 1 thread
  against the others at 4.
- Each process ran in a memory-capped cgroup (6 GB for IO, 8 GB for the
  GFF3 runs). The shared heavy-job queue was not used for most runs because
  it was held by other jobs for hours; load average was 2 to 3 during the IO
  batch.
- Page cache was neither dropped nor deliberately warmed: inputs (14 to
  109 MB) stay cached after the first read, so all times are warm-cache.
  Process wall times of the first repeat are unreliable and are not reported.
- Versions: duckdb 1.5.5 and 2.0.0, nanoparquet 0.5.2, arrow 26.0.1, polars
  2.1.0, data.table 1.18.6.1.

### 5.2 Read a whole file into R

Seconds / op peak MB. duckdb, arrow, polars at 4 threads; nanoparquet at 1.

| Table | duckdb 1.5.5 | duckdb 2.0.0 | nanoparquet (1 th) | arrow | polars |
|---|---|---|---|---|---|
| exons | 4.48 / 1098 | 5.27 / 1125 | 4.13 / 828 | 3.09 / 1004 | 10.09 / 1364 |
| sequences | 8.16 / 3006 | 9.69 / 3014 | 9.39 / 3243 | 4.96 / 3948 | 8.11 / 3103 |
| orfs | 1.94 / 617 | 2.26 / 631 | 1.97 / 453 | 1.29 / 638 | 4.46 / 719 |
| transcripts | 1.88 / 451 | 2.05 / 465 | 1.73 / 346 | 1.14 / 480 | 3.49 / 594 |
| draws | 1.80 / 1517 | 2.41 / 1517 | 1.77 / 1293 | 0.42 / 1590 | 2.87 / 1526 |
| experiment dir (9 files) | 1.54 / 403 | 1.83 / 423 | 1.34 / 275 | 0.93 / 621 | 2.31 / 498 |

At 1 thread (duckdb 1.5.5 / duckdb 2.0.0 / arrow / polars): exons 5.05 /
5.25 / 3.66 / 10.47 s; draws 2.31 / 2.34 / 0.84 / 3.22 s.

Arrow is the fastest reader, 1.45x (exons) to 4.3x (draws) faster than
duckdb 1.5.5; in absolute terms the gain is at most about 3 s per file.
nanoparquet at 1 thread matches duckdb 1.5.5 at 4 threads with 20 to 30% less
peak memory. polars is the slowest on narrow string tables.

### 5.3 Write zstd

Only DuckDB and nanoparquet numbers are valid as published. In the original
script each library wrote a data frame it had read itself; arrow's reader
returns ALTREP vectors backed by Arrow memory, which arrow then writes without
conversion. isoformic writes data frames computed in R, so that shortcut does
not occur in the package. The polars results may have the same artefact and
were not rechecked.

Seconds / file MB. duckdb at 4 threads (its writes do not speed up with
threads: exons 2.24 s at 1 thread), nanoparquet at 1 thread with
`compression_level = 3`.

| Table | duckdb 1.5.5 | duckdb 2.0.0 | nanoparquet (1 th) |
|---|---|---|---|
| exons | 2.26 / 39.3 | 2.43 / 39.3 | 1.81 / 38.0 |
| sequences | 3.23 / 87.4 | 3.26 / 87.4 | 1.99 / 89.7 |
| orfs | 1.28 / 31.9 | 1.39 / 31.9 | 1.36 / 35.0 |
| transcripts | 0.73 / 13.8 | 0.79 / 13.8 | 0.53 / 13.1 |
| draws | 3.90 / 108.7 | 5.06 / 108.7 | 4.59 / 105.1 |
| experiment dir | 1.20 / 54.4 | 1.41 / 54.4 | 1.19 / 54.4 |
| draws as float32 | 3.68 / 109.7 | 4.85 / 109.7 | 4.56 / 106.2 |

Extra memory over the input while writing: duckdb 1.5.5 needs up to about
1 GB for the draws (+1062 MB) and +396 MB for sequences; nanoparquet adds
almost nothing (+1 MB).

Rerun of arrow on the exons write, 4 threads, same batch (higher load, so
only compare within this table):

| Exons write, 4 threads | Time | Peak MB |
|---|---|---|
| arrow, input it read itself (the original method) | 2.65 s | 892 |
| arrow, plain R vectors read by DuckDB | 5.42 s | 971 |
| arrow, input copied out of ALTREP | 5.30 s | 1307 |
| duckdb 1.5.5 | 2.52 s | 809 |

On realistic input arrow writes about 2x slower than DuckDB. The claim of the
first draft that arrow and polars write fastest does not stand.

### 5.4 Projection plus filter

Seconds / op peak MB at 4 threads (1-thread time in brackets). isoformic has
no row-filter site, so this measures a capability, not current use.
Projection alone (the 14 real sites) was not measured.

| Query | duckdb 1.5.5 | duckdb 2.0.0 | nanoparquet (read, then filter in R) | arrow (dplyr) | polars (lazy) |
|---|---|---|---|---|---|
| transcripts: 2 columns, 2000 genes (14,843 rows) | 0.036 / 119 (0.057) | 0.066 / 137 (0.062) | 0.73 / 161 | 0.10 / 283 (0.16) | 0.13 / 152 (0.16) |
| exons: 2 columns, 5000 transcripts (38,920 rows) | 0.064 / 126 (0.155) | 0.164 / 144 (0.144) | 1.13 / 398 | 0.16 / 399 (0.38) | 0.17 / 158 (0.22) |

### 5.5 Cold start and GFF3 build

Cold start (load the namespace, read a 16-row file; base `Rscript` is 58 MB
and 0.21 s): duckdb 1.5.5 0.36 s / 100 MB; duckdb 2.0.0 0.39 s / 113 MB;
nanoparquet 0.019 s / 61 MB; arrow 0.37 s / 122 MB; polars 0.28 s / 108 MB.

GFF3 to the 4 annotation tables (DuckDB runs execute `build_annotation_sql()`
from HEAD; the alternative is a 160-line data.table + stringi + nanoparquet
rewrite that reads in 2.5M-row chunks after a `zcat | grep -v '^#'` step):

| | 4 threads: s / peak MB | 1 thread: s / peak MB |
|---|---|---|
| duckdb 1.5.5 | 86.1 / 4888 | 330.9 / 4507 |
| duckdb 2.0.0 | 83.3 / 4249 | 318.8 / 3870 |
| data.table rewrite | 243.5 / 7130 (hit the 8 GB cgroup cap in all 3 runs) | 251.1 / 7028 (n = 2) |

- DuckDB is 2.8x faster at 4 threads and peaks 2.2 GB lower (2.5 GB at
  1 thread). The data.table peak includes the 4.4 GB scratch file that
  `fread` maps into memory, so the peaks are not like for like. At 4 threads
  the data.table run was likely slowed by memory reclaim at the cap, so the
  absence of a gain from 4 workers is not conclusive.
- At 1 thread the data.table rewrite is faster (251 s against 331 s).
- On a 25% subset without chunking: DuckDB 22.7 s / 1.44 GB, data.table
  44.1 s / 2.72 GB.
- Both produce the same tables (86,411 genes, 670,670 transcripts, 5,227,964
  exons, 528 sequence regions, row order included). The rewrite covers GFF3
  only, not the GTF path, and its `zcat | grep` step does not exist on
  Windows.
- One 1-thread repeat was disturbed by other jobs (load average 24): a
  duckdb 2.0.0 run took 492 s (outside the median) and a data.table run was
  killed by the machine's memory watchdog (hence n = 2).

### 5.6 Equivalence and corrections

- Files written by every library were read back with duckdb 1.5.5 and match
  the sources with `all.equal(tolerance = 0)`; column names were checked
  separately (82 files, 0 mismatches). The other libraries' own readers were
  checked only on 5 DuckDB-written files (transcripts, orfs, ambiguity,
  abundance, col_data), where they return the same types and values. Exons,
  sequences, draws and float32 files were not read back with nanoparquet,
  arrow or polars.
- Metadata differs, values do not: DuckDB writes plain BYTE_ARRAY strings,
  the others add the STRING logical type; arrow writes one row group where
  DuckDB writes several. polars draws files are 43% larger.
- Float32: every writer gives a largest relative error of 5.95e-08 on
  non-integer TPM, and reads back as double.
- **nanoparquet 0.5.2 trap:** `write_parquet(compression = "zstd")` with the
  default `compression_level = NA` marks pages as ZSTD but leaves them
  uncompressed (sequences 1.46 GB instead of 87 MB). Passing
  `parquet_options(compression_level = 3)` fixes it.
- **duckdb 2.0.0 does not scale reads into R with threads.** Draws read,
  rerun: 1.5.5 3.00 s at 1 thread and 2.42 s at 4; 2.0.0 3.23 s and 3.21 s
  (high load; compare within the rerun only). 2.0.0 is 20 to 30% slower than
  1.5.5 reading wide tables at 4 threads; on the GFF3 build it is 3% faster
  with 13% less peak memory.
- Corrected from the first draft: the GFF3 memory gap (2.24 GB, not 2.7 GB);
  "about as fast for plain IO" (arrow reads are faster); "clearly faster for
  projected reads" (only projection plus filter was measured, by 0.04 to
  0.1 s); the write ranking (section 5.3); the mixing of an untouched-read
  thread series with the touched-read tables.

Not benchmarked: the 10 SQL queries on registered frames, the `read_salmon`
PIVOT, the other text parsers, the GTF build, a Rust crate.

## 6. Decision matrix

Scale: ++ clearly better than today, + better, = same, - worse, -- clearly
worse. "Today" is DuckDB for everything.

| Option | Performance | Footprint | Portability (CRAN, conda-forge, Windows) | Features | Migration cost | Risk |
|---|---|---|---|---|---|---|
| Keep DuckDB (1.5.x now, 2.x later) | = | = (56 MB, 0.36 s start) | = (no win-64 conda build; CRAN binary outside the lock) | = (SQL, CSV, out-of-core) | none | = |
| Partial: nanoparquet for plain IO, DuckDB for SQL | = reads, = to + writes (1 thread), -- filters | - (adds 8 MB, removes nothing) | = | = | medium (48 IO sites, zstd level trap) | - (two Parquet writers, metadata differs) |
| Partial: arrow for plain IO | + reads, - writes on real input | -- (adds 25 MB and 11 Imports) | = | = | medium | - |
| Replace with arrow (+ dplyr) | + reads, - writes, unknown for SQL sites | + (drops 56 MB, adds 25 MB + Imports) | + (win-64 on conda-forge) | - (no SQL; range joins and window queries rewritten; out-of-core GFF3 unproven) | high | - |
| Replace with polars | - reads of string tables, writes unverified | - (77 MB) | -- (not on CRAN, not on conda-forge) | = (SQL context, lazy scans) | high | -- (CRAN/Bioconductor cannot depend on it) |
| Replace with nanoparquet + data.table (+ stringi) | -- GFF3 build (2.8x slower at 4 threads), -- filters | ++ (about 14 MB) | + (all on conda-forge win-64) but the GFF3 rewrite needs `zcat`/`grep` | - (no out-of-core; every SQL site rewritten in R) | very high (about 317 lines of annotation SQL, 10 queries, 16 parsers) | -- (behaviour changes in core annotation) |
| Custom Rust crate | unmeasured | + to - (small vendor tarball, Rust toolchain needed) | - (Rust toolchain on CRAN, conda-forge and Windows) | whatever is written | very high (all IO and transforms in-house) | -- (maintenance, CRAN Rust policy) |
| nanoarrow | - | - | - | -- (no Parquet) | - | - |

## 7. Recommendation

**Keep DuckDB as the single storage and query engine. Do not replace it, and
do not add a second Parquet library.**

- The cost that would be removed is small: about 56 MB installed and
  0.35 s at startup. The measured IO gains of the alternatives are seconds per
  call (arrow reads) and disappear or reverse for writes on real input.
- The work that would be added is large: about 317 lines of annotation SQL
  that runs out of R memory, 10 range/window queries on registered frames, and
  16 text-parsing sites. The only full replacement measured (GFF3 in
  data.table) is slower at 4 threads, needs more memory and does not run on
  Windows as written.
- A partial split (nanoparquet or arrow for plain IO) adds a dependency
  without removing DuckDB, introduces a second writer with different Parquet
  metadata, and in nanoparquet's case a silent no-compression default.
- polars is not on CRAN or conda-forge; nanoarrow cannot read Parquet; a Rust
  crate would make isoformic maintain its own IO layer and a Rust build on
  every platform.

Revisit only if one of these happens: DuckDB loses CRAN or conda-forge support
on a platform isoformic targets; a workload shows plain Parquet reads as a
real bottleneck (then consider arrow for that path alone, after rerunning
the write and round-trip checks of section 5); or a dependency-light install
becomes a requirement. In that case the order would be: plain IO first
(`read_parquet()`/`write_parquet()` behind the existing helpers, so no call
site changes), then the text parsers, and the annotation SQL last.

### 7.1 The duckdb 2.x pin

- Keep `duckdb (>= 1.5.4)` in `DESCRIPTION`. conda-forge still ships 1.5.5
  (the 2.0.0 feedstock PR is open), and a higher minimum makes the pixi
  environment unsolvable.
- There is no performance reason to require 2.0: it is slightly better on the
  GFF3 build and 20 to 30% slower reading wide tables into R. Its CRAN build
  also fails on macOS x86_64.
- isoformic already works with 2.0.0 (no single-arrow lambdas, `PIVOT ...
  IN (...)` in place, storage settings set by `duckdb_connect()`, same GFF3
  output). Keep testing against the CRAN 2.x release where CRAN packages are
  used (Windows, R CMD check outside pixi).
- Raise the minimum when conda-forge publishes r-duckdb 2.x, together with a
  pixi lock update. Watch the 2.x read-into-R thread regression and report it
  upstream if it persists.

### 7.2 Threads and memory limits (backlog item)

Implement it; it is the useful DuckDB change.

- Add `isoformic.duckdb_threads` and `isoformic.duckdb_memory_limit`, passed
  in the `duckdb_connect()` config.
- Threads are the real gap: DuckDB uses all cores even under `taskset` or a
  Slurm CPU allocation. Default: the option if set, else
  `SLURM_CPUS_PER_TASK` if set, else DuckDB's own default.
- Memory: DuckDB already takes 80% of a cgroup memory cap, so a Slurm job with
  cgroup enforcement is covered. Default the option to DuckDB's value and let
  users set it where only `ulimit -v` or no cgroup applies; optionally derive
  it from `SLURM_MEM_PER_NODE`.
- The reported connection failure under `ulimit -v` did not reproduce at 4 GB;
  update the backlog text accordingly.
- Document both options and the Slurm variables in `R/isoformic-options.R`
  (and the option count), and add a test that the settings reach
  `current_setting()`.
