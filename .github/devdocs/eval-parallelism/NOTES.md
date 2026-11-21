# Shared-memory parallelism: mirai, mizu and mori

Evaluation for the PLAN.md backlog item "Evaluate shared-memory parallelism
with `mirai`, `mizu` and `mori`, and whether to adopt one" (2026-11-18).

## Question

Today the `fgsea` calls of `run_enrichment()`, `run_productivity_enrichment()`
and `run_apa_enrichment()` run serially: `fgsea_multilevel()` forces
`BiocParallel::SerialParam()` unless the user passes `nproc` or `BPPARAM`.
Those two fork the session (`MulticoreParam`) on Unix. Forked workers share
the session's pages, but tools that sum resident memory (Slurm accounting
without `UsePss`, `top`) count the session once per worker, the likely
source of a reported 160 GB peak. Questions:

1. Is there a faster way to run the many fgsea runs (contrasts x gene set
   collections x transcript types) from one `productive_fraction()` table?
2. What does each option cost in physical memory (PSS, cgroup peak) and in
   summed RSS, against serial and `MulticoreParam`?
3. Do worker counts follow the scheduler's allocation?
4. What does each package add to `Imports`, `Suggests` and `pixi.toml`?

## Status of the evidence

Read this before the numbers.

- The full benchmark ran on 2026-11-20: every strategy of *Method* on both
  workloads with N = 4, the package's own `workers` path in both call
  shapes, N = 2 and N = 8 cells for W1 and an N = 2 `mc` cell for W2; 30
  cells x 3 repeats = 90 runs in shuffled order, one process per run, the
  whole batch under the shared heavy-job lock. The 1-minute load average
  was 1.2 to 2.5 on 16 cores before and after every run. Runs were not
  pinned to CPUs (the *Implementation check* batches below were, with
  `taskset -c 0-3`), so N = 8 could use 8 cores; compare ratios within a
  batch, not times across batches.
- Code: the working tree of 2026-11-20 (kept daemons with the health check
  below), installed into a separate library; the `productive_fraction()`
  tables and both serial references were rebuilt with it. R 4.5.3, mirai
  2.7.4, mizu 0.1.1, mori 0.2.2, fgsea 1.36.2, BiocParallel 1.44.0.
- All 90 runs exited normally and passed their result checks (below). The
  3 repeats agree closely: wall-time ranges are within 0.6 s (most within
  0.1 s), memory ranges within 0.03 GB in W1 (one exception, marked) and
  within 0.6 GB in W2. Tables give the median and the range of 3.
- The single smoke runs of 2026-11-18 that these tables replace agree with
  the new medians within about 10% wherever both exist.
- Package semantics (copying, sharing, RNG, cleanup, errors) were verified
  with small probes; facts from documentation or source only are marked
  `(docs)`.

## Packages compared

All three are by the same author (Posit). Versions: mirai 2.7.4 (nanonext
1.10.4), mizu 0.1.1, mori 0.2.2; parallelly 1.50.0 for core detection.

| | mirai 2.7.4 | mizu 0.1.1 | mori 0.2.2 |
|---|---|---|---|
| What it is | Async evaluation on `Rscript` daemons over NNG sockets; `mirai_map()`; `make_cluster()` | Shared-memory channels, work-stealing pools, parallel map; R workers (`Rscript`) or Python workers (`pymizu`) | Puts R objects in shared memory (ALTREP); starts no workers |
| Data to workers | Serialized and copied into each daemon (verified: a 57 MB frame cost the daemon its full size in private memory). Zero-copy when wrapped in `mori::share()` (verified) | `mizu_submit()` stages arguments in shared memory: a plain frame is mapped, not copied (verified), but staged again on every submit. A list `.x` in `mizu_map()` gave the worker a private copy (verified) | One copy per `share()`. Shares numeric, integer, character, factor, matrix, data.frame, nested lists (verified). Not environments, closures, language objects or S7 objects (returned unchanged, verified) |
| S7 objects | `identical()` FALSE after the trip, `all.equal()` TRUE (verified): the `S7_class` attribute is a closure | Same (verified) | Not shared (verified) |
| Library paths | Daemons do NOT get the session's `.libPaths()`; without `R_LIBS` they hung at start-up (verified) | Workers get the session's `.libPaths()` (verified) | n/a |
| Worker count | `n` must be given; no detection (docs) | `n_workers = 1` by default, max 64; no detection (docs) | n/a |
| RNG | `daemons(seed =)`: one reproducible stream per mirai (verified); unseeded L'Ecuyer streams per daemon are not reproducible (verified) | `mizu_map(.seed =)`: one stream per element, same for any worker count (verified); `mizu_submit()` has no seed (verified) | n/a |
| Errors | A `miraiError` value (returned, not raised) with message, call, condition class, stack trace (verified) | Original condition re-raised on collect, but the call shows as `doTryCatch(...)` (verified); classed `mizu_error_*` conditions; map errors carry the index | n/a |
| Interrupt / host killed | Ctrl-C stops the wait, the task resolves later. Daemons exited within 5 s of the host dying by KILL, TERM or INT, even mid-task (verified) | Ctrl-C stops the wait, the pool survives (verified). Host killed: idle workers exit and shm is removed, a busy worker finishes its current task first (verified) | Region freed on GC or clean exit; left in `/dev/shm` after SIGKILL and SIGTERM, removed by `prune_shared()` (verified) |
| fgsea / BiocParallel inside | `bpparam()` in a daemon is `MulticoreParam(14)`; a `miraiCluster` cannot be a `BPPARAM` (verified) | `fgseaMultilevel()` runs in a worker with `SerialParam` and nested `MulticoreParam(2)` (verified, tiny data); `bpparam()` there is also 14 workers | n/a |
| Dependencies | Imports nanonext (compiled, bundles NNG and mbedtls); R >= 3.6 | No Imports; compiled C with pthread; R >= 4.3; needs `/dev/shm` large enough (64 MB default in containers) | No Imports; compiled C; R >= 4.3 |
| Platforms | All; CRAN checks OK | Linux (kernel >= 5.3), macOS, Windows; CRAN OK on 13 flavours (docs) | Same as mizu; CRAN OK (docs) |
| Maturity | On CRAN since 2022, 52 releases | First CRAN release 2026-09-30; API marked experimental | On CRAN since 2026-04, 4 releases |
| conda-forge (checked 2026-11-18) | `r-mirai` 2.7.4 on linux-64, osx-arm64, win-64. `r-nanonext` 1.10.4 on linux-64 and osx-arm64, but only 1.4.0 on win-64 (too old for mirai 2.7.4, which needs >= 1.10.4) | Not available | Not available |

`parallelly` 1.50.0 (Imports only base packages) is on conda-forge for all
three platforms.

## condathis findings

condathis documents mizu with Python workers (`pymizu` from PyPI, workers
started by `mizu_py_launcher()` from a micromamba environment) and keeps a
mirai + mori + reticulate section for one Python interpreter per daemon.

Findings that transfer to isoformic:

- mizu and mirai + mori move large R data to workers without a copy; mizu
  even does so for unshared frames passed to `mizu_submit()` (confirmed here
  for R workers).
- Workers do not see the session's activation by default: condathis had to
  pass the activated `PATH` to Python workers on Windows. The R analogue
  here is that mirai daemons do not inherit `.libPaths()`.
- Packages missing from conda-forge need special handling in pixi
  (`pymizu` there, `r-mizu` and `r-mori` here).

Findings that do not transfer:

- Python workers, `pymizu`, `mizu_py_launcher()`, the Arrow serializer and the
  Windows DLL and `PATH` issues: isoformic needs R workers only, which mizu
  starts with its default `mizu_launcher()`.
- The motivation in condathis (no `reticulate` in workers, one interpreter
  per environment) has no counterpart: isoformic's tasks are plain R calls to
  `fgsea`.
- Type limits of Python workers (no functions, formulas, S4) do not apply to
  R workers; the S7 limit above does.

## Worker counts and the scheduler

Measured on the 16-core node:

| Setting | `parallelly::availableCores()` | `BiocParallel::multicoreWorkers()` | DuckDB `threads` |
|---|---|---|---|
| none | 16 | 14 | 16 |
| `SLURM_CPUS_PER_TASK=4` | 4 | 14 | 16 |
| `taskset -c 0-3` | 4 | 14 | 16 |
| cgroup `cpu.max` 400% | 4 | 14 | 4 |
| `MC_CORES=5` and `BIOCPARALLEL_WORKER_NUMBER=3` | 3 | 5 | - |
| `BIOCPARALLEL_WORKER_NUMBER=3` alone | 3 | 3 | - |

Only `parallelly::availableCores()` followed the allocation in every case.
BiocParallel's defaults ignore Slurm, CPU affinity and `cpu.max`; DuckDB
follows only `cpu.max`. mirai and mizu do no detection at all.

## Method

- Data: a GENCODE 50 experiment of 295k transcripts x 12 samples (4 groups
  of 3) with 25 inferential draws (855 MB object, held in the session as a
  user session would hold it); 3 synthetic gene set collections of 500 sets.
- W1 (the backlog's starting point): 18 fgsea runs, 6 contrasts (all pairs
  of 4 groups) x 3 collections, `run_productivity_enrichment(rank_by =
  "delta_z")` on precomputed `productive_fraction()` tables. Constants sent
  to workers: 39.4 MB.
- W2 (shared-memory stress test): the per-sample loop of
  `inferential_correlation()`, 12 tasks over the 727 MB draws view. Constants:
  762.8 MB.
- Strategies: serial; `fgsea_nproc` (fgsea's own `nproc`, today's user
  option); `mc` (`BiocParallel::MulticoreParam`); `snow` (`SnowParam`);
  `mirai` (constants sent once per daemon with `everywhere()`); `mirai_mori`;
  `mizu` (`mizu_map()`); `mizu_mori`. N = 4 workers unless stated.
- The package's own path (added 2026-11-20, after `R/parallel.R` was
  implemented): `iso_loop` and `iso_one` run the same 18 fgsea runs through
  `run_productivity_enrichment(workers = N)` (kept daemons of the
  `"isoformic"` profile, per-task seeds), in the two shapes of
  *Implementation check*: `iso_loop` is 6 calls (one per contrast) of 3
  collections, so at most 3 workers per call; `iso_one` is 1 call on the
  first contrast with 18 collections (the 3, six copies each). `iso_tasks`
  sends the W2 loop through the internal `run_tasks()`, as
  `inferential_correlation()` would if it took `workers`: the constants
  travel with every task. Daemon start-up is inside the wall time.
- One process per run, each in its own cgroup scope with a hard memory cap.
  An external sampler read `cgroup.procs` of that scope every 0.2 s (so
  daemons are counted even after reparenting) and summed `Pss` and `Rss` from
  `/proc/<pid>/smaps_rollup`; `memory.peak` of the scope was read at the end
  (cgroup v2; it also covers the `pixi` and `bash` wrappers and page cache).
  Load average recorded before and after each run.
- Each run compared its result with a serial reference using `identical()`.
  Every W1 task calls `set.seed()` itself and W2 results are folded in sample
  order in the session, so `identical()` checks the harness; it says nothing
  about the reproducibility of each backend's own RNG streams.

## Benchmark results

Wall time is the in-process time of the parallel section, including pool
or daemon start-up. Memory columns are peaks over the whole job. Speed-up
is serial median / cell median. Load = range of the 1-minute load average
before and after the 3 runs. All values: median of 3 (range).

### W1: 18 fgsea runs (constants 39.4 MB)

Memory ranges were within 0.03 GB except for `fgsea_nproc` (given).

| Strategy | N | Wall s | Speed-up | PSS GB | cgroup peak GB | Summed RSS GB | Load |
|---|---|---|---|---|---|---|---|
| serial | 1 | 13.51 (13.46-13.56) | 1.00 | 0.95 | 0.91 | 0.96 | 1.2-2.2 |
| `iso_loop` | 1 | 13.20 (13.09-13.38) | 1.02 | 0.95 | 0.90 | 0.95 | 1.5-1.8 |
| `iso_one` | 1 | 13.81 (13.79-13.83) | 0.98 | 0.95 | 0.90 | 0.95 | 1.5-1.9 |
| `mizu_mori` | 2 | 8.17 (8.06-8.18) | 1.65 | 1.20 | 1.15 | 1.28 | 1.7-2.0 |
| `mc` | 2 | 8.23 (8.23-8.29) | 1.64 | 1.60 | 1.55 | 2.61 | 1.5-2.0 |
| `mirai` | 2 | 8.59 (8.54-8.83) | 1.57 | 1.24 | 1.20 | 1.28 | 1.8-2.2 |
| `iso_one` | 2 | 9.99 (9.92-10.11) | 1.35 | 1.63 | 1.59 | 1.67 | 1.7-2.0 |
| `iso_loop` | 2 | 11.69 (11.66-11.74) | 1.16 | 1.63 | 1.59 | 1.67 | 1.6-2.0 |
| `mizu_mori` | 4 | 5.41 (5.40-5.42) | 2.50 | 1.70 | 1.66 | 1.86 | 1.4-2.0 |
| `mizu` | 4 | 5.49 (5.42-5.49) | 2.46 | 1.81 | 1.77 | 1.98 | 1.6-2.0 |
| `mc` | 4 | 5.71 (5.69-5.75) | 2.37 | 2.41 | 2.41 | 4.43 | 1.7-2.0 |
| `mirai_mori` | 4 | 5.77 (5.73-5.81) | 2.34 | 1.80 | 1.77 | 1.93 | 1.8-2.0 |
| `mirai` | 4 | 5.97 (5.96-6.03) | 2.26 | 1.79 | 1.75 | 1.87 | 1.7-2.1 |
| `snow` | 4 | 6.76 (6.75-6.81) | 2.00 | 2.04 | 2.07 | 2.13 | 2.0-2.2 |
| `iso_one` | 4 | 7.28 (7.28-7.30) | 1.86 | 2.22 | 2.18 | 2.30 | 1.9-2.2 |
| `iso_loop` | 4 (3 used) | 8.31 (8.27-8.40) | 1.63 | 1.93 | 1.89 | 1.99 | 1.8 |
| `fgsea_nproc` | 4 | 20.48 (20.31-20.90) | 0.66 | 0.92 (0.92-1.21) | 0.94 (0.94-1.21) | 1.76 | 1.6-2.2 |
| `mizu_mori` | 8 | 4.16 (4.16-4.18) | 3.25 | 2.70 | 2.66 | 3.00 | 1.5-2.3 |
| `mc` | 8 | 4.42 (4.40-4.45) | 3.05 | 3.29 | 3.28 | 6.32 | 1.4-2.5 |
| `mirai` | 8 | 4.92 (4.90-4.93) | 2.75 | 2.90 | 2.87 | 3.06 | 1.5-2.4 |
| `iso_one` | 8 | 6.00 (5.97-6.04) | 2.25 | 3.38 | 3.35 | 3.54 | 1.4-2.3 |

Checks: every prototype run was `identical()` to the serial reference.
For `iso_loop` and `iso_one`, `pathway`, `ES`, `size` and `leadingEdge`
were `identical()` to the serial reference in every run (`iso_one`: each
of the 6 copies against the first contrast); `workers = 1` was
`identical()` across repeats, and every run with 2, 4 or 8 workers was
`identical()` to every other. Against `workers = 1`, the sampled columns
differ as documented in *Implementation status*: largest absolute p-value
difference 0.105 (`iso_loop`) and 0.106 (`iso_one`), correlation of log10
p-values 0.992 and 0.993, no gene set with `padj < 0.05` in either.

Pool start-up inside the wall time (median): `mirai` 0.47, 0.55 and 0.74 s
for 2, 4 and 8 daemons (constants sent with `everywhere()`), `mirai_mori`
0.33 s, `mizu` and `mizu_mori` 0.22 to 0.29 s, `snow` 1.24 s.

### W2: per-sample correlation loop (constants 762.8 MB)

| Strategy | N | Wall s | Speed-up | PSS GB | cgroup peak GB | Summed RSS GB | Load |
|---|---|---|---|---|---|---|---|
| serial | 1 | 8.45 (8.42-8.56) | 1.00 | 1.84 (1.84-1.84) | 1.78 (1.78-1.78) | 1.85 (1.84-1.85) | 1.4-2.3 |
| `mc` | 2 | 5.67 (5.66-5.79) | 1.49 | 3.00 (2.91-3.03) | 3.09 (3.06-3.10) | 4.30 (4.21-4.33) | 1.4-2.1 |
| `mc` | 4 | 3.42 (3.40-3.48) | 2.47 | 3.52 (3.49-3.53) | 3.62 (3.50-3.63) | 6.13 (6.10-6.13) | 1.5-1.9 |
| `mizu` | 4 | 5.38 (5.31-5.42) | 1.57 | 5.56 (5.50-5.75) | 5.70 (5.69-5.70) | 7.63 (7.57-7.82) | 1.5-1.9 |
| `mizu_mori` | 4 | 6.56 (6.42-6.72) | 1.29 | 5.43 (5.42-5.67) | 5.56 (5.56-5.57) | 7.50 (7.49-7.69) | 1.4-1.8 |
| `mirai` | 4 | 8.75 (8.68-8.86) | 0.97 | 5.75 (5.58-5.88) | 5.82 (5.82-5.85) | 5.82 (5.65-5.96) | 1.9-2.0 |
| `mirai_mori` | 4 | 10.47 (10.35-10.50) | 0.81 | 5.69 (5.68-5.77) | 5.69 (5.68-5.70) | 9.19 (9.14-9.40) | 1.4-2.2 |
| `iso_tasks` | 4 | 12.94 (12.85-13.08) | 0.65 | 9.04 (8.72-9.30) | 9.28 (8.98-9.29) | 9.11 (8.79-9.38) | 1.6-2.4 |
| `snow` | 4 | 17.00 (16.96-17.15) | 0.50 | 10.17 (10.10-10.18) | 10.13 (10.13-10.13) | 10.25 (10.18-10.25) | 1.2-2.4 |

Every run was `identical()` to the serial reference. Pool start-up inside
the wall time (median): `mirai` 5.71 s (sending the 763 MB of constants to
4 daemons with `everywhere()`), `mirai_mori` and `mizu_mori` 1.04 to 1.08 s
(`mori::share()` copies the constants into shared memory), `snow` 1.40 s,
`mizu` 0.22 s. The session's own peak (VmHWM) was 1.78 GB serial and
4.31 GB under `iso_tasks`. No `mizu` or `mori` region was left in
`/dev/shm` and no daemon outlived its run.

A probe after the batch (one `mirai` and one `mizu_mori` run, each daemon
reading its own `/proc/self/smaps_rollup` at the end of every task) found
0.92 to 0.93 GB of private anonymous memory per `mizu_mori` worker against
0.97 GB per `mirai` daemon; a fresh R process with isoformic loaded holds
0.08 GB. Subsetting a matrix shared with `mori::share()` (`m[ord, s]`,
also inside a shared list) added the whole matrix (100 MB) to the
process's private anonymous memory (verified in the sharing process, mori
0.2.2, R 4.5.3).

### What the results show

- W1, N = 4: every out-of-process backend and `MulticoreParam` gave 2.0x
  to 2.5x over serial. The order is now established (ranges do not
  overlap): `mizu_mori` 5.41 s, `mizu` 5.49, `mc` 5.71, `mirai_mori` 5.77,
  `mirai` 5.97, `snow` 6.76. mizu was 8% (`mizu`) to 9% (`mizu_mori`)
  faster than mirai with 4 workers, and `mizu_mori` 15% faster with 8
  (4.16 vs 4.92 s); in absolute terms the spread is under 0.8 s on 13.5 s
  of serial work.
- Scaling over 2, 4 and 8 workers: `mizu_mori` 1.65x, 2.50x, 3.25x;
  `mc` 1.64x, 2.37x, 3.05x; `mirai` 1.57x, 2.26x, 2.75x. Eighteen tasks of
  unequal size leave little to gain past 8 workers.
- fgsea's own `nproc = 4` is confirmed slower than serial: 0.66x (20.48 vs
  13.51 s). Parallelism belongs across fgsea runs, not inside one.
- W1 memory: each mirai or mizu worker added 0.13 to 0.24 GB of PSS, each
  forked `mc` worker 0.29 to 0.37 GB. Summed RSS overcounts forks only:
  1.8x to 1.9x PSS under `mc` (4.43 vs 2.41 GB with 4, 6.32 vs 3.29 GB with
  8), within 11% of PSS for every other strategy. PSS and the cgroup peak
  agreed within 0.05 GB in every cell. The (N + 1) x session-RSS
  arithmetic for the reported 160 GB peak still holds (0.96 GB x 5 = 4.8
  GB vs 4.43 GB measured).
- mori adds nothing to W1: the constants are 39 MB (`mirai_mori` 5.77 s and
  1.80 GB vs `mirai` 5.97 s and 1.79 GB).
- The package's `workers` path: one call with 18 tasks (`iso_one`) gained
  1.35x, 1.86x and 2.25x with 2, 4 and 8 workers; the loop of 6 small calls
  (`iso_loop`) 1.16x with 2 and 1.63x with 4 (3 used), so the kept daemons
  keep small calls faster than serial. `workers = 1` runs at serial speed
  (13.20 and 13.81 vs 13.51 s). With 4 workers it is 1.3 s (22%) slower
  than the prototype `mirai` strategy (7.28 vs 5.97 s) and uses 0.43 GB
  more PSS (2.22 vs 1.79 GB; 0.32 vs 0.21 GB per daemon). Each call waits
  for its daemons to connect, every task carries its own ranks, collection
  and seed and collects its warnings, and `iso_one` runs 18 collections on
  one contrast instead of 6 contrasts x 3 collections; the gap was not
  split between these causes.
- W2 is the opposite case. Forking won on both time and physical memory:
  `mc` gave 2.47x with 4 workers for 0.42 GB of PSS per worker, since the
  children read the session's 763 MB of draws without copying them; only
  summed RSS (6.13 GB, 1.7x PSS) is inflated. Every out-of-process backend
  paid about 0.9 to 1.0 GB of PSS per worker (5.4 to 5.8 GB with 4), one
  copy of the inputs each.
- Shared memory saved nothing in W2: `mizu_mori` used 5.43 GB of PSS
  against 5.75 GB for `mirai`. The probe explains it: `draws_sample()`
  reads each draw matrix with `[`, and subsetting a `mori`-shared matrix
  materializes a private copy of it, so every worker ends up holding the
  draws again. The zero-copy shown in *Packages compared* holds for
  passing data, not for this access pattern. Summed RSS made it look worse
  (`mirai_mori` 9.19 GB summed RSS for 5.69 GB PSS: the shared region is
  counted in every process).
- Out-of-process W2 time is dominated by moving the inputs: `mirai` spent
  5.71 of 8.75 s sending them to 4 daemons (no gain over serial),
  `mirai_mori` and `mizu_mori` about 1 s sharing them, and `snow` (0.50x,
  10.2 GB) sent them with the tasks.
- `run_tasks()` is not built for large constants: on W2 (`iso_tasks`) it
  was 0.65x of serial with 9.04 GB of PSS (4.9x serial), and the session's
  own peak rose from 1.78 to 4.31 GB. Loops over the draws should not take
  `workers` through this helper.

## Recommendation

1. Adopt `mirai` for the outer loop over fgsea runs (contrasts x collections
   x transcript types) in `run_enrichment()`, `run_productivity_enrichment()`
   and `run_apa_enrichment()`. Each task receives only its ranks vector and
   one collection, never the experiment object or the session, so neither
   physical memory nor summed RSS grows with the session. The full
   benchmark confirms this choice: within about 10% (N = 4) of the fastest
   backend on W1, with summed RSS within 6% of PSS. Implemented; its path
   is 22% slower than the prototype with 4 workers (*What the results
   show*), worth a profile of `run_tasks()` before more functions take
   `workers`.
2. Put `mirai` in `Suggests`, not `Imports`, used through
   `rlang::check_installed("mirai")` only when more than one worker is asked
   for.
3. Control it with an argument `workers` on those functions (the name
   `test_usage()` already uses), defaulting to a new package option
   `isoformic.workers`, whose default is `1L`: serial, as today. No
   parallelism by default, so CRAN checks, Windows and small inputs stay as
   they are.
4. Worker counts: the user gives a number; `workers = "auto"` resolves to
   `parallelly::availableCores()` (follows Slurm, `taskset`, `cpu.max`).
   Never use `parallel::detectCores()` or BiocParallel's defaults. Cap the
   count at the number of tasks. `parallelly` goes to `Suggests`; without
   it, `"auto"` falls back to 1 with a message.
5. Do not adopt mizu or mori now. Neither is on conda-forge, mizu's API is
   experimental and about seven weeks on CRAN, and the measured gains do
   not pay for that: on W1 mizu was 8% (N = 4) to 15% (N = 8) faster than
   mirai, under a second on these 18 runs, and mori saved no memory there
   (39 MB of inputs). On W2, the large-input case where shared memory
   should help, it saved none either: subsetting a mori-shared matrix makes
   a private copy, so each worker held the draws again (5.43 GB PSS for
   `mizu_mori` vs 5.75 GB for `mirai`). Revisit only if mori keeps subsets
   shared and both packages reach conda-forge.
6. Keep `inferential_correlation()` (W2) and other loops over the draws
   serial, and do not route them through `run_tasks()`: there it was 0.65x
   of serial with 4.9x the memory, and every out-of-process backend copied
   the inputs into each worker (0.9 to 1.0 GB each) for at most 1.57x. The
   only option that paid off on W2 was forking (`mc`: 2.47x with 4 workers
   for 0.42 GB of PSS per worker, Unix only, summed RSS 1.7x PSS). If a
   draws loop ever needs parallelism, that is a separate decision with the
   accounting caveat of *Question*, not an extension of `workers`.
7. Keep `nproc` and `BPPARAM` pass-through for compatibility, but stop
   recommending `nproc` in the documentation (0.66x of serial with
   `nproc = 4`, 3 repeats) and say that `MulticoreParam` forks the session.
8. `pixi.toml`: add `r-mirai` and `r-parallelly` to
   `[target.unix.dependencies]` (next to `bioconductor-biocparallel`). On
   Windows, conda-forge's `r-nanonext` is 1.4.0, too old for mirai 2.7.4:
   do not add `r-mirai` under `[target.win-64.dependencies]`; let
   `pixi run r-packages` (`.github/pixi/windows-packages.R`, which installs
   the `Suggests` from CRAN binaries) install `mirai` and `nanonext`.
   `r-parallelly` 1.50.0 is on conda-forge for win-64 and can go there.

## Risks

- Speed-up is modest at small task counts: daemon start-up (0.3 to 0.6 s,
  plus loading isoformic and fgsea in each daemon) is paid per call. With
  one contrast and one collection there is nothing to parallelise.
- Library paths: mirai daemons do not inherit `.libPaths()` and hung without
  `R_LIBS` in the probe. The implementation must pass the session's library
  paths to the daemons, or daemons will fail for users with a non-default
  library (renv, pak, `R_LIBS_USER` set in the session).
- Oversubscription: inside a daemon `bpparam()` is `MulticoreParam(14)`.
  `fgsea_multilevel()` already forces `SerialParam`, and the daemons must
  keep it that way; BLAS and OpenMP threads should be limited to 1 per
  daemon.
- RNG: fgsea results depend on the RNG. Unless seeds are fixed per task,
  results change with the worker count. Fixing them changes the serial
  results relative to today's (the session stream no longer runs through
  the calls in order), which may move existing test expectations.
- Errors come back as `miraiError` values, not conditions: they must be
  checked and re-raised with a message that names the failing task.
- User-owned daemons: a user who already set `mirai::daemons()` must not
  have them reset; use a dedicated compute profile.
- DuckDB and BiocParallel defaults still ignore Slurm and affinity
  (16 threads / 14 workers on a 4-core allocation). Out of scope here, but
  the same `workers` value could cap DuckDB `threads` later.
- Unmeasured: macOS and Windows behaviour, Slurm (`/dev/shm` per job, cgroup
  accounting of shared pages), and machines with fewer than 16 cores.

## Implementation plan

1. (Done 2026-11-20, see *Benchmark results*.) Rerun the benchmark plan
   (W1 and W2, all strategies, N = 4, 3 repeats) under the heavy-job lock
   on a quiet machine, add N = 2 and N = 8 cells
   for W1, and replace the tables above with medians. Do this before the
   code change is merged.
2. New internal helper in `R/parallel.R`, for example
   `run_tasks(tasks, fun, workers)`:
   - `workers <= 1` or one task: `lapply()`, as today.
   - Otherwise: `rlang::check_installed("mirai")`; start
     `mirai::daemons(n, .compute = "isoformic")` unless that profile is
     already set by the user, and stop it with `on.exit()` only if started
     here; send the session's `.libPaths()` to the daemons before the first
     task; set `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS` to 1 in them.
   - Draw one integer seed per task in the session before dispatch
     (`sample.int(.Machine$integer.max, length(tasks))`); each task calls
     `set.seed()` with its own seed. Results are then `identical()` for any
     worker count, serial included.
   - Run with `mirai::mirai_map(..., .compute = "isoformic")[]`; on any
     `miraiError`, abort with `cli::cli_abort()` naming the task (contrast,
     collection, transcript type) and the original message.
3. `resolve_workers(workers)`: integer >= 1, or `"auto"` via
   `parallelly::availableCores()` when installed (else 1 with a
   `cli::cli_inform()`); capped at the number of tasks.
4. Flatten the loops of `enrich_collections()` and `run_fgsea_by_type()`
   into one task list (contrast x collection x transcript type), each task
   holding a named numeric ranks vector and one collection; keep
   `fgsea_multilevel()` with `SerialParam` inside the task.
5. Add `workers = getOption("isoformic.workers", 1L)` to
   `run_enrichment()`, `run_productivity_enrichment()`,
   `run_apa_enrichment()`; document it in the generics' topics (S7
   signatures explicit).
6. Update `R/isoformic-options.R`: the new option `isoformic.workers` (the
   count of package options goes from two to three), the BiocParallel
   entry, and the environment variables that `parallelly::availableCores()`
   reads (for example `SLURM_CPUS_PER_TASK`, `MC_CORES`), listing them as
   settings of a dependency.
7. DESCRIPTION: `mirai` and `parallelly` in `Suggests`. pixi.toml as in the
   recommendation; check `pixi run r-packages` installs mirai on Windows.
8. Tests (default mode, bundled data): `workers = 2` gives a result
   `identical()` to `workers = 1`; an error inside a task is re-raised with
   the task named; `"auto"` respects `SLURM_CPUS_PER_TASK` (via
   `withr::with_envvar()`); all with `testthat::skip_if_not_installed()` and
   `testthat::skip_on_cran()`.
9. Articles: replace the `nproc = 4` advice with `workers`, and explain
   summed RSS vs PSS for forked workers.

## Implementation status (2026-11-18)

Implemented in `R/parallel.R` (`resolve_workers()`, `run_tasks()`,
`start_workers()`) and the three enrichment functions, before the rerun of
step 1 (it was still waiting for the heavy-job lock). Differences from the
plan above:

- Seeds: with one worker (the default) the tasks run in the session with
  its random stream, in the order of the old loops, so results are
  `identical()` to earlier versions (checked against the code of HEAD
  4a6e983, including the session's RNG state afterwards). Per-task seeds
  apply only when more than one worker runs: the result is then
  `identical()` for any worker count above one, but its sampled p-values
  differ from a serial run (`ES`, `size`, `leadingEdge` do not). A single
  task always runs in the session.
- The task function is sent with the base environment and no srcref, so
  daemons load only `fgsea` and `BiocParallel`, never isoformic (works under
  `devtools::load_all()`).
- `mirai::daemons()` and `mirai_map()` draw from the session's stream; the
  calls are wrapped in `withr::with_preserve_seed()` so the session's stream
  moves by the drawn seeds only.
- Daemons get `R_LIBS` (the session's `.libPaths()`) and
  `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`,
  `VECLIB_MAXIMUM_THREADS`, `R_DATATABLE_NUM_THREADS` set to 1 at spawn.
- Tasks are sent whole (ranks and gene sets per task); no `everywhere()`
  constants. For `run_enrichment()` the transcript-level gene sets of a
  collection are cut down to the IDs of each category's ranking before the
  task is built (`sets_in_ranking()`; fgsea drops the other IDs, so the
  result is identical). Before this, every task carried the whole
  collection: with all GENCODE v50 transcripts, Reactome given twice and 5
  categories (10 tasks), the session's peak RSS went from 898 MB
  (workers = 1) to 1838 MB (workers = 4), since `mirai_map()` serializes
  all tasks at dispatch.
- Daemons do not send warnings back. Each task collects its warning
  messages (`withCallingHandlers()`, muffled) and `run_tasks()` raises them
  again in the session, once per distinct message and task, naming it
  (class `isoformic_worker_warning`). The seeded in-session path does the
  same, so `workers` does not change which warnings are seen.
- Daemons that the user started for the profile get the session's
  `.libPaths()` first for the duration of each task (restored after); their
  thread settings cannot be changed after start-up, so `?run_enrichment`
  shows the `withr::with_envvar()` call that starts them like
  `start_workers()`.
- `mirai (>= 2.3.0)` in Suggests and in `check_installed()`:
  `daemons_set()` first appeared in 2.3.0. Tested with 2.7.4.
- `"auto"` without parallelly falls back to one worker with a message
  (`rlang::is_installed()`, not `check_installed()`).

## Implementation check (2026-11-18)

The implemented `workers` argument (code as above, mirai 2.7.4, parallelly
1.50.0) on W1, all runs pinned with `taskset -c 0-3`, 3 repeats per cell in
shuffled order, one process per run in its own capped cgroup scope, the
whole batch under the heavy-job lock, memory sampled from outside as in
*Method*. The session holds the 855 MB experiment as before. Two shapes of
the same 18 fgsea runs (`run_productivity_enrichment(rank_by = "delta_z")`,
`set.seed()` once before the first call):

- `loop`: 6 calls, one per `productive_fraction()` table (contrast), each
  with the 3 collections: 3 tasks per call, so `workers = 4` and `"auto"`
  are capped at 3, and daemons start and stop in every call.
- `one`: 1 call on one table with 18 collections (the 3 collections, 6
  copies each under different names): 18 tasks, daemons started once.

Wall time is the in-process time of the enrichment calls, including daemon
start-up. Two batches ran: wall times are from the first (1-minute load 6.1
to 9.3 on 16 cores, tight ranges); memory from the second, whose sampler
counted the scope's processes correctly (in the first it latched onto the
batch's own scope; its cgroup peaks agree with the second within 0.01 GB).
The second batch's wall times are not used: the load rose to 20 during it
(ranges up to 2x).

| Shape | `workers` | Used | Wall s, median (range) | PSS GB | cgroup peak GB | Summed RSS GB |
|---|---|---|---|---|---|---|
| loop | 1 | 1 | 14.04 (13.46-14.51) | 1.00 | 0.97 | 1.03 |
| loop | 2 | 2 | 23.02 (22.74-23.11) | 1.60 | 1.57 | 1.67 |
| loop | 4 | 3 | 19.97 (18.97-22.38) | 1.89 | 1.87 | 1.98 |
| loop | `"auto"` | 3 | 21.07 (20.12-24.93) | 1.89-1.92 | 1.87 | 1.98 |
| loop, daemons kept | 4 daemons | 3 | 8.51 (8.33-9.31) | 1.98-2.01 | 1.96 | 2.08 |
| one | 1 | 1 | 14.39 (13.95-14.50) | 1.00 | 0.98 | 1.04 |
| one | 2 | 2 | 10.15 (9.77-10.90) | 1.68 | 1.66 | 1.75 |
| one | 4 | 4 | 7.67 (7.22-8.17) | 2.27 | 2.25 | 2.38 |
| one | `"auto"` | 4 | 7.78 (7.10-7.94) | 2.27 | 2.24 | 2.37 |

Memory ranges across the 3 repeats were within 0.03 GB. "Daemons kept":
`mirai::daemons(4, .compute = "isoformic")` started once by the script
with the settings that `?run_enrichment` shows (start-up included in the
time), then the 6 calls with `workers = 4`; a separate batch of 3 runs at
load 5.3 to 6.8.

- `"auto"`: `parallelly::availableCores()` and `resolve_workers("auto")`
  gave 4 under `taskset -c 0-3` in all 27 runs (affinity read back as
  `0-3`), and `"auto"` used 4 daemons in `one` and the task cap of 3 in
  `loop`; its times and memory match `workers = 4`.
- Results: `workers = 1` was `identical()` across repeats in both shapes;
  `workers = 2`, `4` and `"auto"` were `identical()` to each other in every
  repeat of both shapes. Against `workers = 1`, `ES`, `size` and
  `leadingEdge` (and the other non-sampled columns) were `identical()`;
  `pval`, `padj`, `log2err` and `NES` differ, as documented: the largest
  absolute p-value difference was 0.12, the correlation of log10 p-values
  0.992 (`loop`) and 0.993 (`one`), and the gene sets with `padj < 0.05`
  were the same (2 in `loop`, 0 in `one`; their p-values moved by at most
  0.55 on the log10 scale). Kept daemons gave the same `ES`, `size` and
  `leadingEdge`, identical across repeats, but other p-values than
  per-call daemons, because the script started them after `set.seed()`
  and `mirai::daemons()` draws from the session's stream (verified).
- Speed: with 18 tasks in one call, 4 workers gave 1.9x and 2 workers 1.4x
  over serial. With 3 tasks per call, starting daemons per call costs more
  than it saves: `loop` was 1.4x to 1.6x slower than serial with any
  worker count. Kept daemons removed that cost (1.65x faster than serial at
  3 tasks per call). The speed-up is below the prototype's 2.3x to 2.6x
  (`mirai`, 4 workers) because each call starts its own daemons, which load
  fgsea before their first task.
- Memory: each daemon added about 0.3 GB of PSS (1.00 GB serial, 1.68 GB
  with 2, 2.27 GB with 4 in `one`). Summed RSS stayed within 5% of PSS,
  so tools that add up RSS no longer overcount, unlike `MulticoreParam`
  (4.44 GB summed RSS for 2.39 GB PSS in W1).
- Consequence for the documentation: `workers` pays off when one call has
  many tasks (several collections and transcript categories); for loops of
  small calls (one contrast and a few collections each) users should start
  the daemons once, as `?run_enrichment` shows, or keep `workers = 1`. A
  sentence saying so in the *Workers* section, and in the articles that
  loop over contrasts, would prevent a slowdown.

## Kept daemons (2026-11-20)

The finding above led to a change of design: `run_tasks()` no longer starts
and stops daemons in every call. The first call with `workers` above 1
starts daemons on the `"isoformic"` compute profile (one per task, up to
`workers`) and keeps them for the session. Later calls reuse them, launch
the missing ones with `mirai::launch_local()` when a call uses more (and
wait for them to connect), and start them again when a call requests fewer
than are running, so that no call runs more tasks at once than its
`workers`. `isoformic_workers(n)` starts `n` daemons ahead (used as they
are by every call) and `isoformic_workers(0)` stops them; daemons that the
user started with `mirai::daemons()` are used as they are and never
stopped by isoformic. `.onUnload` stops isoformic's own daemons; at the end
of a session they exit with it (mirai's `autoexit`; an `onexit` finalizer
that called `mirai::status()` segfaulted at exit and was dropped). An
interrupted or failed wait cancels the remaining tasks with
`mirai::stop_mirai()` (verified: without it, a call after an interrupt
queued 8.8 s behind the abandoned tasks; with it, 0.002 s). The per-task
seeds are drawn as before and daemon start-up runs under
`withr::with_preserve_seed()`, so results do not depend on whether the
daemons were already running. Note: a SIGINT sent with `kill` to a
background `Rscript` did not stop the wait on mirai results (it did for a
foreground one), so the test of the interrupt path raises the interrupt
in R; the test of session end uses child sessions that end normally, on an
error and by SIGKILL, and checks that their daemons are gone.

Same harness, data and shapes as *Implementation check* (W1, `taskset -c
0-3`, 3 repeats per cell in shuffled order, own capped cgroup scope per
run, the batch under the heavy-job lock; 1-minute load 1.4 to 5.3 on 16
cores). Wall time includes daemon start-up.

| Shape | `workers` | Used | Wall s, median (range) | Before (per-call daemons) | PSS GB | cgroup peak GB | Summed RSS GB |
|---|---|---|---|---|---|---|---|
| loop | 1 | 1 | 12.81 (12.77-13.05) | 14.04 | 1.03 | 0.97 | 1.03 |
| loop | 2 | 2 | 11.59 (11.28-11.87) | 23.02 | 1.65 | 1.60 | 1.69 |
| loop | 4 | 3 | 7.92 (7.79-8.02) | 19.97 | 1.94 | 1.89 | 2.00 |
| one | 1 | 1 | 13.44 (13.30-13.61) | 14.39 | 1.03 | 0.97 | 1.03 |
| one | 2 | 2 | 9.30 (9.15-9.49) | 10.15 | 1.69 | 1.64 | 1.73 |
| one | 4 | 4 | 6.50 (6.40-6.72) | 7.67 | 2.27 | 2.23 | 2.36 |

- The loop of small calls is now faster than serial with any worker
  count: 1.11x with 2 workers (was 0.61x) and 1.62x with 4 (3 used; was
  0.70x), about the 8.51 s the script-started daemons gave. One call with
  18 tasks gained 1.4x (2) and 2.07x (4) over serial in this batch; the
  serial times are lower than in the first batch (lower load), so compare
  ratios within a batch.
- Memory is unchanged: about 0.3 GB of PSS per daemon once it has run
  `fgsea`. Kept daemons hold that memory while idle, until
  `isoformic_workers(0)` or the end of the session.
- Results: every repeat of every cell was `identical()` to the output of
  the per-call implementation (same inputs, `set.seed()` once), and
  `workers = 2` and `4` were `identical()` to each other.

### Health check of kept daemons (2026-11-20)

Review of the kept daemons found that they were never health-checked: the
profile's dispatcher and address outlive its daemons, so after all daemons
were killed between two calls (SIGKILL, as an OOM kill or a memory cap
does), `ensure_workers()` saw "own daemons, enough of them" and the next
call queued its tasks for no daemon and never returned; when one of two
died, later calls ran on the survivor without saying so. Changes:

- `ensure_workers()` counts the connected daemons
  (`mirai::status()$connections`), not the ones it launched. With none
  connected, isoformic's daemons are started again; with fewer than the
  call needs, the missing ones are launched (`isoformic_workers(n)` daemons
  are kept at `n`). Every start or launch waits until the daemons connect
  (at most 30 s), so the next call never counts daemons still starting as
  dead.
- Daemons the user started are still used as they are, but a call fails
  with `isoformic_worker_error` (naming `isoformic_workers(0)`) when none
  connects within 30 s, instead of waiting forever.
- `wait_tasks()` polls the tasks (every 10 ms, asking the dispatcher for
  the connections once a second) instead of blocking on `m[]`: when tasks
  are still waiting and no daemon has been connected for 30 s (all died
  during the call), the call fails, its remaining tasks are cancelled, and
  the next call starts the daemons again.
- A call with a smaller `workers` than the daemons running restarts them
  with one per task, up to its `workers` (as before); the documentation
  now says so instead of "with its count".

The reviewers' repro scripts now finish: after killing both daemons the
next call returned in 0.5 s on two new daemons; after one died, the next
call ran on two PIDs again. A child-session test covers the three cases
(all killed between calls, one killed, all killed during a call) with a
time limit, checks `identical()` results and that no daemon outlives the
session.

Same harness, plan and data as *Kept daemons* (3 repeats, shuffled,
`taskset -c 0-3`, under the heavy-job lock; 1-minute load 0.8 to 1.5).

| Shape | `workers` | Used | Wall s, median (range) | Before the check | PSS GB | cgroup peak GB | Summed RSS GB |
|---|---|---|---|---|---|---|---|
| loop | 1 | 1 | 12.73 (12.70-12.76) | 12.81 | 1.03 | 0.97 | 1.03 |
| loop | 2 | 2 | 11.31 (11.30-11.33) | 11.59 | 1.66 | 1.61 | 1.70 |
| loop | 4 | 3 | 7.86 (7.78-7.88) | 7.92 | 1.95 | 1.90 | 2.02 |
| one | 1 | 1 | 13.13 (13.12-13.19) | 13.44 | 1.02 | 0.97 | 1.03 |
| one | 2 | 2 | 9.19 (9.17-9.21) | 9.30 | 1.69 | 1.65 | 1.74 |
| one | 4 | 4 | 6.49 (6.44-6.50) | 6.50 | 2.28 | 2.23 | 2.36 |

- The polling wait and the connection checks cost nothing measurable: all
  cells are within 0.3 s of the previous batch (lower load here). The loop
  of small calls stays faster than serial: 1.13x with 2 workers and 1.62x
  with 4 (3 used); one call with 18 tasks 1.43x (2) and 2.02x (4).
- Results: every repeat of every cell was `identical()` to the same cell of
  the previous kept-daemons batch, and `workers = 2` and `4` were
  `identical()` to each other.
- Seeds and `set.seed()`: unchanged. The per-task seeds are drawn in the
  session before dispatch and every daemon start, launch and connection
  wait runs under `withr::with_preserve_seed()`, so restarting dead daemons
  does not move the stream either (the child-session test checks a result
  after a restart against the in-session run). Only daemons that the user
  starts with `mirai::daemons()` after `set.seed()` change the p-values
  that follow, because `mirai::daemons()` draws from the session's stream;
  `?run_enrichment` says so and recommends `isoformic_workers()`.
