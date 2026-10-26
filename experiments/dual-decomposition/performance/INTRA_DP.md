# Parallelism inside one Nussinov DP

This follow-up uses `1b12107` on `codex/decoder-performance` as its baseline,
including the earlier exact sparsification and independent-level OpenMP loop.
The original [SparseRNAfolD applicability discussion](README.md#applicability-of-sparsernafold)
still applies: these are additive pair-weight DD oracles, with optional stacking
support, rather than Turner-energy MFE or partition-function recurrences.

The [official 2022 RNA dataset study](dataset-2022/results/RESULTS.md) evaluates
these backends with CONTRAfold and LPC, including posterior computation,
automatic thresholds and refinement. It measures the usual beam width 100,
checks whether column workers actually activate, and compares serial and
parallel settings under the normal OpenMP waiting policy.

For LPC-DD, the recommended backend is `--dd-dp improved-beam` at posterior
and decoder beam width 100. Use `--dd-threads 2 --dd-dp-threads 1` for
single-RNA latency, or one thread per process for batches of independent RNAs.
The standard beam cases below remain comparisons of the original algorithm.

## Exact wavefront

`--dd-dp wavefront --dd-dp-threads 4` selects the new exact interval backend.
Library callers set `DDOptions::wavefront_dp=true` and `dp_threads=4`.
The configured beam is ignored. `--dd-dp nussinov` keeps the earlier sparse,
left-to-right exact chart; it can still use independent-level threads.

Intervals of length `span` depend only on shorter intervals. A persistent team
distributes those intervals, with a barrier after each length. Within a length,
each right endpoint belongs to exactly one worker. Its descending-left candidate
cursor and retained closed-component list can therefore advance without locks.
The [OpenMP worksharing-loop specification](https://www.openmp.org/spec-html/5.2/openmpse66.html)
defines the team work distribution and implicit barrier used here.

This backend preserves closed-component sparsification. At `(i,j)`, first try
an unpaired `j` and previously retained components with left endpoints `k>i`.
Retain `C(i,j)` only when it improves that interval. A dominated component cannot
improve any longer prefix split. Parent-supported `P(i,j)` states are retained
even if negative or dominated, because an enclosing pair may need that stacking
support. Iterative traceback distinguishes these states from complete intervals.

After endpoint sorting, work remains `O(n^2+n*C+M)` for `C` retained components;
worst-case time is cubic. Unlike the sparse chart, this backend always allocates
the triangular interval table: `6*n*(n+1)` bytes for double scores and integer
choices, plus candidate/stack buffers. It is explicitly selected because this
storage can be wasteful for long sparse inputs. The exact optimal score is
preserved; ties can differ between backends. Thread counts within the wavefront
backend preserve both score and structure.

## LinearNussinov columns

`--dd-dp beam --dd-dp-threads 4` or `--dd-dp improved-beam --dd-dp-threads 4`
distributes closed-component calculations and prefix splits within one column.
Every worker reads completed earlier columns and writes private proposals.
Paired-state trace slots are reserved before the workers run, so there are no
worker allocations, shared vector resizes, or critical sections in the split
loop. A leader merges proposals, applies dominance/beam pruning, and materializes
winning traces before the next endpoint. One team persists across all columns.

Carry states win score ties; among pair derivations the smaller original pair ID
wins. This reproduces the ordinary serial beam's first-encountered derivation
and the improved beam's explicit tie rule, including shuffled input IDs.
Thread counts preserve structures, scores and numbers of beam-pruned states.
Finite-beam decoding does not apply the exact component-sparsification rule.

The initial prototype parallelized all columns. On the sparse n=4000,
degree=12, beam=100 fixture, it took 28.27 ms with two workers and 46.26 ms with
four, versus 21.82 ms serial. Synchronization and merging outweighed the split
work. The final implementation estimates average prefix offers using
`sum_allowed min(beam,left+1)/n`. Below 8192 it stays serial; otherwise worker
count is capped by `average/4096` and the request. These are conservative
heuristics, not a guarantee for every mask, input or machine. Private buffers
add `O(T*n)` space for `T` workers, alongside `O(n*beam+M)` chart/trace storage.

## Thread policy

Both thread settings default to 1. Zero uses the OpenMP runtime limit, with
level-count or beam-work caps. `dp_threads != 1` takes precedence over
`--dd-threads` / `DDOptions::threads`; levels run in order while each DP can
create its own team. Nested teams and global runtime changes are avoided.
Recovery levels retain their dependencies, but their individual decodes use
the same settings, including when recovery workspaces are not shared.
Sparse exact remains serial within one level; select `wavefront` for inner
exact parallelism. Without OpenMP, both settings accept 0/1 and reject values
above 1. Wavefront also works serially in those builds.

`DDNussinov::used_workers()` and `parallel_steps()` expose the actual team and
processed columns/diagonals from the last decode, including beam fallback.
Workspaces are reusable sequentially; concurrent decodes need separate objects.

## Measurements

GCC 13.3.0, Release `-O3 -DNDEBUG`, solver-free Linux aarch64 build on the same
Cortex-X925/A725 host. These DP measurements were restricted to the four
Cortex-X925 CPUs 15–18 (`taskset -c 15-18`) to avoid core-type changes. Synthetic graphs use fixed seed
38433200, signed weights, one warmup and five timed repetitions; construction is
excluded and the workspace is reused. No build/test jobs run during timings.
Raw data and per-mode structure hashes are in [intra_results.json](intra_results.json).

| Case | Baseline 1 | New sparse/beam 1 | Wavefront 1 | Inner 2 | Inner 4 | Speedup vs baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Exact, sparse n=1500, degree=12 | 14.91 ms | 15.02 ms | 13.66 ms | 6.99 ms | 3.78 ms | 3.94x |
| Exact, dense n=1000 | 33.59 ms | 37.58 ms | 41.46 ms | 18.20 ms | 8.33 ms | 4.03x |
| Exact, stacked dense n=700 | 11.96 ms | 12.87 ms | 12.86 ms | 6.01 ms | 3.10 ms | 3.86x |
| Beam=100, sparse n=4000, degree=12 | 20.76 ms | 20.92 ms | — | 21.00 ms | 20.93 ms | 0.99x |
| Beam=512, sparse n=2000, degree=128 | 101.53 ms | 102.49 ms | — | 94.72 ms | 78.48 ms | 1.29x |
| Beam=512, dense n=1000 | 123.24 ms | 123.84 ms | — | 94.72 ms | 65.45 ms | 1.88x |
| Improved beam=512, dense n=1000 | 129.09 ms | 148.01 ms | — | 100.13 ms | 70.46 ms | 1.83x |
| Beam=512, stacked dense n=700 | 47.95 ms | 48.25 ms | — | 40.19 ms | 28.82 ms | 1.66x |

The exact sparse/dense rows compare the baseline sparse backend with the new
wavefront backend. The single-thread wavefront control separates scheduling
speedup from changing chart layout. Beam rows use the same algorithm and beam
on both revisions. The short/default sparse beam selects serial fallback.
These timings cover one DP, excluding BPP calculation and other DD factors.
Objectives agree with the baseline on every fixture. Within each backend,
structures, objectives and pruning counts match across thread counts; finite
beams also match the baseline structure hashes.

The ten saved local BPP files (lengths 51–394) are much shorter and include
process startup, DD factors, recovery and file I/O. They use `-x -r 0 -t .2,.1`
with the default no-lonely rule and 50-iteration budget. Affinity was not fixed
for this end-to-end check. One warmup and five timed batches give:

| Mode | Baseline 1 | New 1 | Inner 4 |
| --- | ---: | ---: | ---: |
| beam | 276.18 ms | 263.57 ms | 252.72 ms |
| improved-beam | 74.01 ms | 67.18 ms | 60.01 ms |
| nussinov | 69.51 ms | 75.06 ms | 73.94 ms |
| wavefront | 75.14 ms | 136.39 ms | 157.22 ms |

The wavefront baseline row uses sparse exact Nussinov. All ten structures,
objectives, bounds and iteration counts match the baseline and between thread
counts in all four modes. These small inputs do not benefit from inner threads;
wavefront is substantially slower than sparse exact here. Prefer the sparse
backend and level parallelism for this workload. Inputs remain local data;
[intra_bpp_results.json](intra_bpp_results.json) records per-file results.

Single-thread dense exact/improved-beam measurements still show overhead
relative to the baseline on this compiler/host, despite retaining the fused
serial kernel. The default sparse beam=100 remains within about 1%. The
measured inner-team benefit is therefore input-dependent, and both threading
settings remain opt-in.

To reproduce, build baseline and candidate with `-DCMAKE_BUILD_TYPE=Release
-DENABLE_ILP=OFF`, then run from this directory with absolute paths:

```sh
c++ -std=c++17 -O3 -fopenmp -I"$BASELINE/src" intra_benchmark.cpp \
  "$BASELINE_BUILD/decoder/libipknot_decoder.a" -o /tmp/intra-before
c++ -std=c++17 -O3 -fopenmp -DIPKNOT_INTRA_THREADS -I"$CANDIDATE/src" \
  intra_benchmark.cpp "$CANDIDATE_BUILD/decoder/libipknot_decoder.a" -o /tmp/intra-after
taskset -c 15-18 python3 run_intra_benchmark.py --baseline /tmp/intra-before \
  --candidate /tmp/intra-after --output intra_results.json
python3 run_bpp_benchmark.py --intra --baseline "$BASELINE_BUILD/ipknot" \
  --candidate "$CANDIDATE_BUILD/ipknot" --input-dir /absolute/path/to/bpp \
  --output intra_bpp_results.json
```

Select an available group of comparable cores on other machines.

## Correctness and build checks

The OpenMP Release build passes all 20 CTest tests. The extended
`dd_sparse_parallel` test compares 640 reused, signed/masked exact fixtures
against an independent dense interval recurrence; both wavefront thread counts
must match in score and traceback, including negative stacking support.
Dense n=384, beam=256 tests cover actual column teams, ordinary/improved beams,
stacking, shuffled IDs, tied weights, changing masks and all-masked fallback.
Trace checks enforce endpoint uniqueness, no crossings, stack support and score.
Plain and constrained DD compare serial/inner/level thread settings 0/1/2/8
across four modes, and shared/unshared recovery. Existing exhaustive short
tests continue to pass. CLI checks cover the new mode, thread parsing, negative
requests and missing-OpenMP rejection.

The OpenMP-disabled build passes five targeted tests: sparse/parallel, both DD
CLI suites, decoder API and installed consumer. Debug ASan/UBSan passes the
extended sparse/parallel and exhaustive-DD tests with
`ASAN_OPTIONS=detect_leaks=0:detect_stack_use_after_return=0`, as required by
this host's ptrace/fake-stack limitations. Leak detection, fake-stack checking
and a race detector are not covered.
