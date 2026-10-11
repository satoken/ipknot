# Nussinov sparsification and DD threading

Baseline: `454db3b` (`dev`), GCC 13.3.0, Release (`-O3 -DNDEBUG`), solver-free
build on Linux aarch64. Measurements use a fixed seed, one warmup, and the median
of five timed runs. The machine exposes Cortex-X925 and Cortex-A725 cores;
affinity was not fixed. The fixtures below are synthetic candidate graphs,
with signed weights and lonely pairs allowed, rather than full RNA prediction
benchmarks. DD cases include graph construction, 20 iterations, and ordinary
primal recovery. DP cases reuse one workspace, excluding construction.

## Applicability of SparseRNAfolD

[Gray, Will and Jabbari (2024)](https://doi.org/10.1186/s13015-024-00256-4)
optimize pseudoknot-free minimum free energy prediction with Turner energies and
dangles. IPknot's DD oracle instead maximizes additive, iteration-dependent pair
weights on each pseudoknot-free level. The paper's decomposition criterion can
be adapted to this max-plus recurrence; its MFE/dangle recurrences are not a
replacement for the oracle or for the partition-function engines supplying BPP.
This implementation independently applies that criterion to the existing DP;
it does not import SparseRNAfolD or implement its trace-arrow garbage collection.

Let `C(k,j)` be the best complete closed component, and `F(i,j)` the best interval
score. An ordinary split has score `F(i,k-1) + C(k,j)`. Process `k` in descending
order. If the closed component cannot improve the interval optimum from an
unpaired right base or already processed later components, omit its prefix-split
updates. Concatenating its dominating decomposition with any prefix is feasible
and no worse. Signed weights do not invalidate this argument. For the no-lonely
grammar, retain the separate `P(k,j)` score and trace even when the complete
component is dominated or negative: a parent may need it for stacking support.

Exact DP also discards suffix-dominated starts. A later start with an equal or
better score can replace an earlier start in every future inside/prefix query.
Neither form of sparsification counts as approximate beam pruning. The exact
oracle score is preserved, with possible changes between tied structures and
therefore potentially different DD trajectories. The decomposition criterion is
disabled for finite beams because its replacement derivations can be pruned.

Per decode, with `C` retained closed candidates, the work is bounded by
`O(n^2 + n*C + M)` after the endpoint rows have been sorted once. Worst-case time
remains cubic; chart/trace space remains `O(n^2 + M)`. This change does not make
DD globally exact or remove crossing/witness candidate limits.

## Parallelism

Both exact Nussinov and the beam/LinearNussinov-style oracle can run concurrently
across independent DD levels. `--dd-threads N` / `DDOptions::threads` controls an
optional OpenMP loop; the default is 1, 0 uses the runtime limit, and worker count
is capped at the number of levels. Each worker owns its decoder workspace and
output; aggregation and multiplier updates keep serial ordering. Worker
exceptions are rethrown after joining. No global OpenMP settings are changed.

Two levels provide at most two-way oracle parallelism. Upper-level primal
recovery depends on recovered lower levels and remains sequential. Single-level
decoding receives no threading speedup from this implementation. Standard exact
interval Nussinov can instead use an interval-length wavefront, but that requires
a different scheduling/chart implementation from this sparse left-to-right DP.
The beam algorithm's endpoint progression and pruning depend on previous
columns, so endpoints cannot simply run in parallel. Fine-grained work within a
column would require merges and a synchronization point at each endpoint; this
change implements the coarser independent-level approach.

Other changes avoid unnecessary paired-state buffers and hash lookup construction
when lonely pairs are allowed, reuse constrained-DD endpoint maxima and
per-level pair indices, reuse recovery masks, query exact inside intervals with
one monotone scan per column, and avoid repeated sorting of dominated columns.

## Measurements

| Case | Size / mode | Baseline, 1 thread | New, 1 thread | New, 2 threads | Best speedup |
| --- | --- | ---: | ---: | ---: | ---: |
| Exact DP, sparse | n=1500, M=17922 | 18.08 ms | 15.22 ms | — | 1.19x |
| Exact DP, dense | n=500, M=124750 | 22.23 ms | 6.16 ms | — | 3.61x |
| Beam DP, sparse | n=4000, M=47922, beam=100 | 21.03 ms | 21.15 ms | — | 0.99x |
| Exact DD | n=1500, M=35844, 2 levels | 1.102 s | 0.693 s | 0.427 s | 2.58x |
| Beam DD | n=4000, M=95844, 2 levels, beam=100 | 1.182 s | 1.188 s | 0.780 s | 1.51x |

All cases above preserved baseline objectives and upper bounds, and serial and
parallel results matched exactly. The standalone beam kernel shows no meaningful
speedup. On the sparse exact-DP fixture, `/usr/bin/time` reported peak process RSS
of 30256 KiB before and 12348 KiB after (about 59% lower). This includes process
and fixture memory, and is not a chart-only measurement. Results depend on input
density, stacking, level count, and available cores; probability engines are
outside these timings.

Raw timings are in [results.json](results.json). To reproduce, build both
revisions with `cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF`, and compile the
same harness against each revision's headers and library:

```sh
c++ -std=c++17 -O3 -I"$BASELINE/src" benchmark.cpp \
  "$BASELINE_BUILD/decoder/libipknot_decoder.a" -o /tmp/bench-before
c++ -std=c++17 -O3 -fopenmp -DIPKNOT_PERFORMANCE_THREADS -I"$CANDIDATE/src" \
  benchmark.cpp "$CANDIDATE_BUILD/decoder/libipknot_decoder.a" -o /tmp/bench-after
python3 run_benchmark.py --baseline /tmp/bench-before \
  --candidate /tmp/bench-after --output results.json
```

Run these commands from this directory and set the four paths to absolute
checkout/build directories. Use GCC/OpenMP for the second harness command.
`--dd-threads 2` requires an OpenMP-enabled build; `-DENABLE_OPENMP=OFF` keeps a
serial-only library. The installed package's CMake export discovers OpenMP for
consumers when enabled.

### Saved RNA BPP files

The local checkout also supplies ten saved `bpp/*.bpp` inputs. With `-x -r 0
-t .2,.1`, the default no-lonely-pair rule, default 50-iteration limit, and one
warmup plus five timed batches, the median wall times for processing all ten
files are:

| Mode | Baseline, 1 thread | New, 1 thread | New, 2 threads | Best speedup |
| --- | ---: | ---: | ---: | ---: |
| Beam | 289.37 ms | 288.13 ms | 225.69 ms | 1.28x |
| Improved beam | 76.37 ms | 66.35 ms | 64.79 ms | 1.18x |
| Exact Nussinov | 321.56 ms | 76.36 ms | 68.82 ms | 4.67x |

All ten structures, objectives, upper bounds, and iteration counts match the
baseline in all three modes and both new thread configurations. Timings include
process startup and writing BPSEQ/trace files, but exclude BPP computation.
Short improved-beam inputs show limited parallel benefit; repeated exploratory
runs sometimes made two threads slower. This is why threading defaults to 1.
These inputs are local data, not added to the branch. Per-file results and
structure hashes are in [bpp_results.json](bpp_results.json). Reproduce using:

```sh
python3 run_bpp_benchmark.py --baseline /absolute/path/to/old/ipknot \
  --candidate /absolute/path/to/new/ipknot --input-dir /absolute/path/to/bpp \
  --output bpp_results.json
```

## Validation

The existing 19-test suite passes, plus `dd_sparse_parallel`, which checks
repeated masked/signed decodes against an independent dense interval recurrence,
confirms decomposable candidates are omitted without beam pruning, compares
serial/parallel plain and constrained DD across all three modes, and exercises
worker exception propagation. Existing short-instance exhaustive enumeration
also checks feasibility and trace scores, including negative stacked pairs.
The CLI tests exercise thread parsing, 0/1/2 thread behavior, and invalid values.
The OpenMP-disabled build passes the sparse/parallel, API, and installed-consumer
tests with its serial fallback and explicit rejection of multiple threads.
The targeted sparse/parallel and exhaustive-DD tests also pass ASan and UBSan
in a Debug build. This host required
`ASAN_OPTIONS=detect_leaks=0:detect_stack_use_after_return=0`: LeakSanitizer is
unsupported under its ptrace environment, and fake-stack checking caused
timeouts. Leak and stack-use-after-return detection are therefore not covered.
