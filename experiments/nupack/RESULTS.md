# NUPACK optimization and linear approximation measurements

Host: Linux AArch64, GCC 13.3, Release (`-O3 -DNDEBUG`), solver-free build. Timings are medians of three separate processes on a shared host; memory is peak RSS. All sequences use the recorded deterministic seed. These are kernel measurements, excluding the IPknot decoder.

## Exact computation

The reference is `b0f965c` with the same default-parameter alignment and single-pair pseudoknot posterior repairs. This isolates optimization performance from those correctness changes. Partition functions and float posterior values matched the repaired reference exactly for every run below. Raw, unrepaired historical predictions can differ.

| Bases | Reference (s) | Optimized (s) | Speedup |
| ---: | ---: | ---: | ---: |
| 20 | 0.0148 | 0.0145 | 1.02× |
| 40 | 0.7679 | 0.1148 | 6.69× |
| 60 | 10.2461 | 1.2204 | 8.40× |

## Beam scaling

Fixed beam width and fixed loop bounds give expected linear time and memory in sequence length. The stored outside forest has constants depending on the beam; larger beams can be expensive. In particular the 2,000-base / beam-100 run used about 1.9 GiB.

| Bases | Beam | Time (s) | Peak RSS (MiB) | Retained states | Retained edges |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 250 | 25 | 0.338 | 43.9 | 37317 | 665161 |
| 250 | 100 | 1.776 | 171.7 | 131852 | 3207578 |
| 500 | 25 | 0.825 | 90.1 | 76793 | 1545920 |
| 500 | 100 | 5.121 | 414.6 | 285411 | 8387631 |
| 1000 | 25 | 1.809 | 168.3 | 156502 | 3042573 |
| 1000 | 100 | 13.402 | 934.4 | 599505 | 19214520 |
| 2000 | 25 | 3.482 | 334.5 | 313091 | 6208328 |
| 2000 | 100 | 30.161 | 1914.2 | 1217895 | 39568257 |

## Approximation error

Errors compare against the full repaired NUPACK DP, with one independent random sequence per length. MAE averages over all possible base-pair coordinates, including zero probabilities. This is a regression and performance sample, not a biological accuracy benchmark. Beam pruning can lose substantial probability mass; increasing the beam does not guarantee monotonic accuracy.

| Bases | Beam | Time (s) | BPP MAE | Maximum BPP error | log Z difference |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 25 | 0.00072 | 1.17158e-05 | 0.00030258 | -0.000468572 |
| 20 | 100 | 0.00109 | 1.39582e-06 | 3.52441e-05 | -5.50864e-05 |
| 40 | 25 | 0.00734 | 0.000198343 | 0.0112279 | -0.0191067 |
| 40 | 100 | 0.01752 | 0.000143494 | 0.00753862 | -0.0137055 |
| 60 | 25 | 0.02811 | 0.0119176 | 0.890848 | -3.15734 |
| 60 | 100 | 0.08956 | 0.000713878 | 0.11313 | -0.147363 |

The approximate grammar permits two crossing bands with arbitrarily long stems and bounded bulges/interior loops, and structured or pseudoknotted gaps. It excludes multiloops inside crossing bands. The retained middle-gap starts and leading multiloop/gap padding are also restricted. Disabling pruning retains these restrictions. Consequently `LinearNUPACK` is opt-in and does not replace full `NUPACK`.

## Verification and reproduction

- CTest: 17/17 tests passed in a Release build with `ENABLE_ILP=OFF` and no ViennaRNA dependency.
- Short-sequence oracle: 64 random sequences at each length 0–12, comparing inside values and all posterior coordinates with full NUPACK.
- Forced hairpin, stacked, bulged, and multiloop weights checked independently by subset partition-function inversion.
- Probability normalization, deterministic pruning, invalid constraints, 1,002-base forced hairpins, and two fully paired nine-pair crossing stems checked.
- AddressSanitizer and UndefinedBehaviorSanitizer run on the NUPACK test executable (leak detection disabled because of the host tracing environment).

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
python3 experiments/nupack/benchmark.py --output experiments/nupack/results.json
```

The benchmark uses only the Python standard library, the local compiler and Git. Raw measurements are in [results.json](results.json). The left-to-right search, prefix-aware pruning and log-space marginal computation were informed by [LinearFold](https://github.com/LinearFold/LinearFold/blob/master/src/LinearFold.cpp) and [LinearPartition](https://github.com/LinearFold/LinearPartition/blob/master/src/LinearPartition.cpp); no additional runtime dependency was added.
