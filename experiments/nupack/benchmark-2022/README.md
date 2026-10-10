# 2022 IPknot dataset: both NUPACK engines vs LPC/LPV

Read [RESULTS.md](RESULTS.md) for the current comparison and completion status. All 11,843 official
single-sequence BPSEQ entries are included; no accuracy-based selection or
parameter tuning is performed. Ambiguous symbols are preserved and cannot pair.

| Engine | Posterior beam | Refinement | Per-sequence wall limit | Address-space limit |
|---|---:|---:|---:|---:|
| NUPACK (nonlinear full DP) | — | 0 | 60 s | 8 GiB |
| LinearNUPACK | 100 | 0 | 1,800 s | 8 GiB |
| LPC | 100 | 1 | 1,800 s | 8 GiB |
| LPV | 100 | 1 | 1,800 s | 8 GiB |

The decoder is shared: DD, beam DP width 100, 50 iterations, patience 0,
crossing width 100, 16 witnesses, automatic thresholds, one thread, no extra
PK score or structure constraints. Each RNA's engines run on the same pinned
Cortex-X925 CPU, in a deterministic shuffled order. Input order is longest
first. Parameters were not chosen using benchmark accuracy.

Exact NUPACK runs that exceed the requested resource limits are recorded as
failures, not assigned empty predicted structures. Accuracy and successful
latencies are reported on two separate common sets: the three linear engines,
and all four engines. The latter set is biased toward short sequences that
full NUPACK can finish. `summary.json` also records success rates and a
diagnostic all-attempt macro F1 with failures assigned zero. It does not replace
the common-set comparison. Timing includes the posterior, decoding and the
specified refinement; summed child wall time is distinct from benchmark elapsed
time.

## Reproduce

From this worktree's root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF -DBUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure

RUN=/tmp/ipknot-nupack-2022-reproduction
SCRIPTS=experiments/nupack/benchmark-2022
python3 "$SCRIPTS/dataset.py" prepare --root "$RUN" --archive /path/to/ipknot_dataset.zip
cc -O2 -Wall -Wextra -Werror "$SCRIPTS/measure.c" -o "$RUN/measure"
NUPACK_BENCHMARK_LAUNCHER="$RUN/measure" python3 "$SCRIPTS/test_evaluation.py"
python3 "$SCRIPTS/evaluate.py" --root "$RUN" --binary "$PWD/build/ipknot" \
  --launcher "$RUN/measure" --phase final --engines lnupack lpc lpv nupack \
  --cpus 5,6,7,8,9,15,16,17,18,19 --exact-timeout 60 --memory-gib 8
python3 "$SCRIPTS/report.py" --root "$RUN" --output "$RUN/report"
```

Choose available homogeneous CPUs on other hosts. `evaluate.py` resumes only
when the binary, input, scripts, options and CPU configuration match the saved
manifest. `report.py` independently rechecks scores, pair validity, ambiguous
endpoints and adjacent stack support; it rejects missing runs by default. Add
`--allow-incomplete --no-bootstrap` for a clearly marked progress report.

The source archive contains the modified engine files; `source.patch` records
tracked changes relative to the base commit in `provenance.json`. The actual
Release executable is saved as `ipknot.gz`. `raw-results.jsonl.gz` retains every
prediction, command, exit code, stderr and child resource measurement;
`per-rna.csv` and `summary.json` include per-sequence and grouped results. The
reference index and inventory identify the exact input set and its checksums.
The official archive checksum is verified during preparation. It contains
11,843 single sequences, 945 fewer than the paper's Table 1; no entries in the
archive are filtered here.

The dataset preparation and resource launcher are adapted from the repository's
existing `experiments/dual-decomposition/paper2022` benchmark. The launcher adds
an explicit `RLIMIT_AS` limit and disables core dumps. Timeout and memory-limit
handling are covered by the standalone harness tests.

## Early full linear accuracy

`complete_linear.py` separately evaluates all RNAs of at most 190 nt on the
otherwise idle Cortex-A725 cores while the original X925 phase continues.
`linear_accuracy.py` combines these predictions with original-phase predictions
for longer RNAs, using this fixed length rule. It independently audits all
35,529 linear-engine predictions and publishes `LINEAR_ACCURACY.md` and
`linear-accuracy.json`. The frozen executable and all model/decoder settings
are identical. Timing and memory from these two CPU classes are not combined.
The original `final` phase remains the sole performance-comparison source.
The auxiliary scripts record the local paths; adjust these constants when
reproducing on another host.
