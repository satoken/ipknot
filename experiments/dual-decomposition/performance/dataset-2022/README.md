# Official 2022 dataset performance study

See [results/RESULTS.md](results/RESULTS.md) for measured speed, accuracy, and limits.
The performance baseline is `454db3b`; the optimized implementation is `eb7cc1c`.
Both use a solver-free GCC 13.3 Release build. The frozen Linux aarch64 binaries
and their hashes are included with the results. No production code changes are
part of this follow-up benchmark.

The full pass runs all 11,843 published single sequences with CONTRAfold/exact
Nussinov and LPC/LinearNussinov, plus the existing improved-beam variant. The
timing passes repeat the established, hash-selected 36 sequences three times.
It compares the old binary, optimized serial decoder, independent-level threads,
and threads inside one DP. Four-core groups must be disjoint and homogeneous;
serial variants use the first core of each group. All three phases run separately.
The first timing pass fixes OpenMP waiting to ACTIVE to compare all backends.
The second uses the normal runtime default, with `OMP_WAIT_POLICY` and
`GOMP_SPINCOUNT` unset, and compares the optimized serial and principal parallel
settings. Recommendations use this second pass's paired serial/parallel timings.

From the repository root on a compatible Linux aarch64 host, use the saved
binaries below. Alternatively build the two source revisions with Release and
`ENABLE_ILP=OFF`, and place the binaries and compiled `measure.c` in `$RUN/bin`.

```sh
HERE="$PWD/experiments/dual-decomposition/performance/dataset-2022"
HELPER="$PWD/experiments/nupack/benchmark-2022/dataset.py"
RUN=/tmp/ipknot-dataset-performance-reproduction
mkdir -p "$RUN/bin"
for file in baseline candidate measure; do
  gzip -dc "$HERE/results/$file.gz" > "$RUN/bin/$file"
  chmod +x "$RUN/bin/$file"
done
unset OMP_WAIT_POLICY GOMP_SPINCOUNT OMP_THREAD_LIMIT
python3 "$HELPER" prepare --root "$RUN" --archive /path/to/ipknot_dataset.zip
cp "$HERE/benchmark_dataset.py" "$RUN/full-runner.py"
python3 "$HERE/benchmark_dataset.py" --root "$RUN" --helper "$HELPER" \
  --baseline "$RUN/bin/baseline" --candidate "$RUN/bin/candidate" \
  --launcher "$RUN/bin/measure" --phase full --wait-policy PASSIVE \
  --modes contra-serial lpc-serial improved-serial
cp "$HERE/benchmark_dataset.py" "$RUN/timing-runner.py"
python3 "$HERE/benchmark_dataset.py" --root "$RUN" --helper "$HELPER" \
  --baseline "$RUN/bin/baseline" --candidate "$RUN/bin/candidate" \
  --launcher "$RUN/bin/measure" --phase timing --wait-policy ACTIVE \
  --sample-per-stratum 3 --repetitions 3 --cpu-groups '5,6,7,8;15,16,17,18' \
  --modes contra-before contra-serial contra-level2 contra-wave1 contra-wave4 \
    lpc-before lpc-serial lpc-level2 lpc-inner4 \
    improved-before improved-serial improved-level2 improved-inner4
python3 "$HERE/benchmark_dataset.py" --root "$RUN" --helper "$HELPER" \
  --baseline "$RUN/bin/baseline" --candidate "$RUN/bin/candidate" \
  --launcher "$RUN/bin/measure" --phase default-policy --wait-policy DEFAULT \
  --sample-per-stratum 3 --repetitions 3 --cpu-groups '5,6,7,8;15,16,17,18' \
  --modes contra-serial contra-level2 contra-wave4 \
    lpc-serial lpc-level2 improved-serial improved-level2
python3 "$HERE/diagnose_columns.py" --root "$RUN" --helper "$HELPER" --cpu 0
cp "$HERE/benchmark_dataset.py" "$RUN/benchmark_dataset.py"
python3 "$HERE/report_dataset.py" --root "$RUN" --helper "$HELPER" \
  --output "$RUN/published"
```

Adjust CPU numbers to the host. The full pass defaults to the ten X925 cores
5–9 and 15–19. Optional historical accuracy archives can be supplied to the
report with `--historical-contra` and `--historical-lpc`; those older binaries
are not used as the matched speed baseline. Resume checks reject changed
inputs, binaries, options, and scripts. Predictions are revalidated on export.
