# NMR-constrained DD recovery (2026-10-20)

The new final recovery pass makes DD substantially closer to ILP on the existing
16-RNA NMR dataset. It retains the existing solution and accepts only feasible
assignments with a higher original objective. Default settings enable it in
both solver-free and ILP-linked application builds.

| Fixed-threshold comparison | Previous DD | Updated DD |
|---|---:|---:|
| Identical base-pair set to ILP, among 14 completed ILP runs | 8/14 | **14/14** |
| Mean reference base-pair F1, all 16 RNAs | 69.88% | **88.51%** |
| Objective improvements relative to previous DD | — | **8 improved, 8 unchanged, 0 worse** |
| Total decoding time, sum of per-RNA medians over 3 runs | 3.181 s | **3.231 s** |

The measured time increase is 1.59%. This includes process startup, model
construction and decoding, using saved BPPs; posterior estimation is excluded.
The complete per-RNA results, objective values, penalties, structures and binary
hashes are in [comparison.csv](results/comparison.csv),
[summary.json](results/summary.json), and [results.json](results/results.json).
The two initial ILP timeouts (2LKR and 2N7M, 120 s each) are excluded from
ILP agreement statistics. Both also timed out at 600 s; these retries are
recorded separately in `additional-results.json`.

## What failed and what changed

The independent DD oracles often select RNA pairs that do not yet form a valid
assignment with the count constraints and NOE witnesses. The existing general
integer repair branches across both RNA and auxiliary variables. Its default
2,048-node budget could retain a very poor incumbent: 9BLM returned the empty
structure with objective -0.5, compared with ILP's 4.837981.

After ordinary recovery, the new pass:

1. Proposes RNA matchings using original scores, averaged DD scores, and
   upper-level-aware scores. It tries both strict stacking and ordinary planar
   matching, completing pseudoknot levels in order.
2. Fixes each distinct proposed RNA matching and searches for consistent
   auxiliary counts, NOE witnesses, topology choices and soft penalties.
3. Checks all variable domains, recorded constraints and implicit planarity
   before adopting a better original objective.

It shares at most 2,048 additional search nodes across six proposals. On these
16 RNAs it used at most 143 additional nodes per solve. The original repair
budget and frontier, dual trajectory, upper bounds and subgradients are
preserved. Completing a restricted search never counts as proving the original
model optimal. The ordinary incumbent remains available if a proposal fails.
This is a feasible-solution heuristic; it does not make finite-budget DD exact.

Simply doubling the old general repair budget to 4,096 nodes did not fix the
main failures. For example, 9BLM still returned -0.5; 2N8V improved only from
1.848585 to 2.130493, versus 5.042596 with the new pass.
Increasing retained witnesses to 128 was also tested and was not adopted as a
default: it costs more and does not consistently improve finite-budget results.
These ablations are recorded in [additional-results.json](results/additional-results.json).

## Remaining difference from ILP

The updated default matches ILP's **base-pair set** in all 14 completed runs.
Its **objective and NOE violation count** match in 10/14. The four exceptions
have identical RNA structures but fewer usable explanations in the retained
NMR candidate graph:

| RNA | Updated bounded DD objective | Full DD / ILP objective | Bounded / full NOE violations |
|---|---:|---:|---:|
| 2N8V | 5.042596 | 5.242596 | 3 / 1 |
| 2N6W | 8.227112 | 8.327112 | 2 / 1 |
| 2MHI | 4.608538 | 4.708538 | 2 / 1 |
| 5IEM | 5.576673 | 5.676673 | 1 / 0 |

Using `--dd-constraints full` with the new pass matched ILP's structures,
objectives and violation counts in all four cases. This diagnostic option
removes candidate restrictions and may be much more expensive on larger inputs.
Neither an open dual gap nor agreement on this small corpus proves equality
with unrestricted ILP on other sequences.

## Additional settings and correctness

| Setting, all 16 RNAs | Previous mean reference F1 | Updated mean reference F1 |
|---|---:|---:|
| NOE penalty 1.0, fixed thresholds | 65.18% | 87.29% |
| NOE and exact HSQC counts, both penalties 0.1, fixed thresholds | 74.09% | 87.58% |
| NOE penalty 0.1, `auto,auto` thresholds | 81.15% | 87.48% |

Both additional fixed-threshold settings had 8 objective improvements and no
decreases. Automatic thresholds compare different models: preserving each
subsolve's incumbent does **not** guarantee better reference accuracy of the
finally selected model. In that experiment 2N6W declined from 95.83% to 92.00%
reference F1. Automatic runs store per-subsolve diagnostics; a single final
objective or violation count is not inferred from the last subsolve.

Validation:

- HiGHS-linked Release build: **86/86 CTests passed**.
- Solver-free Release build: **20/20 CTests passed**, including the installed
  external decoder consumer.
- The two builds produced identical structures and objectives on all 16 RNAs.
- Disabling the pass reproduces the previous DD structures and objectives.
- Synthetic regression models have independently known optima and test competing
  NOE observations, integer count slacks, tiny budgets, both constraint modes,
  preservation of the ordinary search budget and dual bound, and bounded
  propagation work.
- `audit.py` independently verifies recorded domains, integer values, row sums,
  matching, planarity, objective values and state budgets for the feasible
  incumbents and final solutions in the fixed, stronger-penalty, count and
  automatic-threshold experiments: **19,074 feasible incumbent checks across
  164 trace files**, including disabled-pass and solver-free comparisons.

## Settings and reproduction

Existing NMR commands automatically use the new pass. Its explicit setting is:

```sh
ipknot --decoder dd --dd-dp improved-beam --dd-nmr-polish-states 2048 \
  ...your existing NMR options... input.fa
```

`--dd-nmr-polish-states 0` disables it. The original
`--dd-constraint-states 2048` budget remains separate. A larger polish budget is
only a cap; searches can finish earlier, and repeated proposals are skipped.
Trace summaries expose `polish_states`, `polish_calls`, and
`polish_improvements`. `DDOptions::nmr_polish_states` has the same default.
The standalone probability-to-structure decoder API does not handle NMR input;
its existing unconstrained behavior is unchanged.

The baseline is commit `fb0071d6629b06027d936130c8d724a740826e40`. Both binaries
were built in Release mode with the same compiler and HiGHS library. The host
is ARM64 (Cortex-X925/Cortex-A725); each RNA was pinned to the same CPU in
5–8 across versions, with one solver/OpenMP thread. Timed repetitions alternated
version order. ILP timeout retries used separate CPUs and are not part of the
DD timing comparison.

The bundled [inputs.json](inputs.json) contains the existing `new_nmr`
sequences, references, NOE patterns and HSQC counts. The 16 saved MXfold2 BPP
files in [posteriors](posteriors) come from the existing
`ipknot-dd-nmr-l2-study/p0.1-w64` study, with comment headers removed for `-x`.
Lengths are 47–111 nt; 13 RNAs have NOE patterns and 3 are controls without them.
No reference structure is supplied to the predictor. The primary experiment
uses fixed thresholds `0.25,0.25`, no refinement, soft NOE penalty 0.1, all
bulges, coaxial stacking, and the no-adjacent-bulge prior; HSQC counts are used
only in the separate count experiment. DD keeps beam 100, 50 iterations,
64 NMR witnesses and the solver-free relaxed NOE factor.

From the repository root:

```sh
study=experiments/dual-decomposition/nmr-quality
python3 "$study/benchmark.py" --binary /path/to/baseline/ipknot \
  --output /tmp/nmr-before --modes dd ilp --trace-state
python3 "$study/benchmark.py" --binary /path/to/updated/ipknot \
  --output /tmp/nmr-after --modes dd --trace-state
python3 "$study/audit.py" /tmp/nmr-before /tmp/nmr-after
```

Choose available cores with `--cpus 0,1,2,3` if needed. For timing, omit
`--trace-state` and run each version three times into distinct output folders.
Use `summarize.py --help` to combine them. Each result records commands, input
and binary hashes, and CPU assignment; incompatible resumes are rejected.

Variants:

```sh
# Stronger NOE penalty:
python3 "$study/benchmark.py" --binary /path/to/updated/ipknot \
  --output /tmp/nmr-strong --modes dd --penalty 1 --trace-state
# Add the original HSQC counts:
python3 "$study/benchmark.py" --binary /path/to/updated/ipknot \
  --output /tmp/nmr-counts --modes dd --counts --trace-state
# Investigate candidate loss in the four discrepant cases:
python3 "$study/benchmark.py" --binary /path/to/updated/ipknot \
  --output /tmp/nmr-full --modes dd-full --pdbs 2N8V 2N6W 2MHI 5IEM --trace-state
# Reproduce the previous recovery behavior:
python3 "$study/benchmark.py" --binary /path/to/updated/ipknot \
  --output /tmp/nmr-disabled --modes dd --extra --dd-nmr-polish-states 0
```
