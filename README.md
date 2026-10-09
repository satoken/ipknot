IPknot for predicting RNA pseudoknot structures
=========================================================================

Requirements
------------

* C++17 compatible compiler (tested on Apple clang version 14.0.0 and GCC version 10.2.1) 
* cmake (>= 3.18)
* pkg-config
* [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/) (>= 2.2.0, optional for ViennaRNA integration)
* For ILP decoding, one of these MIP solvers (optional for dual decomposition)
	* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
	* [Gurobi Optimizer](http://www.gurobi.com/) (>=8.0)
	* [ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (>=12.0)
	* [SCIP](https://scipopt.org/) (>= 8.0.3)
	* [HiGHS](https://highs.dev/) (>= 1.5.0)

Install
-------
For solver-free dual decomposition:

```sh
cmake -S . -B build-dd -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF
cmake --build build-dd -j
ctest --test-dir build-dd --output-on-failure
build-dd/ipknot --decoder dd examples/drz_Ppac_1_1.fa
```

DD is the default decoder, including in solver-enabled builds. The defaults are
`--decoder dd --dd-beam 100 --dd-max-iter 50 --dd-patience 0
--dd-crossing-beam 100 --dd-witnesses 16`, with LPC beam 100,
`-t auto,auto`, and `-r 1`. The iteration limit applies separately to each
threshold/refinement solve.
`--decoder auto` uses a linked ILP solver when available and DD otherwise.
A solver-enabled build can select either `--decoder ilp` or `--decoder dd`.
DD uses sparse beam Nussinov per level, projected subgradient updates for
cross-level constraints, and feasible structure recovery. Fixed DP and crossing
budgets keep decoding linear in sparse input size. Its crossing witness graph
and automatic PK threshold metric are approximations; predictions may differ
from ILP. Fixed-structure (`-c`) and NMR constraints also work with DD,
including exact/lower-bound counts, soft penalties, bulges and coaxial stacking.
The default `--dd-constraints linear` uses bounded NMR candidates/witnesses,
implicit per-level planarity, and bounded propagation and primal repair. With
fixed widths, iteration/repair budgets, and fixed NMR observation count/length,
constraint construction and decoding take expected O(n+M) time and memory for
sequence length n and sparse input size M. Every returned structure is checked
against all hard constraints. Candidate limits may omit feasible structures;
`--dd-constraints full` retains the original full model for small diagnostics
and has no linear complexity guarantee.
Use `--dd-dp nussinov` to replace beam search with exact per-level
Nussinov interval DP. It ignores `--dd-beam`, retains every interval start,
and works with fixed and NMR constraints in both constraint modes. It has
no linear complexity guarantee: worst-case time is cubic and DP memory is
quadratic. Candidate/witness limits still apply to the retained graph.
The default `--dd-dp beam` retains the configurable beam decoder.
Use `--dd-dp improved-beam --dd-beam 100` for the fixed-width improved
beam used in the DD/LPC comparison. It removes dominated interval states
before beam selection and uses deterministic traceback ties. The crossing
and iteration budgets still apply. Heap and adaptive-beam experiments are
not enabled by this option.

Use `--dd-constraint-recovery-every 10` to resume constrained integer
recovery every 10 iterations (default 0: final-only recovery). Initial,
periodic and final recovery share `--dd-constraint-states`, a total visited
DFS-node limit per threshold/refinement solve. The search frontier and best
feasible solution persist between checkpoints; this does not multiply the
state budget by the number of DD iterations. `--dd-recovery-target baseline`
(default) preserves the original Polyak target; `best` also uses checkpoint
incumbents in that target. Recovery is opt-in because accuracy effects vary
with the input and state budget.
See [constraint decoding details](docs/dual-decomposition-constraints.md) for
supported options, bound interpretation and repair-budget diagnostics.
Learned/energy PK scores support `crossing/blocks` and `projected/blocks` in DD.
Projected scoring adds signed unary pair coefficients while retaining the
original crossing witness constraints, without adding contact product factors.
Its partner pool is sampled with the DD crossing beam before witness trimming.
See [DD projected measurements](experiments/pk-score/dd_projected.md) for the
coefficient definition and the comparison with crossing scoring.
See [PK score integration measurements](experiments/pk-score/dd_integration.md)
for the frozen shape/BPP, Dirks–Pierce and Cao–Chen comparisons with
improved beam 100 and at most 50 DD iterations.
The [sequence-free score evaluation](experiments/pk-score/sequence_free_evaluation.md)
implements a four-coefficient bounded model (`IPKNOT_PK_BOUNDED_V1`) for both
formulations. Its refinement-off evaluation on 900 RNAs does not establish an
accuracy gain; the model remains opt-in.
See [the design and validation notes](docs/dual-decomposition.md) for equations,
budget controls, exact diagnostic settings and benchmark results.

Optional joint integer certificates group distant PK endpoints into bounded
clusters. Add `--dd-unpruned-bound --dd-global-bound --dd-joint-bound 12
--dd-joint-clusters` for stronger stopping certificates. For additional feasible
structure recovery, use `--dd-exchange 12 --dd-exchange-passes 2` and optionally
`--dd-recovery-every 1000000` (first and before stopping). These options keep
fixed-budget sparse decoding linear. See [integer-gap analysis and measured
tradeoffs](docs/dual-decomposition-integer-gap.md) before enabling extra recovery.

For GLPK,

	export PKG_CONFIG_PATH=/path/to/viennarna/lib/pkgconfig:$PKG_CONFIG_PATH
	mkdir build && cd build
	cmake -DCMAKE_BUILD_TYPE=Release ..  # configure
	cmake --build . # build
	cmake --install . # install (optional)

For Gurobi, add ``-DENABLE_GUROBI`` to the configure step:

	cmake -DENABLE_GUROBI -DCMAKE_BUILD_TYPE=Release ..  # configure

For CPLEX, add ``-DENABLE_CPLEX`` to the configure step:

	cmake -DENABLE_CPLEX -DCMAKE_BUILD_TYPE=Release ..  # configure

For SCIP, add ``-DENABLE_SCIP`` to the configure step:

	cmake -DENABLE_SCIP -DCMAKE_BUILD_TYPE=Release ..  # configure

For HiGHS, add ``-DENABLE_HIGHS`` to the configure step:

	cmake -DENABLE_HIGHS -DCMAKE_BUILD_TYPE=Release ..  # configure

To use the optional MXfold2 integration (disabled by default), configure with ``-DMXFOLD2=ON`` and ensure that Python (interpreter and development headers), [pybind11](https://github.com/pybind/pybind11), and the MXfold2 Python package are available in your environment.


Usage
-----

### Single sequences

IPknot can take FASTA formatted RNA sequences as input, then predicts their secondary structures including pseudoknots.

	% ipknot: [options] fasta
	 -h:       show this message
	 -t th:    threshold of base-pairing probabilities for each level (default: auto,auto)
	 -e model: probabilistic model (default: LinearPartition-C)
	 -b:       output the prediction via BPSEQ format

	% ipknot drz_Ppac_1_1.fa
	>drz_Ppac_1_1
	GACUCGCUUGACUGUUCACCUCCCCGUGGUGCGAGUUGGACACCCACCACUCGCAUUCUUCACCUAUUGUUUAAUUGUGCUUGUGGUGGGUGACUGAGAAACAGUC
	.[[[[[...(((((((((((.......))))]]]]].((.((((((((((..((((....................))))..)))))))))).))....)))))))

#### Model

IPknot can calculate the base pairing probability using the following probability models:

* LinearPartition model with CONTRAfold parameters (`LinearPartition-C` or `lpc`) (default)
* LinearPartition model with ViennaRNA parameters (`LinearPartition-V` or `lpv`)
* McCaskill model with Boltzmann likelihood parameters (`Boltzmann`)
* McCaskill model with ViennaRNA parameters (`ViennaRNA`)
* CONTRAfold model (`CONTRAfold`)
* NUPACK model (`NUPACK`)

You can specify the model using `-e` option.

#### Thresholds

IPknot predicts the pseudoknot structure hierarchically. The `-t` option is used to specify the base pairing probability threshold for each level. For example, if you run `ipknot -t 0.25 -t 0.125 seq.fa`, the threshold for the first level is 0.25 and the threshold for the second level is 0.125.

IPknot can search for the best thresholds from multiple combinations of thresholds using *pseudo-expected accuracy*. For example, `ipknot -t 0.5_0.25 -t 0.25_0.125 seq.fa` searches for the combination of 0.5 and 0.125 for the first layer and 0.25 and 0.125 for the second layer, and outputs the secondary structure with the maximum pseudo-expected accuracy as the final prediction. By default, `-t auto -t auto` is specified, where `auto` means `0.5_0.25_0.125_0.0625`.
### Aligned sequences

IPknot can also take CLUSTAL formatted RNA alignments produced by CLUSTALW and MAFFT, then predicts their common secondary structures.

	% clustalw RF00005.fa
	% ipknot -e Boltzmann RF00005.aln
	>J01390-1/6861-6
	--------CAGGUUAGAGCCAGGUGGUU--AGGCGUCUUGUUUGGGUCAAGAAAUU-GUUAUGUUCGAAUCAUAAUAACCUGA-
	........(((((((..(((...........))).(((((.......)))))......(((((.......))))))))))))..

### Folding with constraints

IPknot can fold a given sequence or alignment with some constraints. The constraint is given by a 2-columned TSV file. The first column indicates the position *i* of the base to be constrained. If the second column is given by the number *j*, this line means a base-pair constraint, that is, *i*th base and *j*th base form a base pair. If the second column is respectively given by a character `x`, `|`, `<`, `>`, the *i*th base should be unpaired, paired with another base, paired with a downstream base, paired with a upstream base, respectively.

	% cat constraint.txt
	16 100
	41 x
	42 x
	% ipknot -c constraint.txt drz_Ppac_1_1.fa
	>drz_Ppac_1_1
	GACUCGCUUGACUGUUCACCUCCCCGUGGUGCGAGUUGGACACCCACCACUCGCAUUCUUCACCUAUUGUUUAAUUGUGCUUGUGGUGGGUGACUGAGAAACAGUC
	.........(((((((............(((((((.........[[[[[)))))))....((((...................]]]]])))).......)))))))

This example shows folding with constraints that 16th base and 100th base are paired, 41st and 42nd bases are unpaired.

NMR-derived base-pair and stacking constraints can be supplied with
`--base-pairs` and `--stack-constraint`, respectively. The
`--nmr-bulge-mode` option controls one-nucleotide-bulge candidates:

- `all` (default): consider direct and one-nucleotide-bulged instances.
- `fallback`: use bulged instances for an observation only when it has no
  structurally valid direct instance.
- `none`: consider directly adjacent base pairs only.

The existing `--without-nmr-bulge` flag remains an alias for
`--nmr-bulge-mode none`. Add `--coaxial-stacking` to also allow a
two-base-pair stacking constraint to match flush coaxial stacking between
helix termini in a multibranch loop:

	% ipknot --stack-constraint "GC AU" --coaxial-stacking sequence.fa

Without `--coaxial-stacking`, multibranch-loop coaxial candidates are not
generated and the previous stacking-constraint behavior is preserved.

NOE witnesses that cannot satisfy the existing ordinary endpoint-stacking
support, base-uniqueness, and coaxial topology constraints are removed before
creating their integer variables. This is a necessary-condition test, not a
posterior-probability cutoff: it preserves the feasible RNA structures and
optimal objective, although a solver can select a different tied optimum.
It uses sparse endpoint lists rather than a dense compatibility matrix and
does not apply the endpoint-stacking test when `--allow-isolated` is enabled.

For a stronger, experimental reduction, add
`--nmr-coaxial-no-adjacent-bulge` together with `--coaxial-stacking`:

	% ipknot --stack-constraint "GC AU" --coaxial-stacking \
	    --nmr-coaxial-no-adjacent-bulge sequence.fa

This requires one directly adjacent **selected** base pair on the stem-facing
side of each of the two observed coaxial terminal pairs. It removes contexts
without those continuations and prevents a bulge in that immediate stem step;
bulges elsewhere, including in the supporting third helix, remain allowed.
This option is off by default and can change predictions or make hard NOE
constraints infeasible. It is a local regular-stem approximation, not a claim
that bulges and coaxial stacking never coexist: single-residue bulges can
retain A-form-like continuity ([Ray et al., 2017](https://pubmed.ncbi.nlm.nih.gov/29032448/)),
and coaxial stacking in bulged RNA is experimentally documented
([Bou-Nader et al., 2025](https://www.nature.com/articles/s41467-025-57519-w)).

NMR observations are hard constraints by default. Add `--nmr-soft` to allow
base-pair count mismatches and unsatisfied stacking observations with finite
objective penalties:

	% ipknot --base-pairs "GC=2,AU=3" --stack-constraint "GC AU" \
	    --nmr-soft --nmr-count-penalty 1.0 --nmr-stack-penalty 1.0 sequence.fa

`--nmr-count-penalty` is charged for each missing or excess base pair, while
`--nmr-stack-penalty` is charged once for each unsatisfied stacking
observation. Both weights must be positive. Omitting `--nmr-soft` retains the
hard-constraint behavior. Base-pair uniqueness, pseudoknot-level topology,
stack/bulge witness validity, and coaxial topology remain hard in both modes.
When thresholds are selected automatically, the NMR violation penalty is also
subtracted from the pseudo-expected-accuracy selection score. Its influence on
that selection can be calibrated separately with
`--nmr-threshold-penalty-scale SCALE` (default `1`, finite and nonnegative).
`0` selects thresholds using model-based pseudo expected accuracy alone;
the NMR count/stack penalties still apply inside every integer program.
This experimental setting does not affect fixed-threshold predictions.

By default, `--base-pairs GC=2` means that the complete predicted structure
contains exactly two GC pairs. If the value is instead the number of observed
or assigned peaks, use `--nmr-count-mode lower-bound`; this interprets it as
at least two GC pairs and penalizes only a shortfall in soft mode:

	% ipknot --base-pairs "GC=2,AU=3" --nmr-count-mode lower-bound \
	    --nmr-soft sequence.fa

Different stacking observations use disjoint witness base pairs by default.
`--nmr-allow-shared-stack-pairs` relaxes this assignment rule for sensitivity
analysis when multiple observations may describe overlapping parts of the
same helix.

The NMR formulation groups pair implications and coaxial topology exclusions
by observation (globally for pair implications when sharing is disabled),
and factors common exclusions across coaxial support choices. These are
equivalent constraints. The default candidate filter removes only witnesses
proved infeasible; the optional no-adjacent-bulge prior also narrows the
allowed structures. `--loglevel info` reports the two pruning counts
separately, as well as witness count, compressed row counts, formulation time, and
optimization time. Difficult instances can still take a long time to solve.
The experimental `--nmr-coaxial-no-adjacent-bulge` option requires one
directly adjacent selected stem pair behind each observed coaxial terminus.
It requires `--coaxial-stacking` and applies to ILP and both DD constraint modes.
Bulges elsewhere, including the supporting third helix, remain allowed.

Changing `--nmr-bulge-mode` or disabling coaxial stacking changes the allowed
structures, so those switches are not equivalent performance optimizations.

NOE candidate construction also uses sparse pair-membership and geometric
range indices, caches type/position lookups per sequence, and groups witnesses
by observation. In `fallback` mode, bulged candidates are not enumerated when
direct instances already exist. These construction optimizations do not
introduce probability cutoffs or candidate-count caps. Candidate pruning now
changes the model's insertion order by omitting rejected witnesses. The HiGHS
backend transfers its matrix buffers instead of copying them
and releases construction columns incrementally.

The sparse, fixed-beam `lpc`/`lpv` posterior path is unchanged; no dense
sequence-length-squared table or NOE index is allocated on the ordinary path.
Sequence type indexing takes linear time/storage for a fixed alphabet;
enumeration time additionally depends on the number of matching pairs.
Geometric range-index storage is linear in candidate pairs, with pruning
rather than a worst-case output-sensitive query-time guarantee. All-pair NOE
patterns can still have quadratically many candidates, and exact integer
optimization has no linear-time guarantee. The posterior algorithm's linear
scaling must therefore not be interpreted as an end-to-end guarantee with NOE.

For result-preserving regression and performance measurements, use
[`noe_exact_benchmark.py`](new_nmr/noe_exact_benchmark.py) (`models`,
`predictions`, `fixtures`, `random`, or `formulations`) and
[`noe_sparse_benchmark.py`](new_nmr/noe_sparse_benchmark.py).
Both accept `--before-build` and `--after-build` for separately retained
baseline/optimized build trees. The model recorder is test-only and does not
perform optimization; `formulations` reports construction timing without
fingerprinting overhead, not prediction accuracy or full solver runtime.
The archived-data modes (`models`, `predictions`, and `formulations`) require
the local experiment archives, which are not distributed in this repository.
The `random` and `fixtures` modes can instead use the included test inputs.
See [the integration validation record](docs/pk-score-merge-20261009.md) for
the latest model and prediction comparisons.

### Experimental pseudoknot scores

Opt-in H-type scores can evaluate crossing stems during optimization,
project motif scores onto existing upper-level pair coefficients, restrict
their existing crossing constraints to compatible witnesses, or rerank the
solutions from automatic threshold search. All feature weights are zero
by default. The current presets are experimental decoder scores rather than
validated thermodynamic parameters. See
[the scoring experiment](experiments/pk-score/README.md) for the formulation,
runtime comparisons, limitations, and reproducible commands.
The [crossing-constraint experiment](experiments/pk-score/supported.md)
describes `--pk-h-formulation supported`: it adds no variables or constraint
rows, but narrows the feasible structure set to match projected bonuses.
The [aggregate crossing-score experiment](experiments/pk-score/crossing.md)
uses `--pk-h-formulation crossing` to score actual selected crossing edges
with at most one continuous auxiliary per existing support row, retaining
the original feasible pair assignments and adding no integer variables.
The [weight calibration experiment](experiments/pk-score/calibration.md)
uses the same H-shape features, tuning on bpRNA references and evaluating
separately on Rfam14.5 references with the ordinary automatic-threshold and
refinement workflow.
The [learned crossing-score experiment](experiments/pk-score/learned_implementation.md)
fits a small signed compatibility model using posterior support, competing
partners, stem lengths and loop spans. `--pk-learned-model FILE` applies its
per-pair correction through the same crossing auxiliaries. It is opt-in and
adds no integer variables, BPP calculations or decoder calls; the report
records separate tuning and held-out accuracy and runtime measurements.
The [refinement-off re-evaluation](experiments/pk-score/refinement_off.md)
compares the previous scores on identical initial BPP matrices with `-r 0`,
including newly fitted models and a separate held-out dataset.
The [shape/posterior integration and auxiliary reduction experiment](experiments/pk-score/hybrid_r0.md)
combines geometry and posterior corrections with `--pk-hybrid-shape` and
tests `--pk-crossing-simplify`, which removes auxiliaries only when the
original crossing support and base uniqueness imply an exact direct term.
The [crossing-score speed experiment](experiments/pk-score/crossing_speed.md)
profiles solver work and tests `--pk-crossing-hypograph`, which removes
lower product links while preserving the maximized score for every integer
structure.

### Run with Docker

	docker build . -t ipknot
	docker run -it --rm -v $(pwd):$(pwd) -w $(pwd) ipknot ipknot examples/RF00005.fa

### Run with the web server

The IPknot web server is available at http://ws.sato-lab.org/rtips/ipknot++/.


References
----------

* Sato, K., Kato, Y.: Prediction of RNA secondary structure including pseudoknots for long sequences, *Brief. Bioinform.*, 23(1):bbab395 (Jan. 2022)
* Sato, K., Kato, Y., Hamada, M., Akutsu, T., Asai, K.: IPknot: fast and accurate prediction of RNA secondary structures with pseudoknots using integer programming, *Bioinformatics*, 27(13):i85-i93 (Jul. 2011)
