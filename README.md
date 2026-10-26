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
build-dd/ipknot -e lpc --decoder dd --dd-dp improved-beam examples/drz_Ppac_1_1.fa
```

DD is the default decoder, including in solver-enabled builds. The defaults are
`--decoder dd --dd-dp improved-beam --dd-beam 100 --dd-max-iter 50 --dd-patience 0
--dd-crossing-beam 100 --dd-witnesses 16`, with LPC beam 100,
`-t auto,auto`, and `-r 1`. The iteration limit applies separately to each
threshold/refinement solve.
`--decoder auto` uses a linked ILP solver when available and DD otherwise.
A solver-enabled build can select either `--decoder ilp` or `--decoder dd`.
DD uses sparse beam Nussinov per level, projected subgradient updates for
cross-level constraints, and feasible structure recovery. Fixed DP and crossing
budgets keep decoding linear in sparse input size. Its crossing witness graph
and automatic PK threshold metric are approximations; predictions may differ
from ILP. Fixed-structure (`-c`) constraints also work with DD.
The default `--dd-constraints linear` uses bounded crossing support,
implicit per-level planarity, and bounded propagation and primal repair.
With fixed widths and iteration/repair budgets, constraint construction and
decoding take expected O(n+M) time and memory for sequence length n and sparse
input size M. Every returned structure is checked against all hard constraints.
Crossing limits may omit feasible structures; `--dd-constraints full` retains
the original full structural model for small diagnostics and has no linear
complexity guarantee.
Use `--dd-dp nussinov` to replace beam search with exact per-level
Nussinov interval DP. It ignores `--dd-beam`, removes dominated interval starts
and decomposable closed split candidates without beam pruning,
and works with fixed constraints in both constraint modes. It has
no linear complexity guarantee: worst-case time is cubic and DP memory is
quadratic. Candidate/witness limits still apply to the retained graph.
The default `--dd-dp improved-beam --dd-beam 100` selects the fixed-width improved
beam recommended for LPC-DD. Use `--dd-dp beam` to select the standard beam
decoder explicitly. Improved beam removes dominated interval states
before beam selection and uses deterministic traceback ties. The crossing
and iteration budgets still apply. Heap and adaptive-beam experiments are
not enabled by this option.

Recommended LPC-DD settings for reducing the latency of one RNA are:

```sh
build-dd/ipknot -e lpc --decoder dd --dd-dp improved-beam \
  --beam-size 100 --dd-beam 100 --dd-threads 2 --dd-dp-threads 1 \
  -t auto,auto -r 1 -n 1 input.fa
```

On all 11,843 sequences in the [2022 dataset study](experiments/dual-decomposition/performance/dataset-2022/results/RESULTS.md),
improved beam increased macro F1 from 0.615594 to 0.617121 and micro F1
from 0.508241 to 0.526274, with a 1.59x speedup in summed prediction time
under the same single-thread conditions. Both variants with two level workers
were also compared on 36 stratified RNAs, each repeated three times: improved
beam was 1.48x faster than standard beam. These are aggregate results; medium-length
macro F1 and crossing-pair F1 decreased slightly.
Two level workers require OpenMP. For batches of independent RNAs, or a build
without OpenMP, use `--dd-threads 1` with the same improved-beam settings.

For CONTRAfold-DD, the recommended sparse exact, serial backend is:

```sh
build-dd/ipknot -e CONTRAfold --decoder dd --dd-dp nussinov \
  --dd-threads 1 --dd-dp-threads 1 -t auto,auto -r 1 -n 1 input.fa
```

Use `--dd-threads 2` to solve independent DD levels concurrently in all
DP modes, including constrained DD. The default is 1; 0 uses the OpenMP runtime
limit, capped at the number of levels. Level scores are accumulated in the same
order as serial execution. Primal recovery remains sequential because upper
levels depend on the recovered lower levels. CMake enables OpenMP when available;
use `-DENABLE_OPENMP=OFF` for a serial build. A serial build accepts 0/1 threads
and rejects larger requests. See the [performance study](experiments/dual-decomposition/performance/README.md)
for sparsification, parallelism limits, and reproducible measurements.

Use `--dd-dp wavefront --dd-dp-threads 4` for exact interval-length wavefront
parallelism within a single level. It retains closed-candidate sparsification
and stacking support, but always allocates a quadratic interval chart (about
`6*n*(n+1)` bytes for scores and choices, plus candidate/stack buffers).
The sparse `nussinov` backend remains useful for long, sparse inputs and uses
level parallelism only. Exact backends may choose different tied structures.

`--dd-dp beam --dd-dp-threads 4` and `improved-beam` distribute closed components
and prefix splits within each column, then merge before pruning. Thread counts
preserve serial structures, scores and pruning decisions. Light inputs fall
back to serial processing; the default sparse beam=100 workload usually gains
more from level parallelism. Within-DP threads default to 1; 0 uses the runtime
limit, subject to work-size caps for beam mode. A request other than 1 overrides
`--dd-threads` to avoid nested teams. Recovery levels remain ordered, but each
recovery DP can use the requested inner workers.

Use `--dd-constraint-recovery-every 10` to resume constrained integer
recovery every 10 iterations (default 0: final-only recovery). Initial,
periodic and final recovery share `--dd-constraint-states`, a total visited
DFS-node limit per threshold/refinement solve. The search frontier and best
feasible solution persist between checkpoints; this does not multiply the
state budget by the number of DD iterations. `--dd-recovery-target baseline`
(default) preserves the original Polyak target; `best` also uses checkpoint
incumbents in that target. Periodic recovery is opt-in because accuracy effects
vary with the input and state budget.

Optional joint integer certificates group distant PK endpoints into bounded
clusters. Add `--dd-unpruned-bound --dd-global-bound --dd-joint-bound 12
--dd-joint-clusters` for stronger stopping certificates. For additional feasible
structure recovery, use `--dd-exchange 12 --dd-exchange-passes 2` and optionally
`--dd-recovery-every 1000000` (first and before stopping). These options keep
fixed-budget sparse decoding linear. These additional recovery options are
opt-in and may change predictions under a finite iteration budget.

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

To use the optional MXfold2 integration (disabled by default), configure with ``-DWITH_MXFOLD2=ON`` and ensure that Python (interpreter and development headers), [pybind11](https://github.com/pybind/pybind11), and the MXfold2 Python package are available in your environment.


Usage
-----

### Decode probabilities as a C++ library

The reusable `IPknot::decoder` target accepts supplied base-pair probabilities
and returns partners, pseudoknot levels, and the objective. It shares the CLI's
improved beam, sparse exact Nussinov, level parallelism, and within-DP parallelism,
and supports ILP without linking probability engines or the CLI.
The library defaults to improved beam, matching the CLI and recommended LPC-DD
settings. Set `options.dd.improved_beam = false` to select standard beam.
See [the decoder API and examples](decoder/README.md) for CMake integration,
installation, input formats, and complexity. Refinement and fixed-structure
constraints remain available through IPknot.

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
* Beam approximation with NUPACK energies and two crossing pseudoknot bands
  (`LinearNUPACK` or `lnupack`)

You can specify the model using `-e` option.

`NUPACK` retains the full pseudoknot dynamic program. Its hot paths now use
precomputed table offsets and loop/dangling energies, constant-time poly-C
checks, and a bounded cache of **exact** `expl` results keyed by the float
argument. These optimizations do not use energy quantization or fast-math.
The default parameter reader also consumes the NINIO upper bound correctly;
previously this value was uninitialized and subsequent parameters were shifted
by one. The posterior pass now includes single-pair pseudoknot bands and the
same band penalty as the inside pass. Predictions can change because of these
correctness repairs. BPSEQ restrictions are validated and stored only in the
upper triangle.

For long sequences, select the approximate model explicitly:

```sh
ipknot -e LinearNUPACK --beam-size 100 -r 0 sequence.fa
```

`LinearNUPACK` implements left-to-right beam search inspired by
[LinearFold](https://github.com/LinearFold/LinearFold/blob/master/src/LinearFold.cpp)
and [LinearPartition](https://github.com/LinearFold/LinearPartition/blob/master/src/LinearPartition.cpp),
using this repository's NUPACK energy functions. It computes inside/outside
values in log space on the retained derivation forest. The sparse posterior
path, including constrained refinement, avoids dense probability and
four-dimensional tables. Sequence and alignment inputs are supported, as are
parameter files through `-P`.

Both NUPACK engines preserve ambiguous nucleotide symbols and their positions,
while treating them as unpaired, like LPC/LPV. Sequence-specific loop bonuses
are omitted when the relevant nucleotide is unknown; small interior loops with
unknown bases use generic length/asymmetry energies.

The full 2022 single-sequence comparison uses refinement 0 for both NUPACK
engines and refinement 1 for LPC/LPV; see
[benchmark results and reproduction](experiments/nupack/benchmark-2022/README.md).

The approximation retains ordinary hairpins, bulges, interior loops and
multiloops, plus pseudoknots with two crossing bands. Both bands can have
arbitrary stem length and bounded bulges/interior loops; all three gaps can
contain secondary structures or further pseudoknots. Multiloops **inside a
crossing band** are excluded. Interior loops, left multiloop padding and
leading padding of structured gap states are limited to 30 unpaired bases.
Middle-gap starts come from up to 30 bases of padding or from retained
structured gap states. These bounds
are separate from the beam and are configurable through the C++ constructor.
There is no overall sequence window or maximum pair distance.

Each state family retains at most the requested beam width per right endpoint,
with deterministic tie-breaking and prefix-plus-inside ranking. For fixed
positive beam width and loop bound, expected time and memory are O(sequence
length). Constants depend on the beam and loop bound; retaining the forest's
edges for outside computation can use substantial memory at larger beams
(up to O(n(bL²+b²)) edges for length n, beam b and loop bound L).
Hairpin selection and posterior coordinate ordering also avoid global
comparison sorts. Dense output requested through the C++ matrix interface
still requires O(length²) space. `--beam-size 0` disables pruning but keeps
the grammar bounds; it is not equivalent to full `NUPACK` and has no linear
complexity guarantee. Larger beams provide a tunable approximation, without
an accuracy or monotonicity guarantee.

The correctness tests compare short inside/outside results with the full
pseudoknot DP, independently extract forced motif weights from exact subset
partition functions, check probability normalization and long-distance pairs,
and exercise sparse output and refinement. Reproduce kernel timings and
posterior errors with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_ILP=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
python3 experiments/nupack/benchmark.py --output experiments/nupack/results.json
```

The benchmark applies the same correctness repairs to its reference version
before checking exact posterior identity and speed. Measurements and their
limits are recorded in [experiments/nupack/RESULTS.md](experiments/nupack/RESULTS.md).

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

### Pseudoknot prediction

Ordinary pseudoknot prediction and automatic threshold selection remain
available in both ILP and DD. Experimental PK-specific geometry, learned,
energy, ranking and ensemble scores have been removed; the former
`--pk-*` options are no longer accepted.

### Run with Docker

	docker build . -t ipknot
	docker run -it --rm -v $(pwd):$(pwd) -w $(pwd) ipknot ipknot examples/RF00005.fa

### Run with the web server

The IPknot web server is available at http://ws.sato-lab.org/rtips/ipknot++/.


References
----------

* Sato, K., Kato, Y.: Prediction of RNA secondary structure including pseudoknots for long sequences, *Brief. Bioinform.*, 23(1):bbab395 (Jan. 2022)
* Sato, K., Kato, Y., Hamada, M., Akutsu, T., Asai, K.: IPknot: fast and accurate prediction of RNA secondary structures with pseudoknots using integer programming, *Bioinformatics*, 27(13):i85-i93 (Jul. 2011)
