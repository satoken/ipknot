IPknot for predicting RNA pseudoknot structures using integer programming
=========================================================================

Requirements
------------

* C++17 compatible compiler (tested on Apple clang version 14.0.0 and GCC version 10.2.1) 
* cmake (>= 3.8)
* pkg-config
* [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/) (>= 2.2.0)
* one of these MIP solvers
	* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
	* [Gurobi Optimizer](http://www.gurobi.com/) (>=8.0)
	* [ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (>=12.0)
	* [SCIP](https://scipopt.org/) (>= 8.0.3)
	* [HiGHS](https://highs.dev/) (>= 1.5.0)

Install
-------
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

### Run with Docker

	docker build . -t ipknot
	docker run -it --rm -v $(pwd):$(pwd) -w $(pwd) ipknot ipknot examples/RF00005.fa

### Run with the web server

The IPknot web server is available at http://ws.sato-lab.org/rtips/ipknot++/.


References
----------

* Sato, K., Kato, Y.: Prediction of RNA secondary structure including pseudoknots for long sequences, *Brief. Bioinform.*, 23(1):bbab395 (Jan. 2022)
* Sato, K., Kato, Y., Hamada, M., Akutsu, T., Asai, K.: IPknot: fast and accurate prediction of RNA secondary structures with pseudoknots using integer programming, *Bioinformatics*, 27(13):i85-i93 (Jul. 2011)
