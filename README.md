IPknot for predicting RNA pseudoknot structures using integer programming
=========================================================================

Requirements
------------

* [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/) (>= 2.2.0)
* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
  <!-- [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0)
  or [ILOG CPLEX](http://www.ibm.com/software/products/ibmilogcple/) (>=12.0) -->

Install
-------

	export PKG_CONFIG_PATH=/path/to/viennarna/lib/pkgconfig:$PKG_CONFIG_PATH
	mkdir build && cd build
	cmake -DCMAKE_BUILD_TYPE=Release .. && make 
	make install # optional

Usage
-----

### Single sequences

IPknot can take FASTA formatted RNA sequences as input, then predicts their secondary structures including pseudoknots.

	% ipknot: [options] fasta
	 -h:       show this message
	 -t th:    threshold of base-pairing probabilities for each level
	 -g gamma: weight for true base-pairs equivalent to -t 1/(gamma+1)
               (default: -g 4 -g 8)
     -m:       use McCaskill model (default: CONTRAfold model)
     -i:       allow isolated base-pairs
     -b:       output the prediction via BPSEQ format

	% ipknot drz_Ppac_1_1.fa
	>drz_Ppac_1_1
	GACUCGCUUGACUGUUCACCUCCCCGUGGUGCGAGUUGGACACCCACCACUCGCAUUCUUCACCUAUUGUUUAAUUGUGCUUGUGGUGGGUGACUGAGAAACAGUC
	.((((((..[[..[[..(((.......)))))))))....((((((((((..((((..((............))..))))..)))))))))).((.....]]))]]

### Aligned sequences

IPknot can also take CLUSTAL formatted RNA alignments produced by CLUSTALW and MAFFT, then predicts their common secondary structures.

	% clustalw RF00005:0.fa
	% ipknot RF00005:0.aln
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
	.........(((((((............((((((((.((.......))))))))))....((((....((........))....))))...........)))))))

This example shows folding with constraints that 16th base and 100th base are paired, 41st and 42nd bases are unpaired.

References
----------

* Sato, K., Kato, Y., Hamada, M., Akutsu, T., Asai, K.: IPknot: fast and accurate prediction of RNA secondary structures with pseudoknots using integer programming, *Bioinformatics*, 27(13):i85-i93 (Jul. 2011)
