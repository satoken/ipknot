IPknot for predicting RNA pseudoknot structures using integer programming
=========================================================================

Requirements
------------

* [Vienna RNA package](http://www.tbi.univie.ac.at/~ivo/RNA/) (>= 1.8)
* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/) (>=4.41)
  or [Gurobi Optimizer](http://www.gurobi.com/) (>=2.0)
  or [ILOG CPLEX](http://http://www-01.ibm.com/software/integration/optimization/cplex/) (>=12.0)

Install
-------

For GLPK,

	./configure --with-vienna-rna=/path/to/vienna-rna --with-glpk

For Gurobi, 

	./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi=/path/to/gurobi

For CPLEX,

	env CPPFLAGS='-I/path/to/cplex/include' LDFLAGS='-L/path/to/cplex/lib' \
	./configure --with-vienna-rna=/path/to/vienna-rna --with-cplex

Usage
-----

IPknot can take FASTA formatted RNA sequences as input, the
predict their secondary structures including pseudoknots.

	% ipknot: [options] fasta
	 -h:       show this message
	 -t th:    threshold of base-pairing probabilities for each level
	 -g gamma: weight for true base-pairs equivalent to -t 1/(gamma+1)
               (default: -g 4 -g 8)
     -m:       use McCaskill model (default: CONTRAfold model)
     -i:       allow isolated base-pairs
     -b:       output the prediction via BPSEQ format

	% ipknot ASE00001.fa
	> ASE_00001
	gaggaaagucccgccUCCAGAUCAAGGGAAGUCCCGCGAGGGACAAGGGUAGUACCCUUGGCAACUGCACAGAAAACUUACCCCUAAAUAUUCAAUGAGGAUUUGAUUCGACUCUUACCUUGGCGACAAGGUAAGAUAGAUGAAGAGAAUAUUUAGGGGUUGAAACGCAGUCCUUCCCGGAGCAAGUAGGGGGGUCAAUGAGAAUGAUCUGAAGACCUCCCUUGACGCAUAGUCGAAUCCCCCAAAUacagaagcgggcuua
	.....(((.(((((.((........[[(([[[[[[..((]]]](((((((...)))))))...(((((...........((((((((((((((..............((...((((((((((....))))))))))..))......))))))))))))))......))))).]]))]]...[[....(((((((((....(((....)))...))))))))).(((]]...)))...))...........)).)))))))).

References
----------

* Sato, K., Kato, Y., Hamada, M., Akutsu, T., Asai, K.: IPknot: fast and accurate prediction of RNA secondary structures with pseudoknots using integer programming, *Bioinformatics*, 27(13):i85-i93 (Jul. 2011)
