= IPknot for predicting RNA pseudoknot structures using integer programming

== Requirements

* Boost C++ Library (>=1.42.0) ((<URL:http://www.boost.org/>))
* Vienna RNA package (>= 1.8) ((<URL:http://www.tbi.univie.ac.at/~ivo/RNA/>))
* GNU Linear Programming Kit (>=4.41) ((<URL:http://www.gnu.org/software/glpk/>))
  or Gurobi Optimizer (>=2.0) ((<URL:http://www.gurobi.com/))
  or ILOG CPLEX (>=12.0) ((<URL:http://http://www-01.ibm.com/software/integration/optimization/cplex/>))

== Install

For GLPK,
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-glpk

For Gurobi, 
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi

For CPLEX,
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-cplex

You may have to specify the include path and the library path by CPPFLAGS and LDFLAGS like
 env CPPFLAGS='-I/path/to/gurobi/include' LDFLAGS='-L/path/to/gurobi/lib' \
 ./configure --with-vienna-rna=/path/to/vienna-rna --with-gurobi

Then,
 make
 make install

== Usage

(({ipknot})) can take FASTA formatted RNA sequences as input, the
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

== References

* Sato, K., Kato, Y., Akutsu, T., Asai, K.: RNA pseudoknot prediction
  based on maximizing expected accuracy. in preparation.
