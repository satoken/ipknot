// Offline diagnostic harness: arbitrary signed sparse instances, no ILP link.
#include "dual_decomposition.h"
#include <fstream>
#include <iostream>
int main(int argc, char** argv) {
  try {
    if (argc != 13) throw std::runtime_error("instance trace iterations beam crossing witnesses step projected patience state unpruned diminishing");
    std::ifstream input(argv[1]);
    int n, levels, lonely; double obsolete_score;
    if (!(input >> n >> levels >> lonely >> obsolete_score)) throw std::runtime_error("Invalid instance header");
    std::vector<DDPair> pairs; DDPair p;
    while (input >> p.left >> p.right >> p.level >> p.weight) pairs.push_back(p);
    DDOptions o; o.trace_file=argv[2]; o.max_iterations=std::stoi(argv[3]); o.beam=std::stoi(argv[4]);
    o.crossing_beam=std::stoi(argv[5]); o.witnesses=std::stoi(argv[6]); o.step=std::stod(argv[7]);
    o.projected_norm=std::stoi(argv[8]); o.patience=std::stoi(argv[9]); o.trace_state=std::stoi(argv[10]);
    o.unpruned_bound=std::stoi(argv[11]); o.diminishing_step=std::stoi(argv[12]);
    std::ofstream(o.trace_file).close();
    if (obsolete_score != 0) throw std::runtime_error("PK-specific scores are no longer supported");
    auto result=solve_dual_decomposition(n,pairs,levels,lonely,o);
    std::cout << result.objective << ' ' << result.upper_bound << ' ' << result.stop_reason << '\n';
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
