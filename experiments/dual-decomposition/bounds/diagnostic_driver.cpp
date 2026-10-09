// Offline diagnostic harness: arbitrary signed sparse instances, no ILP link.
#include "dual_decomposition.h"
#include <fstream>
#include <iostream>
int main(int argc, char** argv) {
  try {
    if (argc != 18 && argc != 19) throw std::runtime_error("instance trace iterations beam crossing witnesses step projected patience state unpruned diminishing block bound_every global recovery_every recovery_mask [recovery_target_best]");
    std::ifstream input(argv[1]);
    int n, levels, lonely; double contact_score;
    if (!(input >> n >> levels >> lonely >> contact_score)) throw std::runtime_error("Invalid instance header");
    std::vector<DDPair> pairs; DDPair p;
    while (input >> p.left >> p.right >> p.level >> p.weight) pairs.push_back(p);
    DDOptions o; o.trace_file=argv[2]; o.max_iterations=std::stoi(argv[3]); o.beam=std::stoi(argv[4]);
    o.crossing_beam=std::stoi(argv[5]); o.witnesses=std::stoi(argv[6]); o.step=std::stod(argv[7]);
    o.projected_norm=std::stoi(argv[8]); o.patience=std::stoi(argv[9]); o.trace_state=std::stoi(argv[10]);
    o.unpruned_bound=std::stoi(argv[11]); o.diminishing_step=std::stoi(argv[12]);
    o.bound_block=std::stoi(argv[13]);o.bound_every=std::stoi(argv[14]);o.global_bound=std::stoi(argv[15]);
    o.recovery_every=std::stoi(argv[16]);o.recovery_mode=std::stoi(argv[17]);
    if (argc==19) o.recovery_target_best=std::stoi(argv[18]);
    std::ofstream(o.trace_file).close();
    PKScoreOptions pk; pk.crossing=pk.fixed_blocks=true; pk.intercept=contact_score;
    auto result=solve_dual_decomposition(n,pairs,levels,lonely,o,pk);
    std::cout << result.objective << ' ' << result.upper_bound << ' ' << result.stop_reason << '\n';
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
