// Solver-free diagnostic harness for joint integer certificates and exchanges.
#include "dual_decomposition.h"
#include <fstream>
#include <iostream>
int main(int argc,char** argv) {
  try {
    if(argc!=24) throw std::runtime_error("instance trace iterations beam crossing witnesses state block shift strict joint_width joint_shift matching joint_states recovery_every cache share exchange_width exchange_passes exchange_every exchange_states global clusters");
    std::ifstream in(argv[1]); int n,levels,lonely;double contact;
    if(!(in>>n>>levels>>lonely>>contact)) throw std::runtime_error("Invalid instance header");
    std::vector<DDPair> pairs;DDPair p;
    while(in>>p.left>>p.right>>p.level>>p.weight) pairs.push_back(p);
    DDOptions o;o.trace_file=argv[2];o.max_iterations=std::stoi(argv[3]);o.beam=std::stoi(argv[4]);
    o.crossing_beam=std::stoi(argv[5]);o.witnesses=std::stoi(argv[6]);o.trace_state=std::stoi(argv[7]);
    o.patience=0;o.unpruned_bound=true;
    o.bound_block=std::stoi(argv[8]);o.bound_shift=std::stoi(argv[9]);o.bound_strict_stack=std::stoi(argv[10]);
    o.joint_bound_width=std::stoi(argv[11]);o.joint_bound_shift=std::stoi(argv[12]);o.joint_bound_matching=std::stoi(argv[13]);
    int js=std::stoi(argv[14]),es=std::stoi(argv[21]);
    if(js<0||es<0)throw std::runtime_error("Negative state budget");
    o.joint_bound_states=js;o.recovery_every=std::stoi(argv[15]);o.recovery_cache=std::stoi(argv[16]);
    o.recovery_share=std::stoi(argv[17]);o.exchange_width=std::stoi(argv[18]);o.exchange_passes=std::stoi(argv[19]);
    o.exchange_every=std::stoi(argv[20]);o.exchange_states=es;o.global_bound=std::stoi(argv[22]);o.joint_bound_clusters=std::stoi(argv[23]);
    std::ofstream(o.trace_file).close();
    if (contact != 0) throw std::runtime_error("PK-specific scores are no longer supported");
    const auto r=solve_dual_decomposition(n,pairs,levels,lonely,o);
    std::cout<<r.objective<<' '<<r.upper_bound<<' '<<r.stop_reason<<'\n';
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
