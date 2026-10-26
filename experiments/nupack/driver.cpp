// Standalone kernel benchmark: read an RNA sequence on stdin.
#include <nupack.h>
#ifndef NUPACK_REFERENCE
#include <linear_nupack.h>
#endif
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <sys/resource.h>

int main(int argc,char** argv)
{
  if (argc<2) return 2;
  std::string sequence;
  std::getline(std::cin,sequence);
  double logz=0;
  std::size_t states=0,edges=0;
  std::vector<std::tuple<int,int,double>> pp;
  const auto start=std::chrono::steady_clock::now();
  if (std::string(argv[1])=="exact")
  {
    Nupack<long double> model;
    model.load_default_parameters(); model.load_sequence(sequence);
    logz=std::log(model.calculate_partition_function());
    model.calculate_posterior();
    std::vector<float> bp; std::vector<int> offset;
    model.get_posterior(bp,offset);
    for (int i=0; i<static_cast<int>(sequence.size()); ++i)
      for (int j=i+1; j<static_cast<int>(sequence.size()); ++j)
        if (bp[offset[i+1]+j+1]>0) pp.emplace_back(i,j,bp[offset[i+1]+j+1]);
  }
#ifndef NUPACK_REFERENCE
  else if (std::string(argv[1])=="linear")
  {
    LinearNupack model(argc>2 ? std::stoul(argv[2]) : 100);
    logz=model.calculate(sequence); pp=model.posterior();
    states=model.statistics().retained_states; edges=model.statistics().retained_edges;
  }
#endif
  else return 2;
  const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  rusage usage{}; getrusage(RUSAGE_SELF,&usage);
  long rss=usage.ru_maxrss;
#ifdef __APPLE__
  rss/=1024;
#endif
  std::cout << std::setprecision(17) << "{\"length\":" << sequence.size()
    << ",\"seconds\":" << seconds << ",\"rss_kib\":" << rss
    << ",\"log_z\":" << logz << ",\"states\":" << states << ",\"edges\":" << edges;
  if (argc>3 && std::string(argv[3])=="dump")
  {
    std::cout << ",\"pairs\":[";
    bool comma=false;
    for (auto [i,j,p]:pp)
    {
      if (comma) std::cout << ',';
      comma=true;
      std::cout << '[' << i << ',' << j << ',' << p << ']';
    }
    std::cout << ']';
  }
  std::cout << "}\n";
}
