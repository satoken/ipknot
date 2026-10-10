#include <iostream>
#include <ipknot/decoder.h>
#include <string>

int main(int argc, char **argv) {
  ipknot::DecoderOptions options;
  options.thresholds = {0.2f, 0.1f};
  if (argc > 1) {
    const std::string mode = argv[1];
    if (mode == "nussinov")
      options.dd.nussinov_dp = true;
    else if (mode == "full-dd") {
      options.dd.nussinov_dp = true;
      options.dd.crossing_beam = options.dd.witnesses = 0;
    } else if (mode == "ilp")
      options.backend = ipknot::Backend::ILP;
    else if (mode != "beam")
      return 2;
  }
  try {
    // Two crossing stems; positions and returned partners are zero-based.
    const std::vector<ipknot::PairProbability> probabilities{
        {0, 7, 0.9f}, {1, 6, 0.9f}, {3, 10, 0.9f}, {4, 9, 0.9f}};
    const auto result = ipknot::Decoder(options).decode(11, probabilities);
    for (int i = 0; i < 11; ++i)
      std::cout << i << '\t' << result.bpseq[i] << '\t' << result.levels[i] << '\n';
    std::cout << "objective: " << result.objective << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
