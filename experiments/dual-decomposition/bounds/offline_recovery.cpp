// Frozen-trajectory primal recovery audit. This executable never updates DD
// multipliers: it only replays the supplied coefficients and selected states.
#include "dd_recovery.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    if (argc != 5)
      throw std::invalid_argument("Usage: dd_offline_recovery INPUT FREQUENCY MASK HISTORY");
    const int frequency = std::stoi(argv[2]);
    const unsigned methods = std::stoul(argv[3]);
    const int history = std::stoi(argv[4]);
    if (frequency < 1) throw std::invalid_argument("Recovery frequency must be positive");
    std::ifstream input(argv[1]);
    std::string version;
    input >> version;
    if (version != "DD_RECOVERY_REPLAY_V1") throw std::invalid_argument("Invalid replay version");
    int n, levels, no_lonely, beam;
    std::size_t count;
    input >> n >> levels >> no_lonely >> beam >> count;
    std::vector<DDPair> pairs(count);
    std::vector<unsigned char> allowed(count);
    for (std::size_t id = 0; id < count; ++id) {
      int permitted;
      input >> pairs[id].left >> pairs[id].right >> pairs[id].level >> pairs[id].weight >> permitted;
      allowed[id] = permitted;
    }
    std::size_t row_count;
    input >> row_count;
    std::vector<DDRecoveryRow> rows(row_count);
    for (auto& row : rows) {
      std::size_t contacts;
      input >> row.upper >> contacts;
      row.contacts.resize(contacts);
      for (auto& c : row.contacts) input >> c.first >> c.second;
    }
    std::size_t states;
    input >> states;
    if (!input) throw std::invalid_argument("Truncated replay problem");
    DDPrimalRecovery recovery(n, pairs, levels, beam, no_lonely, allowed, rows, history, methods);
    std::vector<double> weights(count);
    std::vector<unsigned char> selected(count);
    double best_proposal = 0;
    std::cout << std::setprecision(17)
              << "iteration\tbaseline_lb\tbest_proposal_lb\tproposal_lb\tmethod\tproposals\n";
    for (std::size_t state = 0; state < states; ++state) {
      int iteration;
      double baseline;
      input >> iteration >> baseline;
      for (auto& weight : weights) input >> weight;
      for (auto& bit : selected) { int value; input >> value; bit = value; }
      if (!input) throw std::invalid_argument("Truncated replay state");
      recovery.observe_adjusted(weights);
      if (state == 0 || iteration % frequency == 0 || state + 1 == states) {
        auto proposal = recovery.propose(selected);
        best_proposal = std::max(best_proposal, proposal.objective);
        std::cout << iteration << '\t' << baseline << '\t' << best_proposal << '\t'
                  << proposal.objective << '\t' << proposal.method << '\t' << proposal.proposals << '\n';
      }
    }
    std::string extra;
    if (input >> extra) throw std::invalid_argument("Extra data after replay states");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
