#include "sparse_posterior.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>

void require(bool v) { if (!v) throw std::runtime_error("Sparse posterior order or float bits changed"); }
int main() {
  using Matrix = SparsePosteriorAccumulator::Matrix;
  Matrix actual(1025), expected(1025);
  std::mt19937 random(271828);
  for (unsigned i=1;i<30;++i) for (unsigned j=0;j<60;++j)
    actual[i].push_back({random()%1024+1,static_cast<float>(random()%100)/11});
  expected=actual;
  SparsePosteriorAccumulator accumulate(actual);
  for (int step=0;step<5000;++step) {
    unsigned row=random()%1024+1, partner=random()%1024+1;
    float value=(static_cast<int>(random()%201)-100)/17.f;
    auto found=std::find_if(expected[row].begin(),expected[row].end(),
        [&](const auto& pair) { return pair.first==partner; });
    if(found==expected[row].end()) expected[row].push_back({partner,value});
    else found->second+=value;
    accumulate.add(row,partner,value);
  }
  for(unsigned row=0;row<actual.size();++row) {
    require(actual[row].size()==expected[row].size());
    for(unsigned k=0;k<actual[row].size();++k) {
      require(actual[row][k].first==expected[row][k].first);
      std::uint32_t a,b;
      std::memcpy(&a,&actual[row][k].second,sizeof(a));
      std::memcpy(&b,&expected[row][k].second,sizeof(b));
      require(a==b);
    }
  }
  Matrix concentrated(20002);
  SparsePosteriorAccumulator long_row(concentrated);
  for(unsigned j=2;j<20002;++j) long_row.add(1,j,.25f);
  for(unsigned j=20002;j-- >2;) long_row.add(1,j,.5f);
  require(concentrated[1].size()==20000);
  for(unsigned k=0;k<20000;++k) require(concentrated[1][k].first==k+2 && concentrated[1][k].second==.75f);
  std::cout << "Sparse posterior merge preserves row order and float bits, including concentrated rows\n";
}
