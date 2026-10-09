// Exercise real solver matrix handoff with many coefficients and fixed integer
// variables. Optional larger arguments also provide a reproducible RSS test.
#include "ip.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  const int columns = argc > 1 ? std::atoi(argv[1]) : 256;
  const int rows = argc > 2 ? std::atoi(argv[2]) : 1024;
  const int degree = argc > 3 ? std::atoi(argv[3]) : 8;
  if (columns <= 0 || rows <= 0 || degree <= 0 || degree > rows) return 2;
  std::vector<double> expected(rows, 0);
  auto row_index = [rows](int column, int element) {
    return static_cast<int>((17LL * column + element) % rows);
  };
  auto coefficient = [](int column, int element) { return (column / 7 + element) % 2 ? -1.0 : 1.0; };
  for (int column = 0; column < columns; ++column)
    for (int element = 0; element < degree; ++element)
      expected[row_index(column, element)] += coefficient(column, element) * (column % 2);
  IP ip(IP::MAX, 1);
  for (int column = 0; column < columns; ++column)
    assert(ip.make_variable(1, column % 2, column % 2) == column);
  for (int row = 0; row < rows; ++row)
    assert(ip.make_constraint(IP::FX, expected[row], expected[row]) == row);
  for (int column = 0; column < columns; ++column)
    for (int element = 0; element < degree; ++element)
      ip.add_constraint(row_index(column, element), column, coefficient(column, element));
  assert(ip.solve() == columns / 2);
  for (int column = 0; column < columns; ++column) assert(ip.get_value(column) == column % 2);
  std::cout << "Verified " << columns << " columns, " << rows << " rows, "
            << 1LL * columns * degree << " coefficients\n";
}
