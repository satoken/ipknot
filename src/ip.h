/*
 * Copyright (C) 2012 Kengo Sato
 *
 * This file is part of DAFS.
 *
 * DAFS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DAFS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DAFS.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INC_IP_H__
#define __INC_IP_H__

#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

class IPimpl;
struct IPModel;
class IPInfeasible : public std::runtime_error {
public:
  explicit IPInfeasible(const std::string& message) : std::runtime_error(message) {}
};

class IP
{
public:
  typedef enum {MIN, MAX} DirType;
  typedef enum {FR, LO, UP, DB, FX} BoundType;

public:
  IP(DirType dir, int n_th, bool exact = false);
  static bool available();
  // Record the shared formulation without constructing a MIP backend.
  explicit IP(IPModel& model);
  ~IP();
  int make_variable(double coef);
  int make_variable(double coef, int lo, int hi);
  int make_continuous_variable(double coef, double lo = 0.0, double hi = 1.0);
  void add_objective_coefficient(int col, double coefficient);
  // Tag only NOE explanation/context/violation columns in a recorded model.
  void mark_noe_variable(int col);
  int make_constraint(BoundType bnd, double l, double u);
  void add_constraint(int row, int col, double val);
  void update();
  double solve();
  double get_value(int col) const;

private:
  IPimpl* impl_;
  IPModel* model_ = nullptr;
};

// A solver-independent formulation, used by constrained dual decomposition.
struct IPModel {
  struct Variable { double coefficient, lower, upper; bool integer; };
  struct Row {
    IP::BoundType bound;
    double lower, upper;
    std::vector<std::pair<int, double>> terms;
  };
  std::vector<Variable> variables;
  std::vector<int> noe_columns;
  std::vector<Row> rows;
  std::vector<double> solution;
  std::function<double()> optimize;
};

#endif  // __INC_IP_H__

// Local Variables:
// mode: C++
// End:
