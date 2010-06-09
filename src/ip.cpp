/*
 * $Id:$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of IPknot.
 *
 * IPknot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * IPknot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IPknot.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"
#include "ip.h"
#include <vector>
#ifdef WITH_GLPK
#include <glpk.h>
#endif
#ifdef WITH_CPLEX
#include <ilcplex/ilocplex.h>
#endif
#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif

#ifdef WITH_GLPK
class IPimpl
{
public:
  IPimpl(IP::DirType dir, int n_th)
    : ip_(NULL), ia_(1), ja_(1), ar_(1)
  {
    ip_ = glp_create_prob();
    switch (dir)
    {
      case IP::MIN: glp_set_obj_dir(ip_, GLP_MIN); break;
      case IP::MAX: glp_set_obj_dir(ip_, GLP_MAX); break;
    }
  }

  ~IPimpl()
  {
    glp_delete_prob(ip_);
  }

  int make_variable(double coef)
  {
    int col = glp_add_cols(ip_, 1);
    glp_set_col_bnds(ip_, col, GLP_DB, 0, 1);
    glp_set_col_kind(ip_, col, GLP_BV);
    glp_set_obj_coef(ip_, col, coef);
    return col;
  }

  int make_constraint(IP::BoundType bnd, double l, double u)
  {
    int row = glp_add_rows(ip_, 1);
    switch (bnd)
    {
      case IP::FR: glp_set_row_bnds(ip_, row, GLP_FR, l, u); break;
      case IP::LO: glp_set_row_bnds(ip_, row, GLP_LO, l, u); break;
      case IP::UP: glp_set_row_bnds(ip_, row, GLP_UP, l, u); break;
      case IP::DB: glp_set_row_bnds(ip_, row, GLP_DB, l, u); break;
      case IP::FX: glp_set_row_bnds(ip_, row, GLP_FX, l, u); break;
    }
    return row;
  }

  void add_constraint(int row, int col, double val)
  {
    ia_.push_back(row);
    ja_.push_back(col);
    ar_.push_back(val);
  }

  void update() {}

  void solve()
  {
    glp_smcp smcp;
    glp_iocp iocp;
    glp_init_smcp(&smcp); smcp.msg_lev = GLP_MSG_ERR;
    glp_init_iocp(&iocp); iocp.msg_lev = GLP_MSG_ERR;
    glp_load_matrix(ip_, ia_.size()-1, &ia_[0], &ja_[0], &ar_[0]);
    glp_simplex(ip_, &smcp);
    glp_intopt(ip_, &iocp);
  }

  double get_value(int col) const
  {
    return glp_mip_col_val(ip_, col);
  }

private:
  glp_prob *ip_;
  std::vector<int> ia_;
  std::vector<int> ja_;
  std::vector<double> ar_;
};
#endif

#ifdef WITH_GUROBI
class IPimpl
{
public:
  IPimpl(IP::DirType dir, int n_th)
    : env_(NULL), model_(NULL), dir_(dir==IP::MIN ? +1 : -1)
  {
    env_ = new GRBEnv;
    env_->set(GRB_IntParam_Threads, n_th); // # of threads
    model_ = new GRBModel(*env_);
  }

  ~IPimpl()
  {
    delete model_;
    delete env_;
  }

  int make_variable(double coef)
  {
    vars_.push_back(model_->addVar(0, 1, dir_*coef, GRB_BINARY));
    return vars_.size()-1;
  }

  int make_constraint(IP::BoundType bnd, double l, double u)
  {
    bnd_.push_back(bnd);
    l_.push_back(l);
    u_.push_back(u);
    m_.resize(m_.size()+1);
    return m_.size()-1;
  }

  void add_constraint(int row, int col, double val)
  {
    m_[row].push_back(std::make_pair(col, val));
  }

  void update()
  {
    model_->update();
  }

  void solve()
  {
    for (unsigned int i=0; i!=m_.size(); ++i)
    {
      GRBLinExpr c;
      for (unsigned int j=0; j!=m_[i].size(); ++j)
        c += vars_[m_[i][j].first] * m_[i][j].second;
      switch (bnd_[i])
      {
        case IP::LO: model_->addConstr(c >= l_[i]); break;
        case IP::UP: model_->addConstr(c <= u_[i]); break;
        case IP::DB: model_->addConstr(c >= l_[i]); model_->addConstr(c <= u_[i]); break;
        case IP::FX: model_->addConstr(c == l_[i]); break;
      }
    }
    bnd_.clear();
    l_.clear();
    u_.clear();
    m_.clear();
    model_->optimize();
  }

  double get_value(int col) const
  {
    return vars_[col].get(GRB_DoubleAttr_X);
  }

private:
  GRBEnv* env_;
  GRBModel* model_;
  int dir_;

  std::vector<GRBVar> vars_;
  std::vector< std::vector< std::pair<int,double> > > m_;
  std::vector<int> bnd_;
  std::vector<double> l_;
  std::vector<double> u_;
};
#endif  // WITH_GUROBI

#ifdef WITH_CPLEX
class IPimpl
{
public:
  IPimpl(IP::DirType dir, int n_th)
    : env_(), model_(env_), obj_(env_), vars_(env_), cplex_(NULL), dir_(dir), n_th_(n_th)
  {
  }

  ~IPimpl()
  {
    delete cplex_;
  }

  int make_variable(double coef)
  {
    int col = vars_.getSize();
    vars_.add(IloBoolVar(env_));
    obj_ += coef * vars_[col];
    return col;
  }

  int make_constraint(IP::BoundType bnd, double l, double u)
  {
    bnd_.push_back(bnd);
    l_.push_back(l);
    u_.push_back(u);
    m_.resize(m_.size()+1);
    return m_.size()-1;
  }

  void add_constraint(int row, int col, double val)
  {
    m_[row].push_back(std::make_pair(col, val));
  }

  void update()
  {
    switch (dir_)
    {
      case IP::MAX: model_.add(IloMaximize(env_, obj_)); break;
      case IP::MIN: model_.add(IloMinimize(env_, obj_)); break;
    }
  }

  void solve()
  {
    for (unsigned int i=0; i!=m_.size(); ++i)
    {
      IloExpr c(env_);
      for (unsigned int j=0; j!=m_[i].size(); ++j)
        c += vars_[m_[i][j].first] * m_[i][j].second;
      switch (bnd_[i])
      {
        case IP::LO: model_.add(c >= l_[i]); break;
        case IP::UP: model_.add(c <= u_[i]); break;
        case IP::DB: model_.add(c >= l_[i]); model_.add(c <= u_[i]); break;
        case IP::FX: model_.add(c == l_[i]); break;
      }
    }
    bnd_.clear();
    l_.clear();
    u_.clear();
    m_.clear();

    cplex_ = new IloCplex(model_);
    cplex_->solve();
  }

  double get_value(int col) const
  {
    return cplex_->getValue(vars_[col]);
  }

private:
  IloEnv env_;
  IloModel model_;
  IloExpr obj_;
  IloBoolVarArray vars_;
  IloCplex* cplex_;
  std::vector< std::vector< std::pair<int,double> > > m_;
  std::vector<int> bnd_;
  std::vector<double> l_;
  std::vector<double> u_;
  IP::DirType dir_;
  int n_th_;
};
#endif  // WITH_CPLEX

IP::
IP(DirType dir, int n_th)
  : impl_(new IPimpl(dir, n_th))
{
}

IP::
~IP()
{
  delete impl_;
}

int
IP::
make_variable(double coef)
{
  return impl_->make_variable(coef);
}

int
IP::
make_constraint(BoundType bnd, double l, double u)
{
  return impl_->make_constraint(bnd, l, u);
}

void
IP::
add_constraint(int row, int col, double val)
{
  return impl_->add_constraint(row, col, val);
}

void
IP::
update()
{
  impl_->update();
}

void
IP::
solve()
{
  impl_->solve();
}

double
IP::
get_value(int col) const
{
  return impl_->get_value(col);
}
