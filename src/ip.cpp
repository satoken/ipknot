/*
 * Copyright (C) 2012 Kengo Sato
 *
 * This file is part of RactIP.
 *
 * RactIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RactIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RactIP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ip.h"
#include <vector>
#include <cassert>
#include <stdexcept>
#ifdef WITH_GLPK
#include <glpk.h>
#endif
#ifdef WITH_CPLEX
extern "C" {
#include <ilcplex/cplex.h>
};
#endif
#ifdef WITH_GUROBI
#include "gurobi_c++.h"
#endif
#ifdef WITH_SCIP
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#endif
#ifdef WITH_HIGHS
#include <Highs.h>
#endif

#include <cfloat>

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

  int make_variable(double coef, int lo, int hi)
  {
    int col = glp_add_cols(ip_, 1);
    glp_set_col_bnds(ip_, col, GLP_DB, lo, hi);
    glp_set_col_kind(ip_, col, GLP_IV);
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
    assert(row>=0);
    ia_.push_back(row);
    assert(col>=0);
    ja_.push_back(col);
    ar_.push_back(val);
  }

  void update() {}

  double solve()
  {
    glp_smcp smcp;
    glp_iocp iocp;
    glp_init_smcp(&smcp); smcp.msg_lev = GLP_MSG_ERR;
    glp_init_iocp(&iocp); iocp.msg_lev = GLP_MSG_ERR;
    glp_load_matrix(ip_, ia_.size()-1, &ia_[0], &ja_[0], &ar_[0]);
    glp_simplex(ip_, &smcp);
    glp_intopt(ip_, &iocp);
    return glp_mip_obj_val(ip_);
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
    env_->set(GRB_IntParam_OutputFlag, 0); // disable solver's outputs
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

  int make_variable(double coef, int lo, int hi)
  {
    vars_.push_back(model_->addVar(lo, hi, dir_*coef, GRB_INTEGER));
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

  double solve()
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
    return model_->get(GRB_DoubleAttr_ObjVal);
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
    : env_(NULL), lp_(NULL), dir_(dir)
  {
    int status;
    env_ = CPXopenCPLEX(&status);
    if (env_==NULL)
    {
      char errmsg[CPXMESSAGEBUFSIZE];
      CPXgeterrorstring (env_, status, errmsg);
      throw std::runtime_error(errmsg);
    }
    status = CPXsetintparam(env_, CPXPARAM_Threads, n_th);
  }

  ~IPimpl()
  {
    if (lp_) CPXfreeprob(env_, &lp_);
    if (env_) CPXcloseCPLEX(&env_);
  }

  int make_variable(double coef)
  {
    int col = vars_.size();
    vars_.push_back('B');
    coef_.push_back(coef);
    vlb_.push_back(0.0);
    vub_.push_back(1.0);
    return col;
  }

  int make_variable(double coef, int lo, int hi)
  {
    int col = vars_.size();
    vars_.push_back('I');
    coef_.push_back(coef);
    vlb_.push_back(lo);
    vub_.push_back(hi);
    return col;
  }
 
  int make_constraint(IP::BoundType bnd, double l, double u)
  {
    int row = bnd_.size();
    bnd_.resize(bnd_.size()+1);
    rhs_.resize(rhs_.size()+1);
    rngval_.resize(rngval_.size()+1);
    switch (bnd)
    {
      case IP::LO: bnd_[row]='G'; rhs_[row]=l; break;
      case IP::UP: bnd_[row]='L'; rhs_[row]=u; break;
      case IP::DB: bnd_[row]='R'; rhs_[row]=l; rngval_[row]=u-l; break;
      case IP::FX: bnd_[row]='E'; rhs_[row]=l; break;
      case IP::FR: bnd_[row]='R'; rhs_[row]=-DBL_MAX; rngval_[row]=DBL_MAX; break;
    }
    return row;
  }

  void add_constraint(int row, int col, double val)
  {
    m_[col].push_back(std::make_pair(row, val));
  }

  void update()
  {
    m_.resize(vars_.size());
  }

  double solve()
  {
    const int numcols = vars_.size();
    const int numrows = bnd_.size();

    int status;
    lp_ = CPXcreateprob(env_, &status, "");
    if (lp_==NULL) 
      throw std::runtime_error("failed to create LP");
    
    unsigned int n_nonzero=0;
    for (unsigned int i=0; i!=m_.size(); ++i) 
      n_nonzero += m_[i].size();
    std::vector<int> matbeg(numcols, 0);
    std::vector<int> matcnt(numcols, 0);
    std::vector<int> matind(n_nonzero);
    std::vector<double> matval(n_nonzero);
    for (unsigned int i=0, k=0; i!=m_.size(); ++i) 
    {
      matbeg[i] = i==0 ? 0 : matbeg[i-1]+matcnt[i-1];
      matcnt[i] = m_[i].size();
      for (unsigned int j=0; j!=m_[i].size(); ++j, ++k)
      {
        matind[k] = m_[i][j].first;
        matval[k] = m_[i][j].second;
      }
    }
    m_.clear();

    status = CPXcopylp(env_, lp_, numcols, numrows,
                        dir_==IP::MIN ? CPX_MIN : CPX_MAX,
                        &coef_[0], &rhs_[0], &bnd_[0], 
                        &matbeg[0], &matcnt[0], &matind[0], &matval[0],
                        &vlb_[0], &vub_[0], &rngval_[0] );
    vlb_.clear();
    vub_.clear();

    status = CPXcopyctype(env_, lp_, &vars_[0]);
    vars_.clear();

    CPXsetintparam(env_, CPXPARAM_MIP_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Barrier_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Tune_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Network_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Sifting_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Simplex_Display, 0);

    status = CPXmipopt(env_, lp_);
    double objval;
    status = CPXgetobjval(env_, lp_, &objval);
    res_cols_.resize(CPXgetnumcols(env_, lp_));
    status = CPXgetx(env_, lp_, &res_cols_[0], 0, res_cols_.size()-1);

    return objval;
  }

  double get_value(int col) const
  {
    return res_cols_[col];
  }

private:
  CPXENVptr env_;
  CPXLPptr lp_;
  IP::DirType dir_;
  std::vector<char> vars_;
  std::vector<double> coef_;
  std::vector<double> vlb_;
  std::vector<double> vub_;
  std::vector<char> bnd_;
  std::vector<double> rhs_;
  std::vector<double> rngval_;
  std::vector< std::vector< std::pair<int,double> > > m_;
  std::vector<double> res_cols_;
};
#endif  // WITH_CPLEX

#ifdef WITH_SCIP
class IPimpl
{
public:
  IPimpl(IP::DirType dir, int n_th)
    : scip_(nullptr), sol_(nullptr)
    //: ip_(NULL), ia_(1), ja_(1), ar_(1)
  {
    SCIPcreate(&scip_);
    SCIPincludeDefaultPlugins(scip_);
    SCIPcreateProbBasic(scip_, "IPknot");
    switch (dir)
    {
      case IP::MIN: SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE); break;
      case IP::MAX: SCIPsetObjsense(scip_, SCIP_OBJSENSE_MAXIMIZE); break;
    }
  }

  ~IPimpl()
  {
    for (auto& v: vars_)
      SCIPreleaseVar(scip_, &v);
    for (auto& c: cons_)
      SCIPreleaseCons(scip_, &c);
    if (scip_) SCIPfree(&scip_);
  }

  int make_variable(double coef, int lo=0, int hi=1)
  {
    int col = vars_.size();
    SCIP_VAR *var = nullptr;
    char buf[16];
    snprintf(buf, sizeof(buf), "var[%d]", col);
    SCIPcreateVarBasic(scip_,                // SCIP environment
                       &var,                 // reference to the variable
                       buf,                  // name of the variable
                       lo,                   // Lower bound of the variable
                       hi,                   // upper bound of the variable
                       coef,                 // Obj. coefficient. 
                       SCIP_VARTYPE_BINARY   // Binary variable
                      );
    SCIPaddVar(scip_, var);
    vars_.push_back(var);
    return col;
  }

  int make_constraint(IP::BoundType bnd, double l, double u)
  {
    int row = cons_.size();
    SCIP_CONS *cons = nullptr;
    char buf[16];
    snprintf(buf, sizeof(buf), "cons[%d]", row);
    double lhs = -SCIPinfinity(scip_); 
    double rhs =  SCIPinfinity(scip_);
    switch (bnd)
    {
      case IP::FR: break;
      case IP::LO: lhs = l; break; 
      case IP::UP: rhs = u; break;
      case IP::DB: lhs = l; rhs = u; break;
      case IP::FX: lhs = rhs = l; break;
    }
    SCIPcreateConsBasicLinear(scip_,        // SCIP pointer
                              &cons,        // reference to the constraint
                              buf,          // name of the constraint
                              0,            // How many variables are you adding  now
                              nullptr,      // an array of pointers to various variables
                              nullptr,      // an array of values of the coefficients of the corresponding vars
                              lhs,          // LHS of the constraint 
                              rhs           // RHS of the constraint
                             );
    SCIPaddCons(scip_, cons);
    cons_.push_back(cons);
    return row;
  }

  void add_constraint(int row, int col, double val)
  {
    assert(row>=0);
    assert(col>=0);
    SCIPaddCoefLinear(scip_, cons_[row], vars_[col], val);
  }

  void update()
  {
    for (auto c: cons_)
      SCIPaddCons(scip_, c);
  }


  double solve()
  {
    // SCIP_CALL((SCIPwriteOrigProblem(scip, "ipknot.lp", nullptr, FALSE)));
    SCIPsetIntParam(scip_, "display/verblevel", 0);   // We use SCIPsetIntParams to turn off the logging. 
    SCIPsolve(scip_);
    SCIP_STATUS soln_status = SCIPgetStatus(scip_); 
    sol_ = SCIPgetBestSol(scip_);
    return SCIPgetSolOrigObj(scip_, sol_);
  }

  double get_value(int col) const
  {
    return SCIPgetSolVal(scip_, sol_, vars_[col]);
  }

private:
  SCIP *scip_;
  std::vector<SCIP_VAR *> vars_;
  std::vector<SCIP_CONS *> cons_;
  SCIP_SOL *sol_;
};
#endif // WITH_SCIP

#ifdef WITH_HIGHS
class IPimpl
{
public:
  IPimpl(IP::DirType dir, int n_th)
    : model_()
  {
    highs_.setOptionValue("output_flag", false);
    highs_.setOptionValue("threads", n_th);
    switch (dir)
    {
      default:
      case IP::MAX: model_.lp_.sense_ = ObjSense::kMaximize; break;
      case IP::MIN: model_.lp_.sense_ = ObjSense::kMinimize; break;
    }
  }

  ~IPimpl()
  {
  }

  int make_variable(double coef, double l=0, double u=1)
  {
    int col = col_cost_.size();
    col_cost_.push_back(coef);
    col_lower_.push_back(l);
    col_upper_.push_back(u);
    m_.resize(col_cost_.size());
    return col;
  }

  int make_constraint(IP::BoundType bnd, double l, double u)
  {
    int row = row_lower_.size();
    switch (bnd)
    {
      case IP::LO: u=DBL_MAX; break;
      case IP::UP: l=-DBL_MAX; break;
      case IP::DB: break;
      case IP::FX: u=l; break;
      case IP::FR: l=-DBL_MAX; u=DBL_MAX; break;
    }
    row_lower_.push_back(l);
    row_upper_.push_back(u);
    return row;
  }

  void add_constraint(int row, int col, double val)
  {
    m_[col].emplace_back(row, val);
  }

  void update()
  {
    m_.resize(col_cost_.size());
  }

  double solve()
  {
    const int numcols = col_cost_.size();
    const int numrows = row_lower_.size();
    std::vector<HighsInt> start(numcols+1);
    start[0] = 0;
    for (auto i=1; i!=start.size(); i++)
      start[i] = start[i-1] + m_[i-1].size();
    auto non_zeros = start[start.size()-1];
    std::vector<HighsInt> index(non_zeros);
    std::vector<double> value(non_zeros);
    for (auto i=0, k=0; i!=numcols; i++)
    {
      for (auto j=0; j!=m_[i].size(); j++, k++)
      {
        index[k] = m_[i][j].first;
        value[k] = m_[i][j].second;
      }
    }
    m_.clear();

    model_.lp_.num_col_ = numcols;
    model_.lp_.num_row_ = numrows;
    model_.lp_.col_cost_ = col_cost_;
    model_.lp_.col_lower_ = col_lower_;
    model_.lp_.col_upper_ = col_upper_;
    model_.lp_.row_lower_ = row_lower_;
    model_.lp_.row_upper_ = row_upper_;
    model_.lp_.a_matrix_.start_ = start;
    model_.lp_.a_matrix_.index_ = index;
    model_.lp_.a_matrix_.value_ = value;

    // To indicate that variables must take integer values use the HighsLp::integrality vector.
    model_.lp_.integrality_.resize(numcols);
    for (int col=0; col < numcols; col++)
      model_.lp_.integrality_[col] = HighsVarType::kInteger;

    HighsStatus return_status;
  
    // Pass the model to HiGHS
    return_status = highs_.passModel(model_);
    assert(return_status==HighsStatus::kOk);

    // Get a const reference to the LP data in HiGHS
    // const HighsLp& lp = highs_.getLp();
  
    // Solve the model
    return_status = highs_.run();
    assert(return_status==HighsStatus::kOk);
  
    // Get the model status
    const HighsModelStatus& model_status = highs_.getModelStatus();
    assert(model_status==HighsModelStatus::kOptimal);

    const HighsInfo& info = highs_.getInfo();
    return info.objective_function_value;
  }

  double get_value(int col) const
  {
    // Get the solution values and basis
    const HighsSolution& solution = highs_.getSolution();
    return solution.col_value[col];
  }

private:
  IP::DirType dir_;
  HighsModel model_;
  Highs highs_;
  std::vector<double> col_cost_;
  std::vector<double> col_lower_;
  std::vector<double> col_upper_;
  std::vector<double> row_lower_;
  std::vector<double> row_upper_;
  std::vector< std::vector< std::pair<int,double> > > m_;
};
#endif


#if defined(WITH_GLPK) || defined(WITH_CPLEX) || defined(WITH_GUROBI) || defined(WITH_SCIP) || defined(WITH_HIGHS)

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
make_variable(double coef, int lo, int hi)
{
  return impl_->make_variable(coef, lo, hi);
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

double
IP::
solve()
{
  return impl_->solve();
}

double
IP::
get_value(int col) const
{
  return impl_->get_value(col);
}

#else

IP::
IP(DirType dir, int n_th)
{
  throw "no IP solver is linked.";
}

IP::
~IP()
{
  //throw "no IP solver is linked.";
}

int
IP::
make_variable(double coef)
{
  throw "no IP solver is linked.";
  return 0;
}

int
make_variable(double coef, int lo, int hi)
{
  throw "no IP solver is linked.";
  return 0;
}

int
IP::
make_constraint(BoundType bnd, double l, double u)
{
  throw "no IP solver is linked.";
  return 0;
}

void
IP::
add_constraint(int row, int col, double val)
{
  throw "no IP solver is linked.";
}

void
IP::
update()
{
  throw "no IP solver is linked.";
}

double
IP::
solve()
{
  throw "no IP solver is linked.";
}

double
IP::
get_value(int col) const
{
  throw "no IP solver is linked.";
  return 0;
}

#endif
