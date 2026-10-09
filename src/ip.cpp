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
#include <memory>
#include <cassert>
#include <stdexcept>
#include <utility>

#include <cstdlib>
#include <cstdio>
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

  void configure_exact() {} // glp_iocp defaults to zero MIP gap

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

  int make_continuous_variable(double coef, double lo, double hi)
  {
    int col = glp_add_cols(ip_, 1);
    glp_set_col_bnds(ip_, col, GLP_DB, lo, hi);
    glp_set_col_kind(ip_, col, GLP_CV);
    glp_set_obj_coef(ip_, col, coef);
    return col;
  }

  void add_objective_coefficient(int col, double coefficient)
  {
    glp_set_obj_coef(ip_, col, glp_get_obj_coef(ip_, col) + coefficient);
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
    const int lp_status = glp_simplex(ip_, &smcp);
    if (glp_get_status(ip_) == GLP_NOFEAS) throw IPInfeasible("GLPK model is infeasible");
    const int mip_status = glp_intopt(ip_, &iocp);
    if (glp_mip_status(ip_) == GLP_NOFEAS) throw IPInfeasible("GLPK model is infeasible");
    if (lp_status != 0 || mip_status != 0 || glp_mip_status(ip_) != GLP_OPT) {
      throw std::runtime_error("GLPK failed to find an optimal solution");
    }
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

  void configure_exact() {
    model_->set(GRB_DoubleParam_MIPGap, 0.0);
    model_->set(GRB_DoubleParam_MIPGapAbs, 0.0);
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

  int make_continuous_variable(double coef, double lo, double hi)
  {
    vars_.push_back(model_->addVar(lo, hi, dir_*coef, GRB_CONTINUOUS));
    return vars_.size()-1;
  }

  void add_objective_coefficient(int col, double coefficient)
  {
    vars_[col].set(GRB_DoubleAttr_Obj, vars_[col].get(GRB_DoubleAttr_Obj) + dir_ * coefficient);
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
    if (model_->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
      throw IPInfeasible("Gurobi model is infeasible");
    if (model_->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
      throw std::runtime_error("Gurobi failed to find an optimal solution");
    }
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

  void configure_exact() {
    if (CPXsetdblparam(env_, CPXPARAM_MIP_Tolerances_MIPGap, 0.0) ||
        CPXsetdblparam(env_, CPXPARAM_MIP_Tolerances_AbsMIPGap, 0.0))
      throw std::runtime_error("Cannot set CPLEX exact MIP tolerances");
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
    // NMR constraints may introduce auxiliary variables after the initial
    // update().  Keep the sparse column storage in step with vars_ so that
    // add_constraint() can address those columns immediately.
    m_.emplace_back();
    return col;
  }

  int make_variable(double coef, int lo, int hi)
  {
    int col = vars_.size();
    vars_.push_back('I');
    coef_.push_back(coef);
    vlb_.push_back(lo);
    vub_.push_back(hi);
    m_.emplace_back();
    return col;
  }
 
  int make_continuous_variable(double coef, double lo, double hi)
  {
    int col = vars_.size();
    vars_.push_back('C');
    coef_.push_back(coef);
    vlb_.push_back(lo);
    vub_.push_back(hi);
    m_.emplace_back();
    return col;
  }

  void add_objective_coefficient(int col, double coefficient)
  {
    coef_[col] += coefficient;
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
                        coef_.data(), rhs_.data(), bnd_.data(),
                        matbeg.data(), matcnt.data(), matind.data(), matval.data(),
                        vlb_.data(), vub_.data(), rngval_.data() );
    if (status != 0)
      throw std::runtime_error("CPLEX failed to copy the optimization model");
    vlb_.clear();
    vub_.clear();

    status = CPXcopyctype(env_, lp_, vars_.data());
    if (status != 0)
      throw std::runtime_error("CPLEX failed to copy variable types");
    vars_.clear();

    CPXsetintparam(env_, CPXPARAM_MIP_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Barrier_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Tune_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Network_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Sifting_Display, 0);
    CPXsetintparam(env_, CPXPARAM_Simplex_Display, 0);

    status = CPXmipopt(env_, lp_);
    if (status != 0) {
      throw std::runtime_error("CPLEX failed while optimizing the model");
    }
    const int solution_status = CPXgetstat(env_, lp_);
    if (solution_status == CPXMIP_INFEASIBLE) throw IPInfeasible("CPLEX model is infeasible");
    if (solution_status != CPXMIP_OPTIMAL &&
        solution_status != CPXMIP_OPTIMAL_TOL) {
      char status_message[CPXMESSAGEBUFSIZE];
      CPXgetstatstring(env_, solution_status, status_message);
      throw std::runtime_error(
          std::string("CPLEX failed to find an optimal solution: ") +
          status_message);
    }
    double objval;
    status = CPXgetobjval(env_, lp_, &objval);
    res_cols_.resize(CPXgetnumcols(env_, lp_));
    status = CPXgetx(env_, lp_, res_cols_.data(), 0, res_cols_.size()-1);

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

  void configure_exact() {
    if (SCIPsetRealParam(scip_, "limits/gap", 0.0) != SCIP_OKAY ||
        SCIPsetRealParam(scip_, "limits/absgap", 0.0) != SCIP_OKAY)
      throw std::runtime_error("Cannot set SCIP exact MIP tolerances");
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
    const SCIP_VARTYPE variable_type =
        (lo == 0 && hi == 1) ? SCIP_VARTYPE_BINARY : SCIP_VARTYPE_INTEGER;
    char buf[16];
    snprintf(buf, sizeof(buf), "var[%d]", col);
    SCIPcreateVarBasic(scip_,                // SCIP environment
                       &var,                 // reference to the variable
                       buf,                  // name of the variable
                       lo,                   // Lower bound of the variable
                       hi,                   // upper bound of the variable
                       coef,                 // Obj. coefficient. 
                       variable_type         // Binary or bounded integer
                      );
    SCIPaddVar(scip_, var);
    vars_.push_back(var);
    return col;
  }

  int make_continuous_variable(double coef, double lo, double hi)
  {
    int col = vars_.size();
    SCIP_VAR* var = nullptr;
    char buf[32];
    snprintf(buf, sizeof(buf), "continuous[%d]", col);
    SCIPcreateVarBasic(scip_, &var, buf, lo, hi, coef, SCIP_VARTYPE_CONTINUOUS);
    SCIPaddVar(scip_, var);
    vars_.push_back(var);
    return col;
  }

  void add_objective_coefficient(int col, double coefficient)
  {
    SCIPchgVarObj(scip_, vars_[col], SCIPvarGetObj(vars_[col]) + coefficient);
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
    if (soln_status == SCIP_STATUS_INFEASIBLE) throw IPInfeasible("SCIP model is infeasible");
    sol_ = SCIPgetBestSol(scip_);
    if (soln_status != SCIP_STATUS_OPTIMAL || sol_ == nullptr) {
      throw std::runtime_error("SCIP failed to find an optimal solution");
    }
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
    highs_.setOptionValue("output_flag", std::getenv("IPKNOT_HIGHS_LOG") != nullptr);
    highs_.setOptionValue("threads", n_th);
    switch (dir)
    {
      default:
      case IP::MAX: model_.lp_.sense_ = ObjSense::kMaximize; break;
      case IP::MIN: model_.lp_.sense_ = ObjSense::kMinimize; break;
    }
  }

  void configure_exact() {
    if (highs_.setOptionValue("mip_rel_gap", 0.0) != HighsStatus::kOk ||
        highs_.setOptionValue("mip_abs_gap", 0.0) != HighsStatus::kOk)
      throw std::runtime_error("Cannot set HiGHS exact MIP tolerances");
  }

  ~IPimpl()
  {
  }

  int make_variable(double coef, double l=0, double u=1)
  {
    int col = col_cost_.size();
    col_cost_.push_back(coef);
    integrality_.push_back(HighsVarType::kInteger);
    col_lower_.push_back(l);
    col_upper_.push_back(u);
    m_.resize(col_cost_.size());
    return col;
  }

  int make_continuous_variable(double coef, double lo, double hi)
  {
    int col = make_variable(coef, lo, hi);
    integrality_[col] = HighsVarType::kContinuous;
    return col;
  }

  void add_objective_coefficient(int col, double coefficient)
  {
    col_cost_[col] += coefficient;
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
    std::vector<HighsInt> index;
    std::vector<double> value;
    index.reserve(non_zeros);
    value.reserve(non_zeros);
    for (auto i=0; i!=numcols; i++)
    {
      for (const auto& entry : m_[i])
      {
        index.push_back(entry.first);
        value.push_back(entry.second);
      }
      // Release each source column as it is transferred. Large NOE matrices
      // need not keep the complete builder and CSC buffers resident together.
      std::vector<std::pair<int, double>>().swap(m_[i]);
    }
    m_.clear();

    model_.lp_.num_col_ = numcols;
    model_.lp_.num_row_ = numrows;
    model_.lp_.col_cost_ = std::move(col_cost_);
    model_.lp_.col_lower_ = std::move(col_lower_);
    model_.lp_.col_upper_ = std::move(col_upper_);
    model_.lp_.row_lower_ = std::move(row_lower_);
    model_.lp_.row_upper_ = std::move(row_upper_);
    model_.lp_.a_matrix_.start_ = std::move(start);
    model_.lp_.a_matrix_.index_ = std::move(index);
    model_.lp_.a_matrix_.value_ = std::move(value);

    // To indicate that variables must take integer values use the HighsLp::integrality vector.
    model_.lp_.integrality_ = integrality_;

    HighsStatus return_status;
  
    // Pass the model to HiGHS
    return_status = highs_.passModel(std::move(model_));
    assert(return_status==HighsStatus::kOk);

    // Get a const reference to the LP data in HiGHS
    // const HighsLp& lp = highs_.getLp();
  
    // Solve the model
    return_status = highs_.run();
    if (return_status != HighsStatus::kOk) {
      throw std::runtime_error("HiGHS solver returned non-OK status");
    }

    // Get the model status
    const HighsModelStatus& model_status = highs_.getModelStatus();
    if (model_status == HighsModelStatus::kInfeasible) throw IPInfeasible("HiGHS model is infeasible");
    if (model_status != HighsModelStatus::kOptimal) {
      std::string status_str;
      switch (model_status) {
        case HighsModelStatus::kNotset: status_str = "Not set"; break;
        case HighsModelStatus::kLoadError: status_str = "Load error"; break;
        case HighsModelStatus::kModelError: status_str = "Model error"; break;
        case HighsModelStatus::kPresolveError: status_str = "Presolve error"; break;
        case HighsModelStatus::kSolveError: status_str = "Solve error"; break;
        case HighsModelStatus::kPostsolveError: status_str = "Postsolve error"; break;
        case HighsModelStatus::kModelEmpty: status_str = "Model empty"; break;
        case HighsModelStatus::kInfeasible: status_str = "Infeasible"; break;
        case HighsModelStatus::kUnboundedOrInfeasible: status_str = "Unbounded or infeasible"; break;
        case HighsModelStatus::kUnbounded: status_str = "Unbounded"; break;
        case HighsModelStatus::kObjectiveBound: status_str = "Objective bound"; break;
        case HighsModelStatus::kObjectiveTarget: status_str = "Objective target"; break;
        case HighsModelStatus::kTimeLimit: status_str = "Time limit"; break;
        case HighsModelStatus::kIterationLimit: status_str = "Iteration limit"; break;
        case HighsModelStatus::kUnknown: status_str = "Unknown"; break;
        default: status_str = "Other"; break;
      }
      throw std::runtime_error("HiGHS failed to find an optimal solution: " + status_str);
    }

    const HighsInfo& info = highs_.getInfo();
    if (std::getenv("IPKNOT_HIGHS_DIAGNOSTICS"))
      std::fprintf(stderr, "IP solver statistics: nodes=%lld lp_iterations=%d cols=%d rows=%d gap=%.12g\n",
          static_cast<long long>(info.mip_node_count), static_cast<int>(info.simplex_iteration_count),
          numcols, numrows, info.mip_gap);
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
  std::vector<HighsVarType> integrality_;
  std::vector<double> col_lower_;
  std::vector<double> col_upper_;
  std::vector<double> row_lower_;
  std::vector<double> row_upper_;
  std::vector< std::vector< std::pair<int,double> > > m_;
};
#endif


#if !defined(WITH_GLPK) && !defined(WITH_CPLEX) && !defined(WITH_GUROBI) && !defined(WITH_SCIP) && !defined(WITH_HIGHS)
class IPimpl {
  [[noreturn]] static void unavailable() {
    throw std::runtime_error("No ILP solver is linked; use --decoder dd");
  }
public:
  IPimpl(IP::DirType, int) { unavailable(); }
  void configure_exact() { unavailable(); }
  int make_variable(double) { unavailable(); }
  int make_variable(double, int, int) { unavailable(); }
  int make_continuous_variable(double, double, double) { unavailable(); }
  void add_objective_coefficient(int, double) { unavailable(); }
  int make_constraint(IP::BoundType, double, double) { unavailable(); }
  void add_constraint(int, int, double) { unavailable(); }
  void update() { unavailable(); }
  double solve() { unavailable(); }
  double get_value(int) const { unavailable(); }
};
#endif

IP::IP(DirType dir, int n_th, bool exact) : impl_(nullptr) {
  auto impl = std::make_unique<IPimpl>(dir, n_th);
  if (exact) impl->configure_exact();
  impl_ = impl.release();
}
bool IP::available() {
#if defined(WITH_GLPK) || defined(WITH_CPLEX) || defined(WITH_GUROBI) || defined(WITH_SCIP) || defined(WITH_HIGHS)
  return true;
#else
  return false;
#endif
}
void IP::mark_noe_variable(int col) {
  if (model_) {
    if (col < 0 || col >= static_cast<int>(model_->variables.size()))
      throw std::invalid_argument("Invalid NOE variable tag");
    model_->noe_columns.push_back(col);
  }
}
IP::IP(IPModel& model) : impl_(nullptr), model_(&model) {}
IP::~IP() { delete impl_; }

int IP::make_variable(double coef) {
  if (!model_) return impl_->make_variable(coef);
  return make_variable(coef, 0, 1);
}
int IP::make_variable(double coef, int lo, int hi) {
  if (!model_) return impl_->make_variable(coef, lo, hi);
  model_->variables.push_back({coef, double(lo), double(hi), true});
  return static_cast<int>(model_->variables.size()) - 1;
}
int IP::make_continuous_variable(double coef, double lo, double hi) {
  if (!model_) return impl_->make_continuous_variable(coef, lo, hi);
  model_->variables.push_back({coef, lo, hi, false});
  return static_cast<int>(model_->variables.size()) - 1;
}
void IP::add_objective_coefficient(int col, double coefficient) {
  if (model_) model_->variables.at(col).coefficient += coefficient;
  else impl_->add_objective_coefficient(col, coefficient);
}
int IP::make_constraint(BoundType bnd, double l, double u) {
  if (!model_) return impl_->make_constraint(bnd, l, u);
  model_->rows.push_back({bnd, l, u, {}});
  return static_cast<int>(model_->rows.size()) - 1;
}
void IP::add_constraint(int row, int col, double val) {
  if (model_) model_->rows.at(row).terms.emplace_back(col, val);
  else impl_->add_constraint(row, col, val);
}
void IP::update() { if (!model_) impl_->update(); }
double IP::solve() {
  if (model_) {
    if (!model_->optimize) throw std::logic_error("Missing recorded-model decoder");
    return model_->optimize();
  }
  return impl_->solve();
}
double IP::get_value(int col) const {
  return model_ ? model_->solution.at(col) : impl_->get_value(col);
}
