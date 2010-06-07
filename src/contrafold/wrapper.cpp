// 
#include <vector>
#include <string>
#include "InferenceEngine.hpp"
#include "ParameterManager.hpp"
#include "SStruct.hpp"
#include "Defaults.ipp"
#include "contrafold.h"

template < class T > 
struct CONTRAfold<T>::Impl
{
  Impl(bool canonical_only, int max_bp_dist);
  ~Impl() {};

  void SetParameters(const std::string& params);
  void SetConstraint(const std::string& paren);
  
  const T* ComputePosterior(const std::string& seq);
  const T* ComputePosterior(const std::string& seq, std::vector<T>& p);
  T ComputeInside(const std::string& seq);
  T ComputeLogPartitionCoefficient() const;
  T ComputeViterbi(const std::string& seq);

  void InitRand(unsigned int seed) { engine_.InitRand(seed); }
  void PrepareStochasticTraceback(const std::string& seq);
  std::vector<int> StochasticTraceback() const;

  ParameterManager<T> pm_;
  InferenceEngine<T> engine_;
  std::vector<T> w_;
  std::string paren_;
};

template < class T >
CONTRAfold<T>::
CONTRAfold(bool canonical_only, int max_bp_dist)
  : impl_(new Impl(canonical_only, max_bp_dist)), max_bp_dist_(max_bp_dist)
{
}

template < class T >
CONTRAfold<T>::
~CONTRAfold()
{
  delete impl_;
}

template < class T >
void
CONTRAfold<T>::
SetParameters(const std::string& params)
{
  impl_->SetParameters(params);
}

template < class T >
void
CONTRAfold<T>::
SetConstraint(const std::string& paren)
{
  impl_->SetConstraint(paren);
}

template < class T >
const T* 
CONTRAfold<T>::
ComputePosterior(const std::string& seq)
{
  return impl_->ComputePosterior(seq);
}

template < class T >
const T* 
CONTRAfold<T>::
ComputePosterior(const std::string& seq, std::vector<T>& p)
{
  return impl_->ComputePosterior(seq, p);
}

template < class T >
T
CONTRAfold<T>::
ComputeInside(const std::string& seq)
{
  return impl_->ComputeInside(seq);
}

template < class T >
T
CONTRAfold<T>::
ComputeLogPartitionCoefficient() const
{
  return impl_->ComputeLogPartitionCoefficient();
}

template < class T >
T
CONTRAfold<T>::
ComputeViterbi(const std::string& seq)
{
  return impl_->ComputeViterbi(seq);
}

// static
template < class T >
void
CONTRAfold<T>::
init_rand(unsigned long seed)
{
  impl_->InitRand(seed);
}

template < class T >
void
CONTRAfold<T>::
PrepareStochasticTraceback(const std::string& seq)
{
  return impl_->PrepareStochasticTraceback(seq);
}

template < class T >
std::vector<int>
CONTRAfold<T>::
StochasticTraceback() const
{
  return impl_->StochasticTraceback();
}

template < class T > 
CONTRAfold<T>::Impl::
Impl(bool canonical_only, int max_bp_dist)
  : pm_(), engine_(!canonical_only, max_bp_dist), paren_()
{
  engine_.RegisterParameters(pm_);
  if (canonical_only)
    w_ = GetDefaultComplementaryValues<float>();
  else
    w_ = GetDefaultNoncomplementaryValues<float>();
}

template < class T >
void
CONTRAfold<T>::Impl::
SetParameters(const std::string& params)
{
  pm_.ReadFromFile(params, w_);
}

template < class T >
void
CONTRAfold<T>::Impl::
SetConstraint(const std::string& paren)
{
  paren_ = paren;
}

template < class T > 
const T*
CONTRAfold<T>::Impl::
ComputePosterior(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  engine_.ComputeOutside();
  engine_.ComputePosterior();
  delete s;
  return engine_.GetPosterior(0.0);
}

template < class T > 
const T*
CONTRAfold<T>::Impl::
ComputePosterior(const std::string& seq, std::vector<T>& p)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  engine_.ComputeOutside();
  engine_.ComputePosterior();
  delete s;
  return engine_.GetPosterior(0.0, p);
}

template < class T > 
T
CONTRAfold<T>::Impl::
ComputeInside(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  delete s;
  return engine_.ComputeLogPartitionCoefficient();
}

template < class T > 
T
CONTRAfold<T>::Impl::
ComputeLogPartitionCoefficient() const
{
  return engine_.ComputeLogPartitionCoefficient();
}

template < class T > 
T
CONTRAfold<T>::Impl::
ComputeViterbi(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeViterbi();

  return engine_.GetViterbiScore();
}

template < class T > 
void
CONTRAfold<T>::Impl::
PrepareStochasticTraceback(const std::string& seq)
{
  SStruct* s=NULL;
  if (paren_.empty())
    s = new SStruct("unknown", seq);
  else
    s = new SStruct("unknown", seq, paren_);
  engine_.LoadSequence(*s);
  if (!paren_.empty())
    engine_.UseConstraints(s->GetMapping());
  engine_.LoadValues(w_);
  engine_.ComputeInside();
  //engine_.ComputeOutside();
  delete s;
}

template < class T > 
std::vector<int>
CONTRAfold<T>::Impl::
StochasticTraceback() const
{
  return engine_.PredictPairingsStochasticTraceback();
}

template < class T > 
struct CONTRAfoldM<T>::Impl
{
  Impl(bool canonical_only, int max_bp_dist);
  ~Impl()
  {
    for (unsigned int i=0; i!=en_.size(); ++i) if (en_[i]) delete en_[i];
  };

  void SetParameters(const std::string& params);
  //void SetConstraint(const std::string& paren);
  
  //const T* ComputePosterior(const std::string& seq);
  //const T* ComputePosterior(const std::string& seq, std::vector<T>& p);
  //T ComputeInside(const std::string& seq);
  //T ComputeLogPartitionCoefficient() const;
  //T ComputeViterbi(const std::string& seq);

  void InitRand(unsigned int seed) { engine_.InitRand(seed); }
  void PrepareStochasticTraceback(const std::vector<std::string>& aln);
  std::vector<int> StochasticTraceback() const;

  ParameterManager<T> pm_;
  InferenceEngine<T> engine_;
  std::vector<T> w_;
  std::string paren_;
  int max_bp_dist_;
  bool canonical_only_;
  std::vector< std::vector<unsigned int> > idx_;
  std::vector< std::vector<int> > rev_;
  std::vector<InferenceEngine<T>*> en_;  
};

template < class T >
CONTRAfoldM<T>::
CONTRAfoldM(bool canonical_only, int max_bp_dist)
  : impl_(new Impl(canonical_only, max_bp_dist)), max_bp_dist_(max_bp_dist)
{
}

template < class T >
CONTRAfoldM<T>::
~CONTRAfoldM()
{
  delete impl_;
}

template < class T >
void
CONTRAfoldM<T>::
SetParameters(const std::string& params)
{
  impl_->SetParameters(params);
}

// static
template < class T >
void
CONTRAfoldM<T>::
init_rand(unsigned long seed)
{
  impl_->InitRand(seed);
}

template < class T >
void
CONTRAfoldM<T>::
PrepareStochasticTraceback(const std::vector<std::string>& aln)
{
  return impl_->PrepareStochasticTraceback(aln);
}

template < class T >
std::vector<int>
CONTRAfoldM<T>::
StochasticTraceback() const
{
  return impl_->StochasticTraceback();
}

template < class T > 
CONTRAfoldM<T>::Impl::
Impl(bool canonical_only, int max_bp_dist)
  : pm_(), engine_(!canonical_only, max_bp_dist), paren_(),
    max_bp_dist_(max_bp_dist), canonical_only_(canonical_only)
{
  engine_.RegisterParameters(pm_);
  if (canonical_only)
    w_ = GetDefaultComplementaryValues<float>();
  else
    w_ = GetDefaultNoncomplementaryValues<float>();
}

template < class T >
void
CONTRAfoldM<T>::Impl::
SetParameters(const std::string& params)
{
  pm_.ReadFromFile(params, w_);
}

template < class T > 
void
CONTRAfoldM<T>::Impl::
PrepareStochasticTraceback(const std::vector<std::string>& aln)
{
  for (unsigned int i=0; i!=en_.size(); ++i) if (en_[i]) delete en_[i];
  en_.resize(aln.size(), NULL);
  idx_.clear(); idx_.resize(aln.size());
  rev_.clear(); rev_.resize(aln.size());
  for (unsigned int i=0; i!=aln.size(); ++i)
  {
    std::string seq;
    idx_[i].resize(aln[i].size()+1);
    idx_[i][0] = -1u;
    rev_[i].resize(aln[i].size()+1);
    for (unsigned int j=0, k=0; j!=aln[i].size(); ++j)
    {
      if (aln[i][j]=='-')
      {
        idx_[i][j+1] =  -1u;
      }
      else
      {
        idx_[i][j+1] = k+1;
        rev_[i][k+1] = j+1;
        seq.push_back(aln[i][j]);
        k++;
      }
    }

    SStruct s("unknown", seq);
    en_[i] = new InferenceEngine<T>(!canonical_only_, max_bp_dist_);
    en_[i]->RegisterParameters(pm_);
    en_[i]->LoadSequence(s);
    en_[i]->LoadValues(w_);
    en_[i]->ComputeInside();
    en_[i]->ComputeOutside();
  }
}

template < class T > 
std::vector<int>
CONTRAfoldM<T>::Impl::
StochasticTraceback() const
{
  return engine_.PredictPairingsStochasticTracebackM(idx_, rev_, en_);
}

// instantiation
template
class CONTRAfold<float>;

template
class CONTRAfoldM<float>;
