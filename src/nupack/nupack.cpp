// $Id:$

#include <cassert>
#include <vector>

template < class PF_TYPE >
class Nupack<PF_TYPE>::DPTable2
{
public:
  DPTable2() : V_(), N_(0)
  {
  }

  void resize(int n)
  {
    N_=n;
    V_.resize(N_*(N_+1)/2+(N_+1));
  }

  void fill(const PF_TYPE& v)
  {
    std::fill(V_.begin(), V_.end(), v);
  }

  PF_TYPE& operator()(int i, int j)
  {
    return V_[index(i, j)];
  }

  const PF_TYPE& operator()() const
  {
    return V_[index(i, j)];
  }

private:
  int index(int i, int j) const
  {
    assert(i<=j);
    assert(j<=N_);

    return j==i-1 ? N_*(N_+1)/2 + i : i*N_+j-i*(1+i)/2;
  }

private:
  std::vector<PF_TYPE> V_;
  int N_;
};

template < class PF_TYPE >
class Nupack<PF_TYPE>::DPTable4
{
public:
  DPTable4() : V_(), N_(0)
  {
  }

  void resize(int n)
  {
    N_=n;
    V_.resize(N_*(N_-1)*(N_-2)*(N_-3)/2/3/4);
  }

  void fill(const PF_TYPE& v)
  {
    std::fill(V_.begin(), V_.end(), v);
  }

  PF_TYPE& operator()(int i, int d, int e, int j)
  {
    return V_[index(i, d, e, j)];
  }

  const PF_TYPE& operator()(int i, int d, int e, int j) const
  {
    return V_[index(i, d, e, j)];
  }

private:
  int index(int h, int r, int m, int s) const
  {
    int h2 = h*h;
    int h3 = h2*h;
    int h4 = h3*h;
    int m2 = m*m;
    int n2= n*n;
    int n3 = n2*n;
    int r2 = r*r;
    int r3 = r2*r;

    assert(h<=r);
    assert(r<=m);
    assert(m<=s);
    assert(s<=N_);

    return (h==r && m==s) ? V_.size()-1 :
      (-24 - 50*h - 35*h2 - 10*h3 - h4 - 36*m -12*m2 +
       12*n + 70*h*n + 30*h2*n + 4*h3*n + 24*m*n - 12*n2 -30*h*n2 -
       6*h2*n2 + 4*h*n3 + 44*r - 48*n*r + 12*n2*r + 
       24*r2 - 12*n*r2 +  4*r3 + 24*s)/24 ;
  }

private:
  std::vector<PF_TYPE> V_;
  int N_;
};

template < class PF_TYPE >
class Nupack<PF_TYPE>::DPTableX
{
public:
  DPTableX() : V_(), N_(0), D_(0)
  {
  }

  void resize(int d, int n)
  {
    assert(d>=16);
    N_=n;
    D_=d;
    int max_sz=0;
    for (int i=d; i<d+3; ++i)
      max_sz = std::max(max_sz, (N_-i)*(i-5)*(i-1)*(i-2)/2);
    V_.resize(max_sz);
  }

  void fill(const PF_TYPE& v)
  {
    std::fill(V_.begin(), V_.end(), v);
  }

  PF_TYPE& operator()(int i, int d, int e, int s)
  {
    return V_[index(i, d, e, s)];
  }

  const PF_TYPE& operator()(int i, int d, int e, int s) const
  {
    return V_[index(i, d, e, s)];
  }

  void swap(DPTableX& x)
  {
    std::swap(V_, x.V_);
    std::swap(N_, x.N_);
    std::swap(D_, x.D_);
  }

private:
  int index(int i, int h1, int m1, int s) const
  {
    int d1d2 = (D_-1)*(D_-2);
    int d5 = D_-5;
    int h1_i_1 = h1-i-1;
    assert(i+D_<N);
    assert(d-6<=s);
    assert(i<h1);
    return i*d5*d1d2/2 + s*d1d2/2 + 
      h1_i_1*(d-1) - h1_i_1*(h1-i)/2 + m1 - h1 - 1;
  }

private:
  std::vector<PF_TYPE> V_;
  int N_;
  int D_;
};


template < class PF_TYPE >
Nupack<PF_TYPE>::pf_type
Nupack<PF_TYPE>::
calculate_partition_function()
{
  Q.resize(N);
  Qb.resize(N);
  Qm.resize(N);
  Qp.resize(N);
  Qz.resize(N);
  Qg.resize(N);
  Qgl.resize(N);
  Qgr.resize(N);
  Qgls.resize(N);
  Qglr.resize(N);
  DPTableX Qx, Qx1, Qx2;

  for (int l=1; l<=N; ++l)
  {
    std::swap(Qx, Qx1);
    std::swap(Qx1, Qx2);
    Qx2.resize(l+1, N);
    std::fill(Qx2.begin(), Qx2.end(), 0.0);

    for (int i=0; i+l<=N; ++i)
    {
      int j=i+l-1;

      // Qb recursion
      if (allow_paired(i,j))
      {
        Qb(i,j) = exp(-score_hairpin(i,j)/RT);
        for (int d=i+1; d<=j-5; ++d) // all possible rightmost pairs d-e
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e))
            {
              Qb(i,j) += exp(-score_interior(i,d,e,j)/RT) * Qb(d,e);
              Qb(i,j) += Qm(i+1,d-1) * Qb(d,e) * exp(-(alpha1 + 2*alpha2 + alpha3*(j-e-1))/RT);
            }
          }
        }
        for (int d=i+1; d<=j-9; ++d) // all possible rightmost pseudoknots filling [d,e]
        {
          for (int e=d+8; e<=j-1; ++e)
          {
            energy_t temp = alpha1 + beta1m + 3*alpha2 + alpha3*(j-e-1);
            Qb(i,j) += exp(-(temp + alpha3*(d-i-1))/RT) * Qp(d,e);
            Qb(i,j) += Qm(i+1,d-1) * Qp(d,e) * exp(-temp/RT);
          }
        }
      }

      // Qg recursion
      for (int d=i+1; d<=j-5; ++d)
      {
        for (int e=d+4; e<=j-1; ++e)
        {
          Qg(i,d,e,j) += exp(-score_interior(i,d,e,j)/RT);
        }
      }
      fastiloops(i, j, l, Qg, Qx, Qx2);
      for (int d=i+6; d<=j-5; ++d)
      {
        for (int e=d+4; e<=j-1; ++e)
        {
          Qg(i,d,e,j) += Qm(i+1,d-1) * exp(-(alpha1 + 2*alpha2 + alpha3*(j-e-1))/RT);
        }
      }
      for (int d=i+1; d<=j-10; ++d)
      {
        for (int e=d+4; e<=j-6; ++e)
        {
          Qg(i,d,e,j) += exp(-(-alpha1 + 2*alpha2 + alpha3*(d-i-1))/RT) * Qm(e+1,j-1);
        }
      }
      for (int d=i+6; d<=j-10; ++d)
      {
        for (int e=d+4; e<=j-6; ++e)
        {
          Qg(i,d,e,j) += Qm(i+1,d-1) * exp(-(alpha1 + 2*alpha2)/RT) * Qm(e+1,j-1);
        }
      }
      for (int d=i+7; d<=j-6; ++d)
      {
        for (int e=d+4; e<=j-2; ++d)
        {
          for (int f=e+1; f<=j-1; ++f)
          {
            Qg(i,d,e,j) += Qgls(i+1,d,e,f) * exp(-(alpha1 + alpha2 + alpha3*(j-f-1))/RT);
          }
        }
      }
      for (int d=i+2; d<=j-11; ++d)
      {
        for (int e=d+4; e<=j-7; ++e)
        {
          for (int c=i+1; c<=d-1; ++c)
          {
            Qg(i,d,e,j) += exp(-(alpha1 + alpha2 + alpha3*(c-i-1))/RT) * Qgrs(c,d,e,j-1);
          }
        }
      }
      for (int d=i+7; d<=j-11; ++d)
      {
        for (int e=d+4; e<=j-7; ++e)
        {
          for (int c=i+6; c<=d-1; ++c)
          {
            Qg(i,d,e,j) += Qm(i+1,c-1) * Qgrs(c,d,e,j-1) * exp(-(alpha1 + alpha2)/RT);
          }
        }
      }

      // Qgls, Qgrs recursions
      for (int c=i+5; c<=j-6; ++c)
      {
        for (int d=c+1; d<=j-5; ++d)
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            Qgls(i,d,e,j) += exp(-alpha2/RT) * Qm(i,c-1) *Qg(c,d,e,j);
          }
        }
      }
      for (int d=i+1; d<=j-10; ++d)
      {
        for (int e=d+4; e<=j-6; ++e)
        {
          for (int f=e+1; f<=j-5; ++f)
          {
            Qgrs(i,d,e,j) += Qg(i,d,e,f) * Qm(f+1,j) * exp(-alpha2/RT);
          }
        }
      }

      // Qgr, Qgr recursions
      for (int d=i+1; d<=j-5; ++d)
      {
        for (int f=d+4; f<=j-1; ++f)
        {
          for (int e=d; e<=f-3; ++e)
          {
            Qgl(i,e,f,j) += Qg(i,d,f,j) * Qz(d+1,e) * exp(-beta2/RT);
          }
        }
      }
      for (int d=i+1; d<=j-4; ++d)
      {
        for (int e=d+3; e<=j-1; ++e)
        {
          for (int f=e; f<=j-1; ++f)
          {
            Qgr(i,d,e,j) += Qgl(i,d,f,j) * Qz(e,f-1);
          }
        }
      }

      // Qp recursion
      for (int d=i+2; d<=j-4; ++d)
      {
        for (int e=std::max(d+2,i+5); e<=j-3; ++e)
        {
          for (int f=e+1; f<=j-2; ++f)
          {
            Qp(i,j) += Qgl(i,d-1,e,f) * Qgr(d,e-1,f+1,j);
          }
        }
      }

      // Q, Qm, Qz recurtions
      Q(i,j) = 1;               // empty recursion
      Qz(i,j) = exp(-(beta3 * (j-i+1))/RT);
      for (int d=i; d<=j-4; ++d)
      {
        for (int e=d+4; e<=j; ++e)
        {
          Q(i,j) += Q(i,d-1) * Qb(d,e);
          Qm(i,j) += exp(-(alpha2 + alpha3*(d-i) + alpha3*(j-e))/RT) * Qb(d,e);
          Qm(i,j) += Qm(i,d-1) * Qb(d,e) * exp(-(alpha2 + alpha3*(j-e))/RT);
          Qz(i,j) += Qz(i,d-1) * Qb(d,e) * exp(-(beta2 + beta3*(j-e))/RT);
        }
      }
      for (int d=i; d<=j-8; ++d)
      {
        for (int e=d+8; e<=j; ++e)
        {
          Q(i,j) += Q(i,d-1) * Qp(d,e) * exp(-beta1/RT);
          Qm(i,j) += exp(-(beta1m + 2*alpha2 + alpha3*(j-e))/RT) * Qp(d,e);
          Qm(i,j) += Qm(i,d-1) * Qp(d,e) * exp(-(beta1m + 2*alpha2 + alpha3*(j-e))/RT);
          Qz(i,j) += Qz(i,d-1) * Qp(d,e) * exp(-(beta1m + 2*beta2 + beta3*(j-e))/RT);
        }
      }
    }
  }
  return Q(1,N);
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
fastiloops(int i, int j, int l, DPTableX& Qx, DPTableX& Qx2)
{
  if (l>=17)   // smallest subsequence not added to Qg as special case
  {
    for (int d=i+6; d<=j-10; ++d)
    {
      for (int e=d+4; e<=j-6; ++e)
      {
        int l1=4;               // explicitly add in terms for l1=4, l2>=4
        int c=i+l1+1;
        for (int l2=4; l2<=j-e-2; ++l2)
        {
          int s=l1+l2;
          int f=j-l2-1;
          energy_t temp = gamma1(s)+gamma2(std::abs(l1-l2))+gamma3(f,c,f+1,c-1);
          Qx(i,d,e,s) += exp(-temp/RT) * Qg(c,d,e,f);
        }
        if (d>=i+7)
        {
          int l2=4;             // explicitly add in terms of l1>=5, l2=4
          int f=j-l2-1;
          for (int l1=5; l1<=d-i-2; ++l1)
          {
            int s=l1+l2;
            int c=i+l1+1;
            energy_t temp = gamma1(s)+gamma2(std::abs(l1-l2))+gamma3(f,c,f+1,c-1);
            Qx(i,d,e,l1) += exp(-temp/RT) * Qg(c,d,e,f);
          }
        }
      }
    }
  }

  for (int d=i+1; d<=j-5; ++d)
  {
    for (int e=d+4; e<=j-1; ++e) 
    {
      // convert Qx into interior loop energies
      if (l>=17 && allow_paired(i,j))
      {
        for (int s=8; s<=l-9; ++s)
        {
          Qg(i,d,e,j) += Qx(i,d,e,s) * exp(-gamma3(i,j,i+1,j-1)/RT);
        }
      }
      // extend loops for future use
      if (i!=1 && j!=N)
      {
        for (int s=8; s<=l-9; ++s)
        {
          Qx2(i-1,d,e,s+2) = Qx(i,d,e,s) * exp(-(gamma1(s+2)-gamma1(s))/RT);
        }
      }
      // Add small inextensible interior loops to Qg as special cases
      for (int l1=0; l1<=std::min(3,j-e-2); ++l1)
      {
        int c=i+l1+1;
        for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
        {
          int f=j-l2-1;
          Qg(i,d,e,j) += exp(-score_interior(i,c,f,j)/RT) * Qg(c,d,e,f);
        }
      }
      // Add bulge loops and large asymmetric loops as special cases
      for (int l1=0; l1<=std::min(3,d-i-2); ++l1) // cases l1=0,1,2,3, l2>=4
      {
        int c=i+l1+1;
        for (int l2=4; l2<=j-e-2; ++l2)
        {
          int f=j-l2-1;
          Qg(i,d,e,j) += exp(-score_interior(i,c,f,j)/RT) * Qg(c,d,e,f);
        }
      }
      for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
      {
        int f=j-l2-1;
        for (int l1=4; l1<=d-i-2; ++l1)
        {
          int c=i+l1+1;
          Qg(i,d,e,j) += exp(-score_interior(i,c,f,j)/RT) * Qg(c,d,e,f);
        }
      }
    }
  }
}