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
    int n2 = n*n;
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
        Qb(i,j) = exp( -score_hairpin(i,j)/RT );

        for (int d=i+1; d<=j-5; ++d) // all possible rightmost pairs d-e
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e))
            {
              Qb(i,j) +=
                exp( -score_interior(i,d,e,j)/RT ) *
                Qb(d,e);

              if (d>=i+6 && wc_pair(d,e) && wc_pair(i,j))
              {
                Qb(i,j) +=
                  Qm(i+1,d-1) *
                  Qb(d,e) *
                  exp( -( score_multiloop() +
                          score_multiloop_paired(2) +
                          score_multiloop_unpaired(j-e-1) +
                          score_at_penalty(i,j) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j-1) )/RT ) ;
              }
            }
          }
        }

        if (wc_pair(i,j))
        {
          for (int d=i+1; d<=j-6; ++d) // all possible rightmost pseudoknots filling [d,e]
          {
            for (int e=d+5; e<=j-1; ++e)
            {
              Qb(i,j) +=
                exp( -( score_multiloop() +
                        score_pk_multiloop() +
                        score_multiloop_paired(3) +
                        score_multiloop_unpaired(j-e-1 + d-i-1) +
                        score_at_penalty(i,j) +
                        score_dangle(e+1,j-1) +
                        score_dangle(i+1,d-1) )/RT ) *
                Qp(d,e);

              Qb(i,j) +=
                Qm(i+1,d-1) *
                Qp(d,e) *
                exp( -( score_multiloop() +
                        score_pk_multiloop() +
                        score_multiloop_paired(3) +
                        score_multiloop_unpaired(j-e-1) +
                        score_at_penalty(i,j) +
                        score_dangle(e+1,j-1) )/RT );
            }
          }
        }
      }

      // Qg recursion
      if (allow_paired(i,j))
      {
        // case 0: only 1 pair
        Qg(i,i,j,j) = 1;

        // case 1: terminal inner pair
        for (int d=i+1; d<=j-5; ++d)
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e))
              Qg(i,d,e,j) += exp( -score_interior(i,d,e,j)/RT );
          }
        }
      }

      fastiloops(i, j, Qg, Qx, Qx2);

      if (allow_paired(i,j) && wc_pair(i,j))
      {

        // case 2: multiloop left
        for (int d=i+6; d<=j-5; ++d)
        {
          for (int e=d+4; e<=j-1; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e))
            {
              Qg(i,d,e,j) +=
                Qm(i+1,d-1) *
                exp( -( score_multiloop() +
                        score_multiloop_paired(2) +
                        score_multiloop_unpaired(j-e-1) + 
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1, j-1) )/RT );
            }
          }
        }

        // case 3: multiloop right
        for (int d=i+1; d<=j-10; ++d)
        {
          for (int e=d+4; e<=j-6; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e))
            {
              Qg(i,d,e,j) +=
                exp( -( score_multiloop() +
                        score_multiloop_paired(2) +
                        score_multiloop_unpaired(d-i-1) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(i+1, d-1) )/RT ) *
                Qm(e+1,j-1);
            }
          }
        }

        // case 4: multiloop both sides
        for (int d=i+6; d<=j-10; ++d)
        {
          for (int e=d+4; e<=j-6; ++e)
          {
            if (allow_paired(d,e) && wc_pair(d,e))
            {
              Qg(i,d,e,j) +=
                Qm(i+1,d-1) *
                exp( -( score_multiloop() +
                        score_multiloop_paired(2) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) )/RT ) *
                Qm(e+1,j-1);
            }
          }
        }

        // case 5: interior loop + multi left
        for (int d=i+7; d<=j-6; ++d)
        {
          for (int e=d+4; e<=j-2; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int f=e+1; f<=j-1; ++f)
              {
                Qg(i,d,e,j) +=
                  Qgls(i+1,d,e,f) *
                  exp( -( score_multiloop() +
                          score_multiloop_paired(1) +
                          score_multiloop_unpaired(j-f-1) +
                          score_at_penalty(i,j) +
                          score_dangle(f+1,j-1) )/RT );
              }
            }
          }
        }

        // case 6: interior loop + multi right
        for (int d=i+2; d<=j-11; ++d)
        {
          for (int e=d+4; e<=j-7; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int c=i+1; c<=d-1; ++c)
              {
                Qg(i,d,e,j) +=
                  exp( -( score_multiloop() +
                          score_multiloop_paired(1) +
                          score_multiloop_unpaired(c-i-1) +
                          score_at_penalty(i,j) +
                          score_dangle(i+1,c-1) )/RT ) *
                  Qgrs(c,d,e,j-1);
              }
            }
          }
        }

        // case 7: interior loop + multi both sides
        for (int d=i+7; d<=j-11; ++d)
        {
          for (int e=d+4; e<=j-7; ++e)
          {
            if (allow_paired(d,e))
            {
              for (int c=i+6; c<=d-1; ++c)
              {
                Qg(i,d,e,j) +=
                  Qm(i+1,c-1) *
                  Qgrs(c,d,e,j-1) *
                  exp( -( score_multiloop() +
                          score_multiloop_paired(1) +
                          score_at_penalty(i,j) )/RT );
              }
            }
          }
        }
      }

      // Qgls recursion
      for (int c=i+5; c<=j-6; ++c)
      {
        if (allow_paired(c,j) && wc_pair(c,j))
        {
          for (int d=c+1; d<=j-5; ++d)
          {
            for (int e=d+4; e<=j-1; ++e)
            {
              if (allow_paired(d,e))
              {
                Qgls(i,d,e,j) +=
                  exp( -( score_multiloop_paired(1) +
                          score_at_penalty(c,j) )/RT ) *
                  Qm(i,c-1) *
                  Qg(c,d,e,j);
              }
            }
          }
        }
      }

      // Qgrs recursion
      for (int d=i+1; d<=j-10; ++d)
      {
        for (int e=d+4; e<=j-6; ++e)
        {
          if (allow_paired(d,e))
          {
            for (int f=e+1; f<=j-5; ++f)
            {
              if (allow_paired(i,f) && wc_pair(i,f))
              {
                Qgrs(i,d,e,j) +=
                  Qg(i,d,e,f) *
                  Qm(f+1,j) *
                  exp( -( score_multiloop_paired(1) +
                          score_at_penalty(i,f) )/RT );
              }
            }
          }
        }
      }

      // Qgl recursions
      for (int d=i+1; d<=j-5; ++d)
      {
        for (int f=d+4; f<=j-1; ++f)
        {
          if (allow_paired(d,f) && wc_pair(d,f))
          {
            for (int e=d; e<=f-2; ++e) // f-3???
            {
              Qgl(i,e,f,j) +=
                Qg(i,d,f,j) *
                Qz(d+1,e) *
                exp( -( score_pk_paired(1) +
                        score_at_penalty(d,f) )/RT );
            }
          }
        }
      }

      // Qgr recursion
      for (int d=i+1; d<=j-3; ++d)
      {
        for (int e=d+2; e<=j-1; ++e)
        {
          for (int f=e; f<=j-1; ++f)
          {
            Qgr(i,d,e,j) +=
              Qgl(i,d,f,j) *
              Qz(e,f-1);
          }
        }
      }

      // Qp recursion
      // case 1: both Qg are exactly 1 pair
      // first case is exactly 1 pair per Og
      if (j-i>4)
      {
        int a=i;
        int f=j;
        for (int b=a+1; b<=j-4; ++b)
        {
          if (allow_paired(b,j) && wc_pair(b,j))
          {
            int c=b;
            for (int d=std::max(c+1,a+4); d<=j-1; ++d)
            {
              if (allow_paired(a,d) && wc_pair(a,d))
              {
                int e=d;
                Qp(i,j) +=
                  Qg(i,a,d,e) *
                  Qg(b,c,f,j) *
                  exp( -( score_pk_paired(2) +
                          score_at_penalty(a,d) +
                          score_at_penalty(c,f) +
                          score_at_penalty(i,e) +
                          score_at_penalty(b,j) )/RT ) *
                  Qz(e+1,f-1) *
                  Qz(c+1,d-1) *
                  Qz(a+1,b-1);
              }
            }
          }
        }
      }

      if (j-i>6)
      {
        // case 2 left Og is exactly 1 pair, right is 2+
        for (int d=i+1; d<=j-6; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+2,i+4); e<=j-2; ++e)
            {
              int f=e;
              if (allow_paired(i,f) && wc_pair(i,f))
              {
                Qp(i,j) +=
                  Qg(i,i,e,f) *
                  Qz(i+1,d-1) *
                  Qgr(d,e-1,f+1,j) *
                  exp( -( score_pk_paired(1) +
                          score_at_penalty(d,j) +
                          score_at_penalty(i,f)*2 )/RT );
              }
            }
          }
        }

        // case 2 left Qg is 2+ pairs, right is 1
        for (int d=i+2; d<=j-4; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+1,i+4); e<=j-2; ++e)
            {
              for (int f=e+1; f<=j-1; ++f)
              {
                if (allow_paired(i,f) && wc_pair(i,f))
                {
                  Qp(i,j) +=
                    Qg(d,d,j,j) *
                    Qz(d+1,e-1) *
                    Qz(f+1,j-1)
                    exp( -( score_pk_paired(1) +
                            score_at_penalty(d,j)*2 +
                            score_at_penalty(i,f) )/RT );
                }
              }
            }
          }
        }
      }

      // otherwise
      if (j-i>7)
      {
        for (int d=i+2; d<=j-4; ++d)
        {
          if (allow_paired(d,j) && wc_pair(d,j))
          {
            for (int e=std::max(d+2,i+5); e<=j-3; ++e)
            {
              for (int f=e+1; f<=j-2; ++f)
              {
                if (allow_paired(i,f) && wc_pair(i,f))
                {
                  Qp(i,j) +=
                    Qgl(i,d-1,e,f) *
                    Qgr(d,e-1,f+1,j) *
                    exp( -( score_at_penalty(d,j) +
                            score_at_penalty(i,j) )/RT );
                }
              }
            }
          }
        }
      }

      // Q, Qm, Qz recurtions
      Q(i,j) = exp( -score_dangle(i,j)/RT ); // empty recursion

      if (i!=0 && j!=N-1)
      {
        Qz(i,j) =
          exp( -( score_dangle(i,j) +
                  score_pk_unpaired(j-i+1) )/RT );
      }

      for (int d=i; d<=j-4; ++d)
      {
        for (int e=d+4; e<=j; ++e)
        {
          if (allow_paired(d,e) && wc_pair(d,e))
          {
            Q(i,j) +=
              Q(i,d-1) *
              Qb(d,e) *
              exp( -( score_at_penalty(d,e) +
                      score_dangle(e+1,j) )/RT );

            if (i!=0 && j!=N-1)
            {
              Qm(i,j) +=
                exp( -( score_multiloop_paired(1) +
                        score_multiloop_unpaired(d-i + j-e) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1,j) +
                        score_dangle(i,d-1) )/RT ) *
                Qb(d,e);

              if (d>=i+5)
              {
                Qm(i,j) +=
                  Qm(i,d-1) *
                  Qb(d,e) *
                  exp( -( score_multiloop_paired(1) +
                          score_multiloop_unpaired(j-e) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j) )/RT );
              }
                
              Qz(i,j) +=
                Qz(i,d-1) *
                Qb(d,e) *
                exp( -( score_pk_paired(1) +
                        score_pk_unpaired(j-e) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1,j) )/RT );
            }
          }
        }
      }

      for (int d=i; d<=j-5; ++d)
      {
        for (int e=d+5; e<=j; ++e)
        {
          Q(i,j) +=
            Q(i,d-1) *
            Qp(d,e) *
            exp( -( score_pk() +
                    score_dangle(e+1,j) )/RT );

          if (i!=0 && j!=N-1)
          {
            Qm(i,j) +=
              exp( -( score_pk_multiloop() +
                      score_multiloop_paired(2) +
                      score_multiloop_unpaired(d-i + j-e) +
                      score_dangle(e+1,j) +
                      score_dangle(i,d-1) )/RT ) *
              Qp(d,e);

            if (d>=i+5)
            {
              Qm(i,j) +=
                Qm(i,d-1) *
                Qp(d,e) *
                exp( -( score_pk_multiloop() +
                        score_multiloop_paired(2) +
                        score_multiloop_unpaired(j-e) +
                        score_dangle(e+1,j) )/RT );
            }

            Qz(i,j) +=
              Qz(i,d-1) *
              Qp(d,e) *
              exp( -( score_pk_pk() +
                      score_pk_paired(2) +
                      score_pk_unpaired(j-e) +
                      score_dangle(e+1,j) )/RT );
          }
        }
      }
    }
  }

  return Q(0,N-1);
}

template < class PF_TYPE >
void
Nupack<PF_TYPE>::
fastiloops(int i, int j, DPTableX& Qx, DPTableX& Qx2)
{
  int l=j-i+1;
  if (l>=17)   // smallest subsequence not added to Qg as special case
  {
    for (int d=i+6; d<=j-10; ++d)
    {
      for (int e=d+4; e<=j-6; ++e)
      {
        if (allow_paired(d,e))
        {
          int l1=4;               // explicitly add in terms for l1=4, l2>=4
          int c=i+l1+1;
          for (int l2=4; l2<=j-e-2; ++l2)
          {
            int s=l1+l2;
            int f=j-l2-1;
            if (allow_paired(c,f))
            {
              Qx(i,d,e,s) +=
                exp( -( score_interior_asymmetry(std::abs(l1-l2)) +
                        score_interior_mismatch(f,c,f+1,c-1) )/RT ) *
                Qg(c,d,e,f);
            }
          }

          if (d>=i+7)
          {
            int l2=4;             // explicitly add in terms of l1>=5, l2=4
            int f=j-l2-1;
            for (int l1=5; l1<=d-i-2; ++l1)
            {
              int s=l1+l2;
              int c=i+l1+1;
              if (allow_paired(c,f))
              {
                Qx(i,d,e,s) +=
                  exp( -( score_interior_asymmetry(std::abs(l1-l2)) +
                          score_interior_mismatch(f,c,f+1,c-1) )/RT ) *
                  Qg(c,d,e,f);
              }
            }
          }
        }
      }
    }
  }

  for (int d=i+1; d<=j-5; ++d)
  {
    for (int e=d+4; e<=j-1; ++e) 
    {
      if (allow_paired(d,e))
      {
        // convert Qx into interior loop energies
        if (l>=17 && allow_paired(i,j))
        {
          for (int s=8; s<=l-9; ++s)
          {
            Qg(i,d,e,j) +=
              Qx(i,d,e,s) *
              exp( -score_interior_mismatch(i,j,i+1,j-1)/RT );
          }
        }

        // extend loops for future use
        if (i!=0 && j!=N-1)
        {
          for (int s=8; s<=l-9; ++s)
          {
            Qx2(i-1,d,e,s+2) =
              Qx(i,d,e,s) *
              exp( -(score_loop(s+2)-score_loop(s))/RT );
          }
        }
      
        if (allow_paired(i,j))
        {
          // Add small inextensible interior loops to Qg as special cases
          for (int l1=0; l1<=std::min(3,d-i-2); ++l1)
          {
            int c=i+l1+1;
            for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
            {
              int f=j-l2-1;
              if (allow_paired(c,f))
              {
                Qg(i,d,e,j) +=
                  exp( -score_interior(i,c,f,j)/RT ) *
                  Qg(c,d,e,f);
              }
            }
          }
          // Add bulge loops and large asymmetric loops as special cases
          for (int l1=0; l1<=std::min(3,d-i-2); ++l1) // cases l1=0,1,2,3, l2>=4
          {
            int c=i+l1+1;
            for (int l2=4; l2<=j-e-2; ++l2)
            {
              int f=j-l2-1;
              if (allow_paired(c,f))
              {
                Qg(i,d,e,j) +=
                  exp( -score_interior(i,c,f,j)/RT ) *
                  Qg(c,d,e,f);
              }
            }
          }
          for (int l2=0; l2<=std::min(3,j-e-2); ++l2)
          {
            int f=j-l2-1;
            for (int l1=4; l1<=d-i-2; ++l1)
            {
              int c=i+l1+1;
              if (allow_paired(c,f))
              {
                Qg(i,d,e,j) +=
                  exp( -score_interior(i,c,f,j)/RT ) *
                  Qg(c,d,e,f);
              }
            }
          }
        }
      }
    }
  }
}
