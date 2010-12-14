// $Id:$

#include <cassert>
#include <vector>

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
              Qb(i,j) += Qb(d,e) *
                exp( -score_interior(i,d,e,j)/RT );

              if (d>=i+6 && wc_pair(d,e) && wc_pair(i,j))
              {
                Qb(i,j) += Qm(i+1,d-1) * Qb(d,e) *
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
              Qb(i,j) += Qp(d,e) *
                exp( -( score_multiloop() +
                        score_pk_multiloop() +
                        score_multiloop_paired(3) +
                        score_multiloop_unpaired(j-e-1 + d-i-1) +
                        score_at_penalty(i,j) +
                        score_dangle(e+1,j-1) +
                        score_dangle(i+1,d-1) )/RT );

              Qb(i,j) += Qm(i+1,d-1) * Qp(d,e) *
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
              Qg(i,d,e,j) += Qm(i+1,d-1) *
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
              Qg(i,d,e,j) += Qm(e+1,j-1) *
                exp( -( score_multiloop() +
                        score_multiloop_paired(2) +
                        score_multiloop_unpaired(d-i-1) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) +
                        score_dangle(i+1, d-1) )/RT );
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
              Qg(i,d,e,j) += Qm(i+1,d-1) * Qm(e+1,j-1) *
                exp( -( score_multiloop() +
                        score_multiloop_paired(2) +
                        score_at_penalty(i,j) +
                        score_at_penalty(d,e) )/RT );
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
                Qg(i,d,e,j) += Qgls(i+1,d,e,f) *
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
                Qg(i,d,e,j) += Qgrs(c,d,e,j-1) *
                  exp( -( score_multiloop() +
                          score_multiloop_paired(1) +
                          score_multiloop_unpaired(c-i-1) +
                          score_at_penalty(i,j) +
                          score_dangle(i+1,c-1) )/RT );
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
                Qg(i,d,e,j) += Qm(i+1,c-1) * Qgrs(c,d,e,j-1) *
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
                Qgls(i,d,e,j) += Qm(i,c-1) * Qg(c,d,e,j) *
                  exp( -( score_multiloop_paired(1) +
                          score_at_penalty(c,j) )/RT );
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
                Qgrs(i,d,e,j) += Qg(i,d,e,f) * Qm(f+1,j) *
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
              Qgl(i,e,f,j) += Qg(i,d,f,j) * Qz(d+1,e) *
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
            Qgr(i,d,e,j) += Qgl(i,d,f,j) * Qz(e,f-1);
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
                Qp(i,j) += Qg(i,a,d,e) * Qg(b,c,f,j) * Qz(e+1,f-1) * Qz(c+1,d-1) * Qz(a+1,b-1) *
                  exp( -( score_pk_paired(2) +
                          score_at_penalty(a,d) +
                          score_at_penalty(c,f) +
                          score_at_penalty(i,e) +
                          score_at_penalty(b,j) )/RT );
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
                Qp(i,j) += Qg(i,i,e,f) * Qz(i+1,d-1) * Qgr(d,e-1,f+1,j) *
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
                  Qp(i,j) += Qg(d,d,j,j) * Qz(d+1,e-1) * Qz(f+1,j-1) *
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
                  Qp(i,j) += Qgl(i,d-1,e,f) * Qgr(d,e-1,f+1,j) *
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
            Q(i,j) += Q(i,d-1) * Qb(d,e) *
              exp( -( score_at_penalty(d,e) +
                      score_dangle(e+1,j) )/RT );

            if (i!=0 && j!=N-1)
            {
              Qm(i,j) += Qb(d,e) *
                exp( -( score_multiloop_paired(1) +
                        score_multiloop_unpaired(d-i + j-e) +
                        score_at_penalty(d,e) +
                        score_dangle(e+1,j) +
                        score_dangle(i,d-1) )/RT );

              if (d>=i+5)
              {
                Qm(i,j) += Qm(i,d-1) * Qb(d,e) *
                  exp( -( score_multiloop_paired(1) +
                          score_multiloop_unpaired(j-e) +
                          score_at_penalty(d,e) +
                          score_dangle(e+1,j) )/RT );
              }
                
              Qz(i,j) += Qz(i,d-1) * Qb(d,e) *
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
          Q(i,j) += Q(i,d-1) * Qp(d,e) *
            exp( -( score_pk() +
                    score_dangle(e+1,j) )/RT );

          if (i!=0 && j!=N-1)
          {
            Qm(i,j) += Qp(d,e) *
              exp( -( score_pk_multiloop() +
                      score_multiloop_paired(2) +
                      score_multiloop_unpaired(d-i + j-e) +
                      score_dangle(e+1,j) +
                      score_dangle(i,d-1) )/RT );

            if (d>=i+5)
            {
              Qm(i,j) += Qm(i,d-1) * Qp(d,e) *
                exp( -( score_pk_multiloop() +
                        score_multiloop_paired(2) +
                        score_multiloop_unpaired(j-e) +
                        score_dangle(e+1,j) )/RT );
            }

            Qz(i,j) += Qz(i,d-1) * Qp(d,e) *
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
              Qx(i,d,e,s) += Qg(c,d,e,f) *
                exp( -( score_interior_asymmetry(std::abs(l1-l2)) +
                        score_interior_mismatch(f,c,f+1,c-1) )/RT );
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
                Qx(i,d,e,s) += Qg(c,d,e,f) *
                  exp( -( score_interior_asymmetry(std::abs(l1-l2)) +
                          score_interior_mismatch(f,c,f+1,c-1) )/RT );
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
            Qg(i,d,e,j) += Qx(i,d,e,s) *
              exp( -score_interior_mismatch(i,j,i+1,j-1)/RT );
          }
        }

        // extend loops for future use
        if (i!=0 && j!=N-1)
        {
          for (int s=8; s<=l-9; ++s)
          {
            Qx2(i-1,d,e,s+2) = Qx(i,d,e,s) *
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
                Qg(i,d,e,j) += Qg(c,d,e,f) *
                  exp( -score_interior(i,c,f,j)/RT );
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
                Qg(i,d,e,j) += Qg(c,d,e,f) *
                  exp( -score_interior(i,c,f,j)/RT );
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
                Qg(i,d,e,j) += Qg(c,d,e,f) *
                  exp( -score_interior(i,c,f,j)/RT );
              }
            }
          }
        }
      }
    }
  }
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_hairpin(int i, int j) const
{
  energy_t e=0.0;
  bool polyC = true;
  for (int k=i+1; k<j; ++k)
  {
    if (seq[k]!=BASE_C)
    {
      polyC=FALSE;
      break;
    }
  }

  int size=j-i-1;
#if 0
  if (size<3) return INF;
  if (!allow_paired(i,j)) return INF;
#else
  assert(size>=3);
  assert(allow_paired(i,j));
#endif

  e += size<=30 ?
    hairpin37[size-1] : 
    hairpin37[30 - 1] + 1.75*RT*log(size/30.0);

  if (size==3)
  {
    e += score_at_penalty(i,j);
    e += triloop37[seq[i]-1][seq[i+1]-1][seq[i+2]-1][seq[j-1]-1][seq[j]-1];
    if (polyC) e += POLYC3;
  }
  else if (size==4)
  {
    e += tetraloop37[seq[i]-1][seq[i+1]-1][seq[i+2]-1][seq[j-2]-1][seq[j-1]-1][seq[j]-1];
    e += mismatch37[seq[i+1]-1][seq[j-1]-1][pair_type(i,j)];
    if (polyC) e += POLYCSLOPE*size + POLYCINT;
  }
  else /*if (size>4)*/
  {
    e += mismatch37[seq[i+1]-1][seq[j-1]-1][pair_type(i,j)];
    if (polyC) e += POLYCSLOPE*size + POLYCINT;
  }
  return e;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_loop(int l) const
{
  return l<=30 ?
    interior37[l-1] :
    interior37[30-1]+1.75*RT*log(l/30.0);
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior(int i, int h, int m, int j) const
{
  int l1 = h - i - 1;
  int l2 = j - m - 1;
  int size = l1 + l2;
  energy_t e = 0;

  // helix
  if (size==0)
  {
    e += stack37[get_type(i,j)][get_type(h,m)];
  }
  
  // bulge
  else if (l1==0 || l2==0)
  {
    e += size<=30 ?
      bulge37[size-1] :
      bulge37[30-1] + 1.75*RT*log(size/30.0);

    if (l1+l2==1)           //single bulge...treat as a stacked region
    {
      e += score_helix(i,h,m,j);
      e -= SALT_CORRECTION;
    }
    else
    {
      e += score_at_penalty(i,j);
      e += score_at_penalty(h,m);
    }
  }

  // interior loop
  else if (l1>0 && l2>0)
  {
    int asymmetry = std::abs(l1-l2);
    if (asymmetry>1 || size>4)
    {
      e += score_interior_asymmetry(asymmetry);
      if (l1>1 && l2>1)
      {
        e += score_interior_mismatch(m, h, m+1, h-1);
        e += score_interior_mismatch(i, j, i+1, j-1);
      }
      else if (l1==1 || l2==1)
      {                         // assume AA terminal mismatch?
        e += score_interior_mismatch(m, h);
        e += score_interior_mismatch(i, j);
      }
      else
      {
        assert(!"unclassified interior loop");
        exit(1);
      }
    }
    else if (l1==1 && l2==1)
      e += int11[pair_type(i,j)][pair_type(h,m)][seq[i+1]-1][seq[j-1]-1];
    else if (l1==2 && l2==2)
      e += int22[pair_type(i,j)][pair_type(h,m)][seq[i+1]-1][seq[j-1]-1][seq[i+2]-1][seq[j-2]-1];
    else if (l1==1 && l2==2)
      e += int12[get_type(i,j)][seq[j-2]-1][seq[i+1]-1][get_type(h,m)][seq[j-1]-1];
    else if (l1==2 && l2==1)
      e += int12[get_type(m,h)][seq[i+1]-1][seq[j-1]-1][get_type(m,h)][seq[i+2]-1];
    else
    {
      assert(!"error in tabulated interior loop");
      exit(1);
    }
  }
  else
  {
    assert(!"improperly classifed interior loop");
    exit(1);
  }
  return e;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior_mismatch(int i, int j, int k, int l) const
{
/*
  Interior Mismatch calculation

  This calculates the Mismatch interaction energies between positions
  1 -> 5' i k 3'
  2 -> 3' j l 5'
  Interactions energies taken from file tstacki2.dgd.
*/
  return mismatch37[seq[k]-1][seq[l]-1][pair_type(i,j)];
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_interior_asymmetry(int l1, int l2) const
{
  energy_t e=0.0;
  int size = l1+l2;
  int asymmetry = std::abs(l1-l2);
  e += size<=30 ?
    interior37[size-1] :
    interior37[30-1] + 1.75*RT*log(size/30.0);

  //Asymmetry rountine copied from efn.f in Zuker's mfold package.
  int asymmetry_index = 4;
  if( l1 < asymmetry_index) asymmetry_index = l1;
  if( l2 < asymmetry_index) asymmetry_index = l2;
  if( asymmetry*asymmetry_penalty[ asymmetry_index - 1] < max_asymmetry )
    e += asymmetry * asymmetry_penalty[ asymmetry_index - 1];
  else
    e += max_asymmetry; // MAX asymmetry penalty

  return e;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop() const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop_paired(int n) const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop_unpaired(int n) const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop_paired(int n) const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_multiloop_unpaired(int n) const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_at_penalty(int i, int j) const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_dangle(int i, int j) const
{
  energy_t d5=0.0, d3=0.0;
  if (j!=N-1)
    d3 = dangle3_27[pair_type(i-1,j+1)][seq[j]-1];
  if (i!=0)
    d5 = dangle5_27[pair_type(i-1,j+1)][seq[i]-1];
  return d3+d5;
}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk() const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_multiloop() const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_pk() const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_paired(int n) const
{

}

template <class PF_TYPE>
energy_t
Nupack<PF_TYPE>::
score_pk_unpaired(int n) const
{

}

