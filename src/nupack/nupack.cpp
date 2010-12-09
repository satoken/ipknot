// $Id:$

template < class PF_TYPE >
Nupack<PF_TYPE>::pf_type
Nupack<PF_TYPE>::
calculate_partition_function()
{
  DPTablex Qx;
  DPTablex Qx1;
  DPTablex Qx2;

  for (int l=1; l<=N; ++l)
  {
    std::swap(Qx, Qx1);
    std::swap(Qx1, Qx2);
    std::fill(Qx2.begin(), Qx2.end(), 0.0);

    for (int i=1; i<=N-l+1; ++i)
    {
      int j=i+l-1;

      // Qb recursion
      Qb(i,j) = exp(-score_hairpin(i,j)/RT);
      for (int d=i+1; d<=j-5; ++d) // all possible rightmost pairs d-e
      {
        for (int e=d+4; e<=j-1; ++e)
        {
          Qb(i,j) += exp(-score_interior(i,d,e,j)/RT) * Qb(d,e);
          Qb(i,j) += Qm(i+1,d-1) * Qb(d,e) * exp(-(alpha1 + 2*alpha2 + alpha3*(j-e-1))/RT);
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

      // Qg recursion
      for (int d=i+1; d<=j-5; ++d)
      {
        for (int e=d+4; e<=j-1; ++e)
        {
          Qg(i,d,e,j) += exp(-score_interior(i,d,e,j)/RT);
        }
      }
      fastiloop(i, j, l, Qg, Qx, Qx2);
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
fastiloops()
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
      if (l>=17 && allow_bp(i,j))
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
