#include "bpseq.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>

BPSEQ::
BPSEQ()
  : seq_(), bp_()
{
}

BPSEQ::
BPSEQ(const std::string& seq, const std::string& paren)
  : seq_(seq), bp_(seq.size()+1,0)
{
  std::stack<unsigned int> st;
  for (unsigned int i=0; i!=paren.size(); ++i)
  {
    switch (paren[i])
    {
      case '(':
        st.push(i); break;
      case ')':
      {
        int j=st.top();
        st.pop();
        bp_[i+1]=j+1;
        bp_[j+1]=i+1;
      }
      break;
      default: break;
    }
  }
}

bool
BPSEQ::
load(const char* fname)
{
  int i, j, l=0;
  char c;
  std::string s;
  std::ifstream is(fname);
  if (!is.is_open()) return false;
  while (std::getline(is, s))
  {
    if (s[0]>='0' || s[0]<=9)
    {
      std::istringstream ss(s);
      ss >> i >> c >> s;
      l = std::max(l, i);
    }
  }
  
  seq_.clear(); seq_.resize(l);
  bp_.clear(); bp_.resize(l+1,0);
  is.close();
  
  is.open(fname);
  if (!is.is_open()) return false;
  while (std::getline(is, s))
  {
    std::istringstream ss(s);
    ss >> i >> c >> s;
    seq_[i-1] = c;
    if (s[0]>='0' && s[0]<='9')
      bp_[i] = std::atoi(s.c_str());
    else if (s[0]=='x')
      bp_[i] = 0;
    else if (s[0]=='.')
      bp_[i] = DOT;
    else if (s[0]=='<')
      bp_[i] = L;
    else if (s[0]=='>')
      bp_[i] = R;
    else if (s[0]=='|')
      bp_[i] = LR;
  }

  return true;
}

bool
BPSEQ::
save(const char* fname) const
{
  std::ofstream os(fname);
  if (!os.is_open()) return false;
  for (unsigned int i=1; i!=bp_.size(); ++i)
    os << i << " " << seq_[i-1] << " " << bp_[i] << std::endl;
  return true;
}


#ifdef TEST
int
main(int argc, char* argv[])
{
  BPSEQ bpseq("GGCAUUCC",
              "((....))");
  bpseq.save(argv[1]);
  return 0;
}
#endif
