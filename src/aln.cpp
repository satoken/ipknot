/*
 * $Id$
 *
 * Copyright (C) 2008-2010 Kengo Sato
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include "aln.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <stack>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cerrno>

#ifdef HAVE_LIBRNA
namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/PS_dot.h>
#include <ViennaRNA/aln_util.h>
  extern int eos_debug;
};
};
#endif

static
void
append(std::list<Aln>& data, const std::list<std::string>& name, 
       const std::map<std::string,std::string>& seq_map)
{
  std::list<std::string> seq;
  for (std::list<std::string>::const_iterator it=name.begin(); it!=name.end(); ++it)
  {
    seq.push_back(seq_map.find(*it)->second);
  }
  Aln aln(name, seq);
  data.push_back(aln);
}

// static
unsigned int
Aln::
load(std::list<Aln>& data, const char* filename)
{
  std::ifstream ifs(filename);
  if (!ifs) 
    throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);    

  std::map<std::string,std::string> seq_map;
  std::list<std::string> name;
  std::string line;
  while (std::getline(ifs, line))
  {
    if (line.find("CLUSTAL")==0 || line.find("PROBCONS")==0)
    {
      if (!name.empty())
      {
        append(data, name, seq_map);
        name.clear();
        seq_map.clear();
      }
    }
    else if (!line.empty() && line[0]!=' ')
    {
      std::istringstream ss(line);
      std::string n, s;
      if (ss >> n >> s)
      {
        std::map<std::string,std::string>::iterator it=seq_map.find(n);
        if (it==seq_map.end())
        {
          name.push_back(n);
          seq_map.insert(std::make_pair(n,s));
        }
        else
        {
          it->second += s;
        }
      }
    }
  }
  if (!name.empty())
  {
    append(data, name, seq_map);
    name.clear();
    seq_map.clear();
  }

  return data.size();
}

std::string
Aln::
consensus() const
{
#ifdef HAVE_LIBRNA
  // prepare an alignment
  unsigned int length = seq_.front().size();
  char **seqs = new char*[seq_.size()+1];
  seqs[seq_.size()] = NULL;
  std::list<std::string>::const_iterator x;
  unsigned int i=0;
  for (x=seq_.begin(); x!=seq_.end(); ++x) {
    assert(x->size()==length);
    seqs[i] = new char[length+1];
    strcpy(seqs[i], x->c_str());
    for (char* p=seqs[i++]; *p!=0; ++p)
      if (*p>='a' && *p<='z') 
        *p+='A'-'a';
  }

  // make a consensus string
  //char *cons = Vienna::consensus((const char**)seqs);
  char *cons = Vienna::consens_mis((const char**)seqs);
  std::string ret(cons);

  // destroy the alignment
  for (unsigned int i=0; seqs[i]!=NULL; ++i) delete[] seqs[i];
  delete[] seqs;
  free(cons);

  return ret;
#else
  return seq_.front();
#endif
}

#ifdef HAVE_LIBRNA
float
Aln::
energy_of_struct(const std::string& paren) const
{
  float e=0.0;
  std::list<std::string>::const_iterator seq;
  for (seq=seq_.begin(); seq!=seq_.end(); ++seq) {
    std::vector<unsigned int> ppos(paren.size(), -1u);
    std::stack<unsigned int> st;
    for (unsigned int i=0; i!=paren.size(); ++i)
    {
      switch (paren[i])
      {
        case '(':
          st.push(i);
          break;
        case ')':
          ppos[i] = st.top();
          ppos[st.top()] = i;
          st.pop();
          break;
        default:
          break;
      }            
    }
    std::string s(*seq);
    std::string p(paren);
    for (unsigned int i=0; i!=s.size(); ++i)
    {
      if (s[i]=='-')
      {
        p[i]='-';
        if (ppos[i]!=-1u) {
          if (ppos[i]>i || p[ppos[i]]!='-') p[ppos[i]]=' ';
        }
      }
    }
    s.erase(std::remove(s.begin(), s.end(), '-'), s.end());
    p.erase(std::remove(p.begin(), p.end(), '-'), p.end());
#ifdef HAVE_VIENNA20
    e += Vienna::energy_of_structure(s.c_str(), p.c_str(), -1);
#else
    Vienna::eos_debug = -1;
    e += Vienna::energy_of_struct(s.c_str(), p.c_str());
#endif
  }
  return e/num_aln();
}
#endif

#ifdef TEST
#include <iostream>

int main(int argc, char* argv[])
{
  std::list<Aln> data;
  Aln::load(data, argv[1]);
  for (std::list<Aln>::const_iterator a=data.begin(); a!=data.end(); ++a) 
  {
    std::list<std::string>::const_iterator n=a->name().begin();
    std::list<std::string>::const_iterator s=a->seq().begin();
    while (n!=a->name().end() && s!=a->seq().end())
    {
      std::cout << *n << " " << *s << std::endl;
      ++n; ++s;
    }
    std::cout << std::endl;
  }

  return 0;
}
#endif
