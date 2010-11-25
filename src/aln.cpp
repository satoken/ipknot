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
#include <iostream>
#include <sstream>
#include <string>
#include <deque>
#include <map>
#include <stack>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <cstring>

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

using namespace BOOST_SPIRIT_CLASSIC_NS;

struct aln_parser : public grammar< aln_parser >
{
  class format_error : public std::logic_error
  {
  public:
    format_error(const std::string& msg) : std::logic_error(msg) {}
  };
  
  struct WA
  {
    WA() : cur_index(0), cur_name(), cur_seq(), names(), seqs() { }

    unsigned int cur_index;
    std::string cur_name;
    std::string cur_seq;
    std::deque<std::string> names;
    std::deque<std::string> seqs;
  };

  WA& wa;
  aln_parser(WA& x) : wa(x) {}

  struct push_seq
  {
    WA& wa;
    push_seq(WA& x) : wa(x) { }

    template < class Ite >
    void operator()(Ite i1, Ite i2) const
    {
      assert(wa.names.size()==wa.seqs.size());
      if (wa.cur_index >= wa.names.size()) {
	wa.names.push_back(wa.cur_name);
	wa.seqs.push_back(wa.cur_seq);
      } else if (wa.names[wa.cur_index] == wa.cur_name) {
	wa.seqs[wa.cur_index] += wa.cur_seq;
      } else {
	throw format_error("format error: broken sequence name consistency");
      }
      wa.cur_index++;
    }
  };

  struct reset_index
  {
    WA& wa;
    reset_index(WA& x) : wa(x) { }

    template < class Ite >
    void operator()(Ite i1, Ite i2) const
    {
      unsigned int l=wa.seqs[0].size();
      for (unsigned int i=1; i!=wa.seqs.size(); ++i) {
	if (l!=wa.seqs[i].size())
	  throw format_error("format error: broken sequence length consistency");
      }
      wa.cur_index = 0;
    }
  };

  template <class ScannerT>
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t aln;
    rule_t header;
    rule_t head_word;
    rule_t body;
    rule_t body_part;
    rule_t empty;
    rule_t seq;
    rule_t status;
    rule_t status_or_empty;
    
    definition(const aln_parser& self)
    {
      aln = header >> +empty >> body >> *empty;
      head_word = str_p("CLUSTAL") | str_p("PROBCONS");
      header =  head_word >> +print_p >> eol_p;
      empty = *blank_p >> eol_p;
      body_part = +seq[push_seq(self.wa)];
      body = +(body_part[reset_index(self.wa)] >> !status >> *empty);
      seq
	= (+graph_p - head_word)[assign_a(self.wa.cur_name)]
	>> +blank_p >> (+graph_p)[assign_a(self.wa.cur_seq)]
	>> *blank_p >> eol_p;
      status = blank_p >> *(chset<>("*:.") | blank_p) >> eol_p;
    }

    const rule_t& start() const { return aln; }
  };
};

unsigned int
Aln::
load(file_iterator<>& fi)
{
  file_iterator<> s = fi;
  aln_parser::WA wa;
  aln_parser parser(wa);
  parse_info<file_iterator<> > info =  parse(fi, fi.make_end(), parser);
  if (!info.hit) {
    fi = s;
    return 0;
  } else {
    fi = info.stop;
    std::copy(wa.names.begin(), wa.names.end(),
	      std::back_insert_iterator<std::list<std::string> >(name_));
    std::copy(wa.seqs.begin(), wa.seqs.end(),
	      std::back_insert_iterator<std::list<std::string> >(seq_));
    return info.length;
  }
}

// static
unsigned int
Aln::
load(std::list<Aln>& data, const char* filename)
{
  unsigned int n=0;
  file_iterator<> fi(filename);
  if (!fi) {
    std::ostringstream os;
    os << filename << ": no such file";
    throw os.str().c_str();
    //return false;
  }
  while (1) {
    Aln aln;
    if (aln.load(fi)>0) {
      n++;
      data.push_back(aln);
    } else {
      break;
    }
  }
  return n;
}

#ifdef TEST
int
main(int argc, char* argv[])
{
  boost::spirit::file_iterator<> fi(argv[1]);
  if (!fi) {
    //perror(argv[1]);
    return 1;
  }
  Aln aln;
#if 0
  return aln.load(fi) ? 0 : 1;
#else
  std::cout << aln.load(fi) << std::endl;
  std::copy(aln.seq().begin(), aln.seq().end(),
	    std::ostream_iterator<std::string>(std::cout, "\n"));
  return 0;
#endif
}
#endif

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
    boost::to_upper(seqs[i++]);
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
    e += Vienna::energy_of_struct(s.c_str(), p.c_str());
  }
  return e/num_aln();
}
#endif
