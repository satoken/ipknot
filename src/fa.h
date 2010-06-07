/*
 * $Id: fa.h 28 2010-04-20 02:25:33Z satoken $
 * 
 * Copyright (C) 2010 Kengo Sato
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

#ifndef __INC_FA_H__
#define __INC_FA_H__

#include <boost/spirit/iterator/file_iterator.hpp>
#include <string>
#include <list>
//#include "rna.h"

class Fasta
{
public:
  Fasta() : name_(), seq_(), str_() { }

  Fasta(const std::string& name,
	const std::string& seq,
	const std::string& str="")
    : name_(name), seq_(seq), str_(str)
  { }

  Fasta(const Fasta& fa)
    : name_(fa.name_), seq_(fa.seq_), str_(fa.str_)
  { }
    
  Fasta&
  operator=(const Fasta& fa)
  {
    if (this != &fa) {
      name_ = fa.name_;
      seq_ = fa.seq_;
      str_ = fa.str_;
    }
    return *this;
  }

  const std::string& name() const { return name_; }
  const std::string& seq() const { return seq_; }
  std::string& seq() { return seq_; }
  const std::string& str() const { return str_; }
  unsigned int size() const { return seq_.size(); }

  static
  unsigned int
  load(std::list<Fasta>& data, const char* file);
  
  bool
  load(boost::spirit::file_iterator<>& fi);

private:
  std::string name_;
  std::string seq_;
  std::string str_;
};

#endif //  __INC_FA_H__

// Local Variables:
// mode: C++
// End:
