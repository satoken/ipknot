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

#ifndef __INC_ALN_H__
#define __INC_ALN_H__

#include <string>
#include <list>
#include <map>

class Aln
{
public:
  Aln() : name_(), seq_()
  {}

  Aln(const std::list<std::string>& name,
      const std::list<std::string>& seq)
    : name_(name), seq_(seq)
  {}

  Aln(const Aln& aln)
    : name_(aln.name_), seq_(aln.seq_)
  {}

  Aln&
  operator=(const Aln& aln)
  {
    if (this != &aln) {
      name_ = aln.name_;
      seq_  = aln.seq_;
    }
    return *this;
  }
    
  const std::list<std::string>& name() const { return name_; }
  const std::list<std::string>& seq() const { return seq_; }
  std::list<std::string>& seq() { return seq_; }
  unsigned int size() const { return seq_.begin()->size(); }
  unsigned int num_aln() const { return name_.size(); }
  std::string consensus() const;
#ifdef HAVE_LIBRNA
  float energy_of_struct(const std::string& paren) const;
#endif

  static
  unsigned int
  load(std::list<Aln>& data, const char* file );

private:
  std::list<std::string> name_;
  std::list<std::string> seq_;
};

#endif	// __INC_ALN_H__

// Local Variables:
// mode: C++
// End:
