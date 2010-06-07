/*
 * $Id:$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of IPknot.
 *
 * IPknot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * IPknot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with IPknot.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INC_IP_H__
#define __INC_IP_H__

class IPimpl;

class IP
{
public:
  typedef enum {MIN, MAX} DirType;
  typedef enum {FR, LO, UP, DB, FX} BoundType;

public:
  IP(DirType dir, int n_th);
  ~IP();
  int make_variable(double coef);
  int make_constraint(BoundType bnd, double l, double u);
  void add_constraint(int row, int col, double val);
  void update();
  void solve();
  double get_value(int col) const;

private:
  IPimpl* impl_;
};

#endif  // __INC_IP_H__
