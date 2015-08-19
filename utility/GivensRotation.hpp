/* 
    pMMF is an open source software library for efficiently computing
    the MMF factorization of large matrices.

    Copyright (C) 2015 Risi Kondor, Nedelina Teneva, Pramod Mudrakarta
    
    This file is part of pMMF.

    pMMF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#ifndef _GivensRotation
#define _GivensRotation

#include "ElementaryRotation.hpp" 
#include "Rstream.hpp"


class GivensRotation: public ElementaryRotation{
public:

  GivensRotation(int _i1, int _i2, double _cos, double _sin){
    if(_i1<_i2){i1=_i1;i2=_i2;cos=_cos;sin=_sin;}
    else{i1=_i2;i2=_i1;cos=_cos;sin=-_sin;}}

  GivensRotation(int _i1, int _i2, double _alpha);
  
  GivensRotation(const Random dummy, const int n);

public:
  
  string str() const;
  void serialize(Rstream& s) const;
  void serialize(Bofstream& ofs) const;
  void save(const char* filename) {Bofstream ofs(filename); serialize(ofs);}
  GivensRotation(Bifstream& ifs);

public:

  int i1;
  int i2;
  double cos;
  double sin;

  // static default_random_engine randomNumberGenerator;

};



#endif
