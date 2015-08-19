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



#ifndef _Vector
#define _Vector

#include "MMFbase.hpp"
#include <functional>

class Vector {
public:

	Vector(const int _n) :
			n(_n) {
	}

public:
	virtual FIELD& operator()(const int n)=0;
	virtual FIELD operator()(const int n) const=0;
	virtual FIELD read(const int n) const {
		return (*this)(n);
	}

	virtual void foreach(std::function<void(INDEX, FIELD&)> lambda)=0;

	virtual bool isFilled(const int i) const =0;
	virtual int nFilled() const=0;

public:
	virtual int nnz() const=0;
	virtual int argmax() const=0;
	virtual int argmax_abs() const=0;
	virtual FIELD norm2() const=0;

public:
	virtual void apply(const GivensRotation& Q)=0;
	virtual void apply(const KpointRotation& Q)=0;
	virtual void applyInverse(const GivensRotation& Q)=0;
	virtual void applyInverse(const KpointRotation& Q)=0;

public:
	virtual string str(const Dense dummy) const;
	virtual string str(const Sparse dummy) const;
	virtual string str() const {
		return str(Dense());
	}
	;

public:
	int n;
};

ostream& operator<<(ostream& stream, const Vector& x);
#endif
