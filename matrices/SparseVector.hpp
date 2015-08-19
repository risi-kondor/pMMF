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



#ifndef _SparseVector
#define _SparseVector
#include "Vector.hpp"

struct SVpair {
	SVpair() = default;
	SVpair(const INDEX& _first, const FIELD& _second) :
			first(_first), second(_second) {
	}
	;
	string str() {
		ostringstream result;
		result << "(" << first << "," << second << ")";
		return result.str();
	}
	INDEX first;
	FIELD second;
};

struct {
	bool operator()(const SVpair& a, const SVpair& b) {
		return a.first < b.first;
	}
} SVpairComparator;

class SparseVector: public Vector {
public:
	SparseVector(const int _n) :
			Vector(_n) {
	}

public:
	virtual FIELD& operator()(const int i)=0;
	virtual FIELD operator()(const int i) const=0;
	virtual FIELD read(const int i) const=0;

	virtual void insert(const int i, const FIELD value)=0;
	virtual void append(const int i, const FIELD value)=0;
	virtual void zero(const int i)=0;

	virtual void sort()=0;
	virtual void tidy()=0;

public:

public:

	virtual void apply(const GivensRotation& Q)=0;
	virtual void apply(const KpointRotation& Q)=0;
	virtual void applyInverse(const GivensRotation& Q)=0;
	virtual void applyInverse(const KpointRotation& Q)=0;

public:

};

#endif
