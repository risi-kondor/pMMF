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



#ifndef _DenseMatrixFile
#define _DenseMatrixFile 

#include "MatrixFile.hpp"

class DenseMatrixFile: public MatrixFile {
public:
	class ASCII;
	class Matlab;

	class iterator;
	class eof {
	};
public:
	virtual iterator begin()=0;
	virtual void get(iterator& it)=0;
	iterator end();
	virtual DenseMatrixFile& operator>>(FIELD& dest)=0;

public:
};

class DenseMatrixFile::iterator {
public:

	iterator(DenseMatrixFile& _owner) :
			owner(_owner) {
		owner.get(*this);
	}
	iterator(DenseMatrixFile& _owner, const eof& dummy) :
			owner(_owner) {
		endflag = true;
	}

	FIELD operator*() {
		return value;
	}
	FIELD* operator->() {
		return &value;
	}
	iterator& operator++() {
		owner.get(*this);
		return *this;
	}
	bool operator!=(const iterator& x) {
		return endflag != x.endflag;
	}

public:
	FIELD value;
	bool endflag;
	DenseMatrixFile& owner;
};
#endif
