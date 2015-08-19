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



#ifndef _Matrix
#define _Matrix

#include "MMFbase.hpp"
#include <functional>

class Matrix {
public:
	virtual ~Matrix() {
	}

public:
	// constructors
	Matrix(const int _nrows, const int _ncols) :
			nrows(_nrows), ncols(_ncols) {
	}

public:
	// member access
	virtual FIELD& operator()(const int i, const int j)=0;
	virtual FIELD operator()(const int i, const int j) const=0;
	virtual FIELD read(const int i, const int j) const=0;
	virtual bool isFilled(const int i, const int j) const=0;
	virtual int nFilled() const=0;
	virtual bool isSparse() const=0;

	virtual void foreach(std::function<void(INDEX, INDEX, FIELD&)> lambda)=0;
	virtual void foreach_in_column(const int j,
			std::function<void(INDEX, FIELD&)> lambda)=0;

public:
	// scalar valued operations
	virtual int nnz() const=0;

public:
	virtual Cmatrix dot(const Matrix& x) const; // {};

public:
	virtual void applyFromLeft(const GivensRotation& Q)=0;
	virtual void applyFromLeft(const KpointRotation& Q)=0;
	virtual void applyFromLeftT(const GivensRotation& Q)=0;
	virtual void applyFromLeftT(const KpointRotation& Q)=0;

	virtual void applyFromRight(const GivensRotation& Q)=0;
	virtual void applyFromRight(const KpointRotation& Q)=0;
	virtual void applyFromRightT(const GivensRotation& Q)=0;
	virtual void applyFromRightT(const KpointRotation& Q)=0;

public:
	virtual string str(const Dense dummy) const;
	virtual string str(const Sparse dummy) const;
	virtual string str() const {
		return str(Dense());
	}

public:
	int nrows;
	int ncols;

};

ostream& operator<<(ostream& stream, const Matrix& x);

class SparseMatrix: public Matrix {
public:
	using Matrix::Matrix;
	bool isSparse() const {
		return true;
	}
};

class DenseMatrix: public Matrix {
public:
	using Matrix::Matrix;
	bool isFilled(const int i, const int j) const {
		return true;
	}
	int nFilled() const {
		return nrows * ncols;
	}
	bool isSparse() const {
		return false;
	}
};
#endif 
