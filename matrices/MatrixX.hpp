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



#ifndef _MatrixX
#define _MatrixX

#include <random>
#include "Matrix.hpp"
#include "Cmatrix.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "SparseMatrixFile.hpp"
#include "DenseMatrixFile.hpp"
extern default_random_engine randomNumberGenerator;

template<class COLUMNTYPE>
class MatrixX: public SparseMatrix {
public:

	MatrixX(const MatrixX<COLUMNTYPE>& x);
	MatrixX(MatrixX<COLUMNTYPE> && x);
	MatrixX<COLUMNTYPE>& operator=(const MatrixX<COLUMNTYPE>& x);
	MatrixX<COLUMNTYPE>& operator=(MatrixX<COLUMNTYPE> && x);
	~MatrixX() {
		for (int i = 0; i < column.size(); i++)
			delete column[i];
	}

public:
	// constructors
	MatrixX() :
			MatrixX(0, 0) {
	}
	;
	MatrixX(const int _nrows, const int _ncols) :
			SparseMatrix(_nrows, _ncols), column(ncols) {
		for (int i = 0; i < ncols; i++)
			column[i] = new COLUMNTYPE(nrows);
	}
	//MatrixX(const int _nrows, const int _ncols, const Zero& dummy): MatrixX(_nrows,_ncols){}
	//MatrixX(const int _nrows, const int _ncols, const Uninitialized& dummy): SparseMatrix(_nrows,_ncols){}
	MatrixX(const int _nrows, const int _ncols, const Nullcols& dummy) :
			SparseMatrix(_nrows, _ncols), column(ncols, nullptr) {
	}

	MatrixX(const int _nrows, const int _ncols, const Random& dummy); // depricated
	MatrixX(const int _nrows, const class Identity& dummy); // obsolete
	MatrixX(const int _nrows, const class Random& dummy); // obsolete

public:
	// named constructors
	static MatrixX Zero(const int _nrows, const int _ncols) {
		return MatrixX(_nrows, _ncols);
	}
	static MatrixX Identity(const int n);
	static MatrixX Random(const int _nrows, const int _ncols, const double p =
			0.5);
	static MatrixX RandomSymmetric(const int n, const double p = 0.5);

public:
	// conversions
	template<class COLUMNTYPE2>
	MatrixX(const MatrixX<COLUMNTYPE2>& x);
	template<class COLUMNTYPE2>
	MatrixX(const MatrixX<COLUMNTYPE2>& x, const int _nrows, const int _ncols);
	template<class COLUMNTYPE2>
	MatrixX(const MatrixX<COLUMNTYPE2>& x, const Transpose& dummy);
	template<class COLUMNTYPE2>
	MatrixX(const MatrixX<COLUMNTYPE2>& x, const Remap& remap1,
			const Remap& remap2) :
			MatrixX(x, remap1, remap2, x.nrows, x.ncols) {
	}
	MatrixX(const Cmatrix& x);

private:
	// remapping
	template<class COLUMNTYPE2>
	MatrixX(const MatrixX<COLUMNTYPE2>& x, const Remap& remap1,
			const Remap& remap2, const int _nrows, const int _ncols);

public:
	// element access and views
	FIELD& operator()(const int i, const int j) {
		assert(i < nrows);
		assert(j < ncols);
		return (*column[j])(i);
	}
	FIELD operator()(const int i, const int j) const {
		assert(i < nrows);
		assert(j < ncols);
		return (*column[j])(i);
	}
	FIELD read(const int i, const int j) const {
		assert(i < nrows);
		assert(j < ncols);
		return column[j]->read(i);
	}

	COLUMNTYPE& vcolumn(const int j) {
		assert(j < ncols);
		return *column[j];
	}

	void foreach(std::function<void(INDEX, INDEX, FIELD&)> lambda) {
		for (int j = 0; j < ncols; j++)
			for (auto& p : *column[j])
				lambda(p.first, j, p.second);
	}

	void foreach_in_column(const int j,
			std::function<void(INDEX, FIELD&)> lambda) {
		for (auto& p : *column[j])
			lambda(p.first, p.second);
	}

	bool isFilled(const int i, const int j) const {
		assert(j < ncols);
		return column[j]->isFilled(i);
	}

	int nFilled() const {
		int t = 0;
		for (int j = 0; j < ncols; j++)
			t += column[j]->nFilled();
		return t;
	}

	COLUMNTYPE diag() const {
		COLUMNTYPE diagonal(min(nrows, ncols));
		for (int j = 0; j < min(nrows, ncols); j++)
			diagonal(j) = (*column[j])(j);
		return (const COLUMNTYPE) diagonal;
	}

public:
	// scalar valued operations
	int nnz() const;
	FIELD norm2() const;
	FIELD diff2(const MatrixX<COLUMNTYPE>& X) const;

public:
	// vector valued operations
	Cvector operator*(const Cvector& x) const {
		Cvector r(nrows, Zero());
		assert(x.n == ncols);
		for (int j = 0; j < ncols; j++) {
			FIELD t = x.array[j];
			for (auto& p : *column[j])
				r(p.first) += p.second * t;
		}
		return r;
	}

	Cvector dot(const COLUMNTYPE& v) const {
		assert(v.n == nrows);
		Cvector r(ncols);
		for (int j = 0; j < ncols; j++)
			r.array[j] = column[j]->dot(v);
		return r;
	}

public:
	// matrix valued operations
	Cmatrix operator*(const MatrixX<COLUMNTYPE>& x) const {
		Cmatrix M(nrows, x.ncols);
		assert(ncols == nrows);
		return M;
	}

	Cmatrix dot(const MatrixX<COLUMNTYPE>& x) const {
		Cmatrix M(ncols, x.ncols);
		assert(x.nrows == nrows);
		if (this == &x)
			for (int i = 0; i < ncols; i++)
				for (int j = 0; j <= i; j++) {
					M(i, j) = column[i]->dot(*x.column[j]);
					M(j, i) = M(i, j);
				}
		else
			for (int i = 0; i < ncols; i++)
				for (int j = 0; j < x.ncols; j++)
					M(i, j) = column[i]->dot(*x.column[j]);
		return M;
	}

public:
	// in place operations
	MatrixX<COLUMNTYPE>& multiplyRowsBy(const Cvector& v) {
		for (int j = 0; j < ncols; j++)
			(*column[j]) *= v;
		return *this;
	}

	MatrixX<COLUMNTYPE>& multiplyColsBy(const Cvector& v) {
		assert(v.n == ncols);
		for (int j = 0; j < ncols; j++)
			(*column[j]) *= v(j);
		return *this;
	}

	MatrixX<COLUMNTYPE>& divideRowsBy(const Cvector& v) {
		for (int j = 0; j < ncols; j++)
			(*column[j]) /= v;
		return *this;
	}

	MatrixX<COLUMNTYPE>& divideColsBy(const Cvector& v) {
		assert(v.n == ncols);
		for (int j = 0; j < ncols; j++)
			(*column[j]) *= (1.0 / v(j));
		return *this;
	}

	void transpose();

	MatrixX<COLUMNTYPE>& tidy() {
		for (int j = 0; j < ncols; j++)
			column[j]->tidy();
		return *this;
	}

public:
	// rotations
	void applyFromLeft(const GivensRotation& Q) {
		for (int i = 0; i < ncols; i++)
			column[i]->apply(Q);
	}

	void applyFromLeft(const KpointRotation& Q) {
		for (int i = 0; i < ncols; i++)
			column[i]->apply(Q);
	}

	void applyFromLeftT(const GivensRotation& Q) {
		for (int i = 0; i < ncols; i++)
			column[i]->applyInverse(Q);
	}

	void applyFromLeftT(const KpointRotation& Q) {
		for (int i = 0; i < ncols; i++)
			column[i]->applyInverse(Q);
	}

	void applyFromRightT(const GivensRotation& Q) {
		COLUMNTYPE* newcol1 = new COLUMNTYPE(nrows);
		COLUMNTYPE* newcol2 = new COLUMNTYPE(nrows);
		assert(Q.i1 < ncols);
		assert(Q.i2 < ncols);
		newcol1->add(*column[Q.i1], Q.cos);
		newcol1->add(*column[Q.i2], -Q.sin);
		newcol2->add(*column[Q.i2], Q.cos);
		newcol2->add(*column[Q.i1], Q.sin);
		delete column[Q.i1];
		column[Q.i1] = newcol1;
		delete column[Q.i2];
		column[Q.i2] = newcol2;
	}

	void applyFromRightT(const KpointRotation& Q) {
		COLUMNTYPE* col[Q.k];
		for (int j = 0; j < Q.k; j++) {
			assert(Q.ix[j] < ncols);
			col[j] = new COLUMNTYPE(nrows);
			for (int i = 0; i < Q.k; i++)
				col[j]->add(*column[Q.ix[i]], Q.q[i * Q.k + j]);
		}
		for (int j = 0; j < Q.k; j++) {
			delete column[Q.ix[j]];
			column[Q.ix[j]] = col[j];
		}
	}

	void applyFromRight(const GivensRotation& Q) {
		COLUMNTYPE* newcol1 = new COLUMNTYPE(nrows);
		COLUMNTYPE* newcol2 = new COLUMNTYPE(nrows);
		assert(Q.i1 < ncols);
		assert(Q.i2 < ncols);
		newcol1->add(*column[Q.i1], Q.cos);
		newcol1->add(*column[Q.i2], Q.sin);
		newcol2->add(*column[Q.i2], Q.cos);
		newcol2->add(*column[Q.i1], -Q.sin);
		delete column[Q.i1];
		column[Q.i1] = newcol1;
		delete column[Q.i2];
		column[Q.i2] = newcol2;
	}

	void applyFromRight(const KpointRotation& Q) {
		COLUMNTYPE* col[Q.k];
		for (int j = 0; j < Q.k; j++) {
			assert(Q.ix[j] < ncols);
			col[j] = new COLUMNTYPE(nrows);
			for (int i = 0; i < Q.k; i++)
				col[j]->add(*column[Q.ix[i]], Q.q[j * Q.k + i]);
		}
		for (int j = 0; j < Q.k; j++) {
			delete column[Q.ix[j]];
			column[Q.ix[j]] = col[j];
		}
	}

public:
	// I/O
	string str(const Sparse dummy) const;
	string str(const Dense dummy) const {
		return Matrix::str(Dense());
	}
	string str() const {
		return str(Sparse());
	}
	MatrixX(DenseMatrixFile& file);
	MatrixX(SparseMatrixFile& file);

public:
	vector<COLUMNTYPE*> column;
};

// ---- Copying ---------------------------------------------------------------------------------------------------
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE>& x) :
		SparseMatrix(x.nrows, x.ncols), column(x.ncols) {
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: MatrixX copied." << endl;
#endif
	for (int j = 0; j < ncols; j++)
		column[j] = new COLUMNTYPE(*x.column[j]);
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(MatrixX<COLUMNTYPE> && x) :
		SparseMatrix(x.nrows, x.ncols), column(x.ncols) {
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: MatrixX moved." << endl;
#endif
	for (int j = 0; j < ncols; j++) {
		column[j] = x.column[j];
		x.column[j] = nullptr;
	}
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>& MatrixX<COLUMNTYPE>::operator=(
		const MatrixX<COLUMNTYPE>& x) {
	nrows = x.nrows;
	ncols = x.ncols;
	for (auto p : column)
		delete p;
	column.resize(ncols);
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: MatrixX assigned." << endl;
#endif
	for (int j = 0; j < ncols; j++)
		column[j] = new COLUMNTYPE(*x.column[j]);
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>& MatrixX<COLUMNTYPE>::operator=(MatrixX<COLUMNTYPE> && x) {
	nrows = x.nrows;
	ncols = x.ncols;
	for (auto p : column)
		delete p;
	column.resize(ncols);
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: MatrixX move-assigned." << endl;
#endif
	for (int j = 0; j < ncols; j++) {
		column[j] = x.column[j];
		x.column[j] = nullptr;
	}
}

// ---- Constructors ---------------------------------------------------------------------------------------------
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::Identity(const int n) {
	MatrixX M(n, n);
	for (int i = 0; i < n; i++)
		M(i, i) = 1;
	return M;
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const int _nrows, const int _ncols,
		const class Random& random) :
		MatrixX(_nrows, _ncols) {
	double p = random.p;
	uniform_int_distribution<int> distri(0, nrows - 1);
	uniform_int_distribution<int> distrj(0, ncols - 1);
	uniform_real_distribution<FIELD> distr(0, 1);
	for (int i = 0; i < p * nrows * ncols; i++)
		(*this)(distri(randomNumberGenerator), distrj(randomNumberGenerator)) =
				distr(randomNumberGenerator);

}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::Random(const int _nrows,
		const int _ncols, const double p) {
	MatrixX M(_nrows, _ncols);
	uniform_int_distribution<int> distri(0, _nrows - 1);
	uniform_int_distribution<int> distrj(0, _ncols - 1);
	uniform_real_distribution<FIELD> distr(0, 1);
	for (int i = 0; i < p * _nrows * _ncols; i++)
		M(distri(randomNumberGenerator), distrj(randomNumberGenerator)) = distr(
				randomNumberGenerator);
	return M;
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::RandomSymmetric(const int n,
		const double p) {
	MatrixX M(n, n);
	uniform_int_distribution<int> distri(0, n - 1);
	uniform_real_distribution<FIELD> distr(0, 1);
	for (int u = 0; u < p * n * n / 2; u++) {
		int i = distri(randomNumberGenerator);
		int j = distri(randomNumberGenerator);
		FIELD v = distr(randomNumberGenerator);
		M(i, j) = v;
		M(j, i) = v;
	}
	return M;
}

// ---- Conversions -----------------------------------------------------------------------------------------------
template<class COLUMNTYPE>
template<class COLUMNTYPE2>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x) :
		SparseMatrix(x.nrows, x.ncols), column(x.ncols) {
	for (int j = 0; j < ncols; j++)
		column[j] = new COLUMNTYPE(*x.column[j]);
	cout << "Warning: MatrixX convert-copied." << endl;
}

template<class COLUMNTYPE>
template<class COLUMNTYPE2>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x, const int _nrows,
		const int _ncols) :
		SparseMatrix(_nrows, _ncols), column(_ncols) {
	for (int j = 0; j < ncols; j++) {
		column[j] = new COLUMNTYPE(nrows);
		for (auto& p : *x.column[j])
			if (p.first < nrows)
				column[j]->insert(p.first, p.second);
	}
}

template<class COLUMNTYPE>
template<class COLUMNTYPE2>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x,
		const Transpose& dummy) :
		SparseMatrix(x.ncols, x.nrows), column(x.nrows) {
	for (int j = 0; j < ncols; j++)
		column[j] = new COLUMNTYPE(ncols);
	for (int i = 0; i < nrows; i++) {
		x.column[i]->sort();
		for (auto& p : *x.column[i])
			if (p.second != 0)
				column[p.first]->append(i, p.second);
	}
}

template<class COLUMNTYPE>
template<class COLUMNTYPE2> // Warning: won't work for mapping unsorted v or h to l 
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x, const Remap& remap1,
		const Remap& remap2, const int _nrows, const int _ncols) :
		MatrixX(_nrows, _ncols) {
	for (int j = 0; j < x.ncols; j++) {
		const COLUMNTYPE& xcol = *x.column[j];
		COLUMNTYPE& col = *column[remap2.forward[j]];
		for (auto& p : xcol)
			col.append(remap1.forward[p.first], p.second);
	}
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const class Cmatrix& x) :
		SparseMatrix(x.nrows, x.ncols), column(x.ncols) {
	for (int j = 0; j < ncols; j++)
		column[j] = new COLUMNTYPE(x.vcolumn(j));
	cout << "WARNING: dense->sparse conversion" << endl;
}

// ---- Operations ----------------------------------------------------------------------------------------------
template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::transpose() {
	vector<COLUMNTYPE*> newcolumn(nrows);
	for (int i = 0; i < nrows; i++)
		newcolumn[i] = new COLUMNTYPE(ncols);
	for (int i = 0; i < ncols; i++) {
		for (auto& p : *column[i])
			if (p.second != 0)
				newcolumn[p.first]->append(i, p.second);
		delete column[i];
	}
	column = newcolumn;
	swapp(nrows, ncols);
}

// ---- Scalar valued operations ---------------------------------------------------------------------------------
template<class COLUMNTYPE>
FIELD MatrixX<COLUMNTYPE>::norm2() const {
	FIELD t = 0;
	for (int j = 0; j < ncols; j++)
		t += column[j]->norm2();
	return t;
}

template<class COLUMNTYPE>
FIELD MatrixX<COLUMNTYPE>::diff2(const MatrixX<COLUMNTYPE>& X) const {
	assert(X.ncols == ncols);
	FIELD t = 0;
	for (int j = 0; j < ncols; j++)
		t += column[j]->diff2(*X.column[j]);
	return t;
}

template<class COLUMNTYPE>
int MatrixX<COLUMNTYPE>::nnz() const {
	int result = 0;
	for (int j = 0; j < ncols; j++)
		result += column[j]->nnz();
	return result;
}

// ---- I/O --------------------------------------------------------------------------------------------------------
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(SparseMatrixFile& file) :
		MatrixX(file.nrows, file.ncols) {
	IndexValueTriple p;
	file >> p;
	while (p.i != -1) {
		column[p.j]->insert(p.i, p.value);
		file >> p;
	}
	tidy();
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(DenseMatrixFile& file) :
		MatrixX(file.nrows, file.ncols) {
	auto it = file.begin();
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++) {
			FIELD v;
			file >> v;
			column[j]->append(i, v);
		}
}

template<class COLUMNTYPE>
string MatrixX<COLUMNTYPE>::str(const Sparse dummy) const {
	ostringstream stream;
	stream.precision(5);
	for (int j = 0; j < ncols; j++)
		for (auto& it : *column[j])
			stream << "(" << it.first << "," << j << ") : " << it.second
					<< endl;
	return stream.str();
}

template<class COLUMNTYPE>
ostream& operator<<(ostream& stream, const MatrixX<COLUMNTYPE>& x) {
	stream << x.str(Sparse());
	return stream;
}

#endif 
