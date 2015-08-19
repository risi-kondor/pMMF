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



#ifndef _Cmatrix
#define _Cmatrix
#include "Matrix.hpp"
#include "Cvector.hpp"
#include "GivensRotation.hpp"
#include "KpointRotation.hpp"

class Cmatrix: public DenseMatrix, public Serializable {
public:
	Cmatrix(const Cmatrix& x);
	Cmatrix(Cmatrix&& x);
	Cmatrix& operator=(const Cmatrix& x);
	Cmatrix& operator=(Cmatrix&& x);
	~Cmatrix() {
		delete[] array;
	}

public:
	Cmatrix() :
			Cmatrix(0, 0) {
	}
	Cmatrix(const int _nrows, const int _ncols) :
			DenseMatrix(_nrows, _ncols) {
		array = new FIELD[nrows * ncols];
	}
	Cmatrix(const int _nrows, const int _ncols, const Zero dummy) : // obsolete
			Cmatrix(_nrows, _ncols) {
		for (int i = 0; i < nrows * ncols; i++)
			array[i] = 0;
	}

public:
	// named constructors
	static Cmatrix Zero(const int _nrows, const int _ncols) {
		Cmatrix M(_nrows, _ncols);
		for (int i = 0; i < _nrows * _ncols; i++)
			M.array[i] = 0;
		return M;
	}

	static Cmatrix Identity(const int n) {
		Cmatrix M = Zero(n, n);
		for (int i = 0; i < n; i++)
			M.array[i * n + i] = 1;
		return M;
	}

	static Cmatrix Kronecker(const Cmatrix& seed, const int nlevels);
	static Cmatrix Kronecker(const Cmatrix& X1, const Cmatrix& X2);

	static Cmatrix Random(const int _nrows, const int _ncols);
	static Cmatrix RandomSymmetric(const int _nrows);
	static Cmatrix RandomPosDef(const int n);

public:
	// converters
	Cmatrix(const Cmatrix& x, const Remap& remap1, const Remap& remap2);
	Cmatrix(const Cmatrix& x, const Transpose& dummy);

	template<class COLUMNTYPE>
	Cmatrix(const MatrixX<COLUMNTYPE>& x);
	Cmatrix(const EigenMatrixXdAdaptor& M);
	template<class TYPE> TYPE convert() const;

public:
	// element access and views
	FIELD& operator()(const int i, const int j) {
		assert(i < nrows);
		assert(j < ncols);
		return array[j * nrows + i];
	}
	FIELD operator()(const int i, const int j) const {
		assert(i < nrows);
		assert(j < ncols);
		return array[j * nrows + i];
	}
	FIELD read(const int i, const int j) const {
		assert(i < nrows);
		assert(j < ncols);
		return array[j * nrows + i];
	}

	void foreach(std::function<void(INDEX, INDEX, FIELD&)> lambda) {
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				lambda(i, j, array[j * nrows + i]);
	}

	void foreach_in_column(const int j,
			std::function<void(INDEX, FIELD&)> lambda) {
		for (int i = 0; i < nrows; i++)
			lambda(i, array[j * nrows + i]);
	}

	Cvector::Virtual vcolumn(const int j) const {
		return Cvector::Virtual(nrows, &array[j * nrows]);
	}

	Cmatrix submatrix(const int _nrows, const int _ncols, const int ioffset = 0,
			const int joffset = 0) const;
	Cmatrix submatrix(const IndexSet& I1, const IndexSet& I2) const {
		const int k1 = I1.k;
		const int k2 = I2.k;
		Cmatrix R(k1, k2);
		for (int i = 0; i < k1; i++)
			for (int j = 0; j < k2; j++)
				R(i, j) = array[I2[j] * nrows + I1[i]];
		return R;
	}

	Cvector diag() const {
		Cvector diagonal(min(nrows, ncols));
		for (int j = 0; j < min(nrows, ncols); j++)
			diagonal(j) = array[j * (nrows + 1)];
		return diagonal;
	}

public:
	// scalar-valued operations
	int nnz() const {
		int t = 0;
		for (int i = 0; i < nrows * ncols; i++)
			if (array[i] != 0)
				t++;
		return t;
	}

	FIELD norm2() const {
		FIELD t = 0;
		for (int i = 0; i < nrows * ncols; i++)
			t += array[i] * array[i];
		return t;
	}

	FIELD diff2(const Cmatrix& X) const {
		assert(X.nrows == nrows);
		assert(X.ncols == ncols);
		FIELD t = 0;
		for (int i = 0; i < nrows * ncols; i++)
			t += (array[i] - X.array[i]) * (array[i] - X.array[i]);
		return t;
	}

	FIELD columnsum(const int j) const {
		FIELD t = 0;
		for (int i = 0; i < nrows; i++)
			t += array[j * nrows + i];
		return t;
	}

public:
	// vector valued operations
	Cvector operator*(const Cvector& x) const {
		Cvector y = Cvector::Zero(nrows);
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				y.array[i] += array[j * nrows + i] * x.array[j];
		return y;
	}

	Cvector dot(const Cvector& v) const {
		assert(v.n == nrows);
		Cvector r(ncols);
		for (int j = 0; j < ncols; j++) {
			FIELD t = 0;
			for (int k = 0; k < nrows; k++)
				t += array[j * nrows + k] * v.array[k];
			r.array[j] = t;
		}
		return r;
	}

public:
	// matrix valued operations
	Cmatrix operator*(const Cmatrix& X) const {
		Cmatrix M = Zero(nrows, X.ncols);
		assert(ncols == X.nrows);
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < X.ncols; j++)
				for (int k = 0; k < ncols; k++)
					M.array[j * nrows + i] += array[k * nrows + i]
							* X.array[j * X.nrows + k];
		return M;
	}

	Cmatrix operator-(const Matrix& x) const {
		Cmatrix R(nrows, ncols);
		assert(x.nrows == nrows);
		assert(x.ncols == ncols);
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				R(i, j) = (*this)(i, j) - x(i, j);
		return R;
	}

	Cmatrix dot(const Cmatrix& x) const {
		Cmatrix M(ncols, x.ncols);
		if (this == &x)
			for (int i = 0; i < ncols; i++)
				for (int j = 0; j <= i; j++) {
					FIELD t = 0;
					for (int k = 0; k < nrows; k++)
						t += array[i * nrows + k] * x.array[j * x.nrows + k];
					M(i, j) = t;
					M(j, i) = t;
				}
		else
			for (int i = 0; i < ncols; i++)
				for (int j = 0; j < x.ncols; j++) {
					FIELD t = 0;
					for (int k = 0; k < nrows; k++)
						t += array[i * nrows + k] * x.array[j * x.nrows + k];
					M(i, j) = t;
				}
		return M;
	}

public:
	// in-place operations
	Cmatrix& operator+=(const Cmatrix& x) {
		assert(x.nrows == nrows);
		assert(x.ncols == ncols);
		for (int i = 0; i < nrows * ncols; i++)
			array[i] += x.array[i];
		return *this;
	}

	Cmatrix& multiplyRowsBy(const Cvector& v) {
		assert(v.n == nrows);
		for (int j = 0; j < ncols; j++)
			for (int i = 0; i < nrows; i++)
				array[j * nrows + i] = v.array[i] * array[j * nrows + i];
		return *this;
	}

	Cmatrix& multiplyColsBy(const Cvector& v) {
		assert(v.n == ncols);
		for (int j = 0; j < ncols; j++)
			for (int i = 0; i < nrows; i++)
				array[j * nrows + i] = v.array[j] * array[j * nrows + i];
		return *this;
	}

	Cmatrix& divideRowsBy(const Cvector& v) {
		assert(v.n == nrows);
		for (int j = 0; j < ncols; j++)
			for (int i = 0; i < nrows; i++)
				array[j * nrows + i] = array[j * nrows + i] / v.array[i];
		return *this;
	}

	Cmatrix& divideColsBy(const Cvector& v) {
		assert(v.n == ncols);
		for (int j = 0; j < ncols; j++)
			for (int i = 0; i < nrows; i++)
				array[j * nrows + i] = array[j * nrows + i] / v.array[j];
		return *this;
	}

public:
	// solvers
	Cmatrix inverse() const;
	Cmatrix inverseSqrt() const;
	Cvector solve(const Cvector& b) const;
	pair<Cmatrix*, Cvector*> symmetricEigensolver() const;

public:
	// rotations
	void applyFromLeft(const GivensRotation& Q) {
		assert(Q.i1 < nrows);
		assert(Q.i2 < nrows);
		for (int j = 0; j < ncols; j++) {
			FIELD v1 = array[j * nrows + Q.i1];
			FIELD v2 = array[j * nrows + Q.i2];
			array[j * nrows + Q.i1] = Q.cos * v1 - Q.sin * v2;
			array[j * nrows + Q.i2] = Q.cos * v2 + Q.sin * v1;
		}
	}

	void applyFromLeft(const KpointRotation& Q) {
		FIELD temp[Q.k];
		for (int i = 0; i < Q.k; i++)
			assert(Q.ix[i] < nrows);
		for (int j = 0; j < ncols; j++) {
			for (int i = 0; i < Q.k; i++)
				temp[i] = array[j * nrows + Q.ix[i]];
			for (int i = 0; i < Q.k; i++) {
				FIELD t = 0;
				for (int l = 0; l < Q.k; l++)
					t += temp[l] * Q.q[l * Q.k + i];
				array[j * nrows + Q.ix[i]] = t;
			}
		}
	}

	void applyFromLeftT(const GivensRotation& Q) {
		assert(Q.i1 < nrows);
		assert(Q.i2 < nrows);
		for (int j = 0; j < ncols; j++) {
			FIELD v1 = array[j * nrows + Q.i1];
			FIELD v2 = array[j * nrows + Q.i2];
			array[j * nrows + Q.i1] = Q.cos * v1 + Q.sin * v2;
			array[j * nrows + Q.i2] = Q.cos * v2 - Q.sin * v1;
		}
	}

	void applyFromLeftT(const KpointRotation& Q) {
		FIELD temp[Q.k];
		for (int i = 0; i < Q.k; i++)
			assert(Q.ix[i] < nrows);
		for (int j = 0; j < ncols; j++) {
			for (int i = 0; i < Q.k; i++)
				temp[i] = array[j * nrows + Q.ix[i]];
			for (int i = 0; i < Q.k; i++) {
				FIELD t = 0;
				for (int l = 0; l < Q.k; l++)
					t += temp[l] * Q.q[i * Q.k + l];
				array[j * nrows + Q.ix[i]] = t;
			}
		}
	}

	void applyFromRightT(const GivensRotation& Q) {
		assert(Q.i1 < ncols);
		assert(Q.i2 < ncols);
		for (int i = 0; i < nrows; i++) {
			FIELD v1 = array[Q.i1 * nrows + i];
			FIELD v2 = array[Q.i2 * nrows + i];
			array[Q.i1 * nrows + i] = Q.cos * v1 - Q.sin * v2;
			array[Q.i2 * nrows + i] = Q.cos * v2 + Q.sin * v1;
		}
	}

	void applyFromRightT(const KpointRotation& Q) {
		FIELD temp[Q.k];
		for (int i = 0; i < Q.k; i++)
			assert(Q.ix[i] < ncols);
		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < Q.k; j++)
				temp[j] = array[Q.ix[j] * nrows + i];
			for (int j = 0; j < Q.k; j++) {
				FIELD t = 0;
				for (int l = 0; l < Q.k; l++)
					t += temp[l] * Q.q[l * Q.k + j];
				array[Q.ix[j] * nrows + i] = t;
			}
		}
	}

	void applyFromRight(const GivensRotation& Q) {
		assert(Q.i1 < ncols);
		assert(Q.i2 < ncols);
		for (int i = 0; i < nrows; i++) {
			FIELD v1 = array[Q.i1 * nrows + i];
			FIELD v2 = array[Q.i2 * nrows + i];
			array[Q.i1 * nrows + i] = Q.cos * v1 + Q.sin * v2;
			array[Q.i2 * nrows + i] = Q.cos * v2 - Q.sin * v1;
		}
	}

	void applyFromRight(const KpointRotation& Q) {
		FIELD temp[Q.k];
		for (int i = 0; i < Q.k; i++)
			assert(Q.ix[i] < ncols);
		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < Q.k; j++)
				temp[j] = array[Q.ix[j] * nrows + i];
			for (int j = 0; j < Q.k; j++) {
				FIELD t = 0;
				for (int l = 0; l < Q.k; l++)
					t += temp[l] * Q.q[j * Q.k + l];
				array[Q.ix[j] * nrows + i] = t;
			}
		}
	}

public:
	// I/O
	string str(const Dense dummy) const;
	string str(const Sparse dummy) const;
	string str() const;

	Cmatrix(DenseMatrixFile& file);
	Cmatrix(SparseMatrixFile& file);

	Cmatrix(Bifstream& ifs);
	void serialize(Bofstream& ofs) const;
	void serialize(Rstream& rstream) const;

public:
	FIELD* array;
};
ostream& operator<<(ostream& stream, const Cmatrix& v);

// ---- Template methods -----------------------------------------------------------------------------------------
template<class COLUMNTYPE>
Cmatrix::Cmatrix(const MatrixX<COLUMNTYPE>& x) :
		Cmatrix(x.nrows, x.ncols) {
	for (int i = 0; i < nrows * ncols; i++)
		array[i] = 0;
	cout << "Warning!!!!" << endl;
	for (int j = 0; j < ncols; j++)
		for (auto& p : *x.column[j])
			array[j * nrows + p.first] = p.second;
}

#endif 

