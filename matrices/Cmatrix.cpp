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



#include "Cmatrix.hpp"
#include "MatrixX.hpp"
#include "DenseMatrixFile.hpp"
#include "SparseMatrixFile.hpp"
#include "EigenAdaptors.hpp"
extern default_random_engine randomNumberGenerator;

// ---- Copying ---------------------------------------------------------------------------------------------------
Cmatrix::Cmatrix(const Cmatrix& x) :
		Cmatrix(x.nrows, x.ncols) {
	//std::copy(x.array,x.array+nrows*ncols,array);
	for (int i = 0; i < nrows * ncols; i++)
		array[i] = x.array[i];
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: Cmatrix copied." << endl;
#endif
}

Cmatrix::Cmatrix(Cmatrix&& x) :
		DenseMatrix(x.nrows, x.ncols) {
	array = x.array;
	x.array = nullptr;
	x.nrows = 0;
	x.ncols = 0;
	cout << "ssss" << endl;
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: Cmatrix moved." << endl;
#endif
}

Cmatrix& Cmatrix::operator=(const Cmatrix& x) {
	delete[] array;
	nrows = x.nrows;
	ncols = x.ncols;
	array = new FIELD[nrows * ncols];
	for (int i = 0; i < nrows * ncols; i++)
		array[i] = x.array[i];
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: Cmatrix assigned." << endl;
#endif
	return *this;
}

Cmatrix& Cmatrix::operator=(Cmatrix&& x) {
	delete[] array;
	nrows = x.nrows;
	ncols = x.ncols;
	array = x.array;
	x.array = nullptr;
	x.nrows = 0;
	x.ncols = 0;
#ifdef _MATRIXCOPYWARNING
	cout << "WARNING: Cmatrix move-assigned." << endl;
#endif
	return *this;
}

// ---- Conversions -----------------------------------------------------------------------------------------------
Cmatrix::Cmatrix(const Cmatrix& x, const Remap& remap1, const Remap& remap2) :
		Cmatrix(x.nrows, x.ncols) {
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++)
			array[j * nrows + i] = x(remap1.forward[i], remap2.forward[j]);
}

Cmatrix::Cmatrix(const Cmatrix& x, const Transpose& dummy) :
		Cmatrix(x.ncols, x.nrows) {
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++)
			array[j * nrows + i] = x.array[i * ncols + j];
}

Cmatrix::Cmatrix(const EigenMatrixXdAdaptor& X) :
		Cmatrix(X.rows(), X.cols()) {
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++)
			(*this)(i, j) = X(i, j);
}

template<>
Eigen::MatrixXd Cmatrix::convert() const {
	Eigen::MatrixXd M(nrows, ncols);
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++)
			M(i, j) = array[j * nrows + i];
	return M;
}

// ---- Specialized Constructors ----------------------------------------------------------------------------------
Cmatrix Cmatrix::Kronecker(const Cmatrix& seed, const int k) {
	assert(k >= 1);
	if (k == 1)
		return seed;
	Cmatrix Msub = Kronecker(seed, k - 1);
	return Kronecker(seed, Msub);
}

Cmatrix Cmatrix::Kronecker(const Cmatrix& X1, const Cmatrix& X2) {
	int _nrows = X2.nrows * X1.nrows;
	int _ncols = X2.ncols * X1.ncols;
	Cmatrix M(_nrows, _ncols);
	for (int J = 0; J < X1.ncols; J++)
		for (int I = 0; I < X1.nrows; I++) {
			FIELD s = X1(I, J);
			for (int j = 0; j < X2.ncols; j++)
				for (int i = 0; i < X2.nrows; i++)
					M(I * X2.nrows + i, J * X2.ncols + j) = s * X2(i, j);
		}
	return M;
}

Cmatrix Cmatrix::Random(const int _nrows, const int _ncols) {
	Cmatrix M(_nrows, _ncols);
	uniform_real_distribution<FIELD> distr;
	for (int i = 0; i < _nrows * _ncols; i++)
		M.array[i] = distr(randomNumberGenerator);
	return M;
}

Cmatrix Cmatrix::RandomSymmetric(const int n) {
	Cmatrix M(n, n);
	uniform_real_distribution<FIELD> distr;
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++) {
			FIELD t = distr(randomNumberGenerator);
			M.array[j * n + i] = t;
			M.array[i * n + j] = t;
		}
	return M;
}

Cmatrix Cmatrix::RandomPosDef(const int n) {
	Cmatrix A = Random(n, n);
	return A.dot(A);
}

// ---- Submatrices -----------------------------------------------------------------------------------------------
Cmatrix Cmatrix::submatrix(const int _nrows, const int _ncols,
		const int ioffset, const int joffset) const {
	Cmatrix R(_nrows, _ncols);
	assert(ioffset + _nrows <= nrows);
	assert(joffset + _ncols <= ncols);
	for (int j = 0; j < _ncols; j++)
		for (int i = 0; i < _nrows; i++)
			R.array[j * _nrows + i] =
					array[(j + joffset) * nrows + i + ioffset];
	return R;
}

// ---- Solvers ----------------------------------------------------------------------------------------------
Cmatrix Cmatrix::inverse() const {
	Eigen::MatrixXd A = convert<Eigen::MatrixXd>();
	return EigenMatrixXdAdaptor(A.inverse());
}

Cmatrix Cmatrix::inverseSqrt() const {
	Eigen::SelfAdjointEigenSolver < Eigen::MatrixXd
			> solver(convert<Eigen::MatrixXd>());
	return EigenMatrixXdAdaptor(solver.operatorInverseSqrt());
}

Cvector Cmatrix::solve(const Cvector& b) const {
	Eigen::MatrixXd A = convert<Eigen::MatrixXd>();
	Eigen::VectorXd be = b.convert<Eigen::VectorXd>();
	return EigenVectorXdAdaptor(A.colPivHouseholderQr().solve(be));
}

pair<Cmatrix*, Cvector*> Cmatrix::symmetricEigensolver() const {
	pair<Cmatrix*, Cvector*> p;
	Eigen::SelfAdjointEigenSolver < Eigen::MatrixXd
			> solver(convert<Eigen::MatrixXd>());
	p.first = new Cmatrix(EigenMatrixXdAdaptor(solver.eigenvectors()));
	p.second = new Cvector(EigenVectorXdAdaptor(solver.eigenvalues()));
	return p;
}

// ---- I/O -------------------------------------------------------------------------------------------------------
string Cmatrix::str(const Dense dummy) const {
	return Matrix::str(Dense());
}
string Cmatrix::str(const Sparse dummy) const {
	return Matrix::str(Sparse());
}
string Cmatrix::str() const {
	return Matrix::str(Dense());
}
ostream& operator<<(ostream& stream, const Cmatrix& x) {
	stream << x.str();
	return stream;
}

void Cmatrix::serialize(Rstream& rstream) const {
	rstream << "CsytleDenseMatrix{" << Rstream::endl;
	rstream.var("nrows", nrows);
	rstream.var("ncols", ncols);
	rstream.var("array", " *OMITTED*");
	rstream << "}" << Rstream::endl;
}

void Cmatrix::serialize(Bofstream& ofs) const {
	ofs.tag("Cmatrix", 0);
	ofs.write(nrows);
	ofs.write(ncols);
	ofs.write_array(array, nrows * ncols);
}

Cmatrix::Cmatrix(Bifstream& ifs) :
		DenseMatrix(0, 0) {
	ifs.check("Cmatrix", 0);
	ifs.read(nrows);
	ifs.read(ncols);
	ifs.read_array(array);
}

Cmatrix::Cmatrix(DenseMatrixFile& file) :
		Cmatrix(file.nrows, file.ncols) {
	for (int i = 0; i < nrows; i++)
		for (int j = 0; j < ncols; j++)
			file >> array[j * nrows + i];
}

Cmatrix::Cmatrix(SparseMatrixFile& file) :
		Cmatrix(file.nrows, file.ncols) {
	for (int i = 0; i < nrows * ncols; i++)
		array[i] = 0;
	for (auto it = file.begin(); it != file.end(); ++it)
		array[(*it).j * nrows + (*it).i] = (*it).value;
}
