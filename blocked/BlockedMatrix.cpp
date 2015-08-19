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



#include "BlockedMatrix.hpp"

template<>
BlockedMatrix<Cmatrix>::BlockedMatrix(const Cmatrix& M,
		const Bstructure& rstruct, const Bstructure& cstruct) :
		BlockedMatrix(rstruct.size(), cstruct.size()) {
	int joffset = 0;
	for (int J = 0; J < ntowers; J++) {
		int ioffset = 0;
		for (int I = 0; I < nstreets; I++) {
			ptr(I, J) = new Cmatrix(
					M.submatrix(rstruct[I], cstruct[J], ioffset, joffset));
			ioffset += rstruct[I];
		}
		joffset += cstruct[J];
	}
}

template<>
template<>
Cmatrix BlockedMatrix<Cmatrix>::convert<Cmatrix>() const {
	int _nrows = nrows();
	int _ncols = ncols();
	Cmatrix M(_nrows, _ncols);
	Bstructure rstruct = rowstructure();
	Bstructure cstruct = colstructure();

	int joffset = 0;
	for (int J = 0; J < ntowers; J++) {
		int ioffset = 0;
		for (int I = 0; I < nstreets; I++) {
			const Cmatrix& B = (*this)(I, J);
			for (int j = 0; j < B.ncols; j++)
				for (int i = 0; i < B.nrows; i++)
					M(ioffset + i, joffset + j) = B(i, j);
			ioffset += rstruct[I];
		}
		joffset += cstruct[J];
	}

	return M;
}

// ---- I/O ------------------------------------------------------------------------------------------------------

template<>
BlockedMatrix<Cmatrix>::BlockedMatrix(DenseMatrixFile& file,
		const BlockStructure& rstruct, const BlockStructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	assert(file.nrows == nrows());
	assert(file.ncols == ncols());
	for (int I = 0; I < nstreets; I++)
		for (int i = 0; i < rstruct[I]; i++)
			for (int J = 0; J < ntowers; J++)
				for (int j = 0; j < cstruct[J]; j++)
					file >> (*this)(I, i, J, j);
}

template<>
BlockedMatrix<Cmatrix>::BlockedMatrix(SparseMatrixFile& file,
		const BlockStructure& rstruct, const BlockStructure& cstruct) :
		BlockedMatrix(rstruct, cstruct, Zero()) {
	assert(file.nrows == nrows());
	assert(file.ncols == ncols());
	for (auto it = file.begin(); it != file.end(); ++it) {
		int i = (*it).i;
		int I = findI(i, rstruct);
		int j = (*it).j;
		int J = findI(j, cstruct);
		(*this)(I, i, J, j) = (*it).value;
	}
}

template<>
BlockedMatrix<MatrixXv>::BlockedMatrix(SparseMatrixFile& file,
		const BlockStructure& rstruct, const BlockStructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	assert(file.nrows == nrows());
	assert(file.ncols == ncols());
	for (auto it = file.begin(); it != file.end(); ++it) {
		int i = (*it).i;
		int I = findI(i, rstruct);
		int j = (*it).j;
		int J = findI(j, cstruct);
		ptr(I, J)->column[j]->insert(i, (*it).value);
	}
	MultiLoop(nstreets, ntowers,
			[this](const int i, const int j) {ptr(i,j)->tidy();});
}

template<>
BlockedMatrix<MatrixXl>::BlockedMatrix(SparseMatrixFile& file,
		const BlockStructure& rstruct, const BlockStructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	assert(file.nrows == nrows());
	assert(file.ncols == ncols());
	for (auto it = file.begin(); it != file.end(); ++it) {
		int i = (*it).i;
		int I = findI(i, rstruct);
		int j = (*it).j;
		int J = findI(j, cstruct);
		ptr(I, J)->column[j]->insert(i, (*it).value);
	}
	MultiLoop(nstreets, ntowers,
			[this](const int i, const int j) {ptr(i,j)->tidy();});
}

template<>
BlockedMatrix<MatrixXh>::BlockedMatrix(SparseMatrixFile& file,
		const BlockStructure& rstruct, const BlockStructure& cstruct) :
		BlockedMatrix(rstruct, cstruct) {
	assert(file.nrows == nrows());
	assert(file.ncols == ncols());
	for (auto it = file.begin(); it != file.end(); ++it) {
		int i = (*it).i;
		int I = findI(i, rstruct);
		int j = (*it).j;
		int J = findI(j, cstruct);
		ptr(I, J)->column[j]->insert(i, (*it).value);
	}
	MultiLoop(nstreets, ntowers,
			[this](const int i, const int j) {ptr(i,j)->tidy();});
}
