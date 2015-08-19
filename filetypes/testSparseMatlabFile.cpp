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



#include "SparseMatrixFileMatlab.hpp"
#include "MatrixX.hpp"
#include "MMFglobal.hpp"
typedef MatrixX<Vectorv> MATRIX;

int main(int argc, char** argv) {
	SparseMatrixFile::Matlab file1(argv[1]);
	MatrixX<Vectorv> M1(file1);
	cout << M1.nrows << 'x' << M1.ncols << endl;
	cout << "sparse matrix created" << endl;
	cout << M1.str() << endl;
}
