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



#include "DenseMatrixFileASCII.hpp"
#include "Cmatrix.hpp"
#include "MatrixX.hpp"
#include "MMFglobal.hpp"

typedef MatrixX<Vectorv> MATRIX;

int main(int argc, char** argv) {

	char* filename = new char[256];
	sprintf(filename, "mx1.dat");
	if (argc >= 2)
		filename = argv[1];

	DenseMatrixFile::ASCII file1(filename);
	Cmatrix M(file1);
	if (M.nrows <= 12 && M.ncols <= 12) {
		cout << M.str() << endl;
	} else {
		Cmatrix M2(M, 8, 8);
		cout << M2.str() << endl;
	}
}
