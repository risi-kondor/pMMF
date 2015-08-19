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
#include "MatrixX.hpp"
#include "Cmatrix.hpp"
#include "MMFglobal.hpp"
typedef MatrixXv BLOCK;

int main(int argc, char** argv) {

	int nsubrows = 3;
	if (argc >= 2)
		sscanf(argv[1], "%d", &nsubrows);
	int nsubcols = 3;
	if (argc >= 3)
		sscanf(argv[2], "%d", &nsubcols);
	int nstreets = 2;
	if (argc >= 4)
		sscanf(argv[3], "%d", &nstreets);
	int ntowers = 2;
	if (argc >= 5)
		sscanf(argv[4], "%d", &ntowers);

	BlockStructure sx;
	for (int i = 0; i < nstreets; i++)
		sx.push_back(nsubrows);
	BlockStructure sy;
	for (int i = 0; i < ntowers; i++)
		sy.push_back(nsubcols);

	MatrixXl M(nstreets * nsubrows, ntowers * nsubcols, Random());
	cout << M.str(Dense()) << endl;

	BlockedMatrix<BLOCK> B(M, sx, sy);
	cout << B.str(Dense()) << endl;

	cout << "Row map:" << endl;
	BlockedRemap r1x(sx, sx, Random());
	cout << r1x << endl;
	cout << "Column map:" << endl;
	BlockedRemap r1y(sy, sy, Random());
	cout << r1y << endl;

	cout << "Rows remapped:" << endl;
	BlockedMatrix<BLOCK> B1v(B, r1x, Identity());
	cout << B1v << endl;

	cout << "Columns remapped:" << endl;
	BlockedMatrix<BLOCK> B1h(B, Identity(), r1y);
	cout << B1h << endl;

	cout << "Rows and columns remapped:" << endl;
	BlockedMatrix<BLOCK> B1hv(B, r1x, r1y);
	cout << B1hv << endl;

	cout << "Inverse map applied to rows and columns:" << endl;
	BlockedMatrix<BLOCK> B1hvi(B1hv, r1x, r1y, true);
	cout << B1hvi << endl;
}
