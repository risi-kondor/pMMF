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



#include "Log.hpp"
#include "MMFmatrix.hpp"

#include "MatrixX.hpp"
#include "Cmatrix.hpp"
#include "MMFglobal.hpp"

typedef MatrixXh BLOCK;

int main(int argc, char** argv) {

	int n = 10;
	if (argc >= 2)
		sscanf(argv[1], "%d", &n);

	MMFmatrix<BLOCK> A(n, Random(), false);
	cout << A << endl;
	cout << *A.process[0]->gram << endl;

	// for the sparse matrix types, by the end, the Gram matrix is supposed to be all zero.
	for (int i = 0; i < n; i++) {
		A.process[0]->deactivate(i);
		cout << *A.process[0]->gram << endl;
	}

}
