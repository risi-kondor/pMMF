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
#include "MMFglobal.hpp"

int main(int argc, char** argv) {
	int k = 3;
	if (argc >= 2)
		sscanf(argv[1], "%d", &k);
	Cmatrix seed(2, 2);
	seed(0, 0) = 1;
	seed(1, 1) = 1;
	seed(0, 1) = 0.6;
	seed(1, 0) = 0.6;

	Cmatrix M = Cmatrix::Kronecker(seed, 3);
	cout << M.str() << endl;
	cout << Cmatrix::Identity(5) << endl;
}
