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



#include "MMF.hpp"
#include "MMFglobal.hpp"

typedef Cmatrix BLOCK;

int main(int argc, char** argv) {

	char* filename = argv[1];
	int n = 12;
	if (argc >= 3)
		sscanf(argv[2], "%d", &n);
	int nclusters = 2;
	if (argc >= 4)
		sscanf(argv[3], "%d", &nclusters);
	mlog.verbosity = 4;
	if (argc >= 5)
		sscanf(argv[4], "%d", &mlog.verbosity);

	MMFmatrix<BLOCK> A = MMFmatrix<BLOCK>::RandomPosDef(n);
	cout << A.str() << endl;
	MMF mmf(A,
			MMFparams::k(2).nsparsestages(1).ndensestages(1).nclusters(
					nclusters).prenormalize(true));
	mmf.save(filename);

	Bifstream ifs(filename);
	MMF mmf2(ifs);

	BlockedMatrix<Cmatrix> R = mmf.reconstruction();
	cout << R.str() << endl;

	BlockedMatrix<Cmatrix> R2 = mmf2.reconstruction();
	cout << R2.str() << endl;

	cout << "Difference = " << R2.diff2(R) << endl << endl;

}
