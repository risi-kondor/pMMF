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
#include "MMF.hpp"
#include "SparseMatrixFileASCII.hpp"
#include "DenseMatrixFileASCII.hpp"
#include "MMFglobal.hpp"

typedef Cmatrix BLOCK;

int main(int argc, char** argv) {

	int n = 25;
	if (argc >= 2)
		sscanf(argv[1], "%d", &n);
	int nclusters = 4;
	if (argc >= 3)
		sscanf(argv[2], "%d", &nclusters);
	mlog.verbosity = 4;
	if (argc >= 4)
		sscanf(argv[3], "%d", &mlog.verbosity);
	int criterion = 0;
	if (argc >= 5)
		sscanf(argv[4], "%d", &criterion);

	MMFmatrix<BLOCK> A(n, Random());
	cout << endl;

	MMF mmf(A,
			MMFparams::k(2).fraction(0.5).ndensestages(0).nsparsestages(3).nclusters(
					nclusters).prenormalize(true).ncoreclusters(1).dcore(10).selection_normalized(
					false).n_eliminate_per_rotation(1).fraction_eliminate_after(
					0).selection_criterion(criterion));

	FIELD spectral2 = mmf.diffSpectral(A);
	cout << "Spectral norm error= " << sqrt(spectral2) << endl;
	cout << mmf.core.str() << endl;
}
