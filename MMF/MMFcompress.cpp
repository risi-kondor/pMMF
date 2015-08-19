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

typedef MatrixXv BLOCK;

int main(int argc, char** argv) {

	char* filename = argv[1];
	int nsparsestages = 1;
	if (argc >= 3)
		sscanf(argv[2], "%d", &nsparsestages);
	int ndensestages = 1;
	if (argc >= 4)
		sscanf(argv[3], "%d", &ndensestages);
	int nclusters = 2;
	if (argc >= 5)
		sscanf(argv[4], "%d", &nclusters);
	float frac = 0.5;
	if (argc >= 6)
		sscanf(argv[5], "%f", &frac);
	mlog.verbosity = 3;
	if (argc >= 7)
		sscanf(argv[6], "%d", &mlog.verbosity);

	SparseMatrixFile::ASCII file(filename);
	MMFmatrix<BLOCK> A(file);

	MMF mmf(A,
			MMFparams::k(2).nsparsestages(nsparsestages).ndensestages(
					ndensestages).fraction(frac).nclusters(nclusters).minclustersize(
					10).maxclustersize(500).maxclusterdepth(3).bypass(true));
}
