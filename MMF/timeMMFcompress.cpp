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

extern Log mlog;
typedef MatrixXl BLOCK;

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
	mlog.verbosity = 0;
	if (argc >= 7)
		sscanf(argv[6], "%d", &mlog.verbosity);

	SparseMatrixFile::ASCII file(filename);
	BLOCK M(file);

	std::vector<float> timingVec(100);
	std::vector<int> nnzVec(100);
	int count = 0;

	for (int i = 2; i < 90; i = i + 2) {
		BLOCK M2(M, i * 100, i * 100);
		MMFmatrix<BLOCK> A(M2, true, true);
		nnzVec[count] = A.nnz();
		cout << "ncols, nnz " << A.ncols() << " " << A.nnz() << endl;

		float timing = 0;
		for (int t = 0; t < 5; t++) {
			mlog.startClock();
			MMF mmf(A,
					MMFparams::k(2).nsparsestages(15).ndensestages(ndensestages).fraction(
							frac).nclusters(nclusters).minclustersize(10).maxclustersize(
							2000).maxclusterdepth(5).bypass(true));
			timing += mlog.clock();
		}
		timingVec[count] = timing / 5.0;
		cout << timingVec[count] << " " << timing / 5.0 << endl;
		count++;

	}
	cout << "count " << count << endl;
	for (int i = 0; i < count; i++) {
		cout << nnzVec[i] << " " << timingVec[i] << endl;
	}
}
