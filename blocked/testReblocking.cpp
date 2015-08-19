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

	if (argc < 9) {
		cout
				<< "./testReblocking <nsubrowsMax> <nsubcolsMax> <nstreetsMax> <ntowersMax> <niter> <blockType> <increment1> <increment2> <verbosity>"
				<< endl;
	}

	int nsubrowsMax = 3;
	if (argc >= 2)
		sscanf(argv[1], "%d", &nsubrowsMax);
	int nsubcolsMax = 3;
	if (argc >= 3)
		sscanf(argv[2], "%d", &nsubcolsMax);
	int nstreetsMax = 2;
	if (argc >= 4)
		sscanf(argv[3], "%d", &nstreetsMax);
	int ntowersMax = 2;
	if (argc >= 5)
		sscanf(argv[4], "%d", &ntowersMax);
	int niter = 10;
	if (argc >= 6)
		sscanf(argv[5], "%d", &niter);
	int blockType = 4;
	if (argc >= 7)
		sscanf(argv[6], "%d", &blockType);
	int increment1 = 100;
	if (argc >= 8)
		sscanf(argv[7], "%d", &increment1);
	int increment2 = 10;
	if (argc >= 9)
		sscanf(argv[8], "%d", &increment2);
	int verbosity = 0;
	if (argc >= 10)
		sscanf(argv[9], "%d", &verbosity);

	int totalError = 0;
	for (int nsubrows = 1; nsubrows <= nsubrowsMax; nsubrows += increment1) {
		for (int nsubcols = 1; nsubcols <= nsubcolsMax; nsubcols +=
				increment1) {
			for (int nstreets = 1; nstreets <= nstreetsMax; nstreets +=
					increment2) {
				for (int ntowers = 1; ntowers <= ntowersMax; ntowers +=
						increment2) {

					BlockStructure sx;
					for (int i = 0; i < nstreets; i++)
						sx.push_back(nsubrows);
					BlockStructure sy;
					for (int i = 0; i < ntowers; i++)
						sy.push_back(nsubcols);

					if (verbosity >= 1)
						cout << nsubrows << " " << nsubcols << " " << nstreets
								<< " " << ntowers << " ";
					int accumulatedError = 0;
					for (int i = 0; i < niter; i++) {
						BlockedMatrix<BLOCK> B = BlockedMatrix<BLOCK>::Random(
								sx, sy);
						if (verbosity >= 2)
							cout << "B=" << endl << B.str(Dense()) << endl;

						BlockedRemap r1x(sx, sx, Random());
						if (verbosity)
							cout << r1x.str() << endl;
						BlockedRemap r1y(sy, sy, Random());
						if (verbosity)
							cout << r1y.str() << endl;

						BlockedMatrix<BLOCK> B1hv(B, r1x, r1y);
						if (verbosity >= 2)
							cout << B1hv.str() << endl;

						BlockedMatrix<BLOCK> B1hvi(B1hv, r1x, r1y, true);
						if (verbosity >= 2)
							cout << B1hvi.str() << endl;

						accumulatedError += B.diff2(B1hvi);
					}
					totalError += accumulatedError;
					if (verbosity >= 1)
						cout << accumulatedError << endl;
				}
			}
		}
	}
	cout << "Total error=" << totalError << endl;
}
