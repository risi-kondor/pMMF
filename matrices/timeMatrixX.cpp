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
#include "MatrixX.hpp"
#include "GivensRotation.hpp"
#include "MMFglobal.hpp"

template<class TYPE>
void timeTranspose(TYPE& M, const int niter = 1000) {
	mlog.startClock();
	cout << niter << " transpositions of " << typeid(M).name() << "(" << M.nrows
			<< "," << M.ncols << "): " << mlog.clock() << " s" << endl;
}

template<class TYPE>
void timeLeftGivensRotation(TYPE& M, const vector<GivensRotation>& Q) {
	mlog.startClock();
	for (int i = 0; i < Q.size(); i++)
		M.applyFromLeft(Q[i]);
	cout << Q.size() << " left Givens rotations of " << typeid(M).name() << "("
			<< M.nrows << "," << M.ncols << "): " << mlog.clock() << " s"
			<< endl;
}

template<class TYPE>
void timeRightGivensRotation(TYPE& M, const vector<GivensRotation>& Q) {
	mlog.startClock();
	for (int i = 0; i < Q.size(); i++)
		M.applyFromRightT(Q[i]);
	cout << Q.size() << " right Givens rotations of " << typeid(M).name() << "("
			<< M.nrows << "," << M.ncols << "): " << mlog.clock() << " s"
			<< endl;
}

template<class TYPE>
void timeRandomAccess(const TYPE& M, const vector<int>& i1,
		const vector<int>& i2) {
	mlog.startClock();
	for (int i = 0; i < i1.size(); i++)
		FIELD t = M(i1[i], i2[i]);
	cout << i1.size() << " random accesses into " << typeid(M).name() << "("
			<< M.nrows << "," << M.ncols << "): " << mlog.clock() << " s"
			<< endl;
}

int main(int argc, char** argv) {
	int nrows = 6;
	if (argc >= 2)
		sscanf(argv[1], "%d", &nrows);
	int ncols = 6;
	if (argc >= 3)
		sscanf(argv[2], "%d", &ncols);
	float p = 0.5;
	if (argc >= 4)
		sscanf(argv[3], "%f", &p);
	Cmatrix Md = Cmatrix::Random(nrows, ncols);

	MatrixXv Mv(nrows, ncols, Random(p));
	MatrixXl Ml(nrows, ncols, Random(p));
	MatrixXh Mh(nrows, ncols, Random(p));

	if (true) {
		cout << Md.str(Dense()) << endl;
		cout << Mv.str(Dense()) << endl;
		cout << Ml.str(Dense()) << endl;
		cout << Mh.str(Dense()) << endl;
	}
	if (true) {
		timeTranspose(Md);
		timeTranspose(Mv);
		timeTranspose(Ml);
		timeTranspose(Mh);
		cout << endl;
	}
	if (true) {
		vector<GivensRotation> Q;
		for (int i = 0; i < 1000; i++)
			Q.push_back(GivensRotation(Random(), nrows));

		timeLeftGivensRotation(Md, Q);
		timeLeftGivensRotation(Mv, Q);
		timeLeftGivensRotation(Ml, Q);
		timeLeftGivensRotation(Mh, Q);
		cout << endl;
	}
	if (true) {
		vector<GivensRotation> Q;
		for (int i = 0; i < 1000; i++)
			Q.push_back(GivensRotation(Random(), nrows));

		timeRightGivensRotation(Md, Q);
		timeRightGivensRotation(Mv, Q);
		timeRightGivensRotation(Ml, Q);
		timeRightGivensRotation(Mh, Q);
		cout << endl;
	}
	if (true) {
		int ntrials = 100000;
		uniform_int_distribution<int> distri(0, ncols - 1);
		vector<int> i1(ntrials);
		for (int i = 0; i < ntrials; i++)
			i1[i] = distri(randomNumberGenerator);
		vector<int> i2(ntrials);
		for (int i = 0; i < ntrials; i++)
			i2[i] = distri(randomNumberGenerator);

		timeRandomAccess(Md, i1, i2);
		timeRandomAccess(Mv, i1, i2);
		timeRandomAccess(Ml, i1, i2);
		timeRandomAccess(Mh, i1, i2);
		cout << endl;
	}
}
