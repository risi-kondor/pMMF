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



#ifndef _Wtransform
#define _Wtransform

#include "BlockedMatrix.hpp"

template<class VECTOR>
class Wtransform {
public:

	Wtransform(const int nlevels) :
			v(0), w(nlevels, nullptr) {
	}

	Wtransform(const Wtransform& X) :
			v(X.v) {
		cout << "Warning: Wavelet transform copied" << endl;
		for (int i = 0; i < X.w.size(); i++)
			w.push_back(new VECTOR(*X.w[i]));
	}

	Wtransform& operator=(const Wtransform& X) {
		cout << "Warning: Wavelet transform assigned" << endl;
		v = X.v;
		for (int i = 0; i < w.size(); i++)
			delete w[i];
		for (int i = 0; i < X.w.size(); i++)
			w.push_back(new VECTOR(X.w[i]));
	}

	~Wtransform() {
		for (auto* p : w)
			delete p;
	}

public:
	BlockedVector<VECTOR> v;
	vector<VECTOR*> w;

};

#endif
