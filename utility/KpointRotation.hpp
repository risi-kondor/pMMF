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



#ifndef _KpointRotation
#define _KpointRotation
#include "ElementaryRotation.hpp"

class KpointRotation: public ElementaryRotation {
public:

	KpointRotation(const KpointRotation& x) :
			k(x.k) {
#ifdef _UTILITYCOPYWARNING
		cout<<"WARNING: KpointRotation copied."<<endl;
#endif
		ix = new int[k];
		for (int i = 0; i < k; i++)
			ix[i] = x.ix[i];
		q = new double[k * k];
		for (int i = 0; i < k * k; i++)
			q[i] = x.q[i];
		temp = new double[k];
		tempp = new FIELD*[k];
	}

	KpointRotation(KpointRotation&& x) :
			k(x.k) {
		ix = x.ix;
		x.ix = nullptr;
		q = x.q;
		x.q = nullptr;
		temp = x.temp;
		x.temp = nullptr;
		tempp = x.tempp;
		x.tempp = nullptr;
		x.k = 0;
	}

public:

	KpointRotation(const int _k) :
			k(_k) {
		ix = new int[k];
		q = new double[k * k];
		temp = new double[k];
		tempp = new FIELD*[k];
	}
	KpointRotation(const int _k, const Identity dummy);
	KpointRotation(const int _k, const Random dummy, const int n);

	KpointRotation(const IndexSet& I, const Cmatrix& M);

	~KpointRotation() {
		delete[] ix;
		delete[] q;
		delete[] temp;
		delete[] tempp;
	}

public:
	string str() const;
	void serialize(Rstream& s) const;
	void serialize(Bofstream& ofs) const;
	void save(const char* filename) {
		Bofstream ofs(filename);
		serialize(ofs);
	}
	KpointRotation(Bifstream& ifs);

public:

	int k;
	int* ix;
	double* q;
	double* temp;
	FIELD** tempp;
};

#endif
