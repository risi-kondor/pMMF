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



#ifndef _Clusters 
#define _Clusters

#include "MMFbase.hpp"

template<class TYPE>
class Clusters: public vector<vector<TYPE>*> {
public:

	Clusters(const Clusters<TYPE>& x) {
#ifdef _UTILCOPYWARNING
		cout << "WARNING: Clusters copied." << endl;
#endif
		for (auto p : x)
			push_back(new vector<TYPE>(*p));
	}

	Clusters(Clusters<TYPE> && x) {
#ifdef _UTILCOPYWARNING
		cout << "WARNING: Clusters moved." << endl;
#endif
		for (int i = 0; i < x.size(); i++) {
			vector<vector<TYPE>*>::push_back(x[i]);
			x[i] = nullptr;
		}
	}

	Clusters<TYPE>& operator=(const Clusters<TYPE>& x) {
#ifdef _UTILCOPYWARNING
		cout << "WARNING: Clusters assigned." << endl;
#endif
		for (auto p : *this)
			delete p;
		vector<vector<TYPE>*>::clear();
		for (auto p : x)
			push_back(new vector<TYPE>(*p));
		return *this;
	}

	~Clusters() {
		for (auto p : *this)
			delete p;
	}

public:

	Clusters() {
	}

	Clusters(const int n) {
		for (int i = 0; i < n; i++)
			vector<vector<TYPE>*>::push_back(new vector<TYPE>());
	}

	Clusters<TYPE>& add(const Clusters<TYPE>& x) {
		assert(x.size() == vector<vector<TYPE>*>::size());
		for (int i = 0; i < vector<vector<TYPE>*>::size(); i++)
			for (int j = 0; j < x[i]->size(); j++)
				vector<vector<TYPE>*>::operator[](i)->push_back((*x[i])[j]);
		return *this;
	}

public:

	string str() const;

};

#endif
