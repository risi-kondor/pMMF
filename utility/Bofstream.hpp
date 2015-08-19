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



#ifndef _Bofstream
#define _Bofstream

#include "MMFbase.hpp"
#include <fstream>
#include<string.h>

class Bofstream: public ofstream {
public:

	Bofstream(const char* filename) {
		open(filename, ios::binary);
	}

	~Bofstream() {
		close();
	}

public:

	Bofstream& write(const int& x) {
		ofstream::write(reinterpret_cast<const char*>(&x), sizeof(int));
		return *this;
	}
	Bofstream& write(const double& x) {
		ofstream::write(reinterpret_cast<const char*>(&x), sizeof(double));
		return *this;
	}
	Bofstream& write(const bool& x) {
		ofstream::write(reinterpret_cast<const char*>(&x), sizeof(bool));
		return *this;
	}
	Bofstream& write(const BlockIndexPair& x) {
		ofstream::write(reinterpret_cast<const char*>(&x.block), sizeof(INDEX));
		ofstream::write(reinterpret_cast<const char*>(&x.index), sizeof(INDEX));
		return *this;
	}

	template<class TYPE>
	Bofstream& write_vector(const vector<TYPE>& x) {
		int n = x.size();
		ofstream::write(reinterpret_cast<const char*>(&n), sizeof(int));
		for (auto p : x)
			write(p);
		return *this;
	}

	template<class TYPE>
	Bofstream& write_array(const TYPE* x, const int n) {
		ofstream::write(reinterpret_cast<const char*>(&n), sizeof(int));
		ofstream::write(reinterpret_cast<const char*>(x), n * sizeof(TYPE));
		//for(int i=0; i<n; i++) write(x[i]);
		return *this;
	}

public:

	template<class TYPE>
	Bofstream& seriate(const TYPE* x) {
		x->serialize(*this);
		return *this;
	}

	template<class TYPE>
	Bofstream& serialize(const TYPE& x) {
		x.serialize(*this);
		return *this;
	}

	template<class TYPE>
	Bofstream& seriate_vector(const vector<TYPE*>& x) {
		int n = x.size();
		ofstream::write(reinterpret_cast<const char*>(&n), sizeof(int));
		for (auto p : x)
			p->serialize(*this);
		return *this;
	}

	template<class TYPE>
	Bofstream& seriate_array(TYPE** x, const int n) { // why did const have to be removed?
		ofstream::write(reinterpret_cast<const char*>(&n), sizeof(int));
		for (int i = 0; i < n; i++)
			x[i]->serialize(*this);
		return *this;
	}

public:

	void tag(const char* name, const int _version) {
		int version = _version;
		ofstream::write(name, strlen(name));
		ofstream::write(reinterpret_cast<const char*>(&version), sizeof(int));
	}

};

#endif
