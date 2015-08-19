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



#ifndef _Bifstream
#define _Bifstream

#include "MMFbase.hpp"
#include <fstream>
#include<string.h>

extern char strbuffer[255];

class Bifstream: public ifstream {
public:

	Bifstream(const char* filename) {
		open(filename, ios::binary);
	}

	~Bifstream() {
		close();
	}

public:

	Bifstream& read(int& x) {
		ifstream::read(reinterpret_cast<char*>(&x), sizeof(int));
		return *this;
	}
	Bifstream& read(double& x) {
		ifstream::read(reinterpret_cast<char*>(&x), sizeof(double));
		return *this;
	}
	Bifstream& read(bool& x) {
		ifstream::read(reinterpret_cast<char*>(&x), sizeof(bool));
		return *this;
	}
	Bifstream& read(BlockIndexPair& x) {
		ifstream::read(reinterpret_cast<char*>(&x.block), sizeof(INDEX));
		ifstream::read(reinterpret_cast<char*>(&x.index), sizeof(INDEX));
		return *this;
	}

	template<class TYPE>
	Bifstream& read_vector(vector<TYPE>& x) {
		int n;
		ifstream::read(reinterpret_cast<char*>(&n), sizeof(int));
		x.resize(n);
		for (int i = 0; i < n; i++)
			read(x[i]);
		return *this;
	}

	template<class TYPE>
	Bifstream& read_array(TYPE*& x) {
		int n;
		ifstream::read(reinterpret_cast<char*>(&n), sizeof(int));
		x = new TYPE[n];
		ifstream::read(reinterpret_cast<char*>(x), n * sizeof(TYPE));
		return *this;
	}

public:

	template<class TYPE>
	Bifstream& unseriate(TYPE*& x) {
		x = new TYPE(*this);
		return *this;
	}

	template<class TYPE>
	Bifstream& unserialize(TYPE& x) {
		x = TYPE(*this);
		return *this;
	}

	template<class TYPE>
	Bifstream& unseriate_vector(vector<TYPE*>& x) {
		int n;
		ifstream::read(reinterpret_cast<char*>(&n), sizeof(int));
		x.resize(n);
		for (int i = 0; i < n; i++)
			x[i] = new TYPE(*this);
		return *this;
	}

	template<class TYPE>
	Bifstream& unseriate_array(TYPE**& x) {
		int n;
		ifstream::read(reinterpret_cast<char*>(&n), sizeof(int));
		x = new TYPE*[n];
		for (int i = 0; i < n; i++)
			x[i] = new TYPE(*this);
		return *this;
	}

public:

	bool check(char const* name, const int version) {

		ifstream::get(strbuffer, 255, '\0');
		if (strcmp(strbuffer, name) != 0) {
			cout << "Error in Bifstream: " << name << " expected, " << strbuffer
					<< " found." << endl;
			return true;
		}

		int readver;
		ifstream::read(reinterpret_cast<char*>(&readver), sizeof(int));
		if (readver != version) {
			cout << "Error in Bifstream: version " << version << " expected, "
					<< readver << " found." << endl;
			return true;
		}

		return false;
	}

public:

};

#endif
