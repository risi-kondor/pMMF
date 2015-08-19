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



#ifndef _Remap
#define _Remap
#include "MMFbase.hpp"

class Remap {
public:
	Remap() {
		forward = NULL;
		backward = NULL;
		n = 0;
	}

	Remap(const int _n) :
			n(_n) {
		forward = new int[n];
		for (int i = 0; i < n; i++)
			forward[i] = i;
		backward = new int[n];
		for (int i = 0; i < n; i++)
			backward[i] = i;
	}

	Remap(const int _n, const Random random);

	~Remap() {
		delete[] forward;
		delete[] backward;
	}

public:
	int operator()(const int i) const {
		return forward[i];
	}

	Remap& swap(const int i, const int j) {
		int t = forward[i];
		forward[i] = forward[j];
		forward[j] = t;
		backward[forward[i]] = i;
		backward[forward[j]] = j;
		return *this;
	}

	void fixBackwardMap() {
		for (int i = 0; i < n; i++)
			backward[forward[i]] = i;
	}

	string str() const;

public:
	int n;
	int* forward;
	int* backward;

};

#endif 
