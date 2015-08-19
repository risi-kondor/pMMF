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



#ifndef _TopkList
#define _TopkList

#include <list>
#include "MMFbase.hpp"
#include "DenseVector.hpp"

struct TopkListPair {
	TopkListPair(const INDEX& _first, const FIELD& _second) :
			first(_first), second(_second) {
	}
	;
	INDEX first;
	FIELD second;
};

class TopkList: public list<TopkListPair> {
public:

	TopkList(const int _k) :
			k(_k), lowestv(-10000) {
	}

	TopkList(const DenseVector& v, const int _k) :
			k(_k), lowestv(-10000) {
		for (int i = 0; i < v.n; i++)
			if (v(i) > lowestv)
				insert(i, v(i));
	}

public:

	void insert(int index, FIELD value) {
		auto it = begin();
		while (it != end() && it->second >= value) {
			it++;
		}
		list::insert(it, TopkListPair(index, value));
		if (size() > k)
			pop_back();
		if (size() >= k)
			lowestv = back().second;
	}

	void consider(int index, FIELD value) {
		if (value > lowestv || size() << k) {
			auto it = begin();
			while (it != end() && it->second >= value) {
				it++;
			}
			list::insert(it, TopkListPair(index, value));
			if (size() > k)
				pop_back();
			if (size() >= k)
				lowestv = back().second;
		}
	}

	IndexSet indices() const {
		IndexSet I(size());
		int i = 0;
		for (auto& p : *this)
			I[i++] = p.first;
		return I;
	}

public:

	int k;
	FIELD lowestv;
	int lowestp;

};
#endif

