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



#ifndef _Vectorv
#define _Vectorv

#include "SparseVector.hpp"
#include "GivensRotation.hpp"
#include "KpointRotation.hpp"
#include "Cvector.hpp"
#include <algorithm>

class Vectorv: public SparseVector, public vector<SVpair> {
public:
public:
	// constructors
	Vectorv(const int _n) :
			SparseVector(_n), sorted(1) {
	}
	Vectorv(const int _n, const Random& dummy); // obsolete

public:
	// named constructors
	static Vectorv Random(const int _n, const FIELD p = 0.5);
	static Vectorv Zero(const int _n) {
		return Vectorv(_n);
	}

public:
	// converters
	Vectorv(const Cvector& x);
	Vectorv(const Vectorl& x);
	Vectorv(const Vectorh& x);
	Vectorv(const Vectorv& x, const Remap& remap);

public:
	// element access
	FIELD& operator()(const int i) {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first != i) {
			it++;
		}
		if (it == end()) {
			push_back(SVpair(i, 0));
			sorted = 0;
			return back().second;
		}
		return it->second;
	}
	FIELD operator()(const int i) const {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first != i) {
			it++;
		}
		if (it == end())
			return 0;
		else
			return it->second;
	}
	FIELD read(const int i) const {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first != i) {
			it++;
		}
		if (it == end())
			return 0;
		else
			return it->second;
	}
	void foreach(std::function<void(INDEX, FIELD&)> lambda) {
		for (auto& p : *this)
			lambda(p.first, p.second);
	}
	bool isFilled(const int i) const {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first != i) {
			it++;
		}
		return (it != end());
	}

	int nFilled() const {
		return size();
	}

public:
	// sparse vector methods
	FIELD* findptr(const int i) {
		auto it = begin();
		if (sorted) {
			while (it != end() && it->first < i) {
				it++;
			}
			if (it != end() && it->first == i)
				return &it->second;
		} else
			while (it != end()) {
				if (it->first == i)
					return &it->second;
				it++;
			}
		return &dummyZero;
	}
	void insert(const int i, const FIELD value) {
		assert(i < n);
		push_back(SVpair(i, value));
		sorted = 0;
	}
	void append(const int i, const FIELD value) {
		assert(i < n);
		push_back(SVpair(i, value));
	}
	void zero(const int i) {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first != i) {
			it++;
		}
		if (it != end())
			it->second = 0;
	}
	void sort() {
		if (!sorted)
			std::sort(begin(), end(), SVpairComparator);
		sorted = 1;
	}
	void tidy() {
		sort();
		auto v(*this);
		clear();
		int i = -1;
		for (auto& p : v)
			if (p.first != i && p.second != 0) {
				push_back(p);
				i = p.first;
			}
	}

public:
	// scalar-valued operations
	int nnz() const {
		const_cast<Vectorv&>(*this).tidy();
		return size();
	}

	int argmax() const {
		if (size() == 0)
			return 0;
		int best = begin()->first;
		FIELD max = begin()->second;
		for (auto p : *this)
			if (p.second > max) {
				best = p.first;
				max = p.second;
			}
		return best;
	}
	int argmax_abs() const {
		if (size() == 0)
			return 0;
		int best = begin()->first;
		FIELD max = fabs(begin()->second);
		for (auto p : *this)
			if (fabs(p.second) > max) {
				best = p.first;
				max = fabs(p.second);
			}
		return best;
	}
	FIELD norm2() const {
		FIELD t = 0;
		for (auto& p : *this)
			t += p.second * p.second;
		return t;
	}
	FIELD diff2(const Vectorv& x) const {
		assert(x.n == n);
		FIELD t = 0;
		const_cast<Vectorv&>(*this).sort();
		const_cast<Vectorv&>(x).sort();
		auto it = begin();
		for (auto& xv : x) {
			while (it != end() && it->first < xv.first) {
				t += (it->second) * (it->second);
				it++;
			}
			if (it != end() && it->first == xv.first) {
				t += (it->second - xv.second) * (it->second - xv.second);
				it++;
			} else {
				t += xv.second * xv.second;
			}
		}
		return t;
	}
	FIELD dot(const Vectorv& x) const {
		FIELD t = 0;
		const_cast<Vectorv&>(*this).sort();
		const_cast<Vectorv&>(x).sort();
		auto xit = x.begin();
		for (auto& v : *this) {
			int i = v.first;
			while (xit != x.end() && xit->first < i)
				xit++;
			if (xit == x.end())
				break;
			if (xit->first == i)
				t += v.second * xit->second;
		}
		return t;
	}

public:
	// in-place operations

	Vectorv& operator*=(const FIELD c) {
		for (auto& p : *this)
			p.second *= c;
		return *this;
	}

	Vectorv& operator/=(const FIELD c) {
		for (auto& p : *this)
			p.second /= c;
		return *this;
	}

	Vectorv& operator*=(const Cvector& x) {
		assert(x.n == n);
		for (auto& p : *this)
			p.second *= x(p.first);
		return *this;
	}

	Vectorv& operator/=(const Cvector& x) {
		assert(x.n == n);
		for (auto& p : *this)
			p.second /= x(p.first);
		return *this;
	}

	Vectorv& add(const Vectorv& x, const FIELD c = 1) {
		sort();
		const_cast<Vectorv&>(x).sort();
		vector<SVpair> newbies;
		auto it = begin();
		for (auto& xv : x) {
			while (it != end() && it->first < xv.first)
				it++;
			if (it != end() && it->first == xv.first)
				it->second += c * xv.second;
			else
				newbies.push_back(SVpair(xv.first, c * xv.second));
		}
		if (newbies.size() != 0) {
			vector < SVpair > ::insert(end(), newbies.begin(), newbies.end());
			sorted = 0;
		}
		return *this;
	}

public:
	// rotations
	void apply(const GivensRotation& Q) {
		assert(Q.i1 < n);
		assert(Q.i2 < n);
		FIELD* p1 = findptr(Q.i1);
		FIELD* p2 = findptr(Q.i2);
		if (p1 == &dummyZero) {
			if (p2 == &dummyZero)
				return;
			FIELD x2 = *p2;
			*p2 = Q.cos * x2;
			insert(Q.i1, -Q.sin * x2);
		} else {
			FIELD x1 = *p1;
			if (p2 == &dummyZero) {
				*p1 = Q.cos * x1;
				insert(Q.i2, Q.sin * x1);
			} else {
				*p1 = Q.cos * x1 - Q.sin * (*p2);
				*p2 = Q.sin * x1 + Q.cos * (*p2);
			}
		}
	}
	void applyInverse(const GivensRotation& Q) {
		assert(Q.i1 < n);
		assert(Q.i2 < n);
		FIELD* p1 = findptr(Q.i1);
		FIELD* p2 = findptr(Q.i2);
		if (p1 == &dummyZero) {
			if (p2 == &dummyZero)
				return;
			FIELD x2 = *p2;
			*p2 = Q.cos * x2;
			insert(Q.i1, Q.sin * x2);
		} else {
			FIELD x1 = *p1;
			if (p2 == &dummyZero) {
				*p1 = Q.cos * x1;
				insert(Q.i2, -Q.sin * x1);
			} else {
				*p1 = Q.cos * x1 + Q.sin * (*p2);
				*p2 = -Q.sin * x1 + Q.cos * (*p2);
			}
		}
	}

	void apply(const KpointRotation& Q) { // it is assumed that the indices in Q are sorted
		sort();
		vector<int> missing;
		auto it = begin();
		for (int i = 0; i < Q.k; i++) {
			auto it2 = find_if(it, end(),
					[&Q,i](SVpair& p)->bool {return p.first==Q.ix[i];});
			if (it2 == end())
				missing.push_back(Q.ix[i]);
			else
				it = it2++;
		}
		if (missing.size() == Q.k)
			return;
		for (int i = 0; i < missing.size(); i++)
			insert(missing[i], 0);
		int coi = 0;
		for (auto it = begin(); coi < Q.k; it++) {
			assert(it != end());
			if (it->first == Q.ix[coi]) {
				Q.tempp[coi] = &it->second;
				Q.temp[coi] = it->second;
				coi++;
			}
		}
		for (int i = 0; i < Q.k; i++) {
			double s = 0;
			for (int j = 0; j < Q.k; j++)
				s += Q.q[j * Q.k + i] * Q.temp[j];
			*Q.tempp[i] = s;
		}
	}

	void applyInverse(const KpointRotation& Q) { // it is assumed that the indices in Q are sorted
		sort();
		vector<int> missing;
		auto it = begin();
		for (int i = 0; i < Q.k; i++) {
			auto it2 = find_if(it, end(),
					[&Q,i](SVpair& p)->bool {return p.first==Q.ix[i];});
			if (it2 == end())
				missing.push_back(Q.ix[i]);
			else
				it = it2++;
		}
		if (missing.size() == Q.k)
			return;
		for (int i = 0; i < missing.size(); i++)
			insert(missing[i], 0);
		int coi = 0;
		for (auto it = begin(); coi < Q.k; it++) {
			assert(it != end());
			if (it->first == Q.ix[coi]) {
				Q.tempp[coi] = &it->second;
				Q.temp[coi] = it->second;
				coi++;
			}
		}
		for (int i = 0; i < Q.k; i++) {
			double s = 0;
			for (int j = 0; j < Q.k; j++)
				s += Q.q[i * Q.k + j] * Q.temp[j];
			*Q.tempp[i] = s;
		}
	}

public:
	string str(const Sparse dummy) const;
	string str(const Dense dummy) const {
		return Vector::str(Dense());
	}
	string str() const {
		return str(Sparse());
	}

public:
	bool sorted;
	static FIELD dummyZero;

};
#endif
