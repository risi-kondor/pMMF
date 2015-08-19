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



#ifndef _Vectorl
#define _Vectorl

#include "SparseVector.hpp"
#include <list>
#include "GivensRotation.hpp"
#include "KpointRotation.hpp"
#include "Cvector.hpp"

class Vectorl: public SparseVector, public list<SVpair> {
public:

	Vectorl(const int _n) :
			SparseVector(_n) {
	}
	Vectorl(const int _n, const Random& dummy);

public:
	// named constructors

	static Vectorl Random(const int _n, const FIELD p = 0.5);
	static Vectorl Zero(const int _n) {
		return Vectorl(_n);
	}

public:
	// converters
	Vectorl(const Cvector& x);
	Vectorl(const Vectorv& x);
	Vectorl(const Vectorh& x);
	Vectorl(const Vectorl& x, const Remap& remap);

public:
	// element access
	FIELD& operator()(const int i) {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first < i) {
			it++;
		} // list is always ordered
		if (it == end() || it->first > i)
			it = list < SVpair > ::insert(it, SVpair(i, 0));
		return it->second;
	}
	FIELD operator()(const int i) const {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first < i) {
			it++;
		}
		if (it == end() || it->first > i)
			return 0;
		else
			return it->second;
	}
	FIELD read(const int i) const {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first < i) {
			it++;
		}
		if (it == end() || it->first > i)
			return 0;
		else
			return it->second;
	}
	bool isFilled(const int i) const {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first < i) {
			it++;
		}
		return (it != end() && it->first == i);
	}

	int nFilled() const {
		return size();
	}

	void foreach(std::function<void(INDEX, FIELD&)> lambda) {
		for (auto& p : *this)
			lambda(p.first, p.second);
	}

public:
	// sparse methods
	FIELD* findptr(const int i) {
		auto it = begin();
		while (it != end() && it->first < i) {
			it++;
		}
		if (it != end() && it->first == i)
			return &it->second;
		return &dummyZero;
	}
	void insert(const int i, const FIELD value) {
		assert(i < n);
		(*this)(i) = value;
	}
	void append(const int i, const FIELD value) {
		assert(i < n);
		push_back(SVpair(i, value));
	}
	void zero(const int i) {
		assert(i < n);
		auto it = begin();
		while (it != end() && it->first < i) {
			it++;
		}
		if (it != end() && it->first == i)
			erase(it);
	}
	void sort() {
	}
	void tidy() {
		int i = -1;
		for (auto it = begin(); it != end(); it++) {
			while (it->first == i || it->second == 0)
				it = erase(it);
			if (it != end())
				i = it->first;
		}
	}

public:
	// scalar-valued operations
	int nnz() const {
		const_cast<Vectorl&>(*this).tidy();
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
	FIELD diff2(const Vectorl& x) const {
		assert(x.n == n);
		FIELD t = 0;
		auto it = begin();
		for (auto& xv : x) {
			while (it != end() && it->first < xv.first) {
				t += (it->second) * (it->second);
				it++;
			}
			if (it != end() && it->first == xv.first) {
				t += (it->second - xv.second) * (it->second - xv.second);
				it++;
			} else
				t += xv.second * xv.second;
		}
		return t;
	}
	FIELD dot(const Vectorl& x) const {
		FIELD t = 0;
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
	FIELD dot(const Cvector& x) const {
		FIELD t = 0;
		for (auto& v : *this) {
			t += v.second * x(v.first);
			;
		}
		return t;
	}

public:
	// in-place operations
	Vectorl& operator*=(const FIELD c) {
		for (auto& p : *this)
			p.second *= c;
		return *this;
	}
	Vectorl& operator/=(const FIELD c) {
		for (auto& p : *this)
			p.second /= c;
		return *this;
	}
	Vectorl& operator*=(const Cvector& x) {
		assert(x.n == n);
		for (auto& p : *this)
			p.second *= x(p.first);
		return *this;
	}
	Vectorl& operator/=(const Cvector& x) {
		assert(x.n == n);
		for (auto& p : *this)
			p.second /= x(p.first);
		return *this;
	}
	Vectorl& add(const Vectorl& x, const FIELD c = 1) {
		auto it = begin();
		for (auto& xv : x) {
			while (it != end() && it->first < xv.first)
				it++;
			if (it != end() && it->first == xv.first)
				it->second += c * xv.second;
			else {
				list < SVpair > ::insert(it, SVpair(xv.first, c * xv.second));
			}
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

	static FIELD dummyZero;

};

#endif
