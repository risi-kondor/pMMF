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



#ifndef _Vectorh
#define _Vectorh

#include "SparseVector.hpp"
#include <unordered_map>
#include "GivensRotation.hpp"
#include "KpointRotation.hpp"
#include "Cvector.hpp"

class Vectorh: public SparseVector, public unordered_map<INDEX, FIELD> {
public:
	Vectorh(const int _n) :
			SparseVector(_n) {
	}
	Vectorh(const int _n, const Random& dummy);

public:
	// named constructors
	static Vectorh Random(const int _n, const FIELD p = 0.5);
	static Vectorh Zero(const int _n) {
		return Vectorh(_n);
	}

public:
	// converters
	Vectorh(const Cvector& x);
	Vectorh(const Vectorv& x);
	Vectorh(const Vectorl& x);
	Vectorh(const Vectorl& x, const Remap& remap);

public:
	// element access
	FIELD& operator()(const int i) {
		assert(i < n);
		return (*this)[i];
	}
	FIELD operator()(const int i) const {
		assert(i < n);
		const_iterator it = find(i);
		if (it == end())
			return 0;
		else
			return it->second;
	}
	FIELD read(const int i) const {
		assert(i < n);
		const_iterator it = find(i);
		if (it == end())
			return 0;
		else
			return it->second;
	}
	bool isFilled(const int i) const {
		assert(i < n);
		const_iterator it = find(i);
		return (it != end());
	}
	int nFilled() const {
		return size();
	}

public:
	// sparse methods
	FIELD* findptr(const int i) {
		iterator it = find(i);
		if (it != end())
			return &it->second;
		return &dummyZero;
	}
	void insert(const int i, const FIELD value) {
		assert(i < n);
		(*this)[i] = value;
	}
	void append(const int i, const FIELD value) {
		assert(i < n);
		(*this)[i] = value;
	}
	void zero(const int i) {
		assert(i < n);
		erase(i);
	}

	void sort() {
	}
	void tidy() {
		for (auto it = begin(); it != end(); it++)
			if (it->second == 0)
				it = erase(it);
	}
	void foreach(std::function<void(INDEX, FIELD&)> lambda) {
		for (auto& p : *this)
			lambda(p.first, p.second);
	}

public:
	// scalar operations
	int nnz() const {
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
	FIELD diff2(const Vectorh& x) const {
		assert(x.n == n);
		FIELD t = 0;
		for (auto& xp : x)
			if (xp.second != 0)
				t += (xp.second - read(xp.first))
						* (xp.second - read(xp.first));
		for (auto& p : *this)
			if (x.read(p.first) == 0)
				t += p.second * p.second;
		return t;
	}
	FIELD dot(const Vectorh& x) const {
		FIELD t = 0;
		for (auto& p : x)
			t += p.second * this->read(p.first);
		return t;
	}

public:
	// in-place operations
	Vectorh& operator*=(const FIELD c) {
		for (auto& p : *this)
			p.second *= c;
		return *this;
	}

	Vectorh& operator/=(const FIELD c) {
		for (auto& p : *this)
			p.second /= c;
		return *this;
	}

	Vectorh& operator*=(const Cvector& x) {
		assert(x.n == n);
		for (auto& p : *this)
			p.second *= x(p.first);
		return *this;
	}

	Vectorh& operator/=(const Cvector& x) {
		assert(x.n == n);
		for (auto& p : *this)
			p.second /= x(p.first);
		return *this;
	}

	Vectorh& add(const Vectorh& x, const FIELD c = 1) {
		for (auto& xp : x)
			(*this)[xp.first] += c * xp.second;
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

	void apply(const KpointRotation& Q) { // it is not assumed that the indices in Q are sorted
		vector<int> missing;
		for (int i = 0; i < Q.k; i++)
			if (!isFilled(Q.ix[i]))
				missing.push_back(Q.ix[i]);
		if (missing.size() == Q.k)
			return;
		for (int i = 0; i < missing.size(); i++)
			insert(missing[i], 0);
		for (int i = 0; i < Q.k; i++) {
			auto it = find(Q.ix[i]);
			Q.tempp[i] = &it->second;
			Q.temp[i] = it->second;
		}
		for (int i = 0; i < Q.k; i++) {
			double s = 0;
			for (int j = 0; j < Q.k; j++)
				s += Q.q[j * Q.k + i] * Q.temp[j];
			*Q.tempp[i] = s;
		}
	}
	void applyInverse(const KpointRotation& Q) { // it is not assumed that the indices in Q are sorted
		vector<int> missing;
		for (int i = 0; i < Q.k; i++)
			if (!isFilled(Q.ix[i]))
				missing.push_back(Q.ix[i]);
		if (missing.size() == Q.k)
			return;
		for (int i = 0; i < missing.size(); i++)
			insert(missing[i], 0);
		for (int i = 0; i < Q.k; i++) {
			auto it = find(Q.ix[i]);
			Q.tempp[i] = &it->second;
			Q.temp[i] = it->second;
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
