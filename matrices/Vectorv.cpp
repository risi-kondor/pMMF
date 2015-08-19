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



#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "Cvector.hpp"
#include <random>

extern default_random_engine randomNumberGenerator;

Vectorv::Vectorv(const int _n, const class Random& dummy) :
		SparseVector(_n), sorted(0) {
	uniform_real_distribution<FIELD> distr(0, 1);
	for (int i = 0; i < n; i++)
		if (distr(randomNumberGenerator) <= dummy.p)
			push_back(SVpair(i, distr(randomNumberGenerator)));
}

Vectorv Vectorv::Random(const int _n, const FIELD p) {
	Vectorv v(_n);
	uniform_real_distribution<FIELD> distr(0, 1);
	for (int i = 0; i < _n; i++)
		if (distr(randomNumberGenerator) <= p)
			v.push_back(SVpair(i, distr(randomNumberGenerator)));
	return v;
}

// ---- Conversions -----------------------------------------------------------------------------------------------
Vectorv::Vectorv(const Cvector& x) :
		SparseVector(x.n), sorted(true) {
	for (int i = 0; i < n; i++)
		if (x(i) != 0)
			push_back(SVpair(i, x.array[i]));
}

Vectorv::Vectorv(const Vectorl& x) :
		SparseVector(x.n), sorted(true) {
	for (auto& p : x)
		push_back(p);
}

Vectorv::Vectorv(const Vectorh& x) :
		SparseVector(x.n), sorted(false) {
	for (auto& p : x)
		push_back(SVpair(p.first, p.second));
}

Vectorv::Vectorv(const Vectorv& x, const Remap& remap) :
		SparseVector(x.n), sorted(false) {
	for (auto& p : x)
		push_back(SVpair(remap.forward[p.first], p.second));
}

// ---- I/O ------------------------------------------------------------------------------------------------------
string Vectorv::str(const Sparse dummy) const {
	ostringstream res;
	for (auto& it : *this)
		res << it.first << " : " << it.second << endl;
	return res.str();
}

