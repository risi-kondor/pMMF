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



#ifndef _MMFchannel
#define _MMFchannel

#include "MMFbase.hpp"
#include "MMFprocess.hpp"

class MMFchannel: public Serializable {
public:

	MMFchannel(const MMFchannel& x) :
			n(x.n), normalizer(x.normalizer), normalized(x.normalized) {
#ifdef _MMFCOPYWARNING
		cout << "WARNING: MMFchannel copied." << endl;
#endif
		for (auto p : x.givens)
			givens.push_back(new GivensRotation(*p));
		for (auto p : x.kpoint)
			kpoint.push_back(new KpointRotation(*p));
	}

	MMFchannel& operator=(const MMFchannel& x) {
#ifdef _MMFCOPYWARNING
		cout << "WARNING: MMFchannel assigned." << endl;
#endif
		n = x.n;
		normalizer = x.normalizer;
		normalized = x.normalized;
		givens.clear();
		for (auto p : x.givens)
			givens.push_back(new GivensRotation(*p));
		kpoint.clear();
		for (auto p : x.kpoint)
			kpoint.push_back(new KpointRotation(*p));
		return *this;
	}

	~MMFchannel() {
		for (auto p : givens)
			delete p;
		for (auto p : kpoint)
			delete p;
	}

public:

	MMFchannel(const int _n) :
			n(_n) {
	}
	;

	MMFchannel(const int _n, const Random& dummy, const int k = 2);

	template<class BLOCK>
	MMFchannel(const MMFprocess<BLOCK>& x) :
			n(x.tower.ncols) {
		if (x.normalizer.n > 0) {
			normalizer = x.normalizer;
			normalized = true;
		}
		for (auto p : x.givens)
			givens.push_back(new GivensRotation(*p));
		for (auto p : x.kpoint)
			kpoint.push_back(new KpointRotation(std::move(*p)));
	}

public:
	template<class VECTOR>
	void applyTo(VECTOR& v) const {
		assert(v.n == n);
		if (normalized)
			v *= normalizer;
		for (auto Qp : givens)
			v.apply(*Qp);
		for (auto Qp : kpoint)
			v.apply(*Qp);
	}

	template<class VECTOR>
	void applyInverseTo(VECTOR& v) const {
		assert(v.n == n);
		for (int i = givens.size() - 1; i >= 0; i--)
			v.applyInverse(*givens[i]);
		for (int i = kpoint.size() - 1; i >= 0; i--)
			v.applyInverse(*kpoint[i]);
		if (normalized)
			v *= normalizer;
	}

	template<class MATRIX>
	void applyFromLeftTo(MATRIX& M) const {
		if (normalized)
			M.multiplyRowsBy(normalizer);
		for (auto Qp : givens)
			M.applyFromLeft(*Qp);
		for (auto Qp : kpoint)
			M.applyFromLeft(*Qp);
	}

	template<class MATRIX>
	void applyFromRightToT(MATRIX& M) const {
		if (normalized)
			M.multiplyColsBy(normalizer);
		for (auto Qp : givens)
			M.applyFromRightT(*Qp);
		for (auto Qp : kpoint)
			M.applyFromRightT(*Qp);
	}

	template<class MATRIX>
	void applyFromLeftToT(MATRIX& M) const {
		for (int i = givens.size() - 1; i >= 0; i--)
			M.applyFromLeftT(*givens[i]);
		for (int i = kpoint.size() - 1; i >= 0; i--)
			M.applyFromLeftT(*kpoint[i]);
		if (normalized)
			M.multiplyRowsBy(normalizer);
	}

	template<class MATRIX>
	void applyFromRightTo(MATRIX& M) const {
		if (normalized)
			M.multiplyColsBy(normalizer);
		for (int i = givens.size() - 1; i >= 0; i--)
			M.applyFromRight(*givens[i]);
		for (int i = kpoint.size() - 1; i >= 0; i--)
			M.applyFromRight(*kpoint[i]);
	}

public:

	MMFchannel(Bifstream& ifs);
	void serialize(Bofstream& ofs) const;
	void serialize(Rstream& stream) const;

public:
	int n;
	vector<GivensRotation*> givens;
	vector<KpointRotation*> kpoint;

	Cvector normalizer;
	bool normalized = false;

};

#endif
