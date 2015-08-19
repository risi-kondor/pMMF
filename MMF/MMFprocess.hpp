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



#ifndef _MMFprocess
#define _MMFprocess

#include "BlockedMatrix.hpp"
#include "TopkList.hpp"
#include <mutex>

extern Log mlog;

template<class BLOCK>
class MMFprocess {
public:

	MMFprocess(const MMFprocess& x);
	MMFprocess& operator=(const MMFprocess& x);
	~MMFprocess() {
		delete gram;
		for (auto p : givens)
			delete p;
		for (auto p : kpoint)
			delete p;
	}

	MMFprocess(const MMFmatrix<BLOCK>* _owner, const Tower<BLOCK>& _tower,
			const int _myindex);
	void computeGram(const bool selection_normalized = false);

public:

	void computeNormalizer();

	void findGivens(const int nrot, const int nelim = 1,
			const int criterion = 0);
	void findKpoint(const int nrot, const int k, const int nelim = 1,
			const int criterion = 0);
	void findGivens(const double frac, const int nelim = 1,
			const int criterion = 0) {
		findGivens(static_cast<int>(frac * tower.ncols), nelim, criterion);
	}
	void findKpoint(const double frac, const int k, const int nelim = 1,
			const int criterion = 0) {
		findKpoint(static_cast<int>(frac * tower.ncols), k, nelim, criterion);
	}
	void winnow(const double frac);

	vector<vector<int> >* clusterColset(const Tower<BLOCK>& anchors,
			const vector<int>& colset);

public:

	void doGivens(const int i1, const int i2);
	void doKpoint(const IndexSet& I);
	void deactivate(const int r);
	//void removeFromActiveSetSpecialized(const int r);

	template<class MATRIX>
	void applyFromLeftTo(MATRIX& M) const {
		for (auto Qp : givens)
			M.applyFromLeft(*Qp);
		for (auto Qp : kpoint)
			M.applyFromLeft(*Qp);
	}

	template<class MATRIX>
	void applyFromRightToT(MATRIX& M) const {
		for (auto Qp : givens)
			M.applyFromRightT(*Qp);
		for (auto Qp : kpoint)
			M.applyFromRightT(*Qp);
	}

	void foreach_active(std::function<void(INDEX)> lambda) {
		for (int i = 0; i < nactive; i++)
			lambda(activemap(i));
	}

	FIELD offdiag(const int i) const {
		return (*gram)(i, i) - myblock(i, i) * myblock(i, i);
	}

	FIELD quality(const int i1, const int i2, const int criterion) {
		if (criterion == 1) {
			FIELD tr = (*gram)(i1, i1) + (*gram)(i2, i2);
			FIELD det = (*gram)(i1, i1) * (*gram)(i2, i2)
					- (*gram)(i1, i2) * (*gram)(i1, i2);
			return -(tr - sqrt(tr * tr - 4 * det)) / 2;
		}
		if (criterion == 2) {
			double diff = (*gram)(i2, i2) - (*gram)(i1, i1);
			double g = (*gram)(i1, i2);
			if ((float) fabs(100.0 * diff + g) == (float) fabs(100.0 * diff))
				return -min(offdiag(i1), offdiag(i2)); //2.0*g/diff;
			double theta = diff / (2.0 * g);
			double t = 1.0 / (fabs(theta) + sqrt(theta * theta + 1));
			if (theta < 0)
				t = -t;
			double c = 1.0 / sqrt(t * t + 1);
			double s = c * t;
			double gii = (*gram)(i1, i1) * c * c - 2 * (*gram)(i1, i2) * c * s
					+ (*gram)(i2, i2) * s * s;
			double gjj = (*gram)(i2, i2) * c * c + 2 * (*gram)(i1, i2) * c * s
					+ (*gram)(i1, i1) * s * s;
			double aii = myblock(i1, i1) * c * c - 2 * myblock(i1, i2) * c * s
					+ myblock(i2, i2) * s * s;
			double ajj = myblock(i2, i2) * c * c + 2 * myblock(i1, i2) * c * s
					+ myblock(i1, i1) * s * s;
			return -min(gii - aii * aii, gjj - ajj * ajj);
		}
		return -1000;
	}

public:

	MMFmatrix<BLOCK>* owner;
	Tower<BLOCK> tower;
	//GenericTower tower;
	int myindex;
	BLOCK& myblock;

	int nactive;
	Remap activemap;
	vector<bool> activeflag;

	Cmatrix* gram = nullptr;
	FIELD cumulativeError = 0;
	Cvector normalizer;
	Cvector cnorm;

	vector<GivensRotation*> givens;
	vector<KpointRotation*> kpoint;

};

// ---- Constructors -------------------------------------------------------------------------------------------
template<class BLOCK>
MMFprocess<BLOCK>::MMFprocess(const MMFmatrix<BLOCK>* _owner,
		const Tower<BLOCK>& _tower, const int _myindex) :
		owner(const_cast<MMFmatrix<BLOCK>*>(_owner)), tower(0, 0), myindex(
				_myindex), myblock(*_tower.block[myindex]), nactive(
				_tower.ncols), activeflag(_tower.ncols, true), activemap(
				_tower.ncols) {
	tower = std::move(_tower); // why this couldn't go in the init sequence is subtle
}

template<class BLOCK>
void MMFprocess<BLOCK>::computeGram(const bool selection_normalized) {
	mlog.log(4, "  Computing Gram matrix for channel ", myindex, " (n=",
			nactive, ") ...");
	gram = new Cmatrix(tower.dot(tower));
	if (selection_normalized) {
		cnorm = Cvector(tower.ncols);
		for (int i = 0; i < tower.ncols; i++)
			cnorm(i) = sqrt((*gram)(i, i));
	}
	mlog.log(4, "  Gram matrix for channel ", myindex, " complete.");
}

// ---- Copying ------------------------------------------------------------------------------------------------
template<class BLOCK>
MMFprocess<BLOCK>::MMFprocess(const MMFprocess<BLOCK>& x) :
		owner(x.owner), tower(x.tower), myindex(x.myindex), myblock(x.myblock), nactive(
				x.nactive), activeflag(x.activeflag), activemap(x.activemap), gram(
				new Cmatrix(*x.gram)), normalizer(x.normalizer), cnorm(x.cnorm), cumulativeError(
				x.cumulativeError) {
	for (auto p : x.givens)
		givens.push_back(new GivensRotation(*p));
	for (auto p : x.kpoint)
		kpoint.push_back(new KpointRotation(*p));
}

template<class BLOCK>
MMFprocess<BLOCK>& MMFprocess<BLOCK>::operator=(const MMFprocess<BLOCK>& x) {
	owner = x.owner;
	tower = x.tower;
	myindex = x.myindex;
	myblock = x.myblock;
	nactive = x.nactive;
	activeflag = x.activeflag;
	activemap = x.activemap;
	delete gram;
	gram = new Cmatrix(*x.gram);
	cumulativeError = x.cumulativeError;
	normalizer = x.normalizer;
	cnorm = x.cnorm;
	for (auto p : givens)
		delete p;
	givens.clear();
	for (auto p : kpoint)
		delete p;
	kpoint.clear();
	for (auto p : x.givens)
		givens.push_back(new GivensRotation(*p));
	for (auto p : x.kpoint)
		kpoint.push_back(new KpointRotation(*p));
	return *this;
}

// ---- Operations -----------------------------------------------------------------------------------------------
template<class BLOCK>
void MMFprocess<BLOCK>::findGivens(const int nrot, const int nelim,
		const int criterion) {
	mlog.log(4, "  Starting process ", myindex, " (n=", nactive, ") ...");

	for (int rot = 0; rot < nrot && nactive >= 2; rot++) {
		if (owner && owner->nactive <= owner->min_nactive)
			break;

		uniform_int_distribution<int> distri(0, nactive - 1);
		int i1 = activemap(distri(randomNumberGenerator));
		int i2 = activemap(0);
		if (i2 == i1)
			i2 = activemap(1);

		if (criterion) {
			FIELD max = quality(i1, i2, criterion);
			for (int prei = 0; prei < nactive; prei++) {
				int i = activemap(prei);
				FIELD t = quality(i1, i, criterion);
				if (t > max && i != i1) {
					i2 = i;
					max = t;
				}
			}
		} else {
			if (cnorm.n > 0) {
				FIELD max = fabs((*gram)(i2, i1) / cnorm(i2));
				for (int prei = 0; prei < nactive; prei++) {
					int i = activemap(prei);
					FIELD t = fabs((*gram)(i, i1) / cnorm(i));
					if (t > max && i != i1) {
						i2 = i;
						max = t;
					}
				}
			} else {
				FIELD max = fabs((*gram)(i2, i1));
				for (int prei = 0; prei < nactive; prei++) {
					int i = activemap(prei);
					FIELD t = fabs((*gram)(i, i1));
					if (t > max && i != i1) {
						i2 = i;
						max = t;
					}
				}
			}
		}

		doGivens(i1, i2);
		if (nelim == 1) {
			if (offdiag(i1) < offdiag(i2))
				deactivate(i1);
			else
				deactivate(i2);
		}
	}

	mlog.log(4, "  Process ", myindex, " finished.");
}

template<class BLOCK>
void MMFprocess<BLOCK>::findKpoint(const int nrot, const int k, const int nelim,
		const int criterion) {
	mlog.log(4, "  Starting process ", myindex, " (n=", nactive, ") ...");

	for (int rot = 0; rot < nrot; rot++) {
		if (nactive < k)
			break;
		if (owner && owner->nactive <= owner->min_nactive)
			break;

		uniform_int_distribution<int> distri(0, nactive - 1);
		int seed = activemap(distri(randomNumberGenerator));

		TopkList topk(k - 1);
		for (int prei = 0; prei < nactive; prei++) {
			int i = activemap(prei);
			if (i == seed)
				continue;
			if (cnorm.n > 0)
				topk.consider(i, fabs((*gram)(i, seed)) / cnorm(i));
			else
				topk.consider(i, fabs((*gram)(i, seed)));
		}
		IndexSet I(k);
		I.ix[0] = seed;
		{
			int i = 1;
			for (auto& p : topk)
				I[i++] = p.first;
		}
		I.sort();
		doKpoint(I);

		if (nelim == 1) {
			int best = 0;
			FIELD minv = offdiag(I[0]);
			for (int i = 1; i < k; i++) {
				FIELD t = offdiag(I[i]);
				if (t < minv) {
					minv = t;
					best = i;
				}
			}
			deactivate(I[best]);
		}
	}
	mlog.log(4, "  Process ", myindex, " finished.");
}

template<class BLOCK>
void MMFprocess<BLOCK>::doGivens(const int i1, const int i2) {
	if (fabs((*gram)(i1, i2)) < 10e-52)
		return;
	double diff = (*gram)(i2, i2) - (*gram)(i1, i1);
	double g = (*gram)(i1, i2);

	if ((float) fabs(100.0 * diff + g) == (float) fabs(100.0 * diff))
		return; //2.0*g/diff;
	double theta = diff / (2.0 * g);
	double t = 1.0 / (fabs(theta) + sqrt(theta * theta + 1));
	if (theta < 0)
		t = -t;
	double c = 1.0 / sqrt(t * t + 1);
	GivensRotation* newQ = new GivensRotation(i1, i2, c, t * c);
	if (mlog.verbosity > 5) {
		CoutLock lock;
		mlog.skip(1);
		cout << "PROCESS " << myindex << endl << newQ->str() << endl;
	}
	givens.push_back(newQ);

	myblock.applyFromLeft(*newQ);
	myblock.applyFromRightT(*newQ);
	(*gram).applyFromLeft(*newQ);
	(*gram).applyFromRightT(*newQ);
	if (cnorm.n > 0) {
		cnorm(i1) = sqrt((*gram)(i1, i1));
		cnorm(i2) = sqrt((*gram)(i2, i2));
	}
}

template<class BLOCK>
void MMFprocess<BLOCK>::doKpoint(const IndexSet& I) {

	Cmatrix T = gram->submatrix(I, I);
	pair<Cmatrix*, Cvector*> p = T.symmetricEigensolver();
	KpointRotation* newQ = new KpointRotation(I, *p.first);
	if (mlog.verbosity > 5) {
		CoutLock lock;
		mlog.skip(1);
		mlog.skippedlines = 1;
		cout << "PROCESS " << myindex << endl << newQ->str() << endl;
	}
	delete p.first;
	delete p.second;
	kpoint.push_back(newQ);

	myblock.applyFromLeft(*newQ);
	myblock.applyFromRightT(*newQ);
	(*gram).applyFromLeft(*newQ);
	(*gram).applyFromRightT(*newQ);
	if (cnorm.n > 0)
		for (int i = 0; i < I.k; i++)
			cnorm(i) = sqrt((*gram)(i, i));
}

template<class BLOCK>
void MMFprocess<BLOCK>::winnow(const double frac) {
	TopkList topk(static_cast<int>(frac * nactive));
	foreach_active([this,&topk](const int i) {topk.consider(i,-offdiag(i));});
	for (auto& p : topk)
		deactivate(p.first);
}

template<class BLOCK>
void MMFprocess<BLOCK>::deactivate(const int r) {

	if (mlog.verbosity > 5) {
		CoutLock lock;
		mlog.skip(1);
		mlog.skippedlines = 1;
		cout << "PROCESS " << myindex << endl << "Eliminating column " << r
				<< endl << endl;
	}

	activeflag[r] = false;
	if (activemap.backward[r] != nactive - 1)
		activemap.swap(activemap.backward[r], nactive - 1);
	nactive--;
	if (owner)
		owner->nactive--;

	myblock.foreach_in_column(r, [this,r](const int i1, FIELD v1) {
		myblock.foreach_in_column(r,[this,i1,v1](const int i2, FIELD v2) {
					(*gram)(i1,i2)-=v1*v2;});
		if(cnorm.n>0) cnorm(i1)=sqrt((*gram)(i1,i1));});
}

template<class BLOCK>
vector<vector<int> >*
MMFprocess<BLOCK>::clusterColset(const Tower<BLOCK>& anchors,
		const vector<int>& colset) {
	assert(anchors.nblocks == tower.nblocks);
	vector < vector<int> > *clusters = new vector<vector<int> >(anchors.ncols);
	int ncolset = colset.size();
	for (int i = 0; i < ncolset; i++) {
		Cvector v = anchors.block[0]->dot(tower.block[0]->vcolumn(colset[i]));
		for (int u = 1; u < anchors.nblocks; u++)
			v.add(anchors.block[u]->dot(tower.block[u]->vcolumn(colset[i])));
		FIELD max = v(v.argmax_abs());
		vector<int> maximizers;
		for (int j = 0; j < v.n; j++)
			if (v(j) == max)
				maximizers.push_back(j);
		uniform_int_distribution<int> distri(0, maximizers.size() - 1);
		(*clusters)[maximizers[distri(randomNumberGenerator)]].push_back(
				colset[i]);
	}
	return clusters;
}

template<class BLOCK>
void MMFprocess<BLOCK>::computeNormalizer() {
	int nblocks = tower.nblocks;
	int ncols = tower.ncols;
	normalizer = Cvector(ncols);
	Cmatrix T(nblocks, ncols);
	MultiLoop(nblocks,
			[this,ncols,&T](int I) {for(int j=0; j<ncols; j++) T(I,j)=tower.block[I]->vcolumn(j).norm2();});
	for (int j = 0; j < ncols; j++) {
		FIELD t = 0;
		for (int I = 0; I < nblocks; I++)
			t += T(I, j);
		if (t == 0) {
			activeflag[j] = false;
			t = 1;
		}
		normalizer(j) = sqrt(sqrt(t));
	}
}

#endif
