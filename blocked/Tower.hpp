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



#ifndef _Tower
#define _Tower

#include "MMFbase.hpp"
#include "BlockedRemap.hpp"
#include "Cmatrix.hpp"
#include "MultiLoop.hpp"

template<class BLOCK>
class Tower {
public:

	Tower(const Tower<BLOCK>& X);
	Tower(Tower<BLOCK> && X);
	Tower<BLOCK>& operator=(const Tower<BLOCK>& X);
	Tower<BLOCK>& operator=(const Tower<BLOCK> && X);
	virtual ~Tower() {
		if (!nodestruct)
			for (int i = 0; i < nblocks; i++)
				delete block[i];
		delete[] block;
	}
	bool nodestruct = false;

public:

	Tower(const int _nblocks, const int _ncols, const bool _nodestruct = false);

	template<class BLOCK2>
	Tower(const Tower<BLOCK2>& X, const BlockedRemap& rmap, const bool inverse =
			false, const bool push = false);

public:
	// element access and virtual views

public:

	Cmatrix dot(const Tower<BLOCK>& x) const;

	template<class VECTOR>
	Cvector dot(const BlockedVector<VECTOR>& v) const;

public:
	// in-place operations

	Tower<BLOCK>& divideColsBy(const Cvector& v) {
		MultiLoop(nblocks, [this,&v](const int I) {block[I]->divideColsBy(v);});
		return *this;
	}

	Tower<BLOCK>& normalizeColumns() {
		Cmatrix norms(nblocks, ncols);
		MultiLoop(nblocks,
				[this,&norms](const int I) {
					for(int j=0; j<ncols; j++) norms(I,j)=block[I]->vcolumn(j).norm2();});
		Cvector normalizers(ncols);
		for (int j = 0; j < ncols; j++) {
			FIELD t = sqrt(norms.columnsum(j));
			if (t != 0)
				normalizers.array[j] = 1.0 / t;
			else
				normalizers.array[j] = 0;
		}
		MultiLoop(nblocks,
				[this,&normalizers](const int I) {
					for(int j=0; j<ncols; j++) block[I]->vcolumn(j)*=normalizers.array[j];});
		return *this;
	}

public:
	// rotations

	template<class ROTATION>
	void applyFromRightT(const ROTATION & Q) {
		MultiLoop(nblocks, [this,&Q](int i) {block[i]->applyFromRightT(Q);});
	}

	void multiplyColsBy(const Cvector& v) {
		for (int u = 0; u < nblocks; u++)
			block[u]->multiplyColsBy(v);
	}

public:

	string str(const Dense dummy) const;
	string str(const Sparse dummy) const;
	string str() const {
		return str(Dense());
	}
	;

public:

	int ncols;
	int nblocks;
	BLOCK** block;

};

template<class BLOCK>
ostream& operator<<(ostream& stream, const Tower<BLOCK>& x) {
	stream << x.str();
	return stream;
}

// ---- Constructors ----------------------------------------------------------------------------------------------
template<class BLOCK>
Tower<BLOCK>::Tower(const int _nblocks, const int _ncols,
		const bool _nodestruct) :
		//GenericTower(_nblocks,_ncols,_nodestruct){
		nblocks(_nblocks), ncols(_ncols), nodestruct(_nodestruct) {
	block = new BLOCK*[nblocks];
	for (int i = 0; i < nblocks; i++)
		block[i] = NULL;
}

// ---- Copying ---------------------------------------------------------------------------------------------------
template<class BLOCK>
Tower<BLOCK>::Tower(const Tower<BLOCK>& X) :
		Tower(X.nblocks, X.ncols, X.nodestruct) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: Tower copied." << endl;
#endif
	for (int i = 0; i < nblocks; i++)
		block[i] = new BLOCK(*X.block[i]);
}

template<class BLOCK>
Tower<BLOCK>::Tower(Tower<BLOCK> && X) :
		Tower(X.nblocks, X.ncols, X.nodestruct) {
#ifdef _BLOCKEDCOPYWARNING
#endif
	for (int i = 0; i < nblocks; i++) {
		block[i] = X.block[i];
		X.block[i] = nullptr;
	}
	nodestruct = X.nodestruct;
}

template<class BLOCK>
Tower<BLOCK>& Tower<BLOCK>::operator=(const Tower<BLOCK>& X) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: Tower assigned." << endl;
#endif
	ncols = X.ncols;
	if (!nodestruct)
		for (int i = 0; i < nblocks; i++)
			delete block[i];
	nblocks = X.nblocks;
	block = new BLOCK*[nblocks];
	for (int i = 0; i < nblocks; i++)
		block[i] = new BLOCK(*X.block[i]);
	return *this;
}

template<class BLOCK>
Tower<BLOCK>& Tower<BLOCK>::operator=(const Tower<BLOCK> && X) {
#ifdef _BLOCKEDCOPYWARNING
#endif
	ncols = X.ncols;
	if (!nodestruct)
		for (int i = 0; i < nblocks; i++)
			delete block[i];
	nblocks = X.nblocks;
	block = new BLOCK*[nblocks];
	for (int i = 0; i < nblocks; i++) {
		block[i] = X.block[i];
		X.block[i] = nullptr;
	}
	nodestruct = X.nodestruct;
	return *this;
}

// ---- Operations ------------------------------------------------------------------------------------------------
template<class BLOCK>
Cmatrix Tower<BLOCK>::dot(const Tower<BLOCK>& x) const {
	Cmatrix M(block[0]->ncols, x.block[0]->ncols, Zero());  // could be improved
	vector<Cmatrix*> T(nblocks);
	MultiLoop(nblocks,
			[this,&T,&x](const int i) {T[i]=new Cmatrix(block[i]->dot(*x.block[i]));});
	for (int u = 0; u < nblocks; u++) {
		M += *T[u];
		delete T[u];
	}
	return M;
}

template<class BLOCK>
template<class VECTOR>
Cvector Tower<BLOCK>::dot(const BlockedVector<VECTOR>& v) const {
	assert(nblocks == v.nblocks);
	Cvector r = block[0]->dot(v);
	for (int u = 1; u < nblocks; u++)
		r.add(block[u]->dot(v));
	return r;
}

// ---- I/O --------------------------------------------------------------------------------------------------------
template<class BLOCK>
string Tower<BLOCK>::str(const Dense dummy) const {
	ostringstream stream;
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	for (int I = 0; I < nblocks; I++) {
		for (int i = 0; i < block[I]->nrows; i++) {
			stream << "[ ";
			for (int j = 0; j < block[I]->ncols; j++)
				stream << block[I]->read(i, j) << " ";
			stream << "]" << endl;
		}
		if (I < nblocks - 1) {
			stream << "[ ";
			for (int j = 0; j < block[I]->ncols; j++)
				stream << "------";
			stream << "]" << endl;
		}
	}
	return stream.str();
}

template<class BLOCK>
string Tower<BLOCK>::str(const Sparse dummy) const {
	ostringstream stream;
	cout << "Sparse output of towers not supported" << endl;
	return stream.str();
}
#endif
