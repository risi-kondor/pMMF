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



#ifndef _BlockedVector
#define _BlockedVector

#include "MMFbase.hpp"
#include "BlockedRemap.hpp"
//#include "GenericBlockedVector.hpp"
#include "MultiLoop.hpp"

template<class VECTOR>
class BlockedVector {
public:

	BlockedVector(const BlockedVector<VECTOR>& x);
	BlockedVector(BlockedVector<VECTOR> && x);
	BlockedVector<VECTOR>& operator=(const BlockedVector<VECTOR>& x);
	BlockedVector<VECTOR>& operator=(BlockedVector<VECTOR> && x);
	BlockedVector clone() const;
	~BlockedVector() {
		if (!nodestruct)
			for (int i = 0; i < nblocks; i++)
				delete block[i];
		delete[] block;
	}
	bool nodestruct = false;

public:
	// constructors
	BlockedVector(const int _nblocks, const bool _nodestruct = false);
	BlockedVector(const Bstructure& structure, const Zero& dummy);
	BlockedVector(const Bstructure& structure, const Random& dummy);

public:
	// converters and remapping
	BlockedVector(const VECTOR& x);
	BlockedVector(VECTOR& x, const bool nocopy);
	template<class VECTOR2>
	BlockedVector(const VECTOR2& x, const Bstructure& structure);
	VECTOR vector();

	BlockedVector(const BlockedVector<VECTOR>& x, const BlockedRemap& remap,
			const bool inverse = false);

public:
	// indexing
	FIELD& operator()(const int I, const int i) {
		return (*block[I])(i);
	}
	FIELD operator()(const int I, const int i) const {
		return (*block[I])(i);
	}
	FIELD read(const int I, const int i) const {
		return (*block[I])(i);
	}

	int nnz() const {
		int t = 0;
		for (int i = 0; i < nblocks; i++)
			t += block[i]->nnz();
		return t;
	}

public:
	// virtual views
public:
	BlockedVector<VECTOR> operator-(const BlockedVector<VECTOR>& x) const {
		assert(nblocks == x.nblocks);
		BlockedVector<VECTOR> v(nblocks);
		for (int i = 0; i < nblocks; i++)
			v.block[i] = new VECTOR((*block[i]) - (*x.block[i]));
		return v;
	}

	FIELD dot(const BlockedVector<VECTOR>& x) const {
		assert(x.nblocks == nblocks);
		FIELD t = 0;
		for (int i = 0; i < nblocks; i++)
			t += block[i]->dot(*x.block[i]);
		return t;
	}

	FIELD norm2() const {
		std::vector<FIELD> sub(nblocks);
		MultiLoop(nblocks,
				[this,&sub](const int i) {sub[i]=block[i]->norm2();});
		FIELD T = 0;
		for (FIELD t : sub)
			T += t;
		return T;
	}

	FIELD diff2(const BlockedVector<VECTOR>& x) const {
		assert(nblocks == x.nblocks);
		std::vector<FIELD> sub(nblocks);
		MultiLoop(nblocks,
				[this,&sub,&x](const int i) {sub[i]=block[i]->diff2(*x.block[i]);});
		FIELD T = 0;
		for (FIELD t : sub)
			T += t;
		return T;
	}

public:
	// in-place operations
	void operator*=(const FIELD c) {
		MultiLoop(nblocks, [this,c](const int i) {block[i]->operator*=(c);});
	}

public:

	string str(const Dense dummy) const;
	string str(const Sparse dummy) const;
	string str() const {
		return str(Dense());
	}

public:

	int nblocks;
	VECTOR** block;

};

template<class BLOCK>
ostream& operator<<(ostream& stream, const BlockedVector<BLOCK>& x) {
	stream << x.str();
	return stream;
}

// ---- Constructors ---------------------------------------------------------------------------------------------
template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(const int _nblocks, const bool _nodestruct) :
		nblocks(_nblocks), nodestruct(_nodestruct) {
	block = new VECTOR*[nblocks];
	for (int i = 0; i < nblocks; i++)
		block[1] = NULL;
}

template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(const Bstructure& size, const Zero& dummy) :
		BlockedVector(size.size()) {
	for (int i = 0; i < nblocks; i++)
		block[i] = new VECTOR(size[i], dummy);
}

template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(const Bstructure& size,
		const Random& dummy) :
		BlockedVector(size.size()) {
	for (int i = 0; i < nblocks; i++)
		block[i] = new VECTOR(size[i], dummy);
}

// ---- Copying -----------------------------------------------------------------------------------------

template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(const BlockedVector<VECTOR>& x) :
		BlockedVector(x.nblocks) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedVector copied." << endl;
#endif
	for (int i = 0; i < nblocks; i++)
		block[i] = new VECTOR(*x.block[i]);
}

template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(BlockedVector<VECTOR> && x) :
		BlockedVector(x.nblocks) {
#ifdef _BLOCKEDCOPYWARNING
#endif
	nodestruct = x.nodestruct;
	for (int i = 0; i < nblocks; i++) {
		block[i] = x.block[i];
		x.block[i] = nullptr;
	}
}

template<class VECTOR>
BlockedVector<VECTOR> BlockedVector<VECTOR>::clone() const {
#ifdef _BLOCKEDCOPYWARNING
#endif
	BlockedVector v(nblocks);
	v.nodestruct = true;
	for (int i = 0; i < nblocks; i++)
		v.block[i] = block[i];
	return v;
}

template<class VECTOR>
BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator=(
		const BlockedVector<VECTOR>& x) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedVector assigned." << endl;
#endif
	if (!nodestruct)
		for (int i = 0; i < nblocks; i++)
			delete block[i];
	delete[] block;
	nblocks = x.nblocks;
	block = new VECTOR*[nblocks];
	for (int i = 0; i < nblocks; i++)
		block[i] = new VECTOR(*x.block[i]);
	return *this;
}

template<class VECTOR>
BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator=(
		BlockedVector<VECTOR> && x) {
#ifdef _BLOCKEDCOPYWARNING
#endif
	if (!nodestruct)
		for (int i = 0; i < nblocks; i++)
			delete block[i];
	delete[] block;
	nblocks = x.nblocks;
	block = new VECTOR*[nblocks];
	nodestruct = x.nodestruct;
	for (int i = 0; i < nblocks; i++) {
		block[i] = x.block[i];
		x.block[i] = nullptr;
	}
	return *this;
}

// ---- Conversions ---------------------------------------------------------------------------------------------

template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(const VECTOR& x) :
		BlockedVector(1) {
	block[0] = new VECTOR(x);
}

template<class VECTOR>
BlockedVector<VECTOR>::BlockedVector(VECTOR& x, const bool nocopy) :
		BlockedVector(1) {
	if (nocopy) {
		block[0] = &x;
		nodestruct = true;
	} else
		block[0] = new VECTOR(x);
}

template<class VECTOR>
VECTOR BlockedVector<VECTOR>::vector() {
	int t = 0;
	for (int i = 0; i < nblocks; i++)
		t += block[i]->n;
	VECTOR v(t);
	for (int i = 0; i < nblocks; i++)
		for (auto& p : block[i])
			v.append(p);
	return v;
}

// ---- Operations ----------------------------------------------------------------------------------------------
// ---- I/O -------------------------------------------------------------------------------------------------------

template<class VECTOR>
string BlockedVector<VECTOR>::str(const Dense dummy) const {
	ostringstream stream;
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	stream << "[ ";
	for (int I = 0; I < nblocks; I++) {
		for (int i = 0; i < block[i]->size(); i++) {
			stream << block[I]->read(i) << " ";
		}
		if (I < nblocks - 1)
			stream << "| ";
	}
	stream << "]\n";
	return stream.str();
}

template<class VECTOR>
string BlockedVector<VECTOR>::str(const Sparse dummy) const {
	ostringstream stream;
	for (int I = 0; I < nblocks; I++)
		for (int i = 0; i < block[I]->n; i++)
			stream << "(" << I << "," << i << "): " << block[I]->read(i)
					<< endl;
	return stream.str();
}

#endif
