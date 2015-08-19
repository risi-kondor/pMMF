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



#ifndef _Street
#define _Street

#include "MMFbase.hpp"
#include "BlockedRemap.hpp"
#include "Cvector.hpp"

template<class BLOCK>
class Street {
public:

	Street(const Street<BLOCK>& X);
	Street(Street<BLOCK> && X);
	Street<BLOCK>& operator=(const Street<BLOCK>& X);
	virtual ~Street() {
		if (!nodestruct)
			for (int i = 0; i < nblocks; i++)
				delete block[i];
		delete[] block;
	}
	bool nodestruct = false;

public:

	Street(const int _nblocks, const int _nrows,
			const bool _nodestruct = false);

	Street(const Street& X, const BlockedRemap& cmap,
			const bool inverse = false, const bool push = false);
	Street<BLOCK> pullColumns(const BlockedRemap& cmap, const bool inverse =
			false);

public:

	BLOCK* operator()(const int J) const {
		assert(J < nblocks);
		return block[J];
	}
	BLOCK*& operator()(const int J) {
		assert(J < nblocks);
		return block[J];
	}

public:

	template<class ROTATION>
	void applyFromLeft(const ROTATION & Q) {
		for (int u = 0; u < nblocks; u++)
			block[u]->applyFromLeft(Q);
	}

	void multiplyRowsBy(const Cvector& v) {
		for (int u = 0; u < nblocks; u++)
			block[u]->multiplyRowsBy(v);
	}

public:

	string str(const Dense dummy) const;
	string str(const Sparse dummy) const;
	string str() const {
		return str(Dense());
	}
	;

public:

	int nrows;
	int nblocks;
	BLOCK** block;

};

template<class BLOCK>
ostream& operator<<(ostream& stream, const Street<BLOCK>& x) {
	stream << x.str();
	return stream;
}

// ---- Constructors ----------------------------------------------------------------------------------------------

template<class BLOCK>
Street<BLOCK>::Street(const int _nblocks, const int _nrows,
		const bool _nodestruct) :
		nblocks(_nblocks), nrows(_nrows), nodestruct(_nodestruct) {
	block = new BLOCK*[nblocks];
	for (int i = 0; i < nblocks; i++)
		block[i] = NULL;
}

// ---- Copying ---------------------------------------------------------------------------------------------------

template<class BLOCK>
Street<BLOCK>::Street(const Street<BLOCK>& X) :
		Street(X.nblocks, X.nrows, X.nodestruct) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: Street copied." << endl;
#endif
	for (int i = 0; i < nblocks; i++)
		block[i] = new BLOCK(*X.block[i]);
}

template<class BLOCK>
Street<BLOCK>::Street(Street<BLOCK> && X) :
		Street(X.nblocks, X.nrows, X.nodestruct) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: Street moved." << endl;
#endif
	for (int i = 0; i < nblocks; i++) {
		block[i] = X.block[i];
		X.block[i] = nullptr;
	}
	nodestruct = X.nodestruct;
}

template<class BLOCK>
Street<BLOCK>& Street<BLOCK>::operator=(const Street<BLOCK>& X) {
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: Street assigned." << endl;
#endif
	nrows = X.nrows;
	if (!nodestruct)
		for (int i = 0; i < nblocks; i++)
			delete block[i];
	nblocks = X.nblocks;
	block = new BLOCK*[nblocks];
	for (int i = 0; i < nblocks; i++)
		block[i] = new BLOCK(*X.block[i]);
	return *this;
}

// ---- I/O ---------------------------------------------------------------------------------------------------------

template<class BLOCK>
string Street<BLOCK>::str(const Dense dummy) const {
	ostringstream stream;
	stream.precision(3);
	stream.setf(ios_base::fixed, ios_base::floatfield);
	for (int i = 0; i < block[0]->nrows; i++) {
		stream << "[ ";
		for (int I = 0; I < nblocks; I++) {
			for (int j = 0; j < block[I]->ncols; j++)
				stream << block[I]->read(i, j) << " ";
			if (I < nblocks - 1)
				stream << " | ";
		}
		stream << "]" << endl;
	}
	return stream.str();
}

template<class BLOCK>
string Street<BLOCK>::str(const Sparse dummy) const {
	ostringstream stream;
	cout << "Sparse output of streets not supported" << endl;
	return stream.str();
}
#endif
