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



#ifndef _MMFremap
#define _MMFremap

#include "BlockedRemap.hpp"
#include "MMFbase.hpp"

// NB: an MMFremap is a special BlockedRemap, where the last destination block is for wavelets

class MMFremap: public BlockedRemap {
public:

	MMFremap(const MMFremap& X);
	MMFremap(MMFremap&& X);
	MMFremap& operator=(const MMFremap& X);
	MMFremap& operator=(MMFremap&& X);

public:

	MMFremap(BlockedRemap&& x);

public:

	MMFremap(const Bstructure& s, const Clusters<BlockIndexPair>& clusters,
			const bool bypass);

	template<class BLOCK>
	BlockedMatrix<BLOCK> hit(const BlockedMatrix<BLOCK>& A, const bool inverse =
			false);

	template<class BLOCK>
	MMFmatrix<BLOCK> hit(const MMFmatrix<BLOCK>& A, const bool inverse = false);

public:

	bool has_bypass;

};

template<class BLOCK>
BlockedMatrix<BLOCK> MMFremap::hit(const BlockedMatrix<BLOCK>& A,
		const bool inverse) {
	ndest--;
	BlockedMatrix<BLOCK> R(A, *this, *this, inverse, inverse);
	ndest++;
	return R;
}
#endif
