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



#include "Tower.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "Cmatrix.hpp"
#include "MatrixX.hpp"

template<>
template<>
Tower<MatrixXv>::Tower(const Tower<MatrixXv>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	for (int u = 0; u < nblocks; u++)
		if (!inverse)
			block[u] = new MatrixXv(cmap.backward[u]->size(), ncols);
		else
			block[u] = new MatrixXv(cmap.forward[u]->size(), ncols);
	vector<BlockIndexPair>* map;
	for (int xu = 0; xu < X.nblocks; xu++) {
		if (!inverse)
			map = cmap.forward[xu];
		else
			map = cmap.backward[xu];
		for (int j = 0; j < ncols; j++) {
			Vectorv& xcol = *X.block[xu]->column[j];
			for (auto& p : xcol) {
				BlockIndexPair& dest = (*map)[p.first];
				if (dest.block >= nblocks)
					continue; // to allow for non-injective maps
				assert(dest.index < block[dest.block]->nrows);
				block[dest.block]->column[j]->insert(dest.index, p.second);
			}
		}
	}
}

template<>
template<>
Tower<MatrixXl>::Tower(const Tower<MatrixXl>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	for (int u = 0; u < nblocks; u++)
		if (!inverse)
			block[u] = new MatrixXl(cmap.backward[u]->size(), ncols);
		else
			block[u] = new MatrixXl(cmap.forward[u]->size(), ncols);
	vector<BlockIndexPair>* map;
	for (int xu = 0; xu < X.nblocks; xu++) {
		if (!inverse)
			map = cmap.forward[xu];
		else
			map = cmap.backward[xu];
		for (int j = 0; j < ncols; j++) {
			Vectorl& xcol = *X.block[xu]->column[j];
			for (auto& p : xcol) {
				BlockIndexPair& dest = (*map)[p.first];
				if (dest.block >= nblocks)
					continue; // to allow for non-injective maps
				assert(dest.index < block[dest.block]->nrows);
				block[dest.block]->column[j]->insert(dest.index, p.second);
			}
		}
	}
}

template<>
template<>
Tower<MatrixXh>::Tower(const Tower<MatrixXh>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	for (int u = 0; u < nblocks; u++)
		if (!inverse)
			block[u] = new MatrixXh(cmap.backward[u]->size(), ncols);
		else
			block[u] = new MatrixXh(cmap.forward[u]->size(), ncols);
	vector<BlockIndexPair>* map;
	for (int xu = 0; xu < X.nblocks; xu++) {
		if (!inverse)
			map = cmap.forward[xu];
		else
			map = cmap.backward[xu];
		for (int j = 0; j < ncols; j++) {
			Vectorh& xcol = *X.block[xu]->column[j];
			for (auto& p : xcol) {
				BlockIndexPair& dest = (*map)[p.first];
				if (dest.block >= nblocks)
					continue; // to allow for non-injective maps
				assert(dest.index < block[dest.block]->nrows);
				block[dest.block]->column[j]->insert(dest.index, p.second);
			}
		}
	}
}

template<> // by default, uses pull method
template<>
Tower<Cmatrix>::Tower(const Tower<Cmatrix>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	vector<BlockIndexPair>* map;
	if (push) {
		for (int u = 0; u < nblocks; u++)
			if (!inverse)
				block[u] = new Cmatrix(cmap.backward[u]->size(), ncols, Zero());
			else
				block[u] = new Cmatrix(cmap.forward[u]->size(), ncols, Zero());
		for (int u = 0; u < X.nblocks; u++) {
			if (!inverse)
				map = cmap.forward[u];
			else
				map = cmap.backward[u];
			for (int i = 0; i < map->size(); i++) {
				BlockIndexPair& p = (*map)[i];
				assert(p.block < nblocks);
				assert(p.index < block[p.block]->nrows);
				for (int j = 0; j < ncols; j++) {
					(*block[p.block])(p.index, j) = (*X.block[u])(i, j);
				}
			}
		}
	} else {
		for (int u = 0; u < nblocks; u++) {
			if (!inverse)
				map = cmap.backward[u];
			else
				map = cmap.forward[u];
			block[u] = new Cmatrix(map->size(), ncols);
			for (int i = 0; i < map->size(); i++) {
				BlockIndexPair& p = (*map)[i];
				assert(p.block < X.nblocks);
				assert(p.index < X.block[p.block]->nrows);
				for (int j = 0; j < ncols; j++) {
					(*block[u])(i, j) = (*X.block[p.block])(p.index, j);
				}
				//block[u]->array[j*nrows+i]=X.block[p.block]->array[j*X.block[p.block]->nrows+p.index];
			}
		}
	}
}

template<>
template<>
Tower<Cmatrix>::Tower(const Tower<MatrixXv>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	for (int u = 0; u < nblocks; u++)
		if (!inverse)
			block[u] = new Cmatrix(cmap.backward[u]->size(), ncols, Zero());
		else
			block[u] = new Cmatrix(cmap.forward[u]->size(), ncols, Zero());
	vector<BlockIndexPair>* map;
	for (int xu = 0; xu < X.nblocks; xu++) {
		if (!inverse)
			map = cmap.forward[xu];
		else
			map = cmap.backward[xu];
		for (int j = 0; j < ncols; j++) {
			Vectorv& xcol = *X.block[xu]->column[j];
			for (auto& p : xcol) {
				BlockIndexPair& dest = (*map)[p.first];
				if (dest.block >= nblocks)
					continue; // to allow for non-injective maps
				assert(dest.index < block[dest.block]->nrows);
				(*block[dest.block])(dest.index, j) = p.second;
			}
		}
	}
}

template<>
template<>
Tower<Cmatrix>::Tower(const Tower<MatrixXl>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	for (int u = 0; u < nblocks; u++)
		if (!inverse)
			block[u] = new Cmatrix(cmap.backward[u]->size(), ncols, Zero());
		else
			block[u] = new Cmatrix(cmap.forward[u]->size(), ncols, Zero());
	vector<BlockIndexPair>* map;
	for (int xu = 0; xu < X.nblocks; xu++) {
		if (!inverse)
			map = cmap.forward[xu];
		else
			map = cmap.backward[xu];
		for (int j = 0; j < ncols; j++) {
			Vectorl& xcol = *X.block[xu]->column[j];
			for (auto& p : xcol) {
				BlockIndexPair& dest = (*map)[p.first];
				if (dest.block >= nblocks)
					continue; // to allow for non-injective maps
				assert(dest.index < block[dest.block]->nrows);
				(*block[dest.block])(dest.index, j) = p.second;
			}
		}
	}
}

template<>
template<>
Tower<Cmatrix>::Tower(const Tower<MatrixXh>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Tower((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.ncols) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	for (int u = 0; u < nblocks; u++)
		if (!inverse)
			block[u] = new Cmatrix(cmap.backward[u]->size(), ncols, Zero());
		else
			block[u] = new Cmatrix(cmap.forward[u]->size(), ncols, Zero());
	vector<BlockIndexPair>* map;
	for (int xu = 0; xu < X.nblocks; xu++) {
		if (!inverse)
			map = cmap.forward[xu];
		else
			map = cmap.backward[xu];
		for (int j = 0; j < ncols; j++) {
			Vectorh& xcol = *X.block[xu]->column[j];
			for (auto& p : xcol) {
				BlockIndexPair& dest = (*map)[p.first];
				if (dest.block >= nblocks)
					continue; // to allow for non-injective maps
				assert(dest.index < block[dest.block]->nrows);
				(*block[dest.block])(dest.index, j) = p.second;
			}
		}
	}
}
