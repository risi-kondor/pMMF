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



#include "Street.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "Cmatrix.hpp"
#include "MatrixX.hpp"

template<>
Street<MatrixXv>::Street(const Street<MatrixXv>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Street((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.nrows) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	vector<BlockIndexPair>* map;
	for (int u = 0; u < nblocks; u++) {
		if (!inverse)
			map = cmap.backward[u];
		else
			map = cmap.forward[u];
		block[u] = new MatrixXv(nrows, map->size(), Nullcols());
		for (int j = 0; j < map->size(); j++) {
			BlockIndexPair& p = (*map)[j];
			assert(p.block < X.nblocks);
			assert(p.index < X.block[p.block]->ncols);
			block[u]->column[j] = new Vectorv(
					*X.block[p.block]->column[p.index]);
		}
	}
}

template<>
Street<MatrixXl>::Street(const Street<MatrixXl>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Street((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.nrows) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	vector<BlockIndexPair>* map;
	for (int u = 0; u < nblocks; u++) {
		if (!inverse)
			map = cmap.backward[u];
		else
			map = cmap.forward[u];
		block[u] = new MatrixXl(nrows, map->size(), Nullcols());
		for (int j = 0; j < map->size(); j++) {
			BlockIndexPair& p = (*map)[j];
			assert(p.block < X.nblocks);
			assert(p.index < X.block[p.block]->ncols);
			block[u]->column[j] = new Vectorl(
					*X.block[p.block]->column[p.index]);
		}
	}
}

template<>
Street<MatrixXh>::Street(const Street<MatrixXh>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Street((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.nrows) {
	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	vector<BlockIndexPair>* map;
	for (int u = 0; u < nblocks; u++) {
		if (!inverse)
			map = cmap.backward[u];
		else
			map = cmap.forward[u];
		block[u] = new MatrixXh(nrows, map->size(), Nullcols());
		for (int j = 0; j < map->size(); j++) {
			BlockIndexPair& p = (*map)[j];
			assert(p.block < X.nblocks);
			assert(p.index < X.block[p.block]->ncols);
			block[u]->column[j] = new Vectorh(
					*X.block[p.block]->column[p.index]);
		}
	}
}

template<>
Street<Cmatrix>::Street(const Street<Cmatrix>& X, const BlockedRemap& cmap,
		const bool inverse, const bool push) :
		Street((!inverse) * cmap.ndest + (inverse) * cmap.nsource, X.nrows) {

	if (!inverse)
		assert(cmap.nsource == X.nblocks);
	else
		assert(cmap.ndest == X.nblocks);
	vector<BlockIndexPair>* map;
	if (push) {
		for (int u = 0; u < nblocks; u++)
			if (!inverse)
				block[u] = new Cmatrix(nrows, cmap.backward[u]->size(), Zero());
			else
				block[u] = new Cmatrix(nrows, cmap.forward[u]->size(), Zero());
		for (int u = 0; u < X.nblocks; u++) {
			if (!inverse)
				map = cmap.forward[u];
			else
				map = cmap.backward[u];
			for (int j = 0; j < map->size(); j++) {
				BlockIndexPair& p = (*map)[j];
				assert(p.block < nblocks);
				assert(p.index < block[p.block]->ncols);
				for (int i = 0; i < nrows; i++)
					block[p.block]->array[p.index * nrows + i] =
							X.block[u]->array[j * nrows + i];
			}
		}
	} else {
		for (int u = 0; u < nblocks; u++) {
			if (!inverse)
				map = cmap.backward[u];
			else
				map = cmap.forward[u];
			block[u] = new Cmatrix(nrows, map->size());
			for (int j = 0; j < map->size(); j++) {
				BlockIndexPair& p = (*map)[j];
				assert(p.block < X.nblocks);
				assert(p.index < X.block[p.block]->ncols);
				for (int i = 0; i < nrows; i++)
					block[u]->array[j * nrows + i] =
							X.block[p.block]->array[p.index * nrows + i];
			}
		}
	}
}

template<>
Street<MatrixXv> Street<MatrixXv>::pullColumns(const BlockedRemap& cmap,
		const bool inverse) {
	int nnewblocks;
	if (!inverse) {
		nnewblocks = cmap.ndest;
		assert(nblocks == cmap.nsource);
	} else {
		nnewblocks = cmap.nsource;
		assert(nblocks == cmap.ndest);
	}
	Street result(nnewblocks, nrows);
	vector<BlockIndexPair>* map;
	for (int u = 0; u < nnewblocks; u++) {
		if (!inverse)
			map = cmap.backward[u];
		else
			map = cmap.forward[u];
		auto* newblock = new MatrixXv(nrows, map->size());
		for (int j = 0; j < map->size(); j++) {
			BlockIndexPair& p = (*map)[j];
			assert(p.block < nblocks);
			assert(p.index < block[p.block]->ncols);
			newblock->column[j] = new Vectorv(*block[p.block]->column[p.index]);
		}
		result.block[u] = newblock;
	}
	return result;
}

template<>
Street<MatrixXl> Street<MatrixXl>::pullColumns(const BlockedRemap& cmap,
		const bool inverse) {
	int nnewblocks;
	if (!inverse) {
		nnewblocks = cmap.ndest;
		assert(nblocks == cmap.nsource);
	} else {
		nnewblocks = cmap.nsource;
		assert(nblocks == cmap.ndest);
	}
	Street result(nnewblocks, nrows);
	vector<BlockIndexPair>* map;
	for (int u = 0; u < nnewblocks; u++) {
		if (!inverse)
			map = cmap.backward[u];
		else
			map = cmap.forward[u];
		auto* newblock = new MatrixXl(nrows, map->size());
		for (int j = 0; j < map->size(); j++) {
			BlockIndexPair& p = (*map)[j];
			assert(p.block < nblocks);
			assert(p.index < block[p.block]->ncols);
			newblock->column[j] = new Vectorl(*block[p.block]->column[p.index]);
		}
		result.block[u] = newblock;
	}
	return result;
}

template<>
Street<MatrixXh> Street<MatrixXh>::pullColumns(const BlockedRemap& cmap,
		const bool inverse) {
	int nnewblocks;
	if (!inverse) {
		nnewblocks = cmap.ndest;
		assert(nblocks == cmap.nsource);
	} else {
		nnewblocks = cmap.nsource;
		assert(nblocks == cmap.ndest);
	}
	Street result(nnewblocks, nrows);
	vector<BlockIndexPair>* map;
	for (int u = 0; u < nnewblocks; u++) {
		if (!inverse)
			map = cmap.backward[u];
		else
			map = cmap.forward[u];
		auto* newblock = new MatrixXh(nrows, map->size());
		for (int j = 0; j < map->size(); j++) {
			BlockIndexPair& p = (*map)[j];
			assert(p.block < nblocks);
			assert(p.index < block[p.block]->ncols);
			newblock->column[j] = new Vectorh(*block[p.block]->column[p.index]);
		}
		result.block[u] = newblock;
	}
	return result;
}

