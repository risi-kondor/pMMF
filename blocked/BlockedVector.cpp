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



#include "BlockedVector.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "Cvector.hpp"
#include "Cmatrix.hpp"
#include "MatrixX.hpp"

template<>
BlockedVector<Vectorv>::BlockedVector(const BlockedVector<Vectorv>& x,
		const BlockedRemap& remap, const bool inverse) :
		BlockedVector((!inverse) * remap.ndest + inverse * remap.nsource) {
	if (!inverse) {
		assert(remap.nsource == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Vectorv(remap.backward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.forward[xu];
			for (auto& p : *x.block[xu]) {
				BlockIndexPair& dest = submap[p.first];
				block[dest.block]->insert(dest.index, p.second);
			}
		}
	} else {
		assert(remap.ndest == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Vectorv(remap.forward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.backward[xu];
			for (auto& p : *x.block[xu]) {
				BlockIndexPair& dest = submap[p.first];
				block[dest.block]->insert(dest.index, p.second);
			}
		}
	}
}

template<>
BlockedVector<Vectorl>::BlockedVector(const BlockedVector<Vectorl>& x,
		const BlockedRemap& remap, const bool inverse) :
		BlockedVector((!inverse) * remap.ndest + inverse * remap.nsource) {
	if (!inverse) {
		assert(remap.nsource == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Vectorl(remap.backward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.forward[xu];
			for (auto& p : *x.block[xu]) {
				BlockIndexPair& dest = submap[p.first];
				block[dest.block]->insert(dest.index, p.second);
			}
		}
	} else {
		assert(remap.ndest == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Vectorl(remap.forward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.backward[xu];
			for (auto& p : *x.block[xu]) {
				BlockIndexPair& dest = submap[p.first];
				block[dest.block]->insert(dest.index, p.second);
			}
		}
	}
}

template<>
BlockedVector<Vectorh>::BlockedVector(const BlockedVector<Vectorh>& x,
		const BlockedRemap& remap, const bool inverse) :
		BlockedVector((!inverse) * remap.ndest + inverse * remap.nsource) {
	if (!inverse) {
		assert(remap.nsource == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Vectorh(remap.backward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.forward[xu];
			for (auto& p : *x.block[xu]) {
				BlockIndexPair& dest = submap[p.first];
				block[dest.block]->insert(dest.index, p.second);
			}
		}
	} else {
		assert(remap.ndest == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Vectorh(remap.forward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.backward[xu];
			for (auto& p : *x.block[xu]) {
				BlockIndexPair& dest = submap[p.first];
				block[dest.block]->insert(dest.index, p.second);
			}
		}
	}
}

template<>
BlockedVector<Cvector>::BlockedVector(const BlockedVector<Cvector>& x,
		const BlockedRemap& remap, const bool inverse) :
		BlockedVector((!inverse) * remap.ndest + inverse * remap.nsource) {
	if (!inverse) {
		assert(remap.nsource == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Cvector(remap.backward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.forward[xu];
			auto arr = x.block[xu]->array;
			int n = x.block[xu]->n;
			for (int i = 0; i < n; i++) {
				BlockIndexPair& dest = submap[i];
				block[dest.block]->array[dest.index] = arr[i];
			}
		}
	} else {
		assert(remap.ndest == x.nblocks);
		for (int u = 0; u < nblocks; u++)
			block[u] = new Cvector(remap.forward[u]->size());
		for (int xu = 0; xu < x.nblocks; xu++) {
			auto& submap = *remap.backward[xu];
			auto arr = x.block[xu]->array;
			int n = x.block[xu]->n;
			for (int i = 0; i < n; i++) {
				BlockIndexPair& dest = submap[i];
				block[dest.block]->array[dest.index] = arr[i];
			}
		}
	}
}
