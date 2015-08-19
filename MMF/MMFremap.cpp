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



#include "MMFremap.hpp"
#include "Clusters.hpp"

MMFremap::MMFremap(const MMFremap& X) {
	nsource = X.nsource;
	forward = new vector<BlockIndexPair>*[nsource];
	for (int i = 0; i < nsource; i++)
		forward[i] = new vector<BlockIndexPair>(*X.forward[i]);
	ndest = X.ndest;
	backward = new vector<BlockIndexPair>*[ndest];
	for (int i = 0; i < ndest; i++)
		backward[i] = new vector<BlockIndexPair>(*X.backward[i]);
	has_bypass = X.has_bypass;
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFremap copied." << endl;
#endif
}

MMFremap::MMFremap(MMFremap&& x) {
	forward = x.forward;
	x.forward = nullptr;
	nsource = x.nsource;
	x.nsource = 0;
	backward = x.backward;
	x.backward = nullptr;
	ndest = x.ndest;
	x.ndest = 0;
	has_bypass = x.has_bypass;
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFremap moved." << endl;
#endif
}

MMFremap& MMFremap::operator=(const MMFremap& X) {
	for (int i = 0; i < nsource; i++)
		delete forward[i];
	delete[] forward;
	for (int i = 0; i < ndest; i++)
		delete backward[i];
	delete[] backward;
	nsource = X.nsource;
	ndest = X.ndest;
	forward = new vector<BlockIndexPair>*[nsource];
	for (int i = 0; i < nsource; i++)
		forward[i] = new vector<BlockIndexPair>(*X.forward[i]);
	backward = new vector<BlockIndexPair>*[ndest];
	for (int i = 0; i < ndest; i++)
		backward[i] = new vector<BlockIndexPair>(*X.backward[i]);
	has_bypass = X.has_bypass;
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFremap assigned." << endl;
#endif
	return *this;
}

MMFremap& MMFremap::operator=(MMFremap&& x) {
	for (int i = 0; i < nsource; i++)
		delete forward[i];
	delete[] forward;
	for (int i = 0; i < ndest; i++)
		delete backward[i];
	delete[] backward;
	forward = x.forward;
	x.forward = nullptr;
	nsource = x.nsource;
	x.nsource = 0;
	backward = x.backward;
	x.backward = nullptr;
	ndest = x.ndest;
	x.ndest = 0;
	has_bypass = x.has_bypass;
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFremap move-assigned." << endl;
#endif
	return *this;
}

// ---- Conversions -------------------------------------------------------------------------------------
MMFremap::MMFremap(BlockedRemap&& x) :
		BlockedRemap(0, 0) {
	forward = x.forward;
	x.forward = nullptr;
	nsource = x.nsource;
	x.nsource = 0;
	backward = x.backward;
	x.backward = nullptr;
	ndest = x.ndest;
	x.ndest = 0;
#ifdef _MMFCOPYWARNING
	cout << "WARNING: BlockedRemap coverted to MMFremap." << endl;
#endif
}

// ---- Other constructors -------------------------------------------------------------------------------
MMFremap::MMFremap(const Bstructure& s,
		const Clusters<BlockIndexPair>& clusters, const bool bypass) :
		BlockedRemap(s, Bstructure(clusters.size()), Uninitialized()) {
	for (int u = 0; u < ndest; u++)
		for (int i = 0; i < clusters[u]->size(); i++) {
			BlockIndexPair p = (*clusters[u])[i];
			backward[u]->push_back(p);
			(*this)(p.block, p.index) = BlockIndexPair(u, i); // check this
		}
	has_bypass = bypass;
}

