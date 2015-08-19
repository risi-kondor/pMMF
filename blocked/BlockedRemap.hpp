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



#ifndef _BlockedRemap
#define _BlockedRemap

#include "MMFbase.hpp"

class BlockedRemap: public Serializable {
public:

	BlockedRemap(const BlockedRemap& X);
	BlockedRemap(BlockedRemap&& X);
	BlockedRemap& operator=(const BlockedRemap& X);
	BlockedRemap& operator=(BlockedRemap&& X);
	~BlockedRemap();

public:

	BlockedRemap() {
	}
	BlockedRemap(const int _nsource, const int _ndest);
	BlockedRemap(const Bstructure& s1, const Uninitialized& dummy);
	BlockedRemap(const Bstructure& s1, const Identity& dummy); // obsolete
	BlockedRemap(const Bstructure& s1, const Bstructure& s2,
			const Uninitialized& dummy);
	BlockedRemap(const Bstructure& s1, const Bstructure& s2,
			const class Identity& dummy); // obsolete
	BlockedRemap(const Bstructure& s1, const Bstructure& s2,
			const class Random& dummy); // obsolete

public:

	static BlockedRemap Identity(const Bstructure& s);
	static BlockedRemap Identity(const Bstructure& s1, const Bstructure& s2);
	static BlockedRemap Random(const Bstructure& s1, const Bstructure& s2);

public:

	BlockIndexPair& operator()(const int b, const int i) {
		return (*forward[b])[i];
	}
	BlockIndexPair operator()(const int b, const int i) const {
		return (*forward[b])[i];
	}
	BlockIndexPair& operator()(const BlockIndexPair& p) {
		return (*forward[p.block])[p.index];
	}
	BlockIndexPair operator()(const BlockIndexPair& p) const {
		return (*forward[p.block])[p.index];
	}

	BlockIndexPair& inv(const int b, const int i) {
		return (*backward[b])[i];
	}
	BlockIndexPair inv(const int b, const int i) const {
		return (*backward[b])[i];
	}
	BlockIndexPair& inv(const BlockIndexPair& p) {
		return (*backward[p.block])[p.index];
	}
	BlockIndexPair inv(const BlockIndexPair& p) const {
		return (*backward[p.block])[p.index];
	}

	void set(const int xb, const int xi, const int yb, const int yi) {
		assert(xb < nsource);
		assert(xi < forward[xb]->size());
		assert(yb < ndest);
		assert(yi < backward[yb]->size());
		(*forward[xb])[xi] = BlockIndexPair(yb, yi);
		(*backward[yb])[yi] = BlockIndexPair(xb, xi);
	}

	void eliminateEmptyBlocks();

public:

	BlockedRemap applyTo(const BlockedRemap& x);

	BlockStructure sourceStructure() const;
	BlockStructure destStructure() const;

	int getNdest(int n = 0) const;

public:

	string str(const bool inverse = 0) const;

	BlockedRemap(Bifstream& ifs);
	void serialize(Bofstream& ofs) const;
	void serialize(Rstream& rstream) const;

public:

	int nsource;
	int ndest;

	vector<BlockIndexPair>** forward;
	vector<BlockIndexPair>** backward;

};

ostream& operator<<(ostream& stream, const BlockedRemap& x);

#endif
