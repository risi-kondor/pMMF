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



#include "BlockedRemap.hpp"
#include "Remap.hpp"
#include "Clusters.hpp"

BlockedRemap::BlockedRemap(const BlockedRemap& X) :
		nsource(X.nsource), ndest(X.ndest) {
	forward = new vector<BlockIndexPair>*[nsource];
	for (int i = 0; i < nsource; i++)
		forward[i] = new vector<BlockIndexPair>(*X.forward[i]);
	backward = new vector<BlockIndexPair>*[ndest];
	for (int i = 0; i < ndest; i++)
		backward[i] = new vector<BlockIndexPair>(*X.backward[i]);
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedRemap copied." << endl;
#endif
}

BlockedRemap::BlockedRemap(BlockedRemap&& x) {
	forward = x.forward;
	x.forward = nullptr;
	nsource = x.nsource;
	x.nsource = 0;
	backward = x.backward;
	x.backward = nullptr;
	ndest = x.ndest;
	x.ndest = 0;
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedRemap moved." << endl;
#endif
}

BlockedRemap& BlockedRemap::operator=(const BlockedRemap& X) {
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
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedRemap assigned." << endl;
#endif
	return *this;
}

BlockedRemap& BlockedRemap::operator=(BlockedRemap&& x) {
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
#ifdef _BLOCKEDCOPYWARNING
	cout << "WARNING: BlockedRemap move-assigned." << endl;
#endif
	return *this;
}

BlockedRemap::~BlockedRemap() {
	for (int i = 0; i < nsource; i++)
		delete forward[i];
	delete[] forward;
	for (int i = 0; i < ndest; i++)
		delete backward[i];
	delete[] backward;
}

// ---- Constructors ---------------------------------------------------------------------------------------

BlockedRemap::BlockedRemap(const int _nsource, const int _ndest) :
		nsource(_nsource), ndest(_ndest) {
	forward = new vector<BlockIndexPair>*[nsource];
	for (int i = 0; i < nsource; i++)
		forward[i] = new vector<BlockIndexPair>();
	backward = new vector<BlockIndexPair>*[ndest];
	for (int i = 0; i < ndest; i++)
		backward[i] = new vector<BlockIndexPair>();
}

BlockedRemap::BlockedRemap(const Bstructure& s, const Uninitialized& dummy) :
		BlockedRemap(1, s.size()) {
	int t = 0;
	for (int u : s)
		t += u;
	forward[0]->resize(t);
	for (int i = 0; i < ndest; i++)
		backward[i]->resize(s[i]);
}

BlockedRemap::BlockedRemap(const Bstructure& s, const class Identity& dummy) :
		BlockedRemap(s, Uninitialized()) {
	int t = 0;
	for (int I = 0; I < ndest; I++)
		for (int i = 0; i < s[I]; i++) {
			(*forward[0])[t] = BlockIndexPair(I, i);
			(*backward[I])[i] = BlockIndexPair(0, t++);
		}
}

BlockedRemap BlockedRemap::Identity(const Bstructure& s) {
	BlockedRemap R = BlockedRemap(s, Uninitialized());
	int t = 0;
	for (int I = 0; I < R.ndest; I++)
		for (int i = 0; i < s[I]; i++) {
			(*R.forward[0])[t] = BlockIndexPair(I, i);
			(*R.backward[I])[i] = BlockIndexPair(0, t++);
		}
	return R;
}

BlockedRemap::BlockedRemap(const Bstructure& s1, const Bstructure& s2,
		const Uninitialized& dummy) :
		BlockedRemap(s1.size(), s2.size()) {
	for (int i = 0; i < nsource; i++)
		forward[i]->resize(s1[i]);
	for (int i = 0; i < ndest; i++)
		backward[i]->resize(s2[i]);
}

BlockedRemap::BlockedRemap(const Bstructure& s1, const Bstructure& s2,
		const class Identity& dummy) :
		BlockedRemap(s1, s2, Uninitialized()) {
	BlockedRemap r1 = BlockedRemap::Identity(s1);
	BlockedRemap r2 = BlockedRemap::Identity(s2);
	for (int I = 0; I < nsource; I++)
		for (int i = 0; i < s1[I]; i++) {
			BlockIndexPair& p = (*r2.forward[0])[(*r1.backward[I])[i].index];
			(*forward[I])[i] = p;
			(*backward[p.block])[p.index] = BlockIndexPair(I, i);
		}
}

BlockedRemap BlockedRemap::Identity(const Bstructure& s1,
		const Bstructure& s2) {
	BlockedRemap R = BlockedRemap(s1, s2, Uninitialized());
	BlockedRemap r1 = BlockedRemap::Identity(s1);
	BlockedRemap r2 = BlockedRemap::Identity(s2);
	for (int I = 0; I < R.nsource; I++)
		for (int i = 0; i < s1[I]; i++) {
			BlockIndexPair& p = (*r2.forward[0])[(*r1.backward[I])[i].index];
			(*R.forward[I])[i] = p;
			(*R.backward[p.block])[p.index] = BlockIndexPair(I, i);
		}
	return R;
}

BlockedRemap::BlockedRemap(const Bstructure& s1, const Bstructure& s2,
		const class Random& dummy) :
		BlockedRemap(s1, s2, Uninitialized()) {
	BlockedRemap r1 = BlockedRemap::Identity(s1);
	BlockedRemap r2 = BlockedRemap::Identity(s2);
	Remap rand(r1.forward[0]->size(), ::Random());
	for (int I = 0; I < nsource; I++)
		for (int i = 0; i < s1[I]; i++) {
			BlockIndexPair& p = (*r2.forward[0])[rand(
					(*r1.backward[I])[i].index)];
			(*forward[I])[i] = p;
			(*backward[p.block])[p.index] = BlockIndexPair(I, i);
		}
}

BlockedRemap BlockedRemap::Random(const Bstructure& s1, const Bstructure& s2) {
	BlockedRemap R = BlockedRemap(s1, s2, Uninitialized());
	BlockedRemap r1 = BlockedRemap::Identity(s1);
	BlockedRemap r2 = BlockedRemap::Identity(s2);
	Remap rand(r1.forward[0]->size(), ::Random());
	for (int I = 0; I < R.nsource; I++)
		for (int i = 0; i < s1[I]; i++) {
			BlockIndexPair& p = (*r2.forward[0])[rand(
					(*r1.backward[I])[i].index)];
			(*R.forward[I])[i] = p;
			(*R.backward[p.block])[p.index] = BlockIndexPair(I, i);
		}
	return R;
}

// ---- Other ----------------------------------------------------------------------------------------------
BlockedRemap BlockedRemap::applyTo(const BlockedRemap& x) {
	assert(nsource == x.ndest);
	BlockedRemap result(x.sourceStructure(), destStructure(), Uninitialized());
	for (int u = 0; u < nsource; u++) {
		vector<BlockIndexPair>& bmap = *x.backward[u];
		vector<BlockIndexPair>& fmap = *forward[u];
		assert(fmap.size() == bmap.size());
		for (int i = 0; i < fmap.size(); i++) {
			BlockIndexPair b = bmap[i];
			BlockIndexPair f = fmap[i];
			result.set(b.block, b.index, f.block, f.index);
		}
	}
	return result;
}

void BlockedRemap::eliminateEmptyBlocks() {
	Remap remap(ndest);
	int newndest = ndest - 1;
	for (int i = 0; i < ndest - 1; i++) // the last destination block is never eliminated!
		if (backward[i]->size() == 0) {
			if (remap.backward[i] != newndest - 1)
				remap.swap(remap.backward[i], newndest - 1);
			newndest--;
			delete backward[i];
		}
	if (newndest == ndest - 1)
		return;
	remap.swap(newndest, ndest - 1);
	newndest++;

	vector<BlockIndexPair>** newbackward = new vector<BlockIndexPair>*[newndest];
	for (int i = 0; i < newndest; i++)
		newbackward[i] = backward[remap.forward[i]];
	delete[] backward;
	backward = newbackward;
	ndest = newndest;

	for (int I = 0; I < nsource; I++)
		for (auto& p : *forward[I]) {
			p.block = remap.backward[p.block];
			assert(p.block < ndest);
		}
}

// ---- Summary methods ----------------------------------------------------------------------------------------

BlockStructure BlockedRemap::sourceStructure() const {
	BlockStructure s(nsource);
	for (int i = 0; i < nsource; i++)
		s[i] = forward[i]->size();
	return s;
}

BlockStructure BlockedRemap::destStructure() const {
	BlockStructure s(ndest);
	for (int i = 0; i < ndest; i++)
		s[i] = backward[i]->size();
	return s;
}

int BlockedRemap::getNdest(int n) const {
	if (n == 0 && n > ndest)
		n = ndest;
	int t = 0;
	for (int i = 0; i < n; i++)
		t += backward[i]->size();
	return t;
}

// ---- I/O --------------------------------------------------------------------------------------------------------

string BlockedRemap::str(const bool inverse) const {
	ostringstream stream;
	if (inverse)
		for (int I = 0; I < ndest; I++)
			for (int i = 0; i < backward[I]->size(); i++) {
				BlockIndexPair& p = (*backward[I])[i];
				stream << "(" << I << "," << i << ") <- (" << p.block << ","
						<< p.index << ")" << endl;
			}
	else
		for (int I = 0; I < nsource; I++)
			for (int i = 0; i < forward[I]->size(); i++) {
				BlockIndexPair& p = (*forward[I])[i];
				stream << "(" << I << "," << i << ") -> (" << p.block << ","
						<< p.index << ")" << endl;
			}
	return stream.str();
}

void BlockedRemap::serialize(Rstream& rstream) const {
	rstream << "BlockedRemap{" << Rstream::endl;
	rstream.var("nsource", nsource);
	rstream.var("ndest", ndest);
	rstream.var("forward", " *OMITTED*");
	rstream.var("backward", " *OMITTED*");
	rstream << "}" << Rstream::endl;
}

void BlockedRemap::serialize(Bofstream& ofs) const {
	ofs.tag("BlockedMatrix", 0);
	ofs.write(nsource);
	ofs.write(ndest);
	for (int i = 0; i < nsource; i++)
		ofs.write_vector(*forward[i]);
	for (int i = 0; i < ndest; i++)
		ofs.write_vector(*backward[i]);
}

BlockedRemap::BlockedRemap(Bifstream& ifs) {
	ifs.check("BlockedMatrix", 0);
	ifs.read(nsource);
	forward = new vector<BlockIndexPair>*[nsource];
	ifs.read(ndest);
	backward = new vector<BlockIndexPair>*[ndest];
	for (int i = 0; i < nsource; i++) {
		forward[i] = new vector<BlockIndexPair>;
		ifs.read_vector(*forward[i]);
	}
	for (int i = 0; i < ndest; i++) {
		backward[i] = new vector<BlockIndexPair>;
		ifs.read_vector(*backward[i]);
	}
}

ostream& operator<<(ostream& stream, const BlockedRemap& x) {
	stream << x.str();
	return stream;
}

