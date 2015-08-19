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



#ifndef _MMFstage
#define _MMFstage

// #include <chrono>

#include "BlockedVector.hpp"
#include "BlockedMatrix.hpp"
#include "MMFchannel.hpp"
#include "MMFmatrix.hpp"
#include "MMFremap.hpp"
#include "Cmatrix.hpp"
#include "MultiLoop.hpp"

class MMFstage: public Serializable {
public:

	MMFstage(const MMFstage& x);
	MMFstage(MMFstage&& x);
	MMFstage& operator=(const MMFstage& x);
	MMFstage& operator=(MMFstage&& x);
	~MMFstage() {
		for (auto p : channel)
			delete p;
		delete remap;
	}

public:

	template<class BLOCK>
	MMFstage(const MMFmatrix<BLOCK>& A, MMFremap* _remap = NULL); //, bool _remapHasBypass=false);

	MMFstage(const BlockStructure& structure, const Random& dummy, const int k =
			2);

public:

	template<class VECTOR> pair<BlockedVector<VECTOR>*, VECTOR*> transform(
			BlockedVector<VECTOR>& v);
	template<class VECTOR> BlockedVector<VECTOR> transform(
			BlockedVector<VECTOR>& v, VECTOR& wret);

	template<class VECTOR> BlockedVector<VECTOR>* itransform(
			const BlockedVector<VECTOR>& v, const VECTOR& w);
	template<class VECTOR> BlockedVector<VECTOR> itransformB(
			const BlockedVector<VECTOR>& v, const VECTOR& w);

	template<class BLOCK> pair<BlockedMatrix<BLOCK>*, BlockedMatrix<BLOCK>*>
	transform(BlockedMatrix<BLOCK>& M, const bool fromRight = false) const;

	template<class BLOCK> BlockedMatrix<BLOCK>*
	itransform(const BlockedMatrix<BLOCK>& Msub,
			const BlockedMatrix<BLOCK>& Msplit,
			const bool fromRight = false) const;

	template<class BLOCK> void conjugate(BlockedMatrix<BLOCK>& M,
			const bool inverse = false);

public:

	BlockedCmatrix matrix();

public:

	MMFstage(Bifstream& ifs);
	void serialize(Bofstream& ofs) const;
	void serialize(Rstream& stream) const;

private:

	MMFstage() {
	}

public:
	int nchannels;
	vector<MMFchannel*> channel;
	MMFremap* remap;
	Cvector freqs;

};

// ---- Constructors --------------------------------------------------------------------------------------------
template<class BLOCK> // Copy over information from MMFmatrix...
MMFstage::MMFstage(const MMFmatrix<BLOCK>& A, MMFremap* _remap) : //, bool _remapHasBypass):
		nchannels(A.nchannels), channel(A.nchannels), remap(_remap), freqs(
				_remap->backward[_remap->ndest - 1]->size()) {
	MultiLoop(nchannels,
			[this, &A](int i) {channel[i]=new MMFchannel(*A.process[i]);});

	vector<BlockIndexPair>& wavelets = *remap->backward[remap->ndest - 1];
	for (int i = 0; i < wavelets.size(); i++) {
		BlockIndexPair& p = wavelets[i];
		freqs(i) = A(p.block, p.index, p.block, p.index);
	}
}

// ---- Vector operations ---------------------------------------------------------------------------------------
template<class VECTOR> // destructive!!!
pair<BlockedVector<VECTOR>*, VECTOR*> MMFstage::transform(
		BlockedVector<VECTOR>& v) {

	assert(v.nblocks == nchannels);
	MultiLoop(nchannels,
			[this,&v](const int i) {channel[i]->applyTo(*v.block[i]);});

	assert(remap->nsource == nchannels);
	BlockedVector<VECTOR> vbar(v, *remap);

	int nblocks = vbar.nblocks - 1;
	auto* vdash = new BlockedVector<VECTOR>(nblocks);
	for (int i = 0; i < nblocks; i++)
		vdash->block[i] = vbar.block[i];
	vbar.nodestruct = true;

	return pair<BlockedVector<VECTOR>*, VECTOR*>(vdash, vbar.block[nblocks]);
}

template<class VECTOR> // destructive!!!
BlockedVector<VECTOR> MMFstage::transform(BlockedVector<VECTOR>& v,
		VECTOR& wret) {

	assert(v.nblocks == nchannels);
	MultiLoop(nchannels,
			[this,&v](const int i) {channel[i]->applyTo(*v.block[i]);});

	BlockedVector<VECTOR> vbar(v, *remap);

	int nblocks = vbar.nblocks - 1;
	BlockedVector<VECTOR> vdash(nblocks);
	for (int i = 0; i < nblocks; i++) {
		vdash.block[i] = vbar.block[i];
		vbar.block[i] = nullptr;
	}
	wret = std::move(*vbar.block[nblocks]);

	return vdash;
}

template<class VECTOR>
BlockedVector<VECTOR>* MMFstage::itransform(const BlockedVector<VECTOR>& v,
		const VECTOR& w) {

	int nblocks = v.nblocks + 1;
	BlockedVector<VECTOR> vbar(nblocks, true);
	for (int i = 0; i < nblocks - 1; i++)
		vbar.block[i] = v.block[i];
	vbar.block[nblocks - 1] = const_cast<VECTOR*>(&w);

	assert(remap->ndest == nblocks);
	auto* vdash = new BlockedVector<VECTOR>(vbar, *remap, true);

	assert(vdash->nblocks == nchannels);
	MultiLoop(nchannels,
			[this,&vdash](const int i) {channel[i]->applyInverseTo(*vdash->block[i]);});

	return vdash;
}

template<class VECTOR>
BlockedVector<VECTOR> MMFstage::itransformB(const BlockedVector<VECTOR>& v,
		const VECTOR& w) {

	int nblocks = v.nblocks + 1;
	BlockedVector<VECTOR> vbar(nblocks, true);
	for (int i = 0; i < nblocks - 1; i++)
		vbar.block[i] = v.block[i];
	vbar.block[nblocks - 1] = const_cast<VECTOR*>(&w);

	assert(remap->ndest == nblocks);
	BlockedVector<VECTOR> vdash(vbar, *remap, true);
	MultiLoop(nchannels,
			[this,&vdash](const int i) {channel[i]->applyInverseTo(*vdash.block[i]);});

	return vdash;
}

// ---- Matrix Operations ------------------------------------------------------------------------------------
template<class BLOCK> // destructive!!!
pair<BlockedMatrix<BLOCK>*, BlockedMatrix<BLOCK>*> MMFstage::transform(
		BlockedMatrix<BLOCK>& M, const bool fromRight) const {
	int nstr = M.nstreets;
	int ntow = M.ntowers;
	if (!fromRight)
		assert(nstr == nchannels);
	else
		assert(ntow == nchannels);

	MultiLoop(nstr, ntow,
			[this,&M,fromRight](const int i, const int j) {
				mlog.log(4,"  Transforming block [",i,",",j,"] (",M(i,j).nrows,"x",M(i,j).ncols,") ...");
				if(!fromRight) channel[i]->applyFromLeftTo(M(i,j));
				else channel[j]->applyFromRightToT(M(i,j));
				mlog.log(4,"  Block [",i,",",j,"] done.");
			});

	assert(remap->nsource == nchannels);
	mlog.log(3, "Remapping ...");
	BlockedMatrix<BLOCK>* Msub;
	BlockedMatrix<BLOCK>* Msplit;

	if (!fromRight) {
		BlockedMatrix<BLOCK> Mbar(M, *remap, Identity());
		nstr = Mbar.nstreets;
		Mbar.nodestruct = true;

		Msub = new BlockedMatrix<BLOCK>(nstr - 1, ntow);
		for (int i = 0; i < nstr - 1; i++)
			for (int j = 0; j < ntow; j++)
				Msub->ptr(i, j) = Mbar.ptr(i, j);

		Msplit = new BlockedMatrix<BLOCK>(1, ntow);
		for (int j = 0; j < ntow; j++)
			Msplit->ptr(0, j) = Mbar.ptr(nstr - 1, j);
	} else {
		BlockedMatrix<BLOCK> Mbar(M, Identity(), *remap);
		ntow = Mbar.ntowers;
		Mbar.nodestruct = true;

		Msub = new BlockedMatrix<BLOCK>(nstr, ntow - 1);
		for (int i = 0; i < nstr; i++)
			for (int j = 0; j < ntow - 1; j++)
				Msub->ptr(i, j) = Mbar.ptr(i, j);

		Msplit = new BlockedMatrix<BLOCK>(nstr, 1);
		for (int i = 0; i < nstr; i++)
			Msplit->ptr(i, 0) = Mbar.ptr(i, ntow - 1);
	}

	return pair<BlockedMatrix<BLOCK>*, BlockedMatrix<BLOCK>*>(Msub, Msplit);
}

template<class BLOCK>
BlockedMatrix<BLOCK>*
MMFstage::itransform(const BlockedMatrix<BLOCK>& Msub,
		const BlockedMatrix<BLOCK>& Msplit, const bool fromRight) const {
	int nstr = Msub.nstreets + (!fromRight);
	int ntow = Msub.ntowers + fromRight;
	if (!fromRight) {
		assert(Msplit.ntowers == ntow);
		assert(Msplit.nstreets == 1);
		assert(nstr == remap->ndest);
	} else {
		assert(Msplit.ntowers == 1);
		assert(Msplit.nstreets == nstr);
		assert(ntow == remap->ndest);
	}

	BlockedMatrix<BLOCK> Mjoin(nstr, ntow, true);
	if (!fromRight) {
		for (int i = 0; i < nstr - 1; i++)
			for (int j = 0; j < ntow; j++)
				Mjoin.ptr(i, j) = const_cast<BLOCK*>(Msub.ptr(i, j));
		for (int j = 0; j < ntow; j++)
			Mjoin.ptr(nstr - 1, j) = const_cast<BLOCK*>(Msplit.ptr(0, j));
	} else {
		for (int i = 0; i < nstr; i++)
			for (int j = 0; j < ntow - 1; j++)
				Mjoin.ptr(i, j) = const_cast<BLOCK*>(Msub.ptr(i, j));
		for (int i = 0; i < nstr; i++)
			Mjoin.ptr(i, ntow - 1) = const_cast<BLOCK*>(Msplit.ptr(i, 0));
	}

	mlog.log(3, "Remapping ...");
	BlockedMatrix<BLOCK>* M;
	if (!fromRight)
		M = new BlockedMatrix<BLOCK>(Mjoin, *remap, Identity(), true);
	else
		M = new BlockedMatrix<BLOCK>(Mjoin, Identity(), *remap, true);
	nstr = M->nstreets;
	ntow = M->ntowers;

	MultiLoop(nstr, ntow,
			[this,M,fromRight](const int i, const int j) {
				auto& B=(*M)(i,j);
				mlog.log(4,"  Transforming block [",i,",",j,"] (",B.nrows,"x",B.ncols,") ...");
				if(!fromRight) channel[i]->applyFromLeftToT(B);
				else channel[j]->applyFromRightTo(B);
				mlog.log(4,"  Block [",i,",",j,"] done.");
			});

	return M;
}

// ---- Operations -------------------------------------------------------------------------------------------
template<class BLOCK>
void MMFstage::conjugate(BlockedMatrix<BLOCK>& M, const bool inverse) {
	if (!inverse)
		MultiLoop(nchannels, nchannels,
				[this,&M](int i, int j) {
					auto& B=*M.block[j*nchannels+i];
					mlog.log(4,"  Conjugating block [",i,",",j,"] (",B.nrows,"x",B.ncols,") ...");
					channel[i]->applyFromLeftTo(*M.block[j*nchannels+i]);
					channel[j]->applyFromRightToT(*M.block[j*nchannels+i]);
					//this_thread::sleep_for(chrono::milliseconds(300));
					mlog.log(4,"  Block [",i,",",j,"] done.");
				});
	else
		MultiLoop(nchannels, nchannels,
				[this,&M](int i, int j) {
					auto& B=*M.block[j*nchannels+i];
					mlog.log(4,"  Conjugating block [",i,",",j,"] (",B.nrows,"x",B.ncols,") ...");
					channel[i]->applyFromLeftToT(*M.block[j*nchannels+i]);
					channel[j]->applyFromRightTo(*M.block[j*nchannels+i]);
					//this_thread::sleep_for(chrono::milliseconds(300));
					mlog.log(4,"  Block [",i,",",j,"] done.");
				});
}

#endif
