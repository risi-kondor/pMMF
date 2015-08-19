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



#ifndef _MMFmatrix
#define _MMFmatrix

#include "BlockedMatrix.hpp"
#include "MMFstageParams.hpp"
#include "MMFprocess.hpp"
#include "MMFremap.hpp"
#include "MultiLoop.hpp"
#include "DenseMatrixFile.hpp"
#include "Clusters.hpp"

extern Log mlog;
extern default_random_engine randomNumberGenerator;

template<class BLOCK>
class MMFmatrix: public BlockedMatrix<BLOCK> {
public:

	MMFmatrix() {
	}
	;
	MMFmatrix(const MMFmatrix<BLOCK>& X);
	MMFmatrix& operator=(const MMFmatrix<BLOCK>& X);
	~MMFmatrix() {
		for (auto* p : process)
			delete p;
	}

public:

	MMFmatrix(const BLOCK& M, const bool nogram = true);
	MMFmatrix(BLOCK& M, const bool nogram = true, const bool nocopy = false);

	MMFmatrix(const int n, const Random dummy, const bool nogram = true);
	MMFmatrix static RandomPosDef(const int n, const bool nogram = true);

	template<class BLOCK2>
	MMFmatrix(const MMFstageParams& params, const MMFmatrix<BLOCK2>& A,
			const MMFremap& remap, const bool nogram);

public:

	MMFmatrix(DenseMatrixFile& file);
	MMFmatrix(SparseMatrixFile& file);
	MMFmatrix(DenseMatrixFile& file, const int _nchannels);
	MMFmatrix(SparseMatrixFile& file, const int _nchannels);

public:

	void normalize();
	void findRotations(const double frac, const int k = 2, const int nelim = 1,
			const int criterion = 0);
	void winnow(const double frac);
	void applyRotations(const bool diags = false);
	void computeResidual();
	MMFremap compressionMap(const int nclusters, const int minsize = 10,
			const int maxsize = 1000, const int maxdepth = 0,
			const bool useBypass = false);

private:

	Clusters<BlockIndexPair>
	refineClustering(const vector<BlockIndexPair>& colset, int nclust,
			int minsize, int maxsize, int maxdepth, bool useBypass = false);
	Clusters<BlockIndexPair>
	clusterColset(const Tower<BLOCK>& anchors,
			const vector<BlockIndexPair>& colset);

public:

	int nchannels;
	FIELD residual = 0;
	vector<MMFprocess<BLOCK>*> process;
	atomic<int> nactive;
	int min_nactive = 0;
	bool bypass = false; // if true, then the last channel will have no rotations

};

// ---- Copying ---------------------------------------------------------------------------------------------
template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(const MMFmatrix<BLOCK>& X) :
		BlockedMatrix<BLOCK>(X), nchannels(X.nchannels), bypass(X.bypass), residual(
				X.residual) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFmatrix copied." << endl;
#endif
	for (auto p : X.process)
		process.push_back(new MMFprocess<BLOCK>(*p));
	for (auto p : process)
		p->owner = this;
}

template<class BLOCK>
MMFmatrix<BLOCK>& MMFmatrix<BLOCK>::operator=(const MMFmatrix<BLOCK>& X) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFmatrix assigned." << endl;
#endif
	nchannels = X.nchannels;
	residual = X.residual;
	bypass = X.bypass;
	BlockedMatrix<BLOCK>::operator=(X);
	for (auto p : process)
		delete p;
	process.clear();
	for (auto p : X.process)
		process.push_back(new MMFprocess<BLOCK>(*p));
	for (auto p : process)
		p->owner = this;
	return *this;
}

// ---- Constructors ----------------------------------------------------------------------------------------
template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(const BLOCK& M, const bool nogram) :
		BlockedMatrix<BLOCK>(1, 1), nchannels(1) {
	BlockedMatrix<BLOCK>::block[0] = new BLOCK(M);
	process.push_back(
			new MMFprocess<BLOCK>(this, BlockedMatrix<BLOCK>::virtual_tower(0),
					0, nogram)); // changed
	nactive = process[0]->nactive;
}

template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(BLOCK& M, const bool nogram, const bool nocopy) :
		BlockedMatrix<BLOCK>(1, 1), nchannels(1) {
	if (!nocopy)
		BlockedMatrix<BLOCK>::block[0] = new BLOCK(M);
	else {
		BlockedMatrix<BLOCK>::block[0] = &M;
		BlockedMatrix<BLOCK>::nodestruct = true;
	}
	process.push_back(
			new MMFprocess<BLOCK>(this, BlockedMatrix<BLOCK>::virtual_tower(0),
					0));
	if (!nogram)
		process[0]->computeGram();
	nactive = process[0]->nactive;
}

template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(const int n, const Random dummy, const bool nogram) :
		BlockedMatrix<BLOCK>(1, 1), nchannels(1) {
	BlockedMatrix<BLOCK>::block[0] = new BLOCK(BLOCK::RandomSymmetric(n));
	process.push_back(
			new MMFprocess<BLOCK>(this, BlockedMatrix<BLOCK>::virtual_tower(0),
					0));
	if (!nogram)
		process[0]->computeGram();
	nactive = process[0]->nactive;
}

template<class BLOCK>
template<class BLOCK2>
MMFmatrix<BLOCK>::MMFmatrix(const MMFstageParams& params,
		const MMFmatrix<BLOCK2>& A, const MMFremap& remap, const bool nogram) :
		BlockedMatrix<BLOCK>(A, remap, remap), // this will be obsoleted by MMFremap::hit
		nchannels(remap.ndest), process(remap.ndest), min_nactive(
				params.min_nactive), bypass(remap.has_bypass), // changed
		nactive(0) {
	mlog.skip(4).log(3, "Initializing ", nchannels, " channels...");
	MultiLoop(nchannels,
			[this](const int i) {
				process[i]=new MMFprocess<BLOCK>(this,BlockedMatrix<BLOCK>::virtual_tower(i),i);});
	if (!nogram)
		MultiLoop(nchannels - bypass, [this,params](const int i) {
			process[i]->computeGram(params.selection_normalize);});
	for (int i = 0; i < nchannels - bypass; i++)
		nactive += process[i]->nactive;
}

// ---- Named constructors -----------------------------------------------------------------------------
template<class BLOCK>
MMFmatrix<BLOCK> MMFmatrix<BLOCK>::RandomPosDef(const int n,
		const bool nogram) {
	// BLOCK B=BLOCK::RandomPosDef(n); B.nodestruct=true;
	MMFmatrix<BLOCK> M(*new BLOCK(BLOCK::RandomPosDef(n)), nogram, false);
	return M;
}

// ---- Operations ------------------------------------------------------------------------------------------------
template<class BLOCK>
void MMFmatrix<BLOCK>::normalize() {
	MultiLoop(nchannels,
			[this](const int I) {process[I]->computeNormalizer();});
	MultiLoop(nchannels, nchannels, [this](const int I, const int J) {
		auto& M=BlockedMatrix<BLOCK>::operator()(I,J);
		M.divideRowsBy(process[I]->normalizer);
		M.divideColsBy(process[I]->normalizer);
	});
}

template<class BLOCK>
void MMFmatrix<BLOCK>::findRotations(const double frac, const int k,
		const int nelim, const int criterion) {
	assert(process.size() == nchannels);
	if (k == 2)
		MultiLoop(nchannels - bypass, [this,frac,nelim,criterion](const int i) {
			process[i]->findGivens(frac,nelim,criterion);});
	else
		MultiLoop(nchannels - bypass,
				[this,frac,k,nelim,criterion](const int i) {
					process[i]->findKpoint(frac,k,nelim,criterion);});
}

template<class BLOCK>
void MMFmatrix<BLOCK>::winnow(const double frac) {
	assert(process.size() == nchannels);
	MultiLoop(nchannels - bypass,
			[this,frac](const int i) {process[i]->winnow(frac);});
}

template<class BLOCK>
void MMFmatrix<BLOCK>::applyRotations(const bool diags) {
	MultiLoop(nchannels - bypass, nchannels - bypass,
			[this,diags](int i, int j) {
				if(i==j && !diags) return;
				auto& M=*BlockedMatrix<BLOCK>::block[j*nchannels+i];
				mlog.log(5,"  Conjugating block [",i,",",j,"] (",M.nrows,"x",M.ncols,") ...");
				process[i]->applyFromLeftTo(M);
				process[j]->applyFromRightToT(M);
				mlog.log(5,"  Block [",i,",",j,"] done.");});
}

template<class BLOCK>
void MMFmatrix<BLOCK>::computeResidual() {
	Cmatrix err(nchannels - bypass, nchannels - bypass, Zero());
	MultiLoop(nchannels - bypass, nchannels - bypass,
			[this,&err](int I, int J) {

				FIELD t=0;
				BLOCK& M=*BlockedMatrix<BLOCK>::block[J*nchannels+I];
				vector<bool>& activeflagI=process[I]->activeflag;
				vector<bool>& activeflagJ=process[J]->activeflag;

				for(int j=0; j<M.ncols; j++) {
					if(activeflagJ[j]) continue;
					if(I==J) M.foreach_in_column(j,[j,&activeflagI,&t](const int i, const FIELD v) {
								if(i!=j) t+=v*v*(1+activeflagI[i]);});
					else M.foreach_in_column(j,[&activeflagI,&t](const int i, const FIELD v) {
								t+=v*v*(1+activeflagI[i]);});
				}

				err(I,J)=t;});

	residual = 0;
	for (int I = 0; I < nchannels - bypass; I++)
		for (int J = 0; J < nchannels - bypass; J++)
			residual += err(I, J);
}

// ---- Compression map routines -------------------------------------------------------------------------------
template<class BLOCK>
MMFremap MMFmatrix<BLOCK>::compressionMap(const int nclust, const int minsize,
		const int maxsize, const int maxdepth, const bool _bypass) {

	vector<BlockIndexPair>* actives = new vector<BlockIndexPair>();
	vector<BlockIndexPair>* inactives = new vector<BlockIndexPair>();

	for (int u = 0; u < nchannels; u++) {
		MMFprocess<BLOCK>& pr = *process[u];
		for (int i = 0; i < pr.activeflag.size(); i++)
			if (pr.activeflag[i])
				actives->push_back(BlockIndexPair(u, i));
			else
				inactives->push_back(BlockIndexPair(u, i));
	}

	if (nclust < 0) { // even clustering
		Clusters<BlockIndexPair> clusters;
		for (int u = 0; u < -nclust; u++)
			clusters.push_back(new vector<BlockIndexPair>);
		double mult = static_cast<double>(-nclust) / actives->size();
		for (int i = 0; i < actives->size(); i++)
			clusters[static_cast<int>(mult * i)]->push_back((*actives)[i]);
		delete actives;
		if (_bypass)
			clusters.push_back(new vector<BlockIndexPair>());
		clusters.push_back(inactives);
		return MMFremap(BlockedMatrix<BLOCK>::colstructure(), clusters, _bypass);
	}

	if (maxdepth == 0) { // should be obsoleted with nclust=-1
		//{CoutLock lock; cout<<"SHOULD NOT GET HERE"<<endl;}
		Clusters<BlockIndexPair> clusters;
		clusters.push_back(actives);
		if (_bypass)
			clusters.push_back(new vector<BlockIndexPair>());
		clusters.push_back(inactives);
		return MMFremap(BlockedMatrix<BLOCK>::colstructure(), clusters, _bypass);
	}

	Clusters<BlockIndexPair> clusters = refineClustering(*actives, nclust,
			minsize, maxsize, maxdepth - 1, _bypass);
	delete actives;
	clusters.push_back(inactives);
	return MMFremap(BlockedMatrix<BLOCK>::colstructure(), clusters, _bypass);
}

template<class BLOCK>
Clusters<BlockIndexPair> MMFmatrix<BLOCK>::refineClustering(
		const vector<BlockIndexPair>& colset, int nclust, int minsize,
		int maxsize, int maxdepth, bool useBypass) {
	int N = colset.size();
	maxsize = max(maxsize, static_cast<int>(1.5 * minsize));
	minsize = min(minsize, N);
	nclust = max(nclust, static_cast<int>(1.2 * N / maxsize + 1));
	nclust = max(nclust, 1);
	nclust = min(nclust, static_cast<int>(N / 2));

	mlog.log(4, "nclusters=", nclust);
	mlog.log(4, "minclustersize=", minsize);
	mlog.log(4, "maxclustersize=", maxsize);

	vector<int> sample(nclust);
	uniform_int_distribution<int> distri(0, N - 1);
	for (int i = 0; i < nclust; i++) {
		int s = distri(randomNumberGenerator);
		for (int j = 0; j < i; j++)
			if (sample[j] == s) {
				s = distri(randomNumberGenerator);
				j = 0;
			}
		sample[i] = s;
	}

	BlockedRemap samplingmap(nchannels, 1);
	for (int j = 0; j < nclust; j++)
		samplingmap.backward[0]->push_back(colset[sample[j]]);

	BlockedMatrix<BLOCK> anchors = this->pullColumns(samplingmap);
	Tower<BLOCK> anchorTower = anchors.virtual_tower(0);
	anchorTower.normalizeColumns();
	Clusters<BlockIndexPair> initialClusters = clusterColset(anchorTower,
			colset);
	if (mlog.verbosity >= 5) {
		cout << "Initial:" << endl;
		cout << initialClusters.str() << endl;
	}

	vector<int> goodix;
	Clusters<BlockIndexPair> goodClusters;
	vector<BlockIndexPair>* colset2 = new vector<BlockIndexPair>; //  combined contents of the clusters that are too small

	for (int i = 0; i < nclust; i++) {
		if (initialClusters[i]->size() >= minsize) {
			goodClusters.push_back(initialClusters[i]);
			initialClusters[i] = nullptr;
			goodix.push_back(i);
		} else {
			for (int j = 0; j < initialClusters[i]->size(); j++)
				colset2->push_back((*initialClusters[i])[j]);
		}
	}
	int ngood = goodix.size();

	if (!useBypass) { //                                          scatter colset2 amongst the good clusters
		BlockedRemap samplingmap2(nchannels, 1);
		for (int j = 0; j < ngood; j++)
			samplingmap2.backward[0]->push_back(colset[sample[goodix[j]]]);

		BlockedMatrix<BLOCK> anchors2 = this->pullColumns(samplingmap2);
		Tower<BLOCK> anchorTower2 = anchors2.virtual_tower(0);
		anchorTower2.normalizeColumns();
		goodClusters.add(clusterColset(anchorTower2, *colset2));
	}

	if (mlog.verbosity >= 5) {
		cout << "Good:" << endl;
		cout << goodClusters.str() << endl;
	}

	Clusters<BlockIndexPair> finalClusters;

	for (int i = 0; i < ngood; i++) {
		if (goodClusters[i]->size() <= maxsize || maxdepth < 1) {
			finalClusters.push_back(goodClusters[i]);
			goodClusters[i] = nullptr;
		} else {
			Clusters<BlockIndexPair> refinedClusters = refineClustering(
					*goodClusters[i], nclust, minsize, maxsize, maxdepth - 1,
					useBypass);
			for (int j = 0; j < refinedClusters.size() - useBypass; j++) {
				finalClusters.push_back(refinedClusters[j]);
				refinedClusters[j] = nullptr;
			}
			if (useBypass) {
				vector<BlockIndexPair>& bypassCluster =
						*refinedClusters[refinedClusters.size() - 1];
				for (int j = 0; j < bypassCluster.size(); j++)
					colset2->push_back(bypassCluster[j]);
			}
		}
	}

	if (useBypass)
		finalClusters.push_back(colset2);
	else
		delete colset2;

	if (mlog.verbosity >= 5) {
		cout << "Final:" << endl;
		cout << finalClusters.str() << endl;
	}
	return finalClusters;
}

template<class BLOCK>
Clusters<BlockIndexPair> MMFmatrix<BLOCK>::clusterColset(
		const Tower<BLOCK>& anchors, const vector<BlockIndexPair>& colset) {
	int nclusters = anchors.ncols;
	vector < vector<int> > subcolsets(nchannels);
	for (int i = 0; i < colset.size(); i++)
		subcolsets[colset[i].block].push_back(colset[i].index);

	vector<vector<vector<int> >*> subclusters(nchannels);
	MultiLoop(nchannels, [this,&subcolsets,&anchors,&subclusters](const int u) {
		subclusters[u]=process[u]->clusterColset(anchors,subcolsets[u]);
	});

	Clusters<BlockIndexPair> clusters(nclusters);
	for (int u = 0; u < nchannels; u++) {
		for (int i = 0; i < nclusters; i++) {
			vector<int>& clus = (*subclusters[u])[i];
			for (int j = 0; j < clus.size(); j++)
				clusters[i]->push_back(BlockIndexPair(u, clus[j]));
		}
		delete subclusters[u];
	}

	return clusters;
}

// ---- I/O ---------------------------------------------------------------------------------------------
template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(DenseMatrixFile& file) :
		BlockedMatrix<BLOCK>(1, 1), nchannels(1) {
	BlockedMatrix<BLOCK>::block[0] = new BLOCK(file);
	process.push_back(
			new MMFprocess<BLOCK>(this, BlockedMatrix<BLOCK>::virtual_tower(0),
					0));
}

template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(SparseMatrixFile& file) :
		BlockedMatrix<BLOCK>(1, 1), nchannels(1) {
	BlockedMatrix<BLOCK>::block[0] = new BLOCK(file);
	process.push_back(
			new MMFprocess<BLOCK>(this, BlockedMatrix<BLOCK>::virtual_tower(0),
					0));
}

template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(DenseMatrixFile& file, const int _nchannels) :
		BlockedMatrix<BLOCK>(file, nchannels, nchannels), nchannels(_nchannels) {
	MultiLoop(nchannels,
			[this](const int i) {
				process[i]=new MMFprocess<BLOCK>(this,BlockedMatrix<BLOCK>::virtual_tower(i),i,false);});
}

template<class BLOCK>
MMFmatrix<BLOCK>::MMFmatrix(SparseMatrixFile& file, const int _nchannels) :
		BlockedMatrix<BLOCK>(file, nchannels, nchannels), nchannels(_nchannels) {
	MultiLoop(nchannels,
			[this](const int i) {
				process[i]=new MMFprocess<BLOCK>(this,BlockedMatrix<BLOCK>::virtual_tower(i),i,false);});
}
#endif
