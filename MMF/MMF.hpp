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



#ifndef _MMF
#define _MMF

#include "MMFmatrix.hpp"
#include "MMFstage.hpp"
#include "MMFstageParams.hpp"
#include "Wtransform.hpp"
#include "MMFparams.hpp"
#include "MMFremap.hpp"

extern Log mlog;

class MMF: public Serializable {
public:

	MMF() {
	}
	MMF(const MMF& x);
	MMF(MMF&& x);
	MMF& operator=(const MMF& x);
	MMF& operator=(MMF&& x);
	~MMF() {
		for (auto p : stage)
			delete p;
	}

public:
	// Computing MMF

	template<class BLOCK>
	MMF(const BLOCK& M, const MMFparams& params);

	template<class BLOCK>
	MMF(const MMFmatrix<BLOCK>& M, const MMFparams& params);

public:
	// Conversions
	MMF collapse();

public:
	// Transforming vectors and matrices
	template<class VECTOR> Wtransform<VECTOR> transform(
			const BlockedVector<VECTOR>& v) const;
	template<class VECTOR> BlockedVector<VECTOR> itransform(
			const Wtransform<VECTOR>& W) const;

	template<class BLOCK> vector<BlockedMatrix<BLOCK>*> transform(
			const BlockedMatrix<BLOCK>& M, const bool fromRight = false) const;
	template<class BLOCK> BlockedMatrix<BLOCK> itransform(
			const vector<BlockedMatrix<BLOCK>*>& M,
			const bool fromRight = false) const;

	template<class VECTOR> VECTOR hit(const VECTOR& v) const;
	template<class VECTOR> BlockedVector<VECTOR> hit(
			const BlockedVector<VECTOR>& v) const;
	template<class VECTOR> BlockedVector<VECTOR> operator*(
			const BlockedVector<VECTOR>& v) const {
		return hit(v);
	}

	template<class VECTOR> VECTOR project(const VECTOR& v,
			const double eps = 0) const;
	template<class VECTOR> BlockedVector<VECTOR> project(
			const BlockedVector<VECTOR>& v, const double eps = 0) const;

public:
	void invert(const FIELD eps = 0);
	void invertSqrt(const FIELD eps = 0);

	BlockedCmatrix reconstruction() const;

	template<class BLOCK> FIELD diffSpectral(
			const BlockedMatrix<BLOCK>& M) const;

private:
	template<class BLOCK1, class BLOCK2> MMFmatrix<BLOCK2>*
	computeStage(MMFmatrix<BLOCK1>* oldA, const MMFparams& params,
			const bool deleteold, const bool final);

public:
	MMF(Bifstream& ifs);
	void serialize(Bofstream& ofs) const;
	void serialize(Rstream& rstream) const;
	void graph(Rstream& rstream) const;

public:
	FIELD cumulativeError = 0;
	BlockedCmatrix core;
	vector<MMFstage*> stage;

};

// --- vector methods ---------------------------------------------------------------------------------
template<class VECTOR>
Wtransform<VECTOR> MMF::transform(const BlockedVector<VECTOR>& v0) const {
	int nstages = stage.size();
	assert(nstages > 0);
	Wtransform<VECTOR> result(nstages);

	auto v = v0.clone(); // do it with clone
	for (int l = 0; l < nstages; l++) {
		result.w[l] = new VECTOR(0);
		v = stage[l]->transform(v, *result.w[l]);
	}
	result.v = std::move(v);
	return result;
}

template<class VECTOR>
BlockedVector<VECTOR> MMF::itransform(const Wtransform<VECTOR>& W) const {
	int nstages = stage.size();
	assert(W.w.size() == nstages);
	BlockedVector<VECTOR> v = W.v.clone();
	for (int l = nstages - 1; l >= 0; l--)
		v = stage[l]->itransformB(v, *W.w[l]);
	return v;
}

template<class VECTOR>
BlockedVector<VECTOR> MMF::hit(const BlockedVector<VECTOR>& v) const {
	Wtransform<VECTOR> W = transform(v);
	W.v = core * W.v;
	for (int i = 0; i < W.w.size(); i++)
		(*W.w[i]) *= stage[i]->freqs;
	return itransform(W);
}

template<class VECTOR>
VECTOR MMF::hit(const VECTOR& v) const {
	BlockedVector<VECTOR> V(const_cast<VECTOR&>(v), true);
	return std::move(*hit(V).block[0]);
}

template<class VECTOR>
BlockedVector<VECTOR> MMF::project(const BlockedVector<VECTOR>& v,
		const double eps) const {
	Wtransform<VECTOR> W = transform(v);
	for (int i = 0; i < W.w.size(); i++) {
		Cvector& freqs = stage[i]->freqs;
		for (int j = 0; j < freqs.n; j++)
			if (fabs(freqs.array[j]) <= eps)
				W.w[i]->array[j] = 0;
	}
	return itransform(W);
}

template<class VECTOR>
VECTOR MMF::project(const VECTOR& v, const double eps) const {
	BlockedVector<VECTOR> V(const_cast<VECTOR&>(v), true);
	return std::move(*project(V, eps).block[0]);
}

// --- MMF computation --------------------------------------------------------------------------------
template<class BLOCK>
MMF::MMF(const BLOCK& M, const MMFparams& params) {
	const MMFmatrix<BLOCK> B(const_cast<BLOCK>(M), true, true);
	*this = MMF(B, params);
}

template<class BLOCK>
MMF::MMF(const MMFmatrix<BLOCK>& M, const MMFparams& params) {

	if (mlog.verbosity >= 3) {
		mlog.skip(3).log(3,
				"---------------- MMF ----------------------------");
		cout << endl << params << endl;
	} else
		mlog.log(1, "Starting MMF ...");

	int nsparse = params._nsparsestages;
	int ndense = params._ndensestages;
	MMFmatrix<BLOCK>* Asparse = const_cast<MMFmatrix<BLOCK> *>(&M);
	MMFmatrix<Cmatrix>* Adense;

	if (params._prenormalize) {
		mlog.skip(3).log(3, "Normalizing rows/columns ...");
		Asparse->normalize();
	}

	mlog.skip(3).log(3, "Computing initial clustering ...");
	stage.push_back(
			new MMFstage(*Asparse,
					new MMFremap(
							Asparse->compressionMap(params._nclusters,
									params._minclustersize,
									params._maxclustersize,
									params._maxclusterdepth, params._bypass))));

	for (int s = 0; s < nsparse + ndense; s++) {
		bool final = (s == nsparse + ndense - 1);
		if (s < nsparse)
			Asparse = computeStage<BLOCK, BLOCK>(Asparse, params, s > 0, final);
		if (s == nsparse)
			Adense = computeStage<BLOCK, Cmatrix>(Asparse, params, s > 0,
					final);
		if (s > nsparse)
			Adense = computeStage<Cmatrix, Cmatrix>(Adense, params, s > 0,
					final);
	}

	MMFremap& map = *stage[stage.size() - 1]->remap;
	if (mlog.verbosity >= 3)
		mlog.skip(3).log(0, "------- Computing core matrix (n=",
				map.getNdest(map.ndest - 1), ") -----------").skip(3);
	else
		mlog.log(0, "Computing core matrix, n= ", map.getNdest(map.ndest - 1));
	map.ndest--;
	if (ndense == 0) {
		core = BlockedCmatrix(*Asparse, map, map);
		if (stage.size() > 1)
			delete Asparse;
	} else {
		core = BlockedCmatrix(*Adense, map, map);
		delete Adense;
	}
	map.ndest++;

	mlog.skip(3).log(0, "Total Error = ", cumulativeError);
	if (mlog.verbosity >= 3)
		mlog.skip(3).log(3, "----------- MMF complete ------------------------").skip(
				1);
	else
		mlog.log(1, "MMF complete.").skip(1);
}

template<class BLOCK1, class BLOCK2>
MMFmatrix<BLOCK2>*
MMF::computeStage(MMFmatrix<BLOCK1>* oldA, const MMFparams& params,
		const bool deleteold, const bool final) {

	MMFremap& map = *stage[stage.size() - 1]->remap;
	int Ndest = map.getNdest(map.ndest - 1);

	if (mlog.verbosity >= 3)
		if (typeid(BLOCK2) == typeid(Cmatrix))
			mlog.skip(1).log(0, "----------- MMF stage ",
					static_cast<int>(stage.size()), " (n=", Ndest,
					",Dense) ------------").skip(1);
		else
			mlog.skip(1).log(0, "----------- MMF stage ",
					static_cast<int>(stage.size()), " (n=", Ndest,
					",Sparse) -----------").skip(1);
	else if (typeid(BLOCK2) == typeid(Cmatrix))
		mlog.log(1, "MMF stage ", static_cast<int>(stage.size()), " (n=", Ndest,
				",Dense) ...");
	else
		mlog.log(1, "MMF stage ", static_cast<int>(stage.size()), " (n=", Ndest,
				",Sparse) ...");

	mlog.log(3, "Creating new MMFmatrix ...");
	MMFstageParams stage_params(params);
	stage_params.selection_normalize = params._selection_normalize;
	stage_params.min_nactive = params._dcore;
	map.ndest--;
	auto* A = new MMFmatrix<BLOCK2>(stage_params, *oldA, map, false);
	map.ndest++;

	mlog.skip(4).log(3, "Computing rotations ...");
	A->findRotations(params._frac, params._k, params._n_eliminate_per_rotation,
			params._selection_criterion);

	if (params._fraction_eliminate_after > 0) {
		mlog.skip(4).log(3, "Eliminating rows/columns ...");
		A->winnow(params._fraction_eliminate_after);
	}

	mlog.skip(4).log(3, "Applying rotations ...");
	A->applyRotations(false);

	mlog.skip(4).log(3, "Computing error ...");
	// A->computeResidualSpecialized(); obsolete
	A->computeResidual();
	cumulativeError += A->residual;
	mlog.log(2, "Error of this stage = ", A->residual);

	mlog.skip(4).log(3, "Computing next clustering ...");
	MMFremap* remap;
	if (final)
		remap = new MMFremap(A->compressionMap(-params._ncoreclusters));
	else
		remap = new MMFremap(
				A->compressionMap(params._nclusters, params._minclustersize,
						params._maxclustersize, params._maxclusterdepth,
						params._bypass));
	stage.push_back(new MMFstage(*A, remap));
	if (deleteold)
		delete oldA;
	return A;
}

// --- other ---------------------------------------------------------------------------------------
template<class BLOCK>
FIELD MMF::diffSpectral(const BlockedMatrix<BLOCK>& M) const {
	mlog.log(2, "Computing spectral norm error ...");
	BlockedVector<Cvector> v(M.colstructure(), Random());
	FIELD norm = 0;
	FIELD oldnorm = 0;
	for (int i = 0; i < 200; i++) {
		v = hit(v) - M * v;
		norm = v.norm2();
		mlog.log(4, "norm2=", norm);
		if (fabs((norm - oldnorm) / norm) < 0.001)
			break;
		oldnorm = norm;
		v *= 1.0 / sqrt(norm);
	}
	return norm;
}

#endif

