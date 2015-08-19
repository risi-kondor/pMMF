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



#include "MMF.hpp"
#include <functional>
#include <algorithm>
#include <map>

// ---- Copying ------------------------------------------------------------------------------------------------
MMF::MMF(const MMF& x) :
		cumulativeError(x.cumulativeError), core(x.core) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMF copied." << endl;
#endif
	for (auto p : x.stage)
		stage.push_back(new MMFstage(*p));
}

MMF::MMF(MMF&& x) :
		cumulativeError(x.cumulativeError), core(std::move(x.core)) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMF moved." << endl;
#endif
	for (auto p : x.stage)
		stage.push_back(p);
	x.stage.clear();
}

MMF& MMF::operator=(const MMF& x) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMF assigned." << endl;
#endif
	cumulativeError = x.cumulativeError;
	core = x.core;
	for (auto p : stage)
		delete p;
	stage.clear();
	for (auto p : x.stage)
		stage.push_back(new MMFstage(*p));
	return *this;
}

MMF& MMF::operator=(MMF&& x) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMF move-assigned." << endl;
#endif
	cumulativeError = x.cumulativeError;
	core = std::move(x.core);
	for (auto p : stage)
		delete p;
	stage.clear();
	for (auto p : x.stage)
		stage.push_back(p);
	x.stage.clear();
	return *this;
}

// ---- Conversions --------------------------------------------------------------------------------------------
MMF MMF::collapse() {
	MMF x;
	x.cumulativeError = cumulativeError;
	x.core = core;
	x.stage.push_back(new MMFstage(*stage[0]));
	MMFstage& xstage = *x.stage[0];
	assert(xstage.nchannels == 1);

	for (int s = 1; s < stage.size(); s++) {

		for (int u = 0; u < stage[s]->nchannels; u++) {
			MMFchannel& ch = *stage[s]->channel[u];
			for (auto p : ch.givens) {
				int i1 = xstage.remap->inv(u, p->i1).index;
				int i2 = xstage.remap->inv(u, p->i2).index;
				xstage.channel[0]->givens.push_back(
						new GivensRotation(i1, i2, p->cos, p->sin));
			}
			for (auto p : ch.kpoint) {
				KpointRotation* newrot = new KpointRotation(*p);
				for (int i = 0; i < newrot->k; i++)
					newrot->ix[i] = xstage.remap->inv(u, newrot->ix[i]).index;
				xstage.channel[0]->kpoint.push_back(newrot);
			}
		}

		xstage.remap->ndest--;
		MMFremap* newremap = new MMFremap(
				stage[s]->remap->applyTo(*xstage.remap));
		xstage.remap->ndest++;

		int newWaveletChannelIndex = newremap->ndest - 1;
		vector<BlockIndexPair>& wavelets =
				*newremap->backward[newWaveletChannelIndex];
		int t = wavelets.size();
		for (auto& p : *xstage.remap->backward[xstage.remap->ndest - 1]) {
			(*newremap)(p.block, p.index) = BlockIndexPair(
					newWaveletChannelIndex, t++);
			wavelets.push_back(p);
		}
		delete xstage.remap;
		xstage.remap = newremap;

		xstage.freqs = Cvector::merge(stage[s]->freqs, xstage.freqs);
	}
	return x;
}

// ---- Other --------------------------------------------------------------------------------------------------
BlockedCmatrix MMF::reconstruction() const {

	BlockedCmatrix* A = const_cast<BlockedCmatrix*>(&core);

	for (int l = stage.size() - 1; l >= 0; l--) {

		if (mlog.verbosity >= 3)
			mlog.skip(3).log(3, "------ Reconstructing stage ", l, " ---------").skip(
					3);
		else
			mlog.log(1, "Reconstructing stage ", l, " ...");

		mlog.log(3, "Remapping ...");
		MMFremap& map = *stage[l]->remap;
		map.ndest--;
		auto* newA = new BlockedCmatrix(*A, map, map, true, true);
		map.ndest++;
		if (l < stage.size() - 1)
			delete A;
		A = newA;

		vector<BlockIndexPair>& wavelets = *map.backward[map.ndest - 1];
		for (int i = 0; i < wavelets.size(); i++) { // more elegant with diag
			BlockIndexPair& p = wavelets[i];
			(*A)(p.block, p.index, p.block, p.index) = stage[l]->freqs(i);
		}

		if (mlog.verbosity > 5) {
			mlog.skip(1);
			cout << A->str() << endl;
		}

		if (l > 0) { // why do we need this condition?
			mlog.skip(4).log(3, "Rotations ...");
			stage[l]->conjugate(*A, true);
		}
	}

	if (mlog.verbosity >= 3)
		mlog.skip(3).log(3, "------ Reconstruction complete --------").skip(1);
	else
		mlog.log(1, "Reconstruction complete.").skip(1);

	return std::move(*A);
}

void MMF::invert(const FIELD eps) {

	if (eps != 0) {
		Cmatrix& A = *core.block[0];
		for (int i = 0; i < A.nrows; i++)
			A(i, i) += eps;
		core = BlockedMatrix<Cmatrix>(core.inverse());
	} else
		core = BlockedMatrix<Cmatrix>(core.inverse());

	cout << "Inversion done" << endl;

	for (auto s : stage) {
		Cvector& freqs = s->freqs;
		for (int i = 0; i < freqs.n; i++) {
			if (freqs(i) < eps)
				freqs(i) = 0;
			else
				freqs(i) = 1.0 / freqs(i);
		}
	}

	for (auto s : stage) {
		MultiLoop(s->nchannels,
				[this,s](int u) {
					MMFchannel& ch=*s->channel[u];
					if(!ch.normalized) return;
					for(int i=0; i<ch.normalizer.n; i++) ch.normalizer.array[i]=1.0/ch.normalizer.array[i];
				});
	}

}

void MMF::invertSqrt(const FIELD eps) {

	if (eps != 0) {
		Cmatrix& A = *core.block[0];
		for (int i = 0; i < A.nrows; i++)
			A(i, i) += eps;
		core = BlockedMatrix<Cmatrix>(core.inverseSqrt());
	} else
		core = BlockedMatrix<Cmatrix>(core.inverseSqrt());

	MultiForeach<MMFstage*>(stage, [eps](MMFstage* s) {
		Cvector& freqs=s->freqs;
		for(int i=0; i<freqs.n; i++)
		freqs(i)=1.0/sqrt(freqs(i)+eps);});

	for (auto s : stage) {
		MultiForeach<MMFchannel*>(s->channel, [](MMFchannel* ch) {
			if(ch->normalized) return;
			for(int i=0; i<ch->normalizer.n; i++)
			ch->normalizer.array[i]=1.0/sqrt(ch->normalizer.array[i]);});
	}

}

// ---- I/O --------------------------------------------------------------------------------------------------------
void MMF::serialize(Rstream& rstream) const {
	rstream << "MMF{" << Rstream::endl;
	rstream.var("cumulativeError", cumulativeError);
	rstream << "  core=";
	rstream.write(core);
	for (int i = 0; i < stage.size(); i++) {
		rstream << "  stage[" << i << "]=";
		rstream.write(*stage[i]);
	}
	rstream << "}" << Rstream::endl;
}

void MMF::serialize(Bofstream& ofs) const {
	ofs.tag("MMF", 0);
	ofs.write(cumulativeError);
	ofs.serialize(core);
	ofs.seriate_vector(stage);
}

MMF::MMF(Bifstream& ifs) {
	ifs.check("MMF", 0);
	ifs.read(cumulativeError);
	ifs.unserialize(core);
	ifs.unseriate_vector(stage);
}

void MMF::graph(Rstream& rstream) const {
	rstream << "digraph G {" << Rstream::endl;
	// assuming the MMF is single-channel
	int qindex = 0;
	std::map<int, int> nonWvlts;
	std::map<int, int> indices, oldIndices;
	for (int i = 0;
			i
					< (*stage.front()->remap->backward[stage.front()->remap->ndest
							- 2]).size(); ++i)
		indices[i] = i;

	for (int sx = 1; sx < stage.size(); sx++) {
		vector<GivensRotation*> rotations =
				((stage[sx]->channel).front())->givens;
		vector<BlockIndexPair> nonwavelets =
				*stage[sx]->remap->backward[stage[sx]->remap->ndest - 2];
		for (int i = 0; i < rotations.size(); i++) {

			if (nonWvlts.find(rotations[i]->i1) != nonWvlts.end()) {
				rstream << "  Q" << qindex << "->" << "Q"
						<< nonWvlts[rotations[i]->i1] << Rstream::endl;
			} else {
				rstream << "  Q" << qindex << "->" << indices[rotations[i]->i1]
						<< Rstream::endl;
			}
			if (nonWvlts.find(rotations[i]->i2) != nonWvlts.end()) {
				rstream << "  Q" << qindex << "->" << "Q"
						<< nonWvlts[rotations[i]->i2] << Rstream::endl;
			} else {
				rstream << "  Q" << qindex << "->" << indices[rotations[i]->i2]
						<< Rstream::endl;
			}
			qindex++;
		}
		nonWvlts.clear();
		oldIndices = indices;
		indices.clear();
		for (int j = 0; j < nonwavelets.size(); ++j) {
			indices[j] = oldIndices[nonwavelets[j].index];
			for (int i = 0; i < rotations.size(); ++i) {
				if (rotations[i]->i1 == nonwavelets[j].index
						|| rotations[i]->i2 == nonwavelets[j].index)
					nonWvlts[j] = qindex - rotations.size() + i;
			}
		}
	}
	rstream << "}" << Rstream::endl;
}
