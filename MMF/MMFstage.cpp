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



#include "MMFstage.hpp"
#include "BlockedMatrix.hpp"
#include "Cmatrix.hpp"

MMFstage::MMFstage(const MMFstage& x) :
		nchannels(x.nchannels), remap(new MMFremap(*x.remap)), freqs(x.freqs) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFstage copied." << endl;
#endif
	for (auto p : x.channel)
		channel.push_back(new MMFchannel(*p));
}

MMFstage::MMFstage(MMFstage&& x) :
		freqs(std::move(x.freqs)), channel(std::move(x.channel)) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFstage moved." << endl;
#endif
	remap = x.remap;
	x.remap = nullptr;
	nchannels = x.nchannels;
	x.nchannels = 0;
}

MMFstage& MMFstage::operator=(const MMFstage& x) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFstage assigned." << endl;
#endif
	nchannels = x.nchannels;
	delete remap;
	remap = new MMFremap(*x.remap);
	freqs = x.freqs;
	for (auto p : channel)
		delete p;
	channel.clear();
	for (auto p : x.channel)
		channel.push_back(new MMFchannel(*p));
	return *this;
}

MMFstage& MMFstage::operator=(MMFstage&& x) {
#ifdef _MMFCOPYWARNING
	cout << "WARNING: MMFstage move-assigned." << endl;
#endif
	nchannels = x.nchannels;
	x.nchannels = 0;
	delete remap;
	remap = x.remap;
	x.remap = nullptr;
	freqs = std::move(x.freqs);
	for (auto p : channel)
		delete p;
	channel.clear();
	channel = std::move(x.channel);
	return *this;
}

MMFstage::MMFstage(const Bstructure& structure, const Random& dummy,
		const int k) :
		nchannels(structure.size()), channel(structure.size()), remap(NULL) {
	for (int i = 0; i < nchannels; i++)
		channel[i] = new MMFchannel(structure[i], dummy, k);
	remap = new MMFremap(BlockedRemap(structure, structure, Identity()));
}

BlockedCmatrix MMFstage::matrix() {
	Bstructure structure(nchannels);
	for (int i = 0; i < nchannels; i++)
		structure[i] = channel[i]->n;
	BlockedCmatrix R = BlockedCmatrix::Identity(structure);
	for (int I = 0; I < nchannels; I++) {
		auto S = R.virtual_street(I);
		channel[I]->applyFromLeftTo(S);
	}
	return R;
}

void MMFstage::serialize(Rstream& stream) const {
	stream << "MMFstage{" << Rstream::endl;
	stream.var("nchannels", nchannels);
	stream << "  freqs=";
	stream.write(freqs);
	stream << "  remap=";
	stream.write(*remap);
	for (int i = 0; i < nchannels; i++) {
		stream << "  channel[" << i << "]=";
		stream.write(*channel[i]);
	}
	stream << "}" << Rstream::endl;
}

void MMFstage::serialize(Bofstream& ofs) const {
	ofs.tag("MMFstage", 0);
	ofs.write(nchannels);
	ofs.seriate(&freqs);
	ofs.seriate(remap);
	ofs.seriate_vector(channel);
}

MMFstage::MMFstage(Bifstream& ifs) {
	ifs.check("MMFstage", 0);
	ifs.read(nchannels);
	freqs = Cvector(ifs);
	ifs.unseriate(remap);
	ifs.unseriate_vector(channel);
}
