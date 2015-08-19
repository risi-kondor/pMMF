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



#include "MMFchannel.hpp"
#include "GivensRotation.hpp"
#include "KpointRotation.hpp"

MMFchannel::MMFchannel(const int _n, const Random& dummy, const int k) :
		n(_n) {
	if (k == 2)
		for (int i = 0; i < dummy.n; i++)
			givens.push_back(new GivensRotation(Random(), n));
	else
		for (int i = 0; i < dummy.n; i++)
			kpoint.push_back(new KpointRotation(k, Random(), n));
}

void MMFchannel::serialize(Rstream& stream) const {
	stream << "MMFchannel{" << Rstream::endl;
	stream.var("n", n);
	for (int i = 0; i < givens.size(); i++) {
		stream << "  givens[" << i << "]=";
		stream.write(*givens[i]);
	}
	for (int i = 0; i < kpoint.size(); i++) {
		stream << "  kpoint[" << i << "]=";
	}
	stream << "}" << Rstream::endl;
}

void MMFchannel::serialize(Bofstream& ofs) const {
	ofs.tag("MMFchannel", 0);
	ofs.write(n);
	ofs.seriate_vector(givens);
	ofs.seriate_vector(kpoint);
	ofs.seriate(&normalizer);
	ofs.write(normalized);
}

MMFchannel::MMFchannel(Bifstream& ifs) {
	ifs.check("MMFchannel", 0);
	ifs.read(n);
	ifs.unseriate_vector(givens);
	ifs.unseriate_vector(kpoint);
	normalizer = Cvector(ifs);
	ifs.read(normalized);
}
