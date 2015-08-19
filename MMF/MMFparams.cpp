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



#include "MMFparams.hpp"

ostream& operator<<(ostream& stream, const MMFparams& params) {
	stream << params.str();
	return stream;
}

string MMFparams::str() const {
	ostringstream stream;
	stream << "k=" << _k << endl;
	stream << "fraction_rotate=" << _frac << endl;
	stream << "n_eliminate_per_rotation=" << _n_eliminate_per_rotation << endl;
	stream << "fraction_eliminate_after=" << _fraction_eliminate_after << endl;
	stream << endl;
	stream << "nsparsestages=" << _nsparsestages << endl;
	stream << "ndensestages=" << _ndensestages << endl;
	stream << "dcore=" << _dcore << endl;
	stream << endl;
	stream << "nclusters=" << _nclusters << endl;
	stream << "minclustersize=" << _minclustersize << endl;
	stream << "maxclustersize=" << _maxclustersize << endl;
	stream << "maxclusterdepth=" << _maxclusterdepth << endl;
	stream << "bypass=" << _bypass << endl;
	stream << endl;
	stream << "prenormalize=" << _prenormalize << endl;
	stream << "selection_normalize=" << _selection_normalize << endl;
	stream << "selection_criterion=" << _selection_criterion << endl;
	stream << "ncoreclusters=" << _ncoreclusters << endl;
	return stream.str();
}
