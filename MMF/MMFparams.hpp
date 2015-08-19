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



#ifndef _MMFparams
#define _MMFparams

#include "MMFbase.hpp"

class MMFparams {
public:
	MMFparams& nsparsestages(const int n) {
		_nsparsestages = n;
		return *this;
	}
	MMFparams& ndensestages(const int n) {
		_ndensestages = n;
		return *this;
	}
	MMFparams& nclusters(const int n) {
		_nclusters = n;
		return *this;
	}
	MMFparams& minclustersize(const int n) {
		_minclustersize = n;
		return *this;
	}
	MMFparams& maxclustersize(const int n) {
		_maxclustersize = n;
		return *this;
	}
	MMFparams& maxclusterdepth(const int n) {
		_maxclusterdepth = n;
		return *this;
	}
	MMFparams& fraction(const double f) {
		_frac = f;
		return *this;
	}
	MMFparams& bypass(const bool n) {
		_bypass = n;
		return *this;
	}
	MMFparams& prenormalize(const bool n) {
		_prenormalize = n;
		return *this;
	}
	MMFparams& ncoreclusters(const int n) {
		_ncoreclusters = n;
		return *this;
	}
	MMFparams& dcore(const int n) {
		_dcore = n;
		return *this;
	}
	MMFparams& selection_normalize(const bool p) {
		_selection_normalize = p;
		return *this;
	}
	MMFparams& selection_normalized(const bool p) {
		_selection_normalize = p;
		return *this;
	} // obsolete
	MMFparams& n_eliminate_per_rotation(const int n) {
		_n_eliminate_per_rotation = n;
		return *this;
	}
	MMFparams& fraction_eliminate_after(const double f) {
		_fraction_eliminate_after = f;
		return *this;
	} // winnow
	MMFparams& selection_criterion(const int n) {
		_selection_criterion = n;
		return *this;
	}
	static MMFparams k(const int n) {
		MMFparams params;
		params._k = n;
		return params;
	}
	string str() const;
public:
	int _k = 2;
	double _frac = 0.5;
	int _n_eliminate_per_rotation = 1;
	double _fraction_eliminate_after = 0;

	int _nsparsestages = 1;
	int _ndensestages = 1;
	int _dcore = 0;

	int _nclusters = 2;
	int _minclustersize = 1;
	int _maxclustersize = 100;
	int _maxclusterdepth = 4;
	bool _bypass = false;

	bool _prenormalize = false;
	bool _selection_normalize = false;
	int _selection_criterion = 0;
	int _ncoreclusters = 1;
};

ostream& operator<<(ostream& stream, const MMFparams& params);

#endif
