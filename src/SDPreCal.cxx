// Copyright (C) 2014 Lucia Garbini <garbini@mpp.mpg.de>

// This is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


#include "SDPreCal.h"

#include <utility>


using namespace std;

SDPreCal::SDPreCal(std::vector<std::pair<double,double>> source_coll, std::vector<std::pair<double,double>> data_coll) :
	m_source_collection(source_coll),
	m_data_collection(data_coll){
		m_source_size = m_source_collection.size();
		m_data_size = m_data_collection.size();
	}

std::pair<Stats, Stats> SDPreCal::match( Line sline_a, Line sline_b, Line dline_a, Line dline_b, Stats prev_rx, Stats prev_ry ){
	
	double rx = ( dline_b.first - dline_a.first )/( sline_b.first - sline_a.first );
	double ry = ( dline_b.second / sline_b.second );
	prev_rx.add(rx);
	prev_ry.add(ry);
	
	return make_pair( prev_rx, prev_ry );
}

std::pair<Mapping, Stats> SDPreCal::genMap( int sline_i, int dline_i, Mapping prevMap, Stats prevStats ) {
	
	std::pair<Stats,Stats> new_Stats = match();//to be continued!!
	double distance_check = new_Stats.first.sigma()/new_Stats.first.mean();
	double intensity_check = new_Stats.second.sigma()/new_Stats.second.mean()
	std::pair<int, int> next_map = make_pair(sline_i, dline_i);
	
	if ( ( distance_check > distance_thr ||  intensity_check > intensity_thr ) || ( prevMap.first+2>=m_source_size && prevMap.second+2>=m_data_size ) ) {
		return (prevMap.push_back(next_map), new_Stats);
	}
	else {
		double best_err = 0;
		std::pair<Mapping, Stats> check_err_11 = genMap(prevMap.first+1, prevMap.second+1, prevMap, prevStats);
		std::pair<Mapping, Stats> check_err_21 = genMap(prevMap.first+2, prevMap.second+1, prevMap, prevStats);
		std::pair<Mapping, Stats> check_err_12 = genMap(prevMap.first+1, prevMap.second+2, prevMap, prevStats);
		
	}
}
void SDPreCal::calcPreCal(std::vector<std::pair<double, double>> source_lines, std::vector<std::pair<double, double>> data_lines) {	
};

SDPreCal::~SDPreCal() {
}
