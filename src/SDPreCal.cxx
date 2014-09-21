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

#include <iostream>
#include <utility>

#include <TGraph.h>

using namespace std;


namespace rspt {


SDPreCal::SDPreCal( double l_mca, double h_mca ) {
	m_low_mca = l_mca;
	m_high_mca = h_mca;
	debug = true;
	m_source_size = m_source_collection.size();
	m_data_size = m_data_collection.size();
	m_dist_thres = 0.1;
}


SDPreCal::~SDPreCal() {
}


SDPreCal::Stats SDPreCal::match(next_line_info next) {
	double sline_a;
	double sline_b;
	double dline_a;
	double dline_b;

	sline_a = m_source_collection[next.s_ind_a];
	sline_b = m_source_collection[next.s_ind_b];
	dline_a = m_data_collection[next.d_ind_a];
	dline_b = m_data_collection[next.d_ind_b];
	double rx = (dline_b - dline_a) / (sline_b - sline_a);

	if ( (sline_b - sline_a != 0) && (dline_b-dline_a != 0) ) {
		double rx = (dline_b - dline_a) / (sline_b - sline_a);
		if (debug) {
			cerr << "dline_b: " << dline_b << "\tdline_a: " << dline_a << "\tsline_b: " << sline_b << "\tsline_a: " << sline_a << endl;
			cerr << "rx: " << rx << endl;
			cerr << " prev rx_mean: " << next.stats.mean() << endl;
		}

		next.stats.add(rx);
		if (debug) cerr << "rx_mean: " << next.stats.mean() << endl;
	} else {
		if (debug) cerr << "!!!!!ERROR!!!!! line are equal: sline_a: " << next.s_ind_a<<endl;
	}

	Stats out=next.stats;
	return out;
}


SDPreCal::next_line_info SDPreCal::genLineInfo(next_line_info prev,int next_s,int next_d) {
	next_line_info out;
	if (next_s == 0) out.s_ind_a = prev.s_ind_a;
	else out.s_ind_a = prev.s_ind_b;

	if (next_d == 0) out.d_ind_a = prev.d_ind_a;
	else out.d_ind_a = prev.d_ind_b;

	out.s_ind_b = out.s_ind_a + next_s;
	out.d_ind_b = out.d_ind_a + next_d;
	out.stats = prev.stats;
	return out;
}


std::pair<SDPreCal::Mapping, SDPreCal::Stats> SDPreCal::genMap(next_line_info next, SDPreCal::Mapping prevMap, int calls) {
	if (debug) cerr<<endl;
	Stats new_Stats = match(next);
	double distance_check = new_Stats.sigma() / new_Stats.mean();

	if (debug) {
		cerr << "s_line_a: " << next.s_ind_a << "\tsline_b: " << next.s_ind_b << "\tdline_a: " << next.d_ind_a << "\tdline_b: " << next.d_ind_b << endl;
		cerr << "new_Stats.first.sigma(): " << new_Stats.sigma() << "\tnew_Stats.first.mean(): " << new_Stats.mean() << endl;
		cerr << "distance_check: " << distance_check << endl;
	}

	pair<int, int> next_map = make_pair(next.s_ind_b, next.d_ind_b);

	if (distance_check < m_dist_thres) {
		prevMap.push_back(next_map);
		next.stats = new_Stats;
		if (debug) cerr << "new map accepted" << endl;
	}

	if (debug) cerr<<endl;

	next_line_info next11 = genLineInfo(next,1,1);
	next_line_info next12 = genLineInfo(next,1,2);
	next_line_info next21 = genLineInfo(next,2,1);

	if (next11.s_ind_b >= m_source_size && next11.d_ind_b >= m_data_size) {
		if (debug) cerr << "return" << endl;
		return make_pair(prevMap, new_Stats);
	} else if ( next21.s_ind_b >= m_source_size && next12.d_ind_b >= m_data_size
		&& next11.s_ind_b < m_source_size && next11.d_ind_b < m_data_size)
	{
		if(debug) cerr << "calling 11" << endl;
		return  genMap(next11, prevMap,++calls);
	} else if (next21.s_ind_b < m_source_size && next12.d_ind_b < m_data_size) {
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_11 = genMap(next11, prevMap, ++calls);
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_21 = genMap(next21, prevMap, ++calls);
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_12 = genMap(next12, prevMap, ++calls);

		double error_11=calcError(check_err_11.second);
		double error_12=calcError(check_err_12.second);
		double error_21=calcError(check_err_21.second);

		if(debug){
			cerr<<"calling 11"<<endl;
			cerr<<"calling 12"<<endl;
			cerr<<"calling 21"<<endl;
		}

		if (error_11 < error_12){
			if (error_11<error_21) return check_err_11;
			else return check_err_21;
		} else if (error_12 < error_21) {
			return check_err_12;
		} else {
			return check_err_21;
		}
	} else if (next12.d_ind_b < m_data_size && next21.s_ind_b >= m_source_size && next11.s_ind_b < m_source_size) {
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_11 = genMap(next11, prevMap, ++calls);
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_12 = genMap(next12, prevMap, ++calls);
		double error_11 = calcError(check_err_11.second);
		double error_12 = calcError(check_err_12.second);

		if (debug) {
			cerr << "calling 11" << endl;
			cerr << "calling 12" << endl;
		}

		if (error_11<error_12) return check_err_11;
		else return check_err_12;
	} else if (next21.s_ind_b < m_source_size && next12.d_ind_b >= m_data_size &&  next11.d_ind_b < m_data_size) {
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_11 = genMap(next11, prevMap, ++calls);
		pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_21 = genMap(next21, prevMap, ++calls);
		double error_11 = calcError(check_err_11.second);
		double error_21 = calcError(check_err_21.second);

		if (debug) {
			cerr << "calling 11" << endl;
			cerr << "calling 21" << endl;
		}

		if (error_11<error_21) return check_err_11;
		else return check_err_21;
	} else {
		return make_pair(prevMap, new_Stats);
	}
}

bool SDPreCal::compare_pair(std::pair<double,std::pair<int,int> > a, std::pair<double,std::pair<int,int> > b){
    return( a.first < b.first );
}

TF1* SDPreCal::calcPreCal(std::vector<double> source_lines, std::vector<double> data_lines) {
	m_source_collection=source_lines;
	m_data_collection=data_lines;
	m_source_size=m_source_collection.size();
	m_data_size=m_data_collection.size();

	Mapping start;
	Stats stats;

	next_line_info start_ind;

	bool found_start_val = false;
	double start_fact_first = (source_lines[1] - source_lines[0]) / source_lines[1];
	double start_fact_second = (source_lines[2] - source_lines[1]) / source_lines[2];
	cerr << "start_fact_first: " << start_fact_first << endl;
	cerr << "start_fact_second: " << start_fact_second << endl;

	vector<pair<double,pair<int,int> > > comp_firstPair_all;
	vector<pair<double,pair<int,int> > > comp_secondPair_all;
	int i_first_data_line;
	int j_first_data_line;
	int i_second_data_line;
	int j_second_data_line;

	for (int i_first = 0; i_first < data_lines.size(); ++i_first){
		if (data_lines.size() > (i_first+1)) for (int j_first = i_first + 1; j_first < data_lines.size(); ++j_first) {
			double data_fact_first = (data_lines[j_first] - data_lines[i_first]) / data_lines[j_first];
			double comp_first = abs(1 - (data_fact_first/start_fact_first));
			pair<int,int> data_ind_firstPair;
			data_ind_firstPair = make_pair( i_first, j_first );
			pair<double, pair<int,int>> comp_firstPair;
			comp_firstPair = make_pair (comp_first, data_ind_firstPair);
			comp_firstPair_all.push_back(comp_firstPair);
			if (debug) {
				cerr << endl;
				cerr << "i_first: " << i_first << "\tj_first: " << j_first << endl;
				cerr << "data_fact_first: " << data_fact_first << endl;
				cerr << "data_lines[" << i_first << "]: " << data_lines[i_first] << "\tdata_lines[" << j_first << "]: " << data_lines[j_first] << endl;
				cerr << "data_fact_first/start_fact_first: " << data_fact_first/start_fact_first << endl;
				cerr << "comp_first: " << comp_first << endl;
				cerr << endl;
			}
		}
		else {
			cerr <<"Number of peaks found with TSpectrum not enough to proceed with Precalibration"<<endl;
			continue;
		}
		sort( comp_firstPair_all.begin(), comp_firstPair_all.end(), compare_pair);
		if (debug) cerr<<"best comp_first = "<<comp_firstPair_all[0].first<<" \t i_first = "<<comp_firstPair_all[0].second.first<<" , j_first = "<<comp_firstPair_all[0].second.second<<endl;
		if ( (comp_firstPair_all[0].first>0) && (comp_firstPair_all[0].first<0.05) ) {
			i_first_data_line = comp_firstPair_all[0].second.first;
			j_first_data_line =  comp_firstPair_all[0].second.second;
			for (int i_second = j_first_data_line; i_second<data_lines.size(); ++i_second) {
				if (data_lines.size() > (j_first_data_line+1)) for (int j_second = j_first_data_line + 1; j_second<data_lines.size(); ++j_second) {
					double data_fact_second = (data_lines[j_second] - data_lines[i_second]) / data_lines[j_second];
					double comp_second = abs(1-(data_fact_second/start_fact_second ));
					pair<int,int> data_ind_secondPair;
					data_ind_secondPair = make_pair( i_second, j_second );
					pair<double, pair<int,int>> comp_secondPair;
					comp_secondPair = make_pair (comp_second, data_ind_secondPair);
					comp_secondPair_all.push_back(comp_secondPair);
					if (debug) {
						cerr << endl;
						cerr << "i_second: " << i_second << "\tj_second: " << j_second << endl;
						cerr << "data_fact_second: " << data_fact_second << endl;
						cerr << "data_lines[" << i_second << "]: " << data_lines[i_second] << "\tdata_lines[" << j_second << "]: " << data_lines[j_second] << endl;
						cerr << "data_fact_second/start_fact_second "<< data_fact_second/start_fact_second<<endl;
						cerr << "comp_second: " << comp_second << endl;
						cerr << endl;
					}
				}
				else {
					cerr <<"Number of peaks found with TSpectrum not enough to proceed with Precalibration"<<endl;
					continue;
				}
				sort( comp_secondPair_all.begin(), comp_secondPair_all.end(), compare_pair);
				if (debug) cerr<<"best comp_second = "<<comp_secondPair_all[0].first<<" \t i_second = "<<comp_secondPair_all[0].second.first<<" , j_first = "<<comp_secondPair_all[0].second.first<<endl;
				if ( comp_secondPair_all[0].first>0 && comp_secondPair_all[0].first<0.05 ) {
					i_second_data_line = comp_secondPair_all[0].second.first;
					j_second_data_line =  comp_secondPair_all[0].second.second;
					start.push_back(make_pair(0, i_first_data_line));
					start.push_back(make_pair(1, i_second_data_line));
					cerr << endl << "first two matches..." << endl;
					stats.add( (data_lines[j_first_data_line] - data_lines[i_first_data_line]) / (source_lines[1] - source_lines[0]) );
					stats.add( (data_lines[j_second_data_line] - data_lines[i_second_data_line]) / (source_lines[2] - source_lines[1]) );

					cerr << "done\n" << endl;

					found_start_val = true;

					start_ind.d_ind_a = i_second_data_line;
					start_ind.s_ind_a = 1;
					start_ind.d_ind_b = j_second_data_line;
					start_ind.s_ind_b = 2;
					start_ind.stats = stats;

					break;
				}
				else continue;
			}
			break;
		}
		else continue;

		if (found_start_val) break;
	}


	if (found_start_val) {
		cerr << "starting match: ch: " << data_lines[start_ind.d_ind_a] << "\tenergy: " << source_lines[start_ind.s_ind_a] << endl;
		pair<SDPreCal::Mapping, SDPreCal::Stats> map_result = genMap(start_ind, start, 1);
		precal_graph = new TGraph();

		for (int i = 0; i < map_result.first.size(); ++i) {
			double channel = m_data_collection[map_result.first[i].second];
			double energy = m_source_collection[map_result.first[i].first];
			if (debug) cerr << "channel: " << channel << "\tenergy: " << energy << endl;
			precal_graph->SetPoint(precal_graph->GetN(), channel, energy);
		}

		fit = new TF1("fit","pol1", m_low_mca, m_high_mca);
		precal_graph->Fit("fit");
		return dynamic_cast<TF1*>(precal_graph->GetFunction("fit"));
	} else {
		cerr << "does not find a starting point" << endl;
	}

	return 0;
}


} // namespace rspt
