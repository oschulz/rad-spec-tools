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
	

	double frac_source=(sline_b-sline_a)/sline_b;
	double frac_data=(dline_b-dline_a)/dline_b;
	double matching=frac_source/frac_data;
	
	next.stats.add(matching);
	
	return next.stats;
}


SDPreCal::next_line_info SDPreCal::genLineInfo(next_line_info prev,int next_s,int next_d) {
	next_line_info out;
	if(next_s==0){
		out.s_ind_a=prev.s_ind_a;
		out.s_ind_b=prev.s_ind_b;
	}else{
		out.s_ind_a=prev.s_ind_b;
		out.s_ind_b=out.s_ind_a+next_s;
	}
	if(next_d==0){
		out.d_ind_a=prev.d_ind_a;
		out.d_ind_b=prev.d_ind_b;
	}else{
		out.d_ind_a=prev.d_ind_b;
		out.d_ind_b=out.d_ind_a+next_d;
	}
	out.stats=prev.stats;
	return out;
}

pair< SDPreCal::Mapping, SDPreCal::Stats > SDPreCal::genMap(SDPreCal::next_line_info next, SDPreCal::Mapping prevMap) {

	std::pair<SDPreCal::Mapping, SDPreCal::Stats> return_val=make_pair(prevMap,next.stats);
	
	
	next.stats=match(next);
	
	double prev_err;
	
	prev_err=calcError(next.stats);
	if(prev_err>0.01&&next.stats.sigma()/next.stats.mean()>0.05&&next.stats.n()>1){
		return return_val;
	}
	prevMap.push_back(make_pair(next.s_ind_b,next.d_ind_b));
	return_val=make_pair(prevMap,next.stats);
	
	std::pair<SDPreCal::Mapping, SDPreCal::Stats> next_mapping;
	
	for(int next_s_line=next.s_ind_b;next_s_line<m_source_size-1;++next_s_line){
		for(int next_d_line=next.d_ind_b;next_d_line<m_data_size-1;++next_d_line){
			next_mapping=genMap(genLineInfo(next,next_s_line-next.s_ind_b+1,next_d_line-next.d_ind_b+1),prevMap);
			double new_err=calcError(next_mapping.second);
			if(new_err<prev_err){
				prev_err=new_err;
				return_val=next_mapping;
			}
		}
	}
	return return_val;
}


TGraph* SDPreCal::calcPreCal(std::vector< double > source_lines, std::vector< double > data_lines){
	Stats new_Stats;

		m_source_collection=source_lines;
		m_data_collection=data_lines;
		m_source_size=m_source_collection.size();
		m_data_size=m_data_collection.size();

		std::cout<<"no. source lines: "<<m_source_size
				<<"\nno. data lines: "<<m_data_size
				<<std::endl;
	
	std::pair<SDPreCal::Mapping, SDPreCal::Stats> best_mapping;
	std::pair<SDPreCal::Mapping, SDPreCal::Stats> curr_mapping;
	Mapping startMap;

	next_line_info next;
	int n_start_vals=0;

	for(int s_line_start=0;s_line_start<m_source_size;++s_line_start){
		for(int d_line_start=0;d_line_start<m_data_size;++d_line_start){
			
			for(int s_line_next=s_line_start+1;s_line_next<m_source_size;++s_line_next){
				for(int d_line_next=d_line_start+1;d_line_next<m_data_size;++d_line_next){
					startMap.clear();
					next.stats.clear();
					next.s_ind_a=s_line_start;
					next.d_ind_a=d_line_start;
					next.s_ind_b=s_line_next;
					next.d_ind_b=d_line_next;
					next.stats=match(next);

					if(calcError(next.stats)<0.05){
						n_start_vals++;
						startMap.push_back(make_pair(s_line_start,d_line_start));
						curr_mapping=genMap(next,startMap);
						std::cout<<"map ready"<<std::endl;
						if( ( calcError(curr_mapping.second)<calcError(best_mapping.second) ) || best_mapping.second.n()==0){
							best_mapping=curr_mapping;
						}
					}
				}
			}
		}
	}
	std::cout<<"n_start_vals: "<<n_start_vals<<std::endl;
	std::cout<<"best match: "<<calcError(best_mapping.second)<<std::endl;
	std::cout<<"rel. Error: "<<best_mapping.second.sigma()/best_mapping.second.mean()<<std::endl;
	precal_graph=new TGraph();
	precal_graph->SetNameTitle("precal_graph","precal_graph");
	
	for(int i=0;i<best_mapping.first.size();++i){
		double channel=m_data_collection[best_mapping.first[i].second];
		double energy=m_source_collection[best_mapping.first[i].first];
		if ( debug ) std::cerr<<"channel: "<<channel<<"\tenergy: "<<energy<<std::endl;
		precal_graph->SetPoint(precal_graph->GetN(),channel,energy);
	}
	std::cout<<"best_mapping.first.size(): "<<best_mapping.first.size()<<std::endl;
		fit=new TF1("fit","pol1", m_low_mca, m_high_mca);
		precal_graph->Fit("fit");
		
		return dynamic_cast<TGraph*>(precal_graph);
}



} // namespace rspt
