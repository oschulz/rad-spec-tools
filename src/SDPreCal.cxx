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
#include<iostream>
#include <utility>
#include <TGraph.h>


using namespace std;
namespace rspt {

SDPreCal::SDPreCal( double l_mca, double h_mca ) {

	m_low_mca = l_mca;
	m_high_mca = h_mca;
	debug=true;
	m_source_size = m_source_collection.size();
	m_data_size = m_data_collection.size();
	m_dist_thres=0.1;

}

SDPreCal::Stats SDPreCal::match(next_line_info next)
{
	double sline_a;
	double sline_b;
	double dline_a;
	double dline_b;

	sline_a=m_source_collection[next.s_ind_a];
	sline_b=m_source_collection[next.s_ind_b];
	dline_a=m_data_collection[next.d_ind_a];
	dline_b=m_data_collection[next.d_ind_b];
// 	double rx=(dline_b-dline_a)/(sline_b-sline_a);


	
	double frac_source=(sline_b-sline_a)/sline_b;
	double frac_data=(dline_b-dline_a)/dline_b;
	double matching=frac_source/frac_data;
	
		next.stats.add(matching);
// 		if(debug){
// 			std::cerr<<"dline_b: "<<dline_b<<"\tdline_a: "<<dline_a<<"\tsline_b: "<<sline_b<<"\tsline_a: "<<sline_a<<std::endl;
// 			std::cerr<<" prev rx_mean: "<<next.stats.mean()<<std::endl;
// 		}

// 		if(debug){
// 			std::cerr<<"rx_mean: "<<next.stats.mean()<<"\trx_sigma: "<<next.stats.sigma()<<"\trx_n: "<<next.stats.n()<<std::endl;
// 		}


	
	return next.stats;
}

SDPreCal::next_line_info SDPreCal::genLineInfo(next_line_info prev,int next_s,int next_d)
{
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

std::pair<SDPreCal::Mapping, SDPreCal::Stats> SDPreCal::genMap(next_line_info next, SDPreCal::Mapping prevMap) {

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
	
// 	if(next.d_ind_b+1< m_data_size){
// 		next_mapping=genMap(genLineInfo(next,0,1),prevMap);
// 		if(calcError(next_mapping.second)<prev_err){
// 			return_val=next_mapping;
// 		}
// 	}
// 	if(next.s_ind_b+1< m_source_size){
// 
// 		next_mapping=genMap(genLineInfo(next,1,0),prevMap);
// 		if(calcError(next_mapping.second)<prev_err){
// 			return_val=next_mapping;
// 		}
// 	}


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

// std::pair<SDPreCal::Mapping, SDPreCal::Stats> SDPreCal::genMap(next_line_info next, SDPreCal::Mapping prevMap) {
// 
// 	if(debug) std::cerr<<std::endl;
// 	Stats new_Stats = match(next);
// 	double distance_check = new_Stats.sigma()/new_Stats.mean();
// 
// 	if(debug){
// 		std::cerr<<"s_line_a: "<<next.s_ind_a<<"\tsline_b: "<<next.s_ind_b<<"\tdline_a: "<<next.d_ind_a<<"\tdline_b: "<<next.d_ind_b<<std::endl;
// 		std::cerr<<"new_Stats.first.sigma(): "<<new_Stats.sigma()<<"\tnew_Stats.first.mean(): "<<new_Stats.mean()<<std::endl;
// 		std::cerr<<"distance_check: "<<distance_check<<std::endl;
// 	}
// 	std::pair<int, int> next_map = make_pair(next.s_ind_b,next.d_ind_b);
// 	if ( ( distance_check < m_dist_thres  )) {
// 		prevMap.push_back(next_map);
// 		next.stats=new_Stats;
// 		if(debug) std::cerr<<"new map accepted"<<std::endl;
// 	}
// 
// 	if(debug) std::cerr<<std::endl;
// 
// 	next_line_info next11=genLineInfo(next,1,1);
// 	next_line_info next12=genLineInfo(next,1,2);
// 	next_line_info next21=genLineInfo(next,2,1);
// 
// 	if((next11.s_ind_b>=m_source_size &&  next11.d_ind_b>=m_data_size )){
// 		if(debug) std::cerr<<"return"<<std::endl;
// 		return make_pair(prevMap, new_Stats);
// 	}
// 	else if(( next21.s_ind_b>=m_source_size &&   next12.d_ind_b>=m_data_size &&next11.s_ind_b<m_source_size &&next11.d_ind_b<m_data_size)){
// 		if(debug) std::cerr<<"calling 11"<<std::endl;
// 		return  genMap(next11, prevMap,++calls);
// 	}
// 	else if(( next21.s_ind_b<m_source_size &&  next12.d_ind_b<m_data_size )){
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_11 = genMap(next11, prevMap,++calls);
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_21 = genMap(next21, prevMap,++calls);
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_12 = genMap(next12, prevMap,++calls);
// 
// 		double error_11=calcError(check_err_11.second);
// 		double error_12=calcError(check_err_12.second);
// 		double error_21=calcError(check_err_21.second);
// 		if(debug){
// 			std::cerr<<"calling 11"<<std::endl;
// 			std::cerr<<"calling 12"<<std::endl;
// 			std::cerr<<"calling 21"<<std::endl;
// 		}
// 		if(error_11<error_12){
// 			if(error_11<error_21){
// 				return check_err_11;
// 			}
// 			else{
// 				return check_err_21;
// 			}
// 		}
// 		else if(error_12<error_21){
// 			return check_err_12;
// 		}
// 		else{
// 			return check_err_21;
// 		}
// 
// 	}else if(next12.d_ind_b<m_data_size&& next21.s_ind_b>=m_source_size&&  next11.s_ind_b<m_source_size){
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_11 = genMap(next11, prevMap,++calls);
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_12 = genMap(next12, prevMap,++calls);
// 		double error_11=calcError(check_err_11.second);
// 		double error_12=calcError(check_err_12.second);
// 		if(debug){
// 			std::cerr<<"calling 11"<<std::endl;
// 			std::cerr<<"calling 12"<<std::endl;
// 		}
// 		if(error_11<error_12){
// 			return check_err_11;
// 		}
// 		else{
// 			return check_err_12;
// 		}
// 	}
// 	else if( next21.s_ind_b<m_source_size && next12.d_ind_b>=m_data_size &&  next11.d_ind_b<m_data_size){
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_11 = genMap(next11, prevMap,++calls);
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> check_err_21 = genMap(next21, prevMap,++calls);
// 		double error_11=calcError(check_err_11.second);
// 		double error_21=calcError(check_err_21.second);
// 
// 		if(debug){
// 			std::cerr<<"calling 11"<<std::endl;
// 			std::cerr<<"calling 21"<<std::endl;
// 		}
// 		if(error_11<error_21){
// 			return check_err_11;
// 		}
// 		else{
// 			return check_err_21;
// 		}
// 	}
// 	else{
// 		return make_pair(prevMap, new_Stats);
// 	}
// 
// }

TGraph* SDPreCal::calcPreCal(std::vector< double > source_lines, std::vector< double > data_lines)
{
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
						std::cout<<"found next start val("<<n_start_vals<<")"<<std::endl;
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


// TF1* SDPreCal::calcPreCal(vector< double > source_lines, vector<  double > data_lines)
// {
// 	m_source_collection=source_lines;
// 	m_data_collection=data_lines;
// 	m_source_size=m_source_collection.size();
// 	m_data_size=m_data_collection.size();
// 	Mapping start;
// 	Stats stats;
// 
// 	next_line_info start_ind;
// 
// 	bool found_start_val=false;
// 	double start_fact_first=(source_lines[1]-source_lines[0])/source_lines[1];
// 	double start_fact_second=(source_lines[2]-source_lines[1])/source_lines[2];
// 	std::cerr<<"start_fact_first: "<<start_fact_first<<std::endl;
// 	std::cerr<<"start_fact_second: "<<start_fact_second<<std::endl;
// 
// 	for( int i_first=0; i_first<data_lines.size(); ++i_first ) {
// 		for( int j_first=i_first+1; j_first<data_lines.size(); ++j_first ){
// 			double data_fact_first = (data_lines[j_first] - data_lines[i_first] )/ data_lines[j_first];
// 			double comp_first=abs(1-(data_fact_first/start_fact_first ));
// 			if (debug) {
// 				std::cout<<std::endl;
// 				std::cout<<"i_first: "<<i_first<<"\tj_first: "<<j_first<<std::endl;
// 				std::cout<<"data_fact_first: "<<data_fact_first<<std::endl;
// 				std::cout<<"data_lines["<<i_first<<"]: "<<data_lines[i_first]<<"\tdata_lines["<<j_first<<"]: "<<data_lines[j_first]<<std::endl;
// 				std::cout<<"data_fact_first/start_fact_first: "<< data_fact_first/start_fact_first <<std::endl;
// 				std::cout<<"comp_first: "<< comp_first<<std::endl;
// 				std::cout<<std::endl;
// 			}
// 			if( comp_first>0 && comp_first<0.05 ){
// 				std::cout<<"i_second = "<<j_first<<std::endl;
// 				for( int i_second=j_first; i_second<data_lines.size(); ++i_second ){
// 					for( int j_second=j_first+1; j_second<data_lines.size(); ++j_second ){
// 						double data_fact_second=(data_lines[j_second]-data_lines[i_second])/data_lines[j_second];
// 						double comp_second=abs(1-(data_fact_second/start_fact_second ));
// 						if (debug) {
// 							std::cout<<std::endl;
// 							std::cout<<"i_second: "<<i_second<<"\tj_second: "<<j_second<<std::endl;
// 							std::cout<<"data_fact_second: "<<data_fact_second<<std::endl;
// 							std::cout<<"data_lines["<<i_second<<"]: "<<data_lines[i_second]<<"\tdata_lines["<<j_second<<"]: "<<data_lines[j_second]<<std::endl;
// 							std::cout<<"data_fact_second/start_fact_second"<< data_fact_second/start_fact_second<<std::endl;
// 							std::cout<<"comp_second: "<< comp_second<<std::endl;
// 							std::cout<<std::endl;
// 						}
// 						if(comp_second>0&&comp_second<0.20){
// 							start.push_back(make_pair(0,i_first));
// 							start.push_back(make_pair(1,i_second));
// 							std::cerr<<std::endl<<"first two matches..."<<std::endl;
// 							stats.add((data_lines[j_first]-data_lines[i_first])/(source_lines[1]-source_lines[0]));
// 							stats.add((data_lines[j_second]-data_lines[i_second])/(source_lines[2]-source_lines[1]));
// 							std::cerr<<"done\n"<<std::endl;
// 							found_start_val=true;
// 
// 							start_ind.d_ind_a=i_second;
// 							start_ind.s_ind_a=1;
// 							start_ind.d_ind_b=j_second;
// 							start_ind.s_ind_b=2;
// 							start_ind.stats=stats;
// 							break;
// 						}	
// 					}
// 					break;
// 				}
// 				break;
// 			}
// 		}
// 		if(found_start_val){
// 			break;
// 		}// 						startMap.push_back(make_pair(s_line_start,d_line_start));
// 						curr_mapping=genMap(next,startMap);
// 						std::cout<<"map ready"<<std::endl;
// 						if( ( calcError(curr_mapping.second)<calcError(best_mapping.second) ) || best_mapping.second.n()==0){
// 							best_mapping=curr_mapping;
// 	}	
// 
//     
// 	if(found_start_val){
// 		std::cerr<<"starting match: ch: "<<data_lines[start_ind.d_ind_a]<<"\tenergy: "<<source_lines[start_ind.s_ind_a]<<std::endl;
// 		std::pair<SDPreCal::Mapping, SDPreCal::Stats> map_result=genMap(start_ind,start,1);
// 		precal_graph=new TGraph();
// 		for(int i=0;i<map_result.first.size();++i){
// 
// 			double channel=m_data_collection[map_result.first[i].second];
// 			double energy=m_source_collection[map_result.first[i].first];
// 			if ( debug ) std::cerr<<"channel: "<<channel<<"\tenergy: "<<energy<<std::endl;
// 			precal_graph->SetPoint(precal_graph->GetN(),channel,energy);
// 		}
// 
// 		fit=new TF1("fit","pol1", m_low_mca, m_high_mca);
// 		precal_graph->Fit("fit");
// 		return dynamic_cast<TF1*>(precal_graph->GetFunction("fit"));
// 	}
// 	else{
// 		std::cerr<<"does not find a starting point"<<std::endl;
// 	}
// 	return NULL;
// }


SDPreCal::~SDPreCal() {


}
} //namespace rspt
