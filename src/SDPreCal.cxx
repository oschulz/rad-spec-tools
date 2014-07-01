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
// namespace rspt {
    
SDPreCal::SDPreCal() {
        m_prev_source=0;
        m_prev_data=0;
		m_source_size = m_source_collection.size();
		m_data_size = m_data_collection.size();
        m_dist_thres=0.10;
        m_int_thres=999999;
}

pair< SDPreCal::Stats, SDPreCal::Stats > SDPreCal::match(SDPreCal::Line sline_a, SDPreCal::Line sline_b, SDPreCal::Line dline_a, SDPreCal::Line dline_b, SDPreCal::Stats prev_rx, SDPreCal::Stats prev_ry)
{
    if( sline_b.first - sline_a.first!=0){
        double rx=(dline_b.first-dline_a.first)/(sline_b.first-sline_a.first);
        double ry = ( (dline_b.second/dline_a.second) - (sline_b.second/sline_a.second) );
//         std::cout<<std::endl;
        std::cout<<"rx: "<<rx<<std::endl;
        std::cout<<" prev rx_mean: "<<prev_rx.mean()<<std::endl;
        prev_rx.add(rx);
        
        std::cout<<"rx_mean: "<<prev_rx.mean()<<std::endl;
        
//         if( (ry< ( prev_ry.mean()+3*prev_ry.sigma() ) ) && ( ry> (prev_ry.mean()-3*prev_ry.sigma() ) ) ){
//             std::cout<<"added ry"<<std::endl;
//             prev_ry.add(ry);
//         }
    }else{
        std::cerr<<"!!!!!ERROR!!!!! line are equal: sline_a: "<<sline_a.first<<std::endl;
    }
    std::pair<Stats,Stats> output(prev_rx, prev_ry);
    return output;
}


std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> SDPreCal::genMap(std::vector<int> line_ind, SDPreCal::Mapping prevMap, SDPreCal::Stats_pair prevStats ) {
    std::vector<int>next_ind11;
    next_ind11.push_back(line_ind[2]);
    next_ind11.push_back(line_ind[3]);
    next_ind11.push_back(line_ind[2]+1);
    next_ind11.push_back(line_ind[3]+1);
    
    std::vector<int>next_ind12;
    next_ind12.push_back(line_ind[2]);
    next_ind12.push_back(line_ind[3]);
    next_ind12.push_back(line_ind[2]+1);
    next_ind12.push_back(line_ind[3]+2);
    std::vector<int>next_ind21;
    next_ind21.push_back(line_ind[2]);
    next_ind21.push_back(line_ind[3]);
    next_ind21.push_back(line_ind[2]+2);
    next_ind21.push_back(line_ind[3]+1);

    std::vector<int>next_ind01;
    next_ind01.push_back(line_ind[2]-1);
    next_ind01.push_back(line_ind[3]);
    next_ind01.push_back(line_ind[2]+0);
    next_ind01.push_back(line_ind[3]+1);

    std::vector<int>next_ind10;
    next_ind10.push_back(line_ind[2]);
    next_ind10.push_back(line_ind[3]-1);
    next_ind10.push_back(line_ind[2]+1);
    next_ind10.push_back(line_ind[3]+0);
    std::cout<<std::endl;
    std::cout<<"sline_a: "<<line_ind[0]<<"\tsline_b: "<<line_ind[2]<<"\tdline_a: "<<line_ind[1]<<"\tdline_b: "<<line_ind[3]<<std::endl;
    Stats_pair new_Stats = match(m_source_collection[line_ind[0]],
                                 m_source_collection[line_ind[2]],
                                 m_data_collection[line_ind[1]],
                                 m_data_collection[line_ind[3]],
                                prevStats.first,
                                prevStats.second);//to be continued!!
    double distance_check = new_Stats.first.sigma()/new_Stats.first.mean();
	double intensity_check = new_Stats.second.sigma()/new_Stats.second.mean();
    
    std::cout<<"new_Stats.first.sigma(): "<<new_Stats.first.sigma()<<"\tnew_Stats.first.mean(): "<<new_Stats.first.mean()<<std::endl;
    std::cout<<"distance_check: "<<distance_check<<std::endl;
    std::pair<int, int> next_map = make_pair(line_ind[2],line_ind[3]);
    if ( ( distance_check < m_dist_thres ||  intensity_check < m_int_thres )) {
        prevMap.push_back(next_map);
        std::cout<<"new map accepted"<<std::endl;
    }else{
        new_Stats=prevStats;
        if( next_ind01[3]<m_data_size){
            return genMap(next_ind01,prevMap,new_Stats);
        }else{
            next_ind21.clear();
            if( prevMap[prevMap.size()-1].second+2<m_source_size){
                next_ind21.push_back(prevMap[prevMap.size()-1].first+1);
                next_ind21.push_back(prevMap[prevMap.size()-1].second);
                next_ind21.push_back(prevMap[prevMap.size()-1].second+2);
                next_ind21.push_back(prevMap[prevMap.size()-1].second+1);
           
                return genMap(next_ind21,prevMap,new_Stats);
            }else{
                return make_pair(prevMap, new_Stats);
            }
        }
    }
    std::cout<<std::endl;
    if((next_ind11[2]>=m_source_size &&  next_ind11[3]>=m_data_size )){
//         std::cout <<std::endl;
//         std::cout<<"genMap end"<<std::endl;
//         std::cout<<"m_source_size: "<<m_source_size<<"\tm_data_size: "<<m_data_size<<std::endl;
//         std::cout<<"next_ind11: "<<next_ind11[2]<<","<<next_ind11[3]<<std::endl;
//         std::cout<<"next_ind12: "<<next_ind12[2]<<","<<next_ind12[3]<<std::endl;
//         std::cout<<"next_ind21: "<<next_ind21[2]<<","<<next_ind21[3]<<std::endl;
//         std::cout <<std::endl;
        
        return make_pair(prevMap, new_Stats);
    }else if(( next_ind21[2]>=m_source_size &&  next_ind12[3]>=m_data_size )){
        return  genMap(next_ind11, prevMap, prevStats);
    }else if(( next_ind21[2]<m_source_size &&  next_ind12[3]<m_data_size )){
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_11 = genMap(next_ind11, prevMap, new_Stats);
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_21 = genMap(next_ind21, prevMap, new_Stats);
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_12 = genMap(next_ind12, prevMap, new_Stats);
        double error_11=calcError(check_err_11.second.first,check_err_11.second.second);
        double error_12=calcError(check_err_12.second.first,check_err_12.second.second);
        double error_21=calcError(check_err_21.second.first,check_err_21.second.second);
        std::cout<<"error_11: "<<error_11<<"\terror_12: "<<error_12<<"\terror_21: "<<error_21<<std::endl;
        if(error_11<error_12){
            if(error_11<error_21){
                return check_err_11;
            }else{
                return check_err_21;
            }
        }else if(error_12<error_21){
            return check_err_12;
        }else{
            return check_err_21;
        }
    }else if(next_ind12[2]<m_source_size&& next_ind21[3]>=m_data_size){
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_11 = genMap(next_ind11, prevMap, new_Stats);
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_12 = genMap(next_ind12, prevMap, new_Stats);
        double error_11=calcError(check_err_11.second.first,check_err_11.second.second);
        double error_12=calcError(check_err_12.second.first,check_err_12.second.second);
        if(error_11<error_12){
            return check_err_11;
        }else{
            return check_err_12;
        }
        
    }else if( next_ind21[2]<m_source_size && next_ind12[3]>=m_data_size ){
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_11 = genMap(next_ind11, prevMap, new_Stats);
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> check_err_21 = genMap(next_ind21, prevMap, new_Stats);
        double error_11=calcError(check_err_11.second.first,check_err_11.second.second);
        double error_21=calcError(check_err_21.second.first,check_err_21.second.second);
        if(error_11<error_21){
            return check_err_11;
        }else{
            return check_err_21;
        }
    }else if( next_ind11[2]>=m_source_size|| next_ind21[2]>=m_source_size){
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair>check_err_01=genMap(next_ind01,prevMap,new_Stats);
        double error_00=calcError(new_Stats.first,new_Stats.second);
        double error_01=calcError(check_err_01.second.first,check_err_01.second.second);
        if(error_01<error_00){
            return check_err_01;
        }else{
            return make_pair(prevMap, new_Stats);
        }
    }else if( next_ind11[3]>=m_data_size||next_ind12[3]>=m_data_size){
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair>check_err_10=genMap(next_ind10,prevMap,new_Stats);
        double error_00=calcError(new_Stats.first,new_Stats.second);
        double error_10=calcError(check_err_10.second.first,check_err_10.second.second);
        if(error_10<error_00){
            return check_err_10;
        }else{
            return make_pair(prevMap, new_Stats);
        }
    }else{
        std::cout <<std::endl;
        std::cout <<"we are fucked up"<<std::endl;
        std::cout<<"m_source_size: "<<m_source_size<<"\tm_data_size: "<<m_data_size<<std::endl;
        std::cout<<"next_ind11: "<<next_ind11[2]<<","<<next_ind11[3]<<std::endl;
        std::cout<<"next_ind12: "<<next_ind12[2]<<","<<next_ind12[3]<<std::endl;
         std::cout<<"next_ind21: "<<next_ind21[2]<<","<<next_ind21[3]<<std::endl;
         std::cout <<std::endl;
        return make_pair(prevMap, new_Stats);
    }
	
}

TF1* SDPreCal::calcPreCal(vector< std::pair< double, double > > source_lines, vector< std::pair< double, double > > data_lines)
{
    m_source_collection=source_lines;
    m_data_collection=data_lines;
    m_source_size=m_source_collection.size();
    m_data_size=m_data_collection.size();
    Mapping start;
    Stats_pair stats;

    std::vector<int > start_ind;
    
    bool found_start_val=false;
    double start_fact_first=(source_lines[1].first-source_lines[0].first)/source_lines[1].first;
    double start_fact_second=(source_lines[2].first-source_lines[1].first)/source_lines[2].first;
    std::cout<<"start_fact_first: "<<start_fact_first<<std::endl;
    std::cout<<"start_fact_second: "<<start_fact_second<<std::endl;
    for(int i_first=0;i_first<data_lines.size()/2;++i_first){
        for(int j_first=i_first+1;j_first<data_lines.size()/2;++j_first){
            double data_fact_first=(data_lines[j_first].first-data_lines[i_first].first)/data_lines[j_first].first;
//             std::cout<<std::endl;
//             std::cout<<"data_fact: "<<data_fact_first<<std::endl;
//             std::cout<<"data_lines[i_first]: "<<data_lines[i_first].first<<"\tdata_lines[j_first]: "<<data_lines[j_first].first<<std::endl;
//             std::cout<<"i: "<<i_first<<"\tj_first: "<<j_first<<std::endl;
//             std::cout<<" ( 1-(data_fact_first/start_fact_first ) )"<< ( 1-(data_fact_first/start_fact_first ) )<<std::endl;
//             std::cout<<std::endl;
            double comp_first=abs(1-(data_fact_first/start_fact_first ));
            if(comp_first>0&&comp_first<0.05){
                for(int i_second=j_first;i_second<data_lines.size()/2;++i_second){
                    for(int j_second=j_first+1;j_second<data_lines.size()/2;++j_second){
                        double data_fact_second=(data_lines[j_second].first-data_lines[i_second].first)/data_lines[j_second].first;
//                         std::cout<<std::endl;
//                         std::cout<<"data_fact_second: "<<data_fact_second<<std::endl;
//                         std::cout<<"data_lines[i_second]: "<<data_lines[i_second].first<<"\tdata_lines[j_second]: "<<data_lines[j_second].first<<std::endl;
//                         std::cout<<"i: "<<i_second<<"\tj_first: "<<j_second<<std::endl;
//                         std::cout<<" ( 1-(data_fact_second/start_fact_second ) )"<< ( 1-(data_fact_second/start_fact_second) )<<std::endl;
//                         std::cout<<std::endl;
                        double comp_second=abs(1-(data_fact_second/start_fact_second ));
                        if(comp_second>0&&comp_second<0.20){
                            start.push_back(make_pair(0,i_first));
                            start.push_back(make_pair(1,i_second));
                            std::cout<<std::endl<<"first two matches..."<<std::endl;
                            stats.first.add((data_lines[j_first].first-data_lines[i_first].first)/(source_lines[1].first-source_lines[0].first));
                            stats.first.add((data_lines[j_second].first-data_lines[i_second].first)/(source_lines[2].first-source_lines[1].first));
                            std::cout<<"done\n"<<std::endl;
                            found_start_val=true;
                            start_ind.push_back(1);
                            start_ind.push_back(i_second);
                            start_ind.push_back(2);
                            start_ind.push_back(j_second);
                        break;
                        }
                    }
                    break;
                }

                break;
            }
        }
        if(found_start_val){
            break;
        }
    }
    
    
    if(found_start_val){
        std::cout<<"starting match: ch: "<<data_lines[start_ind[1]].first<<"\tenergy: "<<source_lines[start_ind[0]].first<<std::endl;
        std::pair<SDPreCal::Mapping, SDPreCal::Stats_pair> map_result=genMap(start_ind,start,stats);
        precal_graph=new TGraph();
        for(int i=0;i<map_result.first.size();++i){
            
            double channel=m_data_collection[map_result.first[i].second].first;
            double energy=m_source_collection[map_result.first[i].first].first;
            std::cout<<"channel: "<<channel<<"\tenergy: "<<energy<<std::endl;
            precal_graph->SetPoint(precal_graph->GetN(),channel,energy);
        }
        fit=new TF1("fit","pol1",0,8192);
        precal_graph->Fit("fit");
        return dynamic_cast<TF1*>(precal_graph->GetFunction("fit"));
    }else{
        std::cerr<<"does not found a starting point"<<std::endl;
    }
    return NULL;
}


SDPreCal::~SDPreCal() {


}
// } //namespace rspt
