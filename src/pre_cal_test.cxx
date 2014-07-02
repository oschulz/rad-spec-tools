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
#include <iostream>
#include<string>
#include<stdexcept>
#include <stdlib.h>
#include"SDPreCal.h"


#include<TFile.h>
#include<TH1.h>
#include<TF1.h>
#include<TSpectrum.h>
#include "TSystem.h"
#include <TROOT.h>
bool compare_pair(float first,float second){
    return(first<second);
}
int main(int argc, char **argv){

    std::string histname="raw_hist";
    std::string inputname;
    
    int opt = 0;
    while ((opt = getopt(argc, argv, "?n:h:")) != -1) {
        switch (opt) {
            case '?': {; return 0; }
            case 'n': {inputname=optarg ; break;}
            case 'h': {histname=optarg;break;}
            default:{std::string error="Unkown command line option ";error+=opt; throw std::invalid_argument(error.c_str());}
        }
    }
    gSystem->Load("libTree");
    gROOT->ProcessLine("#include <vector>");
    
    
    TFile *input=new TFile(inputname.c_str(),"READONLY");
    TH1* raw_hist=dynamic_cast<TH1*>(input->FindObjectAny(histname.c_str()));
    TSpectrum *spec=new TSpectrum(100);
    if(raw_hist!=0){
        rspt::SDPreCal precal;
//         raw_hist->GetXaxis()->SetRangeUser(500,8196);
        int npeaks=spec->Search(raw_hist,1,"",0.01);
        float* x_pos=spec->GetPositionX();
        float* y_pos=spec->GetPositionY();
        std::vector< double > peakinfo;
        for(int i_peak=0;i_peak<npeaks;++i_peak){
            peakinfo.push_back(x_pos[i_peak]);
        }
        
        std::sort(peakinfo.begin(),peakinfo.end(),compare_pair);
        for(int i_peak=0;i_peak<npeaks;++i_peak){
            std::cout<<"index: "<<i_peak<<"\tpos: "<<peakinfo[i_peak]<<std::endl;
        }
    
        double energy_iso1_spec[] = {209.253, 238.632, 240.986, 300.087, 328.000, 338.320, 463.004, 583.191, 727.330, 860.564, 911.204, 964.766, 968.971, 2614.533}; // here: Th-232
//         double Ig_iso[]={0.0389,0.43,0.041,0.0328,0.0295,0.1127,0.0440,0.845,0.0658,0.1242,0.258,0.0499,0.158,0.99};
//         double energy_iso1_spec[] = {58,1086,1112,1408,2614.5};
        int nlines=14;
        std::vector<double> energy;
        std::vector<double> ig;
        std::vector< double > line_info;
        for(int i=0;i<nlines;i++){
            line_info.push_back(energy_iso1_spec[i]);
        }
        TF1* precalibration=precal.calcPreCal(line_info,peakinfo);
        if(precalibration!=NULL){
            std::cout<<"precal: "<<precalibration->GetParameter(1)<<"x+"<<precalibration->GetParameter(0)<<std::endl;
        }else{
            std::cout<<"precalibration failed"<<std::endl;
        }
       
    }else{
        std::cerr<<"Histogram "<<histname<<" was not found in TFile"<<std::endl;
    }
//     delete spec;
    delete input;
    std::cout<<"end of pre_cal_test"<<std::endl;
    
    return 1;
}