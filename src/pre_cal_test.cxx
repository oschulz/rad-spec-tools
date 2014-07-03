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
#include <stdlib.h>
#include<string>
#include<stdexcept>

#include<TFile.h>
#include<TH1.h>
#include<TF1.h>
#include<TSpectrum.h>
#include "TSystem.h"
#include <TROOT.h>

#include"SDPreCal.h"
#include"SDMultiLineFitter.h"
#include"SDCalibrator.h"
#include"SDFitData.h"


bool compare_pair(float first,float second){
    return(first<second);
}
int main(int argc, char **argv){

    std::string histname="raw_hist";
    std::string inputname;
    std::string rad_source="th";
    float threshold=0.01;
    float sigma=1;
    int opt = 0;
    while ((opt = getopt(argc, argv, "?n:h:r:t:s:")) != -1) {
        switch (opt) {
            case '?': {; return 0; }
            case 'n': {inputname=optarg ; break;}
            case 'h': {histname=optarg;break;}
            case 'r': {rad_source=optarg;break;}
            case 't': {threshold=atof(optarg);break;}
            case 's': {sigma=atof(optarg);break;}
            default:{std::string error="Unkown command line option ";error+=opt; throw std::invalid_argument(error.c_str());}
        }
    }
    gSystem->Load("libTree");
    gROOT->ProcessLine("#include <vector>");
    
    
    TFile *input=new TFile(inputname.c_str(),"READONLY");
    TH1* raw_hist=dynamic_cast<TH1*>(input->FindObjectAny(histname.c_str()));
    TSpectrum *spec=new TSpectrum(100);
    double width;
    if(raw_hist!=0){
        rspt::SDPreCal precal;
        if(rad_source=="gs"){
             raw_hist->GetXaxis()->SetRangeUser(7000,60000);
        }
        int npeaks=spec->Search(raw_hist,sigma,"",threshold);
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
        std::vector<double> energy;
        std::vector<double> ig;
        std::vector< double > line_info;
        if(rad_source=="th232"){
            double energy_iso1_spec[] = {209.253, 238.632, 240.986, 300.087, 328.000, 338.320, 463.004, 583.191, 727.330, 860.564, 911.204, 964.766, 968.971, 2614.533}; // here: Th-232
            int nlines=14;
            width=0.02;
            for(int i=0;i<nlines;i++){
                line_info.push_back(energy_iso1_spec[i]);
            }
        }else if(rad_source=="th228"){
            double energy_iso1_spec[] = {238.632, 510.77, 583.191, 727.330, 860.564,1593, 2614.533}; // here: Th-232
            int nlines=7;
            width=0.02;
            std::cout<<"nlines: "<<nlines<<std::endl;
            for(int i=0;i<nlines;i++){
                line_info.push_back(energy_iso1_spec[i]);
            }
        }else if(rad_source=="gs"){
           
            precal.setDistThres(0.05);
            double energy_iso1_spec[] = {609.3,778.9,964.08,1085.9,1112.1,1408,1764.5,2614.533};
            int nlines=8;
            width=0.08;
            for(int i=0;i<nlines;i++){
                line_info.push_back(energy_iso1_spec[i]);
            }
        }
//         double Ig_iso[]={0.0389,0.43,0.041,0.0328,0.0295,0.1127,0.0440,0.845,0.0658,0.1242,0.258,0.0499,0.158,0.99};
// 
        


        TF1* precalibration=precal.calcPreCal(line_info,peakinfo);
        if(precalibration!=NULL){
            std::cout<<"precal: "<<precalibration->GetParameter(1)<<"x+"<<precalibration->GetParameter(0)<<std::endl;
            TFile *output=new TFile("test_calibration.root","recreate");
            rspt::SDMultiLineFitter mfitter;
            rspt::SDCalibrator cal;
            mfitter.setPreCal(precalibration);
            mfitter.setWidth(width);
            mfitter.setSigma(sigma);
            mfitter.setThreshold(threshold);
            std::vector<rspt::SDFitData*> fits=mfitter.makeCalFits(raw_hist,line_info);
            std::cout<<"fits[0] npeaks: "<<fits[0]->getNPeaks()<<std::endl;
            for(int fit=0;fit<fits.size();++fit){
                cal.addResult(fits[fit]);
            }
            cal.calibrate();
            TList* results=cal.getCalObjects();
            results->Write();
            raw_hist->Write();
            output->Write();
            fits.clear();
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