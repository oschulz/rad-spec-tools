/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2014  <copyright holder> <email>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#include<TMath.h>

#include "SDCalibrator.h"
#include "rsptutils.h"
namespace rspt{

double SQRTQuadFunct(double *x, double *par) {
    double y=0.0;
    y=sqrt( pow(par[0],2) + pow(par[1],2)*x[0] + pow(par[2],2)*pow(x[0],2) );
    return y;
}
    
SDCalibrator::SDCalibrator()
{
    init();
}

SDCalibrator::~SDCalibrator(){

}

void SDCalibrator::init()
{


    if(rescal_graph==0||cal_graph==0){
        setupCalGraphs();
    }
    m_objects=new TList();

    rescal_ch2fch=NULL;
    rescal_e2fe=NULL;
    cal_ch2e=NULL;
    cal_e2ch=NULL;
    cal_graph=NULL;
    rescal_graph=NULL;
    calCanv=NULL;
}

void SDCalibrator::setupCalGraphs()
{
    cal_graph=new TGraphErrors();
    cal_graph->SetNameTitle("cal_graph","calibration graph");
    rescal_graph=new TGraphErrors();
    rescal_graph->SetNameTitle("rescal_graph","resolution graph");
    m_objects->Add(cal_graph);
    m_objects->Add(rescal_graph);
    
}

void SDCalibrator::addResult(SDFitData* data)
{
    if(cal_graph==0&&rescal_graph==0) {
        setupCalGraphs();
    }
        
    int nAdds=0;
    for(unsigned int i=1; i<=data->getNPeaks(); i++)
    {
        if(data->getUsage(i))
        {
            nAdds++;
            // add new calibration point to old graph
            cal_graph->SetPoint(cal_graph->GetN(),data->getEnergy(i),data->getMean(i));
            cal_graph->SetPointError(cal_graph->GetN()-1,0,data->getMeanError(i));
            // add new resolution point to old graph
            if(data->getResUsage(i)){
                rescal_graph->SetPoint(rescal_graph->GetN(),data->getMean(i),TMath::Sqrt(8*TMath::Log(2))*data->getSigma(i));
                rescal_graph->SetPointError(rescal_graph->GetN()-1,data->getMeanError(i),TMath::Sqrt(8*TMath::Log(2))*data->getSigmaError(i));
            }
        }
    }
}

int SDCalibrator::calibrate()
{
    if (cal_graph->GetN() >= 2)
    {
        bool point_removed=false;
        cal_e2ch = new TF1("cal_e2ch","pol1",0,3000); // channel(energy)
        cal_e2ch->SetTitle("Calibration Ch(E)");
        cal_e2ch->SetLineColor(3);
        cal_e2ch->SetLineWidth(1);
        m_objects->Add(cal_e2ch);
        // Fitting resolution with option  "W": Set all weights to 1 for non empty bins; ignore error bars;
        // R: Use the Range specified in the function range
        //Q: Quiet
        cal_graph->Fit("cal_e2ch","WRQ");
        cal_graph->GetXaxis()->SetTitle("Energy");
        cal_graph->GetYaxis()->SetTitle("Channels");
        cal_graph->GetFunction("cal_e2ch")->ResetBit(512);
        
//         for(int point=0;point<cal_graph->GetN();++point){
//             double x;
//             double y;
//             if(cal_graph->GetPoint(point,x,y)){
//                 if(abs(y-cal_e2ch->Eval(x))>0.01*y){
//                     cal_graph->RemovePoint(point);
//                     std::cout <<"remove point x="<<x<<",y=<<"<<y<<std::endl<<std::endl<<std::endl;
//                     point_removed=true;
//                     point--;
//                 }
//             }
//         }
//         if(point_removed==true){
//             return calibrate();
//         }
        // Copy "cal_fit" as "calEqn" (as calibration equation) and invert it (to x=channel, y=energy)
        cal_ch2e = (TF1*)cal_e2ch->Clone("cal_ch2e");
        rspt::transposePol1(&cal_ch2e);
        m_objects->Add(cal_ch2e);
        cal_ch2e->SetRange(0,8192);
        cal_ch2e->SetNameTitle("cal_ch2e","Calibration E(Ch)");
        cal_ch2e->GetXaxis()->SetTitle("Channels");
        cal_ch2e->GetYaxis()->SetTitle("Energy");

        std::cout<<"Calibration function: "<<cal_ch2e->GetParameter(1)<<"*x + "<<cal_ch2e->GetParameter(0)<<std::endl;
    }
    else std::cout << "Less than 2 point in calibration graph. Skipping fit for calibration equation.\n";


    if (rescal_graph->GetN() >= 2)
    {
        TF1 *rescal_ch2fch_lin=new TF1("rescal_ch2fch_lin","pol1",1,10000);

        // Fitting resolution with option  "W": Set all weights to 1 for non empty bins; ignore error bars;
        // because the lower sigmas have a very much smaller error, and therefore the important high
        // energy lines have nearly no influence on the fit!
        // R: Use the Range specified in the function range
        // Q: Quiet

        rescal_graph->SetTitle("Resolution calibration FWHM_{Ch}(Ch)");
        rescal_graph->GetXaxis()->SetTitle("Channels");
        rescal_graph->GetYaxis()->SetTitle("FWHM_{Ch} / channels");
        rescal_graph->Fit("rescal_ch2fch_lin","0Q");

        rescal_ch2fch=new TF1("rescal_ch2fch",SQRTQuadFunct,1,10000,3);
        rescal_ch2fch->SetTitle("Resolution calibration FWHM_{Ch}(Ch)");
        rescal_ch2fch->SetLineColor(4);
        rescal_ch2fch->SetLineWidth(1);
        m_objects->Add(rescal_ch2fch);

        rescal_ch2fch->SetParameter(0,rescal_ch2fch_lin->GetParameter(0));
        rescal_ch2fch->SetParameter(1,0);
        rescal_ch2fch->SetParameter(2,rescal_ch2fch_lin->GetParameter(1));
        delete rescal_ch2fch_lin;


        std::cout << "Fit results for calibration equation:" << std::endl;
        rescal_graph->Fit("rescal_ch2fch");

        rescal_graph->GetFunction("rescal_ch2fch")->ResetBit(512);


        rescal_e2fe=rspt::rescalFCh2Fe(rescal_ch2fch,cal_ch2e);
        m_objects->Add(rescal_e2fe);

        rescal_e2fe->SetLineColor(3);
        rescal_e2fe->SetLineWidth(1);
    }
    else std::cout << "Less than 2 point in resolution graph. Skipping fit for resolution equation.\n";
    return 1;
}



} //namespace rspt