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
#include <cmath>

#include<TAxis.h>

#include "rsptutils.h"
#include "SDCalibrator.h"

namespace rspt{
void transposePol1(TF1 **input) {
    char xTitle[100], yTitle[100];
    strcpy(xTitle,(*input)->GetXaxis()->GetTitle());
    strcpy(yTitle,(*input)->GetYaxis()->GetTitle());
    
    Double_t m, merr, c, cerr;
    c =(*input)->GetParameter(0)+1e-100;
    m = (*input)->GetParameter(1)+1e-100;
    cerr = (*input)->GetParError(0)+1e-100;
    merr = (*input)->GetParError(1)+1e-100;
    (*input)->SetParameter(0,-1.*c/m);
    (*input)->SetParameter(1,1./m);
    (*input)->SetParError(0,c/m*(sqrt(pow((merr/m),2) + pow((cerr/c),2))));
    (*input)->SetParError(1,pow(merr/m,2));
    (*input)->GetXaxis()->SetTitle(yTitle);
    (*input)->GetYaxis()->SetTitle(xTitle);
}

/// @brief Checks if the tested mean is compatible with precalibration
int desiredPeak(int iter,int fitted_lines, std::vector< double > energy, SDFitData *fit, TF1 *cal_ch2e) {
    int peakdesired=0;
    double fit_residual=1000;
    for(unsigned int j=iter; j<iter+fitted_lines; j++) {
        peakdesired=0;
        for(int h=1; h<=fit->getNPeaks(); h++) {
            if (TMath::Abs(energy[j] - cal_ch2e->Eval(fit->getMean(h))) < fit_residual) {
                peakdesired=h;
                fit_residual=TMath::Abs(energy[j] - cal_ch2e->Eval(fit->getMean(h)));
            }
        }
        
        fit->setUsage(peakdesired);
        fit->setEnergy(peakdesired,energy[j]);
        fit_residual=1000;
    }
    return 0;
}

TF1* rescalFCh2Fe(const TF1* rescal_ch2fch, const TF1* cal_ch2e)
{
    TF1 *rescalfct_e2fe = new TF1("rescalfct_e2fe",SQRTQuadFunct,cal_ch2e->Eval(rescal_ch2fch->GetXmin()),cal_ch2e->Eval(rescal_ch2fch->GetXmax()),3);
    /*
     *   double rescalpar[] = {
     *       sqrt(pow(calfct_ch2e->GetParameter(1)*rescalfct_ch2fch->GetParameter(0),2)-pow(rescalfct_ch2fch->GetParameter(1),2)*calfct_ch2e->GetParameter(1)*calfct_ch2e->GetParameter(0)+pow(rescalfct_ch2fch->GetParameter(2)*calfct_ch2e->GetParameter(0),2)),
     *       sqrt(pow(rescalfct_ch2fch->GetParameter(1),2)*calfct_ch2e->GetParameter(1)-2*pow(rescalfct_ch2fch->GetParameter(2),2)*calfct_ch2e->GetParameter(0)),
     *       rescalfct_ch2fch->GetParameter(2)
     };*/
    float p0 = rescal_ch2fch->GetParameter(0);
    float p1 = rescal_ch2fch->GetParameter(1);
    float p2 = rescal_ch2fch->GetParameter(2);
    
    float p0_err = rescal_ch2fch->GetParError(0);
    float p1_err = rescal_ch2fch->GetParError(1);
    float p2_err = rescal_ch2fch->GetParError(2);
    
    float c0 = cal_ch2e->GetParameter(0);
    float c1 = cal_ch2e->GetParameter(1);
    
    float c0_err = cal_ch2e->GetParError(0);
    float c1_err = cal_ch2e->GetParError(1);
    
    // dragon parametrization for conversion of FWHM(Ch) into FWHM(E)
    double rescalpar[3] = {
        sqrt( pow(c1*p0,2) - c0*c1*pow(p1,2) + pow(c0*p2,2) ),
        sqrt( c1*pow(p1,2) - 2*c0*pow(p2,2) ),
        p2
    };
    
    // results of gaussian error propagation
    double par_resE_err[3] = {
        sqrt( pow(c1*pow(p1,2)+2*c0*pow(p2,2),2)*pow(c0_err,2) + pow(2*c1*pow(p0,2)+c0*pow(p1,2),2)*pow(c1,2) + pow(2*c1*pow(p0,2)*p0_err,2) + pow(2*c0*c1*p1*p1_err,2) + pow(2*pow(c0,2)*p2*p2_err,2) ),
        sqrt( pow(2*pow(p2,2)*c0_err,2) + pow(pow(p1,2)*c1_err,2) + pow(2*c1*p1*p1_err,2) + (4*c0*p2*p2_err) ),
        p2_err
    };
    
    rescalfct_e2fe->SetParameter(0,rescalpar[0]);
    rescalfct_e2fe->SetParError(0,par_resE_err[0]);
    rescalfct_e2fe->SetParameter(1,rescalpar[1]);
    rescalfct_e2fe->SetParError(1,par_resE_err[1]);
    rescalfct_e2fe->SetParameter(2,rescalpar[2]);
    rescalfct_e2fe->SetParError(2,par_resE_err[2]);
    rescalfct_e2fe->SetName("rescal_e2fe");
    rescalfct_e2fe->SetTitle("Resolution calibration FWHM_{E}(E)");
    return rescalfct_e2fe;
     }

}// namespace rspt