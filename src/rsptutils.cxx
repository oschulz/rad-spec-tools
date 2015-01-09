// Copyright (C) 2013 Thomas Quante <thomas.quante@tu-dortmund.de>
//               2014 Lucia Garbini <garbini@mpp.mpg.de>
//               2014 Oliver Schulz <oschulz@mpp.mpg.de>

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


#include "rsptutils.h"

#include<iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include <TAxis.h>

using namespace std;

namespace rspt{

double SQRTQuadFunct(double *x, double *par) {
	return sqrt(pow(par[0],2) + pow(par[1],2) * x[0] + pow(par[2],2) * pow(x[0],2));
}


void transposePol1(TF1 *hist) {
	string xTitle(hist->GetXaxis()->GetTitle());
	string yTitle(hist->GetYaxis()->GetTitle());

	Double_t m, mErr, c, cErr;
	c = hist->GetParameter(0) + 1e-100;
	m = hist->GetParameter(1) + 1e-100;
	cErr = hist->GetParError(0) + 1e-100;
	mErr = hist->GetParError(1) + 1e-100;

	hist->SetParameter(0, -1.*c/m);
	hist->SetParameter(1, 1./m);
	hist->SetParError(0, c/m * sqrt(pow( (mErr/m), 2) + pow((cErr/c), 2)) );
	hist->SetParError(1, pow(mErr/m, 2));
	hist->GetXaxis()->SetTitle(yTitle.c_str());
	hist->GetYaxis()->SetTitle(xTitle.c_str());
}

bool compare_pair(std::pair<double,int> a, std::pair<double,int> b){
    return(a.first < b.first);
}

int desiredPeak(int iter, int fitted_lines, std::vector<double> energy, SDFitData *fit, TF1 *cal_ch2e) {
	int peakdesired=0;
	double fit_residual=9999999999;

	for (unsigned int j = iter; j < iter+fitted_lines; j++) {
		cerr << "energy[" << j << "] = " << energy[j] << endl;
		std::vector<std::pair<double, int> > peak_des_all;

		for (int h = 1; h <= fit->getNPeaks(); h++) {
			std::pair<double,int> peak_des;
			peak_des = make_pair(TMath::Abs(energy[j] - cal_ch2e->Eval(fit->getMean(h))), h);
			cerr << "h = "<<h<<", ADC= " << fit->getMean(h) << "\t energy= " << cal_ch2e->Eval(fit->getMean(h)) << "\t residual = " << TMath::Abs(energy[j] - cal_ch2e->Eval(fit->getMean(h))) << endl;
			peak_des_all.push_back(peak_des);
		}

		std::sort(peak_des_all.begin(), peak_des_all.end(), compare_pair);
		if ( peak_des_all[0].first < fit_residual ) {
			peakdesired = peak_des_all[0].second;
			cerr<<"desired peak = "<<peakdesired<<endl;
		}
		else continue;

		fit->setUsage(peakdesired);
		fit->setEnergy(peakdesired,energy[j]);
	}
	return 0;
}


TF1* rescalFCh2Fe(const TF1* rescal_ch2fch, const TF1* cal_ch2e) {
	TF1 *rescalfct_e2fe = new TF1("rescalfct_e2fe", SQRTQuadFunct, cal_ch2e->Eval(rescal_ch2fch->GetXmin()),
		cal_ch2e->Eval(rescal_ch2fch->GetXmax()),3);

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
		sqrt( pow(c1 * p0, 2) - c0*c1*pow(p1, 2) + pow(c0 * p2, 2) ),
		sqrt( c1*pow(p1, 2) - 2*c0*pow(p2, 2) ),
		p2
	};

	// results of gaussian error propagation
	double par_resE_err[3] = {
		sqrt( pow(c1 * pow(p1, 2) + 2 * c0 * pow(p2, 2), 2) * pow(c0_err, 2)
			+ pow(2 * c1 * pow(p0, 2) + c0 * pow(p1, 2), 2) * pow(c1, 2)
			+ pow(2 * c1 * pow(p0, 2) * p0_err, 2)
			+ pow(2 * c0 * c1 * p1 * p1_err, 2)
			+ pow(2 * pow(c0, 2) * p2 * p2_err, 2)
		),
		sqrt( pow( 2 * pow(p2, 2) * c0_err, 2)
			+ pow(pow(p1, 2) * c1_err, 2)
			+ pow( 2 * c1 * p1 * p1_err, 2)
			+ (4 * c0 * p2 * p2_err)
		),
		p2_err
	};

	rescalfct_e2fe->SetParameter(0, rescalpar[0]);
	rescalfct_e2fe->SetParError(0, par_resE_err[0]);
	rescalfct_e2fe->SetParameter(1, rescalpar[1]);
	rescalfct_e2fe->SetParError(1, par_resE_err[1]);
	rescalfct_e2fe->SetParameter(2, rescalpar[2]);
	rescalfct_e2fe->SetParError(2, par_resE_err[2]);
	rescalfct_e2fe->SetName("rescal_e2fe");
	rescalfct_e2fe->SetTitle("Resolution calibration FWHM_{E}(E)");
	return rescalfct_e2fe;
}


} // namespace rspt
