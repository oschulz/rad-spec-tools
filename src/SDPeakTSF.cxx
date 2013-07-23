// Copyright (C) 2013 Oliver Schulz <oliver.schulz@tu-dortmund.de>

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


#include "SDPeakTSF.h"

#include <stdexcept>
#include <cmath>


using namespace std;


namespace rspt {


// Method contains code taken class TSpectrumFit of the CERN ROOT framework,
// version 5.34.09. Code originally written by Miroslav Morhac and repackaged
// for ROOT by Rene Brun.

double SDPeakTSF::shape(double x, double pos, double ampl, double sigma, double t, double b, double step) {  
	double p = 0;
	if (sigma > 0.0001) p = (x - pos) / sigma;
	else {
		if (x == pos) p = 0;
		else p = 10;
	}
	
	double gaussY = 0;
	if (fabs(p) < 3) {
		if ((p * p) < 700) gaussY = exp(-p * p);
		else gaussY = 0;
	}
	
	double skewedGaussY = 0;
	if (t != 0) {
		double c = p + 1. / (2. * b);
		double e = p / b;
		if (e > 700) e = 700;
		skewedGaussY = t * exp(e) * erfc(c) / 2.;
	}

	double stepY = 0;
	if (step != 0) stepY = step * erfc(p) / 2.;

	return ampl * (gaussY + skewedGaussY + stepY);
}


// Method contains code taken class TSpectrumFit of the CERN ROOT framework,
// version 5.34.09. Code originally written by Miroslav Morhac and repackaged
// for ROOT by Rene Brun.

double SDPeakTSF::area(double ampl, double sigma, double t, double b) {  
   double odm_pi = 1.7724538;

   double result = 0;
   if (b != 0) result = 0.5 / b;
   result = - result * result;
   if (fabs(result) < 700) result = ampl * sigma * (odm_pi + t * b * exp(result));
   else result = ampl * sigma * odm_pi;

   return result;
}


SDPeakTSF::SDPeakTSF() {
}


SDPeakTSF::~SDPeakTSF() {
}



double MultiPeakShapeTSF::operator()(double *x, double* p) {
	double result = 0;

	// Quadradic background:
	double bg0 = *p++, bg1 = *p++, bg2 = *p++;
	result += bg0 + (*x) * bg1 + (*x)*(*x) * bg2;

	// Peak sigma and tail step:
	double sigma = *p++, step = *p++;

	// Optional skewed-gaussian tail parameters t and b:
	double t = 0, b = 0;
	if (m_skewEnabled) { t = *p++; b = *p++; }

	for (int i = 0; i < m_nPeaks; ++i) {
		double pos = *p++, ampl = *p++;
		result += SDPeakTSF::shape(*x, pos, ampl, sigma, t, b, step);
	}
	return result;
}


TF1* MultiPeakShapeTSF::newTF1(const char* name, const TAxis *xAxis, TSpectrumFit *tsf) const {
	// double maxPos = numeric_limits<double>::max();
	// double maxPos = 1e10;

	Int_t fromBin = xAxis->GetFirst();
	Int_t toBin = xAxis->GetLast();
	Double_t from = xAxis->GetBinLowEdge(fromBin);
	Double_t to = xAxis->GetBinLowEdge(toBin);
	Double_t scale = (from - to) / (fromBin - toBin);
	
	Double_t bg0 = 0, bg0Err = 0, bg1 = 0, bg1Err = 0, bg2 = 0, bg2Err = 0;
	tsf->GetBackgroundParameters(bg0, bg0Err, bg1, bg1Err, bg2, bg2Err);

	Double_t sigma = 0, sigmaErr = 0;
	tsf->GetSigma(sigma, sigmaErr);

	Double_t t = 0, tErr = 0, b = 0, bErr= 0, step = 0, stepErr = 0;
	tsf->GetTailParameters(t, tErr, b, bErr, step, stepErr);

	Int_t nBgPar = 3;
	Int_t nPeakCommonPar = m_skewEnabled ? 4 : 2;
	Int_t nPeakPar = 2;

	TF1 *tf = new TF1(name, *this, from, to, nBgPar + nPeakCommonPar + m_nPeaks * nPeakPar);
	tf->SetNpx(10000);

	Int_t p = 0;

	Double_t fromSq = from * from;
	Double_t scaleSq = scale * scale;

	//!! TODO: Apply from and scale to bg parms
	tf->SetParName(p, "bg_0"); tf->SetParameter(p, bg0 - bg1 * from / scale + bg2 * fromSq / scaleSq); ++p;
	tf->SetParName(p, "bg_1"); tf->SetParameter(p, bg1 / scale -  2 * bg2 * from / scaleSq); ++p;
	tf->SetParName(p, "bg_2"); tf->SetParameter(p, bg2 / scaleSq); ++p;

	tf->SetParName(p, "sigma"); tf->SetParameter(p, scale * sigma); ++p;
	tf->SetParName(p, "step"); tf->SetParameter(p, step); ++p;
	
	if (m_skewEnabled) {
		tf->SetParName(p, "t"); tf->SetParameter(p, t); ++p;
		tf->SetParName(p, "b"); tf->SetParameter(p, b); ++p;
	}

	for (Int_t i = 0; i < m_nPeaks; ++i) {
		tf->SetParName(p, TString::Format("peak%i_pos",i+1));
		tf->SetParameter(p, from + scale * tsf->GetPositions()[i]);
		++p;
		tf->SetParName(p, TString::Format("peak%i_ampl",i+1));
		tf->SetParameter(p, tsf->GetAmplitudes()[i]);
		++p;
	}

	return tf;
}


MultiPeakShapeTSF::MultiPeakShapeTSF(Int_t n, bool enableSkew)
	: m_nPeaks(n), m_skewEnabled(enableSkew)
{
	if (m_nPeaks < 0) throw invalid_argument("Number of peaks must be >= 0");
}


} // namespace rspt
