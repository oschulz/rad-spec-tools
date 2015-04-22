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


#include "SDPeak.h"

#include <stdexcept>
#include <limits>
#include <cmath>
#include <cassert>


using namespace std;


namespace rspt {


double SDPeak::gauss(double x, double sigma)
	{ return (1. / (sigma * sqrt(2 * M_PI)) ) * exp( - pow(x / (sqrt(2) * sigma), 2.) ); }


double SDPeak::skewedGauss(double x, double sigma, double skew)
	{ return ( exp( x/skew + sigma*sigma/(2.*skew*skew) ) * erfc( x / (sqrt(2)*sigma) + sigma / (sqrt(2)*skew)) ) / (2.*skew); }


double SDPeak::stepWithSigma(double x, double sigma)
	{ return 0.5 * erfc( x / (sqrt(2) * sigma) ); }


double SDPeak::peakShape(double x, double center, double area, double sigma, double stepAmpl, double skewFraction, double skewWidth) {
	double x_c = x - center;
	double skew = skewWidth * center;
	return area * (1 - skewFraction) * gauss(x_c, sigma)
	     + (skewFraction != 0 ? area * skewFraction * skewedGauss(x_c, sigma, skew) : 0)
		 + stepAmpl * stepWithSigma(x_c, sigma);
}


SDPeak::SDPeak() {
}


SDPeak::~SDPeak() {
}



double MultiPeakShape::operator()(double *x, double* p) {
	double result = 0;
	if (m_bg != 0) {
		result = (*m_bg)(x, p);
		p += m_bg->GetNpar();
	}

	double skewFraction = 0, skewWidth = 0;
	if (m_skewEnabled) { skewFraction = *p++; skewWidth = *p++; }

	for (int i = 0; i < m_nPeaks; ++i) {
		double center = *p++, area = *p++, sigma = *p++, stepAmpl = *p++;
		result += SDPeak::peakShape(*x, center, area, sigma, stepAmpl, skewFraction, skewWidth);
	}
	return result;
}


TF1* MultiPeakShape::newTF1(const char* name, TSpectrum *spectrum, double sigma) {
	if (spectrum == 0) throw std::invalid_argument("Input spectrum cannot be a null pointer");
	if (m_nPeaks != spectrum->GetNPeaks()) throw std::invalid_argument("Number of peaks to fit doesn't match number of peaks in input spektrum");

	// double maxPos = numeric_limits<double>::max();
	// double maxPos = 1e10;

	static const int nPeakPar = 4;
	static const int nSkewPar = m_skewEnabled ? 2 : 0;

	Int_t nBgPar = (m_bg != 0) ? m_bg->GetNpar() : 0;
	TF1 *tf = new TF1(name, *this, 0, 16.0 + 16.0 * m_nPeaks , nBgPar + nSkewPar + m_nPeaks * nPeakPar);
	tf->SetNpx(1000);

	for (Int_t i = 0; i < nBgPar; ++i) {
		tf->SetParName(i, TString::Format("bg_%s", m_bg->GetParName(i)));
		tf->SetParameter(i, m_bg->GetParameter(i));
		double a, b;
		m_bg->GetParLimits(i, a, b);
		tf->SetParLimits(i, a, b);
	}

	Int_t p = nBgPar;

	if (m_skewEnabled) {
		tf->SetParName(p, "peak_skewFraction");
		tf->SetParLimits(p, 0, 1);
		tf->SetParameter(p, 0.2);
		++p;
		tf->SetParName(p, "peak_skewWidth");
		tf->SetParLimits(p, 1e-4, 1e-1);
		tf->SetParameter(p, 1e-3);
		++p;
	}

	for (Int_t i = 0; i < m_nPeaks; ++i) {
		double center = spectrum->GetPositionX()[i];
		double height = spectrum->GetPositionY()[i];
		double area = height / SDPeak::gauss(0, sigma);
		double areaLimit = 1000 * area;

		tf->SetParName(p, TString::Format("peak%i_center",i+1));
		tf->SetParameter(p, center);
		++p;
		tf->SetParName(p, TString::Format("peak%i_area",i+1));
		if (areaLimit > 0) tf->SetParLimits(p, 0, areaLimit);
		tf->SetParameter(p, area);
		++p;
		tf->SetParName(p, TString::Format("peak%i_sigma",i+1));
		tf->SetParLimits(p, 0, 10000);
		tf->SetParameter(p, sigma);
		++p;
		tf->SetParName(p, TString::Format("peak%i_stepAmpl",i+1));
		tf->SetParLimits(p, 0, height);
		tf->SetParameter(p, height/double(100));
		++p;
	}
	return tf;
}


MultiPeakShape::MultiPeakShape(Int_t n, bool enableSkew, TF1* bgModel)
	: m_nPeaks(n), m_skewEnabled(enableSkew), m_bg(bgModel)
{
	if (m_nPeaks < 0) throw invalid_argument("Number of peaks must be >= 0");
}


} // namespace rspt
