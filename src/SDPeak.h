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


#ifndef RSPT_SDPEAK_H
#define RSPT_SDPEAK_H

#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>


namespace rspt {


class SDPeak {
public:
	///	@brief	Standard gauss function
	/// Function is normalized to integral of one.
	static double gauss(double x, double sigma);
	static double gauss(double x, double *p) { return gauss(x, p[0]); }

	///	@brief	Fits gauss + polynomial of degree 1
	/// @par    skew  should be > 0.2 * sigma for numerical stability
	/// Function is normalized to integral of one.
	static double skewedGauss(double x, double sigma, double skew);
	static double skewedGauss(double *x, double *p) { return skewedGauss(*x, p[0], p[1]); }

	///	@brief	Step function convolved with gaussian
	/// Converges to 1 for negative x and to 0 for positive x.
	static double stepWithSigma(double x, double sigma);
	static double stepWithSigma(const double *x, const double *p) { return gauss(*x, p[0]); }

	static double peakShape(double x, double center, double area, double sigma, double stepAmpl = 0, double skewFraction = 0, double skewWidth = 0);
	static double peakShape(const double *x, const double *p) { return peakShape(*x, p[0], p[1], p[2], p[3], p[4], p[5]); }

	SDPeak();
	virtual ~SDPeak();
};


class MultiPeakShape: public ROOT::Math::ParamFunctor {
protected:
	Int_t m_nPeaks;
	bool m_skewEnabled;
	TF1 *m_bg;

public:
	Int_t nPeaks() const { return m_nPeaks; }

	double operator()(double* x, double* p);

	TF1* newTF1(const char* name, TSpectrum *spectrum = 0);

	MultiPeakShape(Int_t n, bool enableSkew = true, TF1* bgModel = 0);
};


} // namespace rspt

#endif // RSPT_SDPEAK_H
