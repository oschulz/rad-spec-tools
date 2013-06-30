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


} // namespace rspt
