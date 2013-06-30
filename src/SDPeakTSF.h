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


#ifndef RSPT_SDPEAKTSF_H
#define RSPT_SDPEAKTSF_H

#include <Rtypes.h>


namespace rspt {


class SDPeakTSF {
public:
	// Peak shape
	//
	// pos: peak position
	// sigma: sigma of peak
	// ampl: amplitude of gaussian
	// t: relative amplitude of skewed gaussian
	// b: exponential skewed gaussian, relative to sigma
	// step: amplitude of step function

	static double shape(double x, double pos, double ampl, double sigma, double t, double b, double step);

	// Peak area

	static double area(double ampl, double sigma, double t, double b);

	SDPeakTSF();
	virtual ~SDPeakTSF();
};


} // namespace rspt


#ifdef __CINT__
#pragma link C++ class rspt::SDPeakTSF+;
#endif

#endif // RSPT_SDPEAKTSF_H
