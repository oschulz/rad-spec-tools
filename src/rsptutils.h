// Copyright (C) 2013 Thomas Quante <thomas.quante@tu-dortmund.de>
//               2014 Lucia Garbini <garbini@mpp.mpg.de>

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


#ifndef RSPT_RSPTUTILS_H
#define UTILS_H

#include <TF1.h>

#include "SDFitData.h"


namespace rspt{


void transposePol1(TF1 **input);

/// @brief Checks if the tested mean is compatible with precalibration
int desiredPeak(int iter, int fitted_lines, std::vector< double > energy, SDFitData *fit, TF1 *cal_ch2e) ;

TF1* rescalFCh2Fe(const TF1* rescal_ch2fch, const TF1* cal_ch2e);


} // namespace rspt

#endif // RSPT_RSPTUTILS_H
