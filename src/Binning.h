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


#ifndef RSPT_BINNING_H
#define RSPT_BINNING_H

#include <cmath>
#include <stdint.h>

#include <TAxis.h>


namespace rspt {


class Binning {
protected:
	double m_from;
	double m_binWidth;
	size_t m_nBins;

public:
	double from() const { return m_from; }
	double to() const { return from() + (nBins() - 1) * binWidth(); }
	double until() const { return from() + nBins() * binWidth(); }

	double binWidth() const { return m_binWidth; }
	double nBins() const { return m_nBins; }
	
	double coord(int32_t bin) { return from() + binWidth() * bin; }
	int32_t bin(double x) { return int32_t( floor( (x + m_binWidth/2) / binWidth() ) ); }

	void run();

	Binning(): m_from(0), m_binWidth(0), m_nBins(0) {}
	Binning(double orig, double bwidth, size_t n = 0): m_from(orig), m_binWidth(bwidth), m_nBins(n) {}
	Binning(const TAxis *axis);
};



} // namespace rspt

#endif // RSPT_BINNING_H
