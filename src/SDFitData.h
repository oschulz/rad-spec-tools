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


#ifndef RSPT_SDFITDATA_H
#define RSPT_SDFITDATA_H

#include <TF1.h>


namespace rspt{


class SDFitData {
public:
	SDFitData(TF1 *fit, int npeaks);
	virtual ~SDFitData();

	double getMean(size_t index);
	double getMeanError(size_t index);
	double getSigma(size_t index);
	double getSigmaError(size_t index);

	bool getUsage(size_t index);
	bool setUsage(size_t index, bool use=true);

	int getNPeaks() {return m_npeaks;}

	double getEnergy(size_t index);
	void setEnergy(size_t index, double energy);

	bool getResUsage(size_t index);
	bool setResUsage(size_t index, bool use=true);

protected:
	bool m_valid;

	int m_npeaks;

	TF1 *m_fit;

	std::vector<bool> m_usage;
	std::vector<bool> m_res_usage;
	std::vector<double> m_energy;
};


} // namespace rspt

#endif // RSPT_SDFITDATA_H
