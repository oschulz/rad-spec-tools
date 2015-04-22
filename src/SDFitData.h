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
	SDFitData(TF1 *fit, size_t npeaks);
	virtual ~SDFitData();

	int getNPeaks() { return m_npeaks; }

	double getSkewFraction() const ;
	double getSkewFractionError() const;

	double getSkewWidth() const;
	double getSkewWidthError() const;

	double getStepAmpl(size_t index) const;
	double getStepAmplError(size_t index) const;

	double getPeakArea(size_t index) const;
	double getPeakAreaError(size_t index) const;

	double getMean(size_t index) const;
	double getMeanError(size_t index) const;

	double getSigma(size_t index) const;
	double getSigmaError(size_t index) const;

	bool getUsage(size_t index) const;
	void setUsage(size_t index, bool use=true);

	bool getResUsage(size_t index) const;
	void setResUsage(size_t index, bool use=true);

	double getEnergy(size_t index) const;
	void setEnergy(size_t index, double energy);

protected:
	bool m_valid;

	int m_npeaks;

	TF1 *m_fit;

	std::vector<bool> m_usage;
	std::vector<bool> m_res_usage;
	std::vector<double> m_energy;

	size_t getParNumber(const char* nameForm, size_t index) const;
	double getParValue(const char* nameForm, size_t index) const;
	double getParError(const char* nameForm, size_t index) const;
};


} // namespace rspt

#endif // RSPT_SDFITDATA_H
