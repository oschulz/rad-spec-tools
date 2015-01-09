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


#ifndef RSPT_SDMULTILINEFITTER_H
#define RSPT_SDMULTILINEFITTER_H

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "SDCalibrator.h"
#include "SDFitData.h"


namespace rspt{


class SDMultiLineFitter {
public:
	SDMultiLineFitter();
	virtual ~SDMultiLineFitter() {}

	void setThreshold(double thresh);
	void setSigma(float sig);
	void setTSpecSigma(double sig);
	void setPreCal(double slope, double intercept);
	void setPreCal(TF1 *precal_ch2e);
	void setWidth(double width) { m_width = width; }
	void setRange(double lowEdge, double highEdge);

	void resetPreCal();

	std::vector<rspt::SDFitData*> makeCalFits(TH1* raw_hist, std::vector<double> energy,
                                              double s_factor=0.0099,const char* opt="Q+" , std::vector<bool> *reject_res_cal=0);

	double m_maxADCch;

protected:
    double m_tspec_sigma;
	float m_sigma;
	bool debug;
	double m_threshold;
	double m_low_limit;
	double m_high_limit;
	double m_width;

	SDCalibrator calibrator;

	TCanvas *m_cal_canv;
	TF1 *m_preCalibration_ch2e;
	TF1 *m_preCalibration_e2ch;

	void init();

	std::pair<double, int> getRange(std::vector<double> energy, int iter, int lines_to_fit);
};


} // namespace rspt

#endif // RSPT_SDMULTILINEFITTER_H
