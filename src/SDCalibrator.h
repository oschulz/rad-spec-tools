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


#ifndef SDCALIBRATOR_H
#define SDCALIBRATOR_H
#include <memory>

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TList.h>

#include"SDFitData.h"
namespace rspt{

double SQRTQuadFunct(double *x, double *par);
    
class SDCalibrator
{
public:
	SDCalibrator();
	virtual ~SDCalibrator();

	int calibrate();
	void addResult(SDFitData* data);
	TList* getCalObjects(){return m_objects;};
	void setupCalGraphs();

	double intercept;
	double slope;

	double getIntercept() {return intercept;};
	double getSlope() {return slope;};

protected:

	TCanvas *calCanv;
	TF1 *rescal_ch2fch;
	TF1 *rescal_e2fe;
	TF1 *cal_e2ch;
	TF1 *cal_ch2e;

	TList* m_objects;

	TGraphErrors *cal_graph;
	TGraphErrors *rescal_graph;

	void init();


};
} //namespace rspt
#endif // SDCALIBRATOR_H
