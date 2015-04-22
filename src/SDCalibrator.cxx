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


#include "SDCalibrator.h"

#include <TMath.h>

#include "rsptutils.h"

using namespace std;


namespace rspt{


SDCalibrator::SDCalibrator()
	: m_objects(0)
{
    rescal_graph == 0;
    cal_graph == 0;
	init();
}


SDCalibrator::~SDCalibrator() {
}

void SDCalibrator::clear()
{
	m_objects->Clear();
	delete m_objects;
	m_objects=0;
	if(rescal_ch2fch!=0){
		delete rescal_ch2fch;
		delete rescal_e2fe;
	}

	if(cal_ch2e!=0){
		delete cal_ch2e;
		delete cal_e2ch;
	}

	if(cal_graph!=0){
		delete cal_graph;
		delete rescal_graph;
	}

	if(calCanv !=0){
		delete calCanv;
	}
	
	init();
}


void SDCalibrator::init() {
	m_objects = new TList();

	if (rescal_graph == 0 || cal_graph == 0) setupCalGraphs();

	rescal_ch2fch = 0;
	rescal_e2fe = 0;
	cal_ch2e = 0;
	cal_e2ch = 0;
	cal_graph = 0;
	rescal_graph = 0;
	calCanv = 0;
}


void SDCalibrator::setupCalGraphs() {
	cal_graph = new TGraphErrors();
	cal_graph->SetNameTitle("cal_graph", "calibration graph");
	rescal_graph = new TGraphErrors();
	rescal_graph->SetNameTitle("rescal_graph", "resolution graph");
	m_objects->Add(cal_graph);
	m_objects->Add(rescal_graph);
}


void SDCalibrator::addResult(SDFitData* data) {
	if (cal_graph == 0 && rescal_graph == 0) setupCalGraphs();

	int nAdds=0;
	for (unsigned int i=1; i<=data->getNPeaks(); i++) {
		if ( data->getUsage(i) ) {
			nAdds++;

			// add new calibration point to old graph
			cerr << "to put in the graph: Energy = "<<data->getEnergy(i)<<"\t ADC channel = "<<data->getMean(i) << endl;
			cal_graph->SetPoint(cal_graph->GetN(),data->getEnergy(i),data->getMean(i));
			cal_graph->SetPointError(cal_graph->GetN()-1,0,data->getMeanError(i));

			// add new resolution point to old graph
			if (data->getResUsage(i)) {
				rescal_graph->SetPoint(rescal_graph->GetN(),data->getMean(i),TMath::Sqrt(8*TMath::Log(2))*data->getSigma(i));
				rescal_graph->SetPointError(rescal_graph->GetN()-1,data->getMeanError(i),TMath::Sqrt(8*TMath::Log(2))*data->getSigmaError(i));
			}
		}
	}
}

int SDCalibrator::calibrate() {
	if (cal_graph!=nullptr&& cal_graph->GetN() >= 2) {
		bool point_removed=false;
		cal_e2ch = new TF1("cal_e2ch","pol1",0,3000); // channel(energy)
		cal_e2ch->SetTitle("Calibration Ch(E)");
		cal_e2ch->SetLineColor(3);
		cal_e2ch->SetLineWidth(1);
		m_objects->Add(cal_e2ch);

		// Fitting resolution with option  "W": Set all weights to 1 for non empty bins; ignore error bars;
		// R: Use the Range specified in the function range
		// Q: Quiet
		cal_graph->Fit("cal_e2ch","WRQ");
		cal_graph->GetXaxis()->SetTitle("Energy");
		cal_graph->GetYaxis()->SetTitle("Channels");
		cal_graph->GetFunction("cal_e2ch")->ResetBit(512);

		// Copy "cal_fit" as "calEqn" (as calibration equation) and invert it (to x=channel, y=energy)
		cal_ch2e = (TF1*)cal_e2ch->Clone("cal_ch2e");
		rspt::transposePol1(cal_ch2e);
		m_objects->Add(cal_ch2e);
		cal_ch2e->SetRange(0,60000);//changed for susie
		cal_ch2e->SetNameTitle("cal_ch2e","Calibration E(Ch)");
		cal_ch2e->GetXaxis()->SetTitle("Channels");
		cal_ch2e->GetYaxis()->SetTitle("Energy");

		m_intercept = cal_ch2e->GetParameter(0);
		m_slope = cal_ch2e->GetParameter(1);

		cerr << "Calibration function: " << cal_ch2e->GetParameter(1) << "*x + " << cal_ch2e->GetParameter(0) << endl;
	} else {
		cerr << "Less than 2 point in calibration graph. Skipping fit for calibration equation.\n";
	}


	if (cal_graph!=nullptr&&rescal_graph->GetN() >= 2) {
		TF1 *rescal_ch2fch_lin = new TF1("rescal_ch2fch_lin", "pol1", 1, 10000);

		// Fitting resolution with option  "W": Set all weights to 1 for non empty bins; ignore error bars;
		// because the lower sigmas have a very much smaller error, and therefore the important high
		// energy lines have nearly no influence on the fit!
		// R: Use the Range specified in the function range
		// Q: Quiet

		rescal_graph->SetTitle("Resolution calibration FWHM_{Ch}(Ch)");
		rescal_graph->GetXaxis()->SetTitle("Channels");
		rescal_graph->GetYaxis()->SetTitle("FWHM_{Ch} / channels");
		rescal_graph->Fit("rescal_ch2fch_lin", "0Q");

		rescal_ch2fch=new TF1("rescal_ch2fch", SQRTQuadFunct, 1, 10000, 3);
		rescal_ch2fch->SetTitle("Resolution calibration FWHM_{Ch}(Ch)");
		rescal_ch2fch->SetLineColor(4);
		rescal_ch2fch->SetLineWidth(1);
		m_objects->Add(rescal_ch2fch);

		rescal_ch2fch->SetParameter(0, rescal_ch2fch_lin->GetParameter(0));
		rescal_ch2fch->SetParameter(1, 0);
		rescal_ch2fch->SetParameter(2, rescal_ch2fch_lin->GetParameter(1));
		delete rescal_ch2fch_lin;


		cerr << "Fit results for calibration equation:" << endl;
		rescal_graph->Fit("rescal_ch2fch");

		rescal_graph->GetFunction("rescal_ch2fch")->ResetBit(512);


		rescal_e2fe=rspt::rescalFCh2Fe(rescal_ch2fch,cal_ch2e);
		m_objects->Add(rescal_e2fe);

		rescal_e2fe->SetLineColor(3);
		rescal_e2fe->SetLineWidth(1);
	} else {
		cerr << "Less than 2 point in resolution graph. Skipping fit for resolution equation." << endl;
	}

	return 1;
}


} //namespace rspt
