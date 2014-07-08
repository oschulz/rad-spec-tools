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


#include "HistAnalysis.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <memory>

#include <TSpectrum.h>
#include <TSpectrumFit.h>

#include "SDPeak.h"

using namespace std;


namespace rspt {


TSpectrum* HistAnalysis::findPeaks(TH1 *hist, Option_t* option, double sigma, double threshold) {
	TSpectrum *spec = new TSpectrum;

	//Copy hist because TSpectrum::Search replaces source spectrum by new spectrum
	TH1 *h = dynamic_cast<TH1*>(hist->Clone());

	spec->Search(h, sigma, option, threshold);
	return spec;
}


TSpectrum* HistAnalysis::findSigPeaks(TH1 *hist, Option_t* option, double sigma, double threshold, Int_t nBgIter, double sigThresh) {
	TSpectrum *spec = new TSpectrum();

	TH1 *fg = dynamic_cast<TH1*>(hist->Clone());

	HistAnalysis::removeBackground(fg, "", nBgIter, sigThresh);
	HistAnalysis::filterMinOf3(fg);

	spec->Search(fg, sigma, option, threshold);
	return spec;
}


TF1* HistAnalysis::fitPeaks(TH1 *hist, TSpectrum *spectrum, Option_t* option, Option_t* goption, bool enableSkew, const char* bgModel) {
	TF1* bgFunc = new TF1("bgModel", bgModel, 0, 10);
	MultiPeakShape peakShape(spectrum->GetNPeaks(), enableSkew, bgFunc);
	TF1* sdpeaks = peakShape.newTF1("sdpeaks", spectrum);
	hist->Fit("sdpeaks", option, goption);
	return sdpeaks;
}


TF1* HistAnalysis::findAndFitPeaks(TH1 *hist, Option_t* option, Option_t* goption, double sigma, double threshold, bool enableSkew, const char* bgModel) {
	TSpectrum *spectrum = findPeaks(hist, "goff", sigma, threshold);
	// for (Int_t pIdx = 0; pIdx < spec->GetNPeaks(); ++pIdx) cout << "Found peak at " << spec->GetPositionX()[pIdx] << " with height " << spec->GetPositionY()[pIdx] << endl;
	TF1 *sdpeaks = fitPeaks(hist, spectrum, option, goption, enableSkew, bgModel);
	delete spectrum;
	return sdpeaks;
}


void HistAnalysis::removeBackground(TH1 *hist, Option_t* option, Int_t nBgIter, double threshold) {
	Int_t n = hist->GetNbinsX();

	auto_ptr<TSpectrum> spec = auto_ptr<TSpectrum>(new TSpectrum);

	auto_ptr<TH1> bg = auto_ptr<TH1>(spec->Background(hist, nBgIter, option));

	// Remove non-significant signal
	for (Int_t i = 1; i <= n; i++) {
		double fv = hist->GetBinContent(i);
		double bv = bg->GetBinContent(i);
		double sv = fv - bv;

		// Transition from 0 to sv within one sigma around threshold:
		sv = sv * (0.5 + TMath::Erf(4 * (sv / sqrt(bv) - threshold) )/2.);

		hist->SetBinContent(i, sv > 0 ? sv : 0);
	}
}


void HistAnalysis::filterMinOf3(TH1 *hist) {
	Int_t n = hist->GetNbinsX();
	if (n >= 2) {
		double last2 = hist->GetBinContent(1);
		double last1 = last2;

		for (Int_t i = 2; i <= n+1; ++i) {
			double current = (i<=n) ? hist->GetBinContent(i) : last1;
			hist->SetBinContent(i-1, std::min(current, std::min(last1,last2)));
			last2 = last1;
			last1 = current;
		}
	}
}


void HistAnalysis::copyBins(std::vector<double> &dest, const TH1 *src) {
	Int_t n = src->GetNbinsX();
	dest.resize(n);
	for (int i=0; i<n; ++i) dest[i] = src->GetBinContent(i+1);
}


void HistAnalysis::copyBins(TH1 *dest, const std::vector<double> &src) {
	Int_t n = src.size();
	if (n != dest->GetNbinsX()) throw invalid_argument("Number of bins in dest histogram must match size of src vector");
	for (int i=0; i<n; ++i) dest->SetBinContent(i+1, src[i]);
}


void HistAnalysis::copyBins(std::vector<float> &dest, const TH1 *src) {
	Int_t n = src->GetNbinsX();
	dest.resize(n);
	for (int i=0; i<n; ++i) dest[i] = src->GetBinContent(i+1);
}


void HistAnalysis::copyBins(TH1 *dest, const std::vector<float> &src) {
	Int_t n = src.size();
	if (n != dest->GetNbinsX()) throw invalid_argument("Number of bins in dest histogram must match size of src vector");
	for (int i=0; i<n; ++i) dest->SetBinContent(i+1, src[i]);
}


} // namespace rspt
