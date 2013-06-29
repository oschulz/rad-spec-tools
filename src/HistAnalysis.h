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


#ifndef RSPT_HISTANALYSIS_H
#define RSPT_HISTANALYSIS_H

#include <TH1.h>
#include <TF1.h>
#include <TSpectrum.h>


namespace rspt {


class HistAnalysis {
public:
	static TSpectrum* findPeaks(TH1 *hist, Option_t* option = "goff", double sigma = 4.0, double threshold = 0.1);

	static TSpectrum* findSigPeaks(TH1 *hist, Option_t* option = "goff", double sigma = 4.0, double threshold = 0.01, Int_t nBgIter = 10);

	static TF1* fitPeaks(TH1 *hist, TSpectrum *peaks, Option_t* option = "", Option_t* goption = "", bool enableSkew = true, const char* bkgModel = "pol2");

	static TF1* findAndFitPeaks(TH1 *hist, Option_t* option = "", Option_t* goption = "", double sigma = 4.0, double threshold = 0.1, bool enableSkew = true, const char* bkgModel = "pol2");

	static void removeBackground(TH1 *hist, Option_t* option = "", Int_t nBgIter = 10, double threshold = 3.5);

	static void filterMinOf3(TH1 *hist);
};


} // namespace rspt


#ifdef __CINT__
#endif

#endif // RSPT_HISTANALYSIS_H
