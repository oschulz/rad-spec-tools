/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2014  <copyright holder> <email>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#ifndef SDMULTILINEFITTER_H
#define SDMULTILINEFITTER_H
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "SDCalibrator.h"
#include "SDFitData.h"
namespace rspt{
class SDMultiLineFitter
{

public:
    SDMultiLineFitter();
    virtual ~SDMultiLineFitter();
    void setThreshold(double thresh);
    void setSigma(float sig);
    void setPreCal(double slope,double intercept);
    void setPreCal(TF1 *precal_ch2e);
    void setWidth(double width){m_width=width;};
    void setRange(double lowEdge,double highEdge);

    void resetPreCal();
    std::vector<rspt::SDFitData*> makeCalFits(TH1* raw_hist,
                     std::vector<double> energy,
                     std::vector<bool> *reject_res_cal=0);
protected:
    int Np;
    
    float m_sigma;

    int m_iteration;
    double m_threshold;
    double m_low_limit;
    double m_high_limit;
    double m_width;
    ///x position of the found peaks
    float *m_specXPeak;
    ///y position of the found peaks
    float *m_specYPeak;
    SDCalibrator calibrator;
    
    TH1 *m_raw_hist;
    TCanvas *m_cal_canv;
    TF1 *m_preCalibration_ch2e;
    TF1 *m_preCalibration_e2ch;
    
    void init();
    std::pair<double,int> getRange(std::vector<double> energy,int iter,int lines_to_fit);
};
}
#endif // SDMULTILINEFITTER_H
