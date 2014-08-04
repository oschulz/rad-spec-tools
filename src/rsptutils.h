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

#ifndef UTILS_H
#define UTILS_H
#include <TF1.h>
#include "SDFitData.h"
namespace rspt{
    void transposePol1(TF1 **input);
    int desiredPeak(int iter,int fitted_lines, std::vector< double > energy, SDFitData *fit, TF1 *cal_ch2e) ;
    TF1* rescalFCh2Fe(const TF1* rescal_ch2fch, const TF1* cal_ch2e);
} // namespace rspt
#endif // UTILS_H
