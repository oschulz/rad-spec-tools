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


#include "SDFitData.h"

#include <stdexcept>

using namespace std;


namespace rspt{


SDFitData::SDFitData(TF1 *fit, size_t npeaks)
	: m_fit(fit), m_valid(false), m_npeaks(npeaks)
{
	if (m_fit == 0) throw runtime_error("fit is invalid");

	m_valid=true;
	for (size_t i = 0; i < npeaks; ++i) {
		m_usage.push_back(false);
		m_energy.push_back(0);
		m_res_usage.push_back(true);
	}
}


SDFitData::~SDFitData() {
}


double SDFitData::getMean(size_t index) const {
	return getParValue("peak%i_center",index);
}


double SDFitData::getMeanError(size_t index) const {
	getParError("peak%i_center",index);
}


double SDFitData::getSigma(size_t index) const {
	getParValue("peak%i_sigma", index);
}


double SDFitData::getSigmaError(size_t index) const {
	getParError("peak%i_sigma",index);
}


bool SDFitData::getUsage(size_t index) const {
	if (!m_valid) throw invalid_argument("Invalid SDFitData");

	if (index < 1) throw range_error("index out of range");
	return m_usage.at(index-1);
}


void SDFitData::setUsage(size_t index, bool use) {
	if (!m_valid) throw invalid_argument("Invalid SDFitData");

	if (index < 1) throw range_error("index out of range");
	m_usage.at(index-1) = use;
}


bool SDFitData::getResUsage(size_t index) const {
	if (!m_valid) throw invalid_argument("Invalid SDFitData");

	if (index < 1) throw range_error("index out of range");
	return m_res_usage.at(index-1);
}


void SDFitData::setResUsage(size_t index, bool use) {
	if (!m_valid) throw invalid_argument("Invalid SDFitData");

	if (index < 1) throw range_error("index out of range");
	m_res_usage.at(index-1) = use;
}


double SDFitData::getEnergy(size_t index) const {
	if (!m_valid) throw invalid_argument("Invalid SDFitData");

	if (index < 1) throw range_error("index out of range");
	return m_energy.at(index - 1);
}


void SDFitData::setEnergy(size_t index, double energy) {
	if (index < 1) throw range_error("index out of range");
	m_energy.at(index-1) = energy;
}


size_t SDFitData::getParNumber(const char* nameForm, size_t index) const {
	if (!m_valid) throw invalid_argument("Invalid SDFitData");

	ssize_t parNo = m_fit->GetParNumber(TString::Format(nameForm, int(index)));
	if (parNo < 0) throw out_of_range("No such fit parameter");
	else return parNo;
}


double SDFitData::getParValue(const char* nameForm, size_t index) const {
	return m_fit->GetParameter(getParNumber(nameForm, index));
}


double SDFitData::getParError(const char* nameForm, size_t index) const {
	return m_fit->GetParError(getParNumber(nameForm, index));
}


} // namespace rspt
