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


#include "SDFitData.h"


namespace rspt{


SDFitData::SDFitData(TF1 *fit,int npeaks)
	: m_fit(fit), m_valid(false), m_npeaks(npeaks)
{
	if (m_fit == 0) {
		std::cerr << "fit is invalid" << std::endl;
	} else {
		m_valid=true;
		for (int i = 0; i < npeaks; ++i) {
			m_usage.push_back(false);
			m_energy.push_back(0);
			m_res_usage.push_back(true);
		}
	}
}


SDFitData::~SDFitData() {
}


double SDFitData::getMean(unsigned int index) {
	if (m_valid) {
		int ind = m_fit->GetParNumber(Form("peak%i_center",index));
		if (ind != -1) {
			return m_fit->GetParameter(ind);
		} else {
			std::cerr<<"getMean index out of range"<<std::endl;
			return 0;
		}
	}
}


double SDFitData::getMeanError(unsigned int index) {
	if (m_valid) {
		int ind=m_fit->GetParNumber(Form("peak%i_center",index));
		if (ind != -1) {
			return m_fit->GetParError(ind);
		} else {
			std::cerr << "getMeanError index out of range" << std::endl;
			return 0;
		}
	}
}


double SDFitData::getSigma(unsigned int index) {
	if (m_valid) {
		int ind = m_fit->GetParNumber(Form("peak%i_sigma", index));
		if (ind != -1) {
			return m_fit->GetParameter(ind);
		} else {
			std::cerr<<"getSigma index out of range" << std::endl;
			return 0;
		}
	}
}

double SDFitData::getSigmaError(unsigned int index) {
	if (m_valid) {
		int ind = m_fit->GetParNumber(Form("peak%i_sigma",index));
		if (ind != -1) {
			return m_fit->GetParError(ind);
		} else {
			std::cerr << "getSigmaError index out of range" << std::endl;
			return 0;
		}
	}
}


bool SDFitData::getUsage(unsigned int index){
	if (m_valid) {
		if (index <= m_npeaks&&index > 0) {
			return m_usage[index-1];
		} else {
			std::cerr<<Form("getUsage(%i) index out of range", index) << std::endl;
			return false;
		}
	}
}


bool SDFitData::getResUsage(unsigned int index) {
	if (m_valid) {
		if (index <= m_npeaks && index > 0) {
			return m_res_usage[index-1];
		} else {
			std::cerr << "getResUsage index out of range" << std::endl;
			return false;
		}
	}
}


double SDFitData::getEnergy(unsigned int index) {
	if (m_valid) {
		if (index <= m_npeaks&&index > 0) {
			return m_energy[index - 1];
		} else {
			std::cerr << "getEnergy index out of range" << std::endl;
			return 0;
		}
	}
}


void SDFitData::setEnergy(unsigned int index, double energy) {
	if (index <= m_npeaks && index > 0){
		m_energy[index-1]=energy;
	} else {
		std::cerr << "setEnergy index out of range" << std::endl;
	}
}


bool SDFitData::setResUsage(unsigned int index, bool use) {
	if (m_valid) {
		if (index <= m_npeaks && index > 0) {
			m_res_usage[index-1]=use;
			return true;
		} else {
			std::cerr << "setResUsage index out of range" << std::endl;
			return false;
		}
	}
}


bool SDFitData::setUsage(unsigned int index,bool use) {
	if (m_valid) {
		if (index <= m_npeaks && index > 0) {
			m_usage[index-1] = use;
			return true;
		} else {
			std::cerr << "setUsage index out of range" << std::endl;
			return false;
		}
	}
}


} // namespace rspt
