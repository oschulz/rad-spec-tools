// Copyright (C) 2013 Oliver Schulz <oschulz@mpp.mpg.de>,
//               2014 Lucia Garbini <garbini@mpp.mpg.de>,
//               2014 Thomas Quante <thomas.quante@tu-dortmund.de>

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


#ifndef RSPT_SDPRECAL_H
#define RSPT_SDPRECAL_H

#include <limits>
#include <cmath>

#include <Rtypes.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>


namespace rspt {


class SDPreCal {
public:
	template<typename tp_Type> class DescriptiveStatistics {
	protected:
		size_t m_n;
		tp_Type m_s0, m_s1, m_s2;

		bool m_first;

	public:
		static tp_Type nan() { return std::numeric_limits<tp_Type>::quiet_NaN(); }

		size_t n() const { return m_n; }
		tp_Type sum() const { return m_s1; }
		tp_Type mean() const { return (n() > 0) ? m_s1 / m_s0 : nan(); }
		tp_Type var() const { return (n() > 0) ? (m_s0 * m_s2 - m_s1 * m_s1) / (m_s0 * m_s0) : nan(); }
		tp_Type sigma() const { return (n() > 0) ? sqrt(m_s0 * m_s2 - m_s1 * m_s1) / m_s0 : nan(); }

		void clear() { m_n = 0; m_s0 = 0; m_s1 = 0; m_s2 = 0; }

		DescriptiveStatistics& add(tp_Type x, tp_Type w = 1) {
			m_n += 1;
			m_s0 += w;
			m_s1 += x * w;
			m_s2 += x*x * w;
			return *this;
		}

		DescriptiveStatistics() { clear(); }
	};

	typedef DescriptiveStatistics<double> Stats;
	typedef std::pair< DescriptiveStatistics<double>, DescriptiveStatistics<double> > Stats_pair;
	typedef std::vector<std::pair< int, int > > Mapping;
	typedef std::pair< double,double > Line;

	struct next_line_info {
		short s_ind_a;
		short s_ind_b;
		short d_ind_a;
		short d_ind_b;
		Stats stats;
	};

	SDPreCal(double l_mca, double h_mca);
	virtual ~SDPreCal();

	TGraph* calcPreCal(std::vector<double> source_lines, std::vector<double> data_lines);
	void setDistThres(double thres) { m_dist_thres = thres; }
	next_line_info genLineInfo(next_line_info prev, int next_s, int next_d);

protected:
	bool debug;
	std::vector<double> m_source_collection;
	std::vector<double> m_data_collection;
	int m_source_size;
	int m_data_size;

	double m_dist_thres;
	double m_adc_max;
	double m_low_mca;
	double m_high_mca;

	TF1* fit;
	TGraph *precal_graph;



	inline double calcError(SDPreCal::Stats x) {return std::abs(double(1) - x.mean());};

	
	SDPreCal::Stats match(next_line_info next);

	std::pair<SDPreCal::Mapping, SDPreCal::Stats> genMap(next_line_info next, SDPreCal::Mapping prevMap);

    
};


} // namespace rspt

#endif // RSPT_SDPRECAL_H
