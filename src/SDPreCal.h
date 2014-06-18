// Copyright (C) 2014 Lucia Garbini <garbini@mpp.mpg.de>

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


#ifndef SDPRECAL_H
#define SDPRECAL_H

#include <Rtypes.h>
#include <TMath.h>


class SDPreCal {
public:
	
	SDPreCal(std::vector<std::pair<double,double>> source_coll, std::vector<std::pair<double,double>> data_coll);
	void calcPreCal(std::vector<std::pair<double, double>> source_lines, std::vector<std::pair<double, double>> data_lines);	
	virtual ~SDPreCal();

protected:

	std::vector<std::pair<double,double>> m_source_collection;
	std::vector<std::pair<double,double>> m_data_collection;
	int m_source_size;
	int m_data_size;

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
	
	typedef DescriptiveStatistics <double> Stats;
	typedef std::vector<std::pair< int, int >> Mapping;
	typedef std::pair< double,double > Line;

	inline double error(Stats x, Stats y) {return sqrt(std::pow(x.sigma(), 2) + std::pow(y.sigma(), 2));};

	std::pair<Stats, Stats> match( Line sline_a, Line sline_b, Line dline_a, Line dline_b, Stats prev_rx, Stats prev_ry );

	std::pair<Mapping, Stats> genMap( int sline_i, int dline_i, Mapping prevMap, Stats prevStats );
	
};



#endif // SDPRECAL_H
