// aperture.cc - implementation of the class Aperture

/* This is the implementation of the Aperture class. An Aperture has
 * private variables corresponding to the wavelnegth lambda, the
 * number of discrete fourier sample regions N, its maximum extent 
 * L, the size of a sample region delta, a boolean saying if the 
 * aperture has been discretized or not, an analytical function 
 * representing its complex transmittance, a complex array containing
 * its discretized representation (if it has been discretized) and a 
 * vector array containing the x values corresponding to the discretized
 * complex transmittance. 
 */

#include "aperture.hh"

using namespace std;

	// constructors
	Aperture::Aperture() : m_lambda(500e-9), m_N(two_to_the_30), m_maximum_extent(1), m_delta(m_maximum_extent / m_N),
	 m_discretized(0), m_discretized_representation(0), m_discretized_x_values(0)
	{
		m_analytical_function = NULL;
	}
	
	Aperture::Aperture(double lambda_in, int N_in, double maximum_extent_in) 
	 : m_lambda(lambda_in), m_N(N_in), m_maximum_extent(maximum_extent_in), m_delta(m_maximum_extent / m_N),
	 m_discretized(0), m_discretized_representation(0), m_discretized_x_values(0)
	{
		m_analytical_function = NULL;
	}

	// setters and getters
	void Aperture::SetLambda ( double lambda_in ) { m_lambda = lambda_in; }
	double Aperture::GetLambda () { return m_lambda; }
	
	void Aperture::SetN ( double N_in) { m_N = N_in; }
	int Aperture::GetN () { return m_N; }
	
	void Aperture::SetMaximumExtent ( double max_ext_in ) { m_maximum_extent = max_ext_in; }
	double Aperture::GetMaximumExtent () { return m_maximum_extent; }
	
	void Aperture::SetDelta () { m_delta = m_maximum_extent / m_N; }
	
	double Aperture::GetDelta () { return m_delta; }
		
	void Aperture::SetAnalyticalFunction ( complex_transmittance a )
	{
		m_analytical_function = a;
	}
		
	cav::Complex_Array1 Aperture::GetDiscretizedRepresentation ()
	{
		return m_discretized_representation;
	}
	
	vector<double> Aperture::GetDiscretized_X_Values ()
	{
		return m_discretized_x_values;
	}
		
	// Methods
	
	// This method will discretize the analytical complex
	// transmittance.
	
	void Aperture::discretize ()
	{
		if (m_discretized) { return; } // discretized only once
		else
		{
			m_discretized = 1;
			
			// Discretize the analytical function into m_N regions.
			// Indicate every region in m_discretized_x_values and
			// m_discretized_representation.
			
			for (int j = 0; j < m_N; j++)
			{
				double x_j = (j - (m_N / 2))*m_delta; // eq 10
				complex_number transmittance_value = m_analytical_function(x_j);
						
				m_discretized_x_values.push_back(x_j);
				m_discretized_representation.data().push_back(transmittance_value.real_part);
				m_discretized_representation.data().push_back(transmittance_value.imag_part);
			}
		}
	}
			
