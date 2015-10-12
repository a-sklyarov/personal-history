// screen.hh - defines the class Screen

/* This is the definition of the Screen class. A screen has private
 * variables corresponding to the distance to the screen D, a boolean 
 * saying if the pattern has been generated or not, a comlex array
 * containing the pattern in the form of complex numbers, a vector array
 * containing the intensity pattern and a vector array containing the
 * y values corresponding to the pattern. */

#ifndef SCREEN_HH
#define SCREEN_HH

#include "aperture.hh"

using namespace std;

class Screen
{
	private:
		double m_D; // distance to screen
		bool pattern_generated;
		
		cav::Complex_Array1 m_pattern;
		vector<double> m_intensity_pattern;
		vector<double> m_y_values;
		
	public:
		// constructor
		Screen(double D);
		
		// Setters and Getters
		void SetD(double D);
		double GetD();
		
		cav::Complex_Array1 GetPattern();
		
		vector<double> GetIntensityPattern();
		
		vector<double> GetY_Values();
		
		// Methods
		
		// This method takes as an input an aperture which it then
		// uses to generate the diffraction pattern from it.
		
		void GeneratePattern ( Aperture aperture );
};

#endif
		
