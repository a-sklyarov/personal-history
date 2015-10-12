// lattice.hh

#ifndef LATTICE_HH
#define LATTICE_HH

#include <vector>

using namespace std;

// assume that Boltzman constant k_B = 1 (natural units)
const long J = 1; // exchange energy
const long mu = 1; // magnetic moment

class Lattice
{
	private: 
		// spin state of the lattice with
		// each spin being either +1 or -1
		vector<long> m_state;
		
		long m_N; // dimendion of the N x N lattice
		long m_E; // energy of the current state
		long m_M; // magnetization of the current state
		double m_T; // temperature
		double m_H; // magnetic field
		
		// Array with all possible values of the acceptance 
		// probability - exp(-deltaE / kT). In the absence of a 
		// magnetic field there are only five possible values of deltaE
		// from which only three are positive. In the presence of a
		// magnetic field, however, the possible values of deltaE
		// become 10. This still allows to create an array with all
		// possible values of the acceptance probability and save the
		// evaluation of an exponent at every flip step.
		double exp_vals[10];
		
	public:
		//constructors
		
		// The constructors also initialize the energy and
		// magnetization for the initial state. The consequential
		// changes are handled by the flip method
		Lattice (long N, double T, double H);
		Lattice (long N, double T, vector<long> state, double H);
		 
		// setters and getters
		
		void SetT (double T); // the only setter we need
		void SetH (double H);
		
		vector<long> GetState ();
		long GetN ();
		long GetE ();
		long GetM ();
		double GetT ();
		double GetH ();
		
		// Methods
		
		// Function giving the initial energy of the lattice
		long initial_energy (vector<long> data, long N);
		
		// This method will be used to carry out a single flip-step
		// from the Metropolis algorithm. The call to this method asks
		// for a specific spin to be flipped and the method flips it
		// only if it agrees with the conditions specified in the 
		// Metropolis algorithm. The method also updates the values
		// of the energy and magnetization once the spin is flipped.
		void flip (long row, long col, double p);
};

#endif
