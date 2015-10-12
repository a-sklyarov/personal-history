// lattice.cc

#include <cmath>
#include "lattice.hh"

using namespace std;

// Function to apply the periodic boundary condition to an index,
// taking the index and the maximum value, which the index can take
long periodicBC (long index, long max_index)
{
	long x = abs(index);
	
	if (index >= 0)
	return index % max_index; 
	else
	return (max_index - (x % max_index));
}

// Function returning the [row][col] element of a matrix
// with periodic b.c.
long& element (long row, long col, vector<long> data)
{
	long N = (long) sqrt(data.size());
	return data.at( periodicBC(row, N) * N + periodicBC(col, N));
}

// Function giving the initial energy of the lattice
long Lattice::initial_energy (vector<long> data, long N)
{
	long E = 0;
	
	// nearest neighbour interactions	
	for (long i = 0; i < N; i++)
	{
		for (long j = 0; j < N; j++)
		{
			E -= J * element(i, j, data) * ( element(i, j+1, data) + element(i+1, j, data) );
		}
	}
	
	// interaction with external magnetic field
	if (m_H != 0)
	{
		for (long k = 0; k < N*N; k++)
		{
			E -= mu * m_H * data.at(k);
		}
	}
	
	return E;
}

// Function giving the initial magnetization of the lattice
long initial_magn (vector<long> data, long N)
{
	long M = 0;
	
	for (long i = 0; i < N*N; i++)
	{
		M += data.at(i);
	}
	
	return M;
}

// constructors

Lattice::Lattice (long N, double T, double H) 
	: m_state(N*N, 1), m_N(N), m_E(0), m_M(0), m_T(T), m_H(H)
{
	// initialize energy
	m_E = initial_energy( m_state, m_N );
	
	// initialize magnetization
	m_M = initial_magn( m_state, m_N );
	
	// Initialize the array of possible values of the
	// acceptance probability. The first 5 values correspond to 
	// the case when deltaM<0 and the next 5 to when deltaM>0
	
	/* Possible values of deltaE
	 * i		deltaE1				deltaE2				Neighbours
	 * 0		-4 * 2*J + 2*mu*H	-4 * 2*J - 2*mu*H	(all spins -1)
	 * 1		-2 * 2*J + 2*mu*H	-2 * 2*J - 2*mu*H	(3 spins -1)
	 * 2		0 * 2*J + 2*mu*H	0 * 2*J - 2*mu*H	(2 spins -1)
	 * 3		+2 * 2*J + 2*mu*H	+2 * 2*J - 2*mu*H	(3 spins +1)
	 * 4		+4 * 2*J + 2*mu*H	+4 * 2*J - 2*mu*H	(all spins +1)
	 * 
	 * ( 2*i - 4 ) corresponds to the possible values of the sum of
	 * the four nearest neighbours
	 * 
	 * */
	 	
	for (long i = 0; i < 5; i++)
	{
		long deltaE1 = ( 2*i - 4 ) * 2*J + 2*mu*m_H;
		long deltaE2 = ( 2*i - 4 ) * 2*J - 2*mu*m_H;
		
		exp_vals[i] = exp( -deltaE1 / m_T );
		exp_vals[i+5] = exp( -deltaE2 / m_T );
	}
}

Lattice::Lattice (long N, double T, vector<long> state, double H) 
	: m_state(state), m_N(N), m_E(0), m_M(0), m_T(T), m_H(H)
{
	// initialize energy
	m_E = initial_energy( m_state, m_N );
	
	// initialize magnetization
	m_M = initial_magn( m_state, m_N );
		
	// Initialize the array of possible values of the
	// acceptance probability. The first 5 values correspond to 
	// the case when deltaM<0 and the next 5 to when deltaM>0
	
	/* Possible values of deltaE
	 * i		deltaE1				deltaE2				Neighbours
	 * 0		-4 * 2*J + 2*mu*H	-4 * 2*J - 2*mu*H	(all spins -1)
	 * 1		-2 * 2*J + 2*mu*H	-2 * 2*J - 2*mu*H	(3 spins -1)
	 * 2		0 * 2*J + 2*mu*H	0 * 2*J - 2*mu*H	(2 spins -1)
	 * 3		+2 * 2*J + 2*mu*H	+2 * 2*J - 2*mu*H	(3 spins +1)
	 * 4		+4 * 2*J + 2*mu*H	+4 * 2*J - 2*mu*H	(all spins +1)
	 * 
	 * ( 2*i - 4 ) corresponds to the possible values of the sum of
	 * the four nearest neighbours
	 * 
	 * */
	 	
	for (long i = 0; i < 5; i++)
	{
		long deltaE1 = ( 2*i - 4 ) * 2*J + 2*mu*m_H;
		long deltaE2 = ( 2*i - 4 ) * 2*J - 2*mu*m_H;
		
		exp_vals[i] = exp( -deltaE1 / m_T );
		exp_vals[i+5] = exp( -deltaE2 / m_T );
	}
}
			
// setters and getters

void Lattice::SetT (double T) { m_T = T; }	

void Lattice::SetH (double H) { m_H = H; }

vector<long> Lattice::GetState () { return m_state; }

long Lattice::GetN () {	return m_N; }

long Lattice::GetE () {	return m_E; }

long Lattice::GetM () { return m_M; }

double Lattice::GetT () { return m_T; }

double Lattice::GetH () { return m_H; }

// Methods
/*
long& operator Lattice::() (long row, long col)
{
	return element(row, col, m_state);
}
* */

void Lattice::flip (long row, long col, double p)
{
	bool flipped = 0;
	
	// Assign a variable with the value of the element, which
	// is being considered for flipping, and the sum of its
	// nearest neighbours in order to save repeated
	// evaluations of the function element().
	long s = element(row, col, m_state);
	long sum_nn = element(row, col - 1, m_state) +
				  element(row, col + 1, m_state) +
				  element(row - 1, col, m_state) +
				  element(row + 1, col, m_state);
	
	// deltaM is equal to -2 times the spin of the element
	// before it gets flipped.
	long deltaM = -2 * s;
	
	// deltaE is equal to 2J times the spin of the element
	// before it gets flipped times the sum of the spins
	// of its nearest neighbours
	long deltaE = 2 * J * s * sum_nn;
	
	// in presence of external magnetic field the value
	// of deltaE is additionally modified			  
	if ( m_H != 0 )
	{
		deltaE -= mu * m_H * deltaM;
	}
	
	// decide whether to flip the spin or not, based on the 
	// conditions of the Metropolis algorithm
	
	if ( deltaE < 0 ) 
	{
		flipped = 1;
		// flip the spin element
		m_state[ periodicBC(row, m_N) * m_N + periodicBC(col, m_N)] *= -1;
	}
	else if ( ((exp_vals[ (s*sum_nn + 4)/2 ] > p) && ( deltaM <= 0 )) || 
			   ((exp_vals[ (s*sum_nn + 4)/2 + 5 ] > p) && ( deltaM >= 0)) )
	{
		flipped = 1;
		// flip the spin element
		m_state[ periodicBC(row, m_N) * m_N + periodicBC(col, m_N)] *= -1;
	}
	
	// update energy and magnetization of the system
	if ( flipped )
	{
		m_E += deltaE;
		m_M += deltaM;
	}
}
	
			
			
