// Here begins file MyBonkers.cc

#include <iostream>
#include <fstream> 
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <ctime>

using namespace std;

#define ranf() \
  ((double)random()/(1.0+(double)RAND_MAX)) // Uniform from interval [0,1)
 
struct particle {
  double x    ; // (x,y) coordinates
  double p    ; // momentum
  double im   ; // inverse mass
  double v    ; // velocity
  double T    ; // kinetic energy
  double a    ; // radius of particle
} ; // Note the definition of a structure ends with a semicolon

struct control {
  int verbose    ; // program verbosity
  int printing   ; // period with which to print
  int commanding ; // period with which to issue commands
  int time_in_ms ; // time between frames (commands)
  clock_t latest_time;
  ofstream fout  ; 
};

// Implement     x += v dt for one particle
void PositionStep ( particle &a , double dt )
{ 
  a.x += dt * a.v ; 
}

void v2p( particle &a )
  // propagate the changed velocity into the momentum vector
{
  a.p =  a.v / a.im ;
}

void p2v ( particle &a ) {
  a.v  = a.p * a.im ;
}

void pv2T ( particle &a ) {
  a.T = 0.5*a.v * a.p ; 
}

void collide ( particle & a , particle & b  )
  // two particles collide elastically
{
  // find the relative velocity 
  double velocity_along_line =  a.v - b.v;
  // find a's mass fraction
  double af = a.im / ( a.im + b.im ) ; 
  double bf = b.im / ( a.im + b.im ) ; 
  // reverse the c.o.m. velocity of each along the line of collision
  a.v -= 2.0 * af * velocity_along_line ;
  b.v += 2.0 * bf * velocity_along_line ;
  v2p( a ) ; 
  v2p( b ) ; 
}

// find next collision time 
double nct( particle *a , int NN , int &whichn ) {
  double dt = 1e100 ; 
  // examine all adjacent pairs from 0,1   to  NN-2,NN-1
  for ( int n = 0 ; n < NN-1 ; n ++ ) {
    double relative_velocity = a[n].v - a[n+1].v ;
    if(relative_velocity > 0.0) {
      double collision_time =  ((a[n+1].x-a[n+1].a) - (a[n].x+a[n].a))
	                        /relative_velocity ;
      if ( collision_time < dt ) {
	dt = collision_time ;
	whichn = n ;
      }
    }
  }
  return dt ;
}

void leapForward(  particle *a , int NN , double dt ) {
  for( int n = 0 ; n < NN ; n ++ ) 
    PositionStep( a[n] , dt ) ; 
}

void showState ( particle *a , int n0 ,  int NN, ostream &fout )
{
  for( int n = n0 ; n < NN ; n ++ ) {
      fout << "\t"<<a[n].x;
      fout << "\t"<<a[n].v;
  }
  fout << endl;
}

double kineticEnergy ( particle *a , int NN )
{
  double E = 0.0 ; 
  for( int n=0 ; n < NN ; n ++ ) {
      if ( a[n].im > 0.0 ) // avoid infinitely massive objects
	E+=0.5*a[n].v*a[n].v/a[n].im;
  }
  return E ;
}

void  simulateBonkers( particle *a , int NN ,
		       double &t , double dt , double T ,
		       control &c ) {
  double next_print_dt = 0.0, next_collision_dt ;
  int whichn ;
  int we_are_printing_this_time = 1 ; 
  for(;t<=T;) {
    if( we_are_printing_this_time ) {
      c.fout << t << "\t" << kineticEnergy(a,NN) ;
      showState ( a , 0 , NN , c.fout ) ;
      next_print_dt = dt ;
    }
    // find the next event
    next_collision_dt = nct( a , NN , whichn ) ;
    // ^^ this returns the time, and sets 'whichn'
    if ( next_collision_dt < next_print_dt ) {
      // advance time to that event, have a collision
      leapForward( a ,  NN , next_collision_dt ) ;
      t += next_collision_dt ;
      next_print_dt -= next_collision_dt ;
      collide( a[whichn] , a[whichn+1] ) ;
      we_are_printing_this_time = 0 ; 
    } else {
      leapForward( a , NN , next_print_dt ) ;
      t += next_print_dt ;
      we_are_printing_this_time = 1 ; 
    }
  }
}

// simulation of gas with identical mass of molecules
void pistonSimulation(int N1, int N2, double piston_position, double piston_velocity)
{

	int N = N1 + N2 + 1; // including the piston
	double T = 1000.0; // target time in seconds
	particle *a;
	control c;
	char filename[50]="file1.dat";

	double dt = 0.03; // time step
	double t = 0.0;
	c.verbose = 1;
	c.time_in_ms = 100;

	c.fout.open(filename);

	srandom(time(NULL));

	a = new particle[N+2]; // including the walls

	// random positions and velocities of all particles on both sides
	for (int n=1; n<=N1; n++) {
		a[n].im = 1.0;
		a[n].x = ranf()*(piston_position);
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}

	for (int n=N1+2; n<=N; n++) {
		a[n].im = 1.0;
		a[n].x = ranf()*(1.0*N - piston_position)+piston_position;
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}
	// piston
	a[N1+1].im = 0.001; // 1000 times heavier
	a[N1+1].x = piston_position;
	a[N1+1].v = piston_velocity;
	a[N1+1].a = 0.0;
	// walls
	a[0].im = a[N+1].im = 0.0; // infinite masses
	a[0].x = 0.0; a[N+1].x = 1.0*N;
	a[0].v = a[N+1].v = 0.0;
	a[0].a = a[N+1].a = 0.0;

	for (int n=0; n<N+2; n++) v2p(a[n]);
// For easier changing of the gnuplot commands
	simulateBonkers( a , N+2 , t , dt , T , c );
	cout<<"Column corresponding to the piston: "<<2*(N1+1)+3<<endl;
}

// simulation of gas with random mass of molecules
void pistonSimulation2(int N1, int N2, double piston_position, double piston_velocity)
{

	int N = N1 + N2 + 1; // including the piston
	double T = 1000.0; // target time in seconds
	particle *a;
	control c;
	char filename[50]="file2.dat";

	double dt = 0.03; // time step
	double t = 0.0;
	c.verbose = 1;
	c.time_in_ms = 100;

	c.fout.open(filename);

	srandom(time(NULL));

	a = new particle[N+2]; // including the walls

	// random positions and velocities of all particles on both sides
	for (int n=1; n<=N1; n++) {
		a[n].im = ranf()*10.0;
		a[n].x = ranf()*(piston_position);
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}

	for (int n=N1+2; n<=N; n++) {
		a[n].im = ranf()*10.0;
		a[n].x = ranf()*(1.0*N - piston_position)+piston_position;
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}
	// piston
	a[N1+1].im = 0.001; // a few 100 times heavier
	a[N1+1].x = piston_position;
	a[N1+1].v = piston_velocity;
	a[N1+1].a = 0.0;
	// walls
	a[0].im = a[N+1].im = 0.0; // infinite masses
	a[0].x = 0.0; a[N+1].x = 1.0*N;
	a[0].v = a[N+1].v = 0.0;
	a[0].a = a[N+1].a = 0.0;

	for (int n=0; n<N+2; n++) v2p(a[n]);
// For easier changing of the gnuplot commands
	simulateBonkers( a , N+2 , t , dt , T , c );
	cout<<"Column corresponding to the piston: "<<2*(N1+1)+3<<endl;
}

void tenParticles()
{

	int N = 10;
	double T = 1000.0; // target time in seconds
	particle *a;
	control c;
	char filename[50]="file3.dat";

	double dt = 0.03; // time step
	double t = 0.0;
	c.verbose = 1;
	c.time_in_ms = 100;

	c.fout.open(filename);

	srandom(time(NULL));

	a = new particle[N+2]; // including the walls
	
	// random positions and velocities of all particles
	for (int n=1; n<=N; n++) {
		a[n].im = ranf()*5.0;
		a[n].x = ranf()*N;
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}

		a[2].im = 1.0; // put these two in order to get the ratio 4:1
		a[6].im = 4.0;

	// walls
	a[0].im = a[N+1].im = 0.0; // infinite masses
	a[0].x = 0.0; a[N+1].x = 1.0*N;
	a[0].v = a[N+1].v = 0.0;
	a[0].a = a[N+1].a = 0.0;

	simulateBonkers( a , N+2 , t , dt , T , c ); // and will plot columns 9 and 17
}

void movingWall(int number, double velocity)
{

	
	int N = number;
	double T = 50.0; // target time in seconds
	particle *a;
	control c;
	char filename[50]="file4.dat";

	double dt = 0.03; // time step
	double t = 0.0;
	c.verbose = 1;
	c.time_in_ms = 100;

	c.fout.open(filename);

	srandom(time(NULL));

	a = new particle[N+2]; // including the walls
	
	// random positions and velocities of all particles
	for (int n=1; n<=N; n++) {
		a[n].im = ranf()*5.0;
		a[n].x = ranf()*N;
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}

	// walls
	a[0].im = a[N+1].im = 0.0; // infinite masses
	a[0].x = 0.0; a[N+1].x = 1.0*N;
	a[0].v = 0.0; a[N+1].v = velocity; // giving the right-hand wall velocity
	a[0].a = a[N+1].a = 0.0;

	simulateBonkers( a , N+2 , t , dt , T , c );

}

void swap_check(double &x, double &y)
{
	int temp;
	if (x>y) {temp=y; y=x; x=temp;}
}

void twoPistons(int N1, int N2, int N3, double piston_position1, double piston_position2)
{
// we asssume that piston_position1 is the position of the left piston
	swap_check(piston_position1, piston_position2);
	
	int N = N1 + N2 + N3 + 2; // including the two pistons
	double T = 500.0; // target time in seconds
	particle *a;
	control c;
	char filename[50]="file5.dat";

	double dt = 0.03; // time step
	double t = 0.0;
	c.verbose = 1;
	c.time_in_ms = 100;

	c.fout.open(filename);

	srandom(time(NULL));

	a = new particle[N+2]; // including the walls

	// random positions and velocities of all particles on both sides
	for (int n=1; n<=N1; n++) {
		a[n].im = ranf()*10.0;
		a[n].x = ranf()*(piston_position1);
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}

	for (int n=N1+2; n<=N1+N2+1; n++) {
		a[n].im = ranf()*10.0;
		a[n].x = ranf()*(piston_position2 - piston_position1)+piston_position1;
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}

	for (int n=N1+N2+3; n<=N; n++) {
		a[n].im = ranf()*10.0;
		a[n].x = ranf()*(1.0*N - piston_position2)+piston_position2;
		a[n].v = ranf()*5.0;
		a[n].a = 0.0;
	}
	// piston1
	a[N1+1].im = 0.001; // a few 100 times heavier
	a[N1+1].x = piston_position1;
	a[N1+1].v = 0.0;
	a[N1+1].a = 0.0;
	// piston2
	a[N1+N2+2].im = 0.001; // a few 100 times heavier
	a[N1+N2+2].x = piston_position2;
	a[N1+N2+2].v = 0.0;
	a[N1+N2+2].a = 0.0;
	// walls
	a[0].im = a[N+1].im = 0.0; // infinite masses
	a[0].x = 0.0; a[N+1].x = 1.0*N;
	a[0].v = a[N+1].v = 0.0;
	a[0].a = a[N+1].a = 0.0;

	for (int n=0; n<N+2; n++) v2p(a[n]);

	simulateBonkers( a , N+2 , t , dt , T , c );
// For easier changing of the gnuplot commands
	cout<<"Columns corresponding to the pistons: "<<2*(N1+1)+3<<" "<<2*(N1+N2+2)+3<<endl;

}

int main(int argc, char* argv[])
{
	// pistonSimulation(5, 5, 1.5, 0.0);

	// pistonSimulation2(5, 5, 1.5, 0.0);

	// tenParticles();

	// movingWall(10, -0.1); // plot energy as function of time

	// twoPistons(10, 10, 10, 12.0, 23.0);


  return 0;
}

