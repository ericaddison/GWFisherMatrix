/*--------------------------------
** Mass Gap Research Simulation
** USU FAST Group
** Eric Addison
*--------------------------------*/ 

// GravWave.h
// This is a header file for the GravWave namespace

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include "cosmo.h"
#include "gslWrappers.h"

#ifndef _GW_
#define _GW_

namespace GravWave {
	using gslWrappers::gslRootFinder;

	// scaling values
	const double DSCALE=(1e6*METERPERPARSEC);	// megaparsec
	const double MSCALE=MSUN;					// msun

	class GW_Bin_Insp {
	/// grav wave signal for an inspiralling binary
	private:
	// private data members
		double m1, m2, M;	// mass parameters: m1, m2, total mass,
		double mu, Mc, eta;	// reduced mass, chirp mass, dimensionless mass ratio
		double dm;		// mass difference
		double inc;		// inclination angle
		double Rl;		// luminosity distance
		double alpha, C1; // constant needed for calculations
		double f0;	// min (starting) frequency
		double flso;	// freq. of last stable orbit
		double c,s;		// cos(inc) and sin(inc)
		double N, dt;	// number of points in full output and time step
		double tc, tmax;		// time to coalescence
		double phi_params[3];	// parameters for phase phi(t) expansion		
		double om_params[3];	// parameters for freq omega(t) expansion
		double Hp[5], Hp_pars[16];	// values required for Hplus 
		double cp[7], sp[6];	// cosines and sines of multiples of psi
		double psi, phi_t, om_t; // phase related values
		double Theta, T_8, phi_c; // more phase related values
		gslRootFinder rooter;	// gsl root finding
	// private member functions
		void calc_phi_params();		// calculate phi parameters
		void calc_om_params();		// calculate om parameters
		void set_other_params();	// calculate derived mass parameters
		void set_Hp_params();		// calc parameters for Hp computation
		double get_omega(double &);	// get current value omega(t)
		double get_phi(double &);	// get current valut phi(t)
		double calc_tc();	// calculate time to coal. tc
		void H_plus(); // calculate the Hplus coefficients for different orders
	public:
	// contructors & destructor
		GW_Bin_Insp(double m = 1.4, double n = 1.4, double i = 0.0, double R = 10.0, double f = 0.1, double T = 0.0);
		~GW_Bin_Insp();	
	// set functions
		void set_f0(double);	// set f0 value
		void set_m1(double);	// set m1 value
		void set_m2(double);	// set m2 value
		void set_Rl(double);	// set Rl value
		void set_inc(double);	// set inclination value
		void set_N(double);		// set number of points N
	// get functions
		double get_hp(double);	//calculate h_p(t)
		double get_m1() {return m1;}
		double get_m2() {return m2;}
		double get_f0() {return f0;}
		double get_Rl() {return Rl;}
		double get_inc() {return inc;}
		double get_M() {return M;}
		double get_eta() {return eta;}
		double get_mu() {return mu;}
		double get_Mc() {return Mc;}
		double get_tc() {return tc;}
		double get_flso() {return flso;}
		double get_tmax() {return tmax;}
	// other functions
		void write_to_file(const char * filename = "wave");
	// friends: root finding functions needed for gsl root finder
		friend double fomega(double, void *);
		friend double dfomega(double, void *);
		friend void fdfomega(double, void *, double *, double *);
		friend class Fisher;
	};


// test function for GW_Bin_Insp class
	void GW_Test();

// root finding function
	double fomega(double, void *);
	double dfomega(double, void *);
	void fdfomega(double, void *, double *, double *);


	

} // END namespace GravWave

#endif
