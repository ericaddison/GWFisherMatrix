/*--------------------------------
** Mass Gap Research Simulation
** USU FAST Group
** Eric Addison
*--------------------------------*/ 

/* header file for fisher.cpp */

#include <iostream>
#include <cmath>
#include <fstream>
#include "GravWave.h"
#include "gslWrappers.h"

#ifndef _FISHER_
#define _FISHER_

#define PI 3.1415926535897932384626433832795


namespace GravWave
{
	using gslWrappers::gslCSpline;
	using gslWrappers::gslMatrix;

// noise curve class
	class noise_curve
	{
	private:
		double * f; 		// frequenc points
		double * S;			// noise curve values (S_n(f))
		int N;				// numberof points
		double D_NSNS;		// NSNS CBC horizon distance
		double D_NSBH;		// NSBH CBC horizon distance
		double D_BHBH;		// BHBH CBC horizon distance
		char filename[80];	// filename where the noise curve came from
		gslCSpline Spline;	// spline to interpolate 
	public:
	// constructor and destructor
		noise_curve(const char * filename);	// constructor
		~noise_curve();						// destructor
	// get functions
		double * get_f() {return f;}
		double * get_S() {return S;}
		int size() {return N;}
		gslCSpline & get_spline();			// return the spline
	// operator overload
		double operator()(double x) {return Spline(x);}
	};


	class Fisher
	{
	private:
	// static member variables
		static double gamma;	// Euler's Constant
		static double theta;	// constant from Arun et al paper
		static double lam;		// another constant
		static double FC[17];	// array of values needed for alpha coefficients
		static double VSCALE;	// scale value for building freqs. out of masses
		static double ASCALE;	// amplitude scaling value
	// private member data
		double a[8];			// alpha coeffs. for 3.5PN phase expansion
		double dadeta[8];		// derivs of alpha wrt eta
		double dadMc[2];		// derivs of alpga wrt Mc
		double vlso;			// v for last stable circ orbit
		double v_base;			// coefficient for simple conversion of f to v
		double f0;				// parameter used in analytic noise curve approx.
		double fmin;			// lowest frequency in Noise Curve
		double &flso;			// f for last stable circ orbit
		double fmax;			// maximum frequency (usually flso)
		double S0;				// scaling factor for analytic noise curve
		double &eta;			// dimensionless mass ratio
		double A;				// amplitude of wave
		double &Mc;				// chirp massc
		double &tc;				// time to coalescence
		gslMatrix F;			// the Fisher matrix
		gslMatrix cov;			// Covariance matrix
		gslCSpline &Sn;			// noise curve spline object
		GW_Bin_Insp & GW;		// GW object reference
	// private methods
		void calc_a(double);		// calculate alpha coefficients
		void calc_dadeta(double);	// calculate derivatives of alpha wrt eta
		double dvdeta(double);		// deriv of v wrt eta
		void calc_dadMc();			// deriv of alpha wrt Mc
		double dvdMc(double);		// deriv of v wrt Mc
	public:
		Fisher(GW_Bin_Insp & gw, noise_curve & N);	// constructor
		~Fisher();
		void calc_fisher();		// compute the fisher matrix
		double calc_SNR();		// compute SNR
		void show() { F.show(); }	// show the matrix
		void invert() { cov = ~F; }
		double get_cov(int i, int j) { return cov(i,j); } // invert F
		void show_cov() { cov.show(); }
		double Sn_analytic(double);	// analytic version of noise curve aLIGO
		double get_fmax() { return fmax; }
		double set_fmax(double);	// calculate fmax
	// operator overloads
		double & operator()(int i, int j) { return F(i,j); }
	// friends: g_i functions for Fisher matrix calculations
		friend double gij(double, void*);
		friend double g_tc(double, void*);
		friend double g_phic(double, void*);
		friend double g_Mc(double, void*);
		friend double g_eta(double,void*);
		friend double SNR_func(double f, void* params);
	};

/* params struct for gij functions in Fisher class*/
typedef struct{
  int i;
  int j;  //tells which two functions in funcs to use, i.e. gi and gj
  double (*funcs[4]) (double,void*);
  Fisher * F;
}gij_params;


// functions required for calculating Fisher matrix

double gij(double, void*);
double g_tc(double, void*);
double g_phic(double, void*);
double g_Mc(double, void*);
double g_eta(double,void*);
double SNR_func(double f, void* params);


/* Detection Class */
class source
{
private:
// static members
	static noise_curve Sn;	// one noise curve for all source
// private data members
	char type[4];			// type identifier, i.e. NSNS, BHNS, or BHBH
	GW_Bin_Insp gw;			// grav wave object
	Fisher fish;			// Fisher matrix object
	double L_angles[2];		// binary angular momentum angles (pol,az)
	double p_angles[2];		// binary position azimuthal angles (pol,az)
	double SNR;				// SNR of the signal
	double fisco;			// frequency of innermost stable circular orbit
// public methods
public:
};



void Ftest();


}	// END NAMESPACE GRAVWAVE
#endif
