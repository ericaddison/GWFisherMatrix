/*
	Eric Addison
	Utah State University
	Summer 2009

	This header file contains definitions for cosmology routines
*/


#include <stdio.h>
#include <math.h>

#ifndef _COSMO_
#define _COSMO_

#define METERPERPARSEC 3.08568025e16
#define LYPERPARSEC 3.26163626
#define METERPERLY 9.4605284e15
#define SECPERYEAR 31556926e0
#define PI 3.1415926535897932384626433832795
#define Omega 1.9924e-7                      //angular velocity of earth around the sun
#define METERPERAU 1.496e11                  //meters per AU
#define MSUN 1.98892e30
#define H0 1/(14e9*SECPERYEAR*C)		//Hubble constant in seconds^-1
#define R_WD 7e6
#define R_NS 10e3
#define R_BH 30e3						//wikipedia info for stellar mass BH
#define RSUN 6.955e8;
#define R_MY 100000						// radius of the milky way in light years
#define D_CG 8						// distance from the sun to the center of the galaxy, in kPc, BOB 882

const double G=6.673e-11;
const double C=299792458e0;

//function prototypes

double m_to_ly(double x);
double ly_to_m(double x);
double m_to_au(double x);
double m_to_par(double x);
double par_to_m(double x);
double par_to_ly(double x);
double ly_to_par(double x);
double sec_to_yr(double x);
double yr_to_sec(double x);

double Rl_to_z(double Rl);
double f1(double z1, double Rl);
double z_to_Rl(double z);

double a_to_tc(double a, double m1, double m2);
double f(double x, double a0, double m1, double m2);
double f_to_tc(double f0, double m1, double m2);
double tc_to_omega(double tc, double t1, double m1, double m2);
double tc_to_a(double tc, double t1, double m1, double m2);
double L_E_to_e(double E, double l, double m1, double m2);
double E_to_a(double E, double m1, double m2);

double bisect_t_c(double a, double b, double f(double x, double a0, double m1, double m2), double tol, double a0, double m1, double m2);
double bisect_z(double a, double b, double f(double z1, double Rl), double Rl, double tol);

#endif
