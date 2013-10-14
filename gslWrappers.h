//****************************
//	Eric Addison
//	USU
//	gslWrappers.h
//	header file for wrapper classes to make
//	gsl easier to use in C++
//****************************



#include <iostream>
#include <iomanip>
#include <stdio.h>

// gsl header files 
#ifndef _GSLINC_
#define _GSLINC_
#include <gsl/gsl_rng.h>	// random number generation header
#include <gsl/gsl_randist.h>	// probability distribution functions
#include <gsl/gsl_integration.h> // for numerical integration
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#endif



#ifndef _GSLWRAP_
#define _GSLWRAP_

// Class Declarations

namespace gslWrappers
{

	const double GSL_WRAP_ERR = -98765.4321;

// gslRootFinder: wrapper class for gsl Brent method 1D root finder
	class gslRootFinder
	{
	private:
		const gsl_root_fsolver_type *T;
		gsl_root_fsolver *s;
		double x_lo, x_hi;	// interval end points
		double epsrel, epsabs;	// error tolerances
		double root;		// the root we will find
		int max_iter;		// max iterations
		bool initd;
	public:
		gsl_function F;		// gsl function for root finding	
		gslRootFinder(double xmin, double xmax, double (*f)(double,void*) = NULL, double Niter = 1000, double reltol = 1e-6, double abstol = 0.0);
		~gslRootFinder();
		double findRoot();	// do the root finding
		void rootFinderInit();	// initialize root finder stuff
	// set methods
		void set_max_iter(double Niter) { max_iter = Niter; }
		void set_x_hi(double x) { x_hi = x; }
		void set_x_lo(double x) { x_lo = x; }
		void set_F(double (*f)(double,void*)) { F.function = f; }
		void set_params(void * p) {F.params = p;}
	// get methods
		double get_max_iter() { return max_iter; }
		double get_x_hi() {return x_hi;}
		double get_x_lo() {return x_lo;}
		double get_root() {return root;}
	};


// gslIntegrator: wrapper class for cquad gsl integration
	class gslIntegrator
	{
	private:
		double result;		// holds integration result
		double error;		// integration error
		gsl_function F;		// integrand function
		size_t nevals;		// number of evaluations used
		gsl_integration_cquad_workspace *w;	// integration workspace
		bool initd;			// flag to see if F has been initialized
		bool got_pars;		// flag to see if params has been provided
		double epsrel;		// relative error tolerance
		double epsabs;		// absolute error tolerance
		bool init_checker(const char *);	// initializaton checker
	public:
	// constructors and destructor
		gslIntegrator();
		gslIntegrator(gsl_function G);	
		~gslIntegrator();
		gslIntegrator(double (*f)(double x, void * params), void * params);
		gslIntegrator(double (*f)(double x, void * params));
	// set and get functions
		void set_integrand(gsl_function G);	// change integrand function
		void set_integrand(double (*f)(double x, void * params));	// change integrand function
		void set_params(void * p);	// set params
		void set_epsabs(double e);	// set the absolute error tolerance
		void set_epsrel(double e);	// set the relative error tolerance
		double get_result();		// get the last integration result 
		double get_error();		// get the last integration error
		double get_nevals();		// get the last value of nevals
		double get_epsabs();		// get the absolute error
		double get_epsrel();		// get the relative error
	// integration call
		double integrate(double, double);	// do the integration
		double integrate(double, double, void *); // integrate and provide params
	// operator overload
		double operator()(double,double,void*);	// overload () to use for integration
		double operator()(double,double);
	};

// gslMatrix: wrapper class for the gsl_matrix
	class gslMatrix
	{
	private:
		gsl_matrix * m;		// the actual matrix
		int rows;			// number of rows
		int cols;			// number of columns
	// variables to work with M(i,j) = blah
		double temp;	// reassign value
		int loc[2];			// location of reassigned value
		void reassign();	// reassign a variable if needed
	public:
	// constructors and destructors
		gslMatrix();					// default constructor
		gslMatrix(int, int);			// consructor for empty MxN matrix
		gslMatrix(const gslMatrix &);// copy constructor
		gslMatrix(const gsl_matrix *, int,int);	// another constructor
		gslMatrix(double*,int,int);	// explicit value struct constructor
		~gslMatrix();
	// set and get functions
		int get_rows() { return rows; }	// get the number of rows
		int get_cols() { return cols; } // get the number of cols
		void show(std::ostream & os = std::cout);// print out the matrix
	// other functions
		gslMatrix eye(int n=0);			// get an identity matrix
		gslMatrix inv();				// invert matrix
	// operator overloads
		gslMatrix & operator=(const gslMatrix &);	// assignment operator
		gslMatrix operator+(const gslMatrix &);	// addition
		gslMatrix operator-(const gslMatrix &);	// subtraction
		gslMatrix operator*(const double);		// scalar mult.
		gslMatrix operator-();					// negation
		gslMatrix operator+(const double);		// adda constant
		gslMatrix operator&(const gslMatrix & );	// elementwise mult
		gslMatrix operator|(const gslMatrix & );	// elementwise div
		gslMatrix operator*(const gslMatrix & );	// mat-mat mult
		double & operator()(const int, const int);	// call values
		gslMatrix operator~();						// invert matrix
		gslMatrix operator*=(const gslMatrix &);	// mult. assignment
		gslMatrix operator+=(const gslMatrix &);	// mult. addition
		gslMatrix operator^(int k);					// matrix power
	// friends
		friend gslMatrix operator*(const double, gslMatrix &);
		friend gslMatrix operator+(const double, gslMatrix &);
	};

// gslCSpline: wrapper class for gsl cubic spline interpolation
	class gslCSpline
	{
	private:
	// GSL Interpolation data
		gsl_interp_accel *acc;
		gsl_spline *spline;
		bool initd;				// initialized flag
	public:
		gslCSpline() {initd = 0;}
		gslCSpline(double *, double *, int);
		~gslCSpline();
		double operator()(double);	// overload () to use for interpolation
		void init_spline(double *, double *, int);
	};


// Test functions
	void MatrixTest();
	double func(double, void*);
	void IntTest();
	double gunc(double, void*);
	void rootTest();
	double func2(double, void*);
}	// END NAMESPACE GSLWRAPPERS
#endif
