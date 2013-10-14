/*--------------------------------
** USU FAST Group
** Eric Addison
** gslWrappers.cpp
** implementation of gslWrappers classes
*--------------------------------*/ 


#include "gslWrappers.h"

namespace gslWrappers
{		

//************************************************************
//		gslRootFinder
//************************************************************
/// Root finder class using Brent's method

	gslRootFinder::gslRootFinder(double xmin, double xmax, double (*f)(double,void*), double Niter, double reltol, double abstol)
	/// constructor
	// inputs: xmin, xmax = interval endpoints
	// *f = pointer to function to find root of
	// Niter = max number of iteratons -- default = 1000;
	// reltol = relative error tolerance
	// abstol = absolute error tolerance
	{
		x_lo = xmin;
		x_hi = xmax;
	    F.function = f;
		F.params = NULL; // params needs to be explicitly set if needed
		max_iter = Niter;
		epsrel = reltol;
		epsabs = abstol;
		rootFinderInit();
	}

	gslRootFinder::~gslRootFinder()
	/// destructor
	{ 
			gsl_root_fsolver_free (s);
    }

	void gslRootFinder::rootFinderInit()
	/// initialize root finder stuff
	{
	    T = gsl_root_fsolver_brent;
	    s = gsl_root_fsolver_alloc (T);
	}

	double gslRootFinder::findRoot()
	{

	    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
		int status, iter=0;
		double r=0;
		do {
		   iter++;
		   status = gsl_root_fsolver_iterate (s);
		   r = gsl_root_fsolver_root (s);
		   x_lo = gsl_root_fsolver_x_lower (s);
		   x_hi = gsl_root_fsolver_x_upper (s);
		   status = gsl_root_test_interval (x_lo, x_hi,epsabs,epsrel);
		} while (status == GSL_CONTINUE && iter < max_iter);
		root = r;
		return r;
	}


//************************************************************
//		gslIntegrator
//************************************************************
		
	gslIntegrator::gslIntegrator()
	/// default constructor
	{
		w = gsl_integration_cquad_workspace_alloc(10000); // set up workspace
		initd = false; 	// integrand F has NOT been set
		got_pars = false;
		nevals = result = error = 0;
		epsabs = 0;
		epsrel = 1e-7;		// initialized error tolerances
	}

	gslIntegrator::gslIntegrator(gsl_function G)
	/// constructor with gsl_function argument
	{
		w = gsl_integration_cquad_workspace_alloc(10000); // set up workspace
		F = G;
		initd = true; 	// integrand F has been set
		got_pars = true;
		nevals = result = error = 0;
		epsabs = 0;
		epsrel = 1e-7;		// initialized error tolerances
	}

	gslIntegrator::gslIntegrator(double (*f)(double x, void * params))
	/// constructor with function pointer input
	{
		w = gsl_integration_cquad_workspace_alloc(10000); // set up workspace
		F.function = f;
		initd = true; 	// integrand F has been set
		got_pars = false;
		nevals = result = error = 0;
		epsabs = 0;
		epsrel = 1e-7;		// initialized error tolerances
	}

	gslIntegrator::gslIntegrator(double (*f)(double x, void * params), void * params)
	/// constructor with function pointer input and params
	{
		w = gsl_integration_cquad_workspace_alloc(10000); // set up workspace
		F.function = f;
		F.params = params;
		initd = true; 	// integrand F has been set
		got_pars = true;
		result = nevals = error = 0;
		epsabs = 0;
		epsrel = 1e-7;		// initialized error tolerances
	}


	gslIntegrator::~gslIntegrator()
	/// destructor
	{
		gsl_integration_cquad_workspace_free (w);
	}

	void gslIntegrator::set_integrand(gsl_function G)
	/// change integrand function
	{ F = G; }
	
	void gslIntegrator::set_integrand(double (*f)(double x, void * params))
	/// another change integrand function
	{ F.function = f; initd = true; got_pars = false; }

	void gslIntegrator::set_params(void * p)
	/// set params pointer
	{
		F.params = p;
		got_pars = true;
	}

	bool gslIntegrator::init_checker(const char * fun)
	/// check if initialization is complete
	{
		using std::cerr;

		if( !initd )
		{
			cerr << "\nError: gslWrappers::gslIntegrator::" << fun << " -- integrand not initialized";
			return false;
		}
		else if( !got_pars )
		{
			cerr << "\nError: gslWrappers::gslIntegrator::" << fun << " -- params pointer not initialized";
			return false;
		}
		return true;
	}

	double gslIntegrator::get_result()
	/// return result
	{
		if( !init_checker("get_result()") )
			return 0;
		return result;
	}

	double gslIntegrator::get_error()
	/// return error
	{
		if( !init_checker("get_error()") )
			return 0;
		return error;
	}
	
	double gslIntegrator::get_nevals()
	/// return error
	{
		if( !init_checker("get_nevals()") )
			return 0;
		return nevals;
	}

	void gslIntegrator::set_epsabs(double e)
	/// set absolute error tolerance
	{ epsabs = e; }
	
	void gslIntegrator::set_epsrel(double e)
	/// set relative error tolerance
	{ epsrel = e; }
	
	double gslIntegrator::get_epsabs()
	/// set absolute error tolerance
	{ return epsabs; }
	
	double gslIntegrator::get_epsrel()
	/// set relative error tolerance
	{ return epsrel; }

	double gslIntegrator::integrate(double xmin, double xmax)	
	/// do the integration
	{
		if ( !init_checker("integrate()") )
			return 0;
		
    	gsl_integration_cquad(&F,xmin,xmax,epsabs,epsrel,w,&result,&error,&nevals);
		return result;
	}

	double gslIntegrator::integrate(double xmin, double xmax, void * p)	
	/// do the integration and provide params
	{
		F.params = p;
		got_pars = true;		
		if ( !init_checker("integrate()") )
			return 0;
		
    	gsl_integration_cquad(&F,xmin,xmax,epsabs,epsrel,w,&result,&error,&nevals);
		return result;
	}

	double gslIntegrator::operator()(double xmin, double xmax, void * p)
	{ return integrate(xmin,xmax,p); }	

	double gslIntegrator::operator()(double xmin, double xmax)
	{ return integrate(xmin,xmax); }	

//************************************************************
//		gslMatrix	
//************************************************************
	gslMatrix::gslMatrix()
	/// default constructor
	{
		m = gsl_matrix_alloc (1, 1);	// allocate space for a 1x1
		gsl_matrix_set (m, 1, 1, 0.0);	// set value to 0
		rows = cols = 1;
		loc[0] = loc[1] = 0;
		temp = 0.0;
	}

	gslMatrix::gslMatrix(int M, int N)
	/// construct empty MxN matrix
	{
		m = gsl_matrix_alloc(M,N);		// allocate space for MxN
		for(int i=0; i<M; i++) {		// set all values to zero
			for(int j=0; j<N; j++) 
				gsl_matrix_set(m,i,j,0.0);
		}
		rows = M;
		cols = N;
		loc[0] = loc[1] = 0;
		temp = 0.0;
	}

	gslMatrix::gslMatrix(const gslMatrix & n)
	/// copy constructor
	{
		rows = n.rows;
		cols = n.cols;				// copy rows and cols
		
		m = gsl_matrix_alloc(rows,cols);// allocate new space
		for(int i=0; i<rows; i++) {		// deep copy
			for(int j=0; j<cols; j++) 
				gsl_matrix_set(m,i,j,gsl_matrix_get(n.m,i,j));
		}
		loc[0] = loc[1] = 0;
	temp = gsl_matrix_get(n.m,0,0);		
	}

	gslMatrix::gslMatrix(const gsl_matrix * n, int i,int j)
	{
		rows = i;
		cols = j;				// copy rows and cols
		
		m = gsl_matrix_alloc(rows,cols);// allocate new space
		for(int i=0; i<rows; i++) {		// deep copy
			for(int j=0; j<cols; j++) 
				gsl_matrix_set(m,i,j,gsl_matrix_get(n,i,j));
		}
		loc[0] = loc[1] = 0;
		temp = gsl_matrix_get(m,0,0);		
	}



	gslMatrix::gslMatrix(double *X, int i, int j)
	/// constructor for double ** (inside of MatInit)
	{
		rows = i;
		cols = j;				// copy rows and cols
	
		fflush(stdout);
		m = gsl_matrix_alloc(rows,cols);// allocate new space

		for(int i=0; i<rows; i++) {		// deep copy
			for(int j=0; j<cols; j++) 
				gsl_matrix_set(m,i,j,X[i*cols+j]);
		}
		loc[0] = loc[1] = 0;
		temp = gsl_matrix_get(m,0,0);		
	}

	gslMatrix::~gslMatrix()
	/// destructor
	{ gsl_matrix_free (m); }

	void gslMatrix::reassign()
	/// reassign value to loc position
	{
	   	gsl_matrix_set(m,loc[0],loc[1],temp);
	}	

	double & gslMatrix::operator()(const int i, const int j)
	/// overload the () operator
	{ 
		if ( i<0 || i>=rows || j<0 || j>=cols )
		{
			std::cerr << "\nError: gslWrappers::gslMatrix::operator() -- Matrix indices out of bounds";
			temp = GSL_WRAP_ERR;
			return temp;
		}
		reassign();
		temp = gsl_matrix_get(m,i,j);
		loc[0] = i; loc[1] = j;	
		return temp;
	 }
	

	void gslMatrix::show(std::ostream & os)
	/// print out matrix with default ostream as cout
	{
		reassign();
		using std::endl;
		using std::setw;
		os << endl;
		for (int i=0; i<rows; i++)
		{
			os << "[";
			for (int j=0; j<cols; j++)
				os << setw(4) << gsl_matrix_get(m,i,j) << ", ";

			os << "\b\b ]"  << endl;
		}	
	}


	gslMatrix & gslMatrix::operator=(const gslMatrix & n)
	/// assignment operator
	{
		reassign();
		gsl_matrix_free (m);		// free old memory
		rows = n.rows;
		cols = n.cols;				// copy rows and cols
		
		m = gsl_matrix_alloc(rows,cols);// allocate new space
		for(int i=0; i<rows; i++) {		// deep copy
			for(int j=0; j<cols; j++) 
				gsl_matrix_set(m,i,j,gsl_matrix_get(n.m,i,j));
		}
		loc[0] = loc[1] = 0;
		temp = gsl_matrix_get(n.m,0,0);		
		return *this;
	}


	gslMatrix gslMatrix::operator+(const gslMatrix & n)
	/// addition operator
	{
		reassign();
		gslMatrix sum(*this);	// new sum matrix
		gsl_matrix_add(sum.m,n.m);	// gsl addition routine
		return sum;
	}
		
	gslMatrix gslMatrix::operator-()
	/// negation operator
	{
		reassign();
		gslMatrix neg(*this);	// new sum matrix
		for(int i=0; i<rows; i++) {		// deep copy
			for(int j=0; j<cols; j++) 
				gsl_matrix_set(neg.m,i,j,-gsl_matrix_get(neg.m,i,j));
		}
		return neg;
	}

	gslMatrix gslMatrix::operator-(const gslMatrix & n)
	/// subtraction operator
	{
		reassign();
		gslMatrix diff(*this);	// new sum matrix
		gsl_matrix_sub(diff.m,n.m);	// gsl addition routine
		return diff;	
	}
	
	gslMatrix gslMatrix::operator*(const double k)
	/// scalar multiplication
	{
		reassign();
		gslMatrix prod(*this);
		gsl_matrix_scale(prod.m,k);
		return prod;
	}

	gslMatrix gslMatrix::operator+(const double k)
	/// add a constant value to the whole matrix
	{
		reassign();
		gslMatrix sum(*this);
		gsl_matrix_add_constant(sum.m,k);
		return sum;
	}

	gslMatrix gslMatrix::operator&(const gslMatrix & n)
	/// elementwise multiplication -- like matlab .*
	{
		reassign();
		gslMatrix prod(*this);
		gsl_matrix_mul_elements(prod.m,n.m);
		return prod;
	}

	gslMatrix gslMatrix::operator|(const gslMatrix & n)
	/// elementwise division -- like matlab ./
	{
		reassign();
		gslMatrix quot(*this);
		gsl_matrix_div_elements(quot.m,n.m);
		return quot;
	}

	gslMatrix gslMatrix::operator*(const gslMatrix & n)
	/// matrix-matrix multiplication
	{
		reassign();
		if (cols != n.rows)
		{
			std::cerr << "\nError: gslWrappers::gslMatrix::operator* -- inner matrix dimensions must agree";
			gslMatrix err(1,1);
			err(1,1) = GSL_WRAP_ERR;
			return err;
		}
		gslMatrix prod(rows,n.cols);
		// call to the gsl interface to BLAS functions:
		// http://www.gnu.org/software/gsl/manual/html_node/Level-3-GSL-BLAS-Interface.html
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m,n.m,0.0,prod.m);
		return prod;
	}		

	gslMatrix gslMatrix::operator~()
	/// inversion operator
	{
		reassign();
		
		if(rows != cols)
		{
			std::cerr << "\nError: gslWrappers::gslMatrix::operator~ -- Matrix must be square";
			gslMatrix inv(rows,cols);
			return inv;
		}

		// the following code is weird, but it works.
		// my original attempt didn't work, so this is what I will use
		// similar to gsl example: http://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html	
		gsl_matrix *M = gsl_matrix_alloc( rows, cols);
		gsl_matrix_memcpy( M,m );
		gsl_matrix *inv1 = gsl_matrix_alloc(rows,cols);
		int s;
        gsl_permutation * p = gsl_permutation_alloc (rows);
     	gsl_linalg_LU_decomp (M, p, &s);
        gsl_linalg_LU_invert (M, p, inv1);
		gslMatrix inv(inv1,rows,cols);		// make new gslMatrix using inv1
		gsl_matrix_free(M);
		gsl_matrix_free(inv1);
		gsl_permutation_free(p);
		return inv;							// return inv
	}

	gslMatrix gslMatrix::inv()
	/// another inversion call
	{ reassign(); return ~(*this); }


	gslMatrix gslMatrix::operator*=( const gslMatrix & n)
	{
		reassign();
		(*this) = (*this)*n;
		return *this;
	}

	gslMatrix gslMatrix::operator+=( const gslMatrix & n)
	{
		reassign();
		(*this) = (*this)+n;
		return *this;
	}




	gslMatrix gslMatrix::operator^(int k)
	/// raise a matrix to a power
	{
		reassign();
		if(rows != cols)
		{
			std::cerr << "\nError: gslWrappers::gslMatrix::operator^ -- Matrix must be square";
			gslMatrix inv(1,1);
			return inv;
		}	
		
		if (k == 0)
			return eye();
		else if(k > 0)
		{
			gslMatrix prod(rows,cols);
			prod = (*this);
			for(int i=0;i<k-1;i++)
				prod *= (*this);
			return prod;
		}
		else
		{
			gslMatrix prod(rows,cols);
			gslMatrix inv(rows,cols);
			inv = ~(*this);
			prod = inv;
			for(int i=0;i<abs(k)-1;i++)
				prod *= inv;
			return prod;
		}	
	}

	gslMatrix gslMatrix::eye(int n)
	/// return an identity matrix of size n
	{
		reassign();
		if (n==0)
			n=rows;
		else
			n = ( n>0 ? n : 1);
		gslMatrix I(n,n);
		for(int i=0;i<n;i++)	
			I(i,i) = 1;
		return I;
	}


	gslMatrix operator*(const double k, gslMatrix & M )
	{ return M*k; }

	gslMatrix operator+(const double k, gslMatrix & M )
	{ return M+k; }



//************************************************************
//		gslCSpline
//************************************************************
	gslCSpline::gslCSpline(double * x, double * y, int N)
	/// constructor
	{
	// define spline working variables
		acc = gsl_interp_accel_alloc ();
		spline = gsl_spline_alloc (gsl_interp_cspline,N);  // allocate spline
	// initialize spline
		gsl_spline_init(spline,x,y,N);  
		initd = 1;
	 }


	gslCSpline::~gslCSpline()
	/// destructor
	{
		  gsl_spline_free (spline);
		  gsl_interp_accel_free(acc);
	}	

	double gslCSpline::operator()(double xi)
	/// overload () for interpolation calls
	{
		if (initd)
			return gsl_spline_eval(spline,xi,acc);
		else
		{
			std::cerr << "\nError: gslWrappers::gslCSpline::operator() -- spline not initialized\n";
			return GSL_WRAP_ERR;
		}
	}

	void gslCSpline::init_spline(double * x, double * y, int N)
	{
	// define spline working variables
		acc = gsl_interp_accel_alloc ();
		spline = gsl_spline_alloc (gsl_interp_cspline,N);  // allocate spline

	// initialize spline
		gsl_spline_init(spline,x,y,N);  
		initd = 1; 
	}


//************************************************************
//		Test Functions
//************************************************************
	
	void MatrixTest()
	{
		using std::cout;
		using std::endl;
		cout << "\nCreating empty matrix...";
		gslMatrix M(3,3);
		cout << "success\n";
		cout << "Testing (i,j) opertor:\n";
		cout << "\tM(2,2) = " << M(2,2) << endl;
		cout << "\tM(1,1) = " << M(1,1) << endl;
		cout << "\nTesting (i,j) on assignment on M(1,1)...";
		fflush(stdout);
		M(1,1) = 12.0;
		cout << "Success\n\tM(1,1) = " << M(1,1) << endl;

		cout << "\nTesting matrix show() routine:\n";
		M.show();

		cout << "\nCreating two new matrices for multiplication test...";
		double f[2][3] = {{1.0, 2.0, 3.0},{3.0,4.0,3.0}};
		double g[3][4] = {{1,2,3,4}, {2,3,4,5}, {4,5,6,7}};
		gslMatrix M1(f[0],2,3), M2(g[0],3,4);
		M1.show();
		M2.show();	
		cout << "Product M1*M2 = ";
		(M1*M2).show();	

		cout << "\nTesting identity .eye() function...";
		M1.eye().show();
		M2.eye(4).show();

		cout << "\nTesting inversion...";
		double h[3][3] = { {1,2,3},{4,5,6},{1,3,7} };
		//double h[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
		gslMatrix M3(h[0],3,3);
		M3.show();
		cout << "inv(M3):\n";
		M3.inv().show();
		cout << "\nM3 * inv(M3):";
		( M3 * ~M3).show();	

		cout << "\nTesting matrix power M^3:\n";
		double ff[3][3] = {{ 2,0,0},{0,2,0},{0,0,2}};
		gslMatrix M4(ff[0],3,3);
		(M4^-3).show();

		cout <<"\nTesting += and M+k:\n";
		M4+=(2*M4);
		M4.show();
	}

	double func(double x, void * params)
	{
		double a = *( (double *) params);
		return a*x;
	}
	
	double gunc(double x, void * params)
	{
		(void)params;	// params not used
		return x*cos(x)+exp(-x*x);
	}
	void IntTest()
	{
		
		using std::cout;
		using std::endl;
	
		cout << "Creating a gslIntegrator object with function func...";
		gslIntegrator inter(&func);
		cout << "Success\n";

		cout << "\nAttemtping integration of f(x)=x from x=[0..10]...";
		double a = 1;
		cout << "\n\tresult = " << inter.integrate(0,10,(void*)&a) << "\n";

		cout << "\nAttemtping integration of g(x)=x*cos(x)+exp(-x^2) from x=[0..10]...";
		inter.set_integrand(&gunc);
		inter.integrate(0,1);
		inter.set_params(NULL);
		cout << "\n\tresult = " << inter(0,10) << "\n";


	}	
	
	void rootTest()
	{
		using namespace std;
		cout << "Creating a gslRootFinder object with function func2 = x^2-3...";
		gslRootFinder finder(0,5,&func2);
		cout << "Success!\n\n";

		cout << "Finding root...";
		finder.findRoot();
		cout << "Success! root = " << finder.get_root() << endl;

	}

	double func2(double x, void * params)
	{
	 	(void)params;	
		return x*x - 3;
	}

}	// END NAMESPACE GSLWRAPPERS

