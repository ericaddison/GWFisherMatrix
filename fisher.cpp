/*--------------------------------
** Mass Gap Research Simulation
** USU FAST Group
** Eric Addison
*--------------------------------*/ 

/* This file will calculate the Fisher infomation matrix for the gravitational wave signal of an inspiralling compact binary. */

/* Following the Fisher construction from the paper:
   Parameter estimation of inspiralling compact binaries using 3.5 post-Newtonian gravitational wave phasing: The nonspinning case
   http://prd.aps.org/abstract/PRD/v71/i8/e084008
*/

#include "fisher.h"

#define MAX(x,y) (x>y?x:y)
#define MIN(x,y) (x<y?x:y)
namespace GravWave
{


Fisher::Fisher(GW_Bin_Insp & gw, noise_curve & N) : flso(gw.flso), eta(gw.eta), Mc(gw.Mc), tc(gw.tc), F(4,4) , cov(4,4), Sn(N.get_spline()), GW(gw)
/// constructor with initializer list
// copies required data from noise curve and binary objects
{

	Mc *= MSCALE*G/(C*C*C);
	A = ASCALE*pow(Mc,5.0/6.0)/(gw.get_Rl());	// FIX THIS
// parameter set up 
	vlso = 1.0/sqrt(6.0);	// v for last stable orbit	
	v_base = pow(PI*Mc*pow(eta,-0.6)*VSCALE,0.333333333);
	f0 = 215;
	//f0 = 1;
	S0 = 1e-49;
	fmin = gw.f0; 
	double tmax = gw.tmax;		// get tmax fro gw object
	//fmin = (N.get_f())[0];	// lowest freq value in noise curve
	std::cout << "tmax = " << tmax 
		<< "\ntc = " << tc <<  std::endl;
	fmax = set_fmax(tmax);
}

Fisher::~Fisher()
{}

// Main fisher function
void Fisher::calc_fisher()
{
	using gslWrappers::gslIntegrator;
	
	gslIntegrator int_gij;	// integrator for gij
	gij_params * gij_pars = new gij_params;	// parameter struct to hold stuff
	gij_pars->funcs[0] = &g_tc;
	gij_pars->funcs[1] = &g_phic;
	gij_pars->funcs[2] = &g_Mc;
	gij_pars->funcs[3] = &g_eta;
	gij_pars->F = this;

	/* Compute Fisher matrix entries */
	int_gij.set_integrand(gij);	// always integrating function gij: only 
								// parameters change
	double temp;
	for(int ii=0;ii<4;ii++)
	{
	  gij_pars->i=ii;
	  for(int jj=0;jj<=ii;jj++)
	  {
	   	gij_pars->j=jj;
	   	int_gij.set_params( (void*)gij_pars );
	   	temp = 4*A*A*int_gij(fmin,fmax);
       	F(ii,jj) = temp;
		F(jj,ii) = temp;
	  }
	}

// SNR calc (REMOVE LATER)
	std::cout << "SNR = " << calc_SNR() << "\n";

	delete gij_pars;
	return;
}


double Fisher::calc_SNR()
{ 
	using gslWrappers::gslIntegrator;
	gslIntegrator int_SNR = gslIntegrator(&SNR_func,(void*)this); // initialize SNR integrator
	return 2*A*sqrt(int_SNR(fmin,fmax)); 
}

double SNR_func(double f, void * params)
{
    Fisher * p = (Fisher *)params;
	return pow(f,-7.0/3.0) / (p->Sn(f));
	//return pow(f,-7.0/3.0) / (p->Sn_analytic(f));
}

double Fisher::Sn_analytic(double f)
{
		double x = f/f0;	
		double r =(pow(x,-4.14) - 5.0*pow(x,-2.0));
		r += 111.0*(1.0 - x*x + x*x*x*x/2.0)/(1.0 + x*x/2.0);
		return S0*r;
}

double gij(double f, void * params)
{

  double y;
// recast params pointer as gij_params pointer
  gij_params * p = (gij_params*) params;

// calculate y = g_i(f)*g_j(f)
  y = (*(p->funcs[p->i]))(f,p->F)*(*(p->funcs[p->j]))(f,p->F);

// return value: f^(-7/3)*g_i*g_j/Sn
  return pow(f,-7.0/3.0)*y / (p->F->Sn(f));
//  return pow(f,-7.0/3.0)*y / (p->F->Sn_analytic(f));
}


double g_tc(double f, void * params)
{
  	//Fisher * fish = (Fisher *)params;
  	(void)params;
	return 2.0*PI*(f);	
}

double g_phic(double f, void * params)
{
  (void)f;
  (void)params;		// f and params not used -- turn off warnings
  return -1.0;
}

double g_Mc(double f, void * params)
{

  Fisher * fish = (Fisher *)params;

  double v = fish->v_base*pow(f,0.33333333);	// p[3] is v_base
 
// get required values for function evaluation
  double dv = fish->dvdMc(v);
  fish->calc_dadMc();
  fish->calc_a(v); 
  double Y = 0;

// first summation in g_Mc
  for(int ii=0; ii<=7; ii++)
    Y += fish->a[ii]*(ii-5)*dv*pow(v,ii-6);
 
  Y += fish->dadMc[0] + v*fish->dadMc[1];	// add on terms with dalpha/dMc
  Y *= 3.0*fish->Mc/(128.0*fish->eta);

 return Y;  
}

double g_eta(double f, void * params)
{

  Fisher * fish = (Fisher *)params;

  double v = fish->v_base*pow(f,0.33333333);	// p[3] is v_base


// get required values for function evaluation
  //std::cout << "g_eta: fish = " << fish << "\n";
  double dv = fish->dvdeta(v);
  double Y = 0;
  double etai = 1.0/fish->eta;
  double vp;

// recalculate a and da for current v
  fish->calc_dadeta(v);
  fish->calc_a(v); 

// perform summation to find g_ln(eta)
  for(int ii=0; ii<=7; ii++)
  {
    vp = pow(v,ii-5);
    Y += fish->a[ii]*vp*(-etai + dv*(ii-5)/v);
    Y += vp*fish->dadeta[ii];
  }
  
  Y *= 3.0/128.0;

 return Y;  
}


double Fisher::set_fmax(double tmax)
/// calculate max frequency for Fisher calculation
{
	if( (tmax != 0 && tmax < tc) )	// if want to observe to some time less than tc
	{	
		double Theta = (GW.alpha/5.0)*GW.eta*(tc-tmax), T_8;
		T_8 = pow(Theta,0.125);
		fmax = GW.get_omega(T_8)/(2.0*PI);
	}
	else
		fmax = flso;
	return fmax<flso?fmax:flso;
}


	double Fisher::dvdeta(double v)
	/// Deriv of v wrt eta
	{
	  return -0.2*v/eta;
	}

	double Fisher::dvdMc(double v)
	/// deriv of v wrt Mc
	{
	  return 1.0/3.0*(v/Mc);
	}

	void Fisher::calc_dadMc()
	/// derivs of alpha wrt Mc (only applies for alpha_5 and alpha_6
	{
		dadMc[0] = FC[7]/Mc/3.0+FC[8]*eta/Mc;
		dadMc[1] = FC[13]/Mc/3.0;
	}


	void Fisher::calc_a(double v)
	/// calculate the a (alpha) coefficients
	{
		a[0] = 1;
		a[1] = 0;
		a[2] = FC[0] + FC[1]*eta;
		a[3] = FC[2];
		a[4] = FC[3]+FC[4]*eta+FC[5]*eta*eta;
		a[5] = FC[6]+FC[7]*log(v/vlso)+FC[8]*eta*(1.0+3.0*log(v/vlso));
		a[6] = FC[9]+eta*FC[10]+FC[11]*eta*eta+FC[12]*eta*eta*eta+FC[13]*log(4.0*v);
		a[7] = FC[14]+FC[15]*eta+FC[16]*eta*eta;
	}

	void Fisher::calc_dadeta(double v)
	/// calculate derivatives of alpha wrt eta
	{
		dadeta[0] = 0;
		dadeta[1] = 0;
		dadeta[2] = FC[1];
		dadeta[3] = 0;
		dadeta[4] = FC[4] + 2*FC[5]*eta;
		dadeta[5] = -FC[7]/eta/5.0+FC[8]*(1.0+3.0*log(v/vlso))-3.0/5.0*FC[8];
		dadeta[6] = FC[10]+2.0*FC[11]*eta+3.0*FC[12]*eta*eta-FC[13]/eta/5.0;
		dadeta[7] = FC[15]+2.0*FC[16]*eta;
	}




// Defining static member variables
	  double Fisher::gamma = 0.5772156649;	   // Euler's constant
	  double Fisher::theta = -11831.0/9240.0;  // constant from Arun et al paper
	  double Fisher::lam = -1987.0/3080.0;     // another constant
	  //double Fisher::VSCALE = (C*C*C) / (MSCALE * G);
	  double Fisher::VSCALE = 1;
	  //double Fisher::ASCALE = 2.0*sqrt(5.0/96.0)*pow(G*MSCALE,5.0/6.0) / ( pow(PI,2.0/3.0)*pow(C,1.5)*DSCALE);
	  double Fisher::ASCALE = 2.0*sqrt(5.0/96.0)/ ( pow(PI,2.0/3.0)*DSCALE) *C;
	double Fisher::FC[17] = {
		// Coefficients needed to compute the alpha values
		// alpha 2
		  20.0/9.0 * 743.0/336.0,
		  20.0/9.0 * 11.0/4.0,

		// alpha 3
		  -16.0*PI,

		// alpha 4
		  10.0 * 3058673.0/1016064.0,
		  10.0 * 5429.0/1008.0,
		  10.0 * 617.0/144.0,

		// alpha 5
		  PI * 38645.0/756.0,
		  PI * 38645.0/252.0,
		  -PI * 65.0/9.0,

		// alpha 6
		  11583231236531.0/4694215680.0 - 640.0*PI*PI/3.0 - 6848.0*Fisher::gamma/21.0,
		  -15335597827.0/3048192.0 + 2255.0*PI*PI/12.0 - 1760.0*Fisher::theta/3.0 + 12320.0*Fisher::lam/9.0,
		  76055.0/1728.0,
		  -127825.0/1296.0,
		  -6848.0/21.0,

		// alpha 7
		  PI * 77096675.0/254016.0,
		  PI * 378515.0/1512.0,
		  -PI * 74045.0/756.0
	};

// NOISE CURVE CLASS METHODS

	noise_curve::noise_curve(const char * file)
	/// noise_curve constructor
	{
		using namespace std;	
		ifstream input;
		input.open(file);

		strncpy(filename,file,80);		
		char temp[80];
		int ii,line_cnt=0;

	// find how many lines are in the file
		while( true )
		{
			input.getline(temp,80);
			if(input.eof())
				break;
			line_cnt++;
		}

	// allocate memory
		f = new double[line_cnt];
		S = new double[line_cnt];
		N = line_cnt;

	// rewind file and read in data
		input.clear();		// clear eof bit (and other errors)
		input.seekg(0);		// put get pointer at beginning of file
		
	  for(ii=0; ii<N; ii++)
	  {
		input.getline(temp,80);	// get line of data
		// strtok string token function 
		f[ii] = atof(strtok(temp,",\t ")); // store frequency point
		S[ii] = pow(atof(strtok(NULL,",\t ")),2); // store noise curve value
	  }

		Spline.init_spline(f,S,N);
	}

	noise_curve::~noise_curve()
	/// destructor to free memory
	{
		delete [] f;
		delete [] S;
	}

	gslCSpline & noise_curve::get_spline()
	{ return Spline; }


void Ftest()
{
		using std::cout;
		using std::endl;
		using GravWave::GW_Bin_Insp;
		using gslWrappers::gslCSpline;
		using gslWrappers::gslMatrix;

		double m1, m2, Rl , inc = 0.0*PI/180.0, f0=9;
	// Test values for matching Arun et al (2005) paper	
	// BHBH SNR 10 for analytic noise curve
		//Rl = 8.0*218.7; m1 = m2 = 10;
	// BHNS SNR 10 for analytic noise curve
		//Rl = 8.0*92; m1 = 10; m2 = 1.4;
	// NSNS SNR 10 for analytic noise curve
		Rl = 8.0*44.2; m1 = m2 = 1.4;
		cout << "Creating a GW_Bin_Insp object with:\n" 
			<< "\tm1 = " << m1 << endl
			<< "\tm2 = " << m2 << endl
			<< "\tRl = " << Rl << endl
			<< "\tf0 = " << f0 << endl
			<< "\tinc = " << inc << endl;

		GW_Bin_Insp gw(m1,m2,inc,Rl,f0,-1);
		cout << "Success\n";
	
		cout << "\nCreating a noise_curve object with file \"noise.dat\"...";
		noise_curve Sn("noise.dat");
		cout << "Success\n";
		fflush(stdout);

		cout << "\nTesting noise curve values:\n";
		cout << "\tf(0) = " << (Sn.get_f())[0];
		cout << "\n\tS(0) = " << (Sn.get_S())[0] << endl;
		fflush(stdout);

		cout << "\nTesting spline:\n";
		cout << "\tSn(10) = " << Sn(10) << endl;
		
		cout << "\nCreating a Fisher Object ... ";
		Fisher fish(gw,Sn);
		cout << "Success\n";
		cout << "fmax = " << fish.get_fmax() << "\n\n";
		cout << "Calculating Fisher Matrix ... ";
		fish.calc_fisher();
		cout << " Success!\n";
		fish.show();
		fish.invert();		
	    fish.show_cov();
		cout << endl;	
		cout << "found Dtc = " << 1000*sqrt(fish.get_cov(0,0)) << " msec"<< endl;
		cout << "found Dphic = " << sqrt(fish.get_cov(1,1)) << endl;
		cout << "found DMc/Mc = " << 100*sqrt(fish.get_cov(2,2)) << "\%" << endl;
		cout << "found Deta/eta = " << 100*sqrt(fish.get_cov(3,3)) << "\%";
		cout << endl << endl;
	}



} 	// END NAMESPACE GRAVWAVE
