/*--------------------------------
** Mass Gap Research Simulation
** USU FAST Group
** Eric Addison
*--------------------------------*/ 

// This file contains implementation of the GravWave namespace, specfically 
// the Gwave class


#include "GravWave.h"

namespace GravWave {

using std::cerr;

GW_Bin_Insp::GW_Bin_Insp(double m, double n, double i, double R, double f,double T) : rooter(1.0,1e10,&fomega,1000,1e-6) 
/// Constructor with default arguments
{
	m1 = m;			// mass 1
	m2 = n;			// mass 2
	set_inc(i);		// inclination angle
	Rl = R;			// luminosity distance
	N = 1e5;		// default number of data points
	set_other_params();	// set other mass parameters
	phi_c = 0;		// phase at coalescence, phi(tc)
	flso = 1.0/(pow(6.0,1.5)*PI*M*MSCALE)*C*C*C/G;
	// set up root finder for calculating tc
	rooter.set_params((void*)this);
    set_f0(f);		// set initial frequency f0 and tc
	// set tmax from T or tc
	if ( T <= 0.0 || T > tc )
		tmax = tc;
	else
		tmax = T;
	dt = tmax/(double)N;	// time step for file output
};

GW_Bin_Insp::~GW_Bin_Insp()
/// desctructor (does nothing)
{ }

void GW_Bin_Insp::set_m1(double m)
/// set value of m1
{ m1 = m; set_other_params(); }


void GW_Bin_Insp::set_m2(double m)
/// set value of m2
{ m2 = m; set_other_params(); }


void GW_Bin_Insp::set_other_params()
/// method to compute all derived mass parameters
{
	M = m1+m2;					// total mass
	mu = m1*m2/M;				// reduced mass
	eta = mu/M;					// dimensionless mass ratio
	Mc = pow(eta,0.6)*M;		// chirp mass	
	dm = m1-m2;					// mass difference
	alpha = C*C*C/(G*M*MSCALE);	// constant needed for calculations
	calc_phi_params();			// calculate phi_params array
	calc_om_params();			// calculate om_params array
	C1 = 2*G*M*MSCALE*eta/(C*C*Rl*DSCALE); // value for h_t calculation
    set_Hp_params();			// calc Hp parameters
}


void GW_Bin_Insp::set_inc(double i)
/// method to set inclination and related variables
{
	inc = i;
	c = cos(i);
	s = sin(i);
}


void GW_Bin_Insp::set_f0(double f)
/// method to set the value of f0
{
	f0 = f;	// set value of f0
	calc_tc();	// calculate time to coal. tc
}


void GW_Bin_Insp::set_Rl(double r)
/// method to set value of Rl
{ Rl = r; }


void GW_Bin_Insp::set_N(double n)
/// set number of points in full output and time step
{ 
	N = n; 
	dt = tmax/N;
}


double GW_Bin_Insp::get_hp(double t)
/// calculate h_p waveform for time t
{

// if f0 has not been input
  if (f0 == 0.0)
  {
	cerr << "Error: GravWave::GW_Bin_Insp::calc_h_p(), minimum frequency " 
		<< "f0 not defined\n";
	return 0;
  }

// if t > tc
  if ( t > tc )
  {
	cerr << "Error: GravWave::GW_Bin_Insp::calc_h_p(), time value t exceeds "
		<< "time to coalescence tc\n";
 	return 0;
  }

// precompute some constant coefficients for wave calculations 
  double om_0=10*PI;

// phase related values  
  double Theta = (alpha/5.0)*eta*(tc-t);
  double T_8 = pow(Theta,0.125);

  phi_t = get_phi(T_8);
  om_t = get_omega(T_8);
  psi = phi_t - 2.0*om_t*log(om_t/om_0)/alpha;

// cosines and sines for the current psi
  cp[0] = cos(psi);
  cp[1] = cos(2*psi);
  cp[2] = cos(3*psi);
  cp[3] = cos(4*psi);
  cp[4] = cos(5*psi);
  cp[5] = cos(6*psi); 
  sp[0] = sin(psi);
  sp[2] = sin(3*psi);

// calculate the order coefficients of H_plus
  H_plus();

// expansion parameter x
  double x, x12;
  x = pow(om_t/alpha,0.666666667);; 
  x12 = sqrt(x);

// return value: h_p(t)
  return C1*(Hp[0]*x + Hp[1]*x*x12 + Hp[2]*x*x + Hp[3]*x*x*x12 + Hp[4]*x*x*x);

}

double GW_Bin_Insp::get_omega(double &T_8)
/// get current value omega(t)
// input is the time parameter THETA^(1/8)
{
// alpha = C^3/(G*m)
  double y;
  y = pow(T_8,-3.0);
  for( int ii=0; ii<3; ii++)
	y += om_params[ii]*pow(T_8,-(5.0+ii));
  
  return (alpha/8.0)*y;
}

double GW_Bin_Insp::get_phi(double &T_8)
/// get current value phi(t)
{
// dt = t_c - t

  double y;
  y = pow(T_8,5.0);
  for( int ii=0; ii<3; ii++)
	y += phi_params[ii]*pow(T_8,(3.0-ii));
  
  return phi_c - y/eta;

}

void GW_Bin_Insp::calc_om_params()
/// calculate the parameters in om_params
{  
  /* om_params = contants in omega(t) */
  om_params[0] = (743.0/2688.0 + 11.0/32.0*eta);
  om_params[1] = -3.0*PI/10.0;
  om_params[2] = 1855099.0/14450688.0 + 56975.0/258048.0*eta + 371.0/2048.0*eta*eta;

}

void GW_Bin_Insp::calc_phi_params()
/// Calculate the parameters in phi_params
{
  /* phi_params = constants in phi(t) */
  phi_params[0] = 3715.0/8064.0 + 55.0/96.0*eta;
  phi_params[1] = -3.0*PI/4.0;
  phi_params[2] = 9275495.0/14450688.0 + 284875.0/258048.0*eta + 1855.0/2048.0*eta*eta;

}


 

/* function to find time to coalescence by Brent method root finder */
double GW_Bin_Insp::calc_tc()
/// method to find time to coalescence tc
{

// if f0 has not been input
  if (f0 == 0.0)
  {
	cerr << "Error: GravWave::GW_Bin_Insp::calc_tc(), minimum frequency " 
		<< "f0 not defined\n";
	return 0;
  }

  if (f0 < 1.0)	// pseudo-monochromatic binary, use quarupole formula
  // using Peters formula
  {
	double beta = 64.0/5.0 * (G*G*G*m1*m2*M*pow(MSCALE,3.0))/(C*C*C*C*C);
	double a04 = pow(G*M*MSCALE/(4*PI*PI*f0*f0),1.33333333333);
	tc = a04/(4*beta);
	return tc;
  }

	tc = rooter.findRoot();
	return tc;

	//return tc;
}

void GW_Bin_Insp::H_plus()
/// calculate the Hplus coefficients for different orders
{
// p is double* of Hp_params
  Hp[0] = Hp_pars[0]*cp[1];
  Hp[1] = Hp_pars[1]*cp[0] + Hp_pars[2]*cp[2];
  Hp[2] = Hp_pars[3]*cp[1] + Hp_pars[4]*cp[3];
  Hp[3] = Hp_pars[5]*cp[0] + Hp_pars[6]*cp[2];
  Hp[3] += Hp_pars[7]*cp[4] + Hp_pars[8]*cp[1];
  Hp[4] = Hp_pars[9]*cp[1] + Hp_pars[10]*cp[3];
  Hp[4] += Hp_pars[11]*cp[5] + Hp_pars[12]*sp[0];
  Hp[4] += Hp_pars[13]*cp[0] + Hp_pars[14]*sp[2] + Hp_pars[15]*cp[2];

}

void GW_Bin_Insp::set_Hp_params()
/// set the parameters needed for Hplus calculation
{
  double c4 = c*c*c*c;
  double c6 = c*c*c*c*c*c;
  double dmp = s*dm/M;

  Hp_pars[0] = -(1+c*c);
  Hp_pars[1] = -dmp/8.0*(5.0+c*c);
  Hp_pars[2] = dmp/8.0*9.0*Hp_pars[0];
  Hp_pars[3] = 1.0/6.0*((19.0+9.0*c*c - 2.0*c4) - eta*(19.0 - 11.0*c*c - 6.0*c4));
  Hp_pars[4] = -4.0/3.0*s*s*Hp_pars[0]*(1-3.0*eta);
  Hp_pars[5] = dmp/192.0*((57.0+60.0*c*c - c4) - 2.0*eta*(49.0 - 12.0*c*c - c4));
  Hp_pars[6] = -dmp/192.0*27.0/2.0 * ( (73.0 + 40.0*c*c - 9.0*c4) - 2*eta*(25.0 - 8.0*c*c - 9.0*c4) );
  Hp_pars[7] = dmp/192.0*625.0/2.0 * (1 - 2.0*eta)*s*s*Hp_pars[0];
  Hp_pars[8] = -2*PI*Hp_pars[0];
  Hp_pars[9] = ((22.0 + 396.0*c*c + 145.0*c4 - 5.0*c6) + 5.0/3.0*eta*(706.0 - 216.0*c*c - 251.0*c4 + 15.0*c6) - 5.0*eta*eta*(98.0 - 108.0*c*c + 7.0*c4 + 5.0*c6))/120.0;
  Hp_pars[10] = 2.0/15.0*s*s*( (59.0+35.0*c*c - 8.0*c4) - 5.0/3.0*eta*(131.0 + 59.0*c*c - 24.0*c4) * 5.0*eta*eta*(21.0-3.0*c*c-8.0*c4) );
  Hp_pars[11] = -81.0/40.0*(1.0-5.0*eta + 5.0*eta*eta)*s*s*s*s*Hp_pars[0];
  Hp_pars[12] = dmp/40.0*(11.0+7.0*c*c + 10.0*(5.0+c*c)*log(2));
  Hp_pars[13] = -dmp/40.0*5*PI*(5.0+c*c);
  Hp_pars[14] = -dmp/40.0*27.0*(7.0-10.0*log(1.5))*Hp_pars[0];
  Hp_pars[15] = dmp/40.0*135.0*PI*Hp_pars[0];
}


// functions used for root finding with tc

double fomega(double x, void * params)
{
  if(x<=0)
	return -1;

  struct GW_Bin_Insp *gw = (GW_Bin_Insp*)params;
 
  
  double Theta = (gw->alpha/5.0)*gw->eta*x, T_8, y;
  T_8 = pow(Theta,0.125);

  y = pow(T_8,-3.0);

  for( int ii=0; ii<3; ii++)
	y += gw->om_params[ii]*pow(T_8,-(5.0+ii));
  
  return (gw->alpha/8.0)*y - 2*PI*gw->f0;

}

double dfomega(double x, void * params)
{
  if(x<=0)
	return -1;
  struct GW_Bin_Insp *gw = (GW_Bin_Insp*)params;

  double Theta = (gw->alpha/5.0)*gw->eta*x, T_8, y;

  T_8 = pow(Theta,0.125);

  y = (-3.0)*pow(T_8,-11.0);

  for( int ii=0; ii<3; ii++)
	y += -(5.0+ii)*gw->om_params[ii]*pow(T_8,-(13.0+ii));
  
  return (gw->alpha/8.0)*y/8.0;


}


void fdfomega(double x, void *params, double *f, double *df)
{
  if(x<=0)
  {
	*f = -1;
	*df = -1;
	return;
  }

  struct GW_Bin_Insp *gw = (GW_Bin_Insp*)params;

  double Theta = (gw->alpha/5.0)*gw->eta*x, T_8;

  T_8 = pow(Theta,0.125);

  *f = pow(T_8,-3.0);
  *df = (-3.0)*pow(T_8,-11.0);

  for( int ii=0; ii<3; ii++)
  {
  	 *f += gw->om_params[ii]*pow(T_8,-(5.0+ii));
	*df += -(5.0+ii)*gw->om_params[ii]*pow(T_8,-(13.0+ii));
  }
  *f *= gw->alpha/8.0;
  *f -= 2*PI*gw->f0;
  *df *= gw->alpha/8.0/8.0;

}

void GW_Bin_Insp::write_to_file(const char * filename)
/// write waveform to file for t = 0..tmax
{
	using namespace std;

	if (f0 == 0.0)
	{
	cerr << "Error: GravWave::GW_Bin_Insp::write_to_file(), minimum frequency " 
		<< "f0 not defined\n";
	return;
  	}

	char data[80] = "", meta[80] = "";
	strcat(data, filename); strcat(data, ".dat");
	strcat(meta, filename); strcat(meta, ".info");

	ofstream dat; dat.open(data);
	ofstream nfo; nfo.open(meta);

	// print meta data file
	nfo << "-------------------------------------------------" << endl;
	nfo << "*** Binary Inspiral Gravitational Wave Output ***" << endl;
	nfo << "	" << meta << endl;
	nfo << "-------------------------------------------------" << endl;
	nfo << "\n\nm1 = " << m1 << " Msun" <<  endl;
	nfo << "m2 = " << m2 << " Msun " << endl;
	nfo << "inclination = " << inc*180.0/PI << " degrees" << endl;
	nfo << "Lum. Dist. = " << Rl << " Mpc" << endl;
	nfo << "Time to Coal. tc = " << tc << " sec" << endl;
	nfo << "tmax = " << tmax << " s" << endl;
	nfo << "f0 = " << f0 << " Hz" << endl;
	nfo << "dt = " << dt << " s" << endl;
	nfo << "\n\n\nGenerated by GravWave::GW_Bin_Insp::write_to_file()";

	// write wave data
	for(int i=0; i<N; i++)
		dat << i*dt << ",   " << get_hp(i*dt) << endl;
 	
	dat.close();
	nfo.close();	
	
}	
	
// test function for GW_Bin_Insp class
	void GW_Test()
	{
		using std::cout;
		using std::endl;
		using GravWave::GW_Bin_Insp;

		double m1 = 1.4, m2 = 1.4, Rl = 10.0, inc = 45.0*PI/180.0, f0=0.1;
		
		cout << "Creating a GW_Bin_Insp object with:\n" 
			<< "\tm1 = " << m1 << endl
			<< "\tm2 = " << m2 << endl
			<< "\tRl = " << Rl << endl
			<< "\tf0 = " << f0 << endl
			<< "\tinc = " << inc << endl;

		double tmax = PI*1e7; // one year obs time
		GW_Bin_Insp gw(m1,m2,inc,Rl,f0,tmax);

		cout << "tmax: " << gw.get_tmax() << endl;
		cout << "tc: " << gw.get_tc() << endl;	
		cout << "flso: " << gw.get_flso() << endl;
		cout << "\nWriting waveform to file ...";
		gw.write_to_file();
		cout << " Success\n";
	}
} // END namespace GravWave

