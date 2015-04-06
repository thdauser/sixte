// This code generates a continuous stream of x-ray pulses given an
// input file containing an array of arrival times and photon
// energies. Depending upon the decimation factor only 1 in x samples
// of the data stream are written to an array. The code uses a 4th
// order runge-kutta method to numerically integrate the system on
// differential equations which describes the TES electro-thermal
// response. Tn this code we use a simple linear transition shape with
// alpha (TES temperature sensitivity) and beta (TES current
// sensitivity) dependence. The non-linearity in the TES response due
// to the dissproportionality between resistance and current is
// intrinsically included. Similarly the noise will be non-stationary.
//
// Code by J. Wilms based on the IDL-program pulses_NL.pr by
// Stephen J. Smith - NASA GSFC - 08/APRIL/2009
//
//
// current; muA
// J/K -> pJ/K
// power flow: pW/K   (differential conductance)
//   specify the K and n, spit out G
//   practical use: G is preferred
//      -> allow either G or K, print out both
//  add flag to invert signal, do mapping with imin,imax only
//
// ASTRONOMERS BEWARE: This code uses MKS units!
//
//
// Mapping to I&H
// - rename T_start to T0_start
// - use T0 instead of T1 for current temperature

//
// TODO:
// * write meta interfaces to simulation code
//     - work on photon event list
// * build interface based on record format (tesrecord.h)
// * write current to file
// * improve TES diagnostics
// * calculate taup/taum using proper complex arithmetic!
// * build in material database (for n et al., see I&H, table 2)
// * take into account proper time of arrival of event
// * use stochastic equation solver
// * implement saturation regime
// * generate bitstreams of varying bit length (rather than assuming 16bit)
// * progress bar
//
//

#include "tessim.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fitsio.h>

#include <assert.h>

gsl_rng *rng;

const double kBoltz=1.3806488e-23; // Boltzmann constant [m^2 kg/s^2]
const double eV=1.602176565e-19 ;  // 1eV [J]
const double keV=1.602176565e-16 ;  // 1keV [J]


// Differential equations for the TES to solve. This is the simplest
// calorimeter model assuming a single heat capacity for the 
// TES and absorber

// time is not used in TES_Tdifferential_NL and TES_JAC. 
// Switch off the related warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
int TES_Tdifferential_NL(double time, const double Y[], double eqn[], void *params) {

  // index 0: electrical equation I0=Y[0]
  // index 1: thermal equation    T1=Y[1]

  double II=Y[0];
  double TT=Y[1];

  tesparams *tes=(tesparams *) params;

  // RT(I0,T1)
  double RT=tes->R0+tes->dRdT*(TT-tes->T_start)+tes->dRdI*(II-tes->I0_start);

  // Electrical circuit equation: dI/dt -- dy_0/dt=f_0(t,y_0(t),y_1(t))
  // note: Vdn, Vexc, Vcn are the Johnson noise terms, Vunk is the unknown noise
  eqn[0]=(tes->V0-II*(tes->RL+RT)+tes->Vdn+tes->Vexc+tes->Vcn+tes->Vunk)/tes->Lin;

  // Thermal equation: dT/dt -- dy_1/dt=f_1(t,y_0(t),y_1(t))
  eqn[1]=(II*II*RT-tes->Pb1
	  -II*(tes->Vdn+tes->Vexc+tes->Vunk)
	  +tes->Pnb1+tes->En1)/tes->Ce1;

  return GSL_SUCCESS;
}

int TES_jac (double time, const double Y[], double *dfdy, double dfdt[], void *params) {

  tesparams *tes=(tesparams *) params;

  const double II=Y[0];
  const double TT=Y[1];

  double RT=tes->R0+tes->dRdT*(TT-tes->T_start)+tes->dRdI*(II-tes->I0_start);

  // J_00 = df_0(t,y)/dy0
  double J00=(tes->RL+RT+II*tes->dRdI)/tes->Lin;

  // J_01 = df_0(t,y)/dy1
  double J01=-tes->dRdT*II/tes->Lin;

  // J_10 = df_1(t,y)/dy0
  double J10=(2.0*I*RT-II*II*tes->dRdT-(tes->Vdn+tes->Vexc))/tes->Ce1;

  // J_11 = dt_1(t,y)/dy1
  double J11=II*II*tes->dRdT/tes->Ce1;

  // setup Jacobian matrix per example at 
  // https://www.gnu.org/software/gsl/manual/html_node/ODE-Example-programs.html#ODE-Example-programs
  gsl_matrix_view dfdy_mat=gsl_matrix_view_array(dfdy,2,2);
  gsl_matrix *m=&dfdy_mat.matrix;

  gsl_matrix_set(m,0,0,J00);
  gsl_matrix_set(m,0,1,J01);
  gsl_matrix_set(m,1,0,J10);
  gsl_matrix_set(m,1,1,J11);

  dfdt[0]=0.;
  dfdt[1]=0.;

  return GSL_SUCCESS;
}
#pragma GCC diagnostic pop


double tnoi(tesparams *tes) {
  // Thermal (phonon) noise terms
 
  double n1=tes->n-1.;
  double G1=tes->Gb1*pow(tes->T1/tes->T_start,n1);
  
  double gamma;

  // Irwin, eq. 93
  if (tes->mech==0) {
    gamma=(pow(tes->Tb/tes->T1,n1+2.0)+1.0)/2.0;
  } else {
    gamma=1.;
  }
 
  return(gsl_ran_gaussian(rng, sqrt(4*kBoltz*tes->T1*tes->T1*G1*gamma*tes->bandwidth) ));
}

double tpow(tesparams *tes) {
  // Calculate the power flow between elements of the pixel
  //

  double beta_th=tes->n-1.0;
  double G1=tes->Gb1*pow(tes->T1/tes->T_start,beta_th);

  return G1*(pow(tes->T1,tes->n)-pow(tes->Tb,tes->n))/(tes->n*pow(tes->T1,beta_th));
}

//
// Add properties of tes to the header of fptr. 
// 
// Note: we write all information using SI units (so no conversion to mA etc.!)
//
void tes_fits_write_params(fitsfile *fptr, tesparams *tes,int *status) {
  // bail out if status is not ok
  CHECK_STATUS_VOID(*status);

  fits_update_key(fptr,TDOUBLE,"DELTAT",&tes->delta_t,"Integration step size [s]",status);
  fits_update_key(fptr,TDOUBLE,"IMIN",&tes->imin,"Current corresponding to 0 ADU [A]",status);
  fits_update_key(fptr,TDOUBLE,"IMAX",&tes->imax,"Current corresponding to 65534 ADU [A]",status);
  double dummy=1.0/tes->aducnv;
  fits_update_key(fptr,TDOUBLE,"ADUCNV",&dummy,"ADU conversion factor [A/ADU]",status);
  fits_update_key(fptr,TDOUBLE,"I0_START",&tes->I0_start,"Initial bias current [A]",status);
  fits_update_key(fptr,TDOUBLE,"R0",&tes->R0,"Operating point resistance [Ohm]",status);
  fits_update_key(fptr,TDOUBLE,"RL",&tes->RL,"Shunt/load resistor value [Ohm]",status);
  fits_update_key(fptr,TDOUBLE,"LIN",&tes->Lin,"Circuit inductance [H]",status);
  fits_update_key(fptr,TDOUBLE,"ALPHA",&tes->alpha,"TES sensitivity T/R*dR/dT (alpha)",status);
  fits_update_key(fptr,TDOUBLE,"BETA",&tes->beta,"TES current dependence I/R*dR/dI (beta)",status);
  fits_update_key(fptr,TDOUBLE,"T_START",&tes->T_start,"Initial operating temperature [K]",status);
  fits_update_key(fptr,TDOUBLE,"TB",&tes->Tb,"Heat sink temperature [K]",status);
  fits_update_key(fptr,TDOUBLE,"N",&tes->n,"Heat sink coupling parameter n",status);
  fits_update_key(fptr,TDOUBLE,"CE1",&tes->Ce1,"Absorber+TES heat capacity at Tc",status);
  fits_update_key(fptr,TDOUBLE,"PB1",&tes->Pb1,"Thermal power flow",status);
  fits_update_key(fptr,TDOUBLE,"GB1",&tes->Gb1,"Heat link thermal conductance at Tc",status);
  fits_update_key(fptr,TLONG,"SEED",&tes->seed,"Seed of random number generator",status);
  if (tes->simnoise) {
    fits_update_key(fptr,TINT,"SIMNOISE",&tes->simnoise,"Simulating noise terms",status);
  } else {
    fits_update_key(fptr,TINT,"SIMNOISE",&tes->simnoise,"Not simulating noise terms",status);
  }
  fits_update_key(fptr,TDOUBLE,"M_UNKNOWN",&tes->m_unknown,"Magnitude of unknown noise",status);
}


void tes_print_params(tesparams *tes) {
  //
  // print TES properties to console. If something changes here, do not forget
  // to also edit tes_fits_write_params!
  
  headas_chat(0,"\nStatus of TES Pixel with ID: %s\n",tes->ID);
  headas_chat(0,"\n");
  headas_chat(0,"Start time of simulation [s]            : %15.6f\n",tes->tstart);
  headas_chat(0,"Current time of simulation [s]          : %15.6f\n",tes->time);
  headas_chat(0,"Sample rate [Hz]                        : %15.6f\n",tes->sample_rate);
  headas_chat(0,"Integration step size [mus]             : %15.6f\n",1e6*tes->delta_t);
  headas_chat(0,"Current corresponding to 0 ADU [muA]    : %15.1f\n",1e6*tes->imin);
  headas_chat(0,"Current corresponding to 65534 ADU [muA]: %15.1f\n",1e6*tes->imax);
  headas_chat(0,"ADU to current conv. factor [muA/ADU]   : %15.8e\n",1e6/(tes->aducnv));
  headas_chat(0,"\n");

  headas_chat(0,"Initial bias current [muA]              : %10.5f\n",1e6*tes->I0_start);
  headas_chat(0,"Operating point resistance R0 [mOhm]    : %10.5f\n",1000.*tes->R0);
  headas_chat(0,"Shunt/load resistor value RL [mOhm]     : %10.5f\n",1000.*tes->RL);
  headas_chat(0,"Circuit inductance [nH]                 : %10.5f\n",1e9*tes->Lin);
  headas_chat(0,"TES sensitivity T/R*dR/dT (alpha)       : %10.3f\n",tes->alpha);
  headas_chat(0,"TES current dependence I/R*dR/dI (beta) : %10.3f\n",tes->beta);
  headas_chat(0,"\n");

  headas_chat(0,"Initial operating temperature [mK]      : %10.3f\n",1000.*tes->T_start);
  headas_chat(0,"Heat sink temperature [mK]              : %10.3f\n",1000.*tes->Tb);
  headas_chat(0,"Heat sink coupling parameter n          : %10.2f\n",tes->n);
  headas_chat(0,"\n");

  // CHECK UNITS!
  headas_chat(0,"Absorber+TES heat capacity at Tc [pJ/K] : %10.5f\n",1e12*tes->Ce1);
  headas_chat(0,"Thermal power flow  [pW/K]              : %10.5f\n",1e12*tes->Pb1);
  headas_chat(0,"Heat link thermal conductance at Tc     : %10.5e\n",tes->Gb1);
  headas_chat(0,"\n");

  headas_chat(0,"Effective bias voltage [muV]            : %10.5f\n",1e6*tes->V0);
  headas_chat(0,"Current [muA]                           : %10.5f\n",1e6*tes->I0);
  headas_chat(0,"Temperature [mK]                        : %10.5f\n",1000.*tes->T1);
  headas_chat(0,"Current Resistivity [mOhm]              : %10.5f\n",1000.*tes->RT);
  headas_chat(0,"\n");
    
  headas_chat(0,"Seed of random number generator         : %10lu\n",tes->seed);
  if (tes->simnoise) {
    headas_chat(0,"Simulating noise terms\n");
  } else {
    headas_chat(0,"NOT simulating noise terms\n");
  }

  headas_chat(0,"\n**********************************************************\n");
  headas_chat(0,"Derived Properties of the TES\n");
  headas_chat(0,"\n");

  //
  // The following equations are from Irwin & Hilton, Table 1
  //
  // NOTE: This part was typed jetlagging on a plane somewhere over Greenland.
  // I have not yet checked the equations for correctness
  //

  // Biasing of the detector (Smith, thesis, p.19)
  // NOTE: the |0.3| is arbitrary
  double biaspar=(tes->R0 - tes->RL)/(tes->R0+tes->RL);
  headas_chat(0,"TES bias parameter (R0-RL)/(R0+RL)     :  %15.5f\n",biaspar);
  if (biaspar > 0.3 ) {
    headas_chat(0,"       --> detector is voltage biased\n");
  } else {
    // RL>>R0: current bias case
    if (biaspar < -0.3 ) {
      headas_chat(0,"     --> detector is current biased\n");
    } else {
      headas_chat(0,"WARNING: RL approx R0: expect small/no electrothermal feedback\n");
    }
  }


  // joule power
  double Pj=tes->I0_start*tes->I0_start*tes->R0; 

  // thermal conductance (I&H, after eq. 6)
  double G=tes->Gb1;  // CHECK!

  // loop gain
  double ell=Pj*tes->alpha/(G*tes->T_start);

  headas_chat(0,"Joule power  [J]                       : %10.5e\n",1000.*Pj);
  headas_chat(0,"Loop gain                              : %10.5f\n",ell);

  // heat capacity
  double C=tes->Ce1;

  //
  // time scales
  //

  // natural
  double tau=C/G;

  // constant current
  double tau_I=tau/(1.-ell);

  // zero-inductance effective thermal
  double RLR0=tes->RL/tes->R0;
  double tau_eff=tau*(1.+tes->beta+RLR0)/(1.+tes->beta+RLR0+(1.-RLR0)*ell);

  // electrical
  double tau_el=tes->Lin/(tes->RL+tes->R0*(1.+tes->beta));

  headas_chat(0,"Characteristic timescales (all in ms)\n");
  headas_chat(0,"    natural                            : %10.5f\n",1000.*tau);
  headas_chat(0,"    constant current                   : %10.5f\n",1000.*tau_I);
  headas_chat(0,"    zero-inductance effective thermal  : %10.5f\n",1000.*tau_eff);
  headas_chat(0,"    electrical                         : %10.5f\n",1000.*tau_el);

  if (tau_I<0.) {
    headas_chat(0,"WARNING: constant current tau is negative!\n");
    headas_chat(0,"         The TES might be instable to thermal runaway\n");
    headas_chat(0,"         (but see below)\n");
  }

  // pulse shape timescales

  // NOTE: the following needs to be checked a little bit more, since even
  // for negative tsqrt we might be able to obtain a stable solution!
  double tsqrt=pow((1./tau_el - 1./tau_I),2.)-4.*tes->R0*ell*(2.+tes->beta)/(tau*tes->Lin);
  if (tsqrt<0.) {
    headas_chat(0,"WARNING: Detector is underdamped\n");
    headas_chat(0,"PROCEED WITH CAUTION AND HOPE JOERN IMPROVES THE DIAGNOSTICS!\n");
  } 

  //TODO: we need to be more careful in the next equations and do the
  //      calculation using complex numbers!
  double taum=1./tau_el + 1./tau_I;
  double taup=2./(taum+tsqrt);
  taum=2./(taum-tsqrt);

  headas_chat(0,"Pulse timescales (in ms)\n");
  headas_chat(0,"    rise                               : %10.5f\n",1000.*taup);
  headas_chat(0,"    fall                               : %10.5f\n",1000.*taum);

  //
  // TO DO: estimate pulse heights

  //
  // stability criteria
  //
  double Lsum=ell*(3.+tes->beta-RLR0)+1.+tes->beta+RLR0;
  double Lsqrt=2.*sqrt(ell*(2.+tes->beta)*(ell*(1.-RLR0)+1.+tes->beta+RLR0));
  double Lfac=tes->R0*tau/((ell-1.)*(ell-1.));

  double Lcritp=(Lsum+Lsqrt)*Lfac;
  double Lcritm=(Lsum-Lsqrt)*Lfac;
  
  headas_chat(0,"Inductance at critical damping:\n");
  headas_chat(0,"    L_plus [H]                         : %10.5e\n",Lcritp);
  headas_chat(0,"    L_minus[H]                         : %10.5e\n",Lcritm);

  if (tsqrt<0.) {
    // criteria for underdamped detector
    if (tau>(ell-1.)*tau_el) {
      headas_chat(0,"tau>(ell-1)*tau_el: TES is underdamped, but stable!\n");
    } else {
      headas_chat(0,"WARNING: TES is underdamped and NOT stable\n");
      double Lthresh=tau*(tes->RL+tes->R0*(1.+tes->beta))/(ell-1.0);
      headas_chat(0,"         To be stable, circuit inductance needs to be < %15.5e",Lthresh);
    }
  }

#ifndef SHOW_NOT_TRUSTWORTHY
  if (taup<taum) {
    headas_chat(0,"Tau(rise)<Tau(fall): System is overdamped\n");
  }
  if (fabs(taup/taum)-1. < 1e-5) {
    headas_chat(0,"Tau(rise)=Tau(fall): System is critically damped\n");
  }
#endif
  if (fabs(tes->Lin/Lcritp)-1. < 1e-6) {
    headas_chat(0,"L=L_plus: System is critically damped\n");
  }
  if (fabs(tes->Lin/Lcritm)-1. < 1e-5) {
    headas_chat(0,"L=L_minus: System is critically damped\n");
  }

  if (tes->Lin < Lcritm) {
    headas_chat(0,"L < L_minus: System is overdamped\n");
  }

  if (tes->Lin > Lcritp) {
    headas_chat(0,"L > L_plus: System is overdamped\n");
  }

  if ( (Lcritm < tes->Lin) && (tes->Lin<Lcritp) ) {
    headas_chat(0,"L_minus < L < L_plus: System is underdamped\n");
  }

  //
  // Noise and predicted energy resolution
  //

  // TES voltage noise
  double SV_tes=4.*kBoltz*tes->T_start*tes->R0*(1.+2.*tes->beta);

  // Load voltage noise
  double SV_L=4.*kBoltz*tes->T_start*tes->RL;

  // TFN power noise
  double SP_tfn=4.*kBoltz*tes->T_start*tes->T_start*tes->Gb1*tpow(tes);
  double ell2=ell*ell;
  double I02=tes->I0*tes->I0;
  double FWHM=2.*sqrt(2.*log(2.))*sqrt(
	      tau/ell2*sqrt((ell2*SP_tfn+I02*SV_tes+pow(ell-1.,2.)*I02*SV_L)*
		       (I02*SV_tes+I02*SV_L)));
  headas_chat(0,"Predicted energy resolution (FWHM, eV): %10.3f\n",FWHM/eV);
}


tesparams *tes_init(tespxlparams *par,int *status) {
  //
  // Initialize a TES pixel

  tesparams *tes=malloc(sizeof(tesparams));
  CHECK_NULL_RET(tes,*status,"Memory allocation failed for TES structure",NULL);

  //
  // random number initialization
  //
  rng=gsl_rng_alloc(gsl_rng_taus);

  // if seed is 0, do NOT use the rng's default seed (per gsl), but
  // initialize from the system clock
  if (par->seed==0) {
#if defined( __APPLE__) && defined(__MACH__)
    struct timeval tv;
    gettimeofday(&tv,NULL);
    tes->seed=1000000*tv.tv_sec+tv.tv_usec;
#else
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC,&tp);
    tes->seed=tp.tv_nsec;
#endif
  } else {
    tes->seed=par->seed;
  }
  gsl_rng_set(rng,tes->seed);

  // ID of this pixel
  tes->ID=strdup(par->ID);

  // The photon provider needs to be initialized outside of this
  tes->photoninfo=NULL;
  tes->get_photon=NULL;

  // we have not yet dealt with events
  tes->Nevts=0;

  tes->sample_rate=par->sample_rate; // Sample rate in Hz (typical 100-200kHz)
  tes->timeres=1./tes->sample_rate; // time resolution

  tes->decimate_factor=1; // step size wrt. sample rate

  // current to ADU conversion
  // note: we require imin<imax!
  // note: we assume 16 bits encoding
  assert(par->imin<par->imax);
  tes->imin=par->imin;
  tes->imax=par->imax;
  tes->aducnv=0xFFFE/(tes->imax-tes->imin);

  tes->tstart=par->tstart; // remember the starting time
  tes->time=tes->tstart;   // current time

  tes->delta_t=1./(tes->sample_rate*tes->decimate_factor); //integration step size

  // bandwidth to calculate the noise
  tes->bandwidth=1./(tes->delta_t*2.);

  tes->T_start=par->T_start;// initial operating temperature of the TES [K]
  tes->Tb=par->Tb;      // Heat sink/bath temperature [K]
  tes->R0=par->R0;        // Operating point resistance [Ohm]
  tes->RL=par->RL;     // Shunt / load resistor value
  tes->alpha=par->alpha;      // TES sensitivity (T/R*dR/dT)
  tes->beta=par->beta;       // TES current dependence (I/R*dR/dI)
  tes->Lin=par->Lin;     // Circuit inductance
  tes->n=par->n;          // Temperature dependence of the power flow to the heat sink
  // use this to turn the noise on or off
  tes->simnoise=par->simnoise;

  // hardcode for the moment!
  tes->mech=0;         // thermal coupling

  //JW STILL NEED TO CHECK THIS EQUATION!!!!!!!!!!!!
  // see discussion in I&H around eq 112
  tes->Gb1=220e-12*pow(tes->T_start/0.1,tes->n-1.); 
  tes->therm=1.0; // absorber thermalization time (in units of the step size h). 
                    // For a delta function input of power into the absorber 
                    // it is just set to 1

  // Calculate initial bias current from thermal power balance
  tes->I0_start=sqrt(tes->Gb1/(tes->n*pow(tes->T_start,tes->n-1.))*
		    (pow(tes->T_start,tes->n)-pow(tes->Tb,tes->n))/tes->R0);

  //JW STILL NEED TO CHECK THIS EQUATION!!!!!!!!!!!!
  //JW see discussion around I&H, eq. 111
  tes->Ce1=0.86e-12*(tes->T_start/0.1); //Absorber+TES heat capacity at Tc
  tes->dRdT=tes->alpha*tes->R0/tes->T_start;
  tes->dRdI=tes->beta*tes->R0/tes->I0_start;
 
  // level of SQUID readout and electronics noise
  //JW STILL NEED TO CHECK THIS EQUATION!!!!!!!!!!!!
  tes->squid_noise=2e-12*sqrt(tes->sample_rate/2*32*M_PI)*tes->simnoise;

  // set-up initial input conditions
  tes->I0=tes->I0_start;
  tes->V0=tes->I0*(tes->R0+tes->RL); // Effective bias voltage
  tes->RT=tes->R0;
  tes->T1=tes->T_start;

  // power flow
  tes->Pb1=tpow(tes); 


  // noise terms
  tes->Pnb1=0.;
  tes->Vdn=0.;
  tes->Vexc=0.;
  tes->Vcn=0.;
  tes->Vunk=0.;

  // ...define ODE system
  tes->odesys=malloc(sizeof(gsl_odeiv2_system));
  CHECK_NULL_RET(tes->odesys,*status,"cannot allocate ODE system",NULL);

  tes->odesys->function=TES_Tdifferential_NL;
  tes->odesys->jacobian=TES_jac;
  tes->odesys->dimension=2;
  tes->odesys->params=tes;

  // setup integrator
  tes->odedriver=gsl_odeiv2_driver_alloc_y_new(tes->odesys,gsl_odeiv2_step_rk4,
					       tes->delta_t,1e-6,1e-6);
  return(tes);
}

void tes_free(tesparams *tes) {
  //
  // clean up a tes structure, freeing all memory
  free(tes->ID);
  tes->ID=NULL;

  gsl_odeiv2_driver_free(tes->odedriver);
  tes->odedriver=NULL;

  free(tes->odesys);
  tes->odesys=NULL;

}


//
// logic problem in multiple calls: a photon read here that is
// after tstop will get lost
int tes_propagate(tesparams *tes, double tstop, int *status) {
  CHECK_STATUS_RET(*status,-1);

  // get first photon to deal with
  PixImpact impact;
  impact.time=tstop+100; // initialize to NO photon 
  if (tes->get_photon != NULL) {
    tes->get_photon(&impact,tes->photoninfo,status);
    CHECK_STATUS_RET(*status,-1);
  }

  int samples=0;
  while (tes->time<=tstop) {

    if (tes->simnoise) {
      // thermal noise
      tes->Pnb1=tnoi(tes);

      // Johnson noise terms
      tes->Vdn =gsl_ran_gaussian(rng,sqrt(4.*kBoltz*tes->T1*tes->RT*tes->bandwidth));
      tes->Vexc=gsl_ran_gaussian(rng,sqrt(4.*kBoltz*tes->T1*tes->RT*tes->bandwidth*2.*tes->beta));
      tes->Vcn =gsl_ran_gaussian(rng,sqrt(4.*kBoltz*tes->Tb*tes->RL*tes->bandwidth));
      tes->Vunk=gsl_ran_gaussian(rng,sqrt(4.*kBoltz*tes->T1*tes->RT*tes->bandwidth)*(1.+2*tes->beta) )*tes->m_unknown*tes->m_unknown;
    }

    // absorb next photon?
    tes->En1=0.;
    if (tes->time>=impact.time) {
      tes->Nevts++;
      tes->En1=impact.energy*keV/(tes->delta_t*tes->therm);

      // get the next photon 
      int success=tes->get_photon(&impact,tes->photoninfo,status);
      CHECK_STATUS_RET(*status,-1);
      if (success==0) {
	// there is no further photon to read. Set next impact time to
	// a time outside much after this
	impact.time=tstop+100.;
      }
    }
    double Y[2];
    Y[0]=tes->I0; // current
    Y[1]=tes->T1; // temperature

    int s=gsl_odeiv2_driver_apply_fixed_step(tes->odedriver,&(tes->time),tes->delta_t,1,Y);
    if (s!=GSL_SUCCESS) {
      fprintf(stderr,"Driver error: %d\n",s);
      return(1);
    }

    samples++;

    tes->I0=Y[0];
    tes->T1=Y[1];

    // put the pulse data and the time data into arrays 
    // if the decimation condition is met
    if (samples==tes->decimate_factor) {

      // write pulse data
      // we also add the SQUID/readout noise and subtract the
      // baseline (the equilibrium bias current) and invert the
      // pulses so they are all +ve
      double pulse=tes->I0_start-tes->I0;
      if (tes->simnoise) {
	pulse += gsl_ran_gaussian(rng,tes->squid_noise);
      }

      // write the sucker
      // (note that we do allow NULL here. This could be used,
      // e.g., to propagate the TES for a while without producing
      // output
      if (tes->write_to_stream != NULL ) {
	tes->write_to_stream(tes->time,pulse,tes->streaminfo,status);
      }

      samples=0;
    }

    // Update system properties

    // New resistance value assuming a simple 
    // linear transition with alpha and beta dependence
    tes->RT=tes->R0+tes->dRdT*(tes->T1-tes->T_start)+tes->dRdI*(tes->I0-tes->I0_start);

    // thermal power flow
    tes->Pb1=tpow(tes); 

  }

  return(0);
}

//int tes_propagate_(tesparams
