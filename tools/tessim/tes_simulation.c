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
// ASTRONOMERS BEWARE: This code uses MKS units!

#include "tessim.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

gsl_rng *rng;

const double kBoltz=1.3806488e-23; // Boltzmann constant [m^2 kg/s^2]
const double eV=1.602176565e-19 ;  // 1eV [J]
const double keV=1.602176565e-16 ;  // 1keV [J]

// Differential equations for the TES to solve. This is the simplest
// calorimeter model assuming a single heat capacity for the 
// TES and absorber
int TES_Tdifferential_NL(double time, const double Y[], double eqn[], void *params) {

  // index 0: electrical equation I0=Y[0]
  // index 1: thermal equation    T1=Y[1]

  double II=Y[0];
  double TT=Y[1];

  tesparams *tes=(tesparams *) params;

  // RT(I0,T1)
  double RT=tes->R0+tes->dRdT*(TT-tes->T_start)+tes->dRdI*(II-tes->I0_start);

  // Electrical circuit equation: dI/dt -- dy_0/dt=f_0(t,y_0(t),y_1(t))
  // note: Vdn, Vexc, Vcn are the Johnson noise terms
  eqn[0]=(tes->V0-II*(tes->RL+RT)+tes->Vdn+tes->Vexc+tes->Vcn)/tes->Lin;

  // Thermal equation: dT/dt -- dy_1/dt=f_1(t,y_0(t),y_1(t))
  eqn[1]=(II*II*RT-tes->Pb1
	  -II*(tes->Vdn+tes->Vexc)
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


double tnoi (double T1, double T0, tesparams *tes) {
  // Thermal (phonon) noise terms

  double beta=tes->n-1.;
  double G1=tes->Gb1*pow(T1/tes->T_start,beta);
  
  double gamma;

  // Irwin, eq. 93
  if (tes->mech==0) {
    gamma=(pow(T0/T1,beta+2.0)+1.0)/2.0;
  } else {
    gamma=1.;
  }
 
  // use gsl_ran_gaussian instead?
  return(gsl_ran_ugaussian(rng)*sqrt(4*kBoltz*T1*T1*G1*gamma*tes->bandwidth));
}

double tpow(double T1,double T0, tesparams *tes) {
  // Calculate the power flow between elements of the pixel
  double beta=tes->n-1.0;
  double G1=tes->Gb1*pow(T1/tes->T_start,beta);

  return G1*(pow(T1,tes->n)-pow(T0,tes->n))/(tes->n*pow(T1,beta));
}

void AddToStream(double time,double pulse, tesparams *tes) {
  // only write sample if we still have space in the stream
  if (tes->streamind < tes->stream->Ntime) {
    tes->stream->time[tes->streamind]=time;
    tes->stream->adc_value[tes->streamind++][0]=(uint16_t) (pulse*tes->adu);
  }
}

void WriteStream(tesparams *tes, int * const status) {

  fitsfile *fptr;
  createTESFitsStreamFile(&fptr,
			  tes->streamfile,
			  tes->keywords->telescop,
			  tes->keywords->instrume,
			  tes->keywords->filter,
			  tes->keywords->ancrfile,
			  tes->keywords->respfile,
			  "none", // xmlfile
			  tes->impfile->fptr->Fptr->filename, // name of impactfile
			  tes->keywords->mjdref,
			  tes->keywords->timezero,
			  tes->keywords->tstart,
			  tes->keywords->tstop,
			  'y', //clobber
			  status);

  TESFitsStream *stream=newTESFitsStream(status);
  snprintf(stream->name,9,"ADC%03d",1);
  allocateTESFitsStream(stream,tes->stream->Ntime,1,status);
  // copy over (this is stupid for a single pixel)
  long ii;
  stream->pixID[0]=0;
  for (ii=0; ii<stream->Ntime; ii++) {
    stream->time[ii]=tes->stream->time[ii];
    stream->adc_value[0][ii]=tes->stream->adc_value[ii][0];
  }
  writeTESFitsStream(fptr,
		     stream,
		     tes->keywords->tstart,
		     tes->keywords->tstop,
		     tes->timeres,
		     &(tes->Nevts), 
		     0,-1,status);
  fits_close_file(fptr,status);
  CHECK_STATUS_VOID(*status);

  destroyTESFitsStream(stream);
  free(stream);

}


void tes_print_params(tesparams *tes) {
  headas_chat(0,"\nStatus of TES Pixel with ID: %s\n",tes->ID);
  headas_chat(0,"\n");
  headas_chat(0,"Start time of simulation [s]           : %15.6f\n",tes->tstart);
  headas_chat(0,"Current time of simulation [s]         : %15.6f\n",tes->time);
  headas_chat(0,"Sample rate [Hz]                       : %15.6f\n",tes->sample_rate);
  headas_chat(0,"Integration step size [ms]             : %15.6f\n",1000.*tes->delta_t);
  headas_chat(0,"ADU conversion factor [ADU/A]          : %15.1f\n",tes->adu);
  headas_chat(0,"\n");

  headas_chat(0,"Initial bias current [mA]              : %10.5f\n",1000.*tes->I0_start);
  headas_chat(0,"Operating point resistance [mOhm]      : %10.5f\n",1000.*tes->R0);
  headas_chat(0,"Shunt/load resistor value [mOhm]       : %10.5f\n",1000.*tes->RL);
  headas_chat(0,"Circuit inductance [H]                 : %10.5f\n",tes->Lin);
  headas_chat(0,"TES sensitivity T/R*dR/dT (alpha)      : %10.3f\n",tes->alpha);
  headas_chat(0,"TES current dependence I/R*dR/dI (beta): %10.3f\n",tes->beta);
  headas_chat(0,"\n");

  headas_chat(0,"Initial operating temperature [mK]     : %10.3f\n",1000.*tes->T_start);
  headas_chat(0,"Heat sink temperature [mK]             : %10.3f\n",1000.*tes->Tb);
  headas_chat(0,"Heat sink coupling parameter n         : %10.2f\n",tes->n);
  headas_chat(0,"\n");

  headas_chat(0,"Absorber+TES heat capacity at Tc       : %10.5f\n",tes->Ce1);
  headas_chat(0,"Thermal power flow                     : %10.5f\n",tes->Pb1);
  headas_chat(0,"Heat link thermal conductance at Tc    : %10.5f\n",tes->Gb1);
  headas_chat(0,"\n");

  headas_chat(0,"Effective bias voltage [mV]            : %10.5f\n",1000.*tes->V0);
  headas_chat(0,"Current [mA]                           : %10.5f\n",1000.*tes->I0);
  headas_chat(0,"Temperature [mK]                       : %10.5f\n",1000.*tes->T1);
  headas_chat(0,"Current Resistivity [Ohm]              : %10.5f\n",tes->RT);
  headas_chat(0,"\n");
    
  headas_chat(0,"Seed of random number generator        : %10lu\n",tes->seed);

  if (tes->simnoise) {
    headas_chat(0,"Simulating noise terms\n");
  } else {
    headas_chat(0,"NOT simulating noise terms\n");
  }

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
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC,&tp);
    tes->seed=tp.tv_nsec;
  } else {
    tes->seed=par->seed;
  }
  gsl_rng_set(rng,tes->seed);

  tes->ID=strdup(par->ID);

  tes->streamfile=strdup(par->streamfile);


  // Open the pixel impact file
  // in this initial version we assume that we have one
  // impact file per pixel. 
  // This will obviously be changed later!!
  tes->impactlist=strdup(par->impactlist);

  tes->impfile=openPixImpFile(tes->impactlist, READONLY,status);
  CHECK_STATUS_RET(*status,NULL);

  // we have not yet dealt with events
  tes->Nevts=0;

  // Read standard SIXTE keywords from input file
  tes->keywords=newSixtStdKeywords(status);
  sixt_read_fits_stdkeywords(tes->impfile->fptr,
			     tes->keywords,
			     status);
  CHECK_STATUS_RET(*status,NULL);

  tes->sample_rate=par->sample_rate; // Sample rate in Hz (typical 100-200kHz)
  tes->timeres=1./tes->sample_rate; // time resolution

  tes->decimate_factor=1; // step size wrt. sample rate

  tes->adu=par->adu;

  tes->tstart=par->tstart; // current time is starting time
  tes->time=tes->tstart;

  tes->tstop=par->tstop; // DUMMY FOR THE MOMENT - MAX. TIME FOR INTEGRATION

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

  tes->Gb1=220e-12*pow(tes->T_start/0.1,tes->n-1.); 
  tes->therm=1.0; // absorber thermalization time (in units of the step size h). 
                    // For a delta function input of power into the absorber 
                    // it is just set to 1

  // Calculate initial bias current from thermal power balance
  tes->I0_start=sqrt(tes->Gb1/(tes->n*pow(tes->T_start,tes->n-1.))*
		    (pow(tes->T_start,tes->n)-pow(tes->Tb,tes->n))/tes->R0);
  tes->Ce1=0.86e-12*(tes->T_start/0.1); //Absorber+TES heat capacity at Tc
  tes->dRdT=tes->alpha*tes->R0/tes->T_start;
  tes->dRdI=tes->beta*tes->R0/tes->I0_start;
 
  // level of SQUID readout and electronics noise
  tes->squid_noise=2e-12*sqrt(tes->sample_rate/2*32*M_PI)*tes->simnoise;

  // set-up initial input conditions
  tes->I0=tes->I0_start;
  tes->V0=tes->I0*(tes->R0+tes->RL); // Effective bias voltage
  tes->RT=tes->R0;
  tes->T1=tes->T_start;

  // power flow
  tes->Pb1=tpow(tes->T1,tes->Tb,tes); 


  // noise terms
  tes->Vdn=0.;
  tes->Vexc=0.;
  tes->Vcn=0.;
  tes->Pnb1=0.;

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

  //
  // allocate memory for output
  // (this should be done in a more intelligent way
  // - we should not have to allocate the stream in memory!)
  //

  tes->Nt=(long) ((tes->tstop-tes->tstart)*tes->sample_rate);

  tes->streamind=0;
  tes->stream=newTESDataStream(status);
  if (*status==EXIT_FAILURE) {
    SIXT_ERROR("memory allocation failed for stream structure");
    CHECK_STATUS_RET(*status,NULL);
  }
  allocateTESDataStream(tes->stream,tes->Nt,1,status);
  CHECK_STATUS_RET(*status,NULL);

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

  free(tes->impactlist);
  tes->impactlist=NULL;

  int status;
  freePixImpFile(&(tes->impfile),&status);
  
  destroyTESDataStream(tes->stream);
  free(tes->stream);
  tes->stream=NULL;

  free(tes->streamfile);
  tes->streamfile=NULL;

  freeSixtStdKeywords(tes->keywords);
}


int tes_propagate(tesparams *tes, double tstop, int *status) {

  // first photon
  PixImpact impact;
  getNextImpactFromPixImpFile(tes->impfile,&impact,status);
  CHECK_STATUS_RET(*status,-1);

  int samples=0;
  while (tes->time<=tstop) {

    if (tes->simnoise) {
      // thermal noise
      tes->Pnb1=tnoi(tes->T1,tes->Tb,tes);

      // Johnson noise terms
      tes->Vdn =gsl_ran_gaussian(rng,sqrt(4*kBoltz*tes->T1*tes->RT*tes->bandwidth));
      tes->Vexc=gsl_ran_gaussian(rng,sqrt(4*kBoltz*tes->T1*tes->RT*tes->bandwidth*2.*tes->beta));
      tes->Vcn =gsl_ran_gaussian(rng,sqrt(4*kBoltz*tes->Tb*tes->RL*tes->bandwidth));
    }

    // absorb next photon?
    tes->En1=0.;
    if (tes->time>=impact.time) {
      tes->Nevts++;
      tes->En1=impact.energy*keV/(tes->delta_t*tes->therm);

      // get the next photon 
      int success=getNextImpactFromPixImpFile(tes->impfile,&impact,status);
      CHECK_STATUS_RET(*status,-1);
      if (success==0) {
	// no further photon to read. Set next impact time to
	// a time outside much after this
	impact.time=tes->keywords->tstop+100.;
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
      AddToStream(tes->time,pulse,tes);

      samples=0;
    }

    // Update system properties

    // New resistance value assuming a simple 
    // linear transition with alpha and beta dependence
    tes->RT=tes->R0+tes->dRdT*(tes->T1-tes->T_start)+tes->dRdI*(tes->I0-tes->I0_start);

    // thermal power flow
    tes->Pb1=tpow(tes->T1,tes->Tb,tes); 

  }

  // Dump results
  WriteStream(tes,status);

  // cleanup
  //  freePixImpFile(&(tes->impfile), status);
  // CHECK_STATUS_RET(*status,-1);

  // gsl_odeiv2_driver_free(d);
  // gsl_rng_free(rng);

  return(0);
}
