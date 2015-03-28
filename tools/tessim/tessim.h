// 
// meta struct containing all physical parameters of the TES pixel
//

#ifndef TESSIM_H
#define TESSIM_H


#include "sixt.h"

#include "sixteconfig.h"
#include "tesdatastream.h"
#include "pixelimpactfile.h"

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

typedef struct {
  char *ID;         // pixel ID 
  char *impactlist; // file to read impacts from
  char *streamfile; // file name of output FITS stream

  double tstart; // start time of simulation
  double tstop;  // stop time of simulation

  double T_start; // initial operating temperature of the TES [K]
  double Tb;     // Heat sink/bath temperature [K]
  double R0;     // Operating point resistance [Ohm]
  double RL;     // Shut/load resistor value
  double alpha;  // TES sensitivity (T/R*dR/dT)
  double beta;   // TES current dependence (I/R*dR/dI)
  double Lin;    // Circuit inductance [H]
  double n;      // Temperature dependence of the power flow to the heat sink

  double sample_rate; // sample rate (Hz)
  double adu;    // conversion factor I1->ADU

  int clobber;  // overwrite output files? -- IGNORED SO FAR
  int simnoise; // simulator noise 

  unsigned long seed; // seed of random number generator
  
} tespxlparams;


typedef struct {
  char *ID;      // string containing the pixel ID

  double time;   // current time
  double tstart; // start of simulation
  double tstop;  // end of integration WILL BE REMOVED!!!!!!
  double delta_t;     // integration step size
  double sample_rate; // sample rate (Hz)
  double timeres;     // time resolution for this stream (1/sample_rate)

  double adu;    // conversion factor I1->ADU

  double bandwidth; // needed for the noise simulation --- EXPAND
  int decimate_factor; // step size vs. sample rate



  double T_start; // initial operating temperature of the TES [K]
  double Tb;     // Heat sink/bath temperature [K]
  double R0;     // Operating point resistance [Ohm]
  double RL;     // Shunt/load resistor value [Ohm]
  double Lin;    // Circuit inductance [H]
  double alpha;  // TES sensitivity (T/R*dR/dT)
  double beta;   // TES current dependence (I/R*dR/dI)
  double n;      // Temperature dependence of the power flow to the heat sink

  double dRdT;   // dR/dT at Tc
  double dRdI;   // dR/dI at Tc

  int simnoise;  // simulate noise?

  double Ce1;    // absorber+TES heat capacity at Tc
  double Pb1;    // thermal power flow 
  double Gb1;    // thermal conductance of the bath heat link at Tc

  double Pnb1;   // thermal noise

  double squid_noise; // SQUID readout and electronics noise

  double I0_start; // initial bias current [A]

  double V0;     // Effective bias voltage
  double I0;     // Effective current
  double T1; 

  double RT;     // resistance(T)

  double Vdn;    // Johnson noise terms
  double Vcn;
  double Vexc;

  double En1;   // energy absorbed by a photon during this step [J]

  int mech;     // 1 for (MCCAMMON or IRWIN/HINTON CHAPTER)
  double therm; // thermalization timescale in units of step size

  gsl_odeiv2_system *odesys; // differential equation system
  gsl_odeiv2_driver *odedriver; // ODS solver

  char *impactlist;    // file name of impactlist
  PixImpFile *impfile; // impact file (this will need to be changed for the case
                       // of multiple pixels)

  long Nevts;          // number of simulated events

  unsigned long seed; // rng seed at start of simulation

  TESDataStream *stream; // output stream
  long streamind;    // index of next sample to write
  char *streamfile;  // file name of output FITS stream

  long Nt;           // number of elements in output buffer
  SixtStdKeywords *keywords;   //Standard SIXT Keywords

} tesparams;

tesparams *tes_init(tespxlparams *par,int *status);
int tes_propagate(tesparams *tes, double tstop, int *status);
void tes_free(tesparams *tes);
void tes_print_params(tesparams *tes);

#endif
