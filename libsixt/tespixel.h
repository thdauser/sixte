#ifndef TESPIXEL_H
#define TESPIXEL_H

#include "sixt.h"

#include "sixteconfig.h"
#include "pixelimpactfile.h"
#include "progressbar.h"
#include <gsl/gsl_complex.h>

#include <gsl/gsl_odeiv2.h>


// photon provider function
// a call to this function returns the next photon to be processed
// providerinfo contains arbitrary information used by the function
// must return 1 if successful, 0 if no further photon to process
typedef int (*tes_photon_provider) (PixImpact *photon, void *providerinfo, int *status);


// photon provider functions and their associated data

//get photon from impact file
void tes_photon_from_impact(void *data);

typedef struct {
  char *impactlist;            // file name of impactlist
  PixImpFile *impfile;         // impact file
  SixtStdKeywords *keywords;   // Standard SIXT Keywords
} tes_impactfile_info;

tes_impactfile_info *tes_init_impactlist(char *impactfile, int *status);
int tes_photon_from_impactlist(PixImpact *photon,void *data, int *status);
void tes_free_impactlist(tes_impactfile_info **data, int *status);

//get photon from impact buffer
typedef struct {
// TODO fill this!
  PixImpact *impacts;         // array of impacts saved in this buffer
  long numimpacts;            // number of impacts in this buffer
  long nextimpact;            // id of the next impact to be read in
                              // if this equals numimpacts, there are no more impacts
} tes_impactbuffer_info;

tes_impactbuffer_info *tes_init_impactbuffer(int nImpacts, int *status);
int tes_photon_from_impactbuffer(PixImpact *photon, void *data, int *status);
void tes_free_impactbuffer(tes_impactbuffer_info **data, int *status);


// forward declare of tesparams
typedef struct tesparams tesparams;

// stream write function
// this function is called whenever pulse data are written for the pixel
typedef void (*tes_stream_writer) (tesparams *tes,double time, double pulse, int *status);

// this function is called whenever photon data are written. Data is the streaminfo of tesparams
typedef void (*tes_photon_writer) (tesparams *tes, double time, long phid, int *status);

// dRdI provider
// a call to this function returns the current value of dR/dI
typedef double (*dRdI_provider) (tesparams *tes, const double Y[]);

double get_dRdI(tesparams *tes, const double Y[]);

// dRdT provider
// a call to this function returns the current value of dR/dT
typedef double (*dRdT_provider) (tesparams *tes, const double Y[]);

double get_dRdT(tesparams *tes, const double Y[]);

typedef double (*RTI_provider) (tesparams *tes, const double Y[]);

double get_RTI(tesparams *tes, const double Y[]);

//Frame hits structure
struct structLinkedFrameImpacts{
	double time; //Time of the impact
	double energy; //Deposited energy in the frame
	struct structLinkedFrameImpacts* next;
};

typedef struct structLinkedFrameImpacts LinkedFrameImpacts;

//Frame hits initialiser
LinkedFrameImpacts* newLinkedFrameImpacts(int* const status);

//Frame hits destroyers
void freeLinkedFrameImpacts(LinkedFrameImpacts** const list);

// Frame hit generator
LinkedFrameImpacts* frameHitsList(tesparams *tes, int* const status);

typedef double (*DeltaTb_provider) (tesparams *tes);

double get_DeltaTb(tesparams *tes);

// 
// meta struct containing all physical parameters of the TES pixel
//
struct tesparams {
  char *type;         // string containing the pixel type
  int id;             // number of the pixel

  double time;        // current simulation time
  double tstart;      // start of simulation

  double delta_t;     // integration step size
  double sample_rate; // sample rate (Hz)
  double timeres;     // time resolution for this stream (1/sample_rate)
  unsigned int decimate_factor; // step size vs. sample rate
  double bandwidth; // needed for the noise simulation

  double imin;    // minimum current to encode [A]
  double imax;    // maximum current to encode [A]
  double aducnv;  // conversion factor current->adu

  int acdc;       // boolean for AC biased (true) or DC biased (false)


  double T_start; // initial operating temperature of the TES [K]
  double Tb_start;     // Heat sink/bath temperature at [K]
  DeltaTb_provider Tb; // Delta Heat sink/bath temperature due to frame hits at time t[K]
  LinkedFrameImpacts* frame_hits; //tes frame hits
  double R0;     // Operating point resistance [Ohm]
  double RL;     // Shunt/load resistor value [Ohm]
  double Rpara;  // Parasitic resistor value [Ohm]

  double TTR;    // Transformer Turns Ratio

  double Reff;   // effective resistance (derived from Rl or Rpara) [Ohm]


  double Lin;    // Circuit inductance [H] [only used if DC biased]
  double Lfilter;// Filter inductance [H] [only used if AC biased]

  double Leff;   // Effective inductance (derived from Lin or Lfilter)

  double alpha;  // TES sensitivity (T/R*dR/dT)
  double beta;   // TES current dependence (I/R*dR/dI)
  double n;      // Temperature dependence of the power flow to the heat sink

  dRdT_provider dRdT;   // dR/dT at current timestep
  dRdI_provider dRdI;   // dR/dI at current timestep

  int mech;     // 1 for XXXX (MCCAMMON or IRWIN/HILTON CHAPTER)
  double therm; // thermalization timescale in units of step size

  double Ce1;    // absorber+TES heat capacity at Tc [J/K]
  double Gb1;    // thermal conductance of the bath heat link at Tc [W/K]
  double Pb1;    // thermal power flow


  double I0_start; // initial bias current [A]

  double V0;     // Effective bias voltage
  double I0;     // Effective current
  double T1;     // temperature

  double RT;     // resistance(T)
  RTI_provider RTI; //
  double bias;   // bias percentage (%)

  unsigned long seed; // rng seed at start of simulation

  int simnoise;  // simulate noise?
  int stochastic_integrator; //use stochastic integrator?

  double Pnb1;   // thermal noise

  double Vdn;    // Johnson noise terms
  double Vcn;    // 
  double Vexc;   // 
  double Vunk;   // excess noise 
  double m_excess; // magnitude of excess noise relative to Johnson (bath) noise

  double squid_noise; // SQUID readout and electronics noise 

  int twofluid; //Do we use the twofluid model?
  //TODO Change once more models become available

  int frame_hit; //Option to use frame hits
  double frame_hit_time; //Time of frame event (s)
  char* frame_hit_file; //File name of frame hit model
  int frame_hit_file_length; //File length

  double* frame_hit_time_dep; //Pulse time dependence
  double* frame_hit_shape; //Pulse shape

  PixImpact *impact;   // next photon impact
  double En1;   // energy absorbed by photons during this step [J]
  int n_absorbed;  // number of photons absorbed during this step

  gsl_odeiv2_system *odesys;    // differential equation system
  gsl_odeiv2_driver *odedriver; // ODS solver

  long Nevts;                      // number of simulated events
  void *photoninfo;                // information for photon provider
  tes_photon_provider get_photon;  // photon provider (function pointer)

  void *streaminfo;              // data for pulse
  tes_stream_writer write_to_stream; // function to write the pulse (function pointer)
  tes_photon_writer write_photon; // function to save a processed photon (function pointer)

  progressbar *progressbar; // progress bar (or NULL if we don't show it)

  // FDM crosstalk parameters
  gsl_complex Ioverlap; // carrier overlap current in the pixel's readout
  gsl_complex Iout_start; // output current (I0 + Ioverlap) for the first sample; includes phase rotation
  double Pcommon; // common impedance power for this pixel
  double theta_Vb; // phase of bias Voltage wrt. on-resonance current [radians]

  // readout mode for crosstalk
  int readoutMode;
};

#endif
