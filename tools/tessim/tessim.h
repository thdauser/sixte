#ifndef TESSIM_H
#define TESSIM_H

#include "sixt.h"

#include "sixteconfig.h"
#include "tesdatastream.h"
#include "pixelimpactfile.h"
#include "tesrecord.h"

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>



// parameters for the initializer
typedef struct {
  char *ID;         // pixel ID 
  char *impactlist; // file to read impacts from
  char *streamfile; // file to write data to

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
  double imin;    // minimum current to encode [A]
  double imax;    // maximum current to encode [A]

  int clobber;  // overwrite output files? -- IGNORED SO FAR
  int simnoise; // simulator noise 

  double m_unknown; // magnitude of unknown noise

  unsigned long seed; // seed of random number generator
  
} tespxlparams;

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

// stream write function
// this function is called whenever pulse data are written for the pixel
typedef void (*tes_stream_writer) (double time, double pulse, void *data, int *status);

// this function is called whenever photon data are written. Data is the streaminfo of tesparams
typedef void (*tes_photon_writer) (double time, long phid, void *data, int *status);


// 
// meta struct containing all physical parameters of the TES pixel
//
typedef struct {
  char *ID;      // string containing the pixel ID

  double time;   // current time
  double tstart; // start of simulation

  double delta_t;     // integration step size
  double sample_rate; // sample rate (Hz)
  double timeres;     // time resolution for this stream (1/sample_rate)
  int decimate_factor; // step size vs. sample rate
  double bandwidth; // needed for the noise simulation

  double imin;    // minimum current to encode [A]
  double imax;    // maximum current to encode [A]
  double aducnv;  // conversion factor current->adu

  double T_start; // initial operating temperature of the TES [K]
  double Tb;     // Heat sink/bath temperature [K]
  double R0;     // Operating point resistance [Ohm]
  double RL;     // Shunt/load resistor value [Ohm]
  double Lin;    // Circuit inductance [H] S
  double alpha;  // TES sensitivity (T/R*dR/dT)
  double beta;   // TES current dependence (I/R*dR/dI)
  double n;      // Temperature dependence of the power flow to the heat sink

  double dRdT;   // dR/dT at Tc
  double dRdI;   // dR/dI at Tc


  int mech;     // 1 for XXXX (MCCAMMON or IRWIN/HINTON CHAPTER) 
  double therm; // thermalization timescale in units of step size

  double Ce1;    // absorber+TES heat capacity at Tc 
  double Pb1;    // thermal power flow
  double Gb1;    // thermal conductance of the bath heat link at Tc


  double I0_start; // initial bias current [A]

  double V0;     // Effective bias voltage
  double I0;     // Effective current
  double T1;     // temperature

  double RT;     // resistance(T)

  unsigned long seed; // rng seed at start of simulation

  int simnoise;  // simulate noise?

  double Pnb1;   // thermal noise

  double Vdn;    // Johnson noise terms
  double Vcn;    // 
  double Vexc;   // 
  double Vunk;   // unknown noise 
  double m_unknown; // magnitude of unknown noise relative to Johnson (bath) noise

  double squid_noise; // SQUID readout and electronics noise 

  double En1;   // energy absorbed by a photon during this step [J]


  gsl_odeiv2_system *odesys;    // differential equation system
  gsl_odeiv2_driver *odedriver; // ODS solver

  long Nevts;                      // number of simulated events
  void *photoninfo;                // information for photon provider
  tes_photon_provider get_photon;  // photon provider (function pointer)

  void *streaminfo;              // data for pulse
  tes_stream_writer write_to_stream; // function to write the pulse (function pointer)
  tes_photon_writer write_photon; // function to save a processed photon (function pointer)

} tesparams;

// stream write function using TESDataStream
typedef struct {
  TESDataStream *stream; // output stream
  long streamind;    // index of next sample to write
  double tstart;     // starting time
  double tstop;      // end time
  double imin;       // minimum current to encode
  double imax;       // maximum current to encode
  double aducnv;     // conversion factor
  long Nt;           // number of elements in output buffer
} tes_datastream_info;

// initialize a TES data stream
tes_datastream_info *tes_init_datastream(double tstart, double tstop, tesparams *tes, int *status);
// append a pulse to the data stream. data points on a tes_datastream_info structure
void tes_append_datastream(double time,double pulse,void *data,int *status);
// write the TES data stream to file
void tes_save_datastream(char *streamfile, char *impactfile,tes_datastream_info *data, 
			 tesparams *tes, SixtStdKeywords *keywords, int *status);
// cleanup the TES data stream
void tes_free_datastream(tes_datastream_info **data, int *status);

// stream write function using TesRecords
typedef struct {
  TesRecord *stream; // output stream
  unsigned long streamind;    // index of next sample to write
  double tstart;     // starting time
  double tstop;      // end time
  double imin;       // minimum current to encode
  double imax;       // maximum current to encode
  double aducnv;     // conversion factor
  unsigned long Nt;  // number of elements in output buffer
} tes_record_info;

// initialize a TES data stream
tes_record_info *tes_init_tesrecord(double tstart, double tstop, tesparams *tes, int *status);
// append a pulse to the data stream. data points on a tes_datastream_info structure
void tes_append_tesrecord(double time,double pulse,void *data,int *status);
// append a photon to the data stream
void tes_append_photon_tesrecord(double time, long phid, void *data, int *status);


// write the TES data stream to file
void tes_save_tesrecord(char *streamfile, char *impactfile,tes_record_info *data, 
			 tesparams *tes, SixtStdKeywords *keywords, int *status);


// cleanup the TES data stream
void tes_free_tesrecord(tes_record_info **data, int *status);


// general functions
tesparams *tes_init(tespxlparams *par,int *status);
int tes_propagate(tesparams *tes, double tstop, int *status);
void tes_free(tesparams *tes);
void tes_print_params(tesparams *tes);
void tes_fits_write_params(fitsfile *fptr,tesparams *tes, int *status);


#endif
