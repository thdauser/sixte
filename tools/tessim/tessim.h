#ifndef TESSIM_H
#define TESSIM_H

#include "sixt.h"

#include "sixteconfig.h"
#include "tesdatastream.h"
#include "pixelimpactfile.h"
#include "testriggerfile.h"
#include "tesrecord.h"

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>



// parameters for the initializer
typedef struct {
  char *type;       // type of this pixel
  int id;           // number of the pixel
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
  double Ce1;    // absorber+TES heat capacity at Tc [J/K]
  double Gb1;    // thermal conductance of the bath heat link at Tc [W/K]
  double n;      // Temperature dependence of the power flow to the heat sink
  double I0;     // current at I0 [A]

  double sample_rate; // sample rate (Hz)
  double imin;    // minimum current to encode [A]
  double imax;    // maximum current to encode [A]

  char *trigger; // string defining the trigger strategy
  unsigned long preBufferSize; // number of samples to keep before trigger
  unsigned long triggerSize; // total number of samples to record

  int clobber;  // overwrite output files? -- IGNORED SO FAR
  int simnoise; // simulator noise 

  double m_excess; // magnitude of excess noise

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

// forward declare of tesparams
typedef struct tesparams tesparams;

// stream write function
// this function is called whenever pulse data are written for the pixel
typedef void (*tes_stream_writer) (tesparams *tes,double time, double pulse, int *status);

// this function is called whenever photon data are written. Data is the streaminfo of tesparams
typedef void (*tes_photon_writer) (tesparams *tes, double time, long phid, int *status);

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

  double Ce1;    // absorber+TES heat capacity at Tc [J/K]
  double Gb1;    // thermal conductance of the bath heat link at Tc [W/K]
  double Pb1;    // thermal power flow


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
  double Vunk;   // excess noise 
  double m_excess; // magnitude of excess noise relative to Johnson (bath) noise

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

};

//
// stream write function using TESDataStream
//
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
void tes_append_datastream(tesparams *tes,double time,double pulse,int *status);
// write the TES data stream to file
void tes_close_datastream(tesparams *tes,char *streamfile, char *impactfile,tes_datastream_info *data, 
			 SixtStdKeywords *keywords, int *status);
// cleanup the TES data stream
void tes_free_datastream(tes_datastream_info **data, int *status);

//
// stream write function using TesRecords
//

// buffer size of the tesrecord buffer
// this is a compromise between memory usage and speed
#define TESRECORD_BUFFERSIZE 5000000

typedef struct {
  TesRecord *stream; // output stream
  LONGLONG streamind;    // index of next sample to write
  double tstart;     // starting time
  double tstop;      // end time
  double imin;       // minimum current to encode
  double imax;       // maximum current to encode
  double aducnv;     // conversion factor
  long Nt;           // number of elements in output buffer
  int clobber;       // clobber?
  char *streamfile;  // file name of file we're going to write
  char *impactfile;  // impact file we're using
  fitsfile *fptr;    // file to write to
  SixtStdKeywords *keywords;   // Standard SIXT Keywords
  LONGLONG row;      // last row in fptr we've written
  int timecol;       // column for the time
  int adccol;        // column for the ADC value
  int curcol;        // column for the current
} tes_record_info;

// initialize a TES datastream built on tesrecords
tes_record_info *tes_init_tesrecord(double tstart, double tstop, tesparams *tes, 
				    char *streamfile, char *impactfile,int clobber,
				    SixtStdKeywords *keywords,
				    int *status);
// append a pulse to the data stream. 
void tes_append_tesrecord(tesparams *tes,double time,double pulse,int *status);
// append a photon to the data stream
void tes_append_photon_tesrecord(tesparams *tes, double time, long phid, int *status);

// finish writing datastream to a file
void tes_close_tesrecord(tesparams *tes, int *status);

// cleanup the TES data stream
void tes_free_tesrecord(tes_record_info **data, int *status);

//
// stream write function using TesRecords and triggers
//

typedef struct {
  TesRecord *fifo; // buffer for our fifo
  unsigned int fifoind; // next element to write in fifo

  TesRecord *stream; // output stream after trigger
  unsigned long preBufferSize; // number of samples to keep before trigger
  unsigned long triggerSize; // total number of samples to record
  unsigned int streamind;    // index of next sample to write

  double tstart;     // starting time
  double tstop;      // end time
  double imin;       // minimum current to encode
  double imax;       // maximum current to encode
  double aducnv;     // conversion factor
  int strategy;      // trigger strategy:
                     //   1: moving average
                     //   2: differentiation
  double threshold;  // threshold above which trigger is valid
  unsigned int npts; // number of previous points to include in trigger calculation
  long helper[3]; // helper variables
  unsigned int CanTrigger; // if 0(!!) we can trigger
  unsigned int SuppressTrigger; // #samples not to trigger after a trigger
  TesTriggerFile *fptr;    // file to write to
} tes_trigger_info;

// initialize a TES trigger built on tesrecords
tes_trigger_info *tes_init_trigger(double tstart, double tstop, tesparams *tes, 
				   int strategy,unsigned long preBufferSize,
				   unsigned long triggerSize,
				   double threshold,unsigned int npts,unsigned int suppress,
				   char *streamfile, char *impactfile,int clobber,
				   SixtStdKeywords *keywords,
				   int *status);
// append a pulse to the trigger data stream. 
// this is the routine that also does the trigger
void tes_append_trigger(tesparams *tes,double time,double pulse,int *status);

// trigger types
#define TRIGGER_MOVAVG 1
#define TRIGGER_DIFF 2

// append a photon to the data stream
// NOT YET IMPLEMENTED!
// void tes_append_photon_trigger(tesparams *tes, double time, long phid, int *status);

// finish writing triggers
void tes_close_trigger(tesparams *tes, int *status);

// cleanup the TES data stream
void tes_free_trigger(tes_trigger_info **data, int *status);




// general functions
tesparams *tes_init(tespxlparams *par,int *status);
int tes_propagate(tesparams *tes, double tstop, int *status);
void tes_free(tesparams *tes);
void tes_print_params(tesparams *tes);
void tes_fits_write_params(fitsfile *fptr,tesparams *tes, int *status);
void tes_fits_read_params(char *file,tespxlparams *tes, int *status);


#endif
