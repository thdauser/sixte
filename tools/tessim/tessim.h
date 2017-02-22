#ifndef TESSIM_H
#define TESSIM_H

#include "sixt.h"

#include "sixteconfig.h"
#include "tesdatastream.h"
#include "pixelimpactfile.h"
#include "testriggerfile.h"
#include "tesrecord.h"
#include "progressbar.h"
#include "tespixel.h"
#include "advdet.h"
#include "crosstalk.h"
#include "tescrosstalk.h"

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

  int acdc;    // ac biased?

  double T_start; // initial operating temperature of the TES [K]
  double Tb;     // Heat sink/bath temperature [K]
  double R0;     // Operating point resistance [Ohm]
  double RL;     // Shut/load resistor value [Ohm]
  double Rpara;  // Parasitic resistor value [Ohm]

  double TTR;    // Transformer turns ratio

  double alpha;  // TES sensitivity (T/R*dR/dT)
  double beta;   // TES current dependence (I/R*dR/dI)
  double Lin;    // Circuit inductance [H]
  double Lfilter;// Filter inductance [H]
  double Ce1;    // absorber+TES heat capacity at Tc [J/K]
  double Gb1;    // thermal conductance of the bath heat link at Tc [W/K]
  double n;      // Temperature dependence of the power flow to the heat sink
  double I0;     // current at I0 [A]
  double V0;     // voltage bias [V] (computed for initial equilibrium if negative)

  double sample_rate; // sample rate (Hz)
  double imin;    // minimum current to encode [A]
  double imax;    // maximum current to encode [A]

  char *trigger; // string defining the trigger strategy
  unsigned long preBufferSize; // number of samples to keep before trigger
  unsigned long triggerSize; // total number of samples to record

  int clobber;  // overwrite output files?
  int simnoise; // simulator noise 

  double m_excess; // magnitude of excess noise

  unsigned long seed; // seed of random number generator

  int showprogress;   // show progressbar?

  int doCrosstalk;   // do crosstalk?
  
} tespxlparams;


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
                     //   3: noise
                     //   4: trigger on impact
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
#define TRIGGER_NOISE 3
#define TRIGGER_IMPACT 4

// append a photon to the data stream
// NOT YET IMPLEMENTED!
// void tes_append_photon_trigger(tesparams *tes, double time, long phid, int *status);

// finish writing triggers
void tes_close_trigger(tesparams *tes, int *status);

// cleanup the TES data stream
void tes_free_trigger(tes_trigger_info **data, int *status);




// general functions
tesparams *tes_init(tespxlparams *par,int *status);
int tes_propagate(AdvDet *det, double tstop, int *status);
void tes_free(tesparams *tes);
void tes_print_params(tesparams *tes);
void tes_fits_write_params(fitsfile *fptr,tesparams *tes, int *status);
void tes_fits_read_params(char *file,tespxlparams *tes, int *status);

#define PROGRESSBAR_FACTOR 10

#endif
