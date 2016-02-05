/************************************************************************************************

   This file is part of SIXTE/SIRENA software.

   SIXTE/SIRENA is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE/SIRENA is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

   Copyright 2014:  This file has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

************************************************************************************************/

#ifndef INTEGRASIRENA_H
#define INTEGRASIRENA_H 1

#include "tesrecord.h"
#include "teseventlist.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct {
  /** Number of rows */
  int matrixRows;

  /** Number of columns */
    int matrixColumns;

  /** Matrix */
  gsl_matrix *matrixBody;

} MatrixStruct;

typedef struct {
  /** Duration of the pulse template in the library */
  int template_duration;

  /** Vector containing the pulse template */
  gsl_vector *ptemplate;

  /** Energy of the template */
  double energy;

  /** Pulse Height of the template */
  double pulse_height;

} PulseTemplate;

typedef struct {
  /** Duration of the matched filter in the library */
  int mfilter_duration;

  /** Vector containing the matched filter */
  gsl_vector *mfilter;

  /** Energy of the mfilter */
  double energy;

  /** Pulse Height of the mfilter */
  double pulse_height;

} MatchedFilter;

typedef struct {
  // Duration of the optimalfilter
  int ofilter_duration;

  // Vector containing the optimal filter
  gsl_vector *ofilter;

  /** Energy of the ofilter */
  double energy;

  // Normalization factor
  double nrmfctr;

} OptimalFilterSIRENA;

typedef struct {
  /** Number of templates & matched filters in the structure. */
  int ntemplates;

  /** Energies of the templates */
  gsl_vector *energies;

  /** Pulse Heights of the templates */
  gsl_vector *pulse_heights;

  /** Structure containing all the pulse templates from the library */
  PulseTemplate* pulse_templates;

  /** Structure containing all the pulse templates filtered & derivated from the library */
  PulseTemplate* pulse_templates_filder;

  /** Maximum of pulse_templates_filder */
  gsl_vector *maxDERs;

  /** Structure containing all the pulse templates from the library */
  PulseTemplate* pulse_templates_B0;
  
  /** Structure containing all the matched filters from the library */
  MatchedFilter* matched_filters;

  /** Structure containing all the matched filters from the library */
  MatchedFilter* matched_filters_B0;

  /** Structure containing all the optimal filters from the library */
  OptimalFilterSIRENA* optimal_filters;

  /** Normalization factors of the optimal filters */
  gsl_vector *nrmfctrs;
  
  /** Structure containing all the optimal filters AB from the library */
  //OptimalFilterSIRENA* optimal_filtersab;

  /** Normalization factors of the optimal filters AB */
  //gsl_vector *nrmfctrsab;

  /** Weight matrix */
  //MatrixStruct* W;
  gsl_matrix *V;
  gsl_matrix *W;

  /** T vector */
  gsl_matrix *T;

  /** t escalar */
  gsl_vector *t;

  /** X matrix */
  gsl_matrix *X;

  /** Y vector */
  gsl_matrix *Y;

  /** Z vector */
  gsl_matrix *Z;

  /** r escalar */
  gsl_vector *r;

  /** PAB vector */
  gsl_matrix *PAB;

  /** DAB vector */
  gsl_matrix *DAB;
  
  /** Structure containing all the optimal filters AB from the library */
  OptimalFilterSIRENA* optimal_filtersab;

  /** Normalization factors of the optimal filters AB */
  gsl_vector *nrmfctrsab;

} LibraryCollection;

typedef struct {

  /** Noise standard deviation */
  double noiseStd;
  
  /** Baseline: BASELINE or BASELINR */
  double baseline;
  
  /** Duration of the noise spectrum */
  int noise_duration;
	
  /** Vector containing the noise spectrum */
  gsl_vector *noisespec;

  /** Vector containing the frequecies of the noise spectrum */
  gsl_vector *noisefreqs;

} NoiseSpec;

typedef struct {

  /** Pulse duration (maximum length) */
  int pulse_duration;

  /** Length of filter used during reconstruction (samples)*/
  int grade1;
 
  /** Distance to previous pulse in record (samples)*/
  /** tstart(i)-tend(i-1)*/
  int grade2;
  
  /** Distance to the begining of the previous pulse in record (samples)*/
  /** tstart(i)-tstart(i-1)*/
  int grade2_1;

  /** PIX_ID of the detected pulse*/
  int pixid;
  
  /** Vector containing the pulse adc values */
  gsl_vector *pulse_adc;

  /** Start time of the Pulse */
  double Tstart;

  /** End time of the Pulse */
  double Tend;
 
  /** Rise time of the Pulse */
  double riseTime;

  /** Fall time of the Pulse */
  double fallTime;

  /** Pulse height of the Pulse */
  double pulse_height;

  /** Maximum of the filtered-derived pulse */
  double maxDER;

  /** Energy (KeV) of the Pulse */
  double energy;

  /** Quality of the Pulse */
  double quality;

} PulseDetected;

typedef struct {
  /** Number of detected pulses in the structure. **/
  int ndetpulses;

  /** Array containing all the pulses detetected in record**/
  PulseDetected* pulses_detected;

} PulsesCollection;

/** Structure containing all the pointers and values to run the reconstruction stage */
typedef struct {
 	/**                   **/
	/** SIRENA parameters **/
	/**                   **/
        /** LibraryCollection structure (pulse templates and matched filters)*/		
	LibraryCollection* library_collection;
	
	/** Threshold of each record **/
	double threshold;

	/** Library file (to be used or to be created) **/
	char library_file[256];
	
	/** records file (if any) **/
	char record_file[256];
	
	/** records file pointer **/
	fitsfile *record_file_fptr;

	/** Noise file (to be used) **/
	char noise_file[256];

	/** Output event file **/
	char event_file[256];

	/** Pulse length */
	int pulse_length;
	
	/** Pulses Fall time (s)**/
	double tauFall;

	/** Detection scaleFactor (0.005 – no filtering) **/
	double scaleFactor;
	
	/** Detection samplesUp (samples to confirm threshold overcoming) **/
	double samplesUp;

	/** Detection nSgms (sigmas to establish a threshold for detection) **/
	double nSgms;

	/** Monochromatic energy for library creation **/
	double monoenergy;
	
	/** Running sum length for the RS raw energy estimation, in seconds (only in CALIBRATION) **/
	double LrsT;
	
	/** Baseline averaging length for the RS raw energy estimation, in seconds (only in CALIBRATION) **/
	double LbT;
	
	/** Baseline (in ADC units) **/
	//double baseline;

	/** Run mode (0: calibration/lib creation  1:energy reconstruction) **/
	int mode;

	/** Noise spectrum **/
	NoiseSpec* noise_spectrum;

	/** Pixel Type: SPA, LPA1, LPA2, LPA3 **/
	char PixelType[5];

	/** Filtering Domain: T (Time) or F (Frequency) **/
	char FilterDomain[2];

	/** Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline) **/
	char FilterMethod[3];

	/** Energy Method: OPTFILT, WEIGHT, WEIGHTN or I2R **/
	char EnergyMethod[8];

	//LagsOrNot: LAGS == True or NOLAGS == False **/
	int LagsOrNot;

	//OFIter: Iterate == 1 or NOTIterate == 0 **/
	int OFIter;

	/** Use a library with optimal filters (1) or calculate the optimal filter to each pulse (0) **/
	int OFLib;
	
	//Optimal Filter by using the Matched Filter (MF) or the DAB as matched filter (MF, DAB) **/
	char OFInterp[4];

	/** Optimal Filter length Strategy: FREE, BASE2, BYGRADE or FIXED **/
	char OFStrategy[8];

	//Optimal Filter length (taken into account if OFStrategy=FIXED) **/
	int OFLength;

	/** Write intermediate files **/
	int intermediate;
	
	/** Intermediate file **/
	char detectFile[256];
	
	/** File with the optimal filter info **/
	char filterFile[256];
	
	/** Overwrite files? **/
	int clobber;
	
	/** PP's parameter **/
	int maxPulsesPerRecord;
	
	/** Saturation level of the ADC curves **/
	double SaturationValue;

	/** Tstart of the pulses (to be used instead of calculating them if tstartPulse1 =! 0) **/
	int tstartPulse1;
	int tstartPulse2;
	int tstartPulse3;

} ReconstructInitSIRENA;

/** Destructor. */
#ifdef __cplusplus
extern "C"
#endif
void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init);

/** Constructor. Returns a pointer to an empty ReconstructInitSIRENA data structure. */
#ifdef __cplusplus
extern "C"
#endif
ReconstructInitSIRENA* newReconstructInitSIRENA(int* const status);

/** Initializes the structure ReconstructInitSIRENA with the variables required for SIRENA reconstruction */
#ifdef __cplusplus
extern "C"
#endif

void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, char* const record_file, fitsfile *fptr, char* const library_file,
	char* const event_file, double tauFall,	int pulse_length, double scaleFactor, double samplesUp, double nSgms, 
	int mode, double LrsT, double LbT, char* const noise_file, char* pixel_type, char* filter_domain,
	char* filter_method, char* energy_method, int lagsornot, int ofiter, char oflib, char* ofinterp, char* oflength_strategy, int oflength,
	double monoenergy, int interm, char* detectFile, char* filterFile, char clobber, int maxPulsesPerRecord, double SaturationValue,
	int tstartPulse1, int tstartPulse2, int tstartPulse3, int* const status);

/** Constructor. Returns a pointer to an empty PulsesCollection data structure. */
#ifdef __cplusplus
extern "C"
#endif
PulsesCollection* newPulsesCollection(int* const status);

/** Destructor. */
#ifdef __cplusplus
extern "C"
#endif
void freePulsesCollection(PulsesCollection* PulsesColl);

/** Constructor. Returns a pointer to an empty OptimalFilterSIRENA data structure. */
#ifdef __cplusplus
extern "C"
#endif
OptimalFilterSIRENA* newOptimalFilterSIRENA(int* const status);

/** Destructor. */
#ifdef __cplusplus
extern "C"
#endif
void freeOptimalFilterSIRENA(OptimalFilterSIRENA* PulsesColl);

/** Run reconstruction method with an option for SIRENA*/
#ifdef __cplusplus
extern "C"
#endif
void reconstructRecordSIRENA(TesRecord* record,TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init, int lastRecord, int nRecord, PulsesCollection **pulsesAll, OptimalFilterSIRENA **optimalFilter, int* const status);

/** Create and retrieve a LibraryCollection from a file. */
LibraryCollection* getLibraryCollection(const char* const filename, int mode, char *energy_method, char *filter_method, char oflib, char **ofinterp, int* const status);

/** LibraryCollection constructor. Returns a pointer to an empty LibraryCollection data structure. */
LibraryCollection* newLibraryCollection(int* const status);

/** LibraryCollection Destructor. */
void freeLibraryCollection(LibraryCollection** const library_collection);

/** Allocates memory for an LibraryCollection structure. */
void allocateLibraryCollection(LibraryCollection* library_collection,int ntemplates,int* const status);

/** Create and retrieve a NoiseSpec from a file. */
NoiseSpec* getNoiseSpec(const char* const filename,int mode,char *energy_method,char *filter_method,int* const status);

/** NoiseSpec constructor. Returns a pointer to an empty NoiseSpec data structure. */
NoiseSpec* newNoiseSpec(int* const status);

/** Allocates memory for a NoiseSpec structure. */
void allocateNoiseSpec(NoiseSpec* noise_spectrum,int noise_duration,int* const status);

/** Allocates memory for a Pulsetemplate structure. */
void allocatePulseTemplate(PulseTemplate* pulse_template, int template_duration,int* const status);
/** PulseTemplate destructor. */
void freePulsetemplate(PulseTemplate* pulse_template);


/** Allocates memory for a MatchedFilter structure. */
void allocateMatchedFilter(MatchedFilter* matched_filter,int mfilter_duration,int* const status);
/** MatchedFilter destructor. */
void freeMatchedFilter(MatchedFilter* matched_filter);

#endif /* INTEGRASIRENA_H */




