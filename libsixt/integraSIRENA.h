/***********************************************************************
   This file is part of SIXTE/SIRENA software.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

   Copyright 2014:  INTEGRASIRENA has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P, ESP2014-53672-C3-1-P and RTI2018-096686-B-C21.

***********************************************************************
*                      INTEGRASIRENA
*
*  File:       integraSIRENA.h
*  Developers: Beatriz Cobo
*              cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef INTEGRASIRENA_H
#define INTEGRASIRENA_H 1

#include <stdio.h>

#include "tesrecord.h"
#include "teseventlist.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <time.h>

typedef struct MatrixStruct
{
  /** Number of rows */
  int matrixRows;
  
  /** Number of columns */
  int matrixColumns;
  
  /** Matrix */
  gsl_matrix *matrixBody;
#ifdef __cplusplus
  MatrixStruct();
  MatrixStruct(const MatrixStruct& other);
  MatrixStruct& operator=(const MatrixStruct& other);
  ~MatrixStruct();
#endif

} MatrixStruct;

typedef struct PulseTemplate
{
  /** Duration of the pulse template in the library */
  int template_duration;
  
  /** Vector containing the pulse template */
  gsl_vector *ptemplate;
  
  /** Energy of the template */
  double energy;
  
  /** Pulse Height of the template */
  double pulse_height;
#ifdef __cplusplus
  PulseTemplate();
  PulseTemplate(const PulseTemplate& other);
  PulseTemplate& operator=(const PulseTemplate& other);
  // FIXME
  ~PulseTemplate();
#endif
  
} PulseTemplate;

typedef struct MatchedFilter
{  
  /** Duration of the matched filter in the library */
  int mfilter_duration;
  
  /** Vector containing the matched filter */
  gsl_vector *mfilter;
  
  /** Energy of the mfilter */
  double energy;
  
  /** Pulse Height of the mfilter */
  double pulse_height;
#ifdef __cplusplus
  MatchedFilter();
  MatchedFilter(const MatchedFilter& other);
  MatchedFilter& operator=(const MatchedFilter& other);
  // FIXME
  ~MatchedFilter();
#endif
  
} MatchedFilter;

typedef struct OptimalFilterSIRENA
{
  // Duration of the optimalfilter
  int ofilter_duration;
  
  // Vector containing the optimal filter
  gsl_vector *ofilter;
  
  // Energy of the ofilter
  double energy;
#ifdef __cplusplus
  OptimalFilterSIRENA();
  OptimalFilterSIRENA(const OptimalFilterSIRENA& other);
  OptimalFilterSIRENA& operator=(const OptimalFilterSIRENA& other);
  // FIXME
  ~OptimalFilterSIRENA();
#endif
  
} OptimalFilterSIRENA;

typedef struct LibraryCollection
{
	/** Number of templates & matched filters in the structure. */
	int ntemplates;
        
    /** BASELINE read from the noise file and propagated to the library file**/
    double baseline;
	
	/** Number of fixed length filters in the structure. */
	int nfixedfilters;

    /** Margin to be applied when several energies in the library to choose the proper filter */
    /** To be used in 'find_matchedfilterDAB', 'find_optimalfilterDAB', 'find_prclwn', 'find_Esboundary' and 'find_prclofwm' */
    /** Also included in 'find_matchedfilter' and 'find_optimalfilter' */
    double margin; // %

	/** Energies of the templates */
	gsl_vector *energies;

	/** Pulse Heights of the templates */
	gsl_vector *pulse_heights;

	/** Structure containing all the pulse templates whose length is largeFilter from the library */
	PulseTemplate* pulse_templatesMaxLengthFixedFilter;
	
	/** Structure containing all the pulse templates from the library */
	PulseTemplate* pulse_templates;

	/** Structure containing all the pulse templates filtered & derivated from the library */
	PulseTemplate* pulse_templates_filder;

	/** Maximum of pulse_templates_filder */
	gsl_vector *maxDERs;

	/** 1st sample of pulse_templates_filder */
	gsl_vector *samp1DERs;
	
	/** Structure containing all the pulse templates from the library */
	PulseTemplate* pulse_templates_B0;
	
	/** Structure containing all the matched filters from the library */
	MatchedFilter* matched_filters;

	/** Structure containing all the matched filters from the library */
	MatchedFilter* matched_filters_B0;

	/** Structure containing all the optimal filters from the library */
	OptimalFilterSIRENA* optimal_filters;
	
	/** Structure containing all the fixed optimal filters from the library (FIXFILTF HDU) */
	OptimalFilterSIRENA* optimal_filtersFREQ;

	/** Structure containing all the fixed optimal filters from the library (FIXFILTT HDU) */
	OptimalFilterSIRENA* optimal_filtersTIME;
	
	//MatrixStruct* W;
	gsl_matrix *V;
	gsl_matrix *W;
	
	/** WAB matrix */
	gsl_matrix *WAB;

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

	/** PABMXLFF vector */
	gsl_matrix *PABMXLFF;
	
	/** DAB vector */
	gsl_matrix *DAB;
	
	/** Structure containing all the optimal filters AB from the library */
	OptimalFilterSIRENA* optimal_filtersab;
	
	/** Structure containing all the fixed optimal filters AB in time domain from the library */
	OptimalFilterSIRENA* optimal_filtersabTIME;
	
	/** Structure containing all the fixed optimal filters AB in frequency domain from the library */
	OptimalFilterSIRENA* optimal_filtersabFREQ;

	/** PRECALWN vector */
	gsl_matrix *PRECALWN;
	
	/** PRECALOFWM vector */
	gsl_matrix *PRCLOFWM;
#ifdef __cplusplus
  LibraryCollection();
  LibraryCollection(const LibraryCollection& other);
  LibraryCollection& operator=(const LibraryCollection& other);
  ~LibraryCollection();
#endif
  
} LibraryCollection;

typedef struct NoiseSpec
{
  /** Noise standard deviation */
  double noiseStd;
  
  /** Baseline */
  double baseline;
  
  /** Duration of the noise spectrum */
  int noise_duration;
  
  /** Vector containing the noise spectrum */
  gsl_vector *noisespec;
  
  /** Vector containing the frequecies of the noise spectrum */
  gsl_vector *noisefreqs;
  
  /** Matrix conatining the weight matrixes of the noise for different lengths */
  gsl_matrix *weightMatrixes;
#ifdef __cplusplus
  NoiseSpec();
  NoiseSpec(const NoiseSpec& other);
  NoiseSpec& operator=(const NoiseSpec& other);
  ~NoiseSpec();
#endif

} NoiseSpec;

typedef struct Grading
{
  // Number of grades 
  int ngrades;
  
  // Structure which contains the grading data
  gsl_vector *value;
  gsl_matrix *gradeData;
#ifdef __cplusplus
  Grading();
  Grading(const Grading& other);
  Grading& operator=(const Grading& other);
  ~Grading();
#endif
  
} Grading;// To store the dara read from the XML File

typedef struct I2RData
{
  double I0_START;
  double IMIN;
  double IMAX;
  
  double ADU_CNV;
  double I_BIAS;
  double ADU_BIAS;
  
  double Ifit;

  double V0;
  double RL;
  double L;
  
#ifdef __cplusplus
  I2RData();
  I2RData(const I2RData& other);
  I2RData& operator=(const I2RData& other);
  ~I2RData();
#endif
  
} I2RData;

typedef struct PulseDetected
{
	/** Pulse duration (maximum length) */
	int pulse_duration;

	/** Length of filter used during reconstruction (samples)*/
	int grade1;
      
    /** Distance to previous pulse (samples)*/
    /** tstart(i)-tstart(i-1)*/
	int grade2;

    int grade2_1;

	/** PIX_ID of the detected pulse*/
	int pixid;
        
    ///** PH_IDs corresponding to the record where is the detected pulse*/
	int phid;
    int phid2;
    int phid3;
	
	/** Vector containing the pulse adc values */
	gsl_vector *pulse_adc;
        
    /** Vector containing the pulse adc values and the preBuffer*/
	gsl_vector *pulse_adc_preBuffer;

	/** Start time of the Pulse */
	double Tstart;
        
    /** Start time of the Pulse in samples */
	double TstartSamples;

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

	/** 1st pulse of the filtered-derived pulse */
	double samp1DER;
	
	/** Energy (KeV) of the Pulse */
	double energy;
        
    /** Pulse grade */
	int grading;

    /** Average of the first 4 samples of the derivative of the Pulse */
	double avg_4samplesDerivative;
        
    /** Low resolution energy estimator (4 samples-long filter) */
	double E_lowres;
        
    /** Offset relative to the central point of the parabola */
    double phi;
        
    /** Number of samples shifted to find the maximum of the parabola  */
    int lagsShift;
        
	/** Quality of the Pulse */
	double quality;
        
    /**Number of lags used in detection*/
    int numLagsUsed;
        
    /** Baseline calculated, in general, by using the previous Lb samples to the pulse */
    double bsln;
        
    /** Rms of the baseline calculated, in general, by using the previous Lb samples to the pulse */
    double rmsbsln;
#ifdef __cplusplus
  PulseDetected();
  PulseDetected(const PulseDetected& other);
  PulseDetected& operator=(const PulseDetected& other);
  // FIXME
  ~PulseDetected();
#endif

} PulseDetected;

typedef struct PulsesCollection
{
  /** Number of detected pulses in the structure. **/
  int ndetpulses;
        
  /** Current size of the array **/
  int size;

  /** Array containing all the pulses detetected in record**/
  PulseDetected* pulses_detected;
#ifdef __cplusplus
  PulsesCollection();
  PulsesCollection(const PulsesCollection& other);
  PulsesCollection& operator=(const PulsesCollection& other);
  // FIXME
  ~PulsesCollection();
#endif

} PulsesCollection;

/** Structure containing all the pointers and values to run the reconstruction stage */
typedef struct ReconstructInitSIRENA
{
  /**                   **/
  /** SIRENA parameters **/
  /**                   **/
  /** SIRENA Version **/
  char sirenaVersion[10];
  
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
  
  /** Detection scaleFactor (0.005 ? no filtering) **/
  double scaleFactor;
  
  /** Detection samplesUp (samples to confirm threshold overcoming) **/
  int samplesUp;
  
  /** STC Detection samplesDown (samples below the threshold to look for other pulse) **/
  int samplesDown;
  
  /** Detection nSgms (sigmas to establish a threshold for detection) **/
  double nSgms;
  
  // Detect secondary pulses or not
  int detectSP;

  /** Monochromatic energy for library creation **/
  double monoenergy;
  
  /** Add the PRECALWN HDU in the library file (1) or not (0) **/
  int hduPRECALWN;
  
  /** Add the PRCLOFWM HDU in the library file (1) or not (0) **/
  int hduPRCLOFWM;
  
  /** Length of the longest fixed filter for library creation **/
  int largeFilter;
  
  /** Running sum length for the RS raw energy estimation, in seconds (only in CALIBRATION) **/
  double LrsT;
  
  /** Baseline averaging length for the RS raw energy estimation, in seconds (only in CALIBRATION) **/
  double LbT;
  
  /** Baseline (in ADC units) **/
  //double baseline;
  
  /** Run mode (0: calibration/lib creation  1:energy reconstruction) **/
  int opmode;
  
  /** Detection Mode: AD (Adjusted Derivative) or A1 (Alternative1) **/
  char detectionMode[4];
  
  /** Noise spectrum **/
  NoiseSpec* noise_spectrum;
  
  /** Filtering Domain: T (Time) or F (Frequency) **/
  char FilterDomain[2];
  
  /** Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline) **/
  char FilterMethod[3];
  
  /** Energy Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED or PCA **/
  char EnergyMethod[10];
  
  /** Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R and I2RFITTED) **/
  double filtEev;
  
  /** Constant to apply the I2RFITTED conversion **/
  double Ifit;
  
  //Noise to use with Optimal Filtering: NSD (Noise Spectral Density) or WEIGHTM (weight matrix) **/
  char OFNoise[8];
  
  //LagsOrNot: LAGS == True or NOLAGS == False **/
  int LagsOrNot;
  
  //Number of lags (odd number **/
  int nLags;
  
  //Using 3 lags to analytically calculate a parabola or using 5 lags to fit **/
  int Fitting35;
  
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
  
  //0-padding filter if 0 (from pulseLength to OFLength filter filled in with 0's) or filter with a filter+preBuffer if different from 0
  int preBuffer;
  int preBuffer_min_value;
  int preBuffer_max_value;
  
  // Values from the grading info in the XML file
  int post_min_value;
  int post_max_value;

  /** Write intermediate files **/
  int intermediate;
  
  /** Intermediate file **/
  char detectFile[256];
  
  //Additional error (in samples) added to the detected time"  (Logically, it changes the reconstructed energies) 
  int errorT;
  
  //Sum0Filt: 0-padding: Subtract the sum of the filter (1) or not (0) **/
  int Sum0Filt;
  
  /** Overwrite files? **/
  int clobber;
  
  /** PP's parameter **/
  int maxPulsesPerRecord;
  
  /** Saturation level of the ADC curves **/
  double SaturationValue;
  
  /** Tstart of the pulses (to be used instead of calculating them if tstartPulse1 =! 0) **/
  //int tstartPulse1;
  char tstartPulse1[256];  // Integer number: Sample where the first pulse starts 
                           // or
                           // _nameFile: File where is the tstart (in seconds) of every pulse
  int tstartPulse2;
  int tstartPulse3;
  
  gsl_vector *tstartPulse1_i; // If tstartPulse1=_nameFile, all the tstart in the file are stored in this matrix 
  
  /** Energies for PCA **/
  double energyPCA1;
  double energyPCA2;
  
  /** XML File with instrument definition **/
  char XMLFile[256];
  
  /** grading info read from the XML File **/
  Grading *grading;
  
  I2RData *i2rdata;
#ifdef __cplusplus
  ReconstructInitSIRENA();
  ReconstructInitSIRENA(const ReconstructInitSIRENA& other);
  ~ReconstructInitSIRENA();
  ReconstructInitSIRENA& operator=(const ReconstructInitSIRENA& other);
  ReconstructInitSIRENA* get_threading_object(int n_record);
#endif

} ReconstructInitSIRENA;

#ifdef __cplusplus
extern "C"
#endif
void th_end(ReconstructInitSIRENA* reconstruct_init,
            PulsesCollection** pulsesAll, 
            OptimalFilterSIRENA** optimalFilter);

#ifdef __cplusplus
extern "C"
#endif
int th_get_event_list(TesEventList** test_event, TesRecord** record);

#ifdef __cplusplus
extern "C"
#endif
int is_threading();

#ifdef __cplusplus
extern "C"
#endif
void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init);

#ifdef __cplusplus
extern "C"
#endif
ReconstructInitSIRENA* newReconstructInitSIRENA(int* const status);

#ifdef __cplusplus
extern "C"
#endif

void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, 
                                    char* const record_file, fitsfile *fptr, 
                                    char* const library_file,
                                    char* const event_file,
                                    int pulse_length, double scaleFactor, 
                                    int samplesUp, int samplesDown, 
                                    double nSgms, int detectSP,
                                    int opmode, char* detectionMode,double LrsT, 
                                    double LbT, char* const noise_file, 
                                    char* filter_domain,
                                    char* filter_method, char* energy_method, 
                                    double filtEev, double Ifit, char* ofnoise, 
                                    int lagsornot, int nLags, int Fitting35, int ofiter, char oflib, 
                                    char *ofinterp, char* oflength_strategy, 
                                    int oflength, char preBuffer,
                                    double monoenergy, char hduPRECALWN, 
                                    char hduPRCLOFWM, int largeFilter, 
                                    int interm, char* detectFile, 
                                    int errorT, int Sum0Filt, char clobber, 
                                    int maxPulsesPerRecord, 
                                    double SaturationValue,
                                    //int tstartPulse1, int tstartPulse2, 
                                    char* const tstartPulse1, int tstartPulse2, 
                                    int tstartPulse3, double energyPCA1, 
                                    double energyPCA2, char * const XMLFile, 
                                    int* const status);

#ifdef __cplusplus
extern "C"
#endif
PulsesCollection* newPulsesCollection(int* const status);

#ifdef __cplusplus
extern "C"
#endif
void freePulsesCollection(PulsesCollection* PulsesColl);

#ifdef __cplusplus
extern "C"
#endif
OptimalFilterSIRENA* newOptimalFilterSIRENA(int* const status);

#ifdef __cplusplus
extern "C"
#endif
void freeOptimalFilterSIRENA(OptimalFilterSIRENA* PulsesColl);

#ifdef __cplusplus
extern "C"
#endif
void reconstructRecordSIRENA(TesRecord* record, int trig_reclength, TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init, int lastRecord, int nRecord, PulsesCollection **pulsesAll, OptimalFilterSIRENA **optimalFilter, int* const status);

#ifdef __cplusplus
extern "C"
#endif
void IntegrafreeTesEventListSIRENA(TesEventList* event_list);

LibraryCollection* getLibraryCollection(const char* const filename, int opmode, int hduPRECALWN, int hduPRCLOFWM, int largeFilter, char *filter_domain, int pulse_length, char *energy_method, char *ofnoise, char *filter_method, char oflib, char **ofinterp, double filtEev, int lagsornot, int preBuffer, gsl_vector *pBi, gsl_vector *posti, int* const status);

NoiseSpec* getNoiseSpec(const char* const filename,int opmode,int hduPRCLOFWM,char *energy_method,char *ofnoise,char *filter_method,int* const status);

#ifdef __cplusplus
extern "C"
#endif
void calculateAverageRecord(TesRecord* record, int lastRecord, int *nrecordOK, gsl_vector **averageRecord, int* const status);

#ifdef __cplusplus
extern "C"
#endif
void calculateRecordsError(TesRecord* record, int nrecord, gsl_vector *averageRecord, int* const status);

#endif /* INTEGRASIRENA_H */
