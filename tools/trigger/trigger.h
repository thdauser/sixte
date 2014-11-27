/********************************************************************************************
* (C) This product has been produced by IFCA with funding from the Spanish
* Ministry of Science and Innovation (MICINN) under project
* ESP2006-13608-C02-01, as an in-kind contribution to the EURECA
* project.  It remains the property of CSIC (Spanish Council for
* Scientific Research) and it cannot be modified, adapted or used for
* purposes other than the original ones except by its owner.
*
*                  INSTITUTO DE FISICA DE CANTABRIA
*
*                      TRIGGERING PULSES
*
*
*
*  File:      trigger.h
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene Gonzélez Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
*  29/02/08    First version
*  12/03/08    Deleting variable "infile"
* 								Adding input parameter "deltaV"
*  01/04/08	   Adding function "truncate"
* 			   Adding function "findUpperK"
* 			   Adding function "findLowK"
*  			   Adding function "getQuality"
* 			   Adding function "setQuality"
* 			   Adding input parameter "numBitsQuality"
*  22/04/08    Moving functions "getQuality" and "set Quality" to miscellaneous module included in Utils package.
*  29/04/08	   Changing documentation of the module
*  26/06/08	   Adding input parameter "ql"
* 		       Adding output column "sigma"
*  10/09/08	   Adding input parametrer "nameLog"
* 			   Adding input parametrer "verbosity" and processing it using RIL libraries.
* 			   Adding new function "readInputKeywords"
* 			   Deleting in/output keyword "channelcount"
*  29/09/08	   Included "CREATE" and "PROCESS" keyword in the output FITS file
*			   Included new input parameter: writePulse
*			   Included new input parameter: selectRow
*  19/11/08    Included obtainTau function.
*  14/01/09    Included delDuplicated function.
* 			   Included oldTstart global variable.
*			   Included new input keywords
*  02/02/09    New function "Bins2Seconds".
* 			   Resolved bug of baseline Extension write.
*  20/02/09    Included, modificated and removed input keywords
*  24/02/09    Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
*  12/03/09    Changed documentation
*			   Included new columns in EUR-TRG extension: "Baseline and Sigma" of each pulses.
* 			   Used the keywords libraries.
*  14/04/09    The input parameter selectRow is not used
*              Added IVCAL
*  15/07/09	   Deleted non-used variables
*  05/08/10    Deleted functions "findMean2" & "findPulses2"
* 			   Deleted input parameters: "ColumnNameI", "W", "wb", "sizePulse", "dt_before", "dt_after", "k", "k2" & "n2"
*			   New input parameters added: "tauFALL" & "ntaus"
*			   New parameters in "writePulses" function
*			   New parameters in "procSegment" function and renamed ("procRow")
*  08/11/10	   A new parameter in "writePulses" function
*  			   A new parameter in "procRow" function
*  18/02/11	   I0_Filtered deleted in the output FITS file (not used in PULSESHAPE)(writePulses)
*  14/03/11    Deleted non used variables
*  22/03/11    Deleted some non necessary libraries
*  23/03/11	   Added new variable: "pi"
*  25/03/11    Deleted "ql" variable
*  29/03/11    Updated .h
*  04/04/11    "energy" added
*  17/05/11    seconds2Bins instead bins2Seconds
*  19/05/11    vectorFIL not included in procRow
*  25/05/11    New input parameters: samplesUp and nSgms
*              Deleted input parameter n (and n_b)
*              New parameter: safetyMarginTstart
*  09/06/11    "plspolar" added
*  23/06/11    New pulses models library input FITS file ("readLib" function and some parameters added)
*              sign not included in procRow
*  18/10/11    New output keyword "chngplrt"
*  27/10/11    New parameters to handle with code hardpoints:
*               stopCriteriaMKC
*               kappaMKC
*               limitCriteriaMAX
*               nsAftrtstart
*               levelPrvPulse
*               primaryThresholdCriteria
*  07/03/12    New input parameter "scaleFactor"
*  12/03/12    New input parameters "mode", "b_cF" and "c_cF"
*  26/02/12    Deleted variables related to Savitsky-Golay method (methodUsed, nL, nR)
*  14/12/12    Changes in 'procRow' and 'findSePrPulses' to iteratively look for pulses
*  22/01/13    find_model function renamed as 'find_modelOLD'
*              New 'find_model' function in order to interpolate between pulse shapes
*              New 'interpolate_model' function
*  30/01/13    getEnergy function renamed as 'getPulseHeight'
*  05/02/13    New 'Energy' column in the output FITS file where the estimated energy is written
*  14/02/13    Improved the modularity to find pulses and moved the functions used to find pulses to Utils
*              Deleted a no necessary library: <gsl/gsl_sort_vector.h>
*  26/04/13    tAftrtstart = 0
*  02/05/13    'Atstart' to calculate the tstart precisely (pulse broadening due to the filtering)
*  03/07/13	   tAftrtstart != 0
*  01/10/13	   Variables related to PRETRIGS have been deleted or modified because PRETRIGS info is
*              already not used in XRAYCHAIN
*  06/11/13    tAftrtstart has changed from 'constant' to 'input parameter'
*  18/12/13	   New input parameter 'getTaus'
*  24/03/14	   Added a new output extension EUR-TEST to IFCA tests
*  07/07/14    Adapted to PIL/RIL/Common removed dependencies & removed some includes (set in Utils)
*  03/10/14    New parameters for 'writePulses'
*              New functions 'createLibrary' and 'writeLibrary'
*  15/10/14    New column TAILBFR, just in case there is a tail before the first pulse of the row
* ********************************************************************************************/

#ifndef TRIGGER_H_
#define TRIGGER_H_

// ISDC module
// DAL included
#include <isdc.h>

// Utils module
#include <inoutUtils.h>
#include <pulseProcess.h>
#include <keywords.h>

// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_vector.h>

// Sólo para sacar tstart y pulsos de un evento a ficheros diferentes
double timefirstevent;
gsl_vector *tstartevents;
double biasvolt;
gsl_vector *event;

//MC int status = EPOK;
//MC  FOR DAL -> CFITSIO migration
int  colnum=0, felem=0;
char extname[20];
char keyname[10];
char keyvalstr[1000];
char *tform[1];
char *ttype[1];
char *tunit[1];
int evtcnt=0, ttpls1=0, modeval=0,eventcntLib1=0, lib_id=0;

// Provisional => To be deleted in future
int indice = 0;		// Provisional to export the event data to a file!!!!!!!!!!!!!!!!!!!!!!!!!

// Constants

double safetyMarginTstart = 50e-6;
double stopCriteriaMKC = 1.0;  			// Used in medianKappaClipping
                               	   	   	// Given in %
double kappaMKC = 3.0;					// Used in medianKappaClipping
double limitCriteriaMAX = 5.0; 			// To establish the alignment between pulse and model
							            // Given in %
//double tAftrtstart = 50e-6;			// Time fixed as 0 after tstart+tstartDER in the difference between pulse and model
//double tAftrtstart = 0;				// Time fixed as 0 after tstart+tstartDER in the difference between pulse and model
//long nsAftrtstart;      	            // Number of samples fixed as 0 after tstart+tstartDER in the difference between pulse and model
double levelPrvPulse = 100.0;  		    // Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse
double primaryThresholdCriteria = 60.0; // To establish the threshold to locate piled-up pulses
		                                // Given in %

// INPUT FILES

	//MC dal_element *inObject=NULL; 		// Object contains information of input FITS file
	//MC dal_element *inExten=NULL;			// Extension that contains data
	fitsfile *inObject=NULL;
	
	char inName[255];					// Name of the input FITS file

	//MC dal_element *inLibObject=NULL; 		// Object contains information of pulses models library FITS file
	//MC dal_element *inLibExten=NULL;		// Extension that contains pulses models
	fitsfile *inLibObject=NULL;
	fitsfile *inLibExten=NULL;

	char inLibName[255];				// Name of the pulses models library FITS file

	FILE *fileRef;						// Pointer for file which contains errors and warnings.

	FILE * temporalFile;
	char temporalFileName[255];
	FILE * temporalFile2;
	char temporalFileName2[255];

// INPUT KEYWORDS

	struct strctKey *kg,*ki;
	int kgNumber,kiNumber;

	long eventsz;
	double samprate;
	long eventcnt;
	double ivcal;
	double asquid;
	double energy;
	double plspolar;

// INPUT PARAMETERS
	
	// Parameter to decide which operation mode is going to be used (and other parameters to put it into practice)
	int mode;				// Operation mode
	double b_cF;			// Calibration factors (used only in case of normal operation)
	double c_cF;

	// Parameters to apply the running sum filter
	double LrsT;			// Running sum length (in the RS filter case): T -> Time => Seconds
	double LbT;				// Baseline averaging length (in the RS filter case): T -> Time => Seconds
	double Lrs;				// LrsT in samples
	double Lb;				// LbT in samples
	gsl_vector *Bgsl;
	gsl_vector *Lbgsl;

	// Parameters to find pulses
	double tauFALL;			// Fall time of the pulses
	double scaleFactor; 	// Scale factor to apply to the fall time of the pulses in order to calculate the box-car length
	int samplesUp;
	double nSgms;
	double tAftrtstart;			// Time after the alignment between a pulse and its model to be fixed as 0 in the adjusted derivative
								// (in order to help building the adjusted derivative)
								// Time fixed as 0 after tstart+tstartDER in the difference between pulse and model
	long nsAftrtstart;      	// Number of samples fixed as 0 after tstart+tstartDER in the difference between pulse and model

	int numBitsQual;		// Number of bits using for Quality

	// Parameters to establish the tend of the pulses (size of the pulses)
	int ntaus;
	double sizePulse; 		// Size of pulse (seconds) (ntaus*tauFALL)
	int sizePulse_b;		// sizePulse in bins

	gsl_vector *pulseaverage;
	gsl_vector *row_aux;
	double pulseheighttemplate;

	// General Parameters
	int writePulse = 1;		// Write pulses in the output FITS file. Default value = true
	int getTaus = 0;		// Calculate the approximate rise and fall times of the pulses. Default value = false
	char nameLog[255];		// Output log file name
	int verbosity;			// Verbosity level of the output log file

// OUTPUT VECTORS

// AUXILIARY VARIABLES

	//To avoid the deprecate conversions
	char *straux = new char[255];

	bool append;				// Library FITS file new (append=false) or not (append=true)
	long eventcntLib;
	//gsl_matrix *library;		// Energy-EstEnergy
	//gsl_matrix *models;			// Matrix where all the pulse models of the pulse models library are going to be stored
	gsl_vector *energygsl = gsl_vector_alloc(1);
	gsl_vector *estenergygsl = gsl_vector_alloc(1);
	gsl_matrix *pulsesgsl;

	int totalpulses = 1; 		// Total number of pulses processed
	int ntotalrows = 1;
	double initialtime = 0;
	int nPulsesRow = 0;			// Number of pulses of each row

	long nummodels;				// Number of pulse models (rows number) included in the pulses models library file
	gsl_vector *model;			// Pulse which is going to used as model
	                            // (selected row of the PULSE column of the EUR-LIB extension from the pulses models library file)
								// It will be overwritten with the first derivative of it
	gsl_vector *modelSGN;		// Sign of the pulse model first derivative
	long index_maxModel;		// Maximum of the pulse model first derivative (in order to align each pulse with the model)
	gsl_matrix *library;		// Energy-Taurise-Taufall
	gsl_matrix *models;			// Matrix where all the pulse models of the pulse models library are going to be stored

	// Parameters used to write output FITS files
	IOData obj;

// OUTPUT FILE

	//MC dal_element *trgObject = NULL;
	//MC dal_element *trgExten=NULL;
	//MC dal_element *trgExten1=NULL;		// To build other extension EUR-TEST with intermediate results
	fitsfile *trgObject=NULL;
	
	char *str_i= new char[255];
	char *cadena = new char(255);

	char trgName[255];

	char *unit=NULL, *comment=NULL;

// OUTPUT KEYWORDS

	long nqual;				// Number of bits of quality column
	long trg_id;			// Trigger status
	const char *create;		// Name and version of the module: name-0.0.0
	long chngplrt = 0; 		// 1 => Polarity changed (pulses multiplied by -1)
	                        // 0 => Polarity not changed

//FUNCTIONS
	// Parameter manipulation
	int initModule(int argc, char **argv, int status);
	
	// Read the library FITS file
	int readLib (int status);
	int inDataIteratorLib(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct);

	// Create Extension
	int createTriggerFile(int status);

	// Create TRIGGERRepresentation
	int seconds2Bins (int status);
	int procRow(gsl_vector *vectorNOTFIL, gsl_vector *vectorDER, int *npulses, int status);

	int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl, int status);
	int writePulses(gsl_vector *invectorNOTFIL, gsl_vector *invectorDER, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *energy, gsl_matrix **pulses, int status);

	int createLibrary(int status);
	int writeLibrary(double estenergy, gsl_vector **pulsetemplate, int status);

	/*struct my_f_params { double alpha; double tmax;};
	double my_f (double x, void * p)
	{
	     struct my_f_params *params= (struct my_f_params *)p;
	     double alpha = (params->alpha);
	     double tmax = (params->tmax);

	     return  tmax-(1/x)*log((alpha+x)/alpha);
	}
	double my_f_deriv (double x, void *p)
	{
	    struct my_f_parpulsegsl = gsl_matrix_alloc(nPulsesRow,sizePulse_b);ams *params= (struct my_f_params *)p;

	    double alpha = params->alpha;
	    double tmax = params->tmax;

	    return (1/pow(x,2.0))*(log(x/alpha+1.0)-x/(x+alpha));
	}
    void my_f_fdf (double x, void *p, double *y, double *dy)
	{
	    struct my_f_params *params= (struct my_f_params *)p;

	    double alpha = params->alpha;
	    double tmax = params->tmax;

	    *y = tmax-(1/x)*log((alpha+x)/alpha);
	    *dy = (1/pow(x,2.0))*(log(x/alpha+1.0)-x/(x+alpha));
	}*/

	//int obtainTauRise (double Tmax, double Imax, double tauFall, double *tauRise, int status);

	using namespace std;

#endif /*TRIGGER_H_*/
