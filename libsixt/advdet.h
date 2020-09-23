/*
   This file is part of SIXTE.

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


   Copyright 2014 Thorsten Brand, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef ADVDET_H
#define ADVDET_H 1

#include "sixt.h"
#include "impact.h"
#include "pixelimpact.h"
#include "point.h"
#include "xmlbuffer.h"
#include "teseventlist.h"
#include "pixelimpactfile.h"
#include "tespixel.h"

// For FDM calculations in tessim
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <unistd.h>

////////////////////////////////////////////////////////////////////////
// Constants.
////////////////////////////////////////////////////////////////////////

/** Initial size of RMFLibrary. */
#define RMFLIBRARYSIZE (10)

#define DEFAULTGOODSAMPLE 32768 // TODO : check whether we do want to leave that as a macro

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


typedef struct MatrixCrossTalk MatrixCrossTalk;
typedef struct MatrixEnerdepCrossTalk MatrixEnerdepCrossTalk;
typedef struct MatrixPropCrossTalk MatrixPropCrossTalk;
typedef struct MatrixDerCrossTalk MatrixDerCrossTalk;
typedef struct Channel Channel;
typedef struct IntermodulationCrossTalk IntermodulationCrossTalk;
typedef struct ReadoutChannels ReadoutChannels;
typedef struct FDMSystem FDMSystem;

/** Data structure describing the noise properties of calorimeter
 pixels */
typedef struct{

  /** White noise RMS value */
  double WhiteRMS;

  /** Normalisation of the filter function */
  double H0;

  /** Number of zeros */
  int Nz;

  /** Zeros */
  double *Zeros;

  /** Number of poles */
  int Np;

  /** Poles */
  double *Poles;

  /** 1/f noise properties */
  /** 1/f rms value */
  double OoFRMS;

  /** 1/f knee frequency (e.g. the pole where the flat power spectrum
      turns into a 1/f curve) */
  double OoFKnee;

}TESNoiseProperties;

/** Data structure defining a grading */

typedef struct{

  /** Response matrix file */
  char* rmffile;

  /** ID of the rmf */
  struct RMF* rmf;

  /** The grade values */
  int value;

  /** Size in samples before and after one pulse defining the grade */
  long gradelim_pre;
  long gradelim_post;

}TESGrade;

/** Data structure describing a pixel with arbitrary geometry.
    Current implementation: only rectangulars, parallel to detector
    coordinate system. */
struct AdvPix{

  /** x-shift of pixel in respect to the detector reference point */
  double sx;

  /** y-shift of pixel in respect to the detector reference point */
  double sy;

  /** Width of the pixel [m], centered around pixel reference point. */
  double width;

  /** Height of the pixel [m], centered around pixel reference point. */
  double height;

  /** Index of pixel in detector structure. */
  int pindex;

  /** Version of pulse template. */
  char version[10];

  /** Version index of pulse template */
  int profVersionID;

  /** Sampling frequency */
  int ADCOffset;

  /** Sampling frequency */
  double calfactor;

  /** Noise properties of pixel noise */
  TESNoiseProperties* TESNoise;

  /** Number of grades for this pixel */
  int ngrades;

  /** Different grades for this pixel */
  TESGrade* grades;

  /** Signal if grading has already been defined globally */
  int global_grading;

  /** Ancillary response file */
  char* arffile;

  /** ID of the arf inside general detector (to avoid loading one arf per pixel) */
  struct ARF* arf;

  /** Frequency of the pixel */
  double freq;

  /** Resonance frequency of the pixel */
  /** equals freq unless specified otherwise (i.e. in case of detuned pixels) */
  double resfreq;

  /** Read-out channel to which it belongs */
  Channel* channel;

  /** Cross-talk structures */
  MatrixCrossTalk* thermal_cross_talk;

  /**Electrical cross_talk */
  MatrixEnerdepCrossTalk* electrical_cross_talk;

  /** Proportional cross talk */
  MatrixPropCrossTalk* prop_cross_talk;

  /** Derivative cross talk */
  MatrixDerCrossTalk* der_cross_talk;

  /** tessim struct */
  char* tes_type;
  tesparams *tes;

  /** TDM parameters */
  int row;

}; typedef struct AdvPix AdvPix;

struct AdvRecons{       //SIRENA

  /** Number of grades for this pixel */
  int ngrades;

  /** Different grades for this pixel */
  TESGrade* grades;

}; typedef struct AdvRecons AdvRecons;

/** Data structure containing a library of different RMFs */
typedef struct{

	/** Current size of the allocated library */
	int size;

	/** Number of RMFs in the library */
	int n_rmf;

	/** Array containing the filenames of the loaded rmfs */
	char** filenames;

	/** Array containing the rmf structures */
	struct RMF** rmf_array;

}RMFLibrary;

/** Data structure containing a library of different ARFs */
typedef struct{

	/** Current size of the allocated library */
	int size;

	/** Number of ARFs in the library */
	int n_arf;

	/** Array containing the filenames of the loaded arfs */
	char** filenames;

	/** Array containing the arf structures */
	struct ARF** arf_array;

}ARFLibrary;

/** structure defining the time dependent weights for the crosstalk*/
typedef struct{
	int length;
	char* name_type; // Useless for the moment
	double* time;
	double* weight;
	double weight_t0;  // weight at time t=0 (simultaneous)
} CrosstalkTimedep;

/** structure containing one of the intermodulation crosstalk table*/
typedef struct{

	double**** matrix; // 4d table containing the frequency weights

	double* ampl;
	double* dt; // in seconds
	double* freq;

	int n_ampl;
	int n_dt;
	int n_freq;

	// an event is defined to extent and cause crosstalk from [dt_min,dt_max]
	double dt_min; // in seconds
	double dt_max; // in seconds

}ImodTab;

/** structure containing one of the electrical crosstalk */
typedef struct{

	double*** matrix; // 3d table containing the frequency weights

	double* freq_s;
	double* freq_p;
	double* ener_p;

	int n_freq_s;
	int n_freq_p;
	int n_ener_p;

}ElecTab;

/** structure containing the TDM cross-talk */
typedef struct{

	double*** matrix; // 3d table containing the weights

	double* ener_p;
	double* ener_v;
	double* samples;

	int n_samples;
	int n_ener_p;
	int n_ener_v;

}TDMTab;

/** Data structure describing the geometry of a pixel detector with
    arbitrary pixel geometry. */
typedef struct{

  /** x-shift of detector in respect to the focal point */
  double sx;

  /** y-shift of detector in respect to the focal point */
  double sy;

  /** Number of pixels. */
  int npix;
  
  /** 1 if 'reconstruction' node exits (0 if not) */    //SIRENA
  int nrecons;

  /** Counter for operations on pixels */
  int cpix;

  /** array of pixels. */
  AdvPix *pix;
  
  /** reconstruction info. */                           //SIRENA
  AdvRecons *recons;

  /** File name (without path contributions) of the FITS file
      containing the XML detector definition. */
  char* filename;

  /** Path to the FITS file containing the XML detector definition. */
  char* filepath;

  /** Name of file of pulse templates. */
  char tesproffilename[MAXFILENAME];

  /** Sampling frequency */
  double SampleFreq;

  /** Signal if being inside the 'tesnoisefilter' tag of xml */
  int tesnoisefilter;

  /** Signal if being inside the 'pixel' tag of xml */
  int inpixel;

  /** Signal if the OofNoise is requested in at least one pixel */
  int oof_activated;

  /** RMF library */
  RMFLibrary* rmf_library;

  /** ARF library */
  ARFLibrary* arf_library;

  /** TDM or not TDM */
  int tdm;
  int max_rows;

  /** File listing for each pixel the channel and frequency */
  char* channel_file;

  /** File listing for each pixel the channel and resonance frequency */
  /** Same format as channel_file */
  char* channel_resfreq_file;

  /** List of all readout channels */
  ReadoutChannels* readout_channels;

  /** Structure containing the intermodulation crosstalk table and the time dependence crosstalk table */

  char* crosstalk_intermod_file;
  ImodTab* crosstalk_imod_table;

  /** Structure containing the thermal crosstalk informations */

  char** crosstalk_thermal_timedep_file;
  CrosstalkTimedep** crosstalk_ther_timedep;
  int xt_num_thermal;
  double* xt_dist_thermal;
  double* xt_weight_thermal;

  /** Structure containing the electrical crosstalk table and the dependence crosstalk table  */
  char* crosstalk_elec_file;
  char* crosstalk_elec_timedep_file;
  ElecTab* crosstalk_elec;
  CrosstalkTimedep* crosstalk_elec_timedep;
  float scaling;

  /** Structure containing the proportional TDM cross-talk tables */
  char* TDM_prop_file;
  TDMTab* crosstalk_TDM_prop;
  float prop_TDM_scaling_1;
  float prop_TDM_scaling_2;

  /** Structure containing the derivative TDM cross-talk tables */
  char* TDM_der_file;
  TDMTab* crosstalk_TDM_der;
  float der_TDM_scaling;

  /** Trigger threshold */
  double threshold_event_lo_keV;

  /** Crosstalk ID (which effects should be included */
  int crosstalk_id;

  /** File containing TES pixel definitions */
  char* tes_type_file;

  /** Common impedance and capacitance for FDM Crosstalk */
  double L_Common;
  double C_Common;

}AdvDet;


/////////////////////////////////////////////////////////////////////
// Structures for the Crosstalk
/////////////////////////////////////////////////////////////////////

/** Structure defining the cross talk between pixels, which can be approximated
    by a simple matrix containing weights */
struct MatrixCrossTalk{
	/** number of cross-talk pixels */
	int num_cross_talk_pixels;

	/** Distance to neighbour type*/
	int* cross_talk_index;

	/** Array containing cross-talk pixels */
	AdvPix** cross_talk_pixels;

	/** Cross-talk weights*/
	double* cross_talk_weights;

};

/** Structure defining the cross talk between pixels, which can be approximated
    by a simple matrix containing weights */
struct MatrixEnerdepCrossTalk{
	/** number of cross-talk pixels */
	int num_cross_talk_pixels;

	/** Array containing cross-talk pixels */
	AdvPix** cross_talk_pixels;

	/** Cross-talk weights*/
	double** cross_talk_weights;

	int n_ener;

};

/** Structure defining the proporitonal cross talk between pixels, which can be approximated
    by a simple matrix containing weights */
struct MatrixPropCrossTalk{

	/** Number of cross-talk pixel for type 1 cross-talk */
	int type_1_pix;

	/** Number of cross-talk pixel for type 2 cross-talk */
	int type_2_pix;

	/** Array containing the cross talk pixels */
	int* cross_talk_pixels_1;
	int* cross_talk_pixels_2;

};

/** Structure defining the proporitonal cross talk between pixels, which can be approximated
    by a simple matrix containing weights */
struct MatrixDerCrossTalk{

	/** Number of cross-talk pixel for type 3 cross-talk */
	int num_cross_talk_pixels;

	/** Array containing the cross talk pixels */
	int* cross_talk_pixels;
};


/** Structure defining the cross talk between pixels, which can be approximated
    by a simple matrix containing weights */
struct IntermodulationCrossTalk{
	/** number of cross-talk pixels */
	int num_cross_talk_pixels;

	/** numbrer of combiniations for each pixel */
	//	int* num_pixel_combinations;

	/** Array containing cross-talk pixels */
	AdvPix*** cross_talk_pixels;

	/** Cross-talk weights*/
	ImodTab** cross_talk_info;

};

/** Structure containing linear equations system for FDM Crosstalk **/
struct FDMSystem{
        int num_pixels;
        double* omega_array; // angular frequencies for pixels, length = num_pixels
        double* res_omega_array; // angular resonance frequencies for pixels, length = num_pixels
        double** Z_array; // complex impedances (saved as doubles, since they are imaginary), size = num_pixels x num_pixels
        double L_Common; // common impedance [H]
        double C_Common; // common capacitance [F]
        double* X_L; // imaginary impedance term that's constant and used a lot
        double* capFac; // extra factor due to common capacitance - only for n = 2!
        };


/** Structure of a single channel, including all its contained pixels*/
struct Channel{
	AdvPix** pixels;
	int num_pixels;
	int channel_id;
    FDMSystem* fdmsys;
};

/** Structure containing a certain amount of Impacts, which are not written
    to the output file yet, as further modifications due to additional
    (crosstalk) events might happen or that these events directly influence
    other events */
typedef struct channelImpacts{
	double time_diff_criteria;
	Impact** impacts_saved; // length = num_pixels
	Channel* channel;
} channelImpacts;

/** Structure combining all readout channels */
struct ReadoutChannels{
	Channel* channels;
	double* df_information_band; //in [Hz];
	int num_channels;
};

/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

/** Constructor. Allocates memory for a new TESNoiseProperties data
    structure. */
TESNoiseProperties* newTESNoise(int* const status);

/** Allocates zeros and poles arrays. */
void allocateZerosPoles(TESNoiseProperties* noise,int nzeros,int npoles,int* const status);

/** Duplicates a TESNoiseProperties structure into another. If given structure is NULL, allocate new one. */
TESNoiseProperties* duplicateTESNoise(TESNoiseProperties* noise,int nzeros,int npoles,int* const status);

/** Destructor. Releases all allocated memory. */
void destroyTESNoiseProperties(TESNoiseProperties* noise);

/** Constructor for an empty AdvPix structure */
AdvPix* newAdvPix(int* const status);

/** Destructor of the AdvPix structure */
void freeAdvPix(AdvPix* pix);

/** Constructor for an empty AdvRecons structure */
AdvRecons* newAdvRecons(int* const status);      //SIRENA

/** Destructor of the AdvRecons structure */
void freeAdvRecons(AdvRecons* recons);   //SIRENA

/** Remove the existing grading scheme from the pixel */
void freeGrading(AdvPix* pix);

/** Remove the existing grading scheme from the pixel */
void freeGradingRecons(AdvRecons* recons);      //SIRENA

/** Frees electrical table*/
void freeElecTab(ElecTab* tab, int gr);

/** Free the Readout Channel Structure */
void freeReadoutChannels(ReadoutChannels* rc);

/** Free the FDM System Structure */
void freeFDMSystem(FDMSystem* fdm_sys);

/** Read the advanced detector syntax from the specified XML */
void parseAdvDetXML(AdvDet* const det,
	       const char* const filename,
	       int* const status);

/** Constructor. Allocates memory for a new AdvDet data structure. */
AdvDet* newAdvDet(int* const status);

/** Loads AdvDet information from XML file */
AdvDet* loadAdvDet(const char* const filename,
		     int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the AdvDet data structure to NULL. */
void destroyAdvDet(AdvDet **det);

/** Function testing if an impact lies inside a pixel. */
int CheckAdvPixImpact(AdvPix pix, Impact *imp);

/** Function calculating the exact impact position and time in
    pixel coordinates. */
void CalcAdvPixImpact(AdvPix pix, Impact *imp, PixImpact *piximp);

/** Function determining the pixel indices which have an impact from one
    event. Gives the number of pixels that were hit.*/
int AdvImpactList(AdvDet *det, Impact *imp, PixImpact **piximp);

/** Iterates the different pixels and loads the necessary RMFLibrary */
void loadRMFLibrary(AdvDet* det, int* const status);

/** Adds an RMF to the RMF library. The RMF will only be added if it is not already in the library */
void addRMF(AdvDet* det,AdvPix* pixel,int rmf_index,int* const status);

/** Destructor of the RMF library structure */
void freeRMFLibrary(RMFLibrary* library);

/** Iterates the different pixels and loads the necessary ARFLibrary */
void loadARFLibrary(AdvDet* det, int* const status);

/** Adds an ARF to the ARF library. The ARF will only be added if it is not already in the library */
void addARF(AdvDet* det,AdvPix* pixel,int* const status);

/** Destructor of the ARF library structure */
void freeARFLibrary(ARFLibrary* library);

/** Function to remove overlapping pixels from the detector */
void removeOverlapping(AdvDet* det,int* const status);

/** Constructor for MatrixCrossTalk structure */
MatrixCrossTalk* newMatrixCrossTalk(int* const status);

/** Constructor for MatrixEnerdepCrossTalk structure */
MatrixEnerdepCrossTalk* newMatrixEnerdepCrossTalk(int grade, int* const status);

/** Constructor for MatrixPropCrossTalk structure */
MatrixPropCrossTalk* newMatrixPropCrossTalk(int* const status);

/** Constructor for MatrixDerCrossTalk structure */
MatrixDerCrossTalk* newMatrixDerCrossTalk(int* const status);

/** Destructor for MatrixCrossTalk structure */
void freeMatrixCrossTalk(MatrixCrossTalk* matrix);

/** Destructor for MatrixEnerdepCrossTalk structure */
void freeMatrixEnerdepCrossTalk(MatrixEnerdepCrossTalk* matrix, int gr);

/** Destructor for MatrixPropCrossTalk structure */
void freeMatrixPropCrossTalk(MatrixPropCrossTalk* matrix);

/** Destructor for MatrixDerCrossTalk structure */
void freeMatrixDerCrossTalk(MatrixDerCrossTalk* matrix);

/** Constructor for IntermodulationCrossTalk structure */
IntermodulationCrossTalk* newImodCrossTalk(int* const status);

/** Destructor for IntermodulationCrossTalk structure */
void freeImodCrossTalk(IntermodulationCrossTalk** matrix);

/** Constructor for CrosstalkTimdep structure */
CrosstalkTimedep* newCrossTalkTimedep(int* const status);

/** Destructor for CrosstalkTimdep structure */
void freeCrosstalkTimedep(CrosstalkTimedep* timedep, int gr);

/** free Intermod Table and Strcture */
void freeImodTab(ImodTab* tab, int gr);

/** free TDM Table and Structures */
void freeTDMTab(TDMTab* tab, int gr);

/** free the crosstalk structures */
void freeCrosstalk(AdvDet* det, int gr);

#endif /* ADVDET_H */
