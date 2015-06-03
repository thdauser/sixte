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

////////////////////////////////////////////////////////////////////////
// Constants.
////////////////////////////////////////////////////////////////////////

/** Initial size of RMFLibrary. */
#define RMFLIBRARYSIZE (10)

#define HIRESGRADE 1
#define MEDRESGRADE 2
#define LORESGRADE 3
#define INVGRADE -1

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////

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

/** Data structure describing a pixel with arbitrary geometry.
    Current implementation: only rectangulars, parallel to detector
    coordinate system. */
typedef struct{
  
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

  /** Response matrix file */
  char* rmffile;

  /** ID of the rmf inside general detector (to avoid loading one rmf per pixel) */
  int rmfID;

  /** Ancillary response file */
  char* arffile;

  /** ID of the arf inside general detector (to avoid loading one arf per pixel) */
  int arfID;

}AdvPix;

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

/** Data structure describing the geometry of a pixel detector with
    arbitrary pixel geometry. */
typedef struct{
  
  /** x-shift of detector in respect to the focal point */
  double sx;
  
  /** y-shift of detector in respect to the focal point */
  double sy;
 
  /** Number of pixels. */
  int npix;
  
  /** Counter for operations on pixels */
  int cpix;

  /** array of pixels. */
  AdvPix *pix;

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

}AdvDet;


typedef struct{
	PixImpact *next;
	PixImpact *current;
	PixImpact *previous;
}pixImpPointer;

typedef struct{
	double next,current,previous;
}gradingTimeStruct;

typedef struct {
	gradingTimeStruct *times;
	long row;
}pixGrade;

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

/** Destructore of the AdvPix structure */
void freeAdvPix(AdvPix* pix);

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
void addRMF(AdvDet* det,AdvPix* pixel,int* const status);

/** Destructor of the RMF library structure */
void freeRMFLibrary(RMFLibrary* library);

/** Iterates the different pixels and loads the necessary ARFLibrary */
void loadARFLibrary(AdvDet* det, int* const status);

/** Adds an ARF to the ARF library. The ARF will only be added if it is not already in the library */
void addARF(AdvDet* det,AdvPix* pixel,int* const status);

/** Destructor of the ARF library structure */
void freeARFLibrary(ARFLibrary* library);

/** given grade1 and grade 2, make a decision about the high/mid/los res events **/
int makeGrading(long grade1,long grade2);

/** calculate the grading in samples from the a given impact, and its previous and next impact **/
void calcGradingTimes(double sample_length, gradingTimeStruct pnt,long *grade1, long *grade2, int* status);

/** Process the impacts contained in the piximpacts file with the RMF method */
void processImpactsWithRMF(AdvDet* det,PixImpFile* piximpacfile,TesEventFile* event_file,int* const status);

/** Function to remove overlapping pixels from the detector */
void removeOverlapping(AdvDet* det,int* const status);

#endif /* ADVDET_H */
