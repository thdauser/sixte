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


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef GENDET_H
#define GENDET_H 1

#include "sixt.h"
#include "badpixmap.h"
#include "clocklist.h"
#include "background.h"
#include "event.h"
#include "eventfile.h"
#include "gendetline.h"
#include "genericdetector.h"
#include "genpixgrid.h"
#include "impact.h"
#include "phabkg.h"
#include "point.h"
#include "rmf.h"
#include "xmlbuffer.h"
#include "mxs.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


typedef enum {
  GENDET_TIME_TRIGGERED =1,
  GENDET_EVENT_TRIGGERED=2
} ReadoutTrigger;


typedef enum {
  GS_NONE       =0,
  GS_GAUSS      =1,
  GS_EXPONENTIAL=2
} GenSplitType;


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic split event generator. */
typedef struct {

  /** Type of the GenSplit model. */
  GenSplitType type;

  /** Split parameter 1. For the Gaussian split model this parameter
      represents the charge cloud size in [m] at a photon energy of 1
      keV. For the exponential model this parameter represents the
      denominator in the exponential distribution term. */
  float par1;

  /** Split parameter 2. This value is only defined for the Gaussian
  charge cloud model. It represents the slope of the function giving
  the charge cloud size vs. photon energy. (ccs(E) = par1 + par2
  *sqrt(E)). */
  float par2;

} GenSplit;

/** characteristics of a DEPFET. If the flag
    depfetflag is set to 1, the detector is treated as
    a DEPFET sensor. If the flag istorageflag is set to
    1 additionally, the DEPFET is treated as an IS-DEPFET. */
typedef struct{

  /** Flag for DEPFET =1 */
  int depfetflag;
  /** Flag for IS-part of DEPFET =1 */
  int istorageflag;

  /** Integration time */
  double t_integration;
  /** Clear time */
  double t_clear;
  /** Settling time */
  double t_settling;

  /** clear function pointer */
  double (*clear_fcn)(double time, double energy, double *constants);
  /** clear function constants */
  double *clear_const;

  /** Transfer time (only for IS-DEPFET) */
  double t_transfer;

}DepfetProp;


/** Generic pixelized X-ray detector. The data structure is designed
    in a generic way. The characteristic properties for a particular
    detector are defined in a detector specific XML file. */
typedef struct {

  /** Detector pixel dimensions. */
  GenPixGrid* pixgrid;

  /** Array of pointers to pixel lines. */
  GenDetLine** line;

  /** Pha2Pi correction filename used to convert PHA to PI accounting
      for signal loss due to lower energy thresholds */
  char* pha2pi_filename;

  /** PI RMF filename of RMF build for pha2pi corrected PI values */
  char* pirmf_filename;

  /** SPEC ARF filename of PSF-Chip coverage calibrated ARF.  */
  char* specarf_filename;


  /** Detector response matrix. The RMF file that is originally loaded
      may also contain ARF contributions. But as these already have
      been taken into account in the generation of the input spectra
      for the X-ray sources, the ARF contributions have to be removed
      by normalizing the RSP matrix. */
  char* rmf_filename;
  struct RMF* rmf;

  /** Lower readout threshold in units of [keV]. This threshold is
      applied in the read-out routine before converting the pixel
      charge to a PHA value. Pixel charges below this threshold will
      be discarded. */
  float threshold_readout_lo_keV;

  /** Lower primary event threshold given in unit of [keV]. This value
      is used in the pattern recognition algorithm to find pixels with
      a sufficient charge for a separate event. */
  float threshold_event_lo_keV;

  /** Lower split threshold given in units of [keV]. This value is
      used in the pattern recognition algorithm. */
  float threshold_split_lo_keV;
  /** Lower split threshold given as a fraction of the charge in the
      main pixel. This value is used in the pattern recognition
      algorithm. */
  float threshold_split_lo_fraction;

  /** Upper threshold for a valid event pattern [keV]. Patterns with a
      total signal above this value are discarded. */
  float threshold_pattern_up_keV;

  /** Split model. */
  GenSplit* split;

  /** Flag, whether a eroBackground model is available for the cosmic
      ray detector background. */
  int auxbackground;

  /** Flag whether to to simulate splits in the AUX background */
  int split_bkg;

  /** Models for detector background based on PHA spectra. */
  PHABkg* phabkg[2];

  /** Flag, whether the detector background models (either
      eROSITA-specific model for the cosmic ray detector background or
      the generic PHA detector background model) should be ignored (if
      they are available). */
  int ignore_bkg;

  /** Flag for detector readout trigger. The readout can be triggered
      either by an incoming photon event (GENDET_EVENT_TRIGGERED) or
      by a timing clock (GENDET_TIME_TRIGGERED). */
  ReadoutTrigger readout_trigger;

  /** List of clock operations for time-triggered detectors. */
  ClockList* clocklist;

  /** Time for one read-out frame [s]. */
  double frametime;

  /** Non-paralyzable dead time that is applied at the readout of each
      pixel. */
  double deadtime;

  /** Flag whether there has been any photon interaction since the
      last new frame. */
  int anyphoton;

  /** Charge transfer efficiency (CTE). In a line shift of the pixel
      array the charges in the shifted pixels are multiplied by this
      value in order to account for losses due to the shift. */
  float cte;

  /** Bad pixel map. */
  BadPixMap* badpixmap;

  /** Output EventFile. Note that this FITS file is not closed when
      the GenDet data struct is destroyed. It has to be closed
      manually. */
  EventFile* elf;

  /** Properties of the DEPFET sensor. */
  DepfetProp depfet;

  /** Properties of the MXS source */
  MXSparams* mxs_params;

  /** MIN and MAX indices of readout indices. */
  int rawymin;
  int rawymax;

} GenDet;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenDet data structure and
    initializes it with the values from the specified XML definition
    file. The second parameter determines, whether the cosmic ray
    detector background model should be activated. */
GenDet* newGenDet(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenDet data structure to NULL. */
void destroyGenDet(GenDet** const det);

/** Function which returns the signal result of a photon impacting
    during the depfet clear. Assumes a linear clear behaviour. */
double depfet_get_linear_clear_signal(double time,
				      double energy,
				      double *constants);

/** Function which returns the signal result of a photon impacting
    during the depfet clear. Assumes an exponential clear behaviour. */
double depfet_get_exponential_clear_signal(double time,
					   double energy,
					   double *constants);

/** Add a new photon impact to the detector. The function return value
    is the number of affected valid detector pixels. */
int addGenDetPhotonImpact(GenDet* const det,
			  const Impact* const impact,
			  int* const status);

/** Operate the time-triggered elements of the GenDet detector up to
    the specified point of time. */
void operateGenDetClock(GenDet* const det,
			const double time,
			int* const status);

/** Set the current detector time in the clocklist to the specified
    value. The default value for the start time is 0. */
void setGenDetStartTime(GenDet* const det, const double t0);

/** Shift the lines of the GenDet detector pixel array by one line
    into the direction of the read-out node in line 0. The charges in
    line 1 are added to the charges in line 0, such that the content
    of line 0 is not lost. The charges of all shifted pixels are
    multiplied with the CTE factor in order to account for possible
    losses. */
void GenDetLineShift(GenDet* const det);

/** Read-out a particular line of the GenDet pixel array and store the
    charges in the output EventFile. After read-out the charges
    in the pixels are deleted. */
void GenDetReadoutLine(GenDet* const det,
		       const int lineindex,
		       const int readoutindex,
		       int* const status);

/** Clear a particular line of the GenDet pixel array. */
void GenDetClearLine(GenDet* const det, const int lineindex);

/** Constructor for GenSplit data structure. */
GenSplit* newGenSplit(int* const status);

/** Destructor for GenSplit data structure. */
void destroyGenSplit(GenSplit** const split);

/** Determine split events for a particular photon impact and add the
    fractional charges to the affected pixels. The function return
    value is the number of valid affected pixels inside the detector,
    where charge has been added to. */
int makeGenSplitEvents(GenDet* const det,
		       const struct Point2d* const impact,
		       const float charge,
		       const long ph_id, const long src_id,
		       const double time,
		       int* const status);

/** Add a charge (photon energy [keV]) to a particular pixel in the
    specified GenDetLine. The routine sets the anycharge flag of the
    affected line. */
void addGenDetCharge2Pixel(GenDet* const line,
			   const int column,
			   const int row,
			   const float signal,
			   const double time,
			   const long ph_id, const long src_id);

/** Parse the GenDet definition from an XML file. */
void parseGenDetXML(GenDet* const det,
		    const char* const filename,
		    int* const status);

/** Assign an output EventFile. */
void setGenDetEventFile(GenDet* const det, EventFile* const elf);

/** Set the ignore_bkg flag. */
void setGenDetIgnoreBkg(GenDet* const det, const int ignore);

/** Adds the signel to a depfet detector */
int addDepfetSignal(GenDet* const det,
		const int colnum,
		const int row,
		const float signal,
		const double time,
		const long ph_id,
		const long src_id);


#endif /* GENDET_H */
