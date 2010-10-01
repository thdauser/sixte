#ifndef PNCCDDETECTOR_H
#define PNCCDDETECTOR_H

#include "sixt.h"
#include "impact.h"
#include "eventfile.h"
#include "pnCCDeventfile.h"
#include "pnCCDevent.h"
#include "genericdetector.h"
#include "squarepixels.h"

/** Number of CCDs of the pnCCD */
//TODO set in the future to 12
#define NumberofCCDs 1

////////////////////////////
// Type Declarations
////////////////////////////

typedef struct {
	GenericDetector generic;
	SquarePixels* pixels;
	
	/** Readout mode of the pnCCD */
	int readout_mode;
	/** Integration time */
	double integration time;

	/** Readout time */
	double readout_time;

	/** Frame time */
	double frame_time;

	/** Number of the current frame */
	long frame;
	
	/** Event list FITS file for the pnCCD-specific events */
	pnCCDEventFile eventlist;

}pnCCDDetector;

struct pnCCDDetectorParameters {
	struct GenericDetectorParameters generic;
	struct SquarePixelsParameters pixels;

	double integration_time;
	double t0;
	char* eventlist_filename;
	char* eventlist_template;
};

//**********************************************************
// Function Declarations
//**********************************************************

/** Set up the configuration of the pnCCDDetector. Also the data 
	structures are initialised. */
int initpnCCDDetector(pnCCDDetector*, struct pnCCDDetectorParameters*);

/** Clean up the pnCCDDetector data structure. Release allocated 
	memory and call clean-up routines of underlying structures. */
int cleanuppnCCDDetector(pnCCDDetector*);

/** A check is performed if a readout has to be performed.
	There are different readout modes of the pnCCD. The current time 
	and the readout time are used to determine if a readout has to be
	performed or not.
	The readout itself is initialised with readoutpnCCDDetector.
	*/
int checkReadoutpnCCDDetector(pnCCDDetector*, double time);

/** The readout of the pnCCD is performed. The different readout 
	modes are implemented (timing mode, full frame mode,...). 
	The measured events are stored in an FITS Eventlist. */
inline int readoutpnCCDDetector(pnCCDDetector*);

/** A photon impact is added to the pnCCDDetector pixels. The generated
	charge is determined according to the detector response. 
	If the charge cloud size is set, split events are calculated
	according to a Gaussian charge cloud model. The new charge is added 
	to the charge already contained in the detector pixel, so pileup
	effects are taken into account. */
int addImpact2pnCCDDdetector(pnCCDDetector*, Impact*);

#endif /* PNCCDDETECTOR_H */ 

