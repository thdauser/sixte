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
	SquarePixels pixels[NumberofCCDs];
	
	/** Integration time */
	double integration time;

	/** Readout time */
	double readout_time;

	/** Number of the current frame */
	long frame;
	
	/** Event list FITS file for the pnCCD-specific events */
	pnCCDEventFile eventlist;

}pnCCDDetector;

