#include "pnccddetectors.h"

int initpnCCDDetector(pnCCDDetector* pn,
		struct pnCCDDetectorParameters* parameters)
{
	int status = EXIT_SUCCESS;

	// Call the initialistation routines of the underlying data structures.
	status = initGenericDetector(&pn->generic, &parameters->generic);
	if (EXIT_SUCCESS!=status) return(status);

	// CCDindex can be used for the initialization of the 12 single CCDs 
	// of the pnCCD (at the current stage just one CCD is initialized, due 
	// to the fact, that the timing mode works wirh just one CCD (64x200)
	int ccdindex;
	for(ccdindex=0; ccdindex<NumberofCCDs; ccdindex++) {
		status = initSquarePixels(&pn->pixels[ccdindex], &parameters->pixels);
		if (EXIT_SUCCESS!=status) return(status);
	}

	// Set up the pnCCD configuration
	pn->integration time = parameters->integration_time;
	// !!NOTE: Depends on the mode!!!

	// Set the first readout time such that the first readout is performed
	// immediately at the beginning of the simulation 
	pn->readout_time = parameters->t0;
	pn->frame = 0;
	// !!NOTE: Should be ok for us

	// Create and open new event list file
	// TODO write openNewpnCCDEventFile 

	status = openNewpnCCDEventFile(pn->eventlist, parameters->eventlist_filename, parameters->eventlist_template);
	if (EXIT_SUCCESS!=status) return(status);

	return(status);
}

int cleanuppnCCDDetector(pnCCDDetector* pn)
{
	int status=EXIT_SUCCESS;

	// Call the cleanup routines of the underlying data structures
	int ccdindex;
	for(ccdindex=0; ccdindex<NumberofCCDs; ccdindex++) {
		cleanupSquarePixels(&pn->pixels[ccdindex]);
	}
	cleanupGenericDetector(&pn->generic);
	// TODO write closepnCCDEventFile
	status+=closepnCCDEventFile(&pn->eventlist);

	return(status);
}

int checkReadoutpnCCDDetector(pnCCDDetector* pn, double time)
{
	int status = EXIT_SUCCESS;

	// Check, if the detector integration time was exceeded.
	// In that case, read out the detector.
	if (time > pn->readout_time) {
		// The detector should be read out now!!

		//Readout the detector and create eventlist entries for the actual time
		status = readoutpnCCDDetector(pn);
		// TODO: Write readoutpnCCDDetector!!!
		if (EXIT_SUCCESS!=status) return(status);

		// Clear the detector array
		int ccdindex;
		for(ccdindex=0; ccdindex<NumberofCCDs; ccdindex++) {
			clearSquarePixels(&pn->pixels[ccdindex]);
		}

		// !!! NOTE: Is this ok for us? Why should this only be
		// !!! 		 valid for the eROSITA Detectors?
		// !!!		 Maybe ask Christian
		while (time > pn->readout_time) {
			pn->readout_time += pn->inegration_time;
			pn->frame += 1;
		}

		// Print the time of the current events
		headas_printf("\rtime: %.3lf s ", pn->readout_time);
		fflush(NULL);
	}

	return(status);
}

inline int readoutpnCCDDetector(pnCCDDetector* pn)
{
	int ccdindex;
	int status = EXIT_SUCCESS;
	
	// TODO: together with Michi implement the readout modes
	// read out just one CCD at the moment (CCD0 of Quadrant 1)
	// TODO: first implement the timing mode, afterwards the full frame mode

	return(status);
}


int addImpact2pnCCDDdetector(pnCCDDetector* pn, Impact* impact)
{
	int status = EXIT_SUCCESS;

	// Check if there has to be a read out before adding the new Impact 
	// to the detector
	// !!! Note: in timing mode the charge is continously shifted 
	// !!! 		 the check of the readout has to be implemented 
	// !!!		 carefully 										  !!!
	status=checkReadoutpnCCDDetector(pn, impact->time);

	// Determine a detector channel (PHA channel) according to the RMF.
	// The channel is obtained from the RMF using the corresponding
	// HEAdas routine which is based on drawing a random number
	long channel;
	ReturnChannel(pn->generic.rmf, impact->energy, 1, &channel);

	// Check if the photon is really measured. If the
	// PHA channel returned by the HEAdas RMF function is '-1', 
	// the photon is not detected.
	// This can happen as the RMF is actually an RSP including, e.g., 
	// the detector quantum efficiency and filter transmission.
	if (-1==channel) {
		return(status); // Break the function (and continue with the next photon).
	}
	assert(channel>=0);

	// Get the corresponding created charge.
	// NOTE: In this simulation the charge is represented by the nominal
	// photon energy which corresponds to the PHA channel according to the
	// EBOUNDS table.
	float charge = getEnergy(channel, pn->generic.rmf);

	if (charge > 0.) {
		int x[4], y[4];
		double fraction[4];

		// get the impact point in detector coordinates, check which ccds 
		// is affected

		// TODO: write getCCD 
		// !!! Note: Do we really need getCCD? 
		// 			 We can do anything with the Detector coordinates
		// 			 Maybe it is nicer and more obvious to work with 
		// 			 single CCDs
		int ccd = getCCD(pn, impact->position, x, y);

		int npixels = getSquarePixelsGaussianSplits(&pn->pixels[ccd], &(pn->generic.gcc), impact->position, x, y, fraction);
		
		// Add the charge created by the photon to the affected detector pixels.
		int count;
		for (count=0; count<npixels; count++) {
			if (x[count] != INVALID_PIXEL) {
				pn->pixels[ccdindex].array[x[count]][y[count]].charge += 
					charge * fraction[count];
			}
		}
	}

	return(status);
}

int getCCD(pnCCDDetector* pn, Impact* impact, int *x, int *y);
{
		int status = EXIT_SUCCESS;

		return(status);
}


int getShiftedRaw(pnCCDDetector* pn, Impact* impact)
{
	int status = EXIT_SUCCESS;

	return status;
}
