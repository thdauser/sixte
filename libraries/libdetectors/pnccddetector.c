#include "pnccddetectors.h"

int initpnCCDDetector(pnCCDDetector* pn,
		struct pnCCDDetectorParameters* parameters)
{

	// TODO: make the pnCCDDetector structure

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
	// !!NOTE: Depends on the mode
	pn->integration time = parameters->integration_time;

	// Set the first readout time such that the first readout is performed
	// immediately at the beginning of the simulation 
	// !!NOTE: Should be ok for us
	pn->readout_time = parameters->t0;
	pn->frame = 0;

	// Create and open new event list file
	// TODO write openNEwpnCCDEventFile 

	status = openNewpnCCDEventFile(pn->eventlist, parameters->eventlist_filename, parameters->eventlist_template);
	if (EXIT_SUCCESS!=status) return(status);

	return(status);
}

int getShiftedRaw(pnccdDetector* pn, Impact* impact){

	return status;
}
