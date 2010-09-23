#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "pnccd_simulation.h"

////////////////////////////
/** Main procedure */
int pnccd_simulation_main() {

	struct Parameters parameters; // containing all program parameters read by PIL

	// Detector data structure (containing the pixel array, width, ...)

	// Background structure

	ImpactListFile impactlistfile;

	int status=EXIT_SUCCESS; // Error status

	// Register HEATOOL:
	set_toolname("pnccd_simulation");
	set_toolversion("0.01");

	if (status = getpar(&parameters)) return(status);

	return(status);
}

static int getpar(struct Parameters* parameters) {

	int status=EXIT_SUCCESS; // error status

	// Get the name of the impact list file (FITS)
	if ((status = PILGetFname("impactlist_filename", parameters->impactlist_filename))) {
		HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
	}

	// Get the readout mode (Timing or FullFrame)
	if ((status = PILGetFname("readout_mode", parameters->impactlist_filename))) {
		HD_ERROR_THROW("Error reading the name of the impact list file!\n", status);
	}
	// Number of readout directions.
	else if ((status = PILGetInt("readout_directions", &parameters->readout_directions))) {
		HD_ERROR_THROW("Error reading the number of readout directions!\n", status);
	}
	// Get the readout time for one detector line.
	else if ((status = PILGetReal("readout_time", &parameters->readout_time))) {
		HD_ERROR_THROW("Error reading the readout time per detector line!\n", status);
	}
	// Get the clear time for one detector line.
	else if ((status = PILGetReal("line_clear_time", &parameters->clear_time))) {
		HD_ERROR_THROW("Error reading the clear time per detector line!\n", status);
	}

	// Detector xwidth [pixel]
	else if ((status = PILGetInt("det_width", &parameters->xwidth))) {
		HD_ERROR_THROW("Error reading the width of the detector!\n", status);
	}
	// Detector ywidth [pixel]
	else if ((status = PILGetInt("det_width", &parameters->xwidth))) {
		HD_ERROR_THROW("Error reading the width of the detector!\n", status);
	}
	// [m]
	else if ((status = PILGetReal("det_pixelwidth", &parameters->pixelwidth))) {
		HD_ERROR_THROW("Error reading the width of the detector pixels!\n", status);
	}

	// [m]
	else if ((status = PILGetReal("ccsigma", &parameters->ccsigma))) {
		HD_ERROR_THROW("Error reading the charge cloud sigma!\n", status);
	}
	if (status) return(status);


	// Read the detector thresholds (either integer PHA or float energy):
	int pha_threshold;
	if ((status = PILGetInt("pha_threshold", &pha_threshold))) {
		HD_ERROR_THROW("Error: could not determine detector PHA threshold!\n", status);
		return(status);
	} else {
		parameters->pha_threshold = (long)pha_threshold;
	}
	if (parameters->pha_threshold==-1) {
		if ((status = PILGetReal4("energy_threshold", &parameters->energy_threshold))) {
			HD_ERROR_THROW("Error: could not determine detector energy threshold!\n", status);
			return(status);
		}
	} else {
		parameters->energy_threshold=0.;
	}

	// Get the name of the detector redistribution file (FITS file)
	if ((status = PILGetFname("rmf_filename", parameters->rmf_filename))) {
		HD_ERROR_THROW("Error reading the name of the detector" 
				"redistribution matrix file (RMF)!\n", status);
	}

	// Get the background count rate
	else if ((status = PILGetReal4("background_rate", &parameters->background_rate))) {
		HD_ERROR_THROW("Error: could not determine the detector background rate!\n", status);
	}

	// Get the name of the output event list (FITS file)
	else if ((status = PILGetFname("eventlist_filename", parameters->eventlist_filename))) {
		HD_ERROR_THROW("Error reading the name of the event list file!\n", status);
	}

	// Get the start time of the simulation
	else if ((status = PILGetReal("t0", &parameters->t0))) {
		HD_ERROR_THROW("Error reading the 't0' parameter!\n", status);
	}

	// Get the timespan for the simulation
	else if ((status = PILGetReal("timespan", &parameters->timespan))) {
		HD_ERROR_THROW("Error reading the 'timespan' parameter!\n", status);
	}

	// Get the name of the FITS template directory.
	// First try to read it from the environment variable.
	// If the variable does not exist, read it from the PIL.
	else { 
		char* buffer;
		if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
			strcpy(parameters->eventlist_template, buffer);
		} else {
			if ((status = PILGetFname("fits_templates", parameters->eventlist_template))) {
				HD_ERROR_THROW("Error reading the path of the FITS templates!\n", status);
			}
		}
	}

	// Set the event list template file for the different WFI modes:
	char template_filename[MAXMSG];
	if (0==strcmp(parameters->readout_mode,"TIMING")) {
		strcpy(template_filename, "pn.timing.eventlist.tpl");
	} else if (0==strcmp(parameters->readout_mode,"FullFrame")) {
		strcpy(template_filename, "pn.fullframe.eventlist.tpl");
	} else {
		status = EXIT_FAILURE;
		char msg[MAXMSG];
		sprintf(msg, "Error: Readout mode %s is not supported!\n", parameters->readout_mode);
		HD_ERROR_THROW(msg, status);
		return(status);
	}
	strcat(strcat(parameters->eventlist_template, "/"), template_filename);

	return(status);
}

