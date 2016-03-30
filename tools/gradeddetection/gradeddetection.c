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


   Copyright 2014 Philippe Peille, IRAP
 */

#include "gradeddetection.h"

////////////////////////////////////
/** Main procedure. */
int gradeddetection_main() {

	// Containing all programm parameters read by PIL.
	struct Parameters par;

	// Advanced detector structure
	AdvDet* det=NULL;

	// Piximpact file
	PixImpFile* piximp_file=NULL;

	// Tes event file
	TesEventFile* event_file=NULL;

	// FITS standard keywords
	SixtStdKeywords* keywords=NULL;

	// Error status.
	int status=EXIT_SUCCESS;

	// Register HEATOOL:
	set_toolname("gradeddetection");
	set_toolversion("0.05");

	do { // Beginning of the ERROR handling loop (will at
		// most be run once).

		// Get program parameters.
		gradeddetection_getpar(&par,&status);
		CHECK_STATUS_BREAK(status);

		headas_chat(3, "initialize ...\n");

		// Initialize the random number generator.
		unsigned int seed=getSeed(par.seed);
		sixt_init_rng(seed, &status);
		CHECK_STATUS_BREAK(status);

		// Open files
		piximp_file=openPixImpFile(par.PixImpList, READONLY,&status);
		CHECK_STATUS_BREAK(status);

		keywords=newSixtStdKeywords(&status);
		sixt_read_fits_stdkeywords(piximp_file->fptr,keywords,&status);
		CHECK_STATUS_BREAK(status);

		event_file = opennewTesEventFile(par.TesEventFile,keywords,par.clobber,&status);
		CHECK_STATUS_BREAK(status);

		// Load advanced detector
		det = loadAdvDet(par.AdvXml,&status);
		CHECK_STATUS_BREAK(status);

		// Load RMF library
		loadRMFLibrary(det,&status);
		CHECK_STATUS_BREAK(status);

		// Load crosstalk if needed
		det->crosstalk_id=par.doCrosstalk;
		if (det->crosstalk_id>0){
			headas_chat(3, "initializing crosstalk ...\n");
			init_crosstalk(det, &status);
			if (status!=EXIT_SUCCESS){
				SIXT_ERROR("failed when initializing crosstalk setup");
				break;
			}
			headas_chat(3, "... done\n");
		}

		// Process impacts
		impactsToEvents(det,piximp_file,event_file,par.saveCrosstalk,&status);


	} while(0); // END of the error handling loop.

	// Free memory
	freePixImpFile(&piximp_file, &status);
	destroyAdvDet(&det);
	freeTesEventFile(event_file,&status);
	freeSixtStdKeywords(keywords);

	// Clean up the random number generator.
	sixt_destroy_rng();

	if (EXIT_SUCCESS==status) {
		headas_chat(3, "finished successfully!\n\n");
		return(EXIT_SUCCESS);
	} else {
		return(EXIT_FAILURE);
	}

}

void gradeddetection_getpar(struct Parameters* const par,int* const status) {
	query_simput_parameter_file_name("PixImpList", &(par->PixImpList), status);
	query_simput_parameter_file_name("AdvXml", &(par->AdvXml), status);
	query_simput_parameter_file_name("TesEventFile", &(par->TesEventFile), status);
	query_simput_parameter_double("tstart", &(par->tstart), status);
	query_simput_parameter_double("tstop", &(par->tstop), status);
	char *buf;
	query_simput_parameter_string("doCrosstalk", &buf, status );
	if (strncmp(buf,"yes",3)==0 ||strncmp(buf,"all",3)==0  ){
		par->doCrosstalk = CROSSTALK_ID_ALL;
	} else if (strncmp(buf,"elec",4)==0){
		par->doCrosstalk = CROSSTALK_ID_ELEC;
	} else if (strncmp(buf,"therm",5)==0){
		par->doCrosstalk = CROSSTALK_ID_THERM;
	} else if (strncmp(buf,"imod",4)==0){
		par->doCrosstalk = CROSSTALK_ID_IMOD;
	} else {
		par->doCrosstalk=CROSSTALK_ID_NONE;
	}
	query_simput_parameter_bool("saveCrosstalk", &(par->saveCrosstalk), status);
	query_simput_parameter_int("seed", &par->seed, status);
	query_simput_parameter_bool("clobber", &par->clobber, status);
	query_simput_parameter_bool("history", &par->history, status);
}


