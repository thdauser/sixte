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
 */

#include "xifupipeline.h"


int xifupipeline_main()
{
	// Program parameters.
	struct Parameters par;

	// Tes general parameters
	TESGeneralParameters genpar;

	// Instrument setup.
	GenInst* inst=NULL;

	// Attitude.
	Attitude* ac=NULL;

	// GTI collection.
	GTI* gti=NULL;

	// Catalog of input X-ray sources.
	SourceCatalog* srccat[MAX_N_SIMPUT];
	unsigned long ii;
	for (ii=0; ii<MAX_N_SIMPUT; ii++) {
		srccat[ii]=NULL;
	}

	// Photon list file.
	PhotonFile* plf=NULL;

	// Impact list file.
	ImpactFile* ilf=NULL;

	// Piximpact file.
	PixImpFile* pixilf=NULL;

	// Event list file
	TesEventFile* event_file=NULL;

	// FITS standard keywords
	SixtStdKeywords* keywords=NULL;

	// Advanced detector structure
	AdvDet* det=NULL;

	// Output file for progress status.
	FILE* progressfile=NULL;

	// Tes simulation initialization structure
	TESInitStruct* init=NULL;

	// Pulse reconstruction initialization structure
	ReconstructInit* reconstruct_init = NULL;

	// Error status.
	int status=EXIT_SUCCESS;


	// Register HEATOOL
	set_toolname("xifupipeline");
	set_toolversion("0.06");


	do { // Beginning of ERROR HANDLING Loop.

		// ---- Initialization ----

		// Read the parameters using PIL.
		status=xifupipeline_getpar(&par);
		CHECK_STATUS_BREAK(status);

		headas_chat(3, "initialize ...\n");

		// Determine the prefix for the output files.
		char ucase_buffer[MAXFILENAME]={""};
		strcpy(ucase_buffer, par.Prefix);
		strtoupper(ucase_buffer);
		if (0==strcmp(ucase_buffer,"NONE")) {
			strcpy(par.Prefix, "");
		}

		// Determine the photon list output file.
		char photonlist_filename[MAXFILENAME]={""};
		strcpy(ucase_buffer, par.PhotonList);
		strtoupper(ucase_buffer);
		if (0==strcmp(ucase_buffer,"NONE")) {
		} else {
			strcpy(photonlist_filename, par.Prefix);
			strcat(photonlist_filename, par.PhotonList);
		}

		// Determine the impact list output file.
		char impactlist_filename[MAXFILENAME]={""};
		strcpy(ucase_buffer, par.ImpactList);
		strtoupper(ucase_buffer);
		if (0==strcmp(ucase_buffer,"NONE")) {
			strcpy(impactlist_filename, "");
		} else {
			strcpy(impactlist_filename, par.Prefix);
			strcat(impactlist_filename, par.ImpactList);
		}

		// Determine the piximpact list output file.
		char piximpactlist_filename[MAXFILENAME]={""};
		strcpy(ucase_buffer, par.PixImpactList);
		strtoupper(ucase_buffer);
		int delete_rawdata=0;
		if (0==strcmp(ucase_buffer,"NONE")) {
			delete_rawdata=1;
			strcpy(piximpactlist_filename, par.Prefix);
			strcat(piximpactlist_filename, "piximpact.fits");
		} else {
			strcpy(piximpactlist_filename, par.Prefix);
			strcat(piximpactlist_filename, par.PixImpactList);
		}
		strcpy(par.PixImpactList,piximpactlist_filename);

		//Determine the tes record output file.
		if (! par.UseRMF) {
			char tesrecord_filename[MAXFILENAME];
			strcpy(ucase_buffer, par.TesTriggerFile);
			strtoupper(ucase_buffer);
			if (0==strcmp(ucase_buffer,"NONE")) {
				strcpy(tesrecord_filename, "");
			} else {
				strcpy(tesrecord_filename, par.Prefix);
				strcat(tesrecord_filename, par.TesTriggerFile);
			}
			strcpy(par.TesTriggerFile,tesrecord_filename);
//			free(tesrecord_filename);
		} else {
			strcpy(par.TesTriggerFile, "");
		}


		//Determine the event output file.
		char evtlist_filename[MAXFILENAME];
		strcpy(ucase_buffer, par.EvtFile);
		strtoupper(ucase_buffer);
		if (0==strcmp(ucase_buffer,"NONE")) {
			strcpy(evtlist_filename, par.Prefix);
			strcat(evtlist_filename, "evt.fits");
		} else {
			strcpy(evtlist_filename, par.Prefix);
			strcat(evtlist_filename, par.EvtFile);
		}
		strcpy(par.EvtFile,evtlist_filename);

		// Initialize the random number generator.
		unsigned int seed=getSeed(par.Seed);
		sixt_init_rng(seed, &status);
		CHECK_STATUS_BREAK(status);

		// Set the progress status output file.
		strcpy(ucase_buffer, par.ProgressFile);
		strtoupper(ucase_buffer);
		if (0!=strcmp(ucase_buffer, "STDOUT")) {
			progressfile=fopen(par.ProgressFile, "w+");
			char msg[MAXMSG];
			sprintf(msg, "could not open file '%s' for output of progress status",
					par.ProgressFile);
			CHECK_NULL_BREAK(progressfile, status, msg);
		}

		// Load the instrument configuration.
		inst=loadGenInst(par.XMLFile, seed, &status);
		CHECK_STATUS_BREAK(status);

		// Set the usage of the detector background according to
		// the respective program parameter.
		setGenDetIgnoreBkg(inst->det, !par.Background);

		// Set up the Attitude.
		strcpy(ucase_buffer, par.Attitude);
		strtoupper(ucase_buffer);
		if ((strlen(par.Attitude)==0)||(0==strcmp(ucase_buffer, "NONE"))) {
			// Set up a pointing attitude.
			ac=getPointingAttitude(par.MJDREF, par.TSTART, par.TSTART+par.Exposure,
					par.RA*M_PI/180., par.Dec*M_PI/180., &status);
			CHECK_STATUS_BREAK(status);

		} else {
			// Load the attitude from the given file.
			ac=loadAttitude(par.Attitude, &status);
			CHECK_STATUS_BREAK(status);

			// Check if the required time interval for the simulation
			// is a subset of the period described by the attitude file.
			checkAttitudeTimeCoverage(ac, par.MJDREF, par.TSTART,
					par.TSTART+par.Exposure, &status);
			CHECK_STATUS_BREAK(status);
		}
		// END of setting up the attitude.

		// Get a GTI.
		gti=getGTIFromFileOrContinuous(par.GTIfile,
				par.TSTART, par.TSTART+par.Exposure,
				par.MJDREF, &status);
		CHECK_STATUS_BREAK(status);

		// Load the SIMPUT X-ray source catalogs.
		srccat[0]=loadSourceCatalog(par.Simput, inst->tel->arf, &status);
		CHECK_STATUS_BREAK(status);

		// Optional 2nd catalog.
		if (strlen(par.Simput2)>0) {
			strcpy(ucase_buffer, par.Simput2);
			strtoupper(ucase_buffer);
			if (0!=strcmp(ucase_buffer, "NONE")) {
				srccat[1]=loadSourceCatalog(par.Simput2, inst->tel->arf, &status);
				CHECK_STATUS_BREAK(status);
			}
		}

		// Optional 3rd catalog.
		if (strlen(par.Simput3)>0) {
			strcpy(ucase_buffer, par.Simput3);
			strtoupper(ucase_buffer);
			if (0!=strcmp(ucase_buffer, "NONE")) {
				srccat[2]=loadSourceCatalog(par.Simput3, inst->tel->arf, &status);
				CHECK_STATUS_BREAK(status);
			}
		}

		// Optional 4th catalog.
		if (strlen(par.Simput4)>0) {
			strcpy(ucase_buffer, par.Simput4);
			strtoupper(ucase_buffer);
			if (0!=strcmp(ucase_buffer, "NONE")) {
				srccat[3]=loadSourceCatalog(par.Simput4, inst->tel->arf, &status);
				CHECK_STATUS_BREAK(status);
			}
		}

		// Optional 5th catalog.
		if (strlen(par.Simput5)>0) {
			strcpy(ucase_buffer, par.Simput5);
			strtoupper(ucase_buffer);
			if (0!=strcmp(ucase_buffer, "NONE")) {
				srccat[4]=loadSourceCatalog(par.Simput5, inst->tel->arf, &status);
				CHECK_STATUS_BREAK(status);
			}
		}

		// Optional 6th catalog.
		if (strlen(par.Simput6)>0) {
			strcpy(ucase_buffer, par.Simput6);
			strtoupper(ucase_buffer);
			if (0!=strcmp(ucase_buffer, "NONE")) {
				srccat[5]=loadSourceCatalog(par.Simput6, inst->tel->arf, &status);
				CHECK_STATUS_BREAK(status);
			}
		}


		// --- End of Initialization ---


		// --- Open and set up files ---

		char telescop[MAXMSG]={""}, instrume[MAXMSG]={""};
		if (NULL!=inst->telescop) {
			strcpy(telescop, inst->telescop);
		}
		if (NULL!=inst->instrume) {
			strcpy(instrume, inst->instrume);
		}
		double tstop=gti->stop[gti->ngti-1];

		// Open the output photon list file.
		if (strlen(photonlist_filename)>0) {
			plf=openNewPhotonFile(photonlist_filename,
					telescop, instrume, "Normal",
					inst->tel->arf_filename, inst->det->rmf_filename,
					par.MJDREF, 0.0, par.TSTART, tstop,
					par.clobber, &status);
			CHECK_STATUS_BREAK(status);
		}

		// Open the output impact list file.
		if (strlen(impactlist_filename)>0) {
			ilf=openNewImpactFile(impactlist_filename,
					telescop, instrume, "Normal",
					inst->tel->arf_filename, inst->det->rmf_filename,
					par.MJDREF, 0.0, par.TSTART, tstop,
					par.clobber, &status);
			CHECK_STATUS_BREAK(status);
		}

		// Open the piximpact file
		pixilf=openNewPixImpFile(piximpactlist_filename,telescop, instrume, "Normal",
				inst->tel->arf_filename, inst->det->rmf_filename,par.XMLFile,impactlist_filename,
				par.MJDREF, 0.0, par.TSTART, tstop,
				par.clobber, &status);
		CHECK_STATUS_BREAK(status);

		// ---- TES initialization ----
		// Not in Christian's initialization part, but for the moment we need an already existing piximpact file
		if (!par.UseRMF){

			// Copy parameters in general parameters structure
			copyParams2GeneralStruct(par,&genpar,par.TSTART,tstop);

			// Build up init structure
			init = newInitStruct(&status);
			CHECK_STATUS_BREAK(status);
			tesinitialization(init,&genpar,&status);
			CHECK_STATUS_BREAK(status);
			// Only one piximpact file should be open at a time
			freePixImpFile(&(init->impfile), &status);
			// Deactivate all pixels because we want to treat one by one
			for (int i=0;i<init->det->npix;i++){
				init->activearray[i]=-1;
			}
			event_file = init->event_file;
			det=init->det;

			//Initialize reconstruction
			reconstruct_init = newReconstructInit(&status);
			CHECK_STATUS_BREAK(status);
			initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.PulseLength,
					par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
					par.DerivateExclusion,par.SaturationValue,&status);
			CHECK_STATUS_BREAK(status);
		} else{
			det = loadAdvDet(par.AdvXml,&status);
			keywords = buildSixtStdKeywords(telescop,instrume,"Normal",inst->tel->arf_filename, inst->det->rmf_filename,"NONE",par.MJDREF, 0.0, par.TSTART, tstop,&status);
			event_file = opennewTesEventFile(par.EvtFile,keywords,par.clobber,&status);
			loadRMFLibrary(det,&status);
		}

		if (status!=EXIT_SUCCESS) {
		  printf(" ERROR initializing the TES setup \n");
		  break;
		}
		// ---- End of TES initialization ----


		// Set FITS header keywords.
		// If this is a pointing attitude, store the direction in the output
		// photon list.
		if (1==ac->nentries) {
			// Determine the telescope pointing direction and roll angle.
		  Vector pointing=getTelescopeNz(ac, par.TSTART, &status);
			CHECK_STATUS_BREAK(status);

			// Direction.
			double ra, dec;
			calculate_ra_dec(pointing, &ra, &dec);

			// Roll angle.
			float rollangle=getRollAngle(ac, par.TSTART, &status);
			CHECK_STATUS_BREAK(status);

			// Store the RA and Dec information in the FITS header.
			ra *=180./M_PI;
			dec*=180./M_PI;
			rollangle*=180./M_PI;

			// Photon list file.
			if (NULL!=plf) {
				fits_update_key(plf->fptr, TDOUBLE, "RA_PNT", &ra,
						"RA of pointing direction [deg]", &status);
				fits_update_key(plf->fptr, TDOUBLE, "DEC_PNT", &dec,
						"Dec of pointing direction [deg]", &status);
				fits_update_key(plf->fptr, TFLOAT, "PA_PNT", &rollangle,
						"Roll angle [deg]", &status);
				CHECK_STATUS_BREAK(status);
			}

			// Impact list file.
			if (NULL!=ilf) {
				fits_update_key(ilf->fptr, TDOUBLE, "RA_PNT", &ra,
						"RA of pointing direction [deg]", &status);
				fits_update_key(ilf->fptr, TDOUBLE, "DEC_PNT", &dec,
						"Dec of pointing direction [deg]", &status);
				fits_update_key(ilf->fptr, TFLOAT, "PA_PNT", &rollangle,
						"Roll angle [deg]", &status);
				CHECK_STATUS_BREAK(status);
			}

			// Piximpact list file.
			if (NULL!=pixilf) {
				fits_update_key(pixilf->fptr, TDOUBLE, "RA_PNT", &ra,
						"RA of pointing direction [deg]", &status);
				fits_update_key(pixilf->fptr, TDOUBLE, "DEC_PNT", &dec,
						"Dec of pointing direction [deg]", &status);
				fits_update_key(pixilf->fptr, TFLOAT, "PA_PNT", &rollangle,
						"Roll angle [deg]", &status);
				CHECK_STATUS_BREAK(status);
			}

			// Record file.
			if (!par.UseRMF && NULL!=init->record_file) {
				fits_update_key(init->record_file->fptr, TDOUBLE, "RA_PNT", &ra,
						"RA of pointing direction [deg]", &status);
				fits_update_key(init->record_file->fptr, TDOUBLE, "DEC_PNT", &dec,
						"Dec of pointing direction [deg]", &status);
				fits_update_key(init->record_file->fptr, TFLOAT, "PA_PNT", &rollangle,
						"Roll angle [deg]", &status);
				CHECK_STATUS_BREAK(status);
			}

 			// Tes event file
			fits_update_key(event_file->fptr, TDOUBLE, "RA_PNT", &ra,
					"RA of pointing direction [deg]", &status);
			fits_update_key(event_file->fptr, TDOUBLE, "DEC_PNT", &dec,
					"Dec of pointing direction [deg]", &status);
			fits_update_key(event_file->fptr, TFLOAT, "PA_PNT", &rollangle,
					"Roll angle [deg]", &status);
			CHECK_STATUS_BREAK(status);

		} else {
			// An explicit attitude file is given.
			if (NULL!=plf) {
				fits_update_key(plf->fptr, TSTRING, "ATTITUDE", par.Attitude,
						"attitude file", &status);
			}
			if (NULL!=ilf) {
				fits_update_key(ilf->fptr, TSTRING, "ATTITUDE", par.Attitude,
						"attitude file", &status);
			}
			if (NULL!=pixilf) {
				fits_update_key(pixilf->fptr, TSTRING, "ATTITUDE", par.Attitude,
						"attitude file", &status);
			}
			if (!par.UseRMF && NULL!=init->record_file) {
				fits_update_key(init->record_file->fptr, TSTRING, "ATTITUDE", par.Attitude,
						"attitude file", &status);
			}
			fits_update_key(event_file->fptr, TSTRING, "ATTITUDE", par.Attitude,
					"attitude file", &status);
			CHECK_STATUS_BREAK(status);
		}

		// --- End of opening files ---

		// --- Initialize Crosstalk Structure ---
		det->crosstalk_id=par.doCrosstalk;
		if (det->crosstalk_id>0){
			headas_chat(3, "initializing crosstalk ...\n");
			init_crosstalk(det, &status);
			if (status!=EXIT_SUCCESS){
				SIXT_ERROR("failed when initializing crosstalk setup");
				return EXIT_FAILURE;
			}
			headas_chat(3, "\n");
		}

		// --- Simulation Process ---

		headas_chat(3, "start simulation ...\n");

		// Simulation progress status (running from 0 to 100).
		unsigned int progress=0;
		if (NULL==progressfile) {
			headas_chat(2, "\r%.0lf %%", 0.);
			fflush(NULL);
		} else {
			rewind(progressfile);
			fprintf(progressfile, "%.2lf", 0.);
			fflush(progressfile);
		}

		// Determine the total length of the time interval to
		// be simulated.
		double totalsimtime=sumGTI(gti);

		// Loop over all intervals in the GTI collection.
		int gtibin=0;
		double simtime=0.;
		long current_impact_row = 0;
		long current_impact_write_row = 0;
		long nimpacts=0;
		do {
			// Currently regarded interval.
			double t0=gti->start[gtibin];
			double t1=gti->stop[gtibin];

			// Set the start time for the instrument model.
			//setGenDetStartTime(inst->det, t0);

			// Loop over photon generation and processing
			// till the time of the photon exceeds the requested
			// time interval.
			do {

				// Photon generation.
				Photon ph;
				int isph=phgen(ac, srccat, MAX_N_SIMPUT, t0, t1, par.MJDREF, par.dt,
						inst->tel->fov_diameter, &ph, &status);
				CHECK_STATUS_BREAK(status);

				// If no photon has been generated, break the loop.
				if (0==isph) break;

				// Check if the photon still is within the requested
				// exposure time.
				assert(ph.time<=t1);

				// If requested, write the photon to the output file.
				if (NULL!=plf) {
					status=addPhoton2File(plf, &ph);
					CHECK_STATUS_BREAK(status);
				}

				// Photon imaging.
				Impact imp;
				int isimg=phimg(inst->tel, ac, &ph, &imp, &status);
				CHECK_STATUS_BREAK(status);

				// If the photon is not imaged but lost in the optical system,
				// continue with the next one.
				if (0==isimg) continue;

				// If requested, write the impact to the output file.
				if (NULL!=ilf) {
					addImpact2File(ilf, &imp, &status);
					CHECK_STATUS_BREAK(status);
				}

				// Piximpacts stage
				PixImpact *piximp=NULL;
				int newPixImpacts=AdvImpactList(det, &imp, &piximp);
				nimpacts+=newPixImpacts;
				if(newPixImpacts>0){
					for(int jj=0; jj<newPixImpacts; jj++){
						addImpact2PixImpFile(pixilf, &(piximp[jj]), &status);
					}
				}

				free(piximp);
				CHECK_STATUS_BREAK(status);

				// Program progress output.
				while((unsigned int)((ph.time-t0+simtime)*100./totalsimtime)>progress) {
					progress++;
					if (NULL==progressfile) {
						headas_chat(2, "\r%.0lf %%", progress*1.);
						fflush(NULL);
					} else {
						rewind(progressfile);
						fprintf(progressfile, "%.2lf", progress*1./100.);
						fflush(progressfile);
					}
				}

			} while(1);
			CHECK_STATUS_BREAK(status);
			// END of photon processing loop for the current interval.

			if (!par.UseRMF){
				headas_chat(3, "\nstart event reconstruction ...\n");
				// Generate the data streams for each pixel that has been hit
				int* list_pixels=NULL;
				getListPixelsHit(pixilf,&list_pixels,init->det->npix,&status);
				CHECK_STATUS_BREAK(status);

				// Close piximpact file
				current_impact_write_row = pixilf->row;
				freePixImpFile(&pixilf, &status);

				init->impfile=openPixImpFile(piximpactlist_filename, READONLY,&status);
				//CHECK_STATUS_BREAK(status);

				// Iterate over the pixels that were hit and run simulation
				for(int i=0;i<det->npix;i++){
					if(list_pixels[i]){
						// Activate corresponding pixel
						init->activearray[i]=0;
						genpar.nlo=i;
						genpar.nhi=i;

						// Reinitialize impact file
						init->impfile->row = current_impact_row;

						// Stream generation
						TESDataStream* stream=newTESDataStream(&status);
						CHECK_STATUS_BREAK(status);
						int ismonoc=0;
						float monoen=0.;

						getTESDataStream(stream,
								init->impfile,
								init->profiles,
								init->det,
								t0,
								t1,
								init->det->npix,
								1,
								init->activearray,
								init->Nevts,
								&ismonoc,
								&monoen,
								genpar.seed, // should modify this, we already have a random generator
								&status);
						CHECK_STATUS_BREAK(status);

						// Trigger and reconstruction
						triggerWithImpact(stream,&genpar,init,monoen,reconstruct_init,par.EventListSize,par.Identify,&status);
						CHECK_STATUS_BREAK(status);

						// Release this stream and deactivate the treated pixel
						destroyTESDataStream(stream);
						init->activearray[i]=-1;

					}
				}
				free(list_pixels);
				CHECK_STATUS_BREAK(status);

				// Close piximpact file in read mode (there should be only one file open)
				// after saving current row of impact file (this saves the row for next GTI)
				current_impact_row = init->impfile->row;
				freePixImpFile(&(init->impfile), &status);

				// Reopen piximpact file for writing for next gti
				pixilf = openPixImpFile(piximpactlist_filename,READWRITE,&status);
				pixilf->row = current_impact_write_row;
			} else{
				headas_chat(3, "\nstart event grading ...\n");
				current_impact_write_row = pixilf->row;
				pixilf->row = current_impact_row; // reboot pixilf to first row of the GTI

				impactsToEvents(det,pixilf,event_file,par.saveCrosstalk,&status);

				pixilf->row=current_impact_write_row;
				current_impact_row=current_impact_write_row;
			}

			// Proceed to the next GTI interval.
			simtime+=gti->stop[gtibin]-gti->start[gtibin];
			gtibin++;
			if (gtibin>=gti->ngti) break;

		} while (1);

		if (status==EXIT_FAILURE){
			SIXT_ERROR("xifupipeline had an internal error");
			break;
		}
		// End of loop over the individual GTI intervals.


		// Progress output.
		if (NULL==progressfile) {
			headas_chat(2, "\r%.0lf %%\n", 100.);
			fflush(NULL);
		} else {
			rewind(progressfile);
			fprintf(progressfile, "%.2lf", 1.);
			fflush(progressfile);
		}

		headas_chat(3, "start sky projection ...\n");
		event_file->row=1;
		phproj_advdet(inst,det,ac,event_file,par.TSTART,par.Exposure,par.ProjCenter,&status);
		CHECK_STATUS_BREAK(status);
		
		//Store number of impacts in event file
		fits_update_key(event_file->fptr, TLONG, "NIMP", &nimpacts,
				"Number of impacts", &status);
		
		// Store the GTI extension in the event file.
		saveGTIExt(event_file->fptr, "STDGTI", gti, &status);
		CHECK_STATUS_BREAK(status);

		// --- End of simulation process ---

		if (delete_rawdata){
			headas_chat(3,"removing unwanted RawData file %s \n",piximpactlist_filename);
			status = remove (piximpactlist_filename);
			CHECK_STATUS_BREAK(status);
		}

		CHECK_STATUS_BREAK(status);


	} while(0); // END of ERROR HANDLING Loop.


	// --- Clean up ---

	headas_chat(3, "\ncleaning up ...\n");

	// Release memory.
	freePhotonFile(&plf, &status);
	freeImpactFile(&ilf, &status);
	freePixImpFile(&pixilf, &status);
	for (ii=0; ii<MAX_N_SIMPUT; ii++) {
		freeSourceCatalog(&(srccat[ii]), &status);
	}
	freeGTI(&gti);
	freeAttitude(&ac);
	destroyGenInst(&inst, &status);
	freeTESInitStruct(&init,&status);
	freeReconstructInit(reconstruct_init);
	if(par.UseRMF){
		destroyAdvDet(&det);
		freeTesEventFile(event_file,&status);
	}


	if (NULL!=progressfile) {
		fclose(progressfile);
		progressfile=NULL;
	}

	// Clean up the random number generator.
	sixt_destroy_rng();

	if (EXIT_SUCCESS==status) {
		headas_chat(3, "finished successfully!\n\n");
		return(EXIT_SUCCESS);
	} else {
		headas_chat(3, " ... ERROR when cleaning up! !\n\n");
		return(EXIT_FAILURE);
	}
}


int xifupipeline_getpar(struct Parameters* const par)
{
	// String input buffer.
	char* sbuffer=NULL;

	// Error status.
	int status=EXIT_SUCCESS;

	// check if any obsolete keywords are given
	sixt_check_obsolete_keyword(&status);
	CHECK_STATUS_RET(status,EXIT_FAILURE);

	// Read all parameters via the ape_trad_ routines.

	status=ape_trad_query_string("Prefix", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the prefix for the output files");
		return(status);
	}
	strcpy(par->Prefix, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("PhotonList", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the photon list");
		return(status);
	}
	strcpy(par->PhotonList, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("ImpactList", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the impact list");
		return(status);
	}
	strcpy(par->ImpactList, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("PixImpactList", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the piximpact list");
		return(status);
	}
	strcpy(par->PixImpactList, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("EvtFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the event list");
		return(status);
	}
	strcpy(par->EvtFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("XMLFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the XML file");
		return(status);
	}
	strcpy(par->XMLFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("AdvXml", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the advanced XML file");
		return(status);
	}
	strcpy(par->AdvXml, sbuffer);
	free(sbuffer);

	status=ape_trad_query_bool("Background", &par->Background);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the background flag");
		return(status);
	}

	status=ape_trad_query_string("Attitude", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the attitude");
		return(status);
	}
	strcpy(par->Attitude, sbuffer);
	free(sbuffer);

	status=ape_trad_query_float("RA", &par->RA);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the right ascension of the telescope "
				"pointing");
		return(status);
	}

	status=ape_trad_query_float("Dec", &par->Dec);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the declination of the telescope "
				"pointing");
		return(status);
	}

	status=ape_trad_query_file_name("Simput", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the SIMPUT file");
		return(status);
	}
	strcpy(par->Simput, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Simput2", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the second SIMPUT file");
		return(status);
	}
	strcpy(par->Simput2, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Simput3", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the third SIMPUT file");
		return(status);
	}
	strcpy(par->Simput3, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Simput4", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the forth SIMPUT file");
		return(status);
	}
	strcpy(par->Simput4, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Simput5", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the fifth SIMPUT file");
		return(status);
	}
	strcpy(par->Simput5, sbuffer);
	free(sbuffer);

	status=ape_trad_query_file_name("Simput6", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the sixth SIMPUT file");
		return(status);
	}
	strcpy(par->Simput6, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("GTIfile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the GTI file");
		return(status);
	}
	strcpy(par->GTIfile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("MJDREF", &par->MJDREF);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading MJDREF");
		return(status);
	}

	status=ape_trad_query_double("TSTART", &par->TSTART);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading TSTART");
		return(status);
	}

	status=ape_trad_query_double("Exposure", &par->Exposure);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the exposure time");
		return(status);
	}

	status=ape_trad_query_double("dt", &par->dt);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading dt");
		return(status);
	}

	status=ape_trad_query_bool("SkipInvalids", &par->SkipInvalids);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the SkipInvalids parameter");
		return(status);
	}

	status=ape_trad_query_int("Seed", &par->Seed);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the seed for the random number generator");
		return(status);
	}

	status=ape_trad_query_string("ProgressFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the progress status file");
		return(status);
	}
	strcpy(par->ProgressFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_bool("UseRMF", &par->UseRMF);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the UseRMF parameter");
		return(status);
	}

	// query_simput_parameter_bool("doCrosstalk", &par->doCrosstalk, &status );
	char *buf;
	query_simput_parameter_string("doCrosstalk", &buf, &status );
	if (strncmp(buf,"yes",3)==0 ||strncmp(buf,"all",3)==0  ){
		par->doCrosstalk = CROSSTALK_ID_ALL;
	} else if (strncmp(buf,"elec",4)==0){
		par->doCrosstalk = CROSSTALK_ID_ELEC;
	} else if (strncmp(buf,"therm",5)==0){
		par->doCrosstalk = CROSSTALK_ID_THERM;
	} else if (strncmp(buf,"nlin",4)==0){
		par->doCrosstalk = CROSSTALK_ID_IMOD;
	} else {
		par->doCrosstalk=CROSSTALK_ID_NONE;
	}

	query_simput_parameter_bool("saveCrosstalk", &par->saveCrosstalk, &status );

	if (!par->UseRMF){
		status=ape_trad_query_string("TesTriggerFile", &sbuffer);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the name of the TES Trigger output file");
			return(status);
		}
		strcpy(par->TesTriggerFile, sbuffer);
		free(sbuffer);

		status=ape_trad_query_int("TriggerSize", &par->triggerSize);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the TriggerSize parameter");
			return(status);
		}

		status=ape_trad_query_int("PreBufferSize", &par->preBufferSize);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the PreBufferSize parameter");
			return(status);
		}

		status=ape_trad_query_string("EvtFile", &sbuffer);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the name of the event file");
			return(status);
		}
		strcpy(par->EvtFile, sbuffer);
		free(sbuffer);

		status=ape_trad_query_string("OptimalFilterFile", &sbuffer);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the name of the optimal filter file");
			return(status);
		}
		strcpy(par->OptimalFilterFile, sbuffer);
		free(sbuffer);

		status=ape_trad_query_string("PulseTemplateFile", &sbuffer);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the name of the pulse template file");
			return(status);
		}
		strcpy(par->PulseTemplateFile, sbuffer);
		free(sbuffer);

		status=ape_trad_query_int("PulseLength", &par->PulseLength);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the PulseLength parameter");
			return(status);
		}

		status=ape_trad_query_double("Threshold", &par->Threshold);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the Threshold parameter");
			return(status);
		}

		status=ape_trad_query_double("Calfac", &par->Calfac);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the Calfac parameter");
			return(status);
		}

		status=ape_trad_query_int("EventListSize", &par->EventListSize);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the EventListSize parameter");
			return(status);
		}

		status=ape_trad_query_int("NormalExclusion", &par->NormalExclusion);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the NormalExclusion parameter");
			return(status);
		}

		status=ape_trad_query_int("DerivateExclusion", &par->DerivateExclusion);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the DerivateExclusion parameter");
			return(status);
		}

		status=ape_trad_query_double("SaturationValue", &par->SaturationValue);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the SaturationValue parameter");
			return(status);
		}

		status=ape_trad_query_bool("Identify", &par->Identify);
		if (EXIT_SUCCESS!=status) {
			SIXT_ERROR("failed reading the Identify parameter");
			return(status);
		}
	}
	status=ape_trad_query_bool("clobber", &par->clobber);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the clobber parameter");
		return(status);
	}

	status=ape_trad_query_bool("ProjCenter", &par->ProjCenter);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the ProjCenter parameter");
		return(status);
	}

	return(status);
}

/** Copies the parameters contained in the local parameter structure into the
    general TES parameters structure */
void copyParams2GeneralStruct(const struct Parameters partmp, TESGeneralParameters* const par,double tstart,double tstop){
	strcpy(par->PixImpList,partmp.PixImpactList);
	strcpy(par->XMLFile,partmp.AdvXml);
	strcpy(par->tesTriggerFile,partmp.TesTriggerFile);
	strcpy(par->TesEventFile,partmp.EvtFile);

	par->Nactive=-1;
	par->nlo=-1;
	par->nhi=-1;
	par->triggerSize=partmp.triggerSize;
	par->preBufferSize=partmp.preBufferSize;
	par->Reconstruct=1;
	if (strlen(par->tesTriggerFile)>0){
		par->WriteRecordFile=1;
	} else{
		par->WriteRecordFile=0;
	}
	par->tstart=tstart;
	par->tstop=tstop;

	//par->writeStreamFile=partmp.writeStreamFile;
	par->clobber=partmp.clobber;
	par->history=partmp.history;
	par->check_times=0;

	par->seed=partmp.Seed;
}

/** Iterates over the piximpact file to see which pixels were actually hit */
void getListPixelsHit(PixImpFile* pixilf,int** list_pixels,int npix,int* const status){
	*list_pixels = malloc(npix*sizeof(int));
	if (NULL==*list_pixels){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for list_pixels failed");
		return;
	}
	for(int i=0;i<npix;i++){
		(*list_pixels)[i]=0;
	}

	long* pixids = malloc(pixilf->nrows*sizeof(*pixids));
	if (NULL==pixids){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for pixids failed");
		return;
	}
	int anynul=0;
	fits_read_col(pixilf->fptr, TLONG, pixilf->cpix_id,1,1,pixilf->nrows,
			NULL, pixids, &anynul, status);

	for(int i=0;i<pixilf->nrows;i++){
		(*list_pixels)[(int)(pixids[i]-1)]=1;
	}

}

