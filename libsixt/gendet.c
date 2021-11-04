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

#include "gendet.h"

////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////

GenDet* newGenDet(int* const status) {
	// Allocate memory.
	GenDet* det = (GenDet*) malloc(sizeof(GenDet));
	if (NULL == det) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("memory allocation for GenDet failed");
		return (det);
	}

	// Initialize all pointers with NULL.
	det->pixgrid = NULL;
	det->split = NULL;
	det->line = NULL;
	det->pha2pi_filename = NULL;
	det->pirmf_filename = NULL;
	det->specarf_filename = NULL;
	det->rmf_filename = NULL;
	det->rmf = NULL;
	det->elf = NULL;
	for (int ii = 0; ii < MAX_PHABKG; ii++) {
		det->phabkg[ii] = NULL;
	}
	det->clocklist = NULL;
	det->badpixmap = NULL;

	// Set initial values.
	det->ignore_bkg = 0;
	det->auxbackground = 0;
	det->split_bkg = 0;
	det->anyphoton = 0;
	det->frametime = 0.;
	det->deadtime = 0.;
	det->cte = 1.;
	det->rawymin = INT_MAX;
	det->rawymax = 0;
	det->threshold_readout_lo_keV = 0.;
	det->threshold_event_lo_keV = 0.;
	det->threshold_split_lo_keV = 0.;
	det->threshold_split_lo_fraction = 0.;
	det->threshold_pattern_up_keV = 0.;
	det->readout_trigger = 0;
	det->depfet.depfetflag = 0;
	det->depfet.istorageflag = 0;
	det->depfet.clear_const = NULL;
	det->depfet.clear_fcn = NULL;
        det->mxs_params = NULL;

	// Get empty GenPixGrid.
	det->pixgrid = newGenPixGrid(status);
	CHECK_STATUS_RET(*status, det);

	// Get empty ClockList.
	det->clocklist = newClockList(status);
	CHECK_STATUS_RET(*status, det);

	// Get empty split model.
	det->split = newGenSplit(status);
	CHECK_STATUS_RET(*status, det);

	return (det);
}

void destroyGenDet(GenDet** const det) {
	if (NULL != *det) {
		if (NULL != (*det)->line) {
			int ii;
			for (ii = 0; ii < (*det)->pixgrid->ywidth; ii++) {
				destroyGenDetLine(&(*det)->line[ii]);
			}
			free((*det)->line);
		}
		if (NULL != (*det)->pha2pi_filename) {
			free((*det)->pha2pi_filename);
		}
		if (NULL != (*det)->pirmf_filename) {
			free((*det)->pirmf_filename);
		}
		if (NULL != (*det)->specarf_filename) {
			free((*det)->specarf_filename);
		}

		if (NULL != (*det)->rmf_filename) {
			free((*det)->rmf_filename);
		}
		freeRMF((*det)->rmf);
		destroyClockList(&(*det)->clocklist);
		destroyGenPixGrid(&(*det)->pixgrid);
		destroyGenSplit(&(*det)->split);
		destroyBadPixMap(&(*det)->badpixmap);
		for (int ii = 0; ii < MAX_PHABKG; ii++) {
			destroyPHABkg(&(*det)->phabkg[ii]);
		}
		if ((*det)->depfet.clear_const != NULL) {
			free((*det)->depfet.clear_const);
			(*det)->depfet.clear_const = NULL;
		}
                if ((*det)->mxs_params != NULL) {
                         freeMXSParams( &((*det)->mxs_params) );
                }

		free(*det);
		*det = NULL;
	}
}

double depfet_get_linear_clear_signal(double time, double energy,
		double *constants) {

	return (energy * time / constants[0]);
}

double depfet_get_exponential_clear_signal(double time, double energy,
		double *constants) {

	return (energy * exp(-time / constants[0]));
}

int addGenDetPhotonImpact(GenDet* const det, const Impact* const impact,
		int* const status) {
	// Determine the detected energy.
	float energy;

	if (NULL != det->rmf) {
		// Determine the measured detector channel (PI channel) according
		// to the RMF.
		// The channel is obtained from the RMF using the corresponding
		// HEAdas routine which is based on drawing a random number.
		long channel;
		returnRMFChannel(det->rmf, impact->energy, &channel);

		// Check if the photon is really measured. If the PI channel
		// returned by the HEAdas RMF function is '-1', the photon is not
		// detected. This can happen, if the RMF includes ARF contributions,
		// e.g., the detector quantum efficiency and filter transmission.
		if (channel < det->rmf->FirstChannel) {
			return (0); // Break the function (photon is not detected).
		}

		// Determine the corresponding detected energy.
		// NOTE: In this simulation the collected charge is represented
		// by the nominal photon energy [keV], which corresponds to the
		// PI channel according to the EBOUNDS table.
		energy = getEBOUNDSEnergy(channel, det->rmf, status);
		CHECK_STATUS_RET(*status, 0);
		assert(energy >= 0.);

	} else {
		// The detector has no particular RMF. Therefore we directly
		// use the energy of the incident photon.
		energy = impact->energy;
	}

	// Create split events.
	int npixels = makeGenSplitEvents(det, &impact->position, energy,
			impact->ph_id, impact->src_id, impact->time, status);
	CHECK_STATUS_RET(*status, npixels);

	// Set the flag that there has been a photon interaction.
	det->anyphoton = 1;

	// Return the number of affected pixels.
	return (npixels);
}

void GenDetLineShift(GenDet* const det) {
	headas_chat(5, "lineshift\n");

	// Check if the detector contains more than 1 line.
	if (2 > det->pixgrid->ywidth)
		return;

	// Apply the Charge Transfer Efficiency.
	int ii;
	if (det->cte != 1.) {
		for (ii = 1; ii < det->pixgrid->ywidth; ii++) {
			if (0 != det->line[ii]->anycharge) {
				int jj;
				for (jj = 0; jj < det->line[ii]->xwidth; jj++) {
					if (det->line[ii]->charge[jj] > 0.) {
						det->line[ii]->charge[jj] *= det->cte;
					}
				}
			}
		}
	}

	// Add the charges in line 1 to line 0.
	addGenDetLine(det->line[0], det->line[1]);

	// Clear the charges in line 1, as they are now contained in line 0.
	clearGenDetLine(det->line[1]);

	// Shift the other lines in increasing order and put the newly cleared
	// original line number 1 at the end as the last line.
	GenDetLine* buffer = det->line[1];
	for (ii = 1; ii < det->pixgrid->ywidth - 1; ii++) {
		det->line[ii] = det->line[ii + 1];
	}
	det->line[det->pixgrid->ywidth - 1] = buffer;
}

void setGenDetEventFile(GenDet* const det, EventFile* const elf) {
	det->elf = elf;
}

void setGenDetIgnoreBkg(GenDet* const det, const int ignore) {
	if (0 == ignore) {
		det->ignore_bkg = 0;
	} else {
		det->ignore_bkg = 1;
	}
}

static inline void GenDetReadoutPixel(GenDet* const det, const int lineindex,
		const int readoutindex, const int xindex, const double time,
		int* const status) {
	headas_chat(5, "read out line %d as %d\n", lineindex, readoutindex);

	GenDetLine* line = det->line[lineindex];

	// Check if an output event file is defined.
	if (NULL == det->elf) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("no event file specified (needed for event detection)");
		return;
	}

	if (line->charge[xindex] != 0. || line->ccarry[xindex] != 0.) {
		Event* event = NULL;

		// Error handling loop.
		do {

			// Determine the properties of a new Event object.
			event = getEvent(status);
			CHECK_STATUS_BREAK(*status);

			// Readout the signal from the pixel array ...
			event->signal = line->charge[xindex];
			// ... overwrite the charge with the carry charge for the next cycle...
			line->charge[xindex] = line->ccarry[xindex];
			// ... and delete the carry charge.
			line->ccarry[xindex] = 0.;

			// Set the dead time of the pixel.
			line->deadtime[xindex] = time + det->deadtime;

			// Copy the information about the original photons.
			int jj;
			for (jj = 0; jj < NEVENTPHOTONS; jj++) {
				event->ph_id[jj] = line->ph_id[xindex][jj];
				event->src_id[jj] = line->src_id[xindex][jj];
				// Set the IDs to the carry IDs
				line->ph_id[xindex][jj] = line->carry_ph_id[xindex][jj];
				line->src_id[xindex][jj] = line->carry_src_id[xindex][jj];
				// Set the carry IDs to 0
				line->carry_ph_id[xindex][jj] = 0;
				line->carry_src_id[xindex][jj] = 0;
			}

			// Apply the lower readout threshold. Note that the upper
			// threshold is only applied after pattern recombination.
			if ((event->signal * event->signal)
					<= (det->threshold_readout_lo_keV
							* det->threshold_readout_lo_keV)) {
				break;
			}

			// Apply the detector response if available.
			if (NULL != det->rmf) {
				event->pha = getEBOUNDSChannel(event->signal, det->rmf);
			} else {
				event->pha = 0;
			}

			// Store remaining information.
			event->rawy = readoutindex;
			event->rawx = xindex;
			event->time = time;  // Time of detection.
			event->frame = det->clocklist->frame; // Frame of detection.
			event->npixels = 1;

			// Store the event in the output event file.
			addEvent2File(det->elf, event, status);
			CHECK_STATUS_BREAK(*status);

		} while (0); // END of error handling loop.

		// Release memory
		freeEvent(&event);
	}
}

void GenDetReadoutLine(GenDet* const det, const int lineindex,
		const int readoutindex, int* const status) {
	headas_chat(5, "read out line %d as %d\n", lineindex, readoutindex);

	GenDetLine* line = det->line[lineindex];
	line->last_readouttime = det->clocklist->time;

	if (0 != line->anycharge) {
		int ii;
		for (ii = 0; ii < line->xwidth; ii++) {
			GenDetReadoutPixel(det, lineindex, readoutindex, ii,
					det->clocklist->readout_time, status);
			CHECK_STATUS_BREAK(*status);
		}
		CHECK_STATUS_VOID(*status);
		// END of loop over all pixels in the line.

		// Reset the anycharge flag of this line.
		line->anycharge = line->anycarry;
		line->anycarry = 0;
	}
}

void GenDetClearLine(GenDet* const det, const int lineindex) {

	clearGenDetLine(det->line[lineindex]);
}

/** Apply the bad pixels of the bad pixel map. */
static void insertBadPix(GenDet* const det, const double timespan) {
	int ii;
	for (ii = 0; ii < det->badpixmap->xwidth; ii++) {
		if (1 == det->badpixmap->anybadpix[ii]) {
			int jj;
			for (jj = 0; jj < det->badpixmap->ywidth; jj++) {
				float diff = det->badpixmap->pixels[ii][jj] * timespan;
				if (det->badpixmap->pixels[ii][jj] > 0.) {
					// If the pixel is a hot one, add charge.
					addGenDetCharge2Pixel(det, ii, jj, diff, -1.0, -1, -1);
				} else if (det->badpixmap->pixels[ii][jj] < 0.) {
					// If the pixel is a cold one, remove charge.
					if (det->line[jj]->charge[ii] < (-1.) * diff) {
						det->line[jj]->charge[ii] = 0.;
					} else {
						det->line[jj]->charge[ii] += diff;
					}
				}
			}
			// END of loop over y-coordinate.
		}
	}
	// END of loop over x-coordinate.
}

static void insert_aux_bkg(GenDet* const det, double time, double dt) {
	// Get background events for the required time interval (has
	// to be given in [s]).
	backgroundOutput* list = bkgGetBackgroundList(dt);
	double cosrota = cos(det->pixgrid->rota);
	double sinrota = sin(det->pixgrid->rota);
	int ii;
	for (ii = 0; ii < list->numhits; ii++) {
		// Please note that the detector response matrix is
		// NOT applied to the particle-induced background events,
		// since the response is not available for the therefore
		// necessary high energies.
		// It is important that the high-energetic particle events
		// are only thrown away after the pattern recognition.
		// Otherwise many invalid particle patterns will be reduced
		// to apparently valid event patterns, such that the overall
		// background is too high.
		double xh = list->hit_xpos[ii] * 0.001 * cosrota
				- list->hit_ypos[ii] * 0.001 * sinrota;
		double yh = list->hit_xpos[ii] * 0.001 * sinrota
				+ list->hit_ypos[ii] * 0.001 * cosrota;
                if (det->split_bkg) {
                  struct Point2d pos = {.x=xh, .y=yh};
                  int status=EXIT_SUCCESS;
                  makeGenSplitEvents(det, &pos, list->hit_energy[ii],
                                     -1, -1, time, &status);
                  CHECK_STATUS_VOID(&status);
                } else {
 		  // Add the signal to the detector without
		  // regarding charge cloud splitting effects,
		  // since this is also not done by Tenzer et al. (2010)
		  // and Boller (2011).
                  int x, y;
                  double xr, yr;
                  getGenDetAffectedPixel(det->pixgrid, xh, yh, &x, &y, &xr, &yr);
                  // Check if the pixel indices are valid or if the
                  // specified position lies outside the pixel area.
                  if ((x < 0) || (y < 0))
                          continue;

                  // Add the signal to the pixel.
                  addGenDetCharge2Pixel(det, x, y, list->hit_energy[ii],
                                  time, -1, -1);
                }
	}
	bkgFree(list);
}

static float calc_offaxis_angle(int xi, int yi, int ii, GenDet* const det) {
	// Determine the off-axis angle.
	return atan(
			sqrt(
					pow((xi - det->pixgrid->xrpix + 1.0) * det->pixgrid->xdelt,
							2.0)
							+ pow(
									(yi - det->pixgrid->yrpix + 1.0)
											* det->pixgrid->ydelt, 2.0))
					/ *det->phabkg[ii]->focal_length);
}


// insert PHA background events for the required time interval (and all bkg models)
void insert_pha_bkg(GenDet* const det, double tstart, double dt, int* const status) {
	if (NULL == det->rmf) {
		SIXT_ERROR(
				"RMF needs to be defined for using the PHA background model");
		*status = EXIT_FAILURE;
		return;
	}
	// Loop over all PHA background models.
	int ii;
	for (ii = 0; ii < MAX_PHABKG; ii++) {
		// Check if the model is defined.
		if (det->phabkg[ii] == NULL)
			break;

		// Get background events for the required time interval.
		double bkg_time;
		long bkg_pha;
		while (getPHABkgEvent(det->phabkg[ii],
				det->pixgrid->xwidth * det->pixgrid->xdelt
						* det->pixgrid->ywidth * det->pixgrid->ydelt,
				tstart, tstart + dt,
				&bkg_time, &bkg_pha, status)) {
			CHECK_STATUS_BREAK(*status);
			// Determine the corresponding signal.
			float energy = getEBOUNDSEnergy(bkg_pha, det->rmf, status);
			CHECK_STATUS_BREAK(*status);
			// Determine the affected pixel.
			int xi = (int) (sixt_get_random_number(status)
					* det->pixgrid->xwidth);
			CHECK_STATUS_BREAK(*status);
			int yi = (int) (sixt_get_random_number(status)
					* det->pixgrid->ywidth);
			CHECK_STATUS_BREAK(*status);
			// If specified, apply vignetting.
			if (NULL != det->phabkg[ii]->vignetting) {
				// Check if vignetting function and focal length are given.
				if (NULL == *det->phabkg[ii]->vignetting) {
					*status = EXIT_FAILURE;
					SIXT_ERROR("vignetting function is need for "
							"vignetting-dependent background model");
					break;
				}
				if (0.0 == *det->phabkg[ii]->focal_length) {
					*status = EXIT_FAILURE;
					SIXT_ERROR("focal length is need for vignetting-dependent "
							"background model");
					break;
				}
				// Determine the off-axis angle.
				float theta = calc_offaxis_angle(xi, yi, ii, det);
				// Apply vignetting.
				double p = sixt_get_random_number(status);
				CHECK_STATUS_BREAK(*status);
				if (p
						> get_Vignetting_Factor(*det->phabkg[ii]->vignetting,
								energy, theta, 0.)) {
					// The background event is discarded due to vignetting.
					continue;
				}
			}
			// Add the signal to the pixel.
			addGenDetCharge2Pixel(det, xi, yi, energy, bkg_time, -1, -1);

			// for the Event Triggered Mode, trigger the read-out and advance the frame
			if (GENDET_EVENT_TRIGGERED == det->readout_trigger) {
				// Call the event trigger routine.
				GenDetReadoutPixel(det, yi, yi, xi, bkg_time, status);
				CHECK_STATUS_VOID(*status);

				// In event-triggered mode each event occupies its own frame.
				det->clocklist->frame++;
			}
		}
	}
}

static void jumpToNextFrame(GenDet* const det, const double time) {
	long nframes = (long) ((time - det->clocklist->readout_time)
			/ det->frametime);
	det->clocklist->time += nframes * det->frametime;
	det->clocklist->frame += nframes;
	det->clocklist->readout_time = det->clocklist->time;
	// Don't forget to add the time difference also to the single lines
	int ii;
	for (ii = 0; ii < det->pixgrid->ywidth; ii++) {
		det->line[ii]->last_readouttime += nframes * det->frametime;
	}
}

static void insert_background_events(GenDet* const det,
		double tstart, double dt, int* const status) {

	// Insert background events, if the appropriate PHA background
	// model is defined and should be used.
	if (NULL != det->phabkg[0]) {
		insert_pha_bkg(det, tstart, dt, status);
	}
	// Insert cosmic ray background events,
	// if the appropriate model is defined and should be used.
	if ((1 == det->auxbackground)) {
		insert_aux_bkg(det, tstart, dt );
	}
}

void operateGenDetClock(GenDet* const det, const double time, int* const status) {

	// Event-triggered mode. In this mode only background
	// events are inserted.
	if ( (GENDET_EVENT_TRIGGERED == det->readout_trigger )
			 && (0 == det->ignore_bkg)  ){

		static double last_time = 0.0; // Time of the last function call.

		// Insert background events (PHA and AUX)
		insert_background_events(det, last_time,time - last_time, status);
		CHECK_STATUS_VOID(*status);

		// Remember the time of the function call.
		last_time = time;

	} else if (GENDET_TIME_TRIGGERED == det->readout_trigger) {
		// Time-triggered mode.

		// Get the next element from the clock list.
		CLType type;
		void* element = NULL;
		do {
			CLReadoutLine* clreadoutline = NULL;
			CLClearLine* clclearline = NULL;
			CLWait* clwait = NULL;

			getClockListElement(det->clocklist, time, &type, &element, status);
			CHECK_STATUS_VOID(*status);

			switch (type) {
			case CL_NONE:
				// No operation has to be performed. The clock list is
				// currently in a wait status.
				break;
			case CL_NEWFRAME:
				// The clock list has internally increased the frame counter and readout
				// time.

				// If there has been no photon interaction during the last frame
				// and if no background model is activated, jump over the next empty frames
				// until there is a new photon impact.
				if ((((0 == det->auxbackground) && (NULL == det->phabkg[0])
						&& (NULL == det->phabkg[1])) || (1 == det->ignore_bkg))
						&& (0 == det->anyphoton)) {
					jumpToNextFrame(det, time);
				}

				// Reset the flag.
				det->anyphoton = 0;

				break;
			case CL_WAIT:
				// A waiting period is finished.
				clwait = (CLWait*) element;

				// Insert background events, if the appropriate PHA background
				// model is defined and should be used.
				if (0 == det->ignore_bkg){
					insert_background_events(det, det->clocklist->time, clwait->time, status);
					CHECK_STATUS_VOID(*status);
				}


				// Apply the hot pixels of the bad pixel map (if available) using
				// the pixel values weighted with the waiting time.
				if (NULL != det->badpixmap) {
					insertBadPix(det, clwait->time);
				}
				break;
			case CL_LINESHIFT:
				GenDetLineShift(det);
				break;
			case CL_READOUTLINE:
				clreadoutline = (CLReadoutLine*) element;
				GenDetReadoutLine(det, clreadoutline->lineindex,
						clreadoutline->readoutindex, status);
				CHECK_STATUS_VOID(*status)
				;
				break;
			case CL_CLEARLINE:
				clclearline = (CLClearLine*) element;
				GenDetClearLine(det, clclearline->lineindex);
				break;
			}
			CHECK_STATUS_VOID(*status);
		} while (type != CL_NONE);
	}
}

GenSplit* newGenSplit(int* const status) {
	// Allocate memory.
	GenSplit* split = (GenSplit*) malloc(sizeof(GenSplit));
	CHECK_NULL(split, *status, "memory allocation for GenSplit failed");

	// Initialize all pointers with NULL.

	// Set default values.
	split->type = GS_NONE;
	split->par1 = 0.;
	split->par2 = 0.;

	return (split);
}

void destroyGenSplit(GenSplit** const split) {
	if (NULL != *split) {
		free(*split);
		*split = NULL;
	}
}

static inline int getMinimumDistance(const double array[]) {
	int count, index = 0;
	double minimum = array[0];

	for (count = 1; count < 4; count++) {
		if ((minimum < 0.)
				|| ((array[count] <= minimum) && (array[count] >= 0.))) {
			minimum = array[count];
			index = count;
		}
	}

	return (index);
}

int makeGenSplitEvents(GenDet* const det, const struct Point2d* const position,
		const float signal, const long ph_id, const long src_id,
		const double time, int* const status) {
	// Number of affected pixels.
	int npixels = 0;
	// x- and y-indices of affected pixels.
	int x[4], y[4];
	// Signal fractions in the individual pixels.
	float fraction[4];

	// The following array entries are used to transform between
	// different array indices for accessing neighboring pixels.
	const int xe[4] = { 1, 0, -1, 0 };
	const int ye[4] = { 0, 1, 0, -1 };

	// Which kind of split model has been selected?
	if (GS_NONE == det->split->type) {
		// No split events => all events are singles.
		npixels = 1;

		// Determine the affected detector line and column.
		double xr, yr;
		getGenDetAffectedPixel(det->pixgrid, position->x, position->y, &(x[0]),
				&(y[0]), &xr, &yr);

		// Check if the returned values are valid line and column indices.
		if ((x[0] < 0) || (y[0] < 0)) {
			return (0);
		}

		// The single pixel receives the total photon energy.
		fraction[0] = 1.;

	} else if (GS_GAUSS == det->split->type) {
		// Gaussian split model.

		// Signal cloud sigma as a function of the photon energy.
		const float ccsigma = det->split->par1
				+ det->split->par2 * sqrt(signal);

		// Signal cloud size (3 sigma).
		const float ccsize = ccsigma * 3.;

		// Calculate pixel indices (integer) of the central affected pixel:
		double xr, yr;
		getGenDetAffectedPixel(det->pixgrid, position->x, position->y, &(x[0]),
				&(y[0]), &xr, &yr);

		// Check if the impact position lies inside the detector pixel array.
		if ((x[0] < 0) || (y[0] < 0)) {
			return (0);
		}

		// Calculate the distances from the impact center position to the
		// borders of the surrounding pixel (in [m]).
		double distances[4] = {
		// Distance to right pixel edge.
				(1.0 - xr) * det->pixgrid->xdelt,
				// Distance to upper edge.
				(1.0 - yr) * det->pixgrid->ydelt,
				// Distance to left pixel edge.
				xr * det->pixgrid->xdelt,
				// distance to lower edge
				yr * det->pixgrid->ydelt };

		int mindist = getMinimumDistance(distances);
		if (distances[mindist] < ccsize) {
			// Not a single event!
			x[1] = x[0] + xe[mindist];
			y[1] = y[0] + ye[mindist];

			double mindistgauss = gaussint(distances[mindist] / ccsigma);

			// Search for the next to minimum distance to an edge.
			double minimum = distances[mindist];
			distances[mindist] = -1.;
			int secmindist = getMinimumDistance(distances);
			distances[mindist] = minimum;

			if (distances[secmindist] < ccsize) {
				// Quadruple!
				npixels = 4;

				x[2] = x[0] + xe[secmindist];
				y[2] = y[0] + ye[secmindist];
				x[3] = x[1] + xe[secmindist];
				y[3] = y[1] + ye[secmindist];

				// Calculate the different signal fractions in the 4 affected pixels.
				double secmindistgauss = gaussint(
						distances[secmindist] / ccsigma);
				fraction[0] = (1. - mindistgauss) * (1. - secmindistgauss);
				fraction[1] = mindistgauss * (1. - secmindistgauss);
				fraction[2] = (1. - mindistgauss) * secmindistgauss;
				fraction[3] = mindistgauss * secmindistgauss;

			} else {
				// Double!
				npixels = 2;

				fraction[0] = 1. - mindistgauss;
				fraction[1] = mindistgauss;

			} // END of Double or Quadruple.

		} else {
			// Single event!
			npixels = 1;
			fraction[0] = 1.;
		}
		// END of check for Single event.

		// END of Gaussian split model.

	} else if (GS_EXPONENTIAL == det->split->type) {
		// Exponential split model.
		// None-Gaussian, exponential signal cloud model
		// (concept proposed by Konrad Dennerl).
		npixels = 4;

		// Calculate pixel indices (integer) of central affected pixel.
		double xr, yr;
		getGenDetAffectedPixel(det->pixgrid, position->x, position->y, &(x[0]),
				&(y[0]), &xr, &yr);

		// Check if the impact position lies inside the detector pixel array.
		if ((x[0] < 0) || (y[0] < 0)) {
			return (0);
		}

		// Calculate the distances from the impact center position to the
		// borders of the surrounding pixel (in units [fraction of a pixel edge]).
		double distances[4] = {
		// Distance to right pixel edge.
				(1.0 - xr),
				// Distance to upper edge.
				(1.0 - yr),
				// Distance to left pixel edge.
				xr,
				// distance to lower edge
				yr };

		// Search for the minimum distance to the edges.
		int mindist = getMinimumDistance(distances);
		x[1] = x[0] + xe[mindist];
		y[1] = y[0] + ye[mindist];

		// Search for the next to minimum distance to the edges.
		double minimum = distances[mindist];
		distances[mindist] = -1.;
		int secmindist = getMinimumDistance(distances);
		distances[mindist] = minimum;
		// Pixel coordinates of the 3rd and 4th split partner.
		x[2] = x[0] + xe[secmindist];
		y[2] = y[0] + ye[secmindist];
		x[3] = x[1] + xe[secmindist];
		y[3] = y[1] + ye[secmindist];

		// Now we know the affected pixels and can determine the
		// signal fractions according to the model exp(-(r/0.355)^2).
		// Remember that the array distances[] contains the distances
		// to the pixel borders, whereas here we need the distances from
		// the pixel center for the parameter r.
		// The value 0.355 is given by the parameter ecc->parameter.
		fraction[0] = exp(
				-(pow(0.5 - distances[mindist], 2.)
						+ pow(0.5 - distances[secmindist], 2.))
						/ pow(det->split->par1, 2.));
		fraction[1] = exp(
				-(pow(0.5 + distances[mindist], 2.)
						+ pow(0.5 - distances[secmindist], 2.))
						/ pow(det->split->par1, 2.));
		fraction[2] = exp(
				-(pow(0.5 - distances[mindist], 2.)
						+ pow(0.5 + distances[secmindist], 2.))
						/ pow(det->split->par1, 2.));
		fraction[3] = exp(
				-(pow(0.5 + distances[mindist], 2.)
						+ pow(0.5 + distances[secmindist], 2.))
						/ pow(det->split->par1, 2.));
		// Normalization to 1.
		double sum = fraction[0] + fraction[1] + fraction[2] + fraction[3];
		fraction[0] /= sum;
		fraction[1] /= sum;
		fraction[2] /= sum;
		fraction[3] /= sum;

		// END of exponential split model.

	} else {
		SIXT_ERROR("split model not supported");
		*status = EXIT_FAILURE;
		return (0);
	}

	// Add signal to all valid pixels of the split event.
	int ii, nvalidpixels = 0;
	for (ii = 0; ii < npixels; ii++) {


		if ((x[ii] >= 0) && (x[ii] < det->pixgrid->xwidth) && (y[ii] >= 0)
				&& (y[ii] < det->pixgrid->ywidth)) {

			addGenDetCharge2Pixel(det, x[ii], y[ii], signal * fraction[ii],
					time, ph_id, src_id);
			nvalidpixels++;

			// Call the event trigger routine.
			if (GENDET_EVENT_TRIGGERED == det->readout_trigger) {
				GenDetReadoutPixel(det, y[ii], y[ii], x[ii], time, status);
				CHECK_STATUS_BREAK(*status);

				// In event-triggered mode each event occupies its own frame.
				det->clocklist->frame++;
			}
		}
	}
	CHECK_STATUS_RET(*status, nvalidpixels);

	// Return the number of affected pixels.
	return (nvalidpixels);
}

void addGenDetCharge2Pixel(GenDet* const det, const int column, const int row,
		const float signal, const double time, const long ph_id,
		const long src_id) {
	GenDetLine* line = det->line[row];

	// Check if the pixel is sensitive right now.
	if ((time < line->deadtime[column]) && (time >= 0.0))
		return;

	// Check if pixel is in a readout area (todo: why do we need this? )
	if ((row >= 0 && row < det->pixgrid->ywidth)) {


		float oldcharge = line->charge[column];

		// Add the signal if the event is inside the readout window.
		int sign = 1;
		if (det->depfet.depfetflag != 1) {
			line->charge[column] += signal;
		} else {
			sign = addDepfetSignal(det, column, row, signal, time, ph_id,
					src_id);
		}
		line->anycharge = 1;

		// Set PH_ID and SRC_ID.
		if (oldcharge < 0.001) {
			// If the charge collect in the pixel up to now is below 1eV,
			// overwrite the old PH_ID and SRC_ID by the new value.
			line->ph_id[column][0] = sign * ph_id;
			line->src_id[column][0] = src_id;

		} else if ((signal > 0.001)) {
			// Only store the PH_ID and SRC_ID of the new contribution
			// if its signal is above 1eV.
			long ii;
			for (ii = 0; ii < NEVENTPHOTONS; ii++) {
				if (0 == line->ph_id[column][ii]) {
					line->ph_id[column][ii] = sign * ph_id;
					line->src_id[column][ii] = src_id;
					break;
				}
			}
		}
	}
}

void setGenDetStartTime(GenDet* const det, const double t0) {
	det->clocklist->time = t0;
	det->clocklist->readout_time = t0;
}

int addDepfetSignal(GenDet* const det, const int colnum, const int row,
		const float signal, const double time, const long ph_id,
		const long src_id) {

	GenDetLine* line = det->line[row];

	int sign = 1;

	double rtime = time - line->last_readouttime;
	//printf("New event: rtime=%lf\n", rtime);
//printf("last readout time of line: %lf\n", line->last_readouttime);
	if (det->depfet.istorageflag == 0) {
		// Normal DEPFET.

		// Determine time since the start of the readout cycle

		// Check if the time makes sense
		if (rtime < 0 || rtime > det->frametime) {
			//puts("time oob: rtime<0 || rtime>det->frametime");
			//printf("line     =%d\nframe    =%ld\ntime     =%lf \nl_readout=%lf \nrtime    =%lf \nf_time   =%lf\n",
			//	row, det->clocklist->frame, time, line->last_readouttime, rtime, det->frametime);
			//SIXT_ERROR("time since the start of the frame out of bounds.");
		}

		// t_readout is the time length of the active readout
		double t_readout = det->depfet.t_settling
				+ 2. * det->depfet.t_integration + det->depfet.t_clear;

		// t_wait is the time interval from the beginning of the
		// cycle to the begin of the read-out.
		double t_wait = det->frametime - t_readout;

		if (rtime <= t_wait) {
			// The photon arrives in the normal exposure interval.
			line->charge[colnum] += signal;

		} else {
			double ti1 = rtime - t_wait;
			if (ti1 <= det->depfet.t_integration) {
				// The photon arrives during the first integration.
				// It is detected incompletely but gets cleared afterwards.
				line->charge[colnum] += signal
						* (det->depfet.t_integration - ti1)
						/ det->depfet.t_integration;

			} else {

				// In the following cases, there is always a carry to the next frame
				line->anycarry = 1;
				sign = -1;
//if(sign<0)printf("sign=%d: row=%d; rtime=%le; frametime=%le; t_wait=%le, \n", sign, row, rtime, det->frametime, t_wait);

				// Set PH_ID and SRC_ID in carry-arrays.
				if (line->ccarry[colnum] < 0.001) {
					// If the charge collect in the pixel up to now is below 1eV,
					// overwrite the old PH_ID and SRC_ID by the new value.
					line->carry_ph_id[colnum][0] = ph_id;
					line->carry_src_id[colnum][0] = src_id;

				} else if (signal > 0.001) {
					// Only store the PH_ID and SRC_ID of the new contribution
					// if its signal is above 1eV.
					long ii;
					for (ii = 0; ii < NEVENTPHOTONS; ii++) {
						if (0 == line->carry_ph_id[colnum][ii]) {
							line->carry_ph_id[colnum][ii] = ph_id;
							line->carry_src_id[colnum][ii] = src_id;
							break;
						}
					}
				}

				double tc = ti1 - det->depfet.t_integration;
				if (tc <= det->depfet.t_clear) {
					// The photon arrives during the clear time.
					// It is cleared incompletely and the remaining signal
					// gets detected negatively in this frame, positively
					// in the next
					float rem = det->depfet.clear_fcn(tc, signal,
							det->depfet.clear_const);
					line->charge[colnum] += rem * (-1.);
					line->ccarry[colnum] += rem;

				} else {
					double ts = tc - det->depfet.t_clear;
					if (ts <= det->depfet.t_settling) {
						// The photon arrives during the second settling time.
						// It is detected negatively in this frame, positively
						// in the next frame.
						line->charge[colnum] += signal * (-1.);
						line->ccarry[colnum] += signal;

					} else {
						// The photon arrives during the second integration.
						// It is partially detected negatively in this frame
						// and fully detected positively in the next one.
						float intsig = signal
								* (det->depfet.t_integration - ts
										+ det->depfet.t_settling)
								/ det->depfet.t_integration;
						line->charge[colnum] += intsig * (-1.);
						line->ccarry[colnum] += signal;
					}
				}
			}
		}
	} else {
		// IS-DEPFET

		// t_frame is the time since the frame started
		double t_frame = time - det->clocklist->readout_time;
		// Determine time since the start of the readout cycle

		//Determine time interval of clear
		double cstart = det->frametime
				- (det->depfet.t_integration + det->depfet.t_settling
						+ det->depfet.t_clear);
		double cstop = det->frametime
				- (det->depfet.t_integration + det->depfet.t_settling);

		// Check if the time makes sense
		if (t_frame < 0 || t_frame > det->frametime) {
			//puts("time oob: t_frame<0 || t_frame>det->frametime");
			//SIXT_ERROR("time since the start of the frame out of bounds.");
		}
		if (rtime < 0 || rtime > det->frametime) {
			//puts("time oob: rtime<0 || rtime>det->frametime");
			//SIXT_ERROR("time since the start of the frame out of bounds.");
		}

		double transfertime = det->frametime - det->depfet.t_transfer;

		if (t_frame > transfertime) {
			// The photon arrives during the transfer time
			line->anycarry = 1;
			sign = -1;

			// Set PH_ID and SRC_ID in carry-arrays.
			if (line->ccarry[colnum] < 0.001) {
				// If the charge collect in the pixel up to now is below 1eV,
				// overwrite the old PH_ID and SRC_ID by the new value.
				line->carry_ph_id[colnum][0] = ph_id;
				line->carry_src_id[colnum][0] = src_id;

			} else if (signal > 0.001) {
				// Only store the PH_ID and SRC_ID of the new contribution
				// if its signal is above 1eV.
				long ii;
				for (ii = 0; ii < NEVENTPHOTONS; ii++) {
					if (0 == line->carry_ph_id[colnum][ii]) {
						line->carry_ph_id[colnum][ii] = ph_id;
						line->carry_src_id[colnum][ii] = src_id;
						break;
					}
				}
			}

			double intransfertime = t_frame - transfertime;

			float s_now = signal * (intransfertime) / det->depfet.t_transfer;
			float s_carry = signal - s_now;
			line->charge[colnum] += s_now;
			line->ccarry[colnum] += s_carry;
		} else {
			// The photon arrives during the normal exposure interval

			// Check if the time is in the clear interval
			if (rtime > cstart && rtime < cstop) {

				line->charge[colnum] += det->depfet.clear_fcn((cstop - rtime),
						signal, det->depfet.clear_const);
			} else {
				if (t_frame <= rtime) {

					line->anycarry = 1;
					//	  sign=0; // (set it to zero, as no charge is transfered in on
					// the current cycle)

					// Set PH_ID and SRC_ID in carry-arrays.
					if (line->ccarry[colnum] < 0.001) {
						// If the charge collect in the pixel up to now is below 1eV,
						// overwrite the old PH_ID and SRC_ID by the new value.
						line->carry_ph_id[colnum][0] = ph_id;
						line->carry_src_id[colnum][0] = src_id;

					} else if (signal > 0.001) {
						// Only store the PH_ID and SRC_ID of the new contribution
						// if its signal is above 1eV.
						long ii;
						for (ii = 0; ii < NEVENTPHOTONS; ii++) {
							if (0 == line->carry_ph_id[colnum][ii]) {
								line->carry_ph_id[colnum][ii] = ph_id;
								line->carry_src_id[colnum][ii] = src_id;
								break;
							}
						}
					}
					line->ccarry[colnum] += signal;
				} else {
					line->charge[colnum] += signal;
				}
			}
		}
	}

	return sign;
}
