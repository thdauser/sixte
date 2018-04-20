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

#include "phpat.h"

struct PatternStatistics {
	/** Number of valid patterns. */
	long nvalids;
	/** Number of valid patterns flagged as pile-up. */
	long npvalids;

	/** Number of invalid patterns. */
	long ninvalids;
	/** Number of invalid patterns flagged as pile-up. */
	long npinvalids;

	/** Number of patterns with a particular grade. */
	long ngrade[13];
	/** NUmber of patterns with a particular grade flagged as
	 pile-up. */
	long npgrade[13];
};

static int isNeighbor(const Event* const e1, const Event* const e2) {
	if (((e1->rawx == e2->rawx + 1) && (e1->rawy == e2->rawy))
			|| ((e1->rawx == e2->rawx - 1) && (e1->rawy == e2->rawy))
			|| ((e1->rawx == e2->rawx) && (e1->rawy == e2->rawy + 1))
			|| ((e1->rawx == e2->rawx) && (e1->rawy == e2->rawy - 1))) {
		// The events are neighbors.
		return (1);
	} else {
		return (0);
	}
}

void phpat(GenDet* const det, const EventFile* const src, EventFile* const dest,
		const char* picorr_file, const unsigned int seed, const char skip_invalids, int* const status) {

	// Pattern / grade statistics.
	struct PatternStatistics statistics = { .nvalids = 0, .npvalids = 0,
			.ninvalids = 0, .npinvalids = 0, };
	long ii;
	for (ii = 0; ii < 13; ii++) {
		statistics.ngrade[ii] = 0;
		statistics.npgrade[ii] = 0;
	}

	// List of all events belonging to the current frame.
	const long maxnframelist = 10000;
	Event** framelist = NULL;
	long nframelist = 0;

	// List of all neighboring events in the current frame.
	const long maxnneighborlist = 2000;
	Event** neighborlist = NULL;
	long nneighborlist = 0;

	// Flag, if we analyse the eROSITA-CCD.
	int iseROSITA;

	// Flag, whether the warning that the split threshold lies above the
	// event threshold has already been printed.
	static int threshold_warning_printed = 0;

	// TODO: LOAD Pha2PI Correction File if existent
	Pha2Pi* p2p = initPha2Pi(picorr_file, seed, status);

	// Error handling loop.
	do {

		// Check if the input file contains single-pixel events.
		char evtype[MAXMSG], comment[MAXMSG];
		fits_read_key(src->fptr, TSTRING, "EVTYPE", evtype, comment, status);
		CHECK_STATUS_BREAK_WITH_FITSERROR(*status);
		strtoupper(evtype);
		if (0 != strcmp(evtype, "PIXEL")) {
			*status = EXIT_FAILURE;
			char msg[MAXMSG];
			sprintf(msg, "event type of input file is '%s' (must be 'PIXEL')",
					evtype);
			SIXT_ERROR(msg);
			break;
		}

		// Set the event type in the output file to 'PATTERN'.
		fits_update_key(dest->fptr, TSTRING, "EVTYPE", "PATTERN", "event type",
				status);
		CHECK_STATUS_BREAK_WITH_FITSERROR(*status);

		// Allocate memory.
		framelist = (Event**) malloc(maxnframelist * sizeof(Event*));
		CHECK_NULL_BREAK(framelist, *status,
				"memory allocation for frame list failed");
		neighborlist = (Event**) malloc(maxnneighborlist * sizeof(Event*));
		CHECK_NULL_BREAK(neighborlist, *status,
				"memory allocation for neighbor list failed");

		// Determine the name of the instrument.
		// Particular instruments require a special pattern
		// recombination scheme (e.g. eROSITA).
		char telescope[MAXMSG];
		fits_read_key(src->fptr, TSTRING, "TELESCOP", telescope, comment, status);
		CHECK_STATUS_BREAK_WITH_FITSERROR(*status);
		strtoupper(telescope);
		if (!strcmp(telescope, "EROSITA")) {
			iseROSITA = 1;
		} else {
			iseROSITA = 0;
		}

		// Loop over all events in the input list.
		for (ii = 0; ii <= src->nrows; ii++) {

			// Read the next event from the file, if the end has not been
			// reached so far.
			Event* event = NULL;
			if (ii < src->nrows) {
				event = getEvent(status);
				CHECK_STATUS_BREAK(*status);
				getEventFromFile(src, ii + 1, event, status);
				CHECK_STATUS_BREAK(*status);
			}


			// Check if the new event belongs to a different frame than
			// the previous ones.
			int newframe = 0;
			if ((nframelist > 0) && (NULL != event)) {
				if (event->frame != framelist[0]->frame) {
					newframe = 1;
				}
			}

			// If the new event belongs to a different frame than the
			// previous ones or the end of the input list has been
			// reached, perform a pattern analysis.
			if ((newframe > 0) || ((ii == src->nrows) && (nframelist > 0))) {

				// Loop over all events in the current frame.
				long jj;
				for (jj = 0; jj < nframelist; jj++) {
					if (NULL != framelist[jj]) {

						// Check if the event is below the threshold.
						if ((framelist[jj]->signal * framelist[jj]->signal)
								< (det->threshold_event_lo_keV * det->threshold_event_lo_keV))
							continue;

						// Start a new neighbor list.
						neighborlist[0] = framelist[jj];
						nneighborlist = 1;
						framelist[jj] = NULL;

						// Find the signal maximum in the neighboring pixels.
						Event* maxsignalev = neighborlist[0];
						int updated = 0;
						do {
							updated = 0;
							long ll;
							for (ll = 0; ll < nframelist; ll++) {
								if (NULL != framelist[ll]) {
									if (isNeighbor(maxsignalev, framelist[ll])) {
										if (framelist[ll]->signal > maxsignalev->signal) {
											maxsignalev = framelist[ll];
											updated = 1;
										}
									}
								}
							}
						} while (updated);

						// Determine the split threshold [keV].
						float split_threshold;
						// For eROSITA we need a special treatment (according to
						// a prescription of K. Dennerl).
						if (1 == iseROSITA) {
							if (det->threshold_split_lo_fraction > 0.) {
								float vertical = 0., horizontal = 0.;
								long ll;
								for (ll = 0; ll < nframelist; ll++) {
									if (NULL != framelist[ll]) {
										if (isNeighbor(maxsignalev, framelist[ll])) {
											if (framelist[ll]->rawx == maxsignalev->rawx) {
												if (framelist[ll]->signal > horizontal) {
													horizontal = framelist[ll]->signal;
												}
											} else {
												if (framelist[ll]->signal > vertical) {
													vertical = framelist[ll]->signal;
												}
											}
										}
									}
								}
								split_threshold = det->threshold_split_lo_fraction
										* (maxsignalev->signal + horizontal + vertical);
							} else {
								SIXT_ERROR("eROSITA requires a fractional split threshold");
								*status = EXIT_FAILURE;
								break;
							}
						} else { // Split threshold for generic instruments.
							if (det->threshold_split_lo_fraction > 0.) {
								split_threshold = det->threshold_split_lo_fraction
										* maxsignalev->signal;
							} else {
								split_threshold = det->threshold_split_lo_keV;
							}
						}
						// END of determine the split threshold.

						// Check if the split threshold is above the event threshold.
						if ((split_threshold > det->threshold_event_lo_keV)
								&& (0 == threshold_warning_printed)) {
							char msg[MAXMSG];
							sprintf(msg,
									"split threshold (%.1feV) is above event threshold (%.1feV) "
											"(message is printed only once)",
									split_threshold * 1000.0,
									det->threshold_event_lo_keV * 1000.0);
							SIXT_WARNING(msg);
							threshold_warning_printed = 1;
						}

						// Find all neighboring events above the split threshold.
						long kk;
						for (kk = 0; kk < nneighborlist; kk++) {
							long ll;
							for (ll = 0; ll < nframelist; ll++) {
								if (NULL != framelist[ll]) {
									if (isNeighbor(neighborlist[kk], framelist[ll])) {

										// Check if its signal is below the split threshold.
										if (framelist[ll]->signal < split_threshold) {
											continue;
										}

										// Add the event to the neighbor list.
										if (nneighborlist >= maxnneighborlist) {
											SIXT_ERROR("too many events in the same pattern");
											*status = EXIT_FAILURE;
											break;
										}
										neighborlist[nneighborlist] = framelist[ll];
										nneighborlist++;
										framelist[ll] = NULL;
									}
								}
							}
							CHECK_STATUS_BREAK(*status);
						}
						CHECK_STATUS_BREAK(*status);
						// END of finding all neighbors.

						// Search the pixel with the maximum signal.
						long maxidx = 0;
						for (kk = 1; kk < nneighborlist; kk++) {
							if (neighborlist[kk]->signal > neighborlist[maxidx]->signal) {
								maxidx = kk;
							}
						}
						// END of searching the pixel with the maximum signal.

						// Get a new event.
						Event* event = getEvent(status);
						CHECK_STATUS_BREAK(*status);

						// Set basic properties.
						event->rawx = neighborlist[maxidx]->rawx;
						event->rawy = neighborlist[maxidx]->rawy;
						event->time = neighborlist[maxidx]->time;
						event->frame = neighborlist[maxidx]->frame;
						event->ra = 0.;
						event->dec = 0.;
						event->npixels = nneighborlist;

						// Set the advanced properties.
						// Total signal.
						event->signal = 0.;
						// Flag whether event touches the border of the detector.
						int border = 0;
						for (kk = 0; kk < nneighborlist; kk++) {

							// Determine the total signal.
							event->signal += neighborlist[kk]->signal;
							// If a contribution was negative, flag as invalid
							// (-2, such that it doesn't collide with definition afterwards.
							// Is changed to -1 at the end of the process.)
							if (neighborlist[kk]->signal < 0.) {
								event->type = -2;
							} else {
								event->type = -1;
							}

							// Determine signals in 3x3 matrix.
							if (neighborlist[kk]->rawx == neighborlist[maxidx]->rawx - 1) {
								if (neighborlist[kk]->rawy == neighborlist[maxidx]->rawy - 1) {
									event->signals[0] = neighborlist[kk]->signal;
									event->phas[0] = neighborlist[kk]->pha;
								} else if (neighborlist[kk]->rawy
										== neighborlist[maxidx]->rawy) {
									event->signals[3] = neighborlist[kk]->signal;
									event->phas[3] = neighborlist[kk]->pha;
								} else if (neighborlist[kk]->rawy
										== neighborlist[maxidx]->rawy + 1) {
									event->signals[6] = neighborlist[kk]->signal;
									event->phas[6] = neighborlist[kk]->pha;
								}
							} else if (neighborlist[kk]->rawx == neighborlist[maxidx]->rawx) {
								if (neighborlist[kk]->rawy == neighborlist[maxidx]->rawy - 1) {
									event->signals[1] = neighborlist[kk]->signal;
									event->phas[1] = neighborlist[kk]->pha;
								} else if (neighborlist[kk]->rawy
										== neighborlist[maxidx]->rawy) {
									event->signals[4] = neighborlist[kk]->signal;
									event->phas[4] = neighborlist[kk]->pha;
								} else if (neighborlist[kk]->rawy
										== neighborlist[maxidx]->rawy + 1) {
									event->signals[7] = neighborlist[kk]->signal;
									event->phas[7] = neighborlist[kk]->pha;
								}
							} else if (neighborlist[kk]->rawx
									== neighborlist[maxidx]->rawx + 1) {
								if (neighborlist[kk]->rawy == neighborlist[maxidx]->rawy - 1) {
									event->signals[2] = neighborlist[kk]->signal;
									event->phas[2] = neighborlist[kk]->pha;
								} else if (neighborlist[kk]->rawy
										== neighborlist[maxidx]->rawy) {
									event->signals[5] = neighborlist[kk]->signal;
									event->phas[5] = neighborlist[kk]->pha;
								} else if (neighborlist[kk]->rawy
										== neighborlist[maxidx]->rawy + 1) {
									event->signals[8] = neighborlist[kk]->signal;
									event->phas[8] = neighborlist[kk]->pha;
								}
							}

							// Set PH_IDs and SRC_IDs.
							long ll;
							for (ll = 0; ll < NEVENTPHOTONS; ll++) {
								if (0 == neighborlist[kk]->ph_id[ll])
									break;
								long mm;
								for (mm = 0; mm < NEVENTPHOTONS; mm++) {
									if (event->ph_id[mm] == neighborlist[kk]->ph_id[ll])
										break;
									if (0 == event->ph_id[mm]) {
										event->ph_id[mm] = neighborlist[kk]->ph_id[ll];
										event->src_id[mm] = neighborlist[kk]->src_id[ll];
										break;
									}
								}
							}

							// Check for border pixels.
							if ((0 == neighborlist[kk]->rawx)
									|| (neighborlist[kk]->rawx == det->pixgrid->xwidth - 1)
									|| (det->rawymin == neighborlist[kk]->rawy)
									|| (neighborlist[kk]->rawy == det->rawymax)) {
								border = 1;
							}
						}
						// END of loop over all entries in the neighbor list.


						// Determine the PHA channel corresponding to the total signal.
						if (NULL != det->rmf) {
							event->pha = getEBOUNDSChannel(event->signal, det->rmf);
						} else {
							event->pha = 0;
						}

						// Check for pile-up.
						if (NEVENTPHOTONS >= 2) {
							if (0 != event->ph_id[1]) {
								event->pileup = 1;
							}
						}

						// Determine the event type.
						if (event->type == -2) {
							//Event had negative contributions, flag as invalid.
							event->type = -1;
						} else {
							// First assume that the event is invalid.
							event->type = -1;
							// Border events are declared as invalid.
							if (0 == border) {
								if (1 == nneighborlist) {
									// Single event.
									event->type = 0;

								} else if (2 == nneighborlist) {
									// Check for double types.
									if (event->signals[1] > 0.) {
										event->type = 3; // bottom
									} else if (event->signals[3] > 0.) {
										event->type = 4; // left
									} else if (event->signals[7] > 0.) {
										event->type = 1; // top
									} else if (event->signals[5] > 0.) {
										event->type = 2; // right
									}

								} else if (3 == nneighborlist) {
									// Check for triple types.
									if (event->signals[1] > 0.) {
										// bottom
										if (event->signals[3] > 0.) {
											event->type = 7; // bottom-left
										} else if (event->signals[5] > 0.) {
											event->type = 6; // bottom-right
										}
									} else if (event->signals[7] > 0.) {
										// top
										if (event->signals[3] > 0.) {
											event->type = 8; // top-left
										} else if (event->signals[5] > 0.) {
											event->type = 5; // top-right
										}
									}

								} else if (4 == nneighborlist) {
									// Check for quadruple types.
									if (event->signals[0] > 0.) { // bottom-left
										if ((event->signals[1] > event->signals[0])
												&& (event->signals[3] > event->signals[0])) {
											event->type = 11;
										}
									} else if (event->signals[2] > 0.) { // bottom-right
										if ((event->signals[1] > event->signals[2])
												&& (event->signals[5] > event->signals[2])) {
											event->type = 10;
										}
									} else if (event->signals[6] > 0.) { // top-left
										if ((event->signals[7] > event->signals[6])
												&& (event->signals[3] > event->signals[6])) {
											event->type = 12;
										}
									} else if (event->signals[8] > 0.) { // top-right
										if ((event->signals[7] > event->signals[8])
												&& (event->signals[5] > event->signals[8])) {
											event->type = 9;
										}
									}
								}
							}
						}
						// END of determine the event type.

						// Remove processed events from neighbor list.
						for (kk = 0; kk < nneighborlist; kk++) {
							freeEvent(&neighborlist[kk]);
							neighborlist[kk] = NULL;
						}
						nneighborlist = 0;

						// TODO: PHA2PI CORRECTION
						pha2picorrect( event, p2p , status);
						CHECK_STATUS_BREAK(*status);

						// Check if the total signal of the event is below
						// the upper event threshold.
						if ((det->threshold_pattern_up_keV == 0.)
								|| (event->signal <= det->threshold_pattern_up_keV)) {
							// Update the event statistics.
							if (event->type < 0) {
								statistics.ninvalids++;
								if (event->pileup > 0) {
									statistics.npinvalids++;
								}
							} else {
								statistics.nvalids++;
								statistics.ngrade[event->type]++;
								if (event->pileup > 0) {
									statistics.npvalids++;
									statistics.npgrade[event->type]++;
								}
							}

							// If the event is invalid, check if it should be
							// added to the output file or not.
							if ((0 == skip_invalids) || (event->type >= 0)) {
								// Add the new event to the output file.
								addEvent2File(dest, event, status);
								CHECK_STATUS_BREAK(*status);
							}
						} // End of application of upper threshold.

						// Release memory.
						freeEvent(&event);
					}
				}
				CHECK_STATUS_BREAK(*status);
				// END of loop over all events in the frame list.

				// Delete all remaining events in the frame list.
				// There might still be some, which are below the
				// thresholds.
				for (jj = 0; jj < nframelist; jj++) {
					if (NULL != framelist[jj]) {
						freeEvent(&framelist[jj]);
					}
				}
				nframelist = 0;
			}
			// END of if new frame.

			// Append the new event to the frame list.
			if (NULL != event) {
				if (nframelist >= maxnframelist) {
					SIXT_ERROR("too many events in the same frame");
					*status = EXIT_FAILURE;
					break;
				}
				framelist[nframelist] = event;
				nframelist++;
			}

		}
		CHECK_STATUS_BREAK(*status);
		// END of loop over all events in the input file.

		// Store pattern statistics in the output file.
		// Valids.
		fits_update_key(dest->fptr, TLONG, "NVALID", &statistics.nvalids,
				"number of valid patterns", status);
		fits_update_key(dest->fptr, TLONG, "NPVALID", &statistics.npvalids,
				"number of piled up valid patterns", status);
		// Invalids.
		fits_update_key(dest->fptr, TLONG, "NINVALID", &statistics.ninvalids,
				"number of invalid patterns", status);
		fits_update_key(dest->fptr, TLONG, "NPINVALI", &statistics.npinvalids,
				"number of piled up invalid patterns", status);
		// Numbered grades.
		for (ii = 0; ii < 13; ii++) {
			char keyword[MAXMSG];
			char comment[MAXMSG];
			sprintf(keyword, "NGRAD%ld", ii);
			sprintf(comment, "number of patterns with grade %ld", ii);
			fits_update_key(dest->fptr, TLONG, keyword, &statistics.ngrade[ii],
					comment, status);
			sprintf(keyword, "NPGRA%ld", ii);
			sprintf(comment, "number of piled up patterns with grade %ld", ii);
			fits_update_key(dest->fptr, TLONG, keyword, &statistics.npgrade[ii],
					comment, status);
		}
		// Set Pha2Pi keyword if correction was performed
		if(p2p != NULL){
			fits_update_key(dest->fptr, TSTRING, "PHA2PI", p2p->p2p_filename,
					"Pha2Pi file - set if correction was performed", status);
		}
		CHECK_STATUS_BREAK(*status);

	} while (0); // End of error handling loop.

	// Release memory.
	freePha2Pi(&p2p);

	if (NULL != framelist) {
		for (ii = 0; ii < nframelist; ii++) {
			if (NULL != framelist[ii]) {
				freeEvent(&framelist[ii]);
			}
		}
		free(framelist);
	}
	if (NULL != neighborlist) {
		for (ii = 0; ii < nneighborlist; ii++) {
			if (NULL != neighborlist[ii]) {
				freeEvent(&neighborlist[ii]);
			}
		}
		free(neighborlist);
	}
}

Pha2Pi* getPha2Pi(int* const status) {
	// Allocate memory.
	Pha2Pi* p2p = (Pha2Pi*) malloc(sizeof(Pha2Pi));
	CHECK_NULL_RET(p2p, *status, "memory allocation for Pha2Pi failed", p2p);

	// Initialize.
	p2p->p2p_filename = NULL;
	p2p->randgen = NULL;
	p2p->nrows = 0;
	p2p->ngrades = 0;
	p2p->pha = NULL;
	p2p->pien = NULL;
	p2p->pilow = NULL;
	p2p->pihigh = NULL;

	return (p2p);
}

void freePha2Pi(Pha2Pi** const p2p) {
	if (NULL != *p2p) {
		if (NULL != (*p2p)->p2p_filename) {
			free((*p2p)->p2p_filename);
		}
		if (NULL != (*p2p)->randgen) {
			gsl_rng_free((*p2p)->randgen);
		}
		if (NULL != (*p2p)->pha) {
			free((*p2p)->pha);
		}
		if (NULL != (*p2p)->pien) {
			for (int ii=0; ii<(*p2p)->nrows; ++ii){
				free((*p2p)->pien[ii]);
			}
			free((*p2p)->pien);
		}
		if (NULL != (*p2p)->pilow) {
			for (int ii=0; ii<(*p2p)->nrows; ++ii){
				free((*p2p)->pilow[ii]);
			}
			free((*p2p)->pilow);
		}
		if (NULL != (*p2p)->pihigh) {
			for (int ii=0; ii<(*p2p)->nrows; ++ii){
				free((*p2p)->pihigh[ii]);
			}
			free((*p2p)->pihigh);
		}
		free(*p2p);
		*p2p = NULL;
	}
}

Pha2Pi* initPha2Pi(const char* const filename, const unsigned int seed, int* const status) {
	if (strlen(filename) == 0) {
		return NULL;
	} else if (access(filename, F_OK) != 0) {
		char msg[MAXMSG];
		sprintf(msg, "Pha2Pi correction File '%s' not accessible!", filename);
		SIXT_WARNING(msg);
		return NULL;
	}

	int colnum, typecode;
	long width;
	int anynul = 0;

	Pha2Pi* p2p = getPha2Pi(status);
	CHECK_STATUS_RET(*status, p2p);

	/** Store filename */
	p2p->p2p_filename=(char*)malloc((strlen(filename)+1)*sizeof(char));
	CHECK_NULL_RET(p2p->p2p_filename, *status,
			"memory allocation for p2p_filename failed", p2p);
	strcpy(p2p->p2p_filename, filename);

	/** INITIALIZE RANDOM NUMBER GENERATOR */
	p2p->randgen=gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(p2p->randgen,seed);

	/** LOAD FILE */
	headas_chat(3, "open Pha2Pi file '%s' ...\n", filename);
	fitsfile* fptr;
	fits_open_table(&fptr, filename, READONLY, status);
	CHECK_STATUS_RET(*status, p2p);

	// Determine the number of rows.
	fits_get_num_rows(fptr, &p2p->nrows, status);

	// Load PHA column
	p2p->pha = (long *) malloc(p2p->nrows * sizeof(long));
	fits_get_colnum(fptr, CASEINSEN, "PHA", &colnum, status);
	fits_read_col(fptr, TINT32BIT, colnum, 1, 1, p2p->nrows, NULL, p2p->pha,
			&anynul, status);
	CHECK_STATUS_RET(*status, p2p);

	// CHECK CONSISTENCY: PHA has to start at zero and must be continues
	if( p2p->pha[0] != 0 || p2p->pha[p2p->nrows-1] != p2p->nrows-1 ){
		char msg[MAXMSG];
		sprintf(msg, "Pha2Pi correction File is inconsistent! PHA values must be continues starting with 0!");
		SIXT_ERROR(msg);
		return NULL;
	}

	// Determine the number of grades element dimension.
	fits_get_colnum(fptr, CASEINSEN, "PIEN", &colnum, status);
	fits_get_coltype(fptr, colnum, &typecode, &p2p->ngrades, &width, status);
	CHECK_STATUS_RET(*status, p2p);

	// Load PIEN, PILOW, PIHIGH columns
	p2p->pien = (double **) malloc(p2p->nrows * sizeof(double*));
	for (int ii=0; ii<p2p->nrows; ++ii){
		p2p->pien[ii] = (double *) malloc(p2p->ngrades * sizeof(double));
		fits_read_col(fptr, TDOUBLE, colnum, ii+1, 1, p2p->ngrades, NULL, p2p->pien[ii],
				&anynul, status);
		CHECK_STATUS_RET(*status, p2p);
	}
	fits_get_colnum(fptr, CASEINSEN, "PILOW", &colnum, status);
	p2p->pilow = (double **) malloc(p2p->nrows * sizeof(double*));
	for (int ii=0; ii<p2p->nrows; ++ii){
		p2p->pilow[ii] = (double *) malloc(p2p->ngrades * sizeof(double));
		fits_read_col(fptr, TDOUBLE, colnum, ii+1, 1, p2p->ngrades, NULL, p2p->pilow[ii],
				&anynul, status);
		CHECK_STATUS_RET(*status, p2p);
	}
	fits_get_colnum(fptr, CASEINSEN, "PIHIGH", &colnum, status);
	p2p->pihigh = (double **) malloc(p2p->nrows * sizeof(double*));
	for (int ii=0; ii<p2p->nrows; ++ii){
		p2p->pihigh[ii] = (double *) malloc(p2p->ngrades * sizeof(double));
		fits_read_col(fptr, TDOUBLE, colnum, ii+1, 1, p2p->ngrades, NULL, p2p->pihigh[ii],
				&anynul, status);
		CHECK_STATUS_RET(*status, p2p);
	}

	return (p2p);
}




void pha2picorrect(Event* const evt, const Pha2Pi* const p2p, int* const status) {

	// Do nothing if the Pha2Pi structure is uninitialized
	if( p2p == NULL ){
		return;
	}

	// Do nothing if event is invalid => pi = -1.
	if( evt->type == -1
			|| evt->pha < 0
			|| evt->pha > p2p->pha[p2p->nrows-1]){
		return;
	}
	// Make sure the requested evt->type is tabulated
	else if( evt->type < 0 || evt->type > p2p->ngrades ){
		char msg[MAXMSG];
		sprintf(msg, "Pha2Pi correction event type '%d' not tabulated!", evt->type);
		SIXT_WARNING(msg);
		return;
	}
	// Determine a PI value for the event's PHA value
	else{
	  double ran = gsl_rng_uniform(p2p->randgen);
		const double emin = p2p->pilow[evt->pha][evt->type];
		const double emax = p2p->pihigh[evt->pha][evt->type];
		evt->pi = emin + ran*(emax-emin);
	}

}

