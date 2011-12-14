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
  if (((e1->rawx==e2->rawx+1)&&(e1->rawy==e2->rawy)) ||
      ((e1->rawx==e2->rawx-1)&&(e1->rawy==e2->rawy)) ||
      ((e1->rawx==e2->rawx)&&(e1->rawy==e2->rawy+1)) ||
      ((e1->rawx==e2->rawx)&&(e1->rawy==e2->rawy-1))) {
    // The events are neighbors.
    return(1);
  } else {
    return(0);
  }
}


void phpat(GenDet* const det,
	   EventListFile* const elf,
	   PatternFile* const plf,
	   int* const status)
{
  // Pattern / grade statistics.
  struct PatternStatistics statistics = {
    .nvalids    = 0,
    .npvalids   = 0,
    .ninvalids  = 0,
    .npinvalids = 0,
  };
  long ii;
  for (ii=0; ii<13; ii++) {
    statistics.ngrade[ii]  = 0;
    statistics.npgrade[ii] = 0;
  }

  // List of all events belonging to the current frame.
  const long maxnframelist=10000;
  Event** framelist=NULL;
  long nframelist=0;

  // List of all neighboring events in the current frame.
  const long maxnneighborlist=1000;
  Event** neighborlist=NULL;
  long nneighborlist=0;


  // Error handling loop.
  do {

    // Allocate memory.
    framelist=(Event**)malloc(maxnframelist*sizeof(Event*));
    CHECK_NULL_BREAK(framelist, *status, 
		     "memory allocation for frame list failed");
    neighborlist=(Event**)malloc(maxnneighborlist*sizeof(Event*));
    CHECK_NULL_BREAK(neighborlist, *status, 
		     "memory allocation for neighbor list failed");

    // Determine the name of the instrument.
    // Particular instruments require a special pattern
    // recombination scheme (e.g. eROSITA).
    char telescope[MAXMSG];
    char comment[MAXMSG];
    fits_read_key(elf->fptr, TSTRING, "TELESCOP", 
		  telescope, comment, status);
    CHECK_STATUS_BREAK(*status);
    strtoupper(telescope);
      

    // Loop over all events in the input list.
    for (ii=0; ii<=elf->nrows; ii++) {
      
      // Read the next event from the file, if the end has not been
      // reached so far.
      Event* event=NULL;
      if (ii<elf->nrows) {
	event=getEvent(status);
	CHECK_STATUS_BREAK(*status);
	getEventFromFile(elf, ii+1, event, status);
	CHECK_STATUS_BREAK(*status);
      }

      // Check if the new event belongs to a different frame than 
      // the previous ones.
      int newframe=0;
      if ((nframelist>0) && (NULL!=event)) {
	if (event->frame!=framelist[0]->frame) {
	  newframe=1;
	}
      }

      // If the new event belongs to a different frame than the
      // previous ones or the end of the input list has been
      // reached, perform a pattern analysis.
      if ((newframe>0) || ((ii==elf->nrows)&&(nframelist>0))) {

	// Loop over all events in the current frame.
	long jj;
	for (jj=0; jj<nframelist; jj++) {
	  if (NULL!=framelist[jj]) {
	    
	    // Check if the event is below the threshold.
	    if (framelist[jj]->signal<det->threshold_event_lo_keV) continue;
  
	    // Start a new neighbor list.
	    neighborlist[0]=framelist[jj];
	    nneighborlist=1;
	    framelist[jj]=NULL;
	    
	    // Find the signal maximum in the neighboring pixels.
	    Event* maxsignalev=neighborlist[0];
	    int updated=0;
	    do {
	      updated=0;
	      long ll;
	      for (ll=0; ll<nframelist; ll++) {
		if (NULL!=framelist[ll]) {
		  if (isNeighbor(maxsignalev, framelist[ll])) {
		    if (framelist[ll]->signal>maxsignalev->signal) {
		      maxsignalev=framelist[ll];
		      updated=1;
		    }
		  }
		}
	      }
	    } while(updated);

	    // Determine the split threshold [keV].
	    float split_threshold;
	    // For eROSITA we need a special treatment (according to
	    // a prescription of K. Dennerl).
	    if (!strcmp(telescope, "EROSITA")) {
	      if (det->threshold_split_lo_fraction > 0.) {
		float vertical=0., horizontal=0.;
		long ll;
		for (ll=0; ll<nframelist; ll++) {
		  if (NULL!=framelist[ll]) {
		    if (isNeighbor(maxsignalev, framelist[ll])) {
		      if (framelist[ll]->rawx==maxsignalev->rawx) {
			if (framelist[ll]->signal>horizontal) {
			  horizontal=framelist[ll]->signal;
			}
		      } else {
			if (framelist[ll]->signal>vertical) {
			  vertical=framelist[ll]->signal;
			}
		      }
		    }
		  }
		}
		split_threshold=det->threshold_split_lo_fraction*
		  (maxsignalev->signal+horizontal+vertical);
	      } else {
		SIXT_ERROR("eROSITA requires a fractional split threshold");
		*status=EXIT_FAILURE;
		break;
	      }

	    } else { // Split threshold for generic instruments.
	      if (det->threshold_split_lo_fraction > 0.) {
		split_threshold=
		  det->threshold_split_lo_fraction*maxsignalev->signal;
	      } else {
		split_threshold=det->threshold_split_lo_keV;
	      }
	    }
	    // END of determine the split threshold.

	    // Check if the split threshold is above the event threshold.
	    if (split_threshold > det->threshold_event_lo_keV) {
	      printf("*** warning: split threshold is above event threshold\n");
	    }

	    // Find all neighboring events above the split threshold.
 	    long kk;
	    for (kk=0; kk<nneighborlist; kk++) {
	      long ll;
	      for (ll=0; ll<nframelist; ll++) {
		if (NULL!=framelist[ll]) {
		  if (isNeighbor(neighborlist[kk], framelist[ll])) {

		    // Check if its signal is below the split threshold.
		    if (framelist[ll]->signal<split_threshold) {
		      continue;
		    }
		    
		    // Add the event to the neighbor list.
		    if (nneighborlist>=maxnneighborlist) {
		      SIXT_ERROR("too many events in the same pattern");
		      *status=EXIT_FAILURE;
		      break;
		    }
		    neighborlist[nneighborlist]=framelist[ll];
		    nneighborlist++;
		    framelist[ll]=NULL;
		  }
		}
	      }
	      CHECK_STATUS_BREAK(*status);
	    }
	    CHECK_STATUS_BREAK(*status);
	    // END of finding all neighbors.
	    
	    // Search the pixel with the maximum signal.
	    long maxidx=0;
	    for (kk=1; kk<nneighborlist; kk++) {
	      if (neighborlist[kk]->signal>neighborlist[maxidx]->signal) {
		maxidx=kk;
	      }
	    }
	    // END of searching the pixel with the maximum signal.

	    // Get a new pattern.
	    Pattern* pattern=getPattern(status);
	    CHECK_STATUS_BREAK(*status);
	    
	    // Set basic properties.
	    pattern->rawx   =neighborlist[maxidx]->rawx;
	    pattern->rawy   =neighborlist[maxidx]->rawy;
	    pattern->time   =neighborlist[maxidx]->time;
	    pattern->frame  =neighborlist[maxidx]->frame;
	    pattern->ra     =0.;
	    pattern->dec    =0.;
	    pattern->npixels=nneighborlist;

	    // Set the advanced properties.
	    // Total signal.
	    pattern->signal=0.;
	    // Flag whether pattern touches the border of the detector.
	    int border=0;
	    for (kk=0; kk<nneighborlist; kk++) {
	      
	      // Determine the total signal.
	      pattern->signal+=neighborlist[kk]->signal;

	      // Determine signals in 3x3 matrix.
	      if (neighborlist[kk]->rawx==neighborlist[maxidx]->rawx-1) {
		if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy-1) {
		  pattern->signals[0] = neighborlist[kk]->signal;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy) {
		  pattern->signals[3] = neighborlist[kk]->signal;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy+1) {
		  pattern->signals[6] = neighborlist[kk]->signal;
		}
	      } else if (neighborlist[kk]->rawx==neighborlist[maxidx]->rawx) {
		if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy-1) {
		  pattern->signals[1] = neighborlist[kk]->signal;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy) {
		  pattern->signals[4] = neighborlist[kk]->signal;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy+1) {
		  pattern->signals[7] = neighborlist[kk]->signal;
		}
	      } else if (neighborlist[kk]->rawx==neighborlist[maxidx]->rawx+1) {
		if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy-1) {
		  pattern->signals[2] = neighborlist[kk]->signal;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy) {
		  pattern->signals[5] = neighborlist[kk]->signal;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy+1) {
		  pattern->signals[8] = neighborlist[kk]->signal;
		}
	      }

	      // Set PH_IDs and SRC_IDs.
	      long ll;
	      for (ll=0; ll<NEVENTPHOTONS; ll++) {
		if (0==neighborlist[kk]->ph_id[ll]) break;
		long mm;
		for (mm=0; mm<NPATTERNPHOTONS; mm++) {
		  if (pattern->ph_id[mm]==neighborlist[kk]->ph_id[ll]) break;
		  if (0==pattern->ph_id[mm]) {
		    pattern->ph_id[mm] =neighborlist[kk]->ph_id[ll];
		    pattern->src_id[mm]=neighborlist[kk]->src_id[ll];
		    break;
		  }
		}
	      }

	      // Check for border pixels.
	      if ((0==neighborlist[kk]->rawx)||
		  (neighborlist[kk]->rawx==det->pixgrid->xwidth-1)||
		  (0==neighborlist[kk]->rawy)||
		  (neighborlist[kk]->rawy==det->pixgrid->ywidth-1)) {
		border=1;
	      }
	    }
	    // END of loop over all entries in the neighbor list.

	    // Determine the PHA corresponding to the total signal.
	    if (NULL!=det->rmf) {
	      pattern->pha=getEBOUNDSChannel(pattern->signal, det->rmf);
	    } else {
	      pattern->pha=0;
	    }

	    // Check for pile-up.
	    if (NPATTERNPHOTONS>=2) {
	      if (0!=pattern->ph_id[1]) {
		pattern->pileup=1;
	      }
	    }

	    // Determine the pattern type.
	    // First assume that the pattern is invalid.
	    pattern->type=-1; 
	    // Border events are declared as invalid.
	    if (0==border) {
	      if (1==nneighborlist) {
		// Single event.
		pattern->type=0;

	      } else if (2==nneighborlist) {
		// Check for double types.
		if (pattern->signals[1]>0.) {
		  pattern->type=1; // bottom
		} else if (pattern->signals[3]>0.) {
		  pattern->type=2; // left
		} else if (pattern->signals[7]>0.) {
		  pattern->type=3; // top
		} else if (pattern->signals[5]>0.) {
		  pattern->type=4; // right
		} 

	      } else if (3==nneighborlist) {
		// Check for triple types.
		if (pattern->signals[1]>0.) {
		  // bottom
		  if (pattern->signals[3]>0.) {
		    pattern->type=5; // bottom-left
		  } else if (pattern->signals[5]>0.) {
		    pattern->type=6; // bottom-right
		  }
		} else if (pattern->signals[7]>0.) {
		  // top
		  if (pattern->signals[3]>0.) {
		    pattern->type=7; // top-left
		  } else if (pattern->signals[5]>0.) {
		    pattern->type=8; // top-right
		  }
		}

	      } else if (4==nneighborlist) {
		// Check for quadruple types.
		if (pattern->signals[0]>0.) { // bottom-left
		  if ((pattern->signals[1]>pattern->signals[0])&&
		      (pattern->signals[3]>pattern->signals[0])) {
		    pattern->type=9; 
		  } 
		} else if (pattern->signals[2]>0.) { // bottom-right
		  if ((pattern->signals[1]>pattern->signals[2])&&
		      (pattern->signals[5]>pattern->signals[2])) {
		    pattern->type=10;
		  }
		} else if (pattern->signals[6]>0.) { // top-left
		  if ((pattern->signals[7]>pattern->signals[6])&&
		      (pattern->signals[3]>pattern->signals[6])) {
		    pattern->type=11; 
		  }
		} else if (pattern->signals[8]>0.) { // top-right
		  if ((pattern->signals[7]>pattern->signals[8])&&
		      (pattern->signals[5]>pattern->signals[8])) {
		    pattern->type=12; 
		  }
		} 
	      }
	    }
	    // END of determine the pattern type.

	    // Update the pattern statistics.
	    if (pattern->type<0) {
	      statistics.ninvalids++;
	      if (pattern->pileup>0) {
		statistics.npinvalids++;
	      }
	    } else {
	      statistics.nvalids++;
	      statistics.ngrade[pattern->type]++;
	      if (pattern->pileup>0) {
		statistics.npvalids++;
		statistics.npgrade[pattern->type]++;
	      }
	    }	      

	    // Remove processed events from neighbor list.
	    for (kk=0; kk<nneighborlist; kk++) {
	      freeEvent(&neighborlist[kk]);
	      neighborlist[kk]=NULL;
	    }
	    nneighborlist=0;

	    // Add the new pattern to the output file.
	    addPattern2File(plf, pattern, status);	  
	    CHECK_STATUS_BREAK(*status);

	    // Release memory.
	    freePattern(&pattern);
	  }
	}
	CHECK_STATUS_BREAK(*status);
	// END of loop over all events in the frame list.

	// Delete all remaining events in the frame list.
	// There might still be some, which are below the
	// thresholds.
	for (jj=0; jj<nframelist; jj++) {
	  if (NULL!=framelist[jj]) {
	    freeEvent(&framelist[jj]);
	  }
	}
	nframelist=0;
      }
      // END of if new frame.

      // Append the new event to the frame list.
      if (NULL!=event) {
	if (nframelist>=maxnframelist) {
	  SIXT_ERROR("too many events in the same frame");
	  *status=EXIT_FAILURE;
	  break;
	}
	framelist[nframelist]=event;
	nframelist++;
      }

    }    
    CHECK_STATUS_BREAK(*status);
    // END of loop over all events in the input file.

    // Store pattern statistics in the output file.
    // Valids.
    fits_update_key(plf->fptr, TLONG, "NVALID", 
		    &statistics.nvalids, 
		    "number of valid patterns", status);
    fits_update_key(plf->fptr, TLONG, "NPVALID", 
		    &statistics.npvalids, 
		    "number of piled up valid patterns", status);
    // Invalids.
    fits_update_key(plf->fptr, TLONG, "NINVALID", 
		    &statistics.ninvalids, 
		    "number of invalid patterns", status);
    fits_update_key(plf->fptr, TLONG, "NPINVALI", 
		    &statistics.npinvalids, 
		    "number of piled up invalid patterns", status);
    // Numbered grades.
    for (ii=0; ii<13; ii++) {
      char keyword[MAXMSG];
      char comment[MAXMSG];
      sprintf(keyword, "NGRAD%ld", ii);
      sprintf(comment, "number of patterns with grade %ld", ii);
      fits_update_key(plf->fptr, TLONG, keyword, 
		      &statistics.ngrade[ii], comment, status);
      sprintf(keyword, "NPGRA%ld", ii);
      sprintf(comment, "number of piled up patterns with grade %ld", ii);
      fits_update_key(plf->fptr, TLONG, keyword, 
		      &statistics.npgrade[ii], comment, status);    
    }
    CHECK_STATUS_BREAK(*status);
    
  } while(0); // End of error handling loop.
  

  // Release memory.
  if (NULL!=framelist) {
    for (ii=0; ii<nframelist; ii++) {
      if (NULL!=framelist[ii]) {
	freeEvent(&framelist[ii]);
      }
    }
    free(framelist);
  }
  if (NULL!=neighborlist) {
    for (ii=0; ii<nneighborlist; ii++) {
      if (NULL!=neighborlist[ii]) {
	freeEvent(&neighborlist[ii]);
      }
    }
    free(neighborlist);
  }
}

