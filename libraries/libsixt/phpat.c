#include "phpat.h"

// TODO Apply thresholds.

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
    CHECK_NULL_BREAK(framelist, *status, "memory allocation for frame list failed");
    neighborlist=(Event**)malloc(maxnneighborlist*sizeof(Event*));
    CHECK_NULL_BREAK(neighborlist, *status, "memory allocation for neighbor list failed");
  
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
	    // Start a new neighbor list.
	    neighborlist[0]=framelist[jj];
	    nneighborlist=1;
	    framelist[jj]=NULL;

	    // Find all neighboring events.
	    long kk;
	    for (kk=0; kk<nneighborlist; kk++) {
	      long ll;
	      for (ll=jj+1; ll<nframelist; ll++) {
		if (NULL!=framelist[ll]) {
		  if (((neighborlist[kk]->rawx==framelist[ll]->rawx+1)&&
		       (neighborlist[kk]->rawy==framelist[ll]->rawy)) ||
		      ((neighborlist[kk]->rawx==framelist[ll]->rawx-1)&&
		       (neighborlist[kk]->rawy==framelist[ll]->rawy)) ||
		      ((neighborlist[kk]->rawx==framelist[ll]->rawx)&&
		       (neighborlist[kk]->rawy==framelist[ll]->rawy+1)) ||
		      ((neighborlist[kk]->rawx==framelist[ll]->rawx)&&
		       (neighborlist[kk]->rawy==framelist[ll]->rawy-1))) {
		    // This is a neighbor.
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
	      if (neighborlist[kk]->charge>neighborlist[maxidx]->charge) {
		maxidx=kk;
	      }
	    }
	    // END of searching the pixel with the maximum signal.

	    // Get a new pattern.
	    Pattern* pattern=getPattern(status);
	    CHECK_STATUS_BREAK(*status);
	    
	    // Set basic properties.
	    pattern->event->rawx =neighborlist[maxidx]->rawx;
	    pattern->event->rawy =neighborlist[maxidx]->rawy;
	    pattern->event->time =neighborlist[maxidx]->time;
	    pattern->event->frame=neighborlist[maxidx]->frame;
	    pattern->event->ra   =neighborlist[maxidx]->ra;
	    pattern->event->dec  =neighborlist[maxidx]->dec;
	    pattern->npixels     =nneighborlist;

	    // Set the advanced properties.
	    // Total signal.
	    pattern->event->charge=0.;
	    // Flag whether pattern touches the border of the detector.
	    int border=0;
	    for (kk=0; kk<nneighborlist; kk++) {
	      
	      // Determine the total signal.
	      pattern->event->charge+=neighborlist[kk]->charge;

	      // Determine signals in 3x3 matrix.
	      if (neighborlist[kk]->rawx==neighborlist[maxidx]->rawx-1) {
		if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy-1) {
		  pattern->signals[0] = neighborlist[kk]->charge;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy) {
		  pattern->signals[3] = neighborlist[kk]->charge;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy+1) {
		  pattern->signals[6] = neighborlist[kk]->charge;
		}
	      } else if (neighborlist[kk]->rawx==neighborlist[maxidx]->rawx) {
		if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy-1) {
		  pattern->signals[1] = neighborlist[kk]->charge;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy) {
		  pattern->signals[4] = neighborlist[kk]->charge;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy+1) {
		  pattern->signals[7] = neighborlist[kk]->charge;
		}
	      } else if (neighborlist[kk]->rawx==neighborlist[maxidx]->rawx+1) {
		if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy-1) {
		  pattern->signals[2] = neighborlist[kk]->charge;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy) {
		  pattern->signals[5] = neighborlist[kk]->charge;
		} else if (neighborlist[kk]->rawy==neighborlist[maxidx]->rawy+1) {
		  pattern->signals[8] = neighborlist[kk]->charge;
		}
	      }

	      // Set PH_IDs and SRC_IDs.
	      long ll;
	      for (ll=0; ll<NEVENTPHOTONS; ll++) {
		if (0==neighborlist[kk]->ph_id[ll]) break;
		long mm;
		for (mm=0; mm<NEVENTPHOTONS; mm++) {
		  if (pattern->event->ph_id[mm]==neighborlist[kk]->ph_id[ll]) break;
		  if (0==pattern->event->ph_id[mm]) {
		    pattern->event->ph_id[mm] =neighborlist[kk]->ph_id[ll];
		    pattern->event->src_id[mm]=neighborlist[kk]->src_id[ll];
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
	    pattern->event->pha=getEBOUNDSChannel(pattern->event->charge, det->rmf);

	    // Check for pile-up.
	    if (NEVENTPHOTONS>=2) {
	      if (0!=pattern->event->ph_id[1]) {
		pattern->pileup=1;
	      }
	    }

	    // Determine the pattern type.
	    // First assume that the pattern is invalid.
	    pattern->type = -1; 
	    // Border events are declared as invalid.
	    if (0==border) {
	      if (1==nneighborlist) {
		// Single event.
		pattern->type = 0;

	      } else if (2==nneighborlist) {
		// Check for double types.
		if (pattern->signals[1]>0.) {
		  pattern->type = 1; // bottom
		} else if (pattern->signals[3]>0.) {
		  pattern->type = 2; // left
		} else if (pattern->signals[7]>0.) {
		  pattern->type = 3; // top
		} else if (pattern->signals[5]>0.) {
		  pattern->type = 4; // right
		} 

	      } else if (3==nneighborlist) {
		// Check for triple types.
		if (pattern->signals[1]>0.) {
		  // bottom
		  if (pattern->signals[3]>0.) {
		    pattern->type = 5; // bottom-left
		  } else if (pattern->signals[5]>0.) {
		    pattern->type = 6; // bottom-right
		  }
		} else if (pattern->signals[7]>0.) {
		  // top
		  if (pattern->signals[3]>0.) {
		    pattern->type = 7; // top-left
		  } else if (pattern->signals[5]>0.) {
		    pattern->type = 8; // top-right
		  }
		}

	      } else if (4==nneighborlist) {
		// Check for quadruple types.
		if (pattern->signals[0]>0.) { // bottom-left
		  if ((pattern->signals[1]>pattern->signals[0])&&
		      (pattern->signals[3]>pattern->signals[0])) {
		    pattern->type = 9; 
		  } 
		} else if (pattern->signals[2]>0.) { // bottom-right
		  if ((pattern->signals[1]>pattern->signals[2])&&
		      (pattern->signals[5]>pattern->signals[2])) {
		    pattern->type = 10;
		  }
		} else if (pattern->signals[6]>0.) { // top-left
		  if ((pattern->signals[7]>pattern->signals[6])&&
		      (pattern->signals[3]>pattern->signals[6])) {
		    pattern->type = 11; 
		  }
		} else if (pattern->signals[8]>0.) { // top-right
		  if ((pattern->signals[7]>pattern->signals[8])&&
		      (pattern->signals[5]>pattern->signals[8])) {
		    pattern->type = 12; 
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

	// Now the frame list is empty, too.
	nframelist=0;
      }

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
    fits_update_key(plf->eventlistfile->fptr, TLONG, "NVALID", 
		    &statistics.nvalids, 
		    "number of valid patterns", status);
    fits_update_key(plf->eventlistfile->fptr, TLONG, "NPVALID", 
		    &statistics.npvalids, 
		    "number of piled up valid patterns", status);
    // Invalids.
    fits_update_key(plf->eventlistfile->fptr, TLONG, "NINVALID", 
		    &statistics.ninvalids, 
		    "number of invalid patterns", status);
    fits_update_key(plf->eventlistfile->fptr, TLONG, "NPINVALI", 
		    &statistics.npinvalids, 
		    "number of piled up invalid patterns", status);
    // Numbered grades.
    for (ii=0; ii<13; ii++) {
      char keyword[MAXMSG];
      char comment[MAXMSG];
      sprintf(keyword, "NGRAD%ld", ii);
      sprintf(comment, "number of patterns with grade %ld", ii);
      fits_update_key(plf->eventlistfile->fptr, TLONG, keyword, 
		      &statistics.ngrade[ii], comment, status);
      sprintf(keyword, "NPGRA%ld", ii);
      sprintf(comment, "number of piled up patterns with grade %ld", ii);
      fits_update_key(plf->eventlistfile->fptr, TLONG, keyword, 
		      &statistics.npgrade[ii], comment, status);    
    }
    CHECK_STATUS_BREAK(*status);
    
  } while(0); // End of error handling loop.
  

  // Release memory.
  if (NULL!=framelist) {
    for (ii=0; ii<nframelist; ii++) {
      if (NULL!=framelist[ii]) {
	free(framelist[ii]);
      }
    }
    free(framelist);
  }
  if (NULL!=neighborlist) {
    for (ii=0; ii<nneighborlist; ii++) {
      if (NULL!=neighborlist[ii]) {
	free(neighborlist[ii]);
      }
    }
    free(neighborlist);
  }
}

