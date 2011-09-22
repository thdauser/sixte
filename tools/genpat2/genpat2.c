#include "genpat2.h"


static void copyEvent(Event* const ev1, const Event* const ev2)
{
  ev1->rawx  =ev2->rawx;
  ev1->rawy  =ev2->rawy;
  ev1->pha   =ev2->pha;
  ev1->charge=ev2->charge;
  ev1->time  =ev2->time;
  ev1->frame =ev2->frame;

  int ii;
  for(ii=0; ii<NEVENTPHOTONS; ii++) {
    ev1->ph_id[ii]  =ev2->ph_id[ii];
    ev1->src_id[ii] =ev2->src_id[ii];
  }
}


static struct PatternStatistics emptyPatternStatistics()
{
  struct PatternStatistics stat = {
    .nvalids=0,
    .npvalids=0,
    .ninvalids=0,
    .npinvalids=0
  };
  int ii;
  for (ii=0; ii<256; ii++) {
    stat.ngrade[ii]=0;
    stat.npgrade[ii]=0;
  }

  return(stat);
}


static void writePatternStatistics2FITSHeader(struct PatternStatistics stat, 
					      fitsfile* const fptr, 
					      int* const status)
{
  // Valids.
  if (fits_update_key(fptr, TLONG, "NVALID", &stat.nvalids, 
		      "number of valid patterns", status)) return;
  if (fits_update_key(fptr, TLONG, "NPVALID", &stat.npvalids, 
		      "number of piled up valid patterns", status)) return;
  // Invalids.
  if (fits_update_key(fptr, TLONG, "NINVALID", &stat.ninvalids, 
		      "number of invalid patterns", status)) return;
  if (fits_update_key(fptr, TLONG, "NPINVALI", &stat.npinvalids, 
		      "number of piled up invalid patterns", status)) return;

  // The different grades.
  int ii;
  for (ii=0; ii<256; ii++) {
    char keyword[MAXMSG];
    char comment[MAXMSG];
    sprintf(keyword, "NGRAD%d", ii);
    sprintf(comment, "number of patterns with grade %d", ii);
    if (fits_update_key(fptr, TLONG, keyword, &stat.ngrade[ii], 
			comment, status)) return;
    sprintf(keyword, "NPGRA%d", ii);
    sprintf(comment, "number of piled up patterns with grade %d", ii);
    if (fits_update_key(fptr, TLONG, keyword, &stat.npgrade[ii], 
			comment, status)) return;    
  }
}


static inline void clearGenPatPixels(GenDet* const det, 
				     Event*** const pixels) 
{
  int ii;
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    int jj;
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      freeEvent(&pixels[ii][jj]);
    }
  }
}


static void add2GenPatList(GenDet* const det, 
			   Event*** const pixels, 
			   const int x, const int y, 
			   const float split_threshold,
			   Event** const list, 
			   int* const nlist)
{
  // Check if the pixel is already contained in the list.
  int ii; 
  for (ii=0; ii<*nlist; ii++) {
    if (list[ii]==pixels[x][y]) {
      // The pixel is already contained in the list.
      return;
    }
  }

  // Add the event to the list.
  assert(NULL!=pixels[x][y]);
  list[*nlist]=pixels[x][y];
  (*nlist)++;
  assert(*nlist<10000);

  // Check the surrounding pixels.
#ifdef DIAGONAL_PATTERN_PILEUP
  // Check the diagonal pixels for pattern pile-up.
  int jj;
  int xmin = MAX(0, x-1);
  int xmax = MIN(det->pixgrid->xwidth-1, x+1);
  int ymin = MAX(0, y-1);
  int ymax = MIN(det->pixgrid->ywidth-1, y+1);
  for (ii=xmin; ii<=xmax; ii++) {
    for (jj=ymin; jj<=ymax; jj++) {
      if (NULL==pixels[ii][jj]) continue;
      if (pixels[ii][jj]->charge > split_threshold) {
	add2GenPatList(det, pixels, ii, jj, split_threshold, list, nlist);
      }
    }
  }
#else
  // Simple Pattern check: do NOT check diagonal pixels.
  int min = MAX(0, x-1);
  int max = MIN(det->pixgrid->xwidth-1, x+1);
  for (ii=min; ii<=max; ii++) {
    if (NULL==pixels[ii][y]) continue;
    if (pixels[ii][y]->charge > split_threshold) {
      add2GenPatList(det, pixels, ii, y, split_threshold, list, nlist);
    }
  }
  min = MAX(0, y-1);
  max = MIN(det->pixgrid->ywidth-1, y+1);
  for (ii=min; ii<=max; ii++) {
    if (NULL==pixels[x][ii]) continue;
    if (pixels[x][ii]->charge > split_threshold) {
      add2GenPatList(det, pixels, x, ii, split_threshold, list, nlist);
    }
  }
#endif // END of neglect diagonal pixels.
}


static void findMaxCharge(GenDet* const det,
			  Event*** const pixels,
			  int* const x, 
			  int* const y)
{
  int xn = *x;
  int yn = *y;

#ifdef DIAGONAL_PATTERN_PILEUP
  // Check the diagonal pixels for pattern pile-up.
  int ii, jj;
  int xmin = MAX(0, *x-1);
  int xmax = MIN(det->pixgrid->xwidth-1, *x+1);
  int ymin = MAX(0, *y-1);
  int ymax = MIN(det->pixgrid->ywidth-1, *y+1);
  for (ii=xmin; ii<=xmax; ii++) {
    for (jj=ymin; jj<=ymax; jj++) {
      if (NULL==pixels[ii][jj]) continue;
      if (pixels[ii][jj]->charge > pixels[xn][yn]->charge) {
	xn = ii;
	yn = jj;
      }
    }
  }
#else
  // Simple Pattern check: do NOT check diagonal pixels.
  int ii;
  int min = MAX(0, *x-1);
  int max = MIN(det->pixgrid->xwidth-1, *x+1);
  for (ii=min; ii<=max; ii++) {
    if (NULL==pixels[ii][*y]) continue;	  
    if (pixels[ii][*y]->charge > pixels[xn][yn]->charge) {
      xn = ii;
      yn = *y;
    }
  }
  min = MAX(0, *y-1);
  max = MIN(det->pixgrid->ywidth-1, *y+1);
  for (ii=min; ii<=max; ii++) {
    if (NULL==pixels[*x][ii]) continue;
    if (pixels[*x][ii]->charge > pixels[xn][yn]->charge) {
      xn = *x;
      yn = ii;
    }
  }
#endif

  // If there is a pixel in the neigborhood with a bigger charge than
  // the current maximum, perform an iterative function call.
  if ((xn!=*x) || (yn!=*y)) {
    findMaxCharge(det, pixels, &xn, &yn);
    // Return the new maximum.
    *x = xn;
    *y = yn;
  }
}


static void GenPatIdentification(GenDet* const det, 
				 Event*** const pixels, 
				 PatternFile* const file, 
				 struct PatternStatistics* const patstat,
				 int* const status)
{
  Event* list[10000];
  int nlist;

  // Loop over all pixels, searching charges/PHA values above 
  // the primary event threshold.
  int ii, jj;
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      if (NULL==pixels[ii][jj]) continue;
      if (pixels[ii][jj]->charge > det->threshold_event_lo_keV) {
	// Found an event above the primary event threshold.

	// Find the local charge maximum.
	int maxx=ii, maxy=jj;
	findMaxCharge(det, pixels, &maxx, &maxy);

	// Create a temporary event list of all pixels in the
	// neighborhood above the split threshold.
	float split_threshold;
	if (det->threshold_split_lo_fraction > 0.) {
	  split_threshold = det->threshold_split_lo_fraction*pixels[maxx][maxy]->charge;
	} else {
	  split_threshold = det->threshold_split_lo_keV;
	}
	nlist=0;
	add2GenPatList(det, pixels, maxx, maxy, split_threshold, list, &nlist);
	// Now 'list' contains all events contributing to this pattern.

	// Check if the pattern lies at the borders of the detectors.
	// In that case it is flagged as border pattern and will be
	// treated as invalid.
	int kk;
	int border=0;
	for (kk=0; kk<nlist; kk++) {
	  if ((0==list[kk]->rawx) || (det->pixgrid->xwidth-1==list[kk]->rawx) ||
	      (0==list[kk]->rawy) || (det->pixgrid->ywidth-1==list[kk]->rawy)) {
	    border=1;
	    break;
	  }
	}

	// Check if the pattern covers a larger area than a 3x3 matrix.
	int large=0;
	for (kk=1; kk<nlist; kk++) {
	  if ((abs(list[kk]->rawx-list[0]->rawx)>1) ||
	      (abs(list[kk]->rawy-list[0]->rawy)>1)) {
	    large=1;
	    break;
	  }
	}

	// Sum up the charges of all contributing events
	// above the split threshold (not only 3x3 matrix).
	float total_charge = 0.;
	for (kk=0; kk<nlist; kk++) {
	  total_charge += list[kk]->charge;
	}

	// Determine the pattern grade.
	Pattern pattern = {
	  .pat_type = 0,
	  .pileup   = 0,
	  .event    = *(pixels[maxx][maxy])
	};

	// Here we do not store the PHA values of the invididual 
	// sub-events. Therefore initialize the array with zero values.
	pattern.event.pha = getEBOUNDSChannel(total_charge, det->rmf);
	for (kk=0; kk<9; kk++) {
	  pattern.phas[kk] = 0;
	}

	// Check for pile-up:
	// - multiple photons with different IDs in the same pixel
	// - multiple photons with different IDs in neighboring pixels
	int ph_id =0;
	for(kk=0; kk<nlist; kk++) {
	  int ll;
	  for (ll=0; ll<NEVENTPHOTONS; ll++) {
	    if ((0==kk)&&(0==ll)) {
	      // No photon ID registered yet.
	      ph_id = list[0]->ph_id[0];
	      assert(ph_id!=0);
	    } else if (list[kk]->ph_id[ll]!=0) {
	      if (ph_id!=list[kk]->ph_id[ll]) {
		// Different photon ID.
		pattern.pileup=1;
		break;
	      }
	    }
	  }
	  if (1==pattern.pileup) break;
	}

	
	// Store the pattern in the output file.
#ifdef ONLY_VALID_PATTERNS
	if (pattern.pat_type != det->grading->invalid) {
#endif
	  addPattern2File(file, &pattern, status);	  
#ifdef ONLY_VALID_PATTERNS
	}
#endif
	// END of storing the pattern data in the output file.


	// Store the information about the pattern in the
	// statistics data structure.
	if (pattern.pat_type==det->grading->invalid) {
	  // Invalid pattern.
	  patstat->ninvalids++;
	  if (pattern.pileup > 0) {
	    patstat->npinvalids++;
	  }
	} else  {
	  // Valid pattern.
	  patstat->nvalids++;
	  patstat->ngrade[pattern.pat_type]++;
	  if (pattern.pileup > 0) {
	    patstat->npvalids++;
	    patstat->npgrade[pattern.pat_type]++;
	  }
	}
	// END of gathering statistical data about the pattern type.

	// Delete the events belonging to this pattern from the pixel array
	// in order to prevent them being used another time.
	for (kk=0; kk<nlist; kk++) {
	  freeEvent(&(pixels[list[kk]->rawx][list[kk]->rawy]));
	  list[kk]=NULL;
	}
	// End of deleting all contributing events.

	// TODO Due to the different option of including diagonal split events
	// some events might have been cleared now that could contribute to
	// the 3x3 matrix of a later close-by event.
      }
      // END of found an event above the primary threshold.
    }
  }
  // END of loop over all pixels searching for events above the
  // primary event threshold.
}


////////////////////////////////////
/** Main procedure. */
int genpat_main() {

  // Containing all programm parameters read by PIL
  struct Parameters par; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;
  // Input event list file.
  EventListFile* elf=NULL;
  // Output event file. 
  PatternFile* plf=NULL;
  // Detector pixel array.
  Event*** pixels=NULL;
  // Pattern statistics. Count the numbers of the individual pattern types
  // and store this information in the output event file.
  struct PatternStatistics patstat=emptyPatternStatistics();

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("genpat2");
  set_toolversion("0.02");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Read parameters using PIL library:
    if ((status=getpar(&par))) break;

    // Determine the random number seed.
    int seed;
    if (-1!=par.Seed) {
      seed = par.Seed;
    } else {
      // Determine the seed from the system clock.
      seed = (int)time(NULL);
    }

    // Initialize HEADAS random number generator.
    HDmtInit(seed);

    // Initialize the detector data structure.
    det=newGenDet(par.XMLFile, &status);
    if (EXIT_SUCCESS!=status) break;

    // Check if the detector data structure contains
    // a pattern identifier. Otherwise it is not reasonable
    // to run the pattern identification algorithm.
    if (NULL==det->grading) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: no event grading specified in detector "
		     "XML definition file!\n", status);
      break;
    }

    // Set the input event file.
    elf=openEventListFile(par.EventList, READWRITE, &status);
    if (EXIT_SUCCESS!=status) break;


    // Create and open a new event file.
    plf=openNewPatternFile(par.PatternList, &status);
    if (EXIT_SUCCESS!=status) break;

    // Copy header keywords from the input to the output event file.
    char comment[MAXMSG]; // Buffer.

    // Total number of simulated photons.
    long n_input_photons=0; 
    if (fits_read_key(elf->fptr, TLONG, "NPHOTONS", 
		      &n_input_photons, comment, &status)) break;
    if (fits_update_key(plf->eventlistfile->fptr, TLONG, "NPHOTONS", 
			&n_input_photons, "number of input photons", 
			&status)) break;

    // Total number of detected photons.
    long n_detected_photons=0; 
    if (fits_read_key(elf->fptr, TLONG, "NDETECTD", 
		      &n_detected_photons, comment, &status)) break;
    if (fits_update_key(plf->eventlistfile->fptr, TLONG, "NDETECTD", 
			&n_detected_photons, "number of detected photons", 
			&status)) break;

    // Number of pixels in x-direction.
    long nxdim=0; 
    if (fits_read_key(elf->fptr, TINT, "NXDIM", 
		      &nxdim, comment, &status)) break;
    if (fits_update_key(plf->eventlistfile->fptr, TINT, "NXDIM", 
			&nxdim, comment, &status)) break;

    // Number of pixels in y-direction.
    long nydim=0; 
    if (fits_read_key(elf->fptr, TINT, "NYDIM", 
		      &nydim, comment, &status)) break;
    if (fits_update_key(plf->eventlistfile->fptr, TINT, "NYDIM", 
			&nydim, comment, &status)) break;    
    // END of copying header keywords.


    // Allocate memory for the pixel array used for the pattern identification.
    pixels=(Event***)malloc(det->pixgrid->xwidth*sizeof(Event**));
    if (NULL==pixels) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation for pixel array failed!\n", status);
      break;
    }
    int ii;
    for (ii=0; ii<det->pixgrid->xwidth; ii++) {
      pixels[ii]=(Event**)malloc(det->pixgrid->ywidth*sizeof(Event*));
      if (NULL==pixels[ii]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: Memory allocation for pixel array failed!\n", status);
	break;
      }
      // Initialize Event pointers with NULL.
      int jj;
      for (jj=0; jj<det->pixgrid->ywidth; jj++) {
	pixels[ii][jj] = NULL;
      }
    }
    if (EXIT_SUCCESS!=status) break;

    // --- END of Initialization ---


    // --- Beginning of Pattern Identification Process ---

    headas_chat(3, "start pattern identification ...\n");

    // Loop over all events in the FITS file. The last detector 
    // frame is NOT neglected.
    long row;
    long frame=0;
    int last_loop=0;
    for (row=0; row<=elf->nrows; row++) {
      
      Event event;

      if (row<elf->nrows) {
	last_loop=0;
	// Read the next event from the file.
	getEventFromFile(elf, row+1, &event, &status);
	CHECK_STATUS_BREAK(status);
      } else {
	last_loop = 1;
      }

      // If the event belongs to a new frame, perform the
      // Pattern Identification on the pixel array before
      // starting a new frame.
      if ((event.frame > frame) || (1==last_loop)) {

	// Run the pattern identification and store the pattern 
	// information in the event file.
	GenPatIdentification(det, pixels, plf, &patstat, &status);
	CHECK_STATUS_BREAK(status);
	
	// Delete the old events in the pixel array.
	clearGenPatPixels(det, pixels);

	// Update the frame counter.
	frame = event.frame;
	if (0==frame%100) {
	  headas_printf("\rframe: %ld ", frame);
	  fflush(NULL);
	}
      }

      if (0==last_loop) {
	// Add the event to the pixel array.
	pixels[event.rawx][event.rawy] = getEvent(&status);
	CHECK_STATUS_BREAK(status);
	copyEvent(pixels[event.rawx][event.rawy], &event);
      }
    }
    CHECK_STATUS_BREAK(status);
    headas_printf("\n");
    // END of loop over all events in the FITS file.

    // Store the pattern statistics in the FITS header.
    writePatternStatistics2FITSHeader(patstat, plf->eventlistfile->fptr, &status);
    CHECK_STATUS_BREAK(status);

  } while(0); // END of the error handling loop.

  // --- END of pattern identification process ---


  // --- Cleaning up ---
  headas_chat(3, "cleaning up ...\n");

  // Release HEADAS random number generator.
  HDmtFree();

  // Release memory from pixel array.
  if (NULL!=pixels) {
    int ii;
    for (ii=0; ii<det->pixgrid->xwidth; ii++) {
      if (NULL!=pixels[ii]) {
	free(pixels[ii]);
      }
    }
    free(pixels);
    pixels=NULL;
  }

  // Destroy the detector data structure.
  destroyGenDet(&det, &status);
  
  // Close the files.
  freeEventListFile(&elf, &status);
  destroyPatternFile(&plf, &status);
  
  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}


int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_file_name("XMLFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the XML file!\n", status);
    return(status);
  } 
  strcpy(par->XMLFile, sbuffer);
  free(sbuffer);


  status=ape_trad_query_file_name("EventList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the event list!\n", status);
    return(status);
  } 
  strcpy(par->EventList, sbuffer);
  free(sbuffer);


  status=ape_trad_query_string("PatternList", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the name of the pattern list!\n", status);
    return(status);
  } 
  strcpy(par->PatternList, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("seed", &par->Seed);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the seed for the random number generator!\n", status);
    return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    HD_ERROR_THROW("Error reading the clobber parameter!\n", status);
    return(status);
  }

  return(status);
}


