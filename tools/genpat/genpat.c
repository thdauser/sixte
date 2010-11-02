#include "genpat.h"


static inline void clearGenPatPixels(GenDet* const det, 
				     GenEvent** const pixels) 
{
  GenEvent empty_event = {.frame=0};
  assert(0==empty_event.rawx);
  assert(0==empty_event.rawy);
  assert(0==empty_event.pileup);
  assert(0==empty_event.pha);
  assert(0.==empty_event.charge);
  assert(0==empty_event.frame);

  int ii, jj;
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      pixels[ii][jj] = empty_event;
    }
  }
}


			     
static void addGenPat2List(GenDet* const det, 
			   GenEvent** const pixels, 
			   const int x, const int y, 
			   GenEvent* const list, 
			   int* const nlist)
{
  // Add the event to the list.
  list[*nlist] = pixels[x][y];
  (*nlist)++;
  assert(*nlist<1000);
  // Delete the pixel, such that it is not used twice.
  pixels[x][y].charge = 0.;

  // Check the surrounding pixels.
#ifdef DIAGONAL_PATTERN_PILEUP
  // Check the diagonal pixels for pattern pile-up.
  int ii, jj;
  int xmin = MAX(0, x-1);
  int xmax = MIN(det->pixgrid->xwidth-1, x+1);
  int ymin = MAX(0, y-1);
  int ymax = MIN(det->pixgrid->ywidth-1, y+1);
  for (ii=xmin; ii<=xmax; ii++) {
    for (jj=ymin; jj<=ymax; jj++) {
      if ((pixels[ii][jj].charge > list[0].charge*det->threshold_split_lo_fraction) &&
	  (pixels[ii][jj].charge > det->threshold_split_lo_keV)) {
	addGenPat2List(det, pixels, ii, jj, list, nlist);
      }
    }
  }
  // END of loop over surrounding pixels.
#else
  // Simple Pattern check: do NOT check diagonal pixels.
  int ii;
  int min = MAX(0, x-1);
  int max = MIN(det->pixgrid->xwidth-1, x+1);
  for (ii=min; ii<=max; ii++) {
    if ((pixels[ii][y].charge > list[0].charge*det->threshold_split_lo_fraction) &&
	(pixels[ii][y].charge > det->threshold_split_lo_keV)) {
      addGenPat2List(det, pixels, ii, y, list, nlist);
    }
  }
  min = MAX(0, y-1);
  max = MIN(det->pixgrid->ywidth-1, y+1);
  for (ii=min; ii<=max; ii++) {
    if ((pixels[x][ii].charge > list[0].charge*det->threshold_split_lo_fraction) &&
	(pixels[x][ii].charge > det->threshold_split_lo_keV)) {
      addGenPat2List(det, pixels, x, ii, list, nlist);
    }
  }
#endif // END of neglect diagonal pixels.
}



static void GenPatId(GenDet* const det, 
		     GenEvent** const pixels, 
		     GenPatternFile* const file, 
		     struct PatternStatistics* const patstat,
		     int* const status)
{
  GenEvent list[1000];
  int nlist;

  // Loop over all pixels, searching charges/PHA values above 
  // the primary event threshold.
  int ii, jj;
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      if (pixels[ii][jj].charge > det->threshold_event_lo_keV) {
	// Found an event above the primary event threshold.
	
	// Add the event to the temporary event list and
	// check the surrounding pixels.
	// The function addGetPat2List deletes the charge in
	// the processed pixel such that it is not taken into
	// account later again.
	nlist=0;
	addGenPat2List(det, pixels, ii, jj, list, &nlist);
	// Now 'list' contains all events contributing to this pattern.

	// Check if the pattern lies at the borders of the detectors.
	// In that case it is flagged as border pattern and will be
	// treated as invalid.
	int border=0;
	int kk;
	for (kk=0; kk<nlist; kk++) {
	  if ((0==list[kk].rawx) || (det->pixgrid->xwidth-1==list[kk].rawx) ||
	      (0==list[kk].rawy) || (det->pixgrid->ywidth-1==list[kk].rawy)) {
	    border=1;
	  }
	}

	// Determine the pattern type and orientation.
	int ll;
	// Indices of maximum charged pixels in descending order.
	int idx[1000] = {};
	assert(idx[0]==0);
	assert(idx[999]==0);
	// Determine indices with maximum event charge at the
	// first position and the weakest partner at the end: idx[0]..idx[3].
	for (kk=0; kk<nlist; kk++) {
	  idx[kk]=kk;
	  for (ll=kk-1; ll>=0; ll--) {
	    if (list[idx[ll]].charge<list[kk].charge) {
	      idx[ll+1] = idx[ll];
	      idx[ll]   = kk;
	    } else {
	      break;
	    }
	  }
	}

	// Determine the pattern type and alignment (orientation).
	int pat_type=0;
	switch (nlist) {
	case 1: // Single event.
	  pat_type = 1;
	  break;
	  
	case 2: // Double event.
	  if (list[idx[1]].rawx==list[idx[0]].rawx) {
	    if ((list[idx[1]].rawy==list[idx[0]].rawy+1) ||
		(list[idx[1]].rawy==list[idx[0]].rawy-1)) {
	      pat_type = 2;
	    } else {
	      pat_type = -1;
	    }
	  } else if (list[idx[1]].rawy==list[idx[0]].rawy) {
	    if ((list[idx[1]].rawx==list[idx[0]].rawx+1) ||
		(list[idx[1]].rawx==list[idx[0]].rawx-1)) {
	      pat_type = 2;
	    } else {
	      pat_type = -1;
	    }
	  } else {
	    pat_type = -1;
	  }
	  break;
	  
	case 3: // Triple event.
	  if (list[idx[1]].rawx==list[idx[0]].rawx) {
	    if (((list[idx[1]].rawy==list[idx[0]].rawy+1) ||
		 (list[idx[1]].rawy==list[idx[0]].rawy-1)) &&

		((list[idx[2]].rawy==list[idx[0]].rawy) &&
		 ((list[idx[2]].rawx==list[idx[0]].rawx+1) ||
		  (list[idx[2]].rawx==list[idx[0]].rawx-1)))) {
		  pat_type = 3;
	    } else {
	      pat_type = -1;
	    }
	  } // END of RAWX(maximum) == RAWX(second largest event).

	  else if (list[idx[1]].rawy==list[idx[0]].rawy) {
	    if (((list[idx[1]].rawx==list[idx[0]].rawx+1) ||
		 (list[idx[1]].rawx==list[idx[0]].rawx-1)) &&

		((list[idx[2]].rawx==list[idx[0]].rawx) &&
		 ((list[idx[2]].rawy==list[idx[0]].rawy+1) ||
		  (list[idx[2]].rawy==list[idx[0]].rawy-1)))) {
		  pat_type = 3;
	    } else {
	      pat_type = -1;
	    }
	  } // END of RAWY(maximum) == RAWY(second largest event).
	  else {
	    pat_type = -1;
	  }
	  break;

	case 4: // Quadruple event.
	  
	  if ((((list[idx[1]].rawy == list[idx[3]].rawy) &&
		(list[idx[1]].rawx == list[idx[0]].rawx)) ||
	       ((list[idx[1]].rawy == list[idx[0]].rawy) &&
		(list[idx[1]].rawx == list[idx[3]].rawx))) && 
	      
	      (((list[idx[2]].rawy == list[idx[3]].rawy) &&
		(list[idx[2]].rawx == list[idx[0]].rawx)) ||
	       ((list[idx[2]].rawy == list[idx[0]].rawy) &&
		(list[idx[2]].rawx == list[idx[3]].rawx)))) {
	    pat_type = 4;
	  } else {
	    pat_type = -1;
	  }
	  break;
	  
	default:
	  pat_type=-1;
	  break;
	}
	// END of switch(nlist).

	// Store the pattern in the output file.
	if ((pat_type > 0) && (0==border)) { 
	  // This is a valid pattern: combine the individual pixel
	  // events and store only one pattern.
	  GenPattern pattern = {.pat_type=0};

	  // Set the event data from the central main event.
	  pattern.event    = list[0];
	  pattern.pat_type = pat_type;

	  // Recombine the PHA values of the individual events.
	  float charge = 0.;
	  for (kk=0; kk<nlist; kk++) {
	    charge += list[kk].charge;
	  }	  
	  pattern.event.pha = getEBOUNDSChannel(charge, det->rmf);

	  // Store the PHA values of the pixels above the threshold in  the 
	  // 3x3 matrix around the central pixel.
	  for (kk=0; kk<9; kk++) {
	    pattern.phas[kk] = 0;
	  }
	  for (kk=0; kk<nlist; kk++) {
	    pattern.phas[(list[kk].rawx-list[0].rawx+1) +
			 (3*(list[kk].rawy-list[0].rawy+1))] = list[kk].pha;
	  }

	  // Store the pattern in the output file.
	  addGenPattern2File(file, &pattern, status);	  

	} else if ((-1==pat_type) || (1==border)) {
#ifndef ONLY_VALID_PATTERNS

	  // This is an INvalid pattern: store each individual event as 
	  // a single pattern with patternID = -1.

	  // Transfer PHA of the current pixel.
	  for (kk=0; kk<nlist; kk++) {
	    GenPattern pattern = {.pat_type=0};
	  
	    pattern.event    = list[kk];
	    pattern.phas[4]  = list[kk].pha;
	    pattern.pat_type = -1;

	    // Store the pattern in the output file.
	    addGenPattern2File(file, &pattern, status);	  
	  }
#endif
	} else {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: Could not determine pattern type!\n", *status);
	  return;
	}
	// END of adding the pattern data to the output file.


	// Store the information about the pattern type in the
	// statistics data structure.
	switch(pat_type) {
	case -1:
	  patstat->ninvalids++;
	  break;
	case 1:
	  patstat->nsingles++;
	  break;
	case 2:
	  patstat->ndoubles++;
	  break;
	case 3:
	  patstat->ntriples++;
	  break;
	case 4:
	  patstat->nquadruples++;
	  break;
	default:
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: Invalid Pattern Type!\n", *status);
	  return;
	}

	// Check if it's a valid pattern.
	if (pat_type > 0) {
	  patstat->nvalids++;
	}

	// Check for pile-up.
	if (list[0].pileup > 0) {
	  patstat->npileup++;
	  if (pat_type>0) {
	    patstat->npileup_valid++;
	  } else {
	    patstat->npileup_invalid++;
	  }
	  if (1==pat_type) {
	    patstat->npileup_singles++;
	  }
	}
	// END of gathering statistical data about the pattern type.
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
  struct Parameters parameters; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;
  // Output event file. 
  GenPatternFile* output_file=NULL;
  // Detector pixel array.
  GenEvent** pixels=NULL;
  // Pattern statistics. Count the numbers of the individual pattern types
  // and store this information in the output event file.
  struct PatternStatistics patstat={
    .nsingles=0,
    .ndoubles=0,
    .ntriples=0,
    .nquadruples=0,
    .nvalids=0,
    .ninvalids=0,
    .npileup=0,
    .npileup_singles=0,
    .npileup_valid=0,
    .npileup_invalid=0
  };

  int status=EXIT_SUCCESS; // Error status.


  // Register HEATOOL:
  set_toolname("genpat");
  set_toolversion("0.01");


  do { // Beginning of the ERROR handling loop (will at most be run once).

    // --- Initialization ---

    headas_chat(3, "initialization ...\n");

    // Initialize HEADAS random number generator.
    HDmtInit(SIXT_HD_RANDOM_SEED);

    // Read parameters using PIL library:
    if ((status=getpar(&parameters))) break;

    // Initialize the detector data structure.
    det = newGenDet(parameters.xml_filename, &status);
    if (EXIT_SUCCESS!=status) break;

    // Set the input event file.
    det->eventfile=openGenEventFile(parameters.eventlist_filename, 
				    READWRITE, &status);
    if (EXIT_SUCCESS!=status) break;


    // Create and open the output event file.
    // Filename of the template file.
    char template[MAXMSG];
    // Get the name of the FITS template directory.
    // Try to read it from the environment variable.
    char* buffer;
    if (NULL!=(buffer=getenv("SIXT_FITS_TEMPLATES"))) {
      strcpy(template, buffer);
    } else {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Could not read environment variable 'SIXT_FITS_TEMPLATES'!\n", 
		     status);
      break;
    }
    // Append the filename of the template file itself.
    strcat(template, "/");
    strcat(template, det->patternfile_template);
    // Open a new event file from the specified template.
    output_file = openNewGenPatternFile(parameters.patternlist_filename, template, &status);
    if (EXIT_SUCCESS!=status) break;

    // Copy header keywords from the input to the output event file.
    char comment[MAXMSG]; // Buffer.

    // Total number of detected photons.
    long n_detected_photons=0; 
    if (fits_read_key(det->eventfile->fptr, TLONG, "NDETECTD", 
		      &n_detected_photons, comment, &status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NDETECTD", 
			&n_detected_photons, "number of detected photons", 
			&status)) break;

    // Number of EBOUNDS channels (DETCHANS).
    long detchans=0; 
    if (fits_read_key(det->eventfile->fptr, TLONG, "DETCHANS", 
		      &detchans, comment, &status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "DETCHANS", 
			&detchans, comment, &status)) break;

    // First EBOUNDS channel.
    long tlmin1=0; 
    if (fits_read_key(det->eventfile->fptr, TLONG, "TLMIN1", 
		      &tlmin1, comment, &status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "TLMIN1", 
			&tlmin1, comment, &status)) break;    

    // Last EBOUNDS channel.
    long tlmax1=0; 
    if (fits_read_key(det->eventfile->fptr, TLONG, "TLMAX1", 
		      &tlmax1, comment, &status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "TLMAX1", 
			&tlmax1, comment, &status)) break;    

    // Number of pixels in x-direction.
    long nxdim=0; 
    if (fits_read_key(det->eventfile->fptr, TINT, "NXDIM", 
		      &nxdim, comment, &status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TINT, "NXDIM", 
			&nxdim, comment, &status)) break;

    // Number of pixels in y-direction.
    long nydim=0; 
    if (fits_read_key(det->eventfile->fptr, TINT, "NYDIM", 
		      &nydim, comment, &status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TINT, "NYDIM", 
			&nydim, comment, &status)) break;    
    // END of copying header keywords.


    // Allocate memory for the pixel array used for the pattern identification.
    pixels=(GenEvent**)malloc(det->pixgrid->xwidth*sizeof(GenEvent*));
    if (NULL==pixels) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: Memory allocation for pixel array failed!\n", status);
      break;
    }
    int ii;
    for (ii=0; ii<det->pixgrid->xwidth; ii++) {
      pixels[ii]=(GenEvent*)malloc(det->pixgrid->ywidth*sizeof(GenEvent));
      if (NULL==pixels[ii]) {
	status = EXIT_FAILURE;
	HD_ERROR_THROW("Error: Memory allocation for pixel array failed!\n", status);
	break;
      }
      // Initialize the event data structures with 0 values.
      int jj;
      for (jj=0; jj<det->pixgrid->ywidth; jj++) {
	pixels[ii][jj].time = 0.;
	pixels[ii][jj].pha  = 0;
	pixels[ii][jj].charge = 0.;
	pixels[ii][jj].rawx = 0;
	pixels[ii][jj].rawy = 0;
	pixels[ii][jj].frame  = 0;
	pixels[ii][jj].pileup = 0;
      }
    }
    if (EXIT_SUCCESS!=status) break;

    // --- END of Initialization ---


    // --- Beginning of Pattern Identification Process ---

    headas_chat(3, "start pattern identification ...\n");

    // Loop over all events in the FITS file. The last detector 
    // frame is NOT neglected.
    GenEvent event;
    long row;
    long frame=0;
    int last_loop=0;
    for (row=0; row<=det->eventfile->nrows; row++) {

      if (row<det->eventfile->nrows) {
	last_loop=0;
	// Read the next event from the file.
	getGenEventFromFile(det->eventfile, row+1, &event, &status);
	if (EXIT_SUCCESS!=status) break;
      } else {
	last_loop = 1;
      }

      // If the event belongs to a new frame, perform the
      // Pattern Identification on the pixel array before
      // starting a new frame.
      if ((event.frame > frame) || (1==last_loop)) {

	// Run the pattern identification and store the pattern 
	// information in the event file.
	GenPatId(det, pixels, output_file, &patstat, &status);
	
	// Delete the old events in the pixel array.
	clearGenPatPixels(det, pixels);

	// Update the frame counter.
	frame = event.frame;
	headas_printf("\rframe: %ld ", frame);
	fflush(NULL);
      }

      if (0==last_loop) {
	// Add the event to the pixel array.
	pixels[event.rawx][event.rawy] = event;
      }
    };
    if (EXIT_SUCCESS!=status) break;
    // END of loop over all events in the FITS file.
    headas_printf("\n");

    // Store the pattern statistics in the FITS header.
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NSINGLES", 
			&patstat.nsingles, "number of single patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NDOUBLES", 
			&patstat.ndoubles, "number of double patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NTRIPLES", 
			&patstat.ntriples, "number of triple patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NQUADRUP", 
			&patstat.nquadruples, "number of quadruple patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NVALIDS", 
			&patstat.nvalids, "number of valid patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NINVALID", 
			&patstat.ninvalids, "number of invalid patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NPILEUP", 
			&patstat.npileup, "number of pile-up patterns", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NPILEUPS", 
			&patstat.npileup_singles, 
			"number of singles marked as pile-up", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NPILEUPV", 
			&patstat.npileup_valid, 
			"number of valid patterns marked as pile-up", 
			&status)) break;
    if (fits_update_key(output_file->geneventfile->fptr, TLONG, "NPILEUPI", 
			&patstat.npileup_invalid, 
			"number of invalid patterns marked as pile-up", 
			&status)) break;

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
  
  // Close the output eventfile.
  destroyGenPatternFile(&output_file, &status);
  
  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* const parameters)
{
  int status=EXIT_SUCCESS; // Error status

  // Get the name of the input event list file (FITS file).
  if ((status = PILGetFname("eventlist_filename", 
			    parameters->eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input event list file!\n", status);
  }

  // Get the name of the output event list file (FITS file).
  else if ((status = PILGetFname("patternlist_filename", 
				 parameters->patternlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output pattern list file!\n", status);
  }

  // Get the name of the detector XML description file (FITS file).
  else if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector definition XML file!\n", status);
  }

  return(status);
}


