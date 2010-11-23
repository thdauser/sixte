#include "genpat.h"


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


static inline GenEvent emptyEvent() 
{
  GenEvent empty_event = {.frame=0};
  assert(0==empty_event.rawx);
  assert(0==empty_event.rawy);
  assert(0==empty_event.pileup);
  assert(0==empty_event.pha);
  assert(0.==empty_event.charge);
  assert(0==empty_event.frame);

  return(empty_event);
}


static inline void clearGenPatPixels(GenDet* const det, 
				     GenEvent** const pixels) 
{
  int ii;
#pragma omp parallel for
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    int jj;
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      pixels[ii][jj] = emptyEvent();
    }
  }
}


			     
static void add2GenPatList(GenDet* const det, 
			   GenEvent** const pixels, 
			   const int x, const int y, 
			   const float split_threshold,
			   GenEvent** const list, 
			   int* const nlist)
{
  // Check if the pixel is already contained in the list.
  int ii; 
  for (ii=0; ii<*nlist; ii++) {
    if (list[ii]==&pixels[x][y]) {
      // The pixel is already contained in the list.
      return;
    }
  }

  // Add the event to the list.
  list[*nlist] = &pixels[x][y];
  (*nlist)++;
  assert(*nlist<1000);

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
      if (pixels[ii][jj].charge > split_threshold) {
	add2GenPatList(det, pixels, ii, jj, split_threshold, list, nlist);
      }
    }
  }
#else
  // Simple Pattern check: do NOT check diagonal pixels.
  int min = MAX(0, x-1);
  int max = MIN(det->pixgrid->xwidth-1, x+1);
  for (ii=min; ii<=max; ii++) {
    if (pixels[ii][y].charge > split_threshold) {
      add2GenPatList(det, pixels, ii, y, split_threshold, list, nlist);
    }
  }
  min = MAX(0, y-1);
  max = MIN(det->pixgrid->ywidth-1, y+1);
  for (ii=min; ii<=max; ii++) {
    if (pixels[x][ii].charge > split_threshold) {
      add2GenPatList(det, pixels, x, ii, split_threshold, list, nlist);
    }
  }
#endif // END of neglect diagonal pixels.
}



static void findMaxCharge(GenDet* const det,
			  GenEvent** const pixels,
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
      if (pixels[ii][jj].charge > pixels[xn][yn].charge) {
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
    if (pixels[ii][*y].charge > pixels[xn][yn].charge) {
      xn = ii;
      yn = *y;
    }
  }
  min = MAX(0, *y-1);
  max = MIN(det->pixgrid->ywidth-1, *y+1);
  for (ii=min; ii<=max; ii++) {
    if (pixels[*x][ii].charge > pixeks[xn][yn].charge) {
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
				 GenEvent** const pixels, 
				 GenPatternFile* const file, 
				 struct PatternStatistics* const patstat,
				 int* const status)
{
  GenEvent* list[1000];
  int nlist;

  // Loop over all pixels, searching charges/PHA values above 
  // the primary event threshold.
  int ii, jj;
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      if (pixels[ii][jj].charge > det->threshold_event_lo_keV) {
	// Found an event above the primary event threshold.

	// Find the local charge maximum.
	int maxx=ii, maxy=jj;
	findMaxCharge(det, pixels, &maxx, &maxy);
	
	// Create a temporary event list of all pixels in the
	// neighborhood above the split threshold.
	float split_threshold;
	if (det->threshold_split_lo_fraction > 0) {
	  split_threshold = det->threshold_split_lo_fraction*pixels[maxx][maxy].charge;
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
	    large = 1;
	    break;
	  }
	}

	// Determine the charges of the contributing events
	// above the split threshold in the 3x3 matrix
	// around the main event with the maximum charge.
	float charges[9] = {0., 0., 0.,  0., 0., 0.,  0., 0., 0.};
	charges[4] = pixels[maxx][maxy].charge; // Main pixel.
	for (kk=1; kk<nlist; kk++) {
	  if (list[kk]->rawy == maxy-1) {
	    if (list[kk]->rawx == maxx-1) {
	      charges[0] = list[kk]->charge;
	    } else if (list[kk]->rawx == maxx) {
	      charges[1] = list[kk]->charge;
	    } else if (list[kk]->rawx == maxx+1) {
	      charges[2] = list[kk]->charge;
	    }
	  } else if (list[kk]->rawy == maxy) {
	    if (list[kk]->rawx == maxx-1) {
	      charges[3] = list[kk]->charge;
	    } else if (list[kk]->rawx == maxx+1) {
	      charges[5] = list[kk]->charge;
	    }
	  } else if (list[kk]->rawy == maxy+1) {
	    if (list[kk]->rawx == maxx-1) {
	      charges[6] = list[kk]->charge;
	    } else if (list[kk]->rawx == maxx) {
	      charges[7] = list[kk]->charge;
	    } else if (list[kk]->rawx == maxx+1) {
	      charges[8] = list[kk]->charge;
	    }
	  }
	}
	// END of determine the charge distribution in the 3x3 matrix
	// around the main event.

	// Determine the total charge.
	float total_charge = 0.;
	for (kk=0; kk<9; kk++) {
	  total_charge += charges[kk];
	}

	// Determine the pattern grade.
	GenPattern pattern = {
	  .pat_type = getGenEventGrade(det->grading, charges, 
				       border, large),
	  .event = pixels[maxx][maxy]
	};

	// Combine the PHA values of the individual events.
	pattern.event.pha = getEBOUNDSChannel(total_charge, det->rmf);

	// Store the PHA values of the pixels above the threshold in  the 
	// 3x3 matrix around the central pixel.
	for (kk=0; kk<9; kk++) {
	  pattern.phas[kk] = 0;
	}
	for (kk=0; kk<nlist; kk++) {
	  if ((abs(list[kk]->rawx-maxx)<2) && (abs(list[kk]->rawy-maxy)<2)) {
	    pattern.phas[(list[kk]->rawx-maxx+1) + (3*(list[kk]->rawy-maxy+1))] = 
	      list[kk]->pha;
	  }
	}

	
	// Store the pattern in the output file.
#ifdef ONLY_VALID_PATTERNS
	if (pattern.pat_type != det->grading->invalid) {
#endif
	  addGenPattern2File(file, &pattern, status);	  
#ifdef ONLY_VALID_PATTERNS
	}
#endif
	// END of storing the pattern data in the output file.


	// Store the information about the pattern in the
	// statistics data structure.
	if (pattern.pat_type==det->grading->invalid) {
	  // Invalid pattern.
	  patstat->ninvalids++;
	  if (pattern.event.pileup > 0) {
	    patstat->npinvalids++;
	  }
	} else  {
	  // Valid pattern.
	  patstat->nvalids++;
	  patstat->ngrade[pattern.pat_type]++;
	  if (pattern.event.pileup > 0) {
	    patstat->npvalids++;
	    patstat->npgrade[pattern.pat_type]++;
	  }
	}
	// END of gathering statistical data about the pattern type.

	// Delete the events belonging to the pattern from the pixel array
	// in order to prevent them being used another time. Therefore first 
	// create a new list with contributing events, also the ones below 
	// the original split threshold.
	// Otherwise there might be some events left below the split threshold.
	nlist=0;
	add2GenPatList(det, pixels, maxx, maxy, 0., list, &nlist);
	for (kk=0; kk<nlist; kk++) {
	  *(list[kk]) = emptyEvent();
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
  struct Parameters parameters; 
  // Detector data structure (containing the pixel array, its width, ...).
  GenDet* det=NULL;
  // Output event file. 
  GenPatternFile* output_file=NULL;
  // Detector pixel array.
  GenEvent** pixels=NULL;
  // Pattern statistics. Count the numbers of the individual pattern types
  // and store this information in the output event file.
  struct PatternStatistics patstat=emptyPatternStatistics();

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

    // Check if the detector data structure contains
    // a pattern identifier. Otherwise it is not reasonable
    // to run the pattern identification algorithm.
    if (NULL==det->grading) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: no event grading specified in detector XML definition file!\n",
		     status);
      break;
    }

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
    // Open a new pattern file from the specified template.
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
	GenPatIdentification(det, pixels, output_file, &patstat, &status);
	
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
    writePatternStatistics2FITSHeader(patstat, output_file->geneventfile->fptr, &status);
    if (EXIT_SUCCESS!=status) break;

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


