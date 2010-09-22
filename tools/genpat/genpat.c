#if HAVE_CONFIG_H
#include <config.h>
#else
#error "Do not compile outside Autotools!"
#endif

#include "genpat.h"


static inline void clearGenPatPixels(GenDet* const det, 
				     GenEvent** const pixels) 
{
  GenEvent empty_event = { .time=0. };
  assert(0==empty_event.rawx);
  assert(0==empty_event.rawy);
  assert(0==empty_event.pileup);
  assert(0==empty_event.pha);
  assert(0.==empty_event.charge);
  assert(0==empty_event.frame);
  assert(0==empty_event.pat_type);
  assert(0==empty_event.pat_id);
  assert(0==empty_event.pat_alig);

  int ii, jj;
  for (ii=0; ii<det->pixgrid->xwidth; ii++) {
    for (jj=0; jj<det->pixgrid->ywidth; jj++) {
      pixels[ii][jj] = empty_event;
    }
  }
}


			     
static void addGenPat2List(GenDet* const det, GenEvent** const pixels, 
			   const int x, const int y, 
			   GenEvent* const list, int* const nlist)
{
  // Add the event to the list.
  list[*nlist] = pixels[x][y];
  (*nlist)++;
  assert(*nlist<1000);
  // Delete the pixel, such that it is not used twice.
  pixels[x][y].charge = 0.;

  // Check the surrounding pixels.
  int ii, jj;
  int xmin = MAX(0, x-1);
  int xmax = MIN(det->pixgrid->xwidth-1, x+1);
  int ymin = MAX(0, y-1);
  int ymax = MIN(det->pixgrid->ywidth-1, y+1);
  for (ii=xmin; ii<=xmax; ii++) {
    for (jj=ymin; jj<=ymax; jj++) {
      if (pixels[ii][jj].charge > list[0].charge*det->threshold_split_lo_fraction) {
	addGenPat2List(det, pixels, ii, jj, list, nlist);
      }
    }
  }
  // END of loop over surrounding pixels.
}



static void GenPatId(GenDet* const det, GenEvent** const pixels, 
		     GenEventFile* const file, int* const status)
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
	// Check the surrounding pixels.
	nlist=0;
	addGenPat2List(det, pixels, ii, jj, list, &nlist);
	// Now 'list' contains all events contributing to this pattern.

	// Determine the pattern type and orientation.
	int kk, ll;
	// Indices of maximum charged pixels in descending order.
	int idx[4] = { 0, 0, 0, 0 };
	switch (nlist) {
	case 1:  // Single event.
	  list[0].pat_type = 1;
	  list[0].pat_id   = 5;
	  list[0].pat_alig = 0;
	  break;

	case 2: // Double, Triple, or Quadruple.
	case 3:
	case 4: 

	  // Determine maxidx, minidx; or better: idx[0]..idx[3].
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

	  // The maximum charge in the central pixel.
	  list[idx[0]].pat_id = 5; 

	  // Determine the pattern type and alignment (orientation).
	  int pat_type=0, pat_alig=0;
	  switch (nlist) {
	  case 2: // Doubles.
	    if (list[idx[1]].rawx==list[idx[0]].rawx) {
	      if (list[idx[1]].rawy==list[idx[0]].rawy+1) {
		pat_type = 2;
		list[idx[1]].pat_id = 2;
		pat_alig = 1;
	      } else if (list[idx[1]].rawy==list[idx[0]].rawy-1) {
		pat_type = 2;
		list[idx[1]].pat_id = 8;
		pat_alig = 5;
	      } else {
		pat_type = -1;
	      }
	    } else if (list[idx[1]].rawy==list[idx[0]].rawy) {
	      if (list[idx[1]].rawx==list[idx[0]].rawx+1) {
		pat_type = 2;
		list[idx[1]].pat_id = 6;
		pat_alig = 3;
	      } else if (list[idx[1]].rawx==list[idx[0]].rawx-1) {
		pat_type = 2;
		list[idx[1]].pat_id = 4;
		pat_alig = 7;
	      } else {
		pat_type = -1;
	      }
	    } else {
	      pat_type = -1;
	    }
	    break;

	  case 3: // Triples.
	    if (list[idx[1]].rawx==list[idx[0]].rawx) {
	      if (list[idx[1]].rawy==list[idx[0]].rawy+1) {
		list[idx[1]].pat_id = 2;
		
		if (list[idx[2]].rawy==list[idx[0]].rawy) {
		  if (list[idx[2]].rawx==list[idx[0]].rawx+1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 6;
		    pat_alig            = 1;
		  } else if (list[idx[2]].rawx==list[idx[0]].rawx-1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 4;
		    pat_alig            = 8;
		  } else {
		    pat_type = -1;
		  }
		} else {
		  pat_type = -1;
		}

	      } else if (list[idx[1]].rawy==list[idx[0]].rawy-1) {
		list[idx[1]].pat_id = 8;

		if (list[idx[2]].rawy==list[idx[0]].rawy) {
		  if (list[idx[2]].rawx==list[idx[0]].rawx+1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 6;
		    pat_alig            = 4;
		  } else if (list[idx[2]].rawx==list[idx[0]].rawx-1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 4;
		    pat_alig            = 5;
		  } else {
		    pat_type = -1;
		  }
		} else {
		  pat_type = -1;
		}
	      } else {
		pat_type = -1;
	      }
	    } // END of RAWX(maximum) == RAWX(second largest event).

	    else if (list[idx[1]].rawy==list[idx[0]].rawy) {
	      if (list[idx[1]].rawx==list[idx[0]].rawx+1) {
		list[idx[1]].pat_id = 6;
		
		if (list[idx[2]].rawx==list[idx[0]].rawx) {
		  if (list[idx[2]].rawy==list[idx[0]].rawy+1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 2;
		    pat_alig            = 2;
		  } else if (list[idx[2]].rawy==list[idx[0]].rawy-1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 8;
		    pat_alig            = 3;
		  } else {
		    pat_type = -1;
		  }
		} else {
		  pat_type = -1;
		}

	      } else if (list[idx[1]].rawx==list[idx[0]].rawx-1) {
		list[idx[1]].pat_id = 4;

		if (list[idx[2]].rawx==list[idx[0]].rawx) {
		  if (list[idx[2]].rawy==list[idx[0]].rawy+1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 2;
		    pat_alig            = 7;
		  } else if (list[idx[2]].rawy==list[idx[0]].rawy-1) {
		    pat_type            = 3;
		    list[idx[2]].pat_id = 8;
		    pat_alig            = 6;
		  } else {
		    pat_type = -1;
		  }
		} else {
		  pat_type = -1;
		}
	      } else {
		pat_type = -1;
	      }
	    } // END of RAWY(maximum) == RAWY(second largest event).
	    else {
	      pat_type = -1;
	    }
	    break;

	  case 4: // Quadruples.

	    if ((((list[idx[1]].rawy == list[idx[3]].rawy) &&
		  (list[idx[1]].rawx == list[idx[0]].rawx)) ||
		 ((list[idx[1]].rawy == list[idx[0]].rawy) &&
		  (list[idx[1]].rawx == list[idx[3]].rawx))) && 
		
		(((list[idx[2]].rawy == list[idx[3]].rawy) &&
		  (list[idx[2]].rawx == list[idx[0]].rawx)) ||
		 ((list[idx[2]].rawy == list[idx[0]].rawy) &&
		  (list[idx[2]].rawx == list[idx[3]].rawx)))) {
	      pat_type = 4;

	      // Check the location of the minimum event.
	      if (list[idx[3]].rawx==list[idx[0]].rawx+1) {
		if (list[idx[3]].rawy==list[idx[0]].rawy+1) {
		  list[idx[3]].pat_id = 3;
		  if (list[idx[1]].rawx==list[idx[0]].rawx) {
		    list[idx[1]].pat_id = 2;
		    list[idx[2]].pat_id = 6;
		    pat_alig            = 1;
		  } else {
		    list[idx[1]].pat_id = 6;
		    list[idx[2]].pat_id = 2;
		    pat_alig            = 2;
		  }
		} else if (list[idx[3]].rawy==list[idx[0]].rawy-1) {
		  list[idx[3]].pat_id = 9;
		  if (list[idx[1]].rawx==list[idx[0]].rawx) {
		    list[idx[1]].pat_id = 8;
		    list[idx[2]].pat_id = 6;
		    pat_alig            = 4;
		  } else {
		    list[idx[1]].pat_id = 6;
		    list[idx[2]].pat_id = 8;
		    pat_alig            = 3;
		  }
		}
	      } else if (list[idx[3]].rawx==list[idx[0]].rawx-1) {
		if (list[idx[3]].rawy==list[idx[0]].rawy+1) {
		  list[idx[3]].pat_id = 1;
		  if (list[idx[1]].rawx==list[idx[0]].rawx) {
		    list[idx[1]].pat_id = 2;
		    list[idx[2]].pat_id = 4;
		    pat_alig            = 8;
		  } else {
		    list[idx[1]].pat_id = 4;
		    list[idx[2]].pat_id = 2;
		    pat_alig            = 7;
		  }
		} else if (list[idx[3]].rawy==list[idx[0]].rawy-1) {
		  list[idx[3]].pat_id = 7;
		  if (list[idx[1]].rawx==list[idx[0]].rawx) {
		    list[idx[1]].pat_id = 8;
		    list[idx[2]].pat_id = 4;
		    pat_alig            = 5;
		  } else {
		    list[idx[1]].pat_id = 4;
		    list[idx[2]].pat_id = 8;
		    pat_alig            = 6;
		  }
		}
	      } else { assert(0==1); }

	    } else {
	      pat_type = -1;
	    }
	    break;

	  default:
	    pat_type=-1;
	    break;
	  }
	  // END of inner switch(nlist).

	  // Set the pattern type and alignment (orientation).
	  if (pat_type > 0) { 
	    // Valid pattern!
	    for (kk=0; kk<nlist; kk++) {
	      list[kk].pat_type = pat_type;
	      list[kk].pat_alig = pat_alig;
	    }
	  } else if (-1==pat_type) {
	    // Invalid pattern!
	    for (kk=0; kk<nlist; kk++) {
	      list[kk].pat_type = -1;
	      list[kk].pat_id   =  0;
	      list[kk].pat_alig =  0;
	    }
	  } else {
	    *status=EXIT_FAILURE;
	    HD_ERROR_THROW("Error: Could not determine pattern type!\n", *status);
	    return;
	  }
	  break; // END of double, triple, or quadruple.

	default: // Invalid pattern with more than 4 split partners.
	  for (kk=0; kk<nlist; kk++) {
	    list[kk].pat_type = -1;
	    list[kk].pat_id   =  0;
	    list[kk].pat_alig =  0;
	  }
	  break;
	}
	// END of outer switch(nlist): how many split partners?

	// Store the split pattern information in the output event file.
	for (kk=0; kk<nlist; kk++) {
	  addGenEvent2File(file, &list[kk], status);
	}
	// END of adding the split data to the output event file.
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
  GenEventFile* file=NULL;
  // Detector pixel array.
  GenEvent** pixels=NULL;

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
    det->eventfile=openGenEventFile(parameters.input_eventlist_filename, 
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
    strcat(template, det->eventfile_template);
    // Open a new event file from the specified template.
    file = openNewGenEventFile(parameters.output_eventlist_filename, template, &status);
    if (EXIT_SUCCESS!=status) break;

    // Allocate memory for the pixel array.
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
	GenPatId(det, pixels, file, &status);
	
	// Delete the old events in the pixel array.
	clearGenPatPixels(det, pixels);

	// Update the frame counter.
	frame = event.frame;
      }

      if (0==last_loop) {
	// Add the event to the pixel array.
	pixels[event.rawx][event.rawy] = event;
      }
    };
    if (EXIT_SUCCESS!=status) break;
    // END of loop over all events in the FITS file.

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
  destroyGenEventFile(&file, &status);
  
  if (status == EXIT_SUCCESS) headas_chat(3, "finished successfully\n\n");
  return(status);
}



////////////////////////////////////////////////////////////////
// This routine reads the program parameters using the PIL.
int getpar(struct Parameters* const parameters)
{
  int status=EXIT_SUCCESS; // Error status

  // Get the name of the input event list file (FITS file).
  if ((status = PILGetFname("input_eventlist_filename", 
			    parameters->input_eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the input event list file!\n", status);
  }

  // Get the name of the output event list file (FITS file).
  else if ((status = PILGetFname("output_eventlist_filename", 
				 parameters->output_eventlist_filename))) {
    HD_ERROR_THROW("Error reading the name of the output event list file!\n", status);
  }

  // Get the name of the detector XML description file (FITS file).
  else if ((status = PILGetFname("xml_filename", parameters->xml_filename))) {
    HD_ERROR_THROW("Error reading the name of the detector definition XML file!\n", status);
  }

  return(status);
}


