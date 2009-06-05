/** 
 * Contains routines/functions needed to handle source catalogs.
 */

#include <stdio.h>
#include <stdlib.h>
//#include <stdarg.h>

#include "astrosources.h"
#include "vector.h"


///////////////////////////////////////////////////////////////////
// Function opens the specified point-source catalog files.
PointSourceFiles* get_PointSourceFiles(int nfiles, char** filename, int* status)
{
  PointSourceFiles* psf;
  char msg[MAXMSG];      // error output buffer

  do { // BEGINNING of outer error handling loop

    // Check if number of FITS files exceeds maximum allowed number:
    if(nfiles > MAX_N_POINTSOURCEFILES) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: too many point-source files for input!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }

    // Get memory for the preselected source catalog.
    psf = (PointSourceFiles*) malloc(sizeof(PointSourceFiles));
    if (psf==NULL) {
      *status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the point-source "
	      "catalog files!\n");
      HD_ERROR_THROW(msg, *status);
      break;
    }
    
    // Open the different source catalogs (FITS files)
    int file_counter;    // counter to access the different source filenames for opening
    int hdunum, hdutype; // needed to access the HDUs in the source FITS file

    for(file_counter=0; file_counter<MAX_N_POINTSOURCEFILES; file_counter++) {
      psf->files[file_counter] = NULL;
    }

    for(file_counter=0; (file_counter<nfiles)&&(*status==EXIT_SUCCESS); file_counter++) {
      do {   // BEGINNING of inner error handling loop	
	headas_chat(5, "open point-source catalog in file '%s' ...\n", 
		    filename[file_counter]);
	psf->files[file_counter]=NULL;

	// Open the source catalogue (FITS-file):
	if(fits_open_table(&psf->files[file_counter], filename[file_counter], 
			   READONLY, status)) break;
	// Get the HDU number in the source catalogue:
	if(fits_get_hdu_num(psf->files[file_counter], &hdunum) == 1) {
	  // This is the primary array,
	  // so try to move to the first extension and see if it is a table:
	  fits_movabs_hdu(psf->files[file_counter], 2, &hdutype, status);
	} else {
	  // Get the HDU type:
	  fits_get_hdu_type(psf->files[file_counter], &hdutype, status);
	}
	// Image HDU results in an error message:
	if (hdutype == IMAGE_HDU) {
	  *status=EXIT_FAILURE;
	  sprintf(msg, "Error: FITS extension in source catalog file '%s' is "
		  "not a table but an image (HDU number: %d)!\n", 
		  filename[file_counter], hdunum);
	  HD_ERROR_THROW(msg, *status);
	  break;
	}

	// Determine the column numbers of the r.a., declination and count rate 
	// columns in the FITS table
	fits_get_colnum(psf->files[file_counter], CASEINSEN, "R_A_", 
			&psf->columns[file_counter][0], status);
	fits_get_colnum(psf->files[file_counter], CASEINSEN, "Dec_", 
			&psf->columns[file_counter][1], status);
	fits_get_colnum(psf->files[file_counter], CASEINSEN, "src_cps", 
			&psf->columns[file_counter][2], status);
      } while (0);  // END of inner error handling loop

    } // END of loop for opening the source catalog files

    psf->nfiles = nfiles;  // Store the number of open FITS files.

  } while (0);   // END of outer error handling loop
    
  return(psf);
}




////////////////////////////////////////////////////////////
// Releases the memory which has been allocated to store 
// the source catalogs (all sources and preselected).
void free_PointSourceFiles(PointSourceFiles* psf, int* status)
{
  int counter;
  
  if (psf!=NULL) {
    // Close the FITS files with the source catalogs:
    for (counter = 0; counter < psf->nfiles; counter++) {
      if(psf->files[counter] != NULL) fits_close_file(psf->files[counter], status);
    }
    free(psf);
    psf=NULL;
  }
}






//////////////////////////////////////////////////////////////////////////
int get_PointSourceCatalog(
			   PointSourceFiles* psf,
			   PointSourceCatalog** psc,
			   Vector normal_vector,
			   const double max_align,
			   struct Spectrum_Store spectrum_store
			   )
{
  char msg[MAXMSG];          // error output buffer
  int status=EXIT_SUCCESS;
  
  do { // beginning of outer ERROR handling loop

    // Read-in the different source catalogs
    // from the FITS files and store the sources in an array:

    // field-indices of r.a., dec., and countrate in the different FITS input file tables
    // const int columns[4][3] = {{2,3,7}, {2,3,7}, {1,2,3}, {1,2,3}};
    //                   |  |-> number of columns (should be "3")
    //                   |-> number of different source catalogs
    
    int file_counter;   // counter to access the different source filenames for opening
    double ra, dec;     // source position
    float countrate;    // source countrate

    // Allocate memory:
    if (*psc==NULL) {
      *psc = (PointSourceCatalog*)malloc(sizeof(PointSourceCatalog));
      if (*psc==NULL) {
	status = EXIT_FAILURE;
	sprintf(msg, "Error: not enough memory available to store the point-source "
		"catalog!\n");
	HD_ERROR_THROW(msg, status);
	break;
      } else {
	(*psc)->sources=NULL;
      }
    }
    
    if ((*psc)->sources==NULL) {
      (*psc)->sources = (PointSource*)malloc(MAX_N_POINTSOURCES*sizeof(PointSource));
      if((*psc)->sources==NULL) {
	status = EXIT_FAILURE;
	sprintf(msg, "Error: not enough memory available to store the point-source "
		"catalog!\n");
	HD_ERROR_THROW(msg, status);
	break;
      } 
    } else {
      // Free the light  curves in the old PSC:
      int count;
      for (count=0; count<(*psc)->nsources; count++) {
	if ((*psc)->sources[count].lightcurve != NULL) {
	  free((*psc)->sources[count].lightcurve);
	  (*psc)->sources[count].lightcurve = NULL;
	}
      }
    }// if (*psc)->sources == NULL
    // END of memory allocation


    // At the moment the selected catalog is empty:
    (*psc)->nsources = 0;

    for(file_counter=0; (file_counter<psf->nfiles)&&(status==EXIT_SUCCESS); 
	file_counter++) {

      do {  // beginning of inner ERROR handling loop

	// Determine the number of rows in the source catalogue 
	// (i.e. number of listed sources):
	long nrows;           // number of rows in source catalog
	fits_get_num_rows(psf->files[file_counter], &nrows, &status);

	// Read-in all sources from the individual table rows, one after another:
	int source_counter; // counter to access the individual sources in the catalog
	for (source_counter=0; (source_counter<nrows)&&(status==EXIT_SUCCESS); 
	     source_counter++) {
	  // Read source data (right asension, declination, countrate, ...):
	  if (get_srctbl_row(psf->files[file_counter], source_counter, 
			     psf->columns[file_counter], 
			     &ra, &dec, &countrate, &status)) 
	    break;

	  // Rescale from [deg] -> [rad]
	  ra  *= M_PI/180.; 
	  dec *= M_PI/180.;
	  // Get a unit vector pointing in the direction of the source:
	  Vector source_direction = unit_vector(ra, dec);
	  
	  // Check whether the source should be added to the preselected catalog:
	  if(fabs(scalar_product(&source_direction, &normal_vector)) < max_align) {
	    if((*psc)->nsources > MAX_N_POINTSOURCES) {
	      // Too many sources !
	      status=EXIT_FAILURE;
	      sprintf(msg, "Error: too many sources (%ld)!\n", (*psc)->nsources+1);
	      HD_ERROR_THROW(msg, status);
	      break;
	    }

	    // Add the current source to the selected catalog:
	    (*psc)->sources[(*psc)->nsources].rate = countrate;

	    // save the source direction in the source catalog-array:
	    (*psc)->sources[(*psc)->nsources].ra = ra;
	    (*psc)->sources[(*psc)->nsources].dec = dec;
	    // set lightcurve pointer to NULL
	    (*psc)->sources[(*psc)->nsources].lightcurve = NULL;
	    // so far there was no photon created for this source
	    (*psc)->sources[(*psc)->nsources].t_last_photon = -1.;
	    // source spectrum
	    (*psc)->sources[(*psc)->nsources].spectrum = 
	      &(spectrum_store.spectrum[0]);

	    // increase number of sources in the selected catalog
	    (*psc)->nsources++;
	  }
	} // END of loop over all sources in the current catalog file
      } while (0); // END of inner error handling loop
    } // END of scanning the source catalogs
  } while (0); // END of outer error handling loop

  return(status);
}


