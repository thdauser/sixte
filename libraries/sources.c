/** 
 * Contains routines/functions needed to handle source catalogs.
 */

#include "astrosources.h"


// Function loads the required source catalogs from the specified FITS files and
// stores the source data in an array. Additionally it allocates the memory for the
// preselected catalog.
int get_source_catalogs(
			// preselected catalog
			struct source_cat_entry **selected_catalog,
			// number of FITS sources files to be loaded
			const int n_sourcefiles,
			// FITS file pointers to the input source catalogs
			fitsfile **sourcefiles,
			// column numbers for r.a., declination, and count rate
			int columns[5][3],
			// filenames of the source files
			char source_filename[MAX_NSOURCEFILES][FILENAME_LENGTH] 
      			)
{
  int status=0;     // error handling variable
  char msg[MAXMSG]; // error output buffer

  do {  // BEGINNING of outer error handling loop

    // Get memory for the preselected source catalog.
    *selected_catalog = (struct source_cat_entry *) 
      malloc(MAX_NSOURCES_PRE*sizeof(struct source_cat_entry));
    if (!(*selected_catalog)) {
      status = EXIT_FAILURE;
      sprintf(msg, "Not enough memory available to store the preselected "
	      "source catalog!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Open the different source catalogs (FITS files)
    int file_counter;    // counter to access the different source filenames for opening
    int hdunum, hdutype; // needed to access the HDUs in the source FITS file

    for(file_counter=0; (file_counter<n_sourcefiles)&&(status==EXIT_SUCCESS); 
	file_counter++) {
      do {   // BEGINNING of inner error handling loop
	headas_chat(5, "open source catalog in file '%s' ...\n", 
		    source_filename[file_counter]);
	sourcefiles[file_counter]=NULL;
	// open the source catalogue (FITS-file)
	if(fits_open_file(&sourcefiles[file_counter], source_filename[file_counter], 
			  READONLY, &status)) break;
	// get the HDU number in the source catalogue
	if (fits_get_hdu_num(sourcefiles[file_counter], &hdunum) == 1) {
	  // This is the primary array,
	  // so try to move to the first extension and see if it is a table
	  fits_movabs_hdu(sourcefiles[file_counter], 2, &hdutype, &status);
	} else {
	  // get the HDU type
	  fits_get_hdu_type(sourcefiles[file_counter], &hdutype, &status);
	}

	// image HDU results in an error message
	if (hdutype == IMAGE_HDU) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: FITS extension in source catalog file '%s' is "
		  "not a table but an image (HDU number: %d)!\n", 
		  source_filename[file_counter], hdunum);
	  HD_ERROR_THROW(msg,status);
	  break;
	}

	// Determine the column numbers of the r.a., declination and count rate 
	// columns in the FITS table
	fits_get_colnum(sourcefiles[file_counter], CASEINSEN, "R_A_", 
			&columns[file_counter][0], &status);
	fits_get_colnum(sourcefiles[file_counter], CASEINSEN, "Dec_", 
			&columns[file_counter][1], &status);
	fits_get_colnum(sourcefiles[file_counter], CASEINSEN, "src_cps", 
			&columns[file_counter][2], &status);
      } while (0);  // END of inner error handling loop

    } // END of reading the source catalogs

  } while (0);   // END of outer error handling loop

  return(status);
}



////////////////////////////////////////////////////////////
// Releases the memory which has been allocated to store 
// the source catalogs (all sources and preselected).
void free_source_catalogs(
			  fitsfile **sourcefiles,  // pointers to source catalog files
			  const int n_sourcefiles, // number of source catalog files
			  // preselected source catalog
			  struct source_cat_entry **selected_catalog, 
			  int *status              // error status variable
			  )
{
  int counter;
  
  // free preselected source catalog
  if (*selected_catalog) {
    free(*selected_catalog);
    *selected_catalog=NULL;
  }

  // close the FITS files to the source catalogs
  for (counter = 0; counter < n_sourcefiles; counter++) {
    if(sourcefiles[counter]) fits_close_file(sourcefiles[counter], status);
  }

}







//////////////////////////////////////////////////////////////////////////
// Get the preselected source catalogs with sources along the path of the 
// telescope axis over the sky.
int get_preselected_catalog(
			    // preselected catalog: pointers to source_catalog
			    struct source_cat_entry *selected_catalog, 
			    // total number of sources in selected catalog
			    long *nsources,                    
			    // number of FITS sources files to be loaded
			    const int n_sourcefiles,           
			    // FITS file pointers to the original source catalogs
			    fitsfile **sourcefiles,            
			    // column numbers for r.a., declination, and count rate
			    int columns[5][3],                 
			    // pointing direction of the telescope axis
			    struct vector telescope_direction, 
			    // maximum alignment of source and telescope direction
			    // for beeing chosen to the selected catalog
			    const double pre_max_align,        
			    // storage for individual source spectra
			    struct Spectrum_Store spectrum_store,  
			    // number of different source spectra
			    const int Nspectra    
			    )
{
  int status=0;              // error handling variable
  char msg[MAXMSG];          // error output buffer

  
  do {  // beginning of outer error handling loop

    // Read-in the different source catalogs
    // from the FITS files and store the sources in an array:

    // field-indices of r.a., dec., and countrate in the different FITS input file tables
    // const int columns[4][3] = {{2,3,7}, {2,3,7}, {1,2,3}, {1,2,3}};
    //                   |  |-> number of columns (should be "3")
    //                   |-> number of different source catalogs
    
    int file_counter;     // counter to access the different source filenames for opening
    long nrows;           // number of rows in source catalog
    int source_counter;   // counter to access the individual sources in the catalog
    double rasc, dec;     // source position
    float countrate;      // source countrate
    struct vector source_direction; // source direction (unit vector)
    *nsources=0;          // at the moment the selected catalog is empty

    for(file_counter=0; (file_counter<n_sourcefiles)&&(status==EXIT_SUCCESS); 
	file_counter++) {

      do {  // beginning of inner ERROR handling loop

	// Determine the number of rows in the source catalogue 
	// (i.e. number of listed sources):
	fits_get_num_rows(sourcefiles[file_counter], &nrows, &status);

	// read-in all sources from the individual table rows, one after another
	for (source_counter=0; (source_counter<nrows)&&(status==EXIT_SUCCESS); 
	     source_counter++) {
	  // read source data (right asension, declination, countrate, ...)
	  if (get_srctbl_row(sourcefiles[file_counter], source_counter, 
			     columns[file_counter], &rasc, &dec, &countrate, &status)) 
	    break;

	  // get a unit vector pointing in the direction of the source
	  source_direction = unit_vector(rasc*M_PI/180., dec*M_PI/180.);
	  
	  // check, whether the source should be added to the preselected catalog:
	  if(fabs(scalar_product(source_direction, telescope_direction))<pre_max_align) {
	    if(*nsources > MAX_NSOURCES_PRE) {
	      // too many sources
	      status=EXIT_FAILURE;
	      sprintf(msg, "Error: too many sources (%ld)!\n", *nsources);
	      HD_ERROR_THROW(msg,status);
	      break;
	    }

	    // add source to the preselected catalog
	    selected_catalog[*nsources].rate = countrate;

	    // save the source direction in the source catalog-array:
	    selected_catalog[*nsources].r = source_direction;  // REMOVE
	    selected_catalog[*nsources].ra = rasc;
	    selected_catalog[*nsources].dec = dec;
	    // set lightcurve pointer to NULL
	    selected_catalog[*nsources].lightcurve = NULL;
	    // so far there was no photon created for this source
	    selected_catalog[*nsources].t_last_photon = -1.;
	    // source spectrum   TODO: different spectra
	    selected_catalog[*nsources].spectrum = &(spectrum_store.spectrum[0]);

	    // increase number of sources in preselected catalog
	    (*nsources)++;
	  }
	}  // end of loop over all sources in the current catalog file

      } while (0);  // end of inner error handling loop

    } // end of scanning the source catalogs

  } while (0);   // end of outer error handling loop

  return(status);
}


