#include "spectrum.h"



///////////////////////////////////////////////////////////////////////////////
// Load spectra from PHA files using the method "get_spectrum".
int get_spectra(
		struct Spectrum_Store *spectrum_store,
		long Nchannels,                       // number of PHA channels
		char filenames[N_SPECTRA_FILES][FILENAME_LENGTH],
		int Nfiles                            // number of spectra
		)
{
  int status=EXIT_SUCCESS;  // error handling variable
  char msg[MAXMSG];         // error description output buffer
  
  do {  // beginning of error handling loop
    // get memory
    spectrum_store->spectrum = (struct Spectrum *) 
      malloc(Nfiles * sizeof(struct Spectrum));
    if (!spectrum_store->spectrum) {
      status = EXIT_FAILURE;
      sprintf(msg, "Not enough memory available to store the source spectra!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }
    
    int count;
    for (count=0; count < Nfiles; count++) {
      if ((status=get_spectrum(&(spectrum_store->spectrum[count]), Nchannels, 
			       filenames[count]))!=EXIT_SUCCESS) break;
    }

    spectrum_store->nspectra = Nfiles;

  } while (0);  // end of error handling loop

  return (status);
}



//////////////////////////////////////////////////////////////////////
int get_spectrum(
		 struct Spectrum *spectrum,
		 long Nchannels,
		 char filename[FILENAME_LENGTH]
		 )
{
  fitsfile *pha_fptr=NULL;

  int status=EXIT_SUCCESS;  // error handling variable
  char msg[MAXMSG];         // error description output buffer
  

  do {    // beginning of error handling loop
    // fill the spectrum array with data from the FITS file
    headas_chat(5, "load spectrum from file '%s' ...\n", filename);

    // first open PHA FITS file
    if (fits_open_table(&pha_fptr, filename, READONLY, &status)) break;
  
    int hdunum, hdutype;
    // after opening the FITS file, get the number of the current HDU
    if (fits_get_hdu_num(pha_fptr, &hdunum) == 1) {
      // This is the primary array, so try to move to the first extension 
      // and see if it is a table.
      if (fits_movabs_hdu(pha_fptr, 2, &hdutype, &status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(pha_fptr, &hdutype, &status)) break;
    }

    // image HDU results in an error message
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in file '%s' is not a table "
	      "but an image (HDU number: %d)\n", filename, hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // get the number of rows in the FITS file (number of given points of time)
    long nrows;
    fits_get_num_rows(pha_fptr, &nrows, &status);

    if (nrows != Nchannels) {
      headas_chat(0, "Warning: number of PHA channels in spectrum file '%s' is not "
		  "equivalent to number of detector channels!\n", filename);
    }
  

    // get memory for the spectrum
    spectrum->data = (float *) malloc(nrows * sizeof(float));
    if (!(spectrum->data)) {
      status = EXIT_FAILURE;
      sprintf(msg, "Not enough memory available to store the source spectra!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Read the spectrum from the PHA FITS file and store it in the array.
    long row;
    long channel=0;
    float probability=0., sum=0., normalization=0.;
    for (row=1, sum=0.; (row <= nrows)&&(status==EXIT_SUCCESS); row++) {
      if ((status=read_spec_fitsrow(&channel, &probability, pha_fptr, row))
	  !=EXIT_SUCCESS) break;

      // Check if the channel number is valid.
      // (PHA channel start at 1 !!)
      if ((channel < 0) || (channel > Nchannels)) {
	status=EXIT_FAILURE;
	sprintf(msg, "Error: Invalid channel number (%ld) in file '%s'!\n", channel, 
		filename);
	HD_ERROR_THROW(msg,status);
	break;
      }
      
      // store the probability distribution for the individual PHA channels
      normalization += probability;
      spectrum->data[row-1] = probability;
    }

    if (status != EXIT_SUCCESS) break;

    // Normalize spectrum to 1;
    for (row=0; row<nrows; row++) {
      // TODO: scaling to be implemented
      sum += spectrum->data[row] / normalization;
      spectrum->data[row] = sum; 
    }
    
  } while (0);  // end of error handling loop  


  // clean up:
  if (pha_fptr) fits_close_file(pha_fptr, &status);
  
  return(status);
}





////////////////////////////////////////////////
// Release memory of spectrum array
void free_spectra(struct Spectrum_Store *spectrum_store, long Nfiles)
{
  int count;

  if (spectrum_store->spectrum!=NULL) {
    for (count=0; count<Nfiles; count++) {
      if (spectrum_store->spectrum[count].data!= NULL) {
	free(spectrum_store->spectrum[count].data);
      }
    }
    free(spectrum_store->spectrum);
    spectrum_store->spectrum = NULL;
  }
}


