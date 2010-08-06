#include "spectrum.h"


int loadSpectra(fitsfile* fptr, SpectrumStore* store)
{
  char comment[MAXMSG]; // String buffers.
  int status=EXIT_SUCCESS;

  // Flag whether the spectral name have to be read from the header (1) 
  // or from a binary extension (0).
  int fheader=0; 
  // Number of the current HDU. This has to be determined before shifting
  // the FITS file pointer to the spectrum binary table extension, such 
  // that the original HDU can be restored afterwards.
  int original_hdu=-1;

  do { // Beginning of Error Handling Loop.

    // Try to find out the number of SPECTRA from the FITS header of the source file.
    // First set an error mark on the FITS error stack in order to be able 
    // to delete error messages, if the keyword "NSPECTRA" does not exist.
    fits_write_errmark();
    fits_read_key(fptr, TLONG, "NSPECTRA", &store->nspectra, comment, 
		  &status);
    if (0==status) {
      // The filenames of the spectra are stored in the header.
      fheader=1;

    } else if (KEY_NO_EXIST==status) { //((VALUE_UNDEFINED==status)||
      // The keyword "NSPECTRA" does not exist.
      fheader=0;
      // Clear the error message that the header keyword "NSPECTRA" does not exist.
      fits_clear_errmsg();
      fits_clear_errmark();

      // Get Number of lines in the binary table extension.
      if (fits_get_hdu_num(fptr, &original_hdu)) break;
      if (fits_movnam_hdu(fptr, BINARY_TBL, "SPECTRA", 0, &status)) break;
      if (fits_get_num_rows(fptr, &store->nspectra, &status)) break;

    } else { // Other error.
      break;
    }

    // Check if the point source file specifies any spectra.
    if (store->nspectra<1) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: No spectra specified in point source catalog!\n", status);
      break;
    }
    headas_chat(5, "load required spectra (%d):\n", store->nspectra);

    // Get memory to store the spectra.
    store->spectrum = (Spectrum *)malloc(store->nspectra*sizeof(Spectrum));
    if (NULL==store->spectrum) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Not enough memory available to store the source spectra!\n", status);
      break;
    }
    // Initialize the data structures.
    long count;
    for (count=0; count<store->nspectra; count++) {
      store->spectrum[count].rate = NULL;
    }
    
    // Load the spectra from the PHA files.
    char filename[MAXMSG], key[MAXMSG];
    for (count=0; (EXIT_SUCCESS==status)&&(count<store->nspectra); count++) {

      // If keywords are stored in header or in binary table extension.
      if (1==fheader) {
	// Determine the name of the source file from the FITS header.
	sprintf(key, "SPEC%04ld", count+1);
	if (fits_read_key(fptr, TSTRING, key, filename, comment, &status)) break;
      } else {
	// Determine the name of the source file from the binary table extension.
	int anynul;
	if (fits_read_col(fptr, TSTRING, 1, count+1, 1, 1, "", key, &anynul, &status)) break;
      }

      // Load the spectrum from the named PHA file and add it to the SpectrumStore.
      status = loadSpectrum(&store->spectrum[count], filename);
      if (EXIT_SUCCESS!=status) break;
    }
    if (EXIT_SUCCESS!=status) break;

  } while(0); // End of Error Handling Loop.

  // Restore the original HDU.
  if (original_hdu>=0) {
    int type, status;
    if (fits_movabs_hdu(fptr, original_hdu, &type, &status)) return(status);
  }

  return(status);
}



int loadSpectrum(Spectrum* spectrum, char* filename)
{
  fitsfile *fptr=NULL;

  int status=EXIT_SUCCESS;  // error handling variable
  char msg[MAXMSG];         // error description output buffer
  
  do { // Beginning of ERROR handling loop.

    // Fill the spectrum array with data from the FITS file.
    headas_chat(5, " load spectrum from file '%s' ...\n", filename);

    // First of all open the PHA FITS file:
    if (fits_open_table(&fptr, filename, READONLY, &status)) break;
  
    int hdunum, hdutype;
    // After opening the FITS file, get the number of the current HDU.
    if (1==fits_get_hdu_num(fptr, &hdunum)) {
      // This is the primary array, so try to move to the first extension 
      // and see if it is a table.
      if (fits_movabs_hdu(fptr, 2, &hdutype, &status)) break;
    } else {
      // Get the HDU type.
      if (fits_get_hdu_type(fptr, &hdutype, &status)) break;
    }
    // If the current HDU is an image extension, throw an error message:
    if (IMAGE_HDU==hdutype) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in file '%s' is not a table "
	      "but an image (HDU number: %d)\n", filename, hdunum);
      HD_ERROR_THROW(msg, status);
      break;
    }

    // Determine the number of PHA channels by reading the corresponding
    // FITS header keyword.
    char comment[MAXMSG]; // Buffer
    if (fits_read_key(fptr, TLONG, "DETCHANS", &spectrum->NumberChannels, comment, &status)) 
      break;

    // Get memory for the spectrum:
    spectrum->rate = (float*)malloc(spectrum->NumberChannels*sizeof(float));
    if (NULL==spectrum->rate) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: not enough memory available to store the source spectrum!\n", 
		     status);
      break;
    }

    // Read the spectrum from the PHA FITS file and store it in the array.
    long row;
    long channel=0;
    float probability=0., sum=0., normalization=0.;
    for (row=1,sum=0.; (row<=spectrum->NumberChannels)&&(EXIT_SUCCESS==status); row++) {
      if ((status=read_spec_fitsrow(&channel, &probability, fptr, row))
	  !=EXIT_SUCCESS) break;

      // Check if the channel number is valid.
      // (PHA channel start at 1 !!)
      if ((channel < 0) || (channel > spectrum->NumberChannels)) {
	status=EXIT_FAILURE;
	sprintf(msg, "Error: Invalid channel number (%ld) in file '%s'!\n", channel, 
		filename);
	HD_ERROR_THROW(msg, status);
	break;
      }
      
      // Store the rates for the individual PHA channels:
      normalization += probability;
      spectrum->rate[row-1] = probability;
    }
    if (status != EXIT_SUCCESS) break;

    // Normalize spectrum to 1, i.e., create probability distribution function:
    for (row=0; row<spectrum->NumberChannels; row++) {
      sum += spectrum->rate[row] / normalization;
      spectrum->rate[row] = sum; 
    }
    // Set the last bin explicitly to one in order to avoid 
    // numercial problem due to a value slightly below 1.
    spectrum->rate[spectrum->NumberChannels-1] = 1.;
    
  } while (0); // END of error handling loop  

  // Clean up:
  if (fptr) fits_close_file(fptr, &status);
  
  return(status);
}


void cleanupSpectrum(Spectrum* spectrum)
{
  if (NULL!=spectrum->rate) {
    free(spectrum->rate);
    spectrum->rate=NULL;
  }
}


void freeSpectrumStore(SpectrumStore* store){
  long count;
  if (NULL!=store->spectrum) {
    for(count=0; count<store->nspectra; count++) {
      cleanupSpectrum(&store->spectrum[count]);
    }
    free(store->spectrum);
    store->spectrum=NULL;
  }
  store->nspectra=0;
}

