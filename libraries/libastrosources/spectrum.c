#include "spectrum.h"


///////////////////////////////////////////////////////////////////////////////
int loadSpectra(fitsfile* source_fptr, SpectrumStore* store)
{
  char comment[MAXMSG]; // String buffers.
  int status=EXIT_SUCCESS;

  do { // Beginning of Error Handling Loop.

    // Try to find out the number of SPECTRA from the FITS header of the source file.
    if (fits_read_key(source_fptr, TLONG, "NSPECTRA", &store->nspectra, comment, 
		      &status)) break;

    if (store->nspectra<1) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Error: No spectra specified in point source catalog!\n", status);
      break;
    }

    // Get memory to store the spectra.
    store->spectrum = (Spectrum *)malloc(store->nspectra*sizeof(Spectrum));
    if (NULL==store->spectrum) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Not enough memory available to store the source spectra!\n", status);
      break;
    }
    
    // Load the spectra from the PHA files.
    long count;
    char filename[MAXMSG], key[MAXMSG];
    for (count=0; (EXIT_SUCCESS==status)&&(count<store->nspectra); count++) {
      // Determine the name of the source file from the FITS header.
      sprintf(key, "SPEC%04ld", count+1);
      if (fits_read_key(source_fptr, TSTRING, key, filename, comment, &status)) break;

      // Load the spectrum from the named PHA file and add it to the SpectrumStore.
      status = loadSpectrum(&store->spectrum[count], filename);
      if (EXIT_SUCCESS!=status) break;
    }
    if (EXIT_SUCCESS!=status) break;

  } while(0); // End of Error Handling Loop.

  return(status);
}



int loadSpectrum(Spectrum* spectrum, char* filename)
{
  fitsfile *fptr=NULL;

  int status=EXIT_SUCCESS;  // error handling variable
  char msg[MAXMSG];         // error description output buffer
  
  do { // Beginning of ERROR handling loop.

    // Fill the spectrum array with data from the FITS file.
    headas_chat(5, "load spectrum from file '%s' ...\n", filename);

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

    /*
    // Plot the spectrum to an output file for testing.
    headas_printf("Output of spectrum to file 'spectrum.dat' ...\n");
    FILE* spectrum_file=fopen("spectrum.dat", "w+");
    if (NULL==spectrum_file) {
      status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: Could not open spectrum output file!\n", status);
      break;
    }
    for (row=0; row<spectrum->NumberChannels; row++) {
      fprintf(spectrum_file, "%ld %lf\n", row, spectrum->rate[row]);
    }
    fclose(spectrum_file);
    */

    // Normalize spectrum to 1, i.e., create probability distribution function:
    for (row=0; row<spectrum->NumberChannels; row++) {
      sum += spectrum->rate[row] / normalization;
      spectrum->rate[row] = sum; 
    }
    
  } while (0);  // END of error handling loop  

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
  for(count=0; count<store->nspectra; count++) {
    cleanupSpectrum(&store->spectrum[count]);
  }
  store->nspectra=0;
}

