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


///////////////////////////////////////////////////////////////////////////////
// Load spectra from PHA files using the method "get_spectrum".
int get_spectra(
		struct Spectrum_Store *spectrum_store,
		long NumberChannels, // number of PHA channels in the detector EBOUNDS
		char filenames[N_SPECTRA_FILES][FILENAME_LENGTH],
		int Nfiles           // number of spectra to be read
		)
{
  int status=EXIT_SUCCESS;
  
  do {  // beginning of ERROR handling loop
    // get memory
    spectrum_store->spectrum = (Spectrum *) 
      malloc(Nfiles * sizeof(Spectrum));
    if (!spectrum_store->spectrum) {
      status = EXIT_FAILURE;
      HD_ERROR_THROW("Not enough memory available to store the source spectra!\n", status);
      break;
    }
    
    int count;
    for (count=0; count < Nfiles; count++) {
      if ((status=get_spectrum(&(spectrum_store->spectrum[count]), NumberChannels, 
			       filenames[count]))!=EXIT_SUCCESS) break;
    }

    spectrum_store->nspectra = Nfiles;

    /*
    // Read a FITS PHA file and assign its content to the array
    // of PHA spectra in the Spectrum_Store.
    if ((status=assign_pha_spectrum(spectrum_store, filenames[0]))!=EXIT_SUCCESS) break;
    
    // Check if the number of PHA channels in the spectrum is equivalent to
    // the number of channels in the spectrum:
    if (spectrum_store->pha_spectrum[0].NumberChannels != NumberChannels) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: number of PHA channels in spectrum '%s' (%ld) is not "
	      "equivalent to the number of channels in the EBOUNDS table (%ld)!\n",
	      filenames[0], spectrum_store->pha_spectrum[0].NumberChannels, 
	      NumberChannels);
      HD_ERROR_THROW(msg, status);
      break;
    }

    float sum=0.;
    for(count=0; count<spectrum_store->pha_spectrum[0].NumberChannels; count++) {
      sum+= spectrum_store->pha_spectrum->Pha[count];
    }
    printf("--- PHA sum: %lf ---\n", sum);
    */
  } while (0);  // END of ERROR handling loop

  return (status);
}



/*
//////////////////////////////////////////////////////////////////////
int assign_pha_spectrum(struct Spectrum_Store* store, char* filename)
{
  fitsfile* fptr;

  int status = EXIT_SUCCESS;
  char msg[MAXMSG];

  // Allocate memory:
  store->pha_spectrum = (struct PHA*)malloc(sizeof(struct PHA));
  if (store->pha_spectrum==NULL) {
    status=EXIT_FAILURE;
    sprintf(msg, "Error: could not allocate memory for RMF!\n");
    HD_ERROR_THROW(msg, status);
    return(status);
  }

  // Read the PHA spectrum from a given FITS file using the corresponding routine
  // of libhdsp.
  fits_open_file(&fptr, filename, READONLY, &status);
  if (status != EXIT_SUCCESS) return(status);
  if ((status=ReadPHAtypeI(fptr, 1, store->pha_spectrum))!=EXIT_SUCCESS) return(status);
  fits_close_file(fptr, &status);
  
  return(status);
}
*/


//////////////////////////////////////////////////////////////////////
int get_spectrum(
		 Spectrum *spectrum,
		 long Nchannels,
		 char filename[FILENAME_LENGTH]
		 )
{
  fitsfile *pha_fptr=NULL;

  int status=EXIT_SUCCESS;  // error handling variable
  char msg[MAXMSG];         // error description output buffer
  

  do { // Beginning of ERROR handling loop.

    // Fill the spectrum array with data from the FITS file.
    headas_chat(5, "load spectrum from file '%s' ...\n", filename);

    // First of all open the PHA FITS file:
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

    // If the current HDU is an image extension, throw an error message:
    if (hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in file '%s' is not a table "
	      "but an image (HDU number: %d)\n", filename, hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Determine the number of rows in the FITS file (number of given points of time):
    fits_get_num_rows(pha_fptr, &spectrum->NumberChannels, &status);

    if (spectrum->NumberChannels != Nchannels) {
      headas_chat(0, "Warning: number of PHA channels in spectrum file '%s' is not "
		  "equivalent to number of detector channels!\n", filename);
    }
  

    // Get memory for the spectrum:
    spectrum->rate = (float *) malloc(spectrum->NumberChannels * sizeof(float));
    if (!(spectrum->rate)) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available to store the source spectrum!\n");
      HD_ERROR_THROW(msg,status);
      break;
    }

    // Read the spectrum from the PHA FITS file and store it in the array.
    long row;
    long channel=0;
    float probability=0., sum=0., normalization=0.;
    for (row=1,sum=0.;(row <= spectrum->NumberChannels)&&(status==EXIT_SUCCESS); row++) {
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
      
      // Store the count rate for the individual PHA channels:
      normalization += probability;
      spectrum->rate[row-1] = probability;
    }
    if (status != EXIT_SUCCESS) break;

    // Normalize spectrum to 1 and create probability distribution function:
    for (row=0; row<spectrum->NumberChannels; row++) {
      sum += spectrum->rate[row] / normalization;
      spectrum->rate[row] = sum; 
    }
    
  } while (0);  // end of error handling loop  


  // clean up:
  if (pha_fptr) fits_close_file(pha_fptr, &status);
  
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
	HD_ERROR_THROW(msg,status);
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
    
  } while (0);  // end of error handling loop  

  // clean up:
  if (fptr) fits_close_file(fptr, &status);
  
  return(status);
}



////////////////////////////////////////////////
// Release memory of spectrum array
void free_spectra(struct Spectrum_Store *spectrum_store, long Nfiles)
{
  int count;

  if (spectrum_store->spectrum!=NULL) {
    for (count=0; count<Nfiles; count++) {
      if (spectrum_store->spectrum[count].rate!= NULL) {
	free(spectrum_store->spectrum[count].rate);
      }
    }
    free(spectrum_store->spectrum);
    spectrum_store->spectrum = NULL;
  }
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

