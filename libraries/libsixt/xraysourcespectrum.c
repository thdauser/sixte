#include "xraysourcespectrum.h"


// Conversion factor from [keV] -> [erg].
const float keV2erg = 1.602e-9;


XRaySourceSpectrum* newXRaySourceSpectrum(int* const status)
{
  XRaySourceSpectrum* spec = (XRaySourceSpectrum*)malloc(sizeof(XRaySourceSpectrum));
  CHECK_NULL(spec, *status,
	     "memory allocation for XRaySourceSpectrum failed");

  // Initalize pointers with NULL.
  spec->emin = NULL;
  spec->emax = NULL;
  spec->flux = NULL;
  spec->ratedistr = NULL;

  // Initialize values.
  strcpy(spec->filename, "");

  return(spec);
}


void freeXRaySourceSpectrum(XRaySourceSpectrum** spec)
{
  if (NULL!=*spec) {
    if (NULL!=(*spec)->emin) {
      free((*spec)->emin);
    }    
    if (NULL!=(*spec)->emax) {
      free((*spec)->emax);
    }    
    if (NULL!=(*spec)->flux) {
      free((*spec)->flux);
    }    
    if (NULL!=(*spec)->ratedistr) {
      free((*spec)->ratedistr);
    }    
    free(*spec);
    *spec=NULL;
  }
}


float getRndSpectrumEnergy(const XRaySourceSpectrum* const spec)
{
  // Get a random energy bin according to the given spectral distribution.
  float rand = (float)sixt_get_random_number() * spec->ratedistr[spec->nbins-1];
  int upper = spec->nbins-1, lower=0, mid;
  
  // Determine the energy of the photon (using binary search).
  while (upper>lower) {
    mid = (lower+upper)/2;
    if (spec->ratedistr[mid] < rand) {
      lower = mid+1;
    } else {
      upper = mid;
    }
  }

  // Return an energy chosen randomly out of the determined PHA bin:
  return(spec->emin[lower] + 
	 sixt_get_random_number()*(spec->emax[lower]-spec->emin[lower]));
}


XRaySourceSpectrum* loadXRaySpectrumFilename(const char* const filename,
					     int* const status)
{
  headas_chat(3, "load spectrum from file '%s' ...\n", filename);
  
  XRaySourceSpectrum* spec = newXRaySourceSpectrum(status);
  CHECK_STATUS_RET(*status, spec);
  
  // Store the filename.
  strcpy(spec->filename, filename);

  // Open the FITS file:
  fitsfile* fptr;
  fits_open_file(&fptr, spec->filename, READONLY, status);
  CHECK_STATUS_RET(*status, spec);
  
  int hdunum, hdutype;
  // After opening the FITS file, get the number of the current HDU.
  if (1==fits_get_hdu_num(fptr, &hdunum)) {
    // This is the primary array, so try to move to the first extension 
    // and see if it is a table.
    fits_movabs_hdu(fptr, 2, &hdutype, status);
  } else {
    // Get the HDU type.
    fits_get_hdu_type(fptr, &hdutype, status);
  }
  CHECK_STATUS_RET(*status, spec);
  
  // If the current HDU is an image extension, throw an error message:
  if (IMAGE_HDU==hdutype) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: FITS extension in file '%s' is not a table "
	    "but an image (HDU number: %d)\n", spec->filename, hdunum);
    SIXT_ERROR(msg);
    return(spec);
  }

  // Determine the number of bins.
  fits_get_num_rows(fptr, &spec->nbins, status);
  CHECK_STATUS_RET(*status, spec);

  // Get memory for the spectrum:
  spec->flux = (float*)malloc(spec->nbins*sizeof(float));
  CHECK_NULL(spec->flux, *status, "memory allocation for spectrum failed");
  spec->emin = (float*)malloc(spec->nbins*sizeof(float));
  CHECK_NULL(spec->flux, *status, "memory allocation for spectrum failed");
  spec->emax = (float*)malloc(spec->nbins*sizeof(float));
  CHECK_NULL(spec->flux, *status, "memory allocation for spectrum failed");

  // Determine the column number in the FITS table.
  int cemin, cemax, cflux;
  fits_get_colnum(fptr, CASEINSEN, "E_MIN", &cemin, status);
  fits_get_colnum(fptr, CASEINSEN, "E_MAX", &cemax, status);
  fits_get_colnum(fptr, CASEINSEN, "FLUX", &cflux, status);
  CHECK_STATUS_RET(*status, spec);

  // Read the spectrum from the FITS file and store it in the array.
  long row;
  for (row=0; row<spec->nbins; row++) {

    // Set default values.
    spec->emin[row] = 0.;
    spec->emax[row] = 0.;
    spec->flux[row] = 0.;

    int anynul;
    fits_read_col(fptr, TFLOAT, cemin, row+1, 1, 1, &spec->emin[row], 
		  &spec->emin[row], &anynul, status);
    fits_read_col(fptr, TFLOAT, cemax, row+1, 1, 1, &spec->emax[row], 
		  &spec->emax[row], &anynul, status);
    fits_read_col(fptr, TFLOAT, cflux, row+1, 1, 1, &spec->flux[row], 
		  &spec->flux[row], &anynul, status);
    CHECK_STATUS_RET(*status, spec);

  }
  // END of loop over all entries in the table.

  // Close the FITS file.
  fits_close_file(fptr, status);
  CHECK_STATUS_RET(*status, spec);

  return(spec);
}


void applyARF2Spectrum(XRaySourceSpectrum* const spec, 
		       const struct ARF* const arf, 
		       int* const status)
{
  // Check if we have to allocate memory.
  if (NULL==spec->ratedistr) {
    spec->ratedistr = (float*)malloc(spec->nbins*sizeof(float));
    CHECK_NULL_VOID(spec->ratedistr, *status,
		    "memory allocation for spectrum failed");
  }

  // Multiply each bin by the ARF and the width of the bin.
  // [photons/s/cm^2/keV] -> [photons/s]
  // The ARF contribution corresponding to a particular spectral bin 
  // might be obtained by interpolation.
  long ii;
  for (ii=0; ii<spec->nbins; ii++) {
    // Determine the ARF contribution by interpolation.
    float arf_contribution = 0.;
    long jj;
    for (jj=0; jj<arf->NumberEnergyBins; jj++) {
      if ((arf->LowEnergy[jj]<spec->emax[ii]) && 
	  (arf->HighEnergy[jj]>spec->emin[ii])) {
	float emin = MAX(arf->LowEnergy[jj], spec->emin[ii]);
	float emax = MIN(arf->HighEnergy[jj], spec->emax[ii]);
	assert(emax>emin);
	arf_contribution += arf->EffArea[jj] * (emax-emin);
      }
    }

    spec->ratedistr[ii] = spec->flux[ii] * arf_contribution;

    // Generate a rate distribution. (Similar to a probability
    // distribution, but not normalized to 1, but the final value
    // represents the total photon rate.)
    if (ii>0) {
      spec->ratedistr[ii] += spec->ratedistr[ii-1];
    }
  }
  // END of loop over all bins.
}


float getSpectralEnergyFlux(const XRaySourceSpectrum* const spec, 
			    const float emin, const float emax)
{
  float flux = 0.;

  long ii;
  for (ii=0; ii<spec->nbins; ii++) {
    if ((spec->emin[ii]<emax) && (spec->emax[ii]>emin)) {
      float min = MAX(spec->emin[ii], emin);
      float max = MIN(spec->emax[ii], emax);
      assert(max>min);
      flux += (max-min) * spec->flux[ii] * 0.5*(spec->emin[ii]+spec->emax[ii]);
    }
  }

  // Convert units of 'flux' from [keV/s/cm^2] -> [erg/s/cm^2].
  flux *= keV2erg;

  return(flux);
}


float getSpectralPhotonRate(const XRaySourceSpectrum* const spec,
			    const float emin, const float emax)
{
  // Photon rate.
  float pr = 0.;

  long ii;
  for (ii=0; ii<spec->nbins; ii++) {
    if ((spec->emin[ii]<emax) && (spec->emax[ii]>emin)) {
      float min = MAX(spec->emin[ii], emin);
      float max = MIN(spec->emax[ii], emax);
      assert(max>min);
      float ratediff = spec->ratedistr[ii];
      ratediff -= (0==ii) ? 0. : spec->ratedistr[ii-1];
      pr += (max-min)/(spec->emax[ii]-spec->emin[ii]) * ratediff;
    }
  }

  return(pr);
}


