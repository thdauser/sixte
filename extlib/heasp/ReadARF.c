/* Functions to read the ARF file */

#include <string.h>
#include <stdio.h>
#include "heasp.h"
#include "heaspio.h"

/* Read the SPECRESP extension and return the result in the ARF structure
   Assumes that the file has been opened but not positioned at the correct
   extension */

int ReadARF(fitsfile *fptr, long ARFnumber, struct ARF *arf)
{
  int status=0;
  int colnum=0;
  int anynul=0;
  long ltemp;

  /* Move to the correct SPECRESP extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "SPECRESP", ARFnumber, &status);
  if (!status) {
     headas_chat(5, "Found SPECRESP extension...\n");
  } else {
     headas_chat(1, "***Cannot find SPECRESP extension...\n");
     headas_chat(1, "   FITS status = %d\n", status);
    return(status);
  }

  /* Read the standard keywords and save the values */

  SP_read_key(fptr, TSTRING, "EXTNAME", arf->ARFExtensionName, "UNKNOWN");

  if (SP_read_key(fptr, TSTRING, "HDUVERS", arf->ARFVersion, "UNKNOWN"))
    SP_read_key(fptr, TSTRING, "HDUVERS2", arf->ARFVersion, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "TELESCOP", arf->Telescope, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "INSTRUME", arf->Instrument, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "DETNAM", arf->Detector, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "FILTER", arf->Filter, "UNKNOWN");

  status = SP_read_key(fptr, TLONG, "NAXIS2", &arf->NumberEnergyBins, &ltemp);
  if (status) return(status);

  /* Get the start and stop energies for the bins and the effective areas */

  headas_chat(5, "Allocating %ld for energy bins\n", arf->NumberEnergyBins);

  arf->LowEnergy = (float *) malloc(arf->NumberEnergyBins*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "ENERG_LO", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, arf->NumberEnergyBins, NULL, arf->LowEnergy, &anynul, &status);

  arf->HighEnergy = (float *) malloc(arf->NumberEnergyBins*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "ENERG_HI", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, arf->NumberEnergyBins, NULL, arf->HighEnergy, &anynul, &status);

  arf->EffArea = (float *) malloc(arf->NumberEnergyBins*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "SPECRESP", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, arf->NumberEnergyBins, NULL, arf->EffArea, &anynul, &status);


  /* Check for a cfitsio error and if it occurred write diagnostic information */
  if (status) {
    fits_report_error(stderr, status);
    printf("\n");
  }

  return(status);

}

