/* Functions to write the ARF file */

#include <string.h>
#include <stdio.h>
#include "heasp.h"
#include "heaspio.h"

/* Write the SPECRESP extension */

int WriteARF(fitsfile *fptr, struct ARF *arf)
{
  int tfields;
  int status=0;
  char *ttype[] = {"ENERG_LO", "ENERG_HI", "SPECRESP"};
  char *tform[] = {"E", "E", "E"};
  char *tunit[] = {"keV", "keV", "cm^2"};

  tfields = 3;

  /* Create the new extension */

  fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, tunit, arf->ARFExtensionName, &status);
  if (status) {
     headas_chat(1, "***Failed to create %s extension\n", arf->ARFExtensionName);
    return(status);
  } else {
     headas_chat(5, "Created %s extension\n", arf->ARFExtensionName);
  }

  /* Write the standard keywords */

  SP_write_key(fptr, TSTRING, "HDUCLASS", "OGIP", NULL);
    
  SP_write_key(fptr, TSTRING, "HDUCLAS1", "RESPONSE", NULL);

  SP_write_key(fptr, TSTRING, "HDUCLAS2", "SPECRESP", NULL);
    
  SP_write_key(fptr, TSTRING, "HDUVERS", arf->ARFVersion, "OGIP version number");

  SP_write_key(fptr, TSTRING, "TELESCOP", arf->Telescope, NULL);

  SP_write_key(fptr, TSTRING, "INSTRUME", arf->Instrument, NULL);

  SP_write_key(fptr, TSTRING, "DETNAM", arf->Detector, NULL);

  SP_write_key(fptr, TSTRING, "FILTER", arf->Filter, NULL);

  /* Write the start and stop energies for the bins and the effective area */

  fits_write_col(fptr, TFLOAT, 1, 1, 1, arf->NumberEnergyBins, arf->LowEnergy, &status);  
  fits_write_col(fptr, TFLOAT, 2, 1, 1, arf->NumberEnergyBins, arf->HighEnergy, &status);
  if (status) {
     headas_chat(1, "Failed to write ENERG_LO or ENERG_HI columns\n");
    return(status);
  }

  fits_write_col(fptr, TFLOAT, 3, 1, 1, arf->NumberEnergyBins, arf->EffArea, &status);
  if (status) {
     headas_chat(1, "Failed to write SPECRESP column\n");
    return(status);
  }

  return(status);

}



