/* Functions to write the RMF file */

#include <string.h>
#include <stdio.h>
#include "heasp.h"
#include "heaspio.h"

/* Write the MATRIX extension */

int WriteRMFMatrix(fitsfile *fptr, struct RMF *rmf)
{
  int tfields, i, j;
  long nelt;
  int status=0;
  char *ttype[] = {"ENERG_LO", "ENERG_HI", "N_GRP", "F_CHAN", 
		   "N_CHAN", "MATRIX", "ORDER"};
  char *tform[] = {"E", "E", "I", "PI", "PI", "PE", "PI"};
  char *tunit[] = {"keV", "keV", " ", " ", " ", " ", " "};

  tfields = 6;
  if (rmf->isOrder) tfields++;

  /* Create the new extension */

  fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, tunit, rmf->RMFExtensionName, &status);
  if (status) {
     headas_chat(1, "***Failed to create %s extension\n", rmf->RMFExtensionName);
    return(status);
  } else {
     headas_chat(5, "Created %s extension\n", rmf->RMFExtensionName);
  }

  /* Write the standard keywords */

  SP_write_key(fptr, TSTRING, "HDUCLASS", "OGIP", NULL);
    
  SP_write_key(fptr, TSTRING, "HDUCLAS1", "RESPONSE", NULL);

  SP_write_key(fptr, TSTRING, "HDUCLAS2", "RSP_MATRIX", NULL);
    
  SP_write_key(fptr, TSTRING, "HDUCLAS3", rmf->RMFType, "Type of response matrix");
    
  SP_write_key(fptr, TSTRING, "CHANTYPE", rmf->ChannelType, "Channel type");

  SP_write_key(fptr, TSTRING, "HDUVERS", rmf->RMFVersion, "OGIP version number");

  SP_write_key(fptr, TSTRING, "TELESCOP", rmf->Telescope, NULL);

  SP_write_key(fptr, TSTRING, "INSTRUME", rmf->Instrument, NULL);

  SP_write_key(fptr, TSTRING, "DETNAM", rmf->Detector, NULL);

  SP_write_key(fptr, TSTRING, "FILTER", rmf->Filter, NULL);

  SP_write_key(fptr, TFLOAT, "EFFAREA", &rmf->AreaScaling, "Area scaling factor"); 

  SP_write_key(fptr, TFLOAT, "LO_THRES", &rmf->ResponseThreshold, "Threshold for response value");

  SP_write_key(fptr, TLONG, "DETCHANS", &rmf->NumberChannels, "Number of channels in spectrum");

  SP_write_key(fptr, TLONG, "NUMGRP", &rmf->NumberTotalGroups, "Total number of response groups");

  SP_write_key(fptr, TLONG, "NUMELT", &rmf->NumberTotalElements, "Total number of response elements");

  SP_write_key(fptr, TLONG, "TLMIN4", &rmf->FirstChannel, "First channel number");

  /* Write the start and stop energies for the bins */

  fits_write_col(fptr, TFLOAT, 1, 1, 1, rmf->NumberEnergyBins, rmf->LowEnergy, &status);  
  fits_write_col(fptr, TFLOAT, 2, 1, 1, rmf->NumberEnergyBins, rmf->HighEnergy, &status);
  if (status) {
     headas_chat(1, "Failed to write ENERG_LO or ENERG_HI columns\n");
    return(status);
  }

  /* Write the number of groups for each energy bin */

  fits_write_col(fptr, TLONG, 3, 1, 1, rmf->NumberEnergyBins, rmf->NumberGroups, &status);
  if (status) {
     headas_chat(1, "Failed to write N_GRP column\n");
    return(status);
  }

  /* Loop round energy bins writing out the remaining columns - since these are variable
     length we have to do them row by row */

  for (i=0; i<rmf->NumberEnergyBins; i++) {

    fits_write_col(fptr, TLONG, 4, i+1, 1, rmf->NumberGroups[i], &rmf->FirstChannelGroup[rmf->FirstGroup[i]], &status);

    fits_write_col(fptr, TLONG, 5, i+1, 1, rmf->NumberGroups[i], &rmf->NumberChannelGroups[rmf->FirstGroup[i]], &status);

    nelt = 0;
    for (j=0; j<rmf->NumberGroups[i]; j++) nelt += rmf->NumberChannelGroups[j+rmf->FirstGroup[i]];
    fits_write_col(fptr, TFLOAT, 6, i+1, 1, nelt, &rmf->Matrix[rmf->FirstElement[rmf->FirstGroup[i]]], &status);

    if (rmf->isOrder) fits_write_col(fptr, TLONG, 7, i+1, 1, rmf->NumberGroups[i], &rmf->OrderGroup[rmf->FirstGroup[i]], &status);

    if (status) {
       headas_chat(1, "Failed to write N_CHAN, F_CHAN, MATRIX, or ORDER column for row %d\n", i+1);
      return(status);
    }

  }

  return(status);

}

/****************************** WriteRMFEbounds ********************************/

/* Write the EBOUNDS extension */

int WriteRMFEbounds(fitsfile *fptr, struct RMF *rmf)
{

  int i, status=0;
  long* Channel;
  char *ttype[] = {"CHANNEL", "E_MIN", "E_MAX"};
  char *tform[] = {"I", "E", "E"};
  char *tunit[] = {" ", "keV", "keV"};

  /* Create the new extension */

  fits_create_tbl(fptr, BINARY_TBL, 0, 3, ttype, tform, tunit, rmf->EBDExtensionName, &status);
  if (status) {
     headas_chat(1, "***Failed to create %s extension\n", rmf->EBDExtensionName);
    return(status);
  } else {
     headas_chat(5, "Created %s extension\n", rmf->RMFExtensionName);
  }

  /* Write the standard keywords */

  SP_write_key(fptr, TSTRING, "HDUCLASS", "OGIP", NULL);
    
  SP_write_key(fptr, TSTRING, "HDUCLAS1", "RESPONSE", NULL);

  SP_write_key(fptr, TSTRING, "HDUCLAS2", "EBOUNDS", NULL);
    
  SP_write_key(fptr, TSTRING, "CHANTYPE", rmf->ChannelType, "Channel type");

  SP_write_key(fptr, TSTRING, "HDUVERS", rmf->EBDVersion, "OGIP version");

  SP_write_key(fptr, TSTRING, "TELESCOP", rmf->Telescope, NULL);

  SP_write_key(fptr, TSTRING, "INSTRUME", rmf->Instrument, NULL);

  SP_write_key(fptr, TSTRING, "DETNAM", rmf->Detector, NULL);

  SP_write_key(fptr, TSTRING, "FILTER", rmf->Filter, NULL);

  SP_write_key(fptr, TLONG, "DETCHANS", &rmf->NumberChannels, "Number of channels in spectrum");

  /* Generate the CHANNEL array and write it */

  Channel = (long *) malloc(rmf->NumberChannels*sizeof(long));
  for (i=0; i<rmf->NumberChannels; i++) Channel[i] = rmf->FirstChannel + i;
  fits_write_col(fptr, TLONG, 1, 1, 1, rmf->NumberChannels, Channel, &status);
  free(Channel);
  if (status) {
     headas_chat(1, "Failed to write CHANNEL column\n");
    return(status);
  }

  /* Write the start and stop energies for the channels */

  fits_write_col(fptr, TFLOAT, 2, 1, 1, rmf->NumberChannels, rmf->ChannelLowEnergy, &status);  
  fits_write_col(fptr, TFLOAT, 3, 1, 1, rmf->NumberChannels, rmf->ChannelHighEnergy, &status);
  if (status) {
     headas_chat(1, "Failed to write E_MIN or E_MAX columns\n");
    return(status);
  }

  return(status);

}


