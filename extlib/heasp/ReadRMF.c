/* Functions to read the RMF file */

#include <string.h>
#include <stdio.h>
#include "heasp.h"
#include "heaspio.h"

/* Read the MATRIX extension and return the result in the RMF structure
   Assumes that the file has been opened but not positioned at the correct
   extension */

int ReadRMFMatrix(fitsfile *fptr, long RMFnumber, struct RMF *rmf)
{
  int status=0;
  int status2=0;
  int colnum=0;
  int anynul=0;
  int nfound=0;
  int i, j, ipt, igrp;
  long nelt;
  float rtemp;
  long ltemp;

  /* Move to the correct RMF extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "MATRIX", RMFnumber, &status);
  if (status) {
    status = 0;
    fits_clear_errmsg();
    fits_movnam_hdu(fptr, BINARY_TBL, "SPECRESP MATRIX", RMFnumber, &status);
  }
  if (!status) {
     headas_chat(5, "Found MATRIX or SPECRESP MATRIX extension...\n");
  } else {
     headas_chat(1, "***Cannot find MATRIX or SPECRESP MATRIX extension...\n");
     headas_chat(1, "   FITS status = %d\n", status);
    return(status);
  }

  /* Read the standard keywords and save the values */

  SP_read_key(fptr, TSTRING, "EXTNAME", rmf->RMFExtensionName, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "HDUCLAS3", rmf->RMFType, "UNKNOWN");
    
  SP_read_key(fptr, TSTRING, "CHANTYPE", rmf->ChannelType, "UNKNOWN");

  if (SP_read_key(fptr, TSTRING, "HDUVERS", rmf->RMFVersion, "UNKNOWN")) {
    if (SP_read_key(fptr, TSTRING, "HDUVERS2", rmf->RMFVersion, "UNKNOWN")) {
      if (!SP_read_key(fptr, TSTRING, "RMFVERSN", rmf->RMFVersion, "UNKNOWN")) {
	strcpy(rmf->RMFVersion,"1.0.0");
      }
    }
  }

  SP_read_key(fptr, TSTRING, "TELESCOP", rmf->Telescope, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "INSTRUME", rmf->Instrument, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "DETNAM", rmf->Detector, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "FILTER", rmf->Filter, "UNKNOWN");

  rtemp = 1.0;
  SP_read_key(fptr, TFLOAT, "EFFAREA", &rmf->AreaScaling, &rtemp); 

  rtemp = 0.0;
  SP_read_key(fptr, TFLOAT, "LO_THRES", &rmf->ResponseThreshold, &rtemp);

  ltemp = 0;
  status = SP_read_key(fptr, TLONG, "DETCHANS", &rmf->NumberChannels, &ltemp);
  if (status) return(status);

  status = SP_read_key(fptr, TLONG, "NAXIS2", &rmf->NumberEnergyBins, &ltemp);
  if (status) return(status);

  /* A couple of optional keywords - if they are present it saves us some calculation later */

  ltemp = 0;
  SP_read_key(fptr, TLONG, "NUMGRP", &rmf->NumberTotalGroups, &ltemp);
  SP_read_key(fptr, TLONG, "NUMELT", &rmf->NumberTotalElements, &ltemp);

  /* Get the start and stop energies for the bins */

  headas_chat(5, "Allocating %ld for energy bins\n", rmf->NumberEnergyBins);

  rmf->LowEnergy = (float *) malloc(rmf->NumberEnergyBins*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "ENERG_LO", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, rmf->NumberEnergyBins, NULL, rmf->LowEnergy, &anynul, &status);

  rmf->HighEnergy = (float *) malloc(rmf->NumberEnergyBins*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "ENERG_HI", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, rmf->NumberEnergyBins, NULL, rmf->HighEnergy, &anynul, &status);

  /* Get the number of groups for each energy bin */

  rmf->NumberGroups = (long *) malloc(rmf->NumberEnergyBins*sizeof(long));
  fits_get_colnum(fptr, CASEINSEN, "N_GRP", &colnum, &status);
  fits_read_col(fptr, TLONG, colnum, 1, 1, rmf->NumberEnergyBins, NULL, rmf->NumberGroups, &anynul, &status);

  if (status) {
    fits_report_error(stderr, status);
    printf("\n");
    return(status);
  }

  /* Set up the FirstGroup array - note counts from 0 */

  rmf->FirstGroup = (long *) malloc(rmf->NumberEnergyBins*sizeof(long));
  igrp = 0;
  for (i=0; i<rmf->NumberEnergyBins; i++) {
    rmf->FirstGroup[i] = igrp;
    igrp += rmf->NumberGroups[i];
  }

  /* If the NUMGRP keyword was not read then sum up this column to calculate it */

  if ( rmf->NumberTotalGroups == 0 ) {
    for (i=0; i<rmf->NumberEnergyBins; i++) {
      rmf->NumberTotalGroups += rmf->NumberGroups[i];
    }
     headas_chat(5, "Setting NumberTotalGroups to %ld\n", rmf->NumberTotalGroups);
  }

  /* Get the first channel for each group */

  headas_chat(5, "Allocating %ld for groups\n", rmf->NumberTotalGroups);

  rmf->FirstChannelGroup = (long *) malloc(rmf->NumberTotalGroups*sizeof(long));
  fits_get_colnum(fptr, CASEINSEN, "F_CHAN", &colnum, &status);
  ipt = 0;
  for (i=0; i<rmf->NumberEnergyBins; i++) {
    fits_read_col(fptr, TLONG, colnum, i+1, 1, rmf->NumberGroups[i], NULL, &rmf->FirstChannelGroup[ipt], &anynul, &status);
    ipt += rmf->NumberGroups[i];
  }

  fits_read_keys_lng(fptr, "TLMIN", colnum, 1, &rmf->FirstChannel, &nfound, &status);
  if (status || nfound==0) {
    rmf->FirstChannel = 1;
     headas_chat(5, "Failed to read TLMIN for F_CHAN column - setting FirstChannel to 1\n");
    status = 0;
    fits_clear_errmsg();
  } else {
     headas_chat(5, "TLMIN for F_CHAN column = %ld\n", rmf->FirstChannel);
  }

  /* Get the number of channels for each group */

  rmf->NumberChannelGroups = (long *) malloc(rmf->NumberTotalGroups*sizeof(long));
  fits_get_colnum(fptr, CASEINSEN, "N_CHAN", &colnum, &status);
  ipt = 0;
  for (i=0; i<rmf->NumberEnergyBins; i++) {
    fits_read_col(fptr, TLONG, colnum, i+1, 1, rmf->NumberGroups[i], NULL, &rmf->NumberChannelGroups[ipt], &anynul, &status);
    ipt += rmf->NumberGroups[i];
  }

  if (status) {
    fits_report_error(stderr, status);
    printf("\n");
    return(status);
  }

  /* If the NUMELT keyword was not read then sum up this column to calculate it */

  if ( rmf->NumberTotalElements == 0 ) {
    for (i=0; i<rmf->NumberTotalGroups; i++) {
      rmf->NumberTotalElements += rmf->NumberChannelGroups[i];
    }
     headas_chat(5, "Setting NumberTotalElements to %ld\n", rmf->NumberTotalElements);
  }

  /* Get the memory for the FirstElement array pointing to the first element in the Matrix
     array for each response group */

  rmf->FirstElement = (long *) malloc(rmf->NumberTotalGroups*sizeof(long));

  /* Read the matrix information */

  headas_chat(5, "Allocating %ld for response elements\n", rmf->NumberTotalElements);

  rmf->Matrix = (float *) malloc(rmf->NumberTotalElements*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "MATRIX", &colnum, &status);
  ipt = 0;
  igrp = 0;
  for (i=0; i<rmf->NumberEnergyBins; i++) {
    nelt = 0;
    for (j=0; j<rmf->NumberGroups[i]; j++) {
      rmf->FirstElement[igrp] = nelt+ipt;
      nelt += rmf->NumberChannelGroups[igrp];
      igrp++;
    }
    fits_read_col(fptr, TFLOAT, colnum, i+1, 1, nelt, NULL, &rmf->Matrix[ipt], &anynul, &status);
    ipt += nelt;
  }

  /* Read the optional order information */

  status2 = 0;
  fits_write_errmark();
  rmf->isOrder = 0;
  fits_get_colnum(fptr, CASEINSEN, "ORDER", &colnum, &status2);
  fits_clear_errmark();
  if (!status2) {
    rmf->isOrder = 1;
    rmf->OrderGroup = (long *) malloc(rmf->NumberTotalGroups*sizeof(long));
    fits_read_col(fptr, TLONG, colnum, 1, 1, rmf->NumberTotalGroups, NULL, rmf->OrderGroup, &anynul, &status);
  }

  /* Check for a cfitsio error and if it occurred write diagnostic information */

  if (status) {
    fits_report_error(stderr, status);
    printf("\n");
  }

  return(status);

}

/* Read the EBOUNDS extension and return the result in the RMF structure
   Assumes that the file has been opened but not positioned at the correct
   extension */

int ReadRMFEbounds(fitsfile *fptr, long EBDnumber, struct RMF *rmf)
{
  int status=0;
  int colnum=0;
  int anynul=0;
  long ltemp;

  /* Move to the correct RMF extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "EBOUNDS", EBDnumber, &status);

  headas_chat(5, "Found EBOUNDS extension...\n");

  /* Read the standard keywords and save the values */

  SP_read_key(fptr, TSTRING, "EXTNAME", rmf->EBDExtensionName, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "CHANTYPE", rmf->ChannelType, "UNKNOWN");

  if (SP_read_key(fptr, TSTRING, "HDUVERS", rmf->RMFVersion, "UNKNOWN")) {
    if (SP_read_key(fptr, TSTRING, "HDUVERS2", rmf->RMFVersion, "UNKNOWN")) {
      if (!SP_read_key(fptr, TSTRING, "RMFVERSN", rmf->RMFVersion, "UNKNOWN")) {
	strcpy(rmf->RMFVersion,"1.0.0");
      }
    }
  }

  SP_read_key(fptr, TSTRING, "TELESCOP", rmf->Telescope, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "INSTRUME", rmf->Instrument, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "FILTER", rmf->Filter, "UNKNOWN");

  ltemp = 0;
  SP_read_key(fptr, TLONG, "DETCHANS", &rmf->NumberChannels, &ltemp);
  if (status) return(status);

  /* Get the start and stop energies for the channels */

  headas_chat(5, "Allocating %ld for channels\n", rmf->NumberChannels);

  rmf->ChannelLowEnergy = (float *) malloc(rmf->NumberChannels*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "E_MIN", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, rmf->NumberChannels, NULL, rmf->ChannelLowEnergy, &anynul, &status);

  rmf->ChannelHighEnergy = (float *) malloc(rmf->NumberChannels*sizeof(float));
  fits_get_colnum(fptr, CASEINSEN, "E_MAX", &colnum, &status);
  fits_read_col(fptr, TFLOAT, colnum, 1, 1, rmf->NumberChannels, NULL, rmf->ChannelHighEnergy, &anynul, &status);

  /* Check for a cfitsio error and if it occurred write diagnostic information */

  if (status) {
    fits_report_error(stderr, status);
    printf("\n");
  }

  return(status);

}

