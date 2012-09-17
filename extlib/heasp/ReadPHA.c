/* Functions to read the PHA file */

#include <string.h>
#include <stdio.h>
#include "heasp.h"
#include "heaspio.h"

/******************************** ReturnPHAtype **************************************/

/* Read the SPECTRUM and return the type (1 or 2). Assumes that the file has 
   been opened but not positioned at the correct extension */

int ReturnPHAtype(fitsfile *fptr, long PHAnumber)
{
  int status=0;
  int colnum=0;
  int typecode;
  long repeat, width;
  char instring[FLEN_KEYWORD];

  /* Move to the correct SPECTRUM extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "SPECTRUM", PHAnumber, &status);
  if (!status) {
    headas_chat(5, "Found SPECTRUM extension...\n");
  } else {
    headas_chat(1, "***Cannot find SPECTRUM extension...\n");
    return(status);
  }

  /* Check for a HDUCLAS4 keyword of "TYPEII" */

  SP_read_key(fptr, TSTRING, "HDUCLAS4", instring, "TYPEI");
  if (!strcmp(instring, "TYPEII")) return(2); 

  /* Check for a type II extension by looking for the SPEC_NUM column */

  fits_get_colnum(fptr, CASEINSEN, "SPEC_NUM", &colnum, &status);
  if (!status && colnum != 0) return(2);

  /* If there is no SPEC_NUM column then look for a RATE or COUNTS column
     and check the TFORM */

  status = 0;
  fits_clear_errmsg();
  fits_get_colnum(fptr, CASEINSEN, "RATE", &colnum, &status);
  if (status) {
    status = 0;
    fits_clear_errmsg();
    fits_get_colnum(fptr, CASEINSEN, "COUNTS", &colnum, &status);
  }

  if (status) {
    headas_chat(1, "***Cannot find RATE or COUNTS column in SPECTRUM extension...\n");
    return(status);
  }

  fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
  if (!status && repeat > 1) return(2);

  /* None of the checks for a type II extension panned out so it must be a type I */

  return(1);

}

/******************************** CheckPHAcounts **************************************/

/* Read the SPECTRUM and if the data is COUNTS check that it is integer. Assumes that 
   the file has been opened but not positioned at the correct extension */

int CheckPHAcounts(fitsfile *fptr, long PHAnumber)
{
  int status=0;
  int colnum=0;
  int typecode;
  long repeat, width;
  char instring[FLEN_KEYWORD];

  /* Move to the correct SPECTRUM extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "SPECTRUM", PHAnumber, &status);
  if (!status) {
    headas_chat(5, "Found SPECTRUM extension...\n");
  } else {
    headas_chat(1, "***Cannot find SPECTRUM extension...\n");
    return(status);
  }

  /* Look for a COUNTS column and get the type */

  status = 0;
  fits_clear_errmsg();
  if (fits_get_colnum(fptr, CASEINSEN, "COUNTS", &colnum, &status)) {
    status = 0;
    fits_clear_errmsg();
  } else {
    if (fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status)) {
      return(status);
    } else {
      if ( typecode != TSHORT && typecode != TLONG ) {
	headas_chat(1, "***COUNTS data is non-integer\n");
	return(1);
      }
    }
  }

  return(0);

}

/******************************** ReturnNumberofSpectra **************************************/

/* Read the SPECTRUM and return the number of spectra in the (assumed) type 2 extension. 
   Assumes that the file has been opened but not positioned at the correct extension */

long ReturnNumberofSpectra(fitsfile *fptr, long PHAnumber)
{
  int status=0;
  long repeat;

  /* Move to the correct SPECTRUM extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "SPECTRUM", PHAnumber, &status);
  if (!status) {
    headas_chat(5, "Found SPECTRUM extension...\n");
  } else {
    headas_chat(1, "***ReturnNumberofSpectra: Cannot find SPECTRUM extension...\n");
    return(status);
  }

  /* Read the number of rows in the table */

  SP_read_key(fptr, TLONG, "NAXIS2", &repeat, 0);

  return(repeat);

}

/******************************** ReadPHAtypeI **************************************/

/* Read the Type I SPECTRUM extension and return the result in the PHA structure
   Assumes that the file has been opened but not positioned at the correct
   extension */

int ReadPHAtypeI(fitsfile *fptr, long PHAnumber, struct PHA *pha)
{
  int status=0;
  int colnum=0;
  int nfound=0;
  int i, itemp, ntotal;
  long ltemp, nrows;
  long start, stop;
  long *channels;
  float rtemp;
  char *ctemp1[10];
  char *ctemp2[90];

  /* Move to the correct SPECTRUM extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "SPECTRUM", PHAnumber, &status);
  if (!status) {
    headas_chat(5, "Found SPECTRUM extension...\n");
  } else {
    headas_chat(1, "***Cannot find SPECTRUM extension...\n");
    return(status);
  }

  /* Read the standard keywords and save the values */

  SP_read_key(fptr, TSTRING, "CHANTYPE", pha->ChannelType, "UNKNOWN");

  if (SP_read_key(fptr, TSTRING, "HDUVERS", pha->PHAVersion, "UNKNOWN"))
    SP_read_key(fptr, TSTRING, "HDUVERS1", pha->PHAVersion, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "TELESCOP", pha->Telescope, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "INSTRUME", pha->Instrument, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "DETNAM", pha->Detector, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "FILTER", pha->Filter, "UNKNOWN");

  SP_read_key(fptr, TSTRING, "DATAMODE", pha->Datamode, "UNKNOWN");

  ltemp = 0;
  status = SP_read_key(fptr, TLONG, "DETCHANS", &pha->NumberChannels, &ltemp);
  if (status) return(status);

  SP_read_key(fptr, TSTRING, "HDUCLAS2", pha->Spectrumtype, "TOTAL");

  SP_read_key(fptr, TSTRING, "HDUCLAS3", pha->Datatype, "COUNT");

  SP_read_key(fptr, TSTRING, "RESPFILE", pha->ResponseFile, "NONE");

  SP_read_key(fptr, TSTRING, "ANCRFILE", pha->AncillaryFile, "NONE");

  SP_read_key(fptr, TSTRING, "BACKFILE", pha->BackgroundFile, "NONE");

  SP_read_key(fptr, TSTRING, "CORRFILE", pha->CorrectionFile, "NONE");

  rtemp = 1.0;
  SP_read_key(fptr, TFLOAT, "CORRSCAL", &pha->CorrectionScaling, &rtemp);

  rtemp = 0.0;
  SP_read_key(fptr, TFLOAT, "EXPOSURE", &pha->Exposure, &rtemp);

  itemp = 0;
  SP_read_key(fptr, TINT, "POISSERR", &pha->Poisserr, &itemp);

  /* Read the XFLT keywords  - this is a bit messy because it can't be done in 
     a single cfitsio call */

  for (i=0; i<9; i++) ctemp1[i] = (char *) malloc(FLEN_KEYWORD*sizeof(char));
  fits_read_keys_str(fptr, "XFLT000", 0, 10, ctemp1, &nfound, &status);
  ntotal = nfound;
  if (nfound == 10) {
    for (i=0; i<90; i++) ctemp2[i] = (char *) malloc(FLEN_KEYWORD*sizeof(char));
    fits_read_keys_str(fptr, "XFLT00", 10, 100, ctemp2, &nfound, &status);
    ntotal += nfound;
  }
  for (i=0; i<ntotal; i++) pha->XSPECFilter[i] = (char *) malloc(FLEN_KEYWORD*sizeof(char));
  for (i=0; i<(ntotal<10?ntotal:10); i++) strcpy(pha->XSPECFilter[i], ctemp1[i]);
  if (ntotal > 10) for (i=10; i<(100>ntotal?100:ntotal); i++) strcpy(pha->XSPECFilter[i], ctemp2[i-10]);
  for (i=0; i<9; i++) free(ctemp1[i]);
  if (ntotal > 10) for (i=0; i<90; i++) free(ctemp2[i]);

  /* Check for TLMIN set for the CHANNEL column */

  pha->FirstChannel = 1;
  fits_get_colnum(fptr, CASEINSEN, "CHANNEL", &colnum, &status);
  if (!status) fits_read_keys_lng(fptr, "TLMIN", colnum, 1, &pha->FirstChannel, &nfound, &status);
  status = 0;
  fits_clear_errmsg();

  /* Read NAXIS2 to get the actual number of rows */

  SP_read_key(fptr, TLONG, "NAXIS2", &nrows, &ltemp);

  /* Read the data */

   headas_chat(5, "Allocating %ld for Pha\n", pha->NumberChannels);
  pha->Pha = (float *) malloc(pha->NumberChannels*sizeof(float));
  if (!strcmp(pha->Datatype,"RATE")) {
    SP_read_col(fptr, TFLOAT, "RATE", nrows, pha->Pha);
  } else {
    SP_read_col(fptr, TFLOAT, "COUNTS", nrows, pha->Pha);
  }

  /* Read the statistical error if poisserr is false */

  if (!pha->Poisserr) {
     headas_chat(5, "Allocating %ld for StatError\n", pha->NumberChannels);
    pha->StatError = (float *) malloc(pha->NumberChannels*sizeof(float));
    SP_read_col(fptr, TFLOAT, "STAT_ERR", nrows, pha->StatError);
  }

  /* Read the systematic error if it is there otherwise set the array to 0 */

   headas_chat(5, "Allocating %ld for SysError\n", pha->NumberChannels);
  pha->SysError = (float *) malloc(pha->NumberChannels*sizeof(float));
  if(SP_read_col(fptr, TFLOAT, "SYS_ERR", nrows, pha->SysError)) {
    for (i=0; i<nrows; i++) pha->SysError[i] = 0.0;
  }

  /* Read the QUALITY - first check for a keyword then for a column */

   headas_chat(5, "Allocating %ld for Quality\n", pha->NumberChannels);
  pha->Quality = (int *) malloc(pha->NumberChannels*sizeof(int));
  SP_read_col(fptr, TINT, "QUALITY", nrows, pha->Quality);
  
  /* Read the GROUPING - first check for a keyword then for a column */

   headas_chat(5, "Allocating %ld for Grouping\n", pha->NumberChannels);
  pha->Grouping = (int *) malloc(pha->NumberChannels*sizeof(int));
  SP_read_col(fptr, TINT, "GROUPING", nrows, pha->Grouping);

  /* Read the AREASCAL - first check for a keyword then for a column */

   headas_chat(5, "Allocating %ld for AreaScaling\n", pha->NumberChannels);
  pha->AreaScaling = (float *) malloc(pha->NumberChannels*sizeof(float));
  SP_read_col(fptr, TFLOAT, "AREASCAL", nrows, pha->AreaScaling);

  /* Read the BACKSCAL - first check for a keyword then for a column */

   headas_chat(5, "Allocating %ld for BackScaling\n", pha->NumberChannels);
  pha->BackScaling = (float *) malloc(pha->NumberChannels*sizeof(float));
  SP_read_col(fptr, TFLOAT, "BACKSCAL", nrows, pha->BackScaling);

  /* If the number of rows is not the same as DETCHANS then we need to pad out
     some values. We need to read the channels column to work out which elements
     to pad */
 
  if (nrows != pha->NumberChannels) {
    channels = (long *) malloc(nrows*sizeof(long));
    SP_read_col(fptr, TLONG, "CHANNEL", nrows, channels);
    start = channels[0] - pha->FirstChannel;
    stop  = channels[nrows-1] - pha->FirstChannel;
    for (i=nrows-1; i>=0; i--) {
      pha->Pha[i+start] = pha->Pha[i];
      if(!pha->Poisserr) pha->StatError[i+start] = pha->StatError[i];
      pha->SysError[i+start] = pha->SysError[i];
      pha->Quality[i+start] = pha->Quality[i];
      pha->Grouping[i+start] = pha->Grouping[i];
      pha->AreaScaling[i+start] = pha->AreaScaling[i];
      pha->BackScaling[i+start] = pha->BackScaling[i];
    }
    for (i=0; i<start; i++) {
      pha->Pha[i] = 0.0;
      if(!pha->Poisserr) pha->StatError[i] = 0.0;
      pha->SysError[i] = 0.0;
      pha->Quality[i] = 1;
      pha->Grouping[i] = 1;
      pha->AreaScaling[i] = pha->AreaScaling[start];
      pha->BackScaling[i] = pha->BackScaling[start];
    }
    for (i=stop+1; i<pha->NumberChannels; i++) { 
      pha->Pha[i] = 0.0;
      if(!pha->Poisserr) pha->StatError[i] = 0.0;
      pha->SysError[i] = 0.0;
      pha->Quality[i] = 1;
      pha->Grouping[i] = 1;
      pha->AreaScaling[i] = pha->AreaScaling[start];
      pha->BackScaling[i] = pha->BackScaling[start];
    }
  }


  return(status);

}

/******************************** ReadPHAtypeII **************************************/

/* Read the Type II SPECTRUM extension and return the result in the PHA structure
   Assumes that the file has been opened but not positioned at the correct
   extension */

int ReadPHAtypeII(fitsfile *fptr, long PHAnumber, long NumberSpectra, long *SpectrumNumber, struct PHA **pha)
{
  int status=0;
  int nfound=0;
  int i, j, itemp, colnum;
  long ispec, ltemp, start, stop, nrows;
  long *channels;
  float rtemp;

  /* Move to the correct SPECTRUM extension */

  fits_movnam_hdu(fptr, BINARY_TBL, "SPECTRUM", PHAnumber, &status);
  if (!status) {
    headas_chat(5, "Found SPECTRUM extension...\n");
  } else {
    headas_chat(1, "***Cannot find SPECTRUM extension...\n");
    return(status);
  }

  /* loop round the spectra to read */

  for (i=0; i<NumberSpectra; i++) {

    ispec = SpectrumNumber[i];

  /* Read the standard keywords and save the values */

    SPII_read_key(fptr, TSTRING, "CHANTYPE", ispec, pha[i]->ChannelType, "UNKNOWN");

    if (SPII_read_key(fptr, TSTRING, "HDUVERS", ispec, pha[i]->PHAVersion, "UNKNOWN"))
      SPII_read_key(fptr, TSTRING, "HDUVERS1", ispec, pha[i]->PHAVersion, "UNKNOWN");

    SPII_read_key(fptr, TSTRING, "TELESCOP", ispec, pha[i]->Telescope, "UNKNOWN");

    SPII_read_key(fptr, TSTRING, "INSTRUME", ispec, pha[i]->Instrument, "UNKNOWN");

    SPII_read_key(fptr, TSTRING, "DETNAM", ispec, pha[i]->Detector, "UNKNOWN");

    SPII_read_key(fptr, TSTRING, "FILTER", ispec, pha[i]->Filter, "UNKNOWN");

    SPII_read_key(fptr, TSTRING, "DATAMODE", ispec, pha[i]->Datamode, "UNKNOWN");

    ltemp = 0;
    status = SPII_read_key(fptr, TLONG, "DETCHANS", ispec, &pha[i]->NumberChannels, &ltemp);
    if (status) return(status);

    SPII_read_key(fptr, TSTRING, "HDUCLAS2", ispec, pha[i]->Spectrumtype, "TOTAL");

    SPII_read_key(fptr, TSTRING, "HDUCLAS3", ispec, pha[i]->Datatype, "COUNT");

    SPII_read_key(fptr, TSTRING, "RESPFILE", ispec, pha[i]->ResponseFile, "NONE");

    SPII_read_key(fptr, TSTRING, "ANCRFILE", ispec, pha[i]->AncillaryFile, "NONE");

    SPII_read_key(fptr, TSTRING, "BACKFILE", ispec, pha[i]->BackgroundFile, "NONE");

    SPII_read_key(fptr, TSTRING, "CORRFILE", ispec, pha[i]->CorrectionFile, "NONE");

    rtemp = 1.0;
    SPII_read_key(fptr, TFLOAT, "CORRSCAL", ispec, &pha[i]->CorrectionScaling, &rtemp);

    rtemp = 0.0;
    SPII_read_key(fptr, TFLOAT, "EXPOSURE", ispec, &pha[i]->Exposure, &rtemp);

    itemp = 0;
    SPII_read_key(fptr, TINT, "POISSERR", ispec, &pha[i]->Poisserr, &itemp);

  /* Check for TLMIN set for the CHANNEL column */

    pha[i]->FirstChannel = 1;
    fits_get_colnum(fptr, CASEINSEN, "CHANNEL", &colnum, &status);
    if (!status) fits_read_keys_lng(fptr, "TLMIN", colnum, 1, &pha[i]->FirstChannel, &nfound, &status);
    status = 0;
    fits_clear_errmsg();

  /* Get the actual number of elements in the spectrum vector and set the Datatype in case
     the optional HDUCLAS3 was not used */

    fits_get_colnum(fptr, CASEINSEN, "RATE", &colnum, &status);
    if (status) {
      status = 0;
      fits_clear_errmsg();
      fits_get_colnum(fptr, CASEINSEN, "COUNTS", &colnum, &status);
      if (status) {
	headas_chat(1, "***Failed to find RATE or COUNTS column\n");
	return(status);
      }
      strcpy(pha[i]->Datatype,"COUNT");
       headas_chat(5, "Datatype is COUNT\n");
    } else {
      strcpy(pha[i]->Datatype,"RATE");
       headas_chat(5, "Datatype is RATE\n");
    }
    fits_read_descript(fptr, colnum, ispec, &nrows, &ltemp, &status);
    if (status) {
      status = 0;
      fits_clear_errmsg();
      fits_get_coltype(fptr, colnum, &itemp, &nrows, &ltemp, &status);
      if (status) {
	headas_chat(1, "***Failed to read descriptor for column %d, status = %d\n", colnum, status);
	return(status);
      }
    }
     headas_chat(5, "Column number %d for spectrum %ld has %ld channels\n", colnum, ispec, nrows);

  /* Read the data */

     headas_chat(5, "Allocating %ld for Pha\n", pha[i]->NumberChannels);
    pha[i]->Pha = (float *) malloc(pha[i]->NumberChannels*sizeof(float));
    if (!strcmp(pha[i]->Datatype,"RATE")) {
      SPII_read_col(fptr, TFLOAT, "RATE", ispec, nrows, pha[i]->Pha);
    } else {
      SPII_read_col(fptr, TFLOAT, "COUNTS", ispec, nrows, pha[i]->Pha);
    }

  /* Read the statistical error if poisserr is false */

    if (!pha[i]->Poisserr) {
       headas_chat(5, "Allocating %ld for StatError\n", pha[i]->NumberChannels);
      pha[i]->StatError = (float *) malloc(pha[i]->NumberChannels*sizeof(float));
      SPII_read_col(fptr, TFLOAT, "STAT_ERR", ispec, nrows, pha[i]->StatError);
    }

  /* Read the systematic error if it is there otherwise set the array to 0 */

     headas_chat(5, "Allocating %ld for SysError\n", pha[i]->NumberChannels);
    pha[i]->SysError = (float *) malloc(pha[i]->NumberChannels*sizeof(float));
    if(SPII_read_col(fptr, TFLOAT, "SYS_ERR", ispec, nrows, pha[i]->SysError)) {
      for (i=0; i<nrows; i++) pha[i]->SysError[i] = 0.0;
    }

  /* Read the QUALITY - first check for a keyword then for a column */

     headas_chat(5, "Allocating %ld for Quality\n", pha[i]->NumberChannels);
    pha[i]->Quality = (int *) malloc(pha[i]->NumberChannels*sizeof(int));
    SPII_read_col(fptr, TINT, "QUALITY", ispec, nrows, pha[i]->Quality);
  
  /* Read the GROUPING - first check for a keyword then for a column */

     headas_chat(5, "Allocating %ld for Grouping\n", pha[i]->NumberChannels);
    pha[i]->Grouping = (int *) malloc(pha[i]->NumberChannels*sizeof(int));
    SPII_read_col(fptr, TINT, "GROUPING", ispec, nrows, pha[i]->Grouping);

  /* Read the AREASCAL - first check for a keyword then for a column */

     headas_chat(5, "Allocating %ld for AreaScaling\n", pha[i]->NumberChannels);
    pha[i]->AreaScaling = (float *) malloc(pha[i]->NumberChannels*sizeof(float));
    SPII_read_col(fptr, TFLOAT, "AREASCAL", ispec, nrows, pha[i]->AreaScaling);

  /* Read the BACKSCAL - first check for a keyword then for a column */

     headas_chat(5, "Allocating %ld for BackScaling\n", pha[i]->NumberChannels);
    pha[i]->BackScaling = (float *) malloc(pha[i]->NumberChannels*sizeof(float));
    SPII_read_col(fptr, TFLOAT, "BACKSCAL", ispec, nrows, pha[i]->BackScaling);

  /* If the number of rows is not the same as DETCHANS then we need to pad out
     some values. We need to read the channels column to work out which elements
     to pad */
 
    if (nrows != pha[i]->NumberChannels) {
      channels = (long *) malloc(nrows*sizeof(long));
      SPII_read_col(fptr, TLONG, "CHANNEL", ispec, nrows, channels);
      start = channels[0] - pha[i]->FirstChannel;
      stop  = channels[nrows-1] - pha[i]->FirstChannel;
      for (j=nrows-1; j>=0; j--) {
	pha[i]->Pha[j+start] = pha[i]->Pha[j];
	if(!pha[i]->Poisserr) pha[i]->StatError[j+start] = pha[i]->StatError[j];
	pha[i]->SysError[j+start] = pha[i]->SysError[j];
	pha[i]->Quality[j+start] = pha[i]->Quality[j];
	pha[i]->Grouping[j+start] = pha[i]->Grouping[j];
	pha[i]->AreaScaling[j+start] = pha[i]->AreaScaling[j];
	pha[i]->BackScaling[j+start] = pha[i]->BackScaling[j];
      }
      for (j=0; j<start; j++) {
	pha[i]->Pha[j] = 0.0;
	if(!pha[i]->Poisserr) pha[i]->StatError[j] = 0.0;
	pha[i]->SysError[j] = 0.0;
	pha[i]->Quality[j] = 1;
	pha[i]->Grouping[j] = 1;
	pha[i]->AreaScaling[j] = 1.0;
	pha[i]->BackScaling[j] = 1.0;
      }
      for (j=stop+1; j<pha[i]->NumberChannels; j++) { 
	pha[i]->Pha[j] = 0.0;
	if(!pha[i]->Poisserr) pha[i]->StatError[j] = 0.0;
	pha[i]->SysError[j] = 0.0;
	pha[i]->Quality[j] = 1;
	pha[i]->Grouping[j] = 1;
	pha[i]->AreaScaling[j] = 1.0;
	pha[i]->BackScaling[j] = 1.0;
      }
    }

    /* end of loop over spectra to read */

  }


  return(status);

}

