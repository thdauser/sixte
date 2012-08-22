/* Functions to write the PHA file */

#include <string.h>
#include <stdio.h>
#include "heasp.h"
#include "heaspio.h"

/******************************** WritePHAtypeI **************************************/

/* Write the type I PHA extension */

int WritePHAtypeI(fitsfile *fptr, struct PHA *pha)
{
  int tfields, colnum, i, n;
  void *pnt;
  long *channels;
  int status=0;
  char *ttype[8];
  char *tform[8];
  char *tunit[8];

  /* set up the column descriptors for those attributes which need to be written as
     columns */

  for (i=0; i<9; i++) {
    ttype[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
    tform[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
    tunit[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
  }

  strcpy(ttype[0], "CHANNEL");
  strcpy(tform[0], "J");
  strcpy(tunit[0], " ");

  if (!strcmp(pha->Datatype, "RATE")) {
    strcpy(ttype[1], "RATE");
    strcpy(tform[1], "E");
    strcpy(tunit[1], "counts/s");
  } else {
    strcpy(ttype[1], "COUNTS");
    strcpy(tform[1], "J");
    strcpy(tunit[1], "counts");
  }
     
  tfields = 1;
  if (!pha->Poisserr) {
    tfields++;
    strcpy(ttype[tfields], "STAT_ERR");
    strcpy(tform[tfields], "E");
    strcpy(tunit[tfields], " ");
  }

  if (SP_need_col(pha->SysError, pha->NumberChannels, TFLOAT)) {
    tfields++;
    strcpy(ttype[tfields], "SYS_ERR");
    strcpy(tform[tfields], "E");
    strcpy(tunit[tfields], " ");
  }

  if (SP_need_col(pha->Quality, pha->NumberChannels, TINT)) {
    tfields++;
    strcpy(ttype[tfields], "QUALITY");
    strcpy(tform[tfields], "I");
    strcpy(tunit[tfields], " ");
  }

  if (SP_need_col(pha->Grouping, pha->NumberChannels, TINT)) {
    tfields++;
    strcpy(ttype[tfields], "GROUPING");
    strcpy(tform[tfields], "I");
    strcpy(tunit[tfields], " ");
  }

  if (SP_need_col(pha->AreaScaling, pha->NumberChannels, TFLOAT)) {
    tfields++;
    strcpy(ttype[tfields], "AREASCAL");
    strcpy(tform[tfields], "E");
    strcpy(tunit[tfields], " ");
  }

  if (SP_need_col(pha->BackScaling, pha->NumberChannels, TFLOAT)) {
    tfields++;
    strcpy(ttype[tfields], "BACKSCAL");
    strcpy(tform[tfields], "E");
    strcpy(tunit[tfields], " ");
  }

  tfields++;

  /* Create the new extension */

  fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, tunit, "SPECTRUM", &status);
  if (status) {
     headas_chat(1, "***Failed to create SPECTRUM extension\n");
    return(status);
  } else {
     headas_chat(5, "Created SPECTRUM extension\n");
  }

  /* Write the standard keywords */

  SP_write_key(fptr, TSTRING, "HDUCLASS", "OGIP", NULL);
    
  SP_write_key(fptr, TSTRING, "HDUCLAS1", "SPECTRUM", NULL);

  SP_write_key(fptr, TSTRING, "HDUCLAS2", pha->Spectrumtype, NULL);
    
  SP_write_key(fptr, TSTRING, "HDUCLAS3", pha->Datatype, NULL);
    
  SP_write_key(fptr, TSTRING, "CHANTYPE", pha->ChannelType, "Channel type");

  SP_write_key(fptr, TSTRING, "HDUVERS", pha->PHAVersion, "OGIP version number");

  SP_write_key(fptr, TSTRING, "TELESCOP", pha->Telescope, NULL);

  SP_write_key(fptr, TSTRING, "INSTRUME", pha->Instrument, NULL);

  SP_write_key(fptr, TSTRING, "DETNAM", pha->Detector, NULL);

  SP_write_key(fptr, TSTRING, "FILTER", pha->Filter, NULL);

  SP_write_key(fptr, TSTRING, "DATAMODE", pha->Datamode, NULL);

  SP_write_key(fptr, TLONG, "DETCHANS", &pha->NumberChannels, "Number of channels in spectrum");

  SP_write_key(fptr, TLONG, "TLMIN1", &pha->FirstChannel, "First channel number");

  SP_write_key(fptr, TFLOAT, "EXPOSURE", &pha->Exposure, "Exposure time");

  SP_write_key(fptr, TFLOAT, "CORRSCAL", &pha->CorrectionScaling, "Scaling for correction file");

  SP_write_key(fptr, TLOGICAL, "POISSERR", &pha->Poisserr, "Is error Poisson ?");

  SP_write_key(fptr, TSTRING, "RESPFILE", pha->ResponseFile, NULL);

  SP_write_key(fptr, TSTRING, "ANCRFILE", pha->AncillaryFile, NULL);

  SP_write_key(fptr, TSTRING, "BACKFILE", pha->BackgroundFile, NULL);

  SP_write_key(fptr, TSTRING, "CORRFILE", pha->CorrectionFile, NULL);

  n = sizeof(pha->XSPECFilter)/FLEN_FILENAME/sizeof(char);
  fits_write_keys_str(fptr, "XFLT000", 1, (n>9?n:9), pha->XSPECFilter, NULL, &status);
  pnt = pha->XSPECFilter + 10*FLEN_FILENAME*sizeof(char);
  if ( n > 9 ) fits_write_keys_str(fptr, "XFLT00", 10, (n>99?n:99), pnt, NULL, &status);

  /* Make the array to write to the CHANNEL column */

  channels = (long *) malloc(pha->NumberChannels*sizeof(long));
  for (i=0; i<pha->NumberChannels; i++) channels[i] = i + pha->FirstChannel;

  /* Write the arrays - if an array is all the same value it will be written
     as a keyword */

  colnum = 1;
  SP_write_col(fptr, TLONG, "CHANNEL", &colnum, pha->NumberChannels, channels);

  if (!strcmp(pha->Datatype, "RATE")) {
    SP_write_col(fptr, TFLOAT, "RATE", &colnum, pha->NumberChannels, pha->Pha);
  } else {
    SP_write_col(fptr, TFLOAT, "COUNTS", &colnum, pha->NumberChannels, pha->Pha);
  }

  if (!pha->Poisserr) SP_write_col(fptr, TFLOAT, "STAT_ERR", &colnum, pha->NumberChannels, pha->StatError);

  SP_write_col(fptr, TFLOAT, "SYS_ERR", &colnum, pha->NumberChannels, pha->SysError);

  SP_write_col(fptr, TINT, "QUALITY", &colnum, pha->NumberChannels, pha->Quality);

  SP_write_col(fptr, TINT, "GROUPING", &colnum, pha->NumberChannels, pha->Grouping);

  SP_write_col(fptr, TFLOAT, "AREASCAL", &colnum, pha->NumberChannels, pha->AreaScaling);

  SP_write_col(fptr, TFLOAT, "BACKSCAL", &colnum, pha->NumberChannels, pha->BackScaling);

  return(status);

}

/******************************** WritePHAtypeII **************************************/

int WritePHAtypeII(fitsfile *fptr, long NumberSpectra, struct PHA **pha)
{
  int tfields, i, needcol, isvector, ispec;
  long MaxElements;
  long *channels;
  int status=0;
  char *ttype[23];
  char *tform[23];
  char *tunit[23];
  char **sarray;

  /* set up the column descriptors for those attributes which need to be written as
     columns */

  for (i=0; i<23; i++) {
    ttype[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
    tform[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
    tunit[i] = (char *) malloc(FLEN_VALUE*sizeof(char));
  }

  /* Hopefully all the spectra we are writing have the same number of channels but just
     in case find the maximum then load the channels array */

  MaxElements = pha[0]->NumberChannels;
  for (ispec=1; ispec<NumberSpectra; ispec++) 
    if (pha[ispec]->NumberChannels > MaxElements) MaxElements = pha[ispec]->NumberChannels;
  channels = (long *) malloc(MaxElements*sizeof(long));
  for (i=0; i<MaxElements; i++) channels[i] = i+1;

  /* We know we need SPEC_NUM, CHANNEL, and RATE/COUNTS columns */

  strcpy(ttype[0], "SPEC_NUM");
  strcpy(tform[0], "J");
  strcpy(tunit[0], " ");

  strcpy(ttype[1], "CHANNEL");
  sprintf(tform[1], "%ldJ", MaxElements);
  strcpy(tunit[1], " ");

  if (!strcmp(pha[0]->Datatype, "RATE")) {
    strcpy(ttype[2], "RATE");
    sprintf(tform[2], "%ldE", MaxElements); 
   strcpy(tunit[2], "counts/s");
  } else {
    strcpy(ttype[2], "COUNTS");
    sprintf(tform[2], "%ldE", MaxElements);
    strcpy(tunit[2], "counts");
  }
  tfields = 2;

  /* Find which quantities need to be in columns and which can stay in keywords. To be in
     a keyword the quantity must be the same for all spectra and all channels. We
     also need to find which columns can be scalar and which need to be vector */

  needcol = 0;
  for (ispec=0; ispec<NumberSpectra; ispec++) needcol = needcol || !pha[ispec]->Poisserr;
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "STAT_ERR");
    sprintf(tform[tfields], "%ldE", MaxElements);
    strcpy(tunit[tfields], " ");
  }

  isvector = 0;
  for (ispec=0; ispec<NumberSpectra; ispec++) 
    isvector = isvector || SP_need_col(pha[ispec]->SysError, pha[ispec]->NumberChannels, TFLOAT);
  needcol = isvector;
  if (!needcol) {
    for (ispec=1; ispec<NumberSpectra; ispec++)
      if (pha[ispec]->SysError[0] != pha[0]->SysError[0]) needcol = 1;
  }
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "SYS_ERR");
    if (isvector) {
      sprintf(tform[tfields], "%ldE", MaxElements);
    } else {
      strcpy(tform[tfields], "E");
    }
    strcpy(tunit[tfields], " ");
  }

  isvector = 0;
  for (ispec=0; ispec<NumberSpectra; ispec++) 
    isvector = isvector || SP_need_col(pha[ispec]->Quality, pha[ispec]->NumberChannels, TINT);
  needcol = isvector;
  if (!needcol) {
    for (ispec=1; ispec<NumberSpectra; ispec++)
      if (pha[ispec]->Quality[0] != pha[0]->Quality[0]) needcol = 1;
  }
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "QUALITY");
    if (isvector) {
      sprintf(tform[tfields], "%ldI", MaxElements);
    } else {
      strcpy(tform[tfields], "I");
    }
    strcpy(tunit[tfields], " ");
  }

  isvector = 0;
  for (ispec=0; ispec<NumberSpectra; ispec++) 
    isvector = isvector || SP_need_col(pha[ispec]->Grouping, pha[ispec]->NumberChannels, TINT);
  needcol = isvector;
  if (!needcol) {
    for (ispec=1; ispec<NumberSpectra; ispec++)
      if (pha[ispec]->Grouping[0] != pha[0]->Grouping[0]) needcol = 1;
  }
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "GROUPING");
    if (isvector) {
      sprintf(tform[tfields], "%ldI", MaxElements);
    } else {
      strcpy(tform[tfields], "I");
    }
    strcpy(tunit[tfields], " ");
  }

  isvector = 0;
  for (ispec=0; ispec<NumberSpectra; ispec++) 
    isvector = isvector || SP_need_col(pha[ispec]->AreaScaling, pha[ispec]->NumberChannels, TFLOAT);
  needcol = isvector;
  if (!needcol) {
    for (ispec=1; ispec<NumberSpectra; ispec++)
      if (pha[ispec]->AreaScaling[0] != pha[0]->AreaScaling[0]) needcol = 1;
  }
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "AREASCAL");
    if (isvector) {
      sprintf(tform[tfields], "%ldE", MaxElements);
    } else {
      strcpy(tform[tfields], "E");
    }
    strcpy(tunit[tfields], " ");
  }

  isvector = 0;
  for (ispec=0; ispec<NumberSpectra; ispec++) 
    isvector = isvector || SP_need_col(pha[ispec]->BackScaling, pha[ispec]->NumberChannels, TFLOAT);
  needcol = isvector;
  if (!needcol) {
    for (ispec=1; ispec<NumberSpectra; ispec++)
      if (pha[ispec]->BackScaling[0] != pha[0]->BackScaling[0]) needcol = 1;
  }
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "BACKSCAL");
    if (isvector) {
      sprintf(tform[tfields], "%ldE", MaxElements);
    } else {
      strcpy(tform[tfields], "E");
    }
    strcpy(tunit[tfields], " ");
  }

  /* now we need to find if any of the standard keywords differ between the spectra so must
     be placed in a column. Quantities to check are Exposure, CorrectionScaling, Datatype,
     SpectrumType, ResponseFile, AncillaryFile, BackgroundFile, CorrectionFile, ChannelType, 
     Telescope, Instrument, Detector, Filter, Datamode. I should really check XSPECFilter as
     well but that would be too horrible.*/

  needcol = 0;
  for (ispec=1; ispec<NumberSpectra; ispec++)
    if (pha[ispec]->Exposure != pha[0]->Exposure) needcol = 1;
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "EXPOSURE");
    strcpy(tform[tfields], "E");
    strcpy(tunit[tfields], " ");
  }

  needcol = 0;
  for (ispec=1; ispec<NumberSpectra; ispec++)
    if (pha[ispec]->CorrectionScaling != pha[0]->CorrectionScaling) needcol = 1;
  if (needcol) {
    tfields++;
    strcpy(ttype[tfields], "CORRSCAL");
    strcpy(tform[tfields], "E");
    strcpy(tunit[tfields], " ");
  }

  /* we will need a temporary string array */

  sarray = (char **) malloc(NumberSpectra*sizeof(char *));
  for (ispec=0; ispec<NumberSpectra; ispec++)
    sarray[ispec] = (char *)malloc(FLEN_FILENAME*sizeof(char));

  /* now check all the character variables */
  
  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Datatype);
  SPII_set_table_str("HDUCLAS2", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Spectrumtype);
  SPII_set_table_str("HDUCLAS3", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->ResponseFile);
  SPII_set_table_str("RESPFILE", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->AncillaryFile);
  SPII_set_table_str("ANCRFILE", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->BackgroundFile);
  SPII_set_table_str("BACKFILE", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->CorrectionFile);
  SPII_set_table_str("CORRFILE", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->ChannelType);
  SPII_set_table_str("CHANTYPE", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Telescope);
  SPII_set_table_str("TELESCOP", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Instrument);
  SPII_set_table_str("INSTRUME", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Detector);
  SPII_set_table_str("DETNAM", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Filter);
  SPII_set_table_str("FILTER", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  for (ispec=0; ispec<NumberSpectra; ispec++)
    strcpy(sarray[ispec], pha[ispec]->Datamode);
  SPII_set_table_str("DATAMODE", sarray, NumberSpectra, &tfields, ttype, tform, tunit);

  /* free up the memory for the temporary array */

  free(sarray);

  /* create the binary table - increment tfield to change it from a pointer counting from 0 
     to a number counting from 1.*/

  tfields++;
  fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, tunit, "SPECTRUM", &status);
  if (status) {
     headas_chat(1, "***Failed to create SPECTRUM extension\n");
    return(status);
  } else {
    headas_chat(5, "Created SPECTRUM extension with columns\n");
    for (i=0;i<tfields;i++) headas_chat(5, " %s ", ttype[i]);
    headas_chat(5, "\n");
  }

  /* write the standard keywords that aren't included in the columns */

  SPII_write_key(fptr, tfields, ttype, TSTRING, "HDUCLASS", "OGIP", NULL);
    
  SPII_write_key(fptr, tfields, ttype, TSTRING, "HDUCLAS1", "SPECTRUM", NULL);

  SPII_write_key(fptr, tfields, ttype, TSTRING, "HDUCLAS2", pha[0]->Spectrumtype, NULL);
    
  SPII_write_key(fptr, tfields, ttype, TSTRING, "HDUCLAS3", pha[0]->Datatype, NULL);
    
  SPII_write_key(fptr, tfields, ttype, TSTRING, "CHANTYPE", pha[0]->ChannelType, "Channel type");

  SPII_write_key(fptr, tfields, ttype, TSTRING, "HDUVERS", pha[0]->PHAVersion, "OGIP version number");

  SPII_write_key(fptr, tfields, ttype, TSTRING, "TELESCOP", pha[0]->Telescope, NULL);

  SPII_write_key(fptr, tfields, ttype, TSTRING, "INSTRUME", pha[0]->Instrument, NULL);

  SPII_write_key(fptr, tfields, ttype, TSTRING, "DETNAM", pha[0]->Detector, NULL);

  SPII_write_key(fptr, tfields, ttype, TSTRING, "FILTER", pha[0]->Filter, NULL);

  SPII_write_key(fptr, tfields, ttype, TSTRING, "DATAMODE", pha[0]->Datamode, NULL);

  /* Loop over the spectra writing out the columns -
     the SPII_write_col routine checks whether the named column is in the ttype array
     and if so writes the column */

  for (ispec=1; ispec<=NumberSpectra; ispec++) {

    SPII_write_col(fptr, tfields, ttype, TINT, "SPEC_NUM", ispec, 1, &ispec);

    SPII_write_col(fptr, tfields, ttype, TLONG, "CHANNEL", ispec, pha[ispec-1]->NumberChannels, channels);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "RATE", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->Pha);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "COUNTS", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->Pha);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "STAT_ERR", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->StatError);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "SYS_ERR", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->SysError);

    SPII_write_col(fptr, tfields, ttype, TINT, "QUALITY", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->Quality);

    SPII_write_col(fptr, tfields, ttype, TINT, "GROUPING", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->Grouping);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "AREASCAL", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->AreaScaling);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "BACKSCAL", ispec, pha[ispec-1]->NumberChannels, pha[ispec-1]->BackScaling);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "EXPOSURE", ispec, 1, &pha[ispec-1]->Exposure);

    SPII_write_col(fptr, tfields, ttype, TFLOAT, "CORRSCAL", ispec, 1, &pha[ispec-1]->CorrectionScaling);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "HDUCLAS2", ispec, 1, pha[ispec-1]->Datatype);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "HDUCLAS3", ispec, 1, pha[ispec-1]->Spectrumtype);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "RESPFILE", ispec, 1, pha[ispec-1]->ResponseFile);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "ANCRFILE", ispec, 1, pha[ispec-1]->AncillaryFile);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "BACKFILE", ispec, 1, pha[ispec-1]->BackgroundFile);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "CORRFILE", ispec, 1, pha[ispec-1]->CorrectionFile);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "CHANTYPE", ispec, 1, pha[ispec-1]->ChannelType);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "TELESCOP", ispec, 1, pha[ispec-1]->Telescope);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "INSTRUME", ispec, 1, pha[ispec-1]->Instrument);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "DETNAM", ispec, 1, pha[ispec-1]->Detector); 

    SPII_write_col(fptr, tfields, ttype, TSTRING, "FILTER", ispec, 1, pha[ispec-1]->Filter);

    SPII_write_col(fptr, tfields, ttype, TSTRING, "DATAMODE", ispec, 1, pha[ispec-1]->Datamode);

  }

  return(0);

}


