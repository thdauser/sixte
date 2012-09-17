
#include <string.h>
#include <stdio.h>
#include "heasp.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

/************************* ReadBinningFile **************************/

/* reads an ascii file of binning factors and loads the BinFactors structure */

int SPReadBinningFile(char *filename, struct BinFactors *binning)
{
  int i;
  FILE *fptr;
  char line[512];

  if ( (fptr = fopen(filename, "r")) == NULL ) {
    headas_chat(1,"***ReadBinningFile: Cannot open %s\n", filename);
    return 1;
  }

  /* run through the file getting the number of lines */

  binning->NumberBinFactors = 0;
  while ( fgets(line, 512, fptr) != NULL ) binning->NumberBinFactors++;
  rewind(fptr);

  /* allocate memory for the arrays */

  binning->StartBin = (long *) malloc(binning->NumberBinFactors*sizeof(long));
  binning->EndBin   = (long *) malloc(binning->NumberBinFactors*sizeof(long));
  binning->Binning  = (long *) malloc(binning->NumberBinFactors*sizeof(long));

  if ( binning->StartBin == NULL || binning->StartBin == NULL || binning->StartBin == NULL ) {
    headas_chat(1,"***ReadBinningFile: Cannot allocate memory for binning arrays\n", filename);
    return 2;
  }
    
  /* loop through file again reading lines and storing in arrays */

  for (i=0; i<binning->NumberBinFactors; i++) fscanf(fptr, "%ld %ld %ld", &binning->StartBin[i], &binning->EndBin[i], &binning->Binning[i]);

  return 0;

}

/************************* SPSetGroupArray **************************/

/* Sets up a grouping array using the information in the BinFactors structure
   Assume that output has been allocated enough memory
   Assume that the binning ranges are contiguous the start channel of
   a range is the channel after the end channel of the previous range
   Assume that the first channel is 1 */

int SPSetGroupArray(int inputSize, struct BinFactors *binning, int *grouping)

{

  int status=0;
  int ir, ich, ibin, nbins;


  /* set grouping for any elements before the start of binning */

  for (ich=1; ich<binning->StartBin[0]; ich++) grouping[ich-1] = 1;


  /* loop through ranges for binning */

  for (ir=0; ir<binning->NumberBinFactors; ir++) {

    if (binning->Binning[ir] == -1 ) {

      for (ich=binning->StartBin[ir]; ich<=binning->EndBin[ir]; ich++) grouping[ich-1] = -1;

    } else {

      for (ich=binning->StartBin[ir]; ich<=binning->EndBin[ir]; ich+=binning->Binning[ir]) {

	grouping[ich-1] = 1;
	nbins = MIN(binning->Binning[ir], inputSize-ich);
	for (ibin=1; ibin<nbins; ibin++) {
	  grouping[ich-1+ibin] = 0;
	}

      }

    }

  }

  /* set grouping for any elements after the end of binning */

  for (ich=binning->EndBin[binning->NumberBinFactors-1]+1; ich<=inputSize; ich++) grouping[ich-1] = 1;

  return(status);

}

/************************* SPBinArray **************************/

/* Bin an array based on a grouping array
   assumes that the output array has been malloc'ed correctly
   mode = 1   Sum
   mode = 2   Mean
   mode = 3   First ie value is that of first channel in bin
   mode = 4   Last  ie value is that of first channel in bin */

int SPBinArray(int inputSize, float *input, int *grouping, int mode, float *output)

{

  int status=0;
  int ich, ibin, n;

  ibin = -1;
  for (ich=0; ich<inputSize; ich++) {

    if (grouping[ich] == 1) {
      if ( ibin > -1 && mode == 2 ) output[ibin] /= n;
      ibin++;
      n = 1;
      output[ibin] = input[ich];
    } else if (grouping[ich] == 0) {
      if ( mode == 1 || mode == 2 ) output[ibin] += input[ich];
      if ( mode == 4 ) output[ibin] = input[ich];
      n++;
    }

  }

  if ( mode == 2 ) output[ibin] /= n;

  return(status);

}


