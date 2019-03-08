/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#include "balancing.h"

BalancingArray* newBalancingArray(int* const status)
{
  BalancingArray* balance=(BalancingArray*)malloc(sizeof(BalancingArray));
  if (NULL==balance) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Balancing Array!\n",
		   *status);
    return(balance);
  }

  //Initialization:
  balance->Bmap=NULL;

  balance->naxis1 = 0;
  balance->naxis2 = 0;

  return(balance);
}


BalancingArray* getBalancingArray(ReconArray* recon, SquarePixels* detector_pixels,
				  EventList* ef,  int* const status)
{

  BalancingArray* balance=NULL;
  int x,y;                          //count for memory allocation
  int xcount, ycount;               //count for getting Rmap in case1:same pixel size


  //Get empty balancing array-object
  balance=newBalancingArray(status);
  if (EXIT_SUCCESS!=*status) return(balance);

  balance->naxis1=recon->naxis1;
  balance->naxis2=recon->naxis2;

  //memory-allocation for Bmap
  balance->Bmap=(double**)malloc(balance->naxis1*sizeof(double*));
  if(NULL!=balance->Bmap){
    for(x=0; x < balance->naxis1; x++){
      balance->Bmap[x]=(double*)malloc(balance->naxis2*sizeof(double));
	if(NULL==balance->Bmap[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "BalancingArray!\n", *status);
	  return(balance);
	}
	//Clear the pixels
	for(y=0; y < balance->naxis2; y++){
	  balance->Bmap[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(balance);
  } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "BalancingArray!\n", *status);
      return(balance);
  }//end of memory-allocation


  for(xcount=0; xcount < balance->naxis1; xcount++){
    for(ycount=0; ycount < balance->naxis2; ycount++){
      if(recon->Rmap[xcount][ycount]!=1)
	{
	  balance->Bmap[xcount][ycount]=(recon->Rmap[xcount][ycount]/(double)(detector_pixels->xwidth*
									  detector_pixels->ywidth))*(double)(ef->nrows);
	}
    }
  }
  return(balance);
}


double* SaveBalancingArray1d(BalancingArray* balance, int* status)
{
 double* BalancingArray1d=NULL;

 //Memory-Allocation for 1d-image
 BalancingArray1d = (double*)malloc(balance->naxis1*balance->naxis2*sizeof(double));
 if (!BalancingArray1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-balancing-array!\n", *status);
    return(BalancingArray1d);
 }

    //Create the 1D-image from BalancingArray
  int x, y;
  for (x=0; x<balance->naxis1; x++) {
    for (y=0; y<balance->naxis2; y++) {
	BalancingArray1d[(x+ balance->naxis1*y)] = balance->Bmap[x][y];
   }
  }
 return(BalancingArray1d);
}




void FreeBalancingArray(BalancingArray* balance)
{
  if (balance!=NULL) {
    if ((balance->naxis1>0)&&(NULL!=balance->Bmap)) {
      int count;
      for(count=0; count< balance->naxis1; count++) {
	if (NULL!=balance->Bmap[count]) {
	  free(balance->Bmap[count]);
	}
      }
      free(balance->Bmap);
    }
    free(balance);
  }
}


void FreeBalancingArray1d(double* BalancingArray1d)
{
  if (BalancingArray1d!=NULL) free(BalancingArray1d);
}
