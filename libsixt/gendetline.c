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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#include "gendetline.h"


GenDetLine* newGenDetLine(const int xwidth, int* const status)
{
  GenDetLine* line=NULL;

  // Check if the requested line length is positive.
  if (xwidth<=0) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("invalid line length for GenDetLine");
    return(line);
  }

  line=(GenDetLine*)malloc(sizeof(GenDetLine));
  if (NULL==line) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenDetLine failed");
    return(line);
  }

  // Set all pointers to NULL.
  line->charge  =NULL;
  line->ccarry  =NULL;
  line->deadtime=NULL;
  line->ph_id   =NULL;
  line->src_id  =NULL;
  line->carry_ph_id   =NULL;
  line->carry_src_id  =NULL;

  line->xwidth=0;
  line->anycharge=0;
  line->anycarry=0;

  // Allocate memory.
  line->charge=(float*)malloc(xwidth*sizeof(float));
  CHECK_NULL_RET(line->charge, *status,
		 "memory allocation for GenDetLine failed", line);
  line->ccarry=(float*)malloc(xwidth*sizeof(float));
  CHECK_NULL_RET(line->ccarry, *status,
		 "memory allocation for GenDetLine failed", line);
  line->deadtime=(double*)malloc(xwidth*sizeof(double));
  CHECK_NULL_RET(line->deadtime, *status,
		 "memory allocation for GenDetLine failed", line);
  line->ph_id=(long**)malloc(xwidth*sizeof(long*));
  CHECK_NULL_RET(line->ph_id, *status,
		 "memory allocation for GenDetLine failed", line);
  line->src_id=(long**)malloc(xwidth*sizeof(long*));
  CHECK_NULL_RET(line->src_id, *status,
		 "memory allocation for GenDetLine failed", line);
  line->carry_ph_id=(long**)malloc(xwidth*sizeof(long*));
  CHECK_NULL_RET(line->carry_ph_id, *status,
		 "memory allocation for GenDetLine failed", line);
  line->carry_src_id=(long**)malloc(xwidth*sizeof(long*));
  CHECK_NULL_RET(line->carry_src_id, *status,
		 "memory allocation for GenDetLine failed", line);
  int ii;
  for(ii=0; ii<xwidth; ii++) {
    line->ph_id[ii]=(long*)malloc(NEVENTPHOTONS*sizeof(long));
    CHECK_NULL_RET(line->ph_id[ii], *status,
		   "memory allocation for GenDetLine failed", line);
    line->src_id[ii]=(long*)malloc(NEVENTPHOTONS*sizeof(long));
    CHECK_NULL_RET(line->src_id[ii], *status,
		   "memory allocation for GenDetLine failed", line);
    line->carry_ph_id[ii]=(long*)malloc(NEVENTPHOTONS*sizeof(long));
    CHECK_NULL_RET(line->carry_ph_id[ii], *status,
		   "memory allocation for GenDetLine failed", line);
    line->carry_src_id[ii]=(long*)malloc(NEVENTPHOTONS*sizeof(long));
    CHECK_NULL_RET(line->carry_src_id[ii], *status,
		   "memory allocation for GenDetLine failed", line);
    // Initialize.
    int jj;
    for (jj=0; jj<NEVENTPHOTONS; jj++) {
      line->ph_id[ii][jj] =0;
      line->src_id[ii][jj]=0;
      line->carry_ph_id[ii][jj] =0;
      line->carry_src_id[ii][jj]=0;
    }
  }

  line->xwidth=xwidth;

  // Clear the pixels.
  for(ii=0; ii<line->xwidth; ii++) {
    line->charge[ii]=0.;
    line->ccarry[ii]=0.;
    line->deadtime[ii]=0.;
    int jj;
    for (jj=0; jj<NEVENTPHOTONS; jj++) {
      line->ph_id[ii][jj] =0;
      line->src_id[ii][jj]=0;
      line->carry_ph_id[ii][jj] =0;
      line->carry_src_id[ii][jj]=0;
    }
  }

  return(line);
}


void destroyGenDetLine(GenDetLine** const line)
{
  if (NULL!=*line) {
    if (NULL!=(*line)->charge) {
      free((*line)->charge);
    }
    if (NULL!=(*line)->ccarry) {
      free((*line)->ccarry);
    }
    if (NULL!=(*line)->deadtime) {
      free((*line)->deadtime);
    }
    if (NULL!=(*line)->ph_id) {
      int ii;
      for (ii=0; ii<(*line)->xwidth; ii++) {
	if (NULL!=(*line)->ph_id[ii]) {
	  free((*line)->ph_id[ii]);
	}
      }
      free((*line)->ph_id);
    }
    if (NULL!=(*line)->src_id) {
      int ii;
      for (ii=0; ii<(*line)->xwidth; ii++) {
	if (NULL!=(*line)->src_id[ii]) {
	  free((*line)->src_id[ii]);
	}
      }
      free((*line)->src_id);
    }
    if (NULL!=(*line)->carry_ph_id) {
      int ii;
      for (ii=0; ii<(*line)->xwidth; ii++) {
	if (NULL!=(*line)->carry_ph_id[ii]) {
	  free((*line)->carry_ph_id[ii]);
	}
      }
      free((*line)->carry_ph_id);
    }
    if (NULL!=(*line)->carry_src_id) {
      int ii;
      for (ii=0; ii<(*line)->xwidth; ii++) {
	if (NULL!=(*line)->carry_src_id[ii]) {
	  free((*line)->carry_src_id[ii]);
	}
      }
      free((*line)->carry_src_id);
    }
    free(*line);
    *line=NULL;
  }
}


void clearGenDetLine(GenDetLine* const line)
{
  // Check if the line contains any charge. If not the clearing
  // is not necessary.
  if (1==line->anycharge) {
    int ii;
    for(ii=0; ii<line->xwidth; ii++) {
      if (line->charge[ii]>0.) {
	line->charge[ii]=0.;
	int jj;
	for (jj=0; jj<NEVENTPHOTONS; jj++) {
	  line->ph_id[ii][jj] =0;
	  line->src_id[ii][jj]=0;
	}
      }
    }
    line->anycharge=0;
  }
  // Check if the line contains carry charges for the next
  // read-out cycle.
  if (1==line->anycarry) {
    int ii;
    for(ii=0; ii<line->xwidth; ii++) {
      if (line->ccarry[ii]>0.) {
	line->charge[ii]=line->ccarry[ii];
	line->ccarry[ii]=0.;
	int jj;
	for (jj=0; jj<NEVENTPHOTONS; jj++) {
	  line->ph_id[ii][jj] =line->carry_ph_id[ii][jj];
	  line->src_id[ii][jj]=line->carry_src_id[ii][jj];
	  line->carry_ph_id[ii][jj] =0;
	  line->carry_src_id[ii][jj]=0;
	}
      }
    }
    line->anycarry=0;
  }
}


void addGenDetLine(GenDetLine* const line0, const GenDetLine* const line1)
{
  // Check if the 1st line, which is added to the 0th, contains any charges.
  if (0==line1->anycharge) return;

  // Set the anycharge flag.
  line0->anycharge=1;

  // Add the charges.
  int ii;
  for(ii=0; ii<line1->xwidth; ii++) {
    if (line1->charge[ii]>0.) {
      line0->charge[ii]+=line1->charge[ii];

      // Copy the photon and source IDs.
      int jj, kk;
      for (jj=0, kk=0; jj<NEVENTPHOTONS; jj++) {
	if (0==line1->ph_id[ii][kk]) break;
	if (0==line0->ph_id[ii][jj]) {
	  line0->ph_id[ii][jj] =line1->ph_id[ii][kk];
	  line0->src_id[ii][jj]=line1->src_id[ii][kk];
	  kk++;
	}
      }
    }
  }
}

