#ifndef FITS_PHA_H
#define FITS_PHA_H 1

#include <malloc.h>
#include <stdlib.h>

#include "fitsio.h"


#define PHA_NFIELDS 3
#define SPECTRUM_NFIELDS 2


////////////////////////////////////////////////////////////////
// Spectra

// creates the necessary data for the FITS-table layout
void spectrum_create_tbl_parameter(char *ftype[PHA_NFIELDS], char *fform[PHA_NFIELDS], char *funit[PHA_NFIELDS]);

// writes detector response data to a binary FITS table
int insert_spectrum_fitsrow(long channel, float p, fitsfile *fptr, long row);

// read spectrum line from a FITS table (PHA channel)
int read_spec_fitsrow(long *channel, float *probability, fitsfile *fptr, long row);



////////////////////////////////////////////////////////////////
// Detector redistribution matrix (RMF)

// creates the necessary data for the FITS-table layout
void rmf_create_tbl_parameter(char *ftype[PHA_NFIELDS], char *fform[PHA_NFIELDS], char *funit[PHA_NFIELDS]);

// writes detector redistribution data to a binary FITS table
int insert_rmf_fitsrow(long channel, float Emin, float Emax, fitsfile *fptr, long row);

// read detector redistribution data line from a FITS table (one energy channel)
int read_rmf_fitsrow(float *Emin, float *Emax, int *Ngrp, int *Fchan, int *Nchan, 
		     float *matrix, fitsfile *fptr, long row);

// read EBOUNDS data line from a FITS table (one PHA channel)
int read_ebounds_fitsrow(long *channel, float *Emin, float *Emax, fitsfile *fptr, long row);


#endif

