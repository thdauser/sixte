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


   Copyright 2014 Thorsten Brand, FAU
*/

#ifndef TESPROFTEMPLFILE_H
#define TESPROFTEMPLFILE_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Structure containing the calorimeter profile templates of one extension. */
typedef struct{

  /** Number of time steps. */
  long Nt;
  
  /** Number of Energy steps. */
  int NE;
  
  /** time array. */
  double *time;
  
  /** Energy array. */
  double *energy;
  
  /** Signal array. */
  double **adc_value;  
  
}TESProfilesEntries;

/** Structure containing the calorimeter profile templates of one file. */
typedef struct{
  
  /** Number of versions (pixels etc.). */
  int Nv;
  
  /** version string array. */
  char **version;
  
  /** Array containing the profiles. */
  TESProfilesEntries *profiles;
  
}TESProfiles;

/** Structure containing the input values for pulse profile template generation */
typedef struct{
  
  /** Version name */
  char **version;
  
  /** Sample Frequency (Hz) */
  double freq;
  
  /** Number of energy steps */
  int ne;
  
  /** Array with energies (keV) */
  double *energies;
  
  /** Number of samples per template (2^x, advice >=4096) */
  long nsamp;
  
  /** Pulse rise time (usually equals tfall/50 us) */
  double trise;
  
  /** Pulse fall time (typically 200 us = 2E-4 seconds) */
  double tfall;
  
}TESTemplateInput;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Function which initializes the data structure TESProfilesEntries. */
void newTESProfilesEntries(TESProfilesEntries *prof);

/** Function which initializes the data structure TESProfiles. */
TESProfiles* newTESProfiles(int* const status);

/** Destructor for the data structure TESProfilesEntries. */
void destroyTESProfilesEntries(TESProfilesEntries* prof);

/** Destructor for the data structure TESProfiles. */
void destroyTESProfiles(TESProfiles* prof);

/** Function which returns 1 if a version is already loaded, otherwise 0*/
int testTESProfilesExist(TESProfiles *prof, char *version);

/** Function which reads one FITS table with pulse profile templates. */
void readTESProfiles(char *filename, 
			     char *version, 
			     TESProfiles *prof, 
			     int* const status);

/** Function which creates a pulse profile template file. */
int createTESProfilesFile(char *filename, 
			   const char clobber,
			   char *comment,
			   int* const status);

/** Function which inserts a pulse profile template in a table pointed by fptr */
int InsertTESProfADCCol(fitsfile *fptr, 
			long nt, 
			double energy, 
			double *adc, 
			int* const status);

/** Function which writes one FITS table with pulse profile templates. */
int writeTESProfiles(char *filename, 
			     char *version, 
			     TESProfilesEntries *prof, 
			     const char clobber,
			     char *comment,
			     int* const status);

/** Function which looks for a specific version and returns the index 
    or -1 if the version is not yet present in the template collection */
int findTESProfileVersionIndex(TESProfiles* prof,
			       char *version);

/** Function which returns the index of the best-matching energy
    index in the proper version array. */
int findTESProfileEnergyIndex(TESProfiles* prof, 
			      int version, 
			      double energy);

/** Generate pulse profiles function */
int genTESProfile(TESTemplateInput* pinp, TESProfiles** ptemp, int* const status);

/** Function to calculate a bare exponential pulse */
double ExponentialPulse(double *t, double *trise, double *tfall);

#endif /* TESPROFTEMPLFILE_H */