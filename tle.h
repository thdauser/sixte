//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// general information:
// file:             tle.h (header file of 'tle.c')
// project:          eROSITA - NRTA - simulation
// author:           Christian Schmid
// date:             01/2008
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// description:
// This code implements some useful routines to access the orbit data contained in
// NORAD TLE-files (see tle.c).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef TLE_H
#define TLE_H 1

// TLE access routines

// This function generates a tle-orbit-file entry according to given orbit data.
void tle_output(char *outstring, const char *full_name, const int epoch_year, const double epoch_time, const double e, const double i, const double Omega, const double omega, const double M, const double n);

// This function reads a NORAD TLE-element (1+2 lines), and extracts satellite data/orbit parameters:
int parse_tle_element(char* instring, int* epoch_year, double* epoch_day, double* e, double* i, double* Omega, double* omega, double* M, double* n);


#endif
