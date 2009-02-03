///////////////////////////////////////
// general information:
// file:             tle.c
// related to:       orbit.h, liborbit.a
// project:          eROSITA - NRTA - simulation
// author:           Christian Schmid
// date:             01/2008
///////////////////////////////////////
// description:
// This file belongs to a library (liborbit.a) out of the NRTA-simulation for the eROSITA mission.
// It implements several functions to access NORAD TLE-files.
/////////////////////////////////////////////////////////////////////////
// TLE format:
// 0th line:
//   * full name of the satellite
//
// 1st line:
//   * NORAD number
//   * classification ('U' = unclassified)
//   * launch year
//   * launch number
//   * launch piece
//   * epoch: - year
//            - day (8 decimal number for the fraction of the day)
//
//   * checksum
//
// 2nd line:
//
//   * n - mean motion (in revolutions per day)
//
//   * checksum
//
//////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "strftcpy.h"

// Function to calculate the checksum for a specified TLE output line.
int calc_checksum(const char* output);


///////////////////////////////////////////////////////////////////////////////////////
// This function generates a tle-orbit-file entry according to given orbit data.
// Some orbit parameters are already predefined.
//
// @params:
// char[] - full name of the satllite
// int    - epoch year
// double - epoch time in days with decimal numbers
// double - eccentrictiy
// double - inclination
// double - right ascension of ascending node
// double - argument of perigee
// double - mean anomaly
// double - mean motion in revolutions per day
////////////////////////////////////////////////////////////////////////////////////////
void tle_output(char *outstring, const char *full_name, const int epoch_year, const double epoch_time, const double e, const double i, const double Omega, const double omega, const double M, const double n) {

  // constants
  const int norad_number = 1;         // NORAD number of the satellite
  const char classification = 'U';    // security classification (U = unclassified)
  const int launch_year = 2011;       // year of the launch
  const int launch_number = 1;        // launch-number in the specified year
  const char launch_piece[3] = "A";   // piece-number of the satellite (launcher is usually "B")
  const double sgp_drag_par1 = 0.0;   // SGP-model: \dot{n} * 1/2 (where n is the mean motion)
  const double sgp_drag_par2 = 0.0;   // \dot{\dot{n}} * 1/6
  const double bstar_drag_par = 0.0;  // Bstar drag parameter
  const int ephemeris_type = 0;       // ephermis type
  const int element_number = 1;       // element number (1 - 9999)
  const int revolution_number = 12;   // revolution number (1 - 99999)


  // variables:
  char output[70];      // output buffer
  char sgp_p2[12];
  char bstar_p[12];
  int checksum = 0;     // checksum


  // pre-formatting:
  // 2nd SGP drag parameter:
  sprintf(sgp_p2, "%.4e", sgp_drag_par2);

  // Bstar drag parameter:
  sprintf(bstar_p, "%.4e", bstar_drag_par);

  // printing TLE data set:
  // print the 0th line (containing the common name of the satellite - entirely 24 characters)
  sprintf(output, "%-24s\n", full_name);
  strcpy(outstring, output);

  // print the 1st line (generate output-string - then calculate checksum and plot it)
  sprintf(output, "1 %05d%c %02d%03d%-3s %02d%012.8lf  .%08d  %c%c%c%c%c-%c  %c%c%c%c%c-%c %1d %4d", norad_number, classification, (launch_year>=2000)?(launch_year-2000):(launch_year-1900), launch_number, launch_piece, (epoch_year>=2000)?(epoch_year-2000):(epoch_year-1900), epoch_time, ((int)(sgp_drag_par1*100000000))%100000000, sgp_p2[0] , sgp_p2[2], sgp_p2[3], sgp_p2[4], sgp_p2[5], (sgp_p2[9]==48)?48:sgp_p2[9]-1, bstar_p[0] , bstar_p[2], bstar_p[3], bstar_p[4], bstar_p[5], (bstar_p[9]==48)?48:bstar_p[9]-1, ephemeris_type, element_number%10000);
  strcat(outstring, output);     // append the formatted line to the output string
  sprintf(output, "%1d\n", calc_checksum(output));   // checksum & new line
  strcat(outstring, output);

  // print the 2nd line
  sprintf(output, "2 %05d %8.4lf %8.4lf %07d %8.4lf %8.4lf %11.8lf%05d", norad_number, i, Omega, ((int)(e*10000000))%10000000, omega, M, n, revolution_number%100000);
  strcat(outstring, output);     // append the formatted line to the output string
  sprintf(output, "%1d\n", calc_checksum(output));   // checksum & new line
  strcat(outstring, output);
}



///////////////////////////////////////////////////////////////////////////////////////////////
// This function reads a NORAD TLE-element (1+2 lines), and extracts interesting data, like: //
//    * epoch year and day (day with decimal numbers giving the fraction of the day)
//    * eccentricity
//    * inclination
//    * right ascension of ascending node
//    * argument of perigee
//    * mean anomaly
//    * mean motion in revolutions per day (!)
//
// The checksums are verified. If an error occurs, the return value is -1, successful parsing
// mean return value 0.
///////////////////////////////////////////////////////////////////////////////////////////////
int parse_tle_element(char* instring, int* epoch_year, double* epoch_day, double* e, double* i, double* Omega, double* omega, double* M, double* n) {
  char cbuffer[100];  // input-buffer
  int chksum;         // buffer for the checksum


  // parse TLE data:

  // 1st line:
  // epoch year:
  strftcpy(cbuffer, instring, 43, 2);
  *epoch_year = atoi(cbuffer);
  if (*epoch_year < 57) {
    *epoch_year += 2000;
  } else {
    *epoch_year += 1900;
  }

  // epoch time (seconds since the beginning of the year):
  strftcpy(cbuffer, instring, 45, 12);
  *epoch_day = atof(cbuffer);

  // verify the checksum of the 1st line
  strftcpy(cbuffer, instring, 25, 68);
  chksum = calc_checksum(cbuffer);
  strftcpy(cbuffer, instring, 93, 1);
  if (chksum != atoi(cbuffer)) {
    return(-1);
  }
  

  // 2nd line:
  // inclination
  strftcpy(cbuffer, instring, 103, 8);
  *i = atof(cbuffer);

  // Omega - right ascension of ascending node
  strftcpy(cbuffer,instring,112,8);
  *Omega = atof(cbuffer);

  // eccentricity
  strftcpy(cbuffer, instring, 121, 7);
  *e = atof(cbuffer)/10000000.0;

  // omega - argument of perigee
  strftcpy(cbuffer, instring, 129, 8);
  *omega = atof(cbuffer);

  // M - mean anomaly
  strftcpy(cbuffer, instring, 138, 8);
  *M = atof(cbuffer);

  // n - mean motion
  strftcpy(cbuffer, instring, 147, 11);
  *n = atof(cbuffer);

  // verify the checksum of the 2nd line
  strftcpy(cbuffer, instring, 95, 68);
  chksum = calc_checksum(cbuffer);
  strftcpy(cbuffer, instring, 163, 1);
  if (chksum != atoi(cbuffer)) {
    return(-1);
  }
  

  // parsing was successful (no errors):
  return(0);
}





// Function to calculate the checksum for a specified TLE output line.
int calc_checksum(const char* line) {
  int i;
  int checksum;

  for(i=0,checksum=0;i<68;i++){
    if (line[i]=='-') {     // "-"-sign
      checksum++;
    } else if ((line[i]>47)&&(line[i]<58)) {    // number 0-9
      checksum += line[i]-48;
    }
  }
  return(checksum%10);
}
