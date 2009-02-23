// This program the output of the "tle_output"-routine in the orbit library.
// "tle_output" generates a tle-file entry according to given orbit data.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tle.h"
#include "earth_constants.h"


// Function to calculate the checksum for a specified output line.
int calc_checksum(const char* output);

// Main function
int main(void) {
  // constants
  const char full_name[24] = "eROSITA";     // catalog/common name of the satellite
  int epoch_year = 2011;              // year of epoch
  int epoch_time_s = 86400;           // time since the beginning of the year [s]
  double i = 30.0;                    // inclination [0 - 180) [degrees]
  double Omega = 45.0;                // right ascension of ascending node [0 - 360) [degrees]
  double e = 0.01;                    // eccentricity of the orbit (0.0 - 1.0)
  double omega = 35.0;                // argument of perigee [0 - 360) [degrees]
  double M = 0.0;                     // mean anomaly [0 - 360) [degrees]
  double a = 6958.140;                // semimajor axis

  double epoch_time;                  // time since the beginning of the year [days]
  double n;                           // mean motion [revolutions per day]


  // calculate the day of the epoch from the number of seconds
  epoch_time = (double)epoch_time_s/(24.*3600.);
  // calculate mean motion in revolutions per day
  n = sqrt(mu/pow(a,3.))*(24.*3600)/(2.*M_PI);

  // generate tle-entry
  char outstring[300];
  tle_output(outstring, full_name, epoch_year, epoch_time, i, Omega, e, omega, M, n);
  printf("\nOutput:\n%s", outstring);



  // test the tle_input-routine with the outstring, that was just generated
  parse_tle_element(outstring, &epoch_year, &epoch_time, &i, &Omega, &e, &omega, &M, &n);
  printf("\n\nInput:\nepoch year: %d\nepoch time [days]: %lf\ninclination: %lf\nright ascension of ascending node: %lf\neccentricity: %lf\nargument of perigee %lf\nmean anomaly: %lf\nmean motion [revolutions per day]: %lf\n\n", epoch_year, epoch_time, i, Omega, e, omega, M, n);

  return(EXIT_SUCCESS);
}
