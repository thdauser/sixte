#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pil.h"
#include "vector.h"

// Function solves the Kepler-equation for given mean anomaly 'M' and
// given eccentricity 'e' with the Newton algorithm.
// returns the eccentric anomaly 'E'
double kepler_equation(double M, double e);


// Program calculates the orbit of a satellite according to given orbit parameters:
// -> a - semimajor axis
// -> e - eccentricity
// -> i - inclination
// -> Omega - right ascension of ascending node
// -> omega - argument of perigee
// -> M0    - mean anomaly at t0 (usually M0=0, i.e. the satellite is in perigee)
//
// output: location \vec{r}(t) and velocity \vec{v}(t)

int main(int argc, char** argv) {
  // general constants:
  const double mue = 398601.3;  //G*Msolar; units: (km^3/s^2)
  const double R_e = 6378.140;  // radius of the earth
//  const double satellite_height = 600.0;  // orbit height of the satellite (a=R_e+satellite_height)

        // physical constants for the earth (Flury, p. 69):
        //      mue = 398601.3 km^3/s^2
        //      M_e = 5.976 * 10^27 g
        //      R_e = 6378.140 km  (mean equatorial radius)


  // variables:
  double Epi = 0.0;         // eccentric anomaly (is zero at perigee)
  struct vector rpi, vpi;   // position of perigee and velocity of the satellite in perigee
  double Omega;             // right ascension of ascending node
  double omega;             // argument of perigee
  double M;                 // mean anomaly M(t)
  double E;                 // eccentric anomaly E(t)
  double Ft,Fpt,Gt,Gpt;     // F(t), d/dt F(t), G(t), d/dt G(t)
  struct vector r, v;       // resulting vectors \vec{r}(t) and \vec{v}(t)
  double rl;                // length of vector r
  double cosi, cosi2;       // cosine of inclination: cos(i) and (cos(i))^2

  double t;                 // time counter
  long plotstep = 0;        // plotstep counter
  char outputfile[100];     // filename of the output-file
  FILE *outfptr;            // filepointer to output-file


  // input parameters:
  double a;             // semimajor axis
  double e;             // eccentricity
  double i;             // inclination
  double Omega0;        // initial right ascension of ascending node
  double omega0;        // initial argument of perigee
  double M0;            // mean anomaly at t0=0 ( M0=0.0 => satellite is in perigee)
  double t0;            // time-offset
  double timespan;      // timespan in which the satellite data is calculated
  double dt;            // width of a timestep
  int Nplotintervals;   // orbit data is only plotted, if (t/dt)%Nplotintervals == 0

  // read parameters using PIL library
  // initialize PIL
  int pilread = PILInit(argc, argv);
  // read the length of the semimajor axis a
  pilread = PILGetReal("a", &a);
  // read the eccentricity of the Kepler orbit
  pilread = PILGetReal("e", &e);
  // read the length of the semimajor axis a
  pilread = PILGetReal("i", &i);
  // read the right ascension of the ascending node
  pilread = PILGetReal("gOmega0", &Omega0);
  // read the argument of perigee
  pilread = PILGetReal("komega0", &omega0);
  // read the mean anomaly @ t0
  pilread = PILGetReal("M0", &M0);
  // read the time offset t0
  pilread = PILGetReal("t0", &t0);
  // read the entire timespan of the calculation
  pilread = PILGetReal("timespan", &timespan);
  // read the width of a timestep
  pilread = PILGetReal("dt", &dt);
  // read the number of timesteps performed per data-output
  pilread = PILGetInt("Nplotintervals", &Nplotintervals);
  // get the filename of the output-file
  pilread = PILGetFname("outputfile", outputfile);
  // close PIL interface
  PILClose(PIL_OK);

  // check, whether all parameters were read correctly:
  if (pilread < 0) {
    printf("PIL - parameter input failed: %s\n", PIL_err_handler(pilread));
    return(-1);
  }

  // convert angles [degrees] -> [radians]:
  i = i*M_PI/180.0;
  Omega0 = Omega0*M_PI/180.0;
  omega0 = omega0*M_PI/180.0;
  M0 = M0*M_PI/180.0;

  // start at t0 with initial data:
  Omega = Omega0;
  omega = omega0;
  M = M0;

  // calculation of basic orbit data:
  // length of perigee distance:
  const double rpil = a*(1-e);
  // satellite velocity in perigee:
  const double vpil = sqrt((mue*(1+e))/(a*(1-e)));
  // calculate the squared length of the vector vpi:
  //const double vpil2 = pow(vpil,2.0);
  // calculate mean motion:
  const double n = sqrt(mue/pow(a,3.0));
  // calculate the factor s, that contributes to J2 orbit corrections:
  const double s = -108263*pow(10.0,-8.0)*1.5*n*pow(R_e/(a*(1.0-pow(e,2.0))),2.0); // Nitschke S. 25
  const double ds = dt*s;


  // Now we have the necessary data to calculate the position and velocity of
  // the satellite on its orbit at any given point of time.
  // So we can determine the position of the satellite at any point of time 't':


  // first delete the old output FITS file
  remove(outputfile);

  // open ASCII file for the output of the position and velocity of the satellite
  outfptr=fopen(outputfile, "w");

  if (outfptr != NULL) {
    // loop over the timeinterval from t0 to t0+timespan
    for (t=t0; t<t0+timespan; t+=dt) {
      // calculate the position of the perigee (Flury, p. 38):
      rpi.x = rpil*(cos(Omega)*cos(omega) - sin(Omega)*cos(i)*sin(omega));
      rpi.y = rpil*(sin(Omega)*cos(omega) + cos(Omega)*cos(i)*sin(omega));
      rpi.z = rpil*sin(i)*sin(omega);

      // get the velocity of the satellite in perigee:
      vpi.x = -vpil*(cos(Omega)*sin(omega) + sin(Omega)*cos(omega)*cos(i));
      vpi.y = -vpil*(sin(Omega)*sin(omega) - cos(Omega)*cos(omega)*cos(i));
      vpi.z = vpil*cos(omega)*sin(i);

      // calculate mean anomaly at time t (M(t)):
      M += n*dt;

      // get the eccentric anomaly (E(t)):
      E = kepler_equation(M,e);

      // Now we can determine the functions F(t), G(t) and their derivatives,
      // and therefore the vectors r(t) and v(t) (Flury, ~p.20):
      Ft = 1.0 - a*(1.0 - cos(E-Epi))/rpil;
      Gt = (M-E+Epi+sin(E-Epi))/n; // t-tpi-(E-Epi-sin(E-Epi))/n;

      r.x = Ft*rpi.x + Gt*vpi.x;
      r.y = Ft*rpi.y + Gt*vpi.y;
      r.z = Ft*rpi.z + Gt*vpi.z;

      rl = sqrt(scalar_product(r,r));
      Fpt = -sqrt(mue*a)/(rl*rpil)*sin(E-Epi);
      Gpt = 1.0 - a*(1.0 - cos(E-Epi))/rl;

      v.x = Fpt*rpi.x + Gpt*vpi.x;
      v.y = Fpt*rpi.y + Gpt*vpi.y;
      v.z = Fpt*rpi.z + Gpt*vpi.z;
      // end of unperturbed calculation

      // perturbation effects (according to Nitschke p. 25)
      cosi = cos(i);
      cosi2 = pow(cosi,2.0);
      // change of the ascending node:
      Omega += ds * cos(i);
      // change of the argument of perigee:
      omega += (ds/2.0) * (1.0 - 5.0*cosi2);
      // change of the mean eccentric anomaly:
      M -= (ds/2.0) * sqrt(1.0-pow(e,2.0)) * (3.0*cosi2 - 1.0);

      if (++plotstep==Nplotintervals) {
        fprintf(outfptr, "%lf\t%lf\t%lf\t%lf\n", t, Omega,omega,M);
//        fprintf(outfptr, "%lf\t%lf\t%lf\n", r.x/R_e, r.y/R_e, r.z/R_e);   // \vec{r}(t)
//        fprintf(outfptr, "%lf\t%lf\t%lf\n", v.x, v.y, v.z);             // \vec{v}(t)
        plotstep = 0;
      }
    }
    // close the output-file
    fclose(outfptr);
  } else {
    // output-file was not opened => print error message
    printf("Error: output file '%s' could not be opened!\n\n", outputfile);
  }

  return(EXIT_SUCCESS);
}




double kepler_equation(double M, double e) {
  double E = M;
  double dE;

  // Newton algorithm
  do {
    dE = (M - E + e*sin(E))/(1 - e*cos(E));
    E += dE;
  } while (fabs(dE) > 0.000001);

  return(E);
}
