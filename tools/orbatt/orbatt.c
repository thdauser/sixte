//////////////////////////////////////////////////////////////////////
//
// This program is part of the SIXT software package and calculates a
// satellite's orbit considering perturbation effects.
//
//////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @data     2008/04
// @param    see below
//
//////////////////////////////////////////////////////////////////////
//
// The program calculates the orbit of a satellite according to given
// Keplerian orbit parameters:
// -> a - semimajor axis
// -> e - eccentricity
// -> i - inclination
// -> Omega - right ascension of ascending node
// -> omega - argument of perigee
// -> M0    - mean anomaly at t0 (usually M0=0, i.e. the satellite is in perigee)
// -> t0
//
// output: time t, position \vec{r}(t) and velocity \vec{v}(t)
//
//////////////////////////////////////////////////////////////////////

#include "sixt.h"
#include "vector.h"

#define TOOLSUB orbatt_main
#include "headas_main.c"


struct Parameters {
  /** Initial semimajor axis. */
  double a0;
  /** Initial eccentricity. */
  double e0;
  /** Initial inclination. */
  double i0;
  /** Initial right ascension of ascending node. */
  double Omega0;
  /** Initial argument of perigee. */
  double omega0;
  /** Initial mean anomaly. */
  double M0;
  /** Initial time. */
  double t0;
  /** Timespan for the calculation [s]. */
  double timespan;
  double MJDREF;
  /** Time difference between two calculation steps [s]. */
  double dt;
  /** Filename of the orbit output file. */
  char orbit_filename[MAXFILENAME];

  char clobber;
};


/** Reads the program parameters using PIL. */
static int orbatt_getpar(struct Parameters* const par);

/** Solves the Kepler equation. Calculates the eccentric anomaly from
    the mean anomaly. */
static double kepler_equation(const double M, const double e);


//////////////////////////////////////////////////////////////////////
// Main procedure.
int orbatt_main()
{
  struct Parameters par;
  fitsfile* fptr=NULL;

  int status=EXIT_SUCCESS;


  // Register HEATOOL.
  set_toolname("orbatt");
  set_toolversion("0.04");


  do { // Beginning of error handling program.

    // Read in the parameters using PIL.
    status=orbatt_getpar(&par);
    CHECK_STATUS_BREAK(status);;

    // Check if the output file already exists.
    int exists;
    fits_file_exists(par.orbit_filename, &exists, &status);
    CHECK_STATUS_BREAK(status);
    if (0!=exists) {
      if (0!=par.clobber) {
	// Delete the file.
	remove(par.orbit_filename);
      } else {
	// Throw an error.
	char msg[MAXMSG];
	sprintf(msg, "file '%s' already exists", par.orbit_filename);
	SIXT_ERROR(msg);
	status=EXIT_FAILURE;
	break;
      }
    }

    // Create the output FITS file.
    fits_create_file(&fptr, par.orbit_filename, &status);
    CHECK_STATUS_BREAK(status);;
    char* ttype[] = { "TIME", "X", "Y", "Z", "VX", "VY", "VZ" };
    char* tform[] = { "D"   , "D", "D", "D", "D" , "D" , "D"  };
    char* tunit[] = { "s"   , "m", "m", "m", "m/s","m/s","m/s"};
    fits_create_tbl(fptr, BINARY_TBL, 0, 7, ttype, tform, tunit, 
		    "ORBIT", &status);
    CHECK_STATUS_BREAK(status);;
    HDpar_stamp(fptr, 2, &status);
    CHECK_STATUS_BREAK(status);;
    long nrows=0;
    
    // Set the timing keywords.
    fits_update_key(fptr, TDOUBLE, "TSTART", &par.t0, 
		    "start time", &status);
    double dbuffer=0.0;
    fits_update_key(fptr, TDOUBLE, "TIMEZERO", &dbuffer, 
		    "zero time", &status);
    fits_update_key(fptr, TDOUBLE, "MJDREF", &par.MJDREF, 
		    "MJD for reference time", &status);
    fits_update_key(fptr, TSTRING, "TIMEUNIT", "s", 
		    "Unit for TSTART, TSTOP, TIMEZERO", &status);
    CHECK_STATUS_BREAK(status);;

    // Setup initial orbit parameters.
    const double mu  = 3.987e14;  // G*M_earth ([m^3/s^2])
    const double R_e = 6378.14e3; // [m]
    double a = par.a0;
    double e = par.e0;
    double p = a*(1-pow(e,2.0));
    double i = par.i0;
    double n = sqrt(mu/pow(a,3.0)); // mean motion.
    const double J2  = 1082.6268e-6;
    const double prefactor = -1.5*J2*n*powl(R_e/p,2.);

    // Determine the orbit period (for a Kepler orbit) an print 
    // it to the screen.
    double T = 2*M_PI/n;
    headas_printf("orbit period: %.1lf min\n", T/60.);
    headas_printf("orbit precession period: %.1lf d\n", 
		  2*M_PI/fabs(prefactor*cos(i))/(3600.*24.));

    // Loop over the requested time interval.
    double time; // The absolute time is: 'par.t0 + time'
    for (time=0.; time<=par.timespan; time+=par.dt) {
      // Iterate the Keplerian orbital elements.

      // Secular perturbation effects (J_2 contributions only).
      // Right ascension of ascending node.
      double Omega = par.Omega0 + prefactor * cos(i) * time;
      while (Omega < 0.) {
	Omega += 2.*M_PI;
      }
      while (Omega > 2.*M_PI) {
	Omega -= 2.*M_PI;
      }
      // Argument of perigee.
      double omega = par.omega0 + prefactor*0.5*(1.-5.*pow(cos(i),2.)) * time;
      while (omega < 0.) {
	omega += 2.*M_PI;
      }
      while (omega > 2.*M_PI) {
	omega -= 2.*M_PI;
      }
      // Mean eccentric anomaly.
      double M = par.M0 + 
	(n - prefactor*0.5*sqrt(1.-pow(e,2.))*(3.*pow(cos(i),2.)-1.)) * time;
      while (M < 0.) {
	M += 2.*M_PI;
      }
      while (M > 2.*M_PI) {
	M -= 2.*M_PI;
      }
      // End of secular perturbation terms.

      // Determine the eccentric anomaly E(t) by solving the Kepler 
      // equation (Flury p. 18).
      double E = kepler_equation(M, e);

      // Calculate the true anomaly from the eccentric anomaly (Flury p. 16).
      double f = 2.*atan(sqrt((1.+e)/(1.-e))*tan(E/2.));

      // Now we know the argument of latitude and can calculate 
      // its cosine and sine.
      double cosu = cos(omega+f);
      double sinu = sin(omega+f);
  
      // Now we can determine the position and velocity of the 
      // satellite (Flury p. 38).
      double rl = p/(1.+e*cos(f));
      double vr = sqrt(mu/p)*e*sin(f);
      double vf = sqrt(mu*p)/rl;
  
      // Position and velocity in Carteesian coordinates.
      Vector position = {
	.x = rl*(cos(Omega)*cosu-sin(Omega)*cos(i)*sinu),
	.y = rl*(sin(Omega)*cosu+cos(Omega)*cos(i)*sinu),
	.z = rl*sin(i)*sinu
      };
      Vector velocity = {
	.x = vr*(cos(Omega)*cosu-sin(Omega)*cos(i)*sinu) -
	vf*(cos(Omega)*sinu+sin(Omega)*cosu*cos(i)),
	.y = vr*(sin(Omega)*cosu+cos(Omega)*cos(i)*sinu) - 
	vf*(sin(Omega)*sinu-cos(Omega)*cosu*cos(i)),
	.z = vr*sin(i)*sinu + vf*cosu*sin(i)
      };
      // END of orbit iteration

      // Store the new parameters in the output file.
      double t_abs = par.t0 + time; // Absolute time.
      fits_insert_rows(fptr, nrows++, 1, &status);
      fits_write_col(fptr, TDOUBLE, 1, nrows, 1, 1, &t_abs, &status);
      fits_write_col(fptr, TDOUBLE, 2, nrows, 1, 1, &position.x, &status);
      fits_write_col(fptr, TDOUBLE, 3, nrows, 1, 1, &position.y, &status);
      fits_write_col(fptr, TDOUBLE, 4, nrows, 1, 1, &position.z, &status);
      fits_write_col(fptr, TDOUBLE, 5, nrows, 1, 1, &velocity.x, &status);
      fits_write_col(fptr, TDOUBLE, 6, nrows, 1, 1, &velocity.y, &status);
      fits_write_col(fptr, TDOUBLE, 7, nrows, 1, 1, &velocity.z, &status);
      CHECK_STATUS_BREAK(status);;
    }
    CHECK_STATUS_BREAK(status);;
    // End of loop over time interval.

    // Update the TSTOP header keyword.
    dbuffer=par.t0+time;
    fits_update_key(fptr, TDOUBLE, "TSTOP", &dbuffer, 
		    "stop time", &status);
    CHECK_STATUS_BREAK(status);;

  } while(0); // END of error handling routine.

  if (NULL!=fptr) fits_close_file(fptr, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(0, "finished successfully!\n");
  }

  return(status);
}



static int orbatt_getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  // Read the length of the semimajor axis a.
  status=ape_trad_query_double("a0", &par->a0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'a0' parameter");
    return(status);
  } 
    
  // Read the eccentricity of the Kepler orbit.
  status=ape_trad_query_double("e0", &par->e0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'e0' parameter");
    return(status);
  } 
   
  // Read the inclination.
  status=ape_trad_query_double("i0", &par->i0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'i0' parameter");
    return(status);
  } 

  status=ape_trad_query_double("gOmega0", &par->Omega0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'Omega0' parameter");
    return(status);
  } 

  status=ape_trad_query_double("komega0", &par->omega0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'omega0' parameter");
    return(status);
  } 

  status=ape_trad_query_double("M0", &par->M0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'M0' parameter");
    return(status);
  } 
   
  status=ape_trad_query_double("t0", &par->t0);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 't0' parameter");
    return(status);
  } 

  status=ape_trad_query_double("MJDREF", &par->MJDREF);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading MJDREF");
    return(status);
  } 

  status=ape_trad_query_double("timespan", &par->timespan);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'timespan' parameter");
    return(status);
  } 

  status=ape_trad_query_double("dt", &par->dt);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the 'dt' parameter");
    return(status);
  } 

  status=ape_trad_query_string("orbit_filename", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the orbit file");
    return(status);
  } 
  strcpy(par->orbit_filename, sbuffer);
  free(sbuffer);

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the clobber parameter");
    return(status);
  }

  // Convert angles [degrees] -> [radians]:
  par->i0     = par->i0     *M_PI/180.0;
  par->Omega0 = par->Omega0 *M_PI/180.0;
  par->omega0 = par->omega0 *M_PI/180.0;
  par->M0     = par->M0     *M_PI/180.0;

  return(status);
}


static double kepler_equation(const double M, const double e) 
{
  double E=M;
  double dE;

  // Newton algorithm to solve the Kepler equation
  do {
    dE =(M - E + e*sin(E))/(1.-e*cos(E));
    E +=dE;
  } while (fabs(dE) > 0.00000001);

  return(E);
}

