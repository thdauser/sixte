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
// The orbit parameters can either be given as program arguments on
// the command line or be specified in a TLE file. In the latter case,
// the program call looks like the following:
// "./calc_orbit -tle [name of TLE file]" ( !! NOT IMPLEMENTED YET !! )
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
  /** Timespan for the calculation. */
  double timespan;
  /** Time difference between two calculation steps. */
  double dt;
  /** Filename of the orbit output file. */
  char orbit_filename[MAXFILENAME];
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
  set_toolversion("0.01");


  do { // Beginning of error handling program.

    // Read in the parameters using PIL.
    status = orbatt_getpar(&par);
    if (EXIT_SUCCESS!=status) break;

    // Create the output FITS file.
    remove(par.orbit_filename);
    if (fits_create_file(&fptr, par.orbit_filename, &status)) break;
    char* ttype[] = { "TIME", "X", "Y", "Z", "VX", "VY", "VZ" };
    char* tform[] = { "D"   , "D", "D", "D", "D" , "D" , "D"  };
    char* tunit[] = { "s"   , "m", "m", "m", "m/s","m/s","m/s"};
    if (fits_create_tbl(fptr, BINARY_TBL, 0, 7, ttype, tform, tunit, 
			"ORBIT", &status)) break;
    HDpar_stamp(fptr, 2, &status);
    long nrows=0;
    

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
    headas_printf("orbit period: %lf min\n", T/60.);
    headas_printf("orbit precession period: %lf d\n", 
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
      // satellite (Flury p. 38)
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
    }
    // End of loop over time interval.

  } while(0); // END of error handling routine.

  if (NULL!=fptr) fits_close_file(fptr, &status);

  if (EXIT_SUCCESS==status) {
    headas_chat(0, "finished successfully!\n");
  }

  return(status);
}



static int orbatt_getpar(struct Parameters* const par)
{
  int status=EXIT_SUCCESS; // Error status

  // Read the length of the semimajor axis a.
  if ((status = PILGetReal("a0", &par->a0))) {
    HD_ERROR_THROW("Error reading the 'a0' parameter", status);
  }
    
  // Read the eccentricity of the Kepler orbit.
  else if ((status = PILGetReal("e0", &par->e0))) {
    HD_ERROR_THROW("Error reading the 'e0' parameter", status);
  }
   
  // Read the inclination.
  else if ((status = PILGetReal("i0", &par->i0))) {
    HD_ERROR_THROW("Error reading the 'i0' parameter", status);
  }

  else if ((status = PILGetReal("gOmega0", &par->Omega0))) {
    HD_ERROR_THROW("Error reading the 'Omega0' parameter", status);
  }

  else if ((status = PILGetReal("komega0", &par->omega0))) {
    HD_ERROR_THROW("Error reading the 'omega0' parameter", status);
  }

  else if ((status = PILGetReal("M0", &par->M0))) {
    HD_ERROR_THROW("Error reading the 'M0' parameter", status);
  }
   
  else if ((status = PILGetReal("t0", &par->t0))) {
    HD_ERROR_THROW("Error reading the 't0' parameter", status);
  }

  else if ((status = PILGetReal("timespan", &par->timespan))) {
    HD_ERROR_THROW("Error reading the 'timespan' parameter", status);
  }

  else if ((status = PILGetReal("dt", &par->dt))) {
    HD_ERROR_THROW("Error reading the 'dt' parameter", status);
  }

  else if ((status = PILGetFname("orbit_filename", par->orbit_filename))) {
    HD_ERROR_THROW("Error reading the name of the orbit file", status);
  }
  if (EXIT_SUCCESS!=status) return(status);

  // Convert angles [degrees] -> [radians]:
  par->i0     = par->i0     *M_PI/180.0;
  par->Omega0 = par->Omega0 *M_PI/180.0;
  par->omega0 = par->omega0 *M_PI/180.0;
  par->M0     = par->M0     *M_PI/180.0;

  return(status);
}



static double kepler_equation(const double M, const double e) 
{
  double E = M;
  double dE;

  // Newton algorithm to solve the Kepler equation
  do {
    dE = (M - E + e*sin(E))/(1.-e*cos(E));
    E += dE;
  } while (fabs(dE) > 0.00000001);

  return(E);
}



/*
//////////////////////////////////////////////////////////////////////////
// Performs the actual work: calculate the orbit and save the output
// positions into files (FITS file and optional debug ASCII file).
int calc_orbit_work(
		    // initial Keplerian orbital elements (@t_0):
		    double a0,            // semimajor axis
		    double e0,            // eccentricity
		    double i0,            // inclination
		    double Omega0,        // initial right ascension of ascending node
		    double omega0,        // initial argument of perigee
		    double M0,            // mean anomaly at t0=0 ( M0=0.0 => satellite is in perigee)
		    double t0,            // time-offset
		    // calculation parameters:
		    double timespan,      // timespan in which the satellite data is calculated
		    double dt,            // width of a timestep between two orbit calculations
		    double plotinterval,  // width of a timestep for the outputs of orbit positions
		    // file names
		    const char orbit_fitsfile[],    // filename of the FITS output file
		    const char orbit_asciifile[]    // filename of the ASCII output file
		    )
{
  //  const long double mu = 398601.3; // G*Msolar; units: (km^3/s^2)


  // variables:
  struct orbit_data odata;         // structure containing the orbit data (mainly Keplerian orbital elements)
  struct orbit_position oposition; // structure containing time,  position and velocity of the satellite
  fitsfile *fits_outfptr=NULL;     // filepointer to FITS output-file (orbit file)
  FILE *ascii_outfptr=NULL;        // filepointer to ASCII output-file
  long row;                        // row counter for FITS file
  double t_last_plot;              // point of time of the last data output (to orbit file)


  // error handling variables:
  char msg[MAXMSG];         // buffer for error output messages
  int status=0;             // error status


  do {    // beginning of the error handling loop (only run once)

    // initialize orbit-data structure:
    odata.dt = dt;

    // set the last_plot value in such a way, that there is a data output at the first calculation step:
    t_last_plot = t0 - plotinterval;

//
// loop over all inclination from 0 to 180°C (needed for program verification)
//int inc;
//for (inc = 1; inc < 20; inc++) {
//  i0 = inc*M_PI/180.0;
//

    // calculation of basic orbit data:
    orbit_init(a0, e0, i0, Omega0, omega0, M0, &odata);


    // Now we have the necessary data to calculate the position and velocity of
    // the satellite on its orbit at any given point of time.
    // So we can determine the position of the satellite at any point of time 't':


    ///////////////////////////////////////////////////////
    // orbit calculation
    // loop over the entire required timeinterval from t0 to t0+timespan
    for (oposition.t=t0, row=0; oposition.t<t0+timespan; oposition.t+=odata.dt) {
      // perform an iteration step on the orbital elements
      orbit_step_J234(&odata, &oposition);
	
      // output of orbit data into files, if the required time interval per output is over
      if (oposition.t - t_last_plot >= plotinterval) {
	// output into FITS file
	add_orbtbl_row(fits_outfptr, row++, oposition.t, oposition.r, oposition.v, &status);

	// output into ASCII file
	//        fprintf(ascii_outfptr, "%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n", oposition.t, odata.a, odata.e, odata.i*180./M_PI, odata.Omega*180./M_PI, odata.omega*180./M_PI, odata.M*180./M_PI);
	//        fprintf(ascii_outfptr, "%lf\t%lf\t%lf\t", oposition.r.x, oposition.r.y, oposition.r.z);
	//        fprintf(ascii_outfptr, "%lf\t%lf\t%lf\n", oposition.v.x, oposition.v.y, oposition.v.z);

	// output of orbital elements to STDOUT
	//fprintf(stdout, "%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n", oposition.t, odata.a, odata.e, odata.i*180./M_PI, odata.Omega*180./M_PI, odata.omega*180./M_PI, odata.M*180./M_PI);
	  
	// update the time of the last data output:
	t_last_plot = oposition.t;
      }
    }  // end of the loop over the entire timespan
      

    //    printf("%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t%Lf\n", oposition.t, odata.a, odata.e, odata.i*180./M_PI, odata.Omega*180./M_PI, odata.omega*180./M_PI, odata.M*180./M_PI, (odata.omega+odata.M)*180./M_PI);
    
    // }  // end of loop over different inclinations

  } while (0);     // end of error handling loop


  // clean up (close files):
  status = calc_orbit_cleanup(fits_outfptr, ascii_outfptr, status);

  return(status);
}





//////////////////////////////////////////////////////////////////////////
// performs the clean-up tasks: closing files, releasing memory, ...
int calc_orbit_cleanup(
		       fitsfile *fits_outfptr,   // FITS file pointer to orbit output file
		       FILE *ascii_outfptr,      // file pointer to ASCII orbit output file
		       int estatus               // error status
		       )
{
  char msg[MAXMSG];         // buffer for error output messages
  int status=estatus;       // error status

  // close the FITS orbit file
  if (fits_outfptr) fits_close_file(fits_outfptr, &status);

  // close the ASCII output-file
  if (ascii_outfptr) {
    if (fclose(ascii_outfptr)) {
      status += EXIT_FAILURE;
      sprintf(msg, "Error reading the 'orbit_asciifile' parameter");
      HD_ERROR_THROW(msg,status);
    }
  }
  
  return(status);
}

*/
