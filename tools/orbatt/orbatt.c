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

// Reads the program parameters using PIL.
static int orbatt_getpar(struct Parameters* const par);


//////////////////////////////////////////////////////////////////////
// Main procedure.
int orbatt_main()
{
  struct Parameters par;

  int status=EXIT_SUCCESS;


  // Register HEATOOL.
  set_toolname("orbatt");
  set_toolversion("0.01");


  do { // Beginning of error handling program.

    // Read in the parameters using PIL.
    status = orbatt_getpar(&par);
    if (EXIT_SUCCESS!=status) break;

    // TODO

  } while(0); // END of error handling routine.

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

    // first delete the old output files (ASCII and FITS)
    remove(orbit_asciifile);
    remove(orbit_fitsfile);
    remove("test.fits");

    // open ASCII file for the output of the position and velocity of the satellite
    if(!(ascii_outfptr=fopen(orbit_asciifile, "w"))) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error opening the ASCII file '%s' for output", orbit_asciifile);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // create new FITS output file (orbit file)
    if (fits_create_file(&fits_outfptr, orbit_fitsfile, &status)) break;

    // variables for the creation of the fits-table (field-type, -format and -units of individual columns)
    char *ftype[N_ORBIT_FIELDS];
    char *fform[N_ORBIT_FIELDS];
    char *funit[N_ORBIT_FIELDS];

    // create a binary table in the FITS file
    create_orbtbl_parameter(ftype, fform, funit);
    if (fits_create_tbl(fits_outfptr, BINARY_TBL, 0, N_ORBIT_FIELDS, ftype, fform, funit, "eROSITA orbit", &status)) break;
    


    // write header into FITS file
    // descriptory comments and mission headers:
    fits_write_key (fits_outfptr, TSTRING, "MISSION", "SpectrumXGamma", "name of the mission", &status);
    fits_write_key(fits_outfptr, TSTRING, "COMMENT", "DESCRIPT","orbit file for eROSITA simulation", &status);
    fits_write_key(fits_outfptr, TSTRING, "COMMENT", "CONTENT","describes the orbit of the eROSITA satellite", &status);
    fits_write_key(fits_outfptr, TSTRING, "COMMENT", "FORMAT","(t,x,y,z,vx,vy,vz)", &status);
    fits_write_key(fits_outfptr, TSTRING, "COMMENT", "FORMAT","(time,position,velocity)", &status);
    // about this orbit-calculation program:
    fits_write_key(fits_outfptr, TSTRING, "PROGRAM", "calc_orbit", "orbit-calculation program", &status);
    // history:
    fits_write_key(fits_outfptr, TSTRING, "HISTORY", "program calc_orbit", "program name", &status);
    fits_write_key(fits_outfptr, TSTRING, "HISTORY", "version 1.0", "version number", &status);

    // date & time headers:
    double dummy = 0.0;
    long lbuffer = 0;

    char creation_date[30];
    int timeref;            // is 0, if returned time is in UTC
    fits_get_system_time(creation_date, &timeref, &status);
    fits_write_key(fits_outfptr, TSTRING, "DATE", "2008-03-10","FITS file creation date (yyyy-mm-dd)", &status);

    fits_write_key(fits_outfptr, TSTRING, "DATE-OBS", "2008-03-10T10:00:00","start time for the orbit", &status);
    fits_write_key(fits_outfptr, TSTRING, "DATE-END", "","end time of the orbit", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "MJDSTART", &dummy,"start time of the orbit in Julian date format", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "MJDEND", &dummy,"end time of the orbit in Julian date format", &status);
    dummy = 0.;
    fits_write_key(fits_outfptr, TDOUBLE, "TIMEZERO", &dummy,"Clock correction", &status);
    fits_write_key(fits_outfptr, TLONG, "MJDREFI", &lbuffer,"integer part of reference time", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "MJDREFF", &dummy,"fractional part of reference time", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "TSTART", &t0,"start time of the orbit", &status);
    dummy = t0 + timespan;
    fits_write_key(fits_outfptr, TDOUBLE, "TEND", &dummy,"start time of the orbit", &status);

    // initial orbit parameters:
    fits_write_key(fits_outfptr, TDOUBLE, "t0", &t0, "t0", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "a0", &a0, "inital semimajor axis at t0", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "e0", &e0, "inital eccentricity at t0", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "i0", &i0, "inital inclination at t0", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "ra_an0", &Omega0, "inital right ascension of ascending node at t0", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "arg_per0", &omega0, "inital argument of perigee at t0", &status);
    fits_write_key(fits_outfptr, TDOUBLE, "M0", &M0, "inital mean anomaly at t0", &status);
    
    // HEADAS parstamp: writes all parameters of this program call into the FITS file header
    HDpar_stamp(fits_outfptr, 2, &status);

    // if an error has occured during writing the headers, break and continue with cleaning up
    if (status) break;



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
