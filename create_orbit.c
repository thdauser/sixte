//////////////////////////////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and calculates a satellite's
// orbit considering perturbation effects.
//
/////////////////////////////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @data     2008/04
// @param    see below
//
/////////////////////////////////////////////////////////////////////////////////////////////////
//
// Program calculates the orbit of a satellite according to given Keplerian orbit parameters:
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
// The orbit parameters can either be given as program arguments on the command line or 
// be specified in a TLE file. In the latter case, the program call looks like the following:
// "./calc_orbit -tle [name of TLE file]" ( !! NOT IMPLEMENTED YET !! )
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pil.h"
#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "vector.h"
#include "fits_ctlg.h"
#include "orbit_calc.h"
#include "orbit_calc.c"
#include "earth_constants.h"

#define TOOLSUB calc_orbit_main
#include "headas_main.c"

#define MAXMSG 256
#define FILENAME_LENGTH 128



// reads the program parameters using PIL:
int calc_orbit_getpar(double *a0, double *e0, double *i0, double *Omega0, double *omega0, double *M0, double *t0, double *timespan, double *dt, double *plotinterval, char orbit_fitsfile[], char orbit_asciifile[]);

// performs the actual work: calculate the orbit and save the output
// positions into files (FITS file and optional debug ASCII file)
int calc_orbit_work(double a0, double e0, double i0, double Omega0, double omega0, double M0, double t0, double timespan, double dt, double plotinterval, const char orbit_fitsfile[], const char orbit_asciifile[]);

// performs the clean-up tasks: closing files, releasing memory, ...
int calc_orbit_cleanup(fitsfile *fits_outfptr, FILE *ASCII_outfptr, int estatus);


////////////////////////////////////////////////////////////////////////////
// main procedure
int calc_orbit_main()
{
  // input parameters:
  double a0;            // semimajor axis
  double e0;            // eccentricity
  double i0;            // inclination
  double Omega0;        // initial right ascension of ascending node
  double omega0;        // initial argument of perigee
  double M0;            // mean anomaly at t0=0 ( M0=0.0 => satellite is in perigee)
  double t0;            // time-offset
  double timespan;      // timespan in which the satellite data is calculated
  double dt;            // width of a timestep for an orbit calculation
  double plotinterval;  // width of a timeinterval for the output of an orbit position
  char orbit_asciifile[FILENAME_LENGTH];  // filename of the output-file
  char orbit_fitsfile[FILENAME_LENGTH];   // filename of the FITS orbit-file (outputfile)

  int status=0;         // error status


  // register HEATOOL
  set_toolname("create_orbit");
  set_toolversion("0.01");


  // read in parameters using PIL
  status = calc_orbit_getpar(&a0, &e0, &i0, &Omega0, &omega0, &M0, &t0, &timespan, &dt, &plotinterval, orbit_fitsfile, orbit_asciifile);


  // call the routine, which does the actual work (initialize orbit, calculate orbit positions
  // and store them to the orbit file)
  if (!status) {
    status = calc_orbit_work(a0, e0, i0, Omega0, omega0, M0, t0, timespan, dt, plotinterval, orbit_fitsfile, orbit_asciifile);

    headas_chat(5, "finished\n");
  }


  return(status);
}



////////////////////////////////////////////////////////////////////////////////
// reads the program parameters using PIL:
int calc_orbit_getpar(
		      double *a0,               // initial semimajor axis
		      double *e0,               // initial eccentricity
		      double *i0,               // initial inclination
		      double *Omega0,           // initial right ascension of ascending node
		      double *omega0,           // initial argument of perigee
		      double *M0,               // initial mean anomaly
		      double *t0,               // initial time
		      double *timespan,         // timespan for the calculation
		      double *dt,               // time difference between two calculation steps
		      double *plotinterval,     // time interval between two position outputs
		      char orbit_fitsfile[],    // filename of the FITS output file
		      char orbit_asciifile[]    // filename of the ASCII output file
		      ) 
{
  char msg[MAXMSG];             // error output buffer
  int status=0;                 // error status


/*
  // check if input from a TLE file is requested
  if ((argc>1)&&(strcmp(argv[1],"-tle")==0)) {
    // read the filename of the TLE file (2nd parameter) using PIL library
    pilread = PILGetString("tleflag", tle_flag);
    pilread = PILGetFname("tlefile", tle_filename);

    // read the TLE element from the file into a string
    // open ASCII file for reading
    tle_inpfile=fopen(tle_filename, "r");
    if (tle_inpfile != NULL) {
      // read TLE element into string
      int row = 0;
      while(fgets(inputbuffer, 100, tle_inpfile)&&(row<3)) {
	// append line to tle_element (consists of 3 lines)
	strcat(tle_element, inputbuffer);
      }
      // close TLE ASCII file
      fclose(tle_inpfile);
    }

    // parse the string, using parsing routine from the orbit package
    int epoch_year;
    double epoch_day;
    double n0;
    if(!parse_tle_element(tle_element, &epoch_year, &epoch_day, &e0, &i0, &Omega0, &omega0, &M0, &n0)) {
      a0 = pow(mu/pow(n0*2.*M_PI/(3600.*24.),2.),1./3.);
      t0 = ((double)(epoch_year-2011)*365.25 + epoch_day)*24.*3600;
    } else {
      printf("Error: reading TLE file was not successful!\n");
      return(-1);
    }

  } else {
*/

  // read the length of the semimajor axis a
  if ((status = PILGetReal("a0", a0))) {
    sprintf(msg, "Error reading the 'a0' parameter");
    HD_ERROR_THROW(msg,status);
  }
    
  // read the eccentricity of the Kepler orbit
  else if ((status = PILGetReal("e0", e0))) {
    sprintf(msg, "Error reading the 'e0' parameter");
    HD_ERROR_THROW(msg,status);
  }
   
  // read the length of the semimajor axis a
  else if ((status = PILGetReal("i0", i0))) {
    sprintf(msg, "Error reading the 'i0' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // read the right ascension of the ascending node
  else if ((status = PILGetReal("gOmega0", Omega0))) {
    sprintf(msg, "Error reading the 'Omega0' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // read the argument of perigee
  else if ((status = PILGetReal("komega0", omega0))) {
    sprintf(msg, "Error reading the 'omega0' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // read the mean anomaly @ t0
  else if ((status = PILGetReal("M0", M0))) {
    sprintf(msg, "Error reading the 'M0' parameter");
    HD_ERROR_THROW(msg,status);
  }
   
  // read the time offset t0
  else if ((status = PILGetReal("t0", t0))) {
    sprintf(msg, "Error reading the 't0' parameter");
    HD_ERROR_THROW(msg,status);
  }

/*  } */

  // read the entire timespan of the calculation
  else if ((status = PILGetReal("timespan", timespan))) {
    sprintf(msg, "Error reading the 'timespan' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // read the width of a timestep
  else if ((status = PILGetReal("dt", dt))) {
    sprintf(msg, "Error reading the 'dt' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // read the number of timesteps performed per data-output
  else if ((status = PILGetReal("plotinterval", plotinterval))) {
    sprintf(msg, "Error reading the 'plotinterval' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // get the filename of the FITS output-file
  else if ((status = PILGetFname("orbit_fitsfile", orbit_fitsfile))) {
    sprintf(msg, "Error reading the 'orbit_fitsfile' parameter");
    HD_ERROR_THROW(msg,status);
  }

  // get the filename of the ASCII output-file
  else if ((status = PILGetFname("orbit_asciifile", orbit_asciifile))) {
    sprintf(msg, "Error reading the 'orbit_asciifile' parameter");
    HD_ERROR_THROW(msg,status);
  }


  // convert angles [degrees] -> [radians]:
  *i0 = *i0*M_PI/180.0;
  *Omega0 = *Omega0*M_PI/180.0;
  *omega0 = *omega0*M_PI/180.0;
  *M0 = *M0*M_PI/180.0;


  return(status);
}






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

  /*  
  char tle_flag[100];       // flag for using TLE files
  char tle_filename[100];   // filename of the TLE input file (only used if tleflag=="-tle")
  FILE *tle_inpfile;        // filepointer to the TLE input file (ASCII file)
  char inputbuffer[100];    // input buffer for the TLE file access
  char tle_element[200];    // buffer for the TLE element (before parsing)
  */

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
    
    /*
    int count;
    for (count=0; count<N_ORBIT_FIELDS; count++) {
      if (ftype[count]) free(ftype[count]);
      if (fform[count]) free(fform[count]);
      if (funit[count]) free(funit[count]);
    }
    */



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
	// To obtain a fixed satellite position, use the following lines:
	/*
	oposition.r.x = 1.;
	oposition.r.y = 0.;
	oposition.r.z = 0.;
	oposition.v.x = 0.;
	oposition.v.y = 0.;
	oposition.v.z = 1.;
	*/
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

