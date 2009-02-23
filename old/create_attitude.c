/////////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and calculates 
// an attitude file to a corresponding orbit file.
//
/////////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @data     2008/04
// @param    orbitfile - filename of the orbit file (FITS input file)
// @param    attitudefile - filename of the attitude file (FITS output file)
//
/////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "fits_ctlg.h"
#include "fits_attitude.h"

#define TOOLSUB create_attitude_main
#include "headas_main.c"


#define FILENAME_LENGTH 128   // maximum length of filenames
#define MAXMSG 256            // maximum length of an output/error message




// reads all parameters of 'conv_rosat2fits' using PIL
int create_attitude_getpar(char orbit_filename[], char attitude_filename[]);

// does the acutal work: create FITS file with table, open ASCII file, 
// transfer sources, close files, clean up
int create_attitude_work(const char orbit_filename[], const char attitude_filename[]);

// calculates right ascension and declination of the sun to the corresponding
// Julian Day
void calc_sun_pos(double JD, double *sun_ra, double *sun_dec);


////////////////////////////////////////////////////////////////////////////////////////
// main procedure
int create_attitude_main()
{
  char orbit_filename[FILENAME_LENGTH];      // name of input file (ASCII file)
  char attitude_filename[FILENAME_LENGTH];   // name of output file (FITS file)

  int status=0;             // error report status


  // HEATOOLs: register program
  set_toolname("create_attitude");
  set_toolversion("0.01");


  // read parameters using PIL library
  status = create_attitude_getpar(orbit_filename, attitude_filename);


  if(!status) {
    // perform the actual work: open orbit file, create attitude file, calculate and write attitude data
    status = create_attitude_work(orbit_filename, attitude_filename);

    headas_chat(5, "finished\n");
  }

  return(status);
}



///////////////////////////////////////////////////////////////////////////////////////////////
// reads all parameters of 'create_attitude' using PIL
int create_attitude_getpar(
			   char orbit_filename[],    // filename of the orbit file (FITS input file)
			   char attitude_filename[]  // filename of the attitude file (FITS output file)
			   ) 
{
  int status=0;         // error status
  char msg[MAXMSG];     // buffer for error output messages

  if ((status = PILGetFname("orbitfile", orbit_filename))) {
    sprintf(msg, "Error reading the 'orbitfile' parameter");
    HD_ERROR_THROW(msg,status);
  }

    else if ((status = PILGetFname("attitudefile", attitude_filename))) {
    sprintf(msg, "Error reading the 'attitudefile' parameter");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}



//////////////////////////////////////////////////////////////////////////////////////////////////
// does the acutal work: create FITS file with table, open ASCII file, transfer sources, close files, clean up
int create_attitude_work(
			 const char orbit_filename[],    // filename of the orbit file (FITS input file)
			 const char attitude_filename[]  // filename of the attitude file (FITS output file)
			 )
{
  fitsfile *orbit_fptr=NULL;     // fitsfile pointer to orbit input file
  fitsfile *attitude_fptr=NULL;  // fitsfile pointer to attitude output file
  int orbit_hdunum;              // HDU number in orbit file
  int orbit_hdutype;             // HDU type in orbit file
  long orbit_nrows;              // number of rows in orbit file

  double time;                   // buffer for the orbit data: time
  struct vector r, v;            // -"-                      : position and velocity of the satellite
  struct vector sz0;             // normalized position vector of the satellite (= z-direction in satellite coordinate system)
  struct vector sx0, sy0;        // normalized vectors in satellite coordinate system
  char valtime[19];              // time in human readable format (yyyy-mm-ddThh:mm:ss)
  double JD;                     // Julian day
  double sun_ra, sun_dec;        // right ascension and declination of the sun (in geocentric coordinates)
  struct vector sun_pos;         // unit vector pointing in the direction of the sun
  double view_ra, view_dec;      // right ascension and declination of the telescope viewing direction
  double rollangle;              // roll angle of the satellite (here: angle around the z-axis)
  double aspangle;               // solar aspect angle

  int status=0;                  // error status
  char msg[MAXMSG];              // buffer for error output messages

  do {  // error handling loop (is only run once)

    // open orbit file for data input
    if (fits_open_file(&orbit_fptr, orbit_filename, READONLY, &status)) break;
    // get the number of the current HDU
    if (fits_get_hdu_num(orbit_fptr, &orbit_hdunum) == 1) {
      // this is the primary array, so try to move to the first extension and see if it is a table
      if (fits_movabs_hdu(orbit_fptr, 2, &orbit_hdutype, &status)) break;
    } else {
      // get the HDU type
      if (fits_get_hdu_type(orbit_fptr, &orbit_hdutype, &status)) break;
    }

    // image HDU results in an error message
    if (orbit_hdutype==IMAGE_HDU) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: FITS extension in orbit file '%s' is not a table but an image (HDU number: %d)\n", orbit_filename, orbit_hdunum);
      HD_ERROR_THROW(msg,status);
      break;
    }

    // determine number of rows in the orbit file
    if (fits_get_num_rows(orbit_fptr, &orbit_nrows, &status)) break;

    headas_chat(5, "orbit file '%s' opened ...\n", orbit_filename);



    // delete old attitude file
    remove(attitude_filename);
 
    // create new FITS output file
    if (fits_create_file(&attitude_fptr, attitude_filename, &status)) break;
    headas_chat(5, "attitude file '%s' created ...\n", attitude_filename);

    // variables for the creation of the fits-table (field-type and -format of columns)
    char *ftype[N_ATTITUDE_FIELDS];
    char *fform[N_ATTITUDE_FIELDS];
    char *funit[N_ATTITUDE_FIELDS];
    // create a binary table in the FITS file
    create_attitudetbl_parameter(ftype, fform, funit);
    if (fits_create_tbl(attitude_fptr, BINARY_TBL, 0, N_ATTITUDE_FIELDS, ftype, fform, funit, "eROSITA attitude" , &status)) break;
    /* int fits_create_tbl(fitsfile *fptr, int tbltype, long nrows, int tfields,
       char *ttype[],char *tform[], char *tunit[], char *extname, int *status)*/

    /*
    int count;
    for (count=0; count<N_ATTITUDE_FIELDS; count++) {
      if (ftype[count]) free(ftype[count]);
      if (fform[count]) free(fform[count]);
      if (funit[count]) free(funit[count]);
    }
    */



    // write headers
    // mission headers
    fits_write_key (attitude_fptr, TSTRING, "TELESCOP", "eROSITA", "name of the telescope", &status);
    fits_write_key (attitude_fptr, TSTRING, "MISSION", "SpectrumXGamma", "name of the mission", &status);
    fits_write_key (attitude_fptr, TSTRING, "COMMENT", "DESCRIPT", "attitude file for the eROSITA telescope", &status);

    // date and time headers
    char creation_date[30];
    int timeref;            // is 0, if returned time is in UTC
    fits_get_system_time(creation_date, &timeref, &status);
    fits_write_key(attitude_fptr, TSTRING, "DATE", "2008-03-10","FITS file creation date (yyyy-mm-dd)", &status);

    // if desired by the user, print all program parameters to HISTORY of FITS file (HDU number 1)
    HDpar_stamp(attitude_fptr, 2, &status);



    // check, if errors have occurred on writing the headers
    if (status) break;
    headas_chat(5, "headers written to attitude file '%s' file ...\n", attitude_filename);
    

    // read orbit data from file and calculate attitude data
    headas_chat(5, "reading orbit and calculating attitude data ...\n");
    long row;
    for (row=0; row < orbit_nrows; row++) {
      // read one line of orbit data
      if (get_orbtbl_row(orbit_fptr, row, &time, &r, &v , &status)) break;
      
      // determine the corresponding attitude data
      // time calculation not implemented (TODO)
      strcpy(valtime, "");
      // telescope looks simply outwards (angles in degrees)
      sz0 = normalize_vector(r);
      view_ra = atan2(sz0.y, sz0.x);
      view_dec = asin(sz0.z);
      // determine the position of the sun (according to Meeus p. 163 ff.)
      JD = 2448908.5;
      calc_sun_pos(JD, &sun_ra, &sun_dec);
      // calculate the unit vectors spanning the satellite's intrinsic coordinate system
      sx0 = normalize_vector(v);
      sy0 = normalize_vector(vector_product(sz0, sx0));
      sun_pos = unit_vector(sun_ra, sun_dec);
      // From the sun's position determine roll angle and solar aspect angle.
      // Here the roll angle is the angle around the satellite's z-axis,
      // the solar aspect angle is 0°, if the sun is in the xy-plane of the satellite
      // and 90° if the sun is along the line of sight of the telescope.
      rollangle=atan2(scalar_product(sy0,sun_pos), scalar_product(sx0,sun_pos));
      aspangle=asin(scalar_product(sz0,sun_pos));

//      printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", time, sx0.x, sx0.y, sx0.z, sy0.x, sy0.y, sy0.z, sz0.x, sz0.y, sz0.z);

      // insert new row with calculated attitude data into FITS table
      if ((status=add_attitudetbl_row(attitude_fptr, row, valtime, time, view_ra*180./M_PI, 
				      view_dec*180./M_PI, rollangle*180./M_PI, 
				      aspangle*180./M_PI, status))) break;
//      if (status=add_attitudetbl_row(attitude_fptr, row, valtime, time, 0., 0., 0., 0., status)) break;
    }
    headas_chat(5, "calculated %ld lines of attitude data ...\n", row);
    
  } while(0);  // end of error loop


  // clean up:
  headas_chat(5, "closing files ...\n");
  // close the orbit file
  if(orbit_fptr) fits_close_file(orbit_fptr, &status);
  // close the attitude file
  if(attitude_fptr) fits_close_file(attitude_fptr, &status);

  return(status);
}





//////////////////////////////////////////////////////////////////////////////
// calculates right ascension and declination of the sun to the corresponding
// Julian Day
void calc_sun_pos(
		  double JD,          // Julian Day
		  double *sun_ra,     // right ascension of the sun
		  double *sun_dec     // declination of the sun
		  )
{
  // variables needed for the calculation of the position of the sun
  double JT;                     // Julian century
  double Ls;                     // geometric mean longitude of the sun
  double Ms;                     // mean anomaly of the sun
  double ee;                     // eccentricity of the earth's orbit
  double Cs;                     // sun's equation of center
  double tls;                    // sun's true longitude
  double eps;                    // obliquity of the ecliptic


  // determine the position of the sun (according to Meeus p. 163 ff.)
  JD = 2448908.5;
  JT = (JD-2451545.0)/36525.;
  Ls = 280.46646 + 36000.76983*JT + 0.0003032*pow(JT,2.);
  Ms = 357.52911 + 35999.05029*JT - 0.0001537*pow(JT,2.);
  ee = 0.016708634 - 0.000042037*JT - 0.0000001267*pow(JT,2.);
  Cs = (1.914602 - 0.004817*JT - 0.000014*pow(JT,2.))*sin(Ms*M_PI/180.)
    + (0.019993 - 0.000101*JT)*sin(2.*Ms*M_PI/180.)
    + 0.000289*sin(3.*Ms*M_PI/180.);
  tls = Ls + Cs;
  eps = 23. + 26./60. + 21.448/3600. - 46.8150/3600.*JT - 0.00059/3600.*pow(JT,2.) + 0.001813/3600.*pow(JT,3.);
  *sun_ra = atan2(cos(eps*M_PI/180.)*sin(tls*M_PI/180.), cos(tls*M_PI/180.));
  *sun_dec = asin(sin(eps*M_PI/180.)*sin(tls*M_PI/180.));

}


