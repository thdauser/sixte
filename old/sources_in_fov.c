#include <stdio.h>
#include <math.h>
#include "fitsio.h"
#include "pil.h"
#include "vector.h"
#include "fits_ctlg.h"
#include "check_fov.h"

// Usage (PIL interface):
// ./sources_in_fov <N_inputfiles> <rosat-fitsfile> <rnd_sources-fitsfile> ... <outputfile> <RA> <Dec> ...
//
// The program scans the input-files and selects sources that are inside the specified FOV.
// These are written to the output file.
// Standard shape for the FOV is a circle.
// In order to use a rectangle, you have to remove the "/*" and
// "*/" comments.


int main(int argc, char** argv) {
  // constants
/*   const struct vector h1 = {0.0, 0.0, 1.0};   // normal vector 1
  const struct vector h2 = {0.0, 1.0, 0.0};   // normal vector 2
  const double dec_max = 0.343333 *M_PI/180;  // half-height of the FOV
  const double rasc_max = 0.343333 *M_PI/180; // half-width of the FOV
  const double sin_dec_max = sin(dec_max);
  const double sin_rasc_max = sin(rasc_max);
        // cos(sqrt(dec_max*dec_max + rasc_max*rasc_max)); */

  // field-indices of r.a., dec. etc. in the different FITS input file tables
  const int columns[2][3] = {{2,3,7}, {1,2,3}};

  // variables
  fitsfile *fptr, *fptr2;   // pointer to the FITS file
  int fitsfile_counter;     // counter for the different FITS input files
  int hdunum,hdutype,fitsstatus=0;  // need to access the HDUs in the FITS file
  char *ftype[N_SOURCE_FIELDS];     // table format
  char *fform[N_SOURCE_FIELDS];     // table format
  char *funit[N_SOURCE_FIELDS];     // units of the table entries
  long counter,nrows,nsources=0,output_row=0;

  // parameters:
  int n_inputfiles;         // number of input files
  int input_counter;
  char parname[12];
  char inputfile[5][100];   // input files (rosat-fsc, rnd_sources,...)
  char outputfile[100];     // output file (FITS file)

  // telescope data (input parameters)
  struct vector x0;         // direction of the telescop axis
  double diameter;          // diameter of the FOV (in arcmin ["])
  double min_align;         // minimum cos of for sources inside the FOV
                            //      angle(x0,source) <= 1/2 * diameter

  // source-data
  double rasc,dec;          // right ascension, declination
  float countrate;          // countrate
  struct vector x;            // vector pointing to the examined source
  long nhits = 0, nmiss = 0;  // number of sources inside and outside the FOV respectively


  // read parameters using PIL library
  // initialize PIL
  int pilread = PILInit(argc, argv);
  // get the number of input-files
  pilread = PILGetInt("n_inputfiles", &n_inputfiles);
  // get the filenames of the individual input-files
  for(input_counter=0; input_counter<n_inputfiles; input_counter++) {
    sprintf(parname,"inputfile%d",input_counter+1);
    pilread = PILGetFname(parname, inputfile[input_counter]);
  }
  // get the filename of the output-file
  pilread = PILGetFname("outputfile", outputfile);
  // read right ascension and declination of the telescop pointing direction
  pilread = PILGetReal("ra", &rasc);
  pilread = PILGetReal("dec", &dec);
  // read the diameter of the FOV
  pilread = PILGetReal("diameter", &diameter);
  // close PIL interface
  PILClose(PIL_OK);
  // check, whether all parameters were read correctly:
  if (pilread < 0) {
    printf("PIL - parameter input failed: %s\n", PIL_err_handler(pilread));
    return(-1);
  }



  // create unit vector pointing in the direction of the telescop:
  unit_vector(rasc*M_PI/180,dec*M_PI/180);

  // calculate the minimum cos-value for sources inside the FOV:
  min_align = cos(diameter/120 * M_PI/180);   // angle(x0,source) <= 1/2 * diameter



  // delete the old output FITS file
  remove(outputfile);

  // create new FITS file for output data
  if (fits_create_file(&fptr2, outputfile, &fitsstatus)) {
    printf("Error: could not create output file %s\n", outputfile);
    return(fitsstatus);
  }

  // create a binary table in the output FITS file
  create_srctbl_parameters(ftype,fform,funit);
  if(fits_create_tbl(fptr2, BINARY_TBL, 0, N_SOURCE_FIELDS, ftype, fform, funit, "RASS_FSC" , &fitsstatus)) {
    printf("Error: could not create data table in output file %s\n", outputfile);
    return(fitsstatus);
  }

  // write description data into the header of the FITS-file
  fits_write_key(fptr2, TSTRING, "COMMENT", "VISIBLE SOURCES", "sources from catalogue, which lie inside FOV",&fitsstatus);
  fits_write_key(fptr2, TSTRING, "COMMENT", "displayed data:", "", &fitsstatus);
  fits_write_key(fptr2, TSTRING, "COMMENT", "r.a., dec, cps", "right ascension, declination, countrate", &fitsstatus);


  for(fitsfile_counter=0; fitsfile_counter < n_inputfiles; fitsfile_counter++) {
    // open the source catalogue (FITS-file)
    if(!fits_open_file(&fptr, inputfile[fitsfile_counter], READONLY, &fitsstatus)) {
      if (fits_get_hdu_num(fptr, &hdunum) == 1) {
        // this is the primary array
        // try to move to the first extension and see if it is a table
        fits_movabs_hdu(fptr, 2, &hdutype, &fitsstatus);
      } else {
        fits_get_hdu_type(fptr, &hdutype, &fitsstatus); // Get the HDU type
      }

      if (hdutype == IMAGE_HDU) {
        printf("Error: extension is not a table (default: second extension)\n");
      } else {
        fits_get_num_rows(fptr, &nrows, &fitsstatus);
        nsources += nrows;
        // read all table rows after each other
        for (counter=0; counter<nrows; counter++) {
          // read table row
          if(get_srctbl_row(fptr, counter, &columns[fitsfile_counter][0], &rasc, &dec, &countrate, &fitsstatus)) {
            printf("Error: reading data\n");
            return(fitsstatus);
          } else {
            x = unit_vector(rasc*M_PI/180, dec*M_PI/180);
            if(check_fov(x, x0, /* h1, h2, sin_dec_max, sin_rasc_max, */ min_align) > 0) {
              // source is outside the FOV
              nmiss++;
            } else {
              // source is inside the FOV
              nhits++;
              // add source to the list of visible sources
              add_srctbl_row(fptr2, output_row++, rasc, dec, countrate, &fitsstatus);
            }
          }
        }
      }

      // close the FITS file (source catalogue)
      fits_close_file(fptr, &fitsstatus);
    } else {
      printf("Error: could not open input file '%s'\n", inputfile[fitsfile_counter]);
    }
  }

  // close the FITS file (outputfile)
  fits_close_file(fptr2, &fitsstatus);

  // print any error messages
  if (fitsstatus) {
    fits_report_error(stderr, fitsstatus);
  }

  // print the number of sources inside the FOV
  printf("number of visible sources: %ld of totally %ld\n", nhits, nsources);

  return(fitsstatus);
}
