/////////////////////////////////////////////////////////////////////////
//
// This program is part of the eROSITA simulation and converts 
// the ROSAT faint source catalog from ASCII to FITS format.
//
/////////////////////////////////////////////////////////////////////////
//
// @author   Christian Schmid
// @data     2008/04
// @param    inputfile - filename of the ASCII input file (RASS-FSC)
// @param    outputfile - filename of the FITS output file
//
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fitsio.h"
#include "pil.h"
#include "headas.h"
#include "headas_error.h"

#include "strftcpy.h"

#define TOOLSUB conv_rosat2fits_main
#include "headas_main.c"

#define LINELENGTH 1024     // maximum linelength in the ASCII file
#define NFIELDS 30          // number of fields in the table (= number of columns)
#define FILENAME_LENGTH 128 // maximum length of filenames
#define MAXMSG 256          // maximum length of an output/error message



// reads all parameters of 'conv_rosat2fits' using PIL
int conv_rosat2fits_getpar(char inputfile[], char outputfile[]);

// does the acutal work: create FITS file with table, open ASCII file, 
// transfer sources, close files, clean up
int conv_rosat2fits_work(const char inputfile[], const char outputfile[]);

// creates the necessary parameters to create the table in the FITS file
void rosat_create_tbl_parameter(char **ftype, char **fform, char **funit);

// parses an ASCII input line and writes the corresponding 
// source data to the FITS table
int conv_rosat2fits_insert_fitsrow(char line[], fitsfile *output_fptr, long row);



////////////////////////////////////////////////////////////////////////////////////////
// main procedure
int conv_rosat2fits_main()
{
  char inputfile[FILENAME_LENGTH];      // name of input file (ASCII file)
  char outputfile[FILENAME_LENGTH];     // name of output file (FITS file)

  int status=0;             // error report status


  // HEATOOLs: register program
  set_toolname("conv_rosat2fits");
  set_toolversion("0.01");


  // read parameters using PIL library
  status = conv_rosat2fits_getpar(inputfile, outputfile);


  if(!status) {
    // perform the actual work: open and create files respectively, transfer source data etc.
    status = conv_rosat2fits_work(inputfile, outputfile);

    headas_chat(5, "finished\n");
  }

  return(status);
}



//////////////////////////////////////////////////////////////////////////////////////////
// reads all parameters of 'conv_rosat2fits' using PIL
int conv_rosat2fits_getpar(
			   char inputfile[],   // filename of the ASCII input file
			   char outputfile[]   // filename of the FITS output file
			   ) 
{
  int status=EXIT_SUCCESS;    // error status
  char msg[MAXMSG];           // buffer for error output messages

  if ((status = PILGetFname("inputfile", inputfile))) {
    sprintf(msg, "Error reading the 'inputfile' parameter!");
    HD_ERROR_THROW(msg,status);
  }

  else if ((status = PILGetFname("outputfile", outputfile))) {
    sprintf(msg, "Error reading the 'outputfile' parameter!");
    HD_ERROR_THROW(msg,status);
  }

  return(status);
}



//////////////////////////////////////////////////////////////////////////////////////////
// does the acutal work: create FITS file with table, open ASCII file, transfer sources, close files, clean up
int conv_rosat2fits_work(
			 const char inputfile[],  // filename of the ASCII input file
			 const char outputfile[]  // filename of the FITS output file
			 )
{
  fitsfile *output_fptr=NULL; // fitsfile pointer to output file
  FILE *input_fptr=NULL;      // ASCII file pointer to input file
  char line[LINELENGTH];      // input buffer for ASCII file

  int count;
  char **ftype=NULL;
  char **fform=NULL;
  char **funit=NULL;

  int status=0;               // error status
  char msg[MAXMSG];           // buffer for error output messages


  do {  // error handling loop (is only run once)

    // open ASCII file with input data (ROSAT all sky survey - faint sources)
    if (!(input_fptr=fopen(inputfile, "r+"))) {
      status = EXIT_FAILURE;
      sprintf(msg, "Error opening the ASCII inputfile '%s'\n", inputfile);
      HD_ERROR_THROW(msg,status);
      break;
    }
    headas_chat(5, "ASCII file '%s' opened ...\n", inputfile);


    // delete old FITS output file
    remove(outputfile);
  
    // create new FITS output file
    if (fits_create_file(&output_fptr, outputfile, &status)) break;
    headas_chat(5, "FITS file '%s' created ...\n", outputfile);

    // variables for the creation of the fits-table (field-type and -format of columns)
    ftype = (char **)malloc(NFIELDS * sizeof(char*));
    fform = (char **)malloc(NFIELDS * sizeof(char*));
    funit = (char **)malloc(NFIELDS * sizeof(char*));
    if ((ftype==NULL) || (fform==NULL) || (funit==NULL)) {
      status=EXIT_FAILURE;
      sprintf(msg, "Error: not enough memory available!\n");
      HD_ERROR_THROW(msg,status);
      return(status);
    } else {
      for(count = 0; count < NFIELDS; count++) {
	ftype[count] = (char *) malloc(6 * sizeof(char));
	fform[count] = (char *) malloc(4 * sizeof(char));
	funit[count] = (char *) malloc(30* sizeof(char));
	if ((ftype[count]==NULL) || (fform[count]==NULL) || (funit[count]==NULL)) {
	  status=EXIT_FAILURE;
	  sprintf(msg, "Error: not enough memory available!\n");
	  HD_ERROR_THROW(msg,status);
	  return(status);
	}
      }
    }

    // create a binary table in the FITS file
    rosat_create_tbl_parameter(ftype, fform, funit);
    if (fits_create_tbl(output_fptr, BINARY_TBL, 0, NFIELDS, ftype, fform, funit, "RASS_FSC" , &status)) break;
    /* int fits_create_tbl(fitsfile *fptr, int tbltype, long nrows, int tfields,
       char *ttype[],char *tform[], char *tunit[], char *extname, int *status)*/
    
    // write headers
    fits_write_key (output_fptr, TSTRING, "TELESCOP", "ROSAT", "name of the telescope", &status);
    fits_write_key (output_fptr, TSTRING, "MISSION", "RASS", "rosat all sky survey", &status);
    fits_write_key (output_fptr, TSTRING, "COMMENT", "DESCRIPT", "faint sources of the ROSAT all sky survey", &status);

    // if desired by the user, print all program parameters to HISTORY of FITS file (HDU number 1)
    HDpar_stamp(output_fptr, 2, &status);

    // check, if errors have occurred on writing the headers
    if (status) break;
    headas_chat(5, "headers written into FITS file '%s' ...\n", outputfile);
    

    // read ASCII-data from file
    headas_chat(5, "transferring source data ...\n");
    long row = 0;
    while (fgets(line, LINELENGTH, input_fptr)) {
      // check whether it is a data line with "1" at the beginning (comment lines start with "!")
      if (line[0] == '1') {
	// insert new row with source data from ASCII file to FITS table
	if ((status=conv_rosat2fits_insert_fitsrow(line, output_fptr, ++row))) break;
      }
    }
    headas_chat(3, "%ld sources transferred ...\n", row);
    
  } while(0);  // end of error loop



  // clean up:
  headas_chat(5, "closing files ...\n");

  // close the ASCII input file
  if(input_fptr) fclose(input_fptr);
  // close the FITS file
  if(output_fptr) fits_close_file(output_fptr, &status);

  // free memory
  for (count=0; count<NFIELDS; count++) {
    if (ftype[count]) free (ftype[count]);
    if (fform[count]) free (fform[count]);
    if (funit[count]) free (funit[count]);
  }

  return(status);
}




////////////////////////////////////////////////////////////////////////////////////////////////
// creates the necessary data for the FITS-table layout
void rosat_create_tbl_parameter(char **ftype, char **fform, char **funit) {
    /* determine field types of the table in the FITS file */

    /* 1. : source name */
    strcpy(ftype[0], "SrcNam");
    strcpy(fform[0], "21A");
    strcpy(funit[0], "");

    /* 2. : right ascension */
    strcpy(ftype[1], "R_A_");
    strcpy(fform[1], "D");
    strcpy(funit[1], "decimal degrees");

    /* 3. : declination */
    strcpy(ftype[2], "Dec_");
    strcpy(fform[2], "D");
    strcpy(funit[2], "decimal degrees");

    /* 4. : total position error */
    strcpy(ftype[3], "_p");
    strcpy(fform[3], "I");
    strcpy(funit[3], "arcsec");


    /* 5. : Flags */
    strcpy(ftype[4], "Flags");
    strcpy(fform[4], "5A");
    strcpy(funit[4], "");

    /* 6. : Flags2 */
    strcpy(ftype[5], "Flg2");
    strcpy(fform[5], "4A");
    strcpy(funit[5], "");

    /* 7. : source count rate */
    strcpy(ftype[6], "src_cps");
    strcpy(fform[6], "E");
    strcpy(funit[6], "counts/sec");

    /* 8. : error in source count rate */
    strcpy(ftype[7], "_src_cps");
    strcpy(fform[7], "E");
    strcpy(funit[7], "counts/sec");

    /* 9. : background count rate */
    strcpy(ftype[8], "bgr_cpsa");
    strcpy(fform[8], "E");
    strcpy(funit[8], "counts/sec/arcmin**2");

    /* 10. : exposure time */
    strcpy(ftype[9], "exp");
    strcpy(fform[9], "J");
    strcpy(funit[9], "sec");

    /* 11. : hardness ratio 1 */
    strcpy(ftype[10], "hr1");
    strcpy(fform[10], "E");
    strcpy(funit[10], "");

    /* 12. : error inn hardness ratio 1 */
    strcpy(ftype[11], "_hr1");
    strcpy(fform[11], "E");
    strcpy(funit[11], "");

    /* 13. : hardness ratio 2 */
    strcpy(ftype[12], "hr2");
    strcpy(fform[12], "E");
    strcpy(funit[12], "");

    /* 14. : error inn hardness ratio 2 */
    strcpy(ftype[13], "_hr2");
    strcpy(fform[13], "E");
    strcpy(funit[13], "");

    /* 15. : source extent */
    strcpy(ftype[14], "ext");
    strcpy(fform[14], "I");
    strcpy(funit[14], "arcsec");

    /* 16. : likelihood of source extent */
    strcpy(ftype[15], "extl");
    strcpy(fform[15], "I");
    strcpy(funit[15], "");

    /* 17. : likelihood of source detection */
    strcpy(ftype[16], "srcl");
    strcpy(fform[16], "I");
    strcpy(funit[16], "");

    /* 18. : extraction radius */
    strcpy(ftype[17], "extr");
    strcpy(fform[17], "I");
    strcpy(funit[17], "arcsec");

    /* 19. : priority flag */
    strcpy(ftype[18], "PriFlg");
    strcpy(fform[18], "6A");
    strcpy(funit[18], "");

    /* 20. : PHA range with highest detection likelihood */
    strcpy(ftype[19], "E");
    strcpy(fform[19], "A");
    strcpy(funit[19], "");

    /* 21. : vignetting factor */
    strcpy(ftype[20], "vigf");
    strcpy(fform[20], "E");
    strcpy(funit[20], "");

    /* 22. : date when source was included */
    strcpy(ftype[21], "orgdat");
    strcpy(fform[21], "6A");
    strcpy(funit[21], "");

    /* 23. : date when source properties were changed */
    strcpy(ftype[22], "moddat");
    strcpy(fform[22], "6A");
    strcpy(funit[22], "");

    /* 24. : number of identification candidates */
    strcpy(ftype[23], "_id");
    strcpy(fform[23], "I");
    strcpy(funit[23], "");

    /* 25. : identification number of SASS field */
    strcpy(ftype[24], "field_id");
    strcpy(fform[24], "J");
    strcpy(funit[24], "");

    /* 26. : SASS source number */
    strcpy(ftype[25], "src_");
    strcpy(fform[25], "I");
    strcpy(funit[25], "");

    /* 27. : number of nearby RASS detections */
    strcpy(ftype[26], "rct");
    strcpy(fform[26], "I");
    strcpy(funit[26], "");

    /* 28. : start time of observation */
    strcpy(ftype[27], "itb");
    strcpy(fform[27], "9A");
    strcpy(funit[27], "yymmdd.ff");

    /* 29. : end time of observation */
    strcpy(ftype[28], "ite");
    strcpy(fform[28], "9A");
    strcpy(funit[28], "yymmdd.ff");

    /* 30. : reliability of source detection */
    strcpy(ftype[29], "rl");
    strcpy(fform[29], "I");
    strcpy(funit[29], "");

}



void rosat_create_tbl_parameter2(char *ftype[NFIELDS], char *fform[NFIELDS], char *funit[NFIELDS]) {
    int counter;

    /* determine field types of the table in the FITS file */
    for(counter = 0; counter < NFIELDS; counter++) {
      ftype[counter] = (char *) malloc(6 * sizeof(char));
      fform[counter] = (char *) malloc(4 * sizeof(char));
      funit[counter] = (char *) malloc(30* sizeof(char));
    }

    /* 1. : source name */
    ftype[0] = "SrcNam";
    fform[0] = "21A";
    funit[0] = "";

    /* 2. : right ascension */
    ftype[1] = "R_A_";
    fform[1] = "D";
    funit[1] = "decimal degrees";

    /* 3. : declination */
    ftype[2] = "Dec_";
    fform[2] = "D";
    funit[2] = "decimal degrees";

    /* 4. : total position error */
    ftype[3] = "_p";
    fform[3] = "I";
    funit[3] = "arcsec";

    /* 5. : Flags */
    ftype[4] = "Flags";
    fform[4] = "5A";
    funit[4] = "";

    /* 6. : Flags2 */
    ftype[5] = "Flg2";
    fform[5] = "4A";
    funit[5] = "";

    /* 7. : source count rate */
    ftype[6] = "src_cps";
    fform[6] = "E";
    funit[6] = "counts/sec";

    /* 8. : error in source count rate */
    ftype[7] = "_src_cps";
    fform[7] = "E";
    funit[7] = "counts/sec";

    /* 9. : background count rate */
    ftype[8] = "bgr_cpsa";
    fform[8] = "E";
    funit[8] = "counts/sec/arcmin**2";

    /* 10. : exposure time */
    ftype[9] = "exp";
    fform[9] = "J";
    funit[9] = "sec";

    /* 11. : hardness ratio 1 */
    ftype[10] = "hr1";
    fform[10] = "E";
    funit[10] = "";

    /* 12. : error inn hardness ratio 1 */
    ftype[11] = "_hr1";
    fform[11] = "E";
    funit[11] = "";

    /* 13. : hardness ratio 2 */
    ftype[12] = "hr2";
    fform[12] = "E";
    funit[12] = "";

    /* 14. : error inn hardness ratio 2 */
    ftype[13] = "_hr2";
    fform[13] = "E";
    funit[13] = "";

    /* 15. : source extent */
    ftype[14] = "ext";
    fform[14] = "I";
    funit[14] = "arcsec";

    /* 16. : likelihood of source extent */
    ftype[15] = "extl";
    fform[15] = "I";
    funit[15] = "";

    /* 17. : likelihood of source detection */
    ftype[16] = "srcl";
    fform[16] = "I";
    funit[16] = "";

    /* 18. : extraction radius */
    ftype[17] = "extr";
    fform[17] = "I";
    funit[17] = "arcsec";

    /* 19. : priority flag */
    ftype[18] = "PriFlg";
    fform[18] = "6A";
    funit[18] = "";

    /* 20. : PHA range with highest detection likelihood */
    ftype[19] = "E";
    fform[19] = "A";
    funit[19] = "";

    /* 21. : vignetting factor */
    ftype[20] = "vigf";
    fform[20] = "E";
    funit[20] = "";

    /* 22. : date when source was included */
    ftype[21] = "orgdat";
    fform[21] = "6A";
    funit[21] = "";

    /* 23. : date when source properties were changed */
    ftype[22] = "moddat";
    fform[22] = "6A";
    funit[22] = "";

    /* 24. : number of identification candidates */
    ftype[23] = "_id";
    fform[23] = "I";
    funit[23] = "";

    /* 25. : identification number of SASS field */
    ftype[24] = "field_id";
    fform[24] = "J";
    funit[24] = "";

    /* 26. : SASS source number */
    ftype[25] = "src_";
    fform[25] = "I";
    funit[25] = "";

    /* 27. : number of nearby RASS detections */
    ftype[26] = "rct";
    fform[26] = "I";
    funit[26] = "";

    /* 28. : start time of observation */
    ftype[27] = "itb";
    fform[27] = "9A";
    funit[27] = "yymmdd.ff";

    /* 29. : end time of observation */
    ftype[28] = "ite";
    fform[28] = "9A";
    funit[28] = "yymmdd.ff";

    /* 30. : reliability of source detection */
    ftype[29] = "rl";
    fform[29] = "I";
    funit[29] = "";

}



////////////////////////////////////////////////////////////////////////////////////////////////
// parses an ASCII input line and writes the corresponding source data to the FITS table
int conv_rosat2fits_insert_fitsrow(
				   char line[],           // ASCII textline containing the source data
				   fitsfile *output_fptr, // FITS file pointer to output file
				   long row               // actual row in the FITS file
				   )
{
  // table field buffers
  char *cbuffer[1];         // string buffer
  double dbuffer[1];        // double buffer
  float fbuffer[1];         // float buffer
  long lbuffer[1];
  int ibuffer[1];

  int status=0;             // error status
  char msg[MAXMSG];         // buffer for error output message

  do {  // beginning of error handling loop (is only run once)
    // get memory
    cbuffer[0]=malloc(21*sizeof(char));
    if(!cbuffer[0]) {
      status = -1;
      sprintf(msg, "Error allocating memory");
      HD_ERROR_THROW(msg,status);
      break;
    }


    // insert new row to binary FITS table
    if (fits_insert_rows(output_fptr, row-1, 1, &status)) break;


    // parse table data
    // source name
    strncpy(cbuffer[0],line,21);
    // write column-entry
    if (fits_write_col(output_fptr, TSTRING, 1, row, 1, 1, cbuffer, &status)) break;
    // fits_write_col(fitsfile *fptr, int datatype, int colnum, long firstrow,
    //         long firstelem, long nelements, void *array, int *status)

    // R.A.
    strftcpy(cbuffer[0],line,22,10);
    dbuffer[0] = atof(cbuffer[0]);    // convert to floating point number
    if (fits_write_col(output_fptr, TDOUBLE, 2, row, 1, 1, dbuffer, &status)) break;

    // declination
    strftcpy(cbuffer[0],line,33,9);
    dbuffer[0] = atof(cbuffer[0]);    // convert to floating point number
    if (fits_write_col(output_fptr, TDOUBLE, 3, row, 1, 1, dbuffer, &status)) break;

    // position error
    strftcpy(cbuffer[0],line,43,3);
    ibuffer[0] = atoi(cbuffer[0]);    // convert to integer number
    if (fits_write_col(output_fptr, TINT, 4, row, 1, 1, ibuffer, &status)) break;

    // flags
    strftcpy(cbuffer[0],line,47,5);
    if (fits_write_col(output_fptr, TSTRING, 5, row, 1, 1, cbuffer, &status)) break;

    // flags2
    strftcpy(cbuffer[0],line,53,4);
    if (fits_write_col(output_fptr, TSTRING, 6, row, 1, 1, cbuffer, &status)) break;



    // source count rate
    strftcpy(cbuffer[0],line,59,9);
    fbuffer[0] = atof(cbuffer[0]);

    // convert ROSAT count rate via ROSAT energy flux to eROSITA photon rate (count rate is lower due to PSF losses)
    fbuffer[0] = fbuffer[0] * 5.6e-12 *       1e-7     / 1.6022e-16 /  log(2.34/0.12) *         2471 *  
    //           ROSAT cps    cps->erg/s/cm²  erg->J     J->keV        ln(E1/E0)                A_eROSITA (cm^2)
    //                                                                 observed energy band     light collecting area of eROSITA
    //                                                                 2.-0.5 ?        
       (1./0.5 - 1./10.)     ;
    // energy band of eROSITA

    if (fits_write_col(output_fptr, TFLOAT, 7, row, 1, 1, fbuffer, &status)) break;



    // error in source count rate
    strftcpy(cbuffer[0],line,69,9);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 8, row, 1, 1, fbuffer, &status)) break;

    // background count rate
    strftcpy(cbuffer[0],line,79,9);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 9, row, 1, 1, fbuffer, &status)) break;

    // exposure time
    strftcpy(cbuffer[0],line,89,6);
    lbuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TLONG, 10, row, 1, 1, lbuffer, &status)) break;

    // hardness ratio 1
    strftcpy(cbuffer[0],line,96,5);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 11, row, 1, 1, fbuffer, &status)) break;

    // error in hardness ratio 1
    strftcpy(cbuffer[0],line,102,4);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 12, row, 1, 1, fbuffer, &status)) break;

    // hardness ratio 2
    strftcpy(cbuffer[0],line,107,5);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 13, row, 1, 1, fbuffer, &status)) break;

    // error in hardness ratio 2
    strftcpy(cbuffer[0],line,113,4);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 14, row, 1, 1, fbuffer, &status)) break;

    strftcpy(cbuffer[0],line,118,5);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 15, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,124,4);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 16, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,129,4);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 17, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,134,5);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 18, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,140,6);
    if (fits_write_col(output_fptr, TSTRING, 19, row, 1, 1, cbuffer, &status)) break;
  
    strftcpy(cbuffer[0],line,146,1);
    if (fits_write_col(output_fptr, TSTRING, 20, row, 1, 1, cbuffer, &status)) break;

    strftcpy(cbuffer[0],line,148,4);
    fbuffer[0] = atof(cbuffer[0]);
    if (fits_write_col(output_fptr, TFLOAT, 21, row, 1, 1, fbuffer, &status)) break;

    strftcpy(cbuffer[0],line,153,6);
    if (fits_write_col(output_fptr, TSTRING, 22, row, 1, 1, cbuffer, &status)) break;

    strftcpy(cbuffer[0],line,160,6);
    if (fits_write_col(output_fptr, TSTRING, 23, row, 1, 1, cbuffer, &status)) break;

    strftcpy(cbuffer[0],line,167,4);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 24, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,172,8);
    lbuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TLONG, 25, row, 1, 1, lbuffer, &status)) break;

    strftcpy(cbuffer[0],line,181,4);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 26, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,186,3);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 27, row, 1, 1, ibuffer, &status)) break;

    strftcpy(cbuffer[0],line,190,9);
    if (fits_write_col(output_fptr, TSTRING, 28, row, 1, 1, cbuffer, &status)) break;

    strftcpy(cbuffer[0],line,200,9);
    if (fits_write_col(output_fptr, TSTRING, 29, row, 1, 1, cbuffer, &status)) break;

    strftcpy(cbuffer[0],line,210,2);
    ibuffer[0] = atoi(cbuffer[0]);
    if (fits_write_col(output_fptr, TINT, 30, row, 1, 1, ibuffer, &status)) break;

  } while (0);  // end of error loop

  //clean up:
  // free memory
  if(cbuffer[0]) free(cbuffer[0]);

  return(status);
}
