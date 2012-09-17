#include <stdio.h>
#include "fitsio.h"
#include "headas_error.h"

#ifndef ERR
#define ERR HD_ERR_MAP_STATIC_ENTRY
#endif

#define MSG NULL

static hd_error_map_entry sCFITSIOerrors[] = {
  ERR(SAME_FILE, MSG),           /* input and output files are the same */
  ERR(TOO_MANY_FILES, MSG),      /* tried to open too many FITS files */
  ERR(FILE_NOT_OPENED, MSG),     /* could not open the named file */
  ERR(FILE_NOT_CREATED, MSG),    /* could not create the named file */
  ERR(WRITE_ERROR, MSG),         /* error writing to FITS file */
  ERR(END_OF_FILE, MSG),         /* tried to move past end of file */
  ERR(READ_ERROR, MSG),          /* error reading from FITS file */
  ERR(FILE_NOT_CLOSED, MSG),     /* could not close the file */
  ERR(ARRAY_TOO_BIG, MSG),       /* array dimensions exceed internal limit */
  ERR(READONLY_FILE, MSG),       /* Cannot write to readonly file */
  ERR(MEMORY_ALLOCATION, MSG),   /* Could not allocate memory */
  ERR(BAD_FILEPTR, MSG),         /* invalid fitsfile pointer */
  ERR(NULL_INPUT_PTR, MSG),      /* NULL input pointer to routine */
  ERR(SEEK_ERROR, MSG),          /* error seeking position in file */

  ERR(BAD_URL_PREFIX, MSG),      /* invalid URL prefix on file name */
  ERR(TOO_MANY_DRIVERS, MSG),    /* tried to register too many IO drivers */
  ERR(DRIVER_INIT_FAILED, MSG),   /* driver initialization failed */
  ERR(NO_MATCHING_DRIVER, MSG),   /* matching driver is not registered */
  ERR(URL_PARSE_ERROR, MSG),      /* failed to parse input file URL */
  ERR(RANGE_PARSE_ERROR, MSG),    /* failed to parse input file URL */

  ERR(SHARED_ERRBASE, MSG),
  ERR(SHARED_BADARG, MSG),
  ERR(SHARED_NULPTR, MSG),
  ERR(SHARED_TABFULL, MSG),
  ERR(SHARED_NOTINIT, MSG),
  ERR(SHARED_IPCERR, MSG),
  ERR(SHARED_NOMEM, MSG),
  ERR(SHARED_AGAIN, MSG),
  ERR(SHARED_NOFILE, MSG),
  ERR(SHARED_NORESIZE, MSG),

  ERR(HEADER_NOT_EMPTY, MSG),    /* header already contains keywords */
  ERR(KEY_NO_EXIST, MSG),        /* keyword not found in header */
  ERR(KEY_OUT_BOUNDS, MSG),      /* keyword record number is out of bounds */
  ERR(VALUE_UNDEFINED, MSG),     /* keyword value field is blank */
  ERR(NO_QUOTE, MSG),            /* string is missing the closing quote */
  ERR(BAD_KEYCHAR, MSG),         /* illegal character in keyword name or card */
  ERR(BAD_ORDER, MSG),           /* required keywords out of order */
  ERR(NOT_POS_INT, MSG),         /* keyword value is not a positive integer */
  ERR(NO_END, MSG),              /* couldn't find END keyword */
  ERR(BAD_BITPIX, MSG),          /* illegal BITPIX keyword value*/
  ERR(BAD_NAXIS, MSG),           /* illegal NAXIS keyword value */
  ERR(BAD_NAXES, MSG),           /* illegal NAXISn keyword value */
  ERR(BAD_PCOUNT, MSG),          /* illegal PCOUNT keyword value */
  ERR(BAD_GCOUNT, MSG),          /* illegal GCOUNT keyword value */
  ERR(BAD_TFIELDS, MSG),         /* illegal TFIELDS keyword value */
  ERR(NEG_WIDTH, MSG),           /* negative table row size */
  ERR(NEG_ROWS, MSG),            /* negative number of rows in table */
  ERR(COL_NOT_FOUND, MSG),       /* column with this name not found in table */
  ERR(BAD_SIMPLE, MSG),          /* illegal value of SIMPLE keyword  */
  ERR(NO_SIMPLE, MSG),           /* Primary array doesn't start with SIMPLE */
  ERR(NO_BITPIX, MSG),           /* Second keyword not BITPIX */
  ERR(NO_NAXIS, MSG),            /* Third keyword not NAXIS */
  ERR(NO_NAXES, MSG),            /* Couldn't find all the NAXISn keywords */
  ERR(NO_XTENSION, MSG),         /* HDU doesn't start with XTENSION keyword */
  ERR(NOT_ATABLE, MSG),          /* the CHDU is not an ASCII table extension */
  ERR(NOT_BTABLE, MSG),          /* the CHDU is not a binary table extension */
  ERR(NO_PCOUNT, MSG),           /* couldn't find PCOUNT keyword */
  ERR(NO_GCOUNT, MSG),           /* couldn't find GCOUNT keyword */
  ERR(NO_TFIELDS, MSG),          /* couldn't find TFIELDS keyword */
  ERR(NO_TBCOL, MSG),            /* couldn't find TBCOLn keyword */
  ERR(NO_TFORM, MSG),            /* couldn't find TFORMn keyword */
  ERR(NOT_IMAGE, MSG),           /* the CHDU is not an IMAGE extension */
  ERR(BAD_TBCOL, MSG),           /* TBCOLn keyword value < 0 or > rowlength */
  ERR(NOT_TABLE, MSG),           /* the CHDU is not a table */
  ERR(COL_TOO_WIDE, MSG),        /* column is too wide to fit in table */
  ERR(COL_NOT_UNIQUE, MSG),      /* more than 1 column name matches template */
  ERR(BAD_ROW_WIDTH, MSG),       /* sum of column widths not = NAXIS1 */
  ERR(UNKNOWN_EXT, MSG),         /* unrecognizable FITS extension type */
  ERR(UNKNOWN_REC, MSG),         /* unrecognizable FITS record */
  ERR(END_JUNK, MSG),            /* END keyword is not blank */
  ERR(BAD_HEADER_FILL, MSG),     /* Header fill area not blank */
  ERR(BAD_DATA_FILL, MSG),       /* Data fill area not blank or zero */
  ERR(BAD_TFORM, MSG),           /* illegal TFORM format code */
  ERR(BAD_TFORM_DTYPE, MSG),     /* unrecognizable TFORM datatype code */
  ERR(BAD_TDIM, MSG),            /* illegal TDIMn keyword value */
  ERR(BAD_HEAP_PTR, MSG),        /* invalid BINTABLE heap address */

  ERR(BAD_HDU_NUM, MSG),         /* HDU number < 1 or > MAXHDU */
  ERR(BAD_COL_NUM, MSG),         /* column number < 1 or > tfields */
  ERR(NEG_FILE_POS, MSG),        /* tried to move before beginning of file  */
  ERR(NEG_BYTES, MSG),           /* tried to read or write negative bytes */
  ERR(BAD_ROW_NUM, MSG),         /* illegal starting row number in table */
  ERR(BAD_ELEM_NUM, MSG),        /* illegal starting element number in vector */
  ERR(NOT_ASCII_COL, MSG),       /* this is not an ASCII string column */
  ERR(NOT_LOGICAL_COL, MSG),     /* this is not a logical datatype column */
  ERR(BAD_ATABLE_FORMAT, MSG),   /* ASCII table column has wrong format */
  ERR(BAD_BTABLE_FORMAT, MSG),   /* Binary table column has wrong format */
  ERR(NO_NULL, MSG),             /* null value has not been defined */
  ERR(NOT_VARI_LEN, MSG),        /* this is not a variable length column */
  ERR(BAD_DIMEN, MSG),           /* illegal number of dimensions in array */
  ERR(BAD_PIX_NUM, MSG),         /* first pixel number greater than last pixel */
  ERR(ZERO_SCALE, MSG),          /* illegal BSCALE or TSCALn keyword = 0 */
  ERR(NEG_AXIS, MSG),            /* illegal axis length < 1 */

  ERR(NOT_GROUP_TABLE, MSG),         
  ERR(HDU_ALREADY_MEMBER, MSG),      
  ERR(MEMBER_NOT_FOUND, MSG),        
  ERR(GROUP_NOT_FOUND, MSG),         
  ERR(BAD_GROUP_ID, MSG),            
  ERR(TOO_MANY_HDUS_TRACKED, MSG),   
  ERR(HDU_ALREADY_TRACKED, MSG),     
  ERR(BAD_OPTION, MSG),              
  ERR(IDENTICAL_POINTERS, MSG),      
  ERR(BAD_GROUP_ATTACH, MSG),        
  ERR(BAD_GROUP_DETACH, MSG),        

  ERR(BAD_I2C, MSG),             /* bad int to formatted string conversion */
  ERR(BAD_F2C, MSG),             /* bad float to formatted string conversion */
  ERR(BAD_INTKEY, MSG),          /* can't interprete keyword value as integer */
  ERR(BAD_LOGICALKEY, MSG),      /* can't interprete keyword value as logical */
  ERR(BAD_FLOATKEY, MSG),        /* can't interprete keyword value as float */
  ERR(BAD_DOUBLEKEY, MSG),       /* can't interprete keyword value as double */
  ERR(BAD_C2I, MSG),             /* bad formatted string to int conversion */
  ERR(BAD_C2F, MSG),             /* bad formatted string to float conversion */
  ERR(BAD_C2D, MSG),             /* bad formatted string to double conversion */
  ERR(BAD_DATATYPE, MSG),        /* bad keyword datatype code */
  ERR(BAD_DECIM, MSG),           /* bad number of decimal places specified */
  ERR(NUM_OVERFLOW, MSG),        /* overflow during datatype conversion */

  ERR(DATA_COMPRESSION_ERR, MSG),   /* error in imcompress routines */
  ERR(DATA_DECOMPRESSION_ERR, MSG),  /* error in imcompress routines */
  ERR(NO_COMPRESSED_TILE, MSG),   /* compressed tile doesn't exist */

  ERR(BAD_DATE, MSG),            /* error in date or time conversion */

  ERR(PARSE_SYNTAX_ERR, MSG),    /* syntax error in parser expression */
  ERR(PARSE_BAD_TYPE, MSG),      /* expression did not evaluate to desired type */
  ERR(PARSE_LRG_VECTOR, MSG),    /* vector result too large to return in array */
  ERR(PARSE_NO_OUTPUT, MSG),     /* data parser failed not sent an out column */
  ERR(PARSE_BAD_COL, MSG),       /* bad data encounter while parsing column */
  ERR(PARSE_BAD_OUTPUT, MSG),    /* Output file not of proper type          */

  ERR(ANGLE_TOO_BIG, MSG),       /* celestial angle too large for projection */
  ERR(BAD_WCS_VAL, MSG),         /* bad celestial coordinate or pixel value */
  ERR(WCS_ERROR, MSG),           /* error in celestial coordinate calculation */
  ERR(BAD_WCS_PROJ, MSG),        /* unsupported type of celestial projection */
  ERR(NO_WCS_KEY, MSG),          /* celestial coordinate keywords not found */
  ERR(APPROX_WCS_KEY, MSG),      /* approximate WCS keywords were calculated */

  ERR(NO_CLOSE_ERROR, MSG),      /* special value used internally to switch off */
                                    /* the error message from ffclos and ffchdu */

/*-------, MSG), following error codes are used in the grparser.c file -----------*/
  ERR(NGP_ERRBASE, MSG), /* base chosen so not to interfere with CFITSIO */

  ERR(NGP_NO_MEMORY, MSG),  /* malloc failed */
  ERR(NGP_READ_ERR, MSG),  /* read error from file */
  ERR(NGP_NUL_PTR, MSG),  /* null pointer passed as argument */
  ERR(NGP_EMPTY_CURLINE, MSG), /* line read seems to be empty */
  ERR(NGP_UNREAD_QUEUE_FULL, MSG), /* cannot unread more then 1 line (or single line twice) */
  ERR(NGP_INC_NESTING, MSG),  /* too deep include file nesting (inf. loop ?) */
  ERR(NGP_ERR_FOPEN, MSG),  /* fopen() failed, cannot open file */
  ERR(NGP_EOF, MSG),   /* end of file encountered */
  ERR(NGP_BAD_ARG, MSG),  /* bad arguments passed */
  ERR(NGP_TOKEN_NOT_EXPECT, MSG) /* token not expected here */
};

static hd_error_map sCFITSIOerrorMap = HD_ERR_MAP_INIT(sCFITSIOerrors);

#include "pil_error.h"

#define PILMSG NULL

static hd_error_map_entry sPILerrors[] = {

  ERR(PIL_OK, PILMSG),

  ERR(PIL_ERR_BASE, PILMSG),
  ERR(PIL_ERR_MAX_IDX, PILMSG),

  ERR(PIL_NUL_PTR, PILMSG),
  ERR(PIL_BAD_ARG, PILMSG),
  ERR(PIL_NO_MEM, PILMSG),
  ERR(PIL_NO_FILE, PILMSG),
  ERR(PIL_ERR_FREAD, PILMSG),
  ERR(PIL_ERR_FWRITE, PILMSG),
  ERR(PIL_EOS, PILMSG),
  ERR(PIL_BAD_NAME, PILMSG),
  ERR(PIL_BAD_TYPE, PILMSG),
  ERR(PIL_BAD_MODE, PILMSG),
  ERR(PIL_BAD_LINE, PILMSG),
  ERR(PIL_NOT_IMPLEMENTED, PILMSG),
  ERR(PIL_FILE_NOT_EXIST, PILMSG),
  ERR(PIL_FILE_EXIST, PILMSG),
  ERR(PIL_FILE_NO_RD, PILMSG),
  ERR(PIL_FILE_NO_WR, PILMSG),
  ERR(PIL_LINE_BLANK, PILMSG),
  ERR(PIL_LINE_COMMENT, PILMSG),
  ERR(PIL_LINE_ERROR, PILMSG),
  ERR(PIL_NOT_FOUND, PILMSG),
  ERR(PIL_PFILES_TOO_LONG, PILMSG),
  ERR(PIL_PFILES_FORMAT, PILMSG),
  ERR(PIL_LOCK_FAILED, PILMSG),
  ERR(PIL_BOGUS_CMDLINE, PILMSG),
  ERR(PIL_NO_LOGGER, PILMSG),
  ERR(PIL_LINE_TOO_MANY, PILMSG),
  ERR(PIL_LINE_TOO_FEW, PILMSG),
  ERR(PIL_LINE_UNMATCHED_QUOTE, PILMSG),
  ERR(PIL_LINE_NO_LF, PILMSG),
  ERR(PIL_LINE_EXTRA_SPACES, PILMSG),
  ERR(PIL_BAD_VAL_BOOL, PILMSG),
  ERR(PIL_BAD_VAL_INT, PILMSG),
  ERR(PIL_BAD_VAL_REAL, PILMSG),
  ERR(PIL_BAD_VAL_INT_VAR_VECTOR, PILMSG),
  ERR(PIL_BAD_VAL_INT_VECTOR, PILMSG),
  ERR(PIL_BAD_VAL_REAL_VAR_VECTOR, PILMSG),
  ERR(PIL_BAD_VAL_REAL_VECTOR, PILMSG),
  ERR(PIL_OFF_RANGE, PILMSG),
  ERR(PIL_BAD_ENUM_VALUE, PILMSG),
  ERR(PIL_BAD_FILE_ACCESS, PILMSG),
  ERR(PIL_BAD_VALUE, PILMSG),
  ERR(PIL_VALUE_UNDEFINED, PILMSG),
  ERR(PIL_UNSPECIFIED_ERROR, PILMSG),

  ERR(PIL_ERR_MIN_IDX, PILMSG)
};

static hd_error_map sPILerrorMap = HD_ERR_MAP_INIT(sPILerrors);
