/******************************************************************************
 *   File name: headas_copykeys.c                                             *
 *                                                                            *
 * Description: Copies all general keywords from one HDU to another           *
 *              Does not copy any specific keywords that describe the         *
 *              structure.                                                    *
 *                                                                            *
 *      Author: kaa for HEASARC/GSFC/NASA                                     *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "fitsio.h"
#include "headas_copykeys.h"

#define EXCSIZE 56

#ifdef __cplusplus
extern "C" {
#endif

int HDcopy_keys(fitsfile *inptr, fitsfile *outptr, int docomments, int dohistory, int *status)
{

  int nexc;

  char *inclist[1] = { "*" };
  char *exclist[EXCSIZE] = {
    "SIMPLE", "BITPIX", "NAXIS", "NAXIS#", "EXTEND",
    "XTENSION", "EXTNAME", "PCOUNT", "GCOUNT", "TFIELDS",
    "TTYPE#", "TBCOL#", "TFORM#", "TSCAL#", "TZERO#",
    "TNULL#", "TUNIT#", "THEAP", "TDIM#", "TDISP#",
    "GROUPS", "BSCALE", "BZERO", "BUNIT", "BLANK",
    "CTYPE#", "CRPIX#", "CROTA#", "CRVAL#", "CDELT#",
    "TLMIN#", "TLMAX#", "TDMIN#", "TDMAX#", "OPTIC#",
    "TCRPX#", "TCRVL#", "TCDLT#", "TCTYP#", "TCUNI#",
    "TCD#", "TCROT#", "PLTSCL#", "CD_*", "PC_*",
    "MTYPE#", "MFORM#", "#METYP#", "#DSTYP#", "#DSVAL#",
    "#DSREF#", "CHECKSUM", "DATASUM", "END", "        ", "        "};

  char card[FLEN_CARD];

  nexc = EXCSIZE - 2;
  if ( !docomments ) {
    nexc++;
    strcpy(exclist[nexc-1],"COMMENT");
  }
  if ( !dohistory ) {
    nexc++;
    strcpy(exclist[nexc-1],"HISTORY");
  }

  /* reset to the start of the input header */

  if ( (*status = fits_read_record(inptr, 0, card, status)) ) {
    return *status;
  }

  /* loop round the records till we hit the end */

  while ( *status == 0 ) {

    fits_find_nextkey(inptr, inclist, 1, exclist, nexc, card, status);

    fits_write_record(outptr, card, status);

  }


  *status = 0;

  return *status;

}
