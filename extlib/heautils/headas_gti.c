#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "headas_error.h"
#include "headas_gti.h"

/*
 * GTILIB - library of utility routines for manipulating Good Time Intervals
 *
 * Based heavily on gtilib.f from the extractor
 *
 * Apr 2003 - Initial development
 * May 2003 - Rename to accord with HEADAS programming guidelines
 *
 *  C. Markwardt
 */



/* 
 * HDgti_init - initialize GTI structure
 *
 * Intializes an already-existing GTI structure.  It is assumed that
 * no memory needs to be deallocated.
 *
 * struct gti_struct *gti - pointer to existing GTI structure.
 *
 * RETURNS: status code
 */
int HDgti_init(struct gti_struct *gti)
{
  if (gti == 0) return NULL_INPUT_PTR;

  gti->mjdref = 0;
  gti->timezero = 0;
  gti->ngti = 0;
  gti->maxgti = 0;
  gti->start = 0; gti->stop = 0; gti->dptr = 0;
  
  return 0;
}

/* 
 * HDgti_free - deallocate memory associated with GTI structure
 *
 * Deallocates memory pointed to by GTI structure.  The structure
 * itself is not deallocated, but any time intervals are.
 * The structure is then initialized to zeroes.
 *
 * struct gti_struct *gti - pointer to existing GTI structure.
 *
 * RETURNS: status code
 */
int HDgti_free(struct gti_struct *gti)
{
  if (gti == 0) return NULL_INPUT_PTR;
  
  if (gti->dptr) free(gti->dptr);
  HDgti_init(gti);
  return 0;
}

/* 
 * HDgti_copy - deep copy GTI from one structure to another
 *
 * Perform deep copy of GTI from source to destination structures.
 * The structures must already have been allocated; the contents are
 * then transferred.
 *
 * struct gti_struct *dest - pointer to existing GTI structure, dest of copy
 * struct gti_struct *src - pointer to existing GTI structure, source of copy
 * int *status - pointer to status variable
 *
 * RETURNS: status code 
 */
int HDgti_copy(struct gti_struct *dest, struct gti_struct *src, int *status)
{
  int i;

  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if ((src == 0) || (dest == 0)) return (*status = NULL_INPUT_PTR);

  HDgti_init(dest);
  if (src->ngti > 0) {
    dest->start = (double *)malloc(sizeof(double)*2*src->ngti);
    if (dest->start == 0) return (*status = MEMORY_ALLOCATION);
  }

  dest->mjdref = src->mjdref;
  dest->timezero = src->timezero;
  dest->ngti = src->ngti;
  dest->stop = dest->start + dest->ngti;
  dest->dptr = (void *) dest->start;
  dest->maxgti = src->ngti;

  for (i=0; i<dest->ngti; i++) {
    dest->start[i] = src->start[i];
    dest->stop[i] = src->stop[i];
  }

  return 0;
}


/*
 * HDget_frac_time - get (potentially) fractional time keyword from FITS header
 *
 * Some OGIP FITS timing values are allowed to be specified either as
 * a single keyword, or a pair of keywords which gives integer and
 * fractional values which are to be added.  This allows maximum time
 * precision.  For example, for the observation starting time, it can
 * either be found by:
 *
 *                                  RETURNS    *vali     *valf
 *    TSTART  - single value        TSTART     TSTART    0
 *
 *    TSTARTI - integer value       TSTARTI+F  TSTARTI   TSTARTF
 *    TSTARTF - fractional value
 *
 * This routine read a FITS timing keyword(s).  The returned values
 * are given above.  If the user passes null values for vali and/or
 * valf, then HDget_frac_time will gracefully skip storing those values.
 *
 * fitsfile *fileptr - open FITS file pointer, open to HDU where
 *                     keywords are to be found
 * char *key - keyword name, "TSTART", "TIMEZERO", etc.
 * double *vali - upon return, and if vali is not null, then the whole
 *                number of seconds is stored in *vali.
 * double *valf - upon return, and if valf is not null, then the fractional
 *                number of seconds is stored in *valf.
 * int *status - pointer to status variable
 *
 * RETURNS: (TSTART) or (TSTARTI+TSTARTF)
 *
 */
double HDget_frac_time(fitsfile *fileptr, char *key, 
		       double *vali, double *valf,
		       int *status)
{
  double value, valuei, valuef;
  char keynamei[10];
  char keynamef[10];

  if (status == 0) return 0;
  if (*status) return (*status);
  if ((fileptr == 0) || (key == 0)) return (*status = NULL_INPUT_PTR);

  /* Construct the "I"nteger and "F"ractional keyword names */
  strncpy(keynamei, key, 7);
  keynamei[7] = 0;
  strcat(keynamei, "I");
  strncpy(keynamef, key, 7);
  keynamef[7] = 0;
  strcat(keynamef, "F");

  /* Read pair of values first ... */
  fits_write_errmark();
  fits_read_key(fileptr, TDOUBLE, keynamei, &valuei, NULL, status);
  fits_read_key(fileptr, TDOUBLE, keynamef, &valuef, NULL, status);
  /* It is never an error for either of these two keywords to be missing. */
  fits_clear_errmark();

  /* ... and if not found, then search for single value */
  if (*status) {
    *status = 0;
    fits_read_key(fileptr, TDOUBLE, key, &value, NULL, status);
    if (*status) return 0;

    if (vali) *vali = value;
    if (valf) *valf = 0;
    return value;
  }

  /* Assign fractional values ... */
  if (vali) *vali = valuei;
  if (valf) *valf = valuef;
  return valuei + valuef;  /* ... and return total */

}


/* 
 * HDgti_shell_sort - internal sort routine
 *
 * Sort an array.  The result is an array of indices tbl_index[],
 * which place the source array in ascending or descending order,
 * depending on the value of "ascend."
 *
 * The user is responsible to allocate and deallocate the arrays.
 *
 * int * tbl_index - an n-element array.  The sort indices are
 *                   returned in this array.
 * int n - number of elements in col and tbl_index arrays.
 * double *col - an n-element array, the data to be sorted.
 * int ascend - sort by ascending order (1=yes, 0=no)
 *
 * RETURNS: nothing
 */
void HDgti_shell_sort(int * tbl_index, int n, double *col, 
                int ascend)
{
     double aln2i = 1.442695022;
     double tiny = 1.0e-5;
     int nn,m,j,lognb2;
     int iib,*idx;
     
     lognb2 = (int)(log((double)n) * aln2i + tiny);

     m = n;
     for (nn = 0; nn < lognb2; nn++)
     {
         m >>= 1;
         for (j = m; j < n; j++)
         {
             idx = tbl_index + (j - m);
             iib = tbl_index[j];
             if( !ascend ) 
             {
             	while ((idx >= tbl_index) && (col[*idx] < col[iib]))
             	{
                	  idx[m] = *idx;
                  	idx -= m;
             	}
             }
             else 
             {
             	while ((idx >= tbl_index) && (col[*idx] > col[iib]) )
             	{
                	  idx[m] = *idx;
                  	idx -= m;
             	}
             }
                 
            idx[m] = iib;
          }
      }
}


/* 
 * HDgti_grow - enlarge the storage of an existing GTI structure
 *
 * Takes an existing GTI structure, and enlarges it to a new size,
 * i.e., increases the maximum number of GTI rows that can be stored
 * in the structure.  
 * 
 * The storage can never shrink.  Thus, users can call this function
 * whenever they anticipate needing a certain number of GTI rows.  The
 * existing values will be preserved.
 *
 * HDgti_grow() handles the case where the GTI structure has been
 * initialized and has no allocated storage associated with it, yet.
 *
 * struct gti_struct *gti - pointer to existing GTI structure, to be
 *                           enlarged.
 * int new - new maximum number of rows of storage which gti should contain.
 * int *status - pointer to status variable
 *
 * RETURNS: status code 
 */
int HDgti_grow(struct gti_struct *gti, int new, int *status)
{
  int i;
  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if (gti == 0) return (*status = NULL_INPUT_PTR);
  
  if (new <= gti->maxgti) return 0;

  if (gti->maxgti == 0) {
    gti->start = (double *)malloc(sizeof(double)*new*2);
    gti->stop  = gti->start + new;
    gti->dptr  = (void *) gti->start;
  } else {
    double *newp;
    newp = (double *)realloc(gti->start, sizeof(double)*new*2);
    if (newp == 0) return (*status = MEMORY_ALLOCATION);

    for (i=gti->ngti-1; i>=0; i--) {
      newp[new+i] = newp[gti->ngti+i];
    }
    gti->start = newp;
    gti->stop  = newp+new;
    gti->dptr  = (void *) newp;
  }

  gti->maxgti= new;
  return 0;
}

/* 
 * HDgti_read - read a GTI extension from a FITS file
 *
 * Reads a GTI extension from a FITS file.  If filename is non-NULL,
 * then HDgti_read will open the specified file name (and optional
 * extension name).  If filename is NULL, then it is assumed that the
 * file is already open, and a CFITSIO file handle is available in
 * (*fptr).  
 * 
 * Upon return, if fptr is non-NULL, the file will remain open and
 * *fptr will be a CFITSIO file handle.  Otherwise, the file will be
 * closed.
 * 
 * HDgti_read() will search for the first extension matching the
 * extname parameter.  Set extename="*" to match the first extension
 * of the open file.
 * 
 * The user can choose the column names to read.  The user can also
 * provide a reference GTI.  If the two GTIs have different timezero
 * values, then the times of the new gti are adjusted to match the
 * timezero value of the refer_to GTI.
 * 
 *
 * char *filename - name of FITS file to read (including CFITSIO syntax),
 *                  or NULL
 * struct gti_struct *gti - pointer to existing GTI structure, to be
 *                          filled with data from FITS file.
 * char *extname - name of extension,    or 0 for default of "GTI"
 *                 (set extname="*" to use first open table extension)
 * char *start   - name of START column, or 0 for default of "*START*"
 * char *stop    - name of STOP column,  or 0 for default of "*STOP*"
 * struct gti_struct *refer_to - pointer to existing GTI structure.
 * fitsfile **fptr - three cases:
 *    fptr == 0: 
 *        HDgti_read will open filename, and close it before returning
 *    fptr != 0 and filename != NULL: 
 *        HDgti_read will open filename, file will remain open, and
 *        (*fptr) contains the CFITSIO handle of open file
 *    fptr != 0 and filename == NULL: 
 *        HDgti_read assumes file is already open, and
 *        (*fptr) contains the CFITSIO handle of open file
 * int *status - pointer to status variable
 *
 * RETURNS: status code */
int HDgti_read(char *filename, struct gti_struct *gti, 
	    char *extname, char *start, char *stop,
	    struct gti_struct *refer_to, 
	    fitsfile **fptr, int *status)
{
  int nhdu = 0;
  int curhdu = 0;
  fitsfile *gtifile = 0;
  char fextname[FLEN_CARD];
  int startcol, stopcol;
  long int nrows = 0;

  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if (gti == 0) return (*status = NULL_INPUT_PTR);
  if ((filename == 0) && ((fptr == 0) || (*fptr == 0))) 
    return (*status = NULL_INPUT_PTR);

  /* Initialize GTI */
  HDgti_init(gti);

  /* Initialize input parameters */
  if ((extname == 0) || (*extname == 0)) extname = "GTI";
  if ((start == 0) || (*start == 0)) start = "*START*";
  if ((stop == 0) || (*stop == 0)) stop = "*STOP*";

  if (filename) {
    /* Open file and get HDU information */
    if (fits_open_file(&gtifile, filename, READONLY, status)) {
      return *status;
    }
  } else {
    /* File is already open, use it */
    gtifile = *fptr;
  }

  if (fits_get_num_hdus(gtifile, &nhdu, status)) goto ERRCLEAN;
  fits_get_hdu_num(gtifile, &curhdu);
  if (curhdu == 1) curhdu = 2;  /* Skip primary array */

  /* Scan through the HDUs starting at the current one, until one with
     EXTNAME = '*GTI*' */
  while (curhdu <= nhdu) {
    int hdutype = BINARY_TBL;
    if (fits_movabs_hdu(gtifile, curhdu, &hdutype, status)) {
      goto ERRCLEAN;
    }

    /* Skip image HDUs */
    if (hdutype == BINARY_TBL) {

      /* Wild card matches current HDU */
      if (extname[0] == '*' && extname[1] == 0) break;

      /* Otherwise match by name */
      if (fits_read_key(gtifile, TSTRING, "EXTNAME", fextname, NULL, status)) {
	goto ERRCLEAN;
      }
      if (strstr(fextname, extname)) break;
    }
      
    curhdu ++;
  }
  if (curhdu > nhdu) {
    *status = END_OF_FILE;
    goto ERRCLEAN;
  }

  /* Read keywords and column information from this GTI extension */
  if (fits_get_colnum(gtifile, 0, start, &startcol, status)) goto ERRCLEAN;
  if (fits_get_colnum(gtifile, 0, stop,  &stopcol,  status)) goto ERRCLEAN;
  if (fits_get_num_rows(gtifile, &nrows, status)) goto ERRCLEAN;
  gti->mjdref = HDget_frac_time(gtifile, "MJDREF", 0, 0, status);
  if (*status) {
    *status = 0;
    gti->mjdref = 0;
  }

  /* Read TIMEZERO-like keywords */
  gti->timezero = HDget_frac_time(gtifile, "TIMEZERO", 0, 0, status);
  if (*status) {
    *status = 0;
    gti->timezero = 0;
  }


  /* Start filling in output GTI structure */
  gti->start = (double *)malloc(2*sizeof(double)*nrows);
  if (gti->start == 0) goto ERRCLEAN;
  gti->stop = gti->start + nrows;
  gti->dptr = (void *) gti->start;

  /* Read data from extension */
  fits_read_col(gtifile, TDOUBLE, startcol, 1, 1, nrows, 0,
		gti->start, 0, status);
  fits_read_col(gtifile, TDOUBLE, stopcol, 1, 1, nrows, 0,
		gti->stop, 0, status);
  if (*status) {
    free(gti->dptr);
    gti->start = gti->stop = 0;
    gti->dptr = 0;
    goto ERRCLEAN;
  }

  /* If the refer_to structure was passed, then refer this GTI to the
     reference one.  Basically, apply an offset if the TIMEZERO values
     are different. */
  if (refer_to && (refer_to->timezero != gti->timezero)) {
    int i;
    double dt = gti->timezero - refer_to->timezero;

    for (i=0; i<(2*nrows); i++) {
      gti->start[i] += dt;
    }
    gti->timezero = refer_to->timezero;
  }
  gti->ngti = nrows;
  gti->maxgti = nrows;

 ERRCLEAN:
  /* Either close or return the file pointer, depending on the value
     of fptr */
  if (gtifile) {
    if (fptr) {
      *fptr = gtifile;
    } else {
      int mystatus = 0;
      fits_close_file(gtifile, &mystatus);
      gtifile = 0;
    }
  }
  return *status;
  
}

/* 
 * HDgti_write - create a GTI extension and write it
 *
 * Writes a GTI extension to a FITS file.  A binary table is created
 * and populated with the gti "gti".  Basic OGIP GTI keywords are
 * written.
 *
 * The file is left open upon return so that users may manipulate the
 * header afterwards.
 *
 * The units are assumed to be seconds.  The user can choose the
 * extname, and column names to write.
 * 
 *
 * fitsfile *fptr - already-open FITS file, open for writing
 * struct gti_struct *gti - populated GTI structure, to be written to fptr
 * char *extname - name of extension, or 0 for default of "STDGTI"
 * char *start - name of START column, or 0 for default of "START"
 * char *stop - name of STOP column, or 0 for default of "STOP"
 * int *status - pointer to status variable
 *
 * RETURNS: status code 
 */
int HDgti_write(fitsfile *fptr, struct gti_struct *gti, 
	     char *extname, char *start, char *stop,
	     int *status)
{
  char *colnames[2];
  char *colforms[2] = {"D", "D"};
  char *colunits[2] = {"s", "s"};
  double tstart = 0, tstop = 0;
  int i;

  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if (gti == 0) return (*status = NULL_INPUT_PTR);

  if ((extname == 0) || (*extname == 0)) extname = "STDGTI";
  if ((start == 0) || (*start == 0)) start = "START";
  if ((stop == 0)  || (*stop == 0))  stop  = "STOP";

  colnames[0] = start;
  colnames[1] = stop;
  
  fits_create_tbl(fptr, BINARY_TBL, gti->ngti, 2, 
		  colnames, colforms, colunits, extname, status);
  /* Apply descriptive column comments */
  fits_modify_comment(fptr, "TTYPE1", "GTI start time", status);
  fits_modify_comment(fptr, "TTYPE2", "GTI stop  time", status);

  if (*status) {
    HD_ERROR_THROW("ERROR: Could not create output GTI", *status);
    return *status;
  }

  fits_update_key(fptr, TSTRING, "HDUCLASS", "OGIP", 
		  "Conforms to OGIP/GSFC standards", status);
  fits_update_key(fptr, TSTRING, "HDUCLAS1", "GTI",
		  "Contains good time intervals", status);
  fits_update_key(fptr, TSTRING, "HDUCLAS2", "STANDARD", 
		  "Contains standard good time intervals", status);
  fits_update_key(fptr, TSTRING, "HDUVERS", "1.0.0",
		  "Version of GTI header", status);

  fits_update_key(fptr, TDOUBLE, "TIMEZERO", &(gti->timezero),
		  "Zero-point offset for TIME column", status);
  /* Use special CFITSIO routine to ensure the maximum number of
     digits of precision. */
  fits_write_key_dbl(fptr, "MJDREF", gti->mjdref, 15, 
		     "MJD Epoch of TIME = 0", status);

  if (gti->ngti > 0) {
    tstart = gti->start[0]; 
    tstop  = gti->stop[0];
    for (i=1; i<gti->ngti; i++) {
      if (gti->start[i] < tstart) tstart = gti->start[i];
      if (gti->stop[i]  > tstop)  tstop  = gti->stop[i];
    }
  }

  /* Save TSTART/TSTOP keywords */
  fits_update_key(fptr, TDOUBLE, "TSTART", &tstart,
		  "Start time of GTI", status);
  fits_update_key(fptr, TDOUBLE, "TSTOP", &tstop,
		  "Stop  time of GTI", status);
  
  /* Write actual data */
  if (gti->ngti > 0) {
    fits_write_col(fptr, TDOUBLE, 1, 1, 1, gti->ngti,
		   gti->start, status);
    fits_write_col(fptr, TDOUBLE, 2, 1, 1, gti->ngti,
		   gti->stop, status);
  }

  return (*status);
}


/* 
 * Internal GTI merging routine
 *
 *   (cribbed from GTIMRG1 fortran routine from extractor/gtilib.f)
 *
 */
void HDgtimrg1(struct gti_struct *gti, double *times, 
	     int *types, int *idx, int ilist, int mode, int *status)
{
  int mode_or = (mode == GTI_OR);
  int mode_and = (mode == GTI_AND);
  int a, b, i;

  HDgti_shell_sort(idx, ilist, times, 1);

  /* 
   * Now produce a merged gti in gtistart/gtistop.
   * This section is a finite-state machine that takes the combined, sorted
   * time list (times) and the type codes and creates an ANDed or ORed GTI.
   * If the two original lists are A and B, A+ means list A is on, A- means
   * it's off, etc, and the numbers are the type codes given above, then the
   * FSM looks like:
   *
   *  (A+B+)  3->  <-2 (A+B-)
   * 
   *  1| ^             1| ^
   *   v |0             v |0
   *
   *  (A-B+)  3->  <-2 (A-B-) <-start
   * 
   * For ANDing, a new GTI is begun on the transitions to state A+B+, and
   * a currently open one is ended on the transtitions away from it.
   * For ORing, a new GTI is begun on transitions away from state A-B-,
   * and an open one ended on transtiions to it.
   * There's a complication: If one of the input lists has overlapping
   * GTIs, we might do the wrong thing, such as starting an extra GTI or
   * ending one too soon.  The counters A and B keep track of how many
   * levels of GTI we've seen, so we only take action at the right times.
   * 
   * This code is derived from the FSM in /ftools/xselect/extractor/
   * extractor.f, subroutine fingti(), which was probably written by
   * Bruce O'Neel.  He presorted the arrays instead of using level
   * counters, however.
   */

  a = 0;
  b = 0;
  gti->ngti = 0;

  for (i=0; i<ilist; i++) {

    /* 
     * A turns on.  Start a new GTI if A is not already on, and:
     * AND case: B is on too;  OR case: B is off (else one is already
     * active).
     */

    if (types[idx[i]] == 0) {
      a++;
      if ( (mode_and && (a == 1) && (b >= 1) ) ||
	   (mode_or  && (a == 1) && (b <= 0)) ) {
	
	if (gti->ngti == gti->maxgti) {
	  int new = gti->maxgti;
	  if (gti->maxgti == 0) new = 16;
	  if (new > 256)   new = 256;

	  HDgti_grow(gti, gti->ngti+new, status);
	  if (*status) return;
	}
	gti->start[gti->ngti] = times[idx[i]];
      }

      /* 
       * A turns off.  End the current GTI if: AND: B is on; OR: B is off too.
       */ 

    }  else if (types[idx[i]] == 1) {
      a --;
      if ( (mode_and && (a == 0) && (b >= 1) ) ||
	   (mode_or  && (a == 0) && (b <= 0)) ) {
	gti->stop[gti->ngti] = times[idx[i]];
	gti->ngti ++;
      }

      /* 
       * B turns on. Start a new GTI if B isn't already on, and:
       * AND: A is on too;  OR: B is off.
       */

    } else if (types[idx[i]] == 2) {
      b ++;
      if ( (mode_and && (b == 1) && (a >= 1) ) ||
	   (mode_or  && (b == 1) && (a <= 0)) ) {
	
	if (gti->ngti == gti->maxgti) {
	  int new = gti->maxgti;
	  if (gti->maxgti == 0) new = 16;
	  if (new > 256)   new = 256;

	  HDgti_grow(gti, gti->ngti+new, status);
	  if (*status) return;
	}
	gti->start[gti->ngti] = times[idx[i]];
      }

      /* 
       * B turns off.  End the current GTI if: AND: A is on; OR: A is off too.
       */

    } else if (types[idx[i]] == 3) {
      b --;
      if ( (mode_and && (b == 0) && (a >= 1) ) ||
	   (mode_or  && (b == 0) && (a <= 0)) ) {
	gti->stop[gti->ngti] = times[idx[i]];
	gti->ngti ++;
      }
    }
  }

  return;
}

/* 
 * HDgti_merge - merge two GTIs either using intersection or union
 *
 * Merges two different good time interval lists.  If the mode is
 * GTI_AND, then the intersection between the two lists is determined.
 * If the mode is GTI_OR, then the union between the two lists is
 * found.
 * 
 * This routine is based heavily on the fortran version, taken from
 * the HEASARC extractor (gtilib.f).
 *
 * int mode - merging mode, either GTI_AND or GTI_OR
 * struct gti_struct *gti - pointer to existing GTI structure, result of merge
 * struct gti_struct *agti - pointer to existing GTI structure, 1st input GTI
 * struct gti_struct *bgti - pointer to existing GTI structure, 2nd input GTI
 * int *status - pointer to status variable
 *
 * RETURNS: status code 
 */
int HDgti_merge(int mode, struct gti_struct *gti, 
	     struct gti_struct *agti, struct gti_struct *bgti, 
	     int *status)
{
  double *times = 0;
  int *types = 0, *idx = 0;
  int anum, bnum;
  int ilist, i, j;

  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if ((gti == 0) || (agti == 0) || (bgti == 0)) 
    return (*status = NULL_INPUT_PTR);
  if ((mode != GTI_AND) && (mode != GTI_OR)) return (*status = -1);

  anum = agti->ngti;
  bnum = bgti->ngti;

  HDgti_init(gti);

  /* Special cases - hmm, this defeats the whole purpose of HDgti_clean() */
#if 0
  if ((anum == 0) && (bnum == 0)) return (*status);
  if (mode == GTI_AND) {
    if ((anum == 0) || (bnum == 0)) return (*status);
  } else {
    struct gti_struct *src = 0;
    if (anum == 0) src = bgti;
    if (bnum == 0) src = agti;

    if (src) {
      HDgti_copy(gti, src, status);
      return (*status);
    }
  }
#endif

  /* 
   * Build a combined list
   * times(i) holds the time
   * types(i) holds the type:  0 for GTI list A start, 1 for list A stop,
   *   2 for list B start, and 3 for list B stop.
   * In the process, also filter out invalid GTIs of negative or zero length.
   */

  times = (double *) malloc(sizeof(double)*2*(anum+bnum));
  types = (int *) malloc(sizeof(int)*2*(anum+bnum));
  idx   = (int *) malloc(sizeof(int)*2*(anum+bnum));

  if ((times == 0) || (types == 0) || (idx == 0)) {
    if (times) free(times); times = 0;
    if (types) free(types); types = 0;
    if (idx)   free(idx);   idx = 0;
    return (*status = MEMORY_ALLOCATION);
  }

  /* Fill the array with the first GTI's times ... */
  ilist = 0;
  for (i=0; i<anum; i++) {
    if (agti->start[i] < agti->stop[i]) {
      times[ilist] = agti->start[i];
      types[ilist] = 0;
      ilist++;
      times[ilist] = agti->stop[i];
      types[ilist] = 1;
      ilist++;
    }
  }

  /* ... and the second GTI's times */
  for (i=0; i<bnum; i++) {
    if (bgti->start[i] < bgti->stop[i]) {
      times[ilist] = bgti->start[i];
      types[ilist] = 2;
      ilist++;
      times[ilist] = bgti->stop[i];
      types[ilist] = 3;
      ilist++;
    }
  }

  /* Case of empty list */
  if (ilist == 0) {
    free(times);
    free(types);
    free(idx);
    return (*status);
  }

  /* Pre-populate sort index array */
  for (i=0; i<ilist; i++) idx[i] = i;
  /* Perform actual merge operation */
  HDgtimrg1(gti, times, types, idx, ilist, mode, status);

  /* 
   * The FSM can leave zero-length GTIs, so clean them up...
   */

  i = 0;
  while (i<gti->ngti) {
    if (gti->start[i] >= gti->stop[i]) {
      for (j=i+1; j<gti->ngti; j++) {
	gti->start[j-1] = gti->start[j];
	gti->stop[j-1]  = gti->stop[j];
      }
      gti->ngti --;
    } else {
      i++;
    }
  }

  /* 
   * ...and combine any contiguous ones.
   */

  i = 1;
  while (i<gti->ngti) {
    if (gti->stop[i-1] == gti->start[i]) {
      gti->stop[i-1] = gti->stop[i];
      for (j=i; j<(gti->ngti-1); j++) {
	gti->start[j] = gti->start[j+1];
	gti->stop[j]  = gti->stop[j+1];
      }
      gti->ngti --;
    } else {
      i++;
    }
  }

  free(times);
  free(types);
  free(idx);

  return (*status);
}

/* 
 * HDgti_clean - clean a GTI by sorting, removing duplicates, overlaps
 *
 * HDgti_clean:  Sort a gti list and remove invalid and overlapping
 *   intervals by calling HDgti_merge to OR it with nothing.
 * 
 * This routine is based heavily on the fortran version, taken from
 * the HEASARC extractor.
 *
 * struct gti_struct *gti - pointer to existing GTI structure, result of clean
 * struct gti_struct *ogti - pointer to existing GTI structure, to be cleaned
 * int *status - pointer to status variable
 *
 * RETURNS: status code 
 */
int HDgti_clean(struct gti_struct *gti, struct gti_struct *ogti,
	    int *status)
{
  struct gti_struct null_gti;

  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if ((gti == 0) || (ogti == 0)) return (*status = NULL_INPUT_PTR);

  HDgti_init(&null_gti);
  return HDgti_merge(GTI_OR, gti, ogti, &null_gti, status);
}

/*
 * HDgti_exp - compute overlap exposure of a time bin with GTI
 * 
 * Exposure: What part of the time interval (t1, t2) is contained
 * within the Good Time Intervals?  (For the fractional exposure, just
 * divide gtiexp by (t2-t1).)
 *
 *  NOTE: It is *VITAL* that the GTI list be ordered and
 *     contain no overlaps!  If in doubt, call HDgti_clean first.
 *     also, we assume that t1, t2, and the gti list are all
 *     in the same units.
 *
 * One difference from the FORTRAN version is that if there are no
 * good time intervals (i.e. gti->ngti == 0), then zero is returned
 * here.  In the FORTRAN version, the whole exposure is returned.
 *
 * double t1 - start of time bin.
 * double t2 - stop of time bin.  (note: t1 < t2)
 * struct gti_struct *gti - pointer to existing GTI structure, whose
 *        overlap with the time bin is to be computed.
 * int *status - pointer to status variable
 * 
 * RETURNS: overlap exposure */
double HDgti_exp(double t1, double t2, struct gti_struct *gti, int *status)
{
  double total = 0;
  int i;

  if (status == 0) return 0;
  if (*status) return 0;
  if (gti == 0) {
    *status = NULL_INPUT_PTR;
    return 0;
  }

  /* Invalid input time bin */
  if (t1 > t2) {
    *status = -1;
    return 0;
  }
  if (t1 == t2) return 0;

  /* Return zero if there are no good time intervals.  NOTE: this is
     different from Guerber's FORTRAN GTICLEAN.  The policy decision
     to change the exposure back to non-zero should be made at a
     higher level. */
  if (gti->ngti == 0) return 0;
  
  for (i=0; i<gti->ngti; i++) {
    /* There are 6 cases: */
    /*     1) g1 g2 t1 t2: cycle to next interval */
    if (gti->stop[i] <= t1) {
      continue;
    }

    /*     2) t1 t2 g1 g2: gti list is assumed ordered, so we're done! */
    else if (t2 <= gti->start[i]) {
      break;
    }

    /*     3) g1 t1 g2 t2  (Note: Only 2 comparisons to reach any branch!)*/
    else if (gti->start[i] <= t1) {
      if (gti->stop[i] <= t2) {
	total += (gti->stop[i] - t1);
      } else {
    /*         4) g1 t1 t2 g2 */
	total += (t2 - t1);
      }
    /*         5) t1 g1 g2 t2 */
    } else {
      if (gti->stop[i] <= t2) {
	total += (gti->stop[i] - gti->start[i]);
	
    /*         6) t1 g1 t2 g2 */
      } else {
	total += (t2 - gti->start[i]);
      }
    }

  } /* end for(i) */

  return total;
}

/*
 * HDgti_where - which good time intervals a set of times falls into
 *
 * HDgti_where examines an array of times, and determines which good time
 * intervals the times fall into.  This routine is most highly
 * optimized for times which are sorted in time order.
 *
 * An interval of -1 indicates a time that does not fall into a good
 * time interval.
 *
 * struct gti_struct *gti - pointer to existing GTI structure.
 * int ntimes - size of times array.
 * double *times - an ntimes-element array, times to be examined.
 * int *segs - an ntimes-element array.  Upon return, the values give 
 *   which good time interval the time falls into:
 *    (times[i] >= gti->start[segs[i]]) && (times[i] < gti->stop[segs[i]])
 *   or, if segs[i] is -1, the time does not fall into a GTI
 * int *status - pointer to status variable
 * 
 * RETURNS: status code */
int HDgti_where(struct gti_struct *gti, int ntimes, 
	     double *times, int *segs, int *status)
{
  int i, j;
  int ngti;
  double t, tmin, tmax;
  int s, s1;

  if (status == 0) return NULL_INPUT_PTR;
  if (*status) return (*status);
  if ((gti == 0) || (times == 0) || (segs == 0)) 
    return (*status = NULL_INPUT_PTR);

  if (ntimes == 0) return 0;
  ngti = gti->ngti;
  for (i=0; i<ntimes; i++) segs[i] = -1;
  if (ngti == 0) return 0;

  tmin = gti->start[0];
  tmax = gti->stop[ngti-1];

  i = 0;
  s = -1;  /* Prevents compiler warning */

  while (i < ntimes) {
    t = times[i];

    /* If segment number is not known, then do a linear search */
    if (t >= tmax) {
      /* Check out-of-bounds interval first */
      s = ngti+ngti;
    } else {
      /* Next, check each good time interval, and the bad time
	 interval that precedes it */
      for (j=0; j<ngti; j++) {
	if (t < gti->start[j]) {       /* Bad time interval */
	  s = 2*j;
	  break;
	} else if (t < gti->stop[j]) { /* Good time interval */
	  s = 2*j+1;
	  break;
	}
      }
    }
    
    /* If this is a good time interval, record it */
    if (s % 2) segs[i] = s/2;
    i++;
    if (i >= ntimes) break;

    /* s1 is the good time interval number, or the bad time which
       precedes it */
    s1 = s/2;

    if (s % 2) {  /* GOOD TIME */
      /* Previous time was within segment s1 */
      while ((i < ntimes) && 
	     (times[i] >= gti->start[s1]) && (times[i] < gti->stop[s1])) {
	/* This time is within segment s1 as well */
	segs[i++] = s1;
      }
    } else {      /* BAD TIME */
      
      if (times[i] < tmin) {
	/* Skip over times before start of good time intervals */
	while ((i < ntimes) && (times[i] < tmin)) i++;
      }
      if (i >= ntimes) break;
      if (times[i] >= tmax) {
	/* Skip over times beyond end of good time intervals */
	while ((i < ntimes) && (times[i] >= tmax)) i++;
      }
      if (i >= ntimes) break;
      if (s1 >= 1) {
	/* Advance over events *between* the previous GTI and the next */
	while ((i < ntimes) && 
	       (times[i] >= gti->stop[s1-1]) && (times[i] < gti->start[s1])) {
	  i++;
	}
      }
    }


    /* End of while loop.  i has already advanced to the next
       unrecorded time, or to the end of the array of times, so there
       is no need for i++ here. */
  }
  
  return 0;
}

