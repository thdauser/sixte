/*
 * headas_parstamp(fitsfile *fptr,    FITS file pointer
 *                 int hdunum         HDU in which to write HISTORY block (0 -> current HDU)
 *                 )
 *
 * This routine writes a block of HISTORY keywords detailing the runtime values of all
 * parameters into the specified FITS file HDU. Any parameter beginning with a '@' is
 * taken to indicate an ascii file containing a list of items. In these cases the
 * "parameter = value" string is enclosed in parentheses and the file contents are expanded
 * within the HISTORY block.
 */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include "pil.h"
#include "headas_utils.h"
#include "headas_error.h"
#include "f77_wrap.h"
#include "ape/ape_trad.h"

#define PAREN_YES 1
#define PAREN_NO 0
#define FLAG_NOTE 2


typedef struct ParameterNote ParameterNote;

struct ParameterNote
{
	char * name;
	char * note;
	ParameterNote * next;
};


static int hdpstmp(fitsfile *fptr, int linenum, char *pname, char *pval, int paren);
static int hdfstmp(fitsfile *fptr, int linenum, char *line);


static ParameterNote * headNote;

static ParameterNote *
query_parameter_note (ParameterNote * from, const char * name)
{
	ParameterNote * p = from;
	while (p) {
		if (!strcasecmp(name, p->name))
			return p;
		p = p->next;
	}
	return 0;
}



/* 23May2003: Replaces headas_parstamp. Now takes input status pointer arg */
int HDpar_stamp(fitsfile *fptr, int hdunum, int *status){
  int parcount;
  int histpar, tref;
  char ** par_name = 0;
  char ** name = 0;
  char * value = 0;
  char taskn[73];
  char datestr[32];
  char msg[144];
  char **results=NULL;
  int nitems=0, listcount;

  if (*status>0) return *status;

  get_history(&histpar);
  switch (histpar){
  case -1:
    *status = HD_ERROR_THROW("Error: HDpar_stamp(): no history parameter present", 1);
    return *status;
    break;
  case 0:
    break;
  case 1:
    if (hdunum > 0) { /* ie, if 0 then operate on the current HDU */
      fits_movabs_hdu(fptr, hdunum, NULL, status);
      if (*status) {
	sprintf(msg, "Error moving to HDU %d", hdunum);
	HD_ERROR_THROW(msg,*status);
	return *status;
      }
    }
    fits_get_system_time(datestr, &tref, status);
    if (*status){
      HD_ERROR_THROW("error finding system time",*status);
      return *status;
    }
    get_toolnamev(taskn);
    sprintf(msg,"START PARAMETER list for %s at %s",taskn, datestr);
    /* Surround START (which may wrap) with blank HISTORY lines */
    fits_write_history(fptr, " ", status);
    fits_write_history(fptr, msg, status);
    fits_write_history(fptr, " ", status); 

    /* Get names of parameters. */
    *status = ape_trad_get_par_names(&par_name);
    parcount = 0;
    for (name = par_name; 0 == *status && 0 != *name; ++name, ++parcount) {
      /* Get value of parameter as a string. */
      *status = ape_trad_get_string(*name, &value);
      if (0 != *status) {
        sprintf(msg,"Error getting parameter \"%s\"",*name);
        HD_ERROR_THROW(msg, *status);
        free(value);
        return *status;
      }

      if (*value != '@') {
        if ((*status = hdpstmp(fptr, parcount+1, *name, value, PAREN_NO)) != 0){
          HD_ERROR_THROW("parstamp error",*status);
          free(value);
          return *status;
        }
        {
          ParameterNote * rp = headNote;
          while (NULL != (rp = query_parameter_note(rp, *name))) {
            if ((*status = hdpstmp(fptr, parcount+1, *name, rp->note, FLAG_NOTE)) != 0){
              HD_ERROR_THROW("parstamp error", *status);
              free(value);
              return *status;
            }
            rp = rp->next;
          }
        }
      } else {
        results = expand_item_list(value, &nitems, ' ', 1, 1, 1, status);
        if (*status){
          /* warn on error importing file contents and proceed without expanding */
          *status=0;
          hdpstmp(fptr, parcount+1, *name, value, PAREN_NO);
          break;
        }
        hdpstmp(fptr, parcount+1, *name, value, PAREN_YES);
        for (listcount=0;listcount<nitems;listcount++){
          if (listcount==0){
            sprintf(msg,"START FILE listing: %s",hdbasename(value+1));
            /* Surround START (which may wrap) with blank HISTORY lines */
            fits_write_history(fptr, " ", status); 
            fits_write_history(fptr, msg, status);
            fits_write_history(fptr, " ", status); 
          }
          hdfstmp(fptr, listcount+1, results[listcount]);
        }
        sprintf(msg,"END FILE listing: %s",hdbasename(value+1));
        fits_write_history(fptr, msg, status);
        /* follow END (which may wrap) with a blank line */
        fits_write_history(fptr, " ", status); 
        if (results) free(results);
      }
      free(value);
    }
    sprintf(msg,"END PARAMETER list for %s",taskn);
    fits_write_history(fptr, msg, status);
    /* follow END (which may wrap) with a blank line */
    fits_write_history(fptr, " ", status); 
    break;
  }
  
  return *status;
}
FCALLSCFUN3(INT, HDpar_stamp,HDPAR_STAMP,hdpar_stamp, FITSUNIT, INT, PINT)
     
static int hdpstmp(fitsfile *fptr, int linenum, char *pname, char *pval, int flag){
  char *histstr, msg[73];
  int status=0, numlen, i;

  numlen = sprintf(msg,"%d",linenum);

  if (!(histstr = (char *) malloc(sizeof(char)*(strlen(pname)+strlen(pval)+numlen+8)))){
    sprintf(msg,"malloc error for %s\n",pname);
    status = MEMORY_ALLOCATION;
    HD_ERROR_THROW(msg, status);
    return status;
  }

  if (flag == PAREN_YES){
    sprintf(histstr,"P%d (%s = %s)",linenum,pname,pval);
  }
  else if (flag == FLAG_NOTE) {
    sprintf(histstr, "    %s", pval);
  }else{
    sprintf(histstr,"P%d %s = %s",linenum,pname,pval);
  }

  if (strlen(histstr) <= 72){
    fits_write_history(fptr, histstr, &status);
    free(histstr);
  } else {
    strncpy(msg,histstr,72);
    msg[72]='\0';
    fits_write_history(fptr, msg, &status);
    for (i=0; i < (int) strlen(histstr)/72; i++){
	  if (flag == FLAG_NOTE)
        sprintf(msg, "    ");
      else
        sprintf(msg,"P%d ",linenum);
      strncat(msg,histstr+(i+1)*72,72-numlen-2);
      msg[72]='\0';
      fits_write_history(fptr, msg, &status);
    }
    free(histstr);
  }

  if (!status)
    HD_ERROR_THROW("hdpstmp error", status);
  return status;
}

static int hdfstmp(fitsfile *fptr, int linenum, char *line){
  char *histstr, msg[73];
  int status=0, numlen, i;

  numlen = sprintf(msg,"%d",linenum);

  if (!(histstr = (char *) malloc(sizeof(char)*(strlen(line)+numlen+3)))){
    sprintf(msg,"malloc error for %s\n",line);
    status = MEMORY_ALLOCATION;
    HD_ERROR_THROW(msg, status);
    return status;
  }

  sprintf(histstr,"F%d %s",linenum,line);

  if (strlen(histstr) <= 72){
    fits_write_history(fptr, histstr, &status);
    free(histstr);
  } else {
    strncpy(msg,histstr,72);
    msg[72]='\0';
    fits_write_history(fptr, msg, &status);
    for (i=0; i < (int) strlen(histstr)/72; i++){
      sprintf(msg,"F%d ",linenum);
      strncat(msg,histstr+(i+1)*72,72-numlen-2);
      msg[72]='\0';
      fits_write_history(fptr, msg, &status);
    }
    free(histstr);
  }

  if (!status)
    HD_ERROR_THROW("hdfstmp error", status);
  return status;
}


static char * string_copy (const char * s)
{
	char * copy = malloc(strlen(s) + 1);
	if (copy)
		strcpy(copy, s);
	return copy;
}


static void store_par_note (const char * name, const char * note)
{
static ParameterNote * tailNote;
	/* this allows the same parameter to have multiple resolutions,
	 * could call query_resolved_parameter and update if it exists */
	ParameterNote * p = calloc(1, sizeof(ParameterNote));
	if (p) {
		p->name = string_copy(name);
		p->note = string_copy(note);
		if (p->name && p->note) {
			if (!headNote)
				headNote = p;
			if (tailNote)
				tailNote->next = p;
			tailNote = p;
		}
	}
}

void HDpar_note (const char * parameter, const char * format, ...)
{
	char buffer[1024];
	va_list args;
	va_start(args, format);
#ifdef WIN32
	_vsnprintf(buffer, sizeof(buffer), format, args);
#else
	vsnprintf(buffer, sizeof(buffer), format, args);
#endif
	va_end(args);

	store_par_note(parameter, buffer);
}


void HDrelease_par_notes ()
{
	ParameterNote * p = headNote;
	headNote = 0;
	while (p) {
		ParameterNote * q = p->next;
		free(p->name);
		free(p->note);
		free(p);
		p = q;
	}
}


