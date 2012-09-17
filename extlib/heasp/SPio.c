
#include <string.h>
#include <stdio.h>
#include "heaspio.h"

/******************************** SP_read_key **************************************/

/* reads a keyword */

int SP_read_key(fitsfile *fptr, int datatype, char *keyword_name, void *keyword_value, void *keyword_default)
{

  return(SPII_read_key(fptr, datatype, keyword_name, 0, keyword_value, keyword_default));

}

/******************************** SPII_read_key **************************************/

/* reads a keyword from a type II spectral file - first checks for the keyword then, if
   ispec > 0, for a value in a column */

int SPII_read_key(fitsfile *fptr, int datatype, char *keyword_name, long ispec, void *keyword_value, void *keyword_default)
{
  int status=0;
  int anynul=0;
  int colnum;
  char comment[FLEN_COMMENT];

  fits_read_key(fptr, datatype, keyword_name, keyword_value, comment, &status);
  if (status && ispec > 0) {
    status = 0;
    fits_clear_errmsg();
    fits_get_colnum(fptr, CASEINSEN, keyword_name, &colnum, &status);
    fits_read_col(fptr, datatype, colnum, ispec, 1, 1, NULL, keyword_value, &anynul, &status);
  }

  if (!status) {
    if (datatype == TFLOAT)  headas_chat(5, " %s = %f\n", keyword_name, *(float *)keyword_value);
    if (datatype == TSTRING) headas_chat(5, " %s = %s\n", keyword_name , (char *)keyword_value);
    if (datatype == TINT)    headas_chat(5, " %s = %d\n", keyword_name, *(int *)keyword_value);
    if (datatype == TLONG)   headas_chat(5, " %s = %ld\n", keyword_name, *(long *)keyword_value);
  } else if (status) {

    if (datatype == TFLOAT)  *(float *)keyword_value = *(float *)keyword_default;
    if (datatype == TSTRING) strcpy((char *)keyword_value, (char *)keyword_default);
    if (datatype == TINT)    *(int *)keyword_value = *(int *)keyword_default;
    if (datatype == TLONG)   *(long *)keyword_value = *(int *)keyword_default;

    if (datatype == TFLOAT)  headas_chat(1, "***SPII_read_key: Cannot find %s keyword - setting it to %f\n", keyword_name, *(float *)keyword_default);
    if (datatype == TSTRING) headas_chat(1, "***SPII_read_key: Cannot find %s keyword - setting it to %s\n", keyword_name, (char *)keyword_default);
    if (datatype == TINT)    headas_chat(1, "***SPII_read_key: Cannot find %s keyword - setting it to %d\n", keyword_name, *(int *)keyword_default);
    if (datatype == TLONG)   headas_chat(1, "***SPII_read_key: Cannot find %s keyword - setting it to %ld\n", keyword_name, *(long *)keyword_default);

  }

  if (datatype != TFLOAT && datatype != TSTRING && datatype != TINT && datatype != TLONG)
    headas_chat(1, "%d is an unsupported datatype in SPII_read_key\n", datatype);

  return(status);
}

/******************************** SP_read_col **************************************/

int SP_read_col(fitsfile *fptr, int datatype, char *column_name, long number_values, void *column_values)
{

  return(SPII_read_col(fptr, datatype, column_name, 1, number_values, column_values));
}

/******************************** SPII_read_col **************************************/

int SPII_read_col(fitsfile *fptr, int datatype, char *column_name, long ispec, long number_values, void *column_values)
{
  int status=0;
  int i, colnum, anynul, typecode;
  long repeat, width;
  int *ivalues;
  long *lvalues;
  float *fvalues;

  /* try to read a keyword of name column_name
     if we succeeded then set the array to the returned value */

  if (!SP_read_key(fptr, datatype, column_name, column_values, column_values)) {
    if (datatype == TINT) {
       headas_chat(5, "%s = %d read from keyword\n", column_name, *(int *)column_values);
      ivalues = (int *)column_values;
      for(i=0; i<number_values; i++) ivalues[i] = *(int *)column_values;
    } else if (datatype == TLONG) {
       headas_chat(5, "%s = %ld read from keyword\n", column_name, *(long *)column_values);
      lvalues = (long *)column_values;
      for(i=0; i<number_values; i++) lvalues[i] = *(long *)column_values;
    } else if (datatype == TFLOAT) {
       headas_chat(5, "%s = %f read from keyword\n", column_name, *(float *)column_values);
      fvalues = (float *)column_values;
      for(i=0; i<number_values; i++) fvalues[i] = *(float *)column_values;
    } else {
      printf("%d is an unsupported datatype in SPII_read_col\n", datatype);
    }
    return(0);
  }

  /* get the column number */

  fits_get_colnum(fptr, CASEINSEN, column_name, &colnum, &status);
  if (status) {
    printf("***SPII_read_col: Failed to find column or keyword named %s\n", column_name);
    fits_report_error(stderr, status);
    printf("\n");
    return(status);
  }

  /* find out whether this is a scalar or vector column */

  fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
  if (status) {
    printf("***SPII_read_col: Failed to read information about the column named %s\n", column_name);
    fits_report_error(stderr, status);
    printf("\n");
    return(status);
  }

  /* read the column and replicate the values if the column is a scalar */

  fits_read_col(fptr, datatype, colnum, ispec, 1, repeat, NULL, column_values, &anynul, &status);
  if (status) {
    printf("***SPII_read_col: Failed to read column or keyword named %s\n", column_name);
    fits_report_error(stderr, status);
    printf("\n");
    return(status);
  }

  if (repeat == 1) {
    if (datatype == TINT) {
       headas_chat(5, "%s = %d read from scalar column %d\n", column_name, *(int *)column_values, colnum);
      ivalues = (int *)column_values;
      for(i=0; i<number_values; i++) ivalues[i] = *(int *)column_values;
    } else if (datatype == TLONG) {
       headas_chat(5, "%s = %ld read from scalar column %d\n", column_name, *(long *)column_values, colnum);
      lvalues = (long *)column_values;
      for(i=0; i<number_values; i++) lvalues[i] = *(long *)column_values;
    } else if (datatype == TFLOAT) {
       headas_chat(5, "%s = %f read from scalar column %d\n", column_name, *(float *)column_values, colnum);
      fvalues = (float *)column_values;
      for(i=0; i<number_values; i++) fvalues[i] = *(float *)column_values;
    } else {
      printf("%d is an unsupported datatype in SPII_read_col\n", datatype);
    }
    return(0);
  }

  headas_chat(5, "%s read from column %d\n", column_name, colnum);

  return(0);

}

/******************************** SP_write_key **************************************/

/* writes a keyword */

int SP_write_key(fitsfile *fptr, int datatype, char *keyword_name, void *keyword_value, char *comment)
{
  int status=0;

  fits_write_key(fptr, datatype, keyword_name, keyword_value, comment, &status);

  if (!status) {
    if (datatype == TFLOAT)  headas_chat(5, "Writing %s = %f\n", keyword_name, *(float *)keyword_value);
    if (datatype == TSTRING) headas_chat(5, "Writing %s = %s\n", keyword_name , (char *)keyword_value);
    if (datatype == TINT)    headas_chat(5, "Writing %s = %d\n", keyword_name, *(int *)keyword_value);
    if (datatype == TLOGICAL)headas_chat(5, "Writing %s = %d\n", keyword_name, *(int *)keyword_value);
    if (datatype == TLONG)   headas_chat(5, "Writing %s = %ld\n", keyword_name, *(long *)keyword_value);
  } else if (status) {
    if (datatype == TFLOAT) headas_chat(1, "***SP_write_key: Cannot write %s keyword as %f\n", keyword_name, *(float *)keyword_value);
    if (datatype == TSTRING) headas_chat(1, "***SP_write_key: Cannot write %s keyword as %s\n", keyword_name, (char *)keyword_value);
    if (datatype == TINT) headas_chat(1, "***SP_write_key: Cannot write %s keyword as %d\n", keyword_name, *(int *)keyword_value);
    if (datatype == TLOGICAL) headas_chat(1, "***SP_write_key: Cannot write %s keyword as %d\n", keyword_name, *(int *)keyword_value);
    if (datatype == TLONG) headas_chat(1, "***SP_write_key: Cannot write %s keyword as %ld\n", keyword_name, *(long *)keyword_value);
  }

  if (datatype != TFLOAT && datatype != TSTRING && datatype != TINT && datatype != TLONG && datatype != TLOGICAL)
    headas_chat(1, "%d is an unsupported datatype in SP_write_key\n", datatype);

  return(status);
}

/******************************** SPII_write_key **************************************/

/* writes a keyword if it is not included as a table model column */

int SPII_write_key(fitsfile *fptr, int tfields, char *ttype[], int datatype, char *keyword_name, void *keyword_value, char *comment)
{
  int icol, found;

  /* check whether this keyword is actually a column */

  found = 0;
  for (icol=0; icol<tfields; icol++) {
    if (strcmp(keyword_name, ttype[icol])) found = 1;
  }

  if (!found) SP_write_key(fptr, datatype, keyword_name, keyword_value, comment);

  return(found);

}

/******************************** SP_need_col **************************************/

/* returns true (1) if the array is not all the same value */

int SP_need_col(void *array, long array_size, int datatype)
{
  long i, lfirst;
  int ifirst;
  float first;
  int *iarray;
  long *larray;
  float *farray;

  if (datatype == TINT) {
    iarray = (int *)array;
    ifirst = iarray[0];
    for (i=1; i<array_size; i++) {
      if (ifirst != iarray[i]) {
	return(1);
      }
    }
  } else if (datatype == TLONG) {
    larray = (long *)array;
    lfirst = larray[0];
    for (i=1; i<array_size; i++) {
      if (lfirst != larray[i]) {
	return(1);
      }
    }
  } else if (datatype == TFLOAT) {
    farray = (float *)array;
    first = farray[0];
    for (i=1; i<array_size; i++) {
      if (first != farray[i]) {
	return(1);
      }
    }
  }

  return(0);

}

/******************************** SPII_set_table_str **************************************/

/* Find out whether a given character string needs to go in a column - if so set the table
   parameters and increment the tfield counter */

int SPII_set_table_str(char *column_name, char **column_values, long number_values, int *tfields, char *ttype[], char *tform[], char *tunit[])
{
  long i;
  int length, needcol=0;

  length = strlen(column_values[0]);
  for (i=1; i<number_values; i++) {
    if (strcmp(column_values[i],column_values[0])) needcol = 1;
    if (strlen(column_values[i]) > length) length = strlen(column_values[i]);
  }

  if (needcol) {
    (*tfields)++;
    strcpy(ttype[*tfields], column_name);
    sprintf(tform[*tfields], "A%d", length);
    strcpy(tunit[*tfields], " ");
  }

  return(needcol);

}  
  
/******************************** SP_write_col **************************************/

/* writes a column unless all the values are the same in which case writes a keyword
   note that the colnum argument is returned incremented by one if a column is written */

int SP_write_col(fitsfile *fptr, int datatype, char *column_name, int *colnum, long nrows, void *column_values)
{
  int status=0;

  if (SP_need_col(column_values, nrows, datatype)) {

    /* write the column */

    fits_write_col(fptr, datatype, *colnum, 1, 1, nrows, column_values, &status);
    if (status) {
       headas_chat(1, "Failed to write column %d named %s\n", *colnum, column_name);
      return(status);
    } else {
       headas_chat(5, "Successfully wrote column %d named %s\n", *colnum, column_name);
      (*colnum)++;
    }

  } else {

    /* just need to write the keyword */

    SP_write_key(fptr, datatype, column_name, column_values, NULL);

  }

  return(0);

}

/******************************** SPII_write_col **************************************/

/* writes the vector element to a column */

int SPII_write_col(fitsfile *fptr, int tfields, char *ttype[], int datatype, 
		   char *column_name, int rownum, long vector_length, void *column_values)
{
  int status=0;
  int colnum=0;
  int typecode;
  int i;
  long repeat, width;

  /* check whether there is a column for this name */

  for (i=0; i<tfields; i++) {
    if (!strcmp(column_name, ttype[i])) colnum = i+1;
  }

  if (colnum == 0) return(0);

  /* find the length of the vector that can be written to this column */

  fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
  if (status) {
     headas_chat(1, "Failed to find required length of column %d %s\n", colnum, column_name);
    return(status);
  }

  if ( repeat > vector_length ) repeat = vector_length;

  /* write the column */

  fits_write_col(fptr, datatype, colnum, rownum, 1, repeat, column_values, &status);
  if (status) {
     headas_chat(1, "Failed to write row %d to column %d %s\n", rownum, colnum, column_name);
    return(status);
  } else {
     headas_chat(5, "Successfully wrote row %d to column %d %s\n", rownum, colnum, column_name);
  }

  return(0);

}


