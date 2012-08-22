/* prototypes for I/O wrap-up routines */

#include "fitsio.h"
#include "headas.h"

int SP_read_key(fitsfile *fptr, int datatype, char *keyword_name, void *keyword_value, void *keyword_default);

int SPII_read_key(fitsfile *fptr, int datatype, char *keyword_name, long ispec, void *keyword_value, void *keyword_default);

int SP_read_col(fitsfile *fptr, int datatype, char *column_name, long number_values, void *column_values);

int SPII_read_col(fitsfile *fptr, int datatype, char *column_name, long ispec, long number_values, void *column_values);

int SP_write_key(fitsfile *fptr, int datatype, char *keyword_name, void *keyword_value, char *comment);

int SPII_write_key(fitsfile *fptr, int tfields, char *ttype[], int datatype, char *keyword_name, void *keyword_value, char *comment);

int SP_need_col(void *array, long array_size, int datatype);

int SPII_set_table_str(char *column_name, char **column_values, long number_values, int *tfields, char *ttype[], char *tform[], char *tunit[]);

int SP_write_col(fitsfile *fptr, int datatype, char *column_name, int *colnum, long nrows, void *column_values);

int SPII_write_col(fitsfile *fptr, int tfields, char *ttype[], int datatype, char *column_name, int rownum, long vector_length, void *column_values);

