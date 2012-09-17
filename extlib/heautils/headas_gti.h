
/* Merging mode for gtimerge() */
#define GTI_AND 1
#define GTI_OR  2

/* Good time interval structure */
struct gti_struct {
  double mjdref;
  double timezero;
  int ngti, maxgti;
  double *start, *stop;
  void *dptr;
};

/* Forward declarations of library functions */
extern int HDgti_init(struct gti_struct *gti);
extern int HDgti_free(struct gti_struct *gti);
extern int HDgti_copy(struct gti_struct *dest, struct gti_struct *src, 
		   int *status);
extern int HDgti_grow(struct gti_struct *gti, int new, int *status);
extern int HDgti_read(char *filename, struct gti_struct *gti, 
	    char *extname, char *start, char *stop,
	    struct gti_struct *refer_to, 
	   fitsfile **fptr, int *status);
extern int HDgti_write(fitsfile *fptr, struct gti_struct *gti, 
	     char *extname, char *start, char *stop,
	    int *status);
extern int HDgti_merge(int mode, struct gti_struct *gti, 
	     struct gti_struct *agti, struct gti_struct *bgti, 
	     int *status);
extern int HDgti_clean(struct gti_struct *gti, struct gti_struct *ogti,
	     int *status);
extern double HDgti_exp(double t1, double t2, struct gti_struct *gti, 
		     int *status);
extern int HDgti_where(struct gti_struct *gti, int ntimes, 
	     double *times, int *segs, int *status);

extern double HDget_frac_time(fitsfile *fileptr, char *key, double *keyi, 
			  double *keyf, int *status);
