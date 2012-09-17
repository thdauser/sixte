int HDcmp(const void * data1, const void * data2);

void HDfuncs( double x, double *afunc, unsigned int ma );

void HDpoly_fit(double * x, double * y, double * c, int n, int degree);

void HDsmooth(float * input, float * output, int num, int width);

void HDsort(float * base, int * index, int n);

void HDsvbksb( double **U, double *W, double **V, unsigned int M,
    unsigned int N, unsigned int MP, unsigned int NP,
    double *B, double *X );

void HDsvdcmp( double **A, unsigned int M, unsigned int N,
    unsigned int MP, unsigned int NP, double *W, double **V );

void HDsvdfit( double *X, double *Y, double *Sig, unsigned int NData,
    double *A, unsigned int MA,
    double **U, double **V, double *W, unsigned int MP, unsigned int NP,
    double *ChiSq, void funcs(double x, double *afunc, unsigned int ma) );

void HDsvdvar( double **V, unsigned int MA, unsigned int NP,
    double *W, double **CVM, unsigned int NCVM );
