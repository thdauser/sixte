#include "fft_array.h"

fftw_complex* FFTOfArray_1d(double* Image1d, int ImageSize1, int ImageSize2)
{
  fftw_complex* Input;
  fftw_complex* Output;
  fftw_plan plan;
  int x;

  //Memory-Allocation
  Input=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ImageSize1*ImageSize2/*+ImageSize1*/));
  Output=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ImageSize1*ImageSize2/*+ImageSize1*/));

  //Copy Image1d of type double to Input-Array of type fftw_complex

  for(x=0; x<(/*ImageSize+1*/ImageSize1*ImageSize2);x++){
    Input[x][0]=Image1d[x];
    Input[x][1]=0.;
  }

  plan=fftw_plan_dft_2d(ImageSize1, ImageSize2, Input, Output, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(plan);

  return(Output);

  fftw_destroy(plan);
  fftw_free(Input); fftw_free(Output);
}

fftw_complex* FFTOfArrayInverse_1d(double* Image1d, int ImageSize1, int ImageSize2)
{
  fftw_complex* Input;
  fftw_complex* Output;
  fftw_plan plan;
  int x;

  //Memory-Allocation
  Input=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ImageSize1*ImageSize2);
  Output=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ImageSize1*ImageSize2);

  //Copy Image1d of type double to Input-Array of type fftw_complex

  for(x=0; x<(/*ImageSize1+*/ImageSize1*ImageSize2);x++){
    Input[x][0]=Image1d[x];
    Input[x][1]=0.;
  }

  plan=fftw_plan_dft_2d(ImageSize1, ImageSize2, Input, Output, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(plan);

  return(Output);

  fftw_destroy(plan);
  fftw_free(Input); fftw_free(Output);
}

fftw_complex* FFTOfArrayInverse_fftwcomplex(fftw_complex* Input, int ImageSize1, int ImageSize2)
{
  fftw_complex* Output;
  fftw_plan plan;
  int x;

  //Memory-Allocation
  Output=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ImageSize1*ImageSize2);

  plan=fftw_plan_dft_2d(ImageSize1, ImageSize2, Input, Output, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(plan);

  return(Output);

  fftw_destroy(plan);
  fftw_free(Input); fftw_free(Output);
}

void testFitsImagefft(fftw_complex* Image,char filename, int Size1, int Size2)
{
  int count;
  int status =0;    
  fitsfile* fptr;
  fits_create_file(&fptr, filename ,&status);
  long fpixel[2]={1,1};
  long naxes[2] = {(long)Size1, (long)Size2};
  fits_create_img(fptr, DOUBLE_IMG,2, naxes, &status);
  double* test;
  test= calloc(Size1*Size2, sizeof(double));
  for(count=0; count<Size1*Size2; count++){
     test[count]=Image[count][0];
   }
  fits_write_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1],test, &status);
  fits_close_file(fptr, &status);
}

void testFitsImage1d(double* Image1d,char filename, int Size1, int Size2)
{
  int count;
  int status =0;    
  fitsfile* fptr;
  fits_create_file(&fptr, filename ,&status);
  long fpixel[2]={1,1};
  long naxes[2] = {(long)Size1, (long)Size2};
  fits_create_img(fptr, DOUBLE_IMG,2, naxes, &status);
  double* test;
  test= calloc(Size1*Size2, sizeof(double));
  for(count=0; count<Size1*Size2; count++){
     test[count]=Image1d[count];
   }
  fits_write_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1],test, &status);
  fits_close_file(fptr, &status);
}
