#include "testimg.h"

void createTestImg(void* arg, int type, int sizeX, int sizeY,
		    int shiftX, int shiftY, char* filename, int* const status)
{
  double** data=NULL;

  if(type==TMASKSHADOW){
    MaskShadow* ms =(MaskShadow*)arg;
    data=ms->shadow;
  }else if(type==TREADEVENT){
    ReadEvent* ea =(ReadEvent*)arg;
    data=ea->EventArray;
  }else if(type==TMASKMAP){
    MaskShadow* ms =(MaskShadow*)arg;
    data=ms->map;
  }else if(type==TPROJMASK){
    ProjectedMask* pm =(ProjectedMask*)arg;
    data=pm->map;
  }

  double* buffer1d=NULL;

  //Memory-Allocation for the buffer of 1d-image
  buffer1d = (double*)malloc(sizeX*sizeY*sizeof(double));
  if (!buffer1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-buffer!\n", *status);
  }

  //Create the 1d-image from buffer 
  int x, y;
  for (x=0; x<sizeX; x++) {
    for (y=0; y<sizeY; y++) {
      buffer1d[(x+sizeX*y)] = data[x+shiftX][y+shiftY];
   }
  }

  //Save 1d-image as FITS-file
  testFitsImage1d(buffer1d, filename, sizeX, sizeY, status);

  // Release memory from image input buffer.
  if (NULL!=buffer1d) free(buffer1d);
} 


void testFitsImage1d(double* Image1d, char* filename, int Size1, int Size2, int* const status)
{
  int count;   
  fitsfile* fptr;

  // Check if the file already exists.
  int exists;
  fits_file_exists(filename, &exists, status);
  if (0!=exists) {
    // Delete the file.
    remove(filename);
    }

  fits_create_file(&fptr, filename, status);
  long fpixel[2]={1,1};
  long naxes[2] = {(long)Size1, (long)Size2};
  fits_create_img(fptr, DOUBLE_IMG,2, naxes, status);
  double* test;
  test= calloc(Size1*Size2, sizeof(double));
  for(count=0; count<Size1*Size2; count++){
     test[count]=Image1d[count];
   }
  fits_write_pix(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], test, status);

  // Close the FITS file.
  if (NULL!=fptr) fits_close_file(fptr, status);
}


void createTestImg_fft_part(fftw_complex* Image, char* filename, int type, int Size1, int Size2, int* const status)
{
  int count;   
  fitsfile* fptr;
 
  long fpixel[2]={1,1};
  long naxes[2] = {(long)Size1, (long)Size2};

  fits_create_file(&fptr,filename,status);

  fits_create_img(fptr,DOUBLE_IMG,2,naxes,status);
  double* test;
  test= calloc(Size1*Size2, sizeof(double));
  for(count=0; count<Size1*Size2; count++){
      test[count]=Image[count][type];    
  }
  fits_write_pix(fptr,TDOUBLE,fpixel,naxes[0]*naxes[1],test,status);
  fits_close_file(fptr,status);
}
