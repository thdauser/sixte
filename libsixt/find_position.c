#include "find_position.h"

PixPosition* getPixPosition()
{
  PixPosition* pp=NULL;

  //memory-allocation:
  pp=(PixPosition*)malloc(sizeof(PixPosition));
  if(pp!=NULL){
    pp->posX=0;
    pp->posY=0;
    pp->posRA=0.;
    pp->posDEC=0.;   
    pp->pixval=0.;
    /* pp->errorRA=0.;
    pp->errorDEC=0.;*/
  }

  return(pp);
}

PixPositionList* getPixPositionList()
{
  PixPositionList* ppl=NULL;
 

  //memory-allocation:
  ppl=(PixPositionList*)malloc(sizeof(PixPositionList));
  if(ppl!=NULL){
    ppl->entry=NULL;
    ppl->entryCount=0;
  }

  //allocate memory for the entry-pointer-array
  // size 100 is arbitrary chosen
 
    ppl->entry=(PixPosition*)malloc(100*sizeof(PixPosition));
 

  return(ppl);
}

double findBrightestPix(int threshold, SourceImage* sky_pixels, double pixval, PixPositionList* ppl,
		        float delta, double telRA, double telDEC)
{
  //get empty PixPosition structure (contains pos in x,y,ra,dec for all found sources,
  //as well as errors for ra,dec and the sources' pixval)
  /* PixPosition* pp=getPixPosition();*/
   ppl->entry[ppl->entryCount]=getPixPosition();

  int ii,jj,x,y;
  double pix=0.;

  //find brightest pixel
  for(ii=0; ii<sky_pixels->naxis1; ii++){
      for(jj=0; jj<sky_pixels->naxis2; jj++){

	//first, brightest source
	if(threshold==0){
	  if(sky_pixels->pixel[ii][jj] > pixval && sky_pixels->pixel[ii][jj] > pix){
	  //reset value of 'pix'
	  pix=sky_pixels->pixel[ii][jj];
	  //indices of brightest pixel
	  x=ii; y=jj;
	  }
	}else{
	//source below threshold
	    if(sky_pixels->pixel[ii][jj] > pix && sky_pixels->pixel[ii][jj] < pixval){
	      //reset value of 'pix'
	      pix=sky_pixels->pixel[ii][jj];
	      //indices of brightest pixel
	      x=ii; y=jj;
	    }
 	  }
      }
  }

  //ensure that pixval is still positive
  if(pix>0.){

  //TODO: check, whether enough memory is left for ppl->entry (if not -> realloc?)

  //save pixel coordinates to PixPositionList at entryCount
  ppl->entry[ppl->entryCount]->posX=x;
  ppl->entry[ppl->entryCount]->posY=y;
  ppl->entry[ppl->entryCount]->pixval=sky_pixels->pixel[x][y];

  ppl->entry[ppl->entryCount]->posRA=getPosRa(x,sky_pixels,delta,telRA); 
  ppl->entry[ppl->entryCount]->posDEC=getPosDec(y,sky_pixels,delta,telDEC);


  //error for found position (compared to real .simput-position)
  /*ppl->entry[ppl->entryCount].errorRA=detErrorRa();  //TODO
  ppl->entry[ppl->entryCount].errorDEC=detErrorDec();*/  //TODO

  //increase number of found sources
  ppl->entryCount++;

  //free PixPosition pointer
  // if(pp!=NULL){
  //free(pp);
  //}

  //set pixval to newly found brightest value
  pixval=sky_pixels->pixel[x][y];
  return(pixval);
  }//end pixval positive
  else{
    //free PixPosition pointer
    //  if(pp!=NULL){
    //free(pp);
    //}
    return(0.);
  }
}


double getPosRa(int pix, SourceImage* sky_pixels, float delta, double telRA)
{
  delta=delta*180/M_PI;
  double RA;
  double offAxis=(double)(((float)(pix)-sky_pixels->crpix1)*delta); //TODO:right now:middle of pix?!

  if((offAxis+telRA) < 0.){
     RA=360.+(offAxis+telRA);
  }
  else{
    if((offAxis+telRA)>360.){
      RA=(offAxis+telRA)-360.;
    }else{
      RA=offAxis+telRA;
    }
  }
 
  return(RA);
}


double getPosDec(int pix, SourceImage* sky_pixels, float delta, double telDEC)
{
  delta=delta*180/M_PI;
  double DEC;
  double offAxis=(double)(((float)(pix)-sky_pixels->crpix2)*delta);

  if((offAxis+telDEC) < -90.){
    DEC=(offAxis+telDEC+180.)*(-1);
  }
  else{
    if((offAxis+telDEC)>90.){
      DEC=fabs(offAxis+telDEC-180.);
    }else{
      DEC=offAxis+telDEC;
    }
  }
 
  return(DEC);
}


int getThresholdForSources(double pix, PixPositionList* ppl, SourceImage* sky_pixels)
{
  int th;
  double mean;

  //mean value of all pixels except already found sources
  mean=getMeanValue(ppl, sky_pixels); 
  //pixval:current candidate for bright pix (e.g. a source)
  //ensure that deviation between pix and all other (background) pixels is still big enough
  if((mean>0.) && (pix!=0.) && ((pix-mean) > 250*mean)) //TODO:value??
    {
      th=1;
    }else{
    th=2;
  };

  return(th);
}


double getMeanValue(PixPositionList* ppl, SourceImage* sky_pixels)
{
  double sum_all; //sum of all pixel-values of the sky image
  double sum_src; //sum of the pixel-values of already found sources
  double mean;
  int all, src;   // number of all/source pixels
  int ii, jj;     //counts

   for(ii=0; ii<sky_pixels->naxis1; ii++){
     for(jj=0; jj<sky_pixels->naxis2; jj++){
       sum_all+=sky_pixels->pixel[ii][jj];
     }
   }

   //'entryCount-1' since the count is incremented in fct. 'findBrightestPix'
   //and therefore PixPositionList has one element less
   for(ii=0; ii<(ppl->entryCount-1) ;ii++){
     sum_src+= ppl->entry[ii]->pixval;
   }

   all=sky_pixels->naxis1*sky_pixels->naxis2;
   src=ppl->entryCount-1;

   //mean value of all pixels except already identified sources
   mean=(sum_all-sum_src)/(all-src);
   return(mean);
}


void savePositionList(PixPositionList* ppl, char* filename, int* status)
{
  fitsfile* fptr=NULL;

  // Print information to STDOUT.
  char msg[MAXMSG];
  sprintf(msg, "Store PositionList in file '%s' ...\n", filename);
  headas_chat(5, msg);


  do{
    int count;
    
    //remove old version
    remove(filename);
    
    //number of columns
    int tfields=6; 
    //names of columns
    char* ttype[]={"src_num", "X", "Y", "RA", "DEC", "val"};
    //datatypes of columns
    char* tform[]={"V", "V", "V", "D", "D", "D"};
    char* tunit[]={" ", " ", " ", " ", " ", " "};
    //extension name
    char extname[]="src_positions";
    long firstrow=1;
    long firstelem=1;
    long nelem=ppl->entryCount;
    //arrays for columns
    int src_num_array[nelem];
    int x_array[nelem];
    int y_array[nelem];
    double ra_array[nelem];
    double dec_array[nelem]; 
    double val_array[nelem];

    //fill in the arrays
    for(count=0; count<nelem; count++){
      src_num_array[count]=count;
      x_array[count]=ppl->entry[count]->posX;
      y_array[count]=ppl->entry[count]->posY;
      ra_array[count]=ppl->entry[count]->posRA;
      dec_array[count]=ppl->entry[count]->posDEC;
      val_array[count]=ppl->entry[count]->pixval;
    }

    //create new FITS-file
    if(fits_create_file(&fptr, filename, status)) break;

    //create a table extension
    if(fits_create_tbl(fptr, BINARY_TBL, 1, tfields, ttype, tform, tunit, extname, status)) break;

    fits_write_col(fptr, TUINT, 1, firstrow, firstelem, nelem, src_num_array, status);
    fits_write_col(fptr, TUINT, 2, firstrow, firstelem, nelem, x_array, status);
    fits_write_col(fptr, TUINT, 3, firstrow, firstelem, nelem, y_array, status);
    fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nelem, ra_array, status);
    fits_write_col(fptr, TDOUBLE, 5, firstrow, firstelem, nelem, dec_array, status);
    fits_write_col(fptr, TDOUBLE, 6, firstrow, firstelem, nelem, val_array, status);

  }while(0);

  //close the FITS-file
  if(NULL!=fptr) fits_close_file(fptr, status);
}


void FreePixPositionList(PixPositionList* ppl)
{
  int count;

  if(ppl!=NULL){
    if(ppl->entry!=NULL){
      for(count=0; count<ppl->entryCount; count++){
	free(ppl->entry[count]);
      }
    }
    free(ppl->entry);
  }
    free(ppl);
}

