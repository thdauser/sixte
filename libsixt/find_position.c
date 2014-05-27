#include "find_position.h"

PixPosition* getPixPosition()
{
  PixPosition* pp=NULL;

  //memory-allocation:
  pp=(PixPosition*)malloc(sizeof(PixPosition));
  if(pp!=NULL){
    pp->midPixX=0;
    pp->midPixY=0;
    pp->posX=0.;
    pp->posY=0.;
    pp->posRA=0.;
    pp->posDEC=0.;   
    pp->pixval=0.;
  }

  return(pp);
}


PixPositionList* getPixPositionList(SourceImage* sky_pixels)
{
  PixPositionList* ppl=NULL;
  int ii, jj;

  //memory-allocation:
  ppl=(PixPositionList*)malloc(sizeof(PixPositionList));
  if(ppl!=NULL){
    ppl->entry=NULL;
    ppl->entryCount=0;
    ppl->neighbors=NULL;
    ppl->found_pos=NULL;
  }

  //allocate memory for the entry-pointer-array
  //size 2 is arbitrary chosen
 
  ppl->entry=(PixPosition**)malloc(2*sizeof(PixPosition*));
  ppl->neighbors=(SourceNeighbors**)malloc(2*sizeof(SourceNeighbors*));

  //get temporary array of all pixels except already found sources (including neighbors)
  //memory-allocation: as big as the sky image
  ppl->found_pos=(int**)malloc(sky_pixels->naxis1*sizeof(int*));
  if (NULL!=ppl->found_pos) {
      for(ii=0; ii<sky_pixels->naxis1; ii++) {
	ppl->found_pos[ii] = (int*)malloc(sky_pixels->naxis2*sizeof(int));
	//Clear the pixels
	for(jj=0; jj < sky_pixels->naxis2; jj++){
	  ppl->found_pos[ii][jj]=0;
	}
      }
  }
 
  return(ppl);
}


double* getMedian_list(SourceImage* sky_pixels, int* const status)
{
  double* ml=NULL;
  //memory-allocation: as big as the sky image without found sources + neighbors
  ml=(double*)malloc((sky_pixels->naxis1*sky_pixels->naxis2)*sizeof(double));
   //median_list=(double*)malloc(((sky_pixels->naxis1*sky_pixels->naxis2)-(ppl->entryCount*9))*sizeof(double));
  if (NULL==ml) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: Allocation of median list failed!\n", *status);
	}

  return(ml);
}


SourceNeighbors* getSourceNeighbors()
{
  SourceNeighbors* sn=NULL;

  //memory-allocation:
  sn=(SourceNeighbors*)malloc(sizeof(SourceNeighbors));
  if(sn!=NULL){
    sn->neighbor_list=NULL;
    sn->neighborAmount=0;
  }
  return(sn);
}


double findBrightestPix(int threshold, SourceImage* sky_pixels, double pixval, PixPositionList* ppl,
		        struct wcsprm* wcs, int* const status)
{
  int ii,jj,x,y;
  double pix=0.; //temporary store pixval for comparison with current sky-image pixel to get the brightest

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
	}else{//source below threshold

	     //find next brightest pix
	     if(sky_pixels->pixel[ii][jj] > pix && sky_pixels->pixel[ii][jj] < pixval){	       
	       if(ppl->found_pos[ii][jj]==0.){  //which is the case if the pixel is not identified yet
	       //reset value of 'pix'
	       pix=sky_pixels->pixel[ii][jj];
	       //indices of brightest pixel
	       x=ii; y=jj;
	     }
	   }

	}//END source below threshold (not first source)

      }//END sky_pixels jj
  }//END sky_pixels ii

      if(ppl->entryCount > 1){
	ppl->entry=(PixPosition**)realloc(ppl->entry, (ppl->entryCount+1)*sizeof(PixPosition*));
	if (NULL==ppl->entry) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error:Reallocation of PixPosition entry failed!\n", *status);
	}else{
	ppl->neighbors=(SourceNeighbors**)realloc(ppl->neighbors, (ppl->entryCount+1)*sizeof(SourceNeighbors*));
	if (NULL==ppl->neighbors) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error:Reallocation of PixPosition entry failed!\n", *status);
	}
	}
      }
      //get empty PixPosition structure (contains pos in x,y,ra,dec for all found sources,
      //as well as errors for ra,dec and the sources' pixval)
      ppl->entry[ppl->entryCount]=getPixPosition();

      ppl->entry[ppl->entryCount]->midPixX=x;
      ppl->entry[ppl->entryCount]->midPixY=y;

      //identify and save neighbors;find pos from significant pix-> save to PixPosList at current entry
      findNeighbors(x, y, ppl, sky_pixels, wcs, &status);

      //increase number of found sources
      ppl->entryCount++;

      //set pixval to newly found brightest value
      pixval=sky_pixels->pixel[x][y];

   return(pixval);
}


void findNeighbors(int x, int y, PixPositionList* ppl, SourceImage* sky_pixels, struct wcsprm* wcs, int* const status)
{ 
  int ii, jj; //counts
  double sum_all=0.; //sum of all pixel-values of the sky image
  double src_x=0.,src_y=0.; //weighted mean in x/y direction of all pix
  int xstart,ystart,xamount,yamount;
  double** snl=NULL;
  ppl->neighbors[ppl->entryCount]=getSourceNeighbors();

  //memory-allocation for the neighbor_list at entryCount (current brightest pix)
  //16 pix are allocated even though some sources only contain 9 or 12 pix per source
   snl=(double**)malloc(4*sizeof(double*));
   if (NULL!= snl) {
      for(ii=0; ii<4; ii++) {
	 snl[ii] = (double*)malloc(4*sizeof(double));

	//Clear the pixels
	for(jj=0; jj<4; jj++){
	   snl[ii][jj]=0.;
	}
      }
    }  

  //determine pixel amount in each direction
   //x-direction
   if((sky_pixels->pixel[x+1][y]/sky_pixels->pixel[x][y]) >= 0.8){ //right hand side from brightest pix is also bright
     xamount=4;
     xstart=1;
   }else{
     if((sky_pixels->pixel[x-1][y]/sky_pixels->pixel[x][y]) >= 0.8){ //left hand side from brightest pix is also bright
       xamount=4;
       xstart=2;
     }else{
       xamount=3;
       xstart=1;
     }
   }

   //y-direction
   if((sky_pixels->pixel[x][y+1]/sky_pixels->pixel[x][y]) >= 0.8){ //top of brightest pix is also bright
     yamount=4;
     ystart=1;
   }else{
     if((sky_pixels->pixel[x][y-1]/sky_pixels->pixel[x][y]) >= 0.8){ //bottom of brightest pix is also bright
       yamount=4;
       ystart=2;
     }else{
       yamount=3;
       ystart=1;
     }
   }


  //find all neighbours of current brightest pix
  //'-xstart' since: start before current pix
   for(ii=0; ii<xamount; ii++){
    for(jj=0; jj<yamount; jj++){
      snl[ii][jj]=sky_pixels->pixel[x-xstart+ii][y-ystart+jj];
      //value of 'found_pos' equals (1D) index of sky image (zero else) to be able to compare both
      ppl->found_pos[x-xstart+ii][y-ystart+jj]=(x-xstart+ii)+sky_pixels->naxis1*(y-ystart+jj);
    }
  }

   ppl->neighbors[ppl->entryCount]->neighborAmount=xamount*yamount;

  //determine sum of all counts within the neigbor pix to weight each pix
  for(ii=0; ii<xamount; ii++){
    for(jj=0; jj<yamount; jj++){
      sum_all+= snl[ii][jj];   //TODO: only significant pix
     }
   }

  //weight each pix and determine weighted position
  for(ii=0; ii<xamount; ii++){
    for(jj=0; jj<yamount; jj++){
      src_x+= snl[ii][jj]/sum_all*(x-xstart+ii);
      src_y+= snl[ii][jj]/sum_all*(y-ystart+jj);   //TODO: only significant pix
     }
   }

  ppl->neighbors[ppl->entryCount]->neighbor_list=snl;

  //save pixel coordinates to PixPositionList at entryCount
  ppl->entry[ppl->entryCount]->posX=src_x;
  ppl->entry[ppl->entryCount]->posY=src_y;
  //set pixval to sum of all counts of contributing pixels //TODO: ??
  ppl->entry[ppl->entryCount]->pixval=sum_all;

  //get RA and DEC value of current source with wcs fct. (tangent projection)
  double pixcrd[2]={src_x+1,src_y+1};
  double imgcrd[2];
  double world[2];
  double phi, theta;
  int stat=0;
  wcsp2s(wcs,1,2,pixcrd,imgcrd,&phi,&theta,world,&stat);
  if(0!=stat){
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("wcs projection failed!\n", *status);
  }

  ppl->entry[ppl->entryCount]->posRA=world[0]; 
  ppl->entry[ppl->entryCount]->posDEC=world[1];
}


int getThresholdForSources(double pix, PixPositionList* ppl, SourceImage* sky_pixels,
			   double* median_list, double factorOfSigma)
{
  int th;
  int ii,jj;
  double median=0., sigma=0.;
  long unsigned int median_list_count=0; 
 
  for(ii=0; ii<sky_pixels->naxis1; ii++){
    for(jj=0; jj<sky_pixels->naxis2; jj++){
      if((ii+sky_pixels->naxis1*jj)!=ppl->found_pos[ii][jj]){  //which is the case if the pixel is not identified yet
	median_list[median_list_count]=sky_pixels->pixel[ii][jj];
	median_list_count++;
      }
    }
    }
  
  //median value of all pixels except already found sources
  //median=torben(median_list, median_list_count); //TODO: use other method (quick select?!)
  median=getMean(median_list, median_list_count);
  //pixval:current candidate for bright pix (e.g. a source)
  //ensure that deviation between pix and all other (background) pixels is still big enough
  sigma=getSdev(median_list, median_list_count); //standard deviation
  if((pix!=0.) && ((pix-median) > factorOfSigma*sigma))
    { 
     th=1;
    }else{
     th=2;
     FreeLists(ppl->found_pos, median_list, sky_pixels);
  };

  return(th);
}

double getMean(double* median_list,long unsigned int n)
{
  long unsigned int ii;
  double sum=0.;
  for(ii=0; ii<n; ii++){
    sum+=median_list[ii];
  }
  return(sum/n);
}


double getSdev(double* median_list,long unsigned int n)
{
  int ii;
  double sum=0., mean;
  double dev=0., sigma;

  for(ii=0; ii<n; ii++){
    sum+=median_list[ii];
  }

  mean=sum/(double)(n);

  for(ii=0; ii<n; ii++){
    dev+=(median_list[ii]-mean)*(median_list[ii]-mean);
  }

  sigma=sqrt(dev/((double)(n)-1));
  return(sigma);
}


 /*******************************************************************
  ** The following code is public domain.                          **
  ** Algorithm by Torben Mogensen, implementation by N. Devillard. **
  ** This code in public domain.                                   **
  ******************************************************************/

//TODO:give array of all pixels except already identified ones
double torben(double* m,long unsigned int n)
{
    int         i, less, greater, equal;
    double  min, max, guess, maxltguess, mingtguess;

    min = max = m[0] ;
    for (i=1 ; i<n ; i++) {
        if (m[i]<min) min=m[i];
        if (m[i]>max) max=m[i];
    }

    while (1) {
        guess = (min+max)/2;
        less = 0; greater = 0; equal = 0;
        maxltguess = min ;
        mingtguess = max ;
        for (i=0; i<n; i++) {
            if (m[i]<guess) {
                less++;
                if (m[i]>maxltguess) maxltguess = m[i] ;
            } else if (m[i]>guess) {
                greater++;
                if (m[i]<mingtguess) mingtguess = m[i] ;
            } else equal++;
        }
        if (less <= (n+1)/2 && greater <= (n+1)/2) break ; 
        else if (less>greater) max = maxltguess ;
        else min = mingtguess;
    }
    if (less >= (n+1)/2) return maxltguess;
    else if (less+equal >= (n+1)/2) return guess;
    else return mingtguess;
}



void savePositionList(PixPositionList* ppl, char* filename, int* status)
{
  fitsfile* fptr=NULL;

  // Print information to STDOUT.
  char msg[MAXMSG];
  sprintf(msg, "Store PositionList in file '%s' ...\n", filename);
  headas_chat(5, msg);


  do{
    long unsigned int count;
    
    //remove old version
    remove(filename);
    
    //number of columns
    int tfields=5; 
    //names of columns
    char* ttype[]={"X", "Y", "RA", "DEC", "val"};
    //datatypes of columns
    char* tform[]={"D", "D", "D", "D", "D"};
    char* tunit[]={" ", " ", " ", " ", " "};
    //extension name
    char extname[]="src_positions";
    long firstrow=1;
    long firstelem=1;
    long unsigned int nelem=ppl->entryCount;
    //arrays for columns
    double x_array[nelem-1];
    double y_array[nelem-1];
    double ra_array[nelem-1];
    double dec_array[nelem-1]; 
    double val_array[nelem-1];

    //fill in the arrays
    for(count=0; count<nelem; count++){
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

    fits_write_col(fptr, TDOUBLE, 1, firstrow, firstelem, nelem, x_array, status);
    fits_write_col(fptr, TDOUBLE, 2, firstrow, firstelem, nelem, y_array, status);
    fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, nelem, ra_array, status);
    fits_write_col(fptr, TDOUBLE, 4, firstrow, firstelem, nelem, dec_array, status);
    fits_write_col(fptr, TDOUBLE, 5, firstrow, firstelem, nelem, val_array, status);

  }while(0);

  //close the FITS-file
  if(NULL!=fptr) fits_close_file(fptr, status);
}


void SaveSkyImage3Columns(SourceImage* si,char* filename,int* status)
{
   fitsfile* fptr=NULL;
   long unsigned int *x_array = NULL;
   long unsigned int *y_array = NULL;
   double *z_array = NULL;

   do{
     long unsigned int count=0;
     int x,y;
     
     //remove old version
    remove(filename);
    
    //number of columns
    int tfields=3; 
    //names of columns
    char* ttype[]={"X", "Y", "val"};
     //datatypes of columns
    char* tform[]={"V", "V", "D"};
    char* tunit[]={" ", " ", " "};
    //extension name
    char extname[]="sky_image";
    long firstrow=1;
    long firstelem=1;
    long unsigned int nelem=si->naxis1*si->naxis2;
    //arrays for columns
    x_array=(long unsigned int*)malloc((nelem)*sizeof(long unsigned int));
    y_array=(long unsigned int*)malloc((nelem)*sizeof(long unsigned int));
    z_array=(double*)malloc((nelem)*sizeof(double));
    for(x=0; x<nelem; x++){
      x_array[x]=0;
      y_array[x]=0;
      z_array[x]=0.;
    }    

    //fill in the arrays
    for(x=0; x<si->naxis1; x++){
      for(y=0; y<si->naxis2; y++){
	x_array[count]=x;
	y_array[count]=y;
	z_array[count]=si->pixel[x][y];
	count++;
      }
    }

    //create new FITS-file
    if(fits_create_file(&fptr, filename, status)) break;

    //create a table extension
    if(fits_create_tbl(fptr, BINARY_TBL, 1, tfields, ttype, tform, tunit, extname, status)) break;

    fits_write_col(fptr, TULONG, 1, firstrow, firstelem, nelem, x_array, status);
    fits_write_col(fptr, TULONG, 2, firstrow, firstelem, nelem, y_array, status);
    fits_write_col(fptr, TDOUBLE, 3, firstrow, firstelem, nelem, z_array, status);

   }while(0);

   //close the FITS-file
  if(NULL!=fptr) fits_close_file(fptr, status);

  if(x_array!=NULL){
    free(x_array);
  }
  if(y_array!=NULL){
    free(y_array);
  }
  if(z_array!=NULL){
    free(z_array);
  }
}


void FreePixPositionList(PixPositionList* ppl)
{
  long unsigned int count;
  int ii;

  if(ppl!=NULL){

    if(ppl->entry!=NULL){
      for(count=0; count<ppl->entryCount; count++){
	free(ppl->entry[count]);
      }
      free(ppl->entry);
    }
    

    if(ppl->neighbors!=NULL){
      for(count=0; count<ppl->entryCount; count++){
	if(ppl->neighbors[count]->neighbor_list!=NULL){

	  for(ii=0; ii<4; ii++){
	  free(ppl->neighbors[count]->neighbor_list[ii]);
	  }
	  free(ppl->neighbors[count]->neighbor_list);

	}
	free(ppl->neighbors[count]);
      }
      free(ppl->neighbors);
    }

    free(ppl);
  }
}


void FreeLists(int** found_pos, double*  median_list, SourceImage* sky_pixels)
{
  int count;

  if(found_pos!=NULL){
    for(count=0; count<sky_pixels->naxis1; count++){
      if(found_pos[count]!=NULL){
       free(found_pos[count]);
      }
    }
    free(found_pos);
  }

  if(median_list!=NULL){
    free(median_list);
  }
}
