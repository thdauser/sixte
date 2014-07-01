#include "maskshadow.h"

MaskShadow* getEmptyMaskShadowElement(int* const status)
{
  //memory-allocation:
  MaskShadow* ms=(MaskShadow*)malloc(sizeof(MaskShadow));
  if(ms!=NULL){
    ms->map=NULL;
    ms->shadow=NULL;
  }else{
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for MaskShadow Array!\n",
		   *status);
    return(ms);
  }

  return(ms);
}

MaskShadow* getMaskShadowElement(int const Size1, int const Size2, int* const status)
{
  MaskShadow* ms=NULL;
  int x,y;

  //get empty MaskShadow-element
  ms=getEmptyMaskShadowElement(status);
  if (EXIT_SUCCESS!=*status) return(ms);

  //memory-allocation map (as big as mask with pixel-size of detector)
  ms->map=(double**)malloc((Size1/2+1)*sizeof(double*));
  if(NULL!=ms->map){
    for(x=0; x <= Size1/2; x++){  //'+1' because of shift of map see below (get MaskShadow2)
      ms->map[x]=(double*)malloc((Size2/2+1)*sizeof(double));
	if(NULL==ms->map[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "MaskShadow map!\n", *status);
	 
	  return(ms);
	}
	//Clear the pixels
	for(y=0; y <= Size2/2; y++){
	  ms->map[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(ms);
  }else{
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "MaskShadow map!\n", *status);
      return(ms);
  }//end of memory-allocation map

   //memory-allocation shadow (as big as EventArray -> for later subtraction)
  ms->shadow=(double**)malloc(Size1*sizeof(double*));
  if(NULL!=ms->shadow){
    for(x=0; x < Size1; x++){
      ms->shadow[x]=(double*)malloc(Size2*sizeof(double));
	if(NULL==ms->shadow[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "MaskShadow shadow!\n", *status);
	 
	  return(ms);
	}
	//Clear the pixels
	for(y=0; y < Size2; y++){
	  ms->shadow[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(ms);
  }else{
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "MaskShadow shadow!\n", *status);
      return(ms);
  }//end of memory-allocation shadow


  return(ms);
}

void getMaskRepix(ReconArray* recon, MaskShadow* ms)
{
  int ii,jj;

  for(ii=1; ii<=recon->naxis1/2; ii++){
    for(jj=1; jj<=recon->naxis2/2; jj++){
      ms->map[ii][jj]=(recon->Rmap[ii][jj]+fabs((recon->open_fraction)/(recon->open_fraction - 1.)))/
	(1.-(recon->open_fraction)/(recon->open_fraction - 1.));
    }
  }
  //value at Rmap (ranging from max neg value defined by OF to 1)
  //divided by: distance from max neg value to 1 -> percentage of transparency (0==none to 1==max)
}


void getMaskShadow(MaskShadow* ms, PixPositionList* ppl, SourceImage* sky_pixels, SquarePixels* det_pix,
		    ReconArray* r, const Vector* const nx, const Vector* const ny, double const distance)
{  
  int ii,jj;
  int x,y;
  double* posD=NULL;
  Vector phodir={0.,0.,0.};

  //memory-allcation
  posD=(double*)malloc(2*sizeof(double));

  //get unit vector in photon direction from ra/dec
  double ra=ppl->entry[ppl->entryCount-1]->posRA;
  double dec=ppl->entry[ppl->entryCount-1]->posDEC;
  phodir=normalize_vector(unit_vector(ra*M_PI/180., dec*M_PI/180.));

  //components of the photon direction with respect to the
  //detector coordinate axes nx, ny (in mask plane)
  double x_comp = scalar_product(&phodir, nx);    //TODO: set up telescope attitude
  double y_comp = scalar_product(&phodir, ny);

  //component of the photon direction within the mask plane
  double radius = sqrt(pow(x_comp,2.)+pow(y_comp,2.));
  //azimuthal angle (with respect to the nx-axis)
  double alpha=atan2(y_comp, x_comp);

  for(ii=1; ii<r->naxis1/2; ii++){
    for(jj=1; jj<r->naxis2/2; jj++){
      
      posD[0]=(ii-sky_pixels->crpix1/2+0.5)*det_pix->xpixelwidth
	-radius*distance*cos(alpha)+det_pix->xwidth/2*det_pix->xpixelwidth;
      posD[1]=(jj-sky_pixels->crpix2/2+0.5)*det_pix->ypixelwidth
	-radius*distance*sin(alpha)+det_pix->ywidth/2*det_pix->ypixelwidth;

      if(posD[0]<0. || posD[0]>(det_pix->xwidth*det_pix->xpixelwidth) || posD[1]<0. 
	 || posD[1]>(det_pix->ywidth*det_pix->ypixelwidth)){
      }else{

	if(det_pix->DCU_length!=0.){//detector with gaps

	  double DCA_length=(2.*det_pix->DCU_length)+det_pix->DCU_gap;
	  double DCA=DCA_length+det_pix->DCA_gap;

	  int ratio_x=(int)(posD[0]/DCA);
	  int ratio_y=(int)(posD[1]/DCA);
	  
	  x=-1;y=-1;

	  if(posD[0] < (ratio_x*DCA+det_pix->DCU_length)) 
	    {x = (int)(posD[0]/det_pix->xpixelwidth +1.)-1;}
	  else{if(posD[0] > (ratio_x*DCA+det_pix->DCU_length+det_pix->DCU_gap)
		  && posD[0] < (ratio_x*DCA+DCA_length))
	      {x = (int)(posD[0]/det_pix->xpixelwidth +1.)-1;}
	    else{ }
	  }
	  if(posD[1] < (ratio_y*DCA+det_pix->DCU_length)) 
	    {y = (int)(posD[1]/det_pix->ypixelwidth +1.)-1;}
	  else{if(posD[1]> (ratio_y*DCA+det_pix->DCU_length+det_pix->DCU_gap)
		  && posD[1]< (ratio_y*DCA+DCA_length))
	      {y = (int)(posD[1]/det_pix->ypixelwidth +1.)-1;}
	    else{ }
	  }
	  
	  if(x>=0 && y>=0){
	    ms->shadow[x][y]=ms->map[ii][jj]; //-1 for testing source at 0/0
	    //TODO: CHECK whether really needed?? (also change appropriately at memory allocation)
	  }

	}else{//detector without gaps
	
	  x=(int)(posD[0]/det_pix->xpixelwidth+1.)-1;
	  y=(int)(posD[1]/det_pix->ypixelwidth+1.)-1;
	  ms->shadow[x][y]=ms->map[ii][jj];

	}//detector with or without gaps
      }//does not hit the walls

    }//jj
  }//ii

}

void getMaskShadow2(MaskShadow* ms, struct wcsprm* wcs2, PixPositionList* ppl, SourceImage* sky_pixels, SquarePixels* det_pix, ReconArray* r, int* const status)
{
  double pixcrd[2];
  double imgcrd[2];
  double world[2]={ppl->entry[ppl->entryCount-1]->posRA, ppl->entry[ppl->entryCount-1]->posDEC};
  double phi, theta;
  int stat=0;

  wcss2p(wcs2,1,2,world,&phi,&theta,imgcrd,pixcrd,&stat);
  if(0!=stat){
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("wcs projection failed!\n", *status);
  }

  int ii,jj;
  int x,y;
  double* posD=NULL;

   //memory-allcation
  posD=(double*)malloc(2*sizeof(double));

  for(ii=1; ii<=r->naxis1/2; ii++){
    for(jj=1; jj<=r->naxis2/2; jj++){
      
      //in meters
      posD[0]=(ii-sky_pixels->crpix1/2+0.5)*det_pix->xpixelwidth+(double)(pixcrd[0]*det_pix->xpixelwidth);
      posD[1]=(jj-sky_pixels->crpix2/2+0.5)*det_pix->ypixelwidth+(double)(pixcrd[1]*det_pix->ypixelwidth);

      if(posD[0]<0. || posD[0]>(det_pix->xwidth*det_pix->xpixelwidth) || posD[1]<0. 
	 || posD[1]>(det_pix->ywidth*det_pix->ypixelwidth)){
      }else{

	if(det_pix->DCU_length!=0.){//detector with gaps

	  double DCA_length=(2.*det_pix->DCU_length)+det_pix->DCU_gap;
	  double DCA=DCA_length+det_pix->DCA_gap;

	  int ratio_x=(int)(posD[0]/DCA);
	  int ratio_y=(int)(posD[1]/DCA);
	  
	  x=-1;y=-1;

	  if(posD[0] < (ratio_x*DCA+det_pix->DCU_length)) 
	    {x = (int)(posD[0]/det_pix->xpixelwidth +1.)-1;}
	  else{if(posD[0] > (ratio_x*DCA+det_pix->DCU_length+det_pix->DCU_gap)
		  && posD[0] < (ratio_x*DCA+DCA_length))
	      {x = (int)(posD[0]/det_pix->xpixelwidth +1.)-1;}
	    else{ }
	  }
	  if(posD[1] < (ratio_y*DCA+det_pix->DCU_length)) 
	    {y = (int)(posD[1]/det_pix->ypixelwidth +1.)-1;}
	  else{if(posD[1]> (ratio_y*DCA+det_pix->DCU_length+det_pix->DCU_gap)
		  && posD[1]< (ratio_y*DCA+DCA_length))
	      {y = (int)(posD[1]/det_pix->ypixelwidth +1.)-1;}
	    else{ }
	  }
	  
	  if(x>=0 && y>=0){
	    ms->shadow[x][y]=ms->map[ii-1][jj-1]; //-1 maybe because of wcs starts counting at 1? or at pixcrd?
	    //TODO: CHECK whether really needed?? (also change appropriately at memory allocation)
	  }

	}else{//detector without gaps
	
	  x=(int)(posD[0]/det_pix->xpixelwidth+1.)-1;
	  y=(int)(posD[1]/det_pix->ypixelwidth+1.)-1;
	  ms->shadow[x][y]=ms->map[ii][jj];

	}//detector with or without gaps
      }//does not hit the walls

    }//jj
  }//ii

  if(posD!=NULL){
    free(posD);
  } 
}

void testImageShadow(MaskShadow* ms, SquarePixels* det_pix, int* status)
{
  double* buffer1d=NULL;

  //Memory-Allocation for 1d-image
 buffer1d = (double*)malloc(det_pix->xwidth*det_pix->ywidth*sizeof(double));
 if (!buffer1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-buffer!\n", *status);
 }

 //Create the 1D-image from buffer 
  int x, y;
  for (x=0; x<det_pix->xwidth; x++) {
    for (y=0; y<det_pix->ywidth; y++) {
      buffer1d[(x+det_pix->xwidth*y)] = ms->shadow[x][y];
   }
  }

  testFitsImage1d(buffer1d, "maskshadow.fits",det_pix->xwidth, det_pix->ywidth, &status);
}

void testImageEventArray(ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff, 
			  char* filename, int* status)
{
  double* buffer1d=NULL;

  //Memory-Allocation for 1d-image
 buffer1d = (double*)malloc(det_pix->xwidth*det_pix->ywidth*sizeof(double));
 if (!buffer1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-buffer!\n", *status);
 }

 //Create the 1D-image from buffer 
  int x, y;
  for (x=0; x<det_pix->xwidth; x++) {
    for (y=0; y<det_pix->ywidth; y++) {
      buffer1d[(x+det_pix->xwidth*y)] = ea->EventArray[x+ea->naxis1/2+xdiff][y+ea->naxis2/2+ydiff];
   }
  }

  testFitsImage1d(buffer1d, filename, det_pix->xwidth, det_pix->ywidth, &status);
}


void testImageMap(MaskShadow* ms, ReconArray* r, int* status)
{
  double* buffer1d=NULL;

  //Memory-Allocation for 1d-image
 buffer1d = (double*)malloc(r->naxis1/2*r->naxis2/2*sizeof(double));
 if (!buffer1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-buffer!\n", *status);
 }

 //Create the 1D-image from buffer 
  int x, y;
  for (x=0; x<r->naxis1/2; x++) {
    for (y=0; y<r->naxis2/2; y++) {
      buffer1d[(x+r->naxis1/2*y)] = ms->map[x][y];
   }
  }

  testFitsImage1d(buffer1d, "maskmap.fits",r->naxis1/2, r->naxis2/2, &status);
}


double getNormalization1(MaskShadow* ms, ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff)
{
  int ii,jj;
  double count=0.,sum=0.;

  for (ii=0; ii<det_pix->xwidth; ii++) {
    for (jj=0; jj<det_pix->ywidth; jj++) {
      if(ms->shadow[ii][jj]==0){ //fully opaque
	count+=1.;
	sum+=ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff];
      }else{ //partially opaque
	if(ms->shadow[ii][jj]<1 && ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff] != 0.){
	  count+=(1-ms->shadow[ii][jj]);
	  sum+=(1-ms->shadow[ii][jj])*(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]);
	}	
      }    
    }
  }

  return(sum/count);
}


double getNormalization2(MaskShadow* ms, ReadEvent* ea, SquarePixels* det_pix, int const xdiff, int const ydiff)
{
  int ii,jj;
  long double sum1=0., sum2=0.;

  for (ii=0; ii<det_pix->xwidth; ii++) {
    for (jj=0; jj<det_pix->ywidth; jj++) {
      if(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]>0.){
      sum1+=ms->shadow[ii][jj];
      }
    }
  }

  for (ii=0; ii<det_pix->xwidth; ii++) {
    for (jj=0; jj<det_pix->ywidth; jj++) {
      if(ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff]>0.){
      sum2+=(ms->shadow[ii][jj]*ms->shadow[ii][jj])/ea->EventArray[ii+ea->naxis1/2+xdiff][jj+ea->naxis2/2+ydiff];
      }
    }
  }

  if(sum2!=0.){
    return(sum1/sum2);
  }else{
    return(0.);
  }
}


void FreeMaskShadow(MaskShadow* ms,int const Size1)
{
  int count;

  if(ms->map!=NULL){
    for(count=0; count<(Size1/2+1); count++){
      if(ms->map[count]!=NULL){
       free(ms->map[count]);
      }
    }
    free(ms->map);
  }

  if(ms->shadow!=NULL){
    for(count=0; count<Size1; count++){
      if(ms->shadow[count]!=NULL){
       free(ms->shadow[count]);
      }
    }
    free(ms->shadow);
  }

  if(ms!=NULL){
    free(ms);
  }
}
