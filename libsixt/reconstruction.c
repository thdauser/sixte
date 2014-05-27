#include "reconstruction.h"

ReconArray* newReconArray(int* const status)
{
  ReconArray* recon=(ReconArray*)malloc(sizeof(ReconArray));
  if (NULL==recon) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Reconstruction Array!\n",
		   *status);
    return(recon);
  }

  //Initialization:
  recon->Rmap=NULL;
  recon->open_fraction=0.;

  recon->naxis1 = 0;
  recon->naxis2 = 0;

  return(recon);
}


ReconArray* getReconArray(const CodedMask* const mask, SquarePixels* detector_pixels, int* const status)
{

  ReconArray* recon=NULL;
  int x,y;                          //count for memory allocation
  int xcount, ycount;               //count for getting Rmap in case1:same pixel size
  int xpixelcount=0, ypixelcount=0; //count for getting Rmap in case2:diff pixel size
  double leftsmall=0.,leftbig=0.;   //left border of small(detector) and big(mask) pixel
  double topsmall=0.,topbig=0.;     //top border of small(detector) and big(mask) pixel
  double w, h;                      //width and height of part of small pixel inside current big pixel
  

  //Get empty reconstruction array-object
  recon=newReconArray(status);
  if (EXIT_SUCCESS!=*status) return(recon);

 //in order to get same pixelsize as detector:
 //[naxis1 in maskpixels]*[pixelsize of mask]->absolut size of mask in meters
 //divided by [detector pixelsize]
    recon->naxis1=2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth);
    recon->naxis2=2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth);

  //memory-allocation for Rmap
  recon->Rmap=(double**)malloc(recon->naxis1*sizeof(double*));
  if(NULL!=recon->Rmap){
    for(x=0; x < recon->naxis1; x++){
      recon->Rmap[x]=(double*)malloc(recon->naxis2*sizeof(double));
	if(NULL==recon->Rmap[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "ReconstructionArray!\n", *status);
	  return(recon);
	}
	//Clear the pixels
	for(y=0; y < recon->naxis2; y++){
	  recon->Rmap[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(recon);
  } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "ReconstructionArray!\n", *status);
      return(recon);
  }//end of memory-allocation


   recon->open_fraction = mask->transparency;


  // decide whether pixelsize is different for det and mask
   if(mask->cdelt1 == detector_pixels->xpixelwidth){

  //Scanning over all mask-elements
  for(xcount=0; xcount < mask->naxis1; xcount++){
    for(ycount=0; ycount < mask->naxis2; ycount++){
      if(mask->map[xcount][ycount]==1){
	recon->Rmap[xcount][ycount]=1.;
      }else if(mask->map[xcount][ycount]==0){
	recon->Rmap[xcount][ycount]=
	  (recon->open_fraction)/(recon->open_fraction - 1.);
      }else{
	*status=EXIT_FAILURE;
	HD_ERROR_THROW("Error while scanning mask-elements!\n", *status);
	}
    }
  }//End of scanning over mask-elements
   }//End equal pixelsize

   else{//begin diff pixelsize
     //initialize all to max neg value (depending on OF)
     for(xcount=0; xcount < (recon->naxis1/2); xcount++){
       for(ycount=0; ycount < (recon->naxis2/2); ycount++){
	 recon->Rmap[xcount][ycount]=(recon->open_fraction)/
	   (recon->open_fraction - 1.);
	}
     }

     //distance between max neg value and 1 (which has to be distributed
     // accordingly to new smaller pixels)
     double diff=1.-(recon->open_fraction)/(recon->open_fraction - 1.);

     //Scanning over all mask-elements to get Rmap with same pixel-size as detector
     for(ycount=0; ycount<mask->naxis2;ycount++){
       for(xcount=0; xcount<mask->naxis1;xcount++){

	 topbig=ycount*mask->cdelt2;   //top of current big pixel
	 ypixelcount=topbig/detector_pixels->ypixelwidth;  //count for small pixel (new pixels in Rmap)
	 //current y-pix: top border of big pix/width of one small pix
	 //->determines 1st small in curent big

	 do{//as long as in current big pixel in y-direction
	   topsmall=ypixelcount*detector_pixels->ypixelwidth; //top border of small pix: 
	                                                      //current small pix*width of one
	   if(topsmall<topbig){//1st small in current big starts with part of it in former pix
	     h=detector_pixels->ypixelwidth-(topbig-topsmall);//height that lies in current big:
	                                 //smallwidth-(part of height that lies in former big pix)
	    }else{
	     if((topsmall+detector_pixels->ypixelwidth) <= (topbig+mask->cdelt2)){
              //small pix lies completely in big
	       h=detector_pixels->ypixelwidth;
	     }else{//small pix is at border of curent big->part of it in next big
	       h=topbig+mask->cdelt2-topsmall;//part of height in current big:
	       //top of next big - top of current small
	      }
	    }

	   leftbig=xcount*mask->cdelt1;
	   xpixelcount=leftbig/detector_pixels->xpixelwidth;
	   do{//as long as in current big pixel in x-direction
	     leftsmall=xpixelcount*detector_pixels->xpixelwidth;
	     if(leftsmall<leftbig){
	       w=detector_pixels->xpixelwidth-(leftbig-leftsmall);
	     }else{
	       if((leftsmall+detector_pixels->xpixelwidth) <= (leftbig+mask->cdelt1)){
		 w=detector_pixels->xpixelwidth;
	       }else{
		 w=leftbig+mask->cdelt1-leftsmall;
	       }
	     }
	     if(mask->map[xcount][ycount]==1){//all small transparent pixel-areas
	       // contribute as percentage
	       recon->Rmap[xpixelcount][ypixelcount]+=//one small pix can have 
		         //contributions from parts lying in diff big pix
		 h*w/(detector_pixels->xpixelwidth*detector_pixels->ypixelwidth)*diff;
	       //percentage of area with respect to area of whole small pix
	     }//multiplied with max diff occuring in Rmap values, in order to distribute them accordingly
	     xpixelcount++;
	   }while(leftsmall+detector_pixels->xpixelwidth <= (leftbig+mask->cdelt1));
	   //end current big pixel x-direction


	     ypixelcount++;
	 }while(topsmall+detector_pixels->ypixelwidth <= (topbig+mask->cdelt2));
	 //end current big pixel y-direction
      }
     }  
   }
   return(recon);
}


double* SaveReconArray1d(ReconArray* recon, int* status)
{
 double* ReconArray1d=NULL;

 //Memory-Allocation for 1d-image
 ReconArray1d = (double*)malloc(recon->naxis1*recon->naxis2*sizeof(double));
 if (!ReconArray1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-recon-array!\n", *status);
    return(ReconArray1d);
 }

    //Create the 1D-image from ReconArray 
  int x, y;
  for (x=0; x<recon->naxis1; x++) {
    for (y=0; y<recon->naxis2; y++) {
	ReconArray1d[(x+ recon->naxis1*y)] = recon->Rmap[x][y];
   }
  }
 return(ReconArray1d);  
}


void FreeReconArray1d(double* ReconArray1d)
{
  if (ReconArray1d!=NULL) free(ReconArray1d);
}


void FreeReconArray(ReconArray** const recon)
{
  if ((*recon)!=NULL) {
    if (((*recon)->naxis1>0)&&(NULL!=(*recon)->Rmap)) {
      int count;
      for(count=0; count<(*recon)->naxis1; count++) {
	if (NULL!=(*recon)->Rmap[count]) {
	  free((*recon)->Rmap[count]);
	}
      }
      free((*recon)->Rmap);
    }
    free(*recon);
    *recon=NULL;
  }
}



 
