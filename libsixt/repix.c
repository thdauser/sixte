#include "repix.h"

void repixNoReminder(void* arg_dataBeforeRepix, void* arg_dataAfterRepix, int type, 
		     int Size1, int Size2, double pixelwidth_big, double pixelwidth_small)
{
  int xcount, ycount;                 //bigcount
  int xpixelcount=0, ypixelcount=0;   //smallcount
  double leftsmall=0.,leftbig=0.;     //left border of small and big pixel
  double topsmall=0.,topbig=0.;       //top border of small and big pixel

  double** dataBR=NULL; //before repix
  double** dataAR=NULL; //after repix

  if(type==TREADEVENT){
    ReadEvent* ea =(ReadEvent*)arg_dataBeforeRepix;
    dataBR=ea->EventArray;
    ReadEvent* ear =(ReadEvent*)arg_dataAfterRepix;
    dataAR=ear->EventArray;
  }

  //Scanning over all Array-elements to get Array with smaller pixel-size
  for(ycount=0; ycount<Size1; ycount++){
    for(xcount=0; xcount<Size2; xcount++){
      
      topbig=ycount*pixelwidth_big;   //top of current big pixel
      ypixelcount=ceil(topbig/pixelwidth_small);  //count for new smaller pixels 
       //current y-pix: top border of big pix/width of one small pix->determines 1st small in current big

      do{//as long as in current big pixel in y-direction
	topsmall=ypixelcount*pixelwidth_small; //top border of small pix: current small pix*width of one
	   
	leftbig=xcount*pixelwidth_big;
	xpixelcount=ceil(leftbig/pixelwidth_small);
	do{//as long as in current big pixel in x-direction
	  leftsmall=xpixelcount*pixelwidth_small;
    
	  dataAR[xpixelcount][ypixelcount]=dataBR[xcount][ycount];	     

	  xpixelcount++;
	}while(leftsmall+pixelwidth_small < (leftbig+pixelwidth_big));
	//end current big pixel x-direction
	ypixelcount++;
      }while(topsmall+pixelwidth_small < (topbig+pixelwidth_big));
      //end current big pixel y-direction

    }
  }

}

void repixWithReminder(void* arg_dataBeforeRepix, void* arg_dataAfterRepix, int type, 
		       int Size1, int Size2, double pixelwidth_big, double pixelwidth_small,
		       double MinVal)
{
  int xcount, ycount;                //bigcount
  int xpixelcount=0, ypixelcount=0;  //smallcount
  double leftsmall=0.,leftbig=0.;    //left border of small and big pixel
  double topsmall=0.,topbig=0.;      //top border of small and big pixel
  double w, h;                       //width and height of part of small pixel inside current big pixel
  double diff;                       //difference between MinVal (to which all BR-pixels are initialzed to)
                                       // and current not zero pix
  double dataBRbuffer;               //buffer for BR-data -> has to be typecasted for some input-types

  double** dataBR_double=NULL; //before repix  //only one BR-array needed, depending on input type
  int** dataBR_int=NULL; //before repix
  double** dataAR=NULL; //after repix

  if(type==TMASKMAP){
    CodedMask* cm =(CodedMask*)arg_dataBeforeRepix;
    dataBR_int=cm->map;
    ReconArray* ra =(ReconArray*)arg_dataAfterRepix;
    dataAR=ra->Rmap;
  }else if(type==TPROJMASK){ 
    ProjectedMask* pmBR =(ProjectedMask*)arg_dataBeforeRepix;
    dataBR_double=pmBR->map;
    ProjectedMask* pmAR =(ProjectedMask*)arg_dataAfterRepix;
    dataAR=pmAR->map;
  }

  //Scanning over all Array-elements to get Array with smaller pixel-size
  for(ycount=0; ycount<Size1; ycount++){
    for(xcount=0; xcount<Size2; xcount++){
      
      topbig=ycount*pixelwidth_big;   //top of current big pixel
      ypixelcount=topbig/pixelwidth_small;  //count for new smaller pixels 
       //current y-pix: top border of big pix/width of one small pix->determines 1st small in current big

      do{//as long as in current big pixel in y-direction
	topsmall=ypixelcount*pixelwidth_small; //top border of small pix: current small pix*width of one
	  
	if(topsmall<topbig){//1st small in current big starts with part of it in former pix
	     h=pixelwidth_small-(topbig-topsmall);//height that lies in current big:
	                                 //smallwidth-(part of height that lies in former big pix)
	    }else{
	     if((topsmall+pixelwidth_small) <= (topbig+pixelwidth_big)){
              //small pix lies completely in big
	       h=pixelwidth_small;
	     }else{//small pix is at border of curent big->part of it in next big
	       h=topbig+pixelwidth_big-topsmall;//part of height in current big:
	       //top of next big - top of current small
	      }
	    }

	leftbig=xcount*pixelwidth_big;
	xpixelcount=leftbig/pixelwidth_small;
	do{//as long as in current big pixel in x-direction
	  leftsmall=xpixelcount*pixelwidth_small;
  
	  if(leftsmall<leftbig){
	       w=pixelwidth_small-(leftbig-leftsmall);
	     }else{
	       if((leftsmall+pixelwidth_small) <= (leftbig+pixelwidth_big)){
		 w=pixelwidth_small;
	       }else{
		 w=leftbig+pixelwidth_big-leftsmall;
	       }
	     }

	  if(type==TMASKMAP){
	    dataBRbuffer=(double)dataBR_int[xcount][ycount];
	  }else{
	    dataBRbuffer=dataBR_double[xcount][ycount];
	  }


	  if(dataBRbuffer!=0.){//all small transparent pixel-areas contribute as percentage
	    //distance between MinVal and MaxVal has to be distributed to new smaller pix
	    diff=dataBRbuffer-MinVal; //dataBR[xcount][ycount]=MaxVal;
	                                        //MinVal might be neg-> greater diff
	    dataAR[xpixelcount][ypixelcount]+=//one small pix can have 
	      //contributions from parts lying in diff big pix
	      h*w/(pixelwidth_small*pixelwidth_small)*diff;
	      //percentage of area with respect to area of whole small pix
	  }//multiplied with max diff occuring in Rmap values, in order to distribute them accordingly


	  xpixelcount++;
	}while(leftsmall+pixelwidth_small < (leftbig+pixelwidth_big));
	//end current big pixel x-direction
	ypixelcount++;
      }while(topsmall+pixelwidth_small < (topbig+pixelwidth_big));
      //end current big pixel y-direction

    }
  }

}

