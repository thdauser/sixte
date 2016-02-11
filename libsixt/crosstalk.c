#include "crosstalk.h"



void init_crosstalk(AdvDet* det, int* status){
	// load the channel list
	if (det->readout_channels==NULL){
	//	init_readout_channels(det,status);
		CHECK_STATUS_VOID(*status);
	}
}


channel_list* init_readout_channels(AdvDet* det, int* status){
	channel_list* chans = load_channel_list(det->channel_file, status);
	CHECK_STATUS_RET_NULL(*status);
	return chans;
}

void free_channel_list(channel_list** chans){
	if (*(chans)!=NULL){
		if ((*chans)->chan != NULL){
			free((*chans)->chan);
		}
		if ((*chans)->pixid != NULL){
			free((*chans)->pixid);
		}
		if ((*chans)->freq != NULL){
			free((*chans)->freq);
		}
		free(*chans);
	}
}

channel_list* load_channel_list(char* fname, int* status){

	channel_list* chans = (channel_list*)malloc(sizeof(channel_list));
	CHECK_MALLOC_RET_NULL(chans);

	// init parameters
	chans->len=0;
	chans->chan=NULL;
	chans->freq=NULL;
	chans->pixid=NULL;

	// open the file
	FILE *file=NULL;
	file=fopen(fname, "r");

	int n=0;

	while(!feof(file)){

	     chans->pixid=realloc(chans->pixid, (n+1)*sizeof(int));
	     chans->chan=realloc(chans->chan, (n+1)*sizeof(int));
	     chans->freq=realloc(chans->freq, (n+1)*sizeof(double));

	     if ((chans->pixid==NULL)||(chans->freq==NULL)||(chans->chan==NULL)){
	    	 	 free_channel_list(&chans);
	    	 	 SIXT_ERROR("error when reallocating arrays");
	    	 	 *status=EXIT_FAILURE;
	    	 	 break;
	     }

	     fscanf(file, "%i %i %f\n",
	    		 &(chans->pixid[n]),
				 &(chans->chan[n]),
				 &(chans->freq[n]));

	     printf("reading channel list (line %i): %i %i %f\n",n+1,chans->pixid[n],chans->chan[n],chans->freq[n]);

	    n++;
	}
	return chans;
}
