#include "crosstalk.h"



void init_crosstalk(AdvDet* det, int* status){
	// load the channel list
	if (det->readout_channels==NULL){
		det->readout_channels = get_readout_channels(det,status);
		CHECK_STATUS_VOID(*status);
	}
}

/**static Channel* newChannel(*ReadoutChannels, status){
	Channel* chan = (Channel*) malloc(sizeof (Channel) );
	CHECK_MALLOC_RET_NULL_STATUS(chan,*status);

	chan->channel_id=-1;
	chan->num_pixels=0;
	chan->pixels=NULL;
	return chan;
}*/

static void add_pixel_to_readout(ReadoutChannels* read_chan, AdvPix* pix, int ic, int* status){

	// check if PixID already belongs to a Channel (duplicated naming not allowed)
	if (pix->channel != NULL){
		printf("*** error: Pixel with ID %i already belongs to a Channel (check channel_file)\n",pix->pindex);
		*status = EXIT_FAILURE;
		return;
	}

	if ( (ic <= 0) || ic > read_chan->num_channels ){
		printf("*** error: Channel with ID %i not a valid channel number\n",ic);
		*status = EXIT_FAILURE;
		return;
	}

	// add another pixel to the channel (ID starts at 1)
	// (note that we use here that realloc behaves like malloc for a NULL pointer)
	read_chan->channels[ic-1].pixels = (AdvPix**)
			realloc( read_chan->channels[ic-1].pixels,
					(read_chan->channels[ic-1].num_pixels+1) * sizeof(AdvPix*)
					);
	CHECK_MALLOC_VOID_STATUS(read_chan->channels[ic-1].pixels,*status);

	read_chan->channels[ic-1].pixels[read_chan->channels[ic-1].num_pixels] = pix;

	read_chan->channels[ic-1].num_pixels++;

	return;
}

static int get_num_chans(channel_list* cl, int* status){

	int num_chans = -1;

	for(int ii=0; ii<cl->len; ii++){
		num_chans = MAX(num_chans,cl->chan[ii]);
	}

	if (num_chans <= 0 ){
		CHECK_STATUS_RET(*status,0);
	}

	return num_chans;
}

static ReadoutChannels* init_newReadout(int num_chan, int* status){

	ReadoutChannels* read_chan = (ReadoutChannels*) malloc(sizeof (ReadoutChannels) );
	CHECK_MALLOC_RET_NULL_STATUS(read_chan,*status);

	read_chan->num_channels = num_chan;

	read_chan->channels = (Channel*) malloc (num_chan*sizeof(Channel));
	CHECK_MALLOC_RET_NULL_STATUS(read_chan->channels,*status);

	for (int ii=0; ii<num_chan; ii++){
		read_chan->channels[ii].channel_id = ii+1; // channel IDs go from 1 to num_chan
		read_chan->channels[ii].num_pixels = 0;
		read_chan->channels[ii].pixels = NULL;
	}

	return read_chan;
}

ReadoutChannels* get_readout_channels(AdvDet* det, int* status){

	// get the complete file name
	char* fullfilename = strdup(det->filepath);
	strncat(fullfilename,det->channel_file,MAXFILENAME);

	// get the channel list
	channel_list* chans = load_channel_list(fullfilename, status);
	CHECK_STATUS_RET(*status,NULL);

	// check if the file agrees with the number of pixels
	if (chans->len != det->npix){
		printf("*** error: number of pixels from channel_file %s (%i) is not equal total number of pixels (%i)!\n",
				fullfilename,chans->len,det->npix);
		*status=EXIT_FAILURE;
		return NULL;
	}


	int num_chan = get_num_chans(chans,status);
	CHECK_STATUS_RET(*status,NULL);
	ReadoutChannels* read_chan = init_newReadout(num_chan, status);

	// sort pixels in the channel array
	AdvPix* pix=NULL;
	for (int ii=0; ii < chans->len; ii++ ){

		// check if PixID makes sense (PixID starts at 0)
		if ( (chans->pixid[ii] < 0) || (chans->pixid[ii] > det->npix-1) ){
			printf("*** error: Pixel-ID %i does not belong to the detector \n",chans->pixid[ii]);
			*status = EXIT_FAILURE;
			return NULL;
		}

		pix = &(det->pix[chans->pixid[ii]]);

		// set frequency of the pixel
		pix->freq = chans->freq[ii];
		if (pix->freq <= 0.0){
			printf("*** warning: assiging unrealistic frequency (%.3e) to pixel %i\n",pix->freq,pix->pindex);
		}

		// assign pixel to readout channel
		add_pixel_to_readout(read_chan, pix, chans->chan[ii], status);
		// make sure the pixel knows its channel (Channel ID starts at 1)
		pix->channel = &(read_chan->channels[chans->chan[ii]-1]);

//		printf("Readout Channel: Assigne PixID %i to Channel %i with Freq %.3e\n",
//				pix->pindex,pix->channel->channel_id,pix->freq);
	}

	return read_chan;
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

	     chans->pixid= realloc(chans->pixid, (n+1)*sizeof(int));
	     chans->chan = realloc(chans->chan,  (n+1)*sizeof(int));
	     chans->freq = realloc(chans->freq,  (n+1)*sizeof(double));

	     if ((chans->pixid==NULL)||(chans->freq==NULL)||(chans->chan==NULL)){
	    	 	 free_channel_list(&chans);
	    	 	 SIXT_ERROR("error when reallocating arrays");
	    	 	 *status=EXIT_FAILURE;
	    	 	 break;
	     }
	     // check that always all three numbers are matched
	     if ( (fscanf(file, "%i %i %lf\n",&(chans->pixid[n]),&(chans->chan[n]),&(chans->freq[n]))) != 3){
	    	 printf("*** formatting error in line %i of the channel file %s: check your input file\n",n+1,fname);
	    	 *status=EXIT_FAILURE;
	    	 break;
	     }
	      // printf("reading channel list (line %i): %i %i %lf\n",n+1,chans->pixid[n],chans->chan[n],chans->freq[n]);
	     n++;
	}
    chans->len = n;

	return chans;
}
