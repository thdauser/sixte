#include "crosstalk.h"

/** Calculates distance between two pixels */
static double distance_two_pixels(AdvPix* pix1,AdvPix*pix2){
	// Note: this assumes no pixel overlap (this was ensured at detector loading stage)

	// pix1 is on the right of pix2
	if(pix1->sx>pix2->sx + .5*pix2->width){
		// pix1 is above pix2
		if (pix1->sy > pix2->sy + .5*pix2->height){
			// distance is distance between bottom left corner of pix1 and top right corner of pix2
			return sqrt(pow(pix1->sx - .5*pix1->width - (pix2->sx + .5*pix2->width),2) +
					pow(pix1->sy - .5*pix1->height - (pix2->sy + .5*pix2->height),2));
		}
		// pix1 is below pix2
		if (pix1->sy < pix2->sy - .5*pix2->height){
			// distance is distance between top left corner of pix1 and bottom right corner of pix2
			return sqrt(pow(pix1->sx - .5*pix1->width - (pix2->sx + .5*pix2->width),2) +
					pow(pix1->sy + .5*pix1->height - (pix2->sy - .5*pix2->height),2));
		}
		// pix1 is at the same level as pix2
		// distance is distance between left edge of pix1 and right edge of pix2
		return pix1->sx - .5*pix1->width - (pix2->sx + .5*pix2->width);
	}

	// pix1 is on the left of pix2
	if(pix1->sx<pix2->sx - .5*pix2->width){
		// pix1 is above pix2
		if (pix1->sy > pix2->sy + .5*pix2->height){
			// distance is distance between bottom right corner of pix1 and top left corner of pix2
			return sqrt(pow(pix1->sx + .5*pix1->width - (pix2->sx - .5*pix2->width),2) +
					pow(pix1->sy - .5*pix1->height - (pix2->sy + .5*pix2->height),2));
		}
		// pix1 is below pix2
		if (pix1->sy < pix2->sy - .5*pix2->height){
			// distance is distance between top right corner of pix1 and bottom left corner of pix2
			return sqrt(pow(pix1->sx + .5*pix1->width - (pix2->sx - .5*pix2->width),2) +
					pow(pix1->sy + .5*pix1->height - (pix2->sy - .5*pix2->height),2));
		}
		// pix1 is at the same level as pix2
		// distance is distance between right edge of pix1 and left edge of pix2
		return (pix2->sx - .5*pix2->width) - (pix1->sx + .5*pix1->width);

	}

	// pix1 is at the same left/right position as pix2

	// pix1 is above pix2
	if (pix1->sy > pix2->sy){
		// distance is distance between bottom edge of pix1 and top edge of pix2
		return (pix1->sy - .5*pix1->height - (pix2->sy + .5*pix2->height));
	}
	// pix1 is below pix2
	// distance is distance between top edge of pix1 and bottom edge of pix2
	return (pix2->sy - .5*pix2->height - (pix1->sy + .5*pix1->height));
}

/** Adds a cross talk pixel to the matrix */
static void add_xt_pixel(MatrixCrossTalk** matrix,AdvPix* pixel,double xt_weigth,int* const status){
	CHECK_STATUS_VOID(*status);

	// Allocate matrix if necessary
	if(*matrix==NULL){
		*matrix = newMatrixCrossTalk(status);
		CHECK_MALLOC_VOID_STATUS(*matrix,*status);
	}

	// Increase matrix size
	(*matrix)->cross_talk_pixels = realloc((*matrix)->cross_talk_pixels,((*matrix)->num_cross_talk_pixels+1)*sizeof(*((*matrix)->cross_talk_pixels)));
	CHECK_MALLOC_VOID_STATUS((*matrix)->cross_talk_pixels,*status);
	(*matrix)->cross_talk_weights = realloc((*matrix)->cross_talk_weights,((*matrix)->num_cross_talk_pixels+1)*sizeof(*((*matrix)->cross_talk_weights)));
	CHECK_MALLOC_VOID_STATUS((*matrix)->cross_talk_weights,*status);

	// Affect new values
	(*matrix)->cross_talk_pixels[(*matrix)->num_cross_talk_pixels] = pixel;
	(*matrix)->cross_talk_weights[(*matrix)->num_cross_talk_pixels] = xt_weigth;

	// Now, we can say that the matrix is effectively bigger
	(*matrix)->num_cross_talk_pixels++;
}

// Loads thermal cross-talk for requested pixel
// Concretely, iterates over all the pixels to find neighbours
static void load_thermal_cross_talk(AdvDet* det,int pixid,int* const status){
	double max_cross_talk_dist = det->xt_dist_thermal[det->xt_num_thermal-1];
	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;
	for (int i=0;i<det->npix;i++){
		// Cross-talk is not with yourself ;)
		if (i==pixid){
			continue;
		}
		current_pixel = &(det->pix[i]);

		// Initial quick distance check to avoid spending time on useless pixels
		if((fabs(current_pixel->sx-concerned_pixel->sx) > .5*(current_pixel->width+concerned_pixel->width)+max_cross_talk_dist) ||
				(fabs(current_pixel->sy-concerned_pixel->sy) > .5*(current_pixel->height+concerned_pixel->height)+max_cross_talk_dist)){
			continue;
		}

		// Get distance between two pixels
		double pixel_distance = distance_two_pixels(current_pixel,concerned_pixel);
		if (pixel_distance<0){ // distance should be positive
			*status = EXIT_FAILURE;
			printf("*** error: Distance between pixels %d and %d is negative\n",current_pixel->pindex,concerned_pixel->pindex);
			return;
		}
		//printf("%d - %d : %f\n",pixid,i,pixel_distance*1e6);

		// Iterate over cross-talk values and look for the first matching one
		for (int xt_index=0;xt_index<det->xt_num_thermal;xt_index++){
			if (pixel_distance<det->xt_dist_thermal[xt_index]){
				add_xt_pixel(&concerned_pixel->thermal_cross_talk,current_pixel,det->xt_weight_thermal[xt_index],status);
				CHECK_STATUS_VOID(*status);
				// If one cross talk was identified, go to next pixel (the future cross-talks should be lower order cases)
				break;
			}
		}
	}
}

// Loads electrical cross-talk for requested pixel
// Concretely, iterates over all the pixels of the channel
static void load_electrical_cross_talk(AdvDet* det,int pixid,int* const status){
	CHECK_STATUS_VOID(*status);
	if (det->elec_xt_par==NULL){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load electrical crosstalk with no corresponding information available at detector level");
		return;
	}

	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;
	double weight=0.;


	// Iterate over the channel
	for (int i=0;i<concerned_pixel->channel->num_pixels;i++){
		current_pixel = concerned_pixel->channel->pixels[i];

		// Cross-talk is not with yourself ;)
		if (current_pixel==concerned_pixel) continue;

		// Carrier overlap
		weight = pow(det->elec_xt_par->R0/(4*M_PI*(concerned_pixel->freq - current_pixel->freq)*det->elec_xt_par->Lfprim),2);

		// Common impedence
		weight+=pow(concerned_pixel->freq*det->elec_xt_par->Lcommon/(2*(concerned_pixel->freq - current_pixel->freq)*det->elec_xt_par->Lfsec),2);

		// Add cross-talk pixel
		add_xt_pixel(&concerned_pixel->electrical_cross_talk,current_pixel,weight,status);
		CHECK_STATUS_VOID(*status);

	}

}

/** Load crosstalk timedependence information */
static void load_crosstalk_timedep(AdvDet* det,int* const status){
	if(det->crosstalk_timedep_file==NULL){
		*status=EXIT_FAILURE;
		SIXT_ERROR("Tried to load crosstalk without timedependence information");
		return;
	}

	char fullfilename[MAXFILENAME];
	strcpy(fullfilename,det->filepath);
	strncat(fullfilename,det->crosstalk_timedep_file,MAXFILENAME);

	// open the file
	FILE *file=NULL;
	file=fopen(fullfilename, "r");
	if (file == NULL){
		printf("*** error reading file %s \n",fullfilename);
		*status=EXIT_FAILURE;
		return;
	}

	// initialize crosstalk_timedep structure
	det->crosstalk_timedep = newCrossTalkTimedep(status);
	CHECK_STATUS_VOID(*status);

	double* time=NULL;
	double* weight=NULL;

	// Ignore first line
	char buffer[1000];
	fgets(buffer, 1000, file);

	// Read time dep info
	while(!feof(file)){
		time = realloc(det->crosstalk_timedep->time, (det->crosstalk_timedep->length+1)*sizeof(double));
		weight = realloc(det->crosstalk_timedep->weight, (det->crosstalk_timedep->length+1)*sizeof(double));

		if ((time==NULL)||(weight==NULL)){
			freeCrosstalkTimedep(&det->crosstalk_timedep);
			SIXT_ERROR("error when reallocating arrays in crosstalk timedep loading");
			*status=EXIT_FAILURE;
			return;
		} else {
			det->crosstalk_timedep->time = time;
			det->crosstalk_timedep->weight = weight;
		}
		// check that always all three numbers are matched
		if ( (fscanf(file, "%lg %lg\n",&(det->crosstalk_timedep->time[det->crosstalk_timedep->length]),&(det->crosstalk_timedep->weight[det->crosstalk_timedep->length]))) != 2){
			printf("*** formatting error in line %i of the channel file %s: check your input file\n",det->crosstalk_timedep->length+1,fullfilename);
			*status=EXIT_FAILURE;
			return;
		}
		// printf("reading channel list (line %i): %i %i %lf\n",n+1,chans->pixid[n],chans->chan[n],chans->freq[n]);
		det->crosstalk_timedep->length++;
	}
	fclose(file);

}

void init_crosstalk(AdvDet* det, int* const status){
	// load time dependence
	load_crosstalk_timedep(det,status);
	CHECK_STATUS_VOID(*status);

	// load the channel list
	if (det->readout_channels==NULL){
		det->readout_channels = get_readout_channels(det,status);
		CHECK_STATUS_VOID(*status);
	}

	// load thermal cross talk
	if (det->xt_num_thermal>0){
		for (int i=0;i<det->npix;i++){
			load_thermal_cross_talk(det,i,status);
			CHECK_STATUS_VOID(*status);
		}
	}

	// load electrical cross talk
	if (det->elec_xt_par!=NULL){
		for (int i=0;i<det->npix;i++){
			load_electrical_cross_talk(det,i,status);
			CHECK_STATUS_VOID(*status);
		}
	}
}

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
	char fullfilename[MAXFILENAME];
	strcpy(fullfilename,det->filepath);
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

	// free the channel list
	free_channel_list(&chans);

	return read_chan;
}

void free_channel_list(channel_list** chans){
	if (*(chans)!=NULL){
		free((*chans)->chan);
		free((*chans)->pixid);
		free((*chans)->freq);
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

	if (file == NULL){
		printf("*** error reading file %s \n",fname);
	 	 *status=EXIT_FAILURE;
	 	 return NULL;
	}

	int n=0;
	int* cha;
	int* pix;
	double* freq;
	while(!feof(file)){

	     pix = realloc(chans->pixid, (n+1)*sizeof(int));
	     cha = realloc(chans->chan,  (n+1)*sizeof(int));
	     freq = realloc(chans->freq,  (n+1)*sizeof(double));

	     if ((cha==NULL)||(pix==NULL)||(freq==NULL)){
	    	 	 free_channel_list(&chans);
	    	 	 SIXT_ERROR("error when reallocating arrays");
	    	 	 *status=EXIT_FAILURE;
	    	 	 return NULL;
	     } else {
	    	 chans->pixid = pix;
	    	 chans->chan  = cha;
	    	 chans->freq  = freq;
	     }
	     // check that always all three numbers are matched
	     if ( (fscanf(file, "%i %i %lf\n",&(chans->pixid[n]),&(chans->chan[n]),&(chans->freq[n]))) != 3){
	    	 printf("*** formatting error in line %i of the channel file %s: check your input file\n",n+1,fname);
	    	 *status=EXIT_FAILURE;
	    	 return NULL;
	     }
	      // printf("reading channel list (line %i): %i %i %lf\n",n+1,chans->pixid[n],chans->chan[n],chans->freq[n]);
	     n++;
	}
    chans->len = n;

	return chans;
}


void freeCrosstalk(AdvDet** det){

	freeReadoutChannels( (*det)->readout_channels);

}
