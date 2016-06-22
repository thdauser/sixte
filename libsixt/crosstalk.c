/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2016 Philippe Peille, IRAP; Thomas Dauser, ECAP
*/

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

/** Adds a cross talk pixel to the matrix */
static void add_xt_enerdep_pixel(MatrixEnerdepCrossTalk** matrix,AdvPix* pixel,double* xt_weigth,int n_ener,int* const status){
	CHECK_STATUS_VOID(*status);

	// Allocate matrix if necessary
	if(*matrix==NULL){
		*matrix = newMatrixEnerdepCrossTalk(status);
		(*matrix)->n_ener = n_ener;
		CHECK_MALLOC_VOID_STATUS(*matrix,*status);
	}

	if ((*matrix)->n_ener != n_ener){
		printf(" *** error : number of energy bins is %i, but should be %i \n",
				(*matrix)->n_ener,n_ener);
		SIXT_ERROR(" adding energy-dependent crosstalk weight failed ");
		*status=EXIT_FAILURE;
		return;
	}

	// Increase matrix size
	(*matrix)->cross_talk_pixels = realloc((*matrix)->cross_talk_pixels,
			((*matrix)->num_cross_talk_pixels+1)*sizeof(*((*matrix)->cross_talk_pixels)));
	CHECK_MALLOC_VOID_STATUS((*matrix)->cross_talk_pixels,*status);
	(*matrix)->cross_talk_weights = realloc((*matrix)->cross_talk_weights,
			((*matrix)->num_cross_talk_pixels+1)*sizeof(*((*matrix)->cross_talk_weights)));
	CHECK_MALLOC_VOID_STATUS((*matrix)->cross_talk_weights,*status);

	// allocate space for the array (n_ener bins)
	(*matrix)->cross_talk_weights[(*matrix)->num_cross_talk_pixels] =
			(double*) malloc( n_ener *sizeof(double) );

	// Affect new values
	(*matrix)->cross_talk_pixels[(*matrix)->num_cross_talk_pixels] = pixel;
	for (int ii=0; ii<n_ener; ii++){
		(*matrix)->cross_talk_weights[(*matrix)->num_cross_talk_pixels][ii] = xt_weigth[ii];
	}

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

static void initElecTab(ElecTab** tab, int n_freq_s, int n_freq_p, int n_ener_p,
		double* freq_s, double* freq_p, double* ener_p, int* status){

	(*tab) = (ElecTab*) malloc (sizeof (ElecTab));
	CHECK_MALLOC_VOID_STATUS(tab,*status);

	// make a short pointer
	ElecTab* t = (*tab);

	t->n_freq_s = n_freq_s;
	t->n_freq_p  = n_freq_p;
	t->n_ener_p = n_ener_p;

	t->freq_s = (double*) malloc(n_freq_s * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->freq_s,*status);
	for (int ii=0; ii<n_freq_s; ii++){
		t->freq_s[ii] = freq_s[ii];
	}

	t->freq_p = (double*) malloc(n_freq_p * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->freq_p,*status);
	for (int ii=0; ii<n_freq_p; ii++){
		t->freq_p[ii] = freq_p[ii];
	}

	t->ener_p = (double*) malloc(n_ener_p * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->ener_p,*status);
	for (int ii=0; ii<n_ener_p; ii++){
		t->ener_p[ii] = ener_p[ii];
	}


	// allocate the 3d matrix (n_freq_s x n_freq_p x n_ener_p)
	t->matrix = (double***) malloc (n_freq_s*sizeof(double**));
	CHECK_MALLOC_VOID_STATUS(t->matrix,*status);

	for (int ll=0; ll<n_freq_s; ll++){               // FREQ_S-LOOP
		t->matrix[ll] = (double**) malloc (n_freq_p*sizeof(double*));
		CHECK_MALLOC_VOID_STATUS(t->matrix[ll],*status);

		for (int ii=0; ii<n_freq_p; ii++){             // FREQ_P-LOOP
			t->matrix[ll][ii] = (double*) malloc (n_ener_p*sizeof(double));
			CHECK_MALLOC_VOID_STATUS(t->matrix[ll][ii],*status);

			for (int jj=0; jj<n_ener_p; jj++){      // ENER_P-LOOP
				t->matrix[ll][ii][jj] = 0.0;
			}
		}
	}
	return;
}

static void get_imodtable_axis(int* nrows, double** val, char* extname, char* colname, fitsfile* fptr, int* status){

	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return;
	}

	// get the column id-number
	int colnum;
	if(fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status)) return;

	// get the number of rows
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return;

	// allocate memory for the array
	*val=(double*)malloc(n*sizeof(double));
	CHECK_MALLOC_VOID_STATUS(*val,*status);

    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) n;
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);

	(*nrows) = (int) n;

	return;
}

/** read one electrical crosstalk matrix from a fits file and store it in the structure (units given in keV) */
static void read_elec_matrix(fitsfile* fptr,int n_freq_s, int n_freq_p, int n_ener_p,
		ElecTab* tab, char* extname, int* status){

	// ampl and dt arrays should be initialized before
	assert(tab->freq_s!=NULL);
	assert(tab->freq_p!=NULL);
	assert(tab->ener_p!=NULL);

	// move to the extension containing the data
	int extver = 0;
	if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, extver ,status)){
		printf(" error: moving to extension %s failed \n",extname);
		return;
	}

	// get number of rows
	int naxis;
	if (fits_get_img_dim(fptr, &naxis, status) && (naxis!=4)){
		printf(" error: getting dimensions of intermodulation data array failed \n");
		return;
	}
	long naxes[naxis];


	if (fits_get_img_size(fptr, naxis, naxes, status) ){
		printf(" error: getting dimensions of intermodulation data array failed \n");
		return;
	}
	// check dimensions
	if (naxes[0]!=n_freq_s || naxes[1]!=n_freq_p || naxes[2]!=n_ener_p ){
		*status=EXIT_FAILURE;
		printf(" error: wrong dimensions of the intermodulation data array [%ld %ld %ld] \n",
				naxes[0],naxes[1],naxes[2]);
	}



	// actually read the table now
	int anynul=0;
	double* buf;
	buf = (double*)malloc(n_freq_s*n_freq_p*n_ener_p*sizeof(double));
	CHECK_MALLOC_VOID(buf);

	double nullval=0.0;
	long nelem = n_freq_s*n_freq_p*n_ener_p;
	long fpixel[naxis];
	for (int ii=0;ii<naxis;ii++){
		fpixel[ii]=1;
	}

	fits_read_pix(fptr, TDOUBLE, fpixel,nelem, &nullval,buf, &anynul, status);
	fits_report_error(stdout, *status);
	CHECK_STATUS_VOID(*status);

	headas_chat(7," start reading the electrical crosstalk matrix %s [%ld %ld %ld]",extname,n_freq_s,n_freq_p, n_ener_p);
	for (int ii=0; ii<n_freq_s ; ii++){  // freq_s loop
		for (int jj=0; jj<n_freq_p ; jj++){
			for (int kk=0; kk<n_ener_p ; kk++){
					int ind = 	( kk * n_freq_p + jj ) * n_freq_s + ii;
					assert(ind < nelem);
					tab->matrix[ii][jj][kk] = buf[ind]*1e-3;  // convert to keV !!!
			}
		}
	} // ------------------------  //  end (freq_s loop)
	headas_chat(7," ... done \n");

	free(buf);

	return;
}

/** load the electrical crosstalk tables */
static void load_elec_table(AdvDet* det, int* status){

	// check if the table exists
	CHECK_NULL_VOID(det->crosstalk_elec_file,*status,"no file for the electrical crosstalk table specified");

	char* EXTNAME_FREQ_SIGNAL = "signal_frequency";
	char* EXTNAME_FREQ_PERTURBER = "perturber_frequency";
	char* EXTNAME_ENER_PERTURBER = "perturber_energy";
	char* COLNAME_FREQ_SIGNAL = "FREQ_S";
	char* COLNAME_FREQ_PERTURBER = "FREQ_P";
	char* COLNAME_ENER_PERTURBER = "EN_P";

	// char* EXTNAME_CROSSTALK_TOTAL = "total_crosstalk";
	char* EXTNAME_CARRIER_OVERLAP = "carrier_overlap";
	char* EXTNAME_COMMON_IMPEDANCE = "common_impedance";

	fitsfile *fptr=NULL;

	do {

		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->crosstalk_elec_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the electrical crosstalk table %s \n",fullfilename);

		// read the extensions specifying the axes of the 3d matrix
		int n_freq_s, n_freq_p, n_ener_p;
		double* freq_s;
		double* freq_p;
		double* ener_p;
		get_imodtable_axis(&n_freq_s,&freq_s,EXTNAME_FREQ_SIGNAL,COLNAME_FREQ_SIGNAL,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_freq_p,&freq_p,EXTNAME_FREQ_PERTURBER,COLNAME_FREQ_PERTURBER,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_ener_p,&ener_p,EXTNAME_ENER_PERTURBER,COLNAME_ENER_PERTURBER,fptr,status);
		CHECK_STATUS_BREAK(*status);
		// convert from eV to keV
		for (int ii=0; ii<n_ener_p; ii++){
			ener_p[ii] *=1e-3;
		}

		// initialize the elec crosstalk tables
		initElecTab(&(det->crosstalk_elec_carrier_olap), n_freq_s, n_freq_p, n_ener_p, freq_s, freq_p, ener_p, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing electrical crosstalk table in memory failed");
			break;
		}

		initElecTab(&(det->crosstalk_elec_common_imp), n_freq_s, n_freq_p, n_ener_p, freq_s, freq_p, ener_p, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing electrical crosstalk table in memory failed");
			break;
		}


		read_elec_matrix(fptr,n_freq_s,n_freq_p, n_ener_p, det->crosstalk_elec_carrier_olap,
				EXTNAME_CARRIER_OVERLAP,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading electrical crosstalk table (carrier overlap) %s  failed\n", fullfilename);
			break;
		}

		read_elec_matrix(fptr,n_freq_s,n_freq_p, n_ener_p, det->crosstalk_elec_common_imp,
				EXTNAME_COMMON_IMPEDANCE,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading electrical crosstalk table (common impedance) %s  failed\n", fullfilename);
			break;
		}


//		print_imod_matrix(det->crosstalk_intermod_table);

	} while(0); // END of Error handling loop

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

static void print_elec_matrix(ElecTab* tab_cl,ElecTab* tab_ci){

	FILE * fp;
	fp = fopen ("elec_matrix_output_tab.txt", "w+");

	for (int kk=0; kk<tab_cl->n_ener_p ; kk++){
		for (int jj=0; jj<tab_cl->n_freq_p ; jj++){
				for (int ii=0; ii<tab_cl->n_freq_s ; ii++){
					fprintf(fp, "%i \t %i \t %i  \t %e \t %e \t %e \t %e \t %e \n",
							ii,jj,kk,tab_cl->freq_s[ii],tab_cl->freq_p[jj],tab_cl->ener_p[kk],
							tab_cl->matrix[ii][jj][kk],tab_ci->matrix[ii][jj][kk]);

			}
		}
	} // ------------------------  //  end (freq loop)

	fclose(fp);
}

// get the energy dependent weight (weight needs to be allocated before, values are added
// to the already given values in weight)
static void get_enerdep_weight(ElecTab* tab, double freq_s, double freq_p, double* weight, int n, int* status){

	assert(tab->n_ener_p==n);

	if ((freq_s > tab->freq_s[tab->n_freq_s-1] ) || (freq_s < tab->freq_s[0] )){
		printf(" *** error: signal pixel frequency %g outside tabulated range [%g,%g] \n",
				freq_s,tab->freq_s[0],tab->freq_s[tab->n_freq_s-1]);
		SIXT_ERROR("failed interpolating electrical crosstalk table");
		*status=EXIT_FAILURE;
		return;
	}
	if ((freq_p > tab->freq_p[tab->n_freq_p-1] ) || (freq_p < tab->freq_p[0] )){
		printf(" *** error: perturber pixel frequency %g outside tabulated range [%g,%g] \n",
				freq_p,tab->freq_p[0],tab->freq_p[tab->n_freq_p-1]);
		SIXT_ERROR("failed interpolating electrical crosstalk table");
		*status=EXIT_FAILURE;
		return;
	}

	// interpolate the arrays for the given frequencies
	int ind_freq_s = binary_search(freq_s,tab->freq_s,tab->n_freq_s);
	assert ( (ind_freq_s) >= 0 && (ind_freq_s < tab->n_freq_s-1) ); // should be valid after the above checks
	double d_freq_s =  (freq_s - tab->freq_s[ind_freq_s]) / ( tab->freq_s[ind_freq_s+1] - tab->freq_s[ind_freq_s] );

	// interpolate the arrays for the given frequencies
	int ind_freq_p = binary_search(freq_p,tab->freq_p,tab->n_freq_p);
	assert ( (ind_freq_p) >= 0 && (ind_freq_p < tab->n_freq_p-1) );  // should be valid after the above checks
	double d_freq_p =  (freq_p - tab->freq_p[ind_freq_p]) / ( tab->freq_p[ind_freq_p+1] - tab->freq_p[ind_freq_p] );

	// interpolate in 2d the weights
	for (int ii=0;ii<tab->n_ener_p;ii++){
		weight[ii] +=
				tab->matrix[ind_freq_s][ind_freq_p][ii]*(1-d_freq_s)*(1-d_freq_p) +
				tab->matrix[ind_freq_s+1][ind_freq_p][ii]*(d_freq_s)*(1-d_freq_p) +
				tab->matrix[ind_freq_s][ind_freq_p+1][ii]*(1-d_freq_s)*(d_freq_p) +
				tab->matrix[ind_freq_s+1][ind_freq_p+1][ii]*(d_freq_s)*(d_freq_p);

	}

}

// Loads electrical cross-talk for requested pixel
// Concretely, iterates over all the pixels of the channel
static void load_electrical_cross_talk(AdvDet* det,int pixid,int* const status){
	CHECK_STATUS_VOID(*status);
	if (det->crosstalk_elec_file==NULL){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Tried to load electrical crosstalk with no corresponding file given at detector level");
		return;
	}

	// load the tables for the electrical crosstalk
	if ((det->crosstalk_elec_carrier_olap == NULL) || (det->crosstalk_elec_common_imp == NULL) ){
		load_elec_table(det, status);
		CHECK_STATUS_VOID(*status);
		assert(det->crosstalk_elec_carrier_olap != NULL);
		assert(det->crosstalk_elec_common_imp != NULL);
	//	print_elec_matrix(det->crosstalk_elec_carrier_olap,det->crosstalk_elec_common_imp);
	}


	AdvPix* concerned_pixel = &(det->pix[pixid]);
	AdvPix* current_pixel = NULL;

	// as the electrical crosstalk depends on the energy of the perturber, we need an array
	// of weights depending on the energy
	assert( det->crosstalk_elec_common_imp->n_ener_p == det->crosstalk_elec_carrier_olap->n_ener_p );
	int n_weight=det->crosstalk_elec_common_imp->n_ener_p;
	double weight[n_weight];


	// Iterate over the channel
	for (int i=0;i<concerned_pixel->channel->num_pixels;i++){
		current_pixel = concerned_pixel->channel->pixels[i];

		// Cross-talk is not with yourself ;)
		if (current_pixel==concerned_pixel) continue;

		double freq_s = concerned_pixel->freq;
		double freq_p = current_pixel->freq;

		// set weights to zero;
		for (int jj=0;jj<n_weight;jj++){
			weight[jj] = 0.0;
		}

		// Add Carrier overlap
		get_enerdep_weight(det->crosstalk_elec_carrier_olap,freq_s,freq_p,weight,n_weight,status);

		// Add Common impedence
		get_enerdep_weight(det->crosstalk_elec_common_imp,freq_s,freq_p,weight,n_weight,status);

		// Add cross-talk pixel
		add_xt_enerdep_pixel(&concerned_pixel->electrical_cross_talk,current_pixel,weight,n_weight,status);
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
	strcat(fullfilename,det->crosstalk_timedep_file);

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
	if(fgets(buffer, 1000, file)==NULL){
		printf(" *** error: reading of the time dependence crosstalk file %s failed\n",fullfilename);
		*status=EXIT_FAILURE;
		return;
	}

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

static void initImodTab(ImodTab** tab, int n_ampl, int n_dt, int n_freq,
		double* ampl, double* dt, double* freq, int* status){

	(*tab) = (ImodTab*) malloc (sizeof (ImodTab));
	CHECK_MALLOC_VOID_STATUS(tab,*status);

	// make a short pointer
	ImodTab* t = (*tab);

	t->n_ampl = n_ampl;
	t->n_dt   = n_dt;
	t->n_freq = n_freq;

	t->ampl = (double*) malloc(n_ampl * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->ampl,*status);
	for (int ii=0; ii<n_ampl; ii++){
		t->ampl[ii] = ampl[ii];
	}


	t->dt = (double*) malloc(n_dt * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->dt,*status);

	for (int ii=0; ii<n_dt; ii++){
		t->dt[ii] = dt[ii];
	}
	t->dt_min=dt[0];
	t->dt_max=dt[n_dt-1];

	t->freq = (double*) malloc(n_freq * sizeof(double));
	CHECK_MALLOC_VOID_STATUS(t->freq,*status);

	for (int ii=0; ii<n_freq; ii++){
		t->freq[ii] = freq[ii];
	}


	// allocate the 4d matrix (n_freq x n_dt x n_ampl x nampl)
	t->matrix = (double****) malloc (n_freq*sizeof(double***));
	CHECK_MALLOC_VOID_STATUS(t->matrix,*status);

	for (int ll=0; ll<n_freq; ll++){               // FREQ-LOOP
		t->matrix[ll] = (double***) malloc (n_dt*sizeof(double**));
		CHECK_MALLOC_VOID_STATUS(t->matrix[ll],*status);

		for (int ii=0; ii<n_dt; ii++){             // DT-LOOP
			t->matrix[ll][ii] = (double**) malloc (n_ampl*sizeof(double*));
			CHECK_MALLOC_VOID_STATUS(t->matrix[ll][ii],*status);

			for (int jj=0; jj<n_ampl; jj++){      // AMPL1-LOOP
				t->matrix[ll][ii][jj] = (double*) malloc (n_ampl*sizeof(double));
				CHECK_MALLOC_VOID_STATUS(t->matrix[ll][ii][jj],*status);

				for (int kk=0; kk<n_ampl; kk++){  // AMPL2-LOOP
					t->matrix[ll][ii][jj][kk] = 0.0;
				}
			}
		}
	}
	return;
}

/** read one intermodulation matrix from a fits file and store it in the structure */
static void read_intermod_matrix(fitsfile* fptr,int n_ampl, int n_dt, int n_freq,
		ImodTab* tab, char* extname, int* status){

	// ampl and dt arrays should be initialized before
	assert(tab->ampl!=NULL);
	assert(tab->dt!=NULL);
	assert(tab->freq!=NULL);

	// move to the extension containing the data
	int extver = 0;
//	if (fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status)){
	if (fits_movnam_hdu(fptr, IMAGE_HDU, extname, extver ,status)){
		printf(" error: moving to extension %s failed \n",extname);
		return;
	}

	// get number of rows
//	int colnum;
//	if(fits_get_colnum(fptr, CASEINSEN,colname, &colnum, status)) return;
	int naxis;
	if (fits_get_img_dim(fptr, &naxis, status) && (naxis!=4)){
		printf(" error: getting dimensions of intermodulation data array failed \n");
		return;
	}
	long naxes[naxis];


	if (fits_get_img_size(fptr, naxis, naxes, status) ){
		printf(" error: getting dimensions of intermodulation data array failed \n");
		return;
	}
	// check dimensions
	if (naxes[0]!=n_freq || naxes[1]!=n_dt || naxes[2]!=n_ampl || naxes[3]!=n_ampl ){
		*status=EXIT_FAILURE;
		printf(" error: wrong dimensions of the intermodulation data array [%ld %ld %ld %ld] \n",
				naxes[0],naxes[1],naxes[2],naxes[3]);
	}



	// actually read the table now
	int anynul=0;
	double* buf;
	buf = (double*)malloc(n_freq*n_dt*n_ampl*n_ampl*sizeof(double));
	CHECK_MALLOC_VOID(buf);

	double nullval=0.0;
	long nelem = n_freq*n_dt*n_ampl*n_ampl;
	long fpixel[naxis];
	for (int ii=0;ii<naxis;ii++){
		fpixel[ii]=1;
	}

	fits_read_pix(fptr, TDOUBLE, fpixel,nelem, &nullval,buf, &anynul, status);
	fits_report_error(stdout, *status);
	CHECK_STATUS_VOID(*status);

	headas_chat(7," start reading the intermodulation crosstalk matrix [%i %i %i %i]",n_freq,n_dt,n_ampl,n_ampl);
	for (int ii=0; ii<n_freq ; ii++){  // freq loop
		for (int jj=0; jj<n_dt ; jj++){
			for (int kk=0; kk<n_ampl ; kk++){
				for (int ll=0; ll<n_ampl ; ll++){
					int ind = 	(( ll*n_ampl + kk ) * n_dt + jj ) * n_freq + ii;
					assert(ind < n_dt*n_ampl*n_ampl*n_freq);
					// amplitude is defined descending -> need to reverse it in output array
					tab->matrix[ii][jj][n_ampl-kk-1][n_ampl-ll-1] = buf[ind];
				}
			}
		}
	} // ------------------------  //  end (freq loop)
	headas_chat(7," ... done \n");

	free(buf);

	return;
}

static void print_imod_matrix(ImodTab* tab){

	FILE * fp;

	fp = fopen ("imod_matrix_output_tab.txt", "w+");

	for (int ll=0; ll<tab->n_ampl ; ll++){
		for (int kk=0; kk<tab->n_ampl ; kk++){
			for (int jj=0; jj<tab->n_dt ; jj++){
				for (int ii=0; ii<tab->n_freq ; ii++){  // freq loop
					fprintf(fp, "%i \t %i \t %i \t %i \t %e \t %e \t %e \t %e \t %e \n",
							ii,jj,kk,ll,tab->freq[ii],tab->dt[jj],tab->ampl[tab->n_ampl-kk-1],
							tab->ampl[tab->n_ampl-ll-1],
							tab->matrix[ii][jj][tab->n_ampl-kk-1][tab->n_ampl-ll-1]);

				}
			}
		}
	} // ------------------------  //  end (freq loop)

	fclose(fp);

}

static void reverse_array_dbl(double* arr, int n){
	double t;
	for (int ii = 0; ii < n/2; ii++) {
		t           = arr[ii];
		arr[ii]     = arr[n-1-ii];
		arr[n-1-ii] = t;
	}
}

/** load the intermodulation frequency table */
static void load_intermod_freq_table(AdvDet* det, int* status){

	// check if the table exists
	CHECK_NULL_VOID(det->crosstalk_intermod_file,*status,"no file for the intermodulation table specified");

	char* EXTNAME_AMPLITUDE = "amplitude";
	char* EXTNAME_TIMEOFFSET = "time_offset";
	char* EXTNAME_FREQOFFSET = "frequency_offset";
	char* COLNAME_AMPLITUDE = "AMPLITUDE";
	char* COLNAME_TIMEOFFSET = "DT_SEC";
	char* COLNAME_FREQOFFSET = "D_FREQ";
	char* EXTNAME_CROSSTALK = "crosstalk";

	// we need the time dependence table here
	if (det->crosstalk_timedep==NULL){
		load_crosstalk_timedep(det,status);
		CHECK_STATUS_VOID(*status);
	}

	fitsfile *fptr=NULL;

	do {

		char fullfilename[MAXFILENAME];
		strcpy(fullfilename,det->filepath);
		strcat(fullfilename,det->crosstalk_intermod_file);

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) break;
		headas_chat(5, "   ... reading the intermodulation table %s\n",fullfilename);

		// read the extensions specifying the axes of the 3d matrix
		int n_ampl,n_dt,n_freq;
		double* ampl;
		double* dt;
		double* freq;
		get_imodtable_axis(&n_ampl,&ampl,EXTNAME_AMPLITUDE,COLNAME_AMPLITUDE,fptr,status);
		CHECK_STATUS_BREAK(*status);
		// amplitude is defined descending -> need to reverse it
		reverse_array_dbl(ampl,n_ampl);
		assert(ampl[0]<ampl[1]);

		get_imodtable_axis(&n_dt,&dt,EXTNAME_TIMEOFFSET,COLNAME_TIMEOFFSET,fptr,status);
		CHECK_STATUS_BREAK(*status);
		get_imodtable_axis(&n_freq,&freq,EXTNAME_FREQOFFSET,COLNAME_FREQOFFSET,fptr,status);
		CHECK_STATUS_BREAK(*status);

		// initialize the intermodulation table
		initImodTab(&(det->crosstalk_intermod_table), n_ampl, n_dt, n_freq, ampl, dt, freq, status);
		if (*status!=EXIT_SUCCESS){
			SIXT_ERROR("initializing intermodulation table in memory failed");
			break;
		}

		// set the minimal and maximal dt value (dt_min<0 and dt_max>0 required in the current implementation)
		assert(det->crosstalk_intermod_table->dt_min<0.0 && det->crosstalk_intermod_table->dt_max>0.0);
		if (det->crosstalk_intermod_table->dt_min < det->crosstalk_timedep->time[0]){
			det->crosstalk_intermod_table->dt_min = det->crosstalk_timedep->time[0];
		}
		if (det->crosstalk_intermod_table->dt_max > det->crosstalk_timedep->time[det->crosstalk_timedep->length-1]){
			det->crosstalk_intermod_table->dt_max = det->crosstalk_timedep->time[det->crosstalk_timedep->length-1];
		}

		read_intermod_matrix(fptr,n_ampl,n_dt, n_freq, det->crosstalk_intermod_table,
				EXTNAME_CROSSTALK,status);
		if (*status != EXIT_SUCCESS){
			printf(" *** error: reading intermodulation table %s  failed\n", fullfilename);
			break;
		}

//		print_imod_matrix(det->crosstalk_intermod_table);

	} while(0); // END of Error handling loop

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

/** Load intermodulation cross talk information into a single AdvPix*/
static void load_intermod_cross_talk(AdvDet* det, int* status){
	/** make sure the information is loaded in the table */
	if (det->crosstalk_intermod_table==NULL){
		load_intermod_freq_table(det,status);
		CHECK_STATUS_VOID(*status);
	}
}

static int doCrosstalk(int id, AdvDet* det){
	if ( (id == det->crosstalk_id) || det->crosstalk_id == CROSSTALK_ID_ALL){
		return 1;
	} else {
		return 0;
	}
}

void init_crosstalk(AdvDet* det, int* const status){

	headas_chat(5,"\n[crosstalk] the modes are switched on (%i): \n",det->crosstalk_id);

	// load time dependence
	load_crosstalk_timedep(det,status);
	CHECK_STATUS_VOID(*status);

	// load the channel list
	if (det->readout_channels==NULL){
		det->readout_channels = get_readout_channels(det,status);
		CHECK_STATUS_VOID(*status);
	}


	// load thermal cross talk
	if (det->xt_num_thermal>0 && doCrosstalk(CROSSTALK_ID_THERM,det)){
		headas_chat(5," - thermal crosstalk\n");
		for (int ii=0;ii<det->npix;ii++){
			load_thermal_cross_talk(det,ii,status);
			CHECK_STATUS_VOID(*status);
		}
	}

	// load electrical cross talk
	if (doCrosstalk(CROSSTALK_ID_ELEC,det)){
		headas_chat(5," - electrical crosstalk\n");
		for (int ii=0;ii<det->npix;ii++){
			load_electrical_cross_talk(det,ii,status);
			CHECK_STATUS_VOID(*status);
		}
	}

	// load intermodulation cross talk
	if (det->crosstalk_intermod_file!=NULL && doCrosstalk(CROSSTALK_ID_IMOD,det)){
		headas_chat(5," - intermodulation crosstalk\n");

		// loop through all pixels
		load_intermod_cross_talk(det,status);
		CHECK_STATUS_VOID(*status);
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

	read_chan->df_information_band = (double*) malloc (num_chan*sizeof(double));
	CHECK_MALLOC_RET_NULL_STATUS(read_chan->df_information_band,*status);

	for (int ii=0; ii<num_chan; ii++){
		read_chan->channels[ii].channel_id = ii+1; // channel IDs go from 1 to num_chan
		read_chan->channels[ii].num_pixels = 0;
		read_chan->channels[ii].pixels = NULL;
		read_chan->df_information_band[ii] = 0.0;
	}

	return read_chan;
}

ReadoutChannels* get_readout_channels(AdvDet* det, int* status){

	// get the complete file name
	char fullfilename[MAXFILENAME];
	strcpy(fullfilename,det->filepath);
	strcat(fullfilename,det->channel_file);

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

	// set the information band for each channel now
//	set_df_information_band(read_chan,status);

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
	fclose(file);

    chans->len = n;

	return chans;
}

/** Compute influence of the crosstalk event on an impact using the timedependence table */
int computeCrosstalkInfluence(AdvDet* det,PixImpact* impact,PixImpact* crosstalk,double* influence){
	assert(impact->pixID == crosstalk->pixID);
	// For the moment, the crosstalk event is supposed to happen afterwards, but this could be modified if needed
	double time_difference = crosstalk->time - impact->time;
	double energy_influence = 0.;
	// if impact is close enough to have an influence
	if ((time_difference>det->crosstalk_timedep->time[0]) && (time_difference<det->crosstalk_timedep->time[det->crosstalk_timedep->length-1])){

		headas_chat(7,"TimeI:%.1e TimeC:%.2e CrosstalkE:%.2e ",impact->time,crosstalk->time,crosstalk->energy);
		// Binary search for to find interpolation interval
		int high=det->crosstalk_timedep->length-1;
		int low=0;
		int mid;
		while (high > low) {
			mid=(low+high)/2;
			if (det->crosstalk_timedep->time[mid] < time_difference) {
				low=mid+1;
			} else {
				high=mid;
			}
		}
		int timedep_index = low-1;
		assert(timedep_index<det->crosstalk_timedep->length-1);

		// influence previous impact
		energy_influence= crosstalk->energy*(det->crosstalk_timedep->weight[timedep_index]+
				(det->crosstalk_timedep->weight[timedep_index+1]-det->crosstalk_timedep->weight[timedep_index])/(det->crosstalk_timedep->time[timedep_index+1]-det->crosstalk_timedep->time[timedep_index])*(time_difference-det->crosstalk_timedep->time[timedep_index]));
		impact->energy+=energy_influence;
		*influence+=energy_influence;
		headas_chat(7,", Influence:%.2e (fraction:%.2e)\n",*influence,*influence/crosstalk->energy);
		return 1;
	}
	return 0;
}

/** Cosntructor of CrosstalkProxy structure */
CrosstalkProxy* newCrosstalkProxy(int* const status){
	CrosstalkProxy* xtalk_proxy= (CrosstalkProxy*) malloc (sizeof (*xtalk_proxy));
	CHECK_MALLOC_RET_NULL_STATUS(xtalk_proxy,*status);

	xtalk_proxy->xtalk_impacts = (PixImpact**) malloc(INITXTALKNB*sizeof(PixImpact*));
	CHECK_MALLOC_RET_NULL_STATUS(xtalk_proxy->xtalk_impacts,*status);
	for (int ii=0;ii<INITXTALKNB;ii++){
		xtalk_proxy->xtalk_impacts[ii] = (PixImpact*) malloc(sizeof(PixImpact));
		CHECK_MALLOC_RET_NULL_STATUS(xtalk_proxy->xtalk_impacts[ii],*status);
	}
	xtalk_proxy->xtalk_proxy_size=INITXTALKNB;
	xtalk_proxy->n_active_crosstalk=0;
	xtalk_proxy->current_crosstalk_index=0;
	xtalk_proxy->n_xtalk_pileup=0;

	return xtalk_proxy;
}

/** Destructor of CrosstalkProxy structure */
void freeCrosstalkProxy(CrosstalkProxy** xtalk_proxy){
	if (*(xtalk_proxy)!=NULL){
		if ((*xtalk_proxy)->xtalk_impacts !=NULL){
			for (int ii=0;ii<((*xtalk_proxy)->xtalk_proxy_size);ii++){
				free((*xtalk_proxy)->xtalk_impacts[ii]);
			}
			free((*xtalk_proxy)->xtalk_impacts);
			free(*xtalk_proxy);
		}
	}
}

/** Add crosstalk to proxy */
void addCrosstalk2Proxy(CrosstalkProxy* xtalk_proxy,PixImpact* impact,int* const status){
	// Check that there is room left for new crosstalk
	// If not, allocate new array of double size and copy old impacts (manual copy to reorder circular buffer, should happen rarely anyway)
	if (xtalk_proxy->n_active_crosstalk==xtalk_proxy->xtalk_proxy_size){
		PixImpact** piximp_list = (PixImpact**) malloc((xtalk_proxy->xtalk_proxy_size*2)*sizeof(PixImpact*));
		CHECK_MALLOC_VOID_STATUS(piximp_list,*status);
		for (int ii=0;ii<xtalk_proxy->xtalk_proxy_size;ii++){
			piximp_list[ii] = xtalk_proxy->xtalk_impacts[(xtalk_proxy->current_crosstalk_index-xtalk_proxy->n_active_crosstalk+ii)%xtalk_proxy->xtalk_proxy_size];
		}
		for (int ii=xtalk_proxy->xtalk_proxy_size;ii<2*xtalk_proxy->xtalk_proxy_size;ii++){
			piximp_list[ii] = (PixImpact*) malloc(sizeof(PixImpact));
			CHECK_MALLOC_VOID_STATUS(piximp_list[ii],*status);
		}
		free(xtalk_proxy->xtalk_impacts);
		xtalk_proxy->xtalk_impacts = piximp_list;
		xtalk_proxy->current_crosstalk_index = xtalk_proxy->xtalk_proxy_size;
		xtalk_proxy->xtalk_proxy_size*=2;
	}
	// Copy crosstalk to proxy
	copyPixImpact(xtalk_proxy->xtalk_impacts[xtalk_proxy->current_crosstalk_index % xtalk_proxy->xtalk_proxy_size],impact);
	xtalk_proxy->xtalk_impacts[xtalk_proxy->current_crosstalk_index % xtalk_proxy->xtalk_proxy_size]->nb_pileup=0;
	xtalk_proxy->n_active_crosstalk+=1;
	xtalk_proxy->current_crosstalk_index+=1;
	// Prevent meaningless overrun (TO CHECK)
	if (xtalk_proxy->current_crosstalk_index>=2*xtalk_proxy->xtalk_proxy_size){
		xtalk_proxy->current_crosstalk_index-=xtalk_proxy->xtalk_proxy_size;
	}
}

/** Compute crosstalk influence */
void computeAllCrosstalkInfluence(AdvDet* det,PixImpact * impact,CrosstalkProxy* xtalk_proxy,double* xtalk_energy,int* nb_influences,
		double current_time,int ignore_influence,int full_proxy){
	// If there is no crosstalk in the proxy get out without any influence
	if (xtalk_proxy==NULL || xtalk_proxy->n_active_crosstalk==0){
		*xtalk_energy=0.;
		*nb_influences=0;
		return;
	}

	double energy_influence = 0.;
	int has_influenced= 0;
	PixImpact* crosstalk=NULL;
	int erased_crosstalks=0;
	for (int ii=0;ii<xtalk_proxy->n_active_crosstalk;ii++){
		crosstalk = xtalk_proxy->xtalk_impacts[(xtalk_proxy->current_crosstalk_index-xtalk_proxy->n_active_crosstalk+ii)%xtalk_proxy->xtalk_proxy_size];
		// Compute influence if asked or if the crosstalk is going to be erased (otherwise multiple counting)
		// or if specifically asked to compute all the influences (typically used just before grading of a detected event)
		if ((!ignore_influence) && (full_proxy || current_time-crosstalk->time > -det->crosstalk_timedep->time[0]) ){
			energy_influence=0.;
			has_influenced=computeCrosstalkInfluence(det,impact, crosstalk,&energy_influence);
			if (has_influenced){
				*nb_influences+=crosstalk->nb_pileup+1;
				*xtalk_energy+=energy_influence;
			}
		}
		// If the current crosstalk has no chance of influencing another event, we can clear it
		if (current_time-crosstalk->time > -det->crosstalk_timedep->time[0]){
			erased_crosstalks+=1;
		}
	}
	xtalk_proxy->n_active_crosstalk-=erased_crosstalks;
	// Reboot index if possible
	if (xtalk_proxy->n_active_crosstalk==0){
		xtalk_proxy->current_crosstalk_index=0;
	}
}

/**  Binary search for to find interpolation interval
 *   - return value is the bin [ind,ind+1]
 *   - assume list is sorted ascending */
int binary_search(double val, double* arr, int n){

	if (val < arr[0] || val > arr[n-1]){
		return -1;
	}

	int high=n-1;
	int low=0;
	int mid;
	while (high > low) {
		mid=(low+high)/2;
		if (arr[mid] <= val) {
			low=mid+1;
		} else {
			high=mid;
		}
	}
	return low-1;
}
