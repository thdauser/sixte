// This code contains all the models for the TES tranistion shapes
// as well as the assumptions for frame hits, readout etc...
// Code by J. Wilms, Stephen J. Smith, T. Dauser, P. Peill
// E. Cucchetti, M. Lorenz
//
// current; muA
// J/K -> pJ/K
// power flow: pW/K   (differential conductance)
//   specify the K and n, spit out G
//   practical use: G is preferred
//      -> allow either G or K, print out both
//  add flag to invert signal, do mapping with imin,imax only
//
// ASTRONOMERS BEWARE: This code uses MKS units!
//
//
// Mapping to I&H
// - rename T_start to T0_start
// - use T0 instead of T1 for current temperature

#include "tessim.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <fitsio.h>
#include <sys/time.h>

#include <assert.h>

// Two fluid model parameters
const double two_fluid_n1 = 1.8;
const double two_fluid_n2 = 0.85;
const double two_fluid_Ci = 0.75;
const double two_fluid_Ic0 = 10.25e-3; //[A]
const double two_fluid_Tsat = 0.09; //[K]

// Effect of a frame impact on the bath temperature
double frameImpactEffects(double time_hit, tesparams *tes, int* shift){

	double current_time=tes->time; //current simulation time
	double val=0;

	//If impact is not there yet return 0
	if ((time_hit-current_time)>0){
		return val;
	}

	// If the impact is too much in the past, skip and put in the list to erase
	if ((current_time-time_hit)>MAX_TIME){
		shift+=1;
		return val;
	}

	//Do a binary search of the closest time between current and file
	int n=binary_search(current_time, tes->frame_hit_time_dep, tes->frame_hit_file_length);

	//Linearly interpolate the actual value of DeltaT
	val=tes->frame_hit_shape[n]*((tes->frame_hit_time_dep[n+1]-time_hit)/(tes->frame_hit_time_dep[n+1]-tes->frame_hit_time_dep[n]))+
			tes->frame_hit_shape[n+1]*((time_hit-tes->frame_hit_time_dep[n])/(tes->frame_hit_time_dep[n+1]-tes->frame_hit_time_dep[n]));
	return val;

}

double get_RTI(tesparams *tes, const double Y[]){

	double II=Y[0];
	double TT=Y[1];
	// Computes the value of RTI
	if (tes->twofluid==0){
		return tes->R0+tes->dRdT(tes, Y)*(TT-tes->T_start)+tes->dRdI(tes, Y)*(II-tes->I0_start);
	}
	if (TT>two_fluid_Tsat){
		return tes->R0/tes->bias*0.8;
	}
	double Rn=tes->R0/tes->bias;
	double Cr=0.7+(TT-0.089)*100;
	double Ic=two_fluid_Ic0*pow(1-TT/two_fluid_Tsat, two_fluid_n1);
	double RT=Rn*Cr*pow(1-two_fluid_Ci*pow(Ic/II,two_fluid_n2),1./two_fluid_n2);
	return RT;
}

double get_dRdI(tesparams *tes, const double Y[]){
	
	double II=Y[0];
	double TT=Y[1];
	//Computes the value of dR/dI as a function of the current TES values
	if (tes->twofluid==0){
		return (tes->beta*tes->R0/tes->I0_start);
	}
	if (TT>two_fluid_Tsat){
		return 0;
	}
	double Rn=tes->R0/tes->bias;
	double Cr=0.7+(TT-0.089)*100;
	double factor1=pow(two_fluid_Ic0/II*pow(1-TT/two_fluid_Tsat, two_fluid_n1), two_fluid_n2);
	double factor2=pow(1-two_fluid_Ci*pow(two_fluid_Ic0/II*pow(1-TT/two_fluid_Tsat,two_fluid_n1), two_fluid_n2), 1./two_fluid_n2-1);
	double dRdI=Rn*Cr*two_fluid_Ci/II*factor1*factor2;
	return dRdI;
}

double get_dRdT(tesparams *tes, const double Y[]){
	
	double II=Y[0];
	double TT=Y[1];
	//Computes the value of dR/dT as a function of the current TES values
	if (tes->twofluid==0){
		return (tes->alpha*tes->R0/tes->T_start);
	}
	if (TT>two_fluid_Tsat){
		return 0;
	}
	double Rn=tes->R0/tes->bias;
	double Cr=0.7+(TT-0.089)*100;
	double factor1=pow(two_fluid_Ic0/II, two_fluid_n2)*pow(1-TT/two_fluid_Tsat, two_fluid_n1*two_fluid_n2-1);
	double factor2=pow(1-two_fluid_Ci*pow(two_fluid_Ic0/II*pow(1-TT/two_fluid_Tsat,two_fluid_n1), two_fluid_n2), 1./two_fluid_n2-1);
	double dRdT=Rn*Cr*two_fluid_n1*two_fluid_Ci/tes->T_start*factor1*factor2;
	return dRdT;
}

double get_DeltaTb(tesparams *tes) {

	// If no impact, we are at equilibrium
	if (NULL==tes->frame_hits){
		return tes->Tb_start;
	}

	LinkedFrameImpacts* start = tes->frame_hits;
	double Tb=tes->Tb_start;
	int shift=0;

	//Go through the impacts and find the actual thing
	while(start!=NULL){
		Tb+=frameImpactEffects(start->time, tes, &shift);
		start=tes->frame_hits->next;
	}

	//Shift the list and erase the impact that are no longer useful
	for (int i=0; i<shift;i++){
		LinkedFrameImpacts* toerase=tes->frame_hits;
		tes->frame_hits=tes->frame_hits->next;
		free(toerase);
	}

	return Tb;
}

LinkedFrameImpacts* newLinkedFrameImpacts(int* const status){
	LinkedFrameImpacts* frame_hits=(LinkedFrameImpacts*)malloc(sizeof(LinkedFrameImpacts));
	CHECK_NULL(frame_hits, *status, "memory allocation for LinkedFrameImpacts failed");

	// Initialize pointers with NULL.
	frame_hits->next=NULL;
	return(frame_hits);
}

void freeLinkedFrameImpacts(LinkedFrameImpacts** const list){
  if (NULL!=*list) {
    if (NULL!=(*list)->next) {
    	freeLinkedFrameImpacts(&((*list)->next));
    }
    free(*list);
    *list=NULL;
  }
}

LinkedFrameImpacts* frameHitsList(tesparams *tes, int* const status) {

	//List of impacts
	LinkedFrameImpacts* newlist=NULL;

	if (tes->frame_hit==0){
		return newlist;
	}

	newlist=newLinkedFrameImpacts(status);
	newlist->time=tes->frame_hit_time;

	return newlist;
}

void read_frame_hit_file(tesparams *tes, int* status){
	const char* mode = "r";
	FILE* filetoread;
	filetoread = fopen(tes->frame_hit_file, mode);
	if (filetoread == NULL) {
		SIXT_ERROR("The requested frame hit file is wrong or not in current directory");
        *status=EXIT_FAILURE;
        return;
	}
	int got=2; //I/O check
	int i=0; //counter

	// Allocating memory
	tes->frame_hit_time_dep=(double*)malloc(1000*sizeof(double));
	tes->frame_hit_shape=(double*)malloc(1000*sizeof(double));

	do{
		//If this runs out of memory reallocate more except for first run
		if ((i/1000==0) && (i!=0)){
			tes->frame_hit_time_dep=(double*)realloc(tes->frame_hit_time_dep, (sizeof(tes->frame_hit_time_dep)+1000)*sizeof(double));
			tes->frame_hit_shape=(double*)realloc(tes->frame_hit_shape, (sizeof(tes->frame_hit_shape)+1000)*sizeof(double));
		}
		got=fscanf(filetoread, "%lf; %lf", &tes->frame_hit_time_dep[i], &tes->frame_hit_shape[i]);
		tes->frame_hit_time_dep[i]=tes->frame_hit_time_dep[i]+tes->frame_hit_time; //Shift by the given time
		if (got!=2){
			break;
		}
		i+=1;
	} while(1);
	fclose(filetoread);

	//Allocate the good size of the matrix
	tes->frame_hit_file_length=i;
	tes->frame_hit_time_dep=(double*)realloc(tes->frame_hit_time_dep, tes->frame_hit_file_length*sizeof(double));
	tes->frame_hit_shape=(double*)realloc(tes->frame_hit_shape, tes->frame_hit_file_length*sizeof(double));
	return;
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
