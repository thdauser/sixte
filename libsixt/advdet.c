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


   Copyright 2014-2016 Thorsten Brand, Thomas Dauser, Philippe Peille, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "advdet.h"

/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
  AdvDet* det;
  int status;
};

////////////////////////////////////////////////////////////////////
// Static function declarations.
////////////////////////////////////////////////////////////////////

/** Handler for the start of an XML element. */
static void AdvDetXMLElementStart(void* parsedata,
				   const char* el,
				   const char** attr);
/** Handler for the end of an XML element. */
static void AdvDetXMLElementEnd(void* parsedata, const char* el);

////////////////////////////////////////////////////////////////////
// Program Code
////////////////////////////////////////////////////////////////////

TESNoiseProperties* newTESNoise(int* const status){

  // Allocate memory
  TESNoiseProperties* noise=(TESNoiseProperties*)malloc(sizeof(TESNoiseProperties));
  if(NULL==noise){
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TESNoiseProperties failed");
    return(noise);
  }

  // Initialize values and pointers
  noise->WhiteRMS=0.;
  noise->H0=1.;
  noise->Nz=0;
  noise->Np=0;
  noise->Zeros=NULL;
  noise->Poles=NULL;
  noise->OoFRMS=0.;
  noise->OoFKnee=0.;

  return(noise);
}

void destroyTESNoiseProperties(TESNoiseProperties* noise){
  if (NULL!=noise){
    if(noise->Zeros!=NULL){
      free(noise->Zeros);
    }
    if(noise->Poles!=NULL){
      free(noise->Poles);
    }
    noise->WhiteRMS=0.;
    noise->H0=1.;
    noise->Nz=0;
    noise->Np=0;
    noise->OoFRMS=0.;
    noise->OoFKnee=0.;
  }
  free(noise);
  noise=NULL;
}

AdvPix* newAdvPix(int* const status){
  // Allocate memory
  AdvPix* pix=(AdvPix*)malloc(sizeof(AdvPix));
  if(pix==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("Unable to allocate memory for detector pixel");
    return(pix);
  }
  // Initialize values and pointers
  pix->sx=0.;
  pix->sy=0.;
  pix->width=0.;
  pix->height=0.;
  pix->pindex=0;
  //pix->version=NULL;
  pix->profVersionID=0;
  pix->ADCOffset=0;
  pix->calfactor=0.;
  pix->TESNoise=NULL;
  pix->ngrades=0;
  pix->grades=NULL;
  pix->global_grading=0;
  pix->arffile=NULL;
  pix->arf=NULL;
  pix->freq=0.;
  pix->resfreq=0.;
  pix->channel=NULL;
  pix->electrical_cross_talk=NULL;
  pix->thermal_cross_talk=NULL;
  //pix->intermodulation_cross_talk=NULL;
  pix->prop_cross_talk=NULL;
  pix->der_cross_talk=NULL;
  pix->tes_type=NULL;
  pix->tes=NULL;
  pix->row=0;

  return(pix);
}

AdvRecons* newAdvRecons(int* const status){     //SIRENA
  // Allocate memory
  AdvRecons* recons=(AdvRecons*)malloc(sizeof(AdvRecons));
  if(recons==NULL){
    *status=EXIT_FAILURE;
    SIXT_ERROR("Unable to allocate memory for reconstruction info");
    return(recons);
  }
  // Initialize values and pointers
  recons->ngrades=0;
  recons->grades=NULL;

  return(recons);
}


void freeAdvPix(AdvPix* pix){
  if(NULL!=pix){
    int g=pix->ngrades;
    destroyTESNoiseProperties(pix->TESNoise);
    freeGrading(pix);
    freeMatrixCrossTalk(pix->thermal_cross_talk);
    freeMatrixEnerdepCrossTalk(pix->electrical_cross_talk, g);
    freeMatrixPropCrossTalk(pix->prop_cross_talk);
    freeMatrixDerCrossTalk(pix->der_cross_talk);
  }
}

void freeAdvRecons(AdvRecons* recons){      //SIRENA
  if(NULL!=recons){
    //int g=recons->ngrades;
    freeGradingRecons(recons);
  }
}

void freeGrading(AdvPix* pix){
    for (int i=0;i<pix->ngrades;i++){
    	free(pix->grades[i].rmffile);
    }
    free(pix->grades);
    pix->ngrades=0;
    pix->grades=NULL;
    pix->global_grading=0;
}

void freeGradingRecons(AdvRecons* recons){      //SIRENA
    free(recons->grades);
    recons->ngrades=0;
    recons->grades=NULL;
}

/** Free the FDM System Structure */
void freeFDMSystem(FDMSystem* fdmsys){
        int ii;
        free(fdmsys->omega_array);
        free(fdmsys->res_omega_array);
        free(fdmsys->X_L);
        free(fdmsys->capFac);
        for (ii=0; ii<fdmsys->num_pixels; ii++){
                free(fdmsys->Z_array[ii]);
        }
        free(fdmsys->Z_array);
        free(fdmsys);
}

void freeReadoutChannels(ReadoutChannels* rc){
	if ( rc != NULL ){
		// Free the FDM systems
       for (int ii=0; ii<rc->num_channels; ii++){
    	   if (rc->channels[ii].fdmsys != NULL){
    		   freeFDMSystem(rc->channels[ii].fdmsys);
           }
    	   if (rc->channels[ii].pixels!=NULL){
    		   //Only release first pointer, the rest was released via destroyAdvDet (<-- Dangerous game here! Edoardo's comment)
    		   free(rc->channels[ii].pixels);
    	   }
       }
       free(rc->channels);
       free(rc->df_information_band);
	}
	free(rc);
}

void freeImodTab(ImodTab* tab, int gr){
	if (tab != NULL){
		for(int g=0; g<gr;g++){
			if (tab[g].matrix!=NULL){
				for(int kk=0; kk < tab[g].n_freq ; kk++){ //  --- freq --- //
					if (tab[g].matrix[kk]!=NULL){
						for(int ii=0; ii < tab[g].n_dt ; ii++){ // --- dt --- //
							if (tab[g].matrix[kk][ii]!=NULL){
								for(int jj=0; jj < tab[g].n_ampl ; jj++){
									free(tab[g].matrix[kk][ii][jj]);
								}
								free(tab[g].matrix[kk][ii]);
							}
						}  // --- END dt --- //
						free(tab[g].matrix[kk]);
					}
				}   //  --- END freq --- //
				free(tab[g].matrix);
			}
			free(tab[g].ampl);
			free(tab[g].dt);
			free(tab[g].freq);
		}
		free(tab);
	}
}

void freeElecTab(ElecTab* tab,int gr){
	if (tab != NULL){
		for(int g=0; g<gr;g++){
			if (tab[g].matrix!=NULL){
				for(int kk=0; kk < tab[g].n_freq_s ; kk++){ //  --- freq_p --- //
					if (tab[g].matrix[kk]!=NULL){
						for(int ii=0; ii < tab[g].n_freq_p ; ii++){ // --- freq_s --- //
							if (tab[g].matrix[kk][ii]!=NULL){
								free(tab[g].matrix[kk][ii]);
							}
						}  // --- END freq_s --- //
						free(tab[g].matrix[kk]);
					}
				}   //  --- END freq_p --- //

			}
			free(tab[g].freq_s);
			free(tab[g].freq_p);
			free(tab[g].ener_p);
		}
		free(tab);
	}
}

//TODO Check if this is good
void freeTDMTab(TDMTab* tab, int gr){
	if (tab!=NULL){
		for(int g=0; g<gr;g++){
			if (tab[g].matrix!=NULL){
				for(int kk=0; kk < tab[g].n_samples ; kk++){ //  --- samples --- //
					if (tab[g].matrix[kk]!=NULL){
						for(int ii=0; ii < tab[g].n_ener_p ; ii++){ // --- ener_p --- //
							if (tab[g].matrix[kk][ii]!=NULL){
								free(tab[g].matrix[kk][ii]);
							}
						}   // --- END ener_p --- //
						free(tab[g].matrix[kk]);
					}
				}   //  --- END samples --- //
			}
			free(tab[g].matrix);
			free(tab[g].samples);
			free(tab[g].ener_p);
			free(tab[g].ener_v);
			tab[g].n_ener_p=0;
			tab[g].n_ener_v=0;
			tab[g].n_samples=0;
		}
		free(tab);
	}
}

void freeCrosstalk(AdvDet* det, int gr){

	if (NULL!=det->channel_file){
		free(det->channel_file);
		det->channel_file=NULL;
	}
	if (NULL!=det->channel_resfreq_file){
		free(det->channel_resfreq_file);
		det->channel_resfreq_file=NULL;
	}

	freeReadoutChannels(det->readout_channels);
	//Freeing tables for every grade
	if(gr!=0){
		int len=det->xt_num_thermal;

		//Freeing electrical cross-talk
		free(det->crosstalk_elec_timedep_file);
		det->crosstalk_elec_timedep_file=NULL;
		free(det->crosstalk_elec_file);
		det->crosstalk_elec_file=NULL;

		freeElecTab(det->crosstalk_elec,gr);
		det->crosstalk_elec=NULL;
		freeCrosstalkTimedep(det->crosstalk_elec_timedep,gr);
		free(det->crosstalk_elec_timedep);
		det->crosstalk_elec_timedep=NULL;

		// Freeing intermodulation cross-talk
		free(det->crosstalk_intermod_file);
		det->crosstalk_intermod_file=NULL;

		freeImodTab(det->crosstalk_imod_table,gr);
		det->crosstalk_imod_table=NULL;

		//Freeing thermal cross-talk
		if(det->crosstalk_ther_timedep!=NULL){
			for(int i=0; i<len; i++){
				freeCrosstalkTimedep(det->crosstalk_ther_timedep[i],gr);
				free(det->crosstalk_ther_timedep[i]);
				det->crosstalk_ther_timedep[i]=NULL;
				free(det->crosstalk_thermal_timedep_file[i]);
				det->crosstalk_thermal_timedep_file[i]=NULL;
			}
			free(det->crosstalk_thermal_timedep_file);
			det->crosstalk_thermal_timedep_file=NULL;
			free(det->crosstalk_ther_timedep);
			det->crosstalk_ther_timedep=NULL;
		}
		free(det->xt_dist_thermal);
		det->xt_dist_thermal=NULL;
		free(det->xt_weight_thermal);
		det->xt_weight_thermal=NULL;

		//Freeing Proportional Cross-talk
		freeTDMTab(det->crosstalk_TDM_prop, gr);
		det->crosstalk_TDM_prop=NULL;

		//Freeing Derivative Cross-talk
		freeTDMTab(det->crosstalk_TDM_der, gr);
		det->crosstalk_TDM_der=NULL;
	}
}

AdvDet* newAdvDet(int* const status){

  // Allocate memory.
  AdvDet* det=(AdvDet*)malloc(sizeof(AdvDet));
  if (NULL==det) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for AdvDet failed");
    return(det);
  }

  // Initialize all pointers with NULL.
  det->pix=NULL;
  det->recons=NULL;     //SIRENA
  det->filename=NULL;
  det->filepath=NULL;
  det->sx=0.;
  det->sy=0.;
  det->npix=0;
  det->nrecons=0;   //SIRENA
  det->cpix=0;
  det->SampleFreq=-1.0;
  det->tesnoisefilter=0;
  det->inpixel=0;
  det->oof_activated=0;
  det->rmf_library=NULL;
  det->arf_library=NULL;
  det->tdm=0;
  det->max_rows=0;

  det->xt_dist_thermal=NULL;
  det->xt_weight_thermal=NULL;
  det->xt_num_thermal=0;

  det->channel_file=NULL;
  det->readout_channels=NULL;
  det->channel_resfreq_file=NULL;

  det->crosstalk_intermod_file=NULL;
  det->crosstalk_imod_table=NULL;

  det->crosstalk_thermal_timedep_file=NULL;
  det->crosstalk_ther_timedep=NULL;

  det->crosstalk_elec_file=NULL;
  det->crosstalk_elec_timedep_file=NULL;
  det->crosstalk_elec_timedep=NULL;
  det->crosstalk_elec=NULL;
  det->scaling=1;

  det->TDM_prop_file=NULL;
  det->crosstalk_TDM_prop=NULL;
  det->prop_TDM_scaling_1=1;
  det->prop_TDM_scaling_2=1;

  det->TDM_der_file=NULL;
  det->crosstalk_TDM_der=NULL;
  det->der_TDM_scaling=1;

  det->threshold_event_lo_keV=0.;
  det->crosstalk_id=0;

  det->tes_type_file=NULL;
  det->L_Common=0.;
  det->C_Common=-1.; // default -1 means no common capacitance

  return(det);
}

void destroyAdvDet(AdvDet **det){

	if(NULL!=(*det)){
		int gr=0;
		if(NULL!=(*det)->pix){
			gr =(*det)->pix[0].ngrades;
			for(int i=0;i<(*det)->npix;i++){
				freeAdvPix(&(*det)->pix[i]);
			}
			free((*det)->pix);
		}
		//int grRecons=0;
		if(NULL!=(*det)->recons){             //SIRENA
			//grRecons =(*det)->recons[0].ngrades;
			for(int i=0;i<(*det)->nrecons;i++){
				freeAdvRecons(&(*det)->recons[i]);
			}
			free((*det)->recons);
		}
		if(NULL!=(*det)->filename){
			free((*det)->filename);
		}
		if(NULL!=(*det)->filepath){
			free((*det)->filepath);
		}
		if(NULL!=(*det)->tes_type_file){
			free((*det)->tes_type_file);
		}
		freeRMFLibrary((*det)->rmf_library);
		freeARFLibrary((*det)->arf_library);

		freeCrosstalk(*(det), gr);
		free(*(det));
	}
}

int CheckAdvPixImpact(AdvPix pix, Impact *imp){

  // Calculate impact coordinates in respect to the
  // pixel coordinate system
  double u, v;

  u = imp->position.x - pix.sx;
  v = imp->position.y - pix.sy;

  // Calculate half width and height of the rectangular pixel
  double deltu, deltv;
  deltu=pix.width/2.;
  deltv=pix.height/2.;

  // Check if the impact lies in the rectangular pixel.
  // Return 1 if yes, 0 if not.
  if((u >= -deltu) && (u <= deltu) && (v >= -deltv) && (v <= deltv)){
    return 1;
  }else{
    return 0;
  }
}

void CalcAdvPixImpact(AdvPix pix, Impact *imp, PixImpact *piximp){

  // Calculate impact coordinates in respect to the
  // pixel coordinate system
  double u, v;

  u = imp->position.x - pix.sx;
  v = imp->position.y - pix.sy;

  // Fill piximp fields
  piximp->time = imp->time;
  piximp->energy = imp->energy;
  piximp->ph_id = imp->ph_id;
  piximp->src_id = imp->src_id;
  piximp->detposition.x = imp->position.x;
  piximp->detposition.y = imp->position.y;
  piximp->pixposition.x = u;
  piximp->pixposition.y = v;
}

int AdvImpactList(AdvDet *det, Impact *imp, PixImpact **piximp){

  // Duplicate the impact but transform the coordinates into
  // the detector coordinate system
  Impact detimp;

  detimp.time = imp->time;
  detimp.energy = imp->energy;
  detimp.ph_id = imp->ph_id;
  detimp.src_id = imp->src_id;
  detimp.position.x = imp->position.x - det->sx;
  detimp.position.y = imp->position.y - det->sy;

  int nimpacts=0;

  // loop over all pixels and check for hit
  int ii;

  for(ii=0; ii<det->npix; ii++){
    if(CheckAdvPixImpact(det->pix[ii], &detimp)!=0){
      nimpacts++;
      *piximp=realloc(*piximp, nimpacts*sizeof(**piximp));
      (*piximp)[nimpacts-1].pixID=(long)ii;
      CalcAdvPixImpact(det->pix[ii], &detimp, &((*piximp)[nimpacts-1]));
    }
  }
  return nimpacts;
}


void parseAdvDetXML(AdvDet* const det,
	       const char* const filename,
	       int* const status){

  headas_chat(5, "read advanced detector setup from XML file '%s' ...\n", filename);

  // Read the XML data from the file.
  // Open the specified file.
  printf("Read file %s\n", filename);
  FILE* xmlfile=fopen(filename, "r");
  if (NULL==xmlfile) {
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "failed opening XML "
	    "file '%s' for read access", filename);
    SIXT_ERROR(msg);
    return;
  }

  // The data is read from the XML file and stored in xmlbuffer
  // without any modifications.
  struct XMLBuffer* xmlbuffer=newXMLBuffer(status);
  CHECK_STATUS_VOID(*status);

  // Input buffer with an additional byte at the end for the
  // termination of the string.
  const int buffer_size=256;
  char buffer[buffer_size+1];
  // Number of chars in buffer.
  int len;

  // Read all data from the file.
  do {
    // Get a piece of input into the buffer.
    len=fread(buffer, 1, buffer_size, xmlfile);
    buffer[len]='\0'; // Terminate the string.
    addString2XMLBuffer(xmlbuffer, buffer, status);
    CHECK_STATUS_VOID(*status);
  } while (!feof(xmlfile));

  // Close the file handler to the XML file.
  fclose(xmlfile);

  // Before expanding loops in the XML file, add the included code to it.
  expandIncludesXML(xmlbuffer, filename, status);
  CHECK_STATUS_VOID(*status);

  // Before actually parsing the XML code, expand the loops and
  // arithmetic operations in the XML description.
  // The expansion algorithm repeatetly scans the XML code and
  // searches for loop tags. It replaces the loop tags by repeating
  // the contained XML code.
  expandXML(xmlbuffer, status);
  CHECK_STATUS_VOID(*status);

  // Expand the eventual hexagonloop structure
  expandHexagon(xmlbuffer,status);
  CHECK_STATUS_VOID(*status);

  // Parse XML code in the xmlbuffer using the expat library.
  // Get an XML_Parser object.
  XML_Parser parser=XML_ParserCreate(NULL);
  if (NULL==parser) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("could not allocate memory for XML parser");
    return;
  }

  // Set data that is passed to the handler functions.
  struct XMLParseData xmlparsedata={
    .det  =det,
    .status=EXIT_SUCCESS
  };
  XML_SetUserData(parser, &xmlparsedata);

  // Set the handler functions.
  XML_SetElementHandler(parser, AdvDetXMLElementStart, AdvDetXMLElementEnd);

  // Parse all the data in the string buffer.
  const int done=1;
  if (!XML_Parse(parser, xmlbuffer->text, strlen(xmlbuffer->text), done)) {
    // Parse error.
    *status=EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "failed parsing XML file '%s':\n%s\n",
	    filename, XML_ErrorString(XML_GetErrorCode(parser)));
    printf("%s", xmlbuffer->text);
    SIXT_ERROR(msg);
    return;
  }

  // Check that the correct number of pixels was found
  if ((det->nrecons == 0) && (det->cpix!=det->npix)){     // SIRENA
  //if(det->cpix!=det->npix){
	  *status = EXIT_FAILURE;
	  SIXT_ERROR("Number of pixels given at pixdetector level does not match the number of pixel subelements");
	  det->npix=det->cpix;
  }

  // Check for errors.
  if (EXIT_SUCCESS!=xmlparsedata.status) {
    *status=xmlparsedata.status;
  }

  // Release memory.
  XML_ParserFree(parser);

  // Remove the XML string buffer.
  freeXMLBuffer(&xmlbuffer);
}

static void AdvDetXMLElementStart(void* parsedata,
				   const char* el,
				   const char** attr)
{
	struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;

	// Check if an error has occurred previously.
	CHECK_STATUS_VOID(xmlparsedata->status);

	// Convert the element to an upper case string.
	char Uelement[MAXMSG];
	strcpy(Uelement, el);
	strtoupper(Uelement);

	// Check for advanced detector definition.
	if (!strcmp(Uelement, "PIXDETECTOR")) {
		// Determine npix
		char npix[MAXMSG];
		char readout[MAXMSG];
		getXMLAttributeString(attr, "NPIX", npix);

		xmlparsedata->det->npix=atoi(npix);
		if(xmlparsedata->det->npix<1){
			SIXT_ERROR("Number of pixels in advanced detector description less than 1.");
			return;
		}
		xmlparsedata->det->pix=(AdvPix*)malloc(xmlparsedata->det->npix*sizeof(AdvPix));
		if(xmlparsedata->det->pix==NULL){
			SIXT_ERROR("Unable to allocate memory for advanced detector pixel array.");
			return;
		}

		xmlparsedata->det->sx=getXMLAttributeDouble(attr, "XOFF");
		xmlparsedata->det->sy=getXMLAttributeDouble(attr, "YOFF");
		getXMLAttributeString(attr, "READOUT", readout);

		if (!strcmp(readout,"fdm") || !strcmp(readout,"")){
			xmlparsedata->det->tdm=0;
			printf("Readout type is FDM\n");
		} else if (!strcmp(readout,"tdm")){
			xmlparsedata->det->tdm=1;
			printf("Readout type is TDM\n");
		} else {
			SIXT_ERROR("Readout type is not understood");
			return;
		}


		for (int ii=0;ii<xmlparsedata->det->npix;ii++){
			xmlparsedata->det->pix[ii].TESNoise=NULL;
			xmlparsedata->det->pix[ii].grades=NULL;
			xmlparsedata->det->pix[ii].ngrades=0;
			xmlparsedata->det->pix[ii].global_grading=0;
			xmlparsedata->det->pix[ii].channel=NULL;
			xmlparsedata->det->pix[ii].thermal_cross_talk=NULL;
			xmlparsedata->det->pix[ii].electrical_cross_talk=NULL;
			xmlparsedata->det->pix[ii].prop_cross_talk=NULL;
			xmlparsedata->det->pix[ii].der_cross_talk=NULL;
		}
    } else if (!strcmp(Uelement, "RECONSTRUCTION")) {   //SIRENA
		xmlparsedata->det->nrecons = 1;
        xmlparsedata->det->recons=(AdvRecons*)malloc(xmlparsedata->det->nrecons*sizeof(AdvRecons));

		for (int ii=0;ii<xmlparsedata->det->nrecons;ii++){
            xmlparsedata->det->recons[ii].grades=NULL;
            xmlparsedata->det->recons[ii].ngrades=0;
		}
	} else if (!strcmp(Uelement, "PIXEL")) {
		//if ((xmlparsedata->det->cpix) >= (xmlparsedata->det->npix)) {
        if ((xmlparsedata->det->cpix) > (xmlparsedata->det->npix)) {    //SIRENA
			xmlparsedata->status=EXIT_FAILURE;
			SIXT_ERROR("Number of pixels given at pixdetector level lower than the number of pixel subelements");
			return;
		}
		xmlparsedata->det->inpixel=1;
	} else if (!strcmp(Uelement, "SHAPE")){
		double posx, posy;
		posx=getXMLAttributeDouble(attr, "POSX");
		posy=getXMLAttributeDouble(attr, "POSY");
		if (xmlparsedata->det->inpixel){
			xmlparsedata->det->pix[xmlparsedata->det->cpix].sx=posx*getXMLAttributeDouble(attr, "DELX");
			xmlparsedata->det->pix[xmlparsedata->det->cpix].sy=posy*getXMLAttributeDouble(attr, "DELY");
			xmlparsedata->det->pix[xmlparsedata->det->cpix].width=getXMLAttributeDouble(attr, "WIDTH");
			xmlparsedata->det->pix[xmlparsedata->det->cpix].height=getXMLAttributeDouble(attr, "HEIGHT");
		} else {
			xmlparsedata->status=EXIT_FAILURE;
			SIXT_ERROR("XML syntax error: shape used outside of pixel");
			return;
		}
        } else if (!strcmp(Uelement,"PIXTES")) {
                if (xmlparsedata->det->inpixel){
                        xmlparsedata->det->pix[xmlparsedata->det->cpix].tes_type=(char*)malloc(MAXFILENAME*sizeof(char));
		        CHECK_MALLOC_VOID(xmlparsedata->det->pix[xmlparsedata->det->cpix].tes_type);
                        getXMLAttributeString(attr, "TYPE", xmlparsedata->det->pix[xmlparsedata->det->cpix].tes_type);
		        if (strlen(xmlparsedata->det->pix[xmlparsedata->det->cpix].tes_type) == 0){
			        free(xmlparsedata->det->pix[xmlparsedata->det->cpix].tes_type);
			        xmlparsedata->det->pix[xmlparsedata->det->cpix].tes_type=NULL;
                        }
                } else {
                        xmlparsedata->status=EXIT_FAILURE;
                        SIXT_ERROR("XML syntax error: pixtes used outside of pixel");
                }
	} else if (!strcmp(Uelement,"PULSESHAPE")){
		if (xmlparsedata->det->inpixel){
			xmlparsedata->det->pix[xmlparsedata->det->cpix].calfactor=getXMLAttributeDouble(attr, "CALFACTOR");
			xmlparsedata->det->pix[xmlparsedata->det->cpix].ADCOffset=getXMLAttributeInt(attr, "ADCOFFSET");
			getXMLAttributeString(attr, "TESPROFVER", xmlparsedata->det->pix[xmlparsedata->det->cpix].version);
			xmlparsedata->det->pix[xmlparsedata->det->cpix].profVersionID=-1;
		} else {
			for (int i=0;i<xmlparsedata->det->npix;i++){
				xmlparsedata->det->pix[i].calfactor=getXMLAttributeDouble(attr, "CALFACTOR");
				xmlparsedata->det->pix[i].ADCOffset=getXMLAttributeInt(attr, "ADCOFFSET");
				getXMLAttributeString(attr, "TESPROFVER", xmlparsedata->det->pix[i].version);
				xmlparsedata->det->pix[i].profVersionID=-1;
			}
		}
	} else if (!strcmp(Uelement, "TESWHITENOISE")) {
		if (xmlparsedata->det->inpixel){
			if(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise==NULL){
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise=newTESNoise(&(xmlparsedata->status));
				CHECK_STATUS_VOID(xmlparsedata->status);
			}
			xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->WhiteRMS=getXMLAttributeDouble(attr, "RMS");
		} else {
			for (int i=0;i<xmlparsedata->det->npix;i++){
				if(xmlparsedata->det->pix[i].TESNoise==NULL){
					xmlparsedata->det->pix[i].TESNoise=newTESNoise(&(xmlparsedata->status));
					CHECK_STATUS_VOID(xmlparsedata->status);
				}
				xmlparsedata->det->pix[i].TESNoise->WhiteRMS=getXMLAttributeDouble(attr, "RMS");
			}
		}
	} else if (!strcmp(Uelement, "TESOOFNOISE")) {
		if (xmlparsedata->det->inpixel){
			if(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise==NULL){
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise=newTESNoise(&(xmlparsedata->status));
				CHECK_STATUS_VOID(xmlparsedata->status);
			}
			xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->OoFRMS=getXMLAttributeDouble(attr, "RMS");
			xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->OoFKnee=getXMLAttributeDouble(attr, "FKNEE");
			if (!xmlparsedata->det->oof_activated){
				xmlparsedata->det->oof_activated=1;
			}
		} else {
			for (int i=0;i<xmlparsedata->det->npix;i++){
				if(xmlparsedata->det->pix[i].TESNoise==NULL){
					xmlparsedata->det->pix[i].TESNoise=newTESNoise(&(xmlparsedata->status));
					CHECK_STATUS_VOID(xmlparsedata->status);
				}
				xmlparsedata->det->pix[i].TESNoise->OoFRMS=getXMLAttributeDouble(attr, "RMS");
				xmlparsedata->det->pix[i].TESNoise->OoFKnee=getXMLAttributeDouble(attr, "FKNEE");
			}
			xmlparsedata->det->oof_activated=1;
		}
	} else if (!strcmp(Uelement, "TESNOISEFILTER")) {
		if (xmlparsedata->det->inpixel){
			if(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise==NULL){
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise=newTESNoise(&(xmlparsedata->status));
				CHECK_STATUS_VOID(xmlparsedata->status);
			}
			xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->H0=getXMLAttributeDouble(attr, "NORM");
			xmlparsedata->det->tesnoisefilter=1;
		} else {
			for (int i=0;i<xmlparsedata->det->npix;i++){
				if(xmlparsedata->det->pix[i].TESNoise==NULL){
					xmlparsedata->det->pix[i].TESNoise=newTESNoise(&(xmlparsedata->status));
					CHECK_STATUS_VOID(xmlparsedata->status);
				}
				xmlparsedata->det->pix[i].TESNoise->H0=getXMLAttributeDouble(attr, "NORM");
				xmlparsedata->det->tesnoisefilter=1;
				xmlparsedata->det->pix[i].TESNoise->Np=0;
				xmlparsedata->det->pix[i].TESNoise->Nz=0;
			}
		}
	} else if (!strcmp(Uelement, "NOISEPOLE")) {
		if(xmlparsedata->det->tesnoisefilter==1){
			if(xmlparsedata->det->inpixel){
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Np++;
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Poles=(double*)realloc(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Poles,
						xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Np*sizeof(double));
				if(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Poles==NULL){
					xmlparsedata->status=EXIT_FAILURE;
					SIXT_ERROR("Realloc of noise poles array failed.");
					return;
				}
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Poles[xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Np-1]=getXMLAttributeDouble(attr, "POLE");
			} else {
				for (int i=0;i<xmlparsedata->det->npix;i++){
					xmlparsedata->det->pix[i].TESNoise->Np++;
					xmlparsedata->det->pix[i].TESNoise->Poles=(double*)realloc(xmlparsedata->det->pix[i].TESNoise->Poles,
							xmlparsedata->det->pix[i].TESNoise->Np*sizeof(double));
					if(xmlparsedata->det->pix[i].TESNoise->Poles==NULL){
						xmlparsedata->status=EXIT_FAILURE;
						SIXT_ERROR("Realloc of noise poles array failed.");
						return;
					}
					xmlparsedata->det->pix[i].TESNoise->Poles[xmlparsedata->det->pix[i].TESNoise->Np-1]=getXMLAttributeDouble(attr, "POLE");
				}
			}
		}else{
			xmlparsedata->status=EXIT_FAILURE;
			SIXT_ERROR("XML syntax error: noisepole used outside of tesnoisefilter.");
			return;
		}
	} else if (!strcmp(Uelement, "NOISEZERO")) {
		if(xmlparsedata->det->tesnoisefilter==1){
			if(xmlparsedata->det->inpixel){
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Nz++;
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Zeros=(double*)realloc(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Zeros,
						xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Nz*sizeof(double));
				if(xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Zeros==NULL){
					xmlparsedata->status=EXIT_FAILURE;
					SIXT_ERROR("Realloc of noise zeros array failed.");
					return;
				}
				xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Zeros[xmlparsedata->det->pix[xmlparsedata->det->cpix].TESNoise->Nz-1]=getXMLAttributeDouble(attr, "ZERO");
			} else {
				for (int i=0;i<xmlparsedata->det->npix;i++){
					xmlparsedata->det->pix[i].TESNoise->Nz++;
					xmlparsedata->det->pix[i].TESNoise->Zeros=(double*)realloc(xmlparsedata->det->pix[i].TESNoise->Zeros,
							xmlparsedata->det->pix[i].TESNoise->Nz*sizeof(double));
					if(xmlparsedata->det->pix[i].TESNoise->Zeros==NULL){
						xmlparsedata->status=EXIT_FAILURE;
						SIXT_ERROR("Realloc of noise zeros array failed.");
						return;
					}
					xmlparsedata->det->pix[i].TESNoise->Zeros[xmlparsedata->det->pix[i].TESNoise->Nz-1]=getXMLAttributeDouble(attr, "ZERO");
				}
			}
		}else{
			xmlparsedata->status=EXIT_FAILURE;
			SIXT_ERROR("XML syntax error: noisezero used outside of tesnoisefilter.");
			return;
		}
	} else if (!strcmp(Uelement, "PIXARF")) {
		if (xmlparsedata->det->inpixel){
			char arffile[MAXFILENAME];
			getXMLAttributeString(attr, "FILENAME", arffile);
			xmlparsedata->det->pix[xmlparsedata->det->cpix].arffile=strndup(arffile,MAXFILENAME);
			xmlparsedata->det->pix[xmlparsedata->det->cpix].arf=NULL;
		} else {
			for (int i=0;i<xmlparsedata->det->npix;i++){
				char arffile[MAXFILENAME];
				getXMLAttributeString(attr, "FILENAME", arffile);
				xmlparsedata->det->pix[i].arffile=strndup(arffile,MAXFILENAME);
				xmlparsedata->det->pix[i].arf=NULL;
			}
		}
	}  else if (!strcmp(Uelement,"GRADING")){
        if (xmlparsedata->det->nrecons == 1)    //SIRENA
        {
            if (xmlparsedata->det->inpixel){
                xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades=realloc(xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades,
                                                                               (xmlparsedata->det->recons[xmlparsedata->det->nrecons].ngrades+1)*sizeof(*(xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades)));
                if (xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades==NULL){
                    xmlparsedata->status=EXIT_FAILURE;
                    SIXT_ERROR("Realloc of grading array failed.");
                    return;
                }
                xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades[xmlparsedata->det->recons[xmlparsedata->det->nrecons].ngrades].value=getXMLAttributeInt(attr, "NUM");
                xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades[xmlparsedata->det->recons[xmlparsedata->det->nrecons].ngrades].gradelim_pre=getXMLAttributeLong(attr, "PRE");
                xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades[xmlparsedata->det->recons[xmlparsedata->det->nrecons].ngrades].gradelim_post=getXMLAttributeLong(attr, "POST");
                xmlparsedata->det->recons[xmlparsedata->det->nrecons].grades[xmlparsedata->det->recons[xmlparsedata->det->nrecons].ngrades].grade_preBuffer=getXMLAttributeLong(attr, "PB");

                xmlparsedata->det->recons[xmlparsedata->det->nrecons].ngrades++;
            } else {
                for (int i=0;i<xmlparsedata->det->nrecons;i++){
                    xmlparsedata->det->recons[i].grades=realloc(xmlparsedata->det->recons[i].grades,
                                                             (xmlparsedata->det->recons[i].ngrades+1)*sizeof(*(xmlparsedata->det->recons[i].grades)));
                    if (xmlparsedata->det->recons[i].grades==NULL){
                        xmlparsedata->status=EXIT_FAILURE;
                        SIXT_ERROR("Realloc of grading array failed.");
                        return;
                    }
                    xmlparsedata->det->recons[i].grades[xmlparsedata->det->recons[i].ngrades].value=getXMLAttributeInt(attr, "NUM");
                    xmlparsedata->det->recons[i].grades[xmlparsedata->det->recons[i].ngrades].gradelim_pre=getXMLAttributeLong(attr, "PRE");
                    xmlparsedata->det->recons[i].grades[xmlparsedata->det->recons[i].ngrades].gradelim_post=getXMLAttributeLong(attr, "POST");
                    xmlparsedata->det->recons[i].grades[xmlparsedata->det->recons[i].ngrades].grade_preBuffer=getXMLAttributeLong(attr, "PB");
                    xmlparsedata->det->recons[i].grades[xmlparsedata->det->recons[i].ngrades].rmf=NULL;

                    xmlparsedata->det->recons[i].ngrades++;
                }
            }
        }
        else
        {
            if (xmlparsedata->det->inpixel){
                if (xmlparsedata->det->pix[xmlparsedata->det->cpix].global_grading) {
                    freeGrading(&(xmlparsedata->det->pix[xmlparsedata->det->cpix]));
                    xmlparsedata->det->pix[xmlparsedata->det->cpix].global_grading=0;
                }
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades=realloc(xmlparsedata->det->pix[xmlparsedata->det->cpix].grades,
                                                                               (xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades+1)*sizeof(*(xmlparsedata->det->pix[xmlparsedata->det->cpix].grades)));
                if (xmlparsedata->det->pix[xmlparsedata->det->cpix].grades==NULL){
                    xmlparsedata->status=EXIT_FAILURE;
                    SIXT_ERROR("Realloc of grading array failed.");
                    return;
                }
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades[xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades].value=getXMLAttributeInt(attr, "NUM");
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades[xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades].gradelim_pre=getXMLAttributeLong(attr, "PRE");
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades[xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades].gradelim_post=getXMLAttributeLong(attr, "POST");
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades[xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades].grade_preBuffer=getXMLAttributeLong(attr, "PB");
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades[xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades].rmf=NULL;

                char rmffile[MAXFILENAME];
                getXMLAttributeString(attr, "RMF", rmffile);
                xmlparsedata->det->pix[xmlparsedata->det->cpix].grades[xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades].rmffile=strndup(rmffile,MAXFILENAME);
                xmlparsedata->det->pix[xmlparsedata->det->cpix].ngrades++;
            } else {
                for (int i=0;i<xmlparsedata->det->npix;i++){
                    xmlparsedata->det->pix[i].grades=realloc(xmlparsedata->det->pix[i].grades,
                                                             (xmlparsedata->det->pix[i].ngrades+1)*sizeof(*(xmlparsedata->det->pix[i].grades)));
                    if (xmlparsedata->det->pix[i].grades==NULL){
                        xmlparsedata->status=EXIT_FAILURE;
                        SIXT_ERROR("Realloc of grading array failed.");
                        return;
                    }
                    xmlparsedata->det->pix[i].grades[xmlparsedata->det->pix[i].ngrades].value=getXMLAttributeInt(attr, "NUM");
                    xmlparsedata->det->pix[i].grades[xmlparsedata->det->pix[i].ngrades].gradelim_pre=getXMLAttributeLong(attr, "PRE");
                    xmlparsedata->det->pix[i].grades[xmlparsedata->det->pix[i].ngrades].gradelim_post=getXMLAttributeLong(attr, "POST");
                    xmlparsedata->det->pix[i].grades[xmlparsedata->det->pix[i].ngrades].grade_preBuffer=getXMLAttributeLong(attr, "PB");
                    xmlparsedata->det->pix[i].grades[xmlparsedata->det->pix[i].ngrades].rmf=NULL;

                    char rmffile[MAXFILENAME];
                    getXMLAttributeString(attr, "RMF", rmffile);
                    xmlparsedata->det->pix[i].grades[xmlparsedata->det->pix[i].ngrades].rmffile=strndup(rmffile,MAXFILENAME);
                    xmlparsedata->det->pix[i].ngrades++;
                    xmlparsedata->det->pix[i].global_grading=1;
                }
            }
        }
	} else if(!strcmp(Uelement, "TESPROFILE")){
		getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->tesproffilename);
		double new_samplefreq=getXMLAttributeDouble(attr, "SAMPLEFREQ");
		if ((xmlparsedata->det->SampleFreq!=-1) && (xmlparsedata->det->SampleFreq!=new_samplefreq)){
			SIXT_ERROR("Incompatible sampling frequency values encountered.");
			return;
		}
		xmlparsedata->det->SampleFreq=new_samplefreq;
	} else if(!strcmp(Uelement, "SAMPLEFREQ")){
		double new_samplefreq=getXMLAttributeDouble(attr, "VALUE");
		if ((xmlparsedata->det->SampleFreq!=-1) && (xmlparsedata->det->SampleFreq!=new_samplefreq)){
			SIXT_ERROR("Incompatible sampling frequency values encountered.");
			return;
		}
		xmlparsedata->det->SampleFreq=new_samplefreq;
        } else if (!strcmp(Uelement,"TESFILE")) {
        	xmlparsedata->det->tes_type_file=(char*)malloc(MAXFILENAME*sizeof(char));
            CHECK_MALLOC_VOID(xmlparsedata->det->tes_type_file);
            getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->tes_type_file);
            if (strlen(xmlparsedata->det->tes_type_file) == 0){
            	free(xmlparsedata->det->tes_type_file);
            	xmlparsedata->det->tes_type_file=NULL;
            }

	} else if(!strcmp(Uelement, "CHANNEL_FREQ_LIST"))  {
		xmlparsedata->det->channel_file=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(xmlparsedata->det->channel_file);
		getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->channel_file);

		xmlparsedata->det->channel_resfreq_file=(char*)malloc(MAXFILENAME*sizeof(char));
                CHECK_MALLOC_VOID(xmlparsedata->det->channel_resfreq_file);
		getXMLAttributeString(attr, "RESFILENAME", xmlparsedata->det->channel_resfreq_file);

        } else if(!strcmp(Uelement, "FDM_PARAMETERS"))  {
        	xmlparsedata->det->L_Common=getXMLAttributeDouble(attr, "LCOMMON");
                // retrieval of common capacitance
                // do it via getXMLAttributeString. If it's an empty string, do nothing (default value is -1)
                // if it's not an empty string, do atof, as is done in getXMLAttributeDouble
                char Cbuffer[MAXMSG]; // string buffer
                getXMLAttributeString(attr,"CCOMMON", Cbuffer);
                // an empty string starts with the terminator; this is faster than strcmp
                if (Cbuffer[0] != '\0'){
                  xmlparsedata->det->C_Common = atof(Cbuffer);
                }

	} else if(!strcmp(Uelement, "CROSSTALK"))  {
		// Need to check if we have channels defined
		if (xmlparsedata->det->channel_file == NULL){
			SIXT_ERROR("Trying to implement crosstalk, but no file defining the channels specified in the XML");
			return;
		}
	} else if(!strcmp(Uelement, "THERMALCROSSTALK"))  {
		xmlparsedata->det->xt_dist_thermal = realloc(xmlparsedata->det->xt_dist_thermal,
				(xmlparsedata->det->xt_num_thermal+1)*sizeof(double));
		CHECK_MALLOC_VOID(xmlparsedata->det->xt_dist_thermal);
		xmlparsedata->det->xt_dist_thermal[xmlparsedata->det->xt_num_thermal] =
				getXMLAttributeDouble(attr, "DISTANCE");

		xmlparsedata->det->crosstalk_thermal_timedep_file=(char**)realloc(xmlparsedata->det->crosstalk_thermal_timedep_file, (xmlparsedata->det->xt_num_thermal+1)*sizeof(char*));
		char* tmp=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp);
		char* tmp_path=(char*)malloc(2*MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp_path);
		getXMLAttributeString(attr, "TIMEDEPFILE", tmp);
		xmlparsedata->det->crosstalk_thermal_timedep_file[xmlparsedata->det->xt_num_thermal]=tmp;
		tmp_path=strcpy(tmp_path, xmlparsedata->det->filepath);
		tmp_path=strcat(tmp_path, tmp);
		if (strlen(tmp) == 0 || access(tmp_path, F_OK)==-1){
			xmlparsedata->status=EXIT_FAILURE;
			printf("Warning: Thermal crosstalk time dependence file %s for neighbor %i is not given or does not exist\n", tmp_path, xmlparsedata->det->xt_num_thermal+1);
			free(xmlparsedata->det->crosstalk_thermal_timedep_file[xmlparsedata->det->xt_num_thermal]);
			xmlparsedata->det->crosstalk_thermal_timedep_file[xmlparsedata->det->xt_num_thermal]=NULL;
			SIXT_ERROR("Themal crosstalk time dependence file not found\n");
		}
		free(tmp_path);

		// check that distances are in increasing order
		if((xmlparsedata->det->xt_num_thermal>0) &&
				(xmlparsedata->det->xt_dist_thermal[xmlparsedata->det->xt_num_thermal]<xmlparsedata->det->xt_dist_thermal[xmlparsedata->det->xt_num_thermal-1])){
			xmlparsedata->status=EXIT_FAILURE;
			free(xmlparsedata->det->crosstalk_thermal_timedep_file);
			SIXT_ERROR("Thermal crosstalk distances are supposed to be given in increasing order\n");
			return;
		}

		xmlparsedata->det->xt_weight_thermal = realloc(xmlparsedata->det->xt_weight_thermal,
				(xmlparsedata->det->xt_num_thermal+1)*sizeof(double) );
		CHECK_MALLOC_VOID(xmlparsedata->det->xt_weight_thermal);
		xmlparsedata->det->xt_weight_thermal[xmlparsedata->det->xt_num_thermal] =
				getXMLAttributeDouble(attr, "WEIGHT");

		xmlparsedata->det->xt_num_thermal++;

	} else if(!strcmp(Uelement, "ELECTRICALCROSSTALK")){
		xmlparsedata->det->crosstalk_elec_file=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(xmlparsedata->det->crosstalk_elec_file);
		char* tmp_path=(char*)malloc(2*MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp_path);
		tmp_path=strcpy(tmp_path, xmlparsedata->det->filepath);
		getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->crosstalk_elec_file);
		tmp_path=strcat(tmp_path, xmlparsedata->det->crosstalk_elec_file);
		if (strlen(xmlparsedata->det->crosstalk_elec_file) == 0 || access(tmp_path, F_OK)==-1){
			xmlparsedata->status=EXIT_FAILURE;
			printf("Warning: Electrical crosstalk file %s is not given or does not exist\n", tmp_path);
			free(xmlparsedata->det->crosstalk_elec_file);
			xmlparsedata->det->crosstalk_elec_file=NULL;
			SIXT_ERROR("Electrical crosstalk file not found");
		}
		free(tmp_path);

		xmlparsedata->det->crosstalk_elec_timedep_file=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(xmlparsedata->det->crosstalk_elec_timedep_file);
		char* tmp_path_2=(char*)malloc(2*MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp_path_2);
		tmp_path_2=strcpy(tmp_path_2, xmlparsedata->det->filepath);
		getXMLAttributeString(attr, "TIMEDEPFILE", xmlparsedata->det->crosstalk_elec_timedep_file);
		tmp_path_2=strcat(tmp_path_2, xmlparsedata->det->crosstalk_elec_timedep_file);
		if (strlen(xmlparsedata->det->crosstalk_elec_timedep_file) == 0 || access(tmp_path_2, F_OK)==-1){
			xmlparsedata->status=EXIT_FAILURE;
			printf("Warning: Electrical cross-talk time dependence file %s is not given or does not exist\n", tmp_path_2);
			free(xmlparsedata->det->crosstalk_elec_timedep_file);
			xmlparsedata->det->crosstalk_elec_timedep_file=NULL;
			SIXT_ERROR("Electrical cross-talk time dependence file not found");
		}
		free(tmp_path_2);

	} else if(!strcmp(Uelement, "INTERMODULATION"))  {
		xmlparsedata->det->crosstalk_intermod_file=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(xmlparsedata->det->crosstalk_intermod_file);
		char* tmp_path=(char*)malloc(2*MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp_path);
		tmp_path=strcpy(tmp_path, xmlparsedata->det->filepath);
		getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->crosstalk_intermod_file);
		tmp_path=strcat(tmp_path,xmlparsedata->det->crosstalk_intermod_file);
		if (strlen(xmlparsedata->det->crosstalk_intermod_file) == 0 || access(tmp_path, F_OK)==-1){
			xmlparsedata->status=EXIT_FAILURE;
			printf("Warning: Intermodulation cross-talk file %s is not given or does not exist\n", tmp_path);
			free(xmlparsedata->det->crosstalk_intermod_file);
			xmlparsedata->det->crosstalk_intermod_file=NULL;
			SIXT_ERROR("Intermodulation cross-talk file not found");
		}
		free(tmp_path);

	} else if(!strcmp(Uelement, "PROPCROSSTALK")){
		xmlparsedata->det->TDM_prop_file=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(xmlparsedata->det->TDM_prop_file);
		char* tmp_path=(char*)malloc(2*MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp_path);
		tmp_path=strcpy(tmp_path, xmlparsedata->det->filepath);
		getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->TDM_prop_file);
		tmp_path=strcat(tmp_path, xmlparsedata->det->TDM_prop_file);
		if (strlen(xmlparsedata->det->TDM_prop_file) == 0 || access(tmp_path, F_OK)==-1){
			xmlparsedata->status=EXIT_FAILURE;
			printf("Warning: Proportional cross-talk file %s is not given or does not exist\n", tmp_path);
			free(xmlparsedata->det->TDM_prop_file);
			xmlparsedata->det->TDM_prop_file=NULL;
			SIXT_ERROR("Proportional cross-talk file not found");
		}
		free(tmp_path);

		xmlparsedata->det->prop_TDM_scaling_1=getXMLAttributeDouble(attr,"SCALING1");
		xmlparsedata->det->prop_TDM_scaling_2=getXMLAttributeDouble(attr,"SCALING2");

	} else if(!strcmp(Uelement, "DERCROSSTALK")){
		xmlparsedata->det->TDM_der_file=(char*)malloc(MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(xmlparsedata->det->TDM_der_file);
		char* tmp_path=(char*)malloc(2*MAXFILENAME*sizeof(char));
		CHECK_MALLOC_VOID(tmp_path);
		tmp_path=strcpy(tmp_path, xmlparsedata->det->filepath);
		getXMLAttributeString(attr, "FILENAME", xmlparsedata->det->TDM_der_file);
		tmp_path=strcat(tmp_path, xmlparsedata->det->TDM_der_file);
		if (strlen(xmlparsedata->det->TDM_der_file) == 0 || access(tmp_path, F_OK)==-1){
			xmlparsedata->status=EXIT_FAILURE;
			printf("Warning: Derivative cross-talk file %s is not given or does not exist\n", tmp_path);
			free(xmlparsedata->det->TDM_der_file);
			xmlparsedata->det->TDM_der_file=NULL;
			SIXT_ERROR("Derivative cross-talk file not found");
		}
		free(tmp_path);

		xmlparsedata->det->der_TDM_scaling=getXMLAttributeDouble(attr,"SCALING");

	} else if (!strcmp(Uelement,"THRESHOLD_EVENT_LO_KEV")){
		xmlparsedata->det->threshold_event_lo_keV = getXMLAttributeDouble(attr,"VALUE");
    } else {
		// Unknown tag, display warning.
		char msg[MAXMSG];
		sprintf(msg, "unknown XML tag: <%s>", el);
		SIXT_WARNING(msg);
	}
}

static void AdvDetXMLElementEnd(void* parsedata, const char* el)
{
	struct XMLParseData* xmlparsedata=(struct XMLParseData*)parsedata;

	// Check if an error has occurred previously.
	CHECK_STATUS_VOID(xmlparsedata->status);

	// Convert the element to an upper case string.
	char Uelement[MAXMSG];
	strcpy(Uelement, el);
	strtoupper(Uelement);

	if (!strcmp(Uelement, "TESNOISEFILTER")) {
		xmlparsedata->det->tesnoisefilter=0;
	} else if(!strcmp(Uelement, "PIXEL")){
		xmlparsedata->det->inpixel=0;
		xmlparsedata->det->cpix++;
	}

	// Check if an error has occurred previously.
	CHECK_STATUS_VOID(xmlparsedata->status);
}

AdvDet* loadAdvDet(const char* const filename,
		     int* const status)
{
  // Get a new and empty data structure.
  AdvDet* det=newAdvDet(status);
  CHECK_STATUS_RET(*status, det);

  // Split the reference to the XML detector definition file
  // into path and filename. This has to be done before
  // calling the parser routine for the XML file.
  char filename2[MAXFILENAME];
  char rootname[MAXFILENAME];
  // Make a local copy of the filename variable in order to avoid
  // compiler warnings due to discarded const qualifier at the
  // subsequent function call.
  strcpy(filename2, filename);
  fits_parse_rootname(filename2, rootname, status);
  CHECK_STATUS_RET(*status, det);

  // Split rootname into the file path and the file name.
  char* lastslash=strrchr(rootname, '/');
  if (NULL==lastslash) {
    det->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(det->filepath, *status,
		   "memory allocation for filepath failed", det);
    det->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(det->filename, *status,
		   "memory allocation for filename failed", det);
    strcpy(det->filepath, "");
    strcpy(det->filename, rootname);
  } else {
    lastslash++;
    det->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(det->filename, *status,
		   "memory allocation for filename failed", det);
    strcpy(det->filename, lastslash);

    *lastslash='\0';
    det->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(det->filepath, *status,
		   "memory allocation for filepath failed", det);
    strcpy(det->filepath, rootname);
  }
  // END of storing the filename and filepath.

  // Read in the XML definition of the detector.
  parseAdvDetXML(det, filename, status);
  CHECK_STATUS_RET(*status, det);

  // Remove overlapping pixels with the rule newest survives
  removeOverlapping(det,status);

  return(det);
}


/** Iterates the different pixels and loads the necessary RMFLibrary */
void loadRMFLibrary(AdvDet* det, int* const status){
	det->rmf_library = malloc(sizeof(*(det->rmf_library)));
	if (NULL == det->rmf_library){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation for rmf library failed");
		return;
	}
	det->rmf_library->rmf_array = malloc(RMFLIBRARYSIZE*sizeof(*(det->rmf_library->rmf_array)));
	if (NULL == det->rmf_library->rmf_array){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation for rmf library failed");
		return;
	}
	det->rmf_library->filenames = malloc(RMFLIBRARYSIZE*sizeof(*(det->rmf_library->filenames)));
	if (NULL == det->rmf_library->filenames){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation for rmf library failed");
		return;
	}

	det->rmf_library->size = RMFLIBRARYSIZE;
	det->rmf_library->n_rmf = 0;

	for (int i=0;i<RMFLIBRARYSIZE;i++){
		det->rmf_library->rmf_array[i]=NULL;
		det->rmf_library->filenames[i]=NULL;
	}

	for (int i=0;i<det->npix;i++){
		for (int j=0;j<det->pix[i].ngrades;j++){
			addRMF(det,&(det->pix[i]),j,status);
			CHECK_STATUS_VOID(*status);
		}
	}
}

/** Adds the RMF corresponding to the jth index in the pixel rmffiles array to the RMF library. The RMF will only be added if it is not already in the library */
void addRMF(AdvDet* det,AdvPix* pixel,int rmf_index,int* const status){
	if(NULL==pixel->grades){
		*status=EXIT_FAILURE;
		SIXT_ERROR("Tried to load pixel RMF whereas no RMF file was given in XML. Abort");
		return;
	}
	for (int i=0;i<det->rmf_library->n_rmf;i++){
		if(!strcmp(det->rmf_library->filenames[i],pixel->grades[rmf_index].rmffile)){
			pixel->grades[rmf_index].rmf=det->rmf_library->rmf_array[i];
			return; //If the rmf is already in there, just update the rmfID and return
		}
	}

	// Update size if necessary
	if (det->rmf_library->n_rmf>=det->rmf_library->size){
	  det->rmf_library->size *=2;
	  struct RMF** new_rmf_array = realloc(det->rmf_library->rmf_array,det->rmf_library->size*sizeof(*(det->rmf_library->rmf_array)));
	  if (NULL==new_rmf_array){
	    *status = EXIT_FAILURE;
	    SIXT_ERROR("Size update of RMF library failed");
	    return;
	  }
	  char** new_filenames = realloc(det->rmf_library->filenames,det->rmf_library->size*sizeof(*(det->rmf_library->filenames)));
	  if (NULL==new_filenames){
	    *status = EXIT_FAILURE;
	    SIXT_ERROR("Size update of RMF library failed");
	    return;
	  }

	  det->rmf_library->rmf_array=new_rmf_array;
	  det->rmf_library->filenames=new_filenames;
	}

	//Add RMF to the library
	char filepathname[MAXFILENAME];
	strncpy(filepathname,det->filepath,MAXFILENAME);
	int maxch=strlen(filepathname);
	if (maxch<MAXFILENAME) {
	  strncat(filepathname,pixel->grades[rmf_index].rmffile,MAXFILENAME-maxch);
	  filepathname[MAXFILENAME-1]='\0';
	} else {
	  *status=EXIT_FAILURE;
	  SIXT_ERROR("RMF path and filename are too long");
	  return;
	}
	det->rmf_library->rmf_array[det->rmf_library->n_rmf] = loadNormalizedRMF(filepathname,status);
	det->rmf_library->filenames[det->rmf_library->n_rmf] = strndup(pixel->grades[rmf_index].rmffile,MAXFILENAME);
	pixel->grades[rmf_index].rmf=det->rmf_library->rmf_array[det->rmf_library->n_rmf];
	det->rmf_library->n_rmf++;
}

/** Destructor of the RMF library structure */
void freeRMFLibrary(RMFLibrary* library){
	if (NULL!=library){
		for(int i=0;i<library->size;i++){
			freeRMF(library->rmf_array[i]);
			free(library->filenames[i]);
		}
		free(library->rmf_array);
		free(library->filenames);
		free(library);
	}
	library=NULL;
}

/** Iterates the different pixels and loads the necessary ARFLibrary */
void loadARFLibrary(AdvDet* det, int* const status){
	det->arf_library = malloc(sizeof(*(det->arf_library)));
	if (NULL == det->arf_library){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation for arf library failed");
		return;
	}
	det->arf_library->arf_array = malloc(RMFLIBRARYSIZE*sizeof(*(det->arf_library->arf_array)));
	if (NULL == det->arf_library->arf_array){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation for arf library failed");
		return;
	}
	det->arf_library->filenames = malloc(RMFLIBRARYSIZE*sizeof(*(det->arf_library->filenames)));
	if (NULL == det->arf_library->filenames){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation for arf library failed");
		return;
	}

	det->arf_library->size = RMFLIBRARYSIZE;
	det->arf_library->n_arf = 0;

	for (int i=0;i<RMFLIBRARYSIZE;i++){
		det->arf_library->arf_array[i]=NULL;
		det->arf_library->filenames[i]=NULL;
	}

	for (int i=0;i<det->npix;i++){
		addARF(det,&(det->pix[i]),status);
		CHECK_STATUS_VOID(*status);
	}
}

/** Adds an ARF to the ARF library. The ARF will only be added if it is not already in the library */
void addARF(AdvDet* det,AdvPix* pixel,int* const status){
	if(NULL==pixel->arffile){
		*status=EXIT_FAILURE;
		SIXT_ERROR("Tried to load pixel ARF whereas no ARF file was given in XML. Abort");
		return;
	}
	for (int i=0;i<det->arf_library->n_arf;i++){
		if(!strcmp(det->arf_library->filenames[i],pixel->arffile)){
			pixel->arf=det->arf_library->arf_array[i];
			return; //If the arf is already in there, just update the arfID and return
		}
	}

	// Update size if necessary
	if (det->arf_library->n_arf>=det->arf_library->size){
		det->arf_library->size = det->arf_library->size*2;
		struct ARF** new_arf_array = realloc(det->arf_library->arf_array,det->arf_library->size*sizeof(*(det->arf_library->arf_array)));
		if (NULL==new_arf_array){
			*status = EXIT_FAILURE;
			SIXT_ERROR("Size update of ARF library failed");
			return;
		}
		char** new_filenames = realloc(det->arf_library->filenames,det->arf_library->size*sizeof(*(det->arf_library->filenames)));
		if (NULL==new_filenames){
			*status = EXIT_FAILURE;
			SIXT_ERROR("Size update of ARF library failed");
			return;
		}

		det->arf_library->arf_array=new_arf_array;
		det->arf_library->filenames=new_filenames;
	}

	//Add ARF to the library
	char filepathname[MAXFILENAME];
	strcpy(filepathname,det->filepath);
	strcat(filepathname,pixel->arffile);
	det->arf_library->arf_array[det->arf_library->n_arf] = loadARF(filepathname,status);
	det->arf_library->filenames[det->arf_library->n_arf] = strdup(pixel->arffile);
	pixel->arf=det->arf_library->arf_array[det->arf_library->n_arf];
	det->arf_library->n_arf++;
}

/** Destructor of the ARF library structure */
void freeARFLibrary(ARFLibrary* library){
	if (NULL!=library){
		for(int i=0;i<library->size;i++){
			freeARF(library->arf_array[i]);
			free(library->filenames[i]);
		}
		free(library->arf_array);
		free(library->filenames);
		free(library);
	}
	library=NULL;
}

/** Function to remove overlapping pixels from the detector */
void removeOverlapping(AdvDet* det,int* const status){
	char * active_pixels = malloc(det->npix*sizeof(*active_pixels));
	if (NULL==active_pixels){
		*status = EXIT_FAILURE;
		SIXT_ERROR("Memory allocation failed for active_pixels in removeOverlapping");
		return;
	}

	AdvPix* current_pixel=NULL;
	AdvPix* pixel_to_compare=NULL;
	int number_active_pixels=0;
	for (int i=0;i<det->npix;i++){
		active_pixels[i]=1;
		current_pixel = &(det->pix[i]);
		number_active_pixels++;
		for (int j=0;j<i;j++){
			pixel_to_compare = &(det->pix[j]);
			if(!(current_pixel->sx-.5*current_pixel->width > pixel_to_compare->sx + .5*pixel_to_compare->width || current_pixel->sx+.5*current_pixel->width < pixel_to_compare->sx-.5*pixel_to_compare->width ||
					current_pixel->sy -.5*current_pixel->height > pixel_to_compare->sy+.5*pixel_to_compare->height || current_pixel->sy+.5*current_pixel->height<pixel_to_compare->sy-.5*pixel_to_compare->height)){
				if (active_pixels[j]) number_active_pixels--;
				active_pixels[j]=0;
			}
		}
	}

	AdvPix* new_pix_array=malloc(number_active_pixels*sizeof(*(new_pix_array)));
	if(NULL==new_pix_array){
		*status=EXIT_FAILURE;
		SIXT_ERROR("Memory allocation failed for new pixel array in removeOverlapping");
		return;
	}

	det->cpix=0;
	for (int i=0;i<det->npix;i++){
		if(active_pixels[i]){
			new_pix_array[det->cpix]=det->pix[i];
			new_pix_array[det->cpix].pindex=det->cpix;
			det->cpix++;
		} else{
			freeAdvPix(&(det->pix[i]));
		}
	}
	free(det->pix);
	det->pix = new_pix_array;
	det->npix=number_active_pixels;
	free(active_pixels);
	headas_chat(0,"Number of pixels after removing overlaps: %d\n \n",number_active_pixels);
}

/** Constructor for MatrixCrossTalk structure */
IntermodulationCrossTalk* newImodCrossTalk(int* const status){
	IntermodulationCrossTalk* matrix = (IntermodulationCrossTalk*) malloc(sizeof(IntermodulationCrossTalk));
	CHECK_MALLOC_RET_NULL_STATUS(matrix,*status);

	matrix->num_cross_talk_pixels=0;

	matrix->cross_talk_pixels = NULL;
	matrix->cross_talk_info = NULL;

	return matrix;
}

/** Destructor for MatrixCrossTalk structure */
void freeMatrixCrossTalk(MatrixCrossTalk* matrix){
	if (matrix!=NULL){
		free(matrix->cross_talk_pixels);
		free(matrix->cross_talk_weights);
		free(matrix->cross_talk_index);
	}
	free(matrix);
	matrix=NULL;
}

/** Destructor for MatrixCrossTalk structure */
void freeMatrixEnerdepCrossTalk(MatrixEnerdepCrossTalk* matrix, int gr){
	if (matrix!=NULL){
		for(int k=0;k<gr;k++){
			if (matrix[k].cross_talk_weights!=NULL ){
				for (int ii=0; ii<matrix[k].num_cross_talk_pixels; ii++){
					free(matrix[k].cross_talk_weights[ii]);
				}
				free(matrix[k].cross_talk_weights);
			}
			free(matrix[k].cross_talk_pixels);
		}
	}
	free(matrix);
	matrix=NULL;
}

/** Destructor for MatrixPropCrossTalk structure */
void freeMatrixPropCrossTalk(MatrixPropCrossTalk* matrix){
	if (matrix!=NULL){
		free(matrix->cross_talk_pixels_1);
		free(matrix->cross_talk_pixels_2);
		matrix->type_1_pix=0;
		matrix->type_2_pix=0;
	}
	free(matrix);
	matrix=NULL;
}

/** Destructor for MatrixDerCrossTalk structure */
void freeMatrixDerCrossTalk(MatrixDerCrossTalk* matrix){
	if (matrix!=NULL){
		free(matrix->cross_talk_pixels);
		matrix->num_cross_talk_pixels=0;
	}
	free(matrix);
	matrix=NULL;
}

/** Destructor for MatrixCrossTalk structure */
void freeImodCrossTalk(IntermodulationCrossTalk** matrix){
	if (*matrix!=NULL){
		for (int ii=0; ii <  (*matrix)->num_cross_talk_pixels; ii++){
			free((*matrix)->cross_talk_pixels[ii]);
		}
		free((*matrix)->cross_talk_info);
/*		for (int ii=0; ii <  (*matrix)->num_cross_talk_pixels; ii++){
			free((*matrix)->cross_talk_pixels[ii]);
			free((*matrix)->cross_talk_info[ii]);
		}
		free((*matrix)->num_pixel_combinations); */
	}
	free(*matrix);
}


/** Constructor for MatrixCrossTalk structure */
MatrixCrossTalk* newMatrixCrossTalk(int* const status){
	MatrixCrossTalk* matrix = (MatrixCrossTalk*) malloc(sizeof(*matrix));
	CHECK_MALLOC_RET_NULL_STATUS(matrix,*status);

	matrix->num_cross_talk_pixels=0;
	matrix->cross_talk_index = NULL;
	matrix->cross_talk_pixels = NULL;
	matrix->cross_talk_weights = NULL;

	return matrix;
}

/** Constructor for MatrixEnerdepCrossTalk structure */
MatrixEnerdepCrossTalk* newMatrixEnerdepCrossTalk(int grade, int* const status){
	MatrixEnerdepCrossTalk* matrix = (MatrixEnerdepCrossTalk*) malloc(grade*sizeof(MatrixEnerdepCrossTalk));
	CHECK_MALLOC_RET_NULL_STATUS(matrix,*status);

	for(int l=0;l<grade;l++){
		matrix[l].num_cross_talk_pixels=0;
		matrix[l].cross_talk_pixels = NULL;
		matrix[l].cross_talk_weights = NULL;
		matrix[l].n_ener=0;
	}
	return matrix;
}

/** Constructor for MatrixPropCrossTalk structure */
MatrixPropCrossTalk* newMatrixPropCrossTalk(int* const status){
	MatrixPropCrossTalk* matrix = (MatrixPropCrossTalk*) malloc(sizeof(MatrixPropCrossTalk));
	CHECK_MALLOC_RET_NULL_STATUS(matrix,*status);

	matrix->type_1_pix=0;
	matrix->type_2_pix=0;
	matrix->cross_talk_pixels_1 = NULL;
	matrix->cross_talk_pixels_2 = NULL;

	return matrix;
}

/** Constructor for MatrixPropCrossTalk structure */
MatrixDerCrossTalk* newMatrixDerCrossTalk(int* const status){
	MatrixDerCrossTalk* matrix = (MatrixDerCrossTalk*) malloc(sizeof(MatrixPropCrossTalk));
	CHECK_MALLOC_RET_NULL_STATUS(matrix,*status);

	matrix->num_cross_talk_pixels=0;
	matrix->cross_talk_pixels = NULL;

	return matrix;
}

/** Constructor for CrosstalkTimdep structure */
CrosstalkTimedep* newCrossTalkTimedep(int* const status){
	CrosstalkTimedep* crosstalk_timedep = (CrosstalkTimedep*)malloc(sizeof(CrosstalkTimedep));
	CHECK_MALLOC_RET_NULL_STATUS(crosstalk_timedep,*status);
	crosstalk_timedep->length=0;
	crosstalk_timedep->name_type=NULL; // Useless for the moment
	crosstalk_timedep->time=NULL;
	crosstalk_timedep->weight=NULL;
	crosstalk_timedep->weight_t0 = 1.0;
	return crosstalk_timedep;
}

/** Destructor for CrosstalkTimdep structure */
void freeCrosstalkTimedep(CrosstalkTimedep* timedep, int gr){
	if(timedep!=NULL){
		for(int k=0; k<2*gr;k++){
			timedep[k].length=0;
			free(timedep[k].name_type);
			free(timedep[k].time);
			free(timedep[k].weight);
			timedep[k].weight_t0=0;
		}
	}
}
