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


   Copyright 2015 Philippe Peille, IRAP
*/

#include "tesreconstruction.h"

////////////////////////////////////
/** Main procedure. */
int tesreconstruction_main() {
  time_t ttstart = time(0);
  
  // Containing all programm parameters read by PIL.
  struct Parameters par;
  
  // Error status.
  int status=EXIT_SUCCESS;

  // Register HEATOOL:
  set_toolname("tesreconstruction");
  set_toolversion("0.05");
  
  do { // Beginning of the ERROR handling loop (will at
       // most be run once).
    headas_chat(3, "initialize ...\n");
    // Get program parameters.
    status=getpar(&par);
    CHECK_STATUS_BREAK(status);
    
    // Read XML info
    //--------------
    AdvDet *det = newAdvDet(&status);
    CHECK_STATUS_BREAK(status);
    det = loadAdvDet(par.XMLFile, &status);
    CHECK_STATUS_BREAK(status);
    // Read the sampling rate from XML file
    double sf = -999.; 
    double div = 1;
    sf = det->SampleFreq;
    
    double sampling_rate;
    //int trig_reclength;
    
    char* firstchar = strndup(par.RecordFile, 1);
    char firstchar2[2];
    strcpy(firstchar2,firstchar);
        
    //printf("%s %s %s","File: ",par.RecordFile,"\n");
    // Check input file header is complete to work with xifusim/tessim simulated files
    // -------------------------------------------------------------------------------
    fitsfile* fptr = NULL;
    int numfits;
    int hdunum; // Number of HDUs (RECORDS-file or TESRECORDS-file)
    if (strcmp(firstchar2,"@") == 0)
    {
            //printf("%s %s %s","File: ",strndup(par.RecordFile+1, strlen(par.RecordFile)-1),"\n");
            FILE *filetxt = fopen(strndup(par.RecordFile+1, strlen(par.RecordFile)-1), "r");
            if (filetxt == NULL)    
            {
                    printf("%s","File given in RecordFile does not exist\n");
                    status = 104;
            }
            CHECK_STATUS_BREAK(status);
            
            char filefits[256];
            numfits = 0;
            while(fscanf(filetxt,"%s",filefits)!=EOF)
            {
                numfits++;
            }
            fclose(filetxt);
            
            filetxt = fopen(strndup(par.RecordFile+1, strlen(par.RecordFile)-1), "r");
        
            for (int j=0;j<numfits;j++)   // For every FITS file
            {
                    fgets(filefits, 256, filetxt);
                    strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)
            
                    fits_open_file(&fptr, filefits, READONLY, &status);
                    if (status != 0)    printf("%s","FITS file read from ASCII file does not exist\n");
                    CHECK_STATUS_BREAK(status);  
                    
                    fits_get_num_hdus(fptr, &hdunum,&status);
    
                    if ((hdunum == 8) || (hdunum == 9)) //xifusim simulated file (with TESRECORDS)
                    {    
                        
                        // Move to "Primary" HDU to obtain SAMPLING_RATE
                        fits_movabs_hdu(fptr, 1, NULL, &status); 
                        CHECK_STATUS_BREAK(status);
                        // and read full Primary HDU and store it in 'headerPrimary'
                        int numberkeywords;
                        char *headerPrimary;
                        fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status); 
                        CHECK_STATUS_BREAK(status);
                            
                        // Pointer to where the text "sample_rate=" is in HISTORY block
                        sampling_rate = -999.0;
                        char * sample_rate_pointer;
                        sample_rate_pointer = strstr (headerPrimary,"sample_rate=");    
                        if(!sample_rate_pointer)
                        {
                            // read it from xml file
                            sampling_rate = sf;
                        }
                        else
                        {
                            // Pointer to the next character to "sample_rate=" (12 characters)   
                            sample_rate_pointer = sample_rate_pointer + 12; 
                            char each_character_after_srate[125];		
                            snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                              
                            char characters_after_srate[125];
                            snprintf(characters_after_srate,125,"%c",*sample_rate_pointer);
                                
                            while (*sample_rate_pointer != ' ')
                            {
                                sample_rate_pointer = sample_rate_pointer + 1;
                                snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                                strcat(characters_after_srate,each_character_after_srate); 
                            }
                                
                            sampling_rate = atof(characters_after_srate);
                        }
                            
                        div = sf/sampling_rate;  // Grading info is unique in XML file -> adjust for different sf
            
                    }//if ((hdunum == 8) || (hdunum == 9)) (xifusim file)
                    else
                    {
                        fits_movnam_hdu(fptr, ANY_HDU,"TESRECORDS", 0, &status);
                        CHECK_STATUS_BREAK(status);
                        double keyvalue_double;
                        fits_read_key(fptr,TDOUBLE,"DELTAT",&keyvalue_double,NULL,&status);
                        sampling_rate = 1/keyvalue_double;
                        div = sf/sampling_rate;  // Grading info is unique in XML file -> adjust for different sf
                    }//(tessim)
                    fits_close_file(fptr,&status);
                    CHECK_STATUS_BREAK(status);
            }
            
            fclose(filetxt);
    }
    else
    {
            numfits = 1;
            fits_open_file(&fptr, par.RecordFile, READONLY, &status);
            if (status != 0)    printf("%s","File given in RecordFile does not exist\n");
            
            fits_get_num_hdus(fptr, &hdunum,&status);
    
            if ((hdunum == 8) || (hdunum == 9)) //xifusim simulated file (with TESRECORDS)
            {    
                // Move to "Primary" HDU to obtain SAMPLING_RATE
                fits_movabs_hdu(fptr, 1, NULL, &status); 
                CHECK_STATUS_BREAK(status);
                // and read full Primary HDU and store it in 'headerPrimary'
                int numberkeywords;
                char *headerPrimary;
                fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status); 
                CHECK_STATUS_BREAK(status);
                
                // Pointer to where the text "sample_rate=" is in HISTORY block
                sampling_rate = -999.0;
                char * sample_rate_pointer;
                sample_rate_pointer = strstr (headerPrimary,"sample_rate=");    
                if(!sample_rate_pointer)
                {
                    // read it from xml file
                    sampling_rate = sf;
                }
                else
                {
                    // Pointer to the next character to "sample_rate=" (12 characters)   
                    sample_rate_pointer = sample_rate_pointer + 12; 
                    char each_character_after_srate[125];		
                    snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                        
                    char characters_after_srate[125];
                    snprintf(characters_after_srate,125,"%c",*sample_rate_pointer);
                        
                    while (*sample_rate_pointer != ' ')
                    {
                        sample_rate_pointer = sample_rate_pointer + 1;
                        snprintf(each_character_after_srate,125,"%c",*sample_rate_pointer);
                        strcat(characters_after_srate,each_character_after_srate); 
                    }
                       
                    sampling_rate = atof(characters_after_srate);
                }
                
                div = sf/sampling_rate;  // Grading info is unique in XML file -> adjust for different sf
                
                /*// Pointer to where the text "sample_rate=" is in HISTORY block
                trig_reclength = -999.0;
                char * trig_reclength_pointer;
                trig_reclength_pointer = strstr (headerPrimary,"trig_reclength=");    
                
                // Pointer to the next character to "trig_reclength=" (15 characters)   
                trig_reclength_pointer = trig_reclength_pointer + 15; 
                char each_character_after_treclength[125];		
                snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                
                char characters_after_treclength[125];
                snprintf(characters_after_treclength,125,"%c",*trig_reclength_pointer);
                
                while (*trig_reclength_pointer != ' ')
                {
                    trig_reclength_pointer = trig_reclength_pointer + 1;
                    snprintf(each_character_after_treclength,125,"%c",*trig_reclength_pointer);
                    strcat(characters_after_treclength,each_character_after_treclength); 
                }
                
                trig_reclength = atoi(characters_after_treclength);
                printf("%s %d %s","trig_reclength: ",trig_reclength,"\n");*/
                
            }//if hdunum==8 (xifusim file)
            else //if hdunum!=8 (sixtefile)
            {
                fits_movnam_hdu(fptr, ANY_HDU,"TESRECORDS", 0, &status);
                if (status != 0)
                {
                    status = 0;
                    fits_movnam_hdu(fptr, ANY_HDU,"RECORDS", 0, &status);
                }
                CHECK_STATUS_BREAK(status);
                double keyvalue_double;
                fits_read_key(fptr,TDOUBLE,"DELTAT",&keyvalue_double,NULL,&status);
                sampling_rate = 1./keyvalue_double;
                div = sf/sampling_rate;  // Grading info is unique in XML file -> adjust for different sf
            } // (tessim file)
            fits_close_file(fptr,&status);
            CHECK_STATUS_BREAK(status);
    }

    // Sixt standard keywords structure
    //----------------------------------
    SixtStdKeywords* keywords = newSixtStdKeywords(&status);
    CHECK_STATUS_BREAK(status);
    
    //Open outfile 
    //------------
    TesEventFile * outfile = opennewTesEventFileSIRENA(par.TesEventFile,
                                                 keywords,
                                                 SIRENA_VERSION,
                                                 par.clobber,
                                                 &status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize PP data structures needed for pulse filtering
    //---------------------------------------------------------
    ReconstructInit* reconstruct_init = newReconstructInit(&status);
    CHECK_STATUS_BREAK(status);
    
    // Initialize SIRENA data structures needed for pulse filtering
    //-------------------------------------------------------------
    ReconstructInitSIRENA* reconstruct_init_sirena = newReconstructInitSIRENA(&status);
    CHECK_STATUS_BREAK(status);
    PulsesCollection* pulsesAll = newPulsesCollection(&status);
    CHECK_STATUS_BREAK(status);  
    OptimalFilterSIRENA* optimalFilter = newOptimalFilterSIRENA(&status);
    CHECK_STATUS_BREAK(status);// define a second structure for calibration
    
    // Read the grading data from the XML file and store it in 'reconstruct_init_sirena->grading'
    reconstruct_init_sirena->grading = NULL;
    reconstruct_init_sirena->grading = (Grading*)malloc(sizeof(Grading));
    
    reconstruct_init_sirena->grading->ngrades = 0;
    reconstruct_init_sirena->grading->value  = NULL;
    reconstruct_init_sirena->grading->gradeData = NULL;
    
    if (det->pix->grades == NULL)
    {
        SIXT_ERROR("The provided XMLFile does not have the grading info");
        return(EXIT_FAILURE);
    }
    reconstruct_init_sirena->grading->ngrades=det->pix->ngrades;
    reconstruct_init_sirena->grading->gradeData = gsl_matrix_alloc(det->pix->ngrades,2);
    for (int i=0;i<det->pix->ngrades;i++)
    {
        gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,0,(int) (det->pix->grades[i].gradelim_pre)/div);
        gsl_matrix_set(reconstruct_init_sirena->grading->gradeData,i,1,(int) (det->pix->grades[i].gradelim_post)/div);
    }
    destroyAdvDet(&det);
    
    // Build up TesEventList to recover the results of the reconstruction
    TesEventList* event_list = newTesEventList(&status);
    allocateTesEventListTrigger(event_list,par.EventListSize,&status);
    CHECK_STATUS_BREAK(status);
            
    TesTriggerFile* record_file;
    TesRecord* record;
    int lastRecord = 0, nrecord = 0, nrecord_filei = 0;    //last record required for SIRENA library creation
    
    if (strcmp(firstchar2,"@") == 0)
    {
            FILE *filetxt = fopen(strndup(par.RecordFile+1, strlen(par.RecordFile)-1), "r");
            if (status != 0)    printf("%s","FITS file read from ASCII file does not exist\n");
            CHECK_STATUS_BREAK(status);  
            
            char filefits[256];
        
            for (int j=0;j<numfits;j++)   // For every FITS file
            {
                    fgets(filefits, 256, filetxt);
                    strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)
                    //printf("%s %s %s","FITS file i: ",filefits,"\n");
            
                    // Open record file
                    // ----------------
                    record_file = openexistingTesTriggerFile(filefits,keywords,&status);
                    CHECK_STATUS_BREAK(status);
                    
                    if(!strcmp(par.Rcmethod,"PP"))
                    {
                            initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.PulseLength,
                                                    par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
                                                    par.DerivateExclusion,par.SaturationValue,&status);
                    }
                    else
                    {
                            initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr, 
                                                    par.LibraryFile, par.TesEventFile, par.PulseLength, par.scaleFactor, par.samplesUp, 
                                                    par.samplesDown, par.nSgms, par.detectSP, par.opmode, par.detectionMode, par.LrsT, 
                                                    par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod, par.EnergyMethod, 
                                                    par.filtEev, par.OFNoise, par.LagsOrNot, par.nLags, par.Fitting35, par.OFIter, 
                                                    par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength, par.preBuffer,par.monoenergy, 
                                                    par.hduPRECALWN, par.hduPRCLOFWM, par.largeFilter, par.intermediate, par.detectFile, 
                                                    par.filterFile, par.errorT, par.Sum0Filt, par.clobber, par.EventListSize, par.SaturationValue, par.tstartPulse1, 
                                                    par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);
                                            
                    }  
                    CHECK_STATUS_BREAK(status);
                    
                    // Build up TesRecord to read the file
                    record = newTesRecord(&status);
                    if (record_file->delta_t == -999) record_file->delta_t = 1./sampling_rate;
                    allocateTesRecord(record,record_file->trigger_size,record_file->delta_t,0,&status);
                    CHECK_STATUS_BREAK(status);
                    
                    // Iterate of records and do the reconstruction
                    //int lastRecord = 0, nrecord = 0;    //last record required for SIRENA library creation
                    nrecord_filei = 0;
                    while(getNextRecord(record_file,record,&status))
                    {
                            if(!strcmp(par.Rcmethod,"PP"))
                            {
                                    reconstructRecord(record,event_list,reconstruct_init,0,&status);
                            }
                            else
                            {
                                    nrecord = nrecord + 1;
                                    nrecord_filei = nrecord_filei + 1;
                                    if ((nrecord_filei == record_file->nrows) && (j == numfits-1)) lastRecord=1;  // lastRecord of all the FITS files
                                   
                                    if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RALL") == 0) 
                                        || (strcmp(par.EnergyMethod,"I2RNOL") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0))
                                    {
                                            strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
                                    }
                                
                                    //printf("%s %d %s","**TESRECONSTRUCTION nrecord = ",nrecord,"\n");
                                    reconstructRecordSIRENA(record,event_list,reconstruct_init_sirena,
                                                            lastRecord, nrecord, &pulsesAll, &optimalFilter, &status);
                            }
                            CHECK_STATUS_BREAK(status);

                            if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
                            {
                                    // In THREADING mode, saveEventListToFile is not called until finishing with calculus 
                                    // (ordering is necessary previously)  
                                    if(!is_threading()){    
                                            //printf("\n %p - %f", outfile, record_file->delta_t);
                                            //printf("\nRecord single");
                                            //printf("\n%f - %ld", record->time, record->pixid);
                                            saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                                            CHECK_STATUS_BREAK(status);
                                            //Reinitialize event list
                                            event_list->index=0;//////////////////////!!!!!!!!!!!!!!!!!!!!!OJO!!!!!!!!!!!!!!!!!!! Igual no hay que inicializarlo a 0
                                    }
                                    //else
                                    //        printf("%s","Not prepared to run in THREADING mode with a input ASCII file (with several FITS files)\n");
                            }
                    } // while getNextRecord
                    if(is_threading()) 
                    {
                            //printf("%s","**Threading...waiting \n");
                            th_end(&reconstruct_init_sirena, &pulsesAll, &optimalFilter);
                            //printf("%s %d %s","**Threading...after th_end: pulsesAll->ndetpulses", pulsesAll->ndetpulses,"\n");
                            //printf("%s %d %s","**Threading...after th_end: pulsesAll->size", pulsesAll->size,"\n");
                            int i = 1;
                            int aux = 1;
                            while((aux = th_get_event_list(&event_list, &record)) == 1)
                            {
                                    //printf("%s %d %s","**Threading...i: ", i,"\n");
                                    //printf("%s %d %s","**Threading...event_list->size_energy: ", event_list->size_energy,"\n"); Always 0
                                    //printf("%s %e %s","**Threading...event_list->energies[0]: ", event_list->energies[0],"\n"); Energy value
                                    //printf("%s %e %s","**Threading...event_list->energies[1]: ", event_list->energies[1],"\n"); Not error but non relevant value
                                    //printf("%s %e %s","**Threading...event_list->energies[100000]: ", event_list->energies[100000],"\n"); Not error but non relevant value
                                    saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                                    //printf("%s","**Threading...after de saveEventListToFile \n");
                                    CHECK_STATUS_BREAK(status);
                                    ++i;
                            }
                    }
                    
                    if ((!strcmp(par.Rcmethod,"SIRENA")) && (pulsesAll->ndetpulses == 0)) 
                            printf("%s %s %s","WARNING: no pulses have been detected in the current FITS file: ", filefits,"\n");
                    
                    if (numfits == 0)   // Only one time
                    {
                            // Copy trigger keywords to event file
                            copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
                            CHECK_STATUS_BREAK(status);
                            
                            // Messages providing info of some columns
                            char keywordvalue[9];
                            char comment[MAXMSG];
                            
                            fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
                            CHECK_STATUS_BREAK(status);
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE1", &keywordvalue, NULL, &status);
                            strcpy(comment, "Starting time");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE1", keywordvalue, comment, &status);
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE2", &keywordvalue, NULL, &status);
                            strcpy(comment, "Reconstructed-uncalibrated energy");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE2", keywordvalue, comment, &status);      
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE3", &keywordvalue, NULL, &status);
                            strcpy(comment, "Average first 4 samples (derivative)");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE3", keywordvalue, comment, &status);      
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE4", &keywordvalue, NULL, &status);
                            strcpy(comment, "Optimal filter length");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE4", keywordvalue, comment, &status);      
                            
                            fits_read_key(outfile->fptr, TSTRING, "TTYPE5", &keywordvalue, NULL, &status);
                            strcpy(comment, "Starting time-starting time previous event");
                            fits_update_key(outfile->fptr, TSTRING, "TTYPE5", keywordvalue, comment, &status);
                    }
                    
                    freeTesTriggerFile(&record_file,&status);   // The record_file (every FITS file) is closed
                    
                    CHECK_STATUS_BREAK(status);
            
            }   // for every FITS file
            
            fclose(filetxt);
    }
    else
    {
            // Open record file
            // ----------------
            record_file = openexistingTesTriggerFile(par.RecordFile,keywords,&status);
            CHECK_STATUS_BREAK(status);

            if(!strcmp(par.Rcmethod,"PP")){
                initializeReconstruction(reconstruct_init,par.OptimalFilterFile,par.PulseLength,
                        par.PulseTemplateFile,par.Threshold,par.Calfac,par.NormalExclusion,
                        par.DerivateExclusion,par.SaturationValue,&status);
            }else{
                initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr, 
                        par.LibraryFile, par.TesEventFile, par.PulseLength, par.scaleFactor, par.samplesUp, 
                        par.samplesDown, par.nSgms, par.detectSP, par.opmode, par.detectionMode, par.LrsT, 
                        par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod, par.EnergyMethod, 
                        par.filtEev, par.OFNoise, par.LagsOrNot, par.nLags, par.Fitting35, par.OFIter, 
                        par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength, par.preBuffer, par.monoenergy, 
                        par.hduPRECALWN, par.hduPRCLOFWM, par.largeFilter, par.intermediate, par.detectFile, 
                        par.filterFile, par.errorT, par.Sum0Filt, par.clobber, par.EventListSize, par.SaturationValue, par.tstartPulse1, 
                        par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);
            }
            CHECK_STATUS_BREAK(status);
            
            /*printf("%s %d %s","record_file->trigger_size0: ", record_file->trigger_size,"\n");
            if (record_file->trigger_size > trig_reclength)
            {
                record_file->trigger_size = trig_reclength;
            }*/
            // Build up TesRecord to read the file
            record = newTesRecord(&status);
            if (record_file->delta_t == -999) record_file->delta_t = 1./sampling_rate;
            allocateTesRecord(record,record_file->trigger_size,record_file->delta_t,0,&status);
            CHECK_STATUS_BREAK(status);
            
            // Iterate of records and do the reconstruction
            lastRecord = 0, nrecord = 0;    //last record required for SIRENA library creation
            while(getNextRecord(record_file,record,&status))
            {
                    if(!strcmp(par.Rcmethod,"PP"))
                    {
                            reconstructRecord(record,event_list,reconstruct_init,0,&status);
                    }
                    else
                    {
                            nrecord = nrecord + 1;
                            if(nrecord == record_file->nrows) lastRecord=1;
                            /*if(nrecord < 7905) 
                            {
                              continue;
                            }
                            else if(nrecord > 7905)
                            {
                              status=1;
                              CHECK_STATUS_BREAK(status);
                            }*/
                            /*if(nrecord > 1)
                            {
                            	status=1;
                                CHECK_STATUS_BREAK(status);
                            }*/
                            if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RALL") == 0) 
                                || (strcmp(par.EnergyMethod,"I2RNOL") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0))
                            {
                                strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
                            }
                        
                            //printf("%s %d %s","**TESRECONSTRUCTION nrecord = ",nrecord,"\n");
                            //printf("%s %d %s","record->trigger_size1: ",record->trigger_size,"\n");
                            //printf("%s %d %s", "pixid: ",record->pixid,"\n");
                            //printf("%s %d %s","ph_id: ",record->phid_list->phid_array[0],"\n");
                            reconstructRecordSIRENA(record,event_list,reconstruct_init_sirena,
                                                    lastRecord, nrecord, &pulsesAll, &optimalFilter, &status);
                    }
                    CHECK_STATUS_BREAK(status);

                    if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
                    {
                            // In THREADING mode, saveEventListToFile is not called until finishing with calculus 
                            // (ordering is necessary previously)  
                            if(!is_threading()){    
                                    //printf("\n %p - %f", outfile, record_file->delta_t);
                                    //printf("\nRecord single");
                                    //printf("\n%f - %ld", record->time, record->pixid);
                                    //printf("%s %d %s","**Before saveEventListToFile \n");
                                    //printf("%s %d %s","status2 = ",status,"\n");
                                    saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                                    //printf("%s %d %s","**After saveEventListToFile \n");
                                    //printf("%s %d %s","status3 = ",status,"\n");
                                    CHECK_STATUS_BREAK(status);
                                    //Reinitialize event list
                                    event_list->index=0;
                            }
                    }
            }
            
            if(is_threading()) 
            {
                    //printf("%s","**Threading...waiting \n");
                    th_end(&reconstruct_init_sirena, &pulsesAll, &optimalFilter);
                    //printf("%s %d %s","**Threading...after th_end: pulsesAll->ndetpulses", pulsesAll->ndetpulses,"\n");
                    //printf("%s %d %s","**Threading...after th_end: pulsesAll->size", pulsesAll->size,"\n");
                    int i = 1;
                    int aux = 1;
                    while((aux = th_get_event_list(&event_list, &record)) == 1)
                    {
                            //printf("%s %d %s","**Threading...event_list->size_energy: ", event_list->size_energy,"\n"); //Always 0
                            //printf("%s %e %s","**Threading...event_list->energies[0]: ", event_list->energies[0],"\n"); //Energy value
                            //printf("%s %e %s","**Threading...event_list->energies[1]: ", event_list->energies[1],"\n"); //Not error but non relevant value
                            //printf("%s %e %s","**Threading...event_list->energies[100000]: ", event_list->energies[100000],"\n"); Not error but non relevant value
                            saveEventListToFile(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                            //printf("%s","**Threading...after saveEventListToFile \n");
                            CHECK_STATUS_BREAK(status);
                            ++i;
                    }
            }
            
            if ((!strcmp(par.Rcmethod,"SIRENA")) && (pulsesAll->ndetpulses == 0)) 
            printf("%s","WARNING: no pulses have been detected\n");
            
            // Copy trigger keywords to event file
            copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
            CHECK_STATUS_BREAK(status);
            
            // Messages providing info of some columns
            char keywordvalue[9];
            char comment[MAXMSG];
            
            fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
            CHECK_STATUS_BREAK(status);
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE1", &keywordvalue, NULL, &status);
            strcpy(comment, "Starting time");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE1", keywordvalue, comment, &status);
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE2", &keywordvalue, NULL, &status);
            strcpy(comment, "Reconstructed-uncalibrated energy");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE2", keywordvalue, comment, &status);      
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE3", &keywordvalue, NULL, &status);
            strcpy(comment, "Average first 4 samples (derivative)");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE3", keywordvalue, comment, &status);      
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE4", &keywordvalue, NULL, &status);
            strcpy(comment, "Optimal filter length");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE4", keywordvalue, comment, &status);      
            
            fits_read_key(outfile->fptr, TSTRING, "TTYPE5", &keywordvalue, NULL, &status);
            strcpy(comment, "Starting time-starting time previous event");
            fits_update_key(outfile->fptr, TSTRING, "TTYPE5", keywordvalue, comment, &status);
            
            freeTesTriggerFile(&record_file,&status);
    }
    
    // Save GTI extension to event file
    GTI* gti=getGTIFromFileOrContinuous("none",keywords->tstart, keywords->tstop,keywords->mjdref, &status);
    saveGTIExt(outfile->fptr, "STDGTI", gti, &status);    
    CHECK_STATUS_BREAK(status);
    
    //Free memory
    freeReconstructInit(reconstruct_init);
    freeReconstructInitSIRENA(reconstruct_init_sirena);
    freePulsesCollection(pulsesAll);
    freeOptimalFilterSIRENA(optimalFilter);
    freeTesEventFile(outfile,&status);
    freeTesEventList(event_list);
    freeTesRecord(&record);
    freeSixtStdKeywords(keywords);
    CHECK_STATUS_BREAK(status);
 
  } while(0); // END of the error handling loop.
  
  if (EXIT_SUCCESS==status) 
  {
	headas_chat(3, "finished successfully!\n\n");
        time_t ttcurrent = time(0);
        printf("Elapsed time: %f\n", ((float)(ttcurrent - ttstart)));
	return(EXIT_SUCCESS);
  } 
  else 
  {
        time_t ttcurrent = time(0);
        printf("Elapsed time: %f\n", ((float)(ttcurrent - ttstart)));
	return(status);
  }
}

int getpar(struct Parameters* const par)
{
  // String input buffer.
  char* sbuffer=NULL;

  // Error status.
  int status=EXIT_SUCCESS;

  status=ape_trad_query_string("Rcmethod", &sbuffer);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the reconstruction method");
	  return(status);
  }
  strcpy(par->Rcmethod, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("RecordFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the optimal filter file");
    return(status);
  }
  strcpy(par->RecordFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_string("TesEventFile", &sbuffer);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the name of the event file");
    return(status);
  }
  strcpy(par->TesEventFile, sbuffer);
  free(sbuffer);

  status=ape_trad_query_int("PulseLength", &par->PulseLength);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the PulseLength parameter");
	  return(status);
  }
  assert(par->PulseLength > 0);

  status=ape_trad_query_int("EventListSize", &par->EventListSize);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the EventListSize parameter");
	  return(status);
  }

  status=ape_trad_query_bool("clobber", &par->clobber);
  if (EXIT_SUCCESS!=status) {
	  SIXT_ERROR("failed reading the clobber parameter");
	  return(status);
  }

  status=ape_trad_query_bool("history", &par->history);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the history parameter");
    return(status);
  }

  status=ape_trad_query_double("SaturationValue", &par->SaturationValue);
  if (EXIT_SUCCESS!=status) {
    SIXT_ERROR("failed reading the SaturationValue parameter");
    return(status);
  }

  if(strcmp(par->Rcmethod,"PP")==0){
	// PP  reconstruction method
	status=ape_trad_query_string("OptimalFilterFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the optimal filter file");
		return(status);
	}
	strcpy(par->OptimalFilterFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("PulseTemplateFile", &sbuffer);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the name of the pulse template file");
		return(status);
	}
	strcpy(par->PulseTemplateFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("Threshold", &par->Threshold);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the Threshold parameter");
		return(status);
	}

	status=ape_trad_query_double("Calfac", &par->Calfac);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the Calfac parameter");
		return(status);
	}

	status=ape_trad_query_int("NormalExclusion", &par->NormalExclusion);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the NormalExclusion parameter");
		return(status);
	}

	status=ape_trad_query_int("DerivateExclusion", &par->DerivateExclusion);
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading the DerivateExclusion parameter");
		return(status);
	}
  }else if(strcmp(par->Rcmethod,"SIRENA")==0){
	
	// SIRENA parameters
	status=ape_trad_query_string("LibraryFile", &sbuffer);
	strcpy(par->LibraryFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("scaleFactor", &par->scaleFactor);
    
	status=ape_trad_query_double("samplesUp", &par->samplesUp);
        
        status=ape_trad_query_double("samplesDown", &par->samplesDown);
  
	status=ape_trad_query_double("nSgms", &par->nSgms);
        
        status=ape_trad_query_int("detectSP", &par->detectSP);
  
	status=ape_trad_query_int("opmode", &par->opmode);
        
        status=ape_trad_query_string("detectionMode", &sbuffer);
	strcpy(par->detectionMode, sbuffer);
	free(sbuffer);
  
	status=ape_trad_query_double("LrsT", &par->LrsT);

	status=ape_trad_query_double("LbT", &par->LbT);

	status=ape_trad_query_int("intermediate", &par->intermediate);

	status=ape_trad_query_string("detectFile", &sbuffer);
	strcpy(par->detectFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("filterFile", &sbuffer);
	strcpy(par->filterFile, sbuffer);
	free(sbuffer);

	status=ape_trad_query_double("monoenergy", &par->monoenergy);
	
	status=ape_trad_query_bool("hduPRECALWN", &par->hduPRECALWN);
	status=ape_trad_query_bool("hduPRCLOFWM", &par->hduPRCLOFWM);
	
	status=ape_trad_query_int("largeFilter", &par->largeFilter);

	status=ape_trad_query_string("NoiseFile", &sbuffer);
	strcpy(par->NoiseFile, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_string("FilterDomain", &sbuffer);
	strcpy(par->FilterDomain, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_string("FilterMethod", &sbuffer);
	strcpy(par->FilterMethod, sbuffer);
	free(sbuffer);

	status=ape_trad_query_string("EnergyMethod", &sbuffer);
	strcpy(par->EnergyMethod, sbuffer);
	free(sbuffer);
        
        status=ape_trad_query_double("filtEev", &par->filtEev);

	status=ape_trad_query_string("OFNoise", &sbuffer);
	strcpy(par->OFNoise, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_int("LagsOrNot", &par->LagsOrNot);
        status=ape_trad_query_int("nLags", &par->nLags);
        status=ape_trad_query_int("Fitting35", &par->Fitting35);

	status=ape_trad_query_int("OFIter", &par->OFIter);

	status=ape_trad_query_bool("OFLib", &par->OFLib);
	
	strcpy(par->OFInterp, "DAB");
	
	status=ape_trad_query_string("OFStrategy", &sbuffer);
	strcpy(par->OFStrategy, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_int("OFLength", &par->OFLength);
        
        status=ape_trad_query_int("preBuffer", &par->preBuffer);
        
        status=ape_trad_query_int("errorT", &par->errorT);
        
        status=ape_trad_query_int("Sum0Filt", &par->Sum0Filt);

	//status=ape_trad_query_int("tstartPulse1", &par->tstartPulse1);
        status=ape_trad_query_string("tstartPulse1", &sbuffer);
	strcpy(par->tstartPulse1, sbuffer);
	free(sbuffer);
	
	status=ape_trad_query_int("tstartPulse2", &par->tstartPulse2);
	
	status=ape_trad_query_int("tstartPulse3", &par->tstartPulse3);
	
	status=ape_trad_query_double("energyPCA1", &par->energyPCA1);
	
	status=ape_trad_query_double("energyPCA2", &par->energyPCA2);
	
	status=ape_trad_query_string("XMLFile", &sbuffer);
	strcpy(par->XMLFile, sbuffer);
	free(sbuffer);
	
	if (EXIT_SUCCESS!=status) {
		SIXT_ERROR("failed reading some SIRENA parameter");
		return(status);
	}
	
	MyAssert((par->opmode == 0) || (par->opmode == 1), "opmode must be 0 or 1");
        int isNumber = 1;
        for (int i = 0; i < strlen(par->tstartPulse1); i++) 
        {
            if (isdigit(par->tstartPulse1[i]) == 0)    
            {
                isNumber = 0;
                break;
            }
        }
        if ((isNumber == 0) && (par->opmode == 0))
        {
            SIXT_ERROR("tstartPulse1 can not be a file if CALIBRATION mode");
            return(EXIT_FAILURE);
        }
        if ((isNumber == 0) && (strcmp(par->FilterDomain,"F") == 0))    // It is only implemented tstartPulse1 as a file for time domain
        {
            SIXT_ERROR("It is not possible to work in FREQUENCY domain if tstartPulse1 is a file => Change FilterDomain to TIME domain (T) ");
            return(EXIT_FAILURE);
        }
	  
	//MyAssert((par->intermediate == 0) || (par->intermediate == 1), "intermediate must be 0 or 1");
	if (par->intermediate == 1)
        {
            SIXT_ERROR("'intermediate = 1' is not available at the moment because the code to write an intermediate file is obsolete");
            return(EXIT_FAILURE);
        }
	
        if (par->opmode == 0) MyAssert(par->monoenergy > 0, "monoenergy must be greater than 0");
	
	MyAssert((strcmp(par->FilterDomain,"T") == 0) || (strcmp(par->FilterDomain,"F") == 0), "FilterDomain must be T or F");
	
	MyAssert((strcmp(par->FilterMethod,"F0") == 0) || (strcmp(par->FilterMethod,"B0") == 0),"FilterMethod must be F0 or B0");
	
	MyAssert((strcmp(par->EnergyMethod,"OPTFILT") == 0) || (strcmp(par->EnergyMethod,"WEIGHT") == 0) || (strcmp(par->EnergyMethod,"WEIGHTN") == 0) ||
		(strcmp(par->EnergyMethod,"I2R") == 0) || (strcmp(par->EnergyMethod,"I2RALL") == 0) || (strcmp(par->EnergyMethod,"I2RNOL") == 0) || 
		(strcmp(par->EnergyMethod,"I2RFITTED") == 0) || (strcmp(par->EnergyMethod,"PCA") == 0), "EnergyMethod must be OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED or PCA");
	
	MyAssert((strcmp(par->OFNoise,"NSD") == 0) || (strcmp(par->OFNoise,"WEIGHTM") == 0), "OFNoise must be NSD or WEIGHTM");
        
        MyAssert((strcmp(par->detectionMode,"AD") == 0) || (strcmp(par->detectionMode,"STC") == 0), "detectionMode must be AD or STC");
	
	MyAssert((par->LagsOrNot ==0) || (par->LagsOrNot ==1), "LagsOrNot must me 0 or 1");
        if ((par->nLags)%2 == 0)
	{
		SIXT_ERROR("parameter error: nLags must be odd");
		return(EXIT_FAILURE);
	}
	MyAssert((par->Fitting35 ==3) || (par->Fitting35 ==5), "Fitting35 must me 3 or 5");
        if ((par->Fitting35 ==3) && (par->nLags<3))
        {
                SIXT_ERROR("parameter error: nLags must be at least 3");
		return(EXIT_FAILURE);
        }
        if ((par->Fitting35 ==5) && (par->nLags<5))
        {
                SIXT_ERROR("parameter error: nLags must be at least 5");
		return(EXIT_FAILURE);
        }
        
        MyAssert((par->Sum0Filt ==0) || (par->Sum0Filt ==1), "Sum0Filt must be 0 or 1");

        if ((strcmp(par->EnergyMethod,"WEIGHT") == 0) && (par->LagsOrNot == 1))
	{
		SIXT_ERROR("parameter error: EnergyMethod=WEIGHT and Lags not implemented yet");
		return(EXIT_FAILURE);
	}
	
	MyAssert((par->OFIter ==0) || (par->OFIter ==1), "OFIter must be 0 or 1");
	
        // It was in order to not ask for the noise file if OFLib=1
	/*if ((par->OFLib == 1) && (strcmp(par->FilterMethod,"F0") != 0))
	{
		SIXT_ERROR("parameter error: If OFLib=yes => FilterMethod must be F0");
		return(EXIT_FAILURE);
	}*/
        
        if ((par->PulseLength < par->OFLength) && (strcmp(par->FilterDomain,"F") == 0))
        {
            SIXT_ERROR("Code is not prepared to run 0-padding in Frequency domain");
            // To run 0-padding in Frequency domain the steps should be:
            //1. Take the 8192-samples-length filter in Time domain
            //2. Cut the 0-padding length first samples (the first 4096 samples, or the first 2048 samples...) => 0-padding filter
            //3. FFT of the 0-padding filter
            //4. FFT of the 0-padding pulse (pulse cut according the 0-padding)
            //5. Scalar product in Frequency domain
            return(EXIT_FAILURE);
        }
        
	if ((strcmp(par->EnergyMethod,"WEIGHT") == 0) && (par->OFLib == 1))
	{
		SIXT_ERROR("parameter error: EnergyMethod=WEIGHT => OFLib should be 'no'");
		return(EXIT_FAILURE);
	}
	
	if ((strcmp(par->EnergyMethod,"OPTFILT") == 0) && (strcmp(par->OFNoise,"WEIGHTM") == 0) && (par->OFLib == 0))
	{
		SIXT_ERROR("parameter error: EnergyMethod=OPTFILT && OFNoise=WEIGHTM => OFLib should be 'yes'");
		return(EXIT_FAILURE);
	}
	
	MyAssert((strcmp(par->OFStrategy,"FREE") == 0) || (strcmp(par->OFStrategy,"BYGRADE") == 0) || (strcmp(par->OFStrategy,"FIXED") == 0), 
		 "OFStrategy must be FREE, BYGRADE or FIXED");
	
        MyAssert(par->OFLength > 0, "OFLength must be greater than 0");
        
        MyAssert(par->preBuffer >= 0, "preBuffer must be 0 or greater than 0");
	
	MyAssert(par->energyPCA1 > 0, "energyPCA1 must be greater than 0");
        MyAssert(par->energyPCA2 > 0, "energyPCA2 must be greater than 0");
        
        MyAssert(par->LbT > 0, "LbT must be greater than 0");
	
  } else {
	SIXT_ERROR("failed reading the Rcmethod parameter");
	return(EXIT_FAILURE);
  }
  return(status);
}

void MyAssert(int expr, char* msg)
{
    if (expr == 0)
    {
        printf("%s %s %s"," Assertion failure: ",msg,"\n");
        abort();
    }
}
