#include "initSIRENA.h"
#include <stdlib.h>

/***** SECTION 1 ************************************************************
* MyAssert function: This function displays an error message if the condition in 'expr' is true
*
* Parameters:
* - expr: Condition to be true in order to display the error message
* - msg: Message to be displayed
******************************************************************************/
void MyAssert(int expr, char* msg)
{
    if (expr == 0)
    {
        printf("%s %s %s"," Assertion failure: ",msg,"\n");
        abort();
    }
}
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* checkxmls function: Check if the XML file used to build the library is the same to be used to recconstruct (by checking the checksums)
*
* Parameters:
* - par: Structure containing the input parameters
******************************************************************************/
int checkXmls(struct Parameters* const par)
{
    // Error status.
    int status=EXIT_SUCCESS;

    // Read XMLCKSM from library file or calculate it
    fitsfile* libptr = NULL;
    // Move to "Primary" HDU of the library file
    fits_open_file(&libptr, par->LibraryFile, READONLY, &status);
    if (status != 0)
    {
        SIXT_ERROR("File given in LibraryFile does not exist");
        return(EXIT_FAILURE);
    }
    if (fits_movabs_hdu(libptr, 1, NULL, &status))
    {
        return(EXIT_FAILURE);
    }
    // Read XMLCKSM from library
    int keywordvalue;
    unsigned checksum_libXMLfile;
    char libXMLfile[1024] = "x";
    fits_read_key(libptr, TINT, "XMLCKSM", &keywordvalue, NULL, &status);
    if (status == EXIT_SUCCESS)
    {
        checksum_libXMLfile = keywordvalue;
    }
    else // or read full Primary HDU and store it in 'libheaderPrimary' and calculate XML checksum
    {
        status = EXIT_SUCCESS;
        int numberkeywords;
        char *libheaderPrimary = NULL;
        if (fits_hdr2str(libptr, 0, NULL, 0,&libheaderPrimary, &numberkeywords, &status))
        {
            free(libheaderPrimary);
            return(EXIT_FAILURE);
        }

        char *xml_pointer = NULL;
        char *XMLFile_pointer = NULL;
        char wholelibXMLfile[1024] = "x";
        xml_pointer = strstr(libheaderPrimary,".xml");
        if(!xml_pointer)
        {
            SIXT_ERROR("XML file info not included in Primary HDU in library file");
            return(EXIT_FAILURE);
        }
        XMLFile_pointer = strstr(libheaderPrimary,"XMLFile = ");
        if(!XMLFile_pointer)
        {
            SIXT_ERROR("XML file info not included in Primary HDU in library file");
            return(EXIT_FAILURE);
        }
        strncpy(wholelibXMLfile,XMLFile_pointer+10,xml_pointer-XMLFile_pointer-10+4);  // 10 -> "XMLFile = ", 4 -> ".xml"

        int HISTORYnum = 0;
        int k;
        int lengthTOTAL = (int)strlen(wholelibXMLfile);
        int spaces[lengthTOTAL];
        k=0;
        do
        {
            if(wholelibXMLfile[k]==' ')
            {
                if (HISTORYnum % 2 == 0)
                {  // Par
                    spaces[HISTORYnum] = k-7;
                }
                else
                {
                    spaces[HISTORYnum] = k;
                }

                HISTORYnum++;
            }
            k++;
        }while(k<=lengthTOTAL);

        if (HISTORYnum != 0) HISTORYnum = HISTORYnum/2;

        char libXMLfile2[1024] = "x";

        if (HISTORYnum == 0)
        {
            strcpy(libXMLfile,wholelibXMLfile);
        }
        else
        {
            for (int i=0;i<=HISTORYnum;i++)
            {
                if (i == 0) // First
                {
                    subString (wholelibXMLfile, 0, spaces[i], libXMLfile);
                }
                else
                {
                    if (i == HISTORYnum) // Last
                    {
                        subString (wholelibXMLfile, spaces[2*i-1]+1, lengthTOTAL-spaces[2*i-1]-1, libXMLfile2);
                        strcat(libXMLfile,libXMLfile2);
                        memset(libXMLfile2,0,1024);
                    }
                    else
                    {
                        subString (wholelibXMLfile, spaces[2*i-1]+1, spaces[2*i]-spaces[2*i-1]-1, libXMLfile2);
                        strcat(libXMLfile,libXMLfile2);
                        memset(libXMLfile2,0,1024);
                    }
                }
            }
        }

        if (libheaderPrimary != NULL)
        {
            free(libheaderPrimary);
        }

        FILE *fp_libXMLfile;
        size_t len_libXMLfile;
        char buf_libXMLfile[4096];
        if (NULL == (fp_libXMLfile = fopen(libXMLfile, "rb")))
        {
            printf("Unable to open XML from library, %s, for reading it and calculate its checksum.\n", libXMLfile);
            return -1;
        }
        len_libXMLfile = fread(buf_libXMLfile, sizeof(char), sizeof(buf_libXMLfile), fp_libXMLfile);
        checksum_libXMLfile = checksum(buf_libXMLfile, len_libXMLfile, 0);
    }
    fits_close_file(libptr,&status);

    // Calculate checksum from the XML input FITS file
    FILE *fp_reconsXMLfile;
    size_t len_reconsXMLfile;
    unsigned checksum_reconsXMLfile;
    char reconsXMLfile[1024] = "x";
    char buf_reconsXMLfile[4096];
    strcpy(reconsXMLfile,par->XMLFile);
    if (NULL == (fp_reconsXMLfile = fopen(reconsXMLfile, "rb")))
    {
          printf("Unable to open provided XML, %s, for reading it and calculate its checksum\n", reconsXMLfile);
          return -1;
    }
    len_reconsXMLfile = fread(buf_reconsXMLfile, sizeof(char), sizeof(buf_reconsXMLfile), fp_reconsXMLfile);
    checksum_reconsXMLfile = checksum(buf_reconsXMLfile, len_reconsXMLfile, 0);

    // Compare checksums
    if (checksum_libXMLfile == checksum_reconsXMLfile) status = 0;
    else status = 1;

    if (status != 0)
    {
        SIXT_ERROR("XML file from library FITS file and from input parameter do not match (checksum)");
        printf("The checksum of XML from library, %s, is %#x\n", libXMLfile, checksum_libXMLfile);
        printf("The checksum of provided XML, %s, is %#x\n", reconsXMLfile, checksum_reconsXMLfile);
        return(EXIT_FAILURE);
    }

    return(status);
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* subString: This function extracts some elements from an array of characters
*
* Parameters:
* - input: Array of characters from which some elements are extracted
* - offset: Offset
* - len: Length (number of elements to extract)
* - dest: Array of characters into which the extracted characters are written
******************************************************************************/
char* subString (const char* input, int offset, int len, char* dest)
{
  int input_len = strlen (input);

  if (offset + len > input_len)
  {
     return NULL;
  }

  strncpy (dest, input + offset, len);
  return dest;
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* getSamplingrate_trigreclength_Filei: This function gets the sampling rate and the trig_reclength
*
* Steps:
* - Open FITS file
* - Check if input FITS file have been simulated with TESSIM or XIFUSIM
* - Check if input XML file and XMl file to build the library to be used to
*   reconstruct are the same
* - Get the sampling rate from the HISTORY keyword from the input FITS file
*   and check with sampling rate from XML file
* - If xifusim file => Get 'trig_reclength' from the HISTORY keyword from the
*   input FITS file 'trig_reclength' is necessary if SIRENA is going to run in
*   THREADING mode
* - Close FITS file
* - Free memory
*
* Parameters:
* - inputFile: Input file name
* - par: Input parameters
* - samplingrate: (In) Sampling rate from XML file => (Out) Sampling rate
* - trigreclength: Necessary if SIRENA is going to run in THREADING mode
******************************************************************************/
int getSamplingrate_trigreclength_Filei (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength)
{
    double samplingrate_XML = *samplingrate; // To work with a more meaningful name

    // Error status.
    int status=EXIT_SUCCESS;

    // Open file
    fitsfile* fptr = NULL;
    fits_open_file(&fptr, inputFile, READONLY, &status);
    if (status != 0)
    {
        SIXT_ERROR("File given in RecordFile does not exist");
        return(EXIT_FAILURE);
    }

    // Check if input FITS file have been simulated with TESSIM or XIFUSIM
    char extname[20];
    int extver = 0;
    int tessimOrxifusim = -999;     // 0: tessim, 1: xifusim
    strcpy(extname,"RECORDS");
    fits_movnam_hdu(fptr, ANY_HDU, extname, extver, &status);
    if (status != 0)
    {
        status = 0;
        strcpy(extname,"TESRECORDS");
        fits_movnam_hdu(fptr, ANY_HDU, extname, extver, &status);
        if (status != 0)
        {
            SIXT_ERROR("Cannot move to TESRECORDS HDU in input FITS file");
            return(EXIT_FAILURE);
        }
        else    // TESRECORDS => xifusim
        {
            tessimOrxifusim = 1;
        }
    }
    else // RECORDS => tessim
    {
        tessimOrxifusim = 0;
    }
    if (tessimOrxifusim == -999)
    {
        SIXT_ERROR("Neither the 'RECORDS' nor 'TESRECORDS' HDUs are in the input FITS file");
        return(EXIT_FAILURE);
    }

    // Check if input XML file and XMl file to build the library to be used to reconstruct are the same
    if (par.opmode == 1)
    {
        status = checkXmls(&par);
    }

    // Get the sampling rate from the HISTORY keyword from the input FITS file
    // and check with sampling rate from XML file
    // Move to "Primary" HDU
    fits_movabs_hdu(fptr, 1, NULL, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    // Read full Primary HDU and store it in 'headerPrimary'
    char *headerPrimary = NULL;
    int numberkeywords;
    if (fits_hdr2str(fptr, 0, NULL, 0,&headerPrimary, &numberkeywords, &status))
    {
        free(headerPrimary);
        return(EXIT_FAILURE);
    }
    // Pointer to where the text "sample_rate=" is in HISTORY block
    char *sample_rate_pointer = NULL;
    sample_rate_pointer = strstr (headerPrimary,"sample_rate=");
    /*if(!sample_rate_pointer)
    {
        // Read it from XML file
        *sampling_rate = *samplingrate; // Input parameter 'samplingrate' is sampling frequency read from XML file
    }
    else*/
    // If no 'sample_rate' in HISTORY => sampling frequency from XML is used
    if (sample_rate_pointer)
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
        double samplingrate_HISTORY = atof(characters_after_srate);
        //if (((*samplingrate != -999.0) && (samplingrate_HISTORY != -999.0)) && (*samplingrate != samplingrate_HISTORY))
        if (((samplingrate_XML != -999.0) && (samplingrate_HISTORY != -999.0)) && (samplingrate_XML != samplingrate_HISTORY))
        {
            SIXT_ERROR("Sampling rate from input FITS file and from XML file do not match");
            return(EXIT_FAILURE);
        }
        *samplingrate = samplingrate_HISTORY;
    }

    // If xifusim file => Get 'trig_reclength' from the HISTORY keyword from the input FITS file
    // 'trig_reclength' is necessary if SIRENA is going to run in THREADING mode
    if (tessimOrxifusim == 1) //xifusim simulated file (with TESRECORDS)
    {
        // Pointer to where the text "trig_reclength=" is in HISTORY block
        char *trig_reclength_pointer = NULL;
        trig_reclength_pointer = strstr (headerPrimary,"trig_reclength=");
        if(!trig_reclength_pointer)
        {
            printf("%s","Attention: 'trig_reclength' is not in the input FITS file (necessary if SIRENA is going to run in THREADING mode)\n");
        }
        else
        {
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
            *trigreclength = atoi(characters_after_treclength);
        }
        // Free memory
        if (trig_reclength_pointer != NULL)
        {
            trig_reclength_pointer = NULL;
            free(trig_reclength_pointer);
        }
    }

    // Close file
    fits_close_file(fptr,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Free memory
    if (sample_rate_pointer != NULL)
    {
        sample_rate_pointer = NULL;
        free(sample_rate_pointer);
    }
    if (headerPrimary != NULL)
    {
        free(headerPrimary);
    }

    return(status);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* getSamplingrate_trigreclength: This function gets the sampling rate and the trig_reclength no
*                                matter if only a FITS file or more (inputFile can start with '@'
*                                or not, input file or files can have been simulated with TESSIM
*                                or XIFUSIM)
*
* Parameters:
* - inputFile: Input file name
* - par: Input parameters
* - samplingrate: (In) Sampling rate from XML file => (Out) Sampling rate
* - trigreclength: Necessary if SIRENA is going to run in THREADING mode
* - numfits: Number of FITS files to work with
******************************************************************************/
int getSamplingrate_trigreclength (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength, int* numfits)
{
    // Error status.
    int status=EXIT_SUCCESS;
    char firstcharAux = inputFile[0];
    char firstchar[2] = {firstcharAux , '\0'};

    if (strcmp(firstchar,"@") != 0) // Only A input FITS file
    {
        *numfits = 1;
        status = getSamplingrate_trigreclength_Filei(inputFile, par, samplingrate, trigreclength);
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    }
    else // More than A input FITS file
    {
        // Open @ file to know the number of FITS files to work with
        int numfitsAux = 0;
        FILE *filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");
        if (filetxt == NULL)
        {
            printf("%s","File given in RecordFile does not exist\n");
            status = 104;
        }
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        char filefits[256];
        while(fscanf(filetxt,"%s",filefits)!=EOF)
        {
            numfitsAux++;
        }
        *numfits = numfitsAux;
        fclose(filetxt);

        filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");

        for (int j=0;j<(*numfits);j++)   // For every FITS file
        {
            fgets(filefits, 256, filetxt);
            strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)

            status = getSamplingrate_trigreclength_Filei(inputFile, par, samplingrate, trigreclength);
            if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        }

        fclose(filetxt);
    }

    return(status);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* fillReconstructInitSIRENAGrading: This function reads the grading data from the XML file
*                            and store it in 'reconstruct_init_sirena->grading'
*
* Parameters:
* - par: Input parameters
* - det: Pixel detector
* - reconstruct_init_sirena: Parameters to run SIRENA
******************************************************************************/
int fillReconstructInitSIRENAGrading (struct Parameters par, AdvDet *det, ReconstructInitSIRENA** reconstruct_init_sirena)
{
    // Error status.
    int status=EXIT_SUCCESS;

    (*reconstruct_init_sirena)->grading = NULL;
    (*reconstruct_init_sirena)->grading = (Grading*)malloc(sizeof(Grading));
    (*reconstruct_init_sirena)->grading->ngrades = 0;
    (*reconstruct_init_sirena)->grading->value  = NULL;
    (*reconstruct_init_sirena)->grading->gradeData = NULL;

    if ((det->nrecons == 0) && (det->npix == 0))
    {
        SIXT_ERROR("The provided XMLFile does not have the grading info");
        return(EXIT_FAILURE);
    }
    else if ((det->nrecons == 0) && (det->npix != 0))
    {
        if (det->pix->grades == NULL)
        {
            SIXT_ERROR("The provided XMLFile does not have the grading info");
            return(EXIT_FAILURE);
        }
        (*reconstruct_init_sirena)->grading->ngrades=det->pix->ngrades;
        (*reconstruct_init_sirena)->grading->gradeData = gsl_matrix_alloc(det->pix->ngrades,3);
        for (int i=0;i<det->pix->ngrades;i++)
        {
            if ((int) (det->pix->grades[i].gradelim_post) != 8)
            {
                if ((int) (det->pix->grades[i].grade_preBuffer) >= (int) (det->pix->grades[i].gradelim_post))
                {
                    printf("%s %d %s %d %s","preBuffer=",(int) (det->pix->grades[i].grade_preBuffer)," for filter length ",(int) (det->pix->grades[i].gradelim_post),"\n");
                    SIXT_ERROR("preBuffer values provided in the XML file should be lower than corresponding filter lengths");
                    return(EXIT_FAILURE);
                }
            }
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,0,(int) (det->pix->grades[i].gradelim_pre));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,1,(int) (det->pix->grades[i].gradelim_post));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,(int) (det->pix->grades[i].grade_preBuffer));
            if (((int) (det->pix->grades[i].gradelim_post) == 8) && ((int) (det->pix->grades[i].grade_preBuffer) != 0))
            {
                printf("%s","Attention: preBuffer!=0 for low resolution (filter length 8) but preBuffer=0 is going to be used\n");
                gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,0);
            }
        }
    }
    else if(((det->nrecons != 0) && (det->npix == 0)) || ((det->nrecons == 1) && (det->npix == 1)))
    {
        if (det->recons->grades == NULL)
        {
            SIXT_ERROR("The provided XMLFile does not have the grading info");
            return(EXIT_FAILURE);
        }
        (*reconstruct_init_sirena)->grading->ngrades=det->recons->ngrades;
        (*reconstruct_init_sirena)->grading->gradeData = gsl_matrix_alloc(det->recons->ngrades,3);
        for (int i=0;i<det->recons->ngrades;i++)
        {
            if ((int) (det->recons->grades[i].gradelim_post) != 8)
            {
                if ((int) (det->recons->grades[i].grade_preBuffer) >= (int) (det->recons->grades[i].gradelim_post))
                {
                    printf("%s %d %s %d %s","preBuffer=",(int) (det->recons->grades[i].grade_preBuffer)," for filter length ",(int) (det->recons->grades[i].gradelim_post),"\n");
                    SIXT_ERROR("preBuffer values provided in the XML file should be lower than corresponding filter lengths");
                    return(EXIT_FAILURE);
                }
            }
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,0,(int) (det->recons->grades[i].gradelim_pre));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,1,(int) (det->recons->grades[i].gradelim_post));
            gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,(int) (det->recons->grades[i].grade_preBuffer));
            if (((int) (det->recons->grades[i].gradelim_post) == 8) && ((int) (det->recons->grades[i].grade_preBuffer) != 0))
            {
                printf("%s","Attention: preBuffer!=0 for low resolution (filter length 8) but preBuffer=0 is going to be used\n");
                gsl_matrix_set((*reconstruct_init_sirena)->grading->gradeData,i,2,0);
            }
        }
    }

    /*int OFlengthvsposti = 0;
    if ((par.preBuffer == 1) && (par.opmode == 1))
    {
        for (int i=0;i<(*reconstruct_init_sirena)->grading->ngrades;i++)
        {
            if (par.OFLength == gsl_matrix_get((*reconstruct_init_sirena)->grading->gradeData,i,1))
            {
                OFlengthvsposti = 1;
                break;
            }
        }
        if (OFlengthvsposti == 0)
        printf("%s %d %s","par.: ",par.Pul,"\n");
        printf("%s %d %s","par.OFLength: ",par.OFLength,"\n");
        if ((OFlengthvsposti == 0) && ((*reconstruct_init_sirena)->pulse_length >= (*reconstruct_init_sirena)->OFLength)) // Not 0-padding
        {
            SIXT_ERROR("The grading/preBuffer info of the XML file does not match the OFLength input parameter");
            return(EXIT_FAILURE);
        }
    }*/

    return(status);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* callSIRENA_Filei: This function calls SIRENA to build a library or reconstruct energies
*
* Steps:
* - Open record file
* - initializeReconstructionSIRENA
* - Build up TesRecord to read the file
* - Iterate of records and run SIRENA:
*   - reconstructRecordSIRENA
*   - Save events to the event_list
*   - Copy trigger keywords to event file
*   - Close file
*
* Parameters:
* - inputFile: Input file name
* - keywords: Sixt standard keywords structure
* - reconstruct_init_sirena: Parameters to run SIRENA
* - par: Input parameters
* - sampling_rate: Sampling rate
* - trig_reclength: Necessary if SIRENA is going to run in THREADING mode
* - pulsesAll: Structure containing the detected pulses
* - outfile: Output events FITS file
******************************************************************************/
int callSIRENA_Filei(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFile * outfile)
{
    // Error status.
    int status=EXIT_SUCCESS;

    // Open record file
    TesTriggerFile* record_file;
    record_file = openexistingTesTriggerFile(inputFile,keywords,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    initializeReconstructionSIRENA(reconstruct_init_sirena, par.RecordFile, record_file->fptr,
        par.LibraryFile, par.TesEventFile, par.OFLengthNotPadded, par.scaleFactor, par.samplesUp,
        par.samplesDown, par.nSgms, par.detectSP, par.opmode, par.detectionMode, par.LrsT,
        par.LbT, par.NoiseFile, par.FilterDomain, par.FilterMethod, par.EnergyMethod,
        par.filtEev, par.Ifit, par.OFNoise, par.LagsOrNot, par.nLags, par.Fitting35, par.OFIter,
        par.OFLib, par.OFInterp, par.OFStrategy, par.OFLength, par.preBuffer, par.monoenergy,
        par.hduPRECALWN, par.hduPRCLOFWM, par.largeFilter, par.intermediate, par.detectFile,
        par.errorT, par.Sum0Filt, par.clobber, par.EventListSize, par.SaturationValue, par.tstartPulse1,
        par.tstartPulse2, par.tstartPulse3, par.energyPCA1, par.energyPCA2, par.XMLFile, &status);

    // Build up TesRecord to read the file
    TesRecord* record;
    record = newTesRecord(&status);
    if ((record_file->delta_t == -999) && (sampling_rate == -999))
    {
        SIXT_ERROR("Cannot read or get the sampling rate neither from the input FITS file nor the XML file. Please, include the DELTAT keyword (inverse of sampling rate) in the input FITS file before running the task again");
        return(EXIT_FAILURE);
    }
    if (record_file->delta_t == -999) record_file->delta_t = 1./sampling_rate;
    if ((*trig_reclength) == -999) *trig_reclength = record_file->trigger_size;
    allocateTesRecord(record,*trig_reclength,record_file->delta_t,0,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Iteration of records and running SIRENA
    // Read numrecords (NAXIS2) to simulate a progress bar
    fitsfile *fptr;
    if (fits_open_file(&fptr, inputFile, READONLY, &status)) {
        SIXT_ERROR("Cannot open the input FITS file to read the number of records (NAXIS2)(in 'callSIRENA_Filei' function).");
        return status;
    }
    fits_movnam_hdu(fptr, ANY_HDU,"TESRECORDS", 0, &status);
    if (status != EXIT_SUCCESS)
    {
        status = EXIT_SUCCESS;
        fits_movnam_hdu(fptr, ANY_HDU,"RECORDS", 0, &status);
        if (status != EXIT_SUCCESS)
        {
            SIXT_ERROR("Cannot find the TESRECORDS (or RECORDS) HDU in the input FITS file to read the number of records (NAXIS2)(in 'callSIRENA_Filei' function).");
            return status;
        }
    }
    int numrecords;
    fits_read_key(fptr, TINT, "NAXIS2", &numrecords, NULL, &status);
    fits_close_file(fptr, &status);

    fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    int lastRecord = 0, nrecord = 0, startRecordGroup = 0;    //last record required for SIRENA library creation
    // Build up TesEventList to recover the results of SIRENA
    TesEventList* event_list = newTesEventListSIRENA(&status);
    allocateTesEventListTrigger(event_list,par.EventListSize,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    while(getNextRecord(record_file,record,&lastRecord,&startRecordGroup,&status))
    {
        // Progress bar
        float progress = (float)nrecord / (numrecords-1);
        int bar_width = 50;
        int pos = bar_width * progress;
        printf("Simulating |");
        for (int j = 0; j < bar_width; j++) {
            if (j < pos)
                printf("=");
            else if (j == pos)
                printf(">");
            else
                printf(" ");
        }
        printf("| %.2f%%\r", progress * 100);
        fflush(stdout);
        usleep(100000); // Sleep for a short time to simulate progress
        if (nrecord == numrecords-1)
        {
            printf("\n"); // New line after ending "Simulating...."
        }

        nrecord = nrecord + 1;
        /*if(nrecord < 5887)
        {
            continue;
        }
        else if(nrecord > 5920)
        {
            status=1;
            CHECK_STATUS_BREAK(status);
        }*/
        /*if(nrecord > 1)
        {
          	status=1;
            CHECK_STATUS_BREAK(status);
        }*/
        if ((strcmp(par.EnergyMethod,"I2R") == 0) || (strcmp(par.EnergyMethod,"I2RFITTED") == 0)
            || (strcmp(par.EnergyMethod,"I2RDER") == 0))
        {
            strcpy(reconstruct_init_sirena->EnergyMethod,par.EnergyMethod);
        }

        //printf("%s %d %s","** nrecord = ",nrecord,"\n");
        reconstructRecordSIRENA(record,*trig_reclength, event_list,reconstruct_init_sirena,
            lastRecord, startRecordGroup, &pulsesAll, &status);
        CHECK_STATUS_BREAK(status);

        if ((strcmp(par.EnergyMethod,"PCA") != 0) || ((strcmp(par.EnergyMethod,"PCA") == 0) && lastRecord == 1))
        {
            // In THREADING mode, saveEventListToFileSIRENA is not called until finishing with calculus
            // (ordering is necessary previously)
            if(!is_threading()){
                saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
                CHECK_STATUS_BREAK(status);
                //Reinitialize event list
                event_list->index=0;
            }

            IntegrafreeTesEventListSIRENA(event_list);
        }
    }

    if(is_threading())
    {
        //printf("%s","**Threading...waiting \n");
        //th_end(&reconstruct_init_sirena, &pulsesAll);
        th_end(reconstruct_init_sirena, &pulsesAll);
        int i = 1;
        int aux = 1;
        while((aux = th_get_event_list(&event_list, &record)) == 1)
        {
            saveEventListToFileSIRENA(outfile,event_list,record->time,record_file->delta_t,record->pixid,&status);
            if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
            ++i;
        }
        IntegrafreeTesEventListSIRENA(event_list);
    }

    if (pulsesAll->ndetpulses == 0)
        printf("%s","WARNING: no pulses have been detected\n");

    // Copy trigger keywords to event file
    copyTriggerKeywords(record_file->fptr,outfile->fptr,&status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

    // Messages providing info of some columns
    char keywordvalue[9];
    char comment[MAXMSG];

    fits_movnam_hdu(outfile->fptr, ANY_HDU,"EVENTS", 0, &status);
    if (status != EXIT_SUCCESS) return(EXIT_FAILURE);

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
    freeTesRecord(&record);

    return (status);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* callSIRENA: This function calls SIRENA to build a library or reconstruct energies no
*             matter if only a FITS file or more (inputFile can start with '@' or not)
*
* Parameters:
* - inputFile: Input file name
* - keywords: Sixt standard keywords structure
* - reconstruct_init_sirena: Parameters to run SIRENA
* - par: Input parameters
* - sampling_rate: Sampling rate
* - trig_reclength: Necessary if SIRENA is going to run in THREADING mode
* - pulsesAll: Structure containing the detected pulses
* - outfile: Output events FITS file
******************************************************************************/
int callSIRENA(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFile * outfile)
{
    // Error status.
    int status=EXIT_SUCCESS;

    char firstcharAux = inputFile[0];
    char firstchar[2] = {firstcharAux , '\0'};

    if (strcmp(firstchar,"@") != 0) // Only A input FITS file
    {
        status = callSIRENA_Filei(inputFile, keywords, reconstruct_init_sirena, par, sampling_rate, trig_reclength, pulsesAll, outfile);
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
    }
    else // More than A input FITS file
    {
        // Open @ file to know the number of FITS files to work with
        int numfits = 0;
        FILE *filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");
        if (filetxt == NULL)
        {
            printf("%s","File given in RecordFile does not exist\n");
            status = 104;
        }
        if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        char filefits[256];
        while(fscanf(filetxt,"%s",filefits)!=EOF)
        {
            numfits++;
        }
        fclose(filetxt);

        filetxt = fopen(strndup(inputFile+1, strlen(inputFile)-1), "r");

        for (int j=0;j<numfits;j++)   // For every FITS file
        {
            fgets(filefits, 256, filetxt);
            strtok(filefits, "\n");     // To delete '/n' from filefits (if not, 'fits_open_file' can not open the file)

            status = callSIRENA_Filei(inputFile, keywords, reconstruct_init_sirena, par, sampling_rate, trig_reclength, pulsesAll, outfile);
            if (status != EXIT_SUCCESS) return(EXIT_FAILURE);
        }

        fclose(filetxt);
    }

    return (status);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
