/* Header file for heasp routines
   Defines structures for RMF, RMFchain, ARF, PHA, BinFactors
   Includes prototypes for functions */

#include "fitsio.h"
#include "headas.h"

/* define the RMF structure */

struct RMF {

  long NumberChannels;                            /* Number of spectrum channels */
  long NumberEnergyBins;                          /* Number of response energies */
  long NumberTotalGroups;                         /* Total number of response groups */
  long NumberTotalElements;                       /* Total number of response elements */
  long FirstChannel;                              /* First channel number */
  long isOrder;                                   /* If true grating order information included */

  long* NumberGroups; /*NumberEnergyBins*/        /* Number of response groups for this energy bin */
  long* FirstGroup;  /*NumberEnergyBins*/         /* First response group for this energy bin (counts from 0)*/

  long* FirstChannelGroup; /*NumberTotalGroups*/  /* First channel number in this group */
  long* NumberChannelGroups; /*NumberTotalGroups*//* Number of channels in this group */
  long* FirstElement; /*NumberTotalGroups*/       /* First response element for this group (counts from 0)*/
  long* OrderGroup; /*NumberTotalGroups*/         /* The grating order of this group */

  float* LowEnergy; /*NumberEnergyBins*/         /* Start energy of bin */
  float* HighEnergy; /*NumberEnergyBins*/        /* End energy of bin */

  float* Matrix; /*NumberTotalElements*/         /* Matrix elements */

  float* ChannelLowEnergy; /*NumberChannels*/    /* Start energy of channel */
  float* ChannelHighEnergy; /*NumberChannels*/   /* End energy of channel */

  float AreaScaling;                             /* Value of EFFAREA keyword */
  float ResponseThreshold;                       /* Minimum value in response */

  char ChannelType[FLEN_KEYWORD];                /* Value of CHANTYPE keyword */
  char RMFVersion[FLEN_KEYWORD];                 /* MATRIX extension format version */
  char EBDVersion[FLEN_KEYWORD];                 /* EBOUNDS extension format version */
  char Telescope[FLEN_KEYWORD];                             
  char Instrument[FLEN_KEYWORD];
  char Detector[FLEN_KEYWORD];
  char Filter[FLEN_KEYWORD];
  char RMFType[FLEN_KEYWORD];                    /* Value of HDUCLAS3 keyword in MATRIX extension */
  char RMFExtensionName[FLEN_VALUE];             /* Value of EXTNAME keyword in MATRIX extension */
  char EBDExtensionName[FLEN_VALUE];             /* Value of EXTNAME keyword in EBOUNDS extension */

};

/* define the RMFchan structure - the RMF transposed with channels as rows */

struct RMFchan {

  long NumberChannels;                            /* Number of spectrum channels */
  long NumberEnergyBins;                          /* Number of response energies */
  long NumberTotalGroups;                         /* Total number of response groups */
  long NumberTotalElements;                       /* Total number of response elements */
  long FirstChannel;                              /* First channel number */
  long isOrder;                                   /* If true grating order information included */

  long* NumberGroups; /*NumberChannels*/          /* Number of response groups for this channel bin */
  long* FirstGroup;  /*NumberChannels*/           /* First response group for this channel bin (counts from 0)*/

  long* FirstEnergyGroup; /*NumberTotalGroups*/   /* First energy bin in this group */
  long* NumberEnergyGroups; /*NumberTotalGroups*/ /* Number of energy bins in this group */
  long* FirstElement; /*NumberTotalGroups*/       /* First response element for this group (counts from 0)*/
  long* OrderGroup; /*NumberTotalGroups*/         /* The grating order of this group */

  float* LowEnergy; /*NumberEnergyBins*/         /* Start energy of bin */
  float* HighEnergy; /*NumberEnergyBins*/        /* End energy of bin */

  float* Matrix; /*NumberTotalElements*/         /* Matrix elements */

  float* ChannelLowEnergy; /*NumberChannels*/    /* Start energy of channel */
  float* ChannelHighEnergy; /*NumberChannels*/   /* End energy of channel */

  float AreaScaling;                             /* Value of EFFAREA keyword */
  float ResponseThreshold;                       /* Minimum value in response */

  char ChannelType[FLEN_KEYWORD];                /* Value of CHANTYPE keyword */
  char RMFVersion[FLEN_KEYWORD];                 /* MATRIX extension format version */
  char EBDVersion[FLEN_KEYWORD];                 /* EBOUNDS extension format version */
  char Telescope[FLEN_KEYWORD];                             
  char Instrument[FLEN_KEYWORD];
  char Detector[FLEN_KEYWORD];
  char Filter[FLEN_KEYWORD];
  char RMFType[FLEN_KEYWORD];                    /* Value of HDUCLAS3 keyword in MATRIX extension */
  char RMFExtensionName[FLEN_VALUE];             /* Value of EXTNAME keyword in MATRIX extension */
  char EBDExtensionName[FLEN_VALUE];             /* Value of EXTNAME keyword in EBOUNDS extension */

};

/* define the ARF structure */

struct ARF {

  long NumberEnergyBins;                         /* Number of response energies */

  float* LowEnergy; /*NumberEnergyBins*/         /* Start energy of bin */
  float* HighEnergy; /*NumberEnergyBins*/        /* End energy of bin */

  float* EffArea;    /*NumberEnergyBins*/        /* Effective areas */

  char ARFVersion[FLEN_KEYWORD];                 /* SPECRESP extension format version */
  char Telescope[FLEN_KEYWORD];                             
  char Instrument[FLEN_KEYWORD];
  char Detector[FLEN_KEYWORD];
  char Filter[FLEN_KEYWORD];
  char ARFExtensionName[FLEN_VALUE];             /* Value of EXTNAME keyword in SPECRESP extension */

};

/* define the PHA structure */

struct PHA {

  long NumberChannels;                           /* Number of spectrum channels */
  long FirstChannel;                             /* First channel number */

  float* Pha; /*NumberChannels*/                 /* PHA data */
  float* StatError; /*NumberChannels*/           /* Statistical error */
  float* SysError; /*NumberChannels*/            /* Statistical error */

  int*   Quality; /*NumberChannels*/             /* Data quality */
  int*   Grouping;  /*NumberChannels*/           /* Data grouping */

  float* AreaScaling; /*NumberChannels*/         /* Area scaling factor */
  float* BackScaling; /*NumberChannels*/         /* Background scaling factor */

  float Exposure;                                /* Exposure time */
  float CorrectionScaling;                       /* Correction file scale factor */

  int Poisserr;                                  /* If true, errors are Poisson */
  char Datatype[FLEN_KEYWORD];                   /* "COUNT" for count data and "RATE" for count/sec */
  char Spectrumtype[FLEN_KEYWORD];               /* "TOTAL", "NET", or "BKG" */

  char ResponseFile[FLEN_FILENAME];              /* Response filename */
  char AncillaryFile[FLEN_FILENAME];             /* Ancillary filename */
  char BackgroundFile[FLEN_FILENAME];            /* Background filename */
  char CorrectionFile[FLEN_FILENAME];            /* Correction filename */

  char ChannelType[FLEN_KEYWORD];                /* Value of CHANTYPE keyword */
  char PHAVersion[FLEN_KEYWORD];                 /* PHA extension format version */
  char Telescope[FLEN_KEYWORD];                             
  char Instrument[FLEN_KEYWORD];
  char Detector[FLEN_KEYWORD];
  char Filter[FLEN_KEYWORD];
  char Datamode[FLEN_KEYWORD];

  char *XSPECFilter[100];                       /* Filter keywords */
};

/* define the BinFactors structure */

struct BinFactors {

  long NumberBinFactors;

  long *StartBin;
  long *EndBin;
  long *Binning;

};

/* function proto-types */

/* read the RMF matrix from an open FITS file - if there are multiple RMF extensions then
   read the one in RMFnumber */

int ReadRMFMatrix(fitsfile *fptr, long RMFnumber, struct RMF *rmf);

/* write the RMF matrix to an opened FITS file */

int WriteRMFMatrix(fitsfile *fptr, struct RMF *rmf);

/* read the RMF ebounds from an open FITS file - if there are multiple EBOUNDS extensions then
   read the one in EBDnumber */

int ReadRMFEbounds(fitsfile *fptr, long EBDnumber, struct RMF *rmf);

/* write the RMF ebounds to an opened FITS file */

int WriteRMFEbounds(fitsfile *fptr, struct RMF *rmf);
 
/* return the channel for a photon of the given input energy - draws random
   numbers to return NumberPhotons entries in the channel array */

void ReturnChannel(struct RMF *rmf, float energy, int NumberPhotons, long *channel);

/* normalize the response to unity in each energy */

void NormalizeRMF(struct RMF *rmf);

/* transpose the matrix */

void TransposeRMF(struct RMF *rmf, struct RMFchan *rmfchan);

/* return a single value from the matrix */

float ReturnRMFElement(struct RMF *rmf, long channel, long energybin);
float ReturnRMFchanElement(struct RMFchan *rmfchan, long channel, long energybin);

/* read the effective areas from an open FITS file - if there are multiple SPECRESP extensions then
   read the one in ARFFnumber */

int ReadARF(fitsfile *fptr, long ARFnumber, struct ARF *arf);

/* write the ARF to an opened FITS file */

int WriteARF(fitsfile *fptr, struct ARF *arf);

/* multiply the ARF into the RMF */

long MergeARFRMF(struct ARF *arf, struct RMF *rmf);

/* read the type I PHA extension from an open FITS file - if there are multiple PHA 
   extensions then read the one in PHAnumber */

int ReadPHAtypeI(fitsfile *fptr, long PHAnumber, struct PHA *pha);
 
/* read the type II PHA extension from an open FITS file - if there are multiple PHA
   extensions then read the one in PHAnumber - within the typeII extension reads the 
   spectra listed in the SpectrumNumber vector */

int ReadPHAtypeII(fitsfile *fptr, long PHAnumber, long NumberSpectra, long *SpectrumNumber, struct PHA **pha);

/* write the type I PHA extension to an open FITS file */

int WritePHAtypeI(fitsfile *fptr, struct PHA *pha);

/* write the type II PHA extension to an open FITS file */

int WritePHAtypeII(fitsfile *fptr, long NumberSpectra, struct PHA **pha);

/* return the type of a PHA extension */

int ReturnPHAtype(fitsfile *fptr, long PHAnumber);

/* return 0 if COUNTS column exists and is integer or COUNTS column does not exist */

int CheckPHAcounts(fitsfile *fptr, long PHAnumber);

/* return the number of spectra in a type II PHA extension */

long ReturnNumberofSpectra(fitsfile *fptr, long PHAnumber);

/* Read an ascii file with binning factors and load the binning array */

int SPReadBinningFile(char *filename, struct BinFactors *binning);

/* Set up a grouping array using the BinFactors structure */

int SPSetGroupArray(int inputSize, struct BinFactors *binning, int *grouping);

/* Bin an array using the information in the grouping array */

int SPBinArray(int inputSize, float *input, int *grouping, int mode, float *output);

