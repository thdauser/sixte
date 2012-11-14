#include "geninst.h"


////////////////////////////////////////////////////////////////////
// Function Declarations.
////////////////////////////////////////////////////////////////////


void parseGenInstXML(GenInst* const inst, 
		     const char* const filename, 
		     int* const status);


////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////


GenInst* newGenInst(int* const status) 
{
  // Allocate memory.
  GenInst* inst=(GenInst*)malloc(sizeof(GenInst));
  if (NULL==inst) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenInst failed");
    return(inst);
  }

  // Initialize all pointers with NULL.
  inst->tel=NULL;
  inst->det=NULL;
  inst->filename=NULL;
  inst->filepath=NULL;

  // Allocate memory for the GenTel and GenDet data structs.
  inst->det=newGenDet(status);
  CHECK_STATUS_RET(*status, inst);
  inst->tel=newGenTel(status);
  CHECK_STATUS_RET(*status, inst);

  return(inst);
}


void destroyGenInst(GenInst** const inst, int* const status)
{
  if (NULL!=*inst) {
    if (NULL!=(*inst)->tel) {
      destroyGenTel(&(*inst)->tel);
    }
    if (NULL!=(*inst)->det) {
      destroyGenDet(&(*inst)->det, status);
    }
    if (NULL!=(*inst)->filename) {
      free((*inst)->filename);
    }
    if (NULL!=(*inst)->filepath) {
      free((*inst)->filepath);
    }
    free(*inst);
    *inst=NULL;
  }
}


GenInst* loadGenInst(const char* const filename, int* const status)
{
  // Get a new and empty data structure.
  GenInst* inst=newGenInst(status);
  CHECK_STATUS_RET(*status, inst);

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
  CHECK_STATUS_RET(*status, inst);

  // Split rootname into the file path and the file name.
  char* lastslash=strrchr(rootname, '/');
  if (NULL==lastslash) {
    inst->filepath=(char*)malloc(sizeof(char));
    CHECK_NULL_RET(inst->filepath, *status, 
		   "memory allocation for filepath failed", inst);
    inst->filename=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(inst->filename, *status, 
		   "memory allocation for filename failed", inst);
    strcpy(inst->filepath, "");
    strcpy(inst->filename, rootname);
  } else {
    lastslash++;
    inst->filename=(char*)malloc((strlen(lastslash)+1)*sizeof(char));
    CHECK_NULL_RET(inst->filename, *status, 
		   "memory allocation for filename failed", inst);
    strcpy(inst->filename, lastslash);
      
    *lastslash='\0';
    inst->filepath=(char*)malloc((strlen(rootname)+1)*sizeof(char));
    CHECK_NULL_RET(inst->filepath, *status, 
		   "memory allocation for filepath failed", inst);
    strcpy(inst->filepath, rootname);
  }
  // END of storing the filename and filepath.


  // Read in the XML definition of the detector.
  parseGenInstXML(inst, filename, status);
  CHECK_STATUS_RET(*status, inst);


  // Allocate memory for the detector pixels.
  inst->det->line=
    (GenDetLine**)malloc(inst->det->pixgrid->ywidth*sizeof(GenDetLine*));
  if (NULL==inst->det->line) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for GenDet pixel array failed");
    return(inst);
  }
  int ii;
  for (ii=0; ii<inst->det->pixgrid->ywidth; ii++) {
    inst->det->line[ii]=newGenDetLine(inst->det->pixgrid->xwidth, status);
    if (EXIT_SUCCESS!=*status) return(inst);
  }


  return(inst);
}

