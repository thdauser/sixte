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


 Copyright 2007-2014 Christian Schmid, FAU
 Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
 Erlangen-Nuernberg
 */

#include "geninst.h"

////////////////////////////////////////////////////////////////////
// Local data type declarations.
////////////////////////////////////////////////////////////////////

/** Data structure given to the XML handler to transfer data. */
struct XMLParseData {
	GenInst* inst;
	unsigned int seed;
	int status;
};

////////////////////////////////////////////////////////////////////
// Static variables.
////////////////////////////////////////////////////////////////////

/** Flag indicating that the erodetbkgrndgen module is initialized and
 operational. */
static int auxBkgInitialized = 0;

////////////////////////////////////////////////////////////////////
// Static function declarations.
////////////////////////////////////////////////////////////////////

/** Handler for the start of an XML element. */
static void GenInstXMLElementStart(void* data, const char* el,
		const char** attr);
/** Handler for the end of an XML element. */
static void GenInstXMLElementEnd(void* data, const char* el);

////////////////////////////////////////////////////////////////////
// Program Code.
////////////////////////////////////////////////////////////////////

GenInst* newGenInst(int* const status) {
	// Allocate memory.
	GenInst* inst = (GenInst*) malloc(sizeof(GenInst));
	if (NULL == inst) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("memory allocation for GenInst failed");
		return (inst);
	}

	// Initialize all pointers with NULL.
	inst->tel = NULL;
	inst->det = NULL;
	inst->filename = NULL;
	inst->filepath = NULL;
	inst->telescop = NULL;
	inst->instrume = NULL;

	// Allocate memory for the GenTel and GenDet data structs.
	inst->det = newGenDet(status);
	CHECK_STATUS_RET(*status, inst);
	inst->tel = newGenTel(status);
	CHECK_STATUS_RET(*status, inst);

	return (inst);
}

void destroyGenInst(GenInst** const inst, int* const status) {
	if (NULL != *inst) {
		if (NULL != (*inst)->tel) {
			destroyGenTel(&(*inst)->tel);
		}
		if (NULL != (*inst)->det) {
			destroyGenDet(&(*inst)->det);
		}
		if (1 == auxBkgInitialized) {
			bkgCleanUp(status);
			auxBkgInitialized = 0;
		}
		if (NULL != (*inst)->filename) {
			free((*inst)->filename);
		}
		if (NULL != (*inst)->filepath) {
			free((*inst)->filepath);
		}
		if (NULL != (*inst)->telescop) {
			free((*inst)->telescop);
		}
		if (NULL != (*inst)->instrume) {
			free((*inst)->instrume);
		}
		free(*inst);
		*inst = NULL;
	}
}

void parseGenInstXML(GenInst* const inst, const char* const filename,
		const unsigned int seed, int* const status) {
	headas_chat(5, "read instrument setup from XML file '%s' ...\n", filename);

	// Read the XML data from the file.
	// Open the specified file.
	FILE* xmlfile = fopen(filename, "r");
	if (NULL == xmlfile) {
		*status = EXIT_FAILURE;
		char msg[MAXMSG];
		sprintf(msg, "failed opening XML "
				"file '%s' for read access", filename);
		SIXT_ERROR(msg);
		return;
	}

	// The data are read from the XML file and stored in xmlbuffer
	// without any modifications.
	struct XMLBuffer* xmlbuffer = newXMLBuffer(status);
	CHECK_STATUS_VOID(*status);

	// Input buffer with an additional byte at the end for the
	// termination of the string.
	const int buffer_size = 256;
	char buffer[buffer_size + 1];
	// Number of chars in buffer.
	int len;

	// Read all data from the file.
	do {
		// Get a piece of input into the buffer.
		len = fread(buffer, 1, buffer_size, xmlfile);
		buffer[len] = '\0'; // Terminate the string.
		addString2XMLBuffer(xmlbuffer, buffer, status);
		CHECK_STATUS_VOID(*status);
	} while (!feof(xmlfile));

	// Close the file handler to the XML file.
	fclose(xmlfile);

	// Before expanding loops in the XML file, add the included code to it.
	expandIncludesXML(xmlbuffer, filename, status);
	CHECK_STATUS_VOID(*status);

	// Before actually parsing the XML code, expand the loops and
	// arithmetic operations in the GenDet XML description.
	// The expansion algorithm repeatedly scans the XML code and
	// searches for loop tags. It replaces the loop tags by repeating
	// the contained XML code.
	expandXML(xmlbuffer, status);
	CHECK_STATUS_VOID(*status);

	// Parse XML code in the xmlbuffer using the expat library.
	// Get an XML_Parser object.
	XML_Parser parser = XML_ParserCreate(NULL);
	if (NULL == parser) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("could not allocate memory for XML parser");
		return;
	}

	// Set data that is passed to the handler functions.
	struct XMLParseData xmlparsedata = { .inst = inst, .seed = seed, .status =
			EXIT_SUCCESS };
	XML_SetUserData(parser, &xmlparsedata);

	// Set the handler functions.
	XML_SetElementHandler(parser, GenInstXMLElementStart, GenInstXMLElementEnd);

	// Parse all the data in the string buffer.
	const int done = 1;
	if (!XML_Parse(parser, xmlbuffer->text, strlen(xmlbuffer->text), done)) {
		// Parse error.
		*status = EXIT_FAILURE;
		char msg[MAXMSG];
		sprintf(msg, "failed parsing XML file '%s':\n%s\n", filename,
				XML_ErrorString(XML_GetErrorCode(parser)));
		printf("%s", xmlbuffer->text);
		SIXT_ERROR(msg);
		return;
	}
	// Check for errors.
	if (EXIT_SUCCESS != xmlparsedata.status) {
		*status = xmlparsedata.status;
		return;
	}

	// Release memory.
	XML_ParserFree(parser);

	// Remove the XML string buffer.
	freeXMLBuffer(&xmlbuffer);

	// Check if all required parameters have been read successfully from
	// the XML file.
	if (INT_MAX == inst->det->pixgrid->xwidth) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("no specification of x-width of pixel array");
		return;
	}
	if (INT_MAX == inst->det->pixgrid->ywidth) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("no specification of y-width of pixel array");
		return;
	}

	if (isnan(inst->det->pixgrid->xrpix)) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("no specification of x reference pixel");
		return;
	}
	if (isnan(inst->det->pixgrid->yrpix)) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("no specification of y reference pixel");
		return;
	}

	if (isnan(inst->det->pixgrid->xdelt)) {
		*status = EXIT_FAILURE;
		SIXT_WARNING("no specification of pixel x-width");
		return;
	}
	if (isnan(inst->det->pixgrid->ydelt)) {
		*status = EXIT_FAILURE;
		SIXT_WARNING("no specification of pixel y-width");
		return;
	}

	if (inst->det->pixgrid->xborder < 0.) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("invalid specification of x-border of pixels");
		return;
	}
	if (inst->det->pixgrid->yborder < 0.) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("invalid specification of y-border of pixels");
		return;
	}

	if (NULL == inst->det->rmf) {
		SIXT_WARNING("no specification of response file (RMF/RSP)");
	}
	if (NULL == inst->tel->arf) {
		SIXT_WARNING("no specification of ARF");
	}

	if (NULL == inst->tel->psf) {
		SIXT_WARNING("no specification of PSF");
	}

	if (0. == inst->tel->focal_length) {
		SIXT_WARNING("no specification of the focal length of the telescope");
	}
	if (0. == inst->tel->fov_diameter) {
		SIXT_WARNING("no specification of the diameter of the telescope FoV");
	}

	if (0 == inst->det->readout_trigger) {
		SIXT_WARNING("no specification of the readout trigger");
	}

	if (GS_EXPONENTIAL == inst->det->split->type) {
		if (inst->det->split->par1 == 0.) {
			*status = EXIT_FAILURE;
			SIXT_ERROR("no valid split model parameters in the XML file");
			return;
		}
	}

	if (GS_GAUSS == inst->det->split->type) {
		if ((inst->det->split->par1 == 0.) && (inst->det->split->par2 == 0.)) {
			*status = EXIT_FAILURE;
			SIXT_ERROR("no valid split model parameters in the XML file");
			return;
		}
	}

	// change borders for event driven detectors

	if (GENDET_TIME_TRIGGERED != inst->det->readout_trigger) {
		inst->det->rawymin = 0;
		inst->det->rawymax = inst->det->pixgrid->ywidth - 1;
	}

	// END of checking, if all detector parameters have successfully been
	// read from the XML file.
}

GenInst* loadGenInst(const char* const filename, const unsigned int seed,
		int* const status) {
	// Get a new and empty data structure.
	GenInst* inst = newGenInst(status);
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
	char* lastslash = strrchr(rootname, '/');
	if (NULL == lastslash) {
		inst->filepath = (char*) malloc(sizeof(char));
		CHECK_NULL_RET(inst->filepath, *status,
				"memory allocation for filepath failed", inst);
		inst->filename = (char*) malloc((strlen(rootname) + 1) * sizeof(char));
		CHECK_NULL_RET(inst->filename, *status,
				"memory allocation for filename failed", inst);
		strcpy(inst->filepath, "");
		strcpy(inst->filename, rootname);
	} else {
		lastslash++;
		inst->filename = (char*) malloc((strlen(lastslash) + 1) * sizeof(char));
		CHECK_NULL_RET(inst->filename, *status,
				"memory allocation for filename failed", inst);
		strcpy(inst->filename, lastslash);

		*lastslash = '\0';
		inst->filepath = (char*) malloc((strlen(rootname) + 1) * sizeof(char));
		CHECK_NULL_RET(inst->filepath, *status,
				"memory allocation for filepath failed", inst);
		strcpy(inst->filepath, rootname);
	}
	// END of storing the filename and filepath.

	// Read in the XML definition of the detector.
	parseGenInstXML(inst, filename, seed, status);
	CHECK_STATUS_RET(*status, inst);

	// Allocate memory for the detector pixels.
	inst->det->line = (GenDetLine**) malloc(
			inst->det->pixgrid->ywidth * sizeof(GenDetLine*));
	if (NULL == inst->det->line) {
		*status = EXIT_FAILURE;
		SIXT_ERROR("memory allocation for GenDet pixel array failed");
		return (inst);
	}
	int ii;
	for (ii = 0; ii < inst->det->pixgrid->ywidth; ii++) {
		inst->det->line[ii] = newGenDetLine(inst->det->pixgrid->xwidth, status);
		if (EXIT_SUCCESS != *status)
			return (inst);
	}

	return (inst);
}

static void GenInstXMLElementStart(void* parsedata, const char* el,
		const char** attr) {
	struct XMLParseData* xmlparsedata = (struct XMLParseData*) parsedata;

	// Check if an error has occurred previously.
	CHECK_STATUS_VOID(xmlparsedata->status);

	// Convert the element to an upper case string.
	char Uelement[MAXMSG];
	strcpy(Uelement, el);
	strtoupper(Uelement);

	// Check for different elements.
	if (!strcmp(Uelement, "INSTRUMENT")) {
		// Determine the values of TELESCOP and INSTRUME.
		char telescop[MAXMSG];
		getXMLAttributeString(attr, "TELESCOP", telescop);
		xmlparsedata->inst->telescop = (char*) malloc(
				(strlen(telescop) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->telescop, xmlparsedata->status,
				"memory allocation for TELESCOP failed");
		strcpy(xmlparsedata->inst->telescop, telescop);

		char instrume[MAXMSG];
		getXMLAttributeString(attr, "INSTRUME", instrume);
		xmlparsedata->inst->instrume = (char*) malloc(
				(strlen(instrume) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->instrume, xmlparsedata->status,
				"memory allocation for INSTRUME failed");
		strcpy(xmlparsedata->inst->instrume, instrume);

	} else if (!strcmp(Uelement, "LINESHIFT")) {
		CLLineShift* cllineshift = newCLLineShift(&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);
		append2ClockList(xmlparsedata->inst->det->clocklist, CL_LINESHIFT,
				cllineshift, &xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "NEWFRAME")) {
		CLNewFrame* clnewframe = newCLNewFrame(&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);
		append2ClockList(xmlparsedata->inst->det->clocklist, CL_NEWFRAME,
				clnewframe, &xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "READOUTLINE")) {

		int lineindex = getXMLAttributeInt(attr, "LINEINDEX");
		if (lineindex < 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("negative index for readout line");
			return;
		}
		int readoutindex = getXMLAttributeInt(attr, "READOUTINDEX");
		if (readoutindex < 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("negative index for readout line");
			return;
		}
		if (readoutindex > xmlparsedata->inst->det->rawymax) {
			xmlparsedata->inst->det->rawymax = readoutindex;
		}
		if (readoutindex < xmlparsedata->inst->det->rawymin) {
			xmlparsedata->inst->det->rawymin = readoutindex;
		}
		CLReadoutLine* clreadoutline = newCLReadoutLine(lineindex, readoutindex,
				&xmlparsedata->status);
		append2ClockList(xmlparsedata->inst->det->clocklist, CL_READOUTLINE,
				clreadoutline, &xmlparsedata->status);

	} else if (!strcmp(Uelement, "DIMENSIONS")) {

		xmlparsedata->inst->det->pixgrid->xwidth = getXMLAttributeInt(attr,
				"XWIDTH");
		xmlparsedata->inst->det->pixgrid->ywidth = getXMLAttributeInt(attr,
				"YWIDTH");

	} else if (!strcmp(Uelement, "WCS")) {

		xmlparsedata->inst->det->pixgrid->xrpix = getXMLAttributeFloat(attr,
				"XRPIX");
		xmlparsedata->inst->det->pixgrid->yrpix = getXMLAttributeFloat(attr,
				"YRPIX");
		xmlparsedata->inst->det->pixgrid->xrval = getXMLAttributeFloat(attr,
				"XRVAL");
		xmlparsedata->inst->det->pixgrid->yrval = getXMLAttributeFloat(attr,
				"YRVAL");
		xmlparsedata->inst->det->pixgrid->xdelt = getXMLAttributeFloat(attr,
				"XDELT");
		xmlparsedata->inst->det->pixgrid->ydelt = getXMLAttributeFloat(attr,
				"YDELT");
		xmlparsedata->inst->det->pixgrid->rota = getXMLAttributeFloat(attr,
				"ROTA") * M_PI / 180.;

	} else if (!strcmp(Uelement, "PIXELBORDER")) {

		xmlparsedata->inst->det->pixgrid->xborder = getXMLAttributeFloat(attr,
				"X");
		xmlparsedata->inst->det->pixgrid->yborder = getXMLAttributeFloat(attr,
				"Y");

	} else if (!strcmp(Uelement, "ARF")) {

		// Check if the ARF has been defined previously.
		if (NULL != xmlparsedata->inst->tel->arf) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("ARF already defined (cannot be loaded twice)");
			return;
		}

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("no file specified for ARF");
			return;
		}

		// Store the file name of the ARF.
		xmlparsedata->inst->tel->arf_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->tel->arf_filename,
				xmlparsedata->status,
				"memory allocation for ARF file name failed");
		strcpy(xmlparsedata->inst->tel->arf_filename, filename);

		// Load the ARF.
		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);
		xmlparsedata->inst->tel->arf = loadARF(filepathname,
				&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "PHA2PI")) {

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Store the file name of the Pha2Pi.
		xmlparsedata->inst->det->pha2pi_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->det->pha2pi_filename,
				xmlparsedata->status,
				"memory allocation for Pha2Pi file name failed");
		strcpy(xmlparsedata->inst->det->pha2pi_filename, filename);

	} else if (!strcmp(Uelement, "PIRMF")) {

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Store the file name of the PI RMF.
		xmlparsedata->inst->det->pirmf_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->det->pirmf_filename,
				xmlparsedata->status,
				"memory allocation for SIXTE RMF file name failed");
		strcpy(xmlparsedata->inst->det->pirmf_filename, filename);

	} else if (!strcmp(Uelement, "SPECARF")) {

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Store the file name of the PI RMF.
		xmlparsedata->inst->det->specarf_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->det->specarf_filename,
				xmlparsedata->status,
				"memory allocation for SPEC ARF file name failed");
		strcpy(xmlparsedata->inst->det->specarf_filename, filename);

	} else if (!strcmp(Uelement, "RMF")) {

		// Check if the RMF has been defined previously.
		if (NULL != xmlparsedata->inst->det->rmf) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("RMF already defined (cannot be loaded twice)");
			return;
		}

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("no file specified for RMF");
			return;
		}

		// Store the file name of the RMF.
		xmlparsedata->inst->det->rmf_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->det->rmf_filename,
				xmlparsedata->status,
				"memory allocation for RMF file name failed");
		strcpy(xmlparsedata->inst->det->rmf_filename, filename);

		// Load the RMF.
		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);
		xmlparsedata->inst->det->rmf = loadNormalizedRMF(filepathname,
				&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "RSP")) {

		// Check if the ARF or RMF have been defined previously.
		if ((NULL != xmlparsedata->inst->tel->arf)
				|| (NULL != xmlparsedata->inst->det->rmf)) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("ARF or RMF already defined (cannot be loaded twice)");
			return;
		}

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			SIXT_ERROR("no file specified for RSP");
			xmlparsedata->status = EXIT_FAILURE;
			return;
		}

		// Store the file name of the RSP.
		xmlparsedata->inst->tel->arf_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->tel->arf_filename,
				xmlparsedata->status,
				"memory allocation for ARF file name failed");
		strcpy(xmlparsedata->inst->tel->arf_filename, filename);
		xmlparsedata->inst->det->rmf_filename = (char*) malloc(
				(strlen(filename) + 1) * sizeof(char));
		CHECK_NULL_VOID(xmlparsedata->inst->det->rmf_filename,
				xmlparsedata->status,
				"memory allocation for RMF file name failed");
		strcpy(xmlparsedata->inst->det->rmf_filename, filename);

		// Load the RSP
		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);
		loadArfRmfFromRsp(filepathname, &xmlparsedata->inst->tel->arf,
				&xmlparsedata->inst->det->rmf, &xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "PSF")) {

		// The focal length must be specified before load the PSF.
		// Check if this is the case.
		if (xmlparsedata->inst->tel->focal_length <= 0.) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("telescope focal length must be specified "
					"before loading the PSF");
			return;
		}
		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("no file specified for PSF");
			return;
		}

		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);
		xmlparsedata->inst->tel->psf = newPSF(filepathname,
				xmlparsedata->inst->tel->focal_length, &xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "VIGNETTING")) {

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("no file specified for vignetting");
			return;
		}

		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);
		xmlparsedata->inst->tel->vignetting = newVignetting(filepathname,
				&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "FOCALLENGTH")) {

		xmlparsedata->inst->tel->focal_length = getXMLAttributeFloat(attr,
				"VALUE");

	} else if (!strcmp(Uelement, "FOV")) {

		xmlparsedata->inst->tel->fov_diameter = getXMLAttributeFloat(attr,
				"DIAMETER") * M_PI / 180.;

	} else if (!strcmp(Uelement, "CTE")) {

		xmlparsedata->inst->det->cte = getXMLAttributeFloat(attr, "VALUE");

	} else if (!strcmp(Uelement, "BADPIXMAP")) {

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("no file specified for bad pixel map");
			return;
		}

		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);
		xmlparsedata->inst->det->badpixmap = loadBadPixMap(filepathname,
				&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "AUXBACKGROUND")) {

		if (0 == auxBkgInitialized) {
			// Determine the file containing the simulated background
			// hits from GEANT4.
			char filename[MAXFILENAME];
			getXMLAttributeString(attr, "FILENAME", filename);

			// Check if a file name has been specified.
			if (strlen(filename) == 0) {
				xmlparsedata->status = EXIT_FAILURE;
				SIXT_ERROR("no file specified for aux detector background");
				return;
			}

			char filepathname[MAXFILENAME];
			strcpy(filepathname, xmlparsedata->inst->filepath);
			strcat(filepathname, filename);

			bkgAux rateinfo;
			rateinfo.rate = getXMLAttributeDouble(attr, "RATE");

			bkgInitializeAux(filepathname, xmlparsedata->seed, &rateinfo,
					&xmlparsedata->status);
			CHECK_STATUS_VOID(xmlparsedata->status);
			auxBkgInitialized = 1;

			// Load the optional light curve if available. The light curve
			// defines the time-variability of the background flux.
			getXMLAttributeString(attr, "LIGHTCURVE", filename);

			// Check if a light curve has been specified.
			if (strlen(filename) > 0) {
				strcpy(filepathname, xmlparsedata->inst->filepath);
				strcat(filepathname, filename);
				bkgSetRateFct(filepathname, &xmlparsedata->status);
				CHECK_STATUS_VOID(xmlparsedata->status);
			}
		}
		xmlparsedata->inst->det->auxbackground = 1;

	} else if (!strcmp(Uelement, "PHABACKGROUND")) {

		char filename[MAXFILENAME];
		getXMLAttributeString(attr, "FILENAME", filename);

		// Check if a file name has been specified.
		if (strlen(filename) == 0) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR("no file specified for PHA detector background");
			return;
		}

		char filepathname[MAXFILENAME];
		strcpy(filepathname, xmlparsedata->inst->filepath);
		strcat(filepathname, filename);

		// There can be up to 2 PHA background models.
		int ii;
		for (ii = 0; ii < 2; ii++) {
			if (NULL == xmlparsedata->inst->det->phabkg[ii]) {
				break;
			}
		}
		if (MAX_PHABKG == ii) {
			xmlparsedata->status = EXIT_FAILURE;
			SIXT_ERROR(
					"cannot use more than the given number of PHA background models");
			return;
		}
		xmlparsedata->inst->det->phabkg[ii] = newPHABkg(filepathname,
                                xmlparsedata->inst->det->rmf,
				&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

		// If needed, link the vignetting function and the focal length
		// of the telescope.
		if (getXMLAttributeInt(attr, "VIGNETTING") != 0) {
			xmlparsedata->inst->det->phabkg[ii]->vignetting =
					&xmlparsedata->inst->tel->vignetting;
			xmlparsedata->inst->det->phabkg[ii]->focal_length =
					&xmlparsedata->inst->tel->focal_length;
		}

		// Load the optional light curve if available. The light curve
		// defines the time-variability of the background flux.
		getXMLAttributeString(attr, "LIGHTCURVE", filename);

		// Check if a light curve has been specified.
		if (strlen(filename) > 0) {
			strcpy(filepathname, xmlparsedata->inst->filepath);
			strcat(filepathname, filename);
			bkgSetRateFct(filepathname, &xmlparsedata->status);
			CHECK_STATUS_VOID(xmlparsedata->status);
		}

	} else if (!strcmp(Uelement, "MXS")) {
                // Get filename
                char filename[MAXFILENAME];
        	getXMLAttributeString(attr, "FILENAME", filename);

                // Check if a file name has been specified.
        	if (strlen(filename) == 0) {
        		xmlparsedata->status = EXIT_FAILURE;
        		SIXT_ERROR("no filename specified for MXS");
        		return;
        	}

                char filepathname[MAXFILENAME];
        	strcpy(filepathname, xmlparsedata->inst->filepath);
        	strcat(filepathname, filename);

                double frequency = getXMLAttributeDouble(attr, "FLASH_FREQUENCY");
                double flash_duration = getXMLAttributeDouble(attr, "FLASH_DURATION");
                double rate_det = getXMLAttributeDouble(attr, "RATE");

                int status = EXIT_SUCCESS;
                xmlparsedata->inst->det->mxs_params = loadMXSparams(frequency,
                                                        flash_duration, rate_det,
                                                        filepathname, &status);

                if (status != EXIT_SUCCESS) {
                    xmlparsedata->status = EXIT_FAILURE;
        			SIXT_ERROR("failed to load MXS parameters");
        			return;
                }
        } else if (!strcmp(Uelement, "SPLIT")) {

		char type[MAXMSG];
		getXMLAttributeString(attr, "TYPE", type);
		strtoupper(type);
		if (!strcmp(type, "NONE")) {
			xmlparsedata->inst->det->split->type = GS_NONE;
		} else if (!strcmp(type, "GAUSS")) {
			xmlparsedata->inst->det->split->type = GS_GAUSS;
		} else if (!strcmp(type, "EXPONENTIAL")) {
			xmlparsedata->inst->det->split->type = GS_EXPONENTIAL;
		}
		xmlparsedata->inst->det->split->par1 = getXMLAttributeFloat(attr,
				"PAR1");
		xmlparsedata->inst->det->split->par2 = getXMLAttributeFloat(attr,
				"PAR2");

	} else if (!strcmp(Uelement, "DEPFET")) {

		xmlparsedata->inst->det->depfet.depfetflag = 1;

		xmlparsedata->inst->det->depfet.t_integration = getXMLAttributeDouble(
				attr, "INTEGRATION");

		xmlparsedata->inst->det->depfet.t_clear = getXMLAttributeDouble(attr,
				"CLEAR");

		xmlparsedata->inst->det->depfet.t_settling = getXMLAttributeDouble(attr,
				"SETTLING");

		char type[MAXMSG];
		getXMLAttributeString(attr, "CLEAR_FCN", type);
		strtoupper(type);

		if (!strcmp(type, "EXPONENTIAL")) {

			xmlparsedata->inst->det->depfet.clear_const = (double*) malloc(
					sizeof(double));
			if (xmlparsedata->inst->det->depfet.clear_const == NULL) {
				xmlparsedata->status = EXIT_FAILURE;
				SIXT_ERROR("memory allocation error for depfet constants.");
				return;
			}
			xmlparsedata->inst->det->depfet.clear_const[0] =
					getXMLAttributeDouble(attr, "CLEAR_TAU");
			xmlparsedata->inst->det->depfet.clear_fcn =
					&depfet_get_exponential_clear_signal;
		} else {

			xmlparsedata->inst->det->depfet.clear_const = (double*) malloc(
					sizeof(double));
			if (xmlparsedata->inst->det->depfet.clear_const == NULL) {
				xmlparsedata->status = EXIT_FAILURE;
				SIXT_ERROR("memory allocation error for depfet constants.");
				return;
			}
			xmlparsedata->inst->det->depfet.clear_const[0] =
					xmlparsedata->inst->det->depfet.t_clear;

			xmlparsedata->inst->det->depfet.clear_fcn =
					&depfet_get_linear_clear_signal;
		}

		getXMLAttributeString(attr, "TYPE", type);
		strtoupper(type);
		if (!strcmp(type, "NORMAL")) {
			xmlparsedata->inst->det->depfet.istorageflag = 0;

		} else if (!strcmp(type, "IS")) {
			xmlparsedata->inst->det->depfet.istorageflag = 1;

			xmlparsedata->inst->det->depfet.t_transfer = getXMLAttributeDouble(
					attr, "TRANSFER");
		}

	} else if (!strcmp(Uelement, "READOUT")) {

		char mode[MAXMSG];
		getXMLAttributeString(attr, "MODE", mode);
		strtoupper(mode);
		if (!strcmp(mode, "TIME")) {
			xmlparsedata->inst->det->readout_trigger = GENDET_TIME_TRIGGERED;
		} else if (!strcmp(mode, "EVENT")) {
			xmlparsedata->inst->det->readout_trigger = GENDET_EVENT_TRIGGERED;
		}
		xmlparsedata->inst->det->deadtime = getXMLAttributeDouble(attr,
				"DEADTIME");

	} else if (!strcmp(Uelement, "WAIT")) {

		float waittime = getXMLAttributeFloat(attr, "TIME");
		CLWait* clwait = newCLWait(waittime, &xmlparsedata->status);
		append2ClockList(xmlparsedata->inst->det->clocklist, CL_WAIT, clwait,
				&xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

		// Accumulate the amount of time required for one read-out frame.
		xmlparsedata->inst->det->frametime += waittime;

	} else if (!strcmp(Uelement, "CLEARLINE")) {

		int lineindex = getXMLAttributeInt(attr, "LINEINDEX");
		CLClearLine* clclearline = newCLClearLine(lineindex,
				&xmlparsedata->status);
		append2ClockList(xmlparsedata->inst->det->clocklist, CL_CLEARLINE,
				clclearline, &xmlparsedata->status);
		CHECK_STATUS_VOID(xmlparsedata->status);

	} else if (!strcmp(Uelement, "THRESHOLD_READOUT_LO_KEV")) {

		xmlparsedata->inst->det->threshold_readout_lo_keV =
				getXMLAttributeFloat(attr, "VALUE");
		headas_chat(5, "lower readout threshold: %.3lf keV\n",
				xmlparsedata->inst->det->threshold_readout_lo_keV);

	} else if (!strcmp(Uelement, "THRESHOLD_PATTERN_UP_KEV")) {

		xmlparsedata->inst->det->threshold_pattern_up_keV =
				getXMLAttributeFloat(attr, "VALUE");
		headas_chat(5, "upper pattern threshold: %.3lf keV\n",
				xmlparsedata->inst->det->threshold_pattern_up_keV);

	} else if (!strcmp(Uelement, "THRESHOLD_EVENT_LO_KEV")) {

		xmlparsedata->inst->det->threshold_event_lo_keV = getXMLAttributeFloat(
				attr, "VALUE");
		headas_chat(5, "lower event threshold: %.3lf keV\n",
				xmlparsedata->inst->det->threshold_event_lo_keV);

	} else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_KEV")) {

		xmlparsedata->inst->det->threshold_split_lo_keV = getXMLAttributeFloat(
				attr, "VALUE");
		headas_chat(5, "lower split threshold: %.3lf keV\n",
				xmlparsedata->inst->det->threshold_split_lo_keV);

	} else if (!strcmp(Uelement, "THRESHOLD_SPLIT_LO_FRACTION")) {

		xmlparsedata->inst->det->threshold_split_lo_fraction =
				getXMLAttributeFloat(attr, "VALUE");
		headas_chat(5, "lower split threshold: %.1lf %%\n",
				xmlparsedata->inst->det->threshold_split_lo_fraction * 100.);

	} else if (!strcmp(Uelement, "TELESCOPE")) {

		// Nothing to do here. Do not remove this selection! Otherwise
		// the tag <telescope> will be regarded as unknown.

	} else if (!strcmp(Uelement, "DETECTOR")) {

		// Nothing to do here. Do not remove this selection! Otherwise
		// the tag <detector> will be regarded as unknown.

	} else {
		// Unknown tag, display warning.
		char msg[MAXMSG];
		sprintf(msg, "unknown XML tag: <%s>", el);
		SIXT_WARNING(msg);
	}
}

static void GenInstXMLElementEnd(void* parsedata, const char* el) {
	struct XMLParseData* xmlparsedata = (struct XMLParseData*) parsedata;

	(void) el;

	// Check if an error has occurred previously.
	CHECK_STATUS_VOID(xmlparsedata->status);
}
