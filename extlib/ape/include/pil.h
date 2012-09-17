/** \file pil.h
    \brief Declaration of Ape's PIL emulating interface.

    Ape is backward compatible with a subset of the PIL library developed by
    ISDC. Replacements are provided for the most commonly used API functions.
    Less common functions and functions which rely on details of pil's internal
    structures are not provided. Ape's implementation of these functions does
    not always behave identically to PIL's implementation, but the behavior should be
    close enough to provide compatibility.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_pil_h
#define ape_pil_h

#include "pil_error.h"

/* The following is needed for FILENAME_MAX. */
#include <stdio.h>

#define PIL_LINESIZE            (2000)
#define PIL_PATH_MAX            (FILENAME_MAX)
#define PIL_QUERY_DEFAULT       (1)
#define PIL_QUERY_OVERRIDE      (0)

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Perform standard initialization of local and/or system parameter
    files, using the PFILES environment variable.
    \param argc The number of command line arguments, including the name of executable.
    \param argv The command line arguments. The element argv[0] gives the name of
    the executable or parameter file.
*/
int PILInit(int argc, char ** argv);

/** \brief Perform standard close of all parameter files. If the status is 0,
    the local parameter file will be saved before closing.
    \param status Inherited status, presumably indicating whether a serious error
    occurred.
*/
int PILClose(int status);

/** \brief Query for the given parameter, and return its value converted to a bool,
    expressed as an int with 1 for true/yes and 0 for false/no.

    If the value from the user/parameter file is invalid, or if it cannot be
    converted to a bool, this function returns a non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable.
*/
int PILGetBool(const char * name, int * result);

/** \brief Query for the given parameter, and return its value converted to an int.

    If the value from the user/parameter file is invalid, or if it cannot be
    converted to a int, this function returns a non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable.
*/
int PILGetInt(const char * name, int * result);

/** \brief Query for the given parameter, and return its value converted to a long.

    If the value from the user/parameter file is invalid, or if it cannot be
    converted to a long, this function returns a non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable.
*/
int PILGetLong(const char * name, long * result);

/** \brief Query for the given parameter, and return its value converted to a double.

    If the value from the user/parameter file is invalid, or if it cannot be
    converted to a double, this function returns a non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable.
*/
int PILGetReal(const char * name, double * result);

/** \brief Query for the given parameter, and return its value converted to a float.

    If the value from the user/parameter file is invalid, or if it cannot be
    converted to a float, this function returns a non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable.
*/
int PILGetReal4(const char * name, float * result);

/** \brief Query for the given parameter, and return its value converted to an string.

    If the value from the user/parameter file is invalid, this function returns a
    non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable. The length of the output buffer
    must be at least PIL_LINESIZE characters. The client is responsible for freeing
    any dynamically allocated memory.
*/
int PILGetString(const char * name, char * result);

/** \brief Query for the given parameter, and return its value converted to a file
    name expressed as a string.

    If the value from the user/parameter file is invalid, this function returns a
    non-0 error code.
    \param name The name of the parameter.
    \param result Pointer to the output variable. The length of the output buffer
    must be at least PIL_LINESIZE characters. The client is responsible for freeing
    any dynamically allocated memory.
*/
int PILGetFname(const char * name, char * result);

/** \brief Query for the given parameter, and return its value converted to an string.

    If the value from the user/parameter file is invalid, this function returns a
    non-0 error code. Note that in Ape's implementation, this is identical to
    PILGetString.
    \param name The name of the parameter.
    \param result Pointer to the output variable. The length of the output buffer
    must be at least PIL_LINESIZE characters. The client is responsible for freeing
    any dynamically allocated memory.
*/
int PILGetAsString(const char * name, char * result);

/** \brief Set the value of the given parameter.
    \param name The name of the parameter.
    \param value The new value for that parameter.
*/
int PILPutBool(const char * name, int value);

/** \brief Set the value of the given parameter.
    \param name The name of the parameter.
    \param value The new value for that parameter.
*/
int PILPutInt(const char * name, int value);

/** \brief Set the value of the given parameter.
    \param name The name of the parameter.
    \param value The new value for that parameter.
*/
int PILPutReal(const char * name, double value);

/** \brief Set the value of the given parameter.
    \param name The name of the parameter.
    \param value The new value for that parameter.
*/
int PILPutString(const char * name, const char * value);

/** \brief Set the value of the given parameter.
    \param name The name of the parameter.
    \param value The new value for that parameter.
*/
int PILPutFname(const char * name, const char * value);

/** \brief Get the current query mode.
    \param query_mode The mode of query, either PIL_QUERY_OVERRIDE or PIL_QUERY_DEFAULT.
*/
int PILGetQueryMode(void);

/** \brief Set query mode.
    \param new_mode The new mode to set. If it is PIL_QUERY_OVERRIDE, no prompts will occur. If it is
    PIL_QUERY_DEFAULT, the usual prompts will occur.
*/
int PILOverrideQueryMode(int new_mode);

/** \brief Set function used for output.
    \param func The function to use for standard output.
*/
int PILSetLoggerFunction(int (*func)(char *));

/** \brief Use the given function to check file access.
    \param func The function to use.
*/
int PILSetFileAccessFunction(int (*func)(const char *, const char *));

/** \brief Set reprompt mode. If enabled (reprompt != 0), the given parameter will be prompted for even if it
    was supplied on the command line, is hidden, or was already prompted for, (i.e. always prompted). */
int PILSetReprompt(const char * par_name, int reprompt);

/** \brief Save current parameters. */
int PILFlushParameters(void);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: pil.h,v $
 * Revision 1.8  2009/12/04 19:15:26  peachey
 * Add PILGetLong function.
 *
 * Revision 1.7  2006/06/06 17:47:12  peachey
 * Add cvs log info.
 *
 */
