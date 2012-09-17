/** \file ape_trad.h
    \brief Declaration of traditional XPI/PIL-compliant interface.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_trad_h
#define ape_ape_trad_h

#include "ape/ape_io.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Perform traditional initialization: open and merge user and/or system parameter files for
    the task named in the first command line argument, and interpret remaining arguments as positional
    and/or explicit parameter values. Return eOK (0) if successful.
    \param argc The number of command line arguments.
    \param argv The command line arguments. argv[0] is the name of the executable.
*/
int ape_trad_init(int argc, char ** argv);

/** \brief Perform shutdown of all traditional facilities.
    \param save_flag Flag indicating whether to save (1) or just forget (0) current parameters.
*/
int ape_trad_close(int save_flag);

/** \brief Get the current parameter file.
    \param par_file Pointer to the output file pointer.
*/
int ape_trad_get_current(ApeParFile ** par_file);

/** \brief Get the system parameter file.
    \param par_file Pointer to the output file pointer.
*/
int ape_trad_get_sys_pfile(ApeParFile ** par_file);

/** \brief Get an array containing the names of all named parameters in the file.
    \param par_name Pointer to output array containing the names. Client must free this using ape_util_free_string_array.
    The array is null terminated.
*/
int ape_trad_get_par_names(char *** par_name);

/** \brief Get the type of the given parameter.
    \param par_name The name of the parameter.
    \param par_type The parameter type, e.g. "i", "fr" etc. Client needs to allocate/manage space for at least
    APE_PAR_TYPE_CODE_LEN characters. The first character is guaranteed to specify the type, either b, d, f, i, r, s.
*/
int ape_trad_get_type(const char * par_name, char * par_type);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a boolean,
    stored as a character holding 1 for true and 0 for false.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_bool(const char * par_name, char * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a double.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_double(const char * par_name, double * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a float.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_float(const char * par_name, float * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a file name.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_file_name(const char * par_name, char ** value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as an int.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_int(const char * par_name, int * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a long int.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_long(const char * par_name, long * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a string.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_trad_query_string(const char * par_name, char ** value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a string.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
    \param case_code Code indicating how to handle the case of the string. Can be eUserCase, eLowerCase, eUpperCase, eEnumCase
*/
int ape_trad_query_string_case(const char * par_name, char ** value, char case_code);

/** \brief Get the named parameter interpreted as a boolean, stored as a character holding 1 for true and 0 for false.
    Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_bool(const char * par_name, char * value);

/** \brief Get the named parameter interpreted as a double. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_double(const char * par_name, double * value);

/** \brief Get the named parameter interpreted as a float. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_float(const char * par_name, float * value);

/** \brief Get the named parameter interpreted as a file name. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_file_name(const char * par_name, char ** value);

/** \brief Get the named parameter as an int. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_int(const char * par_name, int * value);

/** \brief Get the named parameter as a long int. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_long(const char * par_name, long * value);

/** \brief Get the named parameter as a string. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_trad_get_string(const char * par_name, char ** value);

/** \brief Get the named parameter as a string. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
    \param case_code Code indicating how to handle the case of the string. Can be eUserCase, eLowerCase, eUpperCase, eEnumCase
*/
int ape_trad_get_string_case(const char * par_name, char ** value, char case_code);

/** \brief Set the value of the named parameter, as a bool.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_bool(const char * par_name, char value);

/** \brief Set the value of the named parameter, as a double.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_double(const char * par_name, double value);

/** \brief Set the value of the named parameter, as a float.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_float(const char * par_name, float value);

/** \brief Set the value of the named parameter, as a string representing a file name.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_file_name(const char * par_name, const char * value);

/** \brief Set the value of the named parameter, as an integer.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_int(const char * par_name, int value);

/** \brief Set the value of the named parameter, as a long integer.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_long(const char * par_name, long value);

/** \brief Set the value of the named parameter, as a short integer.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_short(const char * par_name, short value);

/** \brief Set the value of the named parameter, as a string.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_trad_set_string(const char * par_name, const char * value);

/** \brief Save currently defined parameter set, respecting modes of all parameters and default mode of file.
*/
int ape_trad_save(void);

/** \brief Find the named parameter object.
    \param par_name The name of the parameter.
    \param par Pointer to output parameter object pointer.
*/
int ape_trad_find_par(const char * par_name, ApePar ** par);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_trad.h,v $
 * Revision 1.19  2010/11/12 20:47:44  irby
 * Add ape_trad_set_short() and cfortran macros for ape_trad_set_double,
 * ape_trad_set_long, and ape_trad_set_short.  Also, add string & filename
 * cfortran wrapper functions (including string handling code).
 *
 * Revision 1.18  2010/09/10 19:08:34  irby
 * Add ape_trad_get_byte() wrapper for ape_trad_get_bool() as cfortran has no
 * char, only byte (signed char).  Add ape_trad_get_int() & ape_trad_set_int()
 * (also ape_session_get_int() & ape_session_set_int() and related tests) for
 * Fortran wrappers.
 *
 * Revision 1.17  2007/01/25 17:49:57  peachey
 * Fixed bug in ape_par_get_string_case that caused it not to behave
 * correctly when called for a non-enumerated parameter. Added functions
 * ape_trad_prompt_string_case and ape_trad_get_string_case for getting
 * parameters with control over case of returned string.
 *
 * Revision 1.16  2006/11/30 20:38:37  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.15  2006/11/24 19:47:27  peachey
 * Add ape_trad_query_float, ape_trad_get_float and ape_trad_set_float.
 *
 * Revision 1.14  2006/06/06 13:30:19  peachey
 * Add ape_trad_find_par, for getting underlying ape parameter object.
 *
 * Revision 1.13  2006/06/05 18:59:33  peachey
 * Add apd_trad_get_type for getting the type code of a parameter.
 *
 * Revision 1.12  2006/06/03 01:27:24  peachey
 * Add comment about null termination of array.
 *
 * Revision 1.11  2006/05/31 03:21:43  peachey
 * Rationalize ape_par_get_* family of functions, and use them to implement
 * ape_trad_get_* and ape_trad_query_*.
 *
 * Revision 1.10  2006/05/26 15:22:10  peachey
 * Add ape_trad_set* functions complete ape_trad_query* functions, for completeness
 * and to support PIL compatibility.
 *
 * Revision 1.9  2006/05/23 19:31:41  peachey
 * Add ape_trad_get_sys_pfile, for getting the system parameter file.
 *
 * Revision 1.8  2006/05/22 17:35:41  peachey
 * Add ape_trad_query_string, for prompting for and getting string parameters.
 *
 * Revision 1.7  2006/05/18 14:00:48  peachey
 * Rename ape_trad_get_bool ape_trad_query_bool. Add ape_trad_get_string
 * which only gets string, without prompting.
 *
 * Revision 1.6  2006/05/17 03:15:29  peachey
 * Add ape_trad_get_current, for getting the current parameter file object.
 *
 * Revision 1.5  2006/05/17 02:20:23  peachey
 * Add ape_trad_get_par_names, for getting all named parameters in a file.
 *
 * Revision 1.4  2006/05/13 18:45:57  peachey
 * Add ape_trad_save and use it to save parameter file after each prompt.
 *
 * Revision 1.3  2006/05/12 03:32:07  peachey
 * Add ape_trad_get_bool, and infrastructure for all the ape_trad_get_... family.
 *
 * Revision 1.2  2006/04/19 16:04:18  peachey
 * Add rational shutdown/cleanup/atexit stuff. Public function
 * ape_trad_close performs all shutdown/cleanup. Static function ape_trad_destroy
 * does all deallocations. Static function ape_trad_atexit is registered with
 * ape_util_atexit, and calls ape_trad_close. System and local par files are now
 * kept in static storage by ape_trad_init.
 *
 * Revision 1.1  2006/04/12 14:18:13  peachey
 * Add beginning of traditional facility for duplicating PIL/XPI behavior.
 *
*/
