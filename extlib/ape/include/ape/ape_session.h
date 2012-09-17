/** \file ape_session.h
    \brief Declaration of high-level interface to an entire interactive session, including system
      and local paramter files, and all other parameter file objects needed to perform the various
      required operations.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_session_h
#define ape_ape_session_h

#include "ape/ape_io.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ApeSession;
typedef struct ApeSession ApeSession;

/** \brief Perform initialization of session object: open and merge user and/or system parameter files for
    the task named in the second command line argument, and interpret remaining arguments as positional
    and/or explicit parameter values. Return eOK (0) if successful.
    \param session The session containing parameters.
    \param argc The number of command line arguments.
    \param argv The command line arguments. argv[0] is the name of the executable.
*/
int ape_session_init(ApeSession ** session, int argc, char ** argv);

/** \brief Perform shutdown of a session.
    \param session The session containing parameters.
    \param save_flag Flag indicating whether to save (1) or just forget (0) current parameters.
*/
int ape_session_close(ApeSession * session, int save_flag);

/** \brief Get the current parameter file.
    \param session The session containing parameters.
    \param par_file Pointer to the output file pointer.
*/
int ape_session_get_current(ApeSession * session, ApeParFile ** par_file);

/** \brief Get the system parameter file.
    \param session The session containing parameters.
    \param par_file Pointer to the output file pointer.
*/
int ape_session_get_sys_pfile(ApeSession * session, ApeParFile ** par_file);

/** \brief Get an array containing the names of all named parameters in the file.
    \param session The session containing parameters.
    \param par_name Pointer to output array containing the names. Client must free this using ape_util_free_string_array.
    The array is null terminated.
*/
int ape_session_get_par_names(ApeSession * session, char *** par_name);

/** \brief Get the type of the given parameter.
    \param session The session containing parameters.
    \param par_name The name of the parameter.
    \param par_type The parameter type, e.g. "i", "fr" etc. Client needs to allocate/manage space for at least
    APE_PAR_TYPE_CODE_LEN characters. The first character is guaranteed to specify the type, either b, d, f, i, r, s.
*/
int ape_session_get_type(ApeSession * session, const char * par_name, char * par_type);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a boolean,
    stored as a character holding 1 for true and 0 for false.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_bool(ApeSession * session, const char * par_name, char * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a double.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_double(ApeSession * session, const char * par_name, double * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a float.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_float(ApeSession * session, const char * par_name, float * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a file name.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_file_name(ApeSession * session, const char * par_name, char ** value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as an int.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_int(ApeSession * session, const char * par_name, int * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a long int.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_long(ApeSession * session, const char * par_name, long * value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a string.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
*/
int ape_session_query_string(ApeSession * session, const char * par_name, char ** value);

/** \brief Prompt for the named parameter if necessary, then interpret the parameter value as a string.
    \param par_name The name of the parameter whose value is sought.
    \param value Pointer to the output variable to hold the converted value.
    \param case_code Code indicating how to handle the case of the string. Can be eUserCase, eLowerCase, eUpperCase, eEnumCase
*/
int ape_session_query_string_case(ApeSession * session, const char * par_name, char ** value, char case_code);

/** \brief Get the named parameter interpreted as a boolean, stored as a character holding 1 for true and 0 for false.
    Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_bool(ApeSession * session, const char * par_name, char * value);

/** \brief Get the named parameter interpreted as a double. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_double(ApeSession * session, const char * par_name, double * value);

/** \brief Get the named parameter interpreted as a float. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_float(ApeSession * session, const char * par_name, float * value);

/** \brief Get the named parameter interpreted as a file name. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_file_name(ApeSession * session, const char * par_name, char ** value);

/** \brief Get the named parameter as an int. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_int(ApeSession * session, const char * par_name, int * value);

/** \brief Get the named parameter as a long int. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_long(ApeSession * session, const char * par_name, long * value);

/** \brief Get the named parameter as a string. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
*/
int ape_session_get_string(ApeSession * session, const char * par_name, char ** value);

/** \brief Get the named parameter as a string. Do not prompt. Error checking is still performed.
    \param par_name The name of the parameter whose value is sought.
    \param value The value of the parameter.
    \param case_code Code indicating how to handle the case of the string. Can be eUserCase, eLowerCase, eUpperCase, eEnumCase
*/
int ape_session_get_string_case(ApeSession * session, const char * par_name, char ** value, char case_code);

/** \brief Set the value of the named parameter, as a bool.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_bool(ApeSession * session, const char * par_name, char value);

/** \brief Set the value of the named parameter, as a double.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_double(ApeSession * session, const char * par_name, double value);

/** \brief Set the value of the named parameter, as a float.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_float(ApeSession * session, const char * par_name, float value);

/** \brief Set the value of the named parameter, as a string representing a file name.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_file_name(ApeSession * session, const char * par_name, const char * value);

/** \brief Set the value of the named parameter, as an integer.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_int(ApeSession * session, const char * par_name, int value);

/** \brief Set the value of the named parameter, as a long integer.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_long(ApeSession * session, const char * par_name, long value);

/** \brief Set the value of the named parameter, as a short integer.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_short(ApeSession * session, const char * par_name, short value);

/** \brief Set the value of the named parameter, as a string.
    \param par_name The name of the parameter whose value is being set.
    \param value The new value for the parameter.
*/
int ape_session_set_string(ApeSession * session, const char * par_name, const char * value);

/** \brief Save currently defined parameter set, respecting modes of all parameters and default mode of file.
*/
int ape_session_save(ApeSession * session);

/** \brief Find the named parameter object.
    \param par_name The name of the parameter.
    \param par Pointer to output parameter object pointer.
*/
int ape_session_find_par(ApeSession * session, const char * par_name, ApePar ** par);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_session.h,v $
 * Revision 1.5  2010/11/12 20:49:15  irby
 * Add ape_session_set_short() (wraps to ape_session_set_int).
 *
 * Revision 1.4  2010/09/10 19:08:34  irby
 * Add ape_trad_get_byte() wrapper for ape_trad_get_bool() as cfortran has no
 * char, only byte (signed char).  Add ape_trad_get_int() & ape_trad_set_int()
 * (also ape_session_get_int() & ape_session_set_int() and related tests) for
 * Fortran wrappers.
 *
 * Revision 1.3  2010/06/02 20:03:57  peachey
 * Remove/replace usages of static global ApeSession. Instead,
 * use ape_session_init and ape_session_close to create sessions dynamically.
 *
 * Revision 1.2  2010/06/02 19:43:12  peachey
 * Add ApeSession argument as first parameter to all ape_session calls.
 *
 * Revision 1.1  2010/06/02 19:22:09  peachey
 * Add initial (partial) implementation of ape_session module.
 *
 */
