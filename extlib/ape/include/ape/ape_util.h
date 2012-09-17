/** \file ape_util.h
    \brief Declaration of internal utilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_util_h
#define ape_ape_util_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Register a function to be called when exit is called. (for internal use).
    \param func The function to call.
*/
int ape_util_atexit(void (*func)(void));

/** \brief Concatenate two strings. Return eOK if successful.
    \param s1 The first input string.
    \param s2 The second input string.
    \param result Pointer to the output string, which will be allocated to the correct size. Client is responsible for freeing it.
*/
int ape_util_cat_string(const char * s1, const char * s2, char ** result);

/** \brief Concatenate a directory name and a file name to produce a full file name.
    \param dir_name The name of the directory.
    \param file_name The name of the file.
    \param full_file_name Pointer to output full file name.
*/
int ape_util_append_file_name(const char * dir_name, const char * file_name, char ** full_file_name);

/** \brief Make a copy of a range of a string. Return eOK if successful.
    \param begin Pointer to the first characater in the original input string.
    \param end One past the the last character in the original input string.
    \param result Pointer to the output string, which will be allocated to the correct size. Client is responsible for freeing it.
*/
int ape_util_copy_range(const char * begin, const char * end, char ** result);

/** \brief Make a copy of a string. Return eOK if successful.
    \param orig The original input string.
    \param result Pointer to the output string, which will be allocated to the correct size. Client is responsible for freeing it.
*/
int ape_util_copy_string(const char * orig, char ** result);

/** \brief Free an array of strings.
    \param string_array The array to free (note: must be null terminated!)
*/
void ape_util_free_string_array(char ** string_array);

/** \brief Get the named environment variable, using supplied default value if the variable is not defined.
    \param name The name of the environment variable to get.
    \param value The output pointer to the returned value.
    \param def_value The default value assigned if variable is not defined.
*/
int ape_util_getenv(const char * name, char ** value, const char * def_value);

/** \brief Interpret environment variables to set up Ape internal data structures.
*/
int ape_util_interpret_env(void);

/** \brief Expand environment variables in the input string.
    \param input Input string which may contain environment variables.
    \param output The output string, consisting of the input with environment variables expanded.
*/
int ape_util_expand_env_var(const char * input, char ** output);

/** \brief Split a PFILES double path into its local and system parts.
    \param pfiles Input PFILES value.
    \param loc_pfiles Output local portion of PFILES. May be 0 in which case local part is not returned. Client must free.
    \param sys_pfiles Output system portion of PFILES. May be 0 in which case system part is not returned. Client must free.
*/
int ape_util_parse_pfiles(const char * pfiles, char ** loc_pfiles, char ** sys_pfiles);

/** \brief Convert a string to a bool value, stored as a char. Leading and trailing white space is ignored. In the
    event of a conversion error, the result of the failed conversion will still be stored in result, but the
    returned status will be non-0.
    \param input The input string to convert.
    \param result Pointer to variable in which to store the converted value.
*/
int ape_util_s2b(const char * input, char * result);

/** \brief Convert a string to a double value. Leading and trailing white space is ignored. In the
    event of a conversion error, the result of the failed conversion will still be stored in result, but the
    returned status will be non-0.
    \param input The input string to convert.
    \param result Pointer to variable in which to store the converted value.
*/
int ape_util_s2d(const char * input, double * result);

/** \brief Convert a string to a float value. Leading and trailing white space is ignored. In the
    event of a conversion error, the result of the failed conversion will still be stored in result, but the
    returned status will be non-0.
    \param input The input string to convert.
    \param result Pointer to variable in which to store the converted value.
*/
int ape_util_s2f(const char * input, float * result);

/** \brief Convert a string to an int value. Leading and trailing white space is ignored. In the
    event of a conversion error, the result of the failed conversion will still be stored in result, but the
    returned status will be non-0.
    \param input The input string to convert.
    \param result Pointer to variable in which to store the converted value.
*/
int ape_util_s2i(const char * input, int * result);

/** \brief Convert a string to a long value. Leading and trailing white space is ignored. In the
    event of a conversion error, the result of the failed conversion will still be stored in result, but the
    returned status will be non-0.
    \param input The input string to convert.
    \param result Pointer to variable in which to store the converted value.
*/
int ape_util_s2l(const char * input, long * result);

/** \brief Convert a string to a short value. Leading and trailing white space is ignored. In the
    event of a conversion error, the result of the failed conversion will still be stored in result, but the
    returned status will be non-0.
    \param input The input string to convert.
    \param result Pointer to variable in which to store the converted value.
*/
int ape_util_s2sh(const char * input, short * result);

/** \brief Perform string comparison, strcmp style, optionally case insensitive, and handling NULLs.
    \param s1 First string
    \param s2 Second string
    \param case_insensitive If 0, case sensitive strcmp is done, otherwise comparision is case insensitive.
*/
int ape_util_strcmp(const char * s1, const char * s2, char case_insensitive);

/** \brief Perform string comparison, strcmp style, optionally case insensitive, handling NULLs, between
    two arrays of strings.
    \param s1 First string array.
    \param s2 Second string array.
    \param case_insensitive If 0, case sensitive strcmp is done, otherwise comparision is case insensitive.
*/
int ape_util_cmp_string_array(const char ** s1, const char ** s2, char case_insensitive);

/** \brief Look for a string in a range of pointers to string, return an iterator pointing to the matching string.
    \param begin The beginning iterator (pointer to string) to search.
    \param end The end (one-past-last) iterator (pointer to string) to search.
    \param input The string for which to search.
    \param found Pointer to output pointer to string.
    \param case_insensitive If 0, case sensitive strcmp is done, otherwise comparision is case insensitive.
*/
int ape_util_find_string(char ** begin, char ** end, const char * input, char *** found, char case_insensitive);

/** \brief Check for file access with the given mode, using a function which may be overridden by client to
    determine if the file is available. By default, standard C fopen is used to test access.
    Returned value is 1 if the access is successful, or 0 otherwise.
    \param file_name The name of the file to check.
    \param access The type of access, may be r, w etc. 
*/
int ape_util_check_file_access(const char * file_name, const char * access);

/** \brief Return pointer to function ape currently uses to check file access.
    \param func Pointer to the function being used.
*/
int ape_util_get_file_check_func(int (**func)(const char *, const char *));

/** \brief Configure ape to use the given function to check file access.
    \param func The function to use.
*/
int ape_util_set_file_check_func(int (*func)(const char *, const char *));

/** \brief Operating system independent sleep function.
    \param sleep_time The time to sleep, in seconds.
*/
void ape_util_sleep(int sleep_time);

/** \brief Get a line of text from the user.
    \param prompt Prompt to display.
    \param text Input text.
*/
int ape_util_get_text(const char * prompt, char ** text);

/** \brief Redirect prompts to the stream.
    \param stream The stream.
*/
int ape_util_set_prompt_stream(FILE * stream);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_util.h,v $
 * Revision 1.22  2011/02/01 18:01:46  jpeachey
 * Add support and tests for functions to convert string values into int
 * and short types, and to convert a parameter value into short type.
 *
 * Revision 1.21  2010/11/12 20:53:03  irby
 * Moved some internal static functions from ape_par to ape_util, and allowed
 * new custom_get_text routine to be used instead of ape's standard getter in
 * order to handle xpi call-back.
 *
 * Revision 1.20  2009/06/11 19:09:55  peachey
 * Add ape_util_sleep function, needed for unit test of ape_io with
 * new behavior for hidden parameters (system time stamp overrules local
 * hidden parameters.)
 *
 * Revision 1.19  2007/11/12 19:45:19  peachey
 * Add ape_util_interpret_env function for setting internal state of Ape
 * based on environment variables. Remove spaces at end of line.
 *
 * Revision 1.18  2007/10/10 20:16:39  peachey
 * Add and test ape_util_get_file_check, for retrieving pointer to current
 * file checking function.
 *
 * Revision 1.17  2006/12/21 21:44:29  peachey
 * Add ape_util_find_string, to find the first match to a string
 * in an array of strings.
 *
 * Revision 1.16  2006/11/30 20:38:37  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.15  2006/11/24 19:48:00  peachey
 * Add ape_util_s2f.
 *
 * Revision 1.14  2006/06/15 13:37:08  peachey
 * Add infrastructure to support expansion of environment variables.
 *
 * Revision 1.13  2006/06/07 05:54:04  peachey
 * Reverse the output of ape_util_file_check_access: 1 means access is OK,
 * 0 means it isn't.
 *
 * Revision 1.12  2006/06/05 01:30:47  peachey
 * Include file access checking facilities, including ability for client
 * to supply a custom access checking function.
 *
 * Revision 1.11  2006/05/18 03:17:06  peachey
 * Add ape_util_append_file_name for constructing full file names
 * from directory name + file name. Also improve handling of indef, nan etc.
 * in numeric conversions.
 *
 * Revision 1.10  2006/04/28 02:50:58  peachey
 * Fix sloppy implementation of ape_util_strcmp and add similar
 * comparison for arrays of strings, ape_util_cmp_string_array.
 *
 * Revision 1.9  2006/04/26 01:30:01  peachey
 * Add ape_util_free_string, a utility for freeing arrays of char *.
 *
 * Revision 1.8  2006/04/22 01:35:06  peachey
 * Add public function ape_util_copy_range.
 *
 * Revision 1.7  2006/04/21 14:26:12  peachey
 * Add ape_util_strcmp, for comparing strings with optional case insensitivity.
 *
 * Revision 1.6  2006/04/19 15:40:32  peachey
 * Add ape_util_atexit for registering multiple functions with
 * atexit while still using only one real atexit function.
 *
 * Revision 1.5  2006/04/13 18:43:37  peachey
 * Add and test ape_util_cat_string, for concatenating strings.
 *
 * Revision 1.4  2006/04/12 17:59:19  peachey
 * Add ape_util_parse_pfiles, for splitting into local and system parts of PFILES.
 *
 * Revision 1.3  2006/04/12 14:17:23  peachey
 * Add and test ape_util_getenv.
 *
 * Revision 1.2  2006/04/11 01:03:20  peachey
 * Add ape_util_s2d and its unit test. Improve handling of overflows to silence osx warning.
 *
 * Revision 1.1  2006/04/10 21:10:55  peachey
 * New module for general utilities, such as string conversions.
 *
*/
