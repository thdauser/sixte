/** \file ape_par.h
    \brief Declaration of parameter facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_par_h
#define ape_ape_par_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Number of characters needed to represent the mode of a parameter. */
#define APE_PAR_MODE_CODE_LEN 4

/* Number of characters needed to represent the type of a parameter. */
#define APE_PAR_TYPE_CODE_LEN 4

typedef enum ParFieldId {
  eName = 0,
  eType,
  eMode,
  eValue,
  eMin,
  eMax,
  ePrompt,
  eEndOfField
} ParFieldId;

typedef enum ParPromptModeEnum {
  eDefaultPrompt = 0,
  eNoPrompt = 1,
  eQueryHidden = 2,
  eMultiQuery = 4 /* Issue a new prompt every time ape_par_query is called, even if parameter was supplied on cmd line. */
} ParPromptModeEnum;

/* Enumerated codes governing the case of string parameters values when read. */
typedef enum ParStringCaseEnum {
  eUnknownCase = 0, /* Do not interpret the case at all. Same behavior as eDefaultCase. */
  eDefaultCase = 1, /* Return the string as-is. Do not alter the case of any characters. */
  eLowerCase = 2, /* Convert string to all lower case letters. */
  eUpperCase = 4, /* Convert string to all upper case letters. */
  eEnumCase = 8 /* For enumerated parameters, return the string that matched, not the user's input. */
} ParStringCaseEnum;

struct ApePar;
typedef struct ApePar ApePar;

/** \brief Type for custom prompting functions.
*/
typedef int (*ApeParGetTextFunc)(const char *, const char *, char **);

/** \brief Create an exact duplicate of a parameter.
    \param orig_par The parameter being cloned.
    \param clone_par The output cloned parameter.
*/
int ape_par_clone(ApePar * orig_par, ApePar ** clone_par);

int ape_par_create(ApePar ** ape_par, const char * line);

void ape_par_destroy(ApePar * ape_par);

char * ape_par_get_line(const ApePar * ape_par);

/** \brief Get the name of the given parameter as a string. The caller must free the string. If the
    parameter is unnamed, the error code eUnnamedPar will be returned.
    \param par Pointer to the parameter object whose name is sought.
    \param name The name of the parameter.
*/
int ape_par_get_name(const ApePar * ape_par, char ** name);

/** \brief Get the indicated parameter as a bool. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_bool(const ApePar * par, char * result);

/** \brief Get the indicated parameter as a double. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_double(const ApePar * par, double * result);

/** \brief Get the indicated parameter as a float. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_float(const ApePar * par, float * result);

/** \brief Get the indicated parameter as a file name. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_file_name(const ApePar * par, char ** result);

/** \brief Get the indicated parameter as an int. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_int(const ApePar * par, int * result);

/** \brief Get the indicated parameter as a long. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_long(const ApePar * par, long * result);

/** \brief Get the indicated parameter as a short integer. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_short(const ApePar * par, short * result);

/** \brief Get the indicated parameter as a string. The conversion will occur with a best effort, and any
    problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
*/
int ape_par_get_string(const ApePar * par, char ** result);

/** \brief Get the indicated parameter as a string, with added control over how the case of the string is
    handled. The conversion will occur with a best effort, and any problems will be indicated through the returned status.
    \param par Pointer to the parameter object whose value is sought.
    \param result The value of the field.
    \param case_code Code indicating how to handle the case of the string. See ParStringCaseEnum for allowed values.
*/
int ape_par_get_string_case(const ApePar * par, char ** result, char case_code);

/** \brief Get the comment associated with the indicated parameter. The caller must always free the string.
    \param par Pointer to the parameter object whose value is sought.
    \param comment The value of the comment.
*/
int ape_par_get_comment(const ApePar * par, char ** comment);

/** \brief Get the indicated field of the parameter. The client must free the character array pointed to by result.
    \param par Pointer to the parameter object whose field is sought.
    \param field_id Field identifier.
    \param result The value of the field as a string.
*/
int ape_par_get_field(const ApePar * par, ParFieldId id, char ** result);

/** \brief Get the indicated field of the parameter as an array of strings. The client must free the array pointed to
    by result using the ape_util_free_string_array. The result array will be null terminated.
    \param par Pointer to the parameter object whose field is sought.
    \param field_id Field identifier.
    \param result The value of the field as an array of strings.
*/
int ape_par_get_field_array(const ApePar * par, ParFieldId id, char *** result);

/** \brief Get the type of the parameter.
    \param par Pointer to the parameter object whose type is sought.
    \param type_string Pointer to a character array (must have at least size APE_PAR_TYPE_CODE_LEN) to hold the output type code.
*/
int ape_par_get_type(const ApePar * par, char * type_string);

/** \brief Get the innate mode of the parameter.
    \param par Pointer to the parameter object whose mode is sought.
    \param mode_string Pointer to a character array (must have at least size APE_PAR_MODE_CODE_LEN) to hold the output mode.
*/
int ape_par_get_mode(const ApePar * par, char * mode_string);

/** \brief Get the effective mode of the parameter, taking into account both the automatic mode and whether the parameter
    was already set on the command line.
    \param par Pointer to the parameter object whose mode is sought.
    \param auto_string Pointer to string array which holds the mode to use for automatic parameters. May be 0.
    \param mode_string Pointer to a character array (must have at least size APE_PAR_MODE_CODE_LEN) to hold the output mode.
*/
int ape_par_get_eff_mode(const ApePar * par, const char * auto_string, char * mode_string);

/** \brief Get the minimum value of a parameter as a string. If the minimum gives an enumerated range this still
    gets the entire enumerated range in one string, but returns a status of eParRangeEnum, to indicate this to
    the caller.
    \param par Pointer to the parameter object whose field is being set.
    \param min Pointer to the output string, which the client must free.
*/
int ape_par_get_min_string(const ApePar * par, char ** min);

/** \brief Get the set of enumerated values of a parameter as a vector of strings. If the minimum field contains
    only a single range, this will return a vector with only one string: the single range, but returns a status of
    eParRangeNoEnum, to indicate this to the caller.
    \param par Pointer to the parameter object whose field is being set.
    \param range Pointer to the null-terminated array of output strings, which the client must free using
    ape_util_free_string_array.
*/
int ape_par_get_enum_string(const ApePar * par, char *** range);

/** \brief Set the indicated field of the parameter, using the text to replace the field text completely.
    \param par Pointer to the parameter object whose field is being set.
    \param field_id Field identifier.
    \param field_text The value of the field to set, as a string.
*/
int ape_par_set_field(ApePar * par, ParFieldId id, const char * field_text);

/** \brief Set the value field of the parameter, using the text to replace the field text completely.
    Also, flag the parameter as being modified. To set the value without flagging it as modified, call
    ape_par_set_field. A parameter which has been modified will have effective mode 'h'.
    \param par Pointer to the parameter object whose value field is being set.
    \param value The value to set, as a string.
*/
int ape_par_set_value_string(ApePar * par, const char * value);

/** \brief Set the comment field of the parameter.
    \param par Pointer to the parameter object whose comment field is being set.
    \param comment The value of the comment.
*/
int ape_par_set_comment(ApePar * par, const char * comment);

/** \brief Get global default prompt style, which takes effect only if the parameter's style is eDefaultPrompt.
    Thus individual parameters may override the default value. An exception to this is eNoPrompt, which
    will suppress prompts for all parameters regardless.
    \param prompt_style Current default prompt style.
*/
int ape_par_get_default_prompt_style(int * prompt_style);

/** \brief Set global default prompt style, which takes effect only if the parameter's style is eDefaultPrompt.
    Thus individual parameters may override the default value. An exception to this is eNoPrompt, which
    will suppress prompts for all parameters regardless.
    \param prompt_style Default prompt style desired, a bitwise "or" of options present in enum ParPromptModeEnum.
*/
int ape_par_set_default_prompt_style(int prompt_style);

/** \brief Get prompt style of this parameter, taking into account the default prompt style.
    \param par Pointer to the parameter object.
    \param prompt_style Current prompt style.
*/
int ape_par_get_prompt_style(const ApePar * par, int * prompt_style);

/** \brief Set prompt style of this parameter to control when and if prompts are displayed.
    \param par Pointer to the parameter object.
    \param prompt_style Prompt style desired, a bitwise "or" of options present in enum ParPromptModeEnum.
*/
int ape_par_set_prompt_style(ApePar * par, int prompt_style);

/** \brief Prompt for, and read a string from the standard input, and use it to set the value of a parameter.
    Note that this is a low level prompter; it does not pay attention to the mode of the parameter or to the
    file's hidden "mode" parameter. It does respect global no-prompt option however.
    \param par The parameter for which to prompt, and whose value will be set.
*/
int ape_par_prompt(ApePar * par);

/** \brief Prompt for, and read a string from the standard input, and use it to set the value of a parameter.
    The prompt will only be issued if the effective mode of the parameter is "query".

    Errors in converting the user's input will be reported, but the value will be assigned to the parameter even
    if there is an error, so that reprompts may be more informative.
    \param par The parameter for which to prompt, and whose value will be set.
    \param auto_string The mode for automatic parameters.
*/
int ape_par_query(ApePar * par, const char * auto_string);

/** \brief Redirect prompts to the given stream.
    \param stream The stream to which to send prompts. If stream is a NULL pointer, this generates
    an eNullPointer error.
*/
int ape_par_redirect_prompt_stream(FILE * stream);

/** \brief Check parameter and report any errors in its fields.
    \param par_file The parameter to check.
    \param check_value Flag specifying whether the parameter value itself should be checked, allowing this
    function to be used in two contexts. Check_value == 1 is useful when all aspects of the parameter need to
    be checked, such as when checking user inputs. Check_value == 0 is useful for doing sanity check of the
    file, in which problems with the value are not show-stoppers, because the client and/or user may have
    opportunities to correct the value.
*/
int ape_par_check(const ApePar * par, char check_value);

/** \brief Low level internal function which flags the given parameter as being defined on the command line.
    Under normal circumstances this inhibits prompts.
    \param par The parameter being flagged.
*/
int ape_par_flag_cmd_line(ApePar * par);

int ape_par_register_get_text(ApeParGetTextFunc client_get_text);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_par.h,v $
 * Revision 1.30  2011/02/01 18:01:46  jpeachey
 * Add support and tests for functions to convert string values into int
 * and short types, and to convert a parameter value into short type.
 *
 * Revision 1.29  2010/11/12 20:53:03  irby
 * Moved some internal static functions from ape_par to ape_util, and allowed
 * new custom_get_text routine to be used instead of ape's standard getter in
 * order to handle xpi call-back.
 *
 * Revision 1.28  2009/07/08 19:02:02  peachey
 * Correct comment about options for case sensitivity.
 *
 * Revision 1.27  2007/10/09 16:43:44  peachey
 * Add function ape_par_redirect_prompt_stream to allow prompts to be
 * redirected. This also encapsulates all readline usage to within the ape_par module.
 *
 * Revision 1.26  2007/08/21 19:33:55  peachey
 * Add ape_par_get_name function.
 *
 * Revision 1.25  2007/02/15 21:00:17  peachey
 * Add ape_par_get/set_comment for getting and setting the
 * comment field from a parameter.
 *
 * Revision 1.24  2006/12/21 21:43:43  peachey
 * Add ape_par_get_string_case, similar to ape_par_get_string but more
 * flexible with regard to converting case of returned string values. Change
 * ape_par_get_string to call ape_par_get_string_case internally.
 *
 * Revision 1.23  2006/11/30 20:38:37  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.22  2006/11/24 19:46:43  peachey
 * Add ape_par_get_float.
 *
 * Revision 1.21  2006/06/06 13:28:06  peachey
 * Replace ape_par_get/set_prompt_mode with ape_par_get/set_prompt_style
 * and ape_par_get/set_default_prompt_style. Support no-prompt, query-for-hidden
 * and multi-query styles.
 *
 * Revision 1.20  2006/05/31 18:10:08  peachey
 * Add ape_par_get_prompt_mode/ape_par_set_prompt_mode.
 *
 * Revision 1.19  2006/05/31 03:21:43  peachey
 * Rationalize ape_par_get_* family of functions, and use them to implement
 * ape_trad_get_* and ape_trad_query_*.
 *
 * Revision 1.18  2006/05/31 01:41:37  peachey
 * Rename ape_par_set_string to ape_par_set_field.
 *
 * Revision 1.17  2006/05/31 01:36:27  peachey
 * Rename ape_par_get_string to ape_par_get_field and ape_par_get_string_array to
 * ape_par_get_field_array.
 *
 * Revision 1.16  2006/05/12 00:20:24  peachey
 * Add and test ape_par_query.
 *
 * Revision 1.15  2006/05/11 20:14:50  peachey
 * Add ape_set_value_string, for setting value and flagging it as modified.
 * Change ape_par_get_eff_mode to use this flag to make the effective mode hidden
 * after a parameter was set either on the command line or by a prompt.
 *
 * Revision 1.14  2006/05/11 19:45:47  peachey
 * Add ape_par_get_eff_mode, for getting effective mode of parameter,
 * taking into account mode parameter and command line arguments.
 *
 * Revision 1.13  2006/05/10 16:48:49  peachey
 * Generalize format tests to allow them to be run with or without
 * value checking, to let them be used both in contexts where the whole parameter
 * needs to be checked (checking user input) and contexts in which a bad value
 * is not a show-stopper, because the user or client may yet correct the problem.
 *
 * Revision 1.12  2006/05/09 19:24:25  peachey
 * Add ape_par_check function, for determining whether parameter format is valid.
 *
 * Revision 1.11  2006/05/03 01:29:12  peachey
 * Add ape_par_prompt, which uses readline on unix and homemade code on windows
 * to get a string from the user.
 *
 * Revision 1.10  2006/05/01 14:48:21  peachey
 * Clarify a comment about ape_par_set_string.
 *
 * Revision 1.9  2006/04/26 14:40:54  peachey
 * Add ape_par_get_enum_string and ape_par_get_min_string for interpreting
 * the min field of the parameter.
 *
 * Revision 1.8  2006/04/26 01:34:24  peachey
 * Add ape_par_get_string_array, for getting arrays of strings from
 * fields in parameters. Add helper functions for this: extract_field
 * (derived from get_field, which was refactored to use the new function
 * as well), and get_field_array, which parallels get_field.
 *
 * Revision 1.7  2006/04/24 21:27:26  peachey
 * Add ape_par_get_mode, reusing most of ape_par_get_type in a static function,
 * encode_msg.
 *
 * Revision 1.6  2006/04/24 17:43:00  peachey
 * Add declaration for ape_par_get_type, and macro APE_PAR_TYPE_CODE_LEN
 * to help clients declare output variable with correct size.
 *
 * Revision 1.5  2006/04/20 04:46:46  peachey
 * Add ape_par_set_string function, for setting arbitrary single fields in
 * a parameter.
 *
 * Revision 1.4  2006/04/20 03:43:58  peachey
 * Add and test ape_par_clone.
 *
 * Revision 1.3  2006/04/11 03:26:59  peachey
 * Add ape_par_get_double and its unit test.
 *
 * Revision 1.2  2006/04/11 01:55:24  peachey
 * Publish field identifier enum as part of the interface. Add and test
 * ape_par_get_string, for getting any field as a string.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
