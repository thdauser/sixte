/** \file ape_error.h
    \brief Declaration of error codes and error handling facilities.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_ape_error_h
#define ape_ape_error_h

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Error codes. */
typedef enum ape_error_enum {
  /* 0 */
  eOK,
  eNullPointer,
  eDynAllocFailed,
  eTooManyAtExitFunc,
  eAtExitError,
  /* 5 */
  eInvalidArgument,
  eFileNotFound,
  eFileReadError,
  eFileWriteError,
  eFileRenameError,
  /* 10 */
  eListError,
  eStringRemainder,
  eTypeMismatch,
  eOverflow,
  eUnderflow,
  /* 15 */
  eConversionError,
  eVarNotSet,
  eLineTooLong,
  eTooManyArguments,
  eParameterDuplicated,
  /* 20 */
  eParNotFound,
  eUnnamedPar,
  eInvalidParName, /* Invalid parameter name on command line. */
  eFieldNotFound,
  eNewParameter,
  /* 25 */
  eUnknownMode,
  eUnknownType,
  eFormatError,
  eRangeEnum,
  eRangeNoEnum,
  /* 30 */
  eInputFailure,
  eTooManyFields,
  eMinConversionError,
  eMaxConversionError,
  eValueBelowMin,
  /* 35 */
  eValueAboveMax,
  eInvalidRange,
  eInvalidChoice,
  eInvalidName, /* Invalid parameter name in parameter file. */
  eNoParLoaded,
  /* 40 */
  eUninitialized,
  eUndefinedValue,
  eNan,
  eFileNotAccessible,
  eNoModePar,
  /* 45 */
  eAmbiguousParName
} ape_error_enum;

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: ape_error.h,v $
 * Revision 1.28  2011/02/18 19:37:19  irby
 * Add new error code eAmbiguousParName and notes explaining the
 * difference between eInvalidName and eInvalidParName.
 *
 * Revision 1.27  2007/11/12 16:52:11  peachey
 * New error code eNoModePar to report absence of the mode
 * parameter when it is required.
 *
 * Revision 1.26  2007/08/21 19:35:15  peachey
 * Add new code eUnnamedPar signifying that the parameter has no name field.
 *
 * Revision 1.25  2006/11/08 22:14:33  peachey
 * Add new error code eInvalidParName.
 *
 * Revision 1.24  2006/06/23 18:59:22  peachey
 * Add error code for signalling not-a-number.
 *
 * Revision 1.23  2006/06/23 01:57:17  peachey
 * Add new code for failure to copy temporary file.
 *
 * Revision 1.22  2006/06/06 17:45:23  peachey
 * Remove unused error code eLineParseError.
 *
 * Revision 1.21  2006/06/05 01:31:15  peachey
 * Add eFileNotAccessible, for the benefit of file checking codes.
 *
 * Revision 1.20  2006/05/31 03:52:14  peachey
 * Add eUndefinedValue, for flagging converions of numeric parameters from
 * indef, undef, none, etc.
 *
 * Revision 1.19  2006/05/17 02:18:21  peachey
 * Add eUninitialized code, for when something is called without
 * calling adequate initializer, such as ape_trad_init.
 *
 * Revision 1.18  2006/05/12 03:28:56  peachey
 * Add code signifying request for parameter without parameter file loaded.
 *
 * Revision 1.17  2006/05/10 00:30:30  peachey
 * Rename eMin/MaxTypeMismatch to the more generic eMin/MaxConversionError.
 *
 * Revision 1.16  2006/05/10 00:25:12  peachey
 * Remove eValueTypeMismatch code -- too vague.
 *
 * Revision 1.15  2006/05/09 19:22:40  peachey
 * Add new error codes for reporting parameter format problems.
 *
 * Revision 1.14  2006/05/03 01:29:42  peachey
 * Add error code for input failure.
 *
 * Revision 1.13  2006/05/02 23:46:26  peachey
 * Add eTooManyArguments and eParameterDuplicated, for reporting
 * errors with ape_io_apply_command_line.
 *
 * Revision 1.12  2006/04/28 02:57:16  peachey
 * Remove unneeded code eFileOpenError and reorder remaining parameters.
 *
 * Revision 1.11  2006/04/26 14:32:51  peachey
 * Add status codes eParRangeEnum and eParRangeNoEnum for reporting
 * when client reads enum range as a simple string, or a simple string as an enum range.
 *
 * Revision 1.10  2006/04/25 15:00:10  peachey
 * Add new status code eParFieldNotFoudn, and change comments
 * which help read the error codes so that in general they will need to change
 * less when a new code is added.
 *
 * Revision 1.9  2006/04/24 21:25:54  peachey
 * Add another error code, eUnknownParMode.
 *
 * Revision 1.8  2006/04/24 17:42:00  peachey
 * Add new codes eUnknownParType and eParFormatError, plus comments to help
 * interpret status codes by eye.
 *
 * Revision 1.7  2006/04/24 13:27:58  peachey
 * Add error code for when system par file has a parameter local par file doesn't.
 *
 * Revision 1.6  2006/04/21 01:31:12  peachey
 * Add error code eParNotFound.
 *
 * Revision 1.5  2006/04/19 15:40:53  peachey
 * Add error codes for atexit related problems.
 *
 * Revision 1.4  2006/04/19 00:37:01  peachey
 * Add eLineTooLong error code.
 *
 * Revision 1.3  2006/04/12 14:19:37  peachey
 * Add error codes needed to support problems with environment variables and
 * parameter files not found.
 *
 * Revision 1.2  2006/04/10 21:11:09  peachey
 * Add error conditions for utilities.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
