/** \file pil_error.h
    \brief Declaration of Ape's PIL emulating error codes and error reporting utilities.

    Ape is backward compatible with a subset of the PIL library developed by
    ISDC. Replacements are provided for the most commonly used API functions.
    Less common functions and functions which rely on details of pil's internal
    structures are not provided. Ape's implementation of these functions does
    not always behave identically to PIL's implementation, but the behavior should be
    close enough to provide compatibility.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#ifndef ape_pil_error_h
#define ape_pil_error_h

#ifdef __cplusplus
extern "C" {
#endif

#define PIL_OK                          (0)

#define PIL_ERR_BASE                    (-3000)
#define PIL_ERR_MAX_IDX                 (PIL_ERR_BASE)

#define PIL_NUL_PTR                     (PIL_ERR_BASE - 0)
#define PIL_BAD_ARG                     (PIL_ERR_BASE - 1)
#define PIL_NO_MEM                      (PIL_ERR_BASE - 2)
#define PIL_NO_FILE                     (PIL_ERR_BASE - 3)
#define PIL_ERR_FREAD                   (PIL_ERR_BASE - 4)
#define PIL_ERR_FWRITE                  (PIL_ERR_BASE - 5)
#define PIL_EOS                         (PIL_ERR_BASE - 6)
#define PIL_BAD_NAME                    (PIL_ERR_BASE - 7)
#define PIL_BAD_TYPE                    (PIL_ERR_BASE - 8)
#define PIL_BAD_MODE                    (PIL_ERR_BASE - 9)
#define PIL_BAD_LINE                    (PIL_ERR_BASE - 10)
#define PIL_NOT_IMPLEMENTED             (PIL_ERR_BASE - 11)
#define PIL_FILE_NOT_EXIST              (PIL_ERR_BASE - 12)
#define PIL_FILE_EXIST                  (PIL_ERR_BASE - 13)
#define PIL_FILE_NO_RD                  (PIL_ERR_BASE - 14)
#define PIL_FILE_NO_WR                  (PIL_ERR_BASE - 15)
#define PIL_LINE_BLANK                  (PIL_ERR_BASE - 16)
#define PIL_LINE_COMMENT                (PIL_ERR_BASE - 17)
#define PIL_LINE_ERROR                  (PIL_ERR_BASE - 18)
#define PIL_NOT_FOUND                   (PIL_ERR_BASE - 19)
#define PIL_PFILES_TOO_LONG             (PIL_ERR_BASE - 20)
#define PIL_PFILES_FORMAT               (PIL_ERR_BASE - 21)
#define PIL_LOCK_FAILED                 (PIL_ERR_BASE - 22)
#define PIL_BOGUS_CMDLINE               (PIL_ERR_BASE - 23)
#define PIL_NO_LOGGER                   (PIL_ERR_BASE - 24)
#define PIL_LINE_TOO_MANY               (PIL_ERR_BASE - 25)
#define PIL_LINE_TOO_FEW                (PIL_ERR_BASE - 26)
#define PIL_LINE_UNMATCHED_QUOTE        (PIL_ERR_BASE - 27)
#define PIL_LINE_NO_LF                  (PIL_ERR_BASE - 28)
#define PIL_LINE_EXTRA_SPACES           (PIL_ERR_BASE - 29)
#define PIL_BAD_VAL_BOOL                (PIL_ERR_BASE - 30)
#define PIL_BAD_VAL_INT                 (PIL_ERR_BASE - 31)
#define PIL_BAD_VAL_REAL                (PIL_ERR_BASE - 32)
#define PIL_BAD_VAL_INT_VAR_VECTOR      (PIL_ERR_BASE - 33)
#define PIL_BAD_VAL_INT_VECTOR          (PIL_ERR_BASE - 34)
#define PIL_BAD_VAL_REAL_VAR_VECTOR     (PIL_ERR_BASE - 35)
#define PIL_BAD_VAL_REAL_VECTOR         (PIL_ERR_BASE - 36)
#define PIL_OFF_RANGE                   (PIL_ERR_BASE - 37)
#define PIL_BAD_ENUM_VALUE              (PIL_ERR_BASE - 38)
#define PIL_BAD_FILE_ACCESS             (PIL_ERR_BASE - 39)
#define PIL_BAD_VALUE                   (PIL_ERR_BASE - 40)
#define PIL_VALUE_UNDEFINED             (PIL_ERR_BASE - 41)
#define PIL_UNSPECIFIED_ERROR           (PIL_ERR_BASE - 42)

#define PIL_ERR_MIN_IDX                 (PIL_ERR_MAX_IDX - 42)

const char * PIL_err_handler(int r);

#ifdef __cplusplus
}
#endif

#endif

/*
 * $Log: pil_error.h,v $
 * Revision 1.3  2006/06/06 17:46:59  peachey
 * Add some error codes to cover essential ape errors which have no
 * counterpart in pil error code space.
 *
 */
