#ifndef SIXT_ERROR_H
#define SIXT_ERROR_H (1)


/////////////////////////////////////////////////////////////////
// Macro Definitions.
/////////////////////////////////////////////////////////////////


#define SIXT_ERROR(msg) (sixt_error(__func__, msg))

#define CHECK_STATUS_BREAK(status) \
  if (EXIT_SUCCESS!=status) break;

#define CHECK_STATUS_RET(status, retval)		\
  if (EXIT_SUCCESS!=status) return(retval);

#define CHECK_NULL(a,status,msg) \
  if (NULL==a) { \
    SIXT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return(NULL);\
  }

#define CHECK_NULL_VOID(a,status,msg) \
  if (NULL==a) { \
    SIXT_ERROR(msg); \
    status=EXIT_FAILURE; \
    return;\
  }

#define CHECK_NULL_BREAK(a,status,msg) \
  if (NULL==a) { \
    SIXT_ERROR(msg); \
    status=EXIT_FAILURE; \
    break;\
  }


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Print the given error message for an error occured in the
    specified function. The function name is also part of the
    output. */
void sixt_error(const char* const func, const char* const msg);


#endif /* SIXT_ERROR_H */
