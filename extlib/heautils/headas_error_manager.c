/******************************************************************************
 *   File name: headas_error_manager.c                                        *
 *                                                                            *
 * Description: Implementation of a HEAdas error manager, which maintains     *
 *     an internal integer error status variable, a list of error maps, and a *
 *     stack of error messages.                                               *
 *                                                                            *
 *      Author: James Peachey, LAC, for HEASARC/GSFC/NASA                     *
 *                                                                            *
 *  Change log: see CVS Change log at the end of the file.                    *
 ******************************************************************************/

/******************************************************************************
 * Header files.                                                              *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "headas_error_internal.h"
/******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
 * Global function definitions.                                               *
 ******************************************************************************/
int HDerror_manager_init(hd_error_manager* manager, hd_error_msg* msgArray,
    unsigned int msgArraySize, char* buf, hd_error_map* map,
    unsigned int bufSize, int status) {
  if(HD_OK != status) return status;

  do {
    /* Check arguments before proceeding. */
    if(NULL == manager) status = HD_ERR_NULL_POINTER;

    if(0u != msgArraySize && NULL == msgArray)
      status = HD_ERR_MNGR_ERROR;

    if(0u != bufSize && NULL == buf) status = HD_ERR_MNGR_ERROR;

    if(HD_OK != status) continue;

    manager->MsgArrayBegin = msgArray;
    manager->MsgLast = msgArray;
    manager->MsgArrayEnd = msgArray + msgArraySize;

    manager->TextArrayBegin = buf;
    manager->TextLast = buf;
    manager->TextArrayEnd = buf + bufSize;

    /* Optional map of error messages. */
    manager->Map = map;

    manager->ErrorNumber = HD_OK;

    manager->SilenceDump = 0;

  } while(0);

  return status;
}

int HDerror_manager_throw(hd_error_manager* manager, const char* text,
    const char* fileName, int line, int errNum) {
  do {
    /* If this manager is invalid, throw this error on the default manager. */
    if(NULL == manager) {
      manager = HDerror_manager_get_default();

      /* The default manager should never ever ever be NULL, but just in
         case... */
      if(NULL == manager) {
        errNum = HD_ERR_MNGR_ERROR;
        continue;
      }
    }

    if(HD_OK == errNum) {
      /* No error is currently being thrown, but if an error was
         thrown previously and was not handled, throw the (previous) error
         anyway. */
      if(HD_OK != manager->ErrorNumber) errNum = manager->ErrorNumber;
      else continue;
    }

    /* Reset the error manager only if there was not a previous error. */
    if(HD_OK == manager->ErrorNumber) {

      /* Reset to the beginning of the error message array and buffer. */
      manager->MsgLast = manager->MsgArrayBegin;
      manager->TextLast = manager->TextArrayBegin;

      /* Set the manager's status. */
      manager->ErrorNumber = errNum;
    }

    /* Put the text etc. into a hint. */
    HDerror_manager_hint(manager, text, fileName, line, errNum);

  } while(0);

  return errNum;
}

int HDerror_manager_hint(hd_error_manager* manager, const char* text,
    const char* fileName, int line, int errNum) {
  do {
    /* If this manager is invalid, throw this error on the default manager. */
    if(NULL == manager) {
      manager = HDerror_manager_get_default();

      /* The default manager should never ever ever be NULL, but just in
         case... */
      if(NULL == manager) {
        errNum = HD_ERR_MNGR_ERROR;
        continue;
      }
    }

    if(HD_OK == errNum) {
      /* No error is currently being thrown, but if an error was
         thrown previously and was not handled, throw the (previous) error
         anyway. */
      if(HD_OK != manager->ErrorNumber) errNum = manager->ErrorNumber;
      else continue;
    }

    /* Make sure there is room for the message. */
    if(NULL != manager->MsgLast && manager->MsgLast < manager->MsgArrayEnd) {
      char* newText;

      /* Don't bother copying the text if there is no need or no space. */
      if(NULL == text || NULL == manager->TextLast ||
        manager->TextLast >= manager->TextArrayEnd)
        newText = NULL;
      else {
        /* Start from the current spot in the buffer. */
        char* tmp_p = manager->TextLast;

        /* Copy the text to the buffer while there is room in that buffer. */
        while('\0' != *text && tmp_p < manager->TextArrayEnd - 1)
          *tmp_p++ = *text++;

        /* Make certain of the terminating NULL. */
        *tmp_p = '\0';

        /* Text now begins at the current spot in buffer. */
        newText = manager->TextLast;

        /* Update the current spot to point to the position after the NULL. */
        manager->TextLast = ++tmp_p;
      }

      manager->MsgLast->Text = newText;
      manager->MsgLast->FileName = fileName;
      manager->MsgLast->LineNumber = line;
      manager->MsgLast->ErrorNumber = errNum;
      ++manager->MsgLast;
    }
  } while(0);

  return errNum;
}

int HDerror_manager_reset(hd_error_manager* manager) {
  if(NULL == manager) return HD_ERR_NULL_POINTER;

  /* Reset to the beginning of the error message array and buffer. */
  manager->MsgLast = manager->MsgArrayBegin;
  manager->TextLast = manager->TextArrayBegin;

  /* Reset the error. */
  return manager->ErrorNumber = HD_OK;
}

void HDerror_manager_pop_msg(hd_error_manager* manager, char* stack_empty,
    const char** text, const char** fileName, int* line) {
  do {
    if(NULL == stack_empty) continue;
    *stack_empty = 1;

    if(NULL == manager) continue;

    if(manager->MsgLast > manager->MsgArrayBegin) {
      /* The current message pointer points to where new messages will
         be added, so the top message is the one before the current
         message pointer. */
      --manager->MsgLast;

      /* Return the requested information. */
      if(NULL != text) *text = manager->MsgLast->Text;
      if(NULL != fileName) *fileName = manager->MsgLast->FileName;
      if(NULL != line) *line = manager->MsgLast->LineNumber;

      /* The current message has now been popped, so mark the current
         buffer position as available. */
      if(NULL != manager->MsgLast->Text)
        manager->TextLast = manager->MsgLast->Text;

      /* Indicate status of stack to the caller. */
      *stack_empty = 0;

    } else {
      /* Return null or invalid values for the requested information. */
      if(NULL != text) *text = NULL;
      if(NULL != fileName) *fileName = NULL;
      if(NULL != line) *line = HD_NO_LINE_NUMBER;
    }
  } while(0);
}

int HDerror_manager_get_err(hd_error_manager* manager) {
  if(NULL == manager) return HD_ERR_NULL_POINTER;
  return manager->ErrorNumber;
}

void HDerror_manager_dump(hd_error_manager* manager, FILE* strm) {
  const char* text = NULL;
  const char* fileName = NULL;
  int line = HD_NO_LINE_NUMBER;
  char stack_empty = 0;

  do {
    if(NULL == manager) continue;

    /* Do not display any error messages if in silent mode. */
    if(manager->SilenceDump) continue;

    while(!stack_empty) {
      /* Get the next message off the stack. */
      HDerror_manager_pop_msg(manager, &stack_empty, &text, &fileName, &line);

      /* Print out the text if it's meaningful. */
      if(NULL != text && '\0' != *text) fprintf(strm, "%s ", text);

      /* Add file name and line number information if it's available. */
      if(NULL != fileName && '\0' != fileName) {
        if(HD_NO_LINE_NUMBER == line) fprintf(strm, "(in %s)", fileName);
        else fprintf(strm, "(at %s: %d)", fileName, line);
      }

      /* Print a newline only if something useful was already printed. */
      if(NULL != text || NULL != fileName) fprintf(strm, "\n");

    }
  } while(0);
}

void HDerror_manager_find_map_entry(hd_error_manager* manager, int errNum,
    char* found, const char** pkgid, const char** code, const char** text) {
  hd_error_map* errMap;
  do {
    if(NULL == found) continue;

    *found = 0;

    if(NULL == manager) continue;

    errMap = manager->Map;
    while(NULL != errMap) {
      HDerror_map_get_entry(errMap, errNum, found, pkgid, code, text);
      if(*found) break;
      errMap = HDerror_map_next_map(errMap);
    }
  } while(0);
}

int HDerror_manager_get_map(hd_error_manager* manager, hd_error_map** errMap,
    int status) {
  if(HD_OK != status) return status;
  do {
    if(NULL == errMap) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    *errMap = NULL;

    if(NULL == manager) {
      status = HD_ERR_NULL_POINTER;
      continue;
    }

    *errMap = manager->Map;

  } while(0);
  return status;
}

void HDerror_manager_dump_silence(hd_error_manager* manager, int silent) {
  if(NULL == manager) return;
  manager->SilenceDump = silent;
}

int HDerror_manager_dump_is_silent(hd_error_manager* manager) {
  if(NULL == manager) return 1;
  return manager->SilenceDump;
}
/******************************************************************************/

#ifdef __cplusplus
}
#endif

/******************************************************************************
 * $Log: headas_error_manager.c,v $
 * Revision 1.5  2004/01/07 14:45:52  peachey
 * Add functions HDerror_dump_silence and HDerror_dump_is_silent which
 * set and get a silent flag. When this flag is set, HDerror_dump is a no-op.
 *
 * Revision 1.4  2002/10/09 16:42:22  peachey
 * Make sure that heaerr, heaout, heaprom are used consistently.
 *
 * Revision 1.3  2002/10/04 21:47:16  peachey
 * Replace HDerror_manager_get* family with HDerror_manager_find_map_entry.
 * This allows all the interesting information to be obtained from the
 * error map in one call.
 *
 * Revision 1.2  2002/10/04 19:10:20  peachey
 * Remove a comment.
 *
 * Revision 1.1  2002/09/16 13:45:41  peachey
 * Add headas error handling facilities.
 *
 ******************************************************************************/
