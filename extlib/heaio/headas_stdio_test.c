#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

/* This test program should use only the hdio library. */
#include "headas_stdio.h"

int main(void) {
  /* Right away test tmpfile and tmpnam. */
  const char *fname = HD_tmpnam(0);
  const char *tmp;
  hd_FILE *fptr = HD_tmpfile();
  struct stat myStat;
  char buf[BUFSIZ];
  char vbuf[64];
  fpos_t fpos;
  int status;

  status = HDIO_init();
  if (status) return status;

  /* Make sure tmpfile worked. */
  if(!fptr)
    HD_fprintf(hd_err, "Problem opening a tmpfile\n");
  else if(10 != HD_fprintf(fptr, "ten %s\n", "chars"))
    HD_fprintf(hd_err, "Problem writing to a tmpfile\n");
  else if(EOF == HD_fclose(fptr))
    HD_fprintf(hd_err, "Problem closing a tmpfile\n");

  /* Test tmpnam. */
  fptr = HD_fopen(fname, "w");
  if(!fptr)
    HD_fprintf(hd_err, "Problem opening a tmpnam file\n");

  else if(10 != HD_fprintf(fptr, "ten %s\n", "chars"))
    HD_fprintf(hd_err, "Problem writing to a tmpnam file\n");

  else if(EOF == HD_fclose(fptr))
    HD_fprintf(hd_err, "Problem closing a tmpnam file\n");

  /* Test creating, writing to and closing a normal file. */
  fname = "test-result1";
  fptr = HD_fopen(fname, "w");
  if(!fptr)
    HD_fprintf(hd_err, "Problem opening a normal file\n");

  else if(10 != HD_fprintf(fptr, "ten %s\n", "chars"))
    HD_fprintf(hd_err, "Problem writing to a normal file\n");

  else if(EOF == HD_fclose(fptr))
    HD_fprintf(hd_err, "Problem closing a normal file\n");

  /* Make sure it really was created. */
  if(stat(fname, &myStat))
    HD_fprintf(hd_err, "Can't find the created file\n");

  /* Test renaming a normal file. */
  tmp = "test-result";
  if(HD_rename(fname, tmp))
    HD_fprintf(hd_err, "Problem renaming %s\n", tmp);

  fname = tmp;

  /* Make sure it really was moved. */
  if(stat(tmp, &myStat))
    HD_fprintf(hd_err, "Can't find %s\n", tmp);

  /* Now remove it. */
  if(HD_remove(tmp))
    HD_fprintf(hd_err, "Problem removing %s\n", tmp);

  /* Make sure it really was moved. */
  else if(!stat(tmp, &myStat))
    HD_fprintf(hd_err, "Can still find %s after removing it\n", tmp);

  /* Test setbuf, and use the (large) buffer to test fflush. */
  HD_setbuf(hd_out, buf);

  tmp = "Hit C/R.";

  /* Test fflush. */
  if((strlen(tmp) + 1) != HD_fprintf(hd_prompt, "%s\n", tmp))
    HD_fprintf(hd_err, "Problem typing %s.\n", tmp);
  else if(HD_fflush(hd_out))
    HD_fprintf(hd_err, "Problem flushing %s.\n", tmp);

  HD_printf("You should only see this after you hit that C/R.\n");

  /* Test getchar, and use it to confirm that fflush is working.*/
  if(EOF == HD_getchar())
    HD_fprintf(hd_err, "Problem getting a character.\n");

  if(HD_fflush(hd_out))
    HD_fprintf(hd_err, "Problem flushing %s.\n", tmp);

  /* Test fgetc. */
  if((strlen(tmp) + 1) != HD_fprintf(hd_prompt, "%s\n", tmp))
    HD_fprintf(hd_err, "Problem typing %s.\n", tmp);
  else if(HD_fflush(hd_out))
    HD_fprintf(hd_err, "Problem flushing %s.\n", tmp);
  else if(EOF == HD_fgetc(hd_in))
    HD_fprintf(hd_err, "Problem getting a character.\n");

  /* Test getc. */
  if((strlen(tmp) + 1) != HD_fprintf(hd_prompt, "%s\n", tmp))
    HD_fprintf(hd_err, "Problem typing %s.\n", tmp);
  else if(HD_fflush(hd_out))
    HD_fprintf(hd_err, "Problem flushing %s.\n", tmp);
  else if(EOF == HD_getc(hd_in))
    HD_fprintf(hd_err, "Problem getting a character.\n");

  tmp = "Type something and hit C/R.";

  /* Test fgets. */
  if((strlen(tmp) + 1) != HD_fprintf(hd_prompt, "%s\n", tmp))
    HD_fprintf(hd_err, "Problem typing %s.\n", tmp);
  else if(HD_fflush(hd_out))
    HD_fprintf(hd_err, "Problem flushing %s.\n", tmp);
  else if(!HD_fgets(vbuf, 64, hd_in))
    HD_fprintf(hd_err, "Problem getting a string.\n");

  if(strchr(vbuf, '\n')) *strchr(vbuf, '\n') = '\0';
  HD_printf("\nI think you typed \"%s\".\n", vbuf);
  HD_fflush(hd_out);

  /* Test getenv/setenv. */
  tmp = HD_getenv("HEADAS_STDIO_TEST");
  if (!tmp && getenv("HEADAS_STDIO_TEST"))
    HD_fprintf(hd_err, "HD_getenv is defective; it claims environment variable "
        "HEADAS_STDIO_TEST is unset.");
  else if (tmp && strcmp(tmp, "headas stdio test string"))
    HD_fprintf(hd_err, "Environment variable HEADAS_STDIO_TEST is \"%s\", not "
        "\"headas stdio test string\".\n", tmp);

  tmp = HD_setenv("HEADAS_STDIO_TEST", "headas stdio test string 2");
  if (!tmp || strcmp(tmp, "headas stdio test string 2"))
    HD_fprintf(hd_err, "HD_setenv is defective; it claims environment variable "
        "HEADAS_STDIO_TEST is %s, not \"headas stdio test string 2\".\n", tmp);

  tmp = HD_getenv("HEADAS_STDIO_TEST");
  if (!tmp || strcmp(tmp, "headas stdio test string 2"))
    HD_fprintf(hd_err, "HD_getenv is defective; it claims environment variable "
        "HEADAS_STDIO_TEST is %s, not \"headas stdio test string 2\".\n", tmp);

  HD_printf(
      "\nWhen this executable finishes, you should have a file in the current\n"
      "directory called %s. It should contain the following text:\n"
      "HD_printf   wrote this in vbuf after redirection.\n"
      "fputs wrote this, then it was modified somewhat.\n", fname);

  /* Test freopen */
  if(!(hd_out = HD_freopen(fname, "w", hd_out))) {
    HD_fprintf(hd_err, "Problem reopening hd_out as %s.\n", fname);
    return 0;
  }

  memset(vbuf, '\0', 64);

  /* Test setvbuf */
  HD_setvbuf(hd_out, vbuf, _IOFBF, 64);

  HD_printf("HD_printf   wrote this in vbuf after redirection.\n");

  /* Make sure sprintf now goes to the same buffer as sprintf: */
  HD_sprintf(vbuf, "HD_sprintf\n");

  if(strcmp(vbuf, "HD_sprintf\n\0wrote this in vbuf after redirection.\n"))
    HD_fprintf(hd_err, "vbuf does not contain what it should.\n");

  HD_fflush(hd_out);

  HD_fputs("fputs wrote this, then it was ", hd_out);

  HD_fgetpos(hd_out, &fpos);

  HD_puts("not changed.\n");

  HD_fsetpos(hd_out, &fpos);

  HD_fputs("modified somewha", hd_out);

  HD_putc('t', hd_out);

  HD_putchar('.');

  HD_fputc('\n', hd_out);

  /* finally, fclose a standard stream. */
  HD_fclose(hd_out);
  return 0;
}
