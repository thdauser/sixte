#include <string.h>
#include "cfortran.h"
#include "headas_utils.h"

#define HD_NAME_SIZE 128
#define HD_VERS_SIZE 128

static char hdtoolname[HD_NAME_SIZE] = " ";
static char hdtoolversion[HD_VERS_SIZE] = " ";

void set_toolname(const char *name){

    strncpy(hdtoolname, name, HD_NAME_SIZE);
    hdtoolname[HD_NAME_SIZE - 1] = 0;
}
FCALLSCSUB1(set_toolname, HDNAMESET, hdnameset, STRING)

void get_toolname(char *name){

    strcpy(name, hdtoolname);
}
FCALLSCSUB1(get_toolname, HDNAMEGET, hdnameget, PSTRING)

void set_toolversion(const char *vers){

    strncpy(hdtoolversion, vers, HD_VERS_SIZE);
    hdtoolversion[HD_VERS_SIZE - 1] = 0;
}
FCALLSCSUB1(set_toolversion, HDVERSET, hdverset, STRING)

void get_toolversion(char *vers){

    strcpy(vers, hdtoolversion);
}
FCALLSCSUB1(get_toolversion, HDVERGET, hdverget, PSTRING)

void get_toolnamev(char *str){

    strcpy(str, hdtoolname);
    strcat(str, "_");
    strcat(str, hdtoolversion);
}
FCALLSCSUB1(get_toolnamev, HDNAMEVGET, hdnamevget, PSTRING)
