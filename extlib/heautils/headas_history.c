#include "cfortran.h"
#include "headas_utils.h"

static int headas_history;

void set_history(int hpar){

    headas_history = hpar;
}
FCALLSCSUB1(set_history, HDPHIS, hdphis, INT)

void get_history(int *hvar){

    *hvar = headas_history;
}
FCALLSCSUB1(get_history, HDGHIS, hdghis, PINT)
