#include <stdlib.h>
#include "headas_rand.h"
#include "mt.h"

unsigned long int HDmt_rand (HDmt_state *state) {
	return mt_get(state);
}
double HDmt_drand (HDmt_state *state) {
	return mt_get_double(state);
}
HDmt_state *HDmt_srand (unsigned long int s) {
	mt_state_t* state = calloc(1, sizeof(mt_state_t));
	mt_set(state, s);
	return (state);
}
void HDmt_destroy_state (HDmt_state *state) {
	free (state);
}

HDmt_state *mtstate;

unsigned long int HDmtRand () {
	return mt_get(mtstate);
}
double HDmtDrand () {
	return mt_get_double(mtstate);
}
void HDmtInit (unsigned long int s) {
	mtstate = (HDmt_state *) calloc(1, sizeof(HDmt_state));
	mt_set(mtstate, s);
}
void HDmtFree () {
	free (mtstate);
}
