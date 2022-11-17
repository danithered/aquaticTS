#ifndef _RANDOMSZAM_
#define _RANDOMSZAM_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
extern "C" {
#include <time.h>
#include <stdio.h>
}
	
extern gsl_rng * r;

/*A main-be a kovetkezo sorokat be kell tenni az elejere es a vegere:*/
/*
 * 
 
//randomszam generator inicializalasa main elejen
	time_t timer;
	randomszam_inic(time(&timer), r);

//randomszam generator lezarasa main vegen
	gsl_rng_free(r);
*/

void randomszam_inic(int seed, gsl_rng * rng);

int randomszam_mentes(const char * filename, gsl_rng * rng);

int randomszam_olvasas(const char * filename, gsl_rng * rng);


#endif

