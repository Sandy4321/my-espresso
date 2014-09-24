/**CFile***********************************************************************

  FileName    [main.c]

  PackageName [esp]

  Synopsis    [My implementation of espresso]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

#include <stdio.h>
#include "ansi.h"
#include "util.h"
#include "plaInt.h"
#include "pla.h"

/*
static char rcsid[] = \"$Id: $\";
USE(rcsid);
*/

/*---------------------------------------------------------------------------*/
/* Constant declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Structure declarations                                                    */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

extern double expand_time, irred_time, reduce_time, gasp_time;
extern double tot_time, essen_time, comp_time;

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/


/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

static void StartTimer();
static void StopTimer();
static double GetRtime();

/**AutomaticEnd***************************************************************/


/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis           [Main function]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
main(
  int argc,
  char **argv)
{
FILE *fp;
pla_t *pla;
short verif = 1;
int essen;

  if(argc < 2) {
    fprintf(stdout,"Usage: my_espresso filename\n");
    exit(1);
  }  
  fp = fopen(argv[1], "r");
  if(fp == NULL) {
    fprintf(stdout,"Could not open file %s for reading.\n", argv[1]);
    exit(1);
  }
  pla = plaRead(fp);
  (void) fclose(fp);

  pla_Espresso(pla, 1);
  pla_Print(stdout, pla, 0);

  if(verif == 1) {
    fp = fopen("verif.pla", "w");
    pla_Print(fp, pla, 0);
    (void) fclose(fp);
  }    
  plaCost(pla->F, pla->num_inputs, pla->num_outputs);
  fprintf(stdout, "Cost = %d  in = %d  out = %d\n", pla->F->ncubes,
    pla->F->cost_in, pla->F->cost_out);

  tot_time = expand_time + irred_time + reduce_time + gasp_time + comp_time
    + essen_time;
  fprintf(stdout,"TOTAL TIME: %f\n\n", tot_time);
  fprintf(stdout,"COMPLEMENT TIME: %f(%2.2f)\n", comp_time,
    100*comp_time/tot_time);
  fprintf(stdout,"ESSENCIAL TIME: %f(%2.2f)\n", essen_time,
    100*essen_time/tot_time);
  fprintf(stdout,"EXPAND TIME: %f(%2.2f)\n", expand_time, 
    100*expand_time/tot_time);
  fprintf(stdout,"REDUCE TIME: %f(%2.2f)\n", reduce_time,
    100*reduce_time/tot_time);
  fprintf(stdout,"IRREDUNDANT TIME: %f(%2.2f)\n", irred_time,
    100*irred_time/tot_time);
  fprintf(stdout,"LAST GASP TIME: %f(%2.2f)\n", gasp_time,
    100*gasp_time/tot_time);

  pla_Free(pla);
}

/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/
