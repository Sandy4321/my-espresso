/**CFile***********************************************************************

  FileName    [plaEssen.c]

  PackageName [pla]

  Synopsis    [Essencial primes]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

#include "ansi.h"
#include "util.h"
#include "plaInt.h"

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

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/


/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

/**AutomaticEnd***************************************************************/


/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Definition of internal functions                                          */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaEssencial(
  cover_t **F,
  cover_t **D,
  int *moved,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i,j;
cover_t *ret, *cv1, *cv2, *ctest;
unsigned int *c;
short *dl;

  dl = (short *) ALLOC(short, (*F)->ncubes);
  for(i = 0;i < (*F)->ncubes;i++) dl[i] = 0;

  ret = pla_CoverAlloc(num_inputs);
  for(i = 0;i < (*F)->ncubes;i++) {
    cv1 = pla_CoverAlloc(num_inputs);
    for(j = 0;j < (*F)->ncubes;j++) {
      if(j == i) continue;
      c = plaCubeConsensus((*F)->cubes[j], (*F)->cubes[i], num_inputs,
        num_outputs, cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(cv1, c);
    }  
    cv2 = plaCoverCubeCofactor(cv1, (*F)->cubes[i], num_inputs, num_outputs,
      cube_size);
    pla_CoverFree(cv1);
    ctest = pla_CoverAlloc(num_inputs);
    for(j = 0;j < cv2->ncubes;j++) {
      pla_CubeAdd(ctest, cv2->cubes[j]); 
    }
    FREE(cv2->flags);
    FREE(cv2->cubes);    
    FREE(cv2->p0);
    FREE(cv2->p1);
    FREE(cv2);    

    cv1 = plaCoverConsensus((*D), (*F)->cubes[i], num_inputs, num_outputs,
      cube_size);
    cv2 = plaCoverCubeCofactor(cv1, (*F)->cubes[i], num_inputs, num_outputs,
      cube_size);
    pla_CoverFree(cv1);    
    for(j = 0;j < cv2->ncubes;j++) {
      pla_CubeAdd(ctest, cv2->cubes[j]); 
    }
    FREE(cv2->flags);
    FREE(cv2->cubes);
    FREE(cv2->p0);
    FREE(cv2->p1);
    FREE(cv2);
    if(plaTautology(ctest, num_inputs, num_outputs, cube_size) == 0) {
      c = pla_CubeNew(cube_size, 0);
      COPY_CUBE(c, (*F)->cubes[i], cube_size, j);
      pla_CubeAdd(ret, c);
      dl[i] = 1;
    }
    pla_CoverFree(ctest);
  }
  *moved = (*D)->ncubes;

  for(i = 0;i < ret->ncubes;i++) {
    pla_CubeAdd((*D), ret->cubes[i]);
  }
  FREE(ret->flags);  
  FREE(ret->cubes);
  FREE(ret->p0);
  FREE(ret->p1);
  FREE(ret);
  
  for(i = 0;i < (*F)->ncubes;i++) {
    if(dl[i] == 1) {
      ZERO_VAR((*F)->cubes[i], 0);
    }
  }
  FREE(dl);
  i = 0;
  while(i < (*F)->ncubes) {
    if(VALUE_VAR((*F)->cubes[i], 0) == 0) {
      pla_CubeDel((*F), i, cube_size);
    }
    else {
      i++;
    }  
  }
}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/



