/**CFile***********************************************************************

  FileName    [plaGasp.c]

  PackageName [pla]

  Synopsis    [Last gasp]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

/*
static char rcsid[] = \"$Id: $\";
USE(rcsid);
*/

#include "ansi.h"
#include "util.h"
#include "plaInt.h"

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

  Synopsis           [Main function]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaGasp(
  cover_t **F,
  cover_t *D,
  cover_t *R,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t *tmp, *red;
unsigned int *c, *c2;
int i,j,k;

  red = pla_CoverAlloc(num_inputs);

  for(i = 0;i < (*F)->ncubes;i++) {
    tmp = pla_CoverAlloc(num_inputs);
    for(j = 0;j < (*F)->ncubes;j++) {
      if(i == j) continue;
      c2 = plaCubeCofactor((*F)->cubes[j], (*F)->cubes[i], num_inputs,
        num_outputs, cube_size);
      if(c2 != NIL(unsigned int)) pla_CubeAdd(tmp, c2);
    }
    for(j = 0;j < D->ncubes;j++) {
      c2 = plaCubeCofactor(D->cubes[j], (*F)->cubes[i], num_inputs,
        num_outputs, cube_size);        
      if(c2 != NIL(unsigned int)) pla_CubeAdd(tmp, c2);
    }
    if(tmp->ncubes > 0) {
      c = plaSCCC(tmp, num_inputs, num_outputs, cube_size);      
      c2 = (unsigned int *) ALLOC(unsigned int, cube_size);
      AND_CUBE(c2, c, (*F)->cubes[i], cube_size, k);
      FREE(c); j = 1;
      j = 0;
      for(k = 0;k < num_inputs;k++) {
        if(VALUE_VAR(c2, k) != VALUE_VAR((*F)->cubes[i], k)) {
          j = 1;
          break;
        }      
      }
      if(j == 0) {
        for(k = 0;k < num_outputs;k++) {
          if((VALUE_OUT(c2, k) && !VALUE_OUT((*F)->cubes[i], k)) ||
             (!VALUE_OUT(c2, k) && VALUE_OUT((*F)->cubes[i], k))) {
            j = 1;
            break;
          }
        }
      }
      if(j == 1) {
        pla_CubeAdd(red, c2);      
      }
    }        
    pla_CoverFree(tmp);
  }
  plaExpand(&red, R, num_inputs, num_outputs, cube_size, MINCUBE, 1);
  
  if(red->ncubes > 0) {
    for(i = 0;i < (*F)->ncubes;i++) {
      c = (unsigned int *) ALLOC(unsigned int, cube_size);
      COPY_CUBE(c, (*F)->cubes[i], cube_size, j);
      pla_CubeAdd(red, c);
    }  
    plaIrredundant(&red, D, num_inputs, num_outputs, cube_size);
    if(red->ncubes < (*F)->ncubes) {
      pla_CoverFree(*F);
      *F = red;
    }
    else {
      pla_CoverFree(red);
    }
  }
  else {
    pla_CoverFree(red);
  }
}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/



