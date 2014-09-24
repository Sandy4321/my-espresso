/**CFile***********************************************************************

  FileName    [plaTau.c]

  PackageName [pla]

  Synopsis    [A multiple output tautology algorithm]

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

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
short
plaTautology(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int var, ret, i, j, ins, maxsum, nact, nunate;
cover_t *c1;
unsigned int *ored, *anded;

begin:

  var = -1;
  maxsum = 0;
  nact = 0;
  nunate = 0;
  
  if(cv->ncubes == 0) return 0;
  
  ored = pla_CubeNew(cube_size, 0);
  anded = pla_CubeNew(cube_size, 1);
  
  for(i = 0;i < num_inputs;i++) {
    cv->p0[i] = 0;
    cv->p1[i] = 0;
  }

  for(i = 0;i < cv->ncubes - 1;i++) {
    ins = 1;
    for(j = 0;j < num_inputs;j++) {
      ret = VALUE_VAR(cv->cubes[i], j);    
      if(ret == 1) {
        cv->p1[j]++;
        ins = 0;
      }             
      else
        if(ret == 2) {
          cv->p0[j]++;
          ins = 0;
        }          
    }
    if(ins == 1) {
      for(j = 0;j < num_outputs;j++) {
        if(!VALUE_OUT(cv->cubes[i], j)) break;
      }
      if(j == num_outputs) {
        FREE(ored);
        FREE(anded);
        return 1;
      }
    }
    OR_CUBE(ored, ored, cv->cubes[i], cube_size, j);
    AND_CUBE(anded, anded, cv->cubes[i], cube_size, j);    
  }

  ins = 1;
  for(j = 0;j < num_inputs;j++) {
    ret = VALUE_VAR(cv->cubes[cv->ncubes - 1], j);
    if(ret == 1) {
      cv->p1[j]++;
      ins = 0;
    }             
    else
      if(ret == 2) {
        cv->p0[j]++;
        ins = 0;
      }          
    if((cv->p0[j] == cv->ncubes) ||
       (cv->p1[j] == cv->ncubes)) {
      FREE(ored);
      FREE(anded);
      return 0;
    }
    if((cv->p0[j] > 0) && (cv->p1[j] > 0)) {
      nact++;
      if((cv->p1[j] + cv->p0[j]) > maxsum) {
        maxsum = cv->p1[j] + cv->p0[j];
        var = j;
      }
    }
    else
      if(((cv->p0[j] > 0) && (cv->p1[j] == 0)) ||
         ((cv->p0[j] == 0) && (cv->p1[j] > 0))) {
        nact++;
        nunate++;   
      }
    OR_CUBE(ored, ored, cv->cubes[cv->ncubes - 1], cube_size, i);
    AND_CUBE(anded, anded, cv->cubes[cv->ncubes - 1], cube_size, i);    
  }    
  if(ins == 1) {
    for(j = 0;j < num_outputs;j++) {
      if(!VALUE_OUT(cv->cubes[cv->ncubes - 1], j)) break;
    }
    if(j == num_outputs) {
      FREE(ored);
      FREE(anded);
      return 1;
    }
  }
  ins = 0;
  for(i = 0;i < num_outputs;i++) {
    if(!VALUE_OUT(ored, i)) {
      FREE(ored);
      FREE(anded);
      return 0;
    }
    if(!VALUE_OUT(anded, i)) {
      ins++;
    }
  }
  FREE(ored);
  FREE(anded);
  
  if(ins == 1) nunate++;
  if(ins > 0) nact++;  
  
  if(nunate == nact) {
    return 0;
  }

  if(nact == 1) {
    return 1;
  }

  nact = -1;
  for(i = 0;i < num_inputs;i++) {
    if(!((cv->p0[i] > 0) && (cv->p1[i] > 0))) {
      j = 0;
      while(j < cv->ncubes) {
        if(VALUE_VAR(cv->cubes[j], i) != 3) {
          pla_CubeDel(cv, j, cube_size);
          nact = 1;
        }
        else {
          j++;
        }        
      }
    }  
  }
  if(nact == 1) goto begin;
  
  plaCoverVarCofactorV2(&cv, var, &c1, num_inputs, cube_size);
  
  ret = plaTautology(cv, num_inputs, num_outputs, cube_size);
  if(ret == 0) {
    pla_CoverFree(c1);
    return 0;
  }
  ret = plaTautology(c1, num_inputs, num_outputs, cube_size);
  pla_CoverFree(c1);
  
  return ret;
}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/



