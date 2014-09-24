/**CFile***********************************************************************

  FileName    [plaComp.c]

  PackageName [pla]

  Synopsis    [Complement]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

#include "ansi.h"
#include "util.h"
#include "plaInt.h"
#include <math.h>

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

static cover_t * Complement(cover_t *cv, int num_inputs, int num_outputs, int cube_size);

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
void 
pla_Complement(
  cover_t *F,
  cover_t *D,
  cover_t **R,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t **outputs, *tmp, *cv;
int i, j, k;

  cv = pla_CoverAlloc(num_inputs);
  for(i = 0;i < F->ncubes;i++) {
    pla_CubeAdd(cv, F->cubes[i]);  
  }
  for(i = 0;i < D->ncubes;i++) {
    pla_CubeAdd(cv, D->cubes[i]);  
  }
  pla_CoverFree(*R);

  if(num_outputs <= 1) {
    *R = Complement(cv, num_inputs, num_outputs, cube_size);
  }
  else {
    outputs = (cover_t **) ALLOC(cover_t *, num_outputs);
    for(i = 0;i < num_outputs;i++) {
      outputs[i] = pla_CoverAlloc(num_inputs);
    }
    for(i = 0;i < cv->ncubes;i++) {
      for(j = 0;j < num_outputs;j++) {
        if(VALUE_OUT(cv->cubes[i], j)) {
          pla_CubeAdd(outputs[j], cv->cubes[i]);
          break;
        }
      }
    }
    *R = pla_CoverAlloc(num_inputs);
    
    for(i = 0;i < num_outputs;i++) {
      tmp = Complement(outputs[i], num_inputs, num_outputs, cube_size);
      for(j = 0;j < tmp->ncubes;j++) {
        for(k = 0;k < num_outputs;k++) {
          if(k != i) {
           NEG_OUT(tmp->cubes[j], k);
          }
        }        
        pla_CubeAdd(*R, tmp->cubes[j]);
      }      
      FREE(tmp->flags);
      FREE(tmp->p0);
      FREE(tmp->p1);
      FREE(tmp->cubes);
      FREE(tmp);
    }        
    for(i = 0;i < num_outputs;i++) {    
      FREE(outputs[i]->flags);
      FREE(outputs[i]->p0);
      FREE(outputs[i]->p1);
      FREE(outputs[i]->cubes);
      FREE(outputs[i]);
    }
    FREE(outputs);
  }
  FREE(cv->flags);
  FREE(cv->p0);
  FREE(cv->p1);
  FREE(cv->cubes);
  FREE(cv);
}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static cover_t *
Complement(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t *ret, *c0, *c1, *retc0, *retc1, *cuberet;
int var = -1, i, j, maxsum = 0, ins, t, go = 0;
unsigned int *c;

  if(cv->ncubes == 0) {
    ret = pla_CoverAlloc(num_inputs);
    c = pla_CubeNew(cube_size, 1);
    pla_CubeAdd(ret, c);
    return ret;
  }

  for(i = 0;i < num_inputs;cv->p0[i] = cv->p1[i] = 0,i++);

  for(i = 0;i < cv->ncubes - 1;i++) {
    ins = 1;
    for(j = 0;j < num_inputs;j++) {
      t = VALUE_VAR(cv->cubes[i], j);
      if(t == 1) {
        cv->p1[j]++;
        ins = 0;
      }
      else
        if(t == 2) {
          cv->p0[j]++;
          ins = 0;
        }
    }
    if(ins == 1) {
      ret = pla_CoverAlloc(num_inputs);
      return ret;
    }
  }  

  c = pla_CubeNew(cube_size, 1);

  ins = 1;
  for(j = 0;j < num_inputs;j++) {
    t = VALUE_VAR(cv->cubes[cv->ncubes - 1], j);
    if(t == 1) {
      cv->p1[j]++;
      ins = 0;
    }
    else
      if(t == 2) {
        cv->p0[j]++;
        ins = 0;
      }
    if(((cv->p1[j] + cv->p0[j]) > maxsum) && (cv->p0[j] > 0) && 
        (cv->p1[j] > 0)) {
      maxsum = cv->p1[j] + cv->p0[j];
      var = j;
    }
    if(cv->p0[j] == cv->ncubes) {
      NEG_VAR(c, j);
      go = 1;
    }
    else
      if(cv->p1[j] == cv->ncubes) {
        POS_VAR(c, j);
        go = 1;
      }
  }

  if(ins == 1) {
    FREE(c);
    ret = pla_CoverAlloc(num_inputs);
    return ret;
  }

  if(var == -1) {
    ret = UnateComplement(cv, num_inputs, cube_size);
    FREE(c);    
    return ret;
  }    

  if(go == 1) {
    c0 = pla_CoverAlloc(num_inputs);
    pla_CubeAdd(c0, c);
    cuberet = UnateComplement(c0, num_inputs, cube_size);
    cv = plaCoverCubeCofactor(cv, c, num_inputs, num_outputs, cube_size);
    pla_CoverFree(c0);
  }
  else {
    FREE(c); 
  }

  plaCoverVarCofactor(cv, var, &c0, &c1, num_inputs, cube_size);

  retc0 = Complement(c0, num_inputs, num_outputs, cube_size);
  retc1 = Complement(c1, num_inputs, num_outputs, cube_size);

  ret = plaMwc(retc0, retc1, var, num_inputs, num_outputs, cube_size, 0);
  if(go == 1) {
    for(i = 0;i < cuberet->ncubes;i++) {
      pla_CubeAdd(ret, cuberet->cubes[i]);
    }
    FREE(cuberet->flags);
    FREE(cuberet->p0);
    FREE(cuberet->p1);      
    FREE(cuberet->cubes);
    FREE(cuberet);
    pla_CoverFree(cv);      
  }

  return ret;
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
cover_t *
UnateComplement(
  cover_t *cv,
  int num_inputs,
  int cube_size)
{
cover_t *ret, *c0, *c1, *retc0, *retc1;
int i,j,t,n2,csel,maxn2 = 0,ins,var;
unsigned int *c;

  if(cv->ncubes == 0) {
    ret = pla_CoverAlloc(num_inputs);
    c = pla_CubeNew(cube_size, 1);
    pla_CubeAdd(ret, c);
    return ret;
  }

  for(i = 0;i < num_inputs;i++) {
    cv->p0[i] = 0;
    cv->p1[i] = 0;
  }

  for(i = 0;i < cv->ncubes;i++) {
    ins = 1;
    n2 = 0;
    for(j = 0;j < num_inputs;j++) {
      t = VALUE_VAR(cv->cubes[i], j);
      if(t == 1) {
        cv->p1[j]++;
        ins = 0;
      }
      else
        if(t == 2) {
          cv->p0[j]++; 
          ins = 0;
        }
        else {
          n2++;
        }
    }
    if(ins) {
      ret = pla_CoverAlloc(num_inputs);
      return ret;
    }
    if(n2 > maxn2) {
      maxn2 = n2;
      csel = i; 
    }
  }
  
  if(cv->ncubes == 1) {
    ret = pla_CoverAlloc(num_inputs);
    for(i = 0;i < num_inputs;i++) {
      t = VALUE_VAR(cv->cubes[0], i);
      if(t == 1) {
        c = pla_CubeNew(cube_size, 1);
        NEG_VAR(c, i);
        pla_CubeAdd(ret, c);
      }
      else
        if(t == 2) {
          c = pla_CubeNew(cube_size, 1);
          POS_VAR(c, i);
          pla_CubeAdd(ret, c); 
        }
    }  
    return ret;    
  }
  
  maxn2 = 0;
  for(i = 0;i < num_inputs;i++) {
    if(VALUE_VAR(cv->cubes[csel], i) != 3) {
      if(cv->p0[i] > maxn2) {
        maxn2 = cv->p0[i];
        var = i;
      }
      else
        if(cv->p1[i] > maxn2) {
          maxn2 = cv->p1[i];
          var = i;
        }
    }  
  }
  
  c0 = pla_CoverAlloc(num_inputs);
  c1 = pla_CoverAlloc(num_inputs);
  
  for(i = 0;i < cv->ncubes;i++) {
    if(VALUE_VAR(cv->cubes[i], var) == 3) {
      c = pla_CubeNew(cube_size, 0);
      COPY_CUBE(c, cv->cubes[i], cube_size, j);
      pla_CubeAdd(c0, c);
    }
    c = pla_CubeNew(cube_size, 0);
    COPY_CUBE(c, cv->cubes[i], cube_size, j);
    FULL_VAR(c, var);
    pla_CubeAdd(c1, c);
  }
     
  retc0 = UnateComplement(c0, num_inputs, cube_size);
  retc1 = UnateComplement(c1, num_inputs, cube_size);
  
  pla_CoverFree(c0);
  pla_CoverFree(c1);
  
  if(cv->p0[var] > 0) {
    for(i = 0;i < retc0->ncubes;i++) {
      POS_VAR(retc0->cubes[i], var);      
    }
  }
  else {
    for(i = 0;i < retc0->ncubes;i++) {
      NEG_VAR(retc0->cubes[i], var);      
    }    
  }
  for(i = 0;i < retc1->ncubes;i++) {
    pla_CubeAdd(retc0, retc1->cubes[i]);
  }
  FREE(retc1->flags);
  FREE(retc1->cubes);
  FREE(retc1->p0);
  FREE(retc1->p1);
  FREE(retc1);

  ret = retc0;
  
  return ret;  
}
