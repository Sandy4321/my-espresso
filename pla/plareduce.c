/**CFile***********************************************************************

  FileName    [plaReduce.c]

  PackageName [pla]

  Synopsis    [Reduce algorithm]

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

typedef struct ro ro_t;

/*---------------------------------------------------------------------------*/
/* Structure declarations                                                    */
/*---------------------------------------------------------------------------*/

struct ew {
  int value;
  unsigned int *cube;
};

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

static unsigned int * SCCCU(cover_t *cv, int num_inputs, int num_outputs, int cube_size);
static void ReduceOrder(cover_t *F, int num_inputs, int cube_size);
static int ReduceCompare(const void *c1, const void *c2);

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
plaReduce(
  cover_t **F,
  cover_t *D,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i, j, k, ins;
unsigned int *c, *c2;
cover_t *tmp;

  ReduceOrder(*F, num_inputs, cube_size);

  for(i = 0;i < (*F)->ncubes;i++) {
    tmp = pla_CoverAlloc(num_inputs);
    for(j = 0;j < (*F)->ncubes;j++) {
      if(i == j) continue;
      if(VALUE_VAR((*F)->cubes[j], 0) == 0) continue;      
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
      FREE(c);
      ins = 0;
      for(k = 0;k < num_inputs;k++) {
        if(VALUE_VAR(c2, k) != VALUE_VAR((*F)->cubes[i], k)) {
          ins = 1;
          break;
        }      
      }
      if(ins == 0) {
        for(k = 0;k < num_outputs;k++) {
          if((VALUE_OUT(c2, k) && !VALUE_OUT((*F)->cubes[i], k)) ||
             (!VALUE_OUT(c2, k) && VALUE_OUT((*F)->cubes[i], k))) {
            ins = 1;
            break;  
          }
        }
      }
      if(ins == 1) {
        UNSET_PRIME((*F), i);
        FREE((*F)->cubes[i]);
        (*F)->cubes[i] = c2;
      }
      else {
        FREE(c2);
      }
    }
    pla_CoverFree(tmp);
  }  
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

/**Function********************************************************************

  Synopsis           [Main function]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
unsigned int *
plaSCCC(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
unsigned int *cube1, *cube0;
int var, i;
cover_t *c1;

  var = plaBinateSelect(cv, num_inputs);
  if(var == -1) {
    cube0 = SCCCU(cv, num_inputs, num_outputs, cube_size);
    return cube0;
  }  
  plaCoverVarCofactorV2(&cv, var, &c1, num_inputs, cube_size);

  cube0 = plaSCCC(cv, num_inputs, num_outputs, cube_size);
  cube1 = plaSCCC(c1, num_inputs, num_outputs, cube_size);

  if((cube0 != NIL(unsigned int)) && (cube1 != NIL(unsigned int))) {
    NEG_VAR(cube0, var);
    POS_VAR(cube1, var);
    OR_CUBE(cube0, cube0, cube1, cube_size, i);    
    FREE(cube1);
  }
  else
  if(cube0 != NIL(unsigned int)) {
    NEG_VAR(cube0, var);    
  }
  else
  if(cube1 != NIL(unsigned int)) {
    POS_VAR(cube1, var);
    cube0 = cube1;
  }
  pla_CoverFree(c1);     
  return cube0;
}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis           [Main function]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static unsigned int *
SCCCU(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i, *N, *nN, j, k, tot1 = 0, tot2 = 0;
unsigned int *ret, *tmp1, *tmp2;

  ret = pla_CubeNew(cube_size, 1);

  if(num_outputs == 0) {
    for(i = 0;i < num_inputs;i++) {
      for(j = 0;j < cv->ncubes;j++) {
        for(k = 0;k < num_inputs;k++) {
          if(k != i && VALUE_VAR(cv->cubes[j], k) != 3) break;
        }
        if(k == num_inputs) {
          if(VALUE_VAR(cv->cubes[j], i) == 1) {
            NEG_VAR(ret, i);
          }
          else {
            POS_VAR(ret, i);
          }
          break;
        }        
      }
    }
    return ret;
  }

  N = (int *) ALLOC(int *, num_outputs);
  nN = (int *) ALLOC(int *, num_outputs);
  for(i = 0;i < num_outputs;i++) {
    N[i] = 0;
    nN[i] = 0;
  }

  tmp1 = pla_CubeNew(cube_size, 1);
  tmp2 = pla_CubeNew(cube_size, 1);
  for(i = 0;i < num_outputs;i++) {
    NEG_OUT(tmp1, i);
    NEG_OUT(tmp2, i);    
  }
  POS_VAR(tmp1, 0);
  NEG_VAR(tmp2, 0);
  for(j = 0;j < cv->ncubes;j++) {
    if(plaCubeContain(cv->cubes[j], tmp1, num_inputs, num_outputs)) {
      for(k = 0;k < num_outputs;k++) {
        if(VALUE_OUT(cv->cubes[j], k) && !N[k]) {
          N[k] = 1;
          tot1++;          
        }
      }
    }
    else
      if(plaCubeContain(cv->cubes[j], tmp2, num_inputs, num_outputs)) {
        for(k = 0;k < num_outputs;k++) {
          if(VALUE_OUT(cv->cubes[j], k) && !nN[k]) {
            nN[k] = 1;
            tot2++;
          }
        }        
      }
  }
  if(tot1 == num_outputs) {
    NEG_VAR(ret, 0);
  }
  else
  if(tot2 == num_outputs) {
    POS_VAR(ret, 0);
  }
  for(i = 1;i < num_inputs;i++) {
    tot1 = 0;
    tot2 = 0;
    FULL_VAR(tmp1, i - 1);
    FULL_VAR(tmp2, i - 1);
    POS_VAR(tmp1, i);
    NEG_VAR(tmp2, i);
    for(j = 0;j < cv->ncubes;j++) {
      if(plaCubeContain(cv->cubes[j], tmp1, num_inputs, num_outputs)) {
        for(k = 0;k < num_outputs;k++) {
          if(VALUE_OUT(cv->cubes[j], k) && !N[k]) {
            N[k] = 1;
            tot1++;          
          }
        }
      }
      else
        if(plaCubeContain(cv->cubes[j], tmp2, num_inputs, num_outputs)) {
          for(k = 0;k < num_outputs;k++) {
            if(VALUE_OUT(cv->cubes[j], k) && !nN[k]) {
              nN[k] = 1;
              tot2++;
            }
          }        
        }
    }
    if(tot1 == num_outputs) {
      NEG_VAR(ret, 0);
    }
    else
    if(tot2 == num_outputs) {
      POS_VAR(ret, 0);
    }    
  }
  FREE(tmp2);
  FULL_VAR(tmp1, i - 1);
  POS_OUT(tmp1, 0);
  for(j = 0;j < cv->ncubes;j++) {
    if(plaCubeContain(cv->cubes[j], tmp1, num_inputs, num_outputs)) {
      NEG_OUT(ret, 0);
      break;  
    }
  }    
  for(i = 1;i < num_outputs;i++) {
    NEG_OUT(tmp1, i - 1);
    POS_OUT(tmp1, i);
    for(j = 0;j < cv->ncubes;j++) {
      if(plaCubeContain(cv->cubes[j], tmp1, num_inputs, num_outputs)) {
        NEG_OUT(ret, i);
        break;  
      }
    }    
  }
  FREE(tmp1);
  FREE(N);
  FREE(nN);

  return ret;
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void
ReduceOrder(
  cover_t *F,
  int num_inputs,
  int cube_size)
{
int *cs, i, j;
struct ew **w;

  cs = (int *) ALLOC(int, 2*num_inputs);
  for(i = 0;i < 2*num_inputs;i++) {
    cs[i] = 0;
  }
  for(i = 0;i < F->ncubes;i++) {
    for(j = 0;j < num_inputs;j++) {
      switch(VALUE_VAR(F->cubes[i], j)) {
        case 3: cs[2*j]++;
        case 1: cs[2*j + 1]++; break;
        case 2: cs[2*j]++; break;
        default:
      }
    }
  }
  w = (struct ew **) ALLOC(struct ew *, F->ncubes);
  for(i = 0;i < F->ncubes;i++) {
    w[i] = (struct ew *) ALLOC(struct ew, 1);
    w[i]->value = 0;
    w[i]->cube = F->cubes[i];
    for(j = 0;j < num_inputs;j++) {
      switch(VALUE_VAR(F->cubes[i], j)) {
        case 3: w[i]->value += cs[2*j];
        case 1: w[i]->value += cs[2*j + 1]; break;
        case 2: w[i]->value += cs[2*j]; break;
      }
    }  
  }    
  FREE(cs);  

  qsort(w, F->ncubes, sizeof(struct ew *), ReduceCompare);
  
  for(i = 0;i < F->ncubes;i++) {
    F->cubes[i] = w[i]->cube;    
    FREE(w[i]);
  }
  FREE(w);
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static int
ReduceCompare(
  const void *c1,
  const void *c2)
{
struct ew **p1, **p2;

  p1 = (struct ew **) c1;
  p2 = (struct ew **) c2;
  
  return ((*p2)->value - (*p1)->value);
}
