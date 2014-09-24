/**CFile***********************************************************************

  FileName    [plaIrred.c]

  PackageName [pla]

  Synopsis    [Irredundant Procedures]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

#include "ansi.h"
#include "util.h"
#include "plaInt.h"
#include "sparse.h"
#include "mincov.h"

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

static void Redundant(cover_t **E, cover_t **R, cover_t *F, cover_t *D, int num_inputs, int num_outputs, int cube_size);
static cover_t * PartiallyRedundant(cover_t *E, cover_t *R, cover_t *D, int num_inputs, int num_outputs, int cube_size);
static cover_t * MinimalIrredundant(cover_t *Rp, cover_t *E, cover_t *D, int num_inputs, int num_outputs, int cube_size);
static sm_matrix * NoCovermat(cover_t *Rp, cover_t *E, cover_t *D, int num_inputs, int num_outputs, int cube_size);
static void Ltaut(sm_matrix **B, cover_t *A, cover_t *Bc, int num_inputs, int num_outputs, int cube_size, int *pos_cube);

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
plaIrredundant(
  cover_t **F,
  cover_t *D,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t *E, *R, *Rp, *Rc;
int i;

  Redundant(&E, &R, *F, D, num_inputs, num_outputs, cube_size);
  if(R->ncubes == 0) {
    pla_CoverFree(R);
    pla_CoverFree(E);
    return;
  }
  Rp = PartiallyRedundant(E, R, D, num_inputs, num_outputs, cube_size);
  pla_CoverFree(R);
  if(Rp->ncubes == 0) {
    pla_CoverFree(E);
    return;
  }
  Rc = MinimalIrredundant(Rp, E, D, num_inputs, num_outputs, cube_size);
  pla_CoverFree(Rp);
  pla_CoverFree(*F);
  *F = Rc;
  for(i = 0;i < E->ncubes;i++) {
    pla_CubeAdd(*F, E->cubes[i]);
    (*F)->flags[(*F)->ncubes - 1] = E->flags[i];
  }
  FREE(E->flags);
  FREE(E->p0);
  FREE(E->p1);
  FREE(E->cubes);
  FREE(E);
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
static void
Redundant(
  cover_t **E,
  cover_t **R,
  cover_t *F,
  cover_t *D,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i, j, k;
unsigned int *c;
cover_t *c1;

  *E = pla_CoverAlloc(num_inputs);
  *R = pla_CoverAlloc(num_inputs);
  for(i = 0;i < F->ncubes;i++) {
    c1 = pla_CoverAlloc(num_inputs);
    for(j = 0;j < F->ncubes;j++) {
      if(j == i) continue;
      c = plaCubeCofactor(F->cubes[j], F->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(c1, c);
    }
    for(j = 0;j < D->ncubes;j++) {
      c = plaCubeCofactor(D->cubes[j], F->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(c1, c);
    }
    c = (unsigned int *) ALLOC(unsigned int, cube_size);
    COPY_CUBE(c, F->cubes[i], cube_size, k);
    if(plaTautology(c1, num_inputs, num_outputs, cube_size) == 1) {
      pla_CubeAdd(*R, c);
      (*R)->flags[(*R)->ncubes - 1] = F->flags[i];
    }
    else {
      pla_CubeAdd(*E, c);
      (*E)->flags[(*E)->ncubes - 1] = F->flags[i];      
    }
    pla_CoverFree(c1);
  }
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static cover_t *
PartiallyRedundant(
  cover_t *E,
  cover_t *R,
  cover_t *D,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i, j, k;
unsigned int *c;
cover_t *c1, *ret;

  ret = pla_CoverAlloc(num_inputs);
  
  for(i = 0;i < R->ncubes;i++) {
    c1 = pla_CoverAlloc(num_inputs);
    for(j = 0;j < E->ncubes;j++) {
      c = plaCubeCofactor(E->cubes[j], R->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(c1, c);
    }
    for(j = 0;j < D->ncubes;j++) {
      c = plaCubeCofactor(D->cubes[j], R->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(c1, c);
    }
    if(plaTautology(c1, num_inputs, num_outputs, cube_size) == 0) {
      c = (unsigned int *) ALLOC(unsigned int, cube_size);
      COPY_CUBE(c, R->cubes[i], cube_size, k);
      pla_CubeAdd(ret, c);
      ret->flags[ret->ncubes - 1] = R->flags[i];
    }
    pla_CoverFree(c1);
  }
  
  return ret;
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static cover_t *
MinimalIrredundant(
  cover_t *Rp,
  cover_t *E,
  cover_t *D,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t *Rc;
sm_matrix *B;
int *W, i;
sm_row *row;
sm_element *element;
unsigned int *c;

  Rc = pla_CoverAlloc(num_inputs);  
  B = NoCovermat(Rp, E, D, num_inputs, num_outputs, cube_size);
    
  W = (int *) ALLOC(int, Rp->ncubes);
  for(i = 0;i < Rp->ncubes;W[i] = 1,i++);
  row = sm_minimum_cover(B, W, 1, 0);
  FREE(W);
  sm_foreach_row_element(row, element) {
    c = (unsigned int *) ALLOC(unsigned int, cube_size);
    COPY_CUBE(c, Rp->cubes[element->col_num], cube_size, i);
    pla_CubeAdd(Rc, c);
    Rc->flags[Rc->ncubes - 1] = Rp->flags[element->col_num];
  }
  sm_free(B);
  
  return Rc;
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static sm_matrix *
NoCovermat(
  cover_t *Rp,
  cover_t *E,
  cover_t *D,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
sm_matrix *B;
int i, j, *pos_cube = NIL(int), pc = 0;
unsigned int *c;
cover_t *Bc, *A;

  B = sm_alloc();

  for(i = 0;i < Rp->ncubes;i++) {      
    Bc = pla_CoverAlloc(num_inputs);
    for(j = 0;j < E->ncubes;j++) {
      c = plaCubeCofactor(E->cubes[j], Rp->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(Bc, c);
    }
    for(j = 0;j < D->ncubes;j++) {
      c = plaCubeCofactor(D->cubes[j], Rp->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) pla_CubeAdd(Bc, c);
    }    
    A = pla_CoverAlloc(num_inputs);
    for(j = 0;j < Rp->ncubes;j++) {
      c = plaCubeCofactor(Rp->cubes[j], Rp->cubes[i], num_inputs, num_outputs,
        cube_size);
      if(c != NIL(unsigned int)) {
        pc++;
        pos_cube = (int *) REALLOC(int, pos_cube, pc);
        pla_CubeAdd(A, c);
        pos_cube[pc - 1] = j;
      }
    }   
    Ltaut(&B, A, Bc, num_inputs, num_outputs, cube_size, pos_cube);
    FREE(pos_cube);
    pos_cube = NIL(int);
    pc = 0;

    pla_CoverFree(Bc);
    pla_CoverFree(A);  
  }
    
  return B;
}


/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void
Ltaut(
  sm_matrix **B,
  cover_t *A,
  cover_t *Bc,
  int num_inputs,
  int num_outputs,
  int cube_size,
  int *pos_cube)
{
int var = -1, i, j, ins, maxsum = 0, nact = 0, nunate = 0, ret;
cover_t *A0, *A1, *Bc1;
int *cv2, *pos_cube0 = NIL(int), *pos_cube1 = NIL(int), np0 = 0, np1 = 0;
unsigned int *c;
  
begin:

  cv2 = (int *) ALLOC(int, A->ncubes);
  for(i = 0;i < num_inputs;i++) {
    A->p0[i] = A->p1[i] = 0;
  }  
  
  for(i = 0;i < Bc->ncubes;i++) {
    ins = 1;
    for(j = 0;j < num_inputs;j++) {
      ret = VALUE_VAR(Bc->cubes[i], j);    
      if(ret == 1) {
        A->p1[j]++;
        ins = 0;
      }             
      else
        if(ret == 2) {
          A->p0[j]++;
          ins = 0;
        }          
    }
    if(ins == 1) {
      for(j = 0;j < num_outputs;j++) {
        if(!VALUE_OUT(Bc->cubes[i], j)) break;
      }
      if(j == num_outputs) {
        FREE(cv2);
        return;
      }
    }
  }

  for(i = 0;i < A->ncubes - 1;i++) {
    cv2[i] = 0;
    ins = 1;
    for(j = 0;j < num_inputs;j++) {
      ret = VALUE_VAR(A->cubes[i], j);    
      if(ret == 1) {
        A->p1[j]++;
        ins = 0;
      }             
      else
        if(ret == 2) {
          A->p0[j]++;
          ins = 0;
        }          
    }
    if(ins == 1) {
      for(j = 0;j < num_outputs;j++) {
        if(!VALUE_OUT(A->cubes[i], j)) break;
      }
      if(j == num_outputs) {
        cv2[i] = 1;
      }
    }
  }
  cv2[A->ncubes - 1] = 0;
  ins = 1;
  for(j = 0;j < num_inputs;j++) {
    ret = VALUE_VAR(A->cubes[A->ncubes - 1], j);
    if(ret == 1) {
      A->p1[j]++;
      ins = 0;
    }             
    else
      if(ret == 2) {
        A->p0[j]++;
        ins = 0;
      }          
    if((A->p0[j] > 0) && (A->p1[j] > 0)) {
      nact++;
      if((A->p1[j] + A->p0[j]) > maxsum) {
        maxsum = A->p1[j] + A->p0[j];
        var = j;
      }
    }
    else
      if(((A->p0[j] > 0) && (A->p1[j] == 0)) ||
         ((A->p0[j] == 0) && (A->p1[j] > 0))) {
        nact++;
        nunate++;   
      }
  }    
  if(ins == 1) {
    for(j = 0;j < num_outputs;j++) {
      if(!VALUE_OUT(A->cubes[A->ncubes - 1], j)) break;
    }
    if(j == num_outputs) {
      cv2[A->ncubes - 1] = 1;
    }
  }  
  if(ins == 1) nunate++;
  if(ins > 0) nact++;  
  
  if(nunate == nact) {
    ret = (*B)->nrows;
    for(i = 0;i < A->ncubes;i++) {
      if(cv2[i] == 1) {
        (void) sm_insert(*B, ret, pos_cube[i]);
      }    
    }
    FREE(cv2);
    return;
  }

/*
  nact = -1;
  for(i = 0;i < num_inputs;i++) {
    if(!((A->p0[i] > 0) && (A->p1[i] > 0))) {
      j = 0;
      while(j < A->ncubes) {
        if(VALUE_VAR(A->cubes[j], i) != 3) {
          pla_CubeDel(A, j, cube_size);
          nact = 1;
        }
        else {
          j++;
        }        
      }
      j = 0;
      while(j < Bc->ncubes) {
        if(VALUE_VAR(Bc->cubes[j], i) != 3) {
          pla_CubeDel(Bc, j, cube_size);
          nact = 1;
        }
        else {
          j++;
        }        
      }
    }  
  }
  if(nact == 1) {
    FREE(cv2);
    goto begin;
  }
*/  

  A0 = pla_CoverAlloc(num_inputs);
  A1 = pla_CoverAlloc(num_inputs);
  for(i = 0;i < A->ncubes;i++) {
    switch(VALUE_VAR(A->cubes[i], var)) {
      case 1: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, A->cubes[i], cube_size, j);
        FULL_VAR(c, var);
        pla_CubeAdd(A1, c);
        np1++;
        pos_cube1 = (int *) REALLOC(int, pos_cube1, np1);
        pos_cube1[np1 - 1] = pos_cube[i];
      } break;

      case 2: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, A->cubes[i], cube_size, j);
        FULL_VAR(c, var);
        pla_CubeAdd(A0, c);
        np0++;
        pos_cube0 = (int *) REALLOC(int, pos_cube0, np0);
        pos_cube0[np0 - 1] = pos_cube[i];
      } break;

      case 3: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, A->cubes[i], cube_size, j);
        FULL_VAR(c, var);
        pla_CubeAdd(A1, c);
        np1++;
        pos_cube1 = (int *) REALLOC(int, pos_cube1, np1);
        pos_cube1[np1 - 1] = pos_cube[i];

        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, A->cubes[i], cube_size, j);
        FULL_VAR(c, var);
        pla_CubeAdd(A0, c);      
        np0++;
        pos_cube0 = (int *) REALLOC(int, pos_cube0, np0);
        pos_cube0[np0 - 1] = pos_cube[i];        
      } break;
    }  
  }
  plaCoverVarCofactorV2(&Bc, var, &Bc1, num_inputs, cube_size);

  Ltaut(B, A0, Bc, num_inputs, num_outputs, cube_size, pos_cube0);
  Ltaut(B, A1, Bc1, num_inputs, num_outputs, cube_size, pos_cube1);
  
  pla_CoverFree(A0);
  pla_CoverFree(A1);
  pla_CoverFree(Bc1);
  FREE(pos_cube0);
  FREE(pos_cube1);
}
