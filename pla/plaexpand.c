/**CFile***********************************************************************

  FileName    [plaExpand.c]

  PackageName [pla]

  Synopsis    [Expand]

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

static unsigned int * Expand(unsigned int *cube, cover_t *ON, cover_t *OFF, int num_inputs, int num_outputs, int cube_size, int below, short type, int *nc);
static void ExpandOrder(cover_t *F, int num_inputs, int cube_size);
static int ExpandCompare(const void *c1, const void *c2);

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
plaExpand(
  cover_t **F,
  cover_t *R,
  int num_inputs,
  int num_outputs,
  int cube_size,
  short type,
  short gasp)
{
cover_t *ret;
unsigned int *c;
int i, nc, k, j;

  ret = pla_CoverAlloc(num_inputs);
  i = 0;
  while(i < (*F)->ncubes) {
    if(IS_PRIME((*F), i)) {
      c = (unsigned int *) ALLOC(unsigned int, cube_size);
      COPY_CUBE(c, (*F)->cubes[i], cube_size, k);
      pla_CubeAdd(ret, c);
      SET_PRIME(ret, ret->ncubes - 1);
      pla_CubeDel((*F), i, cube_size);
    }
    else {
      i++;
    }
  }
  ExpandOrder(*F, num_inputs, cube_size);
  
  nc = (*F)->ncubes;
  for(i = 0;i < nc;i++) {
    if(VALUE_VAR((*F)->cubes[i], 0) == 0) continue;
    c = Expand((*F)->cubes[i], (*F), R, num_inputs, num_outputs, cube_size, i,
      type, &j);
    if(!gasp || (gasp && (j > 1))) {
      pla_CubeAdd(ret, c);
      SET_PRIME(ret, ret->ncubes - 1);
    }
  }
  pla_CoverFree(*F);
  *F = ret;
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
static unsigned int *
Expand(
  unsigned int *cube,
  cover_t *ON,
  cover_t *OFF,
  int num_inputs,
  int num_outputs,
  int cube_size,
  int below,
  short type,
  int *nc)
{
sm_matrix *B, *C;
int i,j,t1,t2,ncols = num_inputs + num_outputs;
int *L,lc=0,rc=0, feaseble_covered, covered, *candidates, ncand;
unsigned int *ret;
sm_col *col, *col2, *col3;
sm_row *row, *row2, *raux, *row3;
sm_element *element, *eaux, *e2, *e3;
int *W;

  *nc = 0;
  L = (int *) ALLOC(int, ncols);
  for(i = 0;i < ncols;i++) {
    L[i] = 0;
  }
 
  B = sm_alloc();
  for(i = 0;i < OFF->ncubes;i++) {
    for(j = 0;j < num_inputs;j++) {
      t1 = VALUE_VAR(cube, j);
      t2 = VALUE_VAR(OFF->cubes[i], j);
      if(((t1 == 1) && (t2 == 2)) ||
         ((t1 == 2) && (t2 == 1))) {
        (void) sm_insert(B, i, j);
      }
    }
    for(;j < ncols;j++) {
      if((!VALUE_OUT(cube, j - num_inputs)) &&
         (VALUE_OUT(OFF->cubes[i], j - num_inputs))) {
        (void) sm_insert(B, i, j);
      }
    }
  }
  
  C = sm_alloc();
  for(i = below;i < ON->ncubes;i++) {
    if(VALUE_VAR(ON->cubes[i], 0) == 0) continue;
    for(j = 0;j < num_inputs;j++) {
      t1 = VALUE_VAR(cube, j);
      t2 = VALUE_VAR(ON->cubes[i], j);
      if(((t1 == 1) && (t2 != 1)) ||
         ((t1 == 2) && (t2 != 2))) {
        (void) sm_insert(C, i, j);
      }
    }
    for(;j < ncols;j++) {
      if((!VALUE_OUT(cube, j - num_inputs)) &&
         (VALUE_OUT(ON->cubes[i], j - num_inputs))) {
        (void) sm_insert(C, i, j);
      }
    }
  }
  
  while(((lc + rc) < ncols) && (sm_num_elements(B) > 0) &&
        (sm_num_elements(C) > 0)) {           
    row = B->first_row;
    while(row != NIL(sm_row)) {
      raux = row->next_row;
      if(row->length == 1) {
        col = sm_get_col(B, row->first_col->col_num);
        t1 = col->col_num;
        sm_foreach_col_element(col, element) {
          if((raux != NIL(sm_row)) && (raux->row_num == element->row_num))
            raux = raux->next_row;
          row2 = sm_get_row(B, element->row_num);
          sm_foreach_row_element(row2, e2) {
            if(e2->col_num != element->col_num) {
              col2 = sm_get_col(B, e2->col_num);
              if(col2->length == 1) {
                rc++;
                col3 = sm_get_col(C, col2->col_num);
                sm_foreach_col_element(col3, e3) {
                  row3 = sm_get_row(C, e3->row_num);
                  if(row3->length == 1) ZERO_VAR(ON->cubes[row3->row_num], 0);
                  (*nc)++;
                }
                sm_delcol(C, col2->col_num);
              }
            }
          }
          sm_delrow(B, element->row_num);
        }
        L[t1] = 1;
        lc++;
        col = sm_get_col(C, t1);
        sm_foreach_col_element(col, element) {
          sm_delrow(C, element->row_num);
        }
        if((sm_num_elements(B) == 0) || (sm_num_elements(C) == 0)) break;
      }
      row = raux;
    }
    if((sm_num_elements(B) == 0) || (sm_num_elements(C) == 0)) break;

    candidates = (int *) ALLOC(int, C->nrows);
    ncand = 0;
        
    sm_foreach_row(C, row) {
      feaseble_covered = 1;
      sm_foreach_row(B, row2) {
        covered = 0;
        sm_foreach_row_element(row2, element) {
          if(sm_row_find(row, element->col_num) == NIL(sm_element)) {
            covered = 1;
          }
          if(covered == 1) break;
        }
        if(covered == 0) {
          feaseble_covered = 0;        
          break;
        }
      }  
      if(feaseble_covered == 1) {
        candidates[ncand] = row->row_num;
        ncand++;
      }
    }

    if(ncand == 0) {
      covered = 0;
      sm_foreach_col(C, col) {
        if(col->length > covered) {
          t1 = col->col_num;
          covered = col->length;
        }
      }
      rc++;
      col = sm_get_col(C, t1);
      sm_foreach_col_element(col, element) {
        row = sm_get_row(C, element->row_num);
        if(row->length == 1) ZERO_VAR(ON->cubes[row->row_num], 0);
        (*nc)++;
      }
      sm_delcol(C, t1);
      sm_delcol(B, t1);
    }
    else {
      feaseble_covered = 0;
      t2 = 10000;
      for(i = 0;i < ncand;i++) {
        covered = 0;
        row = sm_get_row(C, candidates[i]);        
        for(j = 0;j < ncand;j++) {
          if(i == j) continue;
          row2 = sm_get_row(C, candidates[j]);
          if(sm_row_contains(row2, row)) {
            covered++;
          }
        }
        if((covered > feaseble_covered) ||
           ((covered == feaseble_covered) && (row->length < t2))) {
          feaseble_covered = covered;
          t1 = candidates[i];
          t2 = row->length;
        }
      }
      
      row = sm_get_row(C, t1);
      element = row->first_col;
      while(element != NIL(sm_element)) {
        eaux = element->next_col;
        t2 = element->col_num;
        rc++;      
        col = sm_get_col(C, t2);
        sm_foreach_col_element(col, e2) {
          row = sm_get_row(C, e2->row_num);
          if(row->length == 1) ZERO_VAR(ON->cubes[row->row_num], 0);
          (*nc)++;
        }
        sm_delcol(C, t2);
        sm_delcol(B, t2);
        element = eaux;
        if((sm_num_elements(B) == 0) || (sm_num_elements(C) == 0)) break;
      }
    }
    FREE(candidates);
  }

  if((sm_num_elements(B) > 0) && (sm_num_elements(C) == 0)) {
    W = (int *) ALLOC(int, ncols);
    for(i = 0;i < ncols;W[i] = 1,i++);
    row = sm_minimum_cover(B, W, 1, 0);
    sm_foreach_row_element(row, element) {
      L[element->col_num] = 1;
    }
    FREE(W);
  }

  ret = pla_CubeNew(cube_size, 1);
  for(i = 0;i < num_inputs;i++) {
    if(L[i]) {
      switch(VALUE_VAR(cube, i)) {
        case 1: POS_VAR(ret, i); break;
        case 2: NEG_VAR(ret, i); break;
      }
    }    
  }
  for(;i < ncols;i++) {
    if(L[i]) {
      if(!VALUE_OUT(cube, i - num_inputs)) {
        NEG_OUT(ret, i - num_inputs);
      }
    }
  }

  if(sm_num_elements(C) > 0) {
    sm_foreach_row(C, row) {
      for(j = 0;j < ncols;j++) {
        if(L[j]) {
          if(sm_row_find(row, j) != NIL(sm_element)) {
            break;
          }
        }
      }
      if(j == ncols) {
        ZERO_VAR(ON->cubes[row->row_num], 0);
        (*nc)++;
      }      
    }    
  }
  FREE(L);
  sm_free(B);   
  sm_free(C);
  
  return ret;
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void
ExpandOrder(
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

  qsort(w, F->ncubes, sizeof(struct ew *), ExpandCompare);
  
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
ExpandCompare(
  const void *c1,
  const void *c2)
{
struct ew **p1, **p2;

  p1 = (struct ew **) c1;
  p2 = (struct ew **) c2;
  
  return ((*p1)->value - (*p2)->value);
}
