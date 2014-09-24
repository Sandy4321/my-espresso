/**CFile***********************************************************************

  FileName    [pla.c]

  PackageName [pla]

  Synopsis    [My implementation of espresso II]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <math.h>
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

double expand_time, irred_time, reduce_time, gasp_time;
double tot_time, comp_time, essen_time;

int ni;
int no;
int cs;

struct timeval  time_start, time_stop;
struct rusage   ru_start, ru_stop;
double start, stop, seconds;

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

static int cube_compare(const void *c1, const void *c2);
static void StartTimer();
static void StopTimer();
static double GetRtime();

/**AutomaticEnd***************************************************************/

/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis           [Allocate memory for a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_Espresso(
  pla_t *pla,
  short stats)
{
int count, bestcost, i, essen, besttotal;
unsigned int *c;

  GLOBAL_NI = pla->num_inputs;

  comp_time = essen_time = expand_time = gasp_time = reduce_time = irred_time = 0;

  plaUnwarpOutput(pla->F, pla->cube_size, pla->num_outputs);
  
  EXECUTE( pla_Complement(pla->F, pla->D, &pla->R, pla->num_inputs, pla->num_outputs,
    pla->cube_size), comp_time);

  plaCost(pla->F, pla->num_inputs, pla->num_outputs);
  if(stats == 1) {
    fprintf(stdout, "Cost = %d  in = %d  out = %d\n", pla->F->ncubes,
      pla->F->cost_in, pla->F->cost_out);
  }

  EXECUTE( plaExpand(&pla->F, pla->R, pla->num_inputs, pla->num_outputs,
    pla->cube_size, MINCUBE, 0), expand_time);
  EXECUTE( plaIrredundant(&pla->F, pla->D, pla->num_inputs, pla->num_outputs,
    pla->cube_size), irred_time); 
  EXECUTE( plaEssencial(&pla->F, &pla->D, &essen, pla->num_inputs, pla->num_outputs,
    pla->cube_size), essen_time);

  do {
    do {
      bestcost = pla->F->ncubes;
      EXECUTE( plaReduce(&pla->F, pla->D, pla->num_inputs, pla->num_outputs,
        pla->cube_size), reduce_time);      
      EXECUTE( plaExpand(&pla->F, pla->R, pla->num_inputs, pla->num_outputs,
        pla->cube_size, MINCUBE, 0), expand_time);           
      EXECUTE( plaIrredundant(&pla->F, pla->D, pla->num_inputs, pla->num_outputs,
        pla->cube_size), irred_time);
    } while(pla->F->ncubes < bestcost);
      
    plaCost(pla->F, pla->num_inputs, pla->num_outputs);

    bestcost = pla->F->ncubes;
    besttotal = bestcost + pla->F->cost_in + pla->F->cost_out;
    
    EXECUTE(plaGasp(&pla->F, pla->D, pla->R, pla->num_inputs, pla->num_outputs,
      pla->cube_size), gasp_time);

    plaCost(pla->F, pla->num_inputs, pla->num_outputs);
    
  } while((pla->F->ncubes < bestcost) || ((pla->F->ncubes == bestcost) &&
    ((pla->F->ncubes + pla->F->cost_in + pla->F->cost_out) < besttotal)));

  for(count = essen;count < pla->D->ncubes;count++) {
    c = (unsigned int *) ALLOC(unsigned int, pla->cube_size);
    COPY_CUBE(c, pla->D->cubes[count], pla->cube_size, i);
    pla_CubeAdd(pla->F, c);
  }   
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_Print(
  FILE *fp,
  pla_t *pla,
  short type)
{
  if(type != 2) {
    fprintf(fp,".i %d\n", pla->num_inputs);
    fprintf(fp,".o %d\n", pla->num_outputs);  
  }
  
  plaCoverPrint(fp, pla->F, pla->num_inputs, pla->num_outputs);
  if(type == 1) {
    fprintf(fp,"DC\n");
    plaCoverPrint(fp, pla->D, pla->num_inputs, pla->num_outputs);
    fprintf(fp,"OFF\n");
    plaCoverPrint(fp,pla->R, pla->num_inputs, pla->num_outputs);
  }
  if(type != 2) {
    fprintf(fp,".e\n");
  }
}

/**Function********************************************************************

  Synopsis           [Allocate memory for a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
pla_t *
pla_Alloc(void)
{
pla_t *ret;

  ret = (pla_t *) ALLOC(pla_t, 1);
  ret->num_inputs = 0;
  ret->num_outputs = 0;
  ret->labels = NIL(char *);
  ret->cube_size = 0;

  ret->F = (cover_t *) ALLOC(cover_t, 1);
  ret->F->nalloc = 0;
  ret->F->ncubes = 0;
  ret->F->cubes = NIL(unsigned int *);
  ret->F->flags = NIL(short);
  ret->F->cost_in = 0;
  ret->F->cost_out = 0;
  ret->F->p0 = NIL(int);
  ret->F->p1 = NIL(int);

  ret->D = (cover_t *) ALLOC(cover_t, 1);
  ret->D->nalloc = 0;
  ret->D->ncubes = 0;
  ret->D->cubes = NIL(unsigned int *); 
  ret->D->flags = NIL(short);
  ret->D->cost_in = 0;
  ret->D->cost_out = 0;
  ret->D->p0 = NIL(int);
  ret->D->p1 = NIL(int);

  ret->R = (cover_t *) ALLOC(cover_t, 1);
  ret->R->nalloc = 0;
  ret->R->ncubes = 0;
  ret->R->cubes = NIL(unsigned int *);
  ret->R->flags = NIL(short);
  ret->R->cost_in = 0;
  ret->R->cost_out = 0;
  ret->R->p0 = NIL(int);
  ret->R->p1 = NIL(int);

  return ret;
}

/**Function********************************************************************

  Synopsis           [Allocate memory for a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
unsigned int *
pla_CubeNew(
  int cube_size,
  short type)
{
unsigned int *cube;
int i;

  cube = (unsigned int *) ALLOC(unsigned int, cube_size); 
  if(type == 0) {
    RESET_CUBE(cube, cube_size, i);
  }
  else {
    FULL_CUBE(cube, cube_size, i);
  }

  return cube;
}

/**Function********************************************************************

  Synopsis           [Add a cube to a pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_CubeAdd(
  cover_t *cv,
  unsigned int *cube)
{
  if(cv->ncubes == cv->nalloc) {
    cv->nalloc++;
    cv->ncubes++;
    cv->cubes = (unsigned int **) REALLOC(unsigned int *, cv->cubes, 
      cv->nalloc);
    cv->flags = (short *) REALLOC(short, cv->flags, cv->nalloc);
  }
  else {
    cv->ncubes++;
  }
  cv->cubes[cv->ncubes - 1] = cube;  
  RESET_FLAG(cv, cv->ncubes - 1);
}

/**Function********************************************************************

  Synopsis           [Allocate memory for a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_Free(pla_t *pla)
{
int i;

  if(pla->labels != NIL(char *)) {
    for(i = 0;i < pla->num_inputs + pla->num_outputs;i++) {
      FREE(pla->labels[i]);
    }
    FREE(pla->labels);
  }
  pla_CoverFree(pla->F);
  pla_CoverFree(pla->D);
  pla_CoverFree(pla->R);
  
  FREE(pla);
}

/**Function********************************************************************

  Synopsis           [Allocate memory for a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_CoverFree(cover_t *cv)
{
int i;

  if(cv->nalloc > 0) {
    for(i = 0;i < cv->nalloc;i++) FREE(cv->cubes[i]);
    FREE(cv->cubes);
    FREE(cv->flags);
  }
  FREE(cv->p0);
  FREE(cv->p1);
  FREE(cv);
}

/**Function********************************************************************

  Synopsis           [Make the cover minimal w.r.t single cube containment]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_Scc(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i, j;

  plaOrder(cv, num_inputs, num_outputs, cube_size);

  for(i = 0;i < cv->ncubes - 1;i++) {
    if(VALUE_VAR(cv->cubes[i], 0) == 0) continue;
    for(j = i+1;j < cv->ncubes;j++) {
      if(plaCubeContain(cv->cubes[i], cv->cubes[j], num_inputs, num_outputs)) {
        ZERO_VAR(cv->cubes[j], 0);
      }        
    }
  }
  i = 0;
  while(i < cv->ncubes) {
    if(VALUE_VAR(cv->cubes[i], 0) == 0) {
      pla_CubeDel(cv, i, cube_size);
    }        
    else {
      i++;
    }
  }
}

/**Function********************************************************************

  Synopsis           [Allocate memory for a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
cover_t *
pla_CoverAlloc(int num_inputs)
{
cover_t *ret;

  ret = (cover_t *) ALLOC(cover_t, 1);
  ret->nalloc = 0;
  ret->ncubes = 0;
  ret->cubes = NIL(unsigned int *);
  ret->flags = NIL(short);
  ret->cost_in = 0;
  ret->cost_out = 0;
  ret->p0 = (int *) ALLOC(int, num_inputs);
  ret->p1 = (int *) ALLOC(int, num_inputs);
  
  return ret;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_CubeDel(
  cover_t *cv,
  int pos,
  int cube_size)
{
int i;
   
  if(cv->ncubes == 1) {
    FREE(cv->cubes[0]);
    cv->ncubes = 0;
  }
  else {
    if(pos != cv->ncubes - 1) {
      COPY_CUBE(cv->cubes[pos], cv->cubes[cv->ncubes - 1], cube_size, i);
      cv->flags[pos] = cv->flags[cv->ncubes - 1];
    }
    FREE(cv->cubes[cv->ncubes - 1]);
    cv->ncubes--;
  }
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_VarDel(
  pla_t *pla,
  int var)
{
int i;

  FREE(pla->labels[var]);
  if(var != pla->num_inputs - 1) {
    for(i = 0;i < pla->F->ncubes;i++) {
      switch(VALUE_VAR(pla->F->cubes[i], pla->num_inputs - 1)) {
        case 1: POS_VAR(pla->F->cubes[i], var); break;
        case 2: NEG_VAR(pla->F->cubes[i], var); break;
        case 3: FULL_VAR(pla->F->cubes[i], var); break;
        default:
      }
    }
    for(i = 0;i < pla->D->ncubes;i++) {
      switch(VALUE_VAR(pla->D->cubes[i], pla->num_inputs - 1)) {
        case 1: POS_VAR(pla->D->cubes[i], var); break;
        case 2: NEG_VAR(pla->D->cubes[i], var); break;
        case 3: FULL_VAR(pla->D->cubes[i], var); break;
        default:
      }
    }
    for(i = 0;i < pla->R->ncubes;i++) {    
      switch(VALUE_VAR(pla->R->cubes[i], pla->num_inputs - 1)) {
        case 1: POS_VAR(pla->R->cubes[i], var); break;
        case 2: NEG_VAR(pla->R->cubes[i], var); break;
        case 3: FULL_VAR(pla->R->cubes[i], var); break;
        default:
      }
    }
    pla->F->p0[var] = pla->F->p0[pla->num_inputs - 1];
    pla->F->p1[var] = pla->F->p1[pla->num_inputs - 1];
    pla->D->p0[var] = pla->D->p0[pla->num_inputs - 1];
    pla->D->p1[var] = pla->D->p1[pla->num_inputs - 1];
    pla->R->p0[var] = pla->R->p0[pla->num_inputs - 1];
    pla->R->p1[var] = pla->R->p1[pla->num_inputs - 1];                            
    pla->labels[var] = pla->labels[pla->num_inputs - 1];
  }
  pla->num_inputs--;
  pla->labels = (char **) REALLOC(char *, pla->labels, pla->num_inputs);

  i = (2*pla->num_inputs + pla->num_outputs)/BPI + 1;
  if(i < pla->cube_size) {
    pla->cube_size = i;
    for(i = 0;i < pla->F->ncubes;i++) {
      pla->F->cubes[i] = (int *) REALLOC(int, pla->F->cubes[i],
        pla->cube_size);
    }
    for(i = 0;i < pla->D->ncubes;i++) {
      pla->D->cubes[i] = (int *) REALLOC(int, pla->D->cubes[i],
        pla->cube_size);
    }
    for(i = 0;i < pla->R->ncubes;i++) {
      pla->R->cubes[i] = (int *) REALLOC(int, pla->R->cubes[i],
        pla->cube_size);
    }
  }
  pla->F->p0 = (int *) REALLOC(int, pla->F->p0, pla->num_inputs);
  pla->F->p1 = (int *) REALLOC(int, pla->F->p1, pla->num_inputs);          	  
  pla->D->p0 = (int *) REALLOC(int, pla->D->p0, pla->num_inputs);
  pla->D->p1 = (int *) REALLOC(int, pla->D->p1, pla->num_inputs);          	  
  pla->R->p0 = (int *) REALLOC(int, pla->R->p0, pla->num_inputs);
  pla->R->p1 = (int *) REALLOC(int, pla->R->p1, pla->num_inputs);
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
pla_CubePrint(
  FILE *fp,
  unsigned int *c,
  int num_inputs,
  int num_outputs)
{
int j, v;

  for(j = 0;j < num_inputs;j++) {
    v = VALUE_VAR(c, j);
    switch(v) {
      case 2 : fprintf(fp,"0"); break;
      case 1 : fprintf(fp,"1"); break;
      case 3 : fprintf(fp,"-"); break;
      case 0 : fprintf(fp,"x"); break;
      default: fprintf(fp,"?");
    }
  }
  fprintf(fp," ");
  for(j = 0;j < num_outputs;j++) {
    if(VALUE_OUT(c, j)) {
      fprintf(fp,"1");
    }
    else {
      fprintf(fp,"0");
    }
  }
  fprintf(fp,"\n");
}

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
plaCoverPrint(
  FILE *fp,
  cover_t *cv,
  int num_inputs,
  int num_outputs)
{
int i;

  for(i = 0;i < cv->ncubes;i++) {  
    pla_CubePrint(fp, cv->cubes[i], num_inputs, num_outputs);
  }
}

/**Function********************************************************************

  Synopsis           [Find a cube]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
int
plaCubeFind(
  cover_t *cv,
  unsigned int *cube,
  int cube_size)
{
int ret, i, j;
  
  for(i = 0;i < cv->ncubes;i++) {
    ret = i;
    for(j = 0;j < cube_size;j++) {
      if(cv->cubes[i][j] != cube[j]) {
        ret = -1;
        break;
      }
    }
    if(ret != -1) return ret;
  }
  
  return -1;
}

/**Function********************************************************************

  Synopsis           [Compute the cost of the F pla from scratch]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaCost(
  cover_t *cv,
  int num_inputs,
  int num_outputs)
{
int i,j,t;

  cv->cost_in = cv->cost_out = 0; 
  for(i = 0;i < cv->ncubes;i++) {
    for(j = 0;j < num_inputs;j++) {
      t = VALUE_VAR(cv->cubes[i], j);
      if((t == 2) || (t == 1)) {
        cv->cost_in++;
      }
    }
    for(j = 0;j < num_outputs;j++) {
      if(VALUE_OUT(cv->cubes[i], j)) {
        cv->cost_out++;
      }
    }
  }
}

/**Function********************************************************************

  Synopsis           [Unwarp a pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaUnwarpOutput(
  cover_t *cv,
  int cube_size,
  int num_outputs)
{
unsigned int *c;
int i, j, k;

  if(num_outputs < 2) return;
  
  for(i = 0;i < cv->ncubes;i++) {
    for(j = 0;j < num_outputs;j++) {
      if(VALUE_OUT(cv->cubes[i], j)) break;
    }
    j++;
    for(;j < num_outputs;j++) {
      if(VALUE_OUT(cv->cubes[i], j)) {
        NEG_OUT(cv->cubes[i], j);
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, cv->cubes[i], cube_size, k);
        for(k = 0;k < num_outputs;k++) {
          if(k == j) {
            POS_OUT(c, k);
          }
          else {
            NEG_OUT(c, k);
          }
        }
        pla_CubeAdd(cv, c);
      }
    }  
  }
}

/**Function********************************************************************

  Synopsis           [Return the distance between two cubes]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
int
plaCubeDist(
  unsigned int *c1,
  unsigned int *c2,
  int ni,
  int no)
{
int i, t1, t2, dist = 0;

  for(i = 0;i < ni;i++) {
    t1 = VALUE_VAR(c1, i);
    t2 = VALUE_VAR(c2, i);
    if(((t1 == 2) && (t2 == 1)) ||
       ((t1 == 1) && (t2 == 2))) {
      dist++;
    }  
  }
  for(i = 0;i < no;i++) {
    if(VALUE_OUT(c1, i) && VALUE_OUT(c2, i)) {
      return dist;
    }
  }
  dist++;

  return dist;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
int
plaCubeContain(
  unsigned int *c1,
  unsigned int *c2,
  int ni,
  int no)
{
int i, t1, t2;

  for(i = 0;i < ni;i++) {
    t1 = VALUE_VAR(c1, i);
    t2 = VALUE_VAR(c2, i);
    if(((t1 == 2) && (t2 != 2)) ||
       ((t1 == 1) && (t2 != 1))) {
      return 0;
    }  
  }
  for(i = 0;i < no;i++) {
    t1 = VALUE_OUT(c1, i);
    t2 = VALUE_OUT(c2, i);
    if(VALUE_OUT(c1, i) && VALUE_OUT(c2, i)) {
      return 1;
    }
  }
  return 0;
}

/**Function********************************************************************

  Synopsis           [Sort a pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaOrder(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
  ni = num_inputs;
  no = num_outputs;
  cs = cube_size;

  qsort(cv->cubes, cv->ncubes, sizeof(unsigned int *), cube_compare);
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
int
plaBinateSelect(
  cover_t *cv,
  int num_inputs)
{
int i,j,t, ret = -1, maxsum = -1;

  for(i = 0;i < num_inputs;cv->p0[i] = cv->p1[i] = 0,i++);

  for(i = 0;i < cv->ncubes;i++) {
    for(j = 0;j < num_inputs;j++) {
      t = VALUE_VAR(cv->cubes[i], j);
      if(t == 1) {
        cv->p1[j]++;
      }
      else
        if(t == 2) {
          cv->p0[j]++;
        }
      if(((cv->p1[j] + cv->p0[j]) > maxsum) && (cv->p0[j] > 0) && 
          (cv->p1[j] > 0)) {
        maxsum = cv->p1[j] + cv->p0[j];
        ret = j;
      }
    }
  }  

  return ret;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
unsigned int *
plaCubeCofactor(
  unsigned int *c1,
  unsigned int *c2,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
unsigned int *ret;
int t1, t2, i;

  ret = (unsigned int *) ALLOC(unsigned int, cube_size);
  for(i = 0;i < num_inputs;i++) {
    t1 = VALUE_VAR(c1, i);
    t2 = VALUE_VAR(c2, i);
    if(((t1 == 2) && (t2 == 1)) ||
       ((t1 == 1) && (t2 == 2))) {
      FREE(ret);
      return NIL(unsigned int);
    }
    else
      if((t2 == 2) || (t2 == 1)) {
        FULL_VAR(ret, i);
      }
      else {
        switch(t1) {
          case 1: POS_VAR(ret, i); break;
          case 2: NEG_VAR(ret, i); break;
          case 3: FULL_VAR(ret, i); break;
          case 0: {
            FREE(ret);
            return NIL(unsigned int);
          } break;
        }
      }
  }  
  for(i = 0;i < num_outputs;i++) {
    if(!VALUE_OUT(c2,i)) {
      POS_OUT(ret, i);
    }  
    else {
      if(!VALUE_OUT(c1,i)) {
        NEG_OUT(ret, i);
      }
      else {
        POS_OUT(ret, i);
      }
    }
  }

  return ret;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
unsigned int *
plaCubeConsensus(
  unsigned int *c1,
  unsigned int *c2,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int i,dist = 0;
unsigned int *ret;

  ret = (unsigned int *) ALLOC(unsigned int, cube_size);  
  AND_CUBE(ret, c1, c2, cube_size, i);

  for(i = 0;i < num_inputs;i++) {
    if(VALUE_VAR(ret, i) == 0) {
      if(dist == 1) {
        FREE(ret);
        return NIL(unsigned int);      
      }
      dist++;
      FULL_VAR(ret, i);
    }
  }
  for(i = 0;i < num_outputs;i++) {
    if(VALUE_OUT(c1,i) && VALUE_OUT(c2,i)) {
      return ret;
    }
  }
  dist++;

  if(dist >= 2) {
    FREE(ret);
    return NIL(unsigned int);
  }
  
  for(i = 0;i < num_outputs;i++) {
    if(VALUE_OUT(c1,i) || VALUE_OUT(c2,i)) {
      POS_OUT(ret, i);    
    }
    else {
      NEG_OUT(ret, i);
    }  
  }

  return ret;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
cover_t *
plaCoverConsensus(
  cover_t *cv,
  unsigned int *cube,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t *ret;
unsigned int *c;
int i;

  ret = pla_CoverAlloc(num_inputs);

  for(i = 0;i < cv->ncubes;i++) {
    c = plaCubeConsensus(cv->cubes[i], cube, num_inputs, num_outputs,
      cube_size);
    if(c != NIL(unsigned int)) {
      pla_CubeAdd(ret, c);
    }
  }

  return ret;
}  

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaCoverVarCofactor(
  cover_t *cv,
  int var,
  cover_t **c0,
  cover_t **c1,
  int num_inputs,
  int cube_size)
{
int i, count;
unsigned int *c;

  if(var > num_inputs) {
    fprintf(stderr,"plaCoverVarCofactor: not a valid input variable!\n");
    return;
  }  
  *c0 = pla_CoverAlloc(num_inputs);
  *c1 = pla_CoverAlloc(num_inputs);  

  for(i = 0;i < cv->ncubes;i++) {
    count = VALUE_VAR(cv->cubes[i], var);
    switch(count) {
      case 1: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, cv->cubes[i], cube_size, count);
        FULL_VAR(c, var);
        pla_CubeAdd(*c1, c);
      } break;

      case 2: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, cv->cubes[i], cube_size, count);
        FULL_VAR(c, var);
        pla_CubeAdd(*c0, c);
      } break;

      case 3: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, cv->cubes[i], cube_size, count);
        FULL_VAR(c, var);
        pla_CubeAdd(*c1, c);

        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, cv->cubes[i], cube_size, count);
        FULL_VAR(c, var);
        pla_CubeAdd(*c0, c);      
      } break;
    }  
  }
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
void
plaCoverVarCofactorV2(
  cover_t **cv,
  int var,
  cover_t **c1,
  int num_inputs,
  int cube_size)
{
int i = 0, count;
unsigned int *c;

  *c1 = pla_CoverAlloc(num_inputs);  

  while(i < (*cv)->ncubes) {
    switch(VALUE_VAR((*cv)->cubes[i], var)) {
      case 1: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, (*cv)->cubes[i], cube_size, count);
        FULL_VAR(c, var);
        pla_CubeAdd(*c1, c);
        pla_CubeDel(*cv, i, cube_size);
      } break;

      case 2: i++; break;

      case 3: {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, (*cv)->cubes[i], cube_size, count);
        FULL_VAR(c, var);
        pla_CubeAdd(*c1, c);
        i++;
      } break;     
    }  
  }
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
cover_t *
plaCoverCubeCofactor(
  cover_t *cv,
  unsigned int *cube,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
cover_t *ret;
unsigned int *c;
int i;

  ret = pla_CoverAlloc(num_inputs);

  for(i = 0;i < cv->ncubes;i++) {
    c = plaCubeCofactor(cv->cubes[i], cube, num_inputs, num_outputs,
      cube_size);
    if(c != NIL(unsigned int)) pla_CubeAdd(ret, c);    
  }

  return ret;
}  

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
cover_t *
plaMwc(
  cover_t *c0,
  cover_t *c1,
  int var,
  int num_inputs,
  int num_outputs,
  int cube_size,
  short contain)
{
cover_t *ret;
int i, j, k, ins;
unsigned int *c;

  ret = pla_CoverAlloc(num_inputs);

  for(i = 0;i < c0->ncubes;i++) {
    j = 0;
    ins = 0;
    while(j < c1->ncubes) {
      for(k = 0;k < cube_size;k++) {
        if(c0->cubes[i][k] != c1->cubes[j][k]) break;
      }
      ins = 0;
      if(k == cube_size) {
        c = (unsigned int *) ALLOC(unsigned int, cube_size);
        COPY_CUBE(c, c0->cubes[i], cube_size, k);
        pla_CubeAdd(ret, c);
        ret->flags[ret->ncubes - 1] = c0->flags[i];
        pla_CubeDel(c1, j, cube_size);
        ins = 1;
        break;
      }
      else
        if(contain == 1) {
          if(plaCubeContain(c0->cubes[i], c1->cubes[j], num_inputs,
            num_outputs)) {
            c = (unsigned int *) ALLOC(unsigned int, cube_size);
            COPY_CUBE(c, c1->cubes[j], cube_size, k);
            pla_CubeAdd(ret, c);
            ret->flags[ret->ncubes - 1] = c1->flags[j];
            pla_CubeDel(c1, j, cube_size);
          }
          else
            if(plaCubeContain(c1->cubes[j], c0->cubes[i], num_inputs,
              num_outputs)) {
              c = (unsigned int *) ALLOC(unsigned int, cube_size);
              COPY_CUBE(c, c0->cubes[i], cube_size, k);
              pla_CubeAdd(ret, c);
              ret->flags[ret->ncubes - 1] = c0->flags[i];
              ins = 1;
              break;
            }
            else {
              j++;
            }
        }
        else {
          j++;
        }
    }
    if(ins == 0) {
      c = (unsigned int *) ALLOC(unsigned int, cube_size);
      COPY_CUBE(c, c0->cubes[i], cube_size, k);
      NEG_VAR(c, var);
      pla_CubeAdd(ret, c);
    }
  }
  for(i = 0;i < c1->ncubes;i++) {
    c = (unsigned int *) ALLOC(unsigned int, cube_size);
    COPY_CUBE(c, c1->cubes[i], cube_size, k);
    POS_VAR(c, var);
    pla_CubeAdd(ret, c);  
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
plaSimplify(
  cover_t *cv,
  int num_inputs,
  int num_outputs,
  int cube_size)
{
int var;
cover_t *c0, *c1, *ret;

  var = plaBinateSelect(cv, num_inputs);
  if(var == -1) {
    pla_Scc(cv, num_inputs, num_outputs, cube_size);
    return cv;
  }
  plaCoverVarCofactor(cv, var, &c0, &c1, num_inputs, cube_size);

  c0 = plaSimplify(c0, num_inputs, num_outputs, cube_size);
  c1 = plaSimplify(c1, num_inputs, num_outputs, cube_size);

  ret = plaMwc(c0, c1, var, num_inputs, num_outputs, cube_size, 1);

  if(ret->ncubes > cv->ncubes) {
    pla_CoverFree(ret);
    ret = cv;
  }
  else {
    pla_CoverFree(cv);
  }  
  pla_CoverFree(c0);
  pla_CoverFree(c1);  
  
  return ret;
}

/**Function********************************************************************

  Synopsis           [Print the pla]

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
cover_t *
plaCoverCopy(
  cover_t *cv,
  int num_inputs,
  int cube_size)
{
cover_t *ret;
int i,j;
unsigned int *c;

  ret = pla_CoverAlloc(num_inputs);
  for(i = 0;i < cv->ncubes;i++) {
    c = (unsigned int *) ALLOC(unsigned int, cube_size);
    COPY_CUBE(c, cv->cubes[i], cube_size, j);
    pla_CubeAdd(ret, c);  
    ret->flags[ret->ncubes - 1] = cv->flags[i];
  }

  return ret;
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
static int
cube_compare(
  const void *c1,
  const void *c2)
{
unsigned int **cube1;
unsigned int **cube2;
int i;

  cube1 = (unsigned int **) c1;
  cube2 = (unsigned int **) c2;

  for(i = 0;i < cs;i++) {
    if((*cube1)[i] != (*cube2)[i]) break;
  }
  if(i == cs) return 0;
   
  if(plaCubeContain(*cube1, *cube2, ni, no)) {
    return -1;
  }
  return 1;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void
StartTimer()
{
  gettimeofday(&time_start, (struct timezone *) 0);    
  getrusage(RUSAGE_SELF, &ru_start);
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void
StopTimer()
{
  getrusage(RUSAGE_SELF, &ru_stop);
  gettimeofday(&time_stop, (struct timezone *) 0);
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static double
GetRtime()
{
  start = ((double) time_start.tv_sec) * 1000000.0 +
    time_start.tv_usec;
  stop = ((double) time_stop.tv_sec) * 1000000.0 +
    time_stop.tv_usec;
  seconds = (stop - start) / 1000000.0;

  return seconds;
}
