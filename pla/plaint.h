/**CHeaderFile*****************************************************************

  FileName    [plaInt.h]

  PackageName [pla]

  Synopsis    [My implementation of espresso]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

  Revision    []

******************************************************************************/

#ifndef _PLAINTH
#define _PLAINTH

#include "pla.h"

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

int GLOBAL_NI;

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */
/*---------------------------------------------------------------------------*/

EXTERN void plaCoverPrint(FILE *fp, cover_t *cv, int num_inputs, int num_outputs);
EXTERN int plaCubeFind(cover_t *cv, unsigned int *cube, int cube_size);
EXTERN void plaCost(cover_t *cv, int num_inputs, int num_outputs);
EXTERN void plaUnwarpOutput(cover_t *cv, int cube_size, int num_outputs);
EXTERN int plaCubeDist(unsigned int *c1, unsigned int *c2, int ni, int no);
EXTERN int plaCubeContain(unsigned int *c1, unsigned int *c2, int ni, int no);
EXTERN void plaOrder(cover_t *cv, int num_inputs, int num_outputs, int cube_size);
EXTERN int plaBinateSelect(cover_t *cv, int num_inputs);
EXTERN unsigned int * plaCubeCofactor(unsigned int *c1, unsigned int *c2, int num_inputs, int num_outputs, int cube_size);
EXTERN unsigned int * plaCubeConsensus(unsigned int *c1, unsigned int *c2, int num_inputs, int num_outputs, int cube_size);
EXTERN cover_t * plaCoverConsensus(cover_t *cv, unsigned int *cube, int num_inputs, int num_outputs, int cube_size);
EXTERN void plaCoverVarCofactor(cover_t *cv, int var, cover_t **c0, cover_t **c1, int num_inputs, int cube_size);
EXTERN void plaCoverVarCofactorV2(cover_t **cv, int var, cover_t **c1, int num_inputs, int cube_size);
EXTERN cover_t * plaCoverCubeCofactor(cover_t *cv, unsigned int *cube, int num_inputs, int num_outputs, int cube_size);
EXTERN cover_t * plaMwc(cover_t *c0, cover_t *c1, int var, int num_inputs, int num_outputs, int cube_size, short contain);
EXTERN cover_t * plaSimplify(cover_t *cv, int num_inputs, int num_outputs, int cube_size);
EXTERN cover_t * plaCoverCopy(cover_t *cv, int num_inputs, int cube_size);
EXTERN cover_t * UnateComplement(cover_t *cv, int num_inputs, int cube_size);
EXTERN void plaEssencial(cover_t **F, cover_t **D, int *moved, int num_inputs, int num_outputs, int cube_size);
EXTERN void plaExpand(cover_t **F, cover_t *R, int num_inputs, int num_outputs, int cube_size, short type, short gasp);
EXTERN void plaGasp(cover_t **F, cover_t *D, cover_t *R, int num_inputs, int num_outputs, int cube_size);
EXTERN void plaIrredundant(cover_t **F, cover_t *D, int num_inputs, int num_outputs, int cube_size);
EXTERN pla_t * plaRead(FILE *fp);
EXTERN void plaReduce(cover_t **F, cover_t *D, int num_inputs, int num_outputs, int cube_size);
EXTERN unsigned int * plaSCCC(cover_t *cv, int num_inputs, int num_outputs, int cube_size);
EXTERN short plaTautology(cover_t *cv, int num_inputs, int num_outputs, int cube_size);

/**AutomaticEnd***************************************************************/

#endif /* _ */
