/**CHeaderFile*****************************************************************

  FileName    [pla.h]

  PackageName [pla]

  Synopsis    [My implementation of espresso]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

  Revision    []

******************************************************************************/

#ifndef _PLAH
#define _PLAH

/*---------------------------------------------------------------------------*/
/* Constant declarations                                                     */
/*---------------------------------------------------------------------------*/

#define MINCUBE 1
#define MINLIT 2

#define BPI 16

#define LOGBPI 4

#define INVBPI 3
#define INVBPI1 4

#define ONE 0xFFFD
#define ZERO 0xFFFE
#define DC 0xFFFC
#define ZEROC 0x0000
#define FULLC 0xFFFF

#define ON_SET 1
#define DC_SET 2
#define OFF_SET 3

#define EXECUTE(oper, time) \
  StartTimer(); oper; StopTimer(); time += GetRtime();

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/

typedef struct cover cover_t;
typedef struct pla pla_t;

/*---------------------------------------------------------------------------*/
/* Structure declarations                                                    */
/*---------------------------------------------------------------------------*/

struct cover {
  unsigned int **cubes;
  short *flags;
  int nalloc;
  int ncubes;
  int cost_in;
  int cost_out;
  int *p0;
  int *p1;
};

struct pla {
  char **labels;
  int num_inputs;
  int num_outputs;
  int cube_size;
  cover_t *F;
  cover_t *D;
  cover_t *R;
};

/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/

#define WHICH_WORD(element) (((element) >> LOGBPI) + 1)
#define WHICH_BIT(element)  ((element) & (BPI-1))

#define RESET_CUBE(cube, csize, count)\
  for(count = 0;count < csize;cube[count] = ZEROC, count++);

#define FULL_CUBE(cube, csize, count)\
  for(count = 0;count < csize;cube[count] = FULLC, count++)

#define SET_INSERT(cube, pos) \
  cube[(pos) >> INVBPI1] |= (1 << WHICH_BIT(pos))

#define SET_REMOVE(cube, pos) \
  cube[(pos) >> INVBPI1] &= ~(1 << WHICH_BIT(pos))

#define IS_IN_SET(cube, pos) \
  (cube[(pos) >> INVBPI1] & (1 << WHICH_BIT(pos)))

#define POS_VAR(cube, pos) \
   {cube[pos >> INVBPI] |= (1 << WHICH_BIT(pos << 1));\
    cube[pos >> INVBPI] &= ~(2 << WHICH_BIT(pos << 1)); }

#define NEG_VAR(cube, pos)\
   {cube[pos >> INVBPI] |= (2 << WHICH_BIT(pos << 1));\
    cube[pos >> INVBPI] &= ~(1 << WHICH_BIT(pos << 1)); }

#define INV_VAR(cube, pos)\
  cube[pos >> INVBPI] ^= (3 << WHICH_BIT(pos << 1))

#define ZERO_VAR(cube, pos)\
  cube[pos >> INVBPI] &= ~(3 << WHICH_BIT(pos << 1))

#define FULL_VAR(cube, pos)\
  cube[(pos) >> INVBPI] |= (3 << WHICH_BIT((pos) << 1))

#define VALUE_VAR(cube, pos)\
  (cube[(pos) >> INVBPI] & (3 << WHICH_BIT((pos) << 1))) >> WHICH_BIT((pos) << 1)

#define POS_OUT(cube,pos) SET_INSERT(cube, pos + (GLOBAL_NI << 1))

#define NEG_OUT(cube,pos) SET_REMOVE(cube, pos + (GLOBAL_NI << 1))

#define VALUE_OUT(cube, pos) IS_IN_SET(cube, pos + (GLOBAL_NI << 1))
  
#define AND_CUBE(c, a, b, csize, count)\
  for(count = 0;count < csize;(c)[count] = (a)[count] & (b)[count],count++)

#define OR_CUBE(c, a, b, csize, count)\
  for(count = 0;count < csize;(c)[count] = (a)[count] | (b)[count],count++)  

#define INV_CUBE(c, a, csize, count)\
  for(count = 0;count < csize;(c)[count] = !(a)[count],count++)  

#define COPY_CUBE(c, a, csize, count)\
  for(count = 0;count < csize;c[count] = a[count],count++)
    
#define RESET_FLAG(cv, cube) cv->flags[cube] = 0x0000

#define SET_PRIME(cv, cube) cv->flags[cube] |= 0x0001
#define UNSET_PRIME(cv, cube) cv->flags[cube] &= 0xFFFE
#define IS_PRIME(cv, cube) (cv->flags[cube] && 0x0001)

/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */
/*---------------------------------------------------------------------------*/

EXTERN void pla_Espresso(pla_t *pla, short stats);
EXTERN void pla_Print(FILE *fp, pla_t *pla, short type);
EXTERN pla_t * pla_Alloc(void);
EXTERN unsigned int * pla_CubeNew(int cube_size, short type);
EXTERN void pla_CubeAdd(cover_t *cv, unsigned int *cube);
EXTERN void pla_Free(pla_t *pla);
EXTERN void pla_CoverFree(cover_t *cv);
EXTERN void pla_Scc(cover_t *cv, int num_inputs, int num_outputs, int cube_size);
EXTERN cover_t * pla_CoverAlloc(int num_inputs);
EXTERN void pla_CubeDel(cover_t *cv, int pos, int cube_size);
EXTERN void pla_VarDel(pla_t *pla, int var);
EXTERN void pla_CubePrint(FILE *fp, unsigned int *c, int num_inputs, int num_outputs);
EXTERN void pla_Complement(cover_t *F, cover_t *D, cover_t **R, int num_inputs, int num_outputs, int cube_size);

/**AutomaticEnd***************************************************************/

#endif /* _ */
