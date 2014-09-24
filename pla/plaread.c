/**CFile***********************************************************************

  FileName    [plaRead.c]

  PackageName [pla]

  Synopsis    [Read a pla description.]

  Description []

  SeeAlso     []

  Author      [Ivan Jeukens]

  Copyright   []

******************************************************************************/

#include <stdio.h>
#include <ctype.h>
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

#define FALSE 0
#define TRUE 1

/*---------------------------------------------------------------------------*/
/* Type declarations                                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Structure declarations                                                    */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Variable declarations                                                     */
/*---------------------------------------------------------------------------*/

static short line_length_error;
static int lineno;

extern int GLOBAL_NI;

/*---------------------------------------------------------------------------*/
/* Macro declarations                                                        */
/*---------------------------------------------------------------------------*/


/**AutomaticStart*************************************************************/

/*---------------------------------------------------------------------------*/
/* Static function prototypes                                                */
/*---------------------------------------------------------------------------*/

static void read_cube(FILE *fp, pla_t *cv);
static void skip_line(FILE *fpin, FILE *fpout, short echo);
static char * get_word(FILE *fp, char *word);

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
pla_t *
plaRead(FILE *fp)
{
pla_t *pla;
int ch;
char word[256];

  lineno = 1;
  line_length_error = FALSE;
  
  pla = pla_Alloc();

loop:
  switch(ch = getc(fp)) {
    case EOF: return pla;

    case '\n': lineno++; break;
    
    case ' ': 
    case '\t': 
    case '\f': 
    case '\r': break;

    case '#': {
      (void) ungetc(ch, fp);
      skip_line(fp, stdout, FALSE);
    } break;

    case '.': {
      if(strcmp(get_word(fp, word), "i") == 0) {
        if(pla->num_inputs != 0) {
          fprintf(stderr, "extra .i ignored\n");
	  skip_line(fp, stdout, FALSE);
	}
	else {
	  (void) fscanf(fp, "%d", &pla->num_inputs);
	  GLOBAL_NI = pla->num_inputs;
          pla->F->p0 = (int *) ALLOC(int, pla->num_inputs);
          pla->F->p1 = (int *) ALLOC(int, pla->num_inputs);          	  
          pla->D->p0 = (int *) ALLOC(int, pla->num_inputs);
          pla->D->p1 = (int *) ALLOC(int, pla->num_inputs);          	  
          pla->R->p0 = (int *) ALLOC(int, pla->num_inputs);
          pla->R->p1 = (int *) ALLOC(int, pla->num_inputs);          	  
        }
      }
      else
      if(strcmp(word,"o") == 0) {       
        if(pla->num_outputs != 0) {
          fprintf(stderr,"extra .o ignored\n");
        }
        else {
          (void) fscanf(fp, "%d", &pla->num_outputs);
          pla->cube_size = (2*pla->num_inputs + pla->num_outputs)/BPI + 1;
        }
      }
      else
      if(strcmp(word,"mv") == 0) {
        fprintf(stderr,".mv not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"p") == 0) {
        (void) fscanf(fp,"%d", &ch);
      }
      else
      if(strcmp(word,"label") == 0) {
         fprintf(stderr,".label not supported\n");
         skip_line(fp, stdout, FALSE); 
      }
      else
      if(strcmp(word,"ilb") == 0) {
        fprintf(stderr,".ilb not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"ob") == 0) {
        fprintf(stderr,".ob no supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if((strcmp(word,"e") == 0) || (strcmp(word,"end") == 0)) {
        return pla;
      }
      else
      if(strcmp(word,"kiss") == 0) {
        fprintf(stderr,".kiss not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"type") == 0) {
        fprintf(stderr,".type not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"symbolic") == 0) {
        fprintf(stderr,".symbolic not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"symbolic-output") == 0) {
        fprintf(stderr,".symbolic-output not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"phase") == 0) {
        fprintf(stderr,".phase not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else
      if(strcmp(word,"pair") == 0) {
        fprintf(stderr,".pair not supported.\n");
        skip_line(fp, stdout, FALSE);
      }
      else {
        fprintf(stderr,"Unkown command %c%s\n", ch, word);
        skip_line(fp, stdout, FALSE);
      }

    } break;
    
    default: {
      (void) ungetc(ch, fp);
      if((pla->num_inputs == 0) && (pla->num_outputs == 0)) {
        fprintf(stderr,"missing .i and/or .o\n");
        skip_line(fp, stdout, FALSE);
      }
      else {
        read_cube(fp, pla);
      }
    }
  }
  goto loop;
    
  return pla;    
}

/*---------------------------------------------------------------------------*/
/* Definition of static functions                                            */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void
read_cube(
  FILE *fp,
  pla_t *pla)
{
int var;
unsigned int *c, *cdc, *con, *coff;
short saveoff = 0, savedc = 0, saveon = 0;

  c = pla_CubeNew(pla->cube_size, 0);

  for(var = 0;var < pla->num_inputs; var++) {
    switch(getc(fp)) {
      case EOF: goto bad_char;

      case '\n': {
        fprintf(stderr,"product term(s) %s\n",
          "span more than one line (warning only)");
        line_length_error = TRUE;
	lineno++;
	var--;
      } break;
      
      case ' ': 
      case '|': 
      case '\t': var--; break;
      
      case '2':
      case '-': {
        FULL_VAR(c, var);
      } break;
      
      case '0': {
        NEG_VAR(c, var);
      } break;
	
      case '1': {
        POS_VAR(c, var);
      } break;

      case '?': break;

      default: goto bad_char;
    }
  }

  con = pla_CubeNew(pla->cube_size, 0);  
  cdc = pla_CubeNew(pla->cube_size, 0);
  coff = pla_CubeNew(pla->cube_size, 0);
  
  for(var = 0;var < pla->cube_size;var++) {
    con[var] = c[var];
    cdc[var] = c[var];
    coff[var] = c[var];
  }
  if(pla->num_outputs == 0) {
    saveon = 1;
  }
  for(var = 0;var < pla->num_outputs;var++) {
    switch(getc(fp)) {
      case EOF: goto bad_char;

      case '\n': {
        fprintf(stderr, "product term(s) %s\n",
	  "span more than one line (warning only)");
        line_length_error = TRUE;
	lineno++;
	var--;
      } break;    

      case ' ': 
      case '|':
      case '\t': var--; break;

      case '4':
      case '1': {
        POS_OUT(con, var);
        saveon = 1;
      } break;
      
      case '3':
      case '0': {
        POS_OUT(coff, var);
        saveoff = 1;
      } break;

      case '2':
      case '-': {
        POS_OUT(cdc, var);
        savedc = 1;
      } break;
    }
  }
  
  if(saveon == 1) {
    pla_CubeAdd(pla->F, con);
  }
  else {
    FREE(con);
  }
  if(savedc == 1) {
    pla_CubeAdd(pla->D, cdc);
  }
  else {
    FREE(cdc);
  }
  if(saveoff == 1) {
    pla_CubeAdd(pla->R, coff);
  }
  else {
    FREE(coff);
  }
  FREE(c);
  
  return;
  
bad_char:
  fprintf(stderr, "(warning): input line #%d ignored\n", lineno);
  skip_line(fp, stdout, TRUE);
  
  return;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static void 
skip_line(
  FILE *fpin,
  FILE *fpout,
  short echo)
{
int ch;

  while((ch = getc(fpin)) != EOF && ch != '\n')
    if(echo) putc(ch, fpout);

  if(echo) putc('\n', fpout);
  lineno++;
}

/**Function********************************************************************

  Synopsis           []

  Description        []

  SideEffects        []

  SeeAlso            []

******************************************************************************/
static char *
get_word(
  FILE *fp,
  char *word)
{
int ch, i = 0;

  while((ch = getc(fp)) != EOF && isspace(ch)) ;

  word[i++] = ch;
  while((ch = getc(fp)) != EOF && ! isspace(ch)) word[i++] = ch;

  word[i++] = '\0';
  return word;
}
