/*
 * File: fileio.c
 * --------------
 * Contains helper functions for opening and reading data from a file.
 *
 * $Id: fileio.c,v 1.1 2002-11-03 00:19:53 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileio.h"

#define BUFFERLENGTH 256

/*
 * Function: getfile
 * Usage: infile = getfile();
 * --------------------------
 * Prompts the user for a valid filename.
 * 
 */
FILE *getfile(void)
{
  char file[BUFFERLENGTH];
  FILE *ifile;

  printf("Please input a file: ");
  scanf("%s",file);
  while(!(ifile = fopen(file,"r"))) {
    printf("Error in opening %s.\n", file);
    printf("Please input another file: ");
    scanf("%s",file);
  }
  return ifile;
}

/*
 * Function: getelement()
 * Usage: n=getelement(ifile);
 * ---------------------------
 * Returns the first numbers in the line up until the first space.
 */
double getelement(FILE *ifile)
{
  char istr[BUFFERLENGTH], ostr[BUFFERLENGTH];
  int i=0;
  getline(ifile,istr,"");
  getchunk(istr,ostr);
  return strtod(ostr,(char **)NULL);
}

/*
 * Function: getchunk()
 * Usage: getchunk(istr,ostr);
 * --------------------------------
 * Gets the first characters in the string istr up until the first
 * space or until the end of line and puts them into ostr.
 */
void getchunk(char *istr, char *ostr)
{
  int i=0;
  while(!isspace(istr[i]) && istr[i] != EOF && istr[i] != '\n' && istr[i]!='\0') {
    ostr[i]=istr[i];
    i++;
  }
  ostr[i]='\0';
}

/*
 * Function: getline()
 * Usage: c = getline(file,str,omit);
 * ----------------------------------
 * Gets the first line of a file that does not contain
 * any one of the characters in the omit string.  Returns 
 * the last character found.
 * 
 */    
char getline(FILE *file, char *line, char *omit)
{
  char c, *start;
  int i, found=0, N=strlen(omit);

  start = line;
  while((c=fgetc(file))!='\n') {
    if(c==EOF) break;
    *(line++)=c;
    for(i=0;i<N;i++)
      if(c==omit[i]) { found=1; break; }
  }
  if(found) {
    getline(file,start,omit);
  }
  *line='\0';
  return c;
}

/*
 * Function: getfield();
 * Usage: x = getfield(file);
 * --------------------------
 * Returns the next floating point number in a file.
 * 
 */
double getfield(FILE *file, char *str)
{
  int i, k, flag;
  char c;
  double field;

  i=0;
  flag = 1;
  c = fgetc(file);

  if(c==EOF) {
    return;
  }

  while(c != ' ' & c != '\n' & c != EOF) {
    str[i++] = c;
    c = fgetc(file);
  }
  if(i==0) return getfield(file,str);

  str[i]='\0';
  for(k=0;k<i;k++) {
    if(!isdigit(str[k]) & str[k]!='.' & str[k]!='-' & str[k]!='+' & str[k]!='e') {
      flag = 0;
      break;
    }
  }
  
  if(flag)
    return strtod(str,(char **)NULL);
  else
    return getfield(file,str);
}

/*
 * Function: getsize();
 * Usage: N = getsize(filename);
 * -----------------------------
 * Returns the number of rows in a file.
 *
 */
int getsize(char *filename)
{
  int N;
  char c;
  FILE *infile = fopen(filename,"r");
  if(!infile)
    printf("Error opening %s\n",filename);

  N=0;
  while((c=fgetc(infile))!=EOF)
    if(c=='\n') N++;
  
  return N;
}

/*
 * Function: GetValue
 * Usage: Nkmax = (int)GetValue("datafile","Nkmax");
 * --------------------------------------------
 * Returns the value of the specified variable defined in datafile.
 *
 */
double GetValue(char *filename, char *str, int *status)
{
  int ispace;
  char c, istr[BUFFERLENGTH], ostr[BUFFERLENGTH];
  FILE *ifile = fopen(filename,"r");
  *status = 0;

  while(1) {
    getline(ifile,istr,"");
    if(strlen(istr)==0)
      break;
    getchunk(istr,ostr);
    if(!strcmp(ostr,str)) {
      for(ispace=strlen(ostr);ispace<strlen(istr);ispace++) 
	if(!isspace(istr[ispace]))
	  break;
      if(ispace==strlen(istr)-1)
	*status=0;
      else
	*status=1;
      getchunk(&(istr[ispace]),ostr);
      break;
    }
  }
  fclose(ifile);
  
  if(*status) 
    return strtod(ostr,(char **)NULL);
  else
    return 0;
}

/*
 * Function: GetString
 * Usage: GetString(string,"file.dat","ufile",&status);
 * ----------------------------------------------------
 * Obtains the string associated with the key "ufile" from the
 * file "file.dat".
 *
 */
void GetString(char *string, char *filename, char *str, int *status)
{
  int ispace;
  char c, istr[BUFFERLENGTH], ostr[BUFFERLENGTH];
  FILE *ifile = fopen(filename,"r");
  *status = 0;

  while(1) {
    getline(ifile,istr,"");
    if(strlen(istr)==0)
      break;
    getchunk(istr,ostr);
    if(!strcmp(ostr,str)) {
      for(ispace=strlen(ostr);ispace<strlen(istr);ispace++) 
	if(!isspace(istr[ispace]))
	  break;
      if(ispace==strlen(istr)-1)
	*status=0;
      else
	*status=1;
      getchunk(&(istr[ispace]),string);
      break;
    }
  }
  fclose(ifile);
}
