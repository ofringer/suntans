/*
 * File: fileio.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains helper functions for opening and reading data from a file.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileio.h"
#include <errno.h>
#include "defaults.h"

#define BUFFERLENGTH 256
#define THISFILE "fileio.c"

/*
 * Function: MyFOpen
 * Usage: fid = MyFOpen(string,"r","GetValue");
 * --------------------------------------------
 * Exits if the requested file does not exist.
 * The third string is useful for determining which
 * function the function was called from.
 *
 */
FILE *MyFOpen(char *file, char *perms, char *caller) {
  extern int errno;
  char str[BUFFERLENGTH];
  FILE *fid = fopen(file,perms);

  if(!fid) {
    sprintf(str,"Error in Function %s while trying to open %s",caller,file);
    perror(str);
    exit(EXIT_FAILURE);
  } else
    return fid;
}

  
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
  if(scanf("%s",file) == 0)
    printf("Error reading from file in getfile\n");
  while(!(ifile = fopen(file,"r"))) {
    printf("Error in opening %s.\n", file);
    printf("Please input another file: ");
    if(scanf("%s",file) == 0)
      printf("Error reading from another file in getfile\n");
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
  mygetline(ifile,istr,"");
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
 * Function: mygetline()
 * Usage: c = mygetline(file,str,omit);
 * ----------------------------------
 * Gets the first line of a file that does not contain
 * any one of the characters in the omit string.  Returns 
 * the last character found.
 * 
 */    
char mygetline(FILE *file, char *line, char *omit)
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
    mygetline(file,start,omit);
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
    if(VERBOSE > 2)
      printf("End of file obtained, if  an error follows it is possible geometry is not consistent.\n");
    return;
  }

  while(c != ' ' & c != '\t' & c != '\n' & c != '\r' & c != EOF) {
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
 * Function: getcolumn();
 * Usage: N = getcolumn(filename);
 * -----------------------------
 * Returns the number of rows in a file.
 *
 */
int getcolumn(char *filename)
{
  int N;
  char c;
  FILE *infile = fopen(filename,"r");
  if(!infile)
    printf("Error opening %s\n",filename);

  N=1;
  while((c=fgetc(infile))!=EOF){
    if(c==' ') 
      N++;
    if(c=='\n')
      break;
  }
  fclose(infile);
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
  FILE *ifile = MyFOpen(filename,"r","GetValue");
  *status = 0;

  while(1) {
    mygetline(ifile,istr,"");
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
 * Function: GetDefaultValue
 * Usage: Nkmax = (int)GetDefaultValue("Nkmax",&status);
 * -----------------------------------------------------
 * Returns the default value of the specified variable.  This is required
 * if the requested value is not in the specified datafile.
 * 
 */
double GetDefaultValue(char *str, int *status) {
  *status=1;

  if(!strcmp(str,"prettyplot")) {

    return prettyplot_DEFAULT;

  }
  else if(!strcmp(str,"linearFS")){

    return linearFS_DEFAULT;

  } else if(!strcmp(str,"grav")){

    return grav_DEFAULT;

  } else if(!strcmp(str,"interp")) {

    return interp_DEFAULT;

  } else if(!strcmp(str,"minimum_depth")) {

    return minimum_depth_DEFAULT;

  } else if(!strcmp(str,"fixdzz")) {

    return fixdzz_DEFAULT;

  } else if(!strcmp(str,"TVDsalt")) {

    return TVDsalt_DEFAULT;

  } else if(!strcmp(str,"TVDtemp")) {

    return TVDtemp_DEFAULT;

  } else if(!strcmp(str,"TVDturb")) {

    return TVDturb_DEFAULT;

  } else if(!strcmp(str,"laxWendroff")) {

    return laxWendroff_DEFAULT;

  } else if(!strcmp(str,"laxWendroff_Vertical")) {

    return laxWendroff_Vertical_DEFAULT;

  } else if(!strcmp(str,"hprecond")) {

    return hprecond_DEFAULT;

  } else if(!strcmp(str,"ntoutStore")) {

    return ntoutStore_DEFAULT;

  } else if(!strcmp(str,"AB")) {

    return AB_DEFAULT;

 } else if(!strcmp(str,"TVDmomentum")) {

    return TVDmomentum_DEFAULT;

 } else if(!strcmp(str,"conserveMomentum")) {

    return conserveMomentum_DEFAULT;

 } else if(!strcmp(str,"wetdry")) {

    return wetdry_DEFAULT;

 } else if(!strcmp(str,"thetaM")) {

    return thetaM_DEFAULT;

 } else if(!strcmp(str,"smoothbot")) {

    return smoothbot_DEFAULT;
  
 } else if(!strcmp(str,"latitude")) {
    
   return latitude_DEFAULT;
   
 } else if(!strcmp(str,"metmodel")) {
    
   return metmodel_DEFAULT;
 
 } else if(!strcmp(str,"varmodel")) {
    
   return varmodel_DEFAULT;
   
 } else if(!strcmp(str,"nugget")) {
    
   return nugget_DEFAULT;
 
 } else if(!strcmp(str,"sill")) {
    
   return sill_DEFAULT;
 
 } else if(!strcmp(str,"range")) {
    
   return range_DEFAULT;
   
} else if(!strcmp(str,"outputNetcdf")) {
    
   return outputNetcdf_DEFAULT;

} else if(!strcmp(str,"Lsw")) {
    
   return Lsw_DEFAULT;
 
} else if(!strcmp(str,"Cda")) {
    
   return Cda_DEFAULT;    
 
} else if(!strcmp(str,"Ce")) {
    
   return Ce_DEFAULT; 
  
} else if(!strcmp(str,"Ch")) {
    
   return Ch_DEFAULT;   

} else if(!strcmp(str,"netcdfBdy")) {
    
   return netcdfBdy_DEFAULT;   

} else if(!strcmp(str,"readinitialnc")) {
    
   return readinitialnc_DEFAULT;   

} else if(!strcmp(str,"calcage")) {
    
   return calcage_DEFAULT;   

} else if(!strcmp(str,"calcaverage")) {
    
   return calcaverage_DEFAULT;   

} else if(!strcmp(str,"maxFaces")) {
    
   return maxFaces_DEFAULT;   


}else {
    *status=0;
    return 0;
  }


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
    mygetline(ifile,istr,"");
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
