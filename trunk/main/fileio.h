/*
 * File: fileio.h
 * --------------
 * Header file for fileio.c.
 *
 * $Id: fileio.h,v 1.1 2002-11-03 00:20:27 fringer Exp $
 * $Log: not supported by cvs2svn $
 *
 */
#ifndef _fileio_h
#define _fileio_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Function: getfile
 * Usage: infile = getfile();
 * --------------------------
 * Prompts the user for a valid filename.
 * 
 */
FILE *getfile(void);

/*
 * Function: getline()
 * Usage: c = getline(file,str,omit);
 * ----------------------------------
 * Gets the first line of a file that does not contain
 * any one of the characters in the omit string.  Returns 
 * the last character found.
 * 
 */    
char getline(FILE *file, char *line, char *omit);

/*
 * Function: getchunk()
 * Usage: getchunk(istri,ostr);
 * ----------------------------
 * Gets the first characters in the string istr up until the first
 * space and puts them into ostr.
 *
 */
void getchunk(char *istr, char *ostr);

/*
 * Function: getelement()
 * Usage: n=getelement(ifile);
 * ---------------------------
 * Returns the first numbers in the line up until the first space.
 */
double getelement(FILE *ifile);

/*
 * Function: getfield();
 * Usage: x = getfield(file);
 * --------------------------
 * Returns the next floating point number in a file.
 * 
 */
double getfield(FILE *file, char *str);

/*
 * Function: getsize();
 * Usage: N = getsize(filename);
 * -----------------------------
 * Returns the number of rows in a file.
 *
 */
int getsize(char *filename);

/*
 * Function: GetValue
 * Usage: Nkmax = (int)GetValue("datafile","Nkmax");
 * --------------------------------------------
 * Returns the value of the specified variable defined in datafile.
 *
 */
double GetValue(char *filename, char *str, int *status);

/*
 * Function: GetString
 * Usage: GetString(string,"file.dat","ufile",&status);
 * ----------------------------------------------------
 * Obtains the string associated with the key "ufile" from the
 * file "file.dat".
 *
 */
void GetString(char *string, char *filename, char *str, int *status);

#endif
