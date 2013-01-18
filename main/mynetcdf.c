/*
* NetCDF IO functions
* -------------------
* All functions either generic or specific involved in netcdf io should go here
* 
*/ 

#include "mynetcdf.h"


/*
* Function: nc_read_3D()
* ----------------------
* Reads a 3D array from a netcdf file and returns the output in an array (not a vector)
*
* Warning: there are no dimension checks performed here so be careful.
* The size of the array dimension should equal 'count' ie [n, k ,j]
*/

void nc_read_3D(int ncid, char *vname, size_t *start, size_t *count, REAL ***tmparray){

    int j, k, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] * (int)count[2] ];
    REAL outdata[(int)count[0]][(int)count[1]][(int)count[2]];

    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &outdata[0][0][0]))) 
	ERR(retval); 

    // Loop through and insert the vector values into an array
    for(n=0;n<(int)count[0];n++){
	for(k=0;k<(int)count[1];k++){
	    for(j=0;j<(int)count[2];j++){
		//Linear index
		//ii = n*(int)count[1]*(int)count[2]+k*(int)count[2]+j;
		//tmparray[n][k][j]=tmpvec[ii];
		tmparray[n][k][j]=outdata[n][k][j];
	    }
	}
    }

}// End function


/*
* Function: nc_read_2D()
* ----------------------
* Reads a 2D array from a netcdf file and returns the output in an array (not a vector)
*
* Warning: there are no dimension checks performed here so be careful.
* The size of the array dimension should equal 'count' ie [n,,j]
*/

void nc_read_2D(int ncid, char *vname, size_t *start, size_t *count, REAL **tmparray, int myproc){

    int j, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] ];
    REAL outdata[(int)count[0]][(int)count[1]];

    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &outdata[0][0]))) 
	ERR(retval); 

    // Loop through and insert the vector values into an array
   
    for(n=0;n<(int)count[0];n++){
	for(j=0;j<(int)count[1];j++){
	    //Linear index
	    //ii = n*(int)count[1]+j;
	    //tmparray[n][j]=tmpvec[ii];
	    tmparray[n][j]=outdata[n][j];
	    //printf("myproc: %d, start[0]: %d, n: %d of %d, j: %d of %d, outdata[n][j]: %f, tmparray[n][j]: %f\n",myproc, (int)start[0], n,(int)count[0],j,(int)count[1],outdata[n][j], tmparray[n][j]);
	}
    }
  
}// End function
