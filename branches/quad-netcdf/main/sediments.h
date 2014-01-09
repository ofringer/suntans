/*
 * File: sediments.h
 * Author : Yun Zhang
 * Institution: Stanford University
 * --------------------------------
 * Header file for sediment.c
 * including all the variables for sediments
 *
 */
#ifndef _sediments_h
#define _sediments_h

#include "suntans.h"
#include "grid.h"
#include "phys.h"

typedef struct _sedimentsT {
REAL ***SediC,// sediment concentration [fraction][cell][Nkmax]
     ***SediCbed, // sediment concentration in bed [fraction][cell][Nlayer]
     ***boundary_sediC, // sediment transport boundary condition [fraction][cell][Nkmax]
     ***Ws, // settling velocity [fraction][cell][Nkmax]
     *Ws0, // constant settling velocity [fraction] -> given in sedi.dat
     *Ds,  // sediment diameter [fraction] -> given in sedi.dat
     *Gsedi, // density of sediment/ density of water [fraction] -> given in sedi.dat
     *Seditb, // store the tb for each cell [cell]
     *Seditbmax, // store the tb max
     //***Erosion, // erosion term in first layer is for suspended sediment transport [fraction][cell][Nlayer]
     //***Erosion_old, //added
     //**Woldsedi, //vertical velocity for sediment particles [cell][Nkmax+1]
     **Wnewsedi, // calculated by SedimentVerticalVelocity
     **SediKappa_tv, // tubulent sediment diffusivity [cell][Nkmax]
     **Erosiontotal, // total erosion for each cell each layer [cell][Nlayer]
     **Erosiontotal_old, // store former step [cell][Nlayer]
     //**Neterosion,  // store the net erosion in Bedinterval steps [cell][Nlayer]
     **Deposition, // deposition term in first layer is for suspended sediment transport [fraction][cell]
     //**Deposition_old, //added
     *Consolid, // consolidation rate for each layer [Nlayer]-> given in sedi.dat
     //*Bedmudratio, // the ratio of mud for each layer [Nlayer]-> given in sedi.dat
     *E0, // constant/basic erosion rate [Nlayer]-> given in sedi.dat
     *Taue, // Erosion critical stress [Nlayer]-> given in sedi.dat
     *Taud, // deposition critical stress [Nlayer]-> given in sedi.dat
     *Prt,  // Prandtl number [fraction]-> now just assume 1
     //**Layermass, // Layer total mass for each layer [cell][Nlayer]
     **Layerthickness, // Layer thickness [cell][Nlayer] //added
     *Thickness, // Layerthickness [Nlayer]-> given in sedi.dat
     *Softhard,  // the layer is soft or hard-> given in sedi.dat [Layer] 
     *Drydensity, // dry density for each layer [Nlayer]
     Chind, // concentration for hindered settling velocity
     Cfloc; // concentration for flocuation for settling velocity
int Nlayer, // number of bed layer -> given in sedi.dat
    Nsize,  // number of fractions -> given in sedi.dat
    WSconstant, // if 1,just use constant settling velocity, 0 otherwise -> given in sedi.dat
    //SETsediment, // because we only give 5 layers properties in sedi.dat, if there is more layers, you can set it as 1, and give the relevant value in SetSediment function. It means you should change the codes 
    ParabolKappa, // whether to use parabolic tubulent diffusivity
    bedInterval, // the interval steps to update bed change
    bedComplex,  // whether consider the possibility to flush away a whole layer
    TBMAX,       // whether to output the tb for each cell
    readSediment; // if 1, we will read sediment file as the IC for sediment Concentration, Now just support Nsizemax=1
  FILE *LayerthickFID, **SedimentFID, *SeditbFID, *SeditbmaxFID;
} sedimentsT;

// Globally allocate the pointer to the sediments structure
sedimentsT *sediments;

void ReadSediProperties(int myproc);
void InitializeSediment(gridT *grid, physT *phys, propT *prop,int myproc);
void AllocateSediment(gridT *grid, int myproc);
//void SetSediment(gridT *grid, int myproc); // used by initializesediment
void SettlingVelocity(gridT *grid, physT *phys, propT *prop, int myproc);
void FreeSediment(gridT *grid, int myproc);
void CalculateErosion(gridT *grid, physT *phys, propT *prop,  int myproc);
void CalculateDeposition(gridT *grid, physT *phys, int myproc); //used by calculateerosion and Bedchange
void BedChange(gridT *grid, physT *phys, propT *prop, int myproc);
void SedimentSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop,int Nosize, REAL theta);
void SedimentVerticalVelocity(gridT *grid, physT *phys,int Nosize,int symbol, int myproc);
void OpenSediFiles(propT *prop, int myproc);
void OutputSediment(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);
void CalculateSediDiffusivity(gridT *grid, physT *phys,int Nosize,int myproc); 
void ISendRecvSediBedData3D(REAL **celldata, gridT *grid, int nlayer, int myproc, MPI_Comm comm);
void ComputeSediments(gridT *grid, physT *phys, propT *prop, int myproc, int numprocs, int blowup, MPI_Comm comm);

#endif


