/*
 * File: sources.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Right-hand sides for momentum and scalars, in the form rhs = dt*SOURCE at time step n.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "phys.h"
#include "sources.h"
#include "memory.h"

void MomentumSource(REAL **usource, gridT *grid, physT *phys, propT *prop) {
  int j, jptr, nc1, nc2, k;
  REAL Coriolis_f, ubar, depth_face;

  /* This is the sponge layer */
  /*
  if(prop->sponge_distance) {
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      ubar = 0;
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
	ubar += grid->dzf[j][k]*phys->u[j][k];
	depth_face += grid->dzf[j][k];
      }
      ubar/=depth_face;
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	usource[j][k]-=prop->dt*exp(-4.0*rSponge[j]/prop->sponge_distance)/
	  prop->sponge_decay*(phys->u[j][k]-ubar);
    }
  }
  */

  /* Coriolis for a 2d problem */
  /*
  if(prop->n==prop->nstart+1) {
    printf("Initializing v coriolis\n");
    v_coriolis = (REAL **)SunMalloc(grid->Ne*sizeof(REAL *),"MomentumSource");
    for(j=0;j<grid->Ne;j++)
      v_coriolis[j] = (REAL *)SunMalloc(grid->Nke[j]*sizeof(REAL),"MomentumSource");
  }

  // Hard-code coriolis here so that it can be zero in the main code
  Coriolis_f=4.9745e-5;

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      
      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	usource[j][k]+=prop->dt*Coriolis_f*(v_coriolis[j][k]*grid->n1[j]-
					    InterpToFace(j,k,phys->uc,phys->u,grid)*grid->n2[j]);
  }

  for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];

      for(k=grid->etop[j];k<grid->Nke[j];k++) 
	v_coriolis[j][k]-=prop->dt*Coriolis_f*InterpToFace(j,k,phys->uc,phys->u,grid);
  }
  */
}
/*
 * Function: SaltSource
 * --------------------
 * Usage: SaltSource(A,B,grid,phys,prop,met)
 *
*/
void SaltSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, metT *met) {
    int i, iptr, k, ktop, gc;
    if(prop->metmodel>0 && prop->metmodel < 4){	
	    REAL L_w, EP, dztop, dzmin_saltflux;
	    L_w = 2.50e6;            // Latent heat of evaporation of water (J kg^{-1})
	    dzmin_saltflux=0.1;	

	    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
		i = grid->cellp[iptr];
		ktop = grid->ctop[i];

		//Zero-term below the surface
		for(k=ktop;k<grid->Nk[i];k++)
		    A[i][k]=B[i][k]=0.0;

		dztop = Max(dzmin_saltflux,grid->dzz[i][ktop]);
		EP =   (met->Hl[i]/L_w - met->rain[i])/RHO0; // m/s

		//B[i][ktop] = EP*RHO0 / (dztop);  
		//A[i][ktop] = EP*phys->s[i][ktop] / (dztop);  
		//A[i][ktop] = -EP*phys->s[i][ktop]/(RHO0*dztop);  
		B[i][ktop] += EP /(dztop);
		met->EP[i]=EP;
	   }
   }else{
	for(i=0;i<grid->Nc;i++)
	  for(k=0;k<grid->Nk[i];k++)
	    A[i][k]=B[i][k]=0;
   }
}
/*
 * Function: HeatSource
 * Usage: HeatSource(A,B,grid,phys,prop);
 * --------------------------------------
 * Source terms for heat equation of the form
 *
 * dT/dt + u dot grad T = d/dz ( kappaT dT/dz) + A + B*T
 *
 * Assuming adiabatic top and bottom.  All non-penetrative fluxes are computed as
 * fluxes through the surface and converted to equivalent sources in the top cell.
 * The shortwave radiation is the only term that penetrates into the water column
 * and acts as a real source term.  For a description  of the linearization of
 * the non-penetrative fluxes please see the latex directory.
 *
 * Terms are based on those described in the paper by Wood et al., 2008, 
 * "Modeling hydrodynamics and heat transport in upper Klamath Lake, Oregon, 
 * and implications for water quality", USGS Scientific Investigations Report 2008-5076.
 *
 * Note that they must be set to zero if there is no heat
 * flux since these are temporary variables with values that may
 * be assigned elsewhere.
 *
 */

void HeatSource(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, metT *met, int myproc, MPI_Comm comm) {
 if(prop->metmodel==1){ //Wood et al., Heat Flux model
  int i, iptr, k, ktop, gc;
  REAL z, dztop, sigma , epsilon_w, epsilon_a, dzmin_heatflux, T_0, T_00, c_p, alpha_0,
    r_LW, alpha_LW, F, e_s, e_a, alpha_E1, alpha_E2, alpha_E3, alpha_E4,
    alpha_E5, L_w;
  REAL H_LE, H_LE0, Delta_H_LE0, H_E, H_E0, Delta_H_E0, 
    H_S, H_S0, Delta_H_S0, C_B, r_p, F_SW, r_SW, k_e;
  REAL fixedT = 5.0; // Temperature about which we might want to linearize, if not using previous timestep value 
  
  //
  REAL CP_WATER = 4186.0;
  REAL R = 0.0;
  REAL RSW = 0.0;
 
  REAL Ta , M , U2, rh ,H_SW, H_LW , lap, mslp;
  
  sigma = 5.67e-8;        // Boltzmann constant (W m^{-2} K^{-4})
  epsilon_w = 0.97;       // Emissivity of water
  dzmin_heatflux = 0.01;  // Minimum allowable depth of top cell for computation of 
                          // non-penetrative fluxes
  T_0 = 273.16;           // Zero C (K)
  T_00 = 237.3;           // Used to calculate latent and sensible heat fluxes
  c_p = 4186.0;             // Specific heat of water at standard conditions (J kg^{-1})
  r_LW = 0.02;            // Fraction of longwave radiation reflected at surface
  L_w = 2.5e6;            // Latent heat of evaporation of water (J kg^{-1})
  C_B = 0.61;             // Bowen ratio (mb C^{-1})
  r_p = 1.0;                // r_p = (p_a/p_sl) = local atmospheric pressure/sea-level atmospheric pressure
  k_e = 2.9;              // Extinction coefficient for shortwave radiation (m^{-1})
  r_SW = 0.05;            // Fraction of shortwave radiation reflected from surface

  alpha_0 = 0.937e-5;     // Proportionality constants ([alpha_0]= K^{-2}, all others dimensionless)
  alpha_E1 = 1.0;           // For consistency, kept this but Wood et al.'s formula is off by 9 orders of magnitude
                          // Wood suggests F = alpha_E1*(alpha_E2 + alpha_E3*U2), but Martin et al. suggest 
                          // F = alpha_E2 + alpha_E3*U2 with much smaller numbers
  alpha_E2 = 2.81e-9;     // Formula from Martin et al. "Hydrodynamics and transport for water quality modeling" 
                          // (in mb^{-1} m s^{-1}
   alpha_E3 = 0.14e-9;     // From Martin et al. 1999
  // Original Values in Wood et al
 // alpha_E1 = 0.26;
 // alpha_E2 = 0.50;
 // alpha_E3 = 0.54;
  //
  alpha_E4 = 6.11;
  alpha_E5 = 17.3;
  alpha_LW = 0.17;
  
  


  /* 
   * Atmospheric properties at the current time step, spatially distributed over
   * the surface in arrays:
   *   Ta = atmospheric temperature
   *   M = percentage cloud cover
   *   U2 = magnitude of wind at 2 m above surface, 
   *   rh = relative humidity (0<rh<1)
   *
   */
 // GetAtmosphericProperties(Ta,M,U2,rh,grid,phys,prop);

  /*
   * Incoming longwave radiation
   * 
   * Horizontal variability depends only on variability of atmosphere
   * and not on water temperature.
   *
   */
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    //fprintf(stderr,"globalcell = %d, u2 = %f, Q0 = %f, H_LW = %f, rh = %f, T_air = %f, lap = %f, mslp = %f\n",grid->mnptr[i],met->u2[i],met->Q0[i],met->longwave[i], met->rh[i],met->T_air[i], met->Pa[i], met->sl[i]);
    ktop = grid->ctop[i];
    M = met->cloud[i];
    Ta = met->Tair[i];
    U2 = sqrt( pow(met->Uwind[i],2.0) + pow(met->Vwind[i],2.0) );
    met->Hsw[i] = shortwave(prop->nctime/86400.0,prop->latitude,met->cloud[i]);
    H_SW = met->Hsw[i];
    rh = met->RH[i]/100.0;
    //H_LW = met->Hlw[i];
    // Emissivity of air
   epsilon_a = alpha_0*(1+alpha_LW*M*M)*pow(Ta + T_0,2.0);
    
    // Incoming longwave
    // H_LW = epsilon_a*sigma*(1-r_LW)*pow(Ta[i]+T_0,2);
    H_LW = epsilon_a*sigma*(1-r_LW)*pow(Ta+T_0,4.0);

    B[i][ktop] = +0;  // No dependence on water temperature
    A[i][ktop] = +H_LW;
  }

  /*
   * Longwave emitted radiation
   *
   */
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    
    ktop = grid->ctop[i];

    // Longwave emitted radiation
    H_LE = -epsilon_w*sigma*pow(T_0 + phys->T[i][ktop],4);
    //H_LE = epsilon_w*sigma*pow(T_0 + fixedT,4);

    // Terms in linearization
    Delta_H_LE0 = 4*H_LE/(T_0 + phys->T[i][ktop]);
    H_LE0 = H_LE - Delta_H_LE0*phys->T[i][ktop];
 
    //Delta_H_LE0 = 4*H_LE/(T_0 + fixedT);
    //H_LE0 = H_LE - Delta_H_LE0*fixedT;
    
    B[i][ktop] += -Delta_H_LE0;
    A[i][ktop] += -H_LE0;
	  
    met->Hlw[i] = H_LW+H_LE;
  }

  /*
   * Sensible (H_S) and latent/evaporative heat flux (H_E)
   *
   */
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    ktop = grid->ctop[i];

    // Evaporative radiation
    F = alpha_E1*(alpha_E2 + alpha_E3*U2);
    e_s = alpha_E4*exp(alpha_E5*phys->T[i][ktop]/(T_00 + phys->T[i][ktop]));

    //e_s = alpha_E4*exp(alpha_E5*fixedT/(T_00 + fixedT));
    e_a = alpha_E4*exp(alpha_E5*Ta/(T_00 + Ta))*rh;
    H_E = -RHO0*L_w*F*(e_s - e_a);

    // Terms in linearization
    Delta_H_E0 = RHO0*L_w*F*alpha_E5*e_s*T_00/pow(phys->T[i][ktop]+T_00,2);
    H_E0 = H_E - Delta_H_E0*phys->T[i][ktop];

    //Delta_H_E0 = RHO0*L_w*F*alpha_E5*e_s*T_00/pow(fixedT+T_00,2);
    //H_E0 = H_E - Delta_H_E0*fixedT;

    // Put flux into top cell
    B[i][ktop] += -Delta_H_E0;
    A[i][ktop] += H_E0;
    met->Hl[i] = H_E;


    // Sensible heat flux
//    if(e_s!=e_a && phys->T[i][ktop]!=Ta && H_E!=0 && fixedT!=Ta) {
//      r_p = met->Pa[i] / met->sl[i];
      r_p=1.0;  	
//      H_S = H_E*C_B*r_p*(phys->T[i][ktop]-Ta);
      H_S = H_E*C_B*r_p*(phys->T[i][ktop]-Ta)/(e_s-e_a);
//      H_S = H_E*C_B*r_p*(fixedT-Ta[i]);
      
      // Terms in linearization
      Delta_H_S0 = H_S*(Delta_H_E0/H_E 
			+ 1/(phys->T[i][ktop]-Ta)
			- alpha_E5*e_s/(e_s-e_a)*T_00/pow(phys->T[i][ktop]+T_00,2));
      H_S0 = H_S - Delta_H_S0*phys->T[i][ktop];

      // Put flux into top cell
      B[i][ktop] += -Delta_H_S0;
      A[i][ktop] += H_S0;

      met->Hs[i]=H_S;
//   }

  }

  // Now divide by rho_0 c_p \Delta z
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    ktop = grid->ctop[i];
    dztop = Max(dzmin_heatflux,grid->dzz[i][ktop]);
     
    A[i][ktop]/=(RHO0*c_p*dztop);
    B[i][ktop]/=(RHO0*c_p*dztop);
  }

  /* 
   * Incoming shortwave radiation
   *
   */
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
    ktop = grid->ctop[i];

    // Make sure values beneath surface are zero
    for(k=ktop+1;k<grid->Nk[i];k++) //originally this should be ktop+1 
      A[i][k]=B[i][k]=0;
     
    REAL topface, botface, depth=0; //These will be used with the heat conservative scheme of solar radiation
    REAL ksw1 = 1/prop->Lsw; //These are extinction coefficients. L1 and L2 are in met.h  
    
    REAL wave1, wave1term1, wave1term2, wave1term3, wave2, wave2term1, wave2term2, wave2term3; //These are the terms from the infinite reflections formula

     for(k=ktop;k<grid->Nk[i];k++) //loop for calculating the local depth
{	depth = depth + grid->dzz[i][k]; 
}
	//fprintf(stderr,"i = %d, Local depth = %f\n", i, depth);
    z=0;
    for(k=ktop;k<grid->Nk[i];k++) {
      z-=grid->dzz[i][k]/2;

      topface = z + grid->dzz[i][k]/2; //Depth of the top face of a cell from the surface
      botface = z - grid->dzz[i][k]/2; //Depth of the bottom face of a cell from the surface

      //F_SW = k_e*H_SW*(1-r_SW)*exp(k_e*z)/(RHO0*c_p);

      wave1term1 = (exp(topface*ksw1) - exp(botface*ksw1))/(topface - botface);
      wave1term2 = (1/(1-exp(-2*depth*ksw1))) * (exp(-(2*depth - topface)*ksw1) - exp(-(2*depth - botface)*ksw1))/(topface - botface);
      wave1term3 = -(1/(1-exp(-2*depth*ksw1))) * (exp(-(2*depth + topface)*ksw1) - exp(-(2*depth + botface)*ksw1))/(topface - botface);  

      wave2 = wave2term1 + wave2term2 + wave2term3;

      F_SW =wave1*H_SW/(RHO0*CP_WATER);

      // Set the coefficients of the sources term here
      B[i][k] += 0;
      A[i][k] += F_SW;
      
      z-=grid->dzz[i][k]/2;
    }
  }
}else if(prop->metmodel==2 || prop->metmodel==3){ // COARE3.0 HeatFlux Algorithm or constant flux coefficients
   int i, k, ktop, iptr;
   REAL dztop;
   REAL *H0=met->Htmp;

   REAL c_p = 4186.0;     // Specific heat of water at standard conditions (J kg^{-1})
   REAL rhocp = RHO0*c_p;
   REAL dHdT;
   REAL eps = 1e-14;
   REAL dzmin_heatflux = 0.01;  // Minimum allowable depth of top cell for computation of   
   // shortwave constants
   REAL ksw1; // light extinction coefficient
   REAL depth, z, topface, botface;
   REAL wave1, wave1term1, wave1term2, wave1term3;
   REAL F_SW;
   REAL dT;
   
   
  // for(j=0;j<grid->Nc;j++){
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
     ktop = grid->ctop[i];
     // Make sure values beneath surface are zero
      for(k=ktop;k<grid->Nk[i];k++)
         A[i][k]=B[i][k]=0.0;
  
     // Set the flux terms in the top cell for the present time step
     H0[i] = ( met->Hs[i] + met->Hl[i] + met->Hlw[i] );
     
     // Calculate T+dt array
    //if(phys->dT[i]<=0.0){
    //   phys->dT[i]=Min(-1e-4,phys->dT[i]);
    // }else{
    //   phys->dT[i]=Max(1e-4,phys->dT[i]);
    // }
    phys->dT[i] = 0.01; //Hard-wire for stability 
    
	
     phys->Ttmp[i][ktop] = phys->T[i][ktop] + phys->dT[i];
   }
   
   //Communicate the temporary Temperature array
   ISendRecvCellData3D(phys->Ttmp,grid,myproc,comm);

   // Evaluate flux terms when the temperature is T + dT
   updateAirSeaFluxes(prop, grid, phys, met, phys->Ttmp);
   
   // Evaluate the implicit source term
  // for(j=0;j<grid->Nc;j++){  
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];
     ktop = grid->ctop[i];
     
     // Put the temperature gradient flux terms into B
     dHdT = ( (met->Hs[i] + met->Hl[i] + met->Hlw[i]) - H0[i]) /(phys->dT[i] + SMALL);
	   
     A[i][ktop] = H0[i] - dHdT * phys->T[i][ktop];
     //A[i][ktop] = -(H0[i] - dHdT * phys->T[i][ktop]);//Testing
     B[i][ktop] = -dHdT;

     dztop = Max(dzmin_heatflux,grid->dzz[i][ktop]);
     A[i][ktop]/=(rhocp*dztop);
     B[i][ktop]/=(rhocp*dztop);

     // Evaluate the shortwave radiation terms and put into A
     for(k=ktop;k<grid->Nk[i];k++) //loop for calculating the local depth (accounts for free surface) 
 	depth = depth + grid->dzz[i][k]; 
     
     //if (depth<dztop)
     //	printf("Warning in sources.c: Depth at cell %d = %6.10f (<%6.10f)\n",j,depth,dztop);
     
     // Set the light extinction coefficient - ensure that thhe extinction depth is less than the water depth
     ksw1 = 1.0 / Min(prop->Lsw, 0.5*depth);

     z=0.0;
     for(k=ktop;k<grid->Nk[i];k++) {
       z-=grid->dzz[i][k]/2.0;
 
       topface = z + grid->dzz[i][k]/2.0; //Depth of the top face of a cell from the surface
       botface = z - grid->dzz[i][k]/2.0; //Depth of the bottom face of a cell from the surface
 
       wave1term1 = (exp(topface*ksw1) - exp(botface*ksw1))/(topface - botface);
       wave1term2 = (1.0/(1.0-exp(-2.0*depth*ksw1))) * (exp(-(2.0*depth - topface)*ksw1) - exp(-(2.0*depth - botface)*ksw1))/(topface - botface);
       wave1term3 = -(1.0/(1.0-exp(-2.0*depth*ksw1))) * (exp(-(2.0*depth + topface)*ksw1) - exp(-(2.0*depth + botface)*ksw1))/(topface - botface);  
 
       wave1 = wave1term1 + wave1term2 + wave1term3;
       
       F_SW = wave1 * met->Hsw[i] / rhocp;
       
       B[i][k] += 0.0;
       A[i][k] += F_SW;
       
       z-=grid->dzz[i][k]/2;

       //Zero for Testing
       //B[i][k] = 0.0;
       //A[i][k] = 0.0;

     }
   }
 }else{ //Set flux terms to zero
  int i,k,iptr;
    //for(i=0;i<grid->Nc;i++)
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      for(k=0;k<grid->Nk[i];k++){
	A[i][k]=B[i][k]=0;
      }
    }
 }// End if
} // End Heatsource

/*
 * Function: InitSponge
 * Usage: InitSponge(grid,myproc);
 * -------------------------------
 * Apply a sponge layer to all type 2 boundaries.
 *
 */
void InitSponge(gridT *grid, int myproc) {
  int Nb, p1, p2, mark, g1, g2;
  int j, n, NeAll, NpAll;
  REAL *xb, *yb, *xp, *yp, r2;
  char str[BUFFERLENGTH];
  FILE *ifile;

  NeAll = MPI_GetSize(EDGEFILE,"InitSponge",myproc);
  NpAll = MPI_GetSize(POINTSFILE,"InitSponge",myproc);

  xp = (REAL *)SunMalloc(NpAll*sizeof(REAL),"InitSponge");
  yp = (REAL *)SunMalloc(NpAll*sizeof(REAL),"InitSponge");
  rSponge = (REAL *)SunMalloc(grid->Ne*sizeof(REAL),"InitSponge");

  // Read in points on entire grid
  ifile = MPI_FOpen(POINTSFILE,"r","InitSponge",myproc);
  for(j=0;j<NpAll;j++) {
    xp[j]=getfield(ifile,str);
    yp[j]=getfield(ifile,str);
    getfield(ifile,str);
  }
  fclose(ifile);

  // Count number of nonzero boundary markers on entire grid
  ifile = MPI_FOpen(EDGEFILE,"r","InitSponge",myproc);
  Nb = 0;
  for(j=0;j<NeAll;j++) {
    fscanf(ifile, "%d %d %d %d %d",&p1,&p2,&mark,&g1,&g2);
    if(mark==2 || mark==3)
      Nb++;
  }
  fclose(ifile);

  xb = (REAL *)SunMalloc(Nb*sizeof(REAL),"InitSponge");
  yb = (REAL *)SunMalloc(Nb*sizeof(REAL),"InitSponge");

  n=0;
  ifile = MPI_FOpen(EDGEFILE,"r","InitSponge",myproc);
  for(j=0;j<NeAll;j++) {
    fscanf(ifile, "%d %d %d %d %d",&p1,&p2,&mark,&g1,&g2);
    if(mark==2 || mark==3) {
      xb[n]=0.5*(xp[p1]+xp[p2]);
      yb[n]=0.5*(yp[p1]+yp[p2]);
      n++;
    }
  }
  fclose(ifile);  

  // Now compute the minimum distance between the edge on the
  // local processor and the boundary and place this in rSponge.
  for(j=0;j<grid->Ne;j++) {
    rSponge[j]=INFTY;

    for(n=0;n<Nb;n++) {
      r2=pow(xb[n]-grid->xe[j],2)+pow(yb[n]-grid->ye[j],2);
      if(r2<rSponge[j])
	rSponge[j]=r2;
    }
    rSponge[j]=sqrt(rSponge[j]);
    //    printf("Processor %d: rSponge[%d]=%f\n",myproc,j,rSponge[j]);
  }
}
  
