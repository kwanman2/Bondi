//Modified by Alexander Tchekhovskoy: MPI+3D
/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"
#include <float.h>


void coord_transform(double *pr,int i, int j,int k) ;
double compute_B_from_A( double (*A)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] );
double normalize_B_by_maxima_ratio(double beta_target, double (*p)[N2M][N3M][NPR], double *norm_value);
double normalize_B_by_beta(double beta_target, double (*p)[N2M][N3M][NPR], double rmax, double *norm_value);
double normalize_B_ratio_bondi(double beta_target, double(*p)[N2M][N3M][NPR], double rmin, double* norm_value);


/////////////////////
//magnetic field geometry and normalization
#define NORMALFIELD (0)
#define VERTICALFIELD (1)
#define MONOPOLE (2)
#define FIELDJONMAD (3)
#define NOFIELD (4)

#define WHICHFIELD FIELDJONMAD

#define NORMALIZE_FIELD_BY_MAX_RATIO (1)
#define NORMALIZE_FIELD_BY_BETAMIN (2)
#define NORMALIZE_FIELD_BONDI (3)
#define WHICH_FIELD_NORMALIZATION NORMALIZE_FIELD_BONDI
//end magnetic field        //tom: don't use NORMALIZE_FIELD_BY_BETAMIN in harmpi
//////////////////////


//double rmax = 0.;
//double rhomax = 1.;

void init()
{
  void init_bondi(void);
  void init_quasi_bondi(void);
  void init_monopole(double Rout_val);

  switch( WHICHPROBLEM ) {
  case MONOPOLE_PROBLEM_2D:
    init_monopole(1e3);
    break;
  case QUASI_BONDI_PROBLEM:
    init_quasi_bondi();
    break;
  case BONDI_PROBLEM:
    init_bondi();
    break;
  }

}

double max(double A, double B) {
    return (A > B) ? A : B;
}

#define M_PI 3.14159265358979323846

//tom: haven't do proper magnetic field normalization
void init_quasi_bondi()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;

  /* for magnetic field */
  double A[N1+D1][N2+D2][N3+D3] ;
  double rho_av,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  
  int iglob, jglob, kglob;
  double rancval;
  
  double amax, aphipow;
  
  double Amin, Amax, cutoff_frac = 0.01;

  double kappa, ssth;
  double rc, mdot, rhor, lb, z1, z2, isco;
  double E, n, f1, f2, rs;
  
  /* some physics parameters */
  gam = 4./3. ;

  /* disk parameters */
  a = 0.9 ;

  //kappa =1.e-3; to be calculated
  beta = 10;

  /* some numerical parameters */
  lim = MC ;
  failed = 0 ;	/* start slow */
  cour = .8 ;
  dt = 1.e-5 ;
  R0 = 0.0 ;
  Rin = 0.8*(1. + sqrt(1. - a*a)) ;  //.98 //0.87
  Rout = 1e5;
  rbr = 200.;                   //tom
  npow2=4.0; //power exponent
  cpow2=1.2; //exponent prefactor (the larger it is, the more hyperexponentiation is)   //tom
  rhor = (1. + sqrt(1. - a * a));
  //mdot = 1e3;   //tom
  mdot = 1000;
  rc = 0.0;

  z1 = 1 + pow(1 - a * a, 1. / 3.) * (pow(1 + a, 1. / 3.) + pow(1 - a, 1. / 3.));
  z2 = pow(3 * a * a + z1 * z1, 1. / 2.);
  isco = 3 + z2 - (a / abs(a)) * pow((3 - z1) * (3 + z1 + 2 * z2), 1. / 2.);

  //angular momentum below bondi radius in the equatorial plane (assuming angular momentum is conserved in the infall phase)
  lb = sqrt(rc);

  t = 0. ;
  hslope = 1.0 ;

  /* notes: BL=1 grid (coord.c)
   if( X[1] > x1br ) {
    theexp += cpow2 * pow(X[1]-x1br,npow2);  //hyperexponential for X[1] > x1br
    
  }
  
  V[0] = X[0];
  V[1] = exp(theexp) + R0 ;
  V[2] = M_PI_2*(1.0+X[2]) + ((1. - hslope)/2.)*sin(M_PI*(1.0+X[2])) ;
  V[3] = X[3] ;
  ...
  x1br = log( rbr - R0 );
  */

  if(N2!=1) {
    //2D/3D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;
  
  set_arrays() ;
  set_grid() ;

  get_phys_coord(4,0,0,&r,&th,&phi) ;
  if(MASTER==mpi_rank) {
    fprintf(stderr,"r[4]: %g\n",r) ;
    fprintf(stderr,"r[4]/rhor: %g",r/(1. + sqrt(1. - a*a))) ;
    if( r > 1. + sqrt(1. - a*a) ) {
      fprintf(stderr, ": INSUFFICIENT RESOLUTION, ADD MORE CELLS INSIDE THE HORIZON\n" );       //tom: it means should put at least 5 cells inside the horizon
    }
    else {
      fprintf(stderr, "\n");
    }
  }

  E = 1.e-7; //tom

  /* output choices */
  tf = 1000000.0 ;

  DTd = 10.; /* dumping frequency, in units of M */
  DTl = 10. ;	/* logfile frequency, in units of M */
  DTi = 10. ; 	/* image file frequ., in units of M */
  DTr = 10. ; /* restart file frequ., in units of M */
  DTr01 = 100. ; /* restart file frequ., in timesteps */

  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;

  umax = 0;

  //ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
  for(iglob=0;iglob<mpi_ntot[1];iglob++) {
    for(jglob=0;jglob<mpi_ntot[2];jglob++) {
      for(kglob=0;kglob<mpi_ntot[3];kglob++) {
        
        rancval = ranc(0);
        i = iglob-mpi_startn[1];
        j = jglob-mpi_startn[2];
        k = kglob-mpi_startn[3];
        if(i<0 ||
           j<0 ||
           k<0 ||
           i>=N1 ||
           j>=N2 ||
           k>=N3){
          continue;
        }
        get_phys_coord(i,j,k,&r,&th,&phi) ;

        ssth = sin(th)* sin(th);   //angular distribution


        if(r > 1.5 * rhor) {        //tom: the inner edge of disk
            
            //ur = -ur0 / sqrt(r / rhor);     //analytic bondi approximate solution
            ur = -sqrt(2 / r);   //free-fall
            uh = 0.;
            if (r < isco && r < rc) {
                up = ssth * sqrt(r / (isco * isco));     //approach to 0 below ISCO
            }
            else if (r < rc) {
                up = ssth / sqrt(r); //keplerian
            }
            else {
                up = lb * ssth / r;  //(assuming angular momentum is conserved in the infall phase)
            }

            rho = mdot / (-4 * M_PI * r * r * ur);
            /*
            rs = 1000; //tom: guess r_s, arbitrary
            double l = lb * ssth;
            for (n = 1; n < 10; n++) {
                f1 = E - pow(l, 2) / (2 * pow(rs, 2)) + 1 / (2 * (rs - 1)) - ((gam + 1) / (2 * (gam - 1))) * (rs / (4 * pow(rs - 1, 2)) - pow(l, 2) / (2 * pow(rs, 2)));
                f2 = pow(l, 2) / pow(rs, 3) - 1 / (2 * pow(rs - 1, 2)) - ((gam + 1) / (2 * (gam - 1))) * (1 / (4 * pow(rs - 1, 2)) - rs / (2 * pow(rs - 1, 3)) + pow(l, 2) / pow(rs, 3));
                rs = rs - f1 / f2;
            }
            rs *= 2;     //tom: fix, based on the values from Sukova 2017, somehow the radius is always underestimated by the factor of 2, I haven't found a problem from the above algorithm yet
            */
            rs = (5 - 3 * gam) / (4 * E * (gam - 1));

            kappa = pow(4 * M_PI * rs * rs * pow(sqrt(2 / rs), (2 / (gam - 1) + 1)) / (mdot * pow(gam, 1 / (gam - 1))), (gam - 1));
            //kappa = 1.e-3;   //temp

            if (MASTER == mpi_rank) {
                if (j < 1 && k < 1 && i < 4) {            //check
                    fprintf(stderr, "circularization radius: %g\n", rc);
                    fprintf(stderr, "Estimated bondi radius with E=%g: %g\n", E, rs);
                    fprintf(stderr, "specific entropy: %g\n", kappa);
                }
            }
            u = kappa * pow(rho, gam - 1.) / (gam - 1.);   //tom
            //u = kappa * pow(rho, gam) / (gam - 1.);
            //u = 0.00002;     //floor value, for no angular momentum flows, check

            if (u > umax) umax = u;

            p[i][j][k][RHO] = rho ;
            p[i][j][k][UU] = u * (1. + 5.e-2 * (rancval - 0.5));
            p[i][j][k][U1] = ur ;
            p[i][j][k][U2] = uh ;
            p[i][j][k][U3] = up ;

            /* convert from 4-vel to 3-vel */
            coord_transform(p[i][j][k], i, j, k);

        }
        else { 
            rho = RHOMINLIMIT;
            u = UUMINLIMIT;
            ur = 0;
            uh = 0. ;
            up = 0. ;


          p[i][j][k][RHO] = rho ;
          p[i][j][k][UU] = u ;
          p[i][j][k][U1] = ur ;
          p[i][j][k][U2] = uh ;

          p[i][j][k][U3] = up ;

        }

        p[i][j][k][B1] = 0. ;
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
      }
    }
  }
  

  bound_prim(p) ;

  if (WHICHFIELD == NORMALFIELD) {
    aphipow = 0.;
  }

  double JONMADHPOW, JONMADR0, JONMADRIN, JONMADROUT, rpow;     //tom: for FIELDJONMAD, 
  JONMADHPOW = 4.0;
  JONMADR0 = -5.0;
  //JONMADR0 = -2.0; 
  JONMADRIN = 4.0; //tom
  JONMADROUT = 1000; //jane: transition to monopolar field  /tom: set larger
  rpow = 1.0;

  /* first find corner-centered vector potential */
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) A[i][j][k] = 0. ;
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
    
      //cannot use get_phys_coords() here because it can only provide coords at CENT
      coord(i,j,k,CORN,X) ;
      bl_coord(X,&r,&th,&phi) ;

      if (WHICHFIELD == VERTICALFIELD) {
          /* vertical field version */
          A[i][j][k] = 0.5 * r * sin(th);
      }
      else if (WHICHFIELD == MONOPOLE) {
          /* radial (monopolar) field version */
          A[i][j][k] = (1 - cos(th));
      }
      else if (WHICHFIELD == NORMALFIELD) {
          /* field-in-disk version */
          /* flux_ct */
          rho_av = 0.25 * (
              p[i][j][k][RHO] +
              p[i - 1][j][k][RHO] +
              p[i][j - 1][k][RHO] +
              p[i - 1][j - 1][k][RHO]);

          q = pow(r, aphipow) * rho_av / rhomax;
          q -= 0.2;
          if (q > 0.) A[i][j][k] = q;
      }
      else if (WHICHFIELD == NOFIELD) {
          A[i][j][k] = 0;
      }
      else if (WHICHFIELD == FIELDJONMAD) {             //subject to be checked
          //double JONMADHPOW, JONMADR0, JONMADRIN, JONMADROUT, rpow;   //tom: declare those variables in the upper level in case there are other functions need one of them 
          //JONMADHPOW = 4.0;
          //JONMADR0 = -5.0;
          //JONMADR0 = -2.0; 
          //JONMADRIN = 4.0; //inner edge I guess?
          //JONMADROUT = 300; //jane: transition to monopolar field
          //rpow = 1.0;

          double jonmadhpow;
          double interp;
          interp = (r - rhor) / (r + 1.e-10);
          if (interp < 0.0) {
              interp = 0.0;     // inside BH
          }
          if (interp > 1.0) {
              interp = 1.0;     //tom why do this
          }

          jonmadhpow = JONMADHPOW * interp + 0.0 * (1.0 - interp);

          double thetafactout, thetafactin, thetafact1;
          thetafactin = 1.0 - cos(th);
          if (th > 0.5 * M_PI) {
              thetafactin = 1.0 - cos(M_PI - th);
          }
          thetafactout = pow(sin(th), 1.0 + jonmadhpow);              // (sin£c)^(1+h)
          thetafact1 = thetafactout * interp + thetafactin * (1.0 - interp);  //tom why do this

          double interp2;
          interp2 = (r - JONMADROUT) / (r + 1.e-10);
          if (interp2 < 0.0) {
              interp2 = 0.0;
          }
          if (interp2 > 1.0) {
              interp2 = 1.0;
          }

          double thetafactout2, thetafactin2, thetafact2;
          thetafactin2 = thetafact1;
          thetafactout2 = 1.0 - cos(th);      //monopolar
          if (th > 0.5 * M_PI) { thetafactout2 = 1.0 - cos(M_PI - th); }
          thetafact2 = thetafactout2 * interp2 + thetafactin2 * (1.0 - interp2);

          double thetafact;
          thetafact = thetafact2; //useless

          double rfact;
          rfact = max(pow(r - JONMADR0, rpow) * 1E40 - 0.02, 0.0);    //MAX(r^£h 10^40 - 0.02, 0)
          if (r >= JONMADROUT) {
              rfact = max(pow(JONMADROUT - JONMADR0, rpow) * 1E40 - 0.02, 0.0); //to monopolar  //MAX(r_0^£h 10^40 - 0.02, 0)
          }
          if (r <= JONMADRIN) {
              rfact = max(pow(JONMADRIN - JONMADR0, rpow) * 1E40 - 0.02, 0.0); //inner edge
          }


          q = rfact * thetafact;

          if (th<1.e-5 || th>M_PI - 1.e-5) {
              q = 0;
          }
          if (q > 0.) { A[i][j][k] = q; }

      }
      else {
          // OLEKFIELD
          rho_av = 0.25 * (
              p[i][j][k][RHO] +
              p[i - 1][j][k][RHO] +
              p[i][j - 1][k][RHO] +
              p[i - 1][j - 1][k][RHO]);

          q = (pow(r, 4.0) * pow(rho_av, 2.0) * 1e40 - 0.02) * pow(sin(th), 4.0);
          if (q > 0.) A[i][j][k] = q;
      }
  }
  
  fixup(p) ;

  //check
  /* now differentiate to find cell-centered B,
     and begin normalization */
  
  bsq_max = compute_B_from_A(A,p);
  
  if(WHICHFIELD == NORMALFIELD) {
    if(MASTER==mpi_rank)
      fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

    /* finally, normalize to set field strength */
    beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;

    if(MASTER==mpi_rank)
      fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
    
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(beta, p, 10*rmax, &norm);
    }
    else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        MPI_Finalize();
        exit(2345);
      }
    }

    if(MASTER==mpi_rank)
      fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;
  }

  if (WHICHFIELD == FIELDJONMAD) {
      /*if (MASTER == mpi_rank)
          fprintf(stderr, "initial b^2_max: %g\n", bsq_max);*/

      /* finally, normalize to set field strength */
      /*beta_act = (gam - 1.) * umax / (0.5 * bsq_max);

      if (MASTER == mpi_rank) {
          fprintf(stderr, "initial beta: %g (target: %g)\n", beta_act, beta);
          fprintf(stderr, "umax: %g \n", umax);
          fprintf(stderr, "bsq_max: %g \n", bsq_max);
      }*/
      if (WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {    //tom: don't use this
          beta_act = normalize_B_by_beta(beta, p, 0.01 * Rout, &norm);
      }
      else if (WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
          beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
      }
      else if (WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BONDI) {
          beta_act = normalize_B_ratio_bondi(beta, p, JONMADRIN, &norm);   //take reference at the density edge
      }
      else {
          if (i_am_the_master) {
              fprintf(stderr, "Unknown magnetic field normalization %d\n",
                  WHICH_FIELD_NORMALIZATION);
              MPI_Finalize();
              exit(2345);
          }
      }

      if (MASTER == mpi_rank)
          fprintf(stderr, "final beta: %g (target: %g)\n", beta_act, beta);
  }
    
  /* enforce boundary conditions */
  fixup(p) ;
  bound_prim(p) ;




#if( DO_FONT_FIX )
  set_Katm();
#endif 


}



//note that only axisymmetric A is supported       
double compute_B_from_A( double (*A)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] )
{
  double bsq_max = 0., bsq_ij ;
  int i, j, k;
  struct of_geom geom;
  
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;
    
    /* flux-ct */
    p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
                       + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
    p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
                      - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;
    
    p[i][j][k][B3] = 0. ;
    
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  return(bsq_max);
}

double normalize_B_by_maxima_ratio(double beta_target, double (*p)[N2M][N3M][NPR], double *norm_value) //tom: old way
{
  double beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  int i, j, k;
  struct of_geom geom;
  double X[NDIM], r, th, ph;
  
  bsq_max = 0;

  ZLOOP{
    get_geometry(i,j,k,CENT,&geom);                 
    bsq_ij = bsq_calc(p[i][j][k],&geom);
    if (bsq_ij > bsq_max) bsq_max = bsq_ij;
    u_ij = p[i][j][k][UU];
    if (u_ij > umax) umax = u_ij;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  /* finally, normalize to set field strength */
  beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;

  if (MASTER == mpi_rank) {
      fprintf(stderr, "initial max b^2: %g \n", bsq_max);
      fprintf(stderr, "initial max u_g: %g \n", umax);
      fprintf(stderr, "initial beta: %g (target: %g)\n", beta_act, beta_target);
  }

  norm = sqrt(beta_act/beta_target); 
  if (MASTER == mpi_rank) {
      fprintf(stderr, "norm: %g\n", norm);
  }
  bsq_max = 0. ; //find the new b^2
  ZLOOP {
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
    
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;

  if(norm_value) {
    *norm_value = norm;
  }
  return(beta_act);
}

double normalize_B_ratio_bondi(double beta_target, double(*p)[N2M][N3M][NPR], double rmin, double* norm_value) //tom: new
{
    double beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
    double norm;
    int i, j, k;
    struct of_geom geom;
    double X[NDIM], r, th, ph;

    bsq_max = 0;
    double closest = 10; //search for the cell with theta closest to pi/2
    ZLOOP {     
      coord(i, j, k, CENT, X);  
      bl_coord(X, &r, &th, &ph);
      get_geometry(i, j, k, CENT, &geom);       //ph: arbitary unless it is not symetric, but have to pick one range to reduce the search time
      if (r > 0.98*rmin && r < 1.02*rmin && ph < 0.2 && ph > -0.2 && th > 1.57 && th < 1.63) { //search near the mid-plane, I am a bit confused with the numerical grid so I can't just limited to only do r
          if (abs(th - M_PI/2 - 0.06) < closest) {
              closest = th;
              fprintf(stderr, "r: %g \n", r);                                              // , and the cell usually doesn't exactly located at the equator I am looking for
              fprintf(stderr, "theta: %g \n", th);
              //fprintf(stderr, "ph: %g \n", ph);
              bsq_ij = bsq_calc(p[i][j][k], &geom);
              if (bsq_ij > bsq_max) bsq_max = bsq_ij;       //tom: measure the things in the certain region above the edge instead, to make sure it is measured closest to the edge
              u_ij = p[i][j][k][UU];
              if (u_ij > umax) umax = u_ij;
          }
      }
      if (r > 3*rmin) { //tom: should be beyond
          continue;
      }
    }
    
#ifdef MPI
        //exchange the info between the MPI processes to get the true max
    MPI_Allreduce(MPI_IN_PLACE, &umax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &bsq_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //
#endif

    /* finally, normalize to set field strength */
    beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
    //beta_act = (gam - 1.) * u_ij / (0.5 * bsq_max);

    if (MASTER == mpi_rank) {
        fprintf(stderr, "initial max b^2: %g \n", bsq_max);
        fprintf(stderr, "initial max u_g: %g \n", umax);
        fprintf(stderr, "initial beta: %g (target: %g)\n", beta_act, beta_target);
    }

    norm = sqrt(beta_act/beta_target); 
    //norm = sqrt(beta_act / beta_target) * 80;
    if (MASTER == mpi_rank) {
        fprintf(stderr, "norm: %g\n", norm);
    }
    bsq_max = 0.; //find the new b^2
    closest = 10.;
    ZLOOP{
      p[i][j][k][B1] *= norm;
      p[i][j][k][B2] *= norm;
      p[i][j][k][B3] *= norm; //

      coord(i, j, k, CENT, X); //
      bl_coord(X, &r, &th, &ph); //
      get_geometry(i, j, k, CENT, &geom);
      if (r > 0.98 * rmin && r < 1.02 * rmin && ph < 0.2 && ph > -0.2 && th > 1.57 && th < 1.63) {      //temp
          if (abs(th - M_PI / 2 - 0.06) < closest) {
              closest = th;
              bsq_ij = bsq_calc(p[i][j][k], &geom);
              if (bsq_ij > bsq_max) bsq_max = bsq_ij; //
          }
      }
      if (r > 3 * rmin) { //tom: should be beyond
          continue;
      }
    }
#ifdef MPI
        //exchange the info between the MPI processes to get the true max
    MPI_Allreduce(MPI_IN_PLACE, &bsq_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); //
#endif

    beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
    //beta_act = (gam - 1.) * u_ij / (0.5 * bsq_max);

    if (norm_value) {
        *norm_value = norm;
    }
    return(beta_act);
}

//normalize the magnetic field using the values inside r < rmax
double normalize_B_by_beta(double beta_target, double (*p)[N2M][N3M][NPR], double rmax, double *norm_value)
{
  double beta_min = 1e100, beta_ij, beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  int i, j, k;
  struct of_geom geom;
  double X[NDIM], r, th, ph;
  
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th, &ph);
    if (r>rmax) {
      continue;
    }
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    u_ij = p[i][j][k][UU];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ; 
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  /* finally, normalize to set field strength */
  beta_act = beta_min;
  
  norm = sqrt(beta_act/beta_target) ; 
  beta_min = 1e100; 
  ZLOOP { 
    p[i][j][k][B1] *= norm ; 
    p[i][j][k][B2] *= norm ;
    p[i][j][k][B3] *= norm ;
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;  
    u_ij = p[i][j][k][UU];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  beta_act = beta_min;

  if(norm_value) {
    *norm_value = norm;
  }

  return(beta_act);
}


// tom: constant density, no initial velocity, it does not support B-field, insensitive to BH spin (maybe the grid setting problem)
void init_bondi()
{
	int i,j,k ;
	double r,th,phi,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1][N3+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax ;

	/* some physics parameters */
	gam = 4./3. ;

	/* black hole parameters */
        a = 0.9375 ;

	kappa = 1.e-3 ;

	/* radius of the inner edge of the initial density distribution */
	rin = 10.;

        /* some numerical parameters */
        lim = MC ;
        failed = 0 ;	/* start slow */
        cour = 0.9 ;
        dt = 1.e-5 ;
	rhor = (1. + sqrt(1. - a*a)) ;
	R0 = -2*rhor ;
        Rin = 0.5*rhor ;
        Rout = 1e3 ;
        rbr = Rout*10.;
        npow2=4.0; //power exponent
        cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)


        t = 0. ;
        hslope = 1.0 ; //uniform angular grid

	if(N2!=1) {
	  //2D problem, use full pi-wedge in theta
	  fractheta = 1.;
	}
	else{
	  //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
	  fractheta = 1.e-2;
	}
        fracphi = 1.;

        set_arrays() ;
        set_grid() ;

	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;
	fprintf(stderr,"rmin: %g\n",r) ;
	fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        /* output choices */
	tf = Rout ;

	DTd = 2. ;	/* dumping frequency, in units of M */
	DTl = 2. ;	/* logfile frequency, in units of M */
	DTi = 2. ; 	/* image file frequ., in units of M */
        DTr = 50 ; /* restart file frequ., in units of M */
        DTr01 = 1000 ; /* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
        rdump01_cnt = 0 ;
	defcon = 1. ;

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
		coord(i,j,k,CENT,X) ;
		bl_coord(X,&r,&th,&phi) ;

		sth = sin(th) ;
		cth = cos(th) ;

		/* regions outside uniform density distribution */
		if(r < rin) {
			rho = 1.e-7*RHOMIN ;
                        u = 1.e-7*UUMIN ;

			/* these values are demonstrably physical
			   for all values of a and r */
			/*
                        ur = -1./(r*r) ;
                        uh = 0. ;
			up = 0. ;
			*/

			ur = 0. ;
			uh = 0. ;
			up = 0. ;

			/*
			get_geometry(i,j,CENT,&geom) ;
                        ur = geom.gcon[0][1]/geom.gcon[0][0] ;
                        uh = geom.gcon[0][2]/geom.gcon[0][0] ;
                        up = geom.gcon[0][3]/geom.gcon[0][0] ;
			*/

			p[i][j][k][RHO] = rho ;
			p[i][j][k][UU] = u ;
			p[i][j][k][U1] = ur ;
			p[i][j][k][U2] = uh ;
			p[i][j][k][U3] = up ;
		}
		/* region inside initial uniform density */
		else { 
		  rho = 1.;
		  u = kappa*pow(rho,gam)/(gam - 1.) ;
		  ur = 0. ;
		  uh = 0. ;


		  p[i][j][k][RHO] = rho ;
		  if(rho > rhomax) rhomax = rho ;
		  p[i][j][k][UU] = u;
		  if(u > umax && r > rin) umax = u ;
		  p[i][j][k][U1] = ur ;
		  p[i][j][k][U2] = uh ;
		  p[i][j][k][U3] = up ;
		  
		  /* convert from 4-vel to 3-vel */
		  coord_transform(p[i][j][k],i,j,k) ;
		}

		p[i][j][k][B1] = 0. ;
		p[i][j][k][B2] = 0. ;
		p[i][j][k][B3] = 0. ;

	}

	fixup(p) ;
	bound_prim(p) ;
    

#if(0) //disable for now
	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2,0,N3) A[i][j][k] = 0. ;
        ZSLOOP(0,N1,0,N2,0,N3) {
                /* vertical field version */
                /*
                coord(i,j,l,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;

                A[i][j][k] = 0.5*r*sin(th) ;
                */

                /* field-in-disk version */
		/* flux_ct */
                rho_av = 0.25*(
                        p[i][j][RHO] +
                        p[i-1][j][RHO] +
                        p[i][j-1][RHO] +
                        p[i-1][j-1][RHO]) ;

                q = rho_av/rhomax - 0.2 ;
                if(q > 0.) A[i][j][k] = q ;

        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,k,CENT,&geom) ;

		/* flux-ct */
		p[i][j][B1] = -(A[i][j][k] - A[i][j+1][k]
				+ A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
		p[i][j][B2] = (A[i][j][k] + A[i][j+1][k]
				- A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;

		p[i][j][B3] = 0. ;

		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
	norm = sqrt(beta_act/beta) ;
	bsq_max = 0. ;
	ZLOOP {
		p[i][j][k][B1] *= norm ;
		p[i][j][k][B2] *= norm ;

		get_geometry(i,j,k,CENT,&geom) ;
		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;

	/* enforce boundary conditions */
	fixup(p) ;
	bound_prim(p) ;
    
#endif

    

    
#if( DO_FONT_FIX ) 
	set_Katm();
#endif 


}

void init_monopole(double Rout_val)
{
	int i,j,k ;
	double r,th,phi,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax ;

	/* some physics parameters */
	gam = 4./3. ;

	/* disk parameters (use fishbone.m to select new solutions) */
        a = 0.9375 ;

	kappa = 1.e-3 ;
	beta = 1.e2 ;

        /* some numerical parameters */
        lim = MC ;
        failed = 0 ;	/* start slow */
        cour = 0.9 ;
        dt = 1.e-5 ;
	rhor = (1. + sqrt(1. - a*a)) ;
	R0 = -4*rhor;
        Rin = 0.7*rhor ;
        Rout = Rout_val ;
        rbr = Rout*10.;
    npow2=4.0; //power exponent
    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)

        t = 0. ;
        hslope = 1. ;

	if(N2!=1) {
	  //2D problem, use full pi-wedge in theta
	  fractheta = 1.;
	}
	else{
	  //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
	  fractheta = 1.e-2;
	}

        fracphi = 1.;

        set_arrays() ;
        set_grid() ;

	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;
	fprintf(stderr,"rmin: %g\n",r) ;
	fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        /* output choices */
	tf = 2*Rout ;

	DTd = 1. ;	/* dumping frequency, in units of M */
	DTl = 50. ;	/* logfile frequency, in units of M */
	DTi = 50. ; 	/* image file frequ., in units of M */
        DTr = 1. ; /* restart file frequ., in units of M */
	DTr01 = 1000 ; 	/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
        rdump01_cnt = 0 ;
	defcon = 1. ;

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
	  coord(i,j,k,CENT,X) ;
	  bl_coord(X,&r,&th,&phi) ;

	  sth = sin(th) ;
	  cth = cos(th) ;

	  /* rho = 1.e-7*RHOMIN ; */
	  /* u = 1.e-7*UUMIN ; */

	  /* rho = pow(r,-4.)/BSQORHOMAX; */
	  /* u = pow(r,-4.*gam)/BSQOUMAX; */

	  rho = RHOMINLIMIT+(r/10./rhor)/pow(r,4)/BSQORHOMAX;
	  u = UUMINLIMIT+(r/10./rhor)/pow(r,4)/BSQORHOMAX;

	  /* these values are demonstrably physical
	     for all values of a and r */
	  /*
	    ur = -1./(r*r) ;
	    uh = 0. ;
	    up = 0. ;
	  */

	  ur = 0. ;
	  uh = 0. ;
	  up = 0. ;

	  /*
	    get_geometry(i,j,CENT,&geom) ;
	    ur = geom.gcon[0][1]/geom.gcon[0][0] ;
	    uh = geom.gcon[0][2]/geom.gcon[0][0] ;
	    up = geom.gcon[0][3]/geom.gcon[0][0] ;
	  */

	  p[i][j][k][RHO] = rho ;
	  p[i][j][k][UU] = u ;
	  p[i][j][k][U1] = ur ;
	  p[i][j][k][U2] = uh ;
	  p[i][j][k][U3] = up ;
	  p[i][j][k][B1] = 0. ;
	  p[i][j][k][B2] = 0. ;
	  p[i][j][k][B3] = 0. ;
	}

	rhomax = 1. ;
	fixup(p) ;
	bound_prim(p) ;

        //leave A[][] a 2D array for now, which means that magnetic field will be axisymmetric
	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2,0,0) A[i][j] = 0. ;
        ZSLOOP(0,N1,0,N2,0,0) {
#if(0)
                /* vertical field version */
                coord(i,j,k,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;
                A[i][j] = 0.5*pow(r*sin(th),2);
#elif(1)
                /* radial (monopolar) field version */
                coord(i,j,k,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;
                A[i][j] = (1-cos(th)) ;
#endif

        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,k,CENT,&geom) ;

		/* flux-ct */
		p[i][j][k][B1] = -(A[i][j] - A[i][j+1]
				+ A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*geom.g) ;
		p[i][j][k][B2] = (A[i][j] + A[i][j+1]
				- A[i+1][j] - A[i+1][j+1])/(2.*dx[1]*geom.g) ;

		p[i][j][k][B3] = 0. ;

		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

	/* enforce boundary conditions */
	fixup(p) ;
	bound_prim(p) ;
    
    



#if( DO_FONT_FIX )
	set_Katm();
#endif 


}


/* this version starts w/ BL 4-velocity and
 * converts to relative 4-velocities in modified
 * Kerr-Schild coordinates */

void coord_transform(double *pr,int ii, int jj, int kk)
{
  double X[NDIM],r,th,phi,ucon[NDIM],uconp[NDIM],trans[NDIM][NDIM],tmp[NDIM] ;
  double AA,BB,CC,discr ;
  double utconp[NDIM], dxdxp[NDIM][NDIM], dxpdx[NDIM][NDIM] ;
  struct of_geom geom ;
  struct of_state q ;
  int i,j,k,m ;

  coord(ii,jj,kk,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  blgset(ii,jj,kk,&geom) ;

  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  AA =     geom.gcov[TT][TT] ;
  BB = 2.*(geom.gcov[TT][1]*ucon[1] +
           geom.gcov[TT][2]*ucon[2] +
           geom.gcov[TT][3]*ucon[3]) ;
  CC = 1. +
          geom.gcov[1][1]*ucon[1]*ucon[1] +
          geom.gcov[2][2]*ucon[2]*ucon[2] +
          geom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(geom.gcov[1][2]*ucon[1]*ucon[2] +
          geom.gcov[1][3]*ucon[1]*ucon[3] +
          geom.gcov[2][3]*ucon[2]*ucon[3]) ;

  discr = BB*BB - 4.*AA*CC ;
  ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;
  /* now we've got ucon in BL coords */

  /* transform to Kerr-Schild */
  /* make transform matrix */
  DLOOP trans[j][k] = 0. ;
  DLOOPA trans[j][j] = 1. ;
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a) ;
  trans[3][1] = a/(r*r - 2.*r + a*a) ;

  /* transform ucon */
  DLOOPA tmp[j] = 0. ;
  DLOOP tmp[j] += trans[j][k]*ucon[k] ;
  DLOOPA ucon[j] = tmp[j] ;
  /* now we've got ucon in KS coords */

  /* transform to KS' coords */
  /* dr^\mu/dx^\nu jacobian, where x^\nu are internal coords */
  dxdxp_func(X, dxdxp);
  /* dx^\mu/dr^\nu jacobian */
  invert_matrix(dxdxp, dxpdx);
  
  for(i=0;i<NDIM;i++) {
    uconp[i] = 0;
    for(j=0;j<NDIM;j++){
      uconp[i] += dxpdx[i][j]*ucon[j];
    }
  }
  //old way of doing things for Gammie coords
  //ucon[1] *= (1./(r - R0)) ;
  //ucon[2] *= (1./(M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]))) ;
  //ucon[3] *= 1.; //!!!ATCH: no need to transform since will use phi = X[3]

  get_geometry(ii, jj, kk, CENT, &geom);
  
  /* now solve for relative 4-velocity that is used internally in the code:
   * we can use the same u^t because it didn't change under KS -> KS' */
  ucon_to_utcon(uconp,&geom,utconp);
  
  pr[U1] = utconp[1] ;
  pr[U2] = utconp[2] ;
  pr[U3] = utconp[3] ;

  /* done! */
}


