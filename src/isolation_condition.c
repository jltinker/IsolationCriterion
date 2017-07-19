#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define OMEGA_M 0.3
#define PI 3.141592741
#define RHO_CRIT 2.775E+11
#define DELTA_HALO 200
#define SPEED_OF_LIGHT 3.0E+5
#define c_on_H0 2997.92
#define BIG_G 4.304E-9 /* BIG G in units of (km/s)^2*Mpc/M_sol */
#define HUBBLE 0.7
#define RT2PI 2.50663

/* External functions
 */
int filesize(FILE *fp);
FILE *openfile(char *ff);
float qromo(float (*func)(float), float a, float b,
	     float (*choose)(float(*)(float), float, float, int));
float midpnt(float (*func)(float), float a, float b, int n);

/* Internal functions
 */
float find_satellites(int i, int ngal, float *ra, float *dec, float *redshift, 
		      float *ang_rad, float *mgal, float *psat, int *ihost, float *mhalo, float *rhalo, float *sigma);
float stellar2halo(float mgal, float redshift, float *mhalo, float *rhalo, float *sigma);
float angular_separation(float a1, float d1, float a2, float d2);
float distance_redshift(float z);
float func_dr1(float z);
float angular_probability(float mass, float dr, float rad, float ang_rad);

/* Global variables
 */
int PHOTOZ=0; // don't use the redshift information, only projected separation
float PHOTOZ_ERROR = 0;

int main(int argc, char **argv)
{
  float *ra, *dec, *z, *mgal, *theta, *psat, *nsat, *mhalo, *rhalo, x1, *sigma;
  int ngal, i, *ihost;
  FILE *fp;
  char aa[1000];

  if(argc<2)
    {
      fprintf(stderr,"isolation_condition filename [PHOTOZ_ERROR] > output\n");
      exit(0);
    }
  fp = openfile(argv[1]);
  ngal = filesize(fp);
  PHOTOZ = 0;
  if(argc>2) 
    {
      PHOTOZ = 1;
      PHOTOZ_ERROR = atof(argv[2]);
      fprintf(stderr,"PHOTOZ_ERROR= (1+z)*%f\n",PHOTOZ_ERROR);
    }

  // declare storage for the galaxies
  ra = vector(1,ngal);
  dec = vector(1,ngal);
  z = vector(1,ngal);
  mgal = vector(1,ngal);
  mhalo = vector(1,ngal);
  sigma = vector(1,ngal);
  rhalo = vector(1,ngal);
  theta = vector(1,ngal);
  psat = vector(1,ngal);
  nsat = vector(1,ngal);
  ihost = ivector(1,ngal);

  // read in the galaxies
  for(i=1;i<=ngal;++i)
    {
      fscanf(fp,"%f %f %f %f",&ra[i], &dec[i], &z[i], &mgal[i]);
      fgets(aa,1000,fp);
      // convert to radians
      ra[i] *= PI/180.;
      dec[i] *= PI/180.;
      psat[i] = 0;
      ihost[i] = -1;
      // check if mgal in right units
      if(mgal[i]>100)mgal[i] = log10(mgal[i]);
    }
  fclose(fp);
  fprintf(stderr,"Done reading %d lines from [%s]\n",ngal,argv[1]);

  // initialize the SHMR conversion tables
  stellar2halo(-1,-1,&x1, &x1, &x1);

  // get the angular radii for all galaxies' halos
  for(i=1;i<=ngal;++i)
    theta[i] = stellar2halo(mgal[i], z[i], &mhalo[i], &rhalo[i], &sigma[i]);

  for(i=1;i<=-ngal;++i)
    {
      printf("BOO %f %f\n", mgal[i], log10(mhalo[i]));
    }

  // loop through each galaxy in the sample, find possible satellites.
  for(i=1;i<=ngal;++i)
    nsat[i] = find_satellites(i,ngal,ra,dec,z,theta,mgal,psat,ihost,mhalo,rhalo,sigma);

  //output results
  for(i=1;i<=ngal;++i)
    {
      printf("%e %d %e %f %e\n",psat[i],ihost[i],mgal[i], theta[i]*180/PI, mhalo[i]);
    }

}


/* This converts galaxy stellar mass to halo mass
 * Using the tabulated stellar-to-halo mass relations of Behroozi+13
 */
float stellar2halo(float mgal, float redshift, float *mhalo, float *rhalo, float *sigma)
{
  static int flag = 1, nx[9];
  static float zz[9], *mgx[9], *mhx[9];
  FILE *fp;
  int i, j, k, iz;
  float mh1, mh2, mh, r, r200, theta;
  char aa[1000];

  if(flag)
    {
      flag = 0;

      zz[0] = 0.1;
      for(i=1;i<=8;++i)
	zz[i] = i;

      for(i=0;i<=8;++i)
	{
	  if(i==0)
	    fp = openfile("shmr_z0.10.dat");
	  else
	    {
	      sprintf(aa,"shmr_z%d.00.dat",i);
	      fp = openfile(aa);
	    }
	  nx[i] = filesize(fp);

	  mgx[i] = vector(1,nx[i]);
	  mhx[i] = vector(1,nx[i]);

	  // get rid of header
	  fgets(aa,1000,fp);
	  for(j=1;j<=nx[i];++j)
	    {
	      fscanf(fp,"%f %f",&mhx[i][j],&mgx[i][j]);
	      fgets(aa,1000,fp);
	      mgx[i][j] = mhx[i][j] + mgx[i][j]; // tables tabulate lg(mg/mh) ratio
	    }
	  fclose(fp);
	}
      fprintf(stderr,"Done initializing SHMR lookup table\n");
      if(mgal<0)return 0;
    }

  // find the redshift bin of the gal
  for(iz=0;iz<9;++iz)
    if(zz[iz]>redshift)break;
  //fprintf(stderr,"%d\n",iz);
  if(iz>0)iz--; // lower redshfit limit
  if(iz>8)iz=8;

  // get the halo mass for lower bin
  for(i=2;i<=nx[iz];++i)
    if(mgx[iz][i]>mgal)break;
  if(i>=nx[iz])i--;
  mh1 = (mgal-mgx[iz][i-1])*(mhx[iz][i]-mhx[iz][i-1])/(mgx[iz][i]-mgx[iz][i-1]) + mhx[iz][i-1];
  //printf("%d %f\n",iz,mh1);
  
  // fprintf(stderr,"%d %f %f\n",i,mgal,mgx[iz][i]);

  // get the halo mass for upper bin
  iz++;
  for(i=2;i<=nx[iz];++i)
    if(mgx[iz][i]>mgal)break;
  if(i>=nx[iz])i--;
  mh2 = (mgal-mgx[iz][i-1])*(mhx[iz][i]-mhx[iz][i-1])/(mgx[iz][i]-mgx[iz][i-1]) + mhx[iz][i-1];
  // printf("%d %f %d %d\n",iz,mh2,i,nx[iz]);


  // interpolate in linear redshift
  mh = (redshift-zz[iz-1])*(mh2-mh1)/(zz[iz]-zz[iz-1]) + mh1;
  if(redshift<zz[0])mh = mh1;

  // convert redshift to distiance
  r = distance_redshift(redshift);
  // convert mass to radius
  mh = pow(10.0,mh);
  r200 = pow(3*mh*HUBBLE/(4*PI*DELTA_HALO*OMEGA_M*RHO_CRIT),0.3333);
  // convert R200 to radians
  theta = r200/r;

  // return halo properties, too
  *mhalo = mh;
  *rhalo = r200;
  *sigma = sqrt(BIG_G*mh/2.0/r200*(1+redshift));

  return theta;

}

/* function(s) for converting redshift to distance
 */
float func_dr1(float z)
{
  return pow(OMEGA_M*(1+z)*(1+z)*(1+z)+(1-OMEGA_M),-0.5);
}
float distance_redshift(float z)
{
  float x;
  if(z<=0)return 0;
  x= c_on_H0*qromo(func_dr1,0.0,z,midpnt);
  return x;
}

/* angular separation between two points in ra/dec [radians]
 */
float angular_separation(float a1, float d1, float a2, float d2)
{
  float cd1,cd2,sd1,sd2,ca1a2,sa1a2;

  return atan((sqrt(cos(d2)*cos(d2)*sin(a2-a1)*sin(a2-a1) + 
		    pow(cos(d1)*sin(d2) - sin(d1)*cos(d2)*cos(a2-a1),2.0)))/
	      (sin(d1)*sin(d2) + cos(d1)*cos(d2)*cos(a2-a1)));
}

/* loop through and find satellites for the given galaxy i
 *
 * NB, i know this is a terribly brute force approach. Feel free to re-write this code if you
 * want to use some fancy KD tree or something. But it runs on typical samples in not a lot of
 * time, and I know it works.
 */
float find_satellites(int i, int ngal, float *ra, float *dec, float *redshift, 
		      float *ang_rad, float *mgal, float *psat, int *ihost, float *mhalo, float *rhalo, float *sigma)
{
  int j, j1;
  float dx, dy, dz, theta, prob_ang, vol_corr, prob_rad, grp_lum, p0, theta_max, x1, photoz_error;

  theta_max = ang_rad[i];
  for(j=1;j<=ngal;++j)
    {
      //if self, skip
      if(j==i)continue;

      //if the galaxy i is smaller than galaxy j, skip it
      if(mgal[i]<mgal[j])continue;

      //if already in group, skip
      //if(psat[j]>0.5)continue;

      /*
      dx = fabs((ra[i]-ra[j])*cos(dec[i]));
      if(dx>3.0*theta_max)continue;
      dy = fabs((dec[i]-dec[j]));
      if(dy>3.0*theta_max)continue;
      */
      if(!PHOTOZ)
	{
	  dz = fabs(redshift[i] - redshift[j])*SPEED_OF_LIGHT;
	  if(dz>6*sigma[i])continue;
	}
      else
	{
	  dz = fabs(redshift[i] - redshift[j])*SPEED_OF_LIGHT;
	  photoz_error = (1+redshift[i])*PHOTOZ_ERROR*1.41*SPEED_OF_LIGHT;
	  if(dz>3*photoz_error)continue;
	}
      theta = angular_separation(ra[i],dec[i],ra[j],dec[j]);
      if(j==5 && i==-1) { fprintf(stdout,"%f %f %f\n",theta,theta_max,theta/theta_max); }
      if(theta>theta_max)continue;
      
      prob_ang = angular_probability(mhalo[i],theta,rhalo[i],theta_max);
      if(j==5 && i==-1) {
	fprintf(stdout,"ACK here %d %e %e %e %e %e %e %e\n",j, dz/sigma[i],
		theta/theta_max,prob_ang, prob_rad, prob_ang*prob_rad, sigma[i], p0);
	fflush(stdout);
	exit(0);
      }

      if(!PHOTOZ)
	prob_rad = exp(-dz*dz/(2*sigma[i]*sigma[i]))/(RT2PI*sigma[i])*SPEED_OF_LIGHT;
      else 
	prob_rad = exp(-dz*dz/(2*photoz_error*photoz_error))/(RT2PI*photoz_error)*SPEED_OF_LIGHT;
      
      p0 = (1 - 1/(1+prob_ang*prob_rad/10));   
      //p0 = (1 - 1/(1+prob_ang*prob_rad/120));   
      if(p0>psat[j]) { psat[j]=p0; ihost[j] = i; }
      //fprintf(stdout,"ACK here %d %e %e %e %e %e %e %eq\n",j, dz/sigma[i],
      //	      theta/theta_max,prob_ang, prob_rad, prob_ang*prob_rad, sigma[i], p0);

    }
  
}

/* This is the projected density of an NFW profile. 
 */
float angular_probability(float mass, float dr, float rad, float ang_rad)
{
  float c, x, rs, delta, f;

  dr = dr*rad/ang_rad;

  c = 10.0*pow(mass/1.0E+14,-0.11);
  rs = rad/c;
  x = dr/rs;

  if(x<1)
    f = 1/(x*x-1)*(1-log((1+sqrt(1-x*x))/x)/(sqrt(1-x*x)));
  if(x==1)
    f = 1.0/3.0;
  if(x>1)
    f = 1/(x*x-1)*(1-atan(sqrt(x*x-1))/sqrt(x*x-1));

  delta = DELTA_HALO/3.0*c*c*c/(log(1+c)-c/(1+c));

  return 1.0/c_on_H0*2*rs*delta*f;
}


