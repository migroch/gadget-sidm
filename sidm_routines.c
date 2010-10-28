#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

#define GSLWORKSIZE 100000

/*! \file sidm_routines.c
 *  \brief Extra fuctions and routines needed for the calculations of dark matter self interactions
 *
 *  This file contains the functions and routines necesary for the computation of 
 *  the self interactions probabilities and the velocity kicks due to the interactios.
 *  Written by Miguel Rocha, rocham@uci.edu. Oct 2010
 */

double g_geo(double r);
double geofactor_integ(double x, void * params);
double geofactor_angle_integ(double u, void * params);
double kernel(double u);

/*! This function calculates the interaction probability between two particles.
 *  It checks if comoving integration is on and does the necesary change of 
 *  variables and units.
 */
double prob_of_interaction(double r, FLOAT Vtarget[3], FLOAT Vno[3], int tbegin, int tend)
  {
    double dT, dloga, dV, dvx, dvy, dvz, prob,h,hubble_a;

    dvx = Vno[0]-Vtarget[0];
    dvy = Vno[1]-Vtarget[1];
    dvz = Vno[2]-Vtarget[2];
    dV = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

    if(All.ComovingIntegrationOn)
      {
        dloga = (tend - tbegin)*All.Timebase_interval;
        hubble_a = (All.Omega0 / (All.Time * All.Time * All.Time)
                    + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) + All.OmegaLambda);
        hubble_a = All.Hubble * sqrt(hubble_a);
        dT = dloga/hubble_a;

        dV = dV/All.Time; /* Convert from internal velocity p = a^2 dx/dt to peculiar velocity v = a dx/dt*/
        dV -= All.Time*hubble_a*r; /* Substract Hubble flow */

        h = All.ForceSoftening[1]*All.Time;
      }
    else
      {
        dT = (tend-tbegin)*All.Timebase_interval;
        h = All.ForceSoftening[1];
      }

    prob = All.InteractionCrossSection*dV*dT*g_geo(r)/(h*h*h*All.UnitLength_in_cm*All.UnitLength_in_cm);

    return prob;

  }

/*! This function returns the value of the geometrical factor needed for 
 *  the calculation of the interaction probability.
 */
double g_geo(double r)
{
  double  f, u;
  int i;

  u = r / 2.0 * GEOFACTOR_TABLE_LENGTH;
  i = (int) u;
  if (i >= GEOFACTOR_TABLE_LENGTH)
    i = GEOFACTOR_TABLE_LENGTH - 1;

  if (i <= 1)
    f = 0.992318  + (GeoFactorTable[0] - 0.992318)*u;
  else
    f = GeoFactorTable[i - 1] + (GeoFactorTable[i] - GeoFactorTable[i - 1]) * (u - i);

  return f;
}

/*! This routine sets the kicks for each particle after it has been decided that they will
 *  interact. It uses an algorithm tha conserves energy and momentum but picks a random direction
 *  so it does not conserves angular momentum.
 */
void  calculate_interact_kick(FLOAT Vtarget[3], FLOAT Vno[3], double kick_target[3], double kick_no[3])
{
  double dV,theta,phi,dvx,dvy,dvz;

  dvx = Vno[0]-Vtarget[0];
  dvy = Vno[1]-Vtarget[1];
  dvz = Vno[2]-Vtarget[2];
  dV = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

  theta = gsl_rng_uniform(random_generator)*M_PI;
  phi = gsl_rng_uniform(random_generator)*2.0*M_PI;

  kick_target[0] = (-Vtarget[0] + Vno[0] + dV*sin(theta)*cos(phi))/2.0;
  kick_no[0] = (Vtarget[0] - Vno[0] - dV*sin(theta)*cos(phi))/2.0;

  kick_target[1] = (-Vtarget[1] + Vno[1] + dV*sin(theta)*sin(phi))/2.0;
  kick_no[1] = (Vtarget[1] - Vno[1] - dV*sin(theta)*sin(phi))/2.0;

  kick_target[2] = (-Vtarget[2] + Vno[2] + dV*cos(theta))/2.0;
  kick_no[2] = (Vtarget[2] - Vno[2] - dV*cos(theta))/2.0;
}

/*! This routine initialazes the table that will be used to get the geometrical factor
 *  as a function of the two particle separations. It populates a table with the results of
 *  the numerical integration.
 */
void init_geofactor_table(void)
{
  int i;
  double result, abserr,r;
  gsl_function F;
  gsl_integration_workspace *workspace;
  workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);

  for(i = 0; i < GEOFACTOR_TABLE_LENGTH; i++)
  {
    r =  2.0/GEOFACTOR_TABLE_LENGTH * (i + 1);
    F.function = &geofactor_integ;
    F.params = &r;
    gsl_integration_qag(&F, 0.0, 1.0, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);
    GeoFactorTable[i] = 2*M_PI*result;
  }
  gsl_integration_workspace_free(workspace);
}

/*! This function returns the integrand of the numerical integration done on init_geofactor_table().
 */
double geofactor_integ(double x, void * params)
{
  double result, abserr,r;
  double newparams[2];

  r = *(double *) params;
  newparams[0] = r;
  newparams[1] = x;

  gsl_function F;
  gsl_integration_workspace *workspace;
  workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);

  F.function = &geofactor_angle_integ;
  F.params = newparams;

  gsl_integration_qag(&F, -1.0, 1.0, 0, 1.0e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41,workspace, &result, &abserr);

  return x*x*kernel(x)*result;
}

/*! This function returns the integrand of the angular part of the integral done on 
 *  init_geofactor_table().
 */
double geofactor_angle_integ(double u, void * params)
{
  double x,r,f;

  r = *(double *) params;
  x = *(double *) (params + sizeof(double));
  f = sqrt(x*x + r*r + 2*x*r*u);

  return kernel(f);
}

/*! This function returns the value W(x). The values of the density kernel as a funtion of x=r/h.
 */
double kernel(double u)
{
  double w;

  if (u < 1.0)
    {
      if (u <= 0.5)
        {
          w = KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u;
          return w;
        }
      else
        {
          w = KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
          return w;
        }
    }
  else
  {
    return 0.0;
  }
}