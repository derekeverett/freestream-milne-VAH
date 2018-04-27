//trigTable is a table with 10 rows for ten combinations or p^(mu)p_(nu) normalized by the energy
#pragma once
#include <stdio.h>
//#include <math.h>
#include "Parameter.h"
#define PI 3.141592654f

void calculateBulkInvReynolds(float *P_T, float *residualBulk, float *R_Pi_Inv, parameters params)
{
  int DIM = params.DIM;
  float TAU = params.TAU;
  #pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {
    R_Pi_Inv[is] = residualBulk[is] / P_T[is];
  }
}

void calculateShearInvReynolds(float *energyDensity, float *P_L, float *P_T, float **residualShear, float *R_pimunu_Inv, parameters params)
{
  int DIM = params.DIM;
  float TAU = params.TAU;
  float tau2 = TAU*TAU;
  #pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {
    float num = residualShear[0][is]*residualShear[0][is] - 2.0 * (residualShear[1][is]*residualShear[1][is] + residualShear[2][is]*residualShear[2][is] + tau2*residualShear[3][is]*residualShear[3][is])
    + (residualShear[4][is]*residualShear[4][is] + residualShear[7][is]*residualShear[7][is] + tau2*tau2*residualShear[9][is]*residualShear[9][is])
    + 2.0 * (residualShear[5][is]*residualShear[5][is] + tau2*residualShear[6][is]*residualShear[6][is] + tau2*residualShear[8][is]*residualShear[8][is]);

    float den = energyDensity[is]*energyDensity[is] + P_L[is]*P_L[is] + 2.0*P_T[is]*P_T[is];
    R_pimunu_Inv[is] = sqrt(num) / sqrt(den);
  }
}
