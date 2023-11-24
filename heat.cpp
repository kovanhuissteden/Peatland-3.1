/***************************************************************************
                   heat.h  -  soil temperature functions PEATLAND
                             -------------------
    begin                : |02-01-2003|
    copyright            : (C) |2003| by |J. van Huissteden|
    email                : |ko.van.huissteden@geo.falw.vu.nl|
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/***************************************************************************
  MODIFICATIONS

  May 2003

  Adaptions for correct modelling of soil freezing

  17 December 2003

  Bug fixes in correction of surface temperature for snow depth
  (These did not affect previous research results)
    Correction for NaN surface temperatures when MaxSnowDepth is set to 0.0
    Correction for bug in reading snowdepths from file
    Surface temperature correction was not done with snowdepths from file

  31 May 2006
  Fixed LatentHeat value replaced by polynomial approximation of temperature
  dependent LatentHeat

  5 September 2008
  Calculation of thermal conductivity from soil constituents adapted to
  Hillel, fixed constants from formula by Mueller removed
  
  June 2011
  Calculation of thermal conductivity cf Balland and Arp, J. Environ. Eng. Sci. 4: 549–558 (2005)
 ***************************************************************************/


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <climits>
#include <cmath>
using namespace std;
#include "matrix.h"
#include "general.h"
#include "initialize.h"
#include "heat.h"

void HeatVar(const int steps, const double tsurf)
/* Solves the heat equation for heat transport into the soil using the classical explicit scheme;
  steps is the number of time steps, tsurf the surface temperature                             */

{
  int i = 1, j;
  double *u, *ut, *r, *rj, *uj, *utj;
  Matrix lm1, lm2, diff, rho;

/* heat capacity: is interpolated from saturated soil heat capacity and dry soil heat capacity
   the thermal conductivity is interpolated using a parabolic function to account for slower rise of tc
   near saturation (see Luckner & Schestakow) */

  ThermalProperties(SoilTemp);                      // Update soil thermal properties for groundwater table changes
  lm1 = Layers.Col(2);                              // auxilary array for interpolation
  lm2 = Layers.Col(2);
  lm2(1, 1) = -0.0;                                 // interpolate thermal diffusivity for heat model layers from soil model layers, assuming that thermal properties below the soil profile are the same as those at its base
  lm2(NrLayers, 1) = -1000.0;
  diff.Resize(NrHeatLayers);
  interp(lm2, ThermDiffVar, HeatLayers, diff, TRUE);  // interpolate diffusivity array
  rho = diff * (TStepHeat / (DStepHeat * DStepHeat));  // initialize numerical solution arrays (classical explicit scheme)
  if (rho.Max() > 0.5)                              // check for stability of the solution
  {
    cout << HEAT_ERROR1 << endl;
    exit(EXIT_FAILURE);
  }
  r = rho.Data();
  u = new double[NrHeatLayers + 2];                 // initialize u and ut solution arrays ofthe heat equation
  *u = SnowVegCorrection(tsurf);                       // first element is upper boundary condition (surface temperature), last is bottom boundary condition (average yearly temperature)
  //cout << HeatCondTop << " " << SnowDepth << endl;
  memcpy(u + 1, TProfile.Data(), NrHeatLayers*sizeof(double));
  *(u + NrHeatLayers + 1) = T_average;
  ut = new double[NrHeatLayers];
  for (i = 1; i <=steps; i++)                       //solve PDE for temperature profile: classic explicit solution
  {                                                 // loop time steps
    rj = r;
    uj = u + 1;
    utj = ut;
    for (j = 1; j <= NrHeatLayers; j++)             // loop depth steps
    {
      *utj = *rj * *(uj + 1) + (1 - 2 * *rj) * *uj + *rj * *(uj - 1);
      rj++;
      uj++;
      utj++;
    }
    memcpy(u + 1, ut, NrHeatLayers*sizeof(double));
    // adjust thermal properties for ice content and latent heat as soon as temperatures anywhere get below 0
    TProfile.PutData(1, NrHeatLayers, u + 1);
    if (TProfile.Min() <= 0)                        // revise thermal properties for freezing
    {
      interp(HeatLayers, TProfile, lm1, SoilTemp, TRUE);
      ThermalProperties(SoilTemp);
      interp(lm2, ThermDiffVar, HeatLayers, diff, TRUE);  // interpolate diffusivity array
    }
  }
  interp(HeatLayers, TProfile, lm1, SoilTemp, TRUE);

  // if (SoilTemp(5) > 0.0) ThermDiffVar.Disp();
}  // end HeatVar


void HeatSimple(double time)
/* Analytical solution of the heat equation based on constant thermal diffusivity and sinusoidal temperature time series
   cf Animo - see Groenendijk and Kroes, 1997                                                             */
{
  int i;
  double omega, dm, d;

  omega = 2* PI / YEAR;
  dm = sqrt(2 * ThermDiff / omega);
  for (i = 1; i <= NrLayers; i++)
  {
    d = -Layers(i, 2);
    SoilTemp(i) = T_average + T_amplitude * exp(-d / dm) * cos(omega * time + PI - d / dm);
  }
} // end HeatSimple

void Temperature()
/* computes temperature soil profile dependent on temperature profile
   t is time step nr                                                          */
{
  switch (ThermModel)
  {
    case 0: {HeatVar((int)(Timestep/TStepHeat), SoilTData(StepNr)); break;} // temperature layer-dependent thermal properties
    case 1: {HeatSimple((StepNr - 0.5) * Timestep ); break;}                // temperature profile simple sinusoidal model
    case 2: SoilTemp = SoilTData.Row(StepNr);                               // temperature from observational data
  }
  if (ProfileOutput.Contains(1)) SoilTemp.Write(output1);                           // write to output file if requested
  if (ProfileOutput.Contains(6)) Ice.Write(output6);                                // write ice content to output
}


void ThermalProperties(Matrix T)
/* computes heat capacity, thermal conductivity and diffusivity based on soil constituents
   approach cf Hillel (1980)
   T: temperature of each soil layer                                                           */
// NB: Check thermal property calculation: for a peat soil the resulting thermal diffusivity is a factor 3 too high
{
  int i, a;
  double fquartz, fi, fw, fo, fm, fa, s, w, lh, fq, Forg, rhop, rhob, sat, fsolid, Ksolid, Kdry, Ksat, Ke, Ksoil, fp, waterfraction;
  double ac = 0.053, alfa = 0.24, beta = 18.3;
  Matrix tp(NrLayers, 2);                               // thermal properties, 1st column conductivity, 2nd column heat capacity

  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i, 4);                              // soil profile horizon number
    // volume fractions of soil constituents
	waterfraction = MoistTheta(i);
	fsolid = 1 - PoreVol(i);
	Forg = PercOrg(a) / 100;
	fquartz = SandFraction(a);
    fi = Ice(i);							// ice volume fraction
    fo = DBD(a) * Forg / DensOrg;         // organic matter volume
    fm = fsolid - fo;						// volume mineral matter
    fw = waterfraction - fi;				// water volume
	fp = PoreVol(i);							// pore volume 
    fa = fp - waterfraction;                    // air volume
	fq = fquartz * fm;				// correct for quartz sand content assuming equal mineral densities
	fm = (1 - fquartz) * fm;
	// thermal conductivity cf Balland & Arp (2005)
	rhob = DBD(a);  // bulk density and density solids, equation 6
	rhop = 1 / (Forg / DensOrg + (1 - Forg) / DensMin); 
	sat = 1.0 - Saturation(i);  // saturation  with water and ice, equation 6
	Ksolid = pow(CondOrg, fo) * pow(CondQuartz, fq) * pow(CondMiner, (1 - fo - fq));  // conductivity of solids, equation 15
	Kdry = ((ac * Ksolid - CondAir) * rhob + CondAir * rhob) / (rhop - (1 - ac) * rhob); // dry conductivity, eq 16
	if (fi > 0.0) // equation 17 and 18 for the Kersten number, 12 and 13 for Ksat depending on presence of ice
	{
		Ksat = pow(Ksolid, (1 - fp)) * pow(CondIce, fi) * pow(CondWater, fp - fi);
		Ke = pow(sat, (1 + fo));
	} else {
		Ksat = pow(Ksolid, (1 - fp)) * pow(CondWater, fp);
		Ke = pow(sat, 0.5 * (1 + fo - alfa * fq - fm)) * pow((pow((1 / (1 + exp(-beta * sat))), 3) - pow(((1 - sat) / 2), 3)), (1 - fo));
	}
	Ksoil = (Ksat - Kdry) * Ke + Kdry;
	tp(i, 1) = Ksoil;
	
    if (i == 1) HeatCondTop = tp(1, 1);
    tp(i, 2) = fo * HCOrg + fm * HCMiner + fi * HCIce + fw * HCWater + fa * HCAir;  // heat capacity
    if (T(i) <= 0)
    {
      s = Freeze(i, 2) / pow((Freeze(i, 3) - T(i)), Freeze(i, 2) + 1);   // slope freezing curve at T; NOTE: this is ice mass per kg dry soil
      lh = LatentHeat(1)*T(i)*T(i) + LatentHeat(2)*T(i) + LatentHeat(3); // determine latent heat lh from temperature using quadratic polynomial approximation
      tp(i, 2) = tp(i, 2) + lh * s * DBD(a);     // add apparent heat capacity from freeze/thaw; multiply with DBD because s is ice mass / kf dry soil
      UnFrozen(i) = Freeze(i, 1) + 1 / pow((Freeze(i, 3) - T(i)), Freeze(i, 2));  // Update unfrozen water conten
      w = MoistTheta(i) - UnFrozen(i) * DBD(a) / DensWater;  // Update ice content: if there is more water than the unfrozen water content according to the unfrozen water function, it is ice
      if (w > 0) Ice(i) = w; else Ice(i) = 0.0;
    } else Ice(i) = 0.0;
    ThermDiffVar(i) = 86400 * tp(i, 1) / tp(i, 2);           // thermal diffusivity in m2d-1, calculated by dividing conductivity by volumetric heat capacity
  }
}  // end ThermalProperties


double SnowVegCorrection(const double tsurf)
/* Corrects soil surface temperature in case of snow cover presence
and also corrects for vegetation effect when there is no snow */
{

  double daycount, tcorrect, a, vcorrect;

  tcorrect = tsurf;
  if (tsurf > 0.0) SnowStartDay = DayOfTheYear;              // this registrates the last day with above zero temperatures, indicating the day at which snowfall may have started
  if (strlen(SnowFile) != 0)                                 // SnowHeight from time series
  {
    SnowDepth = SnowData(StepNr);
    if (SnowDepth > 0.0)                                     // compute snowcover temperature correction
    {
      if (tsurf <= 0.0)                                      // below zero temperatures: isolation effect
      {
        a = (CondSnow * 0.5 * DStepHeat) / (SnowDepth * HeatCondTop);
        tcorrect = (TProfile(1) + a * tsurf) / (1 + a);
        if (tcorrect > 0.0) tcorrect = 0.0;                   // prevent above zero temperatures in early winter
      } else {
        tcorrect = 0.0;                                       // above zero temperature: melting snow causes surface temperatureto be equal to zero
      }
    } 
  } else {                                                    // snowheight from simple linear accumulation -melt model
    if (MaxSnowdepth > 0.0)                                   // skip correction if maximum snowdepth is zero
    {
      if (tsurf <= 0.0)                                        // compute snowcover at below zero
      {
        daycount = (YEAR - SnowStartDay + DayMaxSnowdepth);
        if (DayOfTheYear > 182) SnowDepth = MaxSnowdepth * (DayOfTheYear - SnowStartDay) / daycount;    // snowcover in last months of the year
        if (DayOfTheYear < DayMaxSnowdepth) SnowDepth = MaxSnowdepth * (YEAR - SnowStartDay + DayOfTheYear) / daycount;  // snowcover before maximimum day in first months of the year
        if ((DayOfTheYear >= DayMaxSnowdepth) & (DayOfTheYear < 182)) SnowDepth = MaxSnowdepth;    // after the maximum snowcover day
/* compute correct soil surface temperature cf Granberg et al., Water Resources Research 35:3771-3782 */
        a = (CondSnow * 0.5 * DStepHeat) / (SnowDepth * HeatCondTop);
        tcorrect = (TProfile(1) + a * tsurf) / (1 + a);
        if (tcorrect > 0.0) tcorrect = 0.0;                     // prevent above zero temperatures in early winter
      } else {                                                  // check melting snow
        if (SnowDepth > 0.0)                                    // melting snow present
        {
          tcorrect = 0.0;                                       // surface temperature is zero
          SnowDepth = SnowDepth - SnowMeltrate * tsurf * Timestep; // melting
          if (SnowDepth < 0.0) SnowDepth = 0.0;                 // all snow has disappeared
        }
      }
    }
  }

  //if (SnowDepth == 0.0) tcorrect = VegTScalingFactor * tsurf;			// correction for vegetation when no snow is present
  
  if (SnowDepth == 0.0) {
      vcorrect = VegTScalingFactor * CurrentLAI;
      if (vcorrect <= tcorrect) tcorrect = tcorrect - vcorrect;			// LAI dependent correction for vegetation when no snow is present
  }

  return tcorrect;
} // end SnowCorrection
