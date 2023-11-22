/***************************************************************************
             water.cpp  -  water submodel functions of PEATLAND
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

  September 2004

  Bug corrected: when using pF curves defined by Van Genuchten parameters
  the soil moisture values in MoistTheta were not calculated from pFCurves
  (containing the calculated pF curve) but from the Van Genuchten parameters
  themselves.

  November/december 2007

  Addition of a correction to put lauers that are saturated by a fraction
  < 0.0001 to complete saturation
  to correct small pf curve/poristy errors and to prevent small negative saturations

  Addition of a time lag for anaerobe development and methane production
  after rapid saturation of a layer
  
  April 2011
  Van Genuchten function adapted fro general use and moved to this file
  Calculation of water table using a bucket model added
  
  May 2011 water table model improved with frozen soil
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
#include "heat.h"
#include "water.h"

double  VanGenuchten(double pF, double theta_r, double theta_s, double alfa, double n, int convert)
/* Expands van Genuchten parameters to complete pF curve
pF is the suction values for which the moisture content theta is computed
theta_r, theta_s, alfa, n, l are the parameters of the Van Genuchten function
Note: these are stored in the Layer_pF matrix  read from the parameter file
theta contain is the output  
convert indicates whether pF is given as true pF (convert > 0) or cm (convert <= 0)                                                      

Note: l (last value in rows of Layer_pF) is not used for the water retention curve,
only for conductivity acc to Van Genuchten and will be implemented later */

{
	double theta, h;

	if (convert > 0) h = pow(10, pF); else h = pF;
	theta = theta_r + (theta_s - theta_r) / pow((1 + pow(fabs(alfa * h), n)),(1-1/n));
      // here the formula as cited by Wosten et 1994 has been used (Staringreeks)
    return(theta);
}     // end of function VanGenuchten


double FindWaterBase()
/* returns the base of the water level movement in the profile,
which is under frozen conditions the current water table
or the lowest possibele water table (MinGW) under non-frozen conditions */

{
	int i, indfrozen;
	double frozentop, base;
	
	base = MinGW;  
	i = 1; // find the index of the frozen layer
	while (SoilTemp(i) >= 0.0)
	{
		i++;
		if (i > NrLayers) break;
	}                    
	indfrozen = i;
	if (indfrozen <= NrLayers)
	{
		// if a frozen layer occurs within the profile interpolate its top from the temperature gradient
		if (indfrozen > 1)
		{
			frozentop = -((indfrozen - 0.5) * LayerThickness + LayerThickness * SoilTemp(indfrozen) / (SoilTemp(indfrozen - 1) - SoilTemp(indfrozen)));
		} else
		{
			if (TData(StepNr) > 0.0)
			{ 
				frozentop = -(0.5 * LayerThickness +  (0.5 * LayerThickness) * SoilTemp(indfrozen) / (TData(StepNr) - SoilTemp(indfrozen)));
			} else frozentop = 0.0;
		}
		if (frozentop >= MinGW) base = frozentop;
	    // assuming that if the water table is situated below frozen top level no water can be added
		CurrentFrozen = frozentop;  // registrate frozen layers
	} else CurrentFrozen = FrozenWatertable = -MaxDepthHeat;
	return (base);
	// return(MinGW);
}

double Drainage(double watertable)
/* calculates drainage or seepage from distance to drainage line and hydraulic conductivity
   Drainage term is negative when seepage occurs, positive when drainage occurs; drainage in m/day*/
{
    double head, W; // water table head and drainage resistance
    
    // determine drainage water level and head
    if (strlen(DrainageFile) != 0) DrainLevel = DrainData(StepNr);
    head = watertable - DrainLevel; // negative head for seepage, positive for drainage
    // determine DrainageDepth
    DrainageDepth = Horizons(NrHorizons);
    // reduce drainage depth to unfrozen part of the soil
    if ((-CurrentFrozen) < DrainageDepth) DrainageDepth = -CurrentFrozen;
    W = DrainageDist * DrainageDist / (Ksat * DrainageDepth);
    return(head / W);
}

int belowWT(double currentbase, double watertable)
/* calculates the amount of water below the water table to the base of water table fluctuation range
currentbase is the current base of the soil water profile, either permafrost top or lowest water table */
{
	int i = 1, a;
	int basehorizon, tophorizon;
	double top, base, wt;
	
	// find horizon in whicb the base of the water column is situated
	if (watertable > 0.0) wt = 0.0; else wt = 0.0 - watertable; // exclude ponded water, to compare with horizon depths the watertable should be positive
	while (-currentbase > Horizons(i))
	{
		i++;
		if (i > NrHorizons) break;
	}
	basehorizon = i;
	
	i = 1;	// find horizon in whicb the current water table is situated
	while (wt > Horizons(i))
	{
		i++;
		if (i > NrHorizons) break;
	}
	tophorizon = i;
	BelowWTStorage = 0.0;    
    //cout << tophorizon << " " << basehorizon << " " << currentbase << " " << wt << endl;
	if (tophorizon <= basehorizon) // compute for each horizon the amount of water stored and add to saturated storage
	{
		for (a = tophorizon; a <= basehorizon; a++)
		{	
			if (a == tophorizon) top = wt; else top = Horizons(a - 1); // top horizon or first horizon has always the current water table at top of the water column
			if (a == basehorizon) base = -currentbase; else base = Horizons(a);
			BelowWTStorage = BelowWTStorage + Porosity(a) * (base - top);
		}
	} else if (currentbase == MinGW) cout <<  WT_ERROR1 << endl;
	return tophorizon;
}

double CalcFrozenStorage(double watertable)
/* Calculates the amount of water that is stored in frozen layers above the minimum water table */
{

	int i = 1, a;
	int basehorizon, tophorizon, wthorizon;
	double currenttop, top, base, wt, frozen;
	double zstep = 1.0e-4; // integration step, millimeters
	double storage = 0.0;
	double z, theta, dz, zfinal, h;
	int b, l;
	
	
	
	if (watertable > 0.0) wt = 0.0; else wt = 0.0 - watertable; // exclude ponded water, to compare with horizon depths the watertable should be positive
	// establish the presence of a frozen water table
	if (FrozenWatertable == (-MaxDepthHeat))
	{
		if (watertable < CurrentFrozen) FrozenWatertable =-wt; else FrozenWatertable = CurrentFrozen;
	}
	wt = -FrozenWatertable;
	
	while (-MinGW > Horizons(i))
	{
		i++;
		if (i > NrHorizons) break;
	}
	basehorizon = i;
	i = 1;	// find horizon in whicb the current water table is situated
	while (wt > Horizons(i))
	{
		i++;
		if (i > NrHorizons) break;
	}
	wthorizon = i;
	currenttop = -CurrentFrozen;
	i = 1;	// find horizon in whicb the frozen top or watertable is situated
	while (currenttop > Horizons(i))
	{
		i++;
		if (i > NrHorizons) break;
	}
	tophorizon = i;
	frozen = 0.0;    
	if (tophorizon <= basehorizon) // compute for each horizon the amount of water stored and add to saturated storage
	{
		for (a = tophorizon; a <= basehorizon; a++)
		{	
			if (a == tophorizon) top = wt; else top = Horizons(a - 1); // top horizon or first horizon has always the current water table at top of the water column
			if (a == basehorizon) base = -MinGW; else base = Horizons(a);
			frozen = frozen + Porosity(a) * (base - top);
		}
	}
// Above water table frozen storage
	l = pFVal.Length();
	if (FrozenWatertable < CurrentFrozen)  // above frozen water table water only calculated if the watertable is below zero!
	{
		for (a = tophorizon; a <= wthorizon; a++)  // integrate for every horizon seperately depending on pF curve
		{
			if (a == 1) top = 0.0; else top = Horizons(a - 1); // set top and base for integration
			if (a == wthorizon) base = wt; else base = Horizons(a);
			z = wt - top - 0.5 * zstep; // z is the position with respect to water table
			zfinal = wt - base;
			while (z >= zfinal) // integration of theta over z, start from base of horizon
			{
/*	Applying the Van Genuchten equation makes the model very slow because of the repeated integration; better use interpolation
theta is computed here by interpolating a pF curve at 1 cm suction intervals, computed in paramcheck.cpp */
				h = 100 * z;
				b = (int)(floor(h));
				theta = pFCurves(a, b + 1) - (h - b) * (pFCurves(a, b + 1) - pFCurves(a, b + 2));
				storage += theta * zstep; // add water in zstep to storage
				z -= zstep;
			}
			// last step, if a small difference < zstep remains
			z += zstep;
			if (z > zfinal)
			{
				dz = z - zfinal;
				z = z - 0.5 * dz;
				h = 100 * z;
				b = (int)(floor(h));
				theta = pFCurves(a, b + 1) - (h - b) * (pFCurves(a, b + 1) - pFCurves(a, b + 2));
				storage += theta * dz;
			}
		}
	}
	frozen = frozen + storage;
	return frozen;
}

double aboveWT(int wthorizon,  double watertable)
/* calculates the above water table unsaturated storage
tophorizon is the index to the horizon which contains the water table*/
{
	double zstep = 1.0e-4; // integration step, meters
	double storage = 0.0;
	double base, top, z, theta, wt, dz, zfinal, h;
	int a, b, l;
	
	if (watertable > 0.0) wt = 0.0; else wt = 0.0 - watertable; // exclude ponded water, to compare with horizon depths the watertable should be positive
	l = pFVal.Length();
	//cout << "   " << wt << " " << wthorizon << endl;
	if (wt > 0.0)  // above water table water only calculated if the watertable is below zero!
	{
		for (a = 1; a <= wthorizon; a++)  // integrate for every horizon seperately depending on pF curve
		{
			if (a == 1) top = 0.0; else top = Horizons(a - 1); // set top and base for integration
			if (a == wthorizon) base = wt; else base = Horizons(a);
			z = wt - top - 0.5 * zstep; // z is the position with respect to water table
			zfinal = wt - base;
			while (z >= zfinal) // integration of theta over z, start from base of horizon
			{
/*	Applying the Van Genuchten equation makes the model very slow because of the repeated integration; better use interpolation
theta is computed here by interpolating a pF curve at 1 cm suction intervals, computed in paramcheck.cpp */
				h = 100 * z;
				b = (int)(floor(h));
				theta = pFCurves(a, b + 1) - (h - b) * (pFCurves(a, b + 1) - pFCurves(a, b + 2));
				storage += theta * zstep; // add water in zstep to storage
				z -= zstep;
			}
			// last step, if a small difference < zstep remains
			z += zstep;
			if (z > zfinal)
			{
				dz = z - zfinal;
				z = z - 0.5 * dz;
				h = 100 * z;
				b = (int)(floor(h));
				theta = pFCurves(a, b + 1) - (h - b) * (pFCurves(a, b + 1) - pFCurves(a, b + 2));
				storage += theta * dz;
			}
		}
	}
	return storage;
}


void Watertable()
/* calculates water table position from precipitation and evaporation*/

{
	Matrix result(2);
	double base, V, dV, Vnew, maxV, Vbelow, runoff = 0.0, z, oldGW, cursign, startsign, p, prev, frozen = 0.0, dfrozen = 0.0, snowmelt = 0.0;// V: water volume in the soil
    double drainage = 0.0;
    double runon = 0.0;
	double dz = 1.0e-02;
	double f = 10.0;
	int tophorizon;
	
	// water balance change
	if (WatertableModel == 2)
	{
		base = FindWaterBase(); // determine base of soil watertable fluction range, being either the top of permafrost or the minimum watertable
		// for partly frozen soils we must keep track of changes in the frozen water content;
		// include conservation of frozen  water by calculating the storage in frozen water
        drainage = Drainage(CurrentGW); // drainage term
        // cout << drainage << endl;
		if (base != CurrentGWBase)			// frost table change; add to or release water from frozen storage
		{
				frozen = CalcFrozenStorage(CurrentGW); // calculate the amount of water that is stored frozen
				dfrozen = frozen - FrozenStorage;
				FrozenStorage = frozen;
				CurrentGWBase = base;
				if (CurrentGWBase >= 0.0) AboveWTStorage = BelowWTStorage = 0.0; else // recalculate above and below water table storage
				{
					if (CurrentGW > CurrentGWBase)  // current water table above frost table
					{
						if (BelowWTStorage > dfrozen) BelowWTStorage -= dfrozen; else // add or substract frozen storge chnge from storage below ater table
						{
							BelowWTStorage = 0.0;
							AboveWTStorage = AboveWTStorage - (dfrozen - BelowWTStorage);
						}
					} else AboveWTStorage -= dfrozen;  // if groundwater table below the frozen top, add or substract frozen storage change to above wt storage
				}
				//cout << CurrentGW << ' ' << CurrentFrozen << ' ' << dfrozen << ' ' << FrozenStorage << ' ' << (AboveWTStorage + BelowWTStorage + PondedWater - FrozenStorage) << endl;
				
		} 
		if (TData(StepNr) < 0.0)  // storage of water as snow, no water table change
		{
			SnowStorage +=  Precipitation(StepNr);
		} else
		{
			if (SnowStorage > 0.0) // empty snow storage to runoff and soil assuming that it all occurs in one day
			{
				snowmelt = SnowStorage;
				SnowStorage = 0.0;
			} else snowmelt = 0.0;
			V = AboveWTStorage + BelowWTStorage + PondedWater; // Water balance
            if (strlen(RunOnFile) != 0) runon = RunOn(StepNr); else runon = 0.0;
			if (CurrentGW > 0.0) // calculate water storage change
        // no evaporatien with below zero temperatures because soil that is frozen at the top can't evoprate soil water!
			{
				dV = Precipitation(StepNr) - OpenWaterFactor * Evaporation(StepNr) * CropFactor + snowmelt - drainage + runon; // open water evaporation
			} else
			{ 
				dV = Precipitation(StepNr) - Evaporation(StepNr) * CropFactor * exp(EvapCorrection * CurrentGW) + snowmelt - drainage + runon; // below surface water table: correct for lower evaporation
			}
		//cout << Precipitation(StepNr) << " " << Evaporation(StepNr) << " " << dV << endl;
			oldGW = CurrentGW;
			Vnew = V + dV;  // new water volume
			Vbelow = BelowWTStorage; // preserve below wt storage because belowWT function is being used to recalculate maximum storage
			tophorizon = belowWT(base, 0.0);
			maxV = BelowWTStorage;  // maximum water storage below surface, is recalculated everytime because the base can vary
			BelowWTStorage = Vbelow;
		// cout << maxV << endl;		
			if (Vnew >= maxV) // new water volume larger than max volume: calculate ponding and runoff
			{
				PondedWater = Vnew - maxV;
				AboveWTStorage = 0.0;
				if (PondedWater > RunoffThreshold)
				{
					runoff += PondedWater - RunoffThreshold;
					PondedWater =  RunoffThreshold;
					Vnew = Vnew - runoff;
				} else runoff = 0.0;
				CurrentGW = PondedWater;
				BelowWTStorage = maxV;
			} else // if the new volume < max volume, we should have a water table below surface
			{
				PondedWater = 0.0;
				runoff = 0.0;
				if (base < 0.0) // recalculation of storage and water table only when the water table can go below the soil surface
				{
					if (dV < 0) f = -1.0; else f = 1.0; // solution depth step factor depending on change direction of water table
					if (oldGW > 0.0) z = 0.0; else z = oldGW; // starting point of iteration
					tophorizon = belowWT(base, z);
					AboveWTStorage = aboveWT(tophorizon, z);
					p = BelowWTStorage + AboveWTStorage - Vnew;
				// iteration to find new water tabel
					startsign = -f;
					cursign = startsign;
					while (cursign == startsign)
					{
						z = z + f* dz;
						if (z > 0.0)
						{
							z = 0.0;
							cout << WT_ERROR2 << endl;
							break;
						}
						if (z < base)
						{
							if (z < MinGW) cout << WT_ERROR3 << MinGW << endl;
							break;
						}
						prev = p;
						tophorizon = belowWT(base, z);
						AboveWTStorage = aboveWT(tophorizon, z);
						p = BelowWTStorage + AboveWTStorage - Vnew;
						cursign = p / abs(p);
					}// Interpolate water table and correct for small errors
					if (z <= base) CurrentGW = base; else
					{
						CurrentGW = z - (f * dz) * abs(p) / abs(prev - p);
						if (CurrentGW <= base) CurrentGW = base;
					}
					tophorizon = belowWT(base, CurrentGW);  // recompute current storage
					AboveWTStorage = aboveWT(tophorizon, CurrentGW);
				}
			}
		}
	}
	result(2) = CurrentGW;
	result(1) = Timer + 0.5 * Timestep;
	if (ProfileOutput.Contains(7))
	{
		result.Write(output7);
	}
}

void Moisture(int initial)
/* calculates soil moisture profile assuming that soilmoisture is at gravitational equillibrium with the
groundwater table */
{
  int a, i, j;
  double h, w, lh, theta, t, s, currenttime;
  Matrix result(2);

  t = StepNr;
  if (initial) MatricPotential.Resize(NrLayers); else MatricPotential.Fill(0.0);  // array for matric potential
  if (WatertableModel < 2) {
      CurrentGW = GwData((int)t);  // groundwater level if not modelled; also write to output
	  currenttime = Timer + 0.5 * Timestep;
      result(1) = currenttime;
      result(2) = CurrentGW;
	  if (ProfileOutput.Contains(7))
	  {
		if (!initial) result.Write(output7);
	  }
  }
  if (MoistProfiles.Length() > 1)     // if present, use observational data or soil moisture model output
  {
    MoistTheta = MoistProfiles.Row((int)t);  // store soil moisture
    for (i = 1; i <= NrLayers; i++)     // interpolate corresponding matric potential from pf curves
    { // matric potential needs to be known for correction of decomposition rate for moisture in SOMdecomposition.cpp
      a = (int)Layers(i, 4);
	  j = 1;
      while ((MoistTheta(i) < pFCurves(a, j)) && (j < pFCurves.Cols())) j++;
	  if ((j > 1) && (j < pFCurves.Cols())) MatricPotential(i) = log10(j); else MatricPotential(i) = log10(pFCurves.Cols());
/* If the soil moisture is read from file, also the water table should be supplied in a file.
 * The soil moisture from file is here corrected with water table, to be sure that all layers below the water table are saturated
 */
      if (Layers(i, 1) <= CurrentGW) MoistTheta(i) = pFCurves(a, j);
    }
  }
  else
  {
    if (initial) MoistTheta.Resize(NrLayers);            // initialize theta array
	if (CurrentGW > 0.0) PondedWater = CurrentGW; else PondedWater = 0.0;
    for (i = 1; i <= NrLayers; i++)
    {
	
      a = (int)Layers(i, 4);                             // soil profile horizon number
      MoistTheta(i) = Porosity(a);                       // basic assumption: layer is saturated
      h = Layers(i, 1) - CurrentGW;                      // position of groundwater table with reference to layer top
/* NOTE: contratrary to the former version of PEATLAND, the layer is not saturated
if the groundwater table is anywhere below the top of the layer. If the water table is within the layer
the volumetric water content is computed based on the water content of the saturated part end the theta for
based on the potential halfway the unsaturated part. */
      if (h > 0.00001)                                   // layer completely or partly above the water table
      {
        if (h >= LayerThickness)                         // layer completely above the water table
        {
		  lh = 100 * (h - 0.5 * LayerThickness);  // water content based on tension halfway layer
          w = 0.0;
        } else
        {
		  lh = 100 * (0.5 * h);								// water content based on tension at top of unsaturated part
          w = ((LayerThickness - h) / LayerThickness) * Porosity(a); // water contribution from saturated part
        }
		if (lh > 1.0) MatricPotential(i) = log10(lh); else MatricPotential(i) = 0; // store matric potential for later use in environmental correction
        j = (int)(floor(lh));
		theta = pFCurves(a, j + 1) + (lh - j) * (pFCurves(a, j + 1) - pFCurves(a, j + 2));
		// theta is interpolated from suction curve lookup table with 1 cm resolution
		if (h < LayerThickness) theta = w + h * theta / LayerThickness; // in a partly saturated layer account for water in saturated part
		MoistTheta(i) = theta;                           // store theta in moisture profile array
      }
    }
  }
  if (initial)
  {
    Saturation.Resize(NrLayers);               // calculate saturation of pore volume with water
    OldSat.Resize(NrLayers);
    LastSatTime.Resize(NrLayers);
  }
  OldSat = Saturation;
  TopSat = 0;                                             // finds first saturated layer
  for (i = 1; i <= NrLayers; i++)
  {
    s = (PoreVol(i) - MoistTheta(i)) / PoreVol(i);
    if (s < 0.0001) s = 0.0;
    Saturation(i) = s;
    if (s > 0) TopSat = i;
    if (s == 0.0)
    {
      if (OldSat(i) > 0.0) LastSatTime(i) = 0.0; else LastSatTime(i) += (1.0 * Timestep);
    } else LastSatTime(i) = 0.0;
  }
  TopSat +=1;
  if ((!initial) && ProfileOutput.Contains(2)){
      Saturation.Write(output2b);
      MoistTheta.Write(output2a);
  }
}

