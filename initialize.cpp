/***************************************************************************
               initialize.cpp  -  initialization functions of PEATLAND
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

  June 2003

  Correction of bug in calculation of sinusoidal groundwater table time series

  31 May 2006
  Correction of bug in opening/closing of vertical profile output files
  17 Aug 2021 bug fix in function RootsInit
  
  7 april 2022
  RootLambda was not used in calculating root distribution; fixed

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
#include "water.h"

void MakeLayers()
/* Initializes model layers and interpolates data from the soil profile description
   Initalizes the pore volume array
   Initializes the soil top reference level                                                        */
{
  int i, j, k = 0;

  Layers.Resize(NrLayers, 4);                         // initialize layer matrix
  PoreVol.Resize(NrLayers);                           // initialize pore volume array, to be used in several other functions
  for (i = 1; i <= NrLayers; i++)
  {
    Layers(i, 1) = - (i - 1) * LayerThickness;        // layer top
    Layers(i, 2) = - (i - 0.5) * LayerThickness;      // layer midpoint
    Layers(i, 3) = - i * LayerThickness;              // layer base
    for (j = NrHorizons; j >= 1; j--)
    {
      if (-Layers(i, 2) < Horizons(j))
      {
        Layers(i, 4) = (double)j;
        k = j;
      }
    }
    // 4th column: reference to soil horizon based on layer midpoint
    PoreVol(i) = Porosity(k);                          // pore volume
  }
  RefLevel = 0;                                        // reference level; set at surface at start of model run
}

void RootsInit()                                          // calculates root distribution
{

  int i = 1;
  double s;

  if (NoRootsBelowGWT) MaxRootDepth = -MinGW;         // set max rootdepth to lowest level if roots do not grow beneath the groundwater table
  RootDistrib.Resize(NrLayers);                       // initialize root distribution function
  while (Layers(i, 3) >= (-MaxRootDepth))             // exponential function
  {
    RootDistrib(i) = exp(RootLambda * Layers(i, 2) / MaxRootDepth);
    i++;
  }
  s = RootDistrib.Sum();                              // normalize root distribution function
  RootDistrib /= s;
  RootMass = RootDistrib;                             // Calculate Rootmass
  RootMass *= InitRoots;
}

void OutputInit()                                     // initializes output arrays
{
  TotalReservoir.Resize(NrLayers, NrReservoirs);      // total C per reservoir per layer
  ReservoirTime.Resize(NrOfSteps, NrReservoirs + 3);
/* storage matrix for CO2 per reservoir per timestep; an extra entry is added for CO2 from Methane oxidation */
  LayerTime.Resize(NrOfSteps, NrLayers + 1);           // storage matrix for CO2 per layer per timestep
                                                      // NOTE: only Walters model is used here, in the Matlab version there was a choice of two models
  BioMassRec.Resize(NrOfSteps, 10);                    // storage of biomass, 1: time (days); 2: total plant biomass including roots (kg C.m-2); 3: primary production; 4: // plant respiration; 5: net CO2 flux (soil respiration + plant repiration - primary production); 6: soil respiration + dark respiration of vegetation (mg.m-2.hr-1); 7: // litter mass (kg C.m-2); 8: biomass removed by harvest and grazing 9: LAI

  LayerAnaerobic.Resize(NrOfSteps, NrLayers);
  CarbonBalance.Resize(NrOfSteps, 26);
}

void SOMResInit()                                     // calculates organic C (kg C per layer) per SOM reservoir per layer
{
  int i, j, a;
  double ol;

  InitSOM.Resize(NrLayers, NrReservoirs);
  CO2.Resize(NrLayers, NrReservoirs);
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i, 4);                              // soil profile horizon number
    ol = DBD(a) * LayerThickness * PercOrg(a) / 100;    // kg organic matter per layer
    for (j = 1; j <= NrReservoirs; j++)
    {
      InitSOM(i, j) = ol * InitRes(a, j) * Cfrac(j);    // multiply (kg organic matter per layer) with (fraction organic matter) in reservoir and (reservoir C fraction)
    }
  }
  NewSOM = InitSOM;
  ResYearSOM = InitSOM;
}

void InitHeat()
/* computes heat capacity and thermal conductivity from organic matter content and pore volume
 cf Hillel (1998) and Luckner & Schestakow (1991)                                           */

//////!!!!!!!!! T_average: voor heat model = 0, calculate Taverage from data!
{
  int i, a, l, j;
  double theta_sat, theta_wilt, f0, t, d, w;
  Matrix tp1, tp2, lm1;

  SoilTemp.Resize(NrLayers);                            // Soil temperature array for easy use in other functions, wil be updated by heat functions
  if (ThermModel == 2) return;                          // skip the rest if soil temperatures are taken from file
  if (ThermModel == 0)									// recompute average temperature and amplitude from data
  {														// the average temperature is the lower boundary condition of the temperature solutution
   T_amplitude = TData.Max() - TData.Min();
  }
  ThermDiffVar.Resize(NrLayers);                        // thermal diffusivity
  Freeze.Resize(NrLayers, 3);                           // Freezing fumction parameters
  UnFrozen.Resize(NrLayers);                            // unfrozen water content (kg water / kg dry soil)
  Ice.Resize(NrLayers);                                 // ice content (kg ice/kg dry soil)
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i, 4);                              // soil profile horizon number
    theta_sat = pFCurves(a, 1);                         // parameters for soil freezing based on volumetric water content at saturation and wilting point, or otherwise driest part of pF curve
    if (ApplyVanGenuchten) theta_wilt = Layer_pF(a, 1); else theta_wilt = Layer_pF(a, pFVal.Length());
    f0 = DensWater * theta_sat / DBD(a);                // weight fraction water content at 0 degr
    Freeze(i, 1) = DensWater * theta_wilt / DBD(a);     // weight fraction water content at maximum freezing
    Freeze(i, 2) = FreezingCurve(a);                    // a parameter of freezing curve
    Freeze(i, 3) = exp(- (1 / FreezingCurve(a)) * log (f0 - Freeze(i, 1))); // b parameter
/* determine unfrozen water content*/
    l = T_init.Cols();                                  // determine actual layer temperature from T_Init
    d = -Layers(i, 2);
    j = 1;
    while ((d > T_init(1, j + 1)) && (j < l)) j++;      // look up position of layer
    t = T_init(2, j) + (T_init(2, j + 1) - T_init(2, j)) * (d - T_init(1, j)) / (T_init(1, j + 1) - T_init(1, j));  // interpolate temperature
    if (t < 0) UnFrozen(i) = Freeze(i, 1) + 1 / pow((Freeze(i, 3) - t), Freeze(i, 2));  // Unfrozen water factor
    w = MoistTheta(i) - UnFrozen(i) * DBD(a) / DensWater;  // ice content: if there is more water than the unfrozen water content acxording to the unfrozen water function, it is ice
    if ((t < 0) && (w > 0)) Ice(i) = w ;
  }
/* generalized thermal diffusivity for wet soil for analytical solution of heat equation (ThermModel == 1) */
  NrHeatLayers = (int)ceil(MaxDepthHeat/DStepHeat);
/* number of layers temperature model, generally larger than NrLayers to avoid numerical problems in solving the heat equation */
  TProfile.Resize(NrHeatLayers);                        // initialize temperature profile arrays
  tp1 = T_init.Row(1);                                  // by interpolation from initial T profile
  tp2 = T_init.Row(2);
  HeatLayers.Resize(NrHeatLayers);
  for (i = 1; i <= NrHeatLayers; i++) HeatLayers(i) = 0.5 * DStepHeat + (i - 1) * DStepHeat;
  interp(tp1, tp2, HeatLayers, TProfile, FALSE);
  -HeatLayers;
  lm1 = Layers.Col(2);
  interp(HeatLayers, TProfile, lm1, SoilTemp, TRUE);
}


double thermcond(double f, double l1, double l2)
/*  computes thermal conductivity for for a fractional mixture of two soil constituents
    assuming absence of layering, according to eq 1.108c in Luckner & Schestakow 1991

    Input:
    f	: fraction of the material with lowest conductivity
    l1	: conductivity corresponding with this material
    l2	: conductivity of the other stuff
    Output:
    lambda: thermal conductivity of the mixture                                                 */
{
  double lambda, f1, f2;

  f1 = pow((1 - f), 2/3);
  f2 = l1 / l2;
  lambda = l1 * (f1 + f2 * pow(f, 2/3))/(f1 - (1 - f) + f2 * (2 - f - f1));

  return lambda;
}

void InitTime()           // initializes time system
{
 
  DayNr = StartDay - 0.5 * Timestep;         // day number since day 1 of the year in which the simulation started
  Year = StartYear;
  DayOfTheYear = DayNr;
}


void InitTseries()
// initializes temperature  time series if they have not been read from file
{

  int i;
  double day;

  if (TData.Length() == 1)                  // sinusoidal temperature time series
  {
    day = StartDay + 0.5 * Timestep;
    TData.Resize(NrOfSteps);
    for (i = 1; i <= NrOfSteps; i++)
    {
      TData(i) = T_average - T_amplitude * cos(2 * PI *((day - 15) / YEAR)); // simple sinusoidal function 
      day += Timestep;
    }
  }
  PeatDecay.Resize(NrOfSteps, 2);         // PeatDecay logs true loss of peat matrix
}

void InitWater()
// initializes water table model  
// NB: The lookup table for the pF curves is calculated in function Paramchck!!!!!!!!!
{
  double base, day, meangw;
  int i, tophorizon;
  
  if (WatertableModel == 1)             // sinusoidal groundwater table time series
  {
    day = StartDay + 0.5 * Timestep;
    meangw = MinGW + 0.5 * AmplitudeGW;
    GwData.Resize(NrOfSteps);
    for (i = 1; i <= NrOfSteps; i++)
    {
      GwData(i) = meangw - 0.5 * AmplitudeGW * cos(2 * PI * (day - DayMinGW) / YEAR); // simple sinusoidal function, offset by the day with minimum groundwater table
      day += Timestep;
    }
  }
  if (WatertableModel == 2)             // calculated water table
  {
	GwData.Resize(NrOfSteps);			//prepare time series to store watertable data
	CurrentGW = WatertableInit;			// initial water table
	if (CurrentGW > 0.0) PondedWater = CurrentGW; else PondedWater = 0.0;
	base = FindWaterBase();				// initialize soil water storage
	CurrentGWBase = base;
	if (base > MinGW) FrozenStorage = CalcFrozenStorage(CurrentGW); else FrozenStorage = 0.0; 
	tophorizon = belowWT(base, CurrentGW);
	if (tophorizon > 0) AboveWTStorage = aboveWT(tophorizon, CurrentGW); else AboveWTStorage = 0.0;
	Precipitation = Precipitation / 1000.0; // scale precipitation and evaporation from mm to m
	Evaporation = Evaporation / 1000.0;
    if (strlen(RunOnFile) != 0) RunOn = RunOn / 1000.0;
  }
}

void InitLogFiles()
/*initializes files for logging model state variables */
{
  char buf[1024];
  int i;

  for (i = 1; i <= ProfileOutput.Length(); i++)
  {
    strcpy(buf, &DataDir[0]);
    strcat(buf, &OutputFilePrefix[0]);
    if (ProfileOutput(i) == 1)
    {
      output1 = new ofstream(strcat(buf, OUTPUT1));
      if (!output1)
      {
        cout  << INIT_ERROR2 << OUTPUT1 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 2)
    {
      output2 = new ofstream(strcat(buf, OUTPUT2));
      if (!output2)
      {
        cout  << INIT_ERROR2 << OUTPUT2 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 3)
    {
      output3 = new ofstream(strcat(buf, OUTPUT3));
      if (!output3)
      {
        cout  << INIT_ERROR2 << OUTPUT3 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 4)
    {
      output4 = new ofstream(strcat(buf, OUTPUT4));
      if (!output4)
      {
        cout  << INIT_ERROR2 << OUTPUT4 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 5)
    {
      output5 = new ofstream(strcat(buf, OUTPUT5));
      if (!output5)
      {
        cout  << INIT_ERROR2 << OUTPUT5 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 6)
    {
      output6 = new ofstream(strcat(buf, OUTPUT6));
      if (!output6)
      {
        cout  << INIT_ERROR2 << OUTPUT6 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 7)
    {
        output7 = new ofstream(strcat(buf, OUTPUT7));
        if (!output7)
        {
            cout  << INIT_ERROR2 << OUTPUT7 << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (ProfileOutput(i) == 8)
    {
        output8 = new ofstream(strcat(buf, OUTPUT8));
        if (!output8)
        {
            cout  << INIT_ERROR2 << OUTPUT8 << endl;
            exit(EXIT_FAILURE);
        }
    }
      if (ProfileOutput(i) == 9)
      {
          output9 = new ofstream(strcat(buf, OUTPUT9));
          if (!output9)
          {
              cout  << INIT_ERROR2 << OUTPUT9 << endl;
              exit(EXIT_FAILURE);
          }
      }
    if (ProfileOutput(i) == 11)
    {
      output11 = new ofstream(strcat(buf, OUTPUT11));
      if (!output11)
      {
        cout  << INIT_ERROR2 << OUTPUT11 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 12)
    {
      output12 = new ofstream(strcat(buf, OUTPUT12));
      if (!output12)
      {
        cout  << INIT_ERROR2 << OUTPUT12 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 13)
    {
      output13 = new ofstream(strcat(buf, OUTPUT13));
      if (!output13)
      {
        cout  << INIT_ERROR2 << OUTPUT13 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 14)
    {
      output14 = new ofstream(strcat(buf, OUTPUT14));
      if (!output14)
      {
        cout  << INIT_ERROR2 << OUTPUT14 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 15)
    {
      output15 = new ofstream(strcat(buf, OUTPUT15));
      if (!output15)
      {
        cout  << INIT_ERROR2 << OUTPUT15 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 16)
    {
      output16 = new ofstream(strcat(buf, OUTPUT16));
      if (!output16)
      {
        cout  << INIT_ERROR2 << OUTPUT16 << endl;
        exit(EXIT_FAILURE);
      }
    }
    if (ProfileOutput(i) == 17)
    {
      output17 = new ofstream(strcat(buf, OUTPUT17));
      if (!output17)
      {
        cout  << INIT_ERROR2 << OUTPUT17 << endl;
        exit(EXIT_FAILURE);
      }
	}
  }
}

void CloseLogFiles()        // closes all log files
{
int i;

  for (i = 1; i <= ProfileOutput.Length(); i++)
  {
    if (ProfileOutput(i) == 1) output1->close();
    if (ProfileOutput(i) == 2) output2->close();
    if (ProfileOutput(i) == 3) output3->close();
    if (ProfileOutput(i) == 4) output4->close();
    if (ProfileOutput(i) == 5) output5->close();
    if (ProfileOutput(i) == 6) output6->close();
    if (ProfileOutput(i) == 7) output7->close();
    if (ProfileOutput(i) == 8) output8->close();
    if (ProfileOutput(i) == 9) output9->close();
    if (ProfileOutput(i) == 11) output11->close();
    if (ProfileOutput(i) == 12) output12->close();
    if (ProfileOutput(i) == 13) output13->close();
    if (ProfileOutput(i) == 14) output14->close();
    if (ProfileOutput(i) == 15) output15->close();
    if (ProfileOutput(i) == 16) output16->close();
    if (ProfileOutput(i) == 17) output17->close();
    if (ProfileOutput(i) == 17) output17->close();
  }
}



void  InitDecomp()
/* Corrects aerobic decomposition constants based on assimilation ratio
Decompositon constant for peat is corrected by C/N ratio of layers */

{
  int i, a;
  double k, ad, c;

  /* Deleted code with old corrections
  ad = (1.0 - ResistFrac)/(1.0 + AssimDissim);              // correction for assimilation:
  // (1 - ResistFrac) is part that is not transferred to humus reservoir; 1/(1+AssimDissim) is fraction of this that goes to microbial population
  c = ad + ResistFrac;
  for (i = 1; i <= (NrReservoirs - 2); i++)
  {
    SplitRes(i, 1) = ad;                                // partitioning coeficients to resistant SOM and microbial biomass
    SplitRes(i, 2) = ResistFrac;
    k = Kdecay(i);                                      // correction of k
    Kdecay(i) = k + k * c;
  } */
  
  /* Simplify:
   * Fixed fraction to humus reservoir, based on ResistFrac
   * Fixed fraction to microbial biomass, based on AssimDissim
   * Rename AssimDissim to DissimAssimRatio for correct terminology ??
   * No correction of k
   * A further simplification could be to dowmsize SplitRes to a two element vector as in the current configuration ad and resistFrac is the same for all reservoirs
   * However, maintaining SplitRes as a 5 x 2 matrix allows to have different dissimilation - assimilation ratios per reservoir, e.g. a lower rate for labile stuff
   */
  // New code:
  ad = 1.0/(1.0 + DissimAssimRatio); // assimilation factor, part of decomposed organic matter carbon that is transferred to microbial bioomass
  for (i = 1; i <= (NrReservoirs - 2); i++)
  {
    SplitRes(i, 1) = ad;                                // partitioning coeficients to resistant SOM and microbial biomass
    SplitRes(i, 2) = ResistFrac;
    SplitRes(i, 3) = 1.0 / (1.0 + AnaerobicDARatio);    // partitioning coefficient to microbial biomass from anaerobic CO2
  }
  KPeat.Resize(NrLayers);                               // C/N ratio dependent k for peat
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i,4);                               // soil horizon index
    k = KPeatCN(1) - KPeatCN(2) * CNRatio(a);           // correct k for CN based on empirical linear relation
    KPeat(i) = k + k * c;                               // apply correction as above
  }
// SplitRes: fractions of the decayed component that are tranferred from primary to secondary reservoir
// transfer to biomass: to be based on dissimilation/assimilation ratio, 2 for fungi and 2.3 for bacteria,
// and death rate of biomass
// if transfer to humus (h) amounts 0.1, 0.9 is transferred to CO2 or biomass;
// using dissimilation/assimilation ratio a = 2, the transfer to biomass should be 0.3.
// Formula: transfer fraction to biomass = (1-h)/(a+1)
// The decay constant k is usually measured from transfer to CO2 only. The total transfer (corrected k)
// is the sum of transfer to CO2, microbial biomass and humus.
  CorrFac.Resize(NrLayers, NrReservoirs);
  AnaerobSum.Resize(NrLayers);          // sum of anaerobic CO2 per layer
}  // end of function CorrectDecomp


void InitMethaneModel()
/* Initializations for methane model */
{
  int i, a;
  double mc;

  MethProfile = (InitMethane * PoreVol) * LayerThickness;  // initialize methane profile and convert from millimol/m3 in porewater to millimol/layer
  MethaneMaxConc = MethaneMaxConc * LayerThickness;  // MethaneMaxConc has to be adapted to porevolume on calculation of ebullition in methane.cpp
  OxconCH4.Resize(NrLayers);                          // initialize oxygen consumption CH4 oxidation array
  TotalMethane.Resize(NrOfSteps, 5);                  // initalize total methane
  MethaneFlux.Resize(3);
  MethaneR0Corr.Resize(NrLayers);                     // pH dependent methane production rate
  MethaneR0 = MethaneR0 * LayerThickness;                     // convert MethaneR0 from micromol/L/hr to millimol/layer/hr  1.0e3 (L-m3) * 1.03-3 (micro -milli) *LayerThickness
  MethaneAir = MethaneAir * 1.218581;
  // concentration in the air is recalculated here from ppmv to micromoles CH4-C in a 10 cm layer above the surface 
  // at standard seal level pressure and the current air temperature
  // MethaneAir = MethaneAir * 0.000001 * 1000 * 0.1 * 101325 / 8.3145;
  // MethaneAir = MethaneAir * (ppmv to partial pressure) * (moles to millimoles) * (m3 to 10 cm layer) * standard pressure / gas constant
  // to be divided by actual air absolute temperature at each time step
  // eventually replace with real pressure and air temperature if air pressure becomes input to the model
  
  // Scale MethaneKm from micromol/L (millimol/m3) to micromol/layer
  MethaneKm = MethaneKm / LayerThickness;
  MethaneVmax = MethaneVmax / LayerThickness;
  MethaneReservoirSum = MethaneReservoirs.SumRow(1);
  // MethaneReservoirs.Transp();                         // transpose methane reservoirs array for easier multiplication with carbon reservoirs
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i,4);                             // soil horizon index
    mc = MethaneR0 + (Layer_pH(a) - 7.0) * MethanepHCorr * MethaneR0;  // correction of methane R0 for pH
    if (mc < 0)
    {
      cout << INIT_ERROR3 << endl;
      exit(EXIT_FAILURE);
    }
    MethaneR0Corr(i) = mc;
  }
}
