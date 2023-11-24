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

#define HEAT_ERROR1 "Thermal model: unstable solution, increase DStepHeat or decrease TStepHeat"

/************************ GLOBAL VARIABLES *********************************/

extern double TStepHeat;                 // time step (days) temperature model
extern double DStepHeat;                 // minimum depth step (m) temperature model
extern Matrix SoilTData;                  // soil temperature data
extern double MaxDepthHeat;              // maximum depth temperature model (m)
extern int NrHeatLayers;                 // number of layers temperature model
extern Matrix Saturation;                // pore volume saturation with water
extern Matrix ThermDiffVar;              // thermal diffusivity layer dependent diffusivity
extern Matrix HeatLayers;                // layer midpoints layers thermal model
extern double ThermDiff;                 // thermal diffusivity, if negative it is estimated from soil properties
extern int StepNr;                       // time step number during iteration
extern ofstream *output1;                // output log files
extern double CondOrg;                   // thermal conductivity organic matter J.m-1.s-1.K-1
extern double CondMiner;                 // thermal conductivity mineral matter
extern double CondAir;                   // thermal conductivity air
extern double CondWater;                 // thermal conductivity water
extern double CondQuartz;                 // thermal conductivity quartz
extern double CondIce;                   // thermal conductivity ice
extern double HCIce;                     // volumetric heat capacity ice
extern Matrix LatentHeat;		 // parameters for approximation of temperature-dependent latent heat of fusion of ice J kg-1
extern double MaxSnowdepth;              // Maximum snow depth
extern double DayMaxSnowdepth;           // Day of maximum snow depth
extern double SnowMeltrate;              // Rate of snowmelt per degree C above zero per day
extern double SnowDepth;                 // Actual snow heigth (m)
extern double DayOfTheYear;              // Julian day number of the midpoint of the simulated time step relative to the current year;
extern double SnowStartDay;              // day at wich snow accumulation has started
extern char SnowFile[];                  // file with snowdepth time series
extern Matrix SnowData;                  // Snow depth data from file Snowfile
extern double SnowMeltrate;              // Rate of snowmelt per degree C above zero per day
extern double HeatCondTop;               // heat conductivity of the top soil
extern double CondSnow;                  // thermal conductivity snow
extern Matrix SandFraction;			    // sand weight fraction of mineral fraction
extern double VegTScalingFactor;			// scaling factor for air to soil surface temperature; put to one if soil surface temperature is input; outherwise a valye between 0.6 and 1.0
extern double CurrentLAI;           // leaf area index
extern int ProductionModel;         // productionmodel


/************************ FUNCTION HEADERS *********************************/


void HeatVar(const int steps, const double tsurf);
/* Solves the heat equation for heat transport into the soil using the classical explicit scheme;
  steps is the number of time steps, tsurf the surface temperature                            */

void HeatSimple(double time);
/* Analytical solution of the heat equation based on constant thermal diffusivity and sinusoidal temperature time series
   cf Animo - see Groenendijk and Kroes, 1997                                                  */

void Temperature();
/* computes temperature soil profile dependent on temperature profile                          */

void ThermalProperties(Matrix T);
/* computes heat capacity, thermal conductivity and diffusivity based on soil constituents
   approach cf Hillel (1980)
   T: temperature of each soil layer                                                            */

double SnowVegCorrection(const double tsurf);
/* Corrects soil surface temperature in case of snow cover presence
and for vegetation when there is no snow */
