/***************************************************************************
                  methane.h  -  methane model PEATLAND cf Walter(2000)
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

//#include "matrix.h"

#define METHANE_ERROR1 "Methane model: unlikely high plant transport"
#define METHANE_ERROR2 "Negative CH4 concentration! Numerical instability!"
#define METHANE_ERROR3 "Soil methane concentration > 5.0e5"
#define MAXSOILCH4 5.0e5
//#define CONVCH4CTOKGC 0.012011 // based on moles CH4
#define CONVCH4CTOKGC 1.2011e-5   // based on millimoles CH4


/************************ GLOBAL VARIABLES ***********************************/

extern Matrix MethaneFlux;              // Methane flux at surface;  1st value: ebullition flux; 2nd: plant mediated flux; 3d:diffusive flux
extern Matrix Layers;                   // Layer boundaries
extern int NrLayers;                    // Number of depth steps
extern double LayerThickness;           // Thickness of depth step
extern Matrix MethProfile;              // methane concentration profile
extern Matrix Saturation;               // pore volume saturation with water
extern Matrix PoreVol;                  // porosity each layer
extern Matrix SoilTemp;                 // soil temperatures model layers, interpolated from TProfile
extern int NrReservoirs;                // Number of SOM reservoirs
extern int TopSat;                      // Index of first saturated layer
extern Matrix NewSOM;                   // SOM reservoirs to be changed in each iteration step
extern Matrix MethaneReservoirs;        // This parameter tells which reservoirs are summed for calculating methane production from easily decomposed
extern double MethaneReservoirSum;      // sum of methane reservoirs
extern double MethaneDiff;              // diffusion of methane in air in m2/d
extern double MethaneDiffWater;         // diffusion of methane in water
extern double Timestep;                 // model timestep
extern double MethaneMaxConc;           // maximum methane concentration in pore water
extern double MethaneERateC;            // Ebullition rate constant (1/hr)
extern double CurrentGW;                // Current water table during time step
extern Matrix MethaneFlux;              // Methane flux at surface;  1st value: ebullition flux; 2nd: plant mediated flux; 3d:diffusive flux
extern double MaxProd;                  // Maximum primary productivity (kgC/m2/day)
extern double GrowFuncConst;            // proportionality constant growth functioen - primary productivity for plant transport in methane model
extern double PrimProd;                 // Primary production per time step
extern double MethanePRateC;            // Rate constant for plant transport of methane (1/hr)
extern double MethanePType;             // Vegetation type factor for gas transport by plants range: 0-15
extern double MethanePlantOx;           // Fraction of methane that is oxidized during transport in plants
extern double PartialAnaerobe;           // Determines the slope of the relation of partial anaerobe soil fraction above the water table to soil saturation, >1
extern Matrix AnaerobSum;                // sum of anaerobic CO2 per layer
extern double AnaerobeLagFactor;	// Determines time lag for development of sufficiently anaerobic conditions after saturation of a layer
extern Matrix LastSatTime;              // Time after last complete saturation of layer in days
extern Matrix RootDistrib;              // root distribution function
extern double MethaneAir;               // Methane concentration in the atmosphere specified in ppmv, for calculations inside peatland in moles in a layer of 10 cm above surface
extern double MethaneQ10;               // Q10 value for temperature correction methane production; range 1.7 - 16 ref. in Walther & Heimann 2000
extern double MethaneTRef;              // Reference temperature for temperature sensitivity methane production
extern Matrix MethaneR0Corr;            // pH dependent methane production rate
extern double MethaneOxQ10;             // Q10 value for temperature correction methane oxidation; range 1.4 - 2.1, ref. in Walther & Heimann 2000
extern double MethaneVmax;              // Vmax Michaelis-Menten eq methane oxidation micrMol/hr range 5-50
extern double MethaneKm;                // Km Michaelis-Menten eq methane oxidation micrMol range 1-5
extern double CO2CH4ratio;              // Molar ratio between CH4 and CO2 production; for acetate splitting this is 1, for CO2 reduction 0
extern Matrix TotalMethane;             // storage matrix for CH4 results
extern int StepNr;                      // time step number during iteration
extern Matrix ProfileOutput;            // determines which vertical profiles are sent to log files
extern ofstream *output3;               // methane profile output file
extern Matrix CO2FromMethane;           // CO2 from methane oxidation
extern double DayNr;                    // midpoint of simulated timestep, relative to day 1 of the year in which the simulation started
extern Matrix UnFrozen;                 // Unfrozen water content (kg water / kg dry soil)
extern double DensWater;                // density of water at 0 degr C
extern Matrix DBD;                      // Dry bulk density for each horizon
extern Matrix Ice;                      // Ice content (kg ice / kg dry soil)
extern int ProductionModel;                  // Production model: 0 for simple sinusoidal function; 1 for production dependent on temperature of upper soil layer
extern double CurrentLAI;               // leaf area index
extern double MinProd;                   // Minimum primary productivity
extern Matrix TData;                      // Air temperature data from file Tfile
/******************** FUNCTION DEFINITIONS ***********************************/

void Methane();
/* methane fluxes modified after the model of Walther (2000), Global Biogeochemical Cycles 14, 745 - 765.*/

void ebull(double &flux, Matrix &rate, Matrix &methconc);
/*  calculates ebullition flux and rate of removal of CH4 (sink) by ebullition
    no entrapment cf Walter et al of bubbles is assumed, entrapment is assumed to occur within the same layer
    Input:
    mc		: methane concentration
    Output:

    flux		: total CH4 bubble flux
    rate		: rate of CH4 removal
*/

void planttrans(double &flux, Matrix &CH4oxidation, Matrix &CO2production, Matrix &rate, Matrix &methconc);
/*  calculates plant transport flux and rate of removal of CH4 (sink) by plant transport
    Input:
    methconc:   methane concentration

    Output:

    flux:       total CH4 bubble flux
    oxidation:  total CH4 oxidized in plant root system
    CO2production: toatal CO2 produced by methane oxidation;
    rate:       rate of CH4 removal
*/

void methaneprod(Matrix &prod, Matrix &labileC, Matrix &mc, Matrix &oxidized, Matrix &anaerob, double ttime);
/*  Below water table methane production from soil organic matter and
    oxidation in unsaturated zone as function of soil CH4 concentration

Input:

c       : total labile C soil in micromol;
pc      : peat C
mc      : methane concentration

Output:

prod    : total production or oxidation
ttime   : time in days since start of model time step for calculation of rapid saturation depression of methane production
oxidized: matrix with moles of methane oxidized per layer
anaerob: anaerobically prduced CO2

*/
