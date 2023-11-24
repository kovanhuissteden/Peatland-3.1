/***************************************************************************
     SOMdecomposition.h  -  soil organic matter decomposition PEATLAND
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

/************************ GLOBAL VARIABLES ***********************************/

extern Matrix CorrFac;            // environmental correction factors for aerobic SOM decomposition constants k
extern Matrix RootMass;           // initial root distribution (kg C/m2 in each layer)
extern Matrix Saturation;         // pore volume saturation with water
extern int NrLayers;              // Number of depth steps
extern Matrix SoilTemp;           // soil temperatures model layers, interpolated from TProfile
extern double T_ref;              // reference temperature for correction of decay constants
//extern double MolAct;             // molecular activation energy aerobic organic matter decomposition
extern double Rgas;               // Gas constant
extern Matrix Layer_pH;           // pH
extern Matrix Layers;             // Layer boundaries
extern Matrix MoistTheta;         // soil moisture profile for timestep
extern Matrix pFpoints;           // Curve for determining environmental correction factor for dryness
extern Matrix MatricPotential;    // matric potential of model soil layer
extern double HalfSatPoint;       // half activity saturation point for correction factor aeration
extern double RootAeration;       // Root mass dependent correction (0 - 1) for improved aeration by root growth if 0, this is switched off
extern int NrReservoirs;          // Number of SOM reservoirs
extern double MaxNPP;             // Maximum primary productivity (kgC/m2/day)
extern double MinNPP;             // Minimum primary productivity
extern double PrimingCorrection;  // Root mass dependent priming effect root exudates on slow C reservoirs; value > 0; if 0, this is switched off
extern double NPP;                // Primary production per time step
extern double Timestep;           // model timestep
extern double SpringFactor;       // instantaneous value of spring correction, declared global for use by environmental correction SOM decomposition for priming effect
extern Matrix Kdecay;             // SOM decomposition constants for each reservoir
extern Matrix CNRatio;                    // CN ratios for eac soil layer; the decomposition of peat can be made dependent on these
extern Matrix KPeatCN;                    // constants linear relation of decomposition rate k of peat with CN ratio cf Vermeulen & Hendriks
extern Matrix AerobicQ10;                // Q10 for each reservoir
extern int Q10orArrhenius;             //Switch between temperature correction of (an)aerobic decomposition as Q10 (0) or Arrhenius (1) equation
extern double AnaerobicAssimDissim;      // Assimiltion/Dissimilation ratio anaeroob
extern double CurrentGW;          // Current water table during time step
extern double LayerThickness;     // Thickness of depth step
extern double LitterLayer;              // organic matter stored in above ground litter layer, in kg C / m2
extern double LitterDecomp;              // Litter decomposition kg C per timestep
extern Matrix TData;                           // Air or soil temperature data from file Tfile
extern Matrix BioMassRec;                    // storage of biomass, primary production and plant respiration
extern Matrix NewSOM;             // SOM reservoirs to be changed in each iteration step
extern Matrix OldSOM;                           // SOM reservoirs to be changed in each iteration step (kg C per layer) for calculation of storage change
extern Matrix PeatDecay;                        // logs true loss of peat matrix
extern double PeatLoss;           // Totalized aerobic loss of peat C over all layers
extern Matrix SplitRes;           // partitions decomposed material between CO2 + microbial biomass (1st column) and resistant SOM
extern Matrix CO2;                // CO2 evolved from each layer and reservoir
extern Matrix ProfileOutput;      // determines which vertical profiles are sent to log files
extern ofstream *output5;         // file handles output
extern ofstream *output9;
extern ofstream *output11;
extern ofstream *output12;
extern ofstream *output13;
extern ofstream *output14;
extern ofstream *output15;
extern ofstream *output16;
extern ofstream *output17;
extern int StepNr;                // time step number during iteration
extern Matrix ReservoirTime;      // storage matrix for CO2 per reservoir per timestep
extern Matrix LayerTime;          // storage matrix for CO2 per layer per timestep

extern Matrix CO2FromMethaneOx;     // CO2 from methane oxidation
extern double DayNr;              // midpoint of simulated timestep, relative to day 1 of the year in which the simulation started
extern int AnaerobicCO2;          // Switch for allowing anaerobic decomposition (sulfate etc) resulting in CO2
extern Matrix KAnaerobic;         // Anaerobic decomposition constants
extern double Q10Anaerobic;       // Q10 anaerobic decomposition
extern Matrix LayerAnaerobic;      // Anaerobic CO2 per layer
extern Matrix AnaerobSum;          // sum of anaerobic CO2 per layer
extern double MethaneTRef;        // Reference temperature for temperature sensitivity methane production
extern double KLitter;            // decomposition constant above-ground litter and standing dead biomass
extern double OldLitter;                 // For calculation of storage change of litter layer
extern Matrix CarbonBalance;             // Carbon balance: primary production, C exported, and change in carbon reservoirs in Mol C
/******************** FUNCTION DEFINITIONS ***********************************/


void EnviCor();
/* environmental correction factors for first order decomposition constants
aerobic decomposition of the SOM reservoirs, includes:
Temperature correction
pH correction
soil moisture and soil dryness corrections
Priming correction                                                         */

void Decompose();
/* aerobic decomposition of SOM above the water table.*/

void WriteSOMReservoirs();
/* writes contents of SOM reservoirs to output file */

void CollectCO2();
/* Collect all CO2 and store in output arrays */
