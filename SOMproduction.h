/***************************************************************************
        SOMproduction.cpp  -  soil organic matter production PEATLAND
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

#define PRODUCTION_ERROR1 "Primary production model: negative root mass, decrease RootSenescence parameter"
#define PRODUCTION_ERROR2 "Primary production model: above-ground Biomass below minimum, check BioMassSenescence, LAICarbonFraction parameter"
#define PRODUCTION_ERROR3 "Primary production model: PARunits parameter specified incorrectly, should be 0, 1 or 2"
extern BOOLEAN Verbose;
extern double DayOfTheYear;                  // Julian day number of the midpoint of the simulated time step relative to the current year;
extern int ProductionModel;                  // Production model: 0 for simple sinusoidal function; 1 for production dependent on temperature of upper soil layer
extern double MaxProd;                       // Maximum primary productivity (kgC/m2/day)
extern double MinProd;                       // Minimum primary productivity
extern double Timestep;                      // model timestep
extern double Timer;                      // starting point of simulated time step
extern double DayNr;                         // midpoint of simulated timestep, relative to day 1 of the year in which the simulation started
extern int StepNr;                           // time step number during iteration
extern Matrix ProdTFunc;                     // determines temperature dependent production rate; 1st number is minimum, 2nd optimum
extern Matrix SoilTemp;                      // soil temperatures model layers, interpolated from TProfile
extern double PrimProd;                      // Primary production per time step
extern double SpringCorrection;              // Correction (0-1) for stronger exudation in spring; influences priming and exudate production, if 0, disables spring correction; correction is a factor of 1+SpringCorrection
extern double SatCorr;                       // correction of production for saturation of topsoil, depresses production at high saturation, switched off when 0
extern Matrix Saturation;                    // pore volume saturation with water
extern double TotalPrimProd;                 // total primary production
extern double Shoots;                        // Shoot production per time step
extern double ShootsFactor;                  // mass fraction root growth against shoot growth
extern Matrix RootDistrib;                   // root distribution function
extern int NoRootsBelowGWT;                  // if 1, no roots will grow below groundwater table (No telmatophytes)
extern int NrLayers;                         // Number of depth steps
extern Matrix Layers;                        // Layer boundaries
extern Matrix GwData;                        // Groundwater table data
extern double SpringCorrection;              // Correction (0-1) for stronger exudation in spring; influences priming and exudate production, if 0, disables spring correction; correction is a factor of 1+SpringCorrection
extern double ExudateFactor;                 // mass fraction of of below-ground production that consists of exudates
extern double SpringFactor;                  // instantaneous value of spring correction, declared global for use by environmental correction SOM decomposition for priming effect
extern double RootSenescence;                // root senescence factor = proportion of root mass that dies during each time step
extern Matrix RootMass;                      // initial root distribution (kg C/m2 in each layer)
extern ofstream *output4;                    // output file rootmass
extern Matrix NewSOM;                        // SOM reservoirs to be changed in each iteration step
extern Matrix ProfileOutput;                 // determines which vertical profiles are sent to log files
extern double BioMass;                       // above ground biomass kg C /m2 (standing crop)
extern Matrix Harvest;                       // harvest dates (1st column) and fraction of biomass harvested (2nd column)
extern int HarvestModel;                     // 1 = harvest amount/dates same every year. 2 = harvest amounts/dates defined in input file per year.
extern Matrix HarvestCorrection;            // Correction of GPP after harvest with reduction factor directly after harvest (first) and period of recovery in days (second)
extern Matrix Grazing;                       // Parts of the year in which garazing occurs, each row is a range of days
extern double BioMassSenescence;             // biomass senescence at each DAY as fraction of above-ground biomass
extern Matrix RespFac;                       // factor of primary production that is respirated during growth
extern double PlantRespiration;              // plant respiration
extern Matrix BioMassRec;                    // storage of biomass, primary production and plant respiration
extern Matrix Manure;                        // Manure application dates (column 1) and quantity (column 2) in kg C/m2/d
extern double ManureFluidFrac;               // fluid fraction of manure
extern Matrix ManureLayers;                  // partitioning of manure among layers; first column: fluids; second column: solids
extern int NrOfSteps;                        // number of time steps
extern Matrix LayerTime;                     // storage matrix for CO2 per layer per timestep
extern Matrix NPPData;                       // net primary production data from file NPPfile
extern char NPPFile[256];                    // file with primary production data for each time step; these data override any selected production model
extern double GrowingDegreeDays;             // growing degree days for phenology and photosynthesis model
extern Matrix TData;                           // Air or soil temperature data from file Tfile
extern Matrix PARData;                         // photosynthetic active radiation or cloud cover data
extern int PARunits; // units in which PARData is given, 0: PAR radiation in umol m-2 s-1; 1:total daily radiation in J cm-2; 2: input is cloud cover
extern double Latitude;                      // Latitude of site in degrees (north positive)
extern Matrix Phenology;                // Phenology data for production model 3 and 4
extern ofstream *output8;                 // output stream npp for production model
extern BOOLEAN LeafSenescence;           // indicates whether leaf senescence may occur for photosynthesis model
extern double DayLength;                  // Daylenght for photosynthesis model
extern double AmbientCO2;                // Ambient CO2 concentration
extern char CO2File[];                 // file with yearly averaged CO2 concentration for multi-year runs
extern Matrix CO2Data;                  // CO2 data for longer model runs with photosynthesis model
extern double LitterLayer;              // organic matter stored in above ground litter layer, in kg C / m2
extern double LitterConversion;         // Conversion factor of daily conversion of above ground to below ground litter at reference temperature Tref; the factor is temperature adjusted such that at 0 degrees the conversion factor is also 0
extern double T_ref;              // reference temperature for correction of decay constants
extern double CurrentLAI;           // leaf area index
extern int StartYear;                 // starting year
extern double Year;                        // current simulation year
extern double KBeer;                    // Beer's law constant for photosynthesis models, values around 0.5
extern Matrix PhotoPar;                 // Parameters for photosynthesis model 5 for tundra, Shaver et al, J. Ecology 2007
extern BOOLEAN Harvested;               // indicates the occurrence of harvest
extern double LAICarbonFraction;         // relates leaf area index to kg C/m2 
extern double PotentialLAI;             // potential LAI based on temperature alone
extern double PreviousLAI;              // difference in LAI between successive time steps to calculate litter production in autumn
extern double HarvestGrazing;             // total of harvest and grazing in one time step
extern double TotalManure;                // total of manure added in one time step
extern double OldLitter;                 // For calculation of storage change of litter layer
extern Matrix CarbonBalance;             // Carbon balance: primary production, C exported, and change in carbon reservoirs in Mol C
extern double GPP;                        // Gros primary production
extern int HdayNr;
extern int HYY;
extern double Harv_height;
extern double HDcount;                  // Counts days after harvest
extern int CalendarYear;
extern double HarvestLitter; // fraction of harvest that is left as litter at each harvest
extern double GreenBiomassRatio;        // Ratio of photosyntesizing biomass to total bioamass for ProductionModel 3
extern double LAIovershoot;                    // Allocates more primary production to below-ground biomass if primary production model causes LAI to overshoot its maximum

void OrgProd();
/* Net Primary Production and ists partitioning among roots and shoots */

double SimpleProd();
/* Simple primary production from yearly sinusoidal function */

double TemperatureProd();
/* Primary production from temperature upper soil layer */

double PARcalc();
/* calculates photosynthetically active radiation or returns data from file*/

double PhotoSynthesis(double LAI, double I);
/* Photosynthesis */


double LAICalc();
/* Calculates LAI from photosynthesis */

double RadProd();
/* Primary production from photosynthetically active radiation */

double TundraProd();
/* Primary production from photosynthetically active radiation, for tundra, Shaver et al, J. Ecology 2007 */

void DoHarvest();
// Harvest according to HarvestModel

//void DoHarvest(double minBiomass);
/* Harvest of biomass at fixed dates */

//void DateHarvest(double minBiomass);
/* Harvest of biomass at selected dates that are different for each year*/

void DoGraze();
/* grazing of biomass */

void AddManure();
/* manure addition    */

void CollectBioMass();
/* collects total Biomass, primary production, respiration, net CO2 flux incl. soil respiration
   in BioMassRec
   Elements of each row:
   1. Day number
   2. Biomass
   3. Primary production
   4. Plant respiration
   5. Net CO2 flux incl. soil respiration
*/

void checkHarvestDate();
