/***************************************************************************
              peatland.h  -  general definitions for PEATLAND 2.0
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

#include <iostream>
#include <fstream>
using namespace std;
//#include "matrix.h"

/******************************DEFINITIONS**********************************/

/**************************** GLOBAL VARIABLES *****************************/

BOOLEAN DoPlot = FALSE;                // Output flags
BOOLEAN Verbose = FALSE;

char ParamFile[256] = "params";         // name of the parameter file
char DataDir[256] = "";                 // data directory
char OutputFilePrefix[256] = "";        // output file prefix

/***************************Model configuration*****************************/
int NrLayers = 15;                      // Number of depth steps
int NrReservoirs = 7;                   // Number of SOM reservoirs
double LayerThickness = 0.1;            // Thickness of depth step
double TStepHeat = 0.1;                 // time step (days) temperature model
double DStepHeat = 0.1;                 // minimum depth step (m) temperature model
double MaxDepthHeat = 6.0;              // maximum depth temperature model (m)
double Timestep = 1.0;                  // model timestep
int NrOfSteps = 365;                    // number of time steps
int StartDay = 0;                       // Julian day nr starting day
int ThermModel = 0;                     // choice of soil thermal model:
                                        // 0 for diffusivity varying with layer physical properties
                                        // 1 for constant diffusivity
                                        // 2 for constant temperature model
int WatertableModel = 1;				// choice for water table model: 0 = read from file, 1 = simple sinusoidal 2 = extended 'Yurova' type model
Matrix Layers;                          // Layer boundaries
double RefLevel;                        // surface datum, set at 0 at the start of the simulation
int StepNr = 1;                         // time step number during iteration
double DayNr;                           // midpoint of simulated timestep, relative to day 1 of the year in which the simulation started
double Timer = 0;                       // starting point of simulated time step
double DayOfTheYear;                    // Julian day number of the midpoint of the simulated time step relative to the current year;
double Year = 0;                        // current simulation year
double StartYear = 0.0;                 // starting year, may be fictive

/***************************Physical properties*****************************/
double DensOrg = 1430;                  // density organic matter in peat (density of wood) kgm-3
double DensMin = 2650;                  // density of mineral matter kg m-3
double DensWater = 1000;                // density of water at 0 degr C
double DensIce = 920;                   // density of ice at 0 degr C
double HCOrg = 2.5e6;                   // volumetric heat capacity organic matter J.m-3.K-1
double HCMiner = 2.0e6;                 // volumetric heat capacity mineral matter
double HCAir = 1.2e3;                   // volumetric heat capacity air (saturated with water vapour)
double HCWater = 4.2e6;                 // volumetric heat capacity water
double HCIce = 1.9257e6;                // volumetric heat capacity ice
double CondOrg = 0.25;                  // thermal conductivity organic matter J.m-1.s-1.K-1
double CondMiner = 2.5;                 // thermal conductivity mineral matter
double CondQuartz = 8.8;                 // thermal conductivity quartz
double CondAir = 0.025;                 // thermal conductivity air
double CondWater = 0.57;                // thermal conductivity water
double CondIce = 2.24;                  // thermal conductivity ice
double CondSnow = 0.3;                  // thermal conductivity snow
Matrix LatentHeat(3);		        // parameters for approximation of temperature-dependent latent heat of fusion of ice J kg-1
double Rgas = 8.314;                    // Gas constant
double MethaneDiff = 1.7280;            // diffusion of methane in air in m2/d
double MethaneDiffWater = 1.7280e-4;    // diffusion of methane in water

/***************************Aerobic SOM decomposition*************************/
double AssimDissim = 2.3;               // Assimiltion/Dissimilation ratio
double ResistFrac = 0.1;                // Fraction of decomposited organic material that is transferred to resistant humus fraction
double MolAct = 74826;                  // molecular activation energy aerobic organic matter decomposition
Matrix Cfrac(7);                        // Carbon fraction (kg/kg) each SOM reservoir
Matrix pFpoints(2,2);                   // Curve for determining environmental correction factor for dryness
double HalfSatPoint = 0.1;              // half activity saturation point for correction factor aeration
double RootAeration = 0;                // Root mass dependent correction (0 - 1) for improved aeration by root growth if 0, this is switched off
double PrimingCorrection = 0;           // Root mass dependent priming effect root exudates on slow C reservoirs; value > 0; if 0, this is switched off
Matrix Kdecay(7);                       // SOM decomposition constants for each reservoir
Matrix KPeatCN(2);                      // constants linear relation of decomposition rate k of peat with CN ratio cf Vermeulen & Hendriks
Matrix SplitRes(5,2);                   // partitions decomposed material between CO2 + microbial biomass (1st column) and resistant SOM
Matrix KPeat;                           // Horizon-C/N dpendent peat decomposition rate
int AnaerobicCO2 = 0;                   // Switch for allowing anaerobic decomposition (sulfate etc) resulting in CO2, if 0 not accounted for
Matrix KAnaerobic(7);                   // Anaerobic decomposition constants, for all SOM reservoirs
Matrix LayerAnaerobic;                  // Anaerobic CO2 per layer
Matrix AnaerobSum;                      // sum of anaerobic CO2 per layer
double Q10Anaerobic = 3.5;              // Q10 anaerobic decomposition
double KLitter = 0.5;                   // decomposition constant above-ground litter and standing dead biomass

/***************************SOM production***********************************/
int ProductionModel = 0;                // Production model:
                                        // 0 for simple sinusoidal function;
                                        // 1 for production dependent on temperature of upper soil layer
                                        // 2 for production data read from file, variable NPPFile has to be defined
                                        // 3 for photosynthesis model, PAR data supplied in PARFile (W m-2)
                                        // 4 for photosynthesis model, cloud cover data supplied in PARFile (fractions 0 - 1), PAR calculated, variable Latitude defined
                                        // 5 for photosynthesis model for tundra, Shaver et al, J. Ecology 2007 (not depending on CO2 concentration) PAR data supplied in PARFile (W m-2)
                                        // 6 for photosynthesis model for tundra, Shaver et al, J. Ecology 2007 (not depending on CO2 concentration) cloud cover data supplied in PARFile
double Latitude = 50.0;                 // Latitude of site in degrees (north positive)
double GrowingDegreeDays = 0.0;         // growing degree days for phenology and photosynthesis models
Matrix Phenology(8);                       // Phenology parameters; first number is type of phenology (0 for evergreen, 1 for summergreen)
                                        // second is the base for calculating the heat sum (growing degree days)
                                        // third the heat sum when maximum leaf area index
                                        // is reached, fourth the maximum LAI
double KBeer = 0.5;                     // Beer's law constant for photosynthesis models, values around 0.5
Matrix PhotoPar(4);                     // Parameters for photosynthesis model 5 for tundra, Shaver et al, J. Ecology 2007
                                        // 1: PlantResp0 Plant respiration at zero degrees 0.4-1.5
                                        // 2: Temperature sensitivity factor plant respiration 0.03 - 0.07
                                        // 3: light-saturated photosynthetic rate per unit leaf area (μmol m–2 leaf s–1) val. 10-20
                                        // 4: initial slope of the light response curve (μmol CO2 μmol–1 photons) val 0.002-0.07
double ShootsFactor = 0.5;              // mass fraction root growth against shoot growth
Matrix RespFac(2);                      // factor of primary production that is respirated during growth
Matrix ProdTFunc(2);                    // determines temperature dependent production rate; 1st number is minimum, 2nd optimum
double SatCorr = 0.0;                   // correction of production for saturation of topsoil, depresses production at high saturation, switched off when 0
double GrowFuncConst = 1.0;               // proportionality constant growth functioen - primary productivity for plant transport in methane model
double SpringCorrection = 0;            // Correction (0-1) for stronger exudation in spring; influences priming and exudate production, if 0, disables spring correction; correction is a factor of 1+SpringCorrection
double MaxProd = 0.0057;                // Maximum primary productivity (kgC/m2/day)
double MinProd = 0.0;                   // Minimum primary productivity
double MaxRootDepth = 0.3;              // Maximum root depth (m)
int NoRootsBelowGWT = 1;                // if 1, no roots will grow below groundwater table (No telmatophytes)
double RootLambda = 20;                 // decay rate exponential root distribution function
double RootSenescence = 0.003;          // root senescence factor = proportion of root mass that dies during each time step
double InitRoots = 0.5;                 // Initial root mass in all layers
Matrix RootMass;                        // initial root distribution (kg C/m2 in each layer)
double ExudateFactor = 0.2;             // mass fraction of of below-ground production that consists of exudates
double BioMass = 0.3;                   // above ground biomass kg C /m2 (standing crop)
double BioMassSenescence = 0.001;       // biomass senescence at each DAY as fraction of above-ground biomass
Matrix Harvest;                         // harvest dates (1st column) and fraction of biomass harvested (2nd column)
BOOLEAN Harvested = TRUE;               // indicates the occurrence of harvest
double LAICarbonFraction = 0.1;         // relates leaf area index to kg C/m2 (kg C/m2 per unit of LAI)
Matrix Grazing;                         // Parts of the year in which garazing occurs, each row is a range of days
                                        // followed by the amount of biomass removed (kg C m2/day and the amount of excretion (kg C m2/day)
Matrix Manure;                          // Manure application dates (column 1) and quantity (column 2) in kg C/m2/d
double ManureFluidFrac = 0.75;          // fluid fraction of manure
Matrix ManureLayers;                    // partitioning of manure among layers; first column: fluids; second column: solids
Matrix RootDistrib;                     // root distribution function
char NPPFile[256] = "";                 // file with primary production data for each time step; these data override any selected production model
Matrix NPPData;                         // net primary production data from file NPPfile or model
char PARFile[256] = "";               // file with PAR (photosynthetic active radiation, (J cm-2) or cloud cover data (cloud fraction, 0-1) for NPP model
Matrix PARData;                         // photosynthetic active radiation or cloud cover data
BOOLEAN LeafSenescence = FALSE;           // indicates whether leaf senescence may occur for photosynthesis model
double AmbientCO2 = 385;                // Ambient CO2 concentration
char CO2File[256] = "";                 // file with yearly averaged CO2 concentration for multi-year runs
Matrix CO2Data;                         // CO2 data for longer model runs with photosynthesis model
double DayLength = 12;                  // Daylenght for photosynthesis model
double LitterLayer = 0.0;               // organic matter stored in above ground litter layer, in kg C / m2
double LitterConversion;                // Conversion factor of daily conversion of above ground to below ground litter at reference temperature Tref; the factor is temperature adjusted such that at 0 degrees the conversion factor is also 0

/************************* Methane model ************************************/

Matrix MethaneReservoirs(7);            // This parameter tells which reservoirs are summed for calculating methane production from easily decomposed
                                        // organic matter reservoirs. It allows to delete or reduce some reservoirs; peat always should be deleted
double MethaneReservoirSum;             // sum of methane reservoirs
double MethaneR0 = 0.4;                 // Methane production rate factor for fresh organic C microM/h
                                        // SITE-SPECIFIC PARAMETER, CAN BE USED TO TUNE THE MODEL
                                        // Walther & Heimann 2000: values between 0.3 and 0.6 at high latitude sites and 2.8 at tropical site
double MethanepHCorr = 0.1;             // for every PH unit lower or higher than neutral, MethanepHCorr*R0 is added to R0
double MethaneQ10 = 6.0;                // Q10 value for temperature correction methane production; range 1.7 - 16 ref. in Walther & Heimann 2000
double MethaneOxQ10 = 1.4;              // Q10 value for temperature correction methane oxidation; range 1.4 - 2.1, ref. in Walther & Heimann 2000
double MethaneVmax = 5;                 // Vmax Michaelis-Menten eq methane oxidation micrMol/hr range 5-50
double MethaneKm = 5;                   // Km Michaelis-Menten eq methane oxidation micrMol range 1-5
double MethaneMaxConc = 500;            // maximum methane concentration in pore water
double MethaneERateC = 1;               // Ebullition rate constant (1/hr)
double MethanePRateC = 0.01;            // Rate constant for plant transport of methane (1/hr)
double MethaneAir = 0.1072;             // Methane concentration in the atmosphere - microMol
double MethaneTRef = 10;                // Reference temperature for temperature sensitivity methane production
                                        // if negative, the average yearly temperature is selected as reference temperature
double MethanePType = 15;               // Vegetation type factor for gas transport by plants range: 0-15
double MethanePlantOx = 0.9;            // Fraction of methane that is oxidized during transport in plants
double PartialAnaerobe;                 // Determines the slope of the relation of partial anaerobe soil fraction above the water table to soil saturation, >1
double AnaerobeLagFactor = 0.01;	    // Determines time lag for development of sufficiently anaerobic conditions after saturation of a layer
double CO2CH4ratio = 0.5;               // Molar ratio between CH4 and CO2 production; for acetate splitting this is 0, for CO2 reduction 0.5
Matrix InitMethane;                     // Initial methane concentration profile

/******************************Water table and temperature*********************/

double DayMinGW = 220;                  // day of minimum groundwater table
double MinGW = 0.25;                    // lowest water table level (m below surface)
double AmplitudeGW = 0.20;              // amplitude of water table movement
char GwFile[256] = "";                  // file where the groundwater table time series is stored; if empty a sinusoidal time series is assumed
char SoilMoisture[256] = "";            // file where soil moisture profile time series is stored; if empty the soil moisture will be calculated using very simplified assumptions
char PrecipFile[256] = "";              // file where the precipitation time series is stored; if empty the water table is assumed to be read from file
char EvapFile[256] = "";                // file where the evaporation time series is stored; if empty the water table is assumed to be read from file
double T_average = 10.5;                // average yearly temperature
double T_amplitude = 7.75;              // amplitude of temperature throughout the year
double T_ref = 11;                      // reference temperature for correction of decay constants
double ThermDiff = -1.0;                // thermal diffusivity, if negative it is estimated from soil properties
                                        // this parameter is not used if ThermModel == 0
char TFile[256] = "";                   // file with temperature time series; if empty string, a sinusoidal time series will be defined
Matrix TData;                           // Air or soil temperature data from file Tfile
Matrix Precipitation;					// Precipitation data (read from file)
Matrix Evaporation;						// Evaporation data (read from file)
Matrix GwData;                          // Groundwater table data from GWFile
Matrix T_init;                          // Initial temperature profile
Matrix MoistProfiles;                   // Matrix with soil moisture profiles from observational data or other model
double MaxSnowdepth = 0.0;              // Maximum snow depth
double DayMaxSnowdepth = 60;            // Day of maximum snow depth
double SnowMeltrate =1.0;               // Rate of snowmelt per degree C above zero per day
char SnowFile[256] = "";                // file with temperature time series; if empty string, a sinusoidal time series will be defined
Matrix SnowData;                        // Snow depth data from file Snowfile
double WatertableInit = 0.0;			// Initial water table in m, has to be specified if the watertable is calculated by the model
double EvapCorrection = 0.5;			// Correction factor to reduce evaporation if water table is below surface, for water table model
double RunoffThreshold = 0.1;			// Threshold above which a ponded water layer produces runoff; for water table model
double OpenWaterFactor = 1.0;			// evaporation correction factor for open water evaporation
double CropFactor = 1.0;			    // Makkink Crop factor to correct evaporation for vegetation properties; for water table model
double PondedWater = 0.0;				// layer of ponded water above the soil
double DrainageDist = 1.0;              // Distance to nearest drainage channel (m)
double Ksat = 0.01;                     // Saturated hydraulic conductivity of soil
double DrainageDepth = 1.5;             // depth of aquifer
double DrainLevel = 0.0;                // reference level of water in the drains/river channel with respect to top of soil surface (m)
Matrix DrainData;                       // Matrix to store variable drain water level data
char DrainageFile[256]="";              // file name for drainage water level data
char RunOnFile[256]="";                 // file name for runon data
Matrix RunOn;                           // Matrix with Run-on water quantity
double VegTScalingFactor = 1.0;			// scaling factor for air to soil surface temperature; put to one if soil surface temperature is input; outherwise a valye between 0.6 and 1.0

/******************************Soil profile data*******************************/

char SoilProfile[256] = "";             // Name of the file where the soilprofile data are stored
int NrHorizons;                         // number of horizons
Matrix Horizons;                        // Horizon base depths with respect to surface
Matrix CNRatio;                         // CN ratios for eac soil layer; the decomposition of peat can be made dependent on these
Matrix DBD;                             // Dry bulk density for each horizon
Matrix PercOrg;                         // Percentage organic matter for each horizon
Matrix SandFraction;					// sand weight fraction of mineral fraction
Matrix ClayFraction;					// clay weight fraction of mineral fraction, influences decomposition rate humus reservoir
Matrix Layer_pH;                        // pH
Matrix InitRes;                         // initial contents of SOM reservoirs in % of total SOM including peat
Matrix Layer_pF;                        // pF curves, 2 possible definitions:
                                        // as phi values for potentials defined in pFVal
                                        // or Van genuchten parameters theta_r, theta_s, alfa, l and n
                                        // In the latter case the rows in Layer_pF have only 5 elements and pFVal only 1 column
Matrix pFVal;                           // suction potentialsfor pF curves
Matrix Porosity;                        // porosity of soil layers, can be specified or calculated from organic matter percentage and dry bulk density
Matrix pFCurves;                        // final pF Curves after conversion from van genuchten parameters
Matrix FreezingCurve;                   // Unfrozen water content curve at below zero temperatures
                                        // Relation: f = f_inf + 1/(b-T)^a, where a is a constant given in FreezingCurve
                                        // and b is determined by a and maximum porosity. a has values between 1.5 and 2 (the larger values for sand, the smaller for clay
                                        // f_inf is the weight % unfrozen water content at wilting point ussuming that water with a suction below wilting point does not freeze
Matrix Freeze;                          // contains for each model layer the f_inf, a and b parameters
Matrix UnFrozen;                        // Unfrozen water content (kg water / kg dry soil)
Matrix Ice;                             // Ice content (kg ice / kg dry soil)

/******************************Output data*************************************/

Matrix TotalReservoir;                  // totals C per reservoir per layer
Matrix ReservoirTime;                   // storage matrix for CO2 per reservoir per timestep ; 1st element: day number
Matrix LayerTime;                       // storage matrix for CO2 per layer per timestep ; 1st element: day number
Matrix TotalMethane;                    // storage matrix for CH4 fluxes ; 1st element: day number
Matrix BioMassRec;                      // storage of biomass, primary production and plant respiration
/* collects total Biomass, primary production, respiration, net CO2 flux incl. soil respiration
   in BioMassRec
    1 DayNr;
    2 BioMass + totalroots;                     
    3 Primary Production;
    4 Plant Respiration;
    7 LitterLayer;
    8 Harvest+Grazing;
    9 Manure;
    10 LAI;
*/
Matrix CarbonBalance;                   // Carbon balance: primary production, C exported, and change in carbon reservoirs in Mol C
/*
 * 1: primary production
 * 2: carbon inputs from manure
 * 3: carbon inputs from livestock excretion
 * 4-10: Changes in each soil carbon reservoir
 * 11: total CO2-C emission fom aerobic decomposition
 * 12: CO2-C from above-ground litter decomposition
 * 13: total CH4-C emission
 * 14: CO2-C from CH4 oxidation
 * 15: CO2-C from plant respiration
 * 16: (not yet implemented) storage change (anaerobic) CO2 in soil water
 * 17: storage change CH4 in soil water 
 * 18: (not yet implemented) CO2 - and CH4-change export by groundwater flow
 * 19: harvest
 * 20: grazing
 * 21: above-ground biomass senescence
 * 22: below-ground biomass senescence
 * 23: storage change above-ground biomass
 * 24: storage change below ground biomass
 * 25: balance sum
 */


/*****************************Output Log files*********************************/

Matrix ProfileOutput;                    // determines which vertical profiles are sent to log files
ofstream *output1;                       // output log files
ofstream *output2;
ofstream *output3;
ofstream *output4;
ofstream *output5;
ofstream *output6;
ofstream *output7;
ofstream *output8;
ofstream *output9;
ofstream *output11;
ofstream *output12;
ofstream *output13;
ofstream *output14;
ofstream *output15;
ofstream *output16;
ofstream *output17;

/*****************************Intermediary variables***************************/

Matrix MethProfile;                      // methane concentration profile (mol/m3 in soil pore volume, including bubbles)
Matrix OxconCH4;                         // oxygen consumption methane oxidation
                                         // NB: output of methane model - presently not used in the model
Matrix InitSOM;                          // Initial carbon content (kg C per layer) in each SOM reservoir
Matrix MoistTheta;                       // soil moisture profile for timestep
Matrix MatricPotential;                  // matric potential of model soil layer
int NrHeatLayers = 0;                    // number of layers temperature model
Matrix TProfile;                         // temperature profile soil thermal submodel
Matrix SoilTemp;                         // soil temperatures model layers, interpolated from TProfile
Matrix PoreVol;                          // porosity each layer
Matrix Saturation;                       // pore volume saturation with water
Matrix OldSat;                           // saturation values at previous time step
Matrix LastSatTime;                      // Time after last complete saturation of layer in days
Matrix ThermDiffVar;                     // thermal diffusivity layer dependent diffusivity
Matrix HeatLayers;                       // layer midpoints layers thermal model
Matrix MethaneR0Corr;                    // pH dependent methane production rate
double TotalPrimProd = 0;                // total primary production
Matrix NewSOM;                           // SOM reservoirs to be changed in each iteration step (kg C per layer)
Matrix PeatDecay;                        // logs true loss of peat matrix
Matrix ResYearSOM;                       // logs yearly change of all reservoirs in all layers
int ManureCount = 0;                     // counter manure additions
double PrimProd;                         // Primary production per time step
double Shoots;                           // Shoot production per time step
double SpringFactor;                     // instantaneous value of spring correction, declared global for use by environmental correction SOM decomposition for priming effect
double PlantRespiration;                 // plant respiration
Matrix CorrFac;                          // environmental correction factors for aerobic SOM decomposition constants k
double CurrentGW;                        // Current water table during time step
double PeatLoss = 0;                     // Totalized loss of peat C over all layers
Matrix CO2;                              // CO2 evolved from each layer and reservoir
Matrix MethaneFlux;                      // Methane flux at surface;  1st value: ebullition flux; 2nd: plant mediated flux; 3d:diffusive flux
int TopSat;                              // Index of first saturated layer
Matrix CO2FromMethane;                   // CO2 from methane oxidation
double SnowDepth;                        // Actual snow heigth (m)
double SnowStorage = 0.0;				 // water storage in snow (m)
double SnowStartDay;                     // day at wich snow accumulation has started
double HeatCondTop;                      // heat conductivity of the top soil
double BelowWTStorage;					 // Below water table water storge for water table model; this excludes ponded water with a positive water table
double AboveWTStorage;					 // Above water table water storge for water table model
double MaxBelowWTStorage;			     // Maximum Below water table water storge for water table model
double MinAboveWTStorage;				 // Minimum Above water table water storge for water table model
BOOLEAN ApplyVanGenuchten = FALSE;		 // indicates application of Van Genuchten equation for soil moisture
double FrozenStorage = 0;				 // Frozen water storage for water table model
double CurrentFrozen;					 // Current top of frozen layer for water table model
double CurrentGWBase;					 // current lowest water table level for water table model
double FrozenWatertable;				 // Water table in frozen part of soil
double CurrentLAI = 0.0;                 // leaf area index at timestep
double PotentialLAI = 0.0;               // potential leaf area index at timestep t, if nothing is removed by grazing or harvest
double PreviousLAI = 0.0;              // difference in LAI between successive time steps to calculate litter production in autumn
double HarvestGrazing = 0.0;             // total of harvest and grazing in one time step
double TotalManure = 0.0;                // total of manure added in one time step


