/***************************************************************************
             initialize.h  -  initialization functions of PEATLAND
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


#define INIT_ERROR1 "Insufficient data in initial temperature profile: "
#define INIT_ERROR2 "Cannot open file: "
#define INIT_ERROR3 "Negative methane production rate, check MethaneR0 and pH correction"
#define OUTPUT1     "temperature.dat"
#define OUTPUT2A    "moisture.dat"
#define OUTPUT2B    "wfps.dat"
#define OUTPUT2C    "aerobfraction.dat"
#define OUTPUT3     "methane.dat"
#define OUTPUT4     "roots.dat"
#define OUTPUT5     "labileSOM.dat"
#define OUTPUT6     "ice.dat"
#define OUTPUT7     "watertable.dat"
#define OUTPUT8     "npp.dat"
#define OUTPUT9     "totalSOM.dat"
#define OUTPUT10    "anaerobicCO2reservoirs.dat"
#define OUTPUT11    "peat.dat"
#define OUTPUT12    "liquid_manure.dat"
#define OUTPUT13    "solid_manure.dat"
#define OUTPUT14    "exudates.dat"
#define OUTPUT15    "litter_roots.dat"
#define OUTPUT16    "microbes.dat"
#define OUTPUT17    "humus.dat"


/************************ GLOBAL VARIABLES ***********************************/

extern char OutputFilePrefix[];           // output file prefix;
extern int NrLayers;                      // Number of depth steps /layers
extern double LayerThickness;             // Thickness of layers
extern Matrix Layers;                     // Layer boundaries
extern double RefLevel;                   // Position of base of soil profile with respect to reference level
                                          // at the start of the simulation
extern int NrHorizons;                    // number of horizons
extern Matrix Horizons;                   // Horizon base depths with respect to surface
extern int NoRootsBelowGWT;               // if 1, no roots will grow below groundwater table (No telmatophytes)
extern double RootLambda;                 // decay rate exponential root distribution function
extern char GwFile[];                     // file where the groundwater table time series is stored; if empty a sinusoidal time series is assumed
extern int WatertableModel;				  // choice for water table model: 0 = read from file, 1 = simple sinusoidal 2 = extended 'Yurova' type model
extern double MaxRootDepth;               // Maximum root depth (m)
extern double MinGW;                      // lowest water table level (m below surface)
extern Matrix RootDistrib;                // root distribution function
extern Matrix RootMass;                   // initial root distribution (kg C/m2 in each layer)
extern double InitRoots;                  // Initial root mass in all layers
extern Matrix TotalReservoir;             // totals C per reservoir per layer
extern Matrix ReservoirTime;              // storage matrix for CO2 per reservoir per timestep
extern Matrix AnaerobReservoirTime;       // storage matrix for anaerobic CO2 per reservoir per timestep ; 1st element: day number
extern Matrix LayerTime;                  // storage matrix for CO2 per layer per timestep
extern int NrReservoirs;                  // Number of SOM reservoirs
extern int NrOfSteps;                     // number of time steps
extern Matrix MethProfile;                // methane concentration profile
extern Matrix InitMethane;                // Initial methane concentration profile
extern Matrix OxconCH4;                   // oxygen consumption methane oxidation
extern Matrix TotalMethane;               // storage matrix for CH4 results
extern Matrix InitSOM;                    // Initial carbon content (kg C per layer) in each SOM reservoir
extern Matrix DBD;                        // Dry bulk density for each horizon
extern Matrix PercOrg;                    // Percentage organic matter for each horizon
extern Matrix InitRes;                    // initial contents of SOM reservoirs in % of total SOM including peat
extern Matrix Cfrac;                      // Carbon fraction (kg/kg) each SOM reservoir
extern double DensOrg;                    // density organic matter in peat (density of wood) kgm-3
extern double DensMin;                    // density of mineral matter kg m-3
extern double DensWater;                  // density of water at 0 degr C
extern double DensIce;                    // density of ice at 0 degr C
extern double HCOrg;                      // volumetric heat capacity organic matter J.m-3.K-1
extern double HCMiner;                    // volumetric heat capacity mineral matter
extern double HCAir;                      // volumetric heat capacity air (saturated with water vapour)
extern double HCWater;                    // volumetric heat capacity water
extern double CondOrg;                    // thermal conductivity organic matter J.m-1.s-1.K-1
extern double CondMiner;                  // thermal conductivity mineral matter
extern double CondAir;                    // thermal conductivity air
extern double CondWater;                  // thermal conductivity water
extern Matrix Porosity;                   // porosity of soil layers
extern double ThermDiff;                  // thermal diffusivity, if not defined it is estimated from soil properties
extern Matrix TData;                      // Air temperature data from file Tfile
extern Matrix SoilTData;                  // Air temperature data from file SoilTFile
extern Matrix GwData;                     // Groundwater table data from GWFile
extern double T_average;                  // average yearly temperature
extern double T_amplitude;                // amplitude of temperature throughout the year
extern int StartDay;                      // Julian day nr starting day
extern double Timestep;                   // model timestep
extern double DayMinGW;                   // day of minimum groundwater table
extern double MinGW;                      // lowest water table level at this day (m below surface)
extern double AmplitudeGW;                // amplitude of water table movement
extern double TStepHeat;                  // time step (days) temperature model
extern double DStepHeat;                  // minimum depth step (m) temperature model
extern double MaxDepthHeat;               // maximum depth temperature model (m)
extern int NrHeatLayers;                  // number of layers temperature model
extern Matrix TProfile;                   // temperature profile soil thermal submodel
extern Matrix SoilTemp;                   // soil temperatures model layers, interpolated from TProfile
extern Matrix T_init;                     // Initial temperature profile
extern Matrix PoreVol;                    // Porosity each layer
extern Matrix ThermDiffVar;               // thermal diffusivity for layer dependent diffusivity
extern Matrix HeatLayers;                 // layer midpoints layers thermal model
extern int ThermModel;                    // choice of soil thermal model
extern char DataDir[];                    // data directory
extern ofstream *MethAll;                 // log file methane profiles
extern double ResistFrac;                 // Fraction of decomposited organic material that is transferred to resistant humus fraction
extern double DissimAssimRatio;           // Assimiltion/Dissimilation ratio
extern double AnaerobicDARatio;          // Assimiltion/Dissimilation ratio anaeroob
extern double LitterLayer;              // organic matter stored in above ground litter layer, in kg C / m2
extern double OldLitter;                 // For calculation of storage change of litter layer
extern Matrix SplitRes;                   // partitions decomposed material between CO2 + microbial biomass (1st column) and resistant SOM
extern Matrix CNRatio;                    // CN ratios for eac soil layer; the decomposition of peat can be made dependent on these
extern Matrix Kdecay;                     // SOM decomposition constants for each reservoir
extern Matrix KPeatCN;                    //constants linear relation of decomposition rate k of peat with CN ratio cf Vermeulen & Hendriks.
                                          // First number is the reference C/N value. Second is the slope of the relative decrease within the range of 10-55 C/N.
extern Matrix Layer_pH;                   // pH
extern Matrix MethaneR0Corr;              // pH dependent methane production rate
extern double MethaneR0;                  // Methane production rate factor for fresh organic C microM/h
extern double MethaneAir;               // Methane concentration in the atmosphere specified in ppmv, for calculations inside peatland in moles in a layer of 10 cm above surface
extern double MethaneTRef;              // Reference temperature for temperature sensitivity methane production
extern double MethanepHCorr;              // for every PH unit lower or higher than neutral, MethanepHCorr*R0 is added to R0
extern Matrix BioMassRec;                 // storage of biomass, primary production and plant respiration
//extern ofstream *LabileSOM;               // logs labile C reservoirs, helpful to detect non-steady state behaviour
//extern ofstream *TotalSOM;                // logs all C reservoirs, helpful to detect non-steady state behaviour
extern Matrix NewSOM;                     // SOM reservoirs to be changed in each iteration step
extern Matrix OldSOM;
extern Matrix PeatDecay;                  // logs true loss of peat matrix
extern Matrix ProfileOutput;              // determines which vertical profiles are sent to log files
extern ofstream *output1;                 // output log files
extern ofstream *output2a;
extern ofstream *output2b;
extern ofstream *output2c;
extern ofstream *output3;
extern ofstream *output4;
extern ofstream *output5;
extern ofstream *output6;
extern ofstream *output7;
extern ofstream *output8;
extern ofstream *output9;
extern ofstream *output10;
extern ofstream *output11;
extern ofstream *output12;
extern ofstream *output13;
extern ofstream *output14;
extern ofstream *output15;
extern ofstream *output16;
extern ofstream *output17;
extern double DayNr;                      // midpoint of simulated timestep
extern int YEARdays;                          // Length of one year in days
extern int YearCounter;                       // current simulation year
extern int CalendarYear;
extern int Month;                       // Current month of simulation

extern int Verbose;                       // output flag for run information to console
extern double DayOfTheYear;               // Julian day number of the midpoint of the simulated time step relative to the current year;
extern double Year;                       // current simulation year
extern double Timer;                      // starting point of simulated time step
extern Matrix CorrFac;                    // environmental correction factors for aerobic SOM decomposition constants k
extern Matrix MatricPotential;            // matric potential of model soil layer
extern Matrix CO2;                        // CO2 evolved from each layer and reservoir

extern double MethaneMaxConc;               // maximum methane concentratio
extern Matrix MethaneFlux;                // Methane flux at surface;  1st value: ebullition flux; 2nd: plant mediated flux; 3d:diffusive flux
extern Matrix MethaneReservoirs;          // This parameter tells which reservoirs are summed for calculating methane production from easily decomposed
extern double MethaneReservoirSum;

extern Matrix FreezingCurve;              // Unfrozen water content curve at below zero temperatures
extern Matrix pFCurves;                   // final pF Curves after conversion from van genuchten parameters
extern Matrix pFVal;                      // suction potentials for pF curves
extern Matrix Freeze;                     // contains for each model layer the f_inf, a and b parameters
extern Matrix UnFrozen;                   // Unfrozen water content (weight basis)
extern Matrix Ice;                        // Ice content (kg ice / kg dry soil)
extern Matrix MoistTheta;                 // soil moisture profile for timestep
extern double CurrentGW;                  // Current water table during time step
extern double WatertableInit;			// Initial water table in m, has to be specified if the watertable is calculated by the model

extern int StartYear;                 // starting year
extern char RunOnFile[];                 // file name for runon data
extern Matrix RunOn;                     // Matrix with Run-on water quantity
extern int AnaerobicCO2;                   // Switch for allowing anaerobic decomposition (sulfate etc) resulting in CO2, if 0 not accounted for
extern Matrix LayerAnaerobic;             // Anaerobic CO2 per layer
extern Matrix AnaerobSum;          // sum of anaerobic CO2 per layer
extern Matrix AnaerobSumRes;      // sum of anaerobic CO2 per reservoir
extern double MethaneVmax;              // Vmax Michaelis-Menten eq methane oxidation micrMol/hr range 5-50
extern double MethaneKm;                // Km Michaelis-Menten eq methane oxidation micrMol range 1-5
extern Matrix CarbonBalance;             // Carbon balance: primary production, C exported, and change in carbon reservoirs in Mol C

/************************ FUNCTION HEADERS ************************************/


void MakeLayers();
/* Initializes model layers and interpolates data from the soil profile description
   Initalizes the pore volume array
   Initializes the soil top reference level                                                        */

void RootsInit();                          // calculates root distribution

void OutputInit();                         // initializes output arrays

void SOMResInit();                         // calculates organic C (kg C per layer) per SOM reservoir per layer

void InitHeat();
/* computes heat capacity and thermal conductivity from organic matter content and pore volume
 cf Hillel (1998) and Luckner & Schestakow (1991)                                           */

double thermcond(double f, double l1, double l2);
/*  computes thermal conductivity for for a fractional mixture of two soil constituents
    assuming absence of layering, according to eq 1.108c in Luckner & Schestakow 1991
    Used in function InitHeat
    Input:
    f	: fraction of the material with lowest conductivity
    l1	: conductivity corresponding with this material
    l2	: conductivity of the other stuff
    Output:
    lambda: thermal conductivity of the mixture                                              */

void InitTime();           // initializes time system

void InitTseries();
// initializes temperature and grondwater table time series if they have not been read from file

void InitWater();
// initializes water table model 

void InitLogFiles();
/*initializes files for logging model state variables  */

void CloseLogFiles();                            // closes all log files

void  InitDecomp();
/* Corrects aerobic decomposition constants based on assimilation ratio
Decompositon constant for peat is corrected by C/N ratio of layers */

void InitMethaneModel();
/* Initializations for methane model */

