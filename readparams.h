/***************************************************************************
          readparams.h  -  header file for parameter read functions
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

/****************************************************************************
Readparams contains functions for reading model parameters from the fixed and
variable parameter files.
Unlike the previous version of PEATLAND, the soil profile parameters are not
defined in the general parameter file that is specified on the command line.
Soil profile descriptions are stored in a seperate file, the name of this file
is specified in the general parameter file

The structure of the parameter files is defined as follows:

lines starting with '%' or anything after '%' within a line are comments
(cf Octave and Matlab)
a parameter identifier starts a section where the data of a parameter is stored
this is followed with one or more lines defining the data

Basic data read functions are:

readsingle  - reads a single value parameter
readarray   - reads an array
readmatrix  - reads a matrix

*/


#define OUTPUT_METHANE    "methanefluxes.dat"
#define OUTPUT_CO2RES     "CO2reservoirs.dat"
#define OUTPUT_CO2LAY     "CO2layers.dat"
#define OUTPUT_BIO        "biomass.dat"
#define OUTPUT_WTABLE     "watertable.dat"
#define OUTPUT_ANAEROB    "anaerobCO2.dat"
#define OUTPUT_BALANCE    "carbonbalance.dat"
#define OUTPUT_PEATDECOMP "peatdecomposition.dat"

#define ERRORMSG1     "Cannot open file"
#define ERRORMSG2     "Array/matrix read: there seem to be more items then specified for variable "
#define ERRORMSG3     "Array/matrix read: there are less items then specified for variable "
#define ERRORMSG4     "Array/matrix read: errors found while reading matrix entries in variable "
#define ERRORMSG5     "pF matrix needs either 5 columns or the same number of columns as variable pFVal"
#define ERRORMSG6     "Memory allocation error data read buffer"
#define ERRORMSG7     "Maxinum number of horizons exceeded: should be < "
#define ERRORMSG8     "Time series read error - less data on line than expected."
#define ERRORMSG9     "Time series read error from file "
#define ERRORMSG10    ": Unknown command line option."
#define ERRORMSG11    "Cannot open file: "
#define ERRORMSG12    "Start or end dates are not correct."
#define ERRORMSG13    "Number of columns in Harvest file must be 4"
#define ERRORMSG14    "First harvest before start date. Please change."
#define MSG1          "?             : displays command line options"
#define MSG2          "v             : verbose - displays run information on console during model run"
#define MSG3          "dir=directory : working directory for input and output files"
#define MSG4          "par=filename  : input parameter file, default: params"
#define MSG5          "out=filename  : output file prefix"
#define MSG6          "Parameter: "
#define MSG7          "\tfrom file named: "


/* Extern declarations of global variables; definiton and initialization in peatland.h */

extern BOOLEAN DoPlot;                   // Output flags
extern BOOLEAN Verbose;
extern char ParamFile[];                 // Parameter file name
extern char DataDir[];                   // Data directory
extern char OutputFilePrefix[];          // output file;

extern int NrLayers;                     // Number of depth steps
extern int NrReservoirs;                 // Number of SOM reservoirs
extern double LayerThickness;            // Thickness of depth step
extern double TStepHeat;                 // time step (days) temperature model
extern double DStepHeat;                 // minimum depth step (m) temperature model
extern double MaxDepthHeat;              // maximum depth temperature model (m)

extern double Timestep;                  // model timestep
extern int NrOfSteps;                    // number of time steps
extern int StartDay;                     // Julian day nr starting day
extern int StartMonth;     
extern int StartYear;                 // starting year
extern int EndDay;
extern int EndMonth;
extern int EndYear;   

extern int HDD;
extern int HMM;
extern int HYY;
extern int Harvest_LOC;
extern char HarvestFile[];
extern Matrix HarvestData;
extern int Harvest_LOC_END;
extern int HarvestModel;
extern int PARunits; // units in which PARData is given, 0: PAR radiation in umol m-2 s-1; 1:total daily radiation in J cm-2; 2: input is cloud cover

extern char StartDate[];                // Start date DD/MM/YYYY
extern char EndDate[];                  // End date DD/MM/YYYY
extern int YEARdays;                    // Length of one year in days
extern int NrOfYears;                   // Number of years in simulation  

extern int ThermModel;                   // choice of soil thermal model:
extern int WatertableModel;				 // choice for water table model: 0 = read from file, 1 = simple sinusoidal 2 = extended 'Yurova' type model
extern double DensOrg;                   // density organic matter in peat (density of wood) kgm-3
extern double DensMin;                   // density of mineral matter kg m-3
extern double DensWater;                 // density of water at 0 degr C
extern double DensIce;                   // density of ice at 0 degr C
extern double HCOrg;                     // volumetric heat capacity organic matter J.m-3.K-1
extern double HCMiner;                   // volumetric heat capacity mineral matter
extern double HCAir;                     // volumetric heat capacity air (saturated with water vapour)
extern double HCWater;                   // volumetric heat capacity water
extern double HCIce;                     // volumetric heat capacity ice
extern double CondOrg;                   // thermal conductivity organic matter J.m-1.s-1.K-1
extern double CondMiner;                 // thermal conductivity mineral matter
extern double CondAir;                   // thermal conductivity air
extern double CondWater;                 // thermal conductivity water
extern double CondIce;                   // thermal conductivity ice
extern double CondSnow;                  // thermal conductivity snow
extern double CondQuartz;                 // thermal conductivity quartz
extern double Rgas;                      // Gas constant
extern double MethaneDiff;               // diffusion of methane in air in m2/d
extern double MethaneDiffWater;          // diffusion of methane in water

extern double DissimAssimRatio;          // Assimiltion/Dissimilation ratio
extern double ResistFrac;                // Fraction of decomposited organic material that is transferred to resistant humus fraction
//extern double MolAct;                    // molecular activation energy aerobic organic matter decomposition
extern Matrix Cfrac;                     // Carbon fraction (kg/kg) each SOM reservoir
extern Matrix pFpoints;                  // Curve for determining environmental correction factor for dryness
extern double HalfSatPoint;              // half activity saturation point for correction factor aeration
extern double RootAeration;              // Root mass dependent correction (0 - 1) for improved aeration by root growth if 0, this is switched off
extern double PrimingCorrection;         // Root mass dependent priming effect root exudates on slow C reservoirs; value > 0; if 0, this is switched off
extern Matrix Kdecay;                    // SOM decomposition constants for each reservoir

extern Matrix AerobicQ10;                // Q10 for each reservoir
extern double AnaerobicDARatio;          // Assimiltion/Dissimilation ratio anaeroob

extern Matrix KPeatCN;                   // constants linear relation of decomposition rate k of peat with CN ratio cf Vermeulen & Hendriks
extern double ShootsFactor;              // mass fraction root growth against shoot growth
extern Matrix RespFac;                   // factor of primary production that is respirated during growth
extern Matrix ProdTFunc;                 // determines temperature dependent production rate; 1st number is minimum, 2nd optimum
extern double SatCorr;                   // correction of production for saturation of topsoil, depresses production at high saturation, switched off when 0
extern double GrowFuncConst;             // proportionality constant growth functioen - primary productivity for plant transport in methane model
extern double SpringCorrection;          // Correction (0-1) for stronger exudation in spring; influences priming and exudate production, if 0, disables spring correction; correction is a factor of 1+SpringCorrection

extern double MaxProd;                   // Maximum primary productivity (kgC/m2/day)
extern double MinProd;                   // Minimum primary productivity
extern double MaxRootDepth;              // Maximum root depth (m)
extern int NoRootsBelowGWT;              // if 1, no roots will grow below groundwater table (No telmatophytes)
extern double RootLambda;                // decay rate exponential root distribution function
extern double RootSenescence;            // root senescence factor = proportion of root mass that dies during each time step
extern double InitRoots;                 // Initial root mass in all layers
//extern Matrix RootMass;                  // initial root distribution (kg C/m2 in each layer)
extern double ExudateFactor;             // mass fraction of of below-ground production that consists of exudates
extern double BioMass;                   // above ground biomass kg C /m2 (standing crop)
extern double BioMassSenescence;         // biomass senescence at each DAY as fraction of above-ground biomass
extern Matrix Harvest;                   // harvest dates (1st column) and fraction of biomass harvested (2nd column)
extern Matrix Grazing;                   // Parts of the year in which garazing occurs, each row is a range of days
                                         // followed by the amount of biomass removed (kg C m2/day and the amount of excretion (kg C m2/day)
extern Matrix Manure;                    // Manure application dates (column 1) and quantity (column 2) in kg C/m2/day
extern double ManureFluidFrac;           // fluid fraction of manure
extern Matrix ManureLayers;              // partitioning of manure among layers; first column: fluids; second column: solids
extern char NPPFile[];                // file with primary production data for each time step; these data override any selected production model
extern Matrix NPPData;                   // net primary production data from file NPPfile
extern Matrix PARData;                   // photosynthetic active radiation or cloud cover data
extern Matrix Precipitation;			 // Precipitation data (read from file)
extern Matrix Evaporation;				 // Evaporation data (read from file)

extern Matrix MethaneReservoirs;         // This parameter tells which reservoirs are summed for calculating methane production from easily decomposed
extern double MethaneR0;                 // Methane production rate factor for fresh organic C microM/h

extern double MethanepHCorr;             // for every PH unit lower or higher than neutral, MethanepHCorr*R0 is added to R0
extern double MethaneQ10;                // Q10 value for temperature correction methane production; range 1.7 - 16 ref. in Walther & Heimann 2000
extern double MethaneOxQ10;              // Q10 value for temperature correction methane oxidation; range 1.4 - 2.1, ref. in Walther & Heimann 2000
extern double MethaneVmax;               // Vmax Michaelis-Menten eq methane oxidation micrMol/hr range 5-50
extern double MethaneKm;                 // Km Michaelis-Menten eq methane oxidation micrMol range 1-5
extern double MethaneMaxConc;            // maximum methane concentration in pore water
extern double MethaneERateC;             // Ebullition rate constant (1/hr)
extern double MethanePRateC;             // Rate constant for plant transport of methane (1/hr)
extern double MethaneAir;                // Methane concentration in the atmosphere - microMol
extern double MethaneTRef;               // Reference temperature for temperature sensitivity methane production
extern double MethanePType;              // Vegetation type factor for gas transport by plants range: 0-15
extern double MethanePlantOx;            // Fraction of methane that is oxidized during transport in plants

extern double CO2CH4ratio;              // Molar ratio between CH4 and CO2 production; for acetate splitting this is 1, for CO2 reduction 0

extern double PartialAnaerobe;           // Determines the slope of the relation of partial anaerobe soil fraction above the water table to soil saturation, >1
extern double AnaerobeLagFactor;	     // Determines time lag for development of sufficiently anaerobic conditions after saturation of a layer
extern Matrix InitMethane;               // Initial methane concentration profile
extern double DayMinGW;                  // day of minimum groundwater table
extern double MinGW;                     // lowest water table level at this day (m below surface)
extern double AmplitudeGW;               // amplitude of water table movement
extern char GwFile[];                    // file where the groundwater table time series is stored; if empty a sinusoidal time series is assumed
extern char PrecipFile[];             // file where the precipitation time series is stored; if empty the water table is assumed to be read from file
extern char EvapFile[];               // file where the evaporation time series is stored; if empty the water table is assumed to be read from file
extern char SoilMoisture[];              // file where soil moisture profile time series is stored; if empty the soil moisture will be calculated using very simplified assumptions
extern Matrix MoistProfiles;             // Matrix with soil moisture profiles from observational data or other model

extern double T_average;                 // average yearly temperature
extern double T_amplitude;               // amplitude of temperature throughout the year
extern double T_ref;                     // reference temperature for correction of decay constants
extern double ThermDiff;                 // thermal diffusivity, if not defined it is estimated from soil properties
extern char TFile[];                     // file with temperature time series
extern char SnowFile[];                  // file with snowdepth time series
extern Matrix SnowData;                  // Snow depth data from file Snowfile
extern Matrix TData;                     // Air temperature data from file Tfile
extern Matrix GwData;                    // Groundwater table data from GWFile
extern Matrix T_init;                    // Initial temperature profile
extern char SoilProfile[];               // Name of the file where the soilprofile data are stored
extern int NrHorizons;                   // number of horizons
extern Matrix Horizons;                  // Horizon base depths with respect to surface
extern Matrix CNRatio;                   // CN ratios for eac soil layer; the decomposition of peat can be made dependent on these
extern Matrix DBD;                       // Dry bulk density for each horizon
extern Matrix PercOrg;                   // Percentage organic matter for each horizon
extern Matrix Layer_pH;                  // pH
extern Matrix InitRes;                   // initial contents of SOM reservoirs in % of total SOM including peat
extern Matrix Layer_pF;                  // pF curves, 2 possible definitions:
                                         // as phi values for potentials defined in pFVal
                                         // or Van genuchten parameters theta_r, theta_s, alfa, l and n
                                         // In the latter case the rows in Layer_pF have only 5 elements and pFVal only 1 column
extern Matrix pFVal;                     // suction potentialsfor pF curves

extern Matrix Porosity;                  // porosity of soil layers, can be specified or calculated from organic matter percentage and dry bulk density
extern Matrix ProfileOutput;             // determines which vertical profiles are sent to log files
extern int ProductionModel;              // Production model: 0 for simple sinusoidal function; 1 for production dependent on temperature of upper soil layer, 2 for data supplied externally, 3 - 4 for photosynthesis models
extern double GreenBiomassRatio;        // Ratio of photosyntesizing biomass to total bioamass for ProductionModel 3
extern Matrix HarvestCorrection;            // Correction of GPP after harvest with reduction factor directly after harvest (first) and period of recovery in days (second)
extern Matrix TotalMethane;              // storage matrix for CH4 fluxes ; 1st element: day number
extern Matrix ReservoirTime;             // storage matrix for CO2 per reservoir per timestep ; 1st element: day number
extern Matrix LayerTime;                 // storage matrix for CO2 per layer per timestep ; 1st element: day number
extern Matrix BioMassRec;                // storage of biomass, primary production and plant respiration
extern Matrix FreezingCurve;             // Unfrozen water content curve at below zero temperatures
extern Matrix LatentHeat;		 // parameters for approximation of temperature-dependent latent heat of fusion of ice J kg-1
extern double MaxSnowdepth;              // Maximum snow depth
extern double DayMaxSnowdepth;           // Day of maximum snow depth
extern double SnowMeltrate;              // Rate of snowmelt per degree C above zero per day
extern double WatertableInit;		  	// Initial water table in m, has to be specified if the watertable is calculated by the model
extern double EvapCorrection;	 		// Correction factor to reduce evaporation if water table is below surface, for water table model
extern double RunoffThreshold;			// Threshold above which a ponded water layer produces runoff; for water table model
extern double OpenWaterFactor;			// evaporation correction factor for open water evaporation
extern double CropFactor;			    // Makkink Crop factor to correct evaporation for vegetation properties; for water table model
extern Matrix SandFraction;			    // sand weight fraction of mineral fraction
extern Matrix ClayFraction;					// clay weight fraction of mineral fraction, influences decomposition rate humus reservoir
extern double VegTScalingFactor;			// scaling factor for air to soil surface temperature; put to one if soil surface temperature is input; outherwise a valye between 0.6 and 1.0
extern char PARFile[];                  // file with PAR (W m-2) or cloud cover data for NPP model
extern double Latitude;                    // site latitude, to be used for photosynthesis model if only cloud cover data is supplied
extern Matrix Phenology;                // Phenology data for production model 3 and 4
extern double AmbientCO2;               // Ambient CO2 concentration
extern char CO2File[];                 // file with yearly averaged CO2 concentration for multi-year runs
extern Matrix CO2Data;                  // CO2 data for longer model runs with photosynthesis model
extern double LitterLayer;              // organic matter stored in above ground litter layer, in kg C / m2
extern double LitterConversion;         // Conversion factor of daily conversion of above ground to below ground litter at reference temperature Tref; the factor is temperature adjusted such that at 0 degrees the conversion factor is also 0
extern double DrainageDist;             // Distance to nearest drainage channel (m)
extern double Ksat;                     // Saturated hydraulic conductivity of soil
extern double DrainLevel;               // reference level of water in the drains/river channel with respect to top of soil surface (m)
extern Matrix DrainData;                // Matrix to store variable drain water level data
extern char DrainageFile[];             // file name with drain water levels
extern char RunOnFile[];                // file name for runon data
extern Matrix RunOn;                    // Matrix with Run-on water quantity
extern double KBeer;                    // Beer's law constant for photosynthesis models, values around 0.5
extern Matrix PhotoPar;                 // Parameters for photosynthesis model 5 for tundra, Shaver et al, J. Ecology 2007
extern double LAICarbonFraction;         // relates leaf area index to kg C/m2 
extern int AnaerobicCO2;                       // Switch for allowing anaerobic decomposition (sulfate etc) resulting in CO2
extern Matrix KAnaerobic;                      // Anaerobic decomposition constants
extern double Q10Anaerobic;                    // Q10 anaerobic decomposition
extern double KLitter;                   // decomposition constant above-ground litter and standing dead biomass
extern int AnaerobicCO2;                   // Switch for allowing anaerobic decomposition (sulfate etc) resulting in CO2, if 0 not accounted for
extern Matrix LayerAnaerobic;             // Anaerobic CO2 per layer
extern double HarvestLitter; // fraction of harvest that is left as litter at each harvest
extern Matrix CarbonBalance;             // Carbon balance: primary production, C exported, and change in carbon reservoirs in Mol C
extern Matrix PeatDecay;                        // logs true loss of peat matrix

void cmdline(int argc, char *argv[]);
/* handles the command line options
argc and argv are the variables passed to main */


int readstring(char x[], const char *ident, const char *filename, const int verbose);
/* reads a string parameter
    x       : variable to be read
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored
    verbose : prints warning messages if TRUE                       */

int readscalar(double *x, const char *ident, const char *filename, const int verbose);
/* reads a single value parameter
    x       : variable to be read
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored
    verbose : prints warning messages if TRUE                       */                        // reads a single value parameter


int readarray(double *x, int *len, const char *ident, const char *filename, const int verbose);           // reads an array
/*  x       : variable to be read
    r       : number of rows
    c       : number of columns
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored                       */                 // reads an array


int readmatrix(double *x, int *r, int *c, const char *ident, const char *filename, const int verbose);           // reads a matrix
/*  x       : variable to be read
    r       : number of rows
    c       : number of columns
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored                       */

int readall();                     // read all parameters
/* THIS FUNCTION HAS TO BE ADAPTED IF NEW PARAMETERS ARE ADDED TO THE MODEL */

int readsoil(char *soilname);         // reads soil profile from the file specified in soilname

int readTseries(double *d, const char *filename, const int maxrows, const int cols);
/*  Reads a time series from file, e.g. groundwater level or temperature data
    This may be a multivariate time series, with the variable stored in cols.
    Each line in the file represents observations on one time point

    Tseries : Matrix that stores the data
    filename: the name of the file where the series is stored
    maxrows : the maximum number of data rows to be read
    cols    : the number of data columns                          */

int readHarvest(char *harvestIN);

void checkHarvestDate();

void WriteOutput();
/* writes methane fluxes, CO2 fluxes and plant production/respiration to output files */


bool Assign_Start(string StartDate);
// Assigns values to StartDay, StartMonth, StartYear.

bool Assign_end(string EndDate);
// Assigns values to EndDay, EndMonth, EndYear.

bool isALeapYear(int year);
// Checks if year is a leap year

void isFebruary(bool LeapYear, int daystring);
// If month is February in Start or End year, number of days in month is checked.

int assignYEARdays(int year);
// Assigns YEARdays depending on whether it is a leap year.

// Counts the number of days between two dates.
int count_days(int day_1 = StartDay, int month_1 = StartMonth, int year_1 = StartYear, int day_2 = EndDay, int month_2 = EndMonth, int year_2 = EndYear);

int daysinmonth(int month, int year);
// Checks how many days should be in a month

int count_years(int year_1 = StartYear, int year_2 = EndYear);
// Counts the years of the model simulation.


