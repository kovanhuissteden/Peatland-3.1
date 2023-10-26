/***************************************************************************
             water.h  -  water submodel functions of PEATLAND
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


#define WT_ERROR1 "Error in water table module"
#define WT_ERROR2 "Warning: water table above surface computed from water balance"
#define WT_ERROR3 "Warning: water table drops below minimum : "



/************************ GLOBAL VARIABLES ***********************************/

extern double Timestep;                   // model timestep
extern double Timer;                      // starting point of simulated time step
extern int NrLayers;                      // Number of depth steps /layers
extern double LayerThickness;             // Thickness of layers
extern Matrix Porosity;                   // porosity of soil layers
extern Matrix GwData;                     // Groundwater table data from GWFile
extern Matrix Layers;                     // Layer boundaries
extern Matrix Layer_pF;                   // pF curves, 2 possible definitions:
extern Matrix pFVal;                      // suction potentials for pF curves
extern Matrix MoistTheta;                 // soil moisture profile for timestep
extern Matrix MoistProfiles;              // Matrix with soil moisture profiles from observational data or other model
extern Matrix Saturation;                 // pore volume saturation with water
extern Matrix OldSat;                     // saturation values at previous time step
extern Matrix LastSatTime;                // Time after last complete saturation of layer in days
extern Matrix PoreVol;                    // Porosity each layer
extern int StepNr;                        // time step number during iteration
extern ofstream *output2a;                 // output stream soil moisture profile
extern ofstream *output2b;
extern ofstream *output7;                 // output stream water table
extern Matrix ProfileOutput;              // determines which vertical profiles are sent to log files
extern Matrix MatricPotential;            // matric potential of model soil layer
extern double CurrentGW;                  // Current water table during time step
extern int TopSat;                        // Index of first saturated layer
extern Matrix pFCurves;                   // final pF Curves after conversion from van genuchten parameters
extern int WatertableModel;				  // choice for water table model:0 = read from file, 1 = simple sinusoidal 2 = extended 'Yurova' type model
extern Matrix Precipitation;			 // Precipitation data (read from file)
extern Matrix Evaporation;				 // Evaporation data (read from file)
extern double WatertableInit;		  	// Initial water table in m, has to be specified if the watertable is calculated by the model
extern double EvapCorrection;	 		// Correction factor to reduce evaporation if water table is below surface, for water table model
extern double RunoffThreshold;			// Threshold above which a ponded water layer produces runoff; for water table model
extern double OpenWaterFactor;			// evaporation correction factor for open water evaporation
extern double CropFactor;			    // Makkink Crop factor to correct evaporation for vegetation properties; for water table model
extern double BelowWTStorage;			// Below water table water storge for water table model
extern double AboveWTStorage;			// Above water table water storge for water table model
extern double MaxBelowWTStorage;	    // Maximum Below water table water storge for water table model
extern double MinAboveWTStorage;		// Minimum Above water table water storge for water table model
extern Matrix SoilTemp;                 // soil temperatures model layers, interpolated from TProfile
extern Matrix TData;                    // Air or soil temperature data from file Tfile
extern double MinGW;                    // lowest water table level (m below surface)
extern int NrHorizons;                         // number of horizons
extern Matrix Horizons;                        // Horizon base depths with respect to surface
extern double PondedWater;				// layer of ponded water above the soil
extern BOOLEAN ApplyVanGenuchten;		 // indicates application of Van Genuchten equation for soil moisture
extern double FrozenStorage;			// Frozen water storage for water table model
extern double CurrentFrozen;		    // Current top of frozen layer
extern double MaxDepthHeat;              // maximum depth temperature model (m)
extern double CurrentGWBase;			 // current lowest water table level for water table model
extern double FrozenWatertable;			 // Water table in frozen part of soil
extern double SnowStorage;				 // water storage in snow (m)
extern double DrainageDist;              // Distance to nearest drainage channel (m)
extern double Ksat;                     // Saturated hydraulic conductivity of soil
extern double DrainLevel;                // reference level of water in the drains/river channel with respect to top of soil surface (m)
extern double DrainageDepth;            // Thickness of unit through which drainage occurs
extern Matrix DrainData;                       // Matrix to store variable drain water level data
extern char DrainageFile[];             // file name with drain water levels
extern char RunOnFile[];                 // file name for runon data
extern Matrix RunOn;                     // Matrix with Run-on water quantity

/************************ FUNCTION HEADERS ************************************/
double  VanGenuchten(double pF, double theta_r, double theta_s, double alfa, double n, int convert);
/* Expands van Genuchten parameters to complete pF curve
pF is the suction values for which the moisture content theta is computed
theta_r, theta_s, alfa, n, l are the parameters of the Van Genuchten function
Note: these are stored in the Layer_pF matrix  read from the parameter file
theta contain is the output  
convert indicates whether pF is given as true pF (convert > 0) or cm (convert <= 0)                                                      

Note: l (last value in rows of Layer_pF) is not used for the water retention curve,
only for conductivity acc to Van Genuchten and will be implemented later */


double FindWaterBase();
/* returns the base of the water level movement in the profile,
which is either the topmost frozen layer or the lowest possibele water table */

double Drainage(double watertable);
/* calculates drainage or seepage from distance to drainage line and hydraulic conductivity */

int belowWT(double currentbase, double watertable);
/* calculates the amount of water below the water table to the base of water table fluctuation range
currentbase is the cuurent base of the soil water profile, either permafrost top or lowest water table
the function returns the index of the topmost partly saturated soil horizon, this is set to 0 if there is ponded water
or the index of the topmost frozen layer if it is used for calculating the frozen water storage*/

double CalcFrozenStorage(double watertable);
/* Calculates the amount of water that is stored in frozen layers above the minimum water table */

double aboveWT(int tophorizon, double watertable);
/* calculates the above water table unsaturated storage
tophorizon is the index to the horizon which contains the water table*/

void Watertable();
/* calculates water table position from precipitation and evaporation*/


void Moisture(int initial);
/* calculates soil moisture profile assuming that soilmoisture is at gravitational equillibrium with the
groundwater table
if initial is TRUE, the soil moisture profile is initialized                                  */
