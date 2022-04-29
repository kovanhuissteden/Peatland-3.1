/***************************************************************************
                          paramcheck.h  -  parameter checking functions
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


#define PARAM_ERROR1 " does not match number of horizons."
#define PARAM_ERROR2 "Check pF curves: water content at saturation higher than pore volume in horizon "
#define PARAM_ERROR3 "Initial reservoir content in InitRes does not sum to 1 in horizon "
#define PARAM_ERROR4 "Model layers exceed depth of soil profile."
#define PARAM_ERROR5 "Initial methane profile and number of layers do not match."
#define PARAM_ERROR6 "Matrices Manure and ManureLayers should have 2 columns and ManureLayers no more layers than model."
#define PARAM_ERROR7 "Columns of Matrix ManureLayers should sum to one. "

extern BOOLEAN Verbose;                  // output flag on-screen output
extern double DensOrg;                   // density organic matter in peat (density of wood) kgm-3
extern double DensMin;                   // density of mineral matter kg m-3
extern Matrix PercOrg;                   // Percentage organic matter for each horizon
extern Matrix Porosity;                  // porosity of soil layers, can be specified or calculated from organic matter percentage and dry bulk density
extern Matrix DBD;                       // Dry bulk density for each horizon
extern int NrHorizons;                   // number of horizons
extern int NrLayers;                     // Number of depth steps
extern double LayerThickness;            // Thickness of depth step
extern Matrix Layer_pF;                  // pF curves
extern Matrix pFCurves;                  // final pF Curves
extern Matrix pFVal;                     // suction potentialsfor pF curves
extern Matrix InitRes;                   // initial contents of SOM reservoirs in % of total SOM including peat
extern int NrReservoirs;                 // Number of SOM reservoirs
extern Matrix Horizons;                  // Horizon base depths with respect to surface
extern Matrix InitMethane;               // Initial methane concentration profile
extern Matrix Manure;                    // Manure application dates (column 1) and quantity (column 2) in kg C/m2/day
extern Matrix ManureLayers;              // partitioning of manure among layers; first column: fluids; second column: solids
extern BOOLEAN ApplyVanGenuchten;		 // indicates application of Van Genuchten equation for soil moisture
extern double StepsizepF;				 // pF step size for Van Genuchten function

/************************ FUNCTION HEADERS ************************************/

void Porevol();
/* Calculates porosity from dry bulk density and percentage organic  matter
if the porosity is not defined in the soil profile file */

int Paramchk();
/* performs a number a parameter checks:
match between density, porosity and pF curves
contents of reservoir initialization InitRes
depth of soil profile vs. model layers
length of initial methane profile
manure matrices                                                 */



