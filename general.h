
/***************************************************************************
         general.h  -  some general functions and declarations for PEATLAND
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
 
 May 2012
 
 Corrected:  year lenght for day of the year calculation
 
 ***************************************************************************/

//#include "matrix.h"

#define FALSE 0                        // definitions for boolean operations
#define TRUE 1

typedef int BOOLEAN;


#define DEFAULTS "defaults"             // file with default parameters
#define MAXHORIZON 100                  // max number of horizons
#define MAXLAYER   200                  // max number of layers
#ifndef PI
#define PI 3.141592653589793            // pi 3.141592653589793
#endif
#define DEG2RAD 0.017453292519943       // conversion factor degrees to radians
#define YEAR 365.25                     // length of one year in days
#define C_CO2 3.6641412039                   // C- CO2 mass conversion factor
#define MOLWEIGHTCO2 44.010             // molecular weight CO2
#define MOLWEIGHTC 12.011               // molecular weight C
#define MOLWEIGHTO 15.994               // molecular weight O
#define MOLWEIGHTCH4 16.048             // molecular weight CH4
#define SOLARCONSTANT 1361.5            // solar constant in W/m2
#define CONVKGCTOMOLC 83.2570144035    // conversion factor of kg C to moles C (multiply by this)
#define CONVKGCO2TOMOLC 22.7221086117  // conversion factor of kg C to moles C (multiply by this)


#define GEN_ERROR1 "Function interp: x and y matrix contain an unequal amount of numbers"
#define GEN_ERROR2 "Function interp: x matrix not monotonically rising"
#define GEN_ERROR3 "Function interp: entries in xi out of range of x"


void interp(Matrix &x, Matrix &y, Matrix &xi, Matrix &yi, int decr);
/* interpolates values in yi from data given in x and y and interpolation points in xi
   decr is TRUE if x and xi are decreasing, FALSE otherwise
   x, y, xi and yi may have larger dimensions than row-only or column-only vectors
   However, x and y should have the same number of entries, xi and yi also
   Furthermore x and xi should be monotonically increasing or decreasing
   If not unpredictable results and programm errors will arise
   NO CHECKS ARE BEING DONE ON INCREASE OR DECREASE                                   */


