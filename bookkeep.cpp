/***************************************************************************
           bookkeep.cpp  -  general bookkeping functions for PEATLAND
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
 
 Corrected: day of the year calculation
 
 ***************************************************************************/
#include <cmath>
#include <iostream>
using namespace std;
#include "matrix.h"
#include "bookkeep.h"
#include "general.h"

void TrackTime()            // Time system, updated at end of each model time step
{
    double t, d;
    int yearlength = 365;

    if (Verbose) cout << "Timestep: " << StepNr << " finished " << endl;
    t = Year / 4.0;
    d = t - floor(t);
    if (d == 0.0) yearlength = 366;  // adjust year length to leap year
    t = Year / 1000.0;
    d = t - floor(t);
    if (d == 0.0) yearlength = 365;
    StepNr++;                 // time step number
    DayNr += Timestep;        // day number since day 1 of the year in which the simulation started
    Timer += Timestep;        // day number since start of simulation
    DayOfTheYear += Timestep;
    if (DayOfTheYear > yearlength)
    {
        DayOfTheYear = DayOfTheYear - yearlength;
        Year += 1;
    }
}
