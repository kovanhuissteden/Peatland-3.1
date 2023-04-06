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
#include "readparams.h"

void TrackTime()            // Time system, updated at end of each model time step
{
    double t, d;
    int yearlength = 365;
    int FinalMonthDay = 0;

    YEARdays = assignYEARdays(CalendarYear); 
    monthdays = daysinmonth(Month, CalendarYear);
    FinalMonthDay = count_days(1, 1, CalendarYear, monthdays, Month, CalendarYear);
    // Number of days in the year. Checks for leap year.

    if (DayOfTheYear > FinalMonthDay)  Month += 1;

    StepNr++;                 // time step number
    DayNr += Timestep;        // day number since day 1 of the year in which the simulation started
    Timer += Timestep;        // day number since start of simulation
    DayOfTheYear += Timestep;

    
    if (DayOfTheYear > YEARdays)
    {
        Year += 1;
        DayOfTheYear = DayOfTheYear - YEARdays;
        CalendarYear += 1;
        YearCounter +=1;file:///home/ko/bedrijf/administratie/Uitgaande_facturen/factuur_Bloemendaal_2.pdf
        Month = 1;
    }
    if (Verbose) cout << "Timestep: " << StepNr << " finished " << endl;
}


void checkHarvestDate()
{

    int HDay = 1;

    if ( Harvest_LOC < Harvest_LOC_END)
    {
        Harvest_LOC += 1; 
        //go to next harvest date
    }
    else Harvest_LOC = Harvest_LOC_END; // no more harvest dates in file

    HYY = HarvestData(Harvest_LOC, 1);
    HMM = HarvestData(Harvest_LOC, 2);
    HDay = HarvestData(Harvest_LOC, 3);
    Harv_height = HarvestData(Harvest_LOC, 4);

    HdayNr = count_days(1, 1, HYY, HDay, HMM, HYY); // used to compare to model simulation day nr.
    HDD = HDay;
}


