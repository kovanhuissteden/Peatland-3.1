/***************************************************************************
           bookkeep.h  -  general bookkeping functions for PEATLAND
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

extern int StepNr;                             // time step number during iteration
extern double DayNr;                           // midpoint of simulated timestep
extern double Timer;                           // starting point of simulated time step
extern double Timestep;                        // model timestep
extern int Verbose;                            // output flag for run information to console
extern double DayOfTheYear;                    // Julian day number of the midpoint of the simulated time step relative to the current year;
extern double Year;                            // current simulation year
extern int Month;
extern int YearCounter;                       // current simulation year
extern int CalendarYear;
extern int monthdays;                       // Number of days in the current month.

extern int Harvest_LOC;
extern int HYY;
extern int HMM;
extern int HYY;
extern int HdayNr;
extern double Harv_height;
extern int Harvest_LOC_END;


void TrackTime();                              // Time system, updated at end of each model time step

void checkHarvestDate();