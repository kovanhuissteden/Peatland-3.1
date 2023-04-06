/***************************************************************************
 *   Copyright (C) 2005 by Ko van Huissteden   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/***************************************************************************
 main.cpp  -  Peatland 3.0 main program
 -------------------
 begin                : |02-01-2003|
 copyright            : (C) |2003| by |J. van Huissteden|
 email                : |j.van.huissteden@vu.nl|
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <climits>
using namespace std;
#include "matrix.h"
#include "general.h"
#include "peatland.h"
#include "readparams.h"
#include "paramcheck.h"
#include "initialize.h"
#include "water.h"
#include "heat.h"
#include "bookkeep.h"
#include "SOMproduction.h"
#include "SOMdecomposition.h"
#include "methane.h"

int main(int argc, char *argv[])
{
    
    int i;
    
    /*************************** Data read section *************************************/
    cmdline(argc, argv);                                    // Process command line options
    if (!readall())                                         // Read model parameters
    {
        cout << "Errors reading parameter files" << endl;
        return EXIT_FAILURE;
    }
    if (!readsoil(SoilProfile))                             // Read soil profile data
    {
        cout << "Errors reading soil profile or time series data" << endl;
        return EXIT_FAILURE;
    }
    /***************************Initialization and parameter checking*******************/
    if (!Paramchk()) return EXIT_FAILURE;                   // Check porosity, pF curves, initial reservoir C
    Porevol();                                              // Check porosity and density if necessary
    MakeLayers();                                           // set up model layers and pointers to soil horizons
    InitTime();                                             // set up time system
    InitDecomp();                                           // convert raw decomposition constants to net decomposition excluding microbialbiomass formation
    RootsInit();                                            // initializes root distribution
    SOMResInit();                                           // initializes SOM reservoirs
    Moisture(TRUE);                                         // initialize soil moisture profile
    InitHeat();                                             // initializes thermal model parameters
    InitMethaneModel();                                     // initializations for methane model
    InitTseries();                                          // initialize temperature and groundwater table time series if necessary
    InitWater();
    OutputInit();                                           // initializes output matrices
    InitLogFiles();                                         // initialize output log files for intermediary variables
    for (i = 1; i <= NrOfSteps; i++)                        // model iteration
    {
        if (WatertableModel == 2) Watertable();				  // calculate watertable
        Temperature();                                        // soil temperature
        Moisture(FALSE);                                      // soil moisture profile

        OldSOM = NewSOM;
        OrgProd();                                            // Net primary production and partitioning among roots and shoots
        //cout << NewSOM(1,4) - OldSOM(1,4) << " " << NewSOM(2,4) - OldSOM(2,4) << " "  << NewSOM(3,4) - OldSOM(3,4) << endl;
        EnviCor();                                            // environmental correction factors for first order decomposition constants aerobic decomposition of the SOM reservoirs
        Decompose();                                          // aerobic and anaerobic decomposion of SOM above the water table (anaerobic excl. methane)
        //cout << NewSOM(1,4) - OldSOM(1,4) << " " << NewSOM(2,4) - OldSOM(2,4) << " "  << NewSOM(3,4) - OldSOM(3,4) << endl;
        Methane();                                            // methane fluxes modified after the model of Walther (2000), Global Biogeochemical Cycles 14, 745 - 765
        //cout << NewSOM(1,4) - OldSOM(1,4) << " " << NewSOM(2,4) - OldSOM(2,4) << " "  << NewSOM(3,4) - OldSOM(3,4) << endl;

        WriteSOMReservoirs();                                 // write SOM reservoirs to output files
        CollectCO2();                                         // Collect all CO2 and store in output arrays
        TrackTime();                                          // update time system
    }
    CloseLogFiles();
    CollectBioMass();

/* To be added:
 * CollectCarbonBalance()
 * Carbon balance - new matrix CarbonBalance
 * 1: DayNr
 * 2: carbon released per reservoir incl litter
 * 3: CH4 carbon emitted
 * 4: CH4 storage change
 * 5: storage in biomass change
 * 6: anaerobically produced carbon besides CH4
 * 7: photosynthesis
 */

    WriteOutput();
    cout << "PEATLAND 2.1 model run finished succesfully" << endl;
    return EXIT_SUCCESS;
}
