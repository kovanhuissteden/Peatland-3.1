/***************************************************************************
                          paramcheck.cpp  -   parameter checking functions
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

  April 2004

  Added: display of calculated porosity values to facilitate correction of
  erros in bulk density, pF curves or organic matter content
  
  April 2011
  pF curve: l is put at the end of the LayerPf matrix for calaculation
  of Van Genuchten (l is not yet being used, meant for calculation of
  unsaturated hydraulic conductivity
  function VanGenuchten adapted for general use and moved to water.cpp
  pore volume chaeck replaced by density check; pore volume
  is now taken from the van Genuchten parameters or pF curve (theta_sat)
  
****************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <climits>
#include <cmath>
using namespace std;
#include "matrix.h"
#include "general.h"
#include "paramcheck.h"
#include "water.h"

void Porevol()
/* Calculates porosity from dry bulk density and percentage organic  matter
if the porosity is not defined in the soil profile file */


/* REPLACE WITH RECALCULATION OF density; POROSITY IS KNOWN FROM
PF CURVES (THETA SAT IN CASE OF VAN GENUCHTEN, FIRST VALUE OF Layer_pF otherwise */
{
  int i;
  double p;

  if (Verbose)
  {
	cout << "Density input data:" << endl;
	DBD.Disp();
  }
  if (Porosity(1) == 0.0)           // it is assumed that the porosity parameter not has been set when the first element is 0
  {
    Porosity.Resize(NrHorizons);
    for (i = 1; i <= NrHorizons; i++)
    {
      if (ApplyVanGenuchten) p = Layer_pF(i, 2); else p = Layer_pF(i, 1);  // take porosity from pF data
	  Porosity(i) = p;
	  DBD(i) = (1 - p) * (DensOrg * PercOrg(i) / 100 + DensMin * (100 - PercOrg(i)) / 100 ); // recalc of density
    }
    if (Verbose)
    {
     
	  cout << "Calculated densities" << endl;  // display for checking the input data
      DBD.Disp();
	  cout << "Porosities:" << endl;
      Porosity.Disp();
    }
  }
}

int Paramchk()
/* performs a number a parameter checks:
contents of reservoir initialization InitRes
depth of soil profile vs. model layers
length of initial methane profile
manure matrices                                                 */
{
  int ok = TRUE, i, j, l, m, n, c, a, nsteps;
  double total, pF, theta;
                                                    // check array size pF matrix
  m = Layer_pF.Rows();
  n = Layer_pF.Cols();
  c = pFVal.Cols();
  if (m < NrHorizons)
  {
    cout << "pF data" << PARAM_ERROR1 << endl;
    ok = FALSE;
  }
  
/* Directly applying the Van Genuchten curves makes the model slow, in particular when water table is calculated
The solution of the water table takes repeated integration of the Van Genuchten equation which is computationally 
intensive.
Therefore a lookup table is constructed to read the Van genuchten curve.
Because pF curves of peat soils tend to be steep with pF values > 1, the intervals are better given in cm rather than pF values */
	nsteps = (int)(100 * NrLayers * LayerThickness) + 2;		   // nr of 1 cm depth steps on profile
	pFCurves.Resize(NrHorizons, nsteps);			               // pFCurves: theta value lookup table
	if ((n == 5) && (n != c))                           // Convert Van Genuchten curve to theta values for more rapid interpolation during moisture calculations
	{
		ApplyVanGenuchten = TRUE;
		for (i = 1; i <= NrHorizons; i++)
		{
			for (j = 0; j <= (nsteps - 1); j++) pFCurves(i, j + 1) = VanGenuchten((double)j, Layer_pF(i,1),  Layer_pF(i,2),  Layer_pF(i,3),  Layer_pF(i,4), 0);
		}
	} else
	{
		ApplyVanGenuchten = FALSE; // enumerated pF curve
		l = pFVal.Cols();
		for (i = 1; i <= NrHorizons; i++)
		{
			for (j = 0; j <= (nsteps - 1); j++)
			{
				if (j < 2) theta = Layer_pF(i, 1); else  // interpolate theta values from enumerated pF curve
				{
					a = 1;
					pF = log10(j);
					while (pF > pFVal(1, a)) a++;
					if (a <= l)
					{
						theta = Layer_pF(i, a) + ((pF - pFVal(a))/ (pFVal(a - 1) - pFVal(a))) * (Layer_pF(i, a - 1) - Layer_pF(i, a));
					}
					else theta = Layer_pF(i,l);
				}
				pFCurves(i, j + 1) = theta;
			}
		}
	}
	// pFCurves.Disp();
  for (i = 1; i <= NrHorizons; i++)                  // check the contents of InitRes
  {
    total = 0;
    for (j = 1; j <= NrReservoirs; j++) total +=InitRes(i,j);
    if (abs(total - 1.0) > 0.00001)
    {
      ok = FALSE;
      cout << PARAM_ERROR3 << i << " Sum: " << total <<endl;
    }
  }
  if (NrLayers * LayerThickness > Horizons(NrHorizons)) // check depth of soil profile against model layers
  {
    ok = FALSE;
    cout << PARAM_ERROR4 << endl;
  }
  if (InitMethane.Length() != NrLayers)           // check length of initial methane profile
  {
    ok = FALSE;
    cout << PARAM_ERROR5 <<endl;
  }
  // check manure matrices
  if ((Manure.Cols() != 2) || (ManureLayers.Cols() != 2) || (ManureLayers.Rows() > NrLayers))
  {
    ok = false;
    cout << PARAM_ERROR6 <<endl;
  }
    if ((ManureLayers.SumCol(1) != 1) || (ManureLayers.SumCol(2) != 1))
  {
    ok = false;
    cout << PARAM_ERROR7 <<endl;
  }
  return ok;
}


