/***************************************************************************
                general.cpp  -  general computing functions for PEATLAND
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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <climits>
#include <cmath>
using namespace std;
#include "matrix.h"
#include "general.h"

void interp(Matrix &x, Matrix &y, Matrix &xi, Matrix &yi, int decr)
/* interpolates values in yi from data given in x and y and interpolation points in xi
   decr is TRUE if x and xi are decreasing, FALSE otherwise
   x, y, xi and yi may have larger dimensions than row-only or column-only vectors
   However, x and y should have the same number of entries, xi and yi also
   Furthermore x and xi should be monotonically increasing or decreasing
   If not unpredictable results and programm errors will arise
   NO CHECKS ARE BEING DONE ON INCREASE OR DECREASE                                   */
{
  int i, indx, nx, ny, nxi, nyi;
  double *xd, *yd, *xid, *yid, x1, y1, y2;

  nx = x.Length();                   // check the size of x and y; must contain the same number of data
  ny = y.Length();
  nxi = xi.Length();
  nyi = yi.Length();
  if ((nx != ny) || (nxi != nyi))
  {
    cout << GEN_ERROR1 << endl;
    exit (EXIT_FAILURE);
  }
  xid = xi.Data();
  xd = x.Data();
  yd = y.Data();
  yid = yi.Data();
  indx = 0;
  for (i = 0; i < nxi; i++)                      // loop along values of xi;
  {
    if (decr)                                    // find values in x between which we interpolate
    {
      while (*xid <= *xd)
      {
        indx++;
        xd++;
      }
    } else
    {
      while (*xid >= *xd)
      {
        indx++;
        xd++;
      }
    }
    if ((indx == 0) || (indx >= nx))              // fatal error: xi out of range of x
    {
      cout << GEN_ERROR3 << endl;
      exit (EXIT_FAILURE);
    }
    y1 = *(yd + indx - 1);                        // simple linear interpolation
    y2 = *(yd + indx);
    x1 = *(xd - 1);
    *yid = y1 + (y2 - y1) * ((*xid - x1) / (*xd - x1));
    if (i < (nxi - 1))                             // prodeed to next element of xi and yi
    {
      xid++;
      yid++;
    }
  }
}
