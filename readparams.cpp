/***************************************************************************
                  readparams.cpp  -  parameter read functions
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

  May 2003

  Adaptions for reading of parameters correct modelling of soil freezing

  17 December 2003

  Bug fix in processing of the "?" command line parameter

  September 15, 2004

  Added: messages indicating the reading of input time series from
  external files

  Added: Provisions for using net primary production data from
  external model by reading data from file

  31 May 2006
  Fixed LatentHeat value replaced by polynomial approximation of temperature
  dependent LatentHeat

  29 November 2008
  Bug fix of read errors on time series 
  Added: dispaly of the file from which a parameter is read with verbose option
 
  May 2010 addition of photosynthesis model
  
  March 2021 parameter added for anaerobic CO2 and above ground litter/standing dead biomass decomposition

 ***************************************************************************/

#include <sys/stat.h> 
#include <sys/types.h> 
#include <string>
#include <vector>
#include <sstream>
using namespace std;


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <climits>
#include <cmath>
#include <cfloat>
using namespace std;
#include "matrix.h"
#include "general.h"
#include "readparams.h"

/*
std::string split implementation using delimeter as a character.
*/
vector<string> split(string strToSplit, char delimeter)
{
    stringstream ss(strToSplit);
    string item;
    vector<string> splittedStrings;

    while (getline(ss, item, delimeter))
        {
        splittedStrings.push_back( item );
        }

  return splittedStrings;
}


void cmdline(int argc, char *argv[])
/* handles the command line options
argc and argv are the variables passed to main */
{
  int i,l;
  BOOLEAN recognized;

  for (i = 1; i < argc; i++)                      // handle command line options
  {
    recognized = FALSE;
    l = strlen(argv[i]);
    if (l<2)
    {
      if (*argv[i] == 'v')
      {
        Verbose = TRUE;                          // the v option indicates text output to the console during execution
        recognized = TRUE;
      }
      if (*argv[i] == '?')                       // display help on command line options
      {
        cout << MSG1 << endl;
        cout << MSG2 << endl;
        cout << MSG3 << endl;
        cout << MSG4 << endl;
        cout << MSG5 << endl;
        exit(EXIT_SUCCESS);
      }
    } else
    {
      if (strncmp(argv[i], "dir=", 4) == 0)       // working directory
      {

         strcpy((DataDir), argv[i] + 4);
         l = strlen(DataDir);
#ifdef MSDOS
         DataDir[l] = '\\';
#else
         DataDir[l] = '/';
#endif
         DataDir[l + 1] = '\0';
         recognized = TRUE;
      }
      if (strncmp(argv[i], "par=", 4) == 0)          // specification of parameterfile
      {
        strcpy((ParamFile), argv[i] + 4);
        recognized = TRUE;
      }
      if (strncmp(argv[i], "out=", 4) == 0)          // specification of output file prefix
      {
        strcpy((OutputFilePrefix), argv[i] + 4);
        recognized = TRUE;
      }
    }
    if (!recognized) cout << argv[i] << ERRORMSG10 << endl;
  }
}


int readstring(char x[], const char *ident, const char *filename, const int verbose)  // reads a string parameter
/*  x       : variable to be read
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored
    verbose : prints warning messages if TRUE                       */
{
  char buf[1024];                                                    // input line buffer
  BOOLEAN cont = TRUE;
  BOOLEAN found = FALSE;
  int l;
  int t = '\"';
  char *p, *q;

  strcpy(buf, DataDir);
  ifstream file(strcat(buf,filename));                                          // open file
  if (!file)
  {
    cout <<  ERRORMSG1 << " " << filename << endl;
    return EXIT_FAILURE;
  } else {
    l = strlen(ident);                                               // if succesful, start reading lines until a line starting with the identifier is found
    while (cont)
    {
      file.getline(buf, 1023, '\n');
      if (file.fail()) cont = FALSE;
      if (strncmp(buf, ident, l) == 0)
      {
        if (Verbose) cout << MSG6 << ident << MSG7 << filename << endl;  // display variable name
        cont = FALSE;                                                 // Found!! Now get the data
        found = TRUE;
        p = strchr(buf, t);                                           // find the '=' sign
        if (p != NULL)
        {
          p++;
          q = strpbrk(p, "\"\n");                                     // find the end of the string
          strncpy(x, p, q - p);                                      // copy the string
          x[q - p] = 0;                                              // add the NULL character
        }
      }
    }
    file.close();
  }

  //if (!found) cout << "PARAMETER MISSING! Variable " << ident << " not found in file " << filename << endl;
  if ((!found) & (strcmp(filename, DEFAULTS) == 1)) cout << "PARAMETER MISSING! Variable " << ident << " not found in file " << filename << endl;
  return found;
}

int readscalar(double *x, const char *ident, const char *filename, const int verbose)  // reads a single value parameter
/*  x       : variable to be read
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored
    verbose : prints warning messages if TRUE                       */
{
  char buf[1024];
  char nr[1024];                                                    // input line buffer
  BOOLEAN cont = TRUE;
  BOOLEAN found = FALSE;
  int l;
  int t = '=';
  char *p, *q;

  strcpy(buf, DataDir);
  ifstream file(strcat(buf,filename));                                          // open file
  if (!file)
  {
    cout <<  ERRORMSG1 << " " << filename  << endl;
    return EXIT_FAILURE;
  } else {
    l = strlen(ident);                                               // if succesful, start reading lines until a line starting with the identifier is found
    while (cont)
    {
      file.getline(buf, 1023, '\n');
      if (file.fail()) cont = FALSE;
      if (strncmp(buf, ident, l) == 0)
      {
        if (Verbose) cout << MSG6 << ident << MSG7 << filename << endl;  // display variable name
        cont = FALSE;                                                 // Found!! Now get the data
        found = TRUE;
        p = strchr(buf, t);                                           // find the '=' sign
        if (p != NULL)
        {
          p++;
          q = strpbrk(p, ";%\n");                                     // find the end of the numerical part
          strncpy(nr, p, q - p);                                      // copy the numerical part
          nr[q - p] = 0;                                              // add the NULL character
          *x = atof(nr);                                              // convert to double
        }
      }
    }
    file.close();
  }

  // if ((!found) & verbose) cout << "Variable " << ident << " not found in file " << filename << endl;
  if ((!found) & (strcmp(filename, DEFAULTS) == 1)) cout << "PARAMETER MISSING! Variable " << ident << " not found in file " << filename << endl;
  return found;
}

int readarray(double *x, int *len, const char *ident, const char *filename,  const int verbose)                    // reads an array
/*  x       : variable to be read
    l       : the length of the array; if > 0 this determines the maximum number of array ellements read from the file
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored                       */

{
  char buf[1024];
  char nr[1024];                                                    // input line buffer
  BOOLEAN cont = TRUE;
  BOOLEAN found = FALSE;
  int l, maxlen, count = 0;
  int t = '[';
  char *p, *q;

  if (*len > 0) maxlen = *len; else maxlen = INT_MAX;               // determine data limits
  strcpy(buf, DataDir);
  ifstream file(strcat(buf,filename));                              // open file
  if (!file)
  {
    cout <<  ERRORMSG1 << " " << filename  << endl;
    return EXIT_FAILURE;
  } else {
    l = strlen(ident);                                               // if succesful, start reading lines until a line starting with the identifier is found
    while (cont)
    {
      file.getline(buf, 1023, '\n');
      if (file.fail()) cont = FALSE;
      if (strncmp(buf, ident, l) == 0)
      {
        if (Verbose) cout << MSG6 << ident << MSG7 << filename << endl;  // display variable name
        cont = FALSE;                                                 // Found!! Now get the data
        found = TRUE;
        p = strchr(buf, t);                                           // find the '[' sign that starts the array
        if (p != NULL)
        {
          p++;
          q = p;
          while ((*q != ']') && (count < maxlen))                     // copy numbers until the end of array sign is found or the maximum number is reached
          {
            q = strpbrk(p, " ,];%\n");                                // find the end of the numerical part
            while (*q == ' ') q++;                                    // skip any spaces
            strncpy(nr, p, q - p);                                    // copy the numerical part
            nr[q - p] = 0;                                            // add the NULL character
            *x++ = atof(nr);
            count++;
            p = q;
            if ((*p == ',') || (*p == ';')) p++;                       // skip delimiter
            while (*p == ' ') p++;                                    // skip any spaces
          }
          if (*q != ']') cout << ERRORMSG2 << ident <<endl;            // error message if there is more or less data then there should be
          if ((*len > 0) && (count < maxlen))
          {
            cout << ERRORMSG3 << ident <<endl;
            found = FALSE;                                              // this error may be fatal!
          } else *len = count;                                          // pass array length
        }
      }
    }
    file.close();
  }
  //if ((!found) & verbose) cout << "Variable " << ident << " not found in file " << filename << endl;
  if ((!found) & (strcmp(filename, DEFAULTS) == 1)) cout << "PARAMETER MISSING! Variable " << ident << " not found in file " << filename << endl;
  return found;
}

int readmatrix(double *x, int *r, int *c, const char *ident, const char *filename,  const int verbose)           // reads a matrix
/*  x       : variable to be read
    r       : number of rows
    c       : number of columns
    ident   : identification label in the file for the variable
    filename: the name of the file where x is stored                       */

{
  char buf[1024];
  char nr[1024];                                                    // input line buffer
  BOOLEAN cont = TRUE;
  BOOLEAN found = FALSE;
  int l, cc = 0, rc = 0, previous_cc, maxrows, maxcols, errcount = 0;
  int t = '[';
  char *p, *q;

  if (*r > 0) maxrows = *r; else maxrows = INT_MAX;                 // determine data limits
  if (*c > 0) maxcols = *c; else maxcols = INT_MAX;
  strcpy(buf, DataDir);
  ifstream file(strcat(buf,filename));                              // open file
  if (!file)
  {
    cout <<  ERRORMSG1 << " " << filename  << endl;
    return EXIT_FAILURE;
  } else {
    l = strlen(ident);                                               // if succesful, start reading lines until a line starting with the identifier is found
    while (cont)
    {
      file.getline(buf, 1023, '\n');
      if (file.fail()) cont = FALSE;
      if (strncmp(buf, ident, l) == 0)
      {
        if (Verbose) cout << MSG6 << ident << MSG7 << filename << endl;  // display variable name
        cont = FALSE;                                                 // Found!! Now get the data
        found = TRUE;
        p = strchr(buf, t);                                           // find the '[' sign that starts the array
        if (p != NULL)
        {
          p++;
          while (*p == ' ') p++;
          q = p;
          while ((*q != ']') && (rc < maxrows))                        // row read loop
          {
            cc = 0;
            if (*q == ';') q++;
            while ((*q != ']') && (*q != ';') && (cc < maxcols))     // column read loop copy numbers until the end of array sign is found
            {                                                          // row ends with ;or ]
              q = strpbrk(p, " ,];%\n");                               // find the end of the numerical part
              while (*q == ' ') q++;                                   // skip any spaces
              strncpy(nr, p, q - p);                                   // copy the numerical part
              nr[q - p] = 0;                                           // add the NULL character
              if (strlen(nr) == 0) errcount++;
              *x++ = atof(nr);
              cc++;                                                    // column counter
              p = q;
              if ((*p == ',') || (*p == ';')) p++;                     // skip delimiter and spaces
              while (*p == ' ') p++;
            }
            previous_cc = cc;
            if ((rc > 0) && (cc != previous_cc))                        // error message for 'ragged' or incomplete matrix
            {
              errcount ++;
              found = FALSE;
              break;
            }
            rc++;                                                      // row counter
            if (*p == '\0')                                            // if the next row starts on a new line, read it
            {
              file.getline(buf, 1023, '\n');
              if (file.fail())
              {
                cout << ERRORMSG3 << ident <<endl;
                found = FALSE;                                              // this error may be fatal!
                break;
              }
              p = strpbrk(buf, "-+0123456789.");
              q = p;
            }
          }
          if (*q != ']') cout << ERRORMSG2 << ident <<endl;            // error message if there is more or less data then there should be
          if ((*r < 0) && (rc < maxrows))
          {
            cout << ERRORMSG3 << ident <<endl;
            found = FALSE;                                              // this error may be fatal!
          }
          *c = cc;
          *r = rc;
          if (errcount)
          {
            cout << ERRORMSG4 << ident <<endl;
            found = FALSE;
          }
        }
      }
    }
    file.close();
  }
  if ((!found) & (strcmp(filename, DEFAULTS) == 1)) cout << "PARAMETER MISSING! Variable " << ident << " not found in file " << filename << endl;

  return found;
}


int readHarvest(char *harvestIN)
{

    int row = 0;
    int count = 0, lines = 0, checkval = 1, r = 0, c = 0;
    string line;
    double *buf;

    ifstream file(harvestIN);
    buf = new double[1000];              
    // maximum number of lines in harvest file
    if (buf == NULL)
    {
      cout << ERRORMSG6 << endl;
      exit(EXIT_FAILURE);
    }

    /*if(file.is_open()){
        while(!file.eof())
        {
            getline(file, line);
            lines++;
        }
    file.close();}
*/
    r = lines; //number of harvest dates (number of lines)
    c = 4; // YY, MM, DD, height 

    count = readmatrix(buf, &r, &c, "HarvestData", harvestIN, TRUE);
    HarvestData.Resize(r, c, buf);
    Harvest_LOC = 0;
    Harvest_LOC_END = r; //number of harvest dates (number of lines)
    if (c != 4) // c must be equal to 4: YY, MM, DD, height
    {
        cout << ERRORMSG13 << endl;
        exit(EXIT_FAILURE);
    }
    checkHarvestDate();

    if ((HYY < StartYear) && (HMM < StartMonth) && (HDD < StartDay)){
        cout << ERRORMSG14 << endl; // First harvest before start date
        cout << "HYY " << HYY << endl;
        cout << "HMM " << HMM << endl;
        cout << "HDD " << HDD << endl;
        exit(EXIT_FAILURE);
    }

    if (count < checkval) {
        cout << "count is less than checkval." << endl;
        cout << "count: " << count << endl;
        cout << "checkval: " << checkval << endl;
    }

    if (count < checkval) return FALSE; else return TRUE;

}


int readall()                     // read all parameters
/* THIS FUNCTION HAS TO BE ADAPTED IF NEW PARAMETERS ARE ADDED TO THE MODEL */
{
  double x;
  double *buf;
  int count = 0, found = FALSE, len, r, c, maxpar = 118;



  buf = new double[10000];
  if (buf == NULL)
  {
    cout << "Memory allocation error data read buffer" << endl;
    exit(EXIT_FAILURE);
  }
  count += readscalar(&x, "NrLayers", DEFAULTS, TRUE);                // read first the fixed parameters of the model configuration and constants
  NrLayers = (int)x;
  count += readscalar(&LayerThickness, "LayerThickness", DEFAULTS, TRUE);
 // read first the model configuration and constants
  count += readscalar(&DensOrg, "DensOrg", DEFAULTS, TRUE);
  count += readscalar(&DensMin, "DensMin", DEFAULTS, TRUE);
  count += readscalar(&HCOrg, "HCOrg", DEFAULTS, TRUE);
  count += readscalar(&DensWater, "DensWater", DEFAULTS, TRUE);
  count += readscalar(&DensIce, "DensIce", DEFAULTS, TRUE);
  count += readscalar(&HCMiner, "HCMiner", DEFAULTS, TRUE);
  count += readscalar(&HCAir, "HCAir", DEFAULTS, TRUE);
  count += readscalar(&HCWater, "HCWater", DEFAULTS, TRUE);
  count += readscalar(&HCIce, "HCIce", DEFAULTS, TRUE);
  count += readscalar(&CondOrg, "CondOrg", DEFAULTS, TRUE);
  count += readscalar(&CondQuartz, "CondQuartz", DEFAULTS, TRUE);
  count += readscalar(&CondMiner, "CondMiner", DEFAULTS, TRUE);
  count += readscalar(&CondAir, "CondAir", DEFAULTS, TRUE);
  count += readscalar(&CondWater, "CondWater", DEFAULTS, TRUE);
  count += readscalar(&CondIce, "CondIce", DEFAULTS, TRUE);
  count += readscalar(&CondSnow, "CondSnow", DEFAULTS, TRUE);
  len = 3;
  count += readarray(LatentHeat.Data(), &len, "LatentHeat", DEFAULTS, TRUE);
  count += readscalar(&Rgas, "Rgas", DEFAULTS, TRUE);
  count += readscalar(&MethaneDiff, "MethaneDiff", DEFAULTS, TRUE);
  count += readscalar(&MethaneDiffWater, "MethaneDiffWater", DEFAULTS, TRUE);
// hereafter all the stuff that may be defined in parameters file is read;
// if not found in the parameters file a second attempt is made in the defaults file
  found = readscalar(&TStepHeat, "TStepHeat", ParamFile, FALSE);
  if (!found) found = readscalar(&TStepHeat, "TStepHeat", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&DStepHeat, "DStepHeat", ParamFile, FALSE);
  if (!found) found = readscalar(&DStepHeat, "DStepHeat", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MaxDepthHeat, "MaxDepthHeat", ParamFile, FALSE);
  if (!found) found = readscalar(&MaxDepthHeat, "MaxDepthHeat", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&VegTScalingFactor, "VegTScalingFactor", ParamFile, FALSE);
  if (!found) found = readscalar(&VegTScalingFactor, "VegTScalingFactor", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&x, "ThermModel", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "ThermModel", DEFAULTS, TRUE);
  count += found;
  ThermModel = (int)x;
  found = readscalar(&x, "ProductionModel", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "ProductionModel", DEFAULTS, TRUE);
  count += found;
  ProductionModel = (int)x;
  if (ProductionModel > 1) maxpar += 1;   // adapt number of parameters to be read depending on productionmodel
  if (ProductionModel > 2)
  {
      maxpar += 4;
      found = readscalar(&AmbientCO2, "AmbientCO2", ParamFile, FALSE);
      if (!found) found = readscalar(&AmbientCO2, "AmbientCO2", DEFAULTS, TRUE);
      len = 8;
      found = readarray(Phenology.Data(), &len, "Phenology", ParamFile, FALSE);
      if (!found) found = readarray(Phenology.Data(), &len, "Phenology", DEFAULTS, TRUE);
      count += found;
      found = readscalar(&Latitude, "Latitude", ParamFile, FALSE);
      if (!found) found = readscalar(&Latitude, "Latitude", DEFAULTS, TRUE);
      count += found;
      found = readscalar(&KBeer, "KBeer", ParamFile, FALSE);
      if (!found) found = readscalar(&KBeer, "KBeer", DEFAULTS, TRUE);
      count += found;
      found = readscalar(&LAICarbonFraction, "LAICarbonFraction", ParamFile, FALSE);
      if (!found) found = readscalar(&LAICarbonFraction, "LAICarbonFraction", DEFAULTS, TRUE);
      count += found;
      found = readscalar(&x, "PARunits", ParamFile, FALSE);
      if (!found) found = readscalar(&x, "PARunits", DEFAULTS, TRUE);
      count += found;
      PARunits = (int)x;
  }
  if (ProductionModel == 4)
  {
      maxpar += 1;
      len = 4;
      found = readarray(PhotoPar.Data(), &len, "PhotoPar", ParamFile, FALSE);
      if (!found) found = readarray(PhotoPar.Data(), &len, "PhotoPar", DEFAULTS, TRUE);
      count += found;
  }
  found = readscalar(&GreenBiomassRatio, "GreenBiomassRatio", ParamFile, FALSE);
  if (!found) found = readscalar(&GreenBiomassRatio, "GreenBiomassRatio", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&Timestep, "Timestep", ParamFile, FALSE);
  if (!found) found = readscalar(&Timestep, "Timestep", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&x, "NrOfSteps", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "NrOfSteps", DEFAULTS, TRUE);
  NrOfSteps = (int)x;
    found = readscalar(&x, "StartYear", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "StartYear", DEFAULTS, TRUE);
  StartYear = x;
  count += found;

  found = readscalar(&x, "StartDay", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "StartDay", DEFAULTS, TRUE);
  StartDay = (int)x;
  count += found;

  found = readscalar(&DissimAssimRatio, "DissimAssimRatio", ParamFile, FALSE);
  if (!found) found = readscalar(&DissimAssimRatio, "DissimAssimRatio", DEFAULTS, TRUE);
  count += found;
    found = readscalar(&AnaerobicDARatio, "AnaerobicDARatio", ParamFile, FALSE);
  if (!found) found = readscalar(&AnaerobicDARatio, "AnaerobicDARatio", DEFAULTS, TRUE);

  found = readstring(StartDate, "StartDate", ParamFile, FALSE);
  if (!found) found = readstring(StartDate, "StartDate", DEFAULTS, TRUE);
  count += found;

  found = readstring(EndDate, "EndDate", ParamFile, FALSE);
  if (!found) found = readstring(EndDate, "EndDate", DEFAULTS, TRUE);
  count += found; 

  if ( (Timestep == 1.0) && (Assign_Start(StartDate)) && (Assign_end(EndDate))) 
  {
    NrOfSteps = count_days(StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear);
  }

  if (!NrOfSteps) cout << "Simulation duration could not be calculated. Default assumed." << endl;
  if (Verbose) cout << "Total number of model timesteps (days): " << NrOfSteps << endl;

  found = readscalar(&ResistFrac, "ResistFrac", ParamFile, FALSE);
  if (!found) found = readscalar(&ResistFrac, "ResistFrac", DEFAULTS, TRUE);
  count += found;

  found = readscalar(&KLitter, "KLitter", ParamFile, FALSE);
  if (!found) found = readscalar(&KLitter, "KLitter", DEFAULTS, TRUE);
  count += found;
  len = 7;
  found = readarray(Cfrac.Data(), &len, "Cfrac", ParamFile, FALSE);
  if (!found) found = readarray(Cfrac.Data(), &len, "Cfrac", DEFAULTS, TRUE);
  count += found;
  r = 2;
  c = 2;
  found = readmatrix(pFpoints.Data(), &r, &c, "pFpoints", ParamFile, FALSE);
  if (!found) found = readmatrix(pFpoints.Data(), &r, &c, "pFpoints", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&HalfSatPoint, "HalfSatPoint", ParamFile, FALSE);
  if (!found) found = readscalar(&HalfSatPoint, "HalfSatPoint", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&RootAeration, "RootAeration", ParamFile, FALSE);
  if (!found) found = readscalar(&RootAeration, "RootAeration", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&PrimingCorrection, "PrimingCorrection", ParamFile, FALSE);
  if (!found) found = readscalar(&PrimingCorrection, "PrimingCorrection", DEFAULTS, TRUE);
  count += found;
  len = 7;
  found = readarray(Kdecay.Data(), &len, "Kdecay", ParamFile, FALSE);
  if (!found) found = readarray(Kdecay.Data(), &len, "Kdecay", DEFAULTS, TRUE);
  count += found;

  len = 7;
  found = readarray(AerobicQ10.Data(), &len, "AerobicQ10", ParamFile, FALSE);
  if (!found) found = readarray(AerobicQ10.Data(), &len, "AerobicQ10", DEFAULTS, TRUE);
  count += found;
  len = 2;
/*  found = readarray(KPeatCN.Data(), &len, "KPeatCN", ParamFile, FALSE);
  if (!found) found = readarray(KPeatCN.Data(), &len, "KPeatCN", DEFAULTS, TRUE);
  count += found; */

  found = readscalar(&ShootsFactor, "ShootsFactor", ParamFile, FALSE);
  if (!found) found = readscalar(&ShootsFactor, "ShootsFactor", DEFAULTS, TRUE);
  count += found;
  len = 2;
  found = readarray(RespFac.Data(), &len, "RespFac", ParamFile, FALSE);
  if (!found) found = readarray(RespFac.Data(), &len, "RespFac", DEFAULTS, TRUE);
  count += found;
  len = 2;
  found = readarray(ProdTFunc.Data(), &len, "ProdTFunc", ParamFile, FALSE);
  if (!found) found = readarray(ProdTFunc.Data(), &len, "ProdTFunc", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&SatCorr, "SatCorr", ParamFile, FALSE);
  if (!found) found = readscalar(&SatCorr, "SatCorr", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&GrowFuncConst, "GrowFuncConst", ParamFile, FALSE);
  if (!found) found = readscalar(&GrowFuncConst, "GrowFuncConst", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&SpringCorrection, "SpringCorrection", ParamFile, FALSE);
  if (!found) found = readscalar(&SpringCorrection, "SpringCorrection", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MaxProd, "MaxProd", ParamFile, FALSE);
  if (!found) found = readscalar(&MaxProd, "MaxProd", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MinProd, "MinProd", ParamFile, FALSE);
  if (!found) found = readscalar(&MinProd, "MinProd", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MaxRootDepth, "MaxRootDepth", ParamFile, FALSE);
  if (!found) found = readscalar(&MaxRootDepth, "MaxRootDepth", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&x, "NoRootsBelowGWT", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "NoRootsBelowGWT", DEFAULTS, TRUE);
  NoRootsBelowGWT = (int)x;
  count += found;
  found = readscalar(&RootLambda, "RootLambda", ParamFile, FALSE);
  if (!found) found = readscalar(&RootLambda, "RootLambda", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&RootSenescence, "RootSenescence", ParamFile, FALSE);
  if (!found) found = readscalar(&RootSenescence, "RootSenescence", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&InitRoots, "InitRoots", ParamFile, FALSE);
  if (!found) found = readscalar(&InitRoots, "InitRoots", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&ExudateFactor, "ExudateFactor", ParamFile, FALSE);
  if (!found) found = readscalar(&ExudateFactor, "ExudateFactor", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&BioMass, "BioMass", ParamFile, FALSE);
  if (!found) found = readscalar(&BioMass, "BioMass", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&BioMassSenescence, "BioMassSenescence", ParamFile, FALSE);
  if (!found) found = readscalar(&BioMassSenescence, "BioMassSenescence", DEFAULTS, TRUE);
  count += found;
  len = 2;
  found = readarray(HarvestCorrection.Data(), &len, "HarvestCorrection", ParamFile, FALSE);
  if (!found) found = readarray(HarvestCorrection.Data(), &len, "HarvestCorrection", DEFAULTS, TRUE);
  count += found;
  r = 0;
  c = 2;
  found = readmatrix(buf, &r, &c, "Harvest", ParamFile, FALSE);
  if (!found) found = readmatrix(buf, &r, &c, "Harvest", DEFAULTS, TRUE);
  Harvest.Resize(r, c, buf);
  count += found;
  found = readstring(HarvestFile, "HarvestFile", ParamFile, FALSE);
  if (!found) found = readstring(HarvestFile, "HarvestFile", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&HarvestLitter, "HarvestLitter", ParamFile, FALSE);
  if (!found) found = readscalar(&HarvestLitter, "HarvestLitter", DEFAULTS, TRUE);
  count += found;
  r = 0;
  c = 4;
  found = readmatrix(buf, &r, &c, "Grazing", ParamFile, FALSE);
  if (!found) found = readmatrix(buf, &r, &c, "Grazing", DEFAULTS, TRUE);
  Grazing.Resize(r, c, buf);
  count += found;
  r = 0;
  c = 2;
  found = readmatrix(buf, &r, &c, "Manure", ParamFile, FALSE);
  if (!found) found = readmatrix(buf, &r, &c, "Manure", DEFAULTS, TRUE);
  Manure.Resize(r, c, buf);
  count += found;
  found = readscalar(&ManureFluidFrac, "ManureFluidFrac", ParamFile, FALSE);
  if (!found) found = readscalar(&ManureFluidFrac, "ManureFluidFrac", DEFAULTS, TRUE);
  count += found;
  r = 0;
  c = 2;
  found = readmatrix(buf, &r, &c, "ManureLayers", ParamFile, FALSE);
  if (!found) found = readmatrix(buf, &r, &c, "ManureLayers", DEFAULTS, TRUE);
  ManureLayers.Resize(r, c, buf);
  count += found;
  len = 7;
  found = readarray(MethaneReservoirs.Data(), &len, "MethaneReservoirs", ParamFile, FALSE);
  if (!found) found = readarray(MethaneReservoirs.Data(), &len, "MethaneReservoirs", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneR0, "MethaneR0", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneR0, "MethaneR0", DEFAULTS, TRUE);
  count += found;

  found = readscalar(&MethanepHCorr, "MethanepHCorr", ParamFile, FALSE);
  if (!found) found = readscalar(&MethanepHCorr, "MethanepHCorr", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneQ10, "MethaneQ10", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneQ10, "MethaneQ10", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneOxQ10, "MethaneOxQ10", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneOxQ10, "MethaneOxQ10", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneVmax, "MethaneVmax", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneVmax, "MethaneVmax", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneKm, "MethaneKm", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneKm, "MethaneKm", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneMaxConc, "MethaneMaxConc", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneMaxConc, "MethaneMaxConc", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneERateC, "MethaneERateC", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneERateC, "MethaneERateC", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethanePRateC, "MethanePRateC", ParamFile, FALSE);
  if (!found) found = readscalar(&MethanePRateC, "MethanePRateC", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneAir, "MethaneAir", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneAir, "MethaneAir", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethaneTRef, "MethaneTRef", ParamFile, FALSE);
  if (!found) found = readscalar(&MethaneTRef, "MethaneTRef", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethanePType, "MethanePType", ParamFile, FALSE);
  if (!found) found = readscalar(&MethanePType, "MethanePType", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MethanePlantOx, "MethanePlantOx", ParamFile, FALSE);
  if (!found) found = readscalar(&MethanePlantOx, "MethanePlantOx", DEFAULTS, TRUE);
  count += found;

  found = readscalar(&CO2CH4ratio, "CO2CH4ratio", ParamFile, FALSE);
  if (!found) found = readscalar(&CO2CH4ratio, "CO2CH4ratio", DEFAULTS, TRUE);
  count += found;

  found = readscalar(&PartialAnaerobe, "PartialAnaerobe", ParamFile, FALSE);
  if (!found) found = readscalar(&PartialAnaerobe, "PartialAnaerobe", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&AnaerobeLagFactor, "AnaerobeLagFactor", ParamFile, FALSE);
  if (!found) found = readscalar(&AnaerobeLagFactor, "AnaerobeLagFactor", DEFAULTS, TRUE);
  count += found;
  len = 0;
  found = readarray(buf, &len, "InitMethane", ParamFile, FALSE);
  if (!found) found = readarray(buf, &len, "InitMethane", DEFAULTS, TRUE);
  InitMethane.Resize(len, buf);
  count += found;
  found = readscalar(&DayMinGW, "DayMinGW", ParamFile, FALSE);
  if (!found) found = readscalar(&DayMinGW, "DayMinGW", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MinGW, "MinGW", ParamFile, FALSE);
  if (!found) found = readscalar(&MinGW, "MinGW", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&AmplitudeGW, "AmplitudeGW", ParamFile, FALSE);
  if (!found) found = readscalar(&AmplitudeGW, "AmplitudeGW", DEFAULTS, TRUE);
  count += found;
  found = readstring(GwFile, "GwFile", ParamFile, FALSE);
  if (!found) found = readstring(GwFile, "GwFile", DEFAULTS, TRUE);
  count += found;
   found = readstring(PrecipFile, "PrecipFile", ParamFile, FALSE);
  if (!found) found = readstring(PrecipFile, "PrecipFile", DEFAULTS, TRUE);
  count += found;
  found = readstring(EvapFile, "EvapFile", ParamFile, FALSE);
  if (!found) found = readstring(EvapFile, "EvapFile", DEFAULTS, TRUE);
  count += found;
  found = readstring(DrainageFile, "DrainageFile", ParamFile, FALSE);
  if (!found) found = readstring(DrainageFile, "DrainageFile", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&T_average, "T_average", ParamFile, FALSE);
  if (!found) found = readscalar(&T_average, "T_average", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&T_amplitude, "T_amplitude", ParamFile, FALSE);
  if (!found) found = readscalar(&T_amplitude, "T_amplitude", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&T_ref, "T_ref", ParamFile, FALSE);
  if (!found) found = readscalar(&T_ref, "T_ref", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&MaxSnowdepth, "MaxSnowdepth", ParamFile, FALSE);
  if (!found) found = readscalar(&MaxSnowdepth, "MaxSnowdepth", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&DayMaxSnowdepth, "DayMaxSnowdepth", ParamFile, FALSE);
  if (!found) found = readscalar(&DayMaxSnowdepth, "DayMaxSnowdepth", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&SnowMeltrate, "SnowMeltrate", ParamFile, FALSE);
  if (!found) found = readscalar(&SnowMeltrate, "SnowMeltrate", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&WatertableInit, "WatertableInit", ParamFile, FALSE);
  if (!found) found = readscalar(&WatertableInit, "WatertableInit", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&EvapCorrection, "EvapCorrection", ParamFile, FALSE);
  if (!found) found = readscalar(&EvapCorrection, "EvapCorrection", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&RunoffThreshold, "RunoffThreshold", ParamFile, FALSE);
  if (!found) found = readscalar(&RunoffThreshold, "RunoffThreshold", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&OpenWaterFactor, "OpenWaterFactor", ParamFile, FALSE);
  if (!found) found = readscalar(&OpenWaterFactor, "OpenWaterFactor", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&CropFactor, "CropFactor", ParamFile, FALSE);
  if (!found) found = readscalar(&CropFactor, "CropFactor", DEFAULTS, TRUE);
  count += found;  
  found = readscalar(&DrainageDist, "DrainageDist", ParamFile, FALSE);
  if (!found) found = readscalar(&DrainageDist, "DrainageDist", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&Ksat, "Ksat", ParamFile, FALSE);
  if (!found) found = readscalar(&Ksat, "Ksat", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&DrainLevel, "DrainLevel", ParamFile, FALSE);
  if (!found) found = readscalar(&DrainLevel, "DrainLevel", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&ThermDiff, "ThermDiff", ParamFile, FALSE);          // thermal diffusion needs not be defined
  if (!found) found = readscalar(&ThermDiff, "ThermDiff", DEFAULTS, TRUE);
  r = 2;
  c = 0;
  found = readmatrix(buf, &r, &c, "T_init", ParamFile, FALSE);
  if (!found) found = readmatrix(buf, &r, &c, "T_init", DEFAULTS, TRUE);
  T_init.Resize(r, c, buf);
  if (c > 0) count += found;
  found = readstring(TFile, "TFile", ParamFile, FALSE);
  if (!found) found = readstring(TFile, "TFile", DEFAULTS, TRUE);
  found = readstring(SnowFile, "SnowFile", ParamFile, FALSE);
  if (!found) found = readstring(SnowFile, "SnowFile", DEFAULTS, TRUE);
  if (ProductionModel == 2)
  {
      found = readstring(NPPFile, "NPPFile", ParamFile, FALSE);
      if (!found) found = readstring(NPPFile, "NPPFile", DEFAULTS, TRUE);
      count += found;
  }
  if (ProductionModel > 2)
  {  
      found = readstring(PARFile, "PARFile", ParamFile, FALSE);
      if (!found) found = readstring(PARFile, "PARFile", DEFAULTS, TRUE);
      count += found;
  }
  found = readstring(CO2File, "CO2File", ParamFile, FALSE);
  if (!found) found = readstring(CO2File, "CO2File", DEFAULTS, TRUE);
  count += found;
  found = readstring(RunOnFile, "RunOnFile", ParamFile, FALSE);
  if (!found) found = readstring(RunOnFile, "RunOnFile", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&LitterLayer, "LitterLayer", ParamFile, FALSE);
  if (!found) found = readscalar(&LitterLayer, "LitterLayer", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&LitterConversion, "LitterConversion", ParamFile, FALSE);
  if (!found) found = readscalar(&LitterConversion, "LitterConversion", DEFAULTS, TRUE);
  count += found;
  found = readscalar(&x, "AnaerobicCO2", ParamFile, FALSE);
  if (!found) found = readscalar(&x, "AnaerobicCO2", DEFAULTS, TRUE);
  count += found;
  AnaerobicCO2 = (int)x;
  found = readscalar(&Q10Anaerobic, "Q10Anaerobic", ParamFile, FALSE);
  if (!found) found = readscalar(&Q10Anaerobic, "Q10Anaerobic", DEFAULTS, TRUE);
  count += found;
  len = 7;
  found = readarray(KAnaerobic.Data(), &len, "KAnaerobic", ParamFile, FALSE);
  if (!found) found = readarray(KAnaerobic.Data(), &len, "KAnaerobic", DEFAULTS, TRUE);
  count += found;
  found = readstring(SoilProfile, "SoilProfile", ParamFile, FALSE);
  if (!found) found = readstring(SoilProfile, "SoilProfile", DEFAULTS, TRUE);
  count += found;
  found = readstring(SoilMoisture, "SoilMoisture", ParamFile, FALSE);
  if (!found) found = readstring(SoilMoisture, "SoilMoisture", DEFAULTS, TRUE);
  len = 0;
  found = readarray(buf, &len, "ProfileOutput", ParamFile, FALSE);
  if (!found) found = readarray(buf, &len, "ProfileOutput", DEFAULTS, TRUE);
  ProfileOutput.Resize(len, buf);
  count += found;
  delete buf;
  if (count >= maxpar) {
      if (Verbose) cout << "Number of parameters read: " << count << "; required " << maxpar << endl << endl;
      return TRUE;
    } else {
	if (Verbose)
	{ 
		cout << "Number of parameters read: " << count << "; should be at least " << maxpar << "!" << endl;
	}
	return FALSE;
  }
}

int readsoil(char *soilname)          // reads soil profile from the file specified in soilname
{
  double x, yr;
  int count = 0, len, r, c, nyears;
  double *buf;


  buf = new double[MAXHORIZON];              // maximum number of horizons
  if (buf == NULL)
  {
    cout << ERRORMSG6 << endl;
    exit(EXIT_FAILURE);
  }
  count += readscalar(&x, "NrHorizons", soilname, TRUE);        // read nr of horizons first, then the rest
  NrHorizons = (int)x;
  if (NrHorizons > MAXHORIZON)
  {
    cout << ERRORMSG6 << MAXHORIZON << endl;
    exit(EXIT_FAILURE);
  }
  len = NrHorizons;
  count += readarray(buf, &len, "Horizons", soilname, TRUE);
  Horizons.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "CNRatio", soilname, TRUE);
  CNRatio.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "DBD", soilname, TRUE);
  DBD.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "PercOrg", soilname, TRUE);
  PercOrg.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "SandFraction", soilname, TRUE);
  SandFraction.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "ClayFraction", soilname, TRUE);
  ClayFraction.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "FreezingCurve", soilname, TRUE);
  FreezingCurve.Resize(NrHorizons, buf);
  count += readarray(buf, &len, "Layer_pH", soilname, TRUE);
  Layer_pH.Resize(NrHorizons, buf);
  r = NrHorizons;
  c = NrReservoirs;
  count += readmatrix(buf, &r, &c, "InitRes", soilname, TRUE);
  InitRes.Resize(r, c, buf);
  len = 0;
  count += readarray(buf, &len, "pFVal", soilname, TRUE);
  pFVal.Resize(len, buf);
  c = 0;
  count += readmatrix(buf, &r, &c, "Layer_pF", soilname, TRUE);
  Layer_pF.Resize(r, c, buf);
  if (!((c == 5) || (c == len)))
  {
    cout << ERRORMSG5 << endl;
    count--;
  }
  len = NrHorizons;
  if (readarray(buf, &len, "Porosity", soilname, FALSE)) Porosity.Resize(len, buf);
  if ((ThermModel == 0) && (strlen(TFile) != 0))            // read temperature time series
  {
    TData.Resize(NrOfSteps);
    if (readTseries(TData.Data(), TFile, NrOfSteps, 1) < NrOfSteps)
    {
      cout << ERRORMSG9 << TFile << " Surface temperature time series" << endl;
    } else
    {
      if (Verbose) cout << "Taking temperature data from file: " << TFile <<"\n";
    }
  }
  if ((ThermModel == 0) && (strlen(SnowFile) != 0))            // read snowdepth time series
  {
    SnowData.Resize(NrOfSteps);
    if (readTseries(SnowData.Data(), SnowFile, NrOfSteps, 1) < NrOfSteps)
    {
      cout << ERRORMSG9 << SnowFile << " Snow depth time series" << endl;
      count--;
    } else
    {
      if (Verbose) cout << "Taking snow data from file: " << SnowFile <<"\n";
    }
  }
  if (ProductionModel == 2)         // read primary production time series
  {
    NPPData.Resize(NrOfSteps);
    if (readTseries(NPPData.Data(), NPPFile, NrOfSteps, 1) < NrOfSteps)
    {
      cout << ERRORMSG9 << NPPFile << endl;
      count--;
    } else
    {
      if (Verbose) cout << "Taking NPP data from file: " << NPPFile  << " NPP time series" <<"\n";
    }
  }
  if (ProductionModel > 2)         // read PAR data time series
  {
        PARData.Resize(NrOfSteps);
        if (readTseries(PARData.Data(), PARFile, NrOfSteps, 1) < NrOfSteps)
        {
            cout << ERRORMSG9 << PARFile << endl;
            count--;
        } else
        {
            if (Verbose) cout << "Taking PAR data from file: " << PARFile  << " PAR time series" <<"\n";
        }
  }
  if ((ThermModel == 2) && (strlen(TFile) != 0))            // read temperature time series in case when temperature is specified for each layer
  {
    TData.Resize(NrOfSteps, NrLayers);
    if (readTseries(TData.Data(), TFile, NrOfSteps, NrLayers) < NrOfSteps)
    {
      cout << ERRORMSG9 << TFile << " Soil temperature time series" << endl;
      count--;  
    } else
    {
      if (Verbose) cout << "Taking soil temperature data from file: " << TFile <<"\n";
    }
  }
  if (strlen(GwFile) != 0)                                    // read groundwater table time series
  {
    GwData.Resize(NrOfSteps);
    if (readTseries(GwData.Data(), GwFile, NrOfSteps, 1) < NrOfSteps)
    {
      cout << ERRORMSG9 << GwFile << endl;
      count--;
    } else
    {
      MinGW = GwData.Min();
	  WatertableModel = 0;
	  if (Verbose) cout << "Taking water table data from file: " << GwFile << " Water table time series" << "\n";
    }
  }

  if (strlen(HarvestFile) != 0)                                 // read harvest time series
  {
      if (!readHarvest(HarvestFile)){
          cout << ERRORMSG9 << HarvestFile << endl;
          count--;
      } else
      {
      HarvestModel = 2;
      if (Verbose) cout << "Taking harvest data from file: " << HarvestFile << "\n" << endl;
      }
      
  }

  

  if (strlen(PrecipFile) != 0)                                    // read precipitation time series, this also initilizes the water table module
  {
    Precipitation.Resize(NrOfSteps);
    if (readTseries(Precipitation.Data(), PrecipFile, NrOfSteps, 1) < NrOfSteps)
    {
      cout << ERRORMSG9 << PrecipFile << endl;
      count--;
    } else
    {
	  WatertableModel = 2;										// Set water table to calculated water table
      if (Verbose)
	  {
		cout << "Calculating water table from evaporation and precipitation\n";
		cout << "Taking precipitation data from file: " << PrecipFile << " for water table module" << "\n";
	  }
    }
        
  }
  if (strlen(EvapFile) != 0)                                    // read evaporation time series for water table module
  {
    Evaporation.Resize(NrOfSteps);
    if (readTseries(Evaporation.Data(), EvapFile, NrOfSteps, 1) < NrOfSteps)
    {
      cout << ERRORMSG9 << EvapFile << endl;
      count--;
    } else
    {
      if (Verbose) cout << "Taking evaporation data from file: " << EvapFile << " for water table module" << "\n";
    }
  }
  if (strlen(SoilMoisture) != 0)                            // read soil moisture time series
  {
    MoistProfiles.Resize(NrOfSteps, NrLayers);
    if (readTseries(MoistProfiles.Data(), SoilMoisture, NrOfSteps, NrLayers) < NrOfSteps)
    {
      cout << ERRORMSG9 << SoilMoisture << " Soil moisture time series" << endl;
      count--;
    } else if (Verbose) cout << "Taking soil moisture data from file: " << SoilMoisture <<"\n";
  }
  if (strlen(DrainageFile) != 0)                            // read soil moisture time series
  {
        DrainData.Resize(NrOfSteps);
        if (readTseries(DrainData.Data(), DrainageFile, NrOfSteps, 1) < NrOfSteps)
        {
            cout << ERRORMSG9 << DrainageFile << " River/ditch water level time series" << endl;
            count--;
        } else if (Verbose) cout << "Taking River/ditch water level data from file: " << DrainageFile <<"\n";
  }
  if (strlen(RunOnFile) != 0)                            // read soil moisture time series
  {
        RunOn.Resize(NrOfSteps);
        if (readTseries(RunOn.Data(), RunOnFile, NrOfSteps, 1) < NrOfSteps)
        {
            cout << ERRORMSG9 << RunOnFile << " Run-on water time series" << endl;
            count--;
        } else if (Verbose) cout << "Taking run-on water data from file: " << RunOnFile <<"\n";
  }
  if ((strlen(CO2File) != 0) && (ProductionModel > 2))     // read CO2 input data
    {
        yr = ceil(NrOfSteps / 364);
        nyears = int(yr) + 1;
        CO2Data.Resize(nyears);
        if (readTseries(CO2Data.Data(), CO2File, nyears, 1) < nyears)
        {
            cout << ERRORMSG9 << CO2File << endl;
            count--;
        } else
        {
            if (Verbose) cout << "Taking CO2 data from file: " << CO2File << " for photosynthesis module" << "\n";
        }
    }
  if (count <11) return FALSE; else return TRUE;
}

int readTseries(double *d, const char *filename, const int maxrows, const int cols)
/*  Reads a time series from file, e.g. groundwater level or temperature data
    This may be a multivariate time series, with the variable stored in cols.
    Each line in the file represents observations on one time point
    Tseries : Matrix that stores the data
    filename: the name of the file where the series is stored
    maxrows : the maximum number of data rows to be read
    cols    : the number of data columns                          */
{
  char buf[1024];
  char nr[1024];
  BOOLEAN cont = TRUE;
  int count = 0, rowcount = 0, i;
  char *p, *q;

  strcpy(buf, DataDir);
  ifstream file(strcat(buf,filename));                              // open file
  if (!file)
  {
    cout <<  ERRORMSG1 << " " << filename  << endl;
    return FALSE;
  } else {
    while (cont)
    {
      file.getline(buf, 1023, '\n');
      if (file.fail())
      {
        cont = FALSE;
      }else
      {
        if (cols == 1)                                            // shorter read procedure with one column
        {
          *d++ = atof(buf);
          count++;
        } else
        {
          p = buf;                                                 // start reading numbers
          while (*p == ' ') p++;
          for (i = 0; i < cols; i++)
          {                                                        // find the end of the number
            q = strpbrk(p, " ';\n");
            if (q == NULL) strcpy(nr, p); else
            {
              strncpy(nr, p, q - p);
              nr[q - p] = 0;                                         // add the NULL character
            }
            if (strlen(nr) > 0)
            {
              *d++ = atof(nr);
              count++;
            }
            if ((count > (maxrows * cols)) || ((q == NULL) && (i < (cols - 1))))
            {
              cont = FALSE;
              cout << ERRORMSG8 << endl;
              break;
            } else
            if (q != NULL)
            {
              p = q;
              while (*p == ' ') p++;
            }
          }
        }
      }
      rowcount++;
      if (rowcount == maxrows) cont = FALSE;
    }
    file.close();
  }
  return count;
}


void WriteOutput()
/* writes methane fluxes, CO2 fluxes and plant production/respiration and carbon balance to output files */

{
  char buf[1024];
  ofstream *output;

  strcpy(buf, &DataDir[0]);
  strcat(buf, &OutputFilePrefix[0]);
  output = new ofstream(strcat(buf, OUTPUT_METHANE));
  if (!output)
  {
    cout  << ERRORMSG11 << OUTPUT_METHANE << endl;
    exit(EXIT_FAILURE);
  }
  TotalMethane.Write(output);
  output->close();
  strcpy(buf, &DataDir[0]);
  strcat(buf, &OutputFilePrefix[0]);
  output = new ofstream(strcat(buf, OUTPUT_CO2RES));
  if (!output)
  {
    cout  << ERRORMSG11 << OUTPUT_CO2RES << endl;
    exit(EXIT_FAILURE);
  }
  ReservoirTime.Write(output);
  output->close();
  strcpy(buf, &DataDir[0]);
  strcat(buf, &OutputFilePrefix[0]);
  output = new ofstream(strcat(buf, OUTPUT_CO2LAY));
  if (!output)
  {
    cout  << ERRORMSG11 << OUTPUT_CO2LAY << endl;
    exit(EXIT_FAILURE);
  }
  LayerTime.Write(output);
  output->close();
  strcpy(buf, &DataDir[0]);
  strcat(buf, &OutputFilePrefix[0]);
  output = new ofstream(strcat(buf, OUTPUT_BIO));
  if (!output)
  {
    cout  << ERRORMSG11 << OUTPUT_BIO << endl;
    exit(EXIT_FAILURE);
  }
  BioMassRec.Write(output);
  output->close();
  if (AnaerobicCO2 > 0) {
    strcpy(buf, &DataDir[0]);
    strcat(buf, &OutputFilePrefix[0]);
    output = new ofstream(strcat(buf, OUTPUT_ANAEROB));
    if (!output)
    {
        cout  << ERRORMSG11 << OUTPUT_ANAEROB << endl;
        exit(EXIT_FAILURE);
    }
    LayerAnaerobic.Write(output);
    output->close();
  }

  strcpy(buf, &DataDir[0]);
  strcat(buf, &OutputFilePrefix[0]);
  output = new ofstream(strcat(buf, OUTPUT_BALANCE));
  if (!output)
  {
    cout  << ERRORMSG11 << OUTPUT_BALANCE << endl;
    exit(EXIT_FAILURE);
  }
  CarbonBalance.Write(output);
  output->close();
  
  strcpy(buf, &DataDir[0]);
  strcat(buf, &OutputFilePrefix[0]);
  output = new ofstream(strcat(buf, OUTPUT_PEATDECOMP));
  if (!output)
  {
    cout  << ERRORMSG11 << OUTPUT_PEATDECOMP << endl;
    exit(EXIT_FAILURE);
  }
  PeatDecay.Write(output);
  output->close();
}



int assignYEARdays(int year)
{
  
  if (isALeapYear(year)) return  366;

  else return YEARdays = 365;

}


int daysinmonth(int month, int year)
{
  switch (month)
  {

  case 1: case 3: case 5: case 7: case 8: case 10: case 12: return 31;

  case 4: case 6: case 9: case 11: return 30;

  case 2: if (isALeapYear(year))
        return 29;
      else
        return 28;

  
  default: cout << "That's not a valid month. " << month << endl;
  }

  return 0;
}


bool isALeapYear(int year)
{
  /* Check if the year is divisible by 4 or 
  is divisible by 400 */
  return ((year % 4 == 0 && year % 100 != 0) || ( year % 400 == 0));
}


void isFebruary(bool LeapYear, int daystring)
{
    if (((LeapYear) && (daystring > 29)) || ((!LeapYear) && (daystring > 28)))
    {
      cout <<  ERRORMSG12 << " " << daystring  << endl;
      exit(EXIT_FAILURE);
    }
}


vector<string> FindDate(string IncomingDate)
{

    int daystring, monthstring, yearlength, yearstring;
    bool LeapYear;

    string str = IncomingDate;
    if ( str.empty() ){
        cout <<  ERRORMSG12 << " " << str  << endl;
        exit(EXIT_FAILURE);
    }

    vector<string> splittedString = split(str, '/');
    daystring = std::stoi(splittedString[0]);
    monthstring = std::stoi(splittedString[1]);
    yearstring = std::stoi(splittedString[2]);
    yearlength = splittedString[2].length();

    // Check days, months <12, year is 4 characters:
    if ((daystring > daysinmonth(monthstring, yearstring)) || (monthstring > 12) || (yearlength != 4))
    {
        cout <<  ERRORMSG12 << " " << str  << endl;
        exit(EXIT_FAILURE);
    }
    
    // test for leap year:
    LeapYear = isALeapYear(yearstring);
    if (monthstring == 2) isFebruary(LeapYear, daystring);
    
    return splittedString;

}

bool Assign_Start(string StartDate)
{
      vector<string> splittedString = FindDate(StartDate);
      StartDay = std::stoi( splittedString[0] );
      StartMonth = std::stoi( splittedString[1] );
      StartYear = std::stoi( splittedString[2] );

    // for(int i = 0; i < splittedString.size() ; i++){
    //     cout << i << " " << splittedString[i] << endl;
    // }
    if (Verbose)
    {
      cout << "StartYear: " << StartYear << endl;
      cout << "StartMonth: " << StartMonth << endl;
      cout << "StartDay: " << StartDay << endl;
    }

  // cout << "StartDay, StartMonth, StartYear: " << StartDay << " " << StartMonth << " " << StartYear << endl;  
      return TRUE;
}


void Find_date(string StrDate)
{

      string str = StrDate;
      vector<string> splittedString = FindDate(StrDate);
      HDD = std::stoi( splittedString[0] );
      HMM = std::stoi( splittedString[1] );
      HYY = std::stoi( splittedString[2] );

      if (Verbose)
      {
        cout << "Year: " << HDD << endl;
        cout << "Month: " << HMM << endl;
        cout << "Day: " << HYY << endl;
      }

}


bool Assign_end(string EndDate)
{
      string str = EndDate;
      vector<string> splittedString = FindDate(EndDate);
      EndDay = std::stoi( splittedString[0] );
      EndMonth = std::stoi( splittedString[1] );
      EndYear = std::stoi( splittedString[2] );

      if (Verbose)
      {
        cout << "EndYear: " << EndYear << endl;
        cout << "EndMonth: " << EndMonth << endl;
        cout << "EndDay: " << EndDay << endl;
        
        
      }

      return TRUE;
}


int count_days(int day_1, int month_1, int year_1, int day_2, int month_2, int year_2)
{

  int daysbetween(0);

    if (year_1 == year_2 && month_1 == month_2) // Same year and month
    {
      daysbetween = (day_2 - day_1) + 1;
    }

    else if (year_1 == year_2) // Same year
    {
      daysbetween = daysinmonth(month_1, year_1) - day_1 + 1;
      // daysbetween = daysinmonth(month_1, year_1) - day_1;
      for (int i = month_1 + 1; i < month_2; ++i)
      {
        daysbetween += daysinmonth(i, year_1);
      }
      daysbetween += day_2;
    }

    else // Different month & different year
    {
      daysbetween = daysinmonth(month_1, year_1) - day_1 +1;
      // daysbetween = daysinmonth(month_1, year_1) - day_1;
      for (int i = month_1 + 1; i <= 12; i++)
      {
        daysbetween += daysinmonth(i, year_1);
      }
      for (int i = (year_1 + 1); i < year_2; i++)
      {
        daysbetween += assignYEARdays(i);
      }
      for (int i = 1; i < month_2; i++)
      {
        daysbetween += daysinmonth(i, year_2);
      }
      daysbetween += day_2;
    }

    return daysbetween;
    
}


int count_years(int year_1, int year_2)
{

  if (year_2 < year_1) 
    {
      cout << "StartYear must be older than EndYear." << endl; 
      exit(EXIT_FAILURE);
    }

  if (year_1 == year_2) return 1;

  else return (year_2 - year_1 +1); 

}



