/***************************************************************************
       SOMproduction.cpp  -  soil organic matter production PEATLAND
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

  Bug correction:incorrect calculation of roods added

  September 15, 2004

  Bug fix temperature dependent primary production in function TemperatureProd()
  tmax was calculated incorrectly, giving zero production at the high end of the
  temperature range

  Added in OrgProd(): Provisions for using net primary production data from
  external model by reading data from file

  February 2007
  Added to organic matter production: temperature relation cf Thornley (1988) cited by Mueller (2000) 
  in 'Modelling Soil-Biosphere interactions'
 
  June-July 2012
  added LPJ-derived photosynthesis;
  bug correction in biomass/total respiration recording
  added recording of above-ground litter/dead biomass
 
  July 2013
  added photosynthesis for tundra, Shaver et al, J. Ecology 2007
  made Beers law parameter variable
  
  January 2022
  Bug correction TundraProd (photosynthesis for tundra, Shaver et al, J. Ecology 2007):
  output corrected to kg C/m2/day instead of kg CO2/m2/day in line with other production models
 ***************************************************************************/

#include <cmath>
#include <cstring>
#include "matrix.h"
#include "general.h"
#include "SOMproduction.h"



void OrgProd()
/* Net Primary Production and its partitioning among roots and shoots */
{
  int i, toproots = 0;
  double gw, belowgwt = 0, abovegwt = 0, rootsadded, litter, f, totalroots, f_senescence, b, litterfac, T, maxLAI, minLAI, minBiomass, oldBiomass;
  Matrix rd, exudates, deadroots;
  Matrix result(2);

  // initalize LAI at the first time step
  maxLAI = Phenology(4);
  minLAI = (1.0 - Phenology(6)) * maxLAI;
  DoHarvest();                                                    // harvest
  DoGraze();                                                      // grazing
  switch (ProductionModel)                                       // primary production modelled internally
  {
        case 0: PrimProd = SimpleProd(); break;
        case 1: PrimProd = TemperatureProd(); break; // soil temperature dependent prduction
            /* case 2: PrimProd = TemperatureThornley(); break; OBSOLETE */ 
        case 2: PrimProd = NPPData(StepNr);  // primary production from imported data
        case 3: PrimProd = RadProd(); break; // Haxeltine and Prentice model, PAR data supplied
        case 4: PrimProd = RadProd(); break; // Haxeltine and Prentice model, PAR calculated from cloud cover
        case 5: PrimProd = TundraProd(); break; // Shaver photosynthesis model for tundra, PAR data supplied
        case 6: PrimProd = TundraProd(); break; // Shaver photosynthesis model, PAR calculated from cloud cover
    }
    
	if (ProfileOutput.Contains(8))
	{
		result(2) = PrimProd;
        result(1) = Timer + 0.5 * Timestep;
        result.Write(output8);
	}
    if ((SatCorr > 0) && (Saturation(1) < SatCorr)) PrimProd = (Saturation(1) / SatCorr) * PrimProd;
  // correction factor for poor aeration due to topsoil waterlogging
  TotalPrimProd += PrimProd;
  Shoots = ShootsFactor * PrimProd;                              // shoots production
  gw = GwData(StepNr);                                           // partition roots according to root distribution function
  if (NoRootsBelowGWT)                                           // no roots below groundwater table flag
  {
    rd.Resize(NrLayers);
    for (i = 1; i <= NrLayers; i++)                              // find which part of the root distribution function is above the groundwater table
    {
      if (Layers(i, 1) > gw)
      {
        abovegwt += RootDistrib(i);
        toproots = i;
      }
      belowgwt = 1 - abovegwt;
    }                                                             // reshape root distribution function, to add roots below gwt to the roots above gwt
    if (toproots > 0)
    {
      for (i = 1; i <= toproots; i++) rd(i) = RootDistrib(i) + RootDistrib(i) * belowgwt / abovegwt;
    } else rd(1) = 1.0;                                           // water table above the profile, all roots added to top of profile
  } else rd = RootDistrib;
  rd *= (1 - ShootsFactor) * PrimProd;                            // calculate the amount of roots added
  SpringFactor = 1 + SpringCorrection * sin(2 * PI * (DayNr + 284) / 365);  // spring factor to account for enhanced exudate production during active growing season
  exudates = rd * (SpringFactor * ExudateFactor);                 // exudates
  // exudates.Disp();
  deadroots = RootMass * RootSenescence;                          // organic material dying roots;
  for (i = 1; i <= NrLayers; i++)
  {
    rootsadded = rd(i) - exudates(i) - deadroots(i);              // net root addition, prevent negative rootmass
    RootMass(i) = RootMass(i) + rootsadded;
    if (RootMass(i) < 0) RootMass(i) = 0;
    NewSOM(i, 4) = NewSOM(i, 4) + exudates(i);                    // add exudates and dead roots to SOM reservoirs
    NewSOM(i, 5) = NewSOM(i, 5) + deadroots(i);
  }
  if (ProfileOutput.Contains(4)) RootMass.Write(output4);         // log root mass to output file
  /*
   * cout << "root distribution, exudates, root mass:\n";
  RootDistrib.Disp();
  rd.Disp();
  exudates.Disp();
  RootMass.Disp();
  */
  BioMass += Shoots;                                              // total above ground biomass
  BioMassRec(StepNr, 8) = 0.0;
// Biomass senescence and litter production
  minBiomass = minLAI * LAICarbonFraction;;
  if (ProductionModel < 3) {  // Litter production by dying off of above-ground biomass
      litter = BioMassSenescence * BioMass * Timestep; // a distinction between production models is maintianed for keeping compatibility with earlier versions of the model
      LitterLayer += litter;
      BioMass -= litter;
      if (BioMass < 0.0) BioMass = 0.0;
      CurrentLAI = BioMass / LAICarbonFraction;
  } else { // fraction of biomass shedded in autumn is dependent on decrease of LAI per time step
      if (LeafSenescence) { // autumn biomass senescence by temperare-driven decrease/increase of LAI;
          f_senescence = (PreviousLAI - CurrentLAI) / PreviousLAI;
      } else { //  normal biomass senescence during the growing season
          f_senescence = BioMassSenescence;
          CurrentLAI = CurrentLAI * (1 - f_senescence);
      }
      //cout << DayOfTheYear << ' ' << f_senescence << endl;
      // CurrentLAI = CurrentLAI - f_senescence * BioMassSenescence;
      if (CurrentLAI < minLAI) CurrentLAI = minLAI;
      // surviving = BioMass * (1 - Phenology(6));  // surviving biomass
      oldBiomass = BioMass;
      BioMass = (1 - f_senescence) * BioMass;
      if (BioMass < minBiomass) BioMass = minBiomass;
      if (oldBiomass > BioMass) litter = oldBiomass - BioMass; else litter = 0.0;
      LitterLayer += litter;
  }
  // convert aboveground litter to belowground litter reservoir
    T = SoilTemp(1);
    if (T > 0.0) { // no conversion if the soil is frozen
        litterfac = T * (LitterConversion / T_ref);
        litter = litterfac * LitterLayer;
        NewSOM(1, 5) = NewSOM(1, 5) + litter;
        LitterLayer -= litter;
    }
/* Plant respiration:
   Respiration is linearly dependent on both primary production and total biomass, converted to CO2 eq
   total biomass has to be multiplied by the time step since the unit of the conversion factor is day-1 */
  f = Timestep * 3.6641;                                           // timestep and C - CO2 conversion factor
  totalroots = RootMass.Sum();
  // if (ProductionModel < 3) PlantRespiration = f * (RespFac(1) * PrimProd / Timestep + RespFac(2) * (BioMass + totalroots)); else PlantRespiration = PlantRespiration + f * (RespFac(2) * (BioMass + totalroots));
  if (ProductionModel < 3) PlantRespiration = f * (RespFac(1) * PrimProd / Timestep + RespFac(2) * (BioMass + totalroots));
  // For production model 3 and 4, the plant respiration is only calculated at leaf level, the root respiration has to be added 
  if ((ProductionModel == 3) || (ProductionModel == 4)) PlantRespiration = PlantRespiration + f * (RespFac(2) * totalroots);

/* plant respiration for production model 3 and 4 is calculated as part of the photosynthesis module
 however, this is only the growth respiration, the maintenance respiration still should be included
 units coming from photosynthesis model:  kg CO2 m2 per timestep */
  BioMassRec(StepNr, 1) = DayNr;
  BioMassRec(StepNr, 2) = BioMass + totalroots;                     // log total biomass, primary production and respiration
  BioMassRec(StepNr, 3) = PrimProd;
  BioMassRec(StepNr, 4) = PlantRespiration;
  BioMassRec(StepNr, 7) = LitterLayer;
  BioMassRec(StepNr, 8) = HarvestGrazing;
  HarvestGrazing = 0.0;
  BioMassRec(StepNr, 9) = TotalManure;
  TotalManure = 0.0;
  BioMassRec(StepNr, 10) = CurrentLAI;
}

void DoHarvest()
/* Harvest of biomass at selected dates */
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! add seperate output in Biomass to keep track of grazing and harvest!
{
  int i, l;
  double harvested;

// Harvest matrix: harvest dates (1st column) and fraction of biomass harvested (2nd column)
  l = Harvest.Rows();
  for (i = 1; i <= l; i++)     // if the current day of the year falls within the current time step, cut the grass
  {
    if ((Harvest(i, 1) >= (DayOfTheYear - 0.5 * Timestep)) && (Harvest(i, 1) < (DayOfTheYear + 0.5* Timestep))) {
        harvested = BioMass * Harvest(i, 2);
        BioMass -= harvested;
        HarvestGrazing += harvested;
        CurrentLAI *= (1 - Harvest(i, 2));
        Harvested = TRUE;
    } else Harvested = FALSE;
  }
}

void DoGraze()
/* grazing of biomass */
{
  int i, l;
  double grazed;

/* Grazing matrix:
each row represents a range of days in which grazing occurs (column 1: starting day, 2: ending day
followed by the amount of biomass removed (kg C m2/day and the amount of excretion (kg C m2/day) */

  l = Grazing.Rows();
  for (i = 1; i <= l; i++)     // eat the grass if the current day of the year is within one of the grazing periods
  {
    if ((DayOfTheYear >= Grazing(i, 1)) && (DayOfTheYear < Grazing(i, 2))) {
        grazed = BioMass * Grazing(i, 3);
        BioMass -= grazed;
        HarvestGrazing += grazed;
        CurrentLAI *= (1 - Grazing(i, 3));
        Harvested = TRUE;
    } else Harvested = FALSE;  
  }
}



double SimpleProd()
/* Simple primary production from yearly sinusoidal function */
{
  double a;
  a = Timestep *((MaxProd-MinProd) * (0.5 + 0.5 * sin(2 * PI * (DayNr + 284) / YEAR)) + MinProd);
  return a;
}

double TemperatureProd()
/* Primary production from temperature upper soil layer */
{
  double a, maxT, T, tfac;

  maxT = ProdTFunc(2) + (ProdTFunc(2) - ProdTFunc(1)); // temperature range primary production
  T = SoilTemp(1);
  if ((T >= ProdTFunc(1)) && (T <= maxT))              // sinusoidal approach of temperature optimum function
  {
    tfac = 0.5 * (sin (PI * (T - ProdTFunc(1)) / (ProdTFunc(2) - ProdTFunc(1)) - 0.5 * PI) + 1);
  } else tfac = 0;
  a = Timestep * (MinProd + tfac * (MaxProd - MinProd));

  return a;
}

double TemperatureThornley()
/* Primary production from temperature upper soil layer, temperature function cf Thornley (1998) with steepness parameter = 2*/
{
  double a, maxf, T, tfac, Tmin, Tmax, refT, q, m;

  Tmin = ProdTFunc(1);				       // minimum temperature from ProdTFunc(1); ProdTFunc(2) is ignored, max T set to 45
  Tmax = 45.0;
  q = 2; 					       // steepness parameter, provisionally fixed
  refT = 20;					       // fixed reference temperature
  m = (refT-Tmin) * (Tmax-refT);
  maxf = ((pow((Tmax - ((q+2)/q)*Tmin) / ((q+1)/q), q)) * ((Tmax+Tmin) / (q+1))) / m;	// maximum of temperature function
  T = SoilTemp(1);
  if ((T >= ProdTFunc(1)) && (T <= Tmax))              // sinusoidal approach of temperature optimum function
  {
    tfac = ((pow(T - Tmin,q))*(Tmax - T)) / (m * maxf);
  } else tfac = 0;
  a = Timestep * (MinProd + tfac * (MaxProd - MinProd));

  return a;
}


double PARcalc()
/* calculates photosynthetically active radiation or returns data from file
 units: joule per square meter per day
 based on Haxeltine & Prentice 1996 with some erros corrected*/
{
    double par = 0.0; // photosynthetically active radiation in joule per square meter per day
    double shortwave = 0.0; // total daily shortwave radiation joule per square cm
    double swfrac = 0.44; // fraction of total solar radiation that is visible light
    double parfrac = 0.5; // fraction of visible light that is par
    // double estar = 0.27; // conversion factor megajoule to mol
    double aa; // solar declination
    double c = 0.45, d = 0.9, ni = 1.0; // coeff eq A3 Haxeltine and Prentice
    double beta = 0.17; // albedo eq A3 Haxeltine and Prentice
    double lrad, u, v, h = 0.0, zz; //latitude and parameters for day length h
    double Q0, Rs; // total, shortwave radiation
    double cloudcover; // cloudcover fraction derived from weather station octants
    double B  = 0.7; // Factor depressing rediation with cloudcover cf Budyko

    // daylenght calculation
    aa = -23.4 * cos (2 * PI * (DayOfTheYear + 10) / YEAR );  // sun declination
    lrad = DEG2RAD * Latitude;  // calculation of daylength
    u = sin(lrad) * sin(DEG2RAD *aa); // eq A6, A7 Haxeltine & Prentice
    v = cos(lrad) * cos(DEG2RAD *aa); 
    if (u <= v) h = 0;
    if ((u > -v) & (u < v)) h = 24.0 * acos(-u / v) / PI;
    if (u >= v) h = 24;
    DayLength = h;
    
    /* Eq A9 in Haxeltine & Prentice is h = 24 * acos(-u / v) / (2 * pi);
     this one results in negative day lengths!
     is incorrect cf Bonan
     Formula Prentice 1993:
     cos z = u + v cos h met h time of the day in angular units
     Derivation correct equation:
     instantaneous radiation is (solar constant / radius vector)* cos z (Bonan)
     solar noon; integration over day
     so integration over 0.5h*2pi/24 to -05h2pi/2  (Bonan: h negative after 12:00)
     u + v cos h)dh -> uh + v*sin(h) ->
     uh + v sin(h)  -0.5h*2pi/24 to 0.5h2pi/24 
     u(0.5h2pi/24) + v sin( 0.5h2pi/24) - u(-0.5h2pi/24) - v sin(- 0.5h2pi/24)
     u(0.5h2pi/24) + v sin( 0.5h2pi/24) + u(0.5h2pi/24) + v sin(0.5h2pi/24)
     2u(0.5h2pi/24)+2v sin( 0.5h2pi/24)
     u(h2pi/24)+ 2v sin( hpi/24) 
     Which differs from Haxeltine and Prentice*/
    if ((ProductionModel == 3) || (ProductionModel == 5))
    { // data from file, convert to joule per square meter from joule per cm2 which is usually recorded by weather stations
        // shortwave = swfrac * parfrac * PARData(StepNr) * 1.0e4;
        shortwave = swfrac * PARData(StepNr) * 1.0e4;
    } else { // calculate shortwave radiation
        cloudcover = PARData(StepNr);
        zz = u * h * 2 * PI / 24.0 + 2 * v * sin(h * PI /24.0); // correct integral 
        Q0 = 3600.0 * SOLARCONSTANT *(1 + 2 * 0.01675 * cos(DEG2RAD * (360.0 * DayOfTheYear / YEAR))); // incoming radiation based on solar constantant and variations in earth's orbit
        Rs = zz * (c + d * ni) * (1 - beta) * Q0;
        shortwave = Rs * (1 - B * cloudcover);  // equation 1a in Hurley and Boers 1996
    }
    // par = parfrac * estar * shortwave / 1.0e6;
    par = parfrac * shortwave;
    return par;
}

double PhotoSynthesis(double LAI, double I)
/* Photosynthesis cf Haxeltine and Prentice
 LAI is leaf area index calculated from phenology
 I is photosynthetic active radiation Joule per square meter per day*/
{
    // lpj photosynthesis
    const double cq = 4.6e-6;   //conversion factor for solar radiation at 550 nm  from J/m2 to E/m2 (E mol quanta)
    const double tau25 = 2600.0; // undefined parameter in eq. 19 Sitch et al; used for CO2 compensation point, haxeltine and Prentice call it the CO2/O2 specifity ratio
    double tau = 0.0;
    double gammastar = 0.0;     // CO2 compensation point,
    const double q10tau = 0.57;  // q10 for temperature-sensitive parameter tau 
    double temp_co2[2] = {-4.0, 45.0}; // lower and upper temperature limit for co2 (deg C)
    double temp_photos[2] = {10.0, 30.0};  // lower and upper limit of temperature optimum for photosynthesis(deg C)
    double alpha[2] = {0.08, 0.053}; // effective ecosystem quantum efficiency CO2 uptake first entry C3 plants second C4
    double tm[2] = {45.0, 55.0};  // maximum temperature for photosynthesis, first entry for C3,seconfnd for C4
    double lambda = 0.7;            // parameter relating internal and external pCO2
    double lambdam[2] = {0.8, 0.4}; // optimal ratio of intercellular to ambient CO2
    double b[2] = {0.015, 0.02};    // leaf respiration as fraction of Vmax 1st value for C3, 2nd for C4 plants
    const double alphaa = 0.5;      // fraction of PAR assimilated at ecosystem level relative to leaf level
    const double theta = 0.7;       // shape parameter of co-limitation by light and Rubisco activity
    const double q10ko = 1.2;       // q10 for temperature-sensitive parameter ko
    const double q10kc = 2.1;       // q10 for temperature-sensitive parameter kc
    const double ko25 = 3.0e4;      // value of ko at 25 deg C inhibition constant of O2
    const double kc25 = 30.0;       // value of kc at 25 deg C Michaelis Menten constant of CO2
    const double cmass = 12.0;      // atomic mass of carbon
    double vm = 0.0;                // Vmax
    double p_int = 0.0;             // internal leaf CO2 partial pressure
    double fpar = 0;                // FPAR (PAR at leaf level) cf Haxeltine and Prentice 1996 eq 1
    double apar = 0;                // absorbed par fraction
    double je = 0.0;                // je is PAR-limited photosynthesis rate molC/m2/h
    double jc = 0.0;                // jc is rubisco-activity-limited photosynthesis rate JC, molC/m2/h
    double agd = 0.0;               // daily gross photosynthesis, Agd, gC/m2/day
    double rd = 0.0;                // Daily leaf respiration, Rd, gC/m2/day 
    double nd = 0.0;                // Daily net photosynthesis (at leaf level), And, gC/m2/day
    double adt = 0.0;               // Total daytime net photosynthesis, Adt, gC/m2/day
    double temp, h, pO2, pa, tstress, k1, k2, k3, low, high, s, c1, c2, ko, kc, fac, sigma, phipi;
    int ptype;
    
    b[0] = RespFac(1);  // The first value of respfac is the leaf respiration factor depending on C3 or C4 photosynthesis
    b[1] = RespFac(1);  // The second value of Respfac is the root and stem respiration and is added outside this function
    temp = TData(StepNr);
    if ((LAI == 0.0) || (I == 0.0) || (temp < 0.0)) {PlantRespiration = 0.0; return (0.0);} // no use calculating when there is no PAR, LAI or low temperatures below the minimum requirement!
    ptype  = (int)Phenology(5);      // photosynthesis type, 1 for C3, 2 for C4
    if (ptype > 2) {
        ptype = 1;
        cout << "Photosynthesis: plant type parameter (C3 or C4) not defined correctly, has been set to C3" << endl;
    }
    h = DayLength;                  // daylength in hours 
    pO2 = 101325 * 0.20946;         // ambient O2 partial pressure [Pa]
    if (strlen(CO2File) != 0) {
        AmbientCO2 = CO2Data((int)(Year - StartYear + 1));
    }
    pa = 101325 * AmbientCO2 / 1.0e6; // ambient partial pressure CO2
    // ftemp (tstress) PFT-specific temperature inhibition function 
    // no proper definition in Sitch et al, taken from LPJ source
    if ((h < 0.01) || (temp > tm[ptype - 1])) tstress = 0.0; else { // short daylength or temperatures exceeding maximum
        k1= 2 * log(1/0.99-1) / (temp_co2[0] - temp_photos[0]);
        k2= (temp_co2[0] + temp_photos[1]) * 0.5;
        low = 1 / (1 + exp(k1 * (k2 - temp)));
        k3 = log(0.99/0.01) / (temp_co2[1] - temp_photos[1]);
        high = 1 - 0.01 * exp(k3 * (temp - temp_photos[1]));
        tstress = low * high;
    }
    fpar = 1 - exp(-KBeer * LAI);  // FPAR cf Haxeltine and Prentice 1996 eq 1
    apar = fpar * I * alphaa;       // absorbed par cf Haxeltine and Prentice 1996
    s = (24 / h) * b[ptype - 1];   // eq 16 beware, Sitch uses a variale 'a', being the same as b here!
    if (ptype == 1) {
        tau = tau25 * exp((log(q10tau)) * (temp-25) * 0.1); // tau - CO2 compensation point
        gammastar = pO2 / (2 * tau); // CO2 compensation point
        p_int = lambdam[0] * pa; // internal CO2 partial pressure
        c1 = alpha[0] * tstress * (p_int - gammastar) / (p_int + 2.0 * gammastar); // eq 17
        ko = ko25 * exp((log(q10ko)) * (temp-25) * 0.1);
        kc = kc25 * exp((log(q10kc)) * (temp-25) * 0.1);
        fac = kc * (1 + pO2 / ko);
        c2 = (p_int - gammastar) / (p_int + fac); // eq 18 
        sigma = 1 - (c2 - s)/(c2 - theta * s); // eq 15
        if (sigma > 0) sigma = sqrt(sigma); else sigma = 0.0;
        vm = (1.0 / b[0]) * (c1 / c2) * ((2.0 * theta - 1.0) * s - (2.0 * theta * s - c2) * sigma) * apar * cmass * cq;
        p_int = lambda * pa; // recalculation c1, c2 cf source code LPJ with actual p_int
        c1 = alpha[0] * tstress * (p_int - gammastar) / (p_int + 2.0 * gammastar);
        c2 = (p_int - gammastar) / (p_int + fac);
    } else {
        c1=tstress * alpha[1];
        c2 = 1.0;
        sigma = 1 - (c2 - s)/(c2 - theta * s); // eq 15
        if (sigma > 0) sigma = sqrt(sigma); else sigma = 0.0;
        vm  =(1.0 / b[1]) * c1 / c2 * ((2.0 * theta - 1.0) * s - (2.0 * theta * s - c2) * sigma) * apar * cmass * cq;
        // Parameter accounting for effect of reduced intercellular CO2 concentration on photosynthesis, Phipi. Eqn 14,16, Haxeltine & Prentice 1996
        phipi = lambda / lambdam[1];
        if(phipi < 1.0) c1 = tstress * phipi * alpha[1];
    }
    // je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    // Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h Eqn 3, Haxeltine & Prentice 1996
    je = c1 * apar * cmass * cq / h;
    // Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h Eqn 5, Haxeltine & Prentice 1996
    jc = c2 * vm / 24.0;
    // Calculation of daily gross photosynthesis, Agd, gC/m2/day Eqn 2, Haxeltine & Prentice 1996
    agd = h * (je + jc - sqrt(pow((je + jc), 2.0) - 4.0 * theta * je * jc)) / (2.0 * theta);
    // Daily leaf respiration, Rd, gC/m2/day Eqn 10, Haxeltine & Prentice 1996
    rd = b[ptype - 1] * vm;
    // Daily net photosynthesis (at leaf level), And, gC/m2/day
    nd = agd - rd;
    // Total daytime net photosynthesis, Adt, gC/m2/day Eqn 19, Haxeltine & Prentice 1996
    adt= nd + (1.0 - h / 24) * rd;
    
    PlantRespiration = Timestep * C_CO2 * rd / 1000.0;  // convert from g C/m2/d to kg CO2 m2 per timestep
    return (Timestep * adt / 1000.0);  // convert from g C/m2/d to kg C m2 per timestep
    
}


double LAICalc()
/* Calculates LAI (leaf area index) from photosynthesis 
 this includes not only Photosynthesis models 5 and 6 but also other models 
 so the Phenology parameter should also be defined for the other models*/
{
    double gdd = 0.0, LAI = 0.0, t, p, kgCadded, maxLAI, minLAI, LAIadded, maxGDD, seasonstart, autumnstart;
    BOOLEAN autumn = FALSE; // indication of spring/summer or atumn season
    
    maxLAI = Phenology(4);
    minLAI = (1.0 - Phenology(6)) * maxLAI;
    if (StepNr == 1) {
        LAI = minLAI;
        CurrentLAI = LAI;
        PreviousLAI = CurrentLAI;
    }
    if (Phenology(1) == 1.0) // summergreen phenology depending on heat sum
    {
        maxGDD = Phenology(3);
        seasonstart = Phenology(7);
        //seasonstart = 2.0;
        autumnstart = Phenology(8);
        if (StepNr < seasonstart) LAI = minLAI; // the first 30 days in the first year of the simulation get the minimumLAI, for the next years it depends on the BiomassSenenescence and weather 
        if (DayOfTheYear != (seasonstart + 0.5 * Timestep))
        {   
            gdd = Timestep * (TData(StepNr) - Phenology(2));  // heat sum
            if (DayOfTheYear == (autumnstart + 0.5 * Timestep)) GrowingDegreeDays = maxGDD;
            if ((DayOfTheYear > autumnstart) || (DayOfTheYear < seasonstart)) autumn = TRUE; // autumn; in autumn GDD can decrease
            // if heat sum positive then calculculate LAI
            if (autumn) { // autumn
                if (gdd < 0.0) {  
                    LeafSenescence = TRUE; 
                    GrowingDegreeDays += gdd; // in autumn GDD may decrease
                    if (GrowingDegreeDays < 0.0) GrowingDegreeDays = 0.0;
                } else LeafSenescence = FALSE;
            } else if (gdd > 0.0) GrowingDegreeDays += gdd;
            if (GrowingDegreeDays > Phenology(3)) PotentialLAI = maxLAI; else PotentialLAI = minLAI + (maxLAI - minLAI) * GrowingDegreeDays / Phenology(3);
            // potential LAI can decrease in autumn, so does LAI
            if (CurrentLAI < PotentialLAI) {   // Harvest or grazing has occurred
                kgCadded = PrimProd * ShootsFactor;
                LAIadded = kgCadded / LAICarbonFraction;
                LAI = CurrentLAI + LAIadded;
                if (LAI > PotentialLAI) LAI = PotentialLAI;
            } else LAI = PotentialLAI;
                //maxLAI = Phenology(4) * GrowingDegreeDays / Phenology(3);
                //kgCadded = PrimProd * ShootsFactor; // primary production based LAI growth to account for possible harvest or grazing
                //LAI = CurrentLAI + kgCadded / LAICarbonFraction;
                //if (LAI > maxLAI) LAI = maxLAI;
        } else {
            LeafSenescence = FALSE;
            GrowingDegreeDays = 0.0;
            PotentialLAI = minLAI;
            CurrentLAI = PotentialLAI;
            LAI = minLAI;
        }
    } else LAI = Phenology(4); // evergreen phenology with constant LAI AND NO HARVEST
    PreviousLAI = CurrentLAI;
    CurrentLAI = LAI;
    //cout << LAI << endl;
    return LAI; 
}

double RadProd()
/* Primary production from photosynthetically active radiation, Haxeltine and Prentice */

{
    double par, NPP, LAI;
    
    /* Calculation of Growing Degree Days; counter is reset at start of the year */
    
    par = PARcalc(); // PAR calculation
    LAI = LAICalc();
//    CurrentLAI = LAI;
    NPP = PhotoSynthesis(LAI, par);
    return NPP;
}

double TundraProd()
/* Primary production from photosynthetically active radiation, for tundra, Shaver et al, J. Ecology 2007 */
{
    double par, NPP, LAI, gpp, er, temp, I, convfac, resp0, respbeta, pmax, pslope;
    /* To be added to parameter files:
     resp0 Plant respiration at zero degrees
     respbeta Temperature sensitivity factor plant respiration
     pmax light-saturated photosynthetic rate per unit leaf area (μmol m–2 leaf s–1)
     pslope is the initial slope of the light response curve (μmol CO2 μmol–1 photons)
     */
    
    /* Calculation of Growing Degree Days; counter is reset at start of the year */

    par = PARcalc(); // PAR calculation
    LAI = LAICalc();
//    CurrentLAI = LAI;
    temp = TData(StepNr);
    resp0 = PhotoPar(1);
    respbeta = PhotoPar(2);
    pmax = PhotoPar(3);
    pslope = PhotoPar(4);
    //cout << par << ' ' << temp << ' ' << LAI << ' ' << endl;
    if ((LAI == 0.0) || (par == 0.0) || (temp < 0.0)) {PlantRespiration = 0.0; return (0.0);} // no use calculating when there is no PAR, LAI or low temperatures below the minimum requirement!
    /* GPP calculated cf Shaver from LAI and incoming radiation
     Plant respiration is their ER1 respiration model */
    er = resp0 * LAI * exp(respbeta * temp); //Plant resp μmol m–2 ground s–1)
    I = (par / (DayLength * 3600)) * 4.6035;
    /* convert from joule per square meter per day to μmol photons m–2 s–1
     first from joule m2 day-1 to W m-2 = J m2 s-1 by dividing by numbers of seconds of daylight
     multiplication factor is derived from https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux
     assuming that 66% is radiation between 470 and 520 nm and the rest is 66 nm */
    gpp = (pmax / KBeer) * log((pmax + pslope * I) / (pmax + pslope * I * exp(-KBeer * LAI))); // Gross Primary Production μmol m–2 ground s–1
    //cout << er << ' ' << gpp << endl;
    // Conversions from micromols CO2 m-2 s-1 to kg CO2 m-2 per timestep
    // convfac = Timestep * DayLength * 3600 * 44.01 / 1.0e9;
    // BUG CORRECTION: convert to kg C m-2 per timestep 
    convfac = Timestep * 24 * 3600 * MOLWEIGHTCO2 * 1.0e-9; // convert from μmol CO2 m–2 ground s–1 to kg CO2 m2 h-1; plant respiration occurs during the whole day
    PlantRespiration = er * convfac;
    convfac = Timestep * DayLength * 3600 * MOLWEIGHTC * 1.0e-9; // convert from μmol CO2 m–2 ground s–1 to kg C m2 h-1;photosynthesis only during the day
    NPP = (gpp - er) * convfac;
    if (NPP < 0.0) NPP = 0.0;
    
    //convfac = Timestep * DayLength * 3600 * 44.01 / (1.0e9 * C_CO2); // klopt deze conversiefactor wel?  
    //PlantRespiration = er * convfac;
    //gpp = gpp * convfac;
    //NPP = gpp - PlantRespiration;
    return NPP;
}

void AddManure()                            // manure addition
// !!!!!!!!!!!!!!!!!!!!!!! This should be added to output to allow carbon balance calculation!!!!!!!!!!!!!!!!!!!!!
{
  int i, j;
  double day;

// Manure addition
  for (i = 1; i <= Manure.Rows(); i++)       // check if manure has to be added at the current day
  {
    day = Manure(i, 1);
    if ((day >= (DayOfTheYear - 0.5 * Timestep)) && (day < (DayOfTheYear + 0.5* Timestep)))       // the day should fall within the time step centred on the current day
    {
      for (j = 1; j <= ManureLayers.Rows(); j++)
      {
        NewSOM(j, 2) = NewSOM(j, 2) + Manure(i, 2) * ManureFluidFrac * ManureLayers(j, 1);        // manure fluids
        NewSOM(j, 3) = NewSOM(j, 3) + Manure(i, 2) * (1 - ManureFluidFrac) * ManureLayers(j, 2);
        TotalManure += Manure(i,2);
// manure solids may be added to deeper layers if manure injection is used
      }
    }
  }
// Cattle excretion addition
  for (i = 1; i <= Grazing.Rows(); i++)     // Check if the current day of the year is within one of the grazing periods
  {
    if ((DayOfTheYear >= Grazing(i, 1)) && (DayOfTheYear < Grazing(i, 2))) {
        NewSOM(1, 3) = NewSOM(1, 3) + Grazing(i, 4);
        TotalManure += Grazing(i,4);
    }
// excretion is assumed to be mostly solid and added to the top layer only
  }
}

void CollectBioMass()
// !!!!!!!!!!!!!!!!!!!! add LAI, manure addition and harvest
/* collects total Biomass, primary production, respiration, net CO2 flux incl. soil respiration
   in BioMassRec
   Elements of each row:
   1. Day number
   2. Biomass
   3. Primary production
   4. Plant respiration
   5. Net CO2 flux incl. soil respiration
   6. Soil repiration + dark respiration vegetation
   7. Litter mass
   8. Biomass removed by harvest and grazing
   9. LAI
*/
{
  int i;
  double f, s;
  Matrix lt;

/* Units in which the biomass and respiration data is stored originally
 BioMassRec(StepNr, 1) = DayNr;
 BioMassRec(StepNr, 2) = BioMass + totalroots;      dit is kg C/m2 
 BioMassRec(StepNr, 3) = PrimProd;  dit is kg C/m2/timestep
 BioMassRec(StepNr, 4) = PlantRespiration; kg CO2 m2 per timestep
 BioMassRec(StepNr, 5) = Soil respiration + plant respiration - NPP = NEE
 BiomassRec(StepNr, 6) = Soil respiration + plantrespiration (ecosystem respiration)
 BioMassRec(StepNr, 7) = LitterLayer;
*/
    
    
  f = C_CO2*1.0e6/(24*Timestep);     // conversion factor from kg C/m-2/timestep into mg CO2 m-2/hr
  lt.Resize(NrLayers);
  for (i = 1; i <=NrOfSteps; i++)
  {
    BioMassRec(i, 3) = f * BioMassRec(i, 3);
    BioMassRec(i ,4) = BioMassRec(i, 4) * 1.0e6 / (24 * Timestep);  // plant respiration: transform from kg CO2/m2/timestep into mg CO2 /m2/hr
    lt.PutData(1, 1, LayerTime, i, 2, NrLayers);
    s = lt.Sum();
    BioMassRec(i, 5)  = s + BioMassRec(i, 4) - BioMassRec(i, 3);  // Net CO2 production/uptake incl. soil respiration
    BioMassRec(i, 6)  = s + BioMassRec(i, 4); // soil + plant respiration for comparing with dark chamber / nighttime measurements
  }
}
