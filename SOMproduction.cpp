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
    double gw, belowgwt = 0, abovegwt = 0, rootsadded, litter, f, totalroots, f_senescence, b, litterfac, T, maxLAI, minLAI, minBiomass, oldBiomass, oldrootmass, oldshoots, plantrespC, rootresp, overshoot;
    BOOLEAN autumn;

  Matrix rd, exudates, deadroots;
  Matrix result(2);

  // initalize LAI at the first time step
  maxLAI = Phenology(4);
  minLAI = (1.0 - Phenology(6)) * maxLAI;
  minBiomass = minLAI * LAICarbonFraction;
  OldLitter = LitterLayer;   // store existing litter C for calculation of storage change in function CollectCO2
  oldrootmass = RootMass.Sum(); // store old rootmass and old biomass
  oldshoots = BioMass;
  TotalManure = 0.0;

  //if (HarvestModel == 2) DateHarvest(minBiomass); else DoHarvest(minBiomass);      // harvest
  DoHarvest();
  DoGraze();                                                      // grazing
  AddManure();

  switch (ProductionModel)                                       // primary production modelled internally
  {
        case 0: NPP = SimpleProd(); break;
        case 1: NPP = TemperatureProd(); break; // soil temperature dependent prduction
        case 2: NPP = NPPData(StepNr);  // primary production from imported data
        case 3: NPP = RadProd(); break; // Haxeltine and Prentice model, PAR calculated from total daily radiation.
        case 4: NPP = TundraProd(); break; // Shaver photosynthesis model for tundra, PAR data supplied
  }

  /* Plant respiration:
   Respiration is linearly dependent on both primary production and total biomass, calculated in photoysnthesis models as kg CO2/m2
   total biomass has to be multiplied by the time step since the unit of the conversion factor is day-1 
   For production model 3 and 4, only the leaf respiration is calculated in the model, the root respiration
   has to be added and NPP corrected for it
   For productionmodel 5 and 6, total plantrespiration is included in the model*/
  f = Timestep * 3.6641;                                           // timestep and C - CO2 conversion factor
  if (ProductionModel < 3) {
      totalroots = RootMass.Sum();  // total root mass kg C
      plantrespC = RespFac(1) * NPP / Timestep + RespFac(2) * (BioMass + totalroots);
      NPP -= plantrespC;
      if (NPP < 0.0) NPP = 0.0;
      PlantRespiration = f * plantrespC;
  }
    /* plant respiration for production model 3 and 4 is calculated as part of the photosynthesis module
     however, this differs between model 3 and 4; in 3, no root respiration is included
     units coming from photosynthesis model:  kg CO2 m2 per timestep */
  if (ProductionModel == 3) { // Model 3 has only plant respiration of above-ground biomass
      totalroots = RootMass.Sum();  // total root mass kg C
      rootresp = RespFac(2) * totalroots;
      NPP -= rootresp;
      if (NPP < 0.0) NPP = 0.0;
      PlantRespiration = PlantRespiration + f * rootresp;
  }
  CarbonBalance(StepNr, 15) = PlantRespiration * CONVKGCO2TOMOLC;
  if (ProfileOutput.Contains(8))
  {
    result(2) = NPP;
    result(1) = Timer + 0.5 * Timestep;
    result.Write(output8);
  }
  CarbonBalance(StepNr, 1) = (NPP + PlantRespiration) * CONVKGCTOMOLC;
  TotalNPP += NPP;
  Shoots = ShootsFactor * NPP;                              // shoots production
  overshoot = Shoots * LAIovershoot; // correction for overshooting max LAI by allocating more carbon to beloe-ground biomass
  Shoots = Shoots - overshoot; // decrease allocation to shoots if overshoot occurred
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
    }   // reshape root distribution function, to add roots below gwt to the roots above gwt
    if (toproots > 0)
    {
      for (i = 1; i <= toproots; i++) rd(i) = RootDistrib(i) + RootDistrib(i) * belowgwt / abovegwt;
    } else rd(1) = 1.0; // water table above the profile, all roots added to top of profile
  } else rd = RootDistrib;
  rd *= (1.0 - ShootsFactor) * NPP;                            // calculates the amount of roots added per layer
  rd = rd + overshoot; // additional allocation to belowground biomass if production overshoots max LAI
  SpringFactor = 1 + SpringCorrection * sin(2 * PI * (DayNr + 284) / 365);  // spring factor to account for enhanced exudate production during active growing season
  exudates = rd * (SpringFactor * ExudateFactor);                 // exudates, are a fraction of the new root material added
  deadroots = RootMass * RootSenescence;          // organic material dying roots;
  for (i = 1; i <= NrLayers; i++)
  {
    rootsadded = rd(i) - exudates(i) - deadroots(i);              // net root addition, prevent negative rootmass
    RootMass(i) = RootMass(i) + rootsadded;

    if (RootMass(i) < 0) {
        RootMass(i) = 0;
        cout << PRODUCTION_ERROR1 << endl;
    }
    NewSOM(i, 4) = NewSOM(i, 4) + exudates(i);                    // add exudates and dead roots to SOM reservoirs
    NewSOM(i, 5) = NewSOM(i, 5) + deadroots(i);
  }
  if (ProfileOutput.Contains(4)) RootMass.Write(output4);         // log root mass to output file
  totalroots = RootMass.Sum();  // total root mass kg C
/* plant respiration for production model 3 and 4 is calculated as part of the photosynthesis module
 however, this is only the growth respiration, the maintenance respiration still should be included
 units coming from photosynthesis model:  kg CO2 m2 per timestep */
    BioMass += Shoots;               // total above ground biomass kg C
    BioMassRec(StepNr, 8) = 0.0;
// Biomass senescence and litter production
    oldBiomass = BioMass;
    minBiomass = minLAI * LAICarbonFraction;
    if ((DayOfTheYear > Phenology(8)) || (DayOfTheYear < Phenology(7))) autumn = TRUE; else autumn = FALSE;
    if (ProductionModel < 3) {  // Litter production by dying off of above-ground biomass
        litter = BioMassSenescence * BioMass * Timestep; // a distinction between production models is maintained for keeping compatibility with earlier versions of the model
        BioMass -= litter;
        if (BioMass < 0.0) BioMass = 0.0;
        CurrentLAI = BioMass / LAICarbonFraction;
    } else { // fraction of biomass shedded in autumn is dependent on decrease of LAI per time step
        if (HDcount > 21){ // normal biomass senescence only occurs more than 3 weeks after harvest
          if (LeafSenescence) { // autumn biomass senescence by temperature-driven decrease/increase of LAI;
              f_senescence = (PreviousLAI - CurrentLAI) / PreviousLAI;
          } else { //  normal biomass senescence during the growing season
              f_senescence = BioMassSenescence * Timestep;
              CurrentLAI = CurrentLAI * (1 - f_senescence);
          }
      } else f_senescence = 0;
      BioMass = (1 - f_senescence) * BioMass;
      if ((HDcount == 22) && (CurrentLAI < minLAI)) cout << StepNr << " LAI is not recovering to minimum LAI after harvest" << endl;
      if (((HDcount > 21) || ((HDcount > 7) && (autumn == TRUE))) && (CurrentLAI < minLAI)) CurrentLAI = minLAI;
  }

  if (BioMass <= minBiomass){
    BioMass = minBiomass;
    CurrentLAI = minLAI;
  }
  if (oldBiomass > BioMass) litter = oldBiomass - BioMass; else litter = 0.0;
  LitterLayer += litter;
  // convert aboveground litter to belowground litter reservoir
  T = SoilTemp(1);
  if (T > 0.0) { // only conversion if the soil is not frozen; otherwise no biologicaal activy to physically transport litter into the top layer

        litterfac = T * (LitterConversion / T_ref);
        litter = litterfac * LitterLayer;
        NewSOM(1, 5) = NewSOM(1, 5) + litter;
      LitterLayer -= litter;
  }

  BioMassRec(StepNr, 1) = DayNr;
  BioMassRec(StepNr, 2) = BioMass + totalroots;                     // log total biomass, kg C m-2
  BioMassRec(StepNr, 3) = NPP;    //  primary production kg C m-2 timestep-1
  BioMassRec(StepNr, 4) = PlantRespiration; // kg CO2 m-2 timestep-1
  // litter is stored after litter decomposition in CollectCO2
  BioMassRec(StepNr, 8) = HarvestGrazing;
  CarbonBalance(StepNr, 19) = HarvestGrazing  * CONVKGCTOMOLC;
  CarbonBalance(StepNr, 20) = (BioMass - oldshoots) * CONVKGCTOMOLC;
  CarbonBalance(StepNr, 21) = (totalroots - oldrootmass) * CONVKGCTOMOLC;
  HarvestGrazing = 0.0;
  BioMassRec(StepNr, 9) = TotalManure;
  CarbonBalance(StepNr, 2) = TotalManure * CONVKGCTOMOLC;
  BioMassRec(StepNr, 10) = CurrentLAI;
  BioMassRec(StepNr, 11) = GPP;
}


//  Harvest function according to Merit; combines DoHarvest and DateHarvest
void DoHarvest()
// Harvest of biomass at selected dates

{
  int i, l;
  double harvested, potentialharvest, remaining, add2litter, harvestfraction;
  HDcount += Timestep;      // Count of days after harvest. Used to correct GPP first 6 days after harvest in Shaver model TundraProd()
  if (HarvestModel == 2){
    // Harvest 
    // HYY = year (e.g. 2010), HMM = month (e.g. 3), HDD = dayofyear (e.g. 201) 
    // Harv_height = fraction of biomass harvested
    //
    // Model 2 werkt niet, HYY en HDayNr blijven constant!
    // Oorzaak: Error in het inlezen van de harvest dates
    // DayOfTheYear is een integer, er hoeft geen 0.5 af te worden getrokken
      
    if ((HYY == CalendarYear) && (HdayNr == DayOfTheYear)){
      remaining = BioMass * (1.0 - Harv_height); 
      harvested = BioMass * Harv_height; 
      CurrentLAI *= (1.0 - (harvested/BioMass));
      BioMass -= harvested;
      add2litter = harvested * 0.05; // 5% of harvest is added to litter layer
      LitterLayer += add2litter;
      harvested = harvested - add2litter;
      HarvestGrazing += harvested;
      checkHarvestDate();
      HDcount = 0;
      Harvested = TRUE;
    } else Harvested = FALSE;
  }
  else{
    // Harvest matrix: harvest dates (1st column) and fraction of biomass harvested (2nd column)
    l = Harvest.Rows();
    for (i = 1; i <= l; i++)     // if the current day of the year falls within the current time step, cut the grass
    {
      if ((Harvest(i, 1) >= (DayOfTheYear - 0.5 * Timestep)) && (Harvest(i, 1) < (DayOfTheYear + 0.5 * Timestep))) {
        remaining = BioMass * (1.0 - Harvest(i, 2)); 
        harvested = BioMass * Harvest(i, 2);
        CurrentLAI *= (1.0 - (harvested/BioMass));
        BioMass -= harvested;
        add2litter = harvested * 0.05; // 5% of harvested biomass is lost to litter layer
        LitterLayer += add2litter;
        harvested = harvested - add2litter;
        HarvestGrazing += harvested;
        HDcount = 0;
        Harvested = TRUE;
      } else Harvested = FALSE;
    }
  }
  if (Harvested & Verbose) cout << "Harvest occurred" << endl;
}

//
    
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
  a = Timestep *((MaxNPP-MinNPP) * (0.5 + 0.5 * sin(2 * PI * (DayNr + 284) / YEAR)) + MinNPP);
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
  a = Timestep * (MinNPP + tfac * (MaxNPP - MinNPP));

  return a;
}


double PARcalc()
//calculates photosynthetically active radiation or returns data from file
//units: joule per square meter per day
//based on Haxeltine & Prentice 1996 with some erros corrected. !!This should be adapted for new input: PAR in umol m-2 s-1
// for other models than shaver!!
{
    double par = 0.0; // photosynthetically active radiation in joule per square meter per day
    double shortwave = 0.0; // total daily shortwave radiation joule per square cm
    double aa; // solar declination
    double c = 0.45, d = 0.9, ni = 1.0; // coeff eq A3 Haxeltine and Prentice
    double beta = 0.17; // albedo eq A3 Haxeltine and Prentice
    double lrad, u, v, h = 0.0, zz; //latitude and parameters for day length h
    double Q0, Rs; // total, shortwave radiation
    double cloudcover; // cloudcover fraction derived from weather station octants
    const double B  = 0.7; // Factor depressing rediation with cloudcover cf Budyko
    const double twentyfour = 24.0;
    const double RE2PHOTONS = 4.56;
    const double TR2PARCONV = 0.372;              // The unitless fraction of total global radiation that is PAR.
    
    // daylength calculation
    aa = -23.4 * cos (2 * PI * (DayOfTheYear + 10) / YEAR );  // sun declination
    lrad = DEG2RAD * Latitude;  // calculation of daylength
    u = sin(lrad) * sin(DEG2RAD *aa); // eq A6, A7 Haxeltine & Prentice
    v = cos(lrad) * cos(DEG2RAD *aa); 
    if (u <= v) h = 0;
    if ((u > -v) & (u < v)) h = 24.0 * acos(-u / v) / PI;
    if (u >= v) h = 24;
    DayLength = h;

    // ProductionModel 3 function uses units: J m-2 day-1
    // ProductionModel 4 function uses units: umol m-2 s-1
    // weather stations usually report J cm-2 day-1, which is shortwave radiation intensity integrated over the day or other time periods;
    // also instantaneous data in W/m2 may be reported
    // shortwave radiation as measured by the commonly used pyraonmeter is the wavelength range 300nm to 3000nm = UV+Visible+Near Infrared (Kipp website)
    //Dye (2004) conversie factor: 4.56 umol/joule, here 4.6035

    if (RADunits == 0) // input file is PAR radiation in umol m-2 s-1
    {
        par = PARData(StepNr); // no conversion for ProductionModel 4
        if (ProductionModel == 3) {     // convert par from umol m-2 s-1 to J m-2 day-1
            par = par / RE2PHOTONS * (3600 * twentyfour);
        }
    // ervan uitgaande dat de PAR het totaal is over de daglengte
    } else if (RADunits == 1) // input file is total daily radiation in J cm-2 day-1
    {
        //  convert total daily radiation in J cm-2 day-1 to PAR J m-2 day-1
        //par = PARData(StepNr) * 1.0e4 * parfrac;   // no further conversion forProductionModel 3
        par = PARData(StepNr) * 1.0e4 * TR2PARCONV;   // no further conversion forProductionModel 3
        if (ProductionModel == 4) {par = (par / (twentyfour * 3600)) * RE2PHOTONS ;}
        // convert PAR from J m-2 day-1 to umol m-2 s-1
        // Divide by seconds of daylight hr-1 (3600)
        // multiplication factor (4.6035) is derived from https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux
        // assumes that 66% is radiation between 470 and 520 nm
    } else if (RADunits == 2) // input is cloud cover
    {
        cloudcover = PARData(StepNr);
        zz = u * h * 2 * PI / 24.0 + 2 * v * sin(h * PI /24.0); // correct integral 
        Q0 = 3600.0 * SOLARCONSTANT *(1 + 2 * 0.01675 * cos(DEG2RAD * (360.0 * DayOfTheYear / YEAR))); // incoming radiation based on solar constantant and variations in earth's orbit
        Rs = zz * (c + d * ni) * (1 - beta) * Q0;
        par = Rs * (1 - B * cloudcover) * TR2PARCONV;  // equation 1a in Hurley and Boers 1996
        // For ProductionModel 3 (PAR J m-2 day-1) no further conversion required
        if (ProductionModel == 4)
            // convert PAR from J m-2 day-1 to umol m-2 s-1
            {par = (par / (twentyfour * 3600)) *  RE2PHOTONS;}
    } else {
        cout << PRODUCTION_ERROR3 << endl;
        exit(EXIT_FAILURE);
    }
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

    const double theta = 0.7;       // shape parameter of co-limitation by light and Rubisco activity
    const double q10ko = 1.2;       // q10 for temperature-sensitive parameter ko
    const double q10kc = 2.1;       // q10 for temperature-sensitive parameter kc
    const double ko25 = 3.0e4;      // value of ko at 25 deg C inhibition constant of O2
    const double kc25 = 30.0;       // value of kc at 25 deg C Michaelis Menten constant of CO2
    const double cmass = 12.0;      // atomic mass of carbon
    double alphaa = 0.5;            // fraction of PAR assimilated at ecosystem level relative to leaf level
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
    double temp, h, pO2, pa, tstress, k1, k2, k3, low, high, s, c1, c2, ko, kc, fac, sigma, phipi, NPP;
    int ptype;
    
    b[0] = RespFac(1);  // The first value of respfac is the leaf respiration factor depending on C3 or C4 photosynthesis
    b[1] = RespFac(1);  // The second value of Respfac is the root and stem respiration and is added outside this function
    alphaa = GreenBiomassRatio;
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
    
    agd = h * (je + jc - sqrt(pow((je + jc), 2.0) - 4.0 * theta * je * jc)) / (2.0 * theta); // Calculation of daily gross photosynthesis, Agd, gC/m2/day Eqn 2, Haxeltine & Prentice 1996
    if ((SatCorr > 0.0) && (Saturation(1) < SatCorr)) agd = (Saturation(1) / SatCorr) * agd;  // correction of gross photosynthesis for dry soil conditions
    rd = b[ptype - 1] * vm; // Daily leaf respiration, Rd, gC/m2/day Eqn 10, Haxeltine & Prentice 1996
    nd = agd - rd;  // Daily net photosynthesis (at leaf level), And, gC/m2/day
    adt= nd + (1.0 - h / 24) * rd; // Total daytime net photosynthesis, Adt, gC/m2/day Eqn 19, Haxeltine & Prentice 1996

    /* Merit
    NPP = Timestep * adt / 1000.0;  IS DIT WEL GOED?? in mijn versie is dit nd ----- NB NPP is lokaal gedefinieerd
    // convert from g C/m2/d to kg C m2/d

    NPP is returned!

    GPP = agd * Timestep / 1000;
    // convert from g C/m2/day to kg C m-2 per timestep

    */
    NPP = Timestep * adt / 1000.0;
    // convert from g C/m2/d to kg C m2/d
    GPP = agd * Timestep / 1000;
    // convert from g C/m2/day to kg C m-2 per timestep

    PlantRespiration = Timestep * C_CO2 * rd / 1000.0;  // convert from g C/m2/d to kg CO2 m2 per timestep
    //return (Timestep * nd / 1000.0);  // convert from g C/m2/d to kg C m2 per timestep
    return(NPP);
}


double LAICalc()
/* Calculates LAI (leaf area index) from photosynthesis 
 this includes not only Photosynthesis models 5 and 6 but also other models 
 so the Phenology parameter should also be defined for the other models*/
{
    double gdd = 0.0, LAI = 0.0, t, p, kgCadded, maxLAI, minLAI, PotLAIadded, LAIadded, maxGDD, seasonstart, autumnstart;
    BOOLEAN autumn = FALSE; // indication of spring/summer or atumn season
    
    LAIovershoot = 0.0;
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
        autumnstart = Phenology(8);
        if (StepNr < seasonstart) LAI = minLAI; // the first 30 days in the first year of the simulation get the minimumLAI, for the next years it depends on the BiomassSenenescence and weather
        if ((StepNr < seasonstart) && (BioMass == 0)) LAI = 0;
        if (DayOfTheYear != (seasonstart + 0.5 * Timestep)) //!!!!!moet dit niet > zijn ipv !=
        {   
            gdd = Timestep * (TData(StepNr) - Phenology(2));  // heat sum
            if (DayOfTheYear == (autumnstart + 0.5 * Timestep)) GrowingDegreeDays = maxGDD;
            if ((DayOfTheYear > autumnstart) || (DayOfTheYear < seasonstart)) autumn = TRUE; // autumn; in autumn GDD can decrease
            // if heat sum positive then calculculate LAI
            if (autumn) { // autumn: decrease growing degree days if daily heat sum < 0 and allow leaf senescence
                if (gdd < 0.0) {  
                    LeafSenescence = TRUE; 
                    GrowingDegreeDays += gdd; // in autumn GDD may decrease
                    if (GrowingDegreeDays < 0.0) GrowingDegreeDays = 0.0; // zero minimum of GDD
                } else LeafSenescence = FALSE; // positive daily heat sum: warm weather, further senescence delayed
            } else {
                if (gdd > 0.0) {
                    GrowingDegreeDays += gdd;
                    LeafSenescence = FALSE;
                }
            }
            if (GrowingDegreeDays > maxGDD) PotentialLAI = maxLAI; else PotentialLAI = minLAI + (maxLAI - minLAI) * GrowingDegreeDays / Phenology(3);
            // potential LAI can decrease in autumn, so does LAI
            kgCadded = NPP * ShootsFactor;
            // adapte by using the GreenBiomassRatio, the ratio over photosynthesizing biomass over total bomass
            PotLAIadded = GreenBiomassRatio * kgCadded / LAICarbonFraction; //LAI added according to primary production and GreenBiomassRatio
// Bug correction: previous coude allowed to exceed the maximum LAI strongly, which causes errors in the methane model.
// This has been corrected by introducing the LAIovershoot variable.
// If a LAI overshoot occurs, the carbon allaocation from primary production to shoot production is decreased, in favour
// of root production, if function OrgProd.

            if (CurrentLAI < PotentialLAI) {
                // Harvest or grazing has occurred, or LAI increase by primprod is higher than LAI increase based on gdd
                LAI = CurrentLAI + PotLAIadded;
                if (LAI > maxLAI) { // Restrict LAI to maximum LAI
                    LAIovershoot = (LAI - maxLAI) / maxLAI;
                    // LAI > maximum LAI should result in lower above-ground biomass allocation in favour of below-ground
                    LAI = maxLAI;
                }
            } else LAI = PotentialLAI;
        } else {
            LeafSenescence = FALSE;
            GrowingDegreeDays = 0.0;
            PotentialLAI = minLAI;
            CurrentLAI = PotentialLAI;
            LAI = minLAI;
            if (BioMass == 0) LAI = 0;
        }
    } else LAI = Phenology(4); // evergreen phenology with constant LAI AND NO HARVEST
    PreviousLAI = CurrentLAI;
    CurrentLAI = LAI;
    //cout << "LAI: " << LAI << " Biomass: " << BioMass << endl;
    return LAI; 
}

double RadProd()
/* Primary production from photosynthetically active radiation, Haxeltine and Prentice */

{
    double par, NPP, LAI;
    
    /* Calculation of Growing Degree Days; counter is reset at start of the year */
    
    par = PARcalc(); // PAR calculation
    LAI = LAICalc();
    NPP = PhotoSynthesis(LAI, par);
    return NPP;
}


double TundraProd()
// Primary production from photosynthetically active radiation, for tundra, Shaver et al, J. Ecology 2007
// datatype: 0 PAR data given in μmol photons m–2 s–1; 1: J/cm2 
{
    double par, NPP, LAI, gpp, er, temp, I, resp0, respbeta, pmax, pslope, CDconvfac, Cconvfac, ER, SatTop, SatC, HarC, RecovT;
    // resp0 Plant respiration at zero degrees
    // respbeta Temperature sensitivity factor plant respiration
    // pmax light-saturated photosynthetic rate per unit leaf area (μmol m–2 leaf s–1)
    // pslope is the initial slope of the light response curve (μmol CO2 μmol–1 photons)
    // NB: this model estimates a higer plant respiration than the Haxeltine and Prentice model, which includes only the leaf respiration;
    // here the entire plant respiration is estimated (including roots and stems)
    // Calculation of Growing Degree Days; counter is reset at start of the year
     

    par = PARcalc(); // PAR calculation
    LAI = LAICalc(); // LAI
    temp = TData(StepNr); // air temperatue
    resp0 = PhotoPar(1);
    respbeta = PhotoPar(2);
    pmax = PhotoPar(3);
    pslope = PhotoPar(4);
    SatTop = Saturation(1);
    SatC = SatCorr;
    HarC = HarvestCorrection(1);
    RecovT = HarvestCorrection(2);
    I = (par * 24 / DayLength); // Correct for day length, to have average PAR during day-time instead of 24 hours
    if ((LAI == 0.0) || (I == 0.0) || (temp < 0.0)) {PlantRespiration = 0.0; return (0.0);} // no use calculating when there is no PAR, LAI or low temperatures below the minimum requirement!
    //GPP calculated cf Shaver from LAI and incoming radiation
    // Plant respiration is their ER1 respiration model 
    er = resp0 * LAI * exp(respbeta * temp); //Plant resp μmol m–2 ground s–1)
    gpp = (pmax / KBeer) * log((pmax + pslope * I) / (pmax + pslope * I * exp(-KBeer * LAI))); // Gross Primary Production μmol m–2 ground s–1
    if (HDcount < RecovT) gpp = gpp * (HarC + ((1 - HarC)/ RecovT) * HDcount);
    // Correct for saturated soil
    if ((SatC > 0.0) && (SatTop < SatC)) gpp = (SatTop / SatC) * gpp;  // dry soil correction
    //convfac = Timestep * DayLength * 3600 * MOLWEIGHTC * 1.0e-9; // convert from μmol CO2 m–2 ground s–1 to kg C m2 per time step; photosynthesis only during the day
    //NPP = (gpp - er * DayLength / 24.0) * convfac; // from gpp, only the fraction of ecosystem respiration that occurs during the photoperiod is subtracted
    /* replaced by Merits version
    convfacC = Timestep * 24 *3600 * MOLWEIGHTC * 1.0e-9; // convert from μmol CO2 m–2 ground s–1 to kg C m2 per time step; photosynthesis only during the day
    NPP = (gpp - er) * convfacC; // from gpp, only the fraction of ecosystem respiration that occurs during the photoperiod is subtracted
    if (NPP < 0.0) NPP = 0.0;
    return NPP;*/
    // Merit version
    CDconvfac = Timestep * DayLength * 3600 * MOLWEIGHTC * 1.0e-9;
    Cconvfac = Timestep * 24 * 3600 * MOLWEIGHTC * 1.0e-9;
    ER = er * Cconvfac;
    GPP = gpp * CDconvfac;
    if (ER > GPP) ER = GPP;
    PlantRespiration = ER * MOLWEIGHTCO2 / MOLWEIGHTC;
    NPP = GPP - ER;
    return NPP;
}

/* Version Merit

double TundraProd()
//Primary production from photosynthetically active radiation, for tundra, Shaver et al, J. Ecology 2007 
{
    double par, NPP, LAI, gpp, er, temp, I, convfac, resp0, respbeta, pmax, pslope, CDconvfac, Cconvfac, ER, SatTop, SatC, HarC, RecovT, Correction;
    //double t_corr, Topt, Tb1;
    // To be added to parameter files:
    // resp0 Plant respiration at zero degrees
    //  respbeta Temperature sensitivity factor plant respiration
    // pmax light-saturated photosynthetic rate per unit leaf area (μmol m–2 leaf s–1)
    // pslope is the initial slope of the light response curve (μmol CO2 μmol–1 photons)

    
    // Calculation of Growing Degree Days; counter is reset at start of the year
*    par = PARcalc();
    LAI = LAICalc();
    temp = TData(StepNr);
    resp0 = PhotoPar(1);
    respbeta = PhotoPar(2);
    pmax = PhotoPar(3);
    pslope = PhotoPar(4);
    SatTop = Saturation(1);
    SatC = SatCorr;
*    HarC = HarCorr(1);
*    RecovT = HarCorr(2); 

    if ((LAI == 0.0) || (par == 0.0) || (temp < 0.0)) {PlantRespiration = 0.0; return (0.0);} // no use calculating when there is no PAR, LAI or low temperatures below the minimum requirement!
    // GPP calculated cf Shaver from LAI and incoming radiation
    // Plant respiration is their ER1 respiration model

    er = resp0 * LAI * exp(respbeta * temp); //Plant resp μmol m–2 ground s–1)
*    I = (par * 24 / DayLength); // Correct for day length, to haver average PAR during day-time    !!!!!!!!!!!!!!!! klopt dit wel met PARCalc?
    gpp = (pmax / KBeer) * log((pmax + pslope * I) / (pmax + pslope * I * exp(-KBeer * LAI))); // Gross Primary Production μmol m–2 ground s–1
*    // Correct gpp for recovery after harvest
*    if (HDcount < RecovT) gpp = gpp * (HarC + ((1 - HarC)/ RecovT) * HDcount);
    // Correct for saturated soil
    if ((SatC > 0.0) && (SatTop < SatC)) gpp = (SatTop / SatC) * gpp;                !!!!!!!!!!!!!!!dit zou voor Production model 3 ook moeten!
    // Conversions from micromols CO2 m-2 s-1 to kg C m-2 per timestep, corrected and uncorrected for DayLength
*    CDconvfac = Timestep * DayLength * 3600 * MOLWEIGHTC * 1.0e-9;
*    Cconvfac = Timestep * 24 * 3600 * MOLWEIGHTC * 1.0e-9;
*    ER = er * Cconvfac;
*    GPP = gpp * CDconvfac;           !!!!!!!!!!!!!!!!!!!! even nagaan of dit klopt, mogelijk beter gedaan dan in mijn versie GPP alleen gedurende de dag, ER 24 uur
                                    !!!!!!!!!!!!!!!!!!!!! MAARRR I wordt al gemiddeld over daglengte wat weer de PAR waarmee gerekend wordt beinvloedt!
*   if (ER > GPP) ER = GPP;         !!!!!!!!!!!!!!! dit is wsch ook correcter dan in mijn versie
*    PlantRespiration = ER * MOLWEIGHTCO2 / MOLWEIGHTC;
*    NPP = GPP - ER;
*    return NPP;
}

} */




void AddManure()                            // manure addition

{
  int i, j;
  double day, liquids, solids, totalliquids = 0.0, totalsolids = 0.0, excretion;

// Manure addition
  for (i = 1; i <= Manure.Rows(); i++)       // check if manure has to be added at the current day
  {
    day = Manure(i, 1);
    if ((day >= (DayOfTheYear - 0.5 * Timestep)) && (day < (DayOfTheYear + 0.5* Timestep)))       // the day should fall within the time step covering the current day
    {
      for (j = 1; j <= ManureLayers.Rows(); j++)
      {
        liquids = Manure(i, 2) * ManureFluidFrac * ManureLayers(j, 1);
        NewSOM(j, 2) = NewSOM(j, 2) + liquids;        // manure fluids
        TotalManure += liquids;
        totalliquids += liquids;
        solids = Manure(i, 2) * (1.0 - ManureFluidFrac) * ManureLayers(j, 2);
        NewSOM(j, 3) = NewSOM(j, 3) + solids;
        TotalManure += solids;
        totalsolids += solids;

// manure solids may be added to deeper layers if manure injection is used
      }
    }
  }
// Cattle excretion addition
  for (i = 1; i <= Grazing.Rows(); i++)     // Check if the current day of the year is within one of the grazing periods
  {
    if ((DayOfTheYear >= Grazing(i, 1)) && (DayOfTheYear < Grazing(i, 2))) {
        excretion = Timestep * Grazing(i, 4);
        NewSOM(1, 3) = NewSOM(1, 3) + excretion;
        TotalManure += excretion;

    }
// excretion is assumed to be mostly solid and added to the top layer only
  }
}

void CollectBioMass()
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
 BioMassRec(StepNr, 2) = BioMass + totalroots; kg C/m2 
 BioMassRec(StepNr, 3) = PrimProd; converted here from C/m2/timestep to mg CO2 m-2/hr
 BioMassRec(StepNr, 4) = PlantRespiration; converted here from kg CO2 m2 per timestep to mg CO2 m-2/hr
 BioMassRec(StepNr, 5) = Soil respiration + plant respiration - NPP = NEE  mg CO2 m-2/hr
 BiomassRec(StepNr, 6) = Soil respiration + plantrespiration = ecosystem respiration mg CO2 m-2/hr
 BioMassRec(StepNr, 7) = LitterLayer;
*/
    
    
  f = C_CO2*1.0e6/(24*Timestep);     // conversion factor from kg C/m-2/timestep into mg CO2 m-2/hr
  lt.Resize(NrLayers);
  for (i = 1; i <=NrOfSteps; i++)
  {
    BioMassRec(i, 3) = f * BioMassRec(i, 3);   // net photosynthesis (gross primary production - plant respiration) mg CO2 /m2/hr
    BioMassRec(i ,4) = BioMassRec(i, 4) * 1.0e6 / (24 * Timestep);  // plant respiration: transform from kg CO2/m2/timestep into mg CO2 /m2/hr
    lt.PutData(1, 1, LayerTime, i, 2, NrLayers);  // soil respiration
    s = lt.Sum();
    BioMassRec(i, 5)  = s - BioMassRec(i, 3);  // Net CO2 production/uptake incl. soil respiration
    BioMassRec(i, 6)  = s + BioMassRec(i, 4); // soil + plant respiration for comparing with dark chamber / nighttime measurements
  }
}
