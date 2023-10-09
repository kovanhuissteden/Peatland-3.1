/***************************************************************************
                methane.cpp  - methane model PEATLAND cf Walter(2000)
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

  Bug corrected: at groundwater levels of 0, no methane production occurred
  because of resulting 0 timestep for diffusion equation solution
  calculation of diffusion equation timestep adjusted

  June 2003

  IMPORTANT BUG CORRECTION
  Methane plant and ebullution fluxes were not multiplied by time step length
  resulting in fluxes sensitive to the time step - Crrected by multiplying
  every flux by dt
  Time step mechanism also adapted cf Appelo & Postma

  October 2003

  5th column added to TotalMethane array with the total methane flux
  (summed ebullition, plantflux and diffusive flux)

  April 2004

  Maximum methane concentration in soil profile increased; concentrations over
  10000 appear to occur in reality!

  December 2006

  Partial anaerobic conditions above water table added

  November/december 2007

  Addition of a time lag for anaerobe development and methane production
  after rapid saturation of a layer
 
  July 2012
  Tested for bugs
  Adaptation of plant CH4 flux to new photosynthesis module
 ***************************************************************************/

#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;
#include "matrix.h"
#include "general.h"
#include "methane.h"


void Methane()
/* methane fluxes modified after the model of Walther (2000), Global Biogeochemical Cycles 14, 745 - 765.*/
{
/*
  Translation old variables in Matlab version to new variables:
  lm, lb: Layers col. 2, 3
  c0: MethProfile
  sat: Saturation  : similar as sat, totally saturated = 0
  por: PoreVol
  tp, T: SoilTemp
  newc: NewSOM
  m, n, len: NrLayers, NrReservoirs
  topsat: TopSat;
  peat: column 1 NewSOM
*/
  int i, j, nsteps, a;
  double dtfac, dt, dth, ebullflux = 0.0, plantflux = 0.0, rho, next, prev, topflux, ttime, methair;
  Matrix CH4oxidation, cres, cresstart, va, diff, cc, ccn, ccstart, dc, totprod, ebullrate, plantrate, prod, CO2production, totalCO2, anaerob, anaerobCO2, Closs, Cremoved, mp;

// ALL CALCULATIONS IN MILLIMOLES C
 
  MethaneFlux.Fill(0.0);                  // methane flux totalled per time step
  Closs.Resize(NrLayers, NrReservoirs);   // carbon loss from reservoirs in millimoles C in layer volume
  Cremoved = Closs;                       // C removed be methanogenesis per methane model time step  
  CH4oxidation.Resize(NrLayers);          // logs oxygen consumption by methane oxidation in plants in millimoles/layer
                                          // for later use if O2 diffusion is added to the model
  CO2production.Resize(NrLayers);         // logs all CO2 derived from methane oxidation in soil millimoles/layer
  totalCO2.Resize(NrLayers);              // logs all CO2 generated by oxidation during transport
  anaerobCO2.Resize(NrLayers);            // logs all anaerobic CO2 generated during methane production in millimol C / layer
  anaerob.Resize(NrLayers);               // per model time step anaerobic CO2 generated during methane production in millimol C / layer
  // cres.Mult(NewSOM, MethaneReservoirs);   // sum and weigh easily decomposeable reservoirs
  cres = NewSOM;                          // cres copies NewSOM to maintain an independent copy during calculations
  cres *= (1.0e6 / MOLWEIGHTC);           // convert C reservoir from kg C to millimoles C in layer volume
  cresstart = cres;                       // cresstart stores carbon reservoir at start of time step
  va = PoreVol * Saturation;              // set up volumetric air content, diffusion coeff and temperature arrays
  diff = va* (0.66 * MethaneDiff);        // set up diffusion in aerated zone and saturated zone
  
  for (i = TopSat; i <= NrLayers; i++)    // correct diffusion coefficients for ice and water content
  {
    a = (int)Layers(i, 4);                // soil profile horizon number
    if (Ice(i) > 0) diff(i) = 0.66 * (DBD(a) * UnFrozen(i) / DensWater) * MethaneDiffWater; else diff(i) = 0.66 * PoreVol(i) * MethaneDiffWater;
// With presence of ice the unfrozen water content is taken here to account for smaller diffussivity in the presence of ice
  }
  dtfac = 1 / (0.25 * pow(LayerThickness, 2.0) / diff.Max());  // determine safe time step for numerical solution of PDE
  nsteps = (int)(Timestep * dtfac);
  if (nsteps < (MINMODELSTEPS * (int)Timestep))        // minimum time step is one minute (Timestep is in days)
  {
    nsteps = MINMODELSTEPS * (int)Timestep;
  }
  if (nsteps > (MAXMODELSTEPS * (int)Timestep))
  {
      cout << METHANE_ERROR4 << endl;
      exit(EXIT_FAILURE);
  }
  methair = MethaneAir / (273.15 + TData(StepNr));  // correction of surface aire CH4 concentration for air temperature
  dt = Timestep/nsteps;                    // PDE solution timestep in days
  dth = dt * 24.0;                         // PDE timestep in hours for ebullition rate. plant flux and production
  ttime = 0.5 * dt;			               // time in days since start of model time step for calculation of depression of methane production after rapid soil saturation
  cc = MethProfile;                        // cc, ccn: intermediate results methane concentration profile during diff equation solution iteration
  ccstart = MethProfile;                   // soil CH4 concentration at start of time step, for calculation of storage change
  mp.Resize(NrLayers);                     // CH4 concentration in millimol CH4 per layer volume
  dc.Resize(NrLayers);                     // dc: intermediate variable to calculate ebullition. plant flux, production per time step
  ccn.Resize(NrLayers);                    //  ccn: summed plant flux, ebullition, production
  totprod.Resize(NrLayers);                // collects total CH4 production millimol per layer per timestep
  prod.Resize(NrLayers);                   // instantaneous production
  ebullrate.Resize(NrLayers);               // ebullition in millimol per hour
  plantrate.Resize(NrLayers);               // plant flux in millimol per hour
  for (i = 1; i <= nsteps; i++)            // PDE solution
  {
    dc.Fill(0.0);                          // dc: change in concentration due to plant flux, ebullition, production, oxidation
    ebull(ebullflux, ebullrate, cc);       // ebullition
    dc -= ebullrate * dth;
    ebullflux = ebullflux * dth;
// bubbles reach the soil surface if groundwater is within the top half of the first layer)
// otherwise the gas reaches the unsaturated zone and is subjected there to oxidation and diffusion
    if (CurrentGW >= (-0.5 * LayerThickness))  // record ebullition flux if top layer is more than 50% saturated, else add it to the first unsaturated layer
    {
      MethaneFlux(1) = MethaneFlux(1) + ebullflux;
    } else if (TopSat > 1) dc(TopSat - 1) = dc(TopSat - 1) + ebullflux; else  dc(1) = dc(1) + ebullflux;
    ccn = cc + dc;
    planttrans (plantflux, CH4oxidation, CO2production, plantrate, ccn); // plant flux
    dc = -plantrate * dth;
    // CH4oxidation *= dth;  // CH4 consumed by plant oxidation in time step (NB CH4oxidation = CO2production in planttrans)
    totalCO2 += CO2production * dth; // total CO2 production from oxidation
    //plantox += CO2production.Sum(); // Total plant oxidized CH4
    // previous: plant oxidation was added to top layer :     // CH4oxidation(1) += plantox * dth; 
    MethaneFlux(2) = MethaneFlux(2) + plantflux * dth; // record plant flux
    ccn += dc;
    methaneprod(prod, cres, ccn, CO2production, anaerob, Cremoved, ttime); // methane production including oxidition by methanotrophs above water table
    //cout << (prod.Sum() + anaerob.Sum() + CO2production.Sum()) << " " << Cremoved.Sum() << endl;
    cres -= Cremoved * dth;
    totprod += prod * dth;  // total CH4 production
    anaerobCO2 += anaerob * dth; //Anaerobically procuced CO2
    dc = prod * dth; // change due to production
    ccn += dc; // add change due to CH4 production / oxidation to concentration profile;
    totalCO2 += CO2production * dth; // add CH4 oxidation in the soil above groundwatertable to CO2 oxidation in the plant system
    CH4oxidation = totalCO2; // CH4oxidation is not used any further here, it is kept for later use to calculate oxygen
    for (j = NrLayers; j >= 1; j--)                               // iteration of diffusion equation solution
    {
      rho = diff(j) * dt / pow(LayerThickness, 2.0);
      if (j == 1) next = methair; else next = ccn(j - 1);      // upper boundary condition: concentration is equal to methane concentration in atmosphere
      if (j == NrLayers) prev = ccn(j); else prev = ccn(j + 1);   // lower boundary condition: dC/dz == 0 no change of concentration gradient
      cc(j) = rho * next + (1 - 2 * rho) * ccn(j) + rho * prev;   // classic explicit solution
    }
    if (cc.Min() < 0.0)    // possible error due to numerical instability -  should not occur - usually generated by anomalous CH4 production
    {
      cc.Disp();
      cout << METHANE_ERROR2 << endl;
      exit(EXIT_FAILURE);
    }
    if (cc.Max() > MAXSOILCH4)
    {
      cout << METHANE_ERROR3 << endl;
    }
    topflux = diff(1) * (cc(1) - methair) / LayerThickness;     // diffusive flux at top of profile (note: flux per day)
    MethaneFlux(3) = MethaneFlux(3) + topflux * dt;                // flux per time step
    cc(1) -= topflux * dt;
    MethProfile = cc;
    ttime += dt;
  }
  MethProfile = cc;
  // handle C loss from reservoirs; all data of CO2 and CH4 production remains in millimoles but for substraction of C reservoirs is converted to kg C 
  Closs =  cresstart - cres;  // substract carbon loss by methanogenic decomposition from carbon reservoirs
  Closs *= CONVCH4CTOKGC;
  NewSOM -= Closs;  // update carbon reservoir size by subtracting C lost by CH4 production
  // anarobicCO2 and CO2production from oxidation also have to be added proberly to total anaerobc CO2 and CO2 from methane
  // this is done in function CollectCO2(); CO2 from oxidation is added to top layer
  CarbonBalance(StepNr, 14) = totalCO2.Sum() / 1000.0;   // add oxidized CH4 to carbon balance tracking
  CO2FromMethaneOx = totalCO2 * CONVCH4CTOKGC;                     // CO2 evolved from methane oxidation, convert from millimol CH4-C to kg C
  // CarbonBalance(StepNr, 16) = anaerobCO2.Sum() / 1000.0;
  anaerobCO2 = anaerobCO2 * CONVCH4CTOKGC;   // convert anaerobically generated CO2 to kg C per timestep
  AnaerobSum += anaerobCO2;
  CarbonBalance(StepNr, 13) = MethaneFlux.Sum() / 1000.0;   // add methane carbon to carbon balance
  MethaneFlux *= MOLWEIGHTCH4 / (24 * Timestep);             // conversion from millimoles per timestep to mg per hr
  //MethaneFlux.Disp();
  //cc.Disp();
  
  TotalMethane.PutData(StepNr, 2, MethaneFlux);                     // store current flux in result array
  TotalMethane(StepNr, 1) = DayNr;
  TotalMethane(StepNr, 5) = MethaneFlux.Sum();
  // CH4 soil storage change; + is increase; conversion from millimol/m3 the summed mol over layers
  CarbonBalance(StepNr, 17) = (cc.Sum() - ccstart.Sum()) /1000.0; 
  // convert methane profile back to millimol/m3 in soil pore volume for registration in output in the same unit as the initial profile
  mp = (MethProfile / PoreVol) / LayerThickness;
  if (ProfileOutput.Contains(3)) mp.Write(output3);        // write methane profile to output file if required
}   // end Methane



void ebull(double &flux, Matrix &rate, Matrix &methconc)
/*  calculates ebullition flux and rate of removal of CH4 (sink) by ebullition
    no entrapment cf Walter et al of bubbles is assumed, entrapment is assumed to occur within the same layer
    Input:
    methconc:   methane concentration, mol CH4 in layer
    Output:
    flux:       total CH4 bubble flux in mol CH4 over all layers summed
    rate:       rate of CH4 removal in millimol per layer per hour
*/
{
  int i;
  double mcmax;


  rate.Fill(0.0);
  for (i = TopSat; i <= NrLayers; i++)
  {
    if (Ice(i) > 0) break;    // stop if ice is encountered, no ebullution from ice or unfrozen layers below!
// ebullition occurs when max pore water concentration is exceeded
// Unit Ebullition rate constant MethaneERateC (m3/hr) therefore quantities of CH4 need to be converted to Peatland model time step (days)
// Unit methconc: millimol ch4 / layer volume
// MethaneMaxConc is recalculated in the same unit in InitMethaneModel()
    mcmax = PoreVol(i) * MethaneMaxConc;
    if (methconc(i) > mcmax) rate(i) = MethaneERateC * (methconc(i) - mcmax);
  }
  flux = rate.Sum();     // integrate over all layers to obtain total flux
}  // end ebull


void planttrans(double &flux, Matrix &CH4oxidation, Matrix &CO2production, Matrix &rate, Matrix &methconc)
/*  calculates plant transport flux and rate of removal of CH4 (sink) by plant transport
    Input:
    methconc:   methane concentration

    Output:

    flux:       total CH4  flux millimol removed from all layers
    CH4oxidation:  total CH4 oxidized in plant root system 
    CO2production: toatal CO2 produced by methane oxidation;
    rate:       rate of CH4 removal millimol removed from layer per hour
*/
{
    int i;
    double fgrow, minLAI =0.05;

    rate.Fill(0.0);
    if (ProductionModel < 3) { 
        if (MaxProd > 0) fgrow = GrowFuncConst * (PrimProd / Timestep) / MaxProd;  else fgrow = 0.0;  // growth function for plant transport - in the Walter Heimann model LAI is used but this is not available with the simple production models, therefore fgrow is scaled according to primary production
    } else { // in the original Walter-Heimann model, fgrow is determined by LAI
        fgrow = GrowFuncConst * CurrentLAI;
        if ((fgrow == 0.0) && (SoilTemp(1) > 0.0)) fgrow = GrowFuncConst * minLAI; // this allows for a minimum flux when LAI is zero and upper soil layer is not frozen
    } 
    rate = RootDistrib * methconc;
    rate *= (MethanePRateC * MethanePType * fgrow); // CH4 removed in millimol CH4 per layer in integration time step
    for (i = 1; i <= NrLayers; i++)           // warning message for strange results - should not occur but can arise from erroneous parameters
    // e.g. MethanePRateC * LAI higher than 1
    {
        if (rate(i) > methconc(i))
        {
            rate(i) = methconc(i) - MethaneAir / (273.15 + TData(StepNr));
            cout << METHANE_ERROR1 << endl;
        }
    }
    CH4oxidation = rate * MethanePlantOx; // oxidation in root system (in millimol O2)
    CO2production = rate * MethanePlantOx; // oxidation in root system (in millimol CO2)
    // rate = rate * (1.0 - MethanePlantOx);
    flux = rate.Sum();   // total flux in miilimol CH4 over all layers summed
    flux *= (1.0 - MethanePlantOx);
}    // end planttrans




void methaneprod(Matrix &prod, Matrix &labileC, Matrix &mc, Matrix &oxidized, Matrix &anaerob, Matrix &Cremoved, double ttime)
/*  Below water table methane production per hour from soil organic matter and
    oxidation in unsaturated zone as function of soil CH4 concentration

Input:

labileC : total labile C soil in micromol;
mc      : methane concentration

Output:

prod    : total production or oxidation
Cremoved: carbon removed from carbon reservoirs
ttime   : time in days since start of model time step for calculation of rapid saturation depression of methane production
oxidized: matrix with moles of methane oxidized per layer
anaerob: anaerobically prduced CO2

*/

{
  int i, j;
  double t, f, anaerobe, fsat, pprod, pcons, resprod, satconst, delayfac, allsteps;

  prod.Fill(0.0);
  oxidized.Fill(0.0);
  anaerob.Fill(0.0);
  Cremoved.Fill(0.0);
  for (i=1; i <= NrLayers; i++)
  {
//  calculate lag factor for rapid saturation of layers
    allsteps = (StepNr + 1) * Timestep;
    if (allsteps > LastSatTime(i))     // only for layers that have been unsaturated for some time
    {
      if (AnaerobeLagFactor > 0.0) delayfac = 1 - exp(-AnaerobeLagFactor * (LastSatTime(i) + ttime) * labileC(i,1)); else delayfac = 1.0;
    } else delayfac = 1.0;
    t = SoilTemp(i);				// soil temperature
    if (PartialAnaerobe > 1.0) 			// determination of anaerobic fraction, fsat value is 0.0 when comletely saturated with water
    {
        fsat = Saturation(i);
        satconst = 1.0/PartialAnaerobe;
        if (fsat < satconst)
        {
            anaerobe = 1.0 - PartialAnaerobe * fsat; 
        } else anaerobe = 0.0;
    } else					// no anaerobe fraction above the water table assumed
    {
        if (Saturation(i) == 0.0) anaerobe = 1.0; else anaerobe = 0.0;
    }
    if (t > 0.0)
    {
      pprod = 0.0;
      pcons = 0.0;
      if (anaerobe > 0.0)              // production below groundwater table
      {
        f = pow(MethaneQ10, ((t - MethaneTRef) / 10));  //temperature correction factor f
        for(j=1; j<=NrReservoirs; j++)  // production calculated by C reservoir
        {
            if (MethaneReservoirs(j) > 0.0)  //only the reservoirs with a labilityfactor (0-1) above 0 are included
            {
                 resprod = anaerobe * f * delayfac * MethaneR0Corr(i) * MethaneReservoirs(j) * labileC(i, j);  // produced CH4 (+CO2)
                 Cremoved(i, j) = resprod;  // decrease of labile C reservoir by produced CH4 (milliMol)
                 pprod += resprod; // add to total production for layer
            }
        }
        // from the consumed carbon fraction, a part of the carbon is produced as CO2, another part as CH4, depending on type of reaction (acetate splitting or CO2 reduction
        //this is indicated by the CO2CH4ratio
        anaerob(i) = CO2CH4ratio * pprod;  // records anaerobically produced CO2 by acetate splitting
        pprod = pprod * (1.0 - CO2CH4ratio); // decrease produced CH4-C with produced CO2
      }
      if (anaerobe < 1.0)                                              // oxidation
      {
        f = pow(MethaneOxQ10, ((t - MethaneTRef) / 10));
        pcons = (1.0 - anaerobe) * f * MethaneVmax * mc(i) / (MethaneKm + mc(i)); // oxidation is specified here as negative 'production'
        if (pcons >= mc(i))                              // all methane oxidized
        {
          pcons = mc(i) - 0.001 * MethaneAir / (273.15 + TData(StepNr));
// a small fraction of CH4 (<< atmosphere concentration) is assumed to remain to prevent PDE solution problems in cases with deep and rapid falling groundwater tables
        }
      }
      oxidized(i) += pcons;
      prod(i) = pprod - pcons;                      // total production minus consumption
    }
  }
}    // end methaneprod

