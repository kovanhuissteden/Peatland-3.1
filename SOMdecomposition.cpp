/***************************************************************************
        SOMdecomposition.cpp  -  soil organic matter decomposition PEATLAND
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

  October 2003

  Total CO2 added as last column to output array Reservoirtime

  May 2012
  Calculation of photosyntheses based on photosynthetically active radiation
  added
  
 June-July
 adaptations to added photosynthesis
 added recording of above-ground litter/dead biomass and CO2 from above-ground litter
 ***************************************************************************/
#include <cmath>
#include "matrix.h"
#include "general.h"
#include "SOMdecomposition.h"

void EnviCor()
/* environmental correction factors for first order decomposition constants
aerobic decomposition of the SOM reservoirs, includes:
Temperature correction
pH correction
soil moisture and soil dryness corrections
Priming correction                                                         */
{
  int i, j, a;
  double t, pHcorr, drynesscorr, wetnesscorr, hs, h, s, growfac, dtcorr, CNcorr;
  Matrix roots, priming, tcorr;
  // initialize everything
  roots = RootMass / RootMass.Sum();
  priming.Resize(NrReservoirs);
  priming.Fill(1.0);
  tcorr.Resize(NrReservoirs);
  tcorr.Fill(0.0);
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i,4); // a is reference to soil horizon
    CNcorr = 1.0 + ((KPeatCN(1) - CNRatio(a)) * KPeatCN(2));   // correct k for CN based on empirical linear relation
    if (CNRatio(a) < 10) {CNcorr = ((KPeatCN(1) - 10) * KPeatCN(2));}
    if (CNRatio(a) > 55) {CNcorr = ((KPeatCN(1) - 55) * KPeatCN(2));}
    t = SoilTemp(i);  // temperature correction
    if (Q10orArrhenius == 0) {
        for (j = 1; j <= NrReservoirs; j++) tcorr(j) = pow(AerobicQ10(j), (t - T_ref) / 10.0);
            // temperature correction based on per reservoir specified Q10
    } else { //Arrhenius equation; Q10 is interpreted as molecular activation rate, specified per reservoir
        for (j = 1; j <= NrReservoirs; j++) tcorr(j) = exp((-AerobicQ10(j) / Rgas) * (1 / (t + 273.0) - 1 / (T_ref + 273.0)));
    }
    pHcorr = 1 / (1 + exp(-2.5 * (Layer_pH(a) - 5))); // pH correction cf ANIMO
    // the correction for soil dryness is linearly dependent on matric potential in a user defined range
    // defined in pF points. At the dry end of the range the factor is a fixed low valuse
    // TODO: replace by equations developed by Saurich et al
    drynesscorr = 1;                                  // correction for soil dryness
    if (MatricPotential(i) > pFpoints(1,1))
    {
      if (MatricPotential(i) < pFpoints(1,2))
      {
        drynesscorr =  1 - (1 - pFpoints(2,2)) * (MatricPotential(i) - pFpoints(1,1))/(pFpoints(1,2) - pFpoints(1,1));
      } else drynesscorr = pFpoints(2,2);
    }
    // correction for wetness, depends on water saturation and root mass
    hs = HalfSatPoint * (1 - roots(i) * RootAeration);     // diminish waterlogging correction for saturation/waterlogging in densely rooted zone
    s = Saturation(i);
    h = 2.0 * hs;
    if (s < h) wetnesscorr = s / h; else wetnesscorr = 1.0;  // correction depends linear on saturation
    // priming effect root exudates on slow C reservoirs (peat and humus)
    // depends on root density, relative growth rate of the vegation (more exudation at high rate), and springtime
    growfac = (NPP / Timestep - MinNPP) / (MaxNPP - MinNPP);    // growfac: the relative growth rate of the vegetation
    if (PrimingCorrection > 0)                             // priming correction based on growth rate, root density and time of the year (in spring very active exudation)
    {
      priming(1) = 1.0 + roots(i) * PrimingCorrection * growfac * SpringFactor;  // correction for peat
      priming(NrReservoirs) =  priming(1);                                 // correction for humus
    }
    // final calculation of environmental correction per reservoir and per layer
    for (j = 1; j <= NrReservoirs; j++) {
      // different calculation for peat reservoir, with CN correction
      if (j == 1) CorrFac(i, j) = priming(j) * tcorr(j) * pHcorr * drynesscorr * wetnesscorr * CNcorr;
      else CorrFac(i, j) = priming(j) * tcorr(j) * pHcorr * drynesscorr * wetnesscorr;
    }
  }
}   // end Envicor


void Decompose()
/* aerobic decomposition of SOM above the water table.*/

{
  int i, j, n, m;
  double h, unsatfraction, dt, transfermicrob, transferresist, soilT, anaertemperaturefact, ka, decomp, anaerobtotal = 0.0, aerobtotal = 0.0, tcorr = 2.0, litterTfact;
  Matrix k, aerob, anaerob, anaerobCO2, decomposedC;

  dt = Timestep / YEAR;                           // timestep in years
  anaerobCO2.Resize(NrLayers, NrReservoirs);     // initialize arrays
  anaerobCO2.Fill(0.0);                             // anaerobic CO2
  decomposedC.Resize(NrLayers, NrReservoirs);     // initialize arrays
  decomposedC.Fill(0.0);                             // anaerobic CO2
  anaerob.Resize(NrReservoirs);                 // anaerobic reservoir per layer
  AnaerobSum.Fill(0.0);                         // total anaerobic CO2 produced in timestep
  AnaerobSumRes.Fill(0.0);                         // total anaerobic CO2 produced in timestep per reservoir
  anaerob.Fill(0.0);
  CO2.Fill(0.0);
  PeatLoss = 0.0;
  //cout << CurrentGW << endl;
  for (i = 1; i <= NrLayers; i++)
  {
    soilT = SoilTemp(i); // temperature sensitivity calculation for anaerobic part
    h = Layers(i, 1) - CurrentGW;               // h is used to check if a layer is partly or entirely above the water table
/* calculate aerob and anaerobic fraction of organic C for layers at the groundwater table
   aerobic decomposition is only computed from the above the water table part */
    if (h <= 0.0) unsatfraction = 0.0; else if (h < LayerThickness) unsatfraction = 1.0 - h / LayerThickness; else unsatfraction = 1.0;  // fraction above water table
    aerob = NewSOM.Row(i) * unsatfraction;  // fraction of C above water table
    anaerob = NewSOM.Row(i) * (1.0 - unsatfraction); // fraction of anaerobic C below water table
    if (unsatfraction > 0.0)                      // aerobic decomposition above the water table
    {
      k = CorrFac.Row(i);                         // decomposition constants per layer, corrected for environmental factors
      k *= Kdecay;
      for (j = 1; j <= NrReservoirs; j++)
      {
        decomposedC(i, j) = aerob(j) * (1.0 - exp(- dt * k(j))); // decomposition of reservoir carbon, calculates what remains after timestep; C(t)= C(t-1) exp (-time*k)
        NewSOM (i,j) = NewSOM(i, j) - decomposedC(i,j); // subtract removed carbon from SOM reservoir
      }
      // PeatLoss += CO2(i, 1);  // total aerobic decomposition of peat carbon
      PeatLoss += decomposedC(i, 1);  // total aerobic decomposition of peat carbon
      for (m = 1; m <= 5; m++)  // correct for transfer ot microbial and resistant reservoir for the  first 5 reservoirs; decomposed C is moved to microbial and resistant fraction
      {
        transfermicrob = decomposedC(i, m) * SplitRes(m, 1);  // split removed carbon in CO2 and transfer of assimilated microbial biomass to microbial biomass reservoir
        NewSOM(i, 6) = NewSOM(i, 6) + transfermicrob;    // add microbially assimilated to microbial biomass reservoir
        decomposedC(i, m) = decomposedC(i, m) - transfermicrob;  // substract from decomposed carbon
        transferresist = decomposedC(i, m) * SplitRes(m, 2);  // transfer of resistant fraction
        NewSOM(i, 7) = NewSOM(i, 7) + transferresist;     // add to resistant reservoir
        decomposedC(i, m) = decomposedC(i, m) - transferresist;   // substract from decomposed carbon, remainder is true CO2 created
      }
      aerobtotal += decomposedC.SumRow(i); // total aerobically produced CO2 for carbon balance
      for (m = 1; m <= NrReservoirs; m++) CO2(i, m) += decomposedC(i, m);  // add produced CO2 to total CO2 array
    }  // end of above-water table calculation
    // anaerobic CO2 for layers that are partly or entirely below water table
    if (unsatfraction < 1.0)  // layer partly or entirely below water table
    {
        if (AnaerobicCO2 > 0.0) // anaerobic calculations only included if AnaerobicCO2 parameter is > 0
        { 
            if (Q10orArrhenius == 0) { // temperature correction for anaerobic decomposition
                anaertemperaturefact = pow(Q10Anaerobic, ((soilT - T_ref) / 10.0)); // assuming that the reference temperature for NO3, Fe/Mn and SO4 reduction is the same as for aerobic oxidation
            } else { // temperature correction cf Arrhenius equation
                anaertemperaturefact = exp((-Q10Anaerobic / Rgas) * (1 / (soilT + 273.0) - 1 / (T_ref + 273.0)));
            }
            for (j = 1; j <= NrReservoirs; j++)
            {
                ka = KAnaerobic(j) * anaertemperaturefact;
                anaerobCO2(i, j) = anaerob(j) * ( 1.0 - (exp(- dt * ka)));   // anaerobically produced CO2 as difference between size of anearobic reservoir before and after time step as above for aerobic CO2
                NewSOM(i, j) -= anaerobCO2(i, j);
            }
            for (m = 1; m <= 5; m++) {  // correct for transfer of microbial and resistant reservoir for the  first 5 reservoirs; decomposed C is moved to microbial and resistant fraction
                transfermicrob = anaerobCO2(i, m) * SplitRes(m, 3);  // transfer of assimilated microbial biomass to microbial biomass reservoir
                anaerobCO2(i, m) -= transfermicrob;  // correct anaerobically produced CO2 for assimilation
                NewSOM(i, 6) = NewSOM(i, 6) + transfermicrob; // add to microbial biomass reservoir
                transferresist = anaerobCO2(i, m) * SplitRes(m, 2);  // transfer of resistant fraction
                anaerobCO2(i, m) -= transferresist;  // correct anaerobically produced CO2 for assimilation
                NewSOM(i, 7) = NewSOM(i, 7) + transferresist; // add to resistant reservoir
            }
            for (m = 1; m <= NrReservoirs; m++) CO2(i, m) += anaerobCO2(i, m);       // add anaerobically produced CO2 to CO2 carbon
            AnaerobSum(i) = anaerobCO2.SumRow(i); // sum over all reservoirs for layer total
        }
    } // end anaerobic CO2 calculation
  }
  for (m = 1; m <= NrReservoirs; m++) AnaerobSumRes(m) += anaerobCO2.SumCol(m); // summation of anaerobic CO2 per reservoir
  if (Q10orArrhenius == 0) tcorr = AerobicQ10(5); // litter temperature decomposition correction
  litterTfact = pow(tcorr, ((TData(StepNr) - T_ref) / 10.0)); // temperature sensitivity of litter decomposition is fixed at Q10 = 2.0
  ka = KLitter * litterTfact;
  LitterDecomp = LitterLayer - LitterLayer * (exp(- dt * ka));  // Decomposition of surface litter
  LitterLayer -= LitterDecomp; // decrease litter layer carbon
  CarbonBalance(StepNr, 10) += aerobtotal* CONVKGCTOMOLC;
  CarbonBalance(StepNr, 12) += LitterDecomp * CONVKGCTOMOLC;
}    // end Decompose



void WriteSOMReservoirs()
/* writes contents of SOM reservoirs to output file */
{
  Matrix res, labile, total;
  int i;

    if (ProfileOutput.Contains(5))                   // sum of labile SOM
    {
      labile.Resize(NrLayers);
      res = NewSOM.Range(1, NrLayers, 2, 6);
      for (i = 1; i <= NrLayers; i++) labile(i) = res.SumRow(i);
      labile.Write(output5);
    }
    if (ProfileOutput.Contains(9))                   // sum of all SOM
    {
        total.Resize(NrLayers);
        res = NewSOM.Range(1, NrLayers, 1, 7);
        for (i = 1; i <= NrLayers; i++) total(i) = res.SumRow(i);
        total.Write(output9);
    }
    if (ProfileOutput.Contains(11))                  //  peat
    {
      res = NewSOM.Col(1);
      res.Transp();
      res.Write(output11);
    }
    if (ProfileOutput.Contains(12))                  // liquid manure
    {
      res = NewSOM.Col(2);
      res.Transp();
      res.Write(output12);
    }
    if (ProfileOutput.Contains(13))                  // solid manure
    {
      res = NewSOM.Col(3);
      res.Transp();
     res.Write(output13);
    }
    if (ProfileOutput.Contains(14))                  // exudates
    {
      res = NewSOM.Col(4);
      res.Transp();
      res.Write(output14);
    }
    if (ProfileOutput.Contains(15))                  // litter and roots
    {
      res = NewSOM.Col(5);
      res.Transp();
      res.Write(output15);
    }
    if (ProfileOutput.Contains(16))                  // microbes
    {
      res = NewSOM.Col(6);
      res.Transp();
      res.Write(output16);
    }
    if (ProfileOutput.Contains(17))                  // humus
    {
      res = NewSOM.Col(7);
      res.Transp();
      res.Write(output17);
    }
}

void CollectCO2()
/* Collect all CO2 and store in output arrays */
// record CO2 from CH4 oxidation with aerobic CO2 from decomposition, and record only the anaerobiclly produced CO2 to LayerAnaerobic 
// Also adds soil organic matter reservoir changes to the carbon balance
{
  int i;
  double f, c, ca, ct;
  Matrix storagechange;

  f = C_CO2*1000000.0/(24*Timestep);                         // conversion factor from kg C/m-2/timestep into mg CO2 m-2/hr
  LayerTime(StepNr, 1) = DayNr;
  ReservoirTime(StepNr, 1) = DayNr;
  AnaerobReservoirTime(StepNr, 1) = DayNr;
  ct = 0.0;
  for (i = 1; i <= NrReservoirs; i++)
  {
    c = f * CO2.SumCol(i);
    ReservoirTime(StepNr, i + 1) = c;                       // CO2 per reservoir
    ct += c;                                                // total CO2
  }
  // CO2 from litter decomposition is included in CO2 from top layer; however, it is included seperately in CarbonBalance
  c = f * CO2FromMethaneOx.Sum();                             // CO2 from methane
  ReservoirTime(StepNr, NrReservoirs + 2) = c;
  ct += c;
  c = f * LitterDecomp;
  ReservoirTime(StepNr, NrReservoirs + 3) = c;
  ct += c;
  ReservoirTime(StepNr, NrReservoirs + 4) = ct;             // total CO2 is last column of Reservoirtime
  for (i = 1; i <= NrReservoirs; i++) // Anaerobic CO2 per reservoir
  {
    ca = f * AnaerobSumRes(i);
    AnaerobReservoirTime(StepNr, i + 1) = ca;                       // CO2 per reservoir
    // summation of anaerob CO2 from non-CH4 anaerobic reactions
    // adding CO2 from methane production requires including CO2 per reservoir recording in methane production functions
  }

  for (i = 1; i <= NrLayers; i++)
  {
    LayerTime(StepNr, i + 1) = f * CO2.SumRow(i);
    LayerAnaerobic(StepNr, i) = f * AnaerobSum(i); // CO2 from anaerob reactions, including methane formation
  }
  LayerTime(StepNr, 2) += f * (CO2FromMethaneOx.Sum() + LitterDecomp);  // add CO2 evolved by methane oxidation
  BioMassRec(StepNr, 7) = LitterLayer; // kg C in Litter layer
  storagechange = NewSOM - OldSOM; // record changes in soil organic matter carbon reservoirs, + = increase, - decrease
  for (i = 1; i <= NrReservoirs; i++) CarbonBalance(StepNr, i + 2) = storagechange.SumCol(i) * CONVKGCTOMOLC;  // storage change carbon reservoirs
  CarbonBalance(StepNr, 11) = AnaerobSum.Sum() * CONVKGCTOMOLC;
  CarbonBalance(StepNr, 22) = (LitterLayer - OldLitter) * CONVKGCTOMOLC;
  OldLitter = LitterLayer;
  PeatDecay(StepNr, 1) = storagechange.SumCol(1);
  PeatDecay(StepNr, 2) = PeatLoss;
  for (i = 10; i <= 15; i++) CarbonBalance(StepNr, 23) += CarbonBalance(StepNr, i);  // sum of all carbon emission
  CarbonBalance(StepNr, 23) += CarbonBalance(StepNr, 19); //Carbon loss by harvest and grazing
  CarbonBalance(StepNr, 24) += CarbonBalance(StepNr, 1) + CarbonBalance(StepNr, 2);  // sum of incoming carbon
  for (i = 3; i <= 9; i++) CarbonBalance(StepNr, 25) += CarbonBalance(StepNr, i);  // sum of all carbon storage changes
  CarbonBalance(StepNr, 25) += CarbonBalance(StepNr, 17); // storage change CH4 in soil water
  for (i = 20; i <= 22; i++) CarbonBalance(StepNr, 25) += CarbonBalance(StepNr, i); // storage change biomass and litter
  CarbonBalance(StepNr, 26) += (CarbonBalance(StepNr, 24) - CarbonBalance(StepNr, 25) - CarbonBalance(StepNr, 23)); // total balance
} // end  CollectCO2
