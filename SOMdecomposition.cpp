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
  double t, pHcorr, drynesscorr, wetnesscorr, hs, h, s, growfac;
  Matrix roots, priming, tcorr;

  roots = RootMass / RootMass.Sum();              // initialize everything
  priming.Resize(NrReservoirs);
  priming.Fill(1.0);
  tcorr.Resize(NrReservoirs);
  tcorr.Fill(0.0);
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i,4); // a is reference to soil horizon
    t = SoilTemp(i);
    for (j = 1; j <= NrReservoirs; j++) tcorr(j) = pow(AerobicQ10(j), (t - T_ref) / 10.0); // temperature correction based on per reservoir specified Q10
    // Old code based on Arrhenius equation     
    /*if (t > 4)                                    // temperature correction cf Arrhenius equation downto 4 degrees
    {
      tcorr = exp((-MolAct / Rgas) * (1 / (t + 273.0) - 1 / (T_ref + 273.0)));
    } else
    {                                             // below +4 degrees straight decline to 0 and below 0 no decomposition
      if (t > 0) tcorr = exp((-MolAct / Rgas) * (1 / (277.0) - 1 / (T_ref + 273.0))) * (t / 4); else tcorr = 0;
    }*/
    pHcorr = 1 / (1 + exp(-2.5 * (Layer_pH(a) - 5))); // pH correction cf ANIMO
    // the correction forsoil dryness is linearly dependent on matric potential in a user defined range
    // defined in pF points. At the dry end of the range the factor is a fixed low valuse
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
    growfac = (PrimProd / Timestep - MinProd) / (MaxProd - MinProd);    // growfac: what is the relative growth rate of the vegetation?
    if (PrimingCorrection > 0)                             // priming correction based on growth rate, root density and time of the year (in spring very active exudation)
    {
      priming(1) = 1.0 + roots(i) * PrimingCorrection * growfac * SpringFactor;  // correction for peat
      priming(NrReservoirs) =  priming(1);                                 // correction for humus
    }
    // final calculation of environmental correction per reservoir and per layer
    for (j = 1; j <= NrReservoirs; j++) CorrFac(i, j) = priming(j) * tcorr(j) * pHcorr * drynesscorr * wetnesscorr;
  }
}   // end Envicor


void Decompose()
/* aerobic decomposition of SOM above the water table.*/

{
  int i, j, n, m;
  double h, unsatfraction, dt, transfermicrob, transferresist, soilT, anaertemperaturefact, ka, decomp, anaerobtotal = 0.0, aerobtotal = 0.0, litterTfact;
  Matrix k, aerob, anaerob, anaerobCO2, decomposedC;

  dt = Timestep / YEAR;                           // timestep in years
  anaerobCO2.Resize(NrLayers, NrReservoirs);     // initialize arrays
  anaerobCO2.Fill(0.0);                             // anaerobic CO2
  decomposedC.Resize(NrLayers, NrReservoirs);     // initialize arrays
  decomposedC.Fill(0.0);                             // anaerobic CO2
  anaerob.Resize(NrReservoirs);                 // anaerobic reservoir per layer
  AnaerobSum.Fill(0.0);                         // total anaerobic CO2 produced in timestep
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
    //cout << i << " layer " << unsatfraction <<  " unsatfraction" << endl;
    //anaerob.Disp();
    if (unsatfraction > 0.0)                      // aerobic decomposition above the water table
    {
      k = CorrFac.Row(i);                         // decomposition constants per layer, corrected for environmental factors
      k *= Kdecay;
      for (j = 1; j <= NrReservoirs; j++)
      {
        decomposedC(i, j) = aerob(j) * (1.0 - exp(- dt * k(j))); // decomposition of reservoir carbon, calculates what remains after timestep; C(t)= C(t-1) exp (-time*k)
        NewSOM (i,j) = NewSOM(i, j) - decomposedC(i,j); // subtract removed carbon from SOM reservoir
        // NewSOM(i, j) = aerob(j) * (exp(- dt * k(j)));  // decomposition of reservoir carbon, calculates what remains after timestep; C(t)= C(t-1) exp (-time*k)
        // CO2(i, j) = aerob(j) - NewSOM(i, j);           // difference between old and new reservoir content is CO2 (carbon) evolved by decomposition;
        // still uncorrected for transfer to secondary reservoirs
      }
      // PeatLoss += CO2(i, 1);  // total aerobic decomposition of peat carbon
      PeatLoss += decomposedC(i, 1);  // total aerobic decomposition of peat carbon
      /* This code is cumbersome, replaced by code below
      for (n = 6; n <= NrReservoirs; n++)            // correct CO2 and NewSOM for carbon transferred to microbial biomass and humus for the layer under consideration (layer i)
      // loop microbial biomass and resistant fraction reservoir
      {
        for (m = 1; m <= 5; m++)  // loop first 5 reservoirs that have only losses to microbial biomass and risistant fraction
        {
          transfer = CO2(i, m) * SplitRes(m, (n - 5));
          NewSOM(i, n) = NewSOM(i, n) + transfer;
          CO2(i, m) = CO2(i, m) - transfer;
        }
      } */
      // correct CO2 and NewSOM for carbon transferred to microbial biomass and humus for the layer under consideration (layer i)  
      for (m = 1; m <= 5; m++)  // correct for transfer ot microbial and resistant reservoir for the  first 5 reservoirs; decomposed C is moved to microbial and resistant fraction
      {
        transfermicrob = decomposedC(i, m) * SplitRes(m, 1);  // split removed carbon in CO2 and transfer of assimilated microbial biomass to microbial biomass reservoir
        NewSOM(i, 6) = NewSOM(i, 6) + transfermicrob;    // add microbially assimilated to microbial biomass reservoir
        decomposedC(i, m) = decomposedC(i, m) - transfermicrob;  // substract from decomposed carbon
        transferresist = decomposedC(i, m) * SplitRes(m, 2);  // transfer of resistant fraction
        NewSOM(i, 7) = NewSOM(i, 7) + transferresist;     // add to resistant reservoir
        decomposedC(i, m) = decomposedC(i, m) - transferresist;   // substract from decomposed carbon, remainder is true CO2 created
        /*transfermicrob = CO2(i, m) * SplitRes(m, 1);  // transfer of assimilated microbial biomass to microbial biomass reservoir
        NewSOM(i, 6) = NewSOM(i, 6) + transfermicrob; // add to microbial biomass reservoir
        CO2(i, m) = CO2(i, m) - transfermicrob;       // substract from CO2 carbon
        transferresist = CO2(i, m) * SplitRes(m, 2);  // transfer of resistant fraction
        NewSOM(i, 7) = NewSOM(i, 7) + transferresist; // add to resistant reservoir
        CO2(i, m) = CO2(i, m) - transferresist;       // substract from CO2 carbon*/
      }
      aerobtotal += decomposedC.SumRow(i); // total aerobically produced CO2 for carbon balance
      for (m = 1; m <= NrReservoirs; m++) CO2(i, m) += decomposedC(i, m);  // add produced CO2 to total CO2 array
    }  // end above water table calculation
    // anaerobic CO2 partly or entirely below water table

    if (unsatfraction < 1.0)  // layer partly or entirely below water table
    {
        if (AnaerobicCO2 > 0.0) 
        { // handling of partial anaerobic layer if calculation of anaerobically produced CO2 is requested
            anaertemperaturefact = pow(Q10Anaerobic, ((soilT - MethaneTRef) / 10.0)); // assuming that the reference temperature for Fe and S is the same as for methane
            for (j = 1; j <= NrReservoirs; j++)
            {
                ka = KAnaerobic(j) * anaertemperaturefact;
                anaerobCO2(i, j) = anaerob(j) * ( 1.0 - (exp(- dt * ka)));   // anaerobically produced CO2 as difference between size of anearobic reservoir before and after time step as above for aerobic CO2
                // cout << i << " " << j << " " << anaertemperaturefact << "|" << KAnaerobic(j) << "|" << ka << "|" << "|" << dt <<  "|" << exp(- dt * ka) << "|" << anaerobCO2(i, j) << endl;
                NewSOM(i, j) -= anaerobCO2(i, j);
                //cout << anaerobCO2(i, j) << "|" << NewSOM(i, j) << endl;
            }
            //cout << "-----------------" << endl;
            for (m = 1; m <= 5; m++) {  // correct for transfer ot microbial and resistant reservoir for the  first 5 reservoirs; decomposed C is moved to microbial and resistant fraction
                transfermicrob = anaerobCO2(i, m) * SplitRes(m, 3);  // transfer of assimilated microbial biomass to microbial biomass reservoir
                anaerobCO2(i, m) -= transfermicrob;  // correct anaerobically produced CO2 for assimilation
                NewSOM(i, 6) = NewSOM(i, 6) + transfermicrob; // add to microbial biomass reservoir
                transferresist = anaerobCO2(i, m) * SplitRes(m, 2);  // transfer of resistant fraction
                anaerobCO2(i, m) -= transferresist;  // correct anaerobically produced CO2 for assimilation
                NewSOM(i, 7) = NewSOM(i, 7) + transferresist; // add to resistant reservoir
                //cout << anaerobCO2(i, m) << "|" << NewSOM(i, m) << endl;
            }
            for (m = 1; m <= NrReservoirs; m++) CO2(i, m) += anaerobCO2(i, m);       // add anaerobically produced CO2 to CO2 carbon
            AnaerobSum(i) = anaerobCO2.SumRow(i); // sum over all reservoirs for layer total
            /* old code replaced for code above for partitioning
            for (n = 6; n <= NrReservoirs; n++) {           // correct CO2 and NewSOM for carbon transferred to microbial biomass and humus
                // NB: must be adapted eventuallyto different dissimilation/assimilation ratio of anaeroic deomposition
                for (m = 1; m <= 5; m++) {
                    transfer = CO2(i, m) * SplitRes(m, (n - 5));
                    NewSOM(i, n) = NewSOM(i, n) + transfer;
                    CO2(i, m) = CO2(i, m) - transfer;
                }
            } */
        }
    } // end anaerobic CO2 calculation
    
  }  
/*      }
    } else {
        if (AnaerobicCO2 > 0) {
            for (j = 1; j <= NrReservoirs; j++) { // anaerobic CO2 from sulfate reduction and Fe/Mn reduction IN COMPLETELY ANAEROBIC LAYERS
                ka = KAnaerobic(j) * anaertemperaturefact;
                anaerobCO2(i, j) = anaerob(i, j) * ( 1.0 - (exp(- dt * ka)));   // anaerobically produced CO2 as difference between size of anearobic reservoir before and after time step as above for aerobic CO2
                NewSOM(i, j) -= anaerobCO2(i, j);
                if (j <= 5) { // transfer of C to microbial reservoir and resistant fraction
                    transfermicrob = anaerobCO2(i, j) * SplitRes(j, 3);  // transfer of assimilated microbial biomass to microbial biomass reservoir
                    anaerobCO2(i, j) -= transfermicrob;  // correct anaerobically produced CO2 for assimilation
                    NewSOM(i, 6) = NewSOM(i, 6) + transfermicrob; // add to microbial biomass reservoir
                    transferresist = anaerobCO2(i, j) * SplitRes(j, 2);  // transfer of resistant fraction
                    anaerobCO2(i, j) -= transferresist;  // correct anaerobically produced CO2 for assimilation
                    NewSOM(i, 7) = NewSOM(i, 7) + transferresist; // add to resistant reservoir
                    CO2(i, j) += anaerobCO2(i, j);       // add anaerobically produced CO2 to CO2 carbon    
                } else {
                    CO2(i, j) += anaerobCO2(i, j);
                }
                AnaerobSum(i) += anaerobCO2.SumRow(j); // sum over all reservoirs for layer total
            }
            PeatLoss += CO2(i, 1);
        } else {
            for (j = 1; j <= NrReservoirs; j++) CO2(i, j) = 0.0;
        }
    }
  }*/
  litterTfact = pow(2.0, ((TData(StepNr) - T_ref) / 10.0)); // temperature sensitivity of litter decomposition is fixed at Q10 = 2.0
  ka = KLitter * litterTfact;
  LitterDecomp = LitterLayer - LitterLayer * (exp(- dt * ka));  // Decomposition of surface litter
  // CO2(1,5) = CO2(1,5) + LitterDecomp; // add CO2 from surface litter to CO2 from top layer
  LitterLayer -= LitterDecomp; // decrease litter layer carbon
  // will be added to collectCO2
  // AnaerobSum.Disp();
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
  double f, c, ct;
  Matrix storagechange;

  f = C_CO2*1000000.0/(24*Timestep);                         // conversion factor from kg C/m-2/timestep into mg CO2 m-2/hr
  LayerTime(StepNr, 1) = DayNr;
  ReservoirTime(StepNr, 1) = DayNr;
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
  for (i = 1; i <= NrLayers; i++)
  {
    LayerTime(StepNr, i + 1) = f * CO2.SumRow(i);
    LayerTime(StepNr, i + 1) += f * (CO2FromMethaneOx(i) + LitterDecomp);  // add CO2 evolved by methane oxidation
    // if (AnaerobicCO2 > 0) LayerAnaerobic(StepNr, i) = f * AnaerobSum(i); 
    LayerAnaerobic(StepNr, i) = f * AnaerobSum(i); // CO2 from anaerob reactions, including methane formation
  }
  BioMassRec(StepNr, 7) = LitterLayer; // kg C in Litter layer
  storagechange = NewSOM - OldSOM; // record changes in soil organic matter carbon reservoirs, + = increase, - decrease
  for (i = 1; i <= NrReservoirs; i++) CarbonBalance(StepNr, i + 2) = storagechange.SumCol(i) * CONVKGCTOMOLC;  // storage change carbon reservoirs
  CarbonBalance(StepNr, 11) = AnaerobSum.Sum() * CONVKGCTOMOLC;
  CarbonBalance(StepNr, 22) = (LitterLayer - OldLitter) * CONVKGCTOMOLC;
  OldLitter = LitterLayer;
  PeatDecay(StepNr, 1) = storagechange.SumCol(1);
  PeatDecay(StepNr, 2) = PeatLoss;
  for (i = 10; i <= 14; i++) CarbonBalance(StepNr, 23) += CarbonBalance(StepNr, i);  // sum of all carbon emission
  CarbonBalance(StepNr, 23) += CarbonBalance(StepNr, 19); //Carbon loss by harvest and grazing
  CarbonBalance(StepNr, 24) += CarbonBalance(StepNr, 1) + CarbonBalance(StepNr, 2);  // sum of incoming carbon
  for (i = 3; i <= 9; i++) CarbonBalance(StepNr, 25) += CarbonBalance(StepNr, i);  // sum of all carbon storage changes
  CarbonBalance(StepNr, 25) += CarbonBalance(StepNr, 17); // storage change CH4 in soil water
  for (i = 20; i <= 22; i++) CarbonBalance(StepNr, 25) += CarbonBalance(StepNr, i); // storage change biomass and litter
  CarbonBalance(StepNr, 26) += (CarbonBalance(StepNr, 24) - CarbonBalance(StepNr, 25) - CarbonBalance(StepNr, 23)); // total balance
} // end  CollectCO2
