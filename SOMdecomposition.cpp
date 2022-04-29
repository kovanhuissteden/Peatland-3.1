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
  double t, tcorr, pHcorr, drynesscorr, wetnesscorr, hs, h, s, growfac;
  Matrix roots, priming;

  roots = RootMass / RootMass.Sum();              // initialize everything
  priming.Resize(NrReservoirs);
  priming.Fill(1.0);
  for (i = 1; i <= NrLayers; i++)
  {
    a = (int)Layers(i,4); // a is reference to soil horizon
    t = SoilTemp(i);
    if (t > 4)                                    // temperature correction cf Arrhenius equation downto 4 degrees
    {
      tcorr = exp((-MolAct / Rgas) * (1 / (t + 273.0) - 1 / (T_ref + 273.0)));
    } else
    {                                             // below +4 degrees straight decline to 0 and below 0 no decomposition
      if (t > 0) tcorr = exp((-MolAct / Rgas) * (1 / (277.0) - 1 / (T_ref + 273.0))) * (t / 4); else tcorr = 0;
    }
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
    h = 2 * hs;
    if (s < h) wetnesscorr = s / h; else wetnesscorr = 1;  // correction depends linear on saturation
    // priming effect root exudates on slow C reservoirs (peat and humus)
    growfac = (PrimProd / Timestep - MinProd) / (MaxProd - MinProd);    // growfac: what is the relative growth rate of the vegetation?
    if (PrimingCorrection > 0)                             // priming correction based on growth rate, root density and time of the year (in spring very active exudation)
    {
      priming(1) = 1.0 + roots(i) * PrimingCorrection * growfac * SpringFactor;  // correction for peat
      priming(NrReservoirs) =  priming(1);                                 // correction for humus
    }
    // final calculation of environmental correction per reservoir and per layer
    for (j = 1; j <= NrReservoirs; j++) CorrFac(i, j) = priming(j) * tcorr * pHcorr * drynesscorr * wetnesscorr;
  }
}   // end Envicor


void Decompose()
/* aerobic decomposition of SOM above the water table.*/
{
  int i, j, n, m;
  double h, f, dt, transfer, litterdecomp, soilT, anaertemperaturefact, ka, decomp;
  Matrix k, active, anaerob, anaerobCO2;

  dt = Timestep / YEAR;                           // timestep in years
  anaerobCO2.Resize(NrReservoirs);
  AnaerobSum.Fill(0.0);
  anaerob.Fill(0.0);

  for (i = 1; i <= NrLayers; i++)
  {
    soilT = SoilTemp(i); // temperature sensitivity calculation for anaerobic part
    anaertemperaturefact = pow(Q10Anaerobic, ((soilT - MethaneTRef) / 10)); // assuming that the reference temperature for Fe and S is the same as for methane
    anaerobCO2.Fill(0.0);
    if (Saturation(i) > 0.0)                      // decomposition only above the water table
    {
      k = CorrFac.Row(i);                         // decomposition constants per layer, corrected for environmental factors
      k *= Kdecay;
      h = Layers(i, 1) - CurrentGW;               // part of layer above water table
/* calculate active and anaerobic fraction of organic C for layers at the groundwater table
   decomposition is only computed from the active part, the anaerobic part is kept unchanged */
      if (h < LayerThickness) f = h / LayerThickness; else f = 1.0;   // fraction above water table
      active = NewSOM.Row(i) * f;
      if (f < 1) anaerob = NewSOM.Row(i) * (1.0 - f); // anaerobic C below water table
      for (j = 1; j <= NrReservoirs; j++)
      {
         // kcorr = k(j);
        NewSOM(i, j) = active(j) * (exp(- dt * k(j)));  // decomposition; C(t)= C(t-1) exp (-time*k)
        CO2(i, j) = active(j) - NewSOM(i, j);           // CO2 evolved - uncorrected for transfer to secondary reservoirs
        // dit is eigenlijk wat erover blijft ipv afbraak?
        //!!!!!!!!!!!!!!!!!!!!!! dit is CO2-C - beter een andere variabele naam geven
      }
      PeatLoss += CO2(i, 1);
      for (n = 6; n <= NrReservoirs; n++)            // correct CO2 and NewSOM for carbon transferred to microbial biomass and humus
      {
        for (m = 1; m <= 5; m++)
        {
          transfer = CO2(i, m) * SplitRes(m, (n - 5));
          NewSOM(i, n) = NewSOM(i, n) + transfer;
          CO2(i, m) = CO2(i, m) - transfer;
        }
      }
      // add decomposition of litter layer to that of upper layer
      litterdecomp = LitterLayer - LitterLayer * (exp(- dt * KLitter));
      CO2(1,5) = CO2(1,5) + litterdecomp; 
      LitterLayer -= litterdecomp;
      if (f < 1)
      {
        if (AnaerobicCO2 > 0) { // handling of partial anaerobic layer 
            for (j = 1; j <= NrReservoirs; j++) {
                ka = KAnaerobic(j) * anaertemperaturefact;
                //!!!!!!!!!!!!!!!!!!!this should be one as above for aerobic decomposition!!!!!
                anaerobCO2(j) = anaerob(j) * ( 1.0 - (exp(- dt * ka)));   // 
                NewSOM(i, j) -= anaerobCO2(j);
                CO2(i, j) += anaerobCO2(j);
                AnaerobSum(i) += anaerobCO2(j);
            }
            PeatLoss += CO2(i, 1);
            for (n = 6; n <= NrReservoirs; n++) {           // correct CO2 and NewSOM for carbon transferred to microbial biomass and humus
                // NB: must be adapted eventuallyto different dissimilation/assimilation ratio of anaeroic deomposition
                for (m = 1; m <= 5; m++) {
                    transfer = CO2(i, m) * SplitRes(m, (n - 5));
                    NewSOM(i, n) = NewSOM(i, n) + transfer;
                    CO2(i, m) = CO2(i, m) - transfer;
                }
            }
        } else for (j = 1; j <= NrReservoirs; j++) NewSOM(i, j) += anaerob(j); // add inert (waterlogged) part back to C reservoir
        // !!!!!!!!!!!! KLOPT DIT WEL??????????
      }
    } else {
        if (AnaerobicCO2 > 0) {
            anaerob = NewSOM.Row(i);
            for (j = 1; j <= NrReservoirs; j++) {
            // anaerobic CO2 from sulfate reduction and Fe/Mn reduction IN COMPLETELY ANAEROBIC LAYERS
                ka = KAnaerobic(j) * anaertemperaturefact;
                anaerobCO2(j) = anaerob(j) * ( 1.0 - (exp(- dt * ka)));
                NewSOM(i, j) -= anaerobCO2(j);
                CO2(i, j) += anaerobCO2(j);
                AnaerobSum(i) += anaerobCO2(j);
            }
            PeatLoss += CO2(i, 1);
            for (n = 6; n <= NrReservoirs; n++) {           // correct CO2 and NewSOM for carbon transferred to microbial biomass and humus
                for (m = 1; m <= 5; m++) {
                    transfer = CO2(i, m) * SplitRes(m, (n - 5));
                    NewSOM(i, n) = NewSOM(i, n) + transfer;
                    CO2(i, m) = CO2(i, m) - transfer;
                }
            }
        } else {
            for (j = 1; j <= NrReservoirs; j++) CO2(i, j) = 0.0;
        }
    }
  }
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

{
  int i;
  double f, c, ct;

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
  c = f * CO2FromMethane.Sum();                             // CO2 from methane
  ReservoirTime(StepNr, NrReservoirs + 2) = c;
  ct += c;
  ReservoirTime(StepNr, NrReservoirs + 3) = ct;             // total CO2 is last column of Reservoirtime
  for (i = 1; i <= NrLayers; i++)
  {
    LayerTime(StepNr, i + 1) = f * CO2.SumRow(i);
    LayerTime(StepNr, i + 1) += (f * CO2FromMethane(i));  // add CO2 evolved by methane oxidation
    // if (AnaerobicCO2 > 0) LayerAnaerobic(StepNr, i) = f * AnaerobSum(i); 
    LayerAnaerobic(StepNr, i) = f * AnaerobSum(i); // CO2 from anaerob reactions, including methane formation
  }
} // end  CollectCO2
