#ifndef TREE_UTILS_H
#define TREE_UTILS_H

#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"



/*** tree variables ***/
struct TreeVars
{
  std::vector<float>* genPho_pt;
  std::vector<float>* genPho_eta;
  std::vector<float>* genPho_phi;
  std::vector<int>* genPho_HardProcFinState;
  std::vector<int>* genPho_isPromptFinState;
  std::vector<float>* reso_pt;
  std::vector<float>* reso_eta;
  std::vector<float>* reso_phi;
  std::vector<float>* reso_energy;
  std::vector<int>* reso_charge;
  std::vector<int>* reso_pdgId;
  
  std::vector<int>* resoDau1_n;
  std::vector<std::vector<float> >* resoDau1_pt;
  std::vector<std::vector<float> >* resoDau1_eta;
  std::vector<std::vector<float> >* resoDau1_phi;
  std::vector<std::vector<float> >* resoDau1_energy;
  std::vector<std::vector<int> >* resoDau1_charge;
  std::vector<std::vector<int> >* resoDau1_pdgId;
  
  std::vector<std::vector<int> >* resoDau2_n;
  std::vector<std::vector<std::vector<float> > >* resoDau2_pt;
  std::vector<std::vector<std::vector<float> > >* resoDau2_eta;
  std::vector<std::vector<std::vector<float> > >* resoDau2_phi;
  std::vector<std::vector<std::vector<float> > >* resoDau2_energy;
  std::vector<std::vector<std::vector<int> > >* resoDau2_charge;
  std::vector<std::vector<std::vector<int> > >* resoDau2_pdgId;
  
  float trueNumInteractions;
  
  int vtxs_n;

  
  float rho_all;
  float rho_central;
  float rho_centralNeutral;
  float rho_centralChargedPU;
  
  std::vector<std::string>* trgs_name;
  std::vector<int>* trgs_pass;
  std::vector<int>* trgs_prescale;
  
  std::vector<float>* muons_pt;
  std::vector<float>* muons_ptErr;
  std::vector<float>* muons_eta;
  std::vector<float>* muons_phi;
  std::vector<float>* muons_energy;
  std::vector<int>* muons_charge;
  std::vector<float>* muons_dxy;
  std::vector<float>* muons_dxyErr;
  std::vector<float>* muons_dz;
  std::vector<float>* muons_dzErr;
  std::vector<float>* muons_isLoose;
  std::vector<float>* muons_isMedium;
  std::vector<float>* muons_isTight;
  std::vector<float>* muons_pfIsoChargedHadron;
  std::vector<float>* muons_pfIsoChargedParticle;
  std::vector<float>* muons_pfIsoNeutralHadron;
  std::vector<float>* muons_pfIsoPhoton;
  std::vector<float>* muons_pfIsoPU;
  std::vector<int>* muons_trackerLayersWithMeasurement;

  std::vector<float>* electrons_pt;
  std::vector<float>* electrons_eta;
  std::vector<float>* electrons_phi;
  std::vector<float>* electrons_energy;
  std::vector<float>* electrons_EnergyErrPostCorr;
  std::vector<float>* electrons_EnergyPostCorr;
  std::vector<int>* electrons_charge;
  std::vector<float>* electrons_dxy;
  std::vector<float>* electrons_dxyErr;
  std::vector<float>* electrons_dz;
  std::vector<float>* electrons_dzErr;
  std::vector<float>* electrons_full5x5_R9;
  std::vector<float>* electrons_full5x5_sieie;
  std::vector<float>* electrons_MVAID;

  std::vector<float>* photons_pt;
  std::vector<float>* photons_eta;
  std::vector<float>* photons_phi;
  std::vector<float>* photons_full5x5_R9;
  std::vector<float>* photons_full5x5_sieie;
  std::vector<float>* photons_EnergyPostCorr;


  std::vector<float>* photons_MVAID;
  
  std::vector<float>* jets_pt;
  std::vector<float>* jets_eta;
  std::vector<float>* jets_phi;
  std::vector<float>* jets_energy;
  std::vector<int>* jets_charge;
  std::vector<float>* jets_NHF;
  std::vector<float>* jets_NEMF;
  std::vector<float>* jets_CHF;
  std::vector<float>* jets_MUF;
  std::vector<float>* jets_CEMF;
  std::vector<int>* jets_CM;
  std::vector<int>* jets_NM;
  std::vector<std::vector<float> >* jets_bTag;
  
  /* std::vector<float>* jets_puppi_pt; */
  /* std::vector<float>* jets_puppi_eta; */
  /* std::vector<float>* jets_puppi_phi; */
  /* std::vector<float>* jets_puppi_energy; */
  /* std::vector<int>* jets_puppi_charge; */
  /* std::vector<float>* jets_puppi_NHF; */
  /* std::vector<float>* jets_puppi_NEMF; */
  /* std::vector<float>* jets_puppi_CHF; */
  /* std::vector<float>* jets_puppi_MUF; */
  /* std::vector<float>* jets_puppi_CEMF; */
  /* std::vector<int>* jets_puppi_CM; */
  /* std::vector<int>* jets_puppi_NM; */
  /* std::vector<std::vector<float> >* jets_puppi_bTag; */
  
  float met_pt;
  float met_phi;
  float met_sig;
  
  /* float met_puppi_pt; */
  /* float met_puppi_phi; */
  /* float met_puppi_sig; */
};

void InitTreeVars(TChain* chain_reco, TChain* chain_gen, TreeVars& treeVars);

#endif
