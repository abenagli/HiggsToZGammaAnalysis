#include "interface/TreeUtils.h"



void InitTreeVars(TChain* chain_reco, TChain* chain_gen, TreeVars& treeVars)
{
  treeVars.reso_pt = new std::vector<float>;
  treeVars.reso_eta = new std::vector<float>;
  treeVars.reso_phi = new std::vector<float>;
  treeVars.reso_energy = new std::vector<float>;
  treeVars.reso_charge = new std::vector<int>;
  treeVars.reso_pdgId = new std::vector<int>;
  
  treeVars.resoDau1_n = new std::vector<int>;
  treeVars.resoDau1_pt = new std::vector<std::vector<float> >;
  treeVars.resoDau1_eta = new std::vector<std::vector<float> >;
  treeVars.resoDau1_phi = new std::vector<std::vector<float> >;
  treeVars.resoDau1_energy = new std::vector<std::vector<float> >;
  treeVars.resoDau1_charge = new std::vector<std::vector<int> >;
  treeVars.resoDau1_pdgId = new std::vector<std::vector<int> >;
  
  treeVars.resoDau2_n = new std::vector<std::vector<int> >;
  treeVars.resoDau2_pt = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_eta = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_phi = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_energy = new std::vector<std::vector<std::vector<float> > >;
  treeVars.resoDau2_charge = new std::vector<std::vector<std::vector<int> > >;  
  treeVars.resoDau2_pdgId = new std::vector<std::vector<std::vector<int> > >;

  treeVars.trgs_name = new std::vector<std::string>;
  treeVars.trgs_pass = new std::vector<int>;
  treeVars.trgs_prescale = new std::vector<int>;
  
  treeVars.muons_pt = new std::vector<float>;
  treeVars.muons_eta = new std::vector<float>;
  treeVars.muons_phi = new std::vector<float>;
  treeVars.muons_energy = new std::vector<float>;
  treeVars.muons_charge = new std::vector<int>;
  treeVars.muons_dxy = new std::vector<float>;
  treeVars.muons_dxyErr = new std::vector<float>;
  treeVars.muons_dz = new std::vector<float>;
  treeVars.muons_dzErr = new std::vector<float>;
  treeVars.muons_isLoose = new std::vector<float>;
  treeVars.muons_isMedium = new std::vector<float>;
  treeVars.muons_isTight = new std::vector<float>;
  treeVars.muons_pfIsoChargedHadron = new std::vector<float>;
  treeVars.muons_pfIsoChargedParticle = new std::vector<float>;
  treeVars.muons_pfIsoNeutralHadron = new std::vector<float>;
  treeVars.muons_pfIsoPhoton = new std::vector<float>;
  treeVars.muons_pfIsoPU = new std::vector<float>;
  treeVars.muons_trackerLayersWithMeasurement = new std::vector<int>;

  treeVars.electrons_pt = new std::vector<float>;
  treeVars.electrons_eta = new std::vector<float>;
  treeVars.electrons_phi = new std::vector<float>;
  treeVars.electrons_EnergyPostCorr = new std::vector<float>;
  treeVars.electrons_charge = new std::vector<int>;
  treeVars.electrons_MVAID = new std::vector<float>;

  treeVars.photons_pt = new std::vector<float>;
  treeVars.photons_eta = new std::vector<float>;
  treeVars.photons_phi = new std::vector<float>;
  treeVars.photons_EnergyPostCorr = new std::vector<float>;
  treeVars.photons_MVAID = new std::vector<float>;
  
  treeVars.jets_pt = new std::vector<float>;
  treeVars.jets_eta = new std::vector<float>;
  treeVars.jets_phi = new std::vector<float>;
  treeVars.jets_energy = new std::vector<float>;
  treeVars.jets_charge = new std::vector<int>;
  treeVars.jets_NHF = new std::vector<float>;
  treeVars.jets_NEMF = new std::vector<float>;
  treeVars.jets_CHF = new std::vector<float>;
  treeVars.jets_MUF = new std::vector<float>;
  treeVars.jets_CEMF = new std::vector<float>;
  treeVars.jets_CM = new std::vector<int>;
  treeVars.jets_NM = new std::vector<int>;
  treeVars.jets_bTag = new std::vector<std::vector<float> >;
  
  // treeVars.jets_puppi_pt = new std::vector<float>;
  // treeVars.jets_puppi_eta = new std::vector<float>;
  // treeVars.jets_puppi_phi = new std::vector<float>;
  // treeVars.jets_puppi_energy = new std::vector<float>;
  // treeVars.jets_puppi_charge = new std::vector<int>;
  // treeVars.jets_puppi_bTag = new std::vector<std::vector<float> >;
  // treeVars.jets_puppi_NHF = new std::vector<float>;
  // treeVars.jets_puppi_NEMF = new std::vector<float>;
  // treeVars.jets_puppi_CHF = new std::vector<float>;
  // treeVars.jets_puppi_MUF = new std::vector<float>;
  // treeVars.jets_puppi_CEMF = new std::vector<float>;
  // treeVars.jets_puppi_CM = new std::vector<int>;
  // treeVars.jets_puppi_NM = new std::vector<int>;
  // treeVars.jets_puppi_bTag = new std::vector<std::vector<float> >;
  
  //tree -> SetBranchStatus("*",0);
  
  if( chain_gen != NULL )
  {
    chain_gen -> SetBranchStatus("reso_pt",    1); chain_gen -> SetBranchAddress("reso_pt",    &treeVars.reso_pt);
    chain_gen -> SetBranchStatus("reso_eta",   1); chain_gen -> SetBranchAddress("reso_eta",   &treeVars.reso_eta);
    chain_gen -> SetBranchStatus("reso_phi",   1); chain_gen -> SetBranchAddress("reso_phi",   &treeVars.reso_phi);
    chain_gen -> SetBranchStatus("reso_energy",1); chain_gen -> SetBranchAddress("reso_energy",&treeVars.reso_energy);
    chain_gen -> SetBranchStatus("reso_charge",1); chain_gen -> SetBranchAddress("reso_charge",&treeVars.reso_charge);
    chain_gen -> SetBranchStatus("reso_pdgId", 1); chain_gen -> SetBranchAddress("reso_pdgId", &treeVars.reso_pdgId);
    
    chain_gen -> SetBranchStatus("resoDau1_n",     1); chain_gen -> SetBranchAddress("resoDau1_n",     &treeVars.resoDau1_n);
    chain_gen -> SetBranchStatus("resoDau1_pt",    1); chain_gen -> SetBranchAddress("resoDau1_pt",    &treeVars.resoDau1_pt);
    chain_gen -> SetBranchStatus("resoDau1_eta",   1); chain_gen -> SetBranchAddress("resoDau1_eta",   &treeVars.resoDau1_eta);
    chain_gen -> SetBranchStatus("resoDau1_phi",   1); chain_gen -> SetBranchAddress("resoDau1_phi",   &treeVars.resoDau1_phi);
    chain_gen -> SetBranchStatus("resoDau1_energy",1); chain_gen -> SetBranchAddress("resoDau1_energy",&treeVars.resoDau1_energy);
    chain_gen -> SetBranchStatus("resoDau1_charge",1); chain_gen -> SetBranchAddress("resoDau1_charge",&treeVars.resoDau1_charge);
    chain_gen -> SetBranchStatus("resoDau1_pdgId", 1); chain_gen -> SetBranchAddress("resoDau1_pdgId", &treeVars.resoDau1_pdgId);
    
    chain_gen -> SetBranchStatus("resoDau2_n",     1); chain_gen -> SetBranchAddress("resoDau2_n",     &treeVars.resoDau2_n);
    chain_gen -> SetBranchStatus("resoDau2_pt",    1); chain_gen -> SetBranchAddress("resoDau2_pt",    &treeVars.resoDau2_pt);
    chain_gen -> SetBranchStatus("resoDau2_eta",   1); chain_gen -> SetBranchAddress("resoDau2_eta",   &treeVars.resoDau2_eta);
    chain_gen -> SetBranchStatus("resoDau2_phi",   1); chain_gen -> SetBranchAddress("resoDau2_phi",   &treeVars.resoDau2_phi);
    chain_gen -> SetBranchStatus("resoDau2_energy",1); chain_gen -> SetBranchAddress("resoDau2_energy",&treeVars.resoDau2_energy);
    chain_gen -> SetBranchStatus("resoDau2_charge",1); chain_gen -> SetBranchAddress("resoDau2_charge",&treeVars.resoDau2_charge);
    chain_gen -> SetBranchStatus("resoDau2_pdgId", 1); chain_gen -> SetBranchAddress("resoDau2_pdgId", &treeVars.resoDau2_pdgId);
    
    chain_gen -> SetBranchStatus("trueNumInteractions", 1); chain_gen -> SetBranchAddress("trueNumInteractions", &treeVars.trueNumInteractions);
  }

  if( chain_reco != NULL )
  {
    chain_reco -> SetBranchStatus("vtxs_n",1); chain_reco -> SetBranchAddress("vtxs_n",&treeVars.vtxs_n);
    
    chain_reco -> SetBranchStatus("rho_all",             1); chain_reco -> SetBranchAddress("rho_all",             &treeVars.rho_all);
    chain_reco -> SetBranchStatus("rho_central",         1); chain_reco -> SetBranchAddress("rho_central",         &treeVars.rho_central);
    chain_reco -> SetBranchStatus("rho_centralNeutral",  1); chain_reco -> SetBranchAddress("rho_centralNeutral",  &treeVars.rho_centralNeutral);
    chain_reco -> SetBranchStatus("rho_centralChargedPU",1); chain_reco -> SetBranchAddress("rho_centralChargedPU",&treeVars.rho_centralChargedPU);
    
    chain_reco -> SetBranchStatus("trgs_name",    1); chain_reco -> SetBranchAddress("trgs_name",    &treeVars.trgs_name);
    chain_reco -> SetBranchStatus("trgs_pass",    1); chain_reco -> SetBranchAddress("trgs_pass",    &treeVars.trgs_pass);
    chain_reco -> SetBranchStatus("trgs_prescale",1); chain_reco -> SetBranchAddress("trgs_prescale",&treeVars.trgs_prescale);
    
    chain_reco -> SetBranchStatus("muons_pt",      1); chain_reco -> SetBranchAddress("muons_pt",      &treeVars.muons_pt);
    chain_reco -> SetBranchStatus("muons_eta",     1); chain_reco -> SetBranchAddress("muons_eta",     &treeVars.muons_eta);
    chain_reco -> SetBranchStatus("muons_phi",     1); chain_reco -> SetBranchAddress("muons_phi",     &treeVars.muons_phi);
    chain_reco -> SetBranchStatus("muons_energy",  1); chain_reco -> SetBranchAddress("muons_energy",  &treeVars.muons_energy);
    chain_reco -> SetBranchStatus("muons_charge",  1); chain_reco -> SetBranchAddress("muons_charge",  &treeVars.muons_charge);
    chain_reco -> SetBranchStatus("muons_dxy",     1); chain_reco -> SetBranchAddress("muons_dxy",     &treeVars.muons_dxy);
    chain_reco -> SetBranchStatus("muons_dxyErr",  1); chain_reco -> SetBranchAddress("muons_dxyErr",  &treeVars.muons_dxyErr);
    chain_reco -> SetBranchStatus("muons_dz",      1); chain_reco -> SetBranchAddress("muons_dz",      &treeVars.muons_dz);
    chain_reco -> SetBranchStatus("muons_dzErr",   1); chain_reco -> SetBranchAddress("muons_dzErr",   &treeVars.muons_dzErr);
    chain_reco -> SetBranchStatus("muons_isLoose", 1); chain_reco -> SetBranchAddress("muons_isLoose", &treeVars.muons_isLoose);
    chain_reco -> SetBranchStatus("muons_isMedium",1); chain_reco -> SetBranchAddress("muons_isMedium",&treeVars.muons_isMedium);
    chain_reco -> SetBranchStatus("muons_isTight", 1); chain_reco -> SetBranchAddress("muons_isTight", &treeVars.muons_isTight);
    chain_reco -> SetBranchStatus("muons_pfIsoChargedHadron",  1); chain_reco -> SetBranchAddress("muons_pfIsoChargedHadron",  &treeVars.muons_pfIsoChargedHadron);
    chain_reco -> SetBranchStatus("muons_pfIsoChargedParticle",1); chain_reco -> SetBranchAddress("muons_pfIsoChargedParticle",&treeVars.muons_pfIsoChargedParticle);
    chain_reco -> SetBranchStatus("muons_pfIsoNeutralHadron",  1); chain_reco -> SetBranchAddress("muons_pfIsoNeutralHadron",  &treeVars.muons_pfIsoNeutralHadron);
    chain_reco -> SetBranchStatus("muons_pfIsoPhoton",         1); chain_reco -> SetBranchAddress("muons_pfIsoPhoton",         &treeVars.muons_pfIsoPhoton);
    chain_reco -> SetBranchStatus("muons_pfIsoPU",             1); chain_reco -> SetBranchAddress("muons_pfIsoPU",             &treeVars.muons_pfIsoPU);
    chain_reco -> SetBranchStatus("muons_trackerLayersWithMeasurement",1); chain_reco -> SetBranchAddress("muons_trackerLayersWithMeasurement",&treeVars.muons_trackerLayersWithMeasurement);

    chain_reco -> SetBranchStatus("electrons_pt",      1); chain_reco -> SetBranchAddress("electrons_pt",      &treeVars.electrons_pt);
    chain_reco -> SetBranchStatus("electrons_eta",     1); chain_reco -> SetBranchAddress("electrons_eta",     &treeVars.electrons_eta);
    chain_reco -> SetBranchStatus("electrons_phi",     1); chain_reco -> SetBranchAddress("electrons_phi",     &treeVars.electrons_phi);
    chain_reco -> SetBranchStatus("electrons_EnergyPostCorr",  1); chain_reco -> SetBranchAddress("electrons_EnergyPostCorr",  &treeVars.electrons_EnergyPostCorr);
    chain_reco -> SetBranchStatus("electrons_charge",  1); chain_reco -> SetBranchAddress("electrons_charge",  &treeVars.electrons_charge);
    chain_reco -> SetBranchStatus("electrons_MVAID", 1); chain_reco -> SetBranchAddress("electrons_MVAID", &treeVars.electrons_MVAID);

    chain_reco -> SetBranchStatus("photons_pt",      1); chain_reco -> SetBranchAddress("photons_pt",      &treeVars.photons_pt);
    chain_reco -> SetBranchStatus("photons_eta",     1); chain_reco -> SetBranchAddress("photons_eta",     &treeVars.photons_eta);
    chain_reco -> SetBranchStatus("photons_phi",     1); chain_reco -> SetBranchAddress("photons_phi",     &treeVars.photons_phi);
    chain_reco -> SetBranchStatus("photons_EnergyPostCorr",  1); chain_reco -> SetBranchAddress("photons_EnergyPostCorr",  &treeVars.photons_EnergyPostCorr);
    chain_reco -> SetBranchStatus("photons_MVAID", 1); chain_reco -> SetBranchAddress("photons_MVAID", &treeVars.photons_MVAID);

    
    chain_reco -> SetBranchStatus("jets_pt",      1); chain_reco -> SetBranchAddress("jets_pt",      &treeVars.jets_pt);
    chain_reco -> SetBranchStatus("jets_eta",     1); chain_reco -> SetBranchAddress("jets_eta",     &treeVars.jets_eta);
    chain_reco -> SetBranchStatus("jets_phi",     1); chain_reco -> SetBranchAddress("jets_phi",     &treeVars.jets_phi);
    chain_reco -> SetBranchStatus("jets_energy",  1); chain_reco -> SetBranchAddress("jets_energy",  &treeVars.jets_energy);
    chain_reco -> SetBranchStatus("jets_charge",  1); chain_reco -> SetBranchAddress("jets_charge",  &treeVars.jets_charge);
    chain_reco -> SetBranchStatus("jets_NHF",     1); chain_reco -> SetBranchAddress("jets_NHF",     &treeVars.jets_NHF);
    chain_reco -> SetBranchStatus("jets_NEMF",    1); chain_reco -> SetBranchAddress("jets_NEMF",    &treeVars.jets_NEMF);
    chain_reco -> SetBranchStatus("jets_CHF",     1); chain_reco -> SetBranchAddress("jets_CHF",     &treeVars.jets_CHF);
    chain_reco -> SetBranchStatus("jets_MUF",     1); chain_reco -> SetBranchAddress("jets_MUF",     &treeVars.jets_MUF);
    chain_reco -> SetBranchStatus("jets_CEMF",    1); chain_reco -> SetBranchAddress("jets_CEMF",    &treeVars.jets_CEMF);
    chain_reco -> SetBranchStatus("jets_CM",      1); chain_reco -> SetBranchAddress("jets_CM",      &treeVars.jets_CM);
    chain_reco -> SetBranchStatus("jets_NM",      1); chain_reco -> SetBranchAddress("jets_NM",      &treeVars.jets_NM);
    chain_reco -> SetBranchStatus("jets_bTag",    1); chain_reco -> SetBranchAddress("jets_bTag",    &treeVars.jets_bTag);
    
    // chain_reco -> SetBranchStatus("jets_puppi_pt",      1); chain_reco -> SetBranchAddress("jets_puppi_pt",      &treeVars.jets_puppi_pt);
    // chain_reco -> SetBranchStatus("jets_puppi_eta",     1); chain_reco -> SetBranchAddress("jets_puppi_eta",     &treeVars.jets_puppi_eta);
    // chain_reco -> SetBranchStatus("jets_puppi_phi",     1); chain_reco -> SetBranchAddress("jets_puppi_phi",     &treeVars.jets_puppi_phi);
    // chain_reco -> SetBranchStatus("jets_puppi_energy",  1); chain_reco -> SetBranchAddress("jets_puppi_energy",  &treeVars.jets_puppi_energy);
    // chain_reco -> SetBranchStatus("jets_puppi_charge",  1); chain_reco -> SetBranchAddress("jets_puppi_charge",  &treeVars.jets_puppi_charge);
    // chain_reco -> SetBranchStatus("jets_puppi_NHF",     1); chain_reco -> SetBranchAddress("jets_puppi_NHF",     &treeVars.jets_puppi_NHF);
    // chain_reco -> SetBranchStatus("jets_puppi_NEMF",    1); chain_reco -> SetBranchAddress("jets_puppi_NEMF",    &treeVars.jets_puppi_NEMF);
    // chain_reco -> SetBranchStatus("jets_puppi_CHF",     1); chain_reco -> SetBranchAddress("jets_puppi_CHF",     &treeVars.jets_puppi_CHF);
    // chain_reco -> SetBranchStatus("jets_puppi_MUF",     1); chain_reco -> SetBranchAddress("jets_puppi_MUF",     &treeVars.jets_puppi_MUF);
    // chain_reco -> SetBranchStatus("jets_puppi_CEMF",    1); chain_reco -> SetBranchAddress("jets_puppi_CEMF",    &treeVars.jets_puppi_CEMF);
    // chain_reco -> SetBranchStatus("jets_puppi_CM",      1); chain_reco -> SetBranchAddress("jets_puppi_CM",      &treeVars.jets_puppi_CM);
    // chain_reco -> SetBranchStatus("jets_puppi_NM",      1); chain_reco -> SetBranchAddress("jets_puppi_NM",      &treeVars.jets_puppi_NM);
    // chain_reco -> SetBranchStatus("jets_puppi_bTag",    1); chain_reco -> SetBranchAddress("jets_puppi_bTag",    &treeVars.jets_puppi_bTag);
    
    chain_reco -> SetBranchStatus("met_pt", 1); chain_reco -> SetBranchAddress("met_pt", &treeVars.met_pt);
    chain_reco -> SetBranchStatus("met_phi",1); chain_reco -> SetBranchAddress("met_phi",&treeVars.met_phi);
    chain_reco -> SetBranchStatus("met_sig",1); chain_reco -> SetBranchAddress("met_sig",&treeVars.met_sig);
    
    // chain_reco -> SetBranchStatus("met_puppi_pt", 1); chain_reco -> SetBranchAddress("met_puppi_pt", &treeVars.met_puppi_pt);
    // chain_reco -> SetBranchStatus("met_puppi_phi",1); chain_reco -> SetBranchAddress("met_puppi_phi",&treeVars.met_puppi_phi);
    // chain_reco -> SetBranchStatus("met_puppi_sig",1); chain_reco -> SetBranchAddress("met_puppi_sig",&treeVars.met_puppi_sig);
  }
}
