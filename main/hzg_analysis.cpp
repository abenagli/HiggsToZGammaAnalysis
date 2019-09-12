#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/PUReweighting.h"
#include "interface/OutTree.h"
#include "interface/StdHistoSet.h"
#include "src/RoccoR.cc"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TRandom3.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> hzg_analysis.cpp::usage:   " << argv[0] << " configFileName   [default=0/debug=1]" << std::endl;
    return -1;
  }
  TDirectory* baseDir = gDirectory;
  TRandom3 r;
  
  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);


  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  //--- open files and get trees
  int isData = opts.GetOpt<int>("Input.isData");
  int isSignal = opts.GetOpt<int>("Input.isSignal");
  std::string inputFileList = opts.GetOpt<std::string>("Input.inputFileList");
  // std::string inputFilePU = opts.GetOpt<std::string>("Input.inputFilePU");
  
  TChain* chain_reco = new TChain("DumpReco/reco_tree","DumpReco/reco_tree");
  TChain* chain_gen = NULL; if( !isData ) chain_gen = new TChain("DumpGenParticles/gen_tree","DumpGenParticles/gen_tree");
  
  std::ifstream list(inputFileList.c_str(),std::ios::in);
  std::string fileName;
  while(1)
  {
    getline(list,fileName,'\n');
    if( !list.good() ) break;
    if( fileName.at(0) == '#' ) continue;
    std::cout << "Reading file " << fileName << std::endl;
    
    chain_reco -> Add(fileName.c_str());
    if( !isData ) chain_gen -> Add(fileName.c_str());
  }

  TreeVars treeVars;
  InitTreeVars(chain_reco,chain_gen,treeVars);
  
  int nEntries_recoTree = chain_reco -> GetEntries();
  int nEntries_genTree = 0; if( !isData ) nEntries_genTree = chain_gen -> GetEntries();
  std::cout << "Read " << nEntries_recoTree << " in tree reco_tree" << std::endl;
  std::cout << "Read " << nEntries_genTree  << " in tree gen_tree"  << std::endl;
  
  
  //--- get other variables
  int printGenEvent = opts.GetOpt<int>("Input.printGenEvent");
  int printRecoGenMatchEvent = opts.GetOpt<int>("Input.printRecoGenMatchEvent");
  int printRecoEvent = opts.GetOpt<int>("Input.printRecoEvent");
  
  // float gen_mu_ptMin = opts.GetOpt<float>("Cuts.gen_mu_ptMin");
  // float gen_mu_etaMax = opts.GetOpt<float>("Cuts.gen_mu_etaMax");
  // int genCut = opts.GetOpt<int>("Cuts.genCut");
  
  int genMatch = opts.GetOpt<int>("Cuts.genMatch");  
  // int genMatchCut = opts.GetOpt<int>("Cuts.genMatchCut");
  
  int reco_HLTCut = opts.GetOpt<int>("Cuts.reco_HLTCut");
  std::vector<std::string> reco_HLT_paths = opts.GetOpt<std::vector<std::string> >("Cuts.reco_HLT_paths");
  int reco_HLT_OR = opts.GetOpt<int>("Cuts.reco_HLT_OR");
  
  int doMuon = opts.GetOpt<int>("Cuts.doMuon");
  int doEle = opts.GetOpt<int>("Cuts.doEle");
  int doFilter = opts.GetOpt<int>("Cuts.doFilter");
  int doTrigEff = opts.GetOpt<int>("Cuts.doTrigEff");

  float reco_ptMin1_mu = opts.GetOpt<float>("Cuts.reco_ptMin1_mu");
  float reco_ptMin2_mu = opts.GetOpt<float>("Cuts.reco_ptMin2_mu");
  float reco_eta_mu = opts.GetOpt<float>("Cuts.reco_eta_mu");

  float reco_ptMin1_ele = opts.GetOpt<float>("Cuts.reco_ptMin1_ele");
  float reco_ptMin2_ele = opts.GetOpt<float>("Cuts.reco_ptMin2_ele");
  float reco_eta_ele = opts.GetOpt<float>("Cuts.reco_eta_ele");
  float reco_ID_ele = opts.GetOpt<float>("Cuts.reco_ID_ele");

  float reco_ptMin_gamma = opts.GetOpt<float>("Cuts.reco_ptMin_gamma");
  float reco_ID_gammaB = opts.GetOpt<float>("Cuts.reco_ID_gammaB");
  float reco_ID_gammaE = opts.GetOpt<float>("Cuts.reco_ID_gammaE");
  float reco_eta_gamma = opts.GetOpt<float>("Cuts.reco_eta_gamma");
  float reco_eta1_gamma = opts.GetOpt<float>("Cuts.reco_eta1_gamma");
  float reco_eta2_gamma = opts.GetOpt<float>("Cuts.reco_eta2_gamma");
  float reco_pToM_gamma = opts.GetOpt<float>("Cuts.reco_pToM_gamma");

  float reco_DR_lep = opts.GetOpt<float>("Cuts.reco_DR_lep");
  float reco_sumM = opts.GetOpt<float>("Cuts.reco_sumM");
  
  if (doMuon) std::cout << "Muon tag " << std::endl;
  if (doEle)  {std::cout << "Electron tag " << std::endl;}
  /*
  //--- get scale factors files
  std::vector<std::string> SFparams_HLT = opts.GetOpt<std::vector<std::string> >("ScaleFactors.HLT");
  std::vector<TFile*> inFiles_SF_HLT;
  std::vector<TH2F*> h2_SF_HLT;
  std::vector<float> frac_SF_HLT;
  frac_SF_HLT.push_back(0.);
  for(int ii = 0; ii < int(SFparams_HLT.size()/3); ++ii)
  {
    inFiles_SF_HLT.push_back( TFile::Open(SFparams_HLT.at(0+ii*3).c_str(),"READ") );
    h2_SF_HLT.push_back( (TH2F*)(inFiles_SF_HLT[ii]->Get(SFparams_HLT.at(1+ii*3).c_str())) );
    frac_SF_HLT.push_back( atof(SFparams_HLT.at(2+ii*3).c_str())+frac_SF_HLT.at(frac_SF_HLT.size()-1) );
  }
  
  std::vector<std::string> SFparams_ID = opts.GetOpt<std::vector<std::string> >("ScaleFactors.ID");
  std::vector<TFile*> inFiles_SF_ID;
  std::vector<TH2F*> h2_SF_ID;
  std::vector<float> frac_SF_ID;
  frac_SF_ID.push_back(0.);
  for(int ii = 0; ii < int(SFparams_ID.size()/3); ++ii)
  {
    inFiles_SF_ID.push_back( TFile::Open(SFparams_ID.at(0+ii*3).c_str(),"READ") );
    h2_SF_ID.push_back( (TH2F*)(inFiles_SF_ID[ii]->Get(SFparams_ID.at(1+ii*3).c_str())) );
    frac_SF_ID.push_back( atof(SFparams_ID.at(2+ii*3).c_str())+frac_SF_ID.at(frac_SF_ID.size()-1) );
  }
  
  std::vector<std::string> SFparams_ISO = opts.GetOpt<std::vector<std::string> >("ScaleFactors.ISO");
  std::vector<TFile*> inFiles_SF_ISO;
  std::vector<TH2F*> h2_SF_ISO;
  std::vector<float> frac_SF_ISO;
  frac_SF_ISO.push_back(0.);
  for(int ii = 0; ii < int(SFparams_ISO.size()/3); ++ii)
  {
    inFiles_SF_ISO.push_back( TFile::Open(SFparams_ISO.at(0+ii*3).c_str(),"READ") );
    h2_SF_ISO.push_back( (TH2F*)(inFiles_SF_ISO[ii]->Get(SFparams_ISO.at(1+ii*3).c_str())) );
    frac_SF_ISO.push_back( atof(SFparams_ISO.at(2+ii*3).c_str())+frac_SF_ISO.at(frac_SF_ISO.size()-1) );
  }
  
  std::vector<std::string> kFactors_pt = opts.GetOpt<std::vector<std::string> >("Weights.kFactor_pt");
  TFile* inFile_kFactor_pt = TFile::Open(kFactors_pt.at(0).c_str(),"READ");
  TF1* f_kFactor_pt = (TF1*)( inFile_kFactor_pt->Get(kFactors_pt.at(1).c_str()) );
  */
  
  
  //--- define histograms
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  TFile* outFile = TFile::Open(Form("%s.root",outputFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  OutTree outTree = OutTree("outTree","Final tree for H > Z gamma studies");
  
  /*
  StdHistoSet* genHistos = isData? NULL : new StdHistoSet("gen",outFile);;
  StdHistoSet* genCutHistos = isData? NULL : new StdHistoSet("genCut",outFile);;
  StdHistoSet* recoGenMatchHistos = isData? NULL : new StdHistoSet("recoGenMatch",outFile);;
  StdHistoSet* recoHistos = new StdHistoSet("reco",outFile);
  */
  
  TH1F* h1_nEvents = new TH1F("h1_nEvents",";;entries",20,0.,20.);
  
  TH1F* h_ptReso_mu = new TH1F("h_ptReso_mu","",1000,-10.,10.);
  TH1F* h_ptReso_ele = new TH1F("h_ptReso_ele","",1000,-10.,10.);
  
  std::vector <double> RangeEtaele;
  RangeEtaele.push_back(0);
  RangeEtaele.push_back(1);
  RangeEtaele.push_back(1.5);
  RangeEtaele.push_back(2.5);
  
  std::vector <double> RangePtele;
  RangePtele.push_back(25);
  RangePtele.push_back(50);
  RangePtele.push_back(100);
  RangePtele.push_back(500);


  std::vector <double> RangeEtamu;
  RangeEtamu.push_back(0);
  RangeEtamu.push_back(0.9);
  RangeEtamu.push_back(1.8);
  RangeEtamu.push_back(2.4);


  std::vector <double> RangePtmu;
  RangePtmu.push_back(20);
  RangePtmu.push_back(50);
  RangePtmu.push_back(100);
  RangePtmu.push_back(500);

  
  double* xAxisPt_ele = new double[RangePtele.size()];
  for(unsigned int kk = 0; kk < RangePtele.size(); ++kk)
    xAxisPt_ele[kk] =  RangePtele.at(kk);

  double* xAxisPt_mu = new double[RangePtmu.size()];
  for(unsigned int kk = 0; kk < RangePtmu.size(); ++kk)
    xAxisPt_mu[kk] =  RangePtmu.at(kk);

  double* xAxisEta_ele = new double[RangeEtaele.size()];
  for(unsigned int kk = 0; kk < RangeEtaele.size(); ++kk)
    xAxisEta_ele[kk] =  RangeEtaele.at(kk);

  double* xAxisEta_mu = new double[RangeEtamu.size()];
  for(unsigned int kk = 0; kk < RangeEtamu.size(); ++kk)
    xAxisEta_mu[kk] =  RangeEtamu.at(kk);
  
  TH1F* h_occupancyPt_ele = new TH1F("h_occupancyPt_ele","",RangePtele.size()-1,xAxisPt_ele);
  TH1F* h_occupancyPt_mu = new TH1F("h_occupancyPt_mu","",RangePtmu.size()-1,xAxisPt_mu);
  TH1F* h_occupancyEta_ele = new TH1F("h_occupancyEta_ele","",RangeEtaele.size()-1,xAxisEta_ele);
  TH1F* h_occupancyEta_mu = new TH1F("h_occupancyEta_mu","",RangeEtamu.size()-1,xAxisEta_mu);
  
  std::map<std::string,TH1F*> histosele;
  
  for(int i=0; i<RangeEtaele.size()-1; ++i)
  {
    for(int j=0; j<RangePtele.size()-1; ++j)
    {
      float etaminele=RangeEtaele.at(i);
      float etamaxele=RangeEtaele.at(i+1);
      
      float ptminele=RangePtele.at(j);
      float ptmaxele=RangePtele.at(j+1);
      
      std::string labelele(Form("Res_eta_%.1f-%.1f_pt%.0f-%.0f_ele",etaminele,etamaxele,ptminele,ptmaxele));
      histosele[labelele] = new TH1F(labelele.c_str(),"",1000,-10.,10.);
    }
  }
  
  std::map<std::string, TH1F*> histosmu;
  
  for(int i=0; i<RangeEtamu.size()-1; ++i)
  {
    for(int j=0; j<RangePtmu.size()-1; ++j)
    {
      float etaminmu=RangeEtamu.at(i);
      float etamaxmu=RangeEtamu.at(i+1);
      
      float ptminmu=RangePtmu.at(j);
      float ptmaxmu=RangePtmu.at(j+1);
      
      std::string labelmu(Form("Res_eta_%.1f-%.1f_pt%.0f-%.0f_mu",etaminmu,etamaxmu,ptminmu,ptmaxmu));
      histosmu[labelmu] = new TH1F(labelmu.c_str(),"",1000,-10., 10.);
    }
  }
  
  
  
  /*
    TH2F* h2_mu_pt2_vs_pt1_total = new TH2F("h2_mu_pt2_vs_pt1_total","",100,-0.5,99.5,100,-0.5,99.5);
    TH2F* h2_mu_pt2_vs_pt1_pass  = new TH2F("h2_mu_pt2_vs_pt1_pass", "",100,-0.5,99.5,100,-0.5,99.5);
    TH2F* h2_mu_cuts_pt2_vs_pt1_total = new TH2F("h2_mu_cuts_pt2_vs_pt1_total","",100,-0.5,99.5,100,-0.5,99.5);
    TH2F* h2_mu_cuts_pt2_vs_pt1_pass  = new TH2F("h2_mu_cuts_pt2_vs_pt1_pass", "",100,-0.5,99.5,100,-0.5,99.5);
    
    TEfficiency* p1_eff_mu_pt = new TEfficiency("p1_eff_mu_pt",";muon p_{T} (GeV);#epsilon",200,0.,500.);
    TEfficiency* p1_eff_mu_eta = new TEfficiency("p1_eff_mu_eta",";muon #eta;#epsilon",200,0.,10.);
    
    TH1F* h1_reco_mu_gen_mu_DR = new TH1F("h1_reco_mu_gen_mu_DR","",10000,0.,10.);
    TH1F* h1_reco_mu_gen_mu_ptRatio = new TH1F("h1_reco_mu_gen_mu_ptRatio","",10000,0.,10.);
  */
  
  std::vector<std::pair<std::string,int> > vec_triggerPass;
  
  int nTriggerEvents_noCuts = 0;
  int nTriggerEvents_genCuts = 0;
  int nTriggerEvents_recoLepSelection = 0;
  int nTriggerEvents_AllSelection = 0;
  std::vector<std::pair<std::string,int> > vec_triggerPass_noCuts;
  std::vector<std::pair<std::string,int> > vec_triggerPass_genCuts;
  std::vector<std::pair<std::string,int> > vec_triggerPass_recoLepSelection;
  std::vector<std::pair<std::string,int> > vec_triggerPass_AllSelection;
  TH1F* h1_triggerEff_noCuts = new TH1F("h1_triggerEff_noCuts","",15,0.,15.);
  TH1F* h1_triggerEff_genCuts = new TH1F("h1_triggerEff_genCuts","",15,0.,15.);
  TH1F* h1_triggerEff_recoLepSelection = new TH1F("h1_triggerEff_recoLepSelection","",15,0.,15.);
  TH1F* h1_triggerEff_AllSelection = new TH1F("h1_triggerEff_AllSelection","",15,0.,15.);
  
  
  
  /*
  //----------------------
  //--- pileup reweighting
  std::vector<float> frac_PU;
  std::vector<PUReweighting*> puReweighting;
  if( !isData )
  {
    TFile* inputFile_pileup = TFile::Open(inputFilePU.c_str(),"READ");
    TH1F* pileup_mc = (TH1F*)( inputFile_pileup->Get("pileup_mc") );
    
    std::vector<std::string> pileupFileName_data = opts.GetOpt<std::vector<std::string> >("Weights.pileupFileName_data");
    frac_PU.push_back(0.);
    for(int jj = 0; jj < int(pileupFileName_data.size()/2); ++jj)
    {
      frac_PU.push_back( atof((pileupFileName_data.at(1+jj*2)).c_str()) );
      
      TFile* pileupFile_data = TFile::Open(pileupFileName_data.at(0+jj*2).c_str(),"READ");
      TH1F* pileup_data = (TH1F*)( pileupFile_data->Get("pileup") );
      pileup_data -> SetName(Form("pileup_data_%d",jj));
      
      puReweighting.push_back( new PUReweighting(pileup_data,pileup_mc) );
      
      TH1F* pileup_weights = puReweighting.at(jj) -> GetWeightsHistogram();
      pileup_weights -> SetName(Form("pileup_weights_%d",jj));
      
      outFile -> cd();
      pileup_data -> Write();
      pileup_weights -> Write();
    }
    
    pileup_mc -> Write();
  }
  */
  
  
  
  //----------------
  //--- event weight
  float mcWeight = 1.;
  if( !isData )
  {
    float lumi = opts.GetOpt<float>("Weights.lumi");
    float xsec = opts.GetOpt<float>("Weights.xsec");
    float nevt = opts.GetOpt<int>("Weights.nevt");
    
    mcWeight = 1. * lumi * xsec / nevt;
    std::cout << "mcWeight:  "<< mcWeight << std::endl;
  }
  
  
  
  
  //------------------------------
  //--- muon Rochester corrections
  //RoccoR rc("data/rcdata.2016.v3");
  
  int no_sel=0; 
  int lep_sel=0; 
  int Z_sel=0; 
  int gamma_sel=0; 
  int H_sel=0; 
  
  
  
  //--------------------
  //--- loop over events
  std::vector<int> nEvents(20,0);
  std::vector<std::string> nEvents_label(20,"");
  int nEvents_it;
  
  //for(int entry = 0; entry < 1000; ++entry)
  for(int entry = 0; entry < chain_reco->GetEntries() ; ++entry)
  {
    if( !debugMode && entry%100 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_recoTree << "\r" << std::flush;
    if(  debugMode && entry%1 == 0 ) std::cout << ">>> loop 1/1: reading entry " << entry << " / " << nEntries_recoTree << std::endl;
    chain_reco -> GetEntry(entry);
    if( !isData ) chain_gen -> GetEntry(entry);
    
    nEvents_it = 0;
    nEvents_label[nEvents_it] = "tot";
    ++nEvents[nEvents_it];
    
    
    float frac = r.Uniform(0.,1.);
    float weight = 1.;
    float puWeight = 1.;
    if( !isData )
    {
    //   for(unsigned int jj = 0; jj < frac_PU.size()-1; ++jj)
    //     if( frac >= frac_PU.at(jj) && frac < frac_PU.at(jj+1) )
    //       puWeight = puReweighting.at(jj) -> GetPUWeight(treeVars.trueNumInteractions);
      weight = puWeight * mcWeight;
    }
    
    
    ++nTriggerEvents_noCuts;
    for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
    {
      std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
      std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_noCuts.begin(),vec_triggerPass_noCuts.end(),FindPair(treeVars.trgs_name->at(ii)));
      if( it != vec_triggerPass_noCuts.end() ) it->second += treeVars.trgs_pass->at(ii);
      else vec_triggerPass_noCuts.push_back( p );
    }
    
    
    //--------------
    //--- candidates
    particle gen_H;
    particle gen_Z;
    particle gen_gamma;
    std::vector<particle> gen_mu;
    std::vector<particle> gen_ele;
    
    std::vector<particle> recoGenMatch_mu;
    std::vector<particle> recoGenMatch_ele;
    
    particle reco_H;
    particle reco_Z;
    std::vector<particle> reco_gamma;
    std::vector<particle> reco_ele;
    std::vector<particle> reco_mu;
    std::vector<particle> reco_jets;
    
    
    //-------------------------
    //--- fill particle vectors
    if( debugMode ) std::cout << ">>>>>> start filling particle vectors" << std::endl;
    std::vector<particle> temp_mu;
    std::vector<particle> temp_ele;
    std::vector<particle> temp_gamma;
    std::vector<particle> temp_jets;
    
    
    for(unsigned int ii = 0; ii < treeVars.photons_pt->size(); ++ii)
    {
      particle temp;
      temp.v.SetPtEtaPhiE(treeVars.photons_EnergyPostCorr->at(ii)*sin(2.*atan(exp(-1.*treeVars.photons_eta->at(ii)))),treeVars.photons_eta->at(ii),treeVars.photons_phi->at(ii),treeVars.photons_EnergyPostCorr->at(ii));
      temp.it = ii;
      temp_gamma.push_back(temp);
     }
     if( debugMode ) std::cout << ">>>>>> end filling photon vectors" << std::endl;

   
   
    for(unsigned int ii = 0; ii < treeVars.muons_pt->size(); ++ii)
    {
      particle temp;
      
      double SF = 1.;
      // if( isData )
      //   SF = rc.kScaleDT(treeVars.muons_charge->at(ii),treeVars.muons_pt->at(ii),treeVars.muons_eta->at(ii),treeVars.muons_phi->at(ii));
      // if( !isData && treeVars.muons_trackerLayersWithMeasurement->at(ii) >  0 )
      //   SF = rc.kScaleAndSmearMC(treeVars.muons_charge->at(ii),treeVars.muons_pt->at(ii),treeVars.muons_eta->at(ii),treeVars.muons_phi->at(ii),treeVars.muons_trackerLayersWithMeasurement->at(ii),gRandom->Rndm(),gRandom->Rndm());
      
      temp.v.SetPtEtaPhiE(SF*treeVars.muons_pt->at(ii),treeVars.muons_eta->at(ii),treeVars.muons_phi->at(ii),SF*treeVars.muons_energy->at(ii));
      temp.charge = treeVars.muons_charge->at(ii);
      temp.pdgId = -13 * temp.charge;
      temp.it = ii;
      temp_mu.push_back(temp);
    }
    std::sort(temp_mu.begin(),temp_mu.end(),PtSort());
    if( debugMode ) std::cout << ">>>>>> end filling muon vectors" << std::endl;
    
    
    for(unsigned int ii = 0; ii < treeVars.electrons_pt->size(); ++ii)
    {
      particle temp;
      
      temp.v.SetPtEtaPhiE(treeVars.electrons_EnergyPostCorr->at(ii)*sin(2.*atan(exp(-1.*treeVars.electrons_eta->at(ii)))),treeVars.electrons_eta->at(ii),treeVars.electrons_phi->at(ii),treeVars.electrons_EnergyPostCorr->at(ii));
      temp.charge = treeVars.electrons_charge->at(ii);
      temp.pdgId = -11 * temp.charge;
      temp.it = ii;
      temp_ele.push_back(temp);
    }
    if( debugMode ) std::cout << ">>>>>> end filling electron vectors" << std::endl;
    
    
    for(unsigned int ii = 0; ii < treeVars.jets_pt->size(); ++ii)
    {
      particle temp;
           
      temp.v.SetPtEtaPhiE(treeVars.jets_pt->at(ii),treeVars.jets_eta->at(ii),treeVars.jets_phi->at(ii),treeVars.jets_energy->at(ii));
      temp.it = ii;
      temp_jets.push_back(temp);
     }
    if( debugMode ) std::cout << ">>>>>> end filling jet vectors" << std::endl;
    
    
    
    //---------------------
    //--- retrieve gen info
    if( isSignal )
    {
      if( debugMode ) std::cout << ">>>>>> start retrieve gen info" << std::endl;
      
      gen_H.v.SetPtEtaPhiE(treeVars.reso_pt->at(0),treeVars.reso_eta->at(0),treeVars.reso_phi->at(0),treeVars.reso_energy->at(0));
      gen_H.charge = treeVars.reso_charge->at(0);
      gen_H.pdgId = treeVars.reso_pdgId->at(0);
      
      for(int ii = 0; ii < (*treeVars.resoDau1_n)[0]; ++ii)
      {
        if( fabs((*treeVars.resoDau1_pdgId)[0][ii]) == 23 )
        { 
          for(int jj = 0; jj < (*treeVars.resoDau2_n)[0][ii]; ++jj)
          {
            if( fabs((*treeVars.resoDau2_pdgId)[0][ii][jj]) == 13 )
            {
              particle temp;
              temp.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj],(*treeVars.resoDau2_eta)[0][ii][jj],(*treeVars.resoDau2_phi)[0][ii][jj],(*treeVars.resoDau2_energy)[0][ii][jj]);
              temp.charge = (*treeVars.resoDau2_charge)[0][ii][jj];
              temp.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj];
              
              gen_mu.push_back( temp );
            }
            if( fabs((*treeVars.resoDau2_pdgId)[0][ii][jj]) == 11 )
            {
              particle temp;
              temp.v.SetPtEtaPhiE((*treeVars.resoDau2_pt)[0][ii][jj],(*treeVars.resoDau2_eta)[0][ii][jj],(*treeVars.resoDau2_phi)[0][ii][jj],(*treeVars.resoDau2_energy)[0][ii][jj]);
              temp.charge = (*treeVars.resoDau2_charge)[0][ii][jj];
              temp.pdgId = (*treeVars.resoDau2_pdgId)[0][ii][jj];
              
              gen_ele.push_back( temp );
            }
          }
        }
        
        if( fabs((*treeVars.resoDau1_pdgId)[0][ii]) == 22 )
        { 
          gen_gamma.v.SetPtEtaPhiM((*treeVars.resoDau1_pt)[0][ii],(*treeVars.resoDau1_eta)[0][ii],(*treeVars.resoDau1_phi)[0][ii],0.);
          gen_gamma.charge = (*treeVars.resoDau1_charge)[0][ii];
          gen_gamma.pdgId = (*treeVars.resoDau1_pdgId)[0][ii];
        } 
      }
      if( debugMode ) std::cout << ">>>>>> end retrieve gen info" << std::endl;
      
      std::sort(gen_mu.begin(),gen_mu.end(),PtSort());
      std::sort(gen_ele.begin(),gen_ele.end(),PtSort());
      
      
      //match mu
      if( genMatch && gen_mu.size() >= 2 )
      {
        if( debugMode ) std::cout << ">>>>>> start gen match" << std::endl;
        
        bool recoGenMatch_mu_matching = false;
        
        //-----------------------
        //--- match reco with gen      
        std::vector<int> muIt_genMatch;
        muIt_genMatch.push_back( GetBestMatch(gen_mu.at(0),temp_mu,&muIt_genMatch) );
        muIt_genMatch.push_back( GetBestMatch(gen_mu.at(1),temp_mu,&muIt_genMatch) );
        if( muIt_genMatch.at(0) >= 0 ) recoGenMatch_mu.push_back( temp_mu.at(muIt_genMatch.at(0)) );
        if( muIt_genMatch.at(1) >= 0 ) recoGenMatch_mu.push_back( temp_mu.at(muIt_genMatch.at(1)) );
        
        // fill resolution plots
        if( recoGenMatch_mu.size() > 0 )
        {
          if( IsMatching(gen_mu.at(0),recoGenMatch_mu.at(0),0.01,0.5) )
          {
            h_ptReso_mu -> Fill( recoGenMatch_mu.at(0).v.Pt() - gen_mu.at(0).v.Pt() );
            
            int binPt  = h_occupancyPt_mu -> Fill(gen_mu.at(0).v.Pt());
            int binEta = h_occupancyEta_mu -> Fill(gen_mu.at(0).v.Eta());
            if( binPt >= 1 && binPt <= RangePtmu.size() &&
                binEta >= 1 && binEta <= RangeEtamu.size() )

            {
              float etaminmu = RangeEtamu.at(binEta-1);
              float etamaxmu = RangeEtamu.at(binEta);
              float ptminmu = RangePtmu.at(binPt-1);
              float ptmaxmu = RangePtmu.at(binPt);
              std::string labelmu(Form("Res_eta_%.1f-%.1f_pt%.0f-%.0f_mu",etaminmu,etamaxmu,ptminmu,ptmaxmu));
              histosmu[labelmu] -> Fill( recoGenMatch_mu.at(0).v.Pt() - gen_mu.at(0).v.Pt() );
            }
          }
        }
        
        if( recoGenMatch_mu.size() > 1 )
        {
          if( IsMatching(gen_mu.at(1),recoGenMatch_mu.at(1),0.01,0.5) )
          {
            h_ptReso_mu -> Fill( recoGenMatch_mu.at(1).v.Pt() - gen_mu.at(1).v.Pt() );

	    int binPt  = h_occupancyPt_mu -> Fill(gen_mu.at(1).v.Pt());
            int binEta = h_occupancyEta_mu -> Fill(gen_mu.at(1).v.Eta());
            if( binPt >= 1 && binPt <= RangePtmu.size() &&
                binEta >= 1 && binEta <= RangeEtamu.size() )
            {
              float etaminmu = RangeEtamu.at(binEta-1);
              float etamaxmu = RangeEtamu.at(binEta);
              float ptminmu = RangePtmu.at(binPt-1);
              float ptmaxmu = RangePtmu.at(binPt);
              std::string labelmu(Form("Res_eta_%.1f-%.1f_pt%.0f-%.0f_mu",etaminmu,etamaxmu,ptminmu,ptmaxmu));
              histosmu[labelmu] -> Fill( recoGenMatch_mu.at(1).v.Pt() - gen_mu.at(1).v.Pt() );
            }
          }
        }
      }
      
      //match ele
      if( genMatch && gen_ele.size() >= 2 )
      {
        if( debugMode ) std::cout << ">>>>>> start gen match" << std::endl;
        
        bool recoGenMatch_ele_matching = false;
        
        //-----------------------
        //--- match reco with gen      
        std::vector<int> eleIt_genMatch;
        eleIt_genMatch.push_back( GetBestMatch(gen_ele.at(0),temp_ele,&eleIt_genMatch) );
        eleIt_genMatch.push_back( GetBestMatch(gen_ele.at(1),temp_ele,&eleIt_genMatch) );
        if( eleIt_genMatch.at(0) >= 0 ) recoGenMatch_ele.push_back( temp_ele.at(eleIt_genMatch.at(0)) );
        if( eleIt_genMatch.at(1) >= 0 ) recoGenMatch_ele.push_back( temp_ele.at(eleIt_genMatch.at(1)) );
        
       // fill resolution plots
       if( recoGenMatch_ele.size() > 0 )
       {
         if( IsMatching(gen_ele.at(0),recoGenMatch_ele.at(0),0.01,0.5) )
         {
           h_ptReso_ele -> Fill( recoGenMatch_ele.at(0).v.Pt() - gen_ele.at(0).v.Pt() );

	   int binPt  = h_occupancyPt_ele -> Fill(gen_ele.at(0).v.Pt());
           int binEta = h_occupancyEta_ele -> Fill(gen_ele.at(0).v.Eta());
           if( binPt >= 1 && binPt > RangePtele.size() &&
               binEta >= 1 && binEta <= RangeEtaele.size() )
           {
             float etaminele = RangeEtaele.at(binEta-1);
             float etamaxele = RangeEtaele.at(binEta);
             float ptminele = RangePtele.at(binPt-1);
             float ptmaxele = RangePtele.at(binPt);
             std::string labelele(Form("Res_eta_%.1f-%.1f_pt%.0f-%.0f_ele",etaminele,etamaxele,ptminele,ptmaxele));
             histosele[labelele] -> Fill(recoGenMatch_ele.at(0).v.Pt() - gen_ele.at(0).v.Pt());
           }
         }
       }
       
       if( recoGenMatch_ele.size() > 1 )
       {
         if( IsMatching(gen_ele.at(1),recoGenMatch_ele.at(1),0.01,0.5) )
         {
           h_ptReso_ele -> Fill( recoGenMatch_ele.at(1).v.Pt() - gen_ele.at(1).v.Pt() );

	   int binPt  = h_occupancyPt_ele -> Fill(gen_ele.at(1).v.Pt());
           int binEta = h_occupancyEta_ele -> Fill(gen_ele.at(1).v.Eta());
           if( binPt >= 1 && binPt <= RangePtele.size() &&
               binEta >= 1 && binEta <= RangeEtaele.size() )
           {
             float etaminele = RangeEtaele.at(binEta-1);
             float etamaxele = RangeEtaele.at(binEta);
             float ptminele = RangePtele.at(binPt-1);
             float ptmaxele = RangePtele.at(binPt);
             std::string labelele(Form("Res_eta_%.1f-%.1f_pt%.0f-%.0f_ele",etaminele,etamaxele,ptminele,ptmaxele));
             histosele[labelele] -> Fill(recoGenMatch_ele.at(1).v.Pt() - gen_ele.at(1).v.Pt());
           }
         }
       }
      }
    }
    
    
    
    //---------------------
    //--- trigger selection
    if( debugMode ) std::cout << ">>>>>> start trigger selection" << std::endl;
    bool trgPassAll = true;
    bool trgPassOne = false;
    for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
    {
       //std::cout << "HLT name: " << treeVars.trgs_name->at(ii) << "   pass: " << treeVars.trgs_pass->at(ii) << "    prescale: " << treeVars.trgs_prescale->at(ii) << std::endl;
      for(unsigned int jj = 0; jj < reco_HLT_paths.size(); ++jj)
      {
        std::size_t found = treeVars.trgs_name->at(ii).find( reco_HLT_paths.at(jj) );
        if (found != std::string::npos )
        {
         // std::cout << "HLT name: " << treeVars.trgs_name->at(ii) << "   pass: " << treeVars.trgs_pass->at(ii) << "    prescale: " << treeVars.trgs_prescale->at(ii) << std::endl;
          if( treeVars.trgs_pass->at(ii) == 1 )
          {
            trgPassOne = true;
            std::pair<std::string,int> p(treeVars.trgs_name->at(ii),1);
            std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass.begin(),vec_triggerPass.end(),FindPair(treeVars.trgs_name->at(ii)));
            if( it != vec_triggerPass.end() ) it->second += 1;
            else vec_triggerPass.push_back( p );
          }
          else
          {
            trgPassAll = false;
          }
        }
      }
    }
    // if( reco_HLTCut && (reco_HLT_OR  && !trgPassOne) ) continue;
    // if( reco_HLTCut && (!reco_HLT_OR && !trgPassAll) ) continue;
    if( debugMode ) std::cout << ">>>>>> end trigger selection" << std::endl;
    
    nEvents_it += 1;
    nEvents_label[nEvents_it] = "trigger selection";
    ++nEvents[nEvents_it];
    
    
    
    // da qui parte la divisione lettroni e muoni
    
    
    
    //------------
    //--- do muons
    if (doMuon)
    {
      if( debugMode ) std::cout << ">>>>>> start reco mu selection" << std::endl;
      
      if( temp_gamma.size() < 1 ) continue;
      if( temp_mu.size() < 2 ) continue;
      no_sel++; 
      
      
      float mu1_DR;
      float mu1_iso;
      for(unsigned int ii = 0; ii < temp_mu.size(); ++ii)
      {
        if( temp_mu.at(ii).v.Pt() < reco_ptMin1_mu ) continue;
        if( treeVars.muons_isTight->at(temp_mu.at(ii).it) != 1 ) continue;
        if( fabs(temp_mu.at(ii).v.Eta()) > reco_eta_mu ) continue;
        if( fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)) > 0.05 ) continue;
        if( fabs(treeVars.muons_dz->at(temp_mu.at(ii).it))  > 0.10 ) continue;
      
        mu1_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
        std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
        
        if( mu1_iso/temp_mu.at(ii).v.Pt() > 0.35 ) continue;
        reco_mu.push_back( temp_mu.at(ii) );
        break;
      }
      if( reco_mu.size() < 1 ) continue;
      
      float mu2_DR;
      float mu2_iso;
      for(unsigned int ii = 0; ii < temp_mu.size(); ++ii)
      {
        if( reco_mu.at(0).charge*temp_mu.at(ii).charge == 1 ) continue;
        if( temp_mu.at(ii).v.Pt() < reco_ptMin2_mu ) continue;
        if( fabs(temp_mu.at(ii).v.Eta()) > reco_eta_mu ) continue;
        if( treeVars.muons_isTight->at(temp_mu.at(ii).it) != 1 ) continue;
        if( fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)) > 0.05 ) continue;
        if( fabs(treeVars.muons_dz->at(temp_mu.at(ii).it))  > 0.10 ) continue;
      
        mu2_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
        std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
        if( mu2_iso/temp_mu.at(ii).v.Pt() > 0.35 ) continue;

        reco_mu.push_back( temp_mu.at(ii) );
        break;
      }
      if( reco_mu.size() < 2 ) continue;
      lep_sel++; 
      
      
      reco_Z.v = reco_mu.at(0).v + reco_mu.at(1).v;
      reco_Z.charge = reco_mu.at(0).charge + reco_mu.at(1).charge; 
      if ( reco_Z.v.M() < 50 ) continue;
      Z_sel++;  
      
      
      ++nTriggerEvents_recoLepSelection;
      for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
      {
        std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
        std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_recoLepSelection.begin(),vec_triggerPass_recoLepSelection.end(),FindPair(treeVars.trgs_name->at(ii)));
        if( it != vec_triggerPass_recoLepSelection.end() ) it->second += treeVars.trgs_pass->at(ii);
        else vec_triggerPass_recoLepSelection.push_back( p );
      }
      
      
      for(unsigned int ii = 0; ii < temp_gamma.size(); ++ii)
      {
        if( temp_gamma.at(ii).v.Pt() < reco_ptMin_gamma ) continue;
        if( fabs(temp_gamma.at(ii).v.Eta()) > reco_eta_gamma ) continue;
        if( fabs(temp_gamma.at(ii).v.Eta()) > reco_eta1_gamma && fabs(temp_gamma.at(ii).v.Eta()) < reco_eta2_gamma) continue;
	if (fabs(temp_gamma.at(ii).v.Eta()) < reco_eta1_gamma && treeVars.photons_MVAID->at(temp_gamma.at(ii).it) < reco_ID_gammaB) continue;
	if (fabs(temp_gamma.at(ii).v.Eta()) > reco_eta2_gamma && treeVars.photons_MVAID->at(temp_gamma.at(ii).it) < reco_ID_gammaE) continue;
        
        mu1_DR = DeltaR(temp_gamma.at(ii).v.Eta(),temp_gamma.at(ii).v.Phi(),reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi());
        if (mu1_DR < reco_DR_lep) continue;
        mu2_DR = DeltaR(temp_gamma.at(ii).v.Eta(),temp_gamma.at(ii).v.Phi(),reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi());
        if (mu2_DR < reco_DR_lep) continue;
        
        reco_gamma.push_back( temp_gamma.at(ii) );
        break;
      }
      
      if( reco_gamma.size() < 1 ) continue;
      gamma_sel++; 
      
      
      
      //--------Filter for overlap removal
      bool skipEvent = false;
      if (doFilter)
      {
        float DRMin = 999999;
        int itMin = -1;
        for(unsigned int jj = 0; jj < treeVars.genPho_HardProcFinState->size(); ++jj)
        {
          float DR = DeltaR(treeVars.genPho_eta->at(jj),treeVars.genPho_phi->at(jj),reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi());
          float ptRatio = treeVars.genPho_pt->at(jj) / reco_gamma.at(0).v.Pt();
          if( DR < DRMin && ptRatio > 0.2 && ptRatio < 5 )
          {
            DRMin = DR;
            itMin = jj;
          }
        }
        
        if( itMin >= 0 && DRMin < 0.1)
        {
          if( (treeVars.genPho_HardProcFinState->at(itMin) == 1) || (treeVars.genPho_isPromptFinState->at(itMin) == 1) )
          {
            skipEvent = true;
            
            // std::cout << "gen_pt: " << treeVars.genPho_pt->at(itMin) << "   gen_eta: " << treeVars.genPho_eta->at(itMin) << "   gen_phi: " << treeVars.genPho_phi->at(itMin) << std::endl;
            // std::cout << "reco_pt: " << reco_gamma.at(0).v.Pt() << "   reco_eta: " << reco_gamma.at(0).v.Eta() << "   reco_phi: " << reco_gamma.at(0).v.Phi() << std::endl;
            // std::cout << "DRMin: " << DRMin << std::endl;
            // std::cout << treeVars.genPho_HardProcFinState->at(itMin) << "   " << treeVars.genPho_isPromptFinState->at(itMin) << std::endl;
            // std::cout << std::endl;
          }
        }
      }
      if( doFilter && skipEvent ) continue;
      
      
      
      reco_H.v = reco_mu.at(0).v + reco_mu.at(1).v + reco_gamma.at(0).v;
      reco_H.charge = reco_mu.at(0).charge + reco_mu.at(1).charge; 
      reco_H.pdgId = 25;
      
      if( reco_H.v.M() < 100 || reco_H.v.M() > 180 ) continue;
      if( ((reco_gamma.at(0).v.Pt()) /(reco_H.v.M())) < reco_pToM_gamma ) continue;
      if( (reco_H.v + reco_Z.v).M() < reco_sumM ) continue;
      H_sel++;
      if( debugMode ) std::cout << ">>>>>> end reco mu selection" << std::endl;
      
      
      ++nTriggerEvents_AllSelection;
      for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
      {
        std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
        std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_AllSelection.begin(),vec_triggerPass_AllSelection.end(),FindPair(treeVars.trgs_name->at(ii)));
        if( it != vec_triggerPass_AllSelection.end() ) it->second += treeVars.trgs_pass->at(ii);
        else vec_triggerPass_AllSelection.push_back( p );
      }
      
      
      //-------------
      // third lepton
      if( debugMode ) std::cout << ">>>>>> start third lepton" << std::endl;
      
      float mu3_iso;
      for(unsigned int ii = 0; ii < temp_mu.size(); ++ii) 
      {
        mu3_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
          std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
        
	if( int(ii)==reco_mu.at(0).it || int(ii)==reco_mu.at(1).it ) continue;
        if( temp_mu.at(ii).v.Pt() < 5 ) continue;
        if( treeVars.muons_isTight->at(temp_mu.at(ii).it) != 1 ) continue;
        if( fabs(temp_mu.at(ii).v.Eta()) > reco_eta_mu ) continue;
        if( fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)) > 0.05 ) continue;
        if( fabs(treeVars.muons_dz->at(temp_mu.at(ii).it))  > 0.10 ) continue;
        if( mu3_iso/temp_mu.at(ii).v.Pt() > 0.35 ) continue;
        if( DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),temp_mu.at(ii).v.Eta(),temp_mu.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),temp_mu.at(ii).v.Eta(),temp_mu.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi(),temp_mu.at(ii).v.Eta(),temp_mu.at(ii).v.Phi()) < reco_DR_lep ) continue;
        
        reco_mu.push_back( temp_mu.at(ii) );
        
        break;
      }
      
      for(unsigned int ii = 0; ii < temp_ele.size(); ++ii)
      {
        if( temp_ele.at(ii).v.Pt() < 7) continue;
        
        if( fabs(temp_ele.at(ii).v.Eta()) > reco_eta_ele ) continue;
      	if( treeVars.electrons_MVAID->at(temp_ele.at(ii).it) < reco_ID_ele ) continue;
        if( fabs(treeVars.electrons_dxy->at(temp_ele.at(ii).it)) > 0.05 ) continue;
        if( fabs(treeVars.electrons_dz->at(temp_ele.at(ii).it))  > 0.10 ) continue;
        
        if( DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),temp_ele.at(ii).v.Eta(),temp_ele.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),temp_ele.at(ii).v.Eta(),temp_ele.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi(),temp_ele.at(ii).v.Eta(),temp_ele.at(ii).v.Phi()) < reco_DR_lep ) continue;
        reco_ele.push_back( temp_ele.at(ii) );
        
        break;
      }
      if( debugMode ) std::cout << "MU size "   << reco_mu.size() << std::endl;
      if( debugMode ) std::cout << "ELE size "   << reco_ele.size() << std::endl;
      
      
      //-----
      // Jets
      if( debugMode ) std::cout << ">>>>>> start jets" << std::endl;
      
      for(unsigned int ii = 0; ii < temp_jets.size(); ++ii)
      {
        if (reco_jets.size() == 2) break; 
        
        if( temp_jets.at(ii).v.Et() < 30. ) continue;
        if( fabs(temp_jets.at(ii).v.Eta()) > 4.7 ) continue;
        
        bool skipJet = false;
        
        for(unsigned int jj = 0; jj < reco_mu.size(); ++jj)
          if( DeltaR(temp_jets.at(ii).v.Eta(),temp_jets.at(ii).v.Phi(),reco_mu.at(jj).v.Eta(),reco_mu.at(jj).v.Phi()) < 0.4 ) skipJet = true;
        if( DeltaR(temp_jets.at(ii).v.Eta(),temp_jets.at(ii).v.Phi(),reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi()) < 0.4 ) skipJet = true;
	if( reco_ele.size() > 0 )
          if (DeltaR(temp_jets.at(ii).v.Eta(),temp_jets.at(ii).v.Phi(),reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi())<0.4 ) skipJet = true;
        
        if( skipJet ) continue;
        
        for(unsigned int kk = ii +1 ; kk < temp_jets.size(); ++kk)
        { 
          particle dijet;
          dijet.v = temp_jets.at(ii).v + temp_jets.at(kk).v;
          
          if( temp_jets.at(kk).v.Et() < 30. ) continue;
          if( fabs(treeVars.jets_eta->at(kk)) > 4.7 ) continue;
          
          bool skipJet2 = false;
          
          for(unsigned int jj = 0; jj < reco_mu.size(); ++jj)
            if( DeltaR(temp_jets.at(kk).v.Eta(),temp_jets.at(kk).v.Phi(),reco_mu.at(jj).v.Eta(),reco_mu.at(jj).v.Phi()) < 0.4 ) skipJet2 = true;
          if( DeltaR(temp_jets.at(kk).v.Eta(),temp_jets.at(kk).v.Phi(),reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi()) < 0.4 ) skipJet2 = true;
          if( reco_ele.size() > 0 )
            if (DeltaR(temp_jets.at(kk).v.Eta(),temp_jets.at(kk).v.Phi(),reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi())<0.4 ) skipJet2 = true;
          
          if( skipJet2 ) continue;
          
          if( fabs(DeltaEta(temp_jets.at(ii).v.Eta(),temp_jets.at(kk).v.Eta())) < 3.5 ) continue;
          if( (reco_H.v.Eta()-(temp_jets.at(ii).v.Eta()+temp_jets.at(kk).v.Eta())/2.) > 2.5 ) continue;
          
          if( dijet.v.M() < 500 ) continue;
          if( DeltaPhi(reco_H.v.Phi(),dijet.v.Phi()) < 2.4 ) continue;
          
          reco_jets.push_back(temp_jets.at(ii));
          reco_jets.push_back(temp_jets.at(kk));
          
          break;
        }
      }
      if( debugMode ) std::cout << ">>>>>> end jets" << std::endl;
      
      
      
      //---------------
      // --- categories
      int cat_n; 
      bool isCat = false;
      if( reco_ele.size() > 0 || reco_mu.size() > 2 ) // leptons
      {
        cat_n = 6 ; 
        isCat = true; 
      }
      
      if( !isCat && reco_jets.size() == 2 ) // dijets
      {
        cat_n = 5 ; 
        isCat = true; 
      }
      
      if( !isCat && reco_H.v.Pt() > 60 ) // boosted
      {
        cat_n = 7; 
        isCat = true; 
      }
      
      if( !isCat && (fabs(reco_gamma.at(0).v.Eta()) > 0 && fabs(reco_gamma.at(0).v.Eta()) < 1.4442) ) // untagged 1,2,3
      {
        if( (fabs(reco_mu.at(0).v.Eta()) > 0 && fabs(reco_mu.at(0).v.Eta()) < 2.1) && 
            (fabs(reco_mu.at(1).v.Eta()) > 0 && fabs(reco_mu.at(1).v.Eta()) < 2.1) &&
            ((fabs(reco_mu.at(0).v.Eta()) > 0 && fabs(reco_mu.at(0).v.Eta()) < 0.9)||(fabs(reco_mu.at(1).v.Eta()) > 0 && fabs(reco_mu.at(1).v.Eta()) < 0.9)) )
        {
          if( treeVars.photons_full5x5_R9-> at(reco_gamma.at(0).it) > 0.94 )
          {
            cat_n = 1 ;
            isCat = true; 
          }
          else 
          {
            cat_n = 2;
            isCat = true; 
          }
        }
        
        else if( (fabs(reco_mu.at(0).v.Eta()) > 0.9 && fabs(reco_mu.at(1).v.Eta()) > 0.9) || 
                 ((fabs(reco_mu.at(0).v.Eta()) > 2.1 && fabs(reco_mu.at(0).v.Eta()) < 2.4) || (fabs(reco_mu.at(1).v.Eta())>2.1 && fabs(reco_mu.at(1).v.Eta())< 2.4)) )
        {
          cat_n = 3; 
          isCat = true; 
        }
      }
      
      if( !isCat && (fabs(reco_gamma.at(0).v.Eta()) > 1.566 && fabs(reco_gamma.at(0).v.Eta()) < 2.5) && //untagged 4
          ((fabs(reco_mu.at(0).v.Eta()) > 0 && fabs(reco_mu.at(0).v.Eta()) < 2.4) && 
           (fabs(reco_mu.at(1).v.Eta()) > 0 && fabs(reco_mu.at(1).v.Eta()) < 2.4 )) )
      {
        cat_n = 4 ;
        isCat = true; 
      }
      
      if( !isCat ) cat_n = -100; 
      
      
      //---------------------------
      //--- fill out tree variables
      
      if( debugMode ) std::cout << ">>>>>> start filling out tree" << std::endl;
      outTree.vtxs_n = treeVars.vtxs_n;
      outTree.rho_all = treeVars.rho_all;
      outTree.cat_n = cat_n;
      
      outTree.gamma_pt = reco_gamma.at(0).v.Pt();
      outTree.gamma_eta = reco_gamma.at(0).v.Eta();
      outTree.gamma_phi = reco_gamma.at(0).v.Phi();    
      outTree.gamma_DR1 = DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi());
      outTree.gamma_DR2 = DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi());
      outTree.gamma_IDMVA = treeVars.photons_MVAID->at(reco_gamma.at(0).it);
      outTree.gamma_full5x5_R9 = treeVars.photons_full5x5_R9->at(reco_gamma.at(0).it);
      outTree.gamma_full5x5_sieie = treeVars.photons_full5x5_sieie->at(reco_gamma.at(0).it);
      
      outTree.mu1_pt = reco_mu.at(0).v.Pt();
      outTree.mu2_pt = reco_mu.at(1).v.Pt();
      // outTree.mu1_ptErr = treeVars.muons_ptErr->at(reco_mu.at(0).it);
      // outTree.mu2_ptErr = treeVars.muons_ptErr->at(reco_mu.at(1).it);
      outTree.mu1_eta = reco_mu.at(0).v.Eta();
      outTree.mu2_eta = reco_mu.at(1).v.Eta();
      outTree.mu1_phi = reco_mu.at(0).v.Phi();
      outTree.mu2_phi = reco_mu.at(1).v.Phi();
      outTree.mu1_dxy = treeVars.muons_dxy->at(reco_mu.at(0).it);
      outTree.mu2_dxy = treeVars.muons_dxy->at(reco_mu.at(1).it);
      outTree.mu1_dxyPull = fabs(treeVars.muons_dxy->at(reco_mu.at(0).it))/treeVars.muons_dxyErr->at(reco_mu.at(0).it);
      outTree.mu2_dxyPull = fabs(treeVars.muons_dxy->at(reco_mu.at(1).it))/treeVars.muons_dxyErr->at(reco_mu.at(1).it);
      outTree.mu1_dz = treeVars.muons_dz->at(reco_mu.at(0).it);
      outTree.mu2_dz = treeVars.muons_dz->at(reco_mu.at(1).it);
      outTree.mu1_dzPull = fabs(treeVars.muons_dz->at(reco_mu.at(0).it))/treeVars.muons_dzErr->at(reco_mu.at(0).it);
      outTree.mu2_dzPull = fabs(treeVars.muons_dz->at(reco_mu.at(1).it))/treeVars.muons_dzErr->at(reco_mu.at(1).it);
      outTree.mu1_relIso = mu1_iso/reco_mu.at(0).v.Pt();
      outTree.mu2_relIso = mu2_iso/reco_mu.at(1).v.Pt();
      outTree.mu1_isL = treeVars.muons_isLoose->at(reco_mu.at(0).it);
      outTree.mu2_isL = treeVars.muons_isLoose->at(reco_mu.at(1).it);
      outTree.mu1_isM = treeVars.muons_isMedium->at(reco_mu.at(0).it);
      outTree.mu2_isM = treeVars.muons_isMedium->at(reco_mu.at(1).it);
      outTree.mu1_isT = treeVars.muons_isTight->at(reco_mu.at(0).it);
      outTree.mu2_isT = treeVars.muons_isTight->at(reco_mu.at(1).it);
      outTree.mu_Deta = DeltaEta(reco_mu.at(0).v.Eta(),reco_mu.at(1).v.Eta());
      outTree.mu_Dphi = DeltaPhi(reco_mu.at(0).v.Phi(),reco_mu.at(1).v.Phi());
      outTree.mu_DR = DeltaR(reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi(),reco_mu.at(1).v.Eta(),reco_mu.at(1).v.Phi());
      
      if( reco_mu.size() > 2 )
      {
        outTree.mu3_pt = reco_mu.at(2).v.Pt();
        outTree.mu3_eta = reco_mu.at(2).v.Eta();
        outTree.mu3_phi = reco_mu.at(2).v.Phi();
        outTree.mu3_dxy = treeVars.muons_dxy->at(reco_mu.at(2).it);
        outTree.mu3_dxyPull = fabs(treeVars.muons_dxy->at(reco_mu.at(2).it))/treeVars.muons_dxyErr->at(reco_mu.at(2).it);
        outTree.mu3_dz = treeVars.muons_dz->at(reco_mu.at(2).it);
        outTree.mu3_dzPull = fabs(treeVars.muons_dz->at(reco_mu.at(2).it))/treeVars.muons_dzErr->at(reco_mu.at(2).it);
        outTree.mu3_relIso = mu3_iso/reco_mu.at(2).v.Pt();
        outTree.mu3_isL = treeVars.muons_isLoose->at(reco_mu.at(2).it);
        outTree.mu3_isM = treeVars.muons_isMedium->at(reco_mu.at(2).it);
        outTree.mu3_isT = treeVars.muons_isTight->at(reco_mu.at(2).it);
      }
      else
      {
        outTree.mu3_pt = -100.;
        outTree.mu3_eta = -100.;
        outTree.mu3_phi = -100.;
        outTree.mu3_dxy = -100.;
        outTree.mu3_dxyPull = -100.;
        outTree.mu3_dz = -100.;
        outTree.mu3_dzPull = -100.;
        outTree.mu3_relIso = -100.;
        outTree.mu3_isL = -100.;
        outTree.mu3_isM = -100.;
        outTree.mu3_isT = -100.;
      }
      
      if( reco_ele.size() > 0 )
      {
        outTree.ele1_pt = reco_ele.at(0).v.Pt();
        outTree.ele1_eta = reco_ele.at(0).v.Eta();
        outTree.ele1_phi = reco_ele.at(0).v.Phi();
        //outTree.ele1_dxy = treeVars.electrons_dxy->at(reco_ele.at(0).it);
        //outTree.ele1_dxyPull = fabs(treeVars.electrons_dxy->at(reco_ele.at(0).it))/treeVars.electrons_dxyErr->at(reco_ele.at(0).it);
        //outTree.ele1_dz = treeVars.electrons_dz->at(reco_ele.at(0).it);
        //outTree.ele1_dzPull = fabs(treeVars.electrons_dz->at(reco_ele.at(0).it))/treeVars.electrons_dzErr->at(reco_ele.at(0).it);
        outTree.ele1_IDMVA = treeVars.electrons_MVAID->at(reco_ele.at(0).it);
        outTree.ele1_full5x5_R9 = treeVars.electrons_full5x5_R9->at(reco_ele.at(0).it);
        outTree.ele1_full5x5_sieie = treeVars.electrons_full5x5_sieie->at(reco_ele.at(0).it);
      }
      else
      {
        outTree.ele1_pt = -100;
        outTree.ele1_eta = -100;
        outTree.ele1_phi = -100;
        //outTree.ele1_dxy = treeVars.electrons_dxy->at(reco_ele.at(0).it);
        //outTree.ele1_dxyPull = fabs(treeVars.electrons_dxy->at(reco_ele.at(0).it))/treeVars.electrons_dxyErr->at(reco_ele.at(0).it);
        //outTree.ele1_dz = treeVars.electrons_dz->at(reco_ele.at(0).it);
        //outTree.ele1_dzPull = fabs(treeVars.electrons_dz->at(reco_ele.at(0).it))/treeVars.electrons_dzErr->at(reco_ele.at(0).it);
        outTree.ele1_IDMVA = -100;
        outTree.ele1_full5x5_R9 = -100;
        outTree.ele1_full5x5_sieie = -100;
      }
      
      outTree.jet1_pt = (reco_jets.size()==2 ? reco_jets.at(0).v.Pt() : -100);
      outTree.jet1_eta = (reco_jets.size()==2 ? reco_jets.at(0).v.Eta() : -100);
      outTree.jet1_phi = (reco_jets.size()==2 ? reco_jets.at(0).v.Phi() : -100);
      outTree.jet2_pt = (reco_jets.size()==2 ? reco_jets.at(1).v.Pt() : -100);
      outTree.jet2_eta = (reco_jets.size()==2 ? reco_jets.at(1).v.Eta() : -100);
      outTree.jet2_phi = (reco_jets.size()==2 ? reco_jets.at(1).v.Phi() : -100);
      
      outTree.Z_pt = reco_Z.v.Pt();
      outTree.Z_eta = reco_Z.v.Eta();
      outTree.Z_phi = reco_Z.v.Phi();
      outTree.Z_mass = reco_Z.v.M();

      outTree.H_pt = reco_H.v.Pt();
      outTree.H_eta = reco_H.v.Eta();
      outTree.H_phi = reco_H.v.Phi();
      outTree.H_mass = reco_H.v.M();
      
      outTree.met_pt = treeVars.met_pt;
      outTree.met_phi = treeVars.met_phi;
      outTree.met_sig = treeVars.met_sig;
    
      outTree.weight = weight;
      outTree.weight_MC = mcWeight;
      outTree.weight_PU = puWeight;
      
      outTree.GetTTreePtr()->Fill();
      if( debugMode )  std::cout << ">>>>>> end filling out tree" << std::endl;
      
      //---------------
      //--- print event
      if( debugMode ) std::cout << ">>>>>> start print event" << std::endl;
      if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
        std::cout << "\n\n Event: " << entry << std::endl;
      if( printGenEvent && isSignal )
      {
        std::cout << "--------------------- GEN ---------------------" << std::endl;
        PrintEvent(gen_mu,gen_gamma);
      }
      if( printRecoGenMatchEvent && isSignal )
      {
        std::cout << "--------------------- GEN MATCH ---------------------" << std::endl;
        PrintEvent(recoGenMatch_mu);
      }
      if( printRecoEvent )
      {
        std::cout << "--------------------- RECO ---------------------" << std::endl;
        PrintEvent(reco_mu,reco_gamma.at(0));
      }
      if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
        std::cout << std::endl;
      if( debugMode ) std::cout << ">>>>>> end print event" << std::endl;
      
    } //fine if
    
    
    
    
    
    
    //----------------
    //--- do electrons
    if( doEle )
    {
      if( debugMode ) std::cout << ">>>>>> start reco ele selection" << std::endl;
      
      nEvents_it += 1;
      nEvents_label[nEvents_it] = "start ele selection";
      ++nEvents[nEvents_it];
    
      if( temp_gamma.size() < 1 ) continue;
      
      nEvents_it += 1;
      nEvents_label[nEvents_it] = ">= 1 photon";
      ++nEvents[nEvents_it];
      
      if( temp_ele.size() < 2 ) continue;
      no_sel++;
      
      nEvents_it += 1;
      nEvents_label[nEvents_it] = ">= 2 electrons";
      ++nEvents[nEvents_it];
      
      
      float ele1_DR;
      
      for(unsigned int ii = 0; ii < temp_ele.size(); ++ii)
      {
        
        if( temp_ele.at(ii).v.Pt() < reco_ptMin1_ele ) continue;
        if( fabs(temp_ele.at(ii).v.Eta()) > reco_eta_ele ) continue;
        
        if( fabs(treeVars.electrons_dxy->at(temp_ele.at(ii).it)) > 0.05 ) continue; //
        if( fabs(treeVars.electrons_dz->at(temp_ele.at(ii).it))  > 0.10 ) continue;
        
      	if (treeVars.electrons_MVAID->at(temp_ele.at(ii).it) < reco_ID_ele) continue;
        
        reco_ele.push_back( temp_ele.at(ii) );
        
        break;
      }
      if( reco_ele.size() < 1 ) continue;
      
      float ele2_DR;
      
      for(unsigned int ii = 0; ii < temp_ele.size(); ++ii)
      {
        if( reco_ele.at(0).charge*temp_ele.at(ii).charge == 1 ) continue;
        if( temp_ele.at(ii).v.Pt() < reco_ptMin2_ele ) continue;
        if( fabs(temp_ele.at(ii).v.Eta()) > reco_eta_ele ) continue;
      	if (treeVars.electrons_MVAID->at(temp_ele.at(ii).it) < reco_ID_ele) continue;
        //if( fabs(treeVars.electrons_dxy->at(temp_ele.at(ii).it)) > 0.05 ) continue;
        //if( fabs(treeVars.electrons_dz->at(temp_ele.at(ii).it))  > 0.10 ) continue;
        
        reco_ele.push_back( temp_ele.at(ii) );
        break;
      }
      if( reco_ele.size() < 2 ) continue;
      
      nEvents_it += 1;
      nEvents_label[nEvents_it] = ">= 2 sel. electrons";
      ++nEvents[nEvents_it];
      lep_sel++;
      
      
      reco_Z.v = reco_ele.at(0).v + reco_ele.at(1).v;
      reco_Z.charge = reco_ele.at(0).charge + reco_ele.at(1).charge; 
      if ( reco_Z.v.M() < 50 ) continue;

      nEvents_it += 1;
      nEvents_label[nEvents_it] = "m_{ee} > 50 GeV";
      ++nEvents[nEvents_it];
      Z_sel++;
      
      
      //Trigger eff after leptons selections
      ++nTriggerEvents_recoLepSelection;
      for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
      {
        std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
        std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_recoLepSelection.begin(),vec_triggerPass_recoLepSelection.end(),FindPair(treeVars.trgs_name->at(ii)));
        if( it != vec_triggerPass_recoLepSelection.end() ) it->second += treeVars.trgs_pass->at(ii);
        else vec_triggerPass_recoLepSelection.push_back( p );
      }
      
      
      for(unsigned int ii = 0; ii < temp_gamma.size(); ++ii)
      {
        if(temp_gamma.at(ii).v.Pt() < reco_ptMin_gamma) continue; 
        if( fabs(temp_gamma.at(ii).v.Eta()) > reco_eta_gamma ) continue;
        if( fabs(temp_gamma.at(ii).v.Eta()) > reco_eta1_gamma && fabs(temp_gamma.at(ii).v.Eta()) < reco_eta2_gamma) continue;
	if (fabs(temp_gamma.at(ii).v.Eta()) < reco_eta1_gamma && treeVars.photons_MVAID->at(temp_gamma.at(ii).it) < reco_ID_gammaB) continue;
	if (fabs(temp_gamma.at(ii).v.Eta()) > reco_eta2_gamma && treeVars.photons_MVAID->at(temp_gamma.at(ii).it) < reco_ID_gammaE) continue;
        
        ele1_DR = DeltaR(temp_gamma.at(ii).v.Eta(),temp_gamma.at(ii).v.Phi(),reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi());
        if (ele1_DR < reco_DR_lep) continue;
        ele2_DR = DeltaR(temp_gamma.at(ii).v.Eta(),temp_gamma.at(ii).v.Phi(),reco_ele.at(1).v.Eta(),reco_ele.at(1).v.Phi());
        if (ele2_DR < reco_DR_lep) continue;
        reco_gamma.push_back( temp_gamma.at(ii) );
        break;
      }
      
      if( reco_gamma.size() < 1 ) continue;
      
      nEvents_it += 1;
      nEvents_label[nEvents_it] = ">= 1 sel. photon";
      ++nEvents[nEvents_it];
      gamma_sel++; 
      
      
      
      //--------Filter for overlap removal
      bool skipEvent = false;
      if (doFilter)
      {
        float DRMin = 999999;
        int itMin = -1;
        for(unsigned int jj = 0; jj < treeVars.genPho_HardProcFinState->size(); ++jj)
        {
          float DR = DeltaR(treeVars.genPho_eta->at(jj),treeVars.genPho_phi->at(jj),reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi());
          float ptRatio = treeVars.genPho_pt->at(jj) / reco_gamma.at(0).v.Pt();
          if( DR < DRMin && ptRatio > 0.2 && ptRatio < 5 )
          {
            DRMin = DR;
            itMin = jj;
          }
        }
        
        if( itMin >= 0 && DRMin < 0.1)
        {
          if( (treeVars.genPho_HardProcFinState->at(itMin) == 1) || (treeVars.genPho_isPromptFinState->at(itMin) == 1) )
          {
            skipEvent = true;
            
            // std::cout << "gen_pt: " << treeVars.genPho_pt->at(itMin) << "   gen_eta: " << treeVars.genPho_eta->at(itMin) << "   gen_phi: " << treeVars.genPho_phi->at(itMin) << std::endl;
            // std::cout << "reco_pt: " << reco_gamma.at(0).v.Pt() << "   reco_eta: " << reco_gamma.at(0).v.Eta() << "   reco_phi: " << reco_gamma.at(0).v.Phi() << std::endl;
            // std::cout << "DRMin: " << DRMin << std::endl;
            // std::cout << treeVars.genPho_HardProcFinState->at(itMin) << "   " << treeVars.genPho_isPromptFinState->at(itMin) << std::endl;
            // std::cout << std::endl;
          }
        }
      }
      if( doFilter && skipEvent ) continue;
      
      
      
      reco_H.v = reco_ele.at(0).v + reco_ele.at(1).v + reco_gamma.at(0).v;
      reco_H.charge = reco_ele.at(0).charge + reco_ele.at(1).charge; 
      reco_H.pdgId = 25;
      if ( ((reco_gamma.at(0).v.Pt()) /(reco_H.v.M())) < reco_pToM_gamma ) continue;
      if (reco_H.v.M() < 100 || reco_H.v.M() > 180 ) continue;
      if ((reco_H.v + reco_Z.v).M() < reco_sumM ) continue;

      nEvents_it += 1;
      nEvents_label[nEvents_it] = "H kin. selection";
      ++nEvents[nEvents_it];
      H_sel++;
      
      if( debugMode ) std::cout << ">>>>>> end reco ele selection" << std::endl;
      
      
      //Trigger eff after All selections
      ++nTriggerEvents_AllSelection;
      for(unsigned int ii = 0; ii < treeVars.trgs_name->size(); ++ii)
      {
        std::pair<std::string,int> p(treeVars.trgs_name->at(ii),treeVars.trgs_pass->at(ii));
        std::vector<std::pair<std::string,int> >::iterator it = std::find_if(vec_triggerPass_AllSelection.begin(),vec_triggerPass_AllSelection.end(),FindPair(treeVars.trgs_name->at(ii)));
        if( it != vec_triggerPass_AllSelection.end() ) it->second += treeVars.trgs_pass->at(ii);
        else vec_triggerPass_AllSelection.push_back( p );
      }
      
      //-------------
      // third lepton
      if( debugMode ) std::cout << ">>>>>> start third lepton" << std::endl;
      
      float mu1_iso;
      
      for(unsigned int ii = 0; ii < temp_mu.size(); ++ii) 
      {
        mu1_iso = treeVars.muons_pfIsoChargedHadron->at(temp_mu.at(ii).it) +
          std::max(0.,treeVars.muons_pfIsoNeutralHadron->at(temp_mu.at(ii).it)+treeVars.muons_pfIsoPhoton->at(temp_mu.at(ii).it)-0.5*treeVars.muons_pfIsoPU->at(temp_mu.at(ii).it));
        
        if( temp_mu.at(ii).v.Pt() < 5 ) continue;
        if( treeVars.muons_isTight->at(temp_mu.at(ii).it) != 1 ) continue;
        if( fabs(temp_mu.at(ii).v.Eta()) > reco_eta_mu ) continue;
        if( fabs(treeVars.muons_dxy->at(temp_mu.at(ii).it)) > 0.05 ) continue;
        if( fabs(treeVars.muons_dz->at(temp_mu.at(ii).it))  > 0.10 ) continue;
        if( mu1_iso/temp_mu.at(ii).v.Pt() > 0.35 ) continue;
        if( DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),temp_mu.at(ii).v.Eta(),temp_mu.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi(),temp_mu.at(ii).v.Eta(),temp_mu.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_ele.at(1).v.Eta(),reco_ele.at(1).v.Phi(),temp_mu.at(ii).v.Eta(),temp_mu.at(ii).v.Phi()) < reco_DR_lep ) continue;
        
        reco_mu.push_back( temp_mu.at(ii) );
        
        break;
      }
      
      for(unsigned int ii = 0; ii < temp_ele.size(); ++ii)
      {
        if( int(ii)==reco_ele.at(0).it || int(ii)==reco_ele.at(1).it ) continue;
        if( temp_ele.at(ii).v.Pt() < 7 ) continue;
        
        if( fabs(temp_ele.at(ii).v.Eta()) > reco_eta_ele ) continue;
        if( treeVars.electrons_MVAID->at(temp_ele.at(ii).it) < reco_ID_ele ) continue;
        if( fabs(treeVars.electrons_dxy->at(temp_ele.at(ii).it)) > 0.05 ) continue;
        if( fabs(treeVars.electrons_dz->at(temp_ele.at(ii).it))  > 0.10 ) continue;
        
        if( DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),temp_ele.at(ii).v.Eta(),temp_ele.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi(),temp_ele.at(ii).v.Eta(),temp_ele.at(ii).v.Phi()) < reco_DR_lep ) continue;
        if( DeltaR(reco_ele.at(1).v.Eta(),reco_ele.at(1).v.Phi(),temp_ele.at(ii).v.Eta(),temp_ele.at(ii).v.Phi()) < reco_DR_lep ) continue;
        reco_ele.push_back( temp_ele.at(ii) );
        
        break;
      }
      if( debugMode ) std::cout << "MU size "   << reco_mu.size() << std::endl;
      if( debugMode ) std::cout << "ELE size "   << reco_ele.size() << std::endl;
      
      
      //-----
      // Jets
      if( debugMode ) std::cout << ">>>>>> start jets" << std::endl;
      
      for(unsigned int ii = 0; ii < temp_jets.size(); ++ii)
      {
        if (reco_jets.size() == 2) break; 
        
        if( temp_jets.at(ii).v.Et() < 30. ) continue;
        if( fabs(temp_jets.at(ii).v.Eta()) > 4.7 ) continue;
        
        bool skipJet = false;
        
        for(unsigned int jj = 0; jj < reco_ele.size(); ++jj)
          if( DeltaR(temp_jets.at(ii).v.Eta(),temp_jets.at(ii).v.Phi(),reco_ele.at(jj).v.Eta(),reco_ele.at(jj).v.Phi()) < 0.4 ) skipJet = true;
        if( DeltaR(temp_jets.at(ii).v.Eta(),temp_jets.at(ii).v.Phi(),reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi()) < 0.4 ) skipJet = true;
        if( reco_mu.size() > 0 )
          if (DeltaR(temp_jets.at(ii).v.Eta(),temp_jets.at(ii).v.Phi(),reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi()) < 0.4 ) skipJet = true;
        
        if( skipJet ) continue;
        
        for(unsigned int kk = ii +1 ; kk < temp_jets.size(); ++kk)
        { 
          particle dijet;
          dijet.v = temp_jets.at(ii).v + temp_jets.at(kk).v;
          
          if( temp_jets.at(kk).v.Et() < 30. ) continue;
          if( fabs(treeVars.jets_eta->at(kk)) > 4.7 ) continue;
          
          bool skipJet2 = false;
          
          for(unsigned int jj = 0; jj < reco_ele.size(); ++jj)
            if( DeltaR(temp_jets.at(kk).v.Eta(),temp_jets.at(kk).v.Phi(),reco_ele.at(jj).v.Eta(),reco_ele.at(jj).v.Phi()) < 0.4 ) skipJet2 = true;
          if( DeltaR(temp_jets.at(kk).v.Eta(),temp_jets.at(kk).v.Phi(),reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi()) < 0.4 ) skipJet2 = true;
          if( reco_mu.size() > 0 )
            if (DeltaR(temp_jets.at(kk).v.Eta(),temp_jets.at(kk).v.Phi(),reco_mu.at(0).v.Eta(),reco_mu.at(0).v.Phi()) < 0.4 ) skipJet = true;
          
          if( skipJet2 ) continue;
        
          if( fabs(DeltaEta(temp_jets.at(ii).v.Eta(),temp_jets.at(kk).v.Eta())) < 3.5 ) continue;
          if( (reco_H.v.Eta()-(temp_jets.at(ii).v.Eta()+temp_jets.at(kk).v.Eta())/2.) > 2.5 ) continue;
          
          if( dijet.v.M() < 500 ) continue;
          if( DeltaPhi(reco_H.v.Phi(),dijet.v.Phi()) < 2.4 ) continue;
          
          reco_jets.push_back(temp_jets.at(ii));
          reco_jets.push_back(temp_jets.at(kk));
          
          break;
        }
      }
      if( debugMode ) std::cout << ">>>>>> end jets" << std::endl;
      
      
      nEvents_it += 1;
      nEvents_label[nEvents_it] = "all categories";
      ++nEvents[nEvents_it];
      
      //---------------
      // --- categories
      int cat_n; 
      bool isCat = false;
      if( reco_ele.size() > 2 || reco_mu.size() > 0 ) // leptons
      {
        cat_n = 6 ; 
        isCat = true; 
        nEvents_label[nEvents_it+6] = "cat. 6";
        ++nEvents[nEvents_it+6];
      }
      
      if( !isCat && reco_jets.size() == 2 ) // dijets
      {
        cat_n = 5;
        isCat = true;
        nEvents_label[nEvents_it+5] = "cat. 5";
        ++nEvents[nEvents_it+5];
      }
      
      if( !isCat && reco_H.v.Pt() > 60 ) // boosted
      {
        cat_n = 7; 
        isCat = true; 
        nEvents_label[nEvents_it+7] = "cat. 7";
        ++nEvents[nEvents_it+7];
      }
      
      if( !isCat && (fabs(reco_gamma.at(0).v.Eta()) > 0 && fabs(reco_gamma.at(0).v.Eta()) < 1.4442) ) // untagged 1,2,3
      {
        if( (fabs(reco_ele.at(0).v.Eta()) > 0 && fabs(reco_ele.at(0).v.Eta()) < 1.4442) && 
            (fabs(reco_ele.at(1).v.Eta()) > 0 && fabs(reco_ele.at(1).v.Eta()) < 1.4442))
        {
          if( treeVars.photons_full5x5_R9->at(reco_gamma.at(0).it) > 0.94 )
          {
            cat_n = 1 ;
            isCat = true; 
            nEvents_label[nEvents_it+1] = "cat. 1";
            ++nEvents[nEvents_it+1];
          }
          else 
          {
            cat_n = 2;
            isCat = true; 
            nEvents_label[nEvents_it+2] = "cat. 2";
            ++nEvents[nEvents_it+2];
          }
        }
        
        else if( (fabs(reco_ele.at(0).v.Eta()) > 1.4442 && fabs(reco_ele.at(0).v.Eta()) < 2.5) || 
                 (fabs(reco_ele.at(1).v.Eta()) > 1.4442 && fabs(reco_ele.at(1).v.Eta())< 2.5) )
        {
          cat_n = 3; 
          isCat = true; 
          nEvents_label[nEvents_it+3] = "cat. 3";
          ++nEvents[nEvents_it+3];
        }
      }
      
      if( !isCat && (fabs(reco_gamma.at(0).v.Eta()) > 1.566 && fabs(reco_gamma.at(0).v.Eta()) < 2.5) && // untagged 4
          ((fabs(reco_ele.at(0).v.Eta()) > 0 && fabs(reco_ele.at(0).v.Eta()) < 2.5) && 
           (fabs(reco_ele.at(1).v.Eta()) > 0 && fabs(reco_ele.at(1).v.Eta()) < 2.5 )) )
      {
        cat_n = 4 ;
        isCat = true; 
        nEvents_label[nEvents_it+4] = "cat. 4";
        ++nEvents[nEvents_it+4];
      }
      
      if (!isCat) cat_n=-100; 
      
      
      //---------------------------
      //--- fill out tree variables
  
      if( debugMode ) std::cout << ">>>>>> start filling out tree" << std::endl;
      outTree.vtxs_n = treeVars.vtxs_n;
      outTree.cat_n = cat_n;
      outTree.rho_all = treeVars.rho_all;
      
      outTree.gamma_pt = reco_gamma.at(0).v.Pt();
      outTree.gamma_eta = reco_gamma.at(0).v.Eta();
      outTree.gamma_phi = reco_gamma.at(0).v.Phi();    
      outTree.gamma_DR1 = DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi());
      outTree.gamma_DR2 = DeltaR(reco_gamma.at(0).v.Eta(),reco_gamma.at(0).v.Phi(),reco_ele.at(1).v.Eta(),reco_ele.at(1).v.Phi());
      outTree.gamma_IDMVA = treeVars.photons_MVAID->at(reco_gamma.at(0).it);
      outTree.gamma_full5x5_R9 = treeVars.photons_full5x5_R9->at(reco_gamma.at(0).it);
      outTree.gamma_full5x5_sieie = treeVars.photons_full5x5_sieie->at(reco_gamma.at(0).it);
      
      outTree.ele1_pt = reco_ele.at(0).v.Pt();
      outTree.ele2_pt = reco_ele.at(1).v.Pt();
      outTree.ele1_eta = reco_ele.at(0).v.Eta();
      outTree.ele2_eta = reco_ele.at(1).v.Eta();
      outTree.ele1_phi = reco_ele.at(0).v.Phi();
      outTree.ele2_phi = reco_ele.at(1).v.Phi();
      outTree.ele1_energy = reco_ele.at(0).v.E();
      outTree.ele2_energy = reco_ele.at(1).v.E();
      outTree.ele1_energyErr = treeVars.electrons_EnergyErrPostCorr->at(reco_ele.at(0).it);
      outTree.ele2_energyErr = treeVars.electrons_EnergyErrPostCorr->at(reco_ele.at(1).it);
      outTree.ele1_IDMVA = treeVars.electrons_MVAID->at(reco_ele.at(0).it);
      outTree.ele2_IDMVA = treeVars.electrons_MVAID->at(reco_ele.at(1).it);
      outTree.ele1_full5x5_R9 = treeVars.electrons_full5x5_R9->at(reco_ele.at(0).it);
      outTree.ele2_full5x5_R9 = treeVars.electrons_full5x5_R9->at(reco_ele.at(1).it);
      outTree.ele1_full5x5_sieie = treeVars.electrons_full5x5_sieie->at(reco_ele.at(0).it);
      outTree.ele2_full5x5_sieie = treeVars.electrons_full5x5_sieie->at(reco_ele.at(1).it);
      outTree.ele_Deta = DeltaEta(reco_ele.at(0).v.Eta(),reco_ele.at(1).v.Eta());
      outTree.ele_Dphi = DeltaPhi(reco_ele.at(0).v.Phi(),reco_ele.at(1).v.Phi());
      outTree.ele_DR = DeltaR(reco_ele.at(0).v.Eta(),reco_ele.at(0).v.Phi(),reco_ele.at(1).v.Eta(),reco_ele.at(1).v.Phi());
      
      if (reco_mu.size() >0 )
      {
        outTree.mu1_pt = reco_mu.at(0).v.Pt();
        outTree.mu1_eta = reco_mu.at(0).v.Eta();
        outTree.mu1_phi = reco_mu.at(0).v.Phi();
        outTree.mu1_dxy = treeVars.muons_dxy->at(reco_mu.at(0).it);
        outTree.mu1_dxyPull = fabs(treeVars.muons_dxy->at(reco_mu.at(0).it))/treeVars.muons_dxyErr->at(reco_mu.at(0).it);
        outTree.mu1_dz = treeVars.muons_dz->at(reco_mu.at(0).it);
        outTree.mu1_dzPull = fabs(treeVars.muons_dz->at(reco_mu.at(0).it))/treeVars.muons_dzErr->at(reco_mu.at(0).it);
        outTree.mu1_relIso = mu1_iso/reco_mu.at(0).v.Pt();
        outTree.mu1_isL = treeVars.muons_isLoose->at(reco_mu.at(0).it);
        outTree.mu1_isM = treeVars.muons_isMedium->at(reco_mu.at(0).it);
        outTree.mu1_isT = treeVars.muons_isTight->at(reco_mu.at(0).it);
      }
      else
      {
        outTree.mu1_pt = -100.;
        outTree.mu1_eta = -100.;
        outTree.mu1_phi = -100.;
        outTree.mu1_dxy = -100.;
        outTree.mu1_dxyPull = -100.;
        outTree.mu1_dz = -100.;
        outTree.mu1_dzPull = -100.;
        outTree.mu1_relIso = -100.;
        outTree.mu1_isL = -100.;
        outTree.mu1_isM = -100.;
        outTree.mu1_isT = -100.;
      }
      
      if (reco_ele.size() >2 )
      {
        outTree.ele3_pt = reco_ele.at(2).v.Pt();
        outTree.ele3_eta = reco_ele.at(2).v.Eta();
        outTree.ele3_phi = reco_ele.at(2).v.Phi();
        //outTree.ele1_dxy = treeVars.electrons_dxy->at(reco_ele.at(0).it);
        //outTree.ele1_dxyPull = fabs(treeVars.electrons_dxy->at(reco_ele.at(0).it))/treeVars.electrons_dxyErr->at(reco_ele.at(0).it);
        //outTree.ele1_dz = treeVars.electrons_dz->at(reco_ele.at(0).it);
        //outTree.ele1_dzPull = fabs(treeVars.electrons_dz->at(reco_ele.at(0).it))/treeVars.electrons_dzErr->at(reco_ele.at(0).it);
        outTree.ele3_IDMVA = treeVars.electrons_MVAID->at(reco_ele.at(2).it);
        outTree.ele3_full5x5_R9 = treeVars.electrons_full5x5_R9->at(reco_ele.at(2).it);
        outTree.ele3_full5x5_sieie = treeVars.electrons_full5x5_sieie->at(reco_ele.at(2).it);
      }
      else
      {
        outTree.ele3_pt = -100;
        outTree.ele3_eta = -100;
        outTree.ele3_phi = -100;
        //outTree.ele1_dxy = treeVars.electrons_dxy->at(reco_ele.at(0).it);
        //outTree.ele1_dxyPull = fabs(treeVars.electrons_dxy->at(reco_ele.at(0).it))/treeVars.electrons_dxyErr->at(reco_ele.at(0).it);
        //outTree.ele1_dz = treeVars.electrons_dz->at(reco_ele.at(0).it);
        //outTree.ele1_dzPull = fabs(treeVars.electrons_dz->at(reco_ele.at(0).it))/treeVars.electrons_dzErr->at(reco_ele.at(0).it);
        outTree.ele3_IDMVA = -100;
        outTree.ele3_full5x5_R9 = -100;
        outTree.ele3_full5x5_sieie = -100;
      }
      
      outTree.jet1_pt = (reco_jets.size()==2 ? reco_jets.at(0).v.Pt() : -100);
      outTree.jet1_eta = (reco_jets.size()==2 ? reco_jets.at(0).v.Eta() : -100);
      outTree.jet1_phi = (reco_jets.size()==2 ? reco_jets.at(0).v.Phi() : -100);
      outTree.jet2_pt = (reco_jets.size()==2 ? reco_jets.at(1).v.Pt() : -100);
      outTree.jet2_eta = (reco_jets.size()==2 ? reco_jets.at(1).v.Eta() : -100);
      outTree.jet2_phi = (reco_jets.size()==2 ? reco_jets.at(1).v.Phi() : -100);
      
      outTree.Z_pt = reco_Z.v.Pt();
      outTree.Z_eta = reco_Z.v.Eta();
      outTree.Z_phi = reco_Z.v.Phi();
      outTree.Z_mass = reco_Z.v.M();
    
      outTree.H_pt = reco_H.v.Pt();
      outTree.H_eta = reco_H.v.Eta();
      outTree.H_phi = reco_H.v.Phi();
      outTree.H_mass = reco_H.v.M();

      outTree.met_pt = treeVars.met_pt;
      outTree.met_phi = treeVars.met_phi;
      outTree.met_sig = treeVars.met_sig;
    
      outTree.weight = weight;
      outTree.weight_MC = mcWeight;
      outTree.weight_PU = puWeight;
      
      outTree.GetTTreePtr()->Fill();
      if( debugMode )  std::cout << ">>>>>> end filling out tree" << std::endl;
      
      
      //---------------
      //--- print event
      if( debugMode ) std::cout << ">>>>>> start print event" << std::endl;
      if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
        std::cout << "\n\n Event: " << entry << std::endl;
      if( printGenEvent && isSignal )
      {
        std::cout << "--------------------- GEN ---------------------" << std::endl;
        PrintEvent(gen_ele,gen_gamma);
      }
      if( printRecoGenMatchEvent && isSignal )
      {
        std::cout << "--------------------- GEN MATCH ---------------------" << std::endl;
        PrintEvent(recoGenMatch_mu);
      }
      if( printRecoEvent )
      {
        std::cout << "--------------------- RECO ---------------------" << std::endl;
        PrintEvent(reco_ele,reco_gamma.at(0));
      }
      if( printGenEvent || printRecoGenMatchEvent || printRecoEvent )
        std::cout << std::endl;
      if( debugMode ) std::cout << ">>>>>> end print event" << std::endl;
    }//fine if ELE
    
    
    
  } // loop over events
  std::cout << std::endl;
  std::cout << "NO selection events " << float(no_sel)*mcWeight*35.9 << std::endl;
  std::cout << "lep selection events " << float(lep_sel)*mcWeight*35.9 << std::endl;
  std::cout << "Z selection events " << float(Z_sel)*mcWeight*35.9 << std::endl;
  std::cout << "gamma selection events " << float(gamma_sel)*mcWeight*35.9 << std::endl;
  std::cout << "H selection events " << float(H_sel)*mcWeight*35.9 << std::endl;

  
  if(doTrigEff)
  {
  // sort triggers
  std::cout<< std::setw(140) <<"------------  Trigger Eff with only preselections ------------- " << std::endl;
  std::sort(vec_triggerPass_noCuts.begin(),vec_triggerPass_noCuts.end(),PairSort());
  for(unsigned int ii = 0; ii < vec_triggerPass_noCuts.size(); ++ii)
  {
    std::cout << std::fixed;
    std::cout << std::setw(100) << vec_triggerPass_noCuts.at(ii).first << "   ";
    std::cout << std::setw(5)  << vec_triggerPass_noCuts.at(ii).second << "   ";
    std::cout << std::setw(5)  << 1.*vec_triggerPass_noCuts.at(ii).second/nTriggerEvents_noCuts << "   ";
    std::cout << std::endl;
  }

  // sort triggers
  std::cout<< std::setw(140) <<"------------  Trigger Eff after leptons selections ------------- " << std::endl;
  std::sort(vec_triggerPass_recoLepSelection.begin(),vec_triggerPass_recoLepSelection.end(),PairSort());
  for(unsigned int ii = 0; ii < vec_triggerPass_recoLepSelection.size(); ++ii)
  {
    std::cout << std::fixed;
    std::cout << std::setw(100) << vec_triggerPass_recoLepSelection.at(ii).first << "   ";
    std::cout << std::setw(5)  << vec_triggerPass_recoLepSelection.at(ii).second << "   ";
    std::cout << std::setw(5)  << 1.*vec_triggerPass_recoLepSelection.at(ii).second/nTriggerEvents_recoLepSelection<< "   ";
    std::cout << std::endl;
  }

  // sort triggers
  std::cout << std::setw(140)<<"------------  Trigger Eff after all selections ------------- " << std::endl;
  std::sort(vec_triggerPass_AllSelection.begin(),vec_triggerPass_AllSelection.end(),PairSort());
  for(unsigned int ii = 0; ii < vec_triggerPass_AllSelection.size(); ++ii)
  {
    std::cout << std::fixed;
    std::cout << std::setw(100) << vec_triggerPass_AllSelection.at(ii).first << "   ";
    std::cout << std::setw(5)  << vec_triggerPass_AllSelection.at(ii).second << "   ";
    std::cout << std::setw(5)  << 1.*vec_triggerPass_AllSelection.at(ii).second/nTriggerEvents_AllSelection << "   ";
    std::cout << std::endl;
  }
  
  }
  
  for(int ii = 0; ii < 20; ++ii)
  {
    h1_nEvents -> SetBinContent(1+ii,nEvents[ii]);
    h1_nEvents -> GetXaxis() -> SetBinLabel(1+ii,nEvents_label[ii].c_str());
  }
  
  // std::cout << std::fixed << std::endl;
  // std::cout << "===================================================="  << std::endl;
  // std::cout << "                              Tot. events: "    << std::setw(6) << nEvents_tot << std::endl;
  // std::cout << std::endl;
  // std::cout << "             Events with 2 gen mu in acceptance: " << std::setw(6) << nEvents_genCut          << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_genCut/nEvents_tot             << "%)" << std::endl;
  // std::cout << ">>>>>>       Events with 2 reco gen. matched mu: " << std::setw(6) << nEvents_2RecoGenMatchMu << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoGenMatchMu/nEvents_genCutMu << "%)" << std::endl;
  // std::cout << std::endl;
  // std::cout << "           Events with 2 reco. muons: " << std::setw(6) << nEvents_2RecoMu << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoMu/nEvents_tot << "%)"
  //           << "   + matching gen.: " << std::setw(6) << nEvents_2RecoMu_genMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_2RecoMu_genMatch/nEvents_2RecoMu << "%)" << std::endl;
  // std::cout << ">>>          Events with 3 reco. tracks: " << std::setw(6) << nEvents_3RecoKepi << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoKepi/nEvents_tot << "%)"
  //           << "   + matching gen.: " << std::setw(6) << nEvents_3RecoKepi_genMatch << " (" << std::setw(5) << std::setprecision(1) << 100.*nEvents_3RecoKepi_genMatch/nEvents_3RecoKepi << "%)" << std::endl;
  // std::cout << "===================================================="  << std::endl;
  // std::cout << std::endl;
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  

  
  return 0;
}
