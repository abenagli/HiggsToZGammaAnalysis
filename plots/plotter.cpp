/*
 g++ -Wall -o plotter `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C plotter.cpp
*/
#ifndef CMS_LUMI_H
#include "CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "tdrstyle.h"
#endif

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>


using namespace std;
float DeltaEta(const float& eta1, const float& eta2);
void MakePlot(TH1F** histos, TString title);

int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_sqrtS = "35.9 fb^{-1} 13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	TString MuonFold = "/afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new/";

	TString EleFold = "/afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new/";

	vector<TString> names = {"plots_mc_ttH", "plots_mc_ggH", "plots_mc_VBF", "plots_mc_WplusH","plots_mc_WminusH","plots_mc_ZH", "plots_mc_ZG", "plots_mc_DY", "plots_data"};

	TChain* ttH;
	TChain* ggH;
	TChain* vbf;
	TChain* wph;
	TChain* wmh;
	TChain* zh;
	
	
	TChain* data;

	TChain* Zgamma;
	TChain* DY;


	float lumiFactor = 35.9; 
	int loop_njet = 0;
	bool isMuon = (argc!=1 ? 0 : 1);

	if(isMuon)
	{	
		cout << "Processing muon tag" << endl;
		ttH = new TChain("outTree");
		ttH -> Add(MuonFold + names[0] + "*.root");		
		ggH = new TChain("outTree");
		ggH -> Add(MuonFold + names[1] + "*.root");		
		vbf = new TChain("outTree");
		vbf -> Add(MuonFold + names[2] + "*.root");		
		wph = new TChain("outTree");
		wph -> Add(MuonFold + names[3] + "*.root");	
		wmh = new TChain("outTree");
		wmh -> Add(MuonFold + names[4] + "*.root");		
		zh = new TChain("outTree");
		zh -> Add(MuonFold + names[5] + "*.root");


		data = new TChain("outTree");
		data -> Add(MuonFold + "plots_data_DoubleMuon_*.root");
		
	
		
		Zgamma = new TChain("outTree"); 
		Zgamma -> Add(MuonFold + names[6] + "*.root");
				
		DY = new TChain("outTree");
		DY -> Add(MuonFold +  names[7] +"*.root");
	

	
	
	}


	else
	{	cout << "Processing electron tag" << endl;
		ttH = new TChain("outTree");
		ttH -> Add(EleFold + names[0] + "*.root");		
		ggH = new TChain("outTree");
		ggH -> Add(EleFold + names[1] + "*.root");		
		vbf = new TChain("outTree");
		vbf -> Add(EleFold + names[2] + "*.root");		
		wph = new TChain("outTree");
		wph -> Add(EleFold + names[3] + "*.root");	
		wmh = new TChain("outTree");
		wmh -> Add(EleFold + names[4] + "*.root");		
		zh = new TChain("outTree");
		zh -> Add(EleFold + names[5] + "*.root");


		data = new TChain("outTree");
		data -> Add(EleFold + "plots_data_DoubleEG_*.root");
		
	
		
		Zgamma = new TChain("outTree"); 
		Zgamma -> Add(EleFold + names[6] + "*.root");
				
		DY = new TChain("outTree");
		DY -> Add(EleFold + names[7] + "*.root");		

	}

	int nhisto=9; //numero di 
	

	float weight_MC = 0;
	int cat_n = 0;

	

	float  Z_pt = 0;
	float Z_eta = 0;
	float Z_phi = 0;
	float  Z_mass = 0;

	float  H_pt = 0;
	float H_eta = 0;
	float H_phi = 0;
	float  H_mass = 0;


	float gamma_pt = 0;
	float  gamma_eta = 0;
	float  gamma_phi = 0;
	float  gamma_IDMVA = 0;

	
	
	float  ele1_pt = 0;
	float  ele1_eta = 0;
	float  ele1_phi = 0;
	float  ele1_IDMVA = 0;

	
	float  mu3_pt = 0;
	float  ele3_pt = 0;
	
	float  mu1_pt = 0;
	float  mu1_eta = 0;
	float  mu1_phi = 0;


	float  jet1_pt = 0;
	float  jet1_eta = 0;
	float  jet2_eta = 0;
	float  jet1_phi = 0;


	float met_pt = 0.;
	float met_phi = 0.;
	int n_lep_extra = 0;

	
	TH1F* H_mass_histo[nhisto][7];
	TH1F* Z_mass_histo[nhisto][7];
	TH1F* H_mass_all_histo[nhisto];
	TH1F* Z_mass_all_histo[nhisto];
	TH1F* photons_pt_histo[nhisto];
	TH1F* photons_eta_histo[nhisto];
	TH1F* photons_phi_histo[nhisto];
	TH1F* photons_MVAID_histo[nhisto];

	TH1F* extra_leptons_n[nhisto];
	TH1F* extra_leptons_pt_histo[nhisto];
	TH1F* electrons_pt_histo[nhisto];
	TH1F* electrons_eta_histo[nhisto];
	TH1F* electrons_phi_histo[nhisto];
	TH1F* electrons_MVAID_histo[nhisto];
	

	TH1F* muons_pt_histo[nhisto];
	TH1F* muons_eta_histo[nhisto];
	TH1F* muons_phi_histo[nhisto];



	TH1F* jets_pt_histo[nhisto];
	TH1F* jets_eta_histo[nhisto];
	TH1F* jets_Deta_histo[nhisto];
	TH1F* jets_phi_histo[nhisto];

	TH1F* H_pt_histo[nhisto];
	TH1F* H_eta_histo[nhisto];
	TH1F* H_phi_histo[nhisto];
	TH1F* Z_pt_histo[nhisto];
	TH1F* Z_eta_histo[nhisto];
	TH1F* Z_phi_histo[nhisto];

	

	TH1F* met_pt_histo[nhisto];
	TH1F* met_phi_histo[nhisto];




	for(int i=0; i<nhisto; i++)

	{	
		photons_pt_histo[i] = new TH1F(("photons_pt_histo"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 80, 0, 300 );
		photons_eta_histo[i] = new TH1F(("photons_eta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 25, -2.5, 2.5  );
		photons_phi_histo[i] = new TH1F(("photons_phi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 25, -3.15, 3.15 );
		photons_MVAID_histo[i] = new TH1F(("photons_MVAID_histo"+ std::to_string(i)).c_str(),"; Photon IDMVA; Counts", 25, -1, 1   );
		
		extra_leptons_n[i] = new TH1F(("extra_leptons_number"+ std::to_string(i)).c_str(), "; # of leptons; Counts", 4, -0.5, 3.5 );
		extra_leptons_pt_histo[i] = new TH1F(("extra_leptons_pt"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 20, 0, 40 );

		electrons_pt_histo[i] = new TH1F(("electrons_pt_histo"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 80, 0, 150 );
		electrons_eta_histo[i] = new TH1F(("electrons_eta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 25, -2.5, 2.5  );
		electrons_phi_histo[i] = new TH1F(("electrons_phi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 25, -3.15, 3.15 );
		electrons_MVAID_histo[i] = new TH1F(("electrons_MVAID_histo"+ std::to_string(i)).c_str(),"; Photon IDMVA; Counts", 25, -1, 1   );
		

		muons_pt_histo[i] = new TH1F(("muons_pt_histo"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 50, 0, 150 );
		muons_eta_histo[i] = new TH1F(("muons_eta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 25, -2.5, 2.5  );
		muons_phi_histo[i] = new TH1F(("muons_phi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 25, -3.15, 3.15 );

		jets_pt_histo[i] = new TH1F(("jets_pt_histo"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 80, 0, 300 );
		jets_eta_histo[i] = new TH1F(("jets_eta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 25, -2.5, 2.5  );
		jets_Deta_histo[i] = new TH1F(("jets_Deta_histo"+ std::to_string(i)).c_str(), "; #Delta#eta; Counts", 20, 0, 3.5  );
		jets_phi_histo[i] = new TH1F(("jets_phi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 25, -3.15, 3.15 );

		H_mass_all_histo[i] = new TH1F(("H_mass_all_histo"+  std::to_string(i)).c_str(), "; Invariant mass H (GeV); Counts", 80, 100, 180 );
		H_pt_histo[i] = new TH1F(("H_pt_histo"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 50, 0, 150 );
		H_eta_histo[i] = new TH1F(("H_eta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 25, -2.5, 2.5  );
		H_phi_histo[i] = new TH1F(("H_phi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 25, -3.15, 3.15 );

		Z_pt_histo[i] = new TH1F(("Z_pt_histo"+ std::to_string(i)).c_str(), "; P_{T}(GeV); Counts", 50, 0, 150 );
		Z_eta_histo[i] = new TH1F(("Z_eta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 25, -2.5, 2.5  );
		Z_phi_histo[i] = new TH1F(("Z_phi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 25, -3.15, 3.15 );
		Z_mass_all_histo[i] = new TH1F(("Z_mass_all__histo"+  std::to_string(i)).c_str(), "; Invariant mass Z (GeV); Counts", 80, 70, 110 );
		
		met_pt_histo[i] = new TH1F(("met_pt_histo"+  std::to_string(i)).c_str(), "; P_{T}^{MET} (GeV); Counts", 50, 0, 200 );
		met_phi_histo[i] = new TH1F(("met_phi_histo"+  std::to_string(i)).c_str(), "; #varphi^{MET}; Counts", 25, -3.15, 3.15 );
		

		for(int j=0; j<7; j++)	
                {
		 H_mass_histo[i][j] = new TH1F(("H_mass_histo"+  std::to_string(i)+  std::to_string(j)).c_str(), "; Invariant mass H (GeV); Counts", 80, 100, 180 );
		 Z_mass_histo[i][j] = new TH1F(("Z_mass_histo"+  std::to_string(i)+  std::to_string(j)).c_str(), "; Invariant mass Z (GeV); Counts", 80, 70, 110 );
                }

	}

	for(int n=0; n<nhisto; n++)
	{
		int nentries;
		TChain* serviceTree;
		switch (n)
		{
			case(0):
				nentries = ttH -> GetEntries();
				serviceTree = (TChain*)ttH -> Clone();
				break;


			case(1):
				nentries = ggH -> GetEntries();
				serviceTree = (TChain*)ggH -> Clone();
				break;


			case(2):
				nentries = vbf -> GetEntries();
				serviceTree = (TChain*)vbf -> Clone();
				break;


			case(3):
				nentries = wph -> GetEntries();
				serviceTree = (TChain*)wph -> Clone();
				break;




			case(4):
				nentries =  wmh-> GetEntries();
				serviceTree = (TChain*)wmh -> Clone();
				break;


			case(5):
				nentries = zh -> GetEntries();
				serviceTree = (TChain*)zh -> Clone();
				break;
			case(6):
				nentries = Zgamma -> GetEntries();
				serviceTree = (TChain*)Zgamma -> Clone();
				break;


			case(7):
				nentries = DY -> GetEntries();
				serviceTree = (TChain*)DY -> Clone();
				break;

			case(8):
				nentries = data -> GetEntries();
				serviceTree = (TChain*)data -> Clone();
				break;
		}


		serviceTree -> SetBranchAddress("weight_MC", &weight_MC);
		serviceTree -> SetBranchAddress("cat_n", &cat_n);
		serviceTree -> SetBranchAddress("gamma_pt", &gamma_pt);
		serviceTree -> SetBranchAddress("gamma_eta", &gamma_eta);
		serviceTree -> SetBranchAddress("gamma_phi", &gamma_phi);		
		serviceTree -> SetBranchAddress("gamma_IDMVA", &gamma_IDMVA);
	
		
		serviceTree -> SetBranchAddress("ele1_pt", &ele1_pt);
		serviceTree -> SetBranchAddress("ele3_pt", &ele3_pt);
		serviceTree -> SetBranchAddress("ele1_eta", &ele1_eta);
		serviceTree -> SetBranchAddress("ele1_phi", &ele1_phi);		
		serviceTree -> SetBranchAddress("ele1_IDMVA", &ele1_IDMVA);
	

		serviceTree -> SetBranchAddress("mu1_pt", &mu1_pt);
		serviceTree -> SetBranchAddress("mu3_pt", &mu3_pt);
		serviceTree -> SetBranchAddress("mu1_eta", &mu1_eta);
		serviceTree -> SetBranchAddress("mu1_phi", &mu1_phi);		

		
		serviceTree -> SetBranchAddress("jet1_pt", &jet1_pt);
		serviceTree -> SetBranchAddress("jet1_eta", &jet1_eta);
		serviceTree -> SetBranchAddress("jet2_eta", &jet2_eta);
		serviceTree -> SetBranchAddress("jet1_phi", &jet1_phi);		


		serviceTree -> SetBranchAddress("met_pt", &met_pt);
		serviceTree -> SetBranchAddress("met_phi", &met_phi);

		serviceTree -> SetBranchAddress("H_pt", &H_pt);
		serviceTree -> SetBranchAddress("H_eta", &H_eta);
		serviceTree -> SetBranchAddress("H_phi", &H_phi);
		serviceTree -> SetBranchAddress("H_mass", &H_mass);
		serviceTree -> SetBranchAddress("Z_mass", &Z_mass);
		serviceTree -> SetBranchAddress("Z_pt", &Z_pt);
		serviceTree -> SetBranchAddress("Z_eta", &Z_eta);
		serviceTree -> SetBranchAddress("Z_phi", &Z_phi);
		



		for(int i=0; i<nentries; i++)
		{	
			serviceTree -> GetEntry(i);

			n_lep_extra = 0;

			if(i%10000==0) cout << "Processing tag " << names[n] << ", event " << i << " out of " << nentries << "\r" << flush;			
			
			if (H_mass < 100 || H_mass > 180)  continue;
			if (n == 8) weight_MC= 1./lumiFactor; 
			
			Z_mass_histo[n][cat_n-1] -> Fill(Z_mass, weight_MC*lumiFactor);
			
			photons_pt_histo[n] -> Fill(gamma_pt, weight_MC*lumiFactor);
			photons_eta_histo[n] -> Fill(gamma_eta, weight_MC*lumiFactor);
			photons_phi_histo[n] -> Fill(gamma_phi, weight_MC*lumiFactor);
			photons_MVAID_histo[n] -> Fill(gamma_IDMVA, weight_MC*lumiFactor);

			if (isMuon)
			{
				muons_pt_histo[n] -> Fill(mu1_pt, weight_MC*lumiFactor);
				muons_eta_histo[n] -> Fill(mu1_eta, weight_MC*lumiFactor);
				muons_phi_histo[n] -> Fill(mu1_phi, weight_MC*lumiFactor);

			if (ele1_pt!=-100) 
			{
			n_lep_extra++;
			extra_leptons_pt_histo[n] -> Fill(ele1_pt, weight_MC*lumiFactor);
			}
			if (mu3_pt!=-100) 
			{
				n_lep_extra++;
				extra_leptons_pt_histo[n] -> Fill(mu3_pt, weight_MC*lumiFactor);
			}
				
			}
			else
			{
			electrons_pt_histo[n] -> Fill(ele1_pt, weight_MC*lumiFactor);
			electrons_eta_histo[n] -> Fill(ele1_eta, weight_MC*lumiFactor);
			electrons_phi_histo[n] -> Fill(ele1_phi, weight_MC*lumiFactor);
			electrons_MVAID_histo[n] -> Fill(ele1_IDMVA, weight_MC*lumiFactor);
			if (ele3_pt!=-100) 
			{
			n_lep_extra++;
			extra_leptons_pt_histo[n] -> Fill(ele3_pt, weight_MC*lumiFactor);
			}
			if (mu1_pt!=-100) 
			{
				n_lep_extra++;
				extra_leptons_pt_histo[n] -> Fill(mu1_pt, weight_MC*lumiFactor);
			}
			}

			if (cat_n==5)
			{
			jets_pt_histo[n] -> Fill(jet1_pt, weight_MC*lumiFactor);
			jets_eta_histo[n] -> Fill(jet1_eta, weight_MC*lumiFactor);
			jets_Deta_histo[n] -> Fill(fabs( jet1_eta - jet2_eta ), weight_MC*lumiFactor);
			jets_phi_histo[n] -> Fill(jet1_phi, weight_MC*lumiFactor);
			}
			Z_mass_histo[n][cat_n-1] -> Fill(Z_mass, weight_MC*lumiFactor);
			Z_mass_all_histo[n] -> Fill(Z_mass, weight_MC*lumiFactor);
			Z_pt_histo[n] -> Fill(Z_pt, weight_MC*lumiFactor);
			Z_eta_histo[n] -> Fill(Z_eta, weight_MC*lumiFactor);
			Z_phi_histo[n] -> Fill(Z_phi, weight_MC*lumiFactor);

			H_pt_histo[n] -> Fill(H_pt, weight_MC*lumiFactor);
			H_eta_histo[n] -> Fill(H_eta, weight_MC*lumiFactor);
			H_phi_histo[n] -> Fill(H_phi, weight_MC*lumiFactor);

			met_pt_histo[n] -> Fill(met_pt, weight_MC*lumiFactor);
			met_phi_histo[n] -> Fill(met_phi, weight_MC*lumiFactor);

                        if (n == 8 && (H_mass >115 && H_mass <135)) continue; 
			
			H_mass_histo[n][cat_n-1] -> Fill(H_mass, weight_MC*lumiFactor);
			H_mass_all_histo[n] -> Fill(H_mass, weight_MC*lumiFactor);

			extra_leptons_n[n] -> Fill(n_lep_extra, weight_MC*lumiFactor);
			

		}

		delete serviceTree;
		cout << "Processed tag " << names[n] << ", " << nentries << " events out of " << nentries << endl;

	}
	MakePlot(extra_leptons_n, "ExtraLeptons");
	MakePlot(extra_leptons_pt_histo, "ExtraLeptonsPt");

	MakePlot(Z_mass_all_histo, "ZetaInvariantMass");
	MakePlot(H_mass_all_histo, "HiggsInvariantMass");
	MakePlot(photons_pt_histo, "PhotonsPt");
	MakePlot(photons_eta_histo, "PhotonsEta");
	MakePlot(photons_phi_histo, "PhotonsPhi");	
	MakePlot(photons_MVAID_histo, "PhotonsMVAID");

	MakePlot(Z_pt_histo, "ZetaPt");
	MakePlot(Z_eta_histo, "ZetaEta");
	MakePlot(Z_phi_histo, "ZetaPhi");

	MakePlot(H_pt_histo, "HiggsPt");
	MakePlot(H_eta_histo, "HiggsEta");
	MakePlot(H_phi_histo, "HiggsPhi");


	if (isMuon)
	{
		MakePlot(muons_pt_histo, "MuonsPt");
		MakePlot(muons_eta_histo, "MuonsEta");
		MakePlot(muons_phi_histo, "MuonsPhi");
	}
	else
	{
		MakePlot(electrons_pt_histo, "ElectronsPt");
		MakePlot(electrons_eta_histo, "ElectronsEta");
		MakePlot(electrons_phi_histo, "ElectronsPhi");
		MakePlot(electrons_MVAID_histo, "ElectronsMVAID");

	}
	
		
		
	MakePlot(jets_pt_histo, "JetsPt");
	MakePlot(jets_eta_histo, "JetsEta");
	MakePlot(jets_Deta_histo, "JetsDeltaEta");
	MakePlot(jets_phi_histo, "JetsPhi");
		
		


	MakePlot(met_pt_histo, "met_pt");
	MakePlot(met_phi_histo, "met_phi");



     
	
		for(int i=0; i<7; i++)//loop cat
		{

			TH1F* tmp_Z[nhisto];
			TH1F* tmp_H[nhisto];
		

			for(int j=0; j<nhisto; j++)
			{	
			
				tmp_Z[j] = Z_mass_histo[j][i];
				tmp_H[j] = H_mass_histo[j][i];
				
			}

			
			MakePlot(tmp_Z, ("ZetaInvariantMass_cat"+ to_string(i)).c_str());
			MakePlot(tmp_H, ("HiggsInvariantMass_cat"+ to_string(i)).c_str());
		
		}
	if(isMuon)
	{	system("mv *.png /eos/user/f/fcetorel/www/Plot_Zgamma/muon/");
		system("mv *.pdf /eos/user/f/fcetorel/www/Plot_Zgamma/muon/");
	}
	else
	{	
		system("mv *.png /eos/user/f/fcetorel/www/Plot_Zgamma/ele/");
		system("mv *.pdf /eos/user/f/fcetorel/www/Plot_Zgamma/ele/");
	}


}
float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

void MakePlot(TH1F** histos, TString title)
{
	//float area=0;
	histos[0] -> SetLineWidth(3);			//ttH
	histos[0] -> SetLineColor(kRed + 1);
	histos[0] -> SetFillStyle(0);

	histos[1] -> SetLineWidth(3);			//ggH
	histos[1] -> SetLineColor(kGreen + 2);
	histos[1] -> SetFillStyle(0);

	histos[2] -> SetLineWidth(3);			//VBF
	histos[2] -> SetLineColor(kAzure);
	histos[2] -> SetFillStyle(0);

	histos[3] -> SetLineWidth(3);			//w+h
	histos[3] -> SetLineColor(kViolet - 2);
	histos[3] -> SetFillStyle(0);

	histos[4] -> SetLineWidth(3);			//w-h
	histos[4] -> SetLineColor(kOrange);
	histos[4] -> SetFillStyle(0);

	histos[5] -> SetLineWidth(3);			//zh
	histos[5] -> SetLineColor(kAzure + 8);
	histos[5] -> SetFillStyle(0);

	histos[6] -> SetLineWidth(3);			//Zgamma
	histos[6] -> SetLineColor(kOrange);
	histos[6] -> SetFillStyle(0);

	histos[7] -> SetLineWidth(3);			//DY
	histos[7] -> SetLineColor(kAzure + 8);
	histos[7] -> SetFillStyle(0);

	histos[8] -> SetMarkerStyle(20);		//Data
	histos[8] -> SetMarkerSize(1);
	histos[8] -> SetMarkerColor(kBlack);
	histos[8] -> SetFillStyle(0);



	TLegend* leg = new TLegend(0.65, 0.70, 0.9, 0.85);
	leg -> AddEntry(histos[0], "ttH", "l");
	leg -> AddEntry(histos[1], "ggH", "l");
	leg -> AddEntry(histos[2], "VBF", "l");
	leg -> AddEntry(histos[3], "W+H", "l");
	leg -> AddEntry(histos[4], "W-H", "l");
	leg -> AddEntry(histos[5], "ZH", "l");
	leg -> AddEntry(histos[8], "Data", "p");

	//leg -> AddEntry(histos[7], "Data sidebands 17", "p");


	TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
	//leg2 -> AddEntry(histos[0], "signal x50", "l");
	leg2 -> AddEntry(histos[6], "Zgamma", "l");
	leg2 -> AddEntry(histos[7], "DY", "l");
	leg2 -> AddEntry(histos[8], "Data ", "p");



	TCanvas* c = new TCanvas();
	c -> cd();

	float m = max(max(histos[0]->GetMaximum()/histos[0]->Integral(), histos[1]->GetMaximum()/histos[1]->Integral()), max(histos[2]->GetMaximum()/histos[2]->Integral(), histos[3]->GetMaximum()/histos[3]->Integral()));
	float m2 = max(histos[4]->GetMaximum()/histos[4]->Integral(), histos[5]->GetMaximum()/histos[5]->Integral());
	m = max(m,m2);

	TH1F* axis = new TH1F(*histos[0]);
	axis -> SetMarkerSize(0);
	axis -> SetLineWidth(0);
	axis -> GetYaxis() -> SetTitleOffset(1.5);
	axis -> GetYaxis() -> SetRangeUser(0, 1.1*m);

	axis -> Draw("histo");
	for (int i=0; i<5; i++) 
	{
		histos[i] -> DrawNormalized("histo SAME");
	}
	histos[8] -> DrawNormalized("E1 SAME");

	leg -> Draw("SAME");
	
	CMS_lumi(c, 0, 0);
	c -> SaveAs(title + "Signal.png");
	c -> SaveAs(title + "Signal.pdf");
	


	TCanvas* c2 = new TCanvas();
	c2->cd();
	/*for (int i=1; i<6; i++)
	{
	histos[0]->Add(histos[i]);
	}
	histos[0]->Scale(50);*/
	histos[7] -> Add(histos[6]);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.38, 1, 0.97);
	pad1->SetBottomMargin(0.035);
	pad1->Draw();
	pad1->cd();

	double max2 = histos[7]->GetBinContent(histos[7]->GetMaximumBin());
	if(histos[8]->GetBinContent(histos[8]->GetMaximumBin()) > max2 )
		max2 = histos[8]->GetBinContent(histos[8]->GetMaximumBin());
	


	histos[7]-> GetYaxis() -> SetRangeUser(0, 1.2*max2);
	histos[7]-> GetYaxis() -> SetTitleSize(0.05);
	histos[7]-> GetYaxis() -> SetTitleFont(42);
	histos[7]-> GetYaxis() -> SetLabelSize(0.045);
	histos[7]-> GetXaxis() -> SetLabelSize(0);
	histos[7]-> GetYaxis() -> SetLabelFont(42);


	histos[7] -> Draw("histo");
	//histos[0] -> Draw("histo same");
	
	histos[6] -> Draw("histo SAME");
	histos[8] -> Draw("E1 SAMe");

	
	leg2 -> Draw("SAME");



	c2->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.38);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->SetGridy(1);
	pad2->Draw();
	pad2->cd();

	
	
	TH1F *h = (TH1F*)histos[8]->Clone("h");
	h->SetLineColor(kBlack);
	h->SetMarkerColor(kBlue);
	h->SetMinimum(0.);  // Define Y ..
	h->SetMaximum(2); // .. range
	h->Sumw2();
	h->SetStats(0);      // No statistics on lower plot
	h->Divide(histos[7]);
	h->SetMarkerStyle(21);
	h -> SetTitle("");
	h-> GetYaxis() -> SetTitle("Data/MC");
	
	TF1* line = new TF1("line", "1", -100,300);
	line->SetLineColor(kRed);
	
	

	// Y axis ratio plot settings
	//h->GetYaxis()->SetNdivisions(-10);
	h->GetYaxis()->SetTitleSize(0.09);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleOffset(0.7);
	h->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
	h->GetYaxis()->SetLabelSize(0.09);

	// X axis ratio plot settings
	h->GetXaxis()->SetTitleSize(0.09);
	h->GetXaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
	h->GetXaxis()->SetLabelSize(0.12);
	h->GetXaxis()->SetTitleOffset(1.);
	h->Draw("EP");
	line->Draw("SAME");

	CMS_lumi(c2, 0, 0);


	c2 -> SaveAs(title + "Bkg.png");
	c2 -> SaveAs(title + "Bkg.pdf");
	

	delete c;
	delete leg;
	delete c2;
	delete leg2;


	

	return;
	}





