/*
 g++ -Wall -o cat `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore  cat.cpp

#ifndef CMS_LUMI_H
#include "CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "tdrstyle.h"
#endif
*/
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


void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity);

int main(int argc, char *argv[])
{
	/*writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_sqrtS = "";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);*/

	TString MuonFold = "/afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/mu/new/";

	TString EleFold = "/afs/cern.ch/work/f/fcetorel/private/work2/HtogZ/HiggsToZGammaAnalysis/plots/ele/new/";

	vector<TString> names = {"plots_mc_ttH", "plots_mc_ggH", "plots_mc_VBF", "plots_mc_WplusH","plots_mc_WminusH","plots_mc_ZH", "plots_mc_Zgamma", "plots_mc_DY", "plots_data"};

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


	/*	data = new TChain("out_tree");
		data -> Add(MuonFold + names[8] + "*.root");
		
	
		
		Zgamma = new TChain("out_tree"); 
		Zgamma -> Add(MuonFold + names[6] + ".root");
				
		DY = new TChain("out_tree");
		DY -> Add(MuonFold + names[7] + ".root");
	*/

	
	
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

/*
		data = new TChain("out_tree");
		data -> Add(EleFold + names[8] + "*.root");
		
	
		
		Zgamma = new TChain("out_tree"); 
		Zgamma -> Add(EleFold + names[6] + ".root");
				
		DY = new TChain("out_tree");
		DY -> Add(EleFold + names[7] + ".root");	*/	

	}

	int nhisto=9; //numero di 
	

	float weight_MC = 0;
	int cat_n = 0;

	
	float H_mass = 0.;
	float ele1_pt = 0.;
	float mu3_pt = 0.;



	
	TH1F* H_mass_histo[nhisto][7];
	TH1F* H_mass_tot_histo;


	H_mass_tot_histo = new TH1F("H_mass_histo_tot", "; Invariant mass H (GeV); Counts", 80, 100, 180 );

	for(int i=0; i<6; i++)
	{	
		for(int j=0; j<7; j++)	
                {
		 H_mass_histo[i][j] = new TH1F(("H_mass_histo"+  std::to_string(i)+  std::to_string(j)).c_str(), "; Invariant mass H (GeV); Counts", 80, 100, 180 );
                }

	}

	for(int n=0; n<6; n++)
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
		serviceTree -> SetBranchAddress("ele1_pt", &ele1_pt);
		serviceTree -> SetBranchAddress("mu3_pt", &mu3_pt);


		serviceTree -> SetBranchAddress("H_mass", &H_mass);
		



		for(int i=0; i<nentries; i++)
		{	
			serviceTree -> GetEntry(i);



			if(i%10000==0) cout << "Processing tag " << names[n] << ", event " << i << " out of " << nentries << "\r" << flush;			
			
			if (H_mass < 100 || H_mass > 180)  continue;
			if (cat_n == 6 && ele1_pt != -100) continue; 
			
			H_mass_histo[n][cat_n-1] -> Fill(H_mass, weight_MC*lumiFactor);
			H_mass_tot_histo -> Fill(H_mass, weight_MC*lumiFactor);




			

		}

		delete serviceTree;
		cout << "Processed tag " << names[n] << ", " << nentries << " events out of " << nentries << endl;

	}

       	for (int i=0; i<6; i++)
        {  
       		for (int j=0; j<7; j++)
		{
			float ret[4];
			const bool verbosity = true;  
			FindSmallestInterval(ret, H_mass_histo[i][j], 0.683, verbosity); 
			float nevents_sigma; 
			float min = H_mass_histo[i][j]->FindBin(ret[2]);
			float max = H_mass_histo[i][j]->FindBin(ret[3]);
			nevents_sigma = H_mass_histo[i][j]-> Integral();
			cout << nevents_sigma*0.683 << endl; 
		}
	}
			float nevents_sigma = H_mass_tot_histo-> Integral();
			cout << "totale  " << nevents_sigma*0.683 << endl;

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
/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction, const bool& verbosity)
{
    float integralMax = fraction * histo->Integral();
  
    int N = histo -> GetNbinsX();
    std::vector<float> binCenters(N);
    std::vector<float> binContents(N);
    std::vector<float> binIntegrals(N);
    for(int bin1 = 0; bin1 < N; ++bin1)
    {
        binCenters[bin1] = histo->GetBinCenter(bin1+1);
        binContents[bin1] = histo->GetBinContent(bin1+1);
    
        for(int bin2 = 0; bin2 <= bin1; ++bin2)
            binIntegrals[bin1] += binContents[bin2];
    }
  
    float min = 0.;
    float max = 0.;
    float delta = 999999.;
    for(int bin1 = 0; bin1 < N; ++bin1)
    {
        for(int bin2 = bin1+1; bin2 < N; ++bin2)
        {
            if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
            float tmpMin = histo -> GetBinCenter(bin1);
            float tmpMax = histo -> GetBinCenter(bin2);
      
            if( (tmpMax-tmpMin) < delta )
            {
                delta = (tmpMax - tmpMin);
                min = tmpMin;
                max = tmpMax;
            }
      
            break;
        }
    }
  
    TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
    for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
    {
        if( smallHisto->GetBinCenter(bin) < min )
            smallHisto -> SetBinContent(bin,0);
    
        if( smallHisto->GetBinCenter(bin) > max )
            smallHisto -> SetBinContent(bin,0);
    }
    smallHisto -> SetFillColor(kYellow);
  
    float mean = smallHisto -> GetMean();
    float meanErr = smallHisto -> GetMeanError();  
  
    ret[0] = mean;
    ret[1] = meanErr;
    ret[2] = min;
    ret[3] = max;
}




