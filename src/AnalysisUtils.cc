#include "interface/AnalysisUtils.h"



float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) + 
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}



bool IsMatching(const particle& p1, const particle p2,
                const float& DRMax, const float& ptRatioMax)
{
  float DR = DeltaR(p1.v.Eta(),p1.v.Phi(),p2.v.Eta(),p2.v.Phi());

  bool isMatching = true;
  if( DR > DRMax ) isMatching = false;
  if( fabs(p1.v.Pt()/p2.v.Pt()-1.) > ptRatioMax ) isMatching = false;  

  return isMatching;
}

bool IsMatching(const std::vector<particle>& vp1, const std::vector<particle>& vp2,
                const float& DRMax, const float& ptRatioMax,
                const bool& verbosity)
{
  if( vp1.size() != vp2.size() )
  {
    std::cerr << "!!! IsMatching()::ERROR: the two vectors should have the same size !!!" << std::endl;
    return false;
  }
  
  std::vector<int> trial;
  std::vector<int> bestMatch;
  for(unsigned int ii = 0; ii < vp1.size(); ++ii)
  {
    trial.push_back( ii );
    bestMatch.push_back(-1);
  }
  
  float DRSumMin = 999999.;
  
  do {
    float DRSum = 0.;
    
    if( verbosity ) std::cout << "trying ";
    for(unsigned int ii = 0; ii < vp2.size(); ++ii)
    {
      float DR = DeltaR(vp1.at(trial.at(ii)).v.Eta(),vp1.at(trial.at(ii)).v.Phi(),vp2.at(ii).v.Eta(),vp2.at(ii).v.Phi());
      DRSum += DR;
      if( verbosity ) std::cout << trial.at(ii) << "-" << ii << " (" << DR << ")   ";
    }
    if( verbosity ) std::cout << std::endl;

    if( DRSum < DRSumMin )
    {
      DRSumMin = DRSum;
      bestMatch = trial;
    }
  } while ( std::next_permutation(trial.begin(),trial.end()) );


  bool isMatching = true;
  for(unsigned int ii = 0; ii < bestMatch.size(); ++ii)
  {
    float DR = DeltaR(vp1.at(bestMatch.at(ii)).v.Eta(),vp1.at(bestMatch.at(ii)).v.Phi(),vp2.at(ii).v.Eta(),vp2.at(ii).v.Phi());
    if( DR > DRMax ) isMatching = false;
    if( fabs(vp1.at(bestMatch.at(ii)).v.Pt()/vp2.at(ii).v.Pt()-1.) > ptRatioMax ) isMatching = false;
  }

  if( verbosity ) std::cout << (isMatching ? "MATCHING" : "NOT MATCHING") << std::endl;
  return isMatching;
}



int GetBestMatch(const particle& p, std::vector<particle>& vec, std::vector<int>* vetoVec)
{
  int bestMatch = -1;
  float DRMin = 999999.;
  
  for(unsigned int ii = 0; ii < vec.size(); ++ii)
  {
    bool skip = false;
    if( vetoVec)
      for(unsigned int jj = 0; jj < vetoVec->size(); ++jj)
        if( int(ii) == vetoVec->at(jj) ) skip = true;
    if( skip ) continue;
    
    float DR = DeltaR(vec.at(ii).v.Eta(),vec.at(ii).v.Phi(),p.v.Eta(),p.v.Phi());
    if( DR < DRMin )
    {
      DRMin = DR;
      bestMatch = int(ii);
    }
  }

  return bestMatch;
}



void PrintEvent(std::vector<particle>& lep)
{
  particle Z;
  Z.charge = 0;
  for(unsigned int ii = 0; ii < lep.size(); ++ii)
  {
    Z.v += lep.at(ii).v;
    Z.charge += lep.at(ii).charge;
  }
  Z.pdgId = 23;
  
  std::cout << "Z: " << Z << std::endl;
  for(unsigned int ii = 0; ii < lep.size(); ++ii)
    std::cout << ">>>>>>   lep: " << lep.at(ii) << std::endl;
}



void PrintEvent(std::vector<particle>& lep, const particle& gamma)
{
  particle Z;
  Z.charge = 0;
  for(unsigned int ii = 0; ii < lep.size(); ++ii)
  {
    Z.v += lep.at(ii).v;
    Z.charge += lep.at(ii).charge;
  }
  Z.pdgId = 23;
  
  particle H;
  H.charge = 0;
  H.v = Z.v + gamma.v;
  H.charge = 0;
  H.pdgId = 25;
  
  std::cout << "Z: " << Z << std::endl;
  for(unsigned int ii = 0; ii < lep.size(); ++ii)
    std::cout << ">>>>>>   lep: " << lep.at(ii) << std::endl;
}



std::ostream& operator<<(std::ostream& os, const TLorentzVector& v)
{
  os << "(pt="    << std::fixed << std::setw(6) << std::setprecision(2) << v.Pt()
     << ", eta="  << std::fixed << std::setw(5) << std::setprecision(2) << v.Eta()
     << ", phi="  << std::fixed << std::setw(5) << std::setprecision(2) << v.Phi()
     << ", mass=" << std::fixed << std::setw(7) << std::setprecision(3) << v.M()
     << ")";

  return os;
}

std::ostream& operator<<(std::ostream& os, const particle& p)
{
  os << p.v
     << "   charge: " << std::fixed << std::setw(2) << p.charge
     << "   pdgId: "  << std::fixed << std::setw(6) << p.pdgId
     << " (" << GetParticleName(p.pdgId) << ")";
  
  return os;
}



double ComputeSignificance(TH1F* h_sig, TH1F* h_bkg, const int& mode)
{
  double significance = 0.;
  
  for(int bin = 1; bin <= h_sig->GetNbinsX(); ++bin)
  {
    float S = h_sig->GetBinContent(bin);
    float B = h_bkg->GetBinContent(bin);
    
    if( B > 0 )
    {
      if( mode == 1 ) significance += S*S / B;
      if( mode == 2 ) significance += 2.*((S+B)*log(1.+S/B)-S);
    }
  }
  
  return sqrt(significance);
}
