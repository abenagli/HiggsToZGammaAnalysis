#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "interface/TreeUtils.h"
#include "interface/ParticleNames.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#define PI 3.14159265359



struct particle
{
  TLorentzVector v;
  int charge;
  int pdgId;
  int it;
};



struct PtSort
{
  inline bool operator()(const particle& p1, const particle& p2)
  {
    return (p1.v.Pt() > p2.v.Pt());
  }
};



struct EtaSort
{
  inline bool operator()(const particle& p1, const particle& p2)
  {
    return (fabs(p1.v.Eta()) > fabs(p2.v.Eta()));
  }
};



struct NameSort
{
  inline bool operator()(const particle& p1, const particle& p2)
  {
    return (GetParticleName(p1.pdgId) < GetParticleName(p2.pdgId));
  }
};



struct FindPair
{
FindPair(const std::string& key)
: m_key(key) {}
  std::string m_key;
  bool operator()
    ( const std::pair<std::string,int>& p )
    {
      return( p.first == m_key);
    }
};



struct PairSort
{
  inline bool operator()(const std::pair<std::string,int>& p1, const std::pair<std::string,int>& p2)
  {
    return (p1.second > p2.second);
  }
};



float DeltaEta(const float& eta1, const float& eta2);
float DeltaPhi(const float& phi1, const float& phi2);
float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);

bool IsMatching(const particle& p1, const particle p2,
                const float& DRMax = 0.05, const float& ptRatioMax = 0.1);
bool IsMatching(const std::vector<particle>& vp1, const std::vector<particle>& vp2,
                const float& DRMax = 0.05, const float& ptRatioMax = 0.1,
                const bool& verbosity = false);

int GetBestMatch(const particle& p, std::vector<particle>& vec, std::vector<int>* vetoVec = NULL);


void PrintEvent(std::vector<particle>& lep);
void PrintEvent(std::vector<particle>& lep, const particle& gamma);

std::ostream& operator<<(std::ostream& os, const TLorentzVector& v);
std::ostream& operator<<(std::ostream& os, const particle& p);

double ComputeSignificance(TH1F* h_sig, TH1F* h_bkg, const int& mode = 1);

#endif
