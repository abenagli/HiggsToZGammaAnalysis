#ifndef _OUT__TREE
#define _OUT__TREE

#include "DynamicTTree/interface/DynamicTTreeBase.h"



//---Define the TTree branches
#define DYNAMIC_TREE_NAME OutTree

#define DATA_TABLE                    \
  DATA(float, weight)                 \
  DATA(float, weight_PU)              \
  DATA(float, weight_MC)              \
  DATA(int,   vtxs_n)                 \
  DATA(int,   cat_n)                 \
  DATA(float, rho_all)                \
  DATA(float, H_pt)                   \
  DATA(float, H_eta)                  \
  DATA(float, H_phi)                  \
  DATA(float, H_mass)                 \
  DATA(float, Z_pt)                   \
  DATA(float, Z_eta)                  \
  DATA(float, Z_phi)                  \
  DATA(float, Z_mass)                 \
  DATA(float, mu1_pt)                 \
  DATA(float, mu2_pt)                 \
  DATA(float, mu3_pt)                 \
  DATA(float, mu1_eta)                \
  DATA(float, mu2_eta)                \
  DATA(float, mu3_eta)                \
  DATA(float, mu1_phi)                \
  DATA(float, mu2_phi)                \
  DATA(float, mu3_phi)                \
  DATA(float, mu1_dxy)                \
  DATA(float, mu2_dxy)                \
  DATA(float, mu3_dxy)                \
  DATA(float, mu1_dxyPull)            \
  DATA(float, mu2_dxyPull)            \
  DATA(float, mu3_dxyPull)            \
  DATA(float, mu1_dz)                 \
  DATA(float, mu2_dz)                 \
  DATA(float, mu3_dz)                 \
  DATA(float, mu1_dzPull)             \
  DATA(float, mu2_dzPull)             \
  DATA(float, mu3_dzPull)             \
  DATA(float, mu1_relIso)             \
  DATA(float, mu2_relIso)             \
  DATA(float, mu3_relIso)             \
  DATA(int, mu1_isL)                  \
  DATA(int, mu2_isL)                  \
  DATA(int, mu3_isL)                  \
  DATA(int, mu1_isM)                  \
  DATA(int, mu2_isM)                  \
  DATA(int, mu3_isM)                  \
  DATA(int, mu1_isT)                  \
  DATA(int, mu2_isT)                  \
  DATA(int, mu3_isT)                  \
  DATA(float, mu_Deta)                \
  DATA(float, mu_Dphi)                \
  DATA(float, mu_DR)                  \
  DATA(float, ele1_pt)                 \
  DATA(float, ele2_pt)                 \
  DATA(float, ele3_pt)                 \
  DATA(float, ele1_eta)                \
  DATA(float, ele2_eta)                \
  DATA(float, ele3_eta)                \
  DATA(float, ele1_phi)                \
  DATA(float, ele2_phi)                \
  DATA(float, ele3_phi)                \
  DATA(float, ele1_IDMVA)             \
  DATA(float, ele2_IDMVA)             \
  DATA(float, ele3_IDMVA)             \
  DATA(float, ele1_full5x5_R9)             \
  DATA(float, ele2_full5x5_R9)             \
  DATA(float, ele3_full5x5_R9)             \
  DATA(float, ele1_full5x5_sieie)             \
  DATA(float, ele2_full5x5_sieie)             \
  DATA(float, ele3_full5x5_sieie)             \
  DATA(float, ele_Deta)                \
  DATA(float, ele_Dphi)                \
  DATA(float, ele_DR)                  \
  DATA(float, gamma_pt)               \
  DATA(float, gamma_eta)               \
  DATA(float, gamma_phi)               \
  DATA(float, gamma_IDMVA)             \
  DATA(float, gamma_full5x5_R9)             \
  DATA(float, gamma_full5x5_sieie)     \
  DATA(float, gamma_DR1)             \
  DATA(float, gamma_DR2)             \
  DATA(float, jet1_pt)            \
  DATA(float, jet1_eta)           \
  DATA(float, jet1_phi)           \
  DATA(float, jet1_energy)        \
  DATA(float, jet2_pt)            \
  DATA(float, jet2_eta)           \
  DATA(float, jet2_phi)           \
  DATA(float, jet2_energy)        \
  DATA(float, met_pt)                 \
  DATA(float, met_phi)                \
  DATA(float, met_sig)

#include "DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
