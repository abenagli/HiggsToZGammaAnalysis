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
  DATA(float, mu1_eta)                \
  DATA(float, mu2_eta)                \
  DATA(float, mu1_phi)                \
  DATA(float, mu2_phi)                \
  DATA(float, mu1_dxy)                \
  DATA(float, mu2_dxy)                \
  DATA(float, mu1_dxyPull)            \
  DATA(float, mu2_dxyPull)            \
  DATA(float, mu1_dz)                 \
  DATA(float, mu2_dz)                 \
  DATA(float, mu1_dzPull)             \
  DATA(float, mu2_dzPull)             \
  DATA(float, mu1_relIso)             \
  DATA(float, mu2_relIso)             \
  DATA(int, mu1_isL)                  \
  DATA(int, mu2_isL)                  \
  DATA(int, mu1_isM)                  \
  DATA(int, mu2_isM)                  \
  DATA(int, mu1_isT)                  \
  DATA(int, mu2_isT)                  \
  DATA(float, mu_Deta)                \
  DATA(float, mu_Dphi)                \
  DATA(float, mu_DR)                  \
  DATA(float, ele1_pt)                 \
  DATA(float, ele2_pt)                 \
  DATA(float, ele1_eta)                \
  DATA(float, ele2_eta)                \
  DATA(float, ele1_phi)                \
  DATA(float, ele2_phi)                \
  DATA(float, ele1_IDMVA)             \
  DATA(float, ele2_IDMVA)             \
  DATA(float, ele_Deta)                \
  DATA(float, ele_Dphi)                \
  DATA(float, ele_DR)                  \
  DATA(float, gamma_pt)               \
  DATA(float, gamma_eta)               \
  DATA(float, gamma_phi)               \
  DATA(float, gamma_IDMVA)             \
  DATA(float, gamma_DR1)             \
  DATA(float, gamma_DR2)             \
  DATA(int, jets_all_n)               \
  DATA(int, jets_all_bTagL_n)         \
  DATA(int, jets_all_bTagM_n)         \
  DATA(int, jets_all_bTagT_n)         \
  DATA(int, jets_cen_n)               \
  DATA(int, jets_cen_bTagL_n)         \
  DATA(int, jets_cen_bTagM_n)         \
  DATA(int, jets_cen_bTagT_n)         \
  DATA(int, jets_fwd_n)               \
  DATA(int, jets_fwd_bTagL_n)         \
  DATA(int, jets_fwd_bTagM_n)         \
  DATA(int, jets_fwd_bTagT_n)         \
  DATA(float, jet1_all_pt)            \
  DATA(float, jet1_all_eta)           \
  DATA(float, jet1_all_phi)           \
  DATA(float, jet1_all_energy)        \
  DATA(float, jet2_all_pt)            \
  DATA(float, jet2_all_eta)           \
  DATA(float, jet2_all_phi)           \
  DATA(float, jet2_all_energy)        \
  DATA(float, jet_all_Deta)           \
  DATA(float, jet_all_Dphi)           \
  DATA(float, jet_all_mass)           \
  DATA(float, jet1_cen_pt)            \
  DATA(float, jet1_cen_eta)           \
  DATA(float, jet1_cen_phi)           \
  DATA(float, jet1_cen_energy)        \
  DATA(float, jet2_cen_pt)            \
  DATA(float, jet2_cen_eta)           \
  DATA(float, jet2_cen_phi)           \
  DATA(float, jet2_cen_energy)        \
  DATA(float, jet_cen_Deta)           \
  DATA(float, jet_cen_Dphi)           \
  DATA(float, jet_cen_mass)           \
  DATA(float, met_pt)                 \
  DATA(float, met_phi)                \
  DATA(float, met_sig)

#include "DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
