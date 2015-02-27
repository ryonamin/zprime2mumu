#ifndef __TREEHANDLER_H__
#define __TREEHANDLER_H__
#include "TChain.h"
#include <string>
#include <vector>

class TreeHandler : public TChain 
{
  public :
  TreeHandler(std::string treename, std::vector<std::string>& filelist);
  
  std::vector<float>* mc_pt;
  std::vector<float>* mc_phi;
  std::vector<float>* mc_eta;
  std::vector<float>* mc_mass; 
  std::vector<float>* mc_px; 
  std::vector<float>* mc_py;
  std::vector<float>* mc_pz;               
  std::vector<int>* mc_pdgId;                 
  std::vector<int>* mc_charge;
  std::vector<int>* mc_status;

  // muon
  std::vector<float>* mu_isolationR03_sumPt;
  std::vector<bool>*  mu_isGlobalMuon;
  std::vector<bool>*  mu_isStandAloneMuon;
  std::vector<bool>*  mu_isTrackerMuon;
  std::vector<bool>*  mu_isPFMuon;
  std::vector<bool>*  mu_isPFIsolationValid;
  std::vector<int>*   mu_numberOfMatchedStations;
  std::vector<int>*   mu_numberOfValidPixelHits;

  // global muon
  std::vector<float>* mu_gt_pt;
  std::vector<float>* mu_gt_ptError;
  std::vector<float>* mu_gt_phi;
  std::vector<float>* mu_gt_phiError;
  std::vector<float>* mu_gt_eta;
  std::vector<float>* mu_gt_etaError;
  std::vector<int>*   mu_gt_charge;
  std::vector<float>* mu_gt_px;
  std::vector<float>* mu_gt_py;           
  std::vector<float>* mu_gt_pz;           
  std::vector<float>* mu_gt_dz;
  std::vector<float>* mu_gt_dz_beamspot;
  std::vector<float>* mu_gt_dz_firstPVtx;
  std::vector<float>* mu_gt_dxy;
  std::vector<float>* mu_gt_dxy_beamspot;
  std::vector<float>* mu_gt_dxy_firstPVtx;

  std::vector<float>* mu_gt_normalizedChi2;

  // tev optimized muon
  std::vector<float>* mu_tevOptimized_pt;
  std::vector<float>* mu_tevOptimized_ptError;
  std::vector<float>* mu_tevOptimized_phi;
  std::vector<float>* mu_tevOptimized_phiError;
  std::vector<float>* mu_tevOptimized_eta;
  std::vector<float>* mu_tevOptimized_etaError;
  std::vector<int>*   mu_tevOptimized_charge;
  std::vector<float>* mu_tevOptimized_px;
  std::vector<float>* mu_tevOptimized_py;           
  std::vector<float>* mu_tevOptimized_pz;           
  std::vector<float>* mu_tevOptimized_dz;
  std::vector<float>* mu_tevOptimized_dz_beamSpot;
  std::vector<float>* mu_tevOptimized_dz_firstPVtx;
  std::vector<float>* mu_tevOptimized_dxy;
  std::vector<float>* mu_tevOptimized_dxy_beamSpot;
  std::vector<float>* mu_tevOptimized_dxy_firstPVtx;
  
  void DumpInputs() const;
  private :
  std::vector<std::string> flist;
  ClassDef(TreeHandler,1)
};
#endif
