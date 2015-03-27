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
  std::vector<float>* mu_gt_pt;
  std::vector<float>* mu_gt_px;
  std::vector<float>* mu_gt_py;
  std::vector<float>* mu_gt_pz;
  std::vector<float>* mu_gt_ptError;
  std::vector<float>* mu_gt_phi;
  std::vector<float>* mu_gt_phiError;
  std::vector<float>* mu_gt_theta;
  std::vector<float>* mu_gt_thetaError;
  std::vector<float>* mu_gt_eta;
  std::vector<float>* mu_gt_etaError;
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
  std::vector<bool>* mu_isGlobalMuon;           
  
  void DumpInputs() const;
  private :
  std::vector<std::string> flist;
  ClassDef(TreeHandler,1)
};
#endif
