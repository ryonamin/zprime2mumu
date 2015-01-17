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
  std::vector<float>* muon_pt;
  std::vector<float>* muon_phi;
  std::vector<float>* muon_eta;
  std::vector<float>* muon_tevOptimized_pt;
  std::vector<float>* muon_tevOptimized_ptError;
  std::vector<float>* muon_tevOptimized_phi;
  std::vector<float>* muon_tevOptimized_phiError;
  std::vector<float>* muon_tevOptimized_eta;
  std::vector<float>* muon_tevOptimized_etaError;
  std::vector<int>*   muon_tevOptimized_charge;
  std::vector<float>* muon_tevOptimized_px;
  std::vector<float>* muon_tevOptimized_py;           
  std::vector<float>* muon_tevOptimized_pz;           
  
  void DumpInputs() const;
  private :
  std::vector<std::string> flist;
  ClassDef(TreeHandler,1)
};
#endif
