#include "TreeHandler.h"
#include <iostream> 

ClassImp(TreeHandler)

TreeHandler::TreeHandler(std::string treename, 
      std::vector<std::string>& in) : TChain(treename.c_str()),
   flist(in),
   mc_pt(0),
   mc_phi(0),
   mc_eta(0),
   mc_mass(0),
   mc_px(0),
   mc_py(0),
   mc_pz(0),
   mc_pdgId(0),
   mc_charge(0),
   mc_status(0),
   mu_gt_pt(0),
   mu_gt_ptError(0),
   mu_gt_phi(0),
   mu_gt_phiError(0),
   mu_gt_theta(0),
   mu_gt_thetaError(0),
   mu_tevOptimized_pt(0),
   mu_tevOptimized_ptError(0),
   mu_tevOptimized_phi(0),
   mu_tevOptimized_phiError(0),
   mu_tevOptimized_eta(0),
   mu_tevOptimized_etaError(0),
   mu_tevOptimized_charge(0),
   mu_tevOptimized_px(0),
   mu_tevOptimized_py(0),
   mu_tevOptimized_pz(0)
{
   for (unsigned int i = 0; i < flist.size(); i++ ) {
      Add(flist[i].c_str());
   }
   SetBranchAddress("mc_pt",&mc_pt);
   SetBranchAddress("mc_phi",&mc_phi);
   SetBranchAddress("mc_eta",&mc_eta);
   SetBranchAddress("mc_mass"           , &mc_mass) ;
   SetBranchAddress("mc_px"             , &mc_px) ;
   SetBranchAddress("mc_py"             , &mc_py) ;
   SetBranchAddress("mc_pz"             , &mc_pz) ;
   SetBranchAddress("mc_pdgId"          , &mc_pdgId) ;
   SetBranchAddress("mc_charge"         , &mc_charge) ;
   SetBranchAddress("mc_status"         , &mc_status) ;
   SetBranchAddress("mu_gt_pt"          , &mu_gt_pt) ;
   SetBranchAddress("mu_gt_ptError"     , &mu_gt_ptError) ;
   SetBranchAddress("mu_gt_phi"         , &mu_gt_phi) ;
   SetBranchAddress("mu_gt_phiError"    , &mu_gt_phiError) ;
   SetBranchAddress("mu_gt_theta"       , &mu_gt_theta) ;
   SetBranchAddress("mu_gt_thetaError"  , &mu_gt_thetaError) ;
   SetBranchAddress("mu_tevOptimized_pt",&mu_tevOptimized_pt) ;
   SetBranchAddress("mu_tevOptimized_ptError",&mu_tevOptimized_ptError) ;
   SetBranchAddress("mu_tevOptimized_phi",&mu_tevOptimized_phi) ;
   SetBranchAddress("mu_tevOptimized_phiError",&mu_tevOptimized_phiError) ;
   SetBranchAddress("mu_tevOptimized_eta",&mu_tevOptimized_eta) ;
   SetBranchAddress("mu_tevOptimized_etaError",&mu_tevOptimized_etaError) ;
   SetBranchAddress("mu_tevOptimized_charge",&mu_tevOptimized_charge) ;
   SetBranchAddress("mu_tevOptimized_px",&mu_tevOptimized_px) ;
   SetBranchAddress("mu_tevOptimized_py",&mu_tevOptimized_py) ;
   SetBranchAddress("mu_tevOptimized_pz",&mu_tevOptimized_pz) ;
}

void TreeHandler::DumpInputs() const {
   std::cout << "======== Input Files =======" << std::endl;
   for (unsigned int i = 0; i < flist.size(); i++ ) {
      std::cout << " " << i+1 << ") " << flist[i] << std::endl;
   }
   std::cout << "============================" << std::endl;
}
