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
			                                 muon_pt(0),
			                                 muon_ptError(0),
			                                 muon_phi(0),
			                                 muon_phiError(0),
			                                 muon_theta(0),
			                                 muon_thetaError(0),
			                                 muon_tevOptimized_pt(0),
			                                 muon_tevOptimized_ptError(0),
			                                 muon_tevOptimized_phi(0),
			                                 muon_tevOptimized_phiError(0),
			                                 muon_tevOptimized_eta(0),
			                                 muon_tevOptimized_etaError(0),
			                                 muon_tevOptimized_charge(0),
			                                 muon_tevOptimized_px(0),
			                                 muon_tevOptimized_py(0),
			                                 muon_tevOptimized_pz(0)
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
  SetBranchAddress("muon_pt"            , &muon_pt) ;
  SetBranchAddress("muon_ptError"       , &muon_ptError) ;
  SetBranchAddress("muon_phi"           , &muon_phi) ;
  SetBranchAddress("muon_phiError"      , &muon_phiError) ;
  SetBranchAddress("muon_theta"         , &muon_theta) ;
  SetBranchAddress("muon_thetaError"    , &muon_thetaError) ;
  SetBranchAddress("muon_tevOptimized_pt",&muon_tevOptimized_pt) ;
  SetBranchAddress("muon_tevOptimized_ptError",&muon_tevOptimized_ptError) ;
  SetBranchAddress("muon_tevOptimized_phi",&muon_tevOptimized_phi) ;
  SetBranchAddress("muon_tevOptimized_phiError",&muon_tevOptimized_phiError) ;
  SetBranchAddress("muon_tevOptimized_eta",&muon_tevOptimized_eta) ;
  SetBranchAddress("muon_tevOptimized_etaError",&muon_tevOptimized_etaError) ;
  SetBranchAddress("muon_tevOptimized_charge",&muon_tevOptimized_charge) ;
  SetBranchAddress("muon_tevOptimized_px",&muon_tevOptimized_px) ;
  SetBranchAddress("muon_tevOptimized_py",&muon_tevOptimized_py) ;
  SetBranchAddress("muon_tevOptimized_pz",&muon_tevOptimized_pz) ;
}

void TreeHandler::DumpInputs() const {
  std::cout << "======== Input Files =======" << std::endl;
  for (unsigned int i = 0; i < flist.size(); i++ ) {
    std::cout << " " << i+1 << ") " << flist[i] << std::endl;
  }
  std::cout << "============================" << std::endl;
}
