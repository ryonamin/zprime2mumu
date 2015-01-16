#include "TreeHandler.h"

ClassImp(TreeHandler)

TreeHandler::TreeHandler(TTree* in) : fT(in),
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
  fT->SetBranchAddress("mc_pt",&mc_pt);
  fT->SetBranchAddress("mc_phi",&mc_phi);
  fT->SetBranchAddress("mc_eta",&mc_eta);
  fT->SetBranchAddress("mc_mass"           , &mc_mass) ;
  fT->SetBranchAddress("mc_px"             , &mc_px) ;
  fT->SetBranchAddress("mc_py"             , &mc_py) ;
  fT->SetBranchAddress("mc_pz"             , &mc_pz) ;
  fT->SetBranchAddress("mc_pdgId"          , &mc_pdgId) ;
  fT->SetBranchAddress("mc_charge"         , &mc_charge) ;
  fT->SetBranchAddress("mc_status"         , &mc_status) ;
  fT->SetBranchAddress("muon_tevOptimized_pt",&muon_tevOptimized_pt) ;
  fT->SetBranchAddress("muon_tevOptimized_ptError",&muon_tevOptimized_ptError) ;
  fT->SetBranchAddress("muon_tevOptimized_phi",&muon_tevOptimized_phi) ;
  fT->SetBranchAddress("muon_tevOptimized_phiError",&muon_tevOptimized_phiError) ;
  fT->SetBranchAddress("muon_tevOptimized_eta",&muon_tevOptimized_eta) ;
  fT->SetBranchAddress("muon_tevOptimized_etaError",&muon_tevOptimized_etaError) ;
  fT->SetBranchAddress("muon_tevOptimized_charge",&muon_tevOptimized_charge) ;
  fT->SetBranchAddress("muon_tevOptimized_px",&muon_tevOptimized_px) ;
  fT->SetBranchAddress("muon_tevOptimized_py",&muon_tevOptimized_py) ;
  fT->SetBranchAddress("muon_tevOptimized_pz",&muon_tevOptimized_pz) ;
}
