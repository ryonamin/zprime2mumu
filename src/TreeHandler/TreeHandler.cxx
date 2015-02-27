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
                                                         mu_isolationR03_sumPt(0),
							 mu_isGlobalMuon(0),
                                                         mu_isStandAloneMuon(0),
                                                         mu_isTrackerMuon(0),
                                                         mu_isPFMuon(0),
                                                         mu_isPFIsolationValid(0),
                                                         mu_numberOfMatchedStations(0),
                                                         mu_numberOfValidPixelHits(0),
			                                 mu_gt_pt(0),
			                                 mu_gt_ptError(0),
			                                 mu_gt_phi(0),
			                                 mu_gt_phiError(0),
			                                 mu_gt_eta(0),
			                                 mu_gt_etaError(0),
			                                 mu_gt_charge(0),
			                                 mu_gt_px(0),
			                                 mu_gt_py(0),
			                                 mu_gt_pz(0),
                                                         mu_gt_dz(0),
                                                         mu_gt_dz_beamspot(0),
                                                         mu_gt_dz_firstPVtx(0),
                                                         mu_gt_dxy(0),
                                                         mu_gt_dxy_beamspot(0),
                                                         mu_gt_dxy_firstPVtx(0),
                                                         mu_gt_normalizedChi2(0),
			                                 mu_tevOptimized_pt(0),
			                                 mu_tevOptimized_ptError(0),
			                                 mu_tevOptimized_phi(0),
			                                 mu_tevOptimized_phiError(0),
			                                 mu_tevOptimized_eta(0),
			                                 mu_tevOptimized_etaError(0),
			                                 mu_tevOptimized_charge(0),
			                                 mu_tevOptimized_px(0),
			                                 mu_tevOptimized_py(0),
			                                 mu_tevOptimized_pz(0),
                                                         mu_tevOptimized_dz(0),
                                                         mu_tevOptimized_dz_beamspot(0),
                                                         mu_tevOptimized_dz_firstPVtx(0),
                                                         mu_tevOptimized_dxy(0),
                                                         mu_tevOptimized_dxy_beamspot(0),
                                                         mu_tevOptimized_dxy_firstPVtx(0)
{
  for (unsigned int i = 0; i < flist.size(); i++ ) {
    Add(flist[i].c_str());
  }
  // mc
  SetBranchAddress("mc_pt"                         , &mc_pt);
  SetBranchAddress("mc_phi"                        , &mc_phi);
  SetBranchAddress("mc_eta"                        , &mc_eta);
  SetBranchAddress("mc_mass"                       , &mc_mass) ;
  SetBranchAddress("mc_px"                         , &mc_px) ;
  SetBranchAddress("mc_py"                         , &mc_py) ;
  SetBranchAddress("mc_pz"                         , &mc_pz) ;
  SetBranchAddress("mc_pdgId"                      , &mc_pdgId) ;
  SetBranchAddress("mc_charge"                     , &mc_charge) ;
  SetBranchAddress("mc_status"                     , &mc_status) ;
  // muon
  SetBranchAddress("mu_isolationR03_sumPt"         , &mu_isolationR03_sumPt) ;
  SetBranchAddress("mu_isGlobalMuon"               , &mu_isGlobalMuon           ) ; 
  SetBranchAddress("mu_isStandAloneMuon"           , &mu_isStandAloneMuon       ) ;  
  SetBranchAddress("mu_isTrackerMuon"              , &mu_isTrackerMuon          ) ; 
  SetBranchAddress("mu_isPFMuon"                   , &mu_isPFMuon               ) ;  
  SetBranchAddress("mu_isPFIsolationValid"         , &mu_isPFIsolationValid     ) ;      
  SetBranchAddress("mu_numberOfMatchedStations"    , &mu_numberOfMatchedStations) ; 
  SetBranchAddress("mu_numberOfValidPixelHits"     , &mu_numberOfValidPixelHits ) ; 
  // global muon
  SetBranchAddress("mu_gt_pt"                      , &mu_gt_pt) ;
  SetBranchAddress("mu_gt_ptError"                 , &mu_gt_ptError) ;
  SetBranchAddress("mu_gt_phi"                     , &mu_gt_phi) ;
  SetBranchAddress("mu_gt_phiError"                , &mu_gt_phiError) ;
  SetBranchAddress("mu_gt_eta"                     , &mu_gt_eta) ;
  SetBranchAddress("mu_gt_etaError"                , &mu_gt_etaError) ;
  SetBranchAddress("mu_gt_charge"                  , &mu_gt_charge) ;
  SetBranchAddress("mu_gt_px"                      , &mu_gt_px) ;
  SetBranchAddress("mu_gt_py"                      , &mu_gt_py) ;
  SetBranchAddress("mu_gt_pz"                      , &mu_gt_pz) ;
  SetBranchAddress("mu_gt_dz"                      , &mu_gt_dz            ) ;  
  SetBranchAddress("mu_gt_dz_beamspot"             , &mu_gt_dz_beamspot   ) ;   
  SetBranchAddress("mu_gt_dz_firstPVtx"            , &mu_gt_dz_firstPVtx  ) ; 
  SetBranchAddress("mu_gt_dxy"                     , &mu_gt_dxy           ) ;  
  SetBranchAddress("mu_gt_dxy_beamspot"            , &mu_gt_dxy_beamspot  ) ; 
  SetBranchAddress("mu_gt_dxy_firstPVtx"           , &mu_gt_dxy_firstPVtx ) ; 

  SetBranchAddress("mu_gt_normalizedChi2"          , &mu_gt_normalizedChi2) ; 
  // tev optimized muon
  SetBranchAddress("mu_tevOptimized_pt"            , &mu_tevOptimized_pt) ;
  SetBranchAddress("mu_tevOptimized_ptError"       , &mu_tevOptimized_ptError) ;
  SetBranchAddress("mu_tevOptimized_phi"           , &mu_tevOptimized_phi) ;
  SetBranchAddress("mu_tevOptimized_phiError"      , &mu_tevOptimized_phiError) ;
  SetBranchAddress("mu_tevOptimized_eta"           , &mu_tevOptimized_eta) ;
  SetBranchAddress("mu_tevOptimized_etaError"      , &mu_tevOptimized_etaError) ;
  SetBranchAddress("mu_tevOptimized_charge"        , &mu_tevOptimized_charge) ;
  SetBranchAddress("mu_tevOptimized_px"            , &mu_tevOptimized_px) ;
  SetBranchAddress("mu_tevOptimized_py"            , &mu_tevOptimized_py) ;
  SetBranchAddress("mu_tevOptimized_pz"            , &mu_tevOptimized_pz) ;
  SetBranchAddress("mu_tevOptimized_dz"            , &mu_tevOptimized_dz            ) ;  
  SetBranchAddress("mu_tevOptimized_dz_beamspot"   , &mu_tevOptimized_dz_beamspot   ) ;   
  SetBranchAddress("mu_tevOptimized_dz_firstPVtx"  , &mu_tevOptimized_dz_firstPVtx  ) ; 
  SetBranchAddress("mu_tevOptimized_dxy"           , &mu_tevOptimized_dxy           ) ;  
  SetBranchAddress("mu_tevOptimized_dxy_beamspot"  , &mu_tevOptimized_dxy_beamspot  ) ; 
  SetBranchAddress("mu_tevOptimized_dxy_firstPVtx" , &mu_tevOptimized_dxy_firstPVtx ) ; 
}

void TreeHandler::DumpInputs() const {
  std::cout << "======== Input Files =======" << std::endl;
  for (unsigned int i = 0; i < flist.size(); i++ ) {
    std::cout << " " << i+1 << ") " << flist[i] << std::endl;
  }
  std::cout << "============================" << std::endl;
}
