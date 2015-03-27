#include "MCParticleFinder.h"
#include "TreeHandler.h"

#include <TLorentzVector.h>
#include <iostream>
#include <cmath>        // std::abs

using namespace std;

ClassImp(MCParticleFinder)

MCParticleFinder::MCParticleFinder(TreeHandler& in) : fT(in),
   fnsig_pt(5.),
   fnsig_phi(5.),
   fnsig_eta(5.)
{
}

MCParticleFinder::~MCParticleFinder(){
}

int MCParticleFinder::getMatchedMCId(int recoId)
{
   std::vector<int> nMatched;
   for ( int ip = 0; ip < fT.mc_pt->size(); ip++ ) {
      if ( (*fT.mc_status)[ip] != 1 ) continue; // only check the final state particles
      if ( TMath::Abs((*fT.mc_pdgId)[ip]) != 13 ) continue;
      TLorentzVector *lv_mc = new TLorentzVector() ;
      lv_mc->SetXYZM((*fT.mc_px)[ip], (*fT.mc_py)[ip], (*fT.mc_pz)[ip], 105.6583715e-3) ;
      TLorentzVector *lv_gt = new TLorentzVector() ;
      lv_gt->SetXYZM((*fT.mu_gt_px)[recoId], (*fT.mu_gt_py)[recoId], (*fT.mu_gt_pz)[recoId], 105.6583715e-3) ;
      //if ( (*fT.mc_charge)[ip] != (*fT.mu_gt_charge)[recoId] ) continue;
      if ( (*fT.mc_pt)[ip] < 500) continue ;
      if ( (*fT.mu_gt_pt)[recoId] < 30) continue ;
      if ( std::abs((*fT.mc_eta)[ip]) > 1) continue ;
      if ( std::abs(lv_gt->Eta()) > 10) continue ;
      if ( (lv_mc->Phi()/lv_gt->Phi() > 1.2 || lv_mc->Phi()/lv_gt->Phi() < 0.8) ) continue ;
      if ( (lv_mc->Eta()/lv_gt->Eta() > 1.2 || lv_mc->Eta()/lv_gt->Eta() < 0.8) ) continue ;
      //if ( (*fT.mc_pt)[ip] < (*fT.mu_gt_pt)[recoId] - fnsig_pt*(*fT.mu_gt_ptError)[recoId] ) continue;
      //if ( (*fT.mc_pt)[ip] > (*fT.mu_gt_pt)[recoId] + fnsig_pt*(*fT.mu_gt_ptError)[recoId] ) continue;
      //if ( (*fT.mc_phi)[ip] < (*fT.mu_gt_px)[recoId]/(*fT.mu_gt_pt)[recoId] - fnsig_phi*(*fT.mu_gt_phiError)[recoId] ) continue;
      //if ( (*fT.mc_phi)[ip] > (*fT.mu_gt_px)[recoId]/(*fT.mu_gt_pt)[recoId] + fnsig_phi*(*fT.mu_gt_phiError)[recoId] ) continue;
      //if ( (*fT.mc_eta)[ip] < (*fT.mu_gt_eta)[recoId] - fnsig_eta*(*fT.mu_gt_etaError)[recoId] ) continue;
      //if ( (*fT.mc_eta)[ip] > (*fT.mu_gt_eta)[recoId] + fnsig_eta*(*fT.mu_gt_etaError)[recoId] ) continue;
      //if ( (*fT.mc_eta)[ip] < (*fT.mu_gt_eta)[recoId] - fnsig_eta*(*fT.mu_gt_etaError)[recoId] ) continue;
      //if ( (*fT.mc_phi)[ip]/(*fT.mu_gt_px)[recoId]/(*fT.mu_gt_py)[recoId] - fnsig_phi*(*fT.mu_gt_phiError)[recoId] ) continue;
      nMatched.push_back(ip);
   }
   if (nMatched.size()==1) return nMatched[0];
   else if (nMatched.size() > 1 ) return -1;
   else  return -2;
}

