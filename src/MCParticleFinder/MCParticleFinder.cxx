#include "MCParticleFinder.h"
#include "TreeHandler.h"
#include <iostream>

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
    if ( (*fT.mc_charge)[ip] != (*fT.mu_tevOptimized_charge)[recoId] ) continue;
    if ( (*fT.mc_pt)[ip] < (*fT.mu_tevOptimized_pt)[recoId] - fnsig_pt*(*fT.mu_tevOptimized_ptError)[recoId] ) continue;
    if ( (*fT.mc_pt)[ip] > (*fT.mu_tevOptimized_pt)[recoId] + fnsig_pt*(*fT.mu_tevOptimized_ptError)[recoId] ) continue;
    if ( (*fT.mc_phi)[ip] < (*fT.mu_tevOptimized_phi)[recoId] - fnsig_phi*(*fT.mu_tevOptimized_phiError)[recoId] ) continue;
    if ( (*fT.mc_phi)[ip] > (*fT.mu_tevOptimized_phi)[recoId] + fnsig_phi*(*fT.mu_tevOptimized_phiError)[recoId] ) continue;
    if ( (*fT.mc_eta)[ip] < (*fT.mu_tevOptimized_eta)[recoId] - fnsig_eta*(*fT.mu_tevOptimized_etaError)[recoId] ) continue;
    if ( (*fT.mc_eta)[ip] > (*fT.mu_tevOptimized_eta)[recoId] + fnsig_eta*(*fT.mu_tevOptimized_etaError)[recoId] ) continue;
    nMatched.push_back(ip);
  }
  if (nMatched.size()==1) return nMatched[0];
  else if (nMatched.size() > 1 ) return -1;
  else  return -2;
}

