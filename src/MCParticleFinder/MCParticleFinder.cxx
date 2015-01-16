#include "MCParticleFinder.h"
#include "TreeHandler.h"
#include <iostream>

ClassImp(MCParticleFinder)

MCParticleFinder::MCParticleFinder(TreeHandler& in) : fT(in),
                                                      fnsig(5.)
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
    if ( (*fT.mc_charge)[ip] != (*fT.muon_tevOptimized_charge)[recoId] ) continue;
    if ( (*fT.mc_pt)[ip] < (*fT.muon_tevOptimized_pt)[recoId] - fnsig*(*fT.muon_tevOptimized_ptError)[recoId] ) continue;
    if ( (*fT.mc_pt)[ip] > (*fT.muon_tevOptimized_pt)[recoId] + fnsig*(*fT.muon_tevOptimized_ptError)[recoId] ) continue;
    if ( (*fT.mc_phi)[ip] < (*fT.muon_tevOptimized_phi)[recoId] - fnsig*(*fT.muon_tevOptimized_phiError)[recoId] ) continue;
    if ( (*fT.mc_phi)[ip] > (*fT.muon_tevOptimized_phi)[recoId] + fnsig*(*fT.muon_tevOptimized_phiError)[recoId] ) continue;
    if ( (*fT.mc_eta)[ip] < (*fT.muon_tevOptimized_eta)[recoId] - fnsig*(*fT.muon_tevOptimized_etaError)[recoId] ) continue;
    if ( (*fT.mc_eta)[ip] > (*fT.muon_tevOptimized_eta)[recoId] + fnsig*(*fT.muon_tevOptimized_etaError)[recoId] ) continue;
    nMatched.push_back(ip);
  }
  if (nMatched.size()==1) return nMatched[0];
  else if (nMatched.size() > 1 ) return -1;
  else  return -2;
}

