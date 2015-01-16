#include "TreeHandler.h"
#include "MCDimuonReco.h"
#include "TLorentzVector.h"
#include <iostream>
using namespace std;

ClassImp(MCDimuonReco)

MCDimuonReco::MCDimuonReco(TreeHandler& in) : fT(in),
                                              fmass(0.)
{
}

bool MCDimuonReco::findHighPtDimuon() 
{
  float leadingPt = 0;
  float subleadingPt = 0;
  int leadingPtId = -1;
  int subleadingPtId = -1;
  for ( int i = 0; i < fT.mc_pt->size(); i++ ) {
    if (fT.mc_status->at(i)!=1) continue; // Select final state particles 
    if (TMath::Abs(fT.mc_pdgId->at(i))!=13) continue; // Select muons
    float pt = fT.mc_pt->at(i);
    if ( leadingPt < pt ) {
      subleadingPt = leadingPt;
      subleadingPtId = leadingPtId;
      leadingPt = pt;
      leadingPtId = i;
    } else if ( leadingPt > pt && subleadingPt < pt ) {
      subleadingPt = pt;
      subleadingPtId = i;
    }
  } 

  if ( leadingPtId > 0 && subleadingPtId > 0 ) { // Find a high Pt dimuon pair
    TLorentzVector p1,p2;
    p1.SetXYZM( fT.mc_px->at(leadingPtId),
                fT.mc_py->at(leadingPtId),
                fT.mc_pz->at(leadingPtId),
                fT.mc_mass->at(leadingPtId) );
    p2.SetXYZM( fT.mc_px->at(subleadingPtId),
                fT.mc_py->at(subleadingPtId),
                fT.mc_pz->at(subleadingPtId),
                fT.mc_mass->at(subleadingPtId) );
    TLorentzVector dimuon(p1+p2);
    fmass = dimuon.M();
    return true;
  } else return false; 
  
}
