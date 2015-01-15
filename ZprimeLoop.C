#define ZprimeLoop_cxx
#include "ZprimeLoop.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std ;

void ZprimeLoop::Loop()
{
   //   In a ROOT session, you can do:
   //      Root > .L ZprimeLoop.C
   //      Root > ZprimeLoop t
   //      Root > t.GetEntry(12); // Fill t data members with entry number 12
   //      Root > t.Show();       // Show values of entry 12
   //      Root > t.Show(16);     // Read and show values of entry 16
   //      Root > t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   //by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   TLorentzVector gen_muon1  ;
   TLorentzVector reco_muon1  ;
   TLorentzVector gen_muon2  ;
   TLorentzVector reco_muon2  ;
   TLorentzVector gen_Zprime ;
   TLorentzVector reco_Zprime ;
   TH1F* h_mc_ZpMass = new TH1F("", "", 500, -100, 100) ;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      if(jentry % 1000 == 0) cout << jentry*100./nentries << "%" << endl ;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // If there's no muon, this loop is useless
      if(mc_n == 0) cout << "No muon in this event." << endl ;
      if(mc_n < 2) continue ;

      // generated level
      int big_pT_index = 0 ;
      int bigger_pT_index = 0 ; 
      double big_pT = 0. ; 
      double bigger_pT = 0. ;

      for(unsigned long i = 0; i < mc_pt->size(); i++) {
	 if(mc_status->at(i) != 1) continue ;
	 if(fabs(mc_pdgId->at(i)) != 13) continue ;
	 if(mc_pt->at(i) > bigger_pT) {
	    bigger_pT = mc_pt->at(i) ;
	    bigger_pT_index = i ;
	 }
      }
      for(unsigned long i = 0; i < mc_pt->size(); i++) {
	 if(mc_status->at(i) != 1) continue ;
	 if(fabs(mc_pdgId->at(i)) != 13) continue ;
	 if(i == bigger_pT_index) continue ;
	 if(mc_pt->at(i) > big_pT) {
	    big_pT = mc_pt->at(i) ;
	    big_pT_index = i ;
	 }
      }
      
      gen_muon1.SetPxPyPzE(mc_px->at(bigger_pT_index), mc_py->at(bigger_pT_index), mc_pz->at(bigger_pT_index), mc_energy->at(bigger_pT_index)) ;
      gen_muon2.SetPxPyPzE(mc_px->at(big_pT_index), mc_py->at(big_pT_index), mc_pz->at(big_pT_index), mc_energy->at(big_pT_index)) ;
      gen_Zprime = gen_muon1 + gen_muon2 ;

      // reconstructed level
      big_pT_index = 0 ;
      bigger_pT_index = 0 ; 
      big_pT = 0. ; 
      bigger_pT = 0. ;

      if(muon_pt->size() < 2) continue ; 

      for(unsigned long i = 0; i < muon_pt->size(); i++) {
	 if(muon_pt->at(i) > bigger_pT) {
	    bigger_pT = muon_pt->at(i) ;
	    bigger_pT_index = i ;
	 }
      }
      for(unsigned long i = 0; i < muon_pt->size(); i++) {
	 if(i == bigger_pT_index) continue ;
	 if(muon_pt->at(i) > big_pT) {
	    big_pT = muon_pt->at(i) ;
	    big_pT_index = i ;
	 }
      }

      //cout << big_pT_index << "," << bigger_pT_index << endl ;
      reco_muon1.SetPtEtaPhiM(muon_pt->at(bigger_pT_index), muon_eta->at(bigger_pT_index), muon_phi->at(bigger_pT_index), 105.6583715e-3) ;
      reco_muon2.SetPtEtaPhiM(muon_pt->at(big_pT_index), muon_eta->at(big_pT_index), muon_phi->at(big_pT_index), 105.6583715e-3) ;

      reco_Zprime = reco_muon1 + reco_muon2 ;


      h_mc_ZpMass->Fill(gen_Zprime.M() - reco_Zprime.M()) ; 
   }

   h_mc_ZpMass->Draw() ;

}
