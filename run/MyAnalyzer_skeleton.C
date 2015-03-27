// Add include file for compilation
// === ROOT === 
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
// === STL ===
#include <iostream>
// === My Package ==
#include "TreeHandler.h" // Hold TTree 
#include "MCParticleFinder.h"

using namespace std;

void MyAnalyzer_skeleton(vector<string>& fpaths, string outfilename) {
   TreeHandler T("IIHEAnalysis",fpaths);
   std::cout << "Processing  ..." << std::endl;
   T.DumpInputs();

   std::cout << "  Total number of events = " << T.GetEntries() << std::endl;

   TFile fout(outfilename.c_str(),"RECREATE");

   // Prepare output objects
   TH1F hPt("hReco_Pt","Reco Pt",1000,0,5000);
   TH1F hPtError("hReco_PtError","Reco PtError",1000,0,5000);
   TH1F hPtRes("hReco_PtRes","Reco PtRes",200,-1,1);

   // for finding corresponding mc particles
   MCParticleFinder mcpfinder(T);
   // define how tightly you require for variable(pt,phi,eta) to be matched between mc and reco. 
   mcpfinder.setNSig_Pt(10);  // require 10 sigma for matching in pt
   mcpfinder.setNSig_Phi(10); // require 10 sigma for matching in phi 
   mcpfinder.setNSig_Eta(10); // require 10 sigma for matching in eta

   // Event Loop
   for ( int ev = 0; ev < T.GetEntries(); ev++ ) {
      T.GetEntry(ev);
      if(ev%1000 == 0) cout << "ev n. " << ev << endl;
      // Write processes here
      // An example for loop of a vector
      for ( unsigned int ip = 0; ip < T.mu_gt_pt->size(); ip++ ) {
	 if(!T.mu_isGlobalMuon->at(ip)) continue ;
	 //if(T.mu_gt_pt->at(ip)<30) continue ;
	 //cout << T.mu_gt_pt->at(ip) << "," << T.mu_gt_eta->at(ip) << "," << T.mu_gt_phi->at(ip) << endl; 
	 hPt.Fill(T.mu_gt_pt->at(ip));
	 hPtError.Fill(T.mu_gt_ptError->at(ip));
	 //std::cout << "Reco. Pt = " << T.mu_gt_pt->at(ip) << std::endl;

	 // find the matched mc particle
	 // mcId > 0 : only one candidate is found 
	 int mcId = mcpfinder.getMatchedMCId(ip); 
	 //if (mcId<0) std::cout << "No matched MC particle." << std::endl;
	 //else std::cout << "MC Pt = " << T.mc_pt->at(mcId) << std::endl;
	 if(mcId >= 0) {
	    //cout << "I exist" << endl ;
	    double ptRes = (-T.mc_pt->at(mcId)+T.mu_gt_pt->at(ip))/T.mc_pt->at(mcId) ;	 
	    hPtRes.Fill(ptRes);
	 }
      }

      //cout << "--" << endl;
      // Dealing with MC particles
      /*
	 for ( unsigned int ip = 0; ip < T.mc_pt->size(); ip++ ) {
	 if(TMath::Abs(T.mc_pdgId->at(ip)) != 13) continue ;
	 cout << T.mc_pt->at(ip) << "," << T.mc_eta->at(ip) << "," << T.mc_phi->at(ip) << endl; 
	 }
	 cout << endl;
	 */
      // Write output file
   }

   hPt.Write();
   fout.Write();

}
