// Add include file for compilation
// === ROOT === 
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
// === STL ===
#include <iostream>
// === My Package ==
#include "TreeHandler.h" // Hold TTree 

using namespace std;

void MyAnalyzer_skeleton(vector<string>& fpaths, string outfilename) {
  TreeHandler T("IIHEAnalysis",fpaths);
  std::cout << "Processing  ..." << std::endl;
  T.DumpInputs();
  
  std::cout << "  Total number of events = " << T.GetEntries() << std::endl;

  TFile fout(outfilename.c_str(),"RECREATE");

  // Prepare output objects
  TH1F hPt("hReco_Pt","Reco Pt",1000,0,5000);

  // Event Loop
  for ( int ev = 0; ev < T.GetEntries(); ev++ ) {
    T.GetEntry(ev);
    // Write processes here
    // An example for loop of a vector
    for ( unsigned int ip = 0; ip < T.muon_pt->size(); ip++ ) {
      hPt.Fill(T.muon_pt->at(ip));
    }
  }

  // Write output file
  hPt.Write();
  fout.Write();
}
