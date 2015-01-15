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

void MyAnalyzer_skeleton(const string fpath = "") {
  TFile fin(fpath.c_str());
  // Get IIHETree object.
  TreeHandler T(static_cast<TTree*>(fin.Get("IIHEAnalysis")));
  
  std::cout << "Processing " << fpath << " ..." << std::endl;
  std::cout << "  Number of events = " << T.GetEntries() << std::endl;

  // Get input file name from the path
  string ifname = fpath.substr(0,fpath.find(".root"));
  ifname = ifname.substr(ifname.find_last_of("/")+1,ifname.length());
  
  // Output root file
  string ofname = "Output_" + ifname + ".root";
  TFile fout(ofname.c_str(),"RECREATE");

  // Prepare output objects
  TH1F hPt("hReco_Pt","Reco Pt",1000,0,5000);

  // Event Loop
  for ( int ev = 0; ev < T.GetEntries(); ev++ ) {
    T.GetEntry(ev);
    // Write processes here
    // An example for loop of a vector
    for ( unsigned int ip = 0; ip < T.muon_tevOptimized_pt->size(); ip++ ) {
      hPt.Fill(T.muon_tevOptimized_pt->at(ip));
    }
  }

  // Write output file
  hPt.Write();
  fout.Write();
}
