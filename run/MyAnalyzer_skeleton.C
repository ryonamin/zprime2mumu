// Add include file for compilation
// === ROOT === 
#include "TFile.h"
#include "TTree.h"
// === STL ===
#include <iostream>
// === My Package ==
#include "TreeHandler.h" // Hold TTree 

using namespace std;

void MyAnalyzer_skeleton(const string fname = "") {
  TFile fin(fname.c_str());
  // Get IIHETree object.
  TreeHandler T(static_cast<TTree*>(fin.Get("IIHEAnalysis")));
  
  std::cout << "Processing " << fname << " ..." << std::endl;
  std::cout << "  Number of events = " << T.GetEntries() << std::endl;

  // Event Loop
  for ( int ev = 0; ev < T.GetEntries(); ev++ ) {
    // Write
  }

  // Release the memory before destruct TFile object.
  // Otherwise we will get a segmentation fault.
  //T = 0;
}
