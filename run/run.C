{
  // Input files
  const string inputs[] = {
    "/Users/ryonamin/tmp/Zptest/root/ZprimeSSMToMuMu_M-2500_14TeV-pythia6_GEM2019Upg14DR-final_phase1_PU50bx25.root",
    "/Users/ryonamin/tmp/Zptest/root/ZprimeToMuMu_M-1000_Tune4C_13TeV-pythia8_tsg_PU40bx25_POSTLS162_V2-v1.root",
    "/Users/ryonamin/tmp/Zptest/root/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8_PU40bx25_POSTLS162_V2-v1.root"
  };

  // Load additional packages
  gSystem->Load("TreeHandler"); // TTree wrapper
  gSystem->Load("Util");        // Helper
  gSystem->Load("MCParticleFinder"); // Matching Reco-Gen 
  gSystem->Load("MCDimuonReco");     // reco dimuon for MC

  // Load your macro with compilation
  //gROOT->LoadMacro("MyAnalyzer_skeleton.C+");
  gROOT->LoadMacro("MyAnalyzer.C+");

  // Run it
  int nentries = sizeof(inputs) / sizeof(string);
  for (int ientry = 0; ientry < nentries; ientry++ ) { 
    MyAnalyzer(inputs[ientry]); 
  }
}
