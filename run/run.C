{
  // Input files
  vector<string> inputs;
  inputs.push_back("/Users/fzenoni/zprime2mumu/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8_tsg_PU40bx25_POSTLS162_V2-v1.root");
  //inputs.push_back("/Users/ryonamin/tmp/Zptest/root/ZprimeToMuMu_M-1000_Tune4C_13TeV-pythia8_tsg_PU40bx25_POSTLS162_V2-v1.root");
  //inputs.push_back("/Users/ryonamin/tmp/Zptest/root/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8_PU40bx25_POSTLS162_V2-v1.root");
  //inputs.push_back("/Users/ryonamin/tmp/Zptest/root/v2/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8_tsg_PU40bx25_POSTLS162_V2-v1.root");

  // Load additional packages
  gSystem->Load("TreeHandler"); // TTree wrapper
  gSystem->Load("Util");        // Helper
  gSystem->Load("MCParticleFinder"); // Matching Reco-Gen 
  gSystem->Load("MCDimuonReco");     // reco dimuon for MC


  // Load your macro with compilation
  gROOT->LoadMacro("MyAnalyzer_skeleton.C+");
  //gROOT->LoadMacro("MyAnalyzer.C+");

  // Run it
  MyAnalyzer_skeleton(inputs,"Output.root"); 
  //MyAnalyzer(inputs,"Output.root"); 
}
