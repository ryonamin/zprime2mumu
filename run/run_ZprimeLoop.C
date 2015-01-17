{
  // Input files
  // Note that if you specify multiple input files,
  // you will merge them and get one output file.
  vector<string> inputs;
  //inputs.push_back("/Users/ryonamin/tmp/Zptest/root/ZprimeSSMToMuMu_M-2500_14TeV-pythia6_GEM2019Upg14DR-final_phase1_PU50bx25.root");
  //inputs.push_back("/Users/ryonamin/tmp/Zptest/root/ZprimeToMuMu_M-1000_Tune4C_13TeV-pythia8_tsg_PU40bx25_POSTLS162_V2-v1.root");
  inputs.push_back("/Users/ryonamin/tmp/Zptest/root/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8_PU40bx25_POSTLS162_V2-v1.root");

  // Load additional packages
  gSystem->Load("Util");          // Helper 

  // Load your macro with compilation
  gROOT->LoadMacro("ZprimeLoop.C++");

  // Run it
  ZprimeLoop m(inputs,"Output_ZprimeLoop.root"); 
  m.Loop();
 
}
