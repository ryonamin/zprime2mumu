// ROOT
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1I.h"
// STL
#include <vector>
#include <iostream>
#include <sstream>
// Local
#include "TreeHandler.h"
#include "Util.h"
#include "MCParticleFinder.h"
#include "MCDimuonReco.h"

using namespace std;

void MyAnalyzer(vector<string>& fpaths, string outfilename) {
  TreeHandler T("IIHEAnalysis",fpaths);
  std::cout << "Processing  ..." << std::endl;
  T.DumpInputs();
  std::cout << "  Total umber of events = " << T.GetEntries() << std::endl;

  // Cut parameters 
  const float highPtThreshold = 200; // 40 GeV [HLT], 200 GeV [tevOptimiaze]

  // Output root file
  TFile fout(outfilename.c_str(),"RECREATE");

  MCDimuonReco mcdimuonreco(T);
  MCParticleFinder mcpfinder(T);
  mcpfinder.setNSig(10.); // require 5 sigma for matching in pt, phi, eta

  int nevents = T.GetEntries();

  TH1F* hMCMuonPhi            = new TH1F("hMCMuonPhi","Phi of MC Muons",100,-5,5);
  TH1F* hMCMuonEta            = new TH1F("hMCMuonEta","Eta of MC Muons",200,-10,10);
  TH1I* hMatchedMC            = new TH1I("hMatchedMC","# of Matched MC particles",2,0,2);
  TH1I* hNMatchedPerEvent     = new TH1I("hNMatchedPerEvent","# of Matched MC particles per Event",5,0,5);
  TH1I* hNHighPtMuonsPerEvent = new TH1I("hNHighPtMuonsPerEvent","# of High Pt Muons per Event",5,0,5);
  TH1F* hDimuonMass           = new TH1F("hDimuonMass","Invariant Mass",1000,0,5000);
  TH1F* hMatchedMCDimuonMass  = new TH1F("hMatchedMCDimuonMass","Invariant Mass (MC matched Reco)",1000,0,5000);
  TH1F* hDimuonDeltaMass      = new TH1F("hDimuonDeltaMass","Delta Mass (Reco-MC)",1000,-1000,1000);
  TH1F* hDimuonDeltaMassNorm  = new TH1F("hDimuonDeltaMassNorm","Delta Mass (Reco-MC)/MC",1000,-0.5,0.5);
  TH1F* hMCDimuonMass         = new TH1F("hMCDimuonMass","Invariant Mass (MC)",1000,0,5000);
  TH1F* hRecoMuonPhi          = new TH1F("hRecoMuonPhi","Phi of Reco Muons",100,-5,5);
  TH1F* hRecoMuonEta          = new TH1F("hRecoMuonEta","Eta of Reco Muons",200,-10,10);
  TH1F* hMatchedMuonPhi       = new TH1F("hMatchedMuonPhi","Phi of Matched-reco Muons",100,-5,5);
  TH1F* hMatchedMuonEta       = new TH1F("hMatchedMuonEta","Eta of Matched-reco Muons",200,-10,10);

  int nHighPtMuonsTotal = 0;
  int nHighPtMuonsWithMatchedMC = 0;

  for ( int ev = 0; ev < nevents; ev++ ) {
    T.GetEntry(ev);

    for ( unsigned int mp = 0; mp < T.mc_phi->size(); mp++) {
      if (TMath::Abs(T.mc_pdgId->at(mp))==13 && T.mc_status->at(mp) == 1 && T.mc_pt->at(mp) > highPtThreshold ) {
        hMCMuonPhi->Fill(T.mc_phi->at(mp));
        hMCMuonEta->Fill(T.mc_eta->at(mp));
      }
    }

    if ( mcdimuonreco.findHighPtDimuon() ) {
      hMCDimuonMass->Fill(mcdimuonreco.getDimuonMass());
    }

    int nMatched = 0;
    int nHighPtMuons = 0;
    float leadingPt = 0;
    float subleadingPt = 0;
    int leadingPtId = -1;
    int leadingPtMCId = -1;
    int subleadingPtId = -1;
    int subleadingPtMCId = -1;

    for ( unsigned int ip = 0; ip < T.muon_tevOptimized_pt->size(); ip++) {
      float pt = T.muon_tevOptimized_pt->at(ip);
      float phi = T.muon_tevOptimized_phi->at(ip);
      float eta = T.muon_tevOptimized_eta->at(ip);
      hRecoMuonPhi->Fill(phi);
      hRecoMuonEta->Fill(eta);
      int mcId = mcpfinder.getMatchedMCId(ip); 
      // mcId < 0 --> no matched MC particle.
      if (mcId > 0) {
	hMatchedMC->Fill(1);
	nMatched++;
        hMatchedMuonPhi->Fill(phi);
        hMatchedMuonEta->Fill(eta);
      } else {
        hMatchedMC->Fill(0);
	//cerr << "This is rejected : " << T.muon_tevOptimized_pt->at(ip) << " retval=" << mcId << endl;
      }
      if ( pt > highPtThreshold ) {
	nHighPtMuons++;
        nHighPtMuonsTotal++;
        if (mcId>0) nHighPtMuonsWithMatchedMC++;
      }

      if (leadingPt < pt ) {
	subleadingPt = leadingPt;
	subleadingPtId = leadingPtId;
	subleadingPtMCId = leadingPtMCId;
	leadingPt = pt;
	leadingPtId = ip;
	leadingPtMCId = mcId;
      } else if (leadingPt >= pt && subleadingPt < pt) {
	subleadingPt = pt;
	subleadingPtId = ip;
	subleadingPtMCId = mcId;
      } 
    }
    hNMatchedPerEvent->Fill(nMatched);
    hNHighPtMuonsPerEvent->Fill(nHighPtMuons);
    //cerr << "Leading Pt = " << leadingPt << " Sub-Leading Pt = " << subleadingPt << endl;
    if (leadingPtId>=0&&subleadingPtId>=0 && subleadingPt > highPtThreshold ) {
      TLorentzVector rp1,rp2;
      rp1.SetXYZM(T.muon_tevOptimized_px->at(leadingPtId),
                 T.muon_tevOptimized_py->at(leadingPtId),
                 T.muon_tevOptimized_pz->at(leadingPtId),
		 0.105); 
      rp2.SetXYZM(T.muon_tevOptimized_px->at(subleadingPtId),
                 T.muon_tevOptimized_py->at(subleadingPtId),
                 T.muon_tevOptimized_pz->at(subleadingPtId),
		 0.105); 
      TLorentzVector dimuon(rp1+rp2);
      hDimuonMass->Fill(dimuon.M());

      //cerr << leadingPtMCId << " " << subleadingPtMCId << endl; 
      if (leadingPtMCId>0&&subleadingPtMCId>0) {
        TLorentzVector mp1,mp2;
        mp1.SetXYZM(T.mc_px->at(leadingPtMCId),
                    T.mc_py->at(leadingPtMCId),
                    T.mc_pz->at(leadingPtMCId),
          	  T.mc_mass->at(leadingPtMCId)); 
        mp2.SetXYZM(T.mc_px->at(subleadingPtMCId),
                    T.mc_py->at(subleadingPtMCId),
                    T.mc_pz->at(subleadingPtMCId),
          	  T.mc_mass->at(subleadingPtMCId)); 
        TLorentzVector dimuonMC(mp1+mp2);
        hMatchedMCDimuonMass->Fill(dimuonMC.M());
	hDimuonDeltaMass->Fill(dimuon.M()-dimuonMC.M());
	hDimuonDeltaMassNorm->Fill( (dimuon.M()-dimuonMC.M())/dimuonMC.M() );
      }

    } else {
#if 0
      cerr << "Rejected reco particle pt1 = " << leadingPt << " pt2=" << subleadingPt  
           << "  pt1Id = " << leadingPtId << " pt2Id = " << subleadingPtId;
      if (leadingPtId>=0&&subleadingPtId>=0) {
      cerr << " eta1 = " << T.muon_tevOptimized_eta->at(leadingPtId)
	   << " eta2 = " << T.muon_tevOptimized_eta->at(subleadingPtId)
	   << " pz1 = " << T.muon_tevOptimized_pz->at(leadingPtId)
	   << " pz2 = " << T.muon_tevOptimized_pz->at(subleadingPtId);
      }
      cerr << endl;
#endif
    }
  }

  cerr << "Matching efficiency for high Pt Muons (> " << highPtThreshold << " GeV) : " 
       << float(nHighPtMuonsWithMatchedMC)/float(nHighPtMuonsTotal) << endl;

  TH1F* hEfficiencyPhi = static_cast<TH1F*>(hMatchedMuonPhi->Clone("hEfficiencyPhi"));
  hEfficiencyPhi->Sumw2();
  hEfficiencyPhi->Divide(hMCMuonPhi);
  hEfficiencyPhi->SetTitle("tevOptimized Muon Detection Efficiency in Phi");
  TH1F* hEfficiencyEta = static_cast<TH1F*>(hMatchedMuonEta->Clone("hEfficiencyEta"));
  hEfficiencyEta->Sumw2();
  hEfficiencyEta->Divide(hMCMuonEta);
  hEfficiencyEta->SetTitle("tevOptimized Muon Detection Efficiency in Eta");
  TH1F* hPurityPhi     = static_cast<TH1F*>(hMatchedMuonPhi->Clone("hPurityPhi"));
  hPurityPhi->Sumw2();
  hPurityPhi->Divide(hRecoMuonPhi);
  hPurityPhi->SetTitle("tevOptimized Muon Detection Purity in Phi");
  TH1F* hPurityEta     = static_cast<TH1F*>(hMatchedMuonEta->Clone("hPurityEta"));
  hPurityEta->Sumw2();
  hPurityEta->Divide(hRecoMuonEta);
  hPurityEta->SetTitle("tevOptimized Muon Detection Purity in Eta");

  // Write histgrams
  hMatchedMC->Write();
  hNMatchedPerEvent->Write();
  hNHighPtMuonsPerEvent->Write();
  hDimuonMass->Write();
  hMatchedMCDimuonMass->Write();
  hMCDimuonMass->Write();
  hDimuonDeltaMass->Write();
  hDimuonDeltaMassNorm->Write();
  hMCMuonPhi->Write();
  hMCMuonEta->Write();
  hRecoMuonPhi->Write();
  hRecoMuonEta->Write();
  hMatchedMuonPhi->Write();
  hMatchedMuonEta->Write();
  hEfficiencyPhi->Write();
  hEfficiencyEta->Write();
  hPurityPhi->Write();
  hPurityEta->Write();
  fout.Close();
}
