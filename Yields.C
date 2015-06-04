#define Yields_cxx
#include "Yields.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TF1.h>
#include <TGraphErrors.h>

using namespace std;

void Yields::Loop()
{
   //   In a ROOT session, you can do:
   //      Root > .L Yields.C
   //      Root > Yields t
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


   TFile * myfile = new TFile("out.root","recreate");

   //TH1D * h1d_pt_barrel = new TH1D("h1d_pt_barrel","",100,-1,1);
   //TH2D * h2d_pt_barrel = new TH2D("h2d_pt_barrel","",100,-1,1,10,0,2600);  
   //TH1D * h1d_pt_endcaps = new TH1D("h1d_pt_endcaps","",100,-1,1);
   //TH2D * h2d_pt_endcaps = new TH2D("h2d_pt_endcaps","",100,-1,1,10,0,2600);  

   //TH1D*  h_inv_genmass = new TH1D("h_inv_genmass","",100,0,5200) ;
   //TH1D*  h_inv_recomass = new TH1D("h_inv_recomass","",100,0,5200) ;

   //TH1D * h1d_mass = new TH1D("h1d_mass","",100,-1,1);
   //TH2D * h2d_mass = new TH2D("h2d_mass","",100,-1,1,10,0,5200);

   TH1D * h_gen1_px = new TH1D("h_gen1_px","",100,-3000,3000);
   TH1D * h_gen1_py = new TH1D("h_gen1_py","",100,-3000,3000);
   TH1D * h_gen1_pz = new TH1D("h_gen1_pz","",100,-3000,3000);
   TH1D * h_gen1_e =  new TH1D("h_gen1_e" ,"",100,0,6000);

   TH1D * h_gen2_px = new TH1D("h_gen2_px","",100,-3000,3000);
   TH1D * h_gen2_py = new TH1D("h_gen2_py","",100,-3000,3000);
   TH1D * h_gen2_pz = new TH1D("h_gen2_pz","",100,-3000,3000);
   TH1D * h_gen2_e =  new TH1D("h_gen2_e" ,"",100,0,6000);

   TH1D * h_reco1_px = new TH1D("h_reco1_px","",100,-6000,6000);
   TH1D * h_reco1_py = new TH1D("h_reco1_py","",100,-6000,6000);
   TH1D * h_reco1_pz = new TH1D("h_reco1_pz","",100,-6000,6000);
   TH1D * h_reco1_e  = new TH1D("h_reco1_e" ,"",100,0,6000);

   TH1D * h_reco2_px = new TH1D("h_reco2_px","",100,-6000,6000);
   TH1D * h_reco2_py = new TH1D("h_reco2_py","",100,-6000,6000);
   TH1D * h_reco2_pz = new TH1D("h_reco2_pz","",100,-6000,6000);
   TH1D * h_reco2_e  = new TH1D("h_reco2_e" ,"",100,0,6000);

   TH1D * h_gen1_eta  = new TH1D("h_gen1_eta" ,"",100,-5,5);
   TH1D * h_gen2_eta  = new TH1D("h_gen2_eta" ,"",100,-5,5);
   TH1D * h_gen1_pt  = new TH1D("h_gen1_pt" ,"",100,0,3000);
   TH1D * h_gen2_pt  = new TH1D("h_gen2_pt" ,"",100,0,3000);

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   int no_reco = 0;
   int counter_gen[4] = {0};
   int counter_reco[5] = {0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if(jentry%1000 == 0) cout << jentry << " over " << nentries << endl;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);  nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      bool gem_reco = false;
      bool gem_gen = false;
      std::vector<TLorentzVector> genmu, recomu;
      TLorentzVector genmu_tmp, recomu_tmp;

      for(int i = 0; i < (int)mc_n; i++) {
	 if((*mc_status)[i] != 1) continue ; // final state (observable)
	 if(abs((*mc_pdgId)[i])!=13) continue; // we want muons
	 genmu_tmp.SetPxPyPzE((*mc_px)[i],(*mc_py)[i],(*mc_pz)[i],(*mc_energy)[i]);
	 genmu.push_back(genmu_tmp) ;
      }

      TLorentzVector genmu_1, genmu_2;
      if(genmu.size() > 1) { // there are at least two muons around
	 counter_gen[0]++ ; // all events
	 genmu_1 = genmu[0] ;
	 genmu_2 = genmu[1] ;
	 for(int i = 0; i < (int)genmu.size(); i++) {
	    if(genmu[i].Pt() > genmu_1.Pt() && genmu[i] != genmu_1) {
	       genmu_2 = genmu_1;
	       genmu_1 = genmu[i];
	    }
	    else if(genmu[i].Pt() > genmu_2.Pt() && genmu[i] != genmu_1) {
	       genmu_2 = genmu[i];
	    }
	 }
	 
	 h_gen1_eta->Fill(genmu_1.Eta());
	 h_gen2_eta->Fill(genmu_2.Eta());      
	 h_gen1_pt->Fill(genmu_1.Pt());
	 h_gen2_pt->Fill(genmu_2.Pt());
	 
	 if(abs(genmu_1.Eta()) < 2.4 && abs(genmu_2.Eta()) < 2.4) {
	    counter_gen[1]++ ; 
	    if(genmu_1.Pt() > 10 && genmu_2.Pt() > 10) {
	       h_gen1_px->Fill(genmu_1.Px());
	       h_gen1_py->Fill(genmu_1.Py());
	       h_gen1_pz->Fill(genmu_1.Pz());
	       h_gen1_e->Fill(genmu_1.P());
	       h_gen2_px->Fill(genmu_2.Px());
	       h_gen2_py->Fill(genmu_2.Py());
	       h_gen2_pz->Fill(genmu_2.Pz());
	       h_gen2_e->Fill(genmu_2.P());
	       counter_gen[2]++ ; 
	       if((abs(genmu_1.Eta()) > 1.5 && abs(genmu_1.Eta() < 2.2)) || (abs(genmu_2.Eta()) > 1.5 && abs(genmu_2.Eta()) < 2.2))
		  counter_gen[3]++ ;
	    }
	 }
      }

      //-------------------------

      for(int i = 0; i < (int)(*mu_gt_pt).size(); i++) {
	 if(!(*mu_isGlobalMuon)[i]) continue;
	 recomu_tmp.SetPxPyPzE((*mu_gt_px)[i],(*mu_gt_py)[i],(*mu_gt_pz)[i],(*mu_gt_p)[i]);
	 recomu.push_back(recomu_tmp) ;
      }

      int index1(0), index2(0);
      TLorentzVector recomu_1, recomu_2;
      if(recomu.size() > 1) { // there are at least two muons around
	 counter_reco[0]++ ; // all events
	 recomu_1 = recomu[0] ;
	 index1 = 0;
	 recomu_2 = recomu[1] ;
	 index2 = 1;
	 for(int i = 0; i < (int)recomu.size(); i++) {
	    if(recomu[i].Pt() > recomu_1.Pt() && recomu[i] != recomu_1) {
	       recomu_2 = recomu_1;
	       index2 = index1;
	       recomu_1 = recomu[i];
	       index1 = i;
	    }
	    else if(recomu[i].Pt() > recomu_2.Pt() && recomu[i] != recomu_1) {
	       recomu_2 = recomu[i];
	       index2 = i;
	    }
	 }
	 if(abs(recomu_1.Eta()) < 2.4 && abs(recomu_2.Eta() < 2.4)) {
	    counter_reco[1]++ ; 
	    if(recomu_1.Pt() > 10 && recomu_2.Pt() > 10) {
	        h_reco1_px->Fill(recomu_1.Px());
	        h_reco1_py->Fill(recomu_1.Py());
	        h_reco1_pz->Fill(recomu_1.Pz());
	        h_reco1_e->Fill(recomu_1.P());
	        h_reco2_px->Fill(recomu_2.Px());
	        h_reco2_py->Fill(recomu_2.Py());
	        h_reco2_pz->Fill(recomu_2.Pz());
	        h_reco2_e->Fill(recomu_2.P());
	        counter_reco[2]++ ; 
		if(PassHighPtSel(index1) && PassHighPtSel(index2)) counter_reco[3]++;
	        if((abs(recomu_1.Eta()) > 1.5 && abs(recomu_1.Eta()) < 2.2) || (abs(recomu_2.Eta()) > 1.5 && abs(recomu_2.Eta()) < 2.2)) {
		  counter_reco[4]++ ;
	       }
	    }
	 }
      }
   }

   cout << "G E N  L E V E L" << endl;
   cout << "All events:                 " << counter_gen[0] << endl;
   cout << "eta 2.4 acceptance:         " << counter_gen[1] << endl;
   cout << "10 GeV pT cut:              " << counter_gen[2] << endl;
   cout << "At least one muon in GE1/1: " << counter_gen[3] << endl;
   cout << "==================" << endl;
   cout << "R E C O  L E V E L" << endl;
   cout << "All events:                 " << counter_reco[0] << endl;
   cout << "eta 2.4 acceptance:         " << counter_reco[1] << endl;
   cout << "10 GeV pT cut:              " << counter_reco[2] << endl;
   cout << "Pass selection:             " << counter_reco[3] << endl;
   cout << "At least one muon in GE1/1: " << counter_reco[4] << endl;

   myfile->Write();
   myfile->Close();
}

bool Yields::PassHighPtSel(int n){
   //This is an old selection that needs to be updated... but better than nothing. Compare with selection in AN-15-061
   //Be careful with it as I am not sure I found the right variables in the IIHE tree. 
   //This tree is really a mess at the moment... Apparently All the muons (regardless of their type) are stored in a same vector and you end up with empty entries. 
   // Also some of the variables are labelled mu_gt_X, but others mu_X, which is confusing... 
   float muon_pt_min = 45.;
   float muon_dptOverPt = 0.3;
   float muon_etaMax = 2.4;

   int muon_nHitsMinPixel = 1;
   int muon_nHitsMinMuon = 1;
   int muon_nLayersMin = 6;
   float muon_impactParamMaxXY = 0.2;   // in cm                                                                                                   

   int muon_nSegMatchMin = 2;
   float muon_relIsoCutMax = 0.1;

   if ((*mu_gt_pt)[n] > muon_pt_min
	 && (*mu_gt_ptError)[n] / (*mu_gt_pt)[n] < muon_dptOverPt
	 //&& fabs((*mu_gt_eta)[n]) < muon_etaMax
	 && (*mu_numberOfValidPixelHits)[n] >= muon_nHitsMinPixel
	 //      && mu_gt_nhitsmuons[n] >= muon_nHitsMinMuon  I had to comment this line because I don't find the variable in the tree... 
	 // but I think this condition (above) is superseeded by mu_numberOfMatchedStations>=muon_nSegMatchMin      
	 //       && mu_gt_nlayerswithhits[n] >= muon_nLayersMin // I had to comment this line because I don't find the variable in the tree...
	 //I don't know what the impact is...
	 && fabs((*mu_gt_dxy_firstPVtx)[n]) < muon_impactParamMaxXY 
	 && (*mu_numberOfMatchedStations)[n] >= muon_nSegMatchMin// mu_numberOfMatchedStations
	 && (*mu_isTrackerMuon)[n]
	 && (*mu_isolationR03_sumPt)[n] / (*mu_gt_pt)[n] < muon_relIsoCutMax
      )
      return true;
   return false;
}


void Yields::Resolvs(TH2D * myh2d){
   TString hname=myh2d->GetName(); 
   const int nbbins = myh2d->GetNbinsY();
   //if(nbbins == 0) cout << "Danger" << endl;
   double resol[nbbins], xaxis[nbbins],errresol[nbbins],errxaxis[nbbins];
   //Loop over the mass/pt bins and fit a gaussian to each slice. 
   for(int i = 1 ; i <= myh2d->GetNbinsY(); i++){

      TH1D *myproj_x = myh2d->ProjectionX("_px",i,i);
      TH1D *myproj_y = myh2d->ProjectionY("_py",1,1);
      TF1 *f1 = new TF1("f1", "gaus", -0.3, 0.3);
      myproj_x->Fit(f1, "R");
      //Parameter 2 is sigma
      resol[i-1] = f1->GetParameter(2);
      errresol[i-1] = f1->GetParError(2);
      xaxis[i-1]= myproj_y->GetBinCenter(i);
      //For the x error, I take half of the bin
      errxaxis[i-1]= (myproj_y->GetBinLowEdge(i+1) -myproj_y->GetBinLowEdge(i))/2;

      TString suffix; 
      suffix.Form("%-.4g", myproj_y->GetBinLowEdge(i));
      myproj_x->Write(hname+"_"+suffix);

   }  
   TGraphErrors * mygph = new TGraphErrors (nbbins,xaxis,resol,errxaxis,errresol);
   mygph->SetName("gph_"+hname);
   mygph->Write();
}
