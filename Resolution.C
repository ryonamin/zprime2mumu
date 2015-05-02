#define Resolution_cxx
#include "Resolution.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TF1.h>
#include <TGraphErrors.h>
using namespace std;

void Resolution::Loop()
{
   //   In a ROOT session, you can do:
   //      Root > .L Resolution.C
   //      Root > Resolution t
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

   TH1D * h1d_pt_barrel = new TH1D("h1d_pt_barrel","",100,-1,1);
   TH2D * h2d_pt_barrel = new TH2D("h2d_pt_barrel","",100,-1,1,10,0,2600);  
   TH1D * h1d_pt_endcaps = new TH1D("h1d_pt_endcaps","",100,-1,1);
   TH2D * h2d_pt_endcaps = new TH2D("h2d_pt_endcaps","",100,-1,1,10,0,2600);  

   TH1D * h1d_mass = new TH1D("h1d_mass","",100,-1,1);
   TH2D * h2d_mass = new TH2D("h2d_mass","",100,-1,1,10,0,5200);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry</*nentries*/3000;jentry++) {

      if(jentry%1000 == 0) cout << jentry << " over " << nentries << endl;
      TLorentzVector recomu_1,recomu_2, genmu_1,genmu_2;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);  nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      for(int i = 0; i < (*mu_gt_pt).size();i++){
	 //The first line just ensures that the global muon exists.
	 //In Aidan's code all muons are stored in a same vector but not all of them are global. 
	 if((*mu_gt_pt)[i] <10) continue;
	 double dR = 0.3; 
	 double mupt = 0;
	 bool matchgen = false;
	 TLorentzVector recomu,genmu; 
	 recomu.SetPxPyPzE((*mu_gt_px)[i],(*mu_gt_py)[i],(*mu_gt_pz)[i],(*mu_gt_p)[i]);
	 if(!PassHighPtSel(i)) continue;

	 for(int j = 0; j < mc_n; j++){
	    //Look only at muons 
	    if(abs((*mc_pdgId)[j])!=13) continue; 
	    //If ptgen<30 we don't care
	    if((*mc_pt)[j] < 30) continue;
	    TLorentzVector genmu_loop;
	    genmu_loop.SetPxPyPzE((*mc_px)[j],(*mc_py)[j],(*mc_pz)[j],(*mc_energy)[j]);
	    //Gen-Reco matching: do a delta R matching and save the gen particle idx
	    //Also require that the charge is correctly measured
	    if(genmu_loop.DeltaR(recomu) < dR && genmu_loop.Pt() > mupt && (*mc_charge)[j] == (*mu_gt_charge)[i]) {
	       matchgen = true;
	       mupt = genmu_loop.Pt();
	       genmu = genmu_loop;
	    }
	 }

	 //if the reco matches the gen, fill histos
	 if(matchgen) {
	    if(recomu_1.E() != 0 && recomu_2.E() == 0) { recomu_2 = recomu; genmu_2 = genmu; }
	    if(recomu_1.E() == 0){ recomu_1 = recomu  ; genmu_1 =genmu; }
	    //Split barrel and endcaps
	    if(fabs( (*mu_gt_eta)[i]) < 1.2) {
	       h2d_pt_barrel->Fill(1-recomu.Pt()/genmu.Pt(), genmu.Pt());
	       h1d_pt_barrel->Fill(1-recomu.Pt()/genmu.Pt());
	    }
	    else {
	       h2d_pt_endcaps->Fill(1-recomu.Pt()/genmu.Pt(), genmu.Pt());
	       h1d_pt_endcaps->Fill(1-recomu.Pt()/genmu.Pt());
	    }
	 }	  
      }

      //Compute the invariant mass     
      double genmass = (genmu_1 + genmu_2).M();
      double recomass = (recomu_1 + recomu_2).M();
      //If two reco muons matching gen muons are found, fill histos
      if(recomu_1.E() != 0 && recomu_2.E() != 0){

	 h1d_mass->Fill(1-recomass/genmass);
	 h2d_mass->Fill(1-recomass/genmass, genmass);
      }


   }
   myfile->Write();  
   Resolvs(h2d_mass);
   Resolvs(h2d_pt_barrel);
   Resolvs(h2d_pt_endcaps);
}




bool Resolution::PassHighPtSel(int n){
   //This is an old selection that needs to be updated... but better than nothing.
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
	 && fabs((*mu_gt_eta)[n]) < muon_etaMax
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


void Resolution::Resolvs(TH2D * myh2d){
   TString hname=myh2d->GetName(); 
   const int nbbins = myh2d->GetNbinsY();
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
