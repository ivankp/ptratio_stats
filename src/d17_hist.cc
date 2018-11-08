#include <iostream>
#include <vector>
#include <stdlib.h>
#include <utility>
#include <chrono>

#include <boost/optional.hpp>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include "Higgs2diphoton.hh"

using namespace std;
bool sort_pt (TLorentzVector i, TLorentzVector j){return (i.Pt()>j.Pt());}
Higgs2diphoton Hdecay;

int main(){
  //Loading in data
  TChain chain_d("mini");

  chain_d.Add("$TMPDIR/data1516.root");
  chain_d.Add("$TMPDIR/data17.root");

  //Setting chains
  constexpr size_t max = 15;
  Int_t photon_n, jet_n;
  Float_t photon_m[max], photon_eta[max], photon_phi[max], photon_pt[max];
  Float_t jet_m[max], jet_eta[max], jet_phi[max], jet_pt[max];

  chain_d.SetBranchAddress("photon_n", &photon_n);
  chain_d.SetBranchAddress("photon_m", photon_m);
  chain_d.SetBranchAddress("photon_eta", photon_eta);
  chain_d.SetBranchAddress("photon_phi", photon_phi);
  chain_d.SetBranchAddress("photon_pt", photon_pt);

  chain_d.SetBranchAddress("jet_n", &jet_n);
  chain_d.SetBranchAddress("jet_m", jet_m);
  chain_d.SetBranchAddress("jet_eta", jet_eta);
  chain_d.SetBranchAddress("jet_phi", jet_phi);
  chain_d.SetBranchAddress("jet_pt", jet_pt);

  //===============================================
  //ALL THE TREES
  Float_t ptratio;
  //===============================================
  vector<TString>num;
  num.emplace_back("300");
  num.emplace_back("350");
  num.emplace_back("400");
  num.emplace_back("400p");

  for (int i=0; i<4; ++i){
    TString outname = "ptratio_data_" + num[i] + ".root";
    TFile* out = TFile::Open(outname, "RECREATE");
    cout<<"Writing "<<outname<<"....."<<endl;

    //===============================================
    //Ouputs
    int nbins = 38.;

    //Signal
    TH1D* h_ptratio_s_data = new TH1D("Data_Signal__ptratio_"+num[i], "p_{T}^{#gamma_{1}} / p_{T}^{#gamma_{2}}", nbins, 1, 20);
    //Sideband
    TH1D* h_ptratio_b_right = new TH1D("Data_Right__ptratio_"+num[i], "p_{T}^{#gamma_{1}} /p_{T}^{#gamma_{2}}", nbins, 1, 20);
    TH1D* h_ptratio_b_left = new TH1D("Data_Left__ptratio_"+num[i], "p_{T}^{#gamma_{1}} /p_{T}^{#gamma_{2}}", nbins, 1, 20);

    //Data
    //Go through each event
    const Long64_t entries_d = chain_d.GetEntries();
    for (long k=0; k<entries_d; ++k){
      chain_d.GetEntry(k);
      //Assign Particles
      TLorentzVector yy, y1, y2, jet, photons, jets;
      vector<TLorentzVector> t_y, t_jet;

      for (int j=0; j<photon_n; ++j){
        photons.SetPtEtaPhiM(photon_pt[j], photon_eta[j], photon_phi[j], photon_m[j]);
        t_y.emplace_back(photons);
      }
      if (jet_n!=0){
        for (int j=0; j<jet_n; ++j){
          jets.SetPtEtaPhiM(jet_pt[j], jet_eta[j], jet_phi[j], jet_m[j]);
          t_jet.emplace_back(jets);
        }
      }
      else{
        jets.SetPtEtaPhiM(0,0,0,0);
        t_jet.emplace_back(jets);
      }

      //Sort by pT
      sort(t_jet.begin(), t_jet.end(), sort_pt);
      sort(t_y.begin(), t_y.end(), sort_pt);

      jet = t_jet[0];
      y1 = t_y[0];
      y2 = t_y[1];
      yy = y1+y2;

      //Cuts
      //Rapidity Cut
      bool select = (abs(y2.Rapidity())<2.4);
      select &= (abs(y1.Rapidity())<2.4);
      //PseudoRapidity Cut
      select &= (abs(y1.PseudoRapidity())<2.37);
      select &= !(1.37<abs(y1.PseudoRapidity()) && abs(y1.PseudoRapidity())<1.52);
      select &= (abs(y2.PseudoRapidity())<2.37);
      select &= !(1.37<abs(y2.PseudoRapidity()) && abs(y2.PseudoRapidity())<1.52);
      //Pt Cut
      select &= (y1.Pt()>.35*yy.M());
      select &= (y2.Pt()>.25*yy.M());
      //Delta R
      select &= (y1.DeltaR(jet)>.4);
      select &= (y2.DeltaR(jet)>.4);
      //Jet Cuts
      for (int i=0; i<jet_n; ++i){
        select &= (t_jet[i].Pt()>30);
        select &= (t_jet[i].Rapidity()<4.4);
      }

      if (i==0){select &= (yy.Pt()>250. && yy.Pt()<300.);}
      if (i==1){select &= (yy.Pt()>300. && yy.Pt()<350);}
      if (i==2){select &= (yy.Pt()>350. && yy.Pt()<400);}
      if (i==3){select &= (yy.Pt()>400.);}

      if (select){
        ptratio = abs(y1.Pt()/y2.Pt());

        bool sig = ((yy.M()>122.)&&(yy.M()<128.));
        bool left = ((yy.M()>105.)&&(yy.M()<121.));
        bool right = ((yy.M()>129.)&&(yy.M()<160.));

        if (sig){h_ptratio_s_data->Fill(ptratio);}
        if (left){h_ptratio_b_left->Fill(ptratio);}
        if (right){h_ptratio_b_right->Fill(ptratio);}
      }
    }

    cout<<outname<<" written"<<endl;

    out->Write();
    out->Close();
  }
}
