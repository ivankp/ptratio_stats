#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include "cat.hh"
#include "timed_counter.hh"

using std::vector;
using std::cout;
using std::endl;

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

const double inf = std::numeric_limits<double>::infinity();

template <typename T, typename F>
void loop1(const T& xs, F&& f) {
  auto it = next(begin(xs));
  const auto _end = end(xs);
  for (; it<_end; ++it) f(it);
}

bool photon_eta_cut(double abs_eta) {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

template <typename T>
void pt_sort(T& xs) {
  sort(begin(xs),end(xs),[](const auto& a, const auto& b) {
    return a.Pt() > b.Pt();
  });
}

int main(int argc, char* argv[]) {
  // Loading in data
  TChain tree("mini");

  tree.Add("data1516.root");
  tree.Add("data17.root");

  // Setting chains
  Int_t photon_n=8, jet_n=32;
  vector<Float_t>
    photon_pt(photon_n), photon_eta(photon_n),
    photon_phi(photon_n), photon_m(photon_n),
    jet_pt(jet_n), jet_eta(jet_n),
    jet_phi(jet_n), jet_m(jet_n);

  tree.SetBranchAddress("photon_n"  ,&photon_n);
  tree.SetBranchAddress("photon_pt" , photon_pt .data());
  tree.SetBranchAddress("photon_eta", photon_eta.data());
  tree.SetBranchAddress("photon_phi", photon_phi.data());
  tree.SetBranchAddress("photon_m"  , photon_m  .data());
  tree.SetBranchAddress("jet_n"  ,&jet_n);
  tree.SetBranchAddress("jet_pt" , jet_pt .data());
  tree.SetBranchAddress("jet_eta", jet_eta.data());
  tree.SetBranchAddress("jet_phi", jet_phi.data());
  tree.SetBranchAddress("jet_m"  , jet_m  .data());

  const vector<double>
    bins_pt_yy { 250, 300, 350, 400, inf },
    bins_m_yy  { 105, 121, 129, 160 };

  TFile fout("ptrat.root","recreate");

  vector<TH1D*> h_ptrat;
  h_ptrat.reserve( bins_pt_yy.size()*bins_m_yy.size() );
  loop1(bins_pt_yy, [&](auto pt_yy){
    loop1(bins_m_yy, [&](auto m_yy){
      h_ptrat.push_back(new TH1D(cat(
        "ptrat"
        " pt_yy:",*prev(pt_yy),'-',*pt_yy,
        " m_yy:",*prev(m_yy),'-',*m_yy
      ).c_str(), "p_{T}^{#gamma_{1}} / p_{T}^{#gamma_{2}}", 38, 1, 20));
    });
  });

  auto fill = [&](double pt_yy, double m_yy, double x) -> void {
    auto index = [](const auto& edges, auto x) -> int {
      const auto _begin = begin(edges);
      const auto _end = end(edges);
      const auto it = std::upper_bound(_begin, _end, x);
      if (it==_begin || it==_end) return -1;
      return std::distance(_begin,it)-1;
    };
    const auto i_pt_yy = index(bins_pt_yy,pt_yy);
    if (i_pt_yy<0) return;
    const auto i_m_yy = index(bins_m_yy,m_yy);
    if (i_m_yy<0) return;

    static int n = bins_m_yy.size()-1;
    h_ptrat[i_pt_yy*n+i_m_yy]->Fill(x);
  };

  for (timed_counter<Long64_t> ent(tree.GetEntries()); !!ent; ++ent) {
    tree.GetEntry(ent);

    if (photon_n<2) continue;

    // Assign Particles
    vector<TLorentzVector> photons(photon_n)/*, jets(jet_n)*/;

    for (int i=0; i<photon_n; ++i) {
      photons[i].SetPtEtaPhiM(
        photon_pt[i], photon_eta[i], photon_phi[i], photon_m[i]);
    }
    /*
    for (int i=0; i<jet_n; ++i) {
      jets[i].SetPtEtaPhiM(
        jet_pt[i], jet_eta[i], jet_phi[i], jet_m[i]);
    }
    */

    // Sort by pT
    pt_sort(photons);
    // pt_sort(jets);

    const TLorentzVector& y1 = photons[0];
    const TLorentzVector& y2 = photons[1];

    const auto y1_eta = std::abs(y1.Eta());
    const auto y2_eta = std::abs(y2.Eta());

    const auto y1_pt = y1.Pt();
    const auto y2_pt = y2.Pt();

    if ( // Apply photon cuts
      (y1_pt < 0.35*125.) or
      (y2_pt < 0.25*125.) or
      photon_eta_cut(y1_eta) or
      photon_eta_cut(y2_eta)
    ) continue;

    const TLorentzVector yy = y1 + y2;
    fill( yy.Pt(), yy.M(), std::abs(y1.Pt()/y2.Pt()) );
  }

  fout.Write();
  cout << "Wrote " << fout.GetName() << endl;
}
