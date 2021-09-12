#include <iostream>
#include <string>
#include <vector>
#include <initializer_list>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <yaml-cpp/yaml.h>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <FileInPath.h>
#include <Options.h>
#include <Utils.h>
using namespace ROOT;
namespace po = boost::program_options;

std::map<double, std::pair<double, double>> TH1toMap(TH1F * h_weight) {
  std::map<double, std::pair<double, double> > myMap; //label = bin low edge; pair = (bin content; bin error)
  int nBinsX = h_weight->GetNbinsX();
  for(int bin_X = 1; bin_X <= nBinsX; bin_X++){
    myMap[h_weight->GetXaxis()->GetBinLowEdge(bin_X)].first = h_weight->GetBinContent(bin_X); 
    myMap[h_weight->GetXaxis()->GetBinLowEdge(bin_X)].second = h_weight->GetBinError(bin_X); 
    
  }
  return myMap;
}

TH1F * GetWeightHisto(TH1F * h_ll, TH1F * h_gamma,  TString obs) {

  h_ll->Scale(1/h_ll->Integral());
  h_gamma->Scale(1/h_gamma->Integral());
  auto h_weight = (TH1F*) h_ll->Clone();
  h_weight->Divide(h_gamma);
  h_weight->SetName(obs + "_weight");
  h_weight->SetTitle(obs + "_weight");
  if (obs.EqualTo("nvtx")){
  auto h_nvtx_weight_smooth = utils::SmoothenHisto(h_weight);
  return h_nvtx_weight_smooth;
  }
  else
    return h_weight;
}
RDF::RNode applyCutsCommon(RDF::RNode df){
  std::function<bool(int)> nGoodJet = [](int ngood_jets) { return ngood_jets >= 2; };
  std::function<bool(int)> nGoodbJet = [](int ngood_bjets) { return ngood_bjets == 0; };
  std::function<bool(int)> nTaus = [](int ntau) { return ntau == 0; };
  std::function<bool(float)> zpt = [](float Z_pt) { return Z_pt > 55.; };
  std::function<bool(float)> met = [](float met) { return met < 60.; };
  std::function<bool(float)> delta_phi_j_met = [](float delta_phi_j_met) { return abs(delta_phi_j_met) > 0.5; };
  std::function<bool(float)> delta_phi_ZMet_bst = [](float delta_phi_ZMet_bst) { return abs(delta_phi_ZMet_bst) > 1.0; };
  std::function<bool(float)> Jet_etas_multiplied = [](float Jet_etas_multiplied) { return Jet_etas_multiplied < 0.; };
  std::function<bool(float)> dijet_Mjj = [](float dijet_Mjj) { return dijet_Mjj > 400.; };
  std::function<bool(float)> dijet_abs_dEta = [](float dijet_abs_dEta) { return dijet_abs_dEta > 2.4; };
  auto tmp = df.
                     Filter(nGoodJet, {"ngood_jets"}, "2-jets").
                     Filter(nGoodbJet, {"ngood_bjets"}, "bveto").
                     Filter(nTaus, {"nhad_taus"}, "tau-veto").
                     Filter(met, {"met_pt"}, "MET<60").
                     Filter(delta_phi_j_met, {"delta_phi_j_met"}, "delta_phi_j_met").
                    // Filter(delta_phi_ZMet_bst, {"delta_phi_ZMet_bst"}).
                     //Filter(Jet_etas_multiplied, {"Jet_etas_multiplied"}).
                     //Filter(dijet_Mjj, {"dijet_Mjj"}).
                     //Filter(dijet_abs_dEta, {"dijet_abs_dEta"}).
                     Filter(zpt, {"Z_pt"}, "Boson_pt");
  return tmp;
}

RDF::RNode applyCutsLL(RDF::RNode df){
  std::function<bool(int)> lep_cat = [](int lep_category) { return lep_category == 1 || lep_category == 3; };
  std::function<bool(float)> deltaR_ll = [](float deltaR_ll) { return deltaR_ll < 2.5; };
  auto tmp = df.Filter(lep_cat, {"lep_category"}).
                Filter(deltaR_ll, {"deltaR_ll"})
;
  return tmp;
}
RDF::RNode applyCutsPhoton(RDF::RNode df){
  std::function<bool(float)> boson_eta = [](float Z_eta) { return abs(Z_eta) < 2.4; };
  auto tmp = df.
             Filter(boson_eta, {"Z_eta"});
  return tmp;
}
  float applyWeight   (TH1F * h_weight, float boson_eta) {
    //auto h_weight = (TH1F*) gROOT->Get("eta_weight");
    auto const bin = h_weight->FindFixBin(boson_eta);
    return h_weight->GetBinContent(bin);
    //std::cout << h_weight->GetBinContent(bin) <<std::endl;
  }
int main(int argc, char **argv){

  Options options(argc, argv);
  YAML::Node const config = options.GetConfig();
  std::string tree = Options::NodeAs<std::string>(config, {"tree_name"});
  std::string dilepton_filename = Options::NodeAs<std::string>(config, {"dilepton_files"});
  std::string photon_filename = Options::NodeAs<std::string>(config, {"photon_files"});
  bool isMC = Options::NodeAs<bool>(config, {"isMC"});
  std::string dileptonPath = FileInPath::Resolve(dilepton_filename);
  std::string photonPath = FileInPath::Resolve(photon_filename);
  std::vector<std::string> dilepton_files;
  std::vector<std::string> photon_files;
  FileInPath::GetFilenames(dileptonPath, dilepton_files);
  FileInPath::GetFilenames(photonPath, photon_files);
  RDataFrame::ColumnNames_t varibles = {"nJet","Muon_pt" "Jet_pt_nom", 
      "lead_jet_pt", "lead_jet_phi", "trail_jet_pt", "trail_jet_eta",
      "trail_jet_phi", "lep_category", "ngood_jets", "ngood_bjets",
      "nhad_taus", "met_pt", "met_phi", "delta_R_ll", "delta_phi_j_met",
      "Jet_etas_multiplied", "dijet_Mjj", "dijet_abs_dEta", 
      "Z_pt", "Z_eta", "Z_phi", "Z_mass",
      "nPhoton", "nMuon", "nElectron", "Pileup_nPU",
      "ngood_leptons", "nextra_leptons",
      "deltaPhiClosestJetMet", "deltaPhiFarthestJetMet",
      "delta_phi_ZMet_bst", "delta_phi_ZMet"};
  RDataFrame::ColumnNames_t storeBranches = {"boson_pt", "boson_eta", "Pileup_nPU", "met_pt"};
  if (isMC) {
    //varibles.insert(0, RDataFrame::ColumnNames_t{"weight"});
    for (auto &br : RDataFrame::ColumnNames_t{"weight", "puWeight", "w_muon_SF", "w_electron_SF",
        "PrefireWeight", "nvtxWeight", "TriggerSFWeight", "btagEventWeight", "Jet_partonFlavour",
        "Jet_qgl"})
      varibles.push_back(br);
  }
  RDataFrame df_dilepton(tree, dilepton_files, varibles);
  RDataFrame df_photon(tree, photon_files, varibles);
  auto df_ll_filtered_tmp = applyCutsCommon(df_dilepton);
  auto df_gamma_filtered_tmp = applyCutsCommon(df_photon);
  auto df_ll_filtered = applyCutsCommon(df_ll_filtered_tmp);
  auto df_gamma_filtered = applyCutsCommon(df_gamma_filtered_tmp);
  auto df_ll_aug = df_ll_filtered.
                  Define("boson_pt", "Z_pt").
                  Define("boson_eta", "abs(Z_eta)").
                  Define("nvtx", "Pileup_nPU").
                  Define("ptmiss", "met_pt").
                  Define("corr", "weight");

  auto df_gamma_aug = df_gamma_filtered.
                  Define("boson_pt", "Z_pt").
                  Define("boson_eta", "abs(Z_eta)").
                  Define("nvtx", "Pileup_nPU").
                  Define("ptmiss", "met_pt").
                  Define("corr", "weight");
  auto h_nvtx_ll = (TH1F*)df_ll_aug.Histo1D({"nvtx_ll","nvtx_ll",100,0,100},"nvtx", "corr").GetValue().Clone();
  auto h_nvtx_gamma = (TH1F*)df_gamma_aug.Histo1D({"nvtx_gamma","nvtx_gamma",100,0,100},"nvtx", "corr").GetValue().Clone();
  auto h_eta_ll = (TH1F*)df_ll_aug.Histo1D({"eta_ll","eta_ll",24,0,2.4},"boson_eta", "corr").GetValue().Clone();
  auto h_eta_gamma = (TH1F*)df_gamma_aug.Histo1D({"eta_gamma","eta_gamma",24,0,2.4},"boson_eta", "corr").GetValue().Clone();
  //auto f1 = new TFile("weights.root", "RECREATE");
  auto h_nvtx_weight = GetWeightHisto(h_nvtx_ll, h_nvtx_gamma, "nvtx");
  auto h_eta_weight = GetWeightHisto(h_eta_ll, h_eta_gamma, "eta");
  //h_nvtx_weight->Write();
  //h_eta_weight->Write();
  //f1->Close();
  //std::cout<<applyEtaWeight(h_eta_weight, 1.1) << std::endl;
  auto applyEtaWeight = [] (float boson_eta) {
    TH1F * h_weight = (TH1F *)gROOT->Get("eta_weight");
    auto const bin = h_weight->FindFixBin(boson_eta);
    return h_weight->GetBinContent(bin);
  };
  auto df_gamma_eta_weighted_tmp = df_gamma_aug.Define("eta_weight", applyEtaWeight, {"boson_eta"});
  auto df_gamma_eta_weighted = df_gamma_eta_weighted_tmp.Define("corr_eta", "corr * eta_weight");
  const float  binning[] = {55,75,100,150,200,250,1500};
  RDF::TH1DModel pt_model("pt_ll","pt_ll",6, binning);
  RDF::TH1DModel pt_model_gamma("pt_gamma","pt_gamma",6, binning);
  auto h_pt_ll = (TH1F*)df_ll_aug.Histo1D(pt_model,"boson_pt", "corr").GetValue().Clone();
  auto h_pt_gamma = (TH1F*)df_gamma_eta_weighted.Histo1D(pt_model_gamma,"boson_pt", "corr_eta").GetValue().Clone();
  auto h_pt_weight = GetWeightHisto(h_pt_ll, h_pt_gamma, "pt");
  
  auto applyPtWeight = [] (float boson_pt) {
    TH1F * h_weight = (TH1F *)gROOT->Get("pt_weight");
    auto const bin = h_weight->FindFixBin(boson_pt);
    return h_weight->GetBinContent(bin);
  };
  auto df_gamma_eta_pt_weighted = df_gamma_eta_weighted.Define("pt_weight", applyPtWeight, {"boson_pt"});
  std::cout << "All stats:" << std::endl;
  auto allCutsReport = df_gamma_filtered.Report();
  std::cout << "Name\tAll\tPass\tEfficiency" << std::endl;
   for (auto &&cutInfo : allCutsReport) {
      std::cout << cutInfo.GetName() << "\t" << cutInfo.GetAll() << "\t" << cutInfo.GetPass() << "\t"
                << cutInfo.GetEff() << " %" << std::endl;
   }
 
  df_gamma_eta_pt_weighted.Snapshot("Events","llll.root");
  auto f1 = new TFile("weights.root", "RECREATE");
  h_nvtx_weight->Write();
  h_eta_weight->Write();
  h_pt_weight->Write();
  h_pt_gamma->Write();
  h_pt_ll->Write();
  //df_gamma_aug.Foreach(applyEtaWeight, {"boson_eta"});
  f1->Close();
  //std::cout << wrapper.getEventNum(lep_cat) << std::endl;
  return 0;
}

