#include <iostream>
#include <string>
#include <vector>
#include <initializer_list>

#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TRatioPlot.h>
#include <yaml-cpp/yaml.h>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <FileInPath.h>
#include <Options.h>
#include <Utils.h>
using namespace ROOT;
namespace po = boost::program_options;


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
                     Filter(delta_phi_j_met, {"delta_phi_j_met"}, "delta_phi_j_met").
                     //Filter(delta_phi_ZMet_bst, {"delta_phi_ZMet_bst"}).
                     Filter(Jet_etas_multiplied, {"Jet_etas_multiplied"}, "opposite_jet_eta").
                     Filter(dijet_Mjj, {"dijet_Mjj"}, "dijet_Mjj<400").
                     Filter(dijet_abs_dEta, {"dijet_abs_dEta"}, "dijet_abs_eta<2.4").
                     Filter(zpt, {"Z_pt"}, "boson_pt>55");
  return tmp;
}
int main(int argc, char **argv){

  Options options(argc, argv);
  YAML::Node const config = options.GetConfig();
  std::string tree = Options::NodeAs<std::string>(config, {"tree_name"});
 // std::string dilepton_filename = Options::NodeAs<std::string>(config, {"dilepton_files"});
  std::vector<std::string> files = Options::GetStrings(config, {"samples"});
  //std::string photon_filename = Options::NodeAs<std::string>(config, {"photon_files"});
  bool isMC = Options::NodeAs<bool>(config, {"isMC"});
  
  RDataFrame::ColumnNames_t varibles_test = {"nJet"};
  RDataFrame::ColumnNames_t varibles = {"nJet","Jet_pt_nom", 
      "lead_jet_pt", "lead_jet_phi", "trail_jet_pt", "trail_jet_eta",
      "trail_jet_phi", "lep_category", "ngood_jets", "ngood_bjets",
      "nhad_taus", "met_pt", "met_phi", "delta_R_ll", "delta_phi_j_met",
      "Jet_etas_multiplied", "dijet_Mjj", "dijet_abs_dEta", 
      "Z_pt", "Z_eta", "Z_phi", "Z_mass",
      "Pileup_nPU",
      //"ngood_leptons", "nextra_leptons",
      //"deltaPhiClosestJetMet", "deltaPhiFarthestJetMet",
      //"delta_phi_ZMet_bst",
      //"nloose_photons",
       "delta_phi_ZMet",};
  for (auto &aName : files)
  {
    std::cout<< aName << std::endl;
    std::vector<std::string> filepaths;
    FileInPath::GetFilenames(FileInPath::Resolve(Options::NodeAs<std::string>(config, {"file_paths", aName})), filepaths);
  //std::string photonPath = FileInPath::Resolve(photon_filename);
    if (isMC) {
      //varibles.insert(0, RDataFrame::ColumnNames_t{"weight"});
      for (auto &br : RDataFrame::ColumnNames_t{"weight", "puWeight", "w_muon_SF", "w_electron_SF",
          "Jet_qgl"})
        varibles.push_back(br);
    }
    if (boost::contains(aName, "GJet")) {
      //varibles.insert(0, RDataFrame::ColumnNames_t{"weight"});
      for (auto &br : RDataFrame::ColumnNames_t{
        "nPhoton", "nMuon", "nElectron",
        "ngood_leptons", "nextra_leptons",
        //"deltaPhiClosestJetMet", "deltaPhiFarthestJetMet",
        //"delta_phi_ZMet_bst",
        //"nloose_photons"
        })
        varibles.push_back(br);
    }
  // enabling Multi-Thread
    ROOT::EnableImplicitMT();
  // data-frame initializing
    RDataFrame df(tree, filepaths, varibles);
    RDataFrame df_run("Runs", filepaths, {"genEventSumw2"});
    std::cout<< df_run.Sum("genEventSumw2").GetValue() <<std::endl;
    long int num = df_run.Sum("genEventSumw2").GetValue();
  // apply common selections
    double lumi = Options::NodeAs<double>(config, {"lumi"});;
    double xsec = Options::NodeAs<double>(config, {"xsec", aName});
    double factor = lumi * xsec / num;
    int log_value = (int)log10(factor);
    int precision = 3;
    if (log_value < 0) precision = abs(log_value) + 6;
    else precision = 6;
    std::string reweight_factor = utils::to_string_with_precision(factor, precision); 
    std::cout<< reweight_factor <<std::endl;
    auto df_tmp = df.Define("xsec_reweight", reweight_factor); 
    auto df_tmp2 = applyCutsCommon(df_tmp);
    varibles.push_back("xsec_reweight");
    std::string dir_name = "";
    if (boost::contains(aName, "DY")) dir_name = "outputs/DY/";
    else if (boost::contains(aName, "GJet")) dir_name = "outputs/GJet/";
    df_tmp2.Snapshot("Events", dir_name + aName+".root", varibles);
  }
  return 0;
}

