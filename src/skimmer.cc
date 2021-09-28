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
                     Filter(Jet_etas_multiplied, {"Jet_etas_multiplied"}, "opposite_jet_eta").
                     Filter(dijet_Mjj, {"dijet_Mjj"}, "dijet_Mjj<400").
                     Filter(dijet_abs_dEta, {"dijet_abs_dEta"}, "dijet_abs_eta<2.4").
                     Filter(zpt, {"Z_pt"}, "boson_pt>55");
  return tmp;
}

RDF::RNode applyCutsStrings(RDF::RNode df, RDataFrame::ColumnNames_t branches) {
  if (branches.size() == 0 ) return df;
  std::string cut = "";
  for (auto &br: branches) {
    if (&br != &branches.back())
      cut = cut + br + " == 1 || ";
    else
      cut = cut + br + " == 1";
  }
  return df.Filter(cut);
}

int main(int argc, char **argv){

  Options options(argc, argv);
  YAML::Node const config = options.GetConfig();
  std::string tree = Options::NodeAs<std::string>(config, {"tree_name"});
  std::vector<std::string> files = Options::GetStrings(config, {"samples"});
  
  RDataFrame::ColumnNames_t varibles; 
  for (auto &branch: Options::GetStrings(config,{"branch", "common"}))
    varibles.push_back(branch);
  for (auto &aName : files)
  {
    std::vector<std::string> filepaths;
    FileInPath::GetFilenames(FileInPath::Resolve(Options::NodeAs<std::string>(config, {"attributes", aName, "path"})), filepaths);
    if (boost::contains(Options::NodeAs<std::string>(config, {"attributes", aName, "type"}), "MC"))
      for (auto &br : Options::GetStrings(config,{"branch", "MC"}))
        varibles.push_back(br);
    if (boost::contains(Options::NodeAs<std::string>(config, {"attributes", aName, "type"}), "photon")) {
      for (auto &br : Options::GetStrings(config,{"branch", "photon"}))
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
    double xsec = Options::NodeAs<double>(config, {"attributes", aName, "xsec"});
    double factor = lumi * xsec / num;
    int log_value = (int)log10(factor);
    int precision = 3;
    if (log_value < 0) precision = abs(log_value) + 6;
    else precision = 6;
    std::string reweight_factor = utils::to_string_with_precision(factor, precision); 
    std::cout<< reweight_factor <<std::endl;
    auto df_tmp = df.Define("xsec_reweight", reweight_factor); 
    auto df_tmp2 = applyCutsCommon(df_tmp);
    auto triggers = Options::GetStrings(config,{"HLT", Options::NodeAs<std::string>(config, {"attributes", aName, "type"})});
    auto df_tmp3 = applyCutsStrings(df_tmp2, triggers);
    varibles.push_back("xsec_reweight");
    std::string dir_name = "";
    if (boost::contains(Options::NodeAs<std::string>(config, {"attributes", aName, "process"}), "HTbinDY")) dir_name = "outputs/HTbinDY/";
    else if (boost::contains(Options::NodeAs<std::string>(config, {"attributes", aName, "process"}), "GJet")) dir_name = "outputs/GJet/";
    else if (boost::contains(Options::NodeAs<std::string>(config, {"attributes", aName, "process"}), "inclusiveDY")) dir_name = "outputs/inclusiveDY/";
    df_tmp3.Snapshot("Events", dir_name + aName+".root", varibles);
  }
  return 0;
}

