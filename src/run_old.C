// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 11 Fev 2012
//
// run_analysis()
// compiles the analysis routines and runs them to do the analysis

void run_analysis()
{
  gROOT -> Reset();
  gROOT->SetStyle("Plain");

 //------------ FWLite libraries ------
  gSystem->Load("libFWCoreFWLite.so");
  gSystem->Load("libKKousourQCDAnalysis.so");
  gSystem->Load("../RooUnfold-1.1.1/libRooUnfold.so");
  AutoLibraryLoader::enable();

//verbose level
  bool detail = true;
  bool show_results = true;
  bool show_steps = true;

//enable testing
  bool test = false; //if set to true will get only some entries from histograms and save outputs with extension _test

//steps to be done
  bool get_trigger_efficiencies = false; //not used anymore
  bool compute_trigger_efficiency = false; //not used anymore
  bool compute_trigger_turn_on = false;
  bool compute_trigger_correction_v1 = false; //not used anymore
  bool check_trigger_correction_v1 = false; //not used anymore
  bool compute_trigger_correction_v2 = false; //not used anymore
  bool check_trigger_correction_v2 = false; //not used anymore
  bool apply_trigger_correction = false; //not used anymore
  bool apply_delta_phi_factor = false; // not needed or recomended
  bool plot_final_trigger_efficiency = false;
  bool counting_combination_statistics = false; //need to have root output files set
  bool get_vertex_distribution_v0 = false; //need to run
  bool compute_vertex_weights_v1 = false; //need to run
  bool get_vertex_distribution_v1 = false; //need to run
  bool compute_vertex_weights_v2 = false; //need to run
  bool get_vertex_distribution_v2 = false; //need to run
  bool compute_vertex_weights_v3 = false; //need to run
  bool get_vertex_distribution_v3 = false; //need to run
  bool compute_vertex_weights_v4 = false; //need to run
  bool get_vertex_distribution_v4 = false; //need to run
  bool compute_vertex_weights_v5 = false; //need to run
  bool plot_status_vertex_reweight = false; //need to run
  bool get_pileup_normalization = false; //not used anymore - i think
  bool read_mc_ntuples_gen = false; //need to run
  bool read_mc_ntuples_det = false; //need to run
  bool read_data_ntuples = false; //need to run
  bool normalize_mc = false; //need to run
  bool plot_control_dist = false; //need to change and run
  bool compute_corrections = false; //need to run
  bool apply_corrections = false; //need to run
  bool merge_data = false; //need to change and run
  bool compute_model_uncertainty = false; //need to change and run
  bool compute_jes_uncertainty = false; //need to change and run
  bool estimate_combination_systematic = false; //need to change and run
  bool compute_total_uncertainty = false; //need to change and run
  bool apply_uncertainty = false; //need to run
  bool convert_predictions = false; //not needed anymore
  bool create_ratio = true; //need to run
  bool do_final_plots = true; //need to run
  bool create_unfolding_response = false; //need to change and run
  bool check_response_matrix = false;
  bool unfold = false;

//directory configuration
  string output_dir = "../output/";

// histogram directories
  string histograms_dir = output_dir + "histograms/";
  string trigger_syst_dir = histograms_dir + "trigger_systematic/";
  string trigger_dir = histograms_dir + "trigger_eff/";
  string trigger_corr_dir = histograms_dir + "trigger_correction/";
  string eff_dir = histograms_dir + "efficiency/";
  string vertex_dir = histograms_dir + "vertex/";
  string vertex_weights_dir = histograms_dir + "vertex_weights/";
  string pileup_norm_dir = histograms_dir + "pileup_norm/";
  string mc_det_dir = histograms_dir + "xsec_mc_det/";
  string mc_gen_dir = histograms_dir + "xsec_mc_gen/";
  string mc_norm_dir = histograms_dir + "normalized_mc/";
  string raw_data_dir = histograms_dir + "raw_xsec_data/";
  string corrections_dir = histograms_dir + "corrections/";
  string corrected_dir = histograms_dir + "corrected_data/";
  string model_unc_dir = histograms_dir + "model_uncertainty/";
  string jes_unc_dir = histograms_dir + "jes_uncertainty/";
  string total_unc_dir = histograms_dir + "total_uncertainty/";
  string mc_pred_dir = histograms_dir + "mc_prediction/";
  string ratio_dir = histograms_dir + "ratio_to_data/";
  string response_dir = histograms_dir + "unfolding_response/";
  string unfolded_dir = histograms_dir + "unfolded/";

//plot directories
  string trigger_syst_plots = output_dir + "trigger_systematic/";
  string trigger_plots = output_dir + "trigger_eff/";
  string trigger_corr_plots = output_dir + "trigger_correction/";
  string vertex_weights_plots = output_dir + "vertex_weights/";
  string mc_norm_plots = output_dir + "normalized_mc/";
  string control_plots = output_dir + "control_dist/";
  string corrections_plots = output_dir + "corrections/";
  string corrected_plots = output_dir + "corrected_data/";
  string merged_data_plots = output_dir + "merged_data/";
  string model_unc_plots = output_dir + "model_uncertainty/";
  string jes_unc_plots = output_dir + "jes_uncertainty/";
  string total_unc_plots = output_dir + "total_uncertainty/";
  string data_unc_plots = output_dir + "data_uncertainty/";
  string ratio_plots = output_dir + "ratios_to_data/";
  string final_plots = output_dir + "final_xsec/";
  string check_response_dir = output_dir + "check_response/";
  string unfolding_plots = output_dir + "unfolding/";

//ntuple directories
  string ntuple_storage_path = "root://eoscms//eos/cms/store/user/";
  string ntuple_mc_path = ntuple_storage_path + "ilknur/7TeVMC/CMSSW_4_2_4/2010_MC_with_LowPU/GlobalTag_V16_RemoteGlidein/";
  string ntuple_pythia6_path = ntuple_mc_path + "pythia6/TuneZ2star_pythia6_START42_V16/";
  string ntuple_pythia8_path = ntuple_mc_path + "pythia8/Tune4C_pythia8_START42_V16/";
  string ntuple_data_path = ntuple_storage_path + "ilknur/7TeVdata/CMSSW_4_2_4/2010_DATA_with_QCD_NTuples_with_pTCut_5GeV/";
  string ntuple_trigeff_path = ntuple_storage_path + "ilknur/Trigger_Eff_Ntuple/";

//temp vars
  string prefix;


//MC Ntuples locations and luminosities

   //Pythia6_TuneZ2star
   string mc_p6_z2[5];
   double lumi_p6_z2[5];
   mc_p6_z2[0] = ntuple_pythia6_path + "QCD_Pt_10to25_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_10to25_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[0] = 3.651e9/5205128.0;
   mc_p6_z2[1] = ntuple_pythia6_path + "QCD_Pt_25to40_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_25to40_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[1] = 1.059e8/5194564.0;
   mc_p6_z2[2] = ntuple_pythia6_path + "QCD_Pt_40to80_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_40to80_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[2] = 1.771e7/5011067.0;
   mc_p6_z2[3] = ntuple_pythia6_path + "QCD_Pt_80to150_fwdJet_TuneZ2star_HFshowerLibrary_7TeV_pythia6/QCD_Pt_80to150_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[3] = 8.756e5/4876756.0;
   mc_p6_z2[4] = "../../../../../../../../work/c/cipriano/public/pythia6_150toinf/QCD_Pt_150toInf_fwdJet_TuneZ2star_7TeV_pythia6_ProcessedTree_data.root";
   lumi_p6_z2[4] = 4.76e4/5710726.0;
   string vertex_p6_z2_allvertex = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex.root";
   string vertex_p6_z2_allvertex_jetmettau_v1 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v1.root";
   string vertex_p6_z2_allvertex_jetmet_v1 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v1.root";
   string vertex_p6_z2_allvertex_jet_v1 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v1.root";
   string vertex_p6_z2_allvertex_jetmettau_v2 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v2.root";
   string vertex_p6_z2_allvertex_jetmet_v2 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v2.root";
   string vertex_p6_z2_allvertex_jet_v2 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v2.root";
   string vertex_p6_z2_allvertex_jetmettau_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v3.root";
   string vertex_p6_z2_allvertex_jetmet_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v3.root";
   string vertex_p6_z2_allvertex_jet_v3 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v3.root";
   string vertex_p6_z2_allvertex_jetmettau_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMETTau_2010A_v4.root";
   string vertex_p6_z2_allvertex_jetmet_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_JetMET_2010A_v4.root";
   string vertex_p6_z2_allvertex_jet_v4 = vertex_dir + "vertex_Pythia6_TuneZ2star_allvertex_Jet_2010B_v4.root";
   string vertex_weights_p6_z2_jetmettau_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v1.root";
   string vertex_weights_p6_z2_jetmet_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v1.root";
   string vertex_weights_p6_z2_jet_v1 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v1.root";
   string vertex_weights_p6_z2_jetmettau_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v2.root";
   string vertex_weights_p6_z2_jetmet_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v2.root";
   string vertex_weights_p6_z2_jet_v2 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v2.root";
   string vertex_weights_p6_z2_jetmettau_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v3.root";
   string vertex_weights_p6_z2_jetmet_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v3.root";
   string vertex_weights_p6_z2_jet_v3 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v3.root";
   string vertex_weights_p6_z2_jetmettau_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v4.root";
   string vertex_weights_p6_z2_jetmet_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v4.root";
   string vertex_weights_p6_z2_jet_v4 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v4.root";
   string vertex_weights_p6_z2_jetmettau_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMETTau2010A_v5.root";
   string vertex_weights_p6_z2_jetmet_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_JetMET2010A_v5.root";
   string vertex_weights_p6_z2_jet_v5 = vertex_weights_dir + "vertex_Pythia6_TuneZ2star_Jet2010B_v5.root";
   string norm_p6_z2 = pileup_norm_dir + "norm_Pythia6_TuneZ2star.root";
   string norm_p6_z2_jetmettau = pileup_norm_dir + "norm_Pythia6_TuneZ2star_JetMETTau_2010A.root";
   string norm_p6_z2_jetmet = pileup_norm_dir + "norm_Pythia6_TuneZ2star_JetMET_2010A.root";
   string norm_p6_z2_jet = pileup_norm_dir + "norm_Pythia6_TuneZ2star_Jet_2010B.root";
   string out_p6_z2_gen_nopileup = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup.root";
   string out_p6_z2_gen_allvertex = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex.root";
   string out_p6_z2_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup_JetMETTau_2010A.root"; 
   string out_p6_z2_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex_JetMETTau_2010A.root";
   string out_p6_z2_gen_nopileup_jetmet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup_JetMET_2010A.root"; 
   string out_p6_z2_gen_allvertex_jetmet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex_JetMET_2010A.root";
   string out_p6_z2_gen_nopileup_jet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_nopileup_Jet_2010B.root"; 
   string out_p6_z2_gen_allvertex_jet = mc_gen_dir + "xsec_Pythia6_TuneZ2star_gen_allvertex_Jet_2010B.root";
   string out_p6_z2_gen_norm = output_dir + "histograms/normalized_mc/xsec_p6_z2_gen.root";
   string out_p6_z2_gen_jetmettau_norm = mc_norm_dir + "xsec_p6_z2_JetMETTau_2010A_gen.root";
   string out_p6_z2_gen_jetmet_norm = mc_norm_dir + "xsec_p6_z2_JetMET_2010A_gen.root";
   string out_p6_z2_gen_jet_norm = mc_norm_dir + "xsec_p6_z2_Jet_2010B_gen.root";
   string out_p6_z2_det_1vertex = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex.root"; 
   string out_p6_z2_det_allvertex = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex.root";
   string out_p6_z2_det_1vertex_jetmettau = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex_JetMETTau_2010A.root"; 
   string out_p6_z2_det_allvertex_jetmettau = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMETTau_2010A.root";
   string out_p6_z2_det_1vertex_jetmet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex_JetMET_2010A.root"; 
   string out_p6_z2_det_allvertex_jetmet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_JetMET_2010A.root";
   string out_p6_z2_det_1vertex_jet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_1vertex_Jet_2010B.root"; 
   string out_p6_z2_det_allvertex_jet = mc_det_dir + "xsec_Pythia6_TuneZ2star_det_allvertex_Jet_2010B.root";
   string out_p6_z2_det_norm = mc_norm_dir + "xsec_p6_z2_det.root";
   string out_p6_z2_det_jetmettau_norm = mc_norm_dir + "xsec_p6_z2_JetMETTau_2010A_det.root";
   string out_p6_z2_det_jetmet_norm = mc_norm_dir + "xsec_p6_z2_JetMET_2010A_det.root";
   string out_p6_z2_det_jet_norm = mc_norm_dir + "xsec_p6_z2_Jet_2010B_det.root";
   string corr_pileup_p6_z2 = corrections_dir + "correction_pileup_p6_z2.root";
   string corr_pileup_p6_z2_jetmettau = corrections_dir + "correction_pileup_p6_z2_JetMETTau_2010A.root";
   string corr_pileup_p6_z2_jetmet = corrections_dir + "correction_pileup_p6_z2_JetMET_2010A.root";
   string corr_pileup_p6_z2_jet = corrections_dir + "correction_pileup_p6_z2_Jet_2010B.root";
   string corr_detector_p6_z2 = corrections_dir + "correction_detector_p6_z2.root";
   string corr_detector_p6_z2_jetmettau = corrections_dir + "correction_detector_p6_z2_JetMETTau_2010A.root";
   string corr_detector_p6_z2_jetmet = corrections_dir + "correction_detector_p6_z2_JetMET_2010A.root";
   string corr_detector_p6_z2_jet = corrections_dir + "correction_detector_p6_z2_Jet_2010B.root";
   string corr_final_p6_z2 = corrections_dir + "correction_final_p6_z2.root";
   string corr_final_p6_z2_jetmettau = corrections_dir + "correction_final_p6_z2_JetMETTau_2010A.root";
   string corr_final_p6_z2_jetmet = corrections_dir + "correction_final_p6_z2_JetMET_2010A.root";
   string corr_final_p6_z2_jet = corrections_dir + "correction_final_p6_z2_Jet_2010B.root";
   string out_p6_z2_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_jetmettau_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_jetmet_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_jet_unf_allvertex = response_dir + "unfresp_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_data = unfolded_dir + "Data_unfolded_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia6_TuneZ2star_allvertex.root";
   string out_p6_z2_unfolded_data_jetmettau = unfolded_dir + "Data_unfolded_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia6_TuneZ2star_JetMETTau_2010A_allvertex.root";
   string out_p6_z2_unfolded_data_jetmet = unfolded_dir + "Data_unfolded_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia6_TuneZ2star_JetMET_2010A_allvertex.root";
   string out_p6_z2_unfolded_data_jet = unfolded_dir + "Data_unfolded_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   string out_p6_z2_unfolded_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia6_TuneZ2star_Jet_2010B_allvertex.root";
   int n_files_p6_z2 = 5;

   //Pythia8_Tune4C
   string mc_p8_4c[5];
   double lumi_p8_4c[5];
   mc_p8_4c[0] = ntuple_pythia8_path + "QCD_Pt_10to25_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_10to25_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[0] = 4.334e9/9695621.0;
   mc_p8_4c[1] = ntuple_pythia8_path + "QCD_Pt_25to40_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_25to40_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[1] = 1.226e8/9949451.0;
   mc_p8_4c[2] = ntuple_pythia8_path + "QCD_Pt_40to80_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_40to80_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[2] = 2.021e7/9658287.0;
   mc_p8_4c[3] = ntuple_pythia8_path + "QCD_Pt_80to150_fwdJet_Tune4C_HFshowerLibrary_7TeV_pythia8/QCD_Pt_80to150_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[3] = 9.855e5/9117322.0;
   mc_p8_4c[4] = "../../../../../../../../work/p/pgunnel/QCD_Pt_150toInf_fwdJet_Tune4C_7TeV_pythia8_ProcessedTree_data.root";
   lumi_p8_4c[4] = 5.283e4/9761171.0;
   string vertex_p8_4c_allvertex = vertex_dir + "vertex_Pythia8_Tune4C_allvertex.root";
   string vertex_p8_4c_allvertex_jetmettau_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v1.root";
   string vertex_p8_4c_allvertex_jetmet_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v1.root";      
   string vertex_p8_4c_allvertex_jet_v1 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v1.root";
   string vertex_p8_4c_allvertex_jetmettau_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v2.root";
   string vertex_p8_4c_allvertex_jetmet_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v2.root";      
   string vertex_p8_4c_allvertex_jet_v2 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v2.root";
   string vertex_p8_4c_allvertex_jetmettau_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v3.root";
   string vertex_p8_4c_allvertex_jetmet_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v3.root";      
   string vertex_p8_4c_allvertex_jet_v3 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v3.root";
   string vertex_p8_4c_allvertex_jetmettau_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMETTau_2010A_v4.root";
   string vertex_p8_4c_allvertex_jetmet_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_JetMET_2010A_v4.root";      
   string vertex_p8_4c_allvertex_jet_v4 = vertex_dir + "vertex_Pythia8_Tune4C_allvertex_Jet_2010B_v4.root";
   string vertex_weights_p8_4c_jetmettau_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v1.root";
   string vertex_weights_p8_4c_jetmet_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v1.root";
   string vertex_weights_p8_4c_jet_v1 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v1.root";
   string vertex_weights_p8_4c_jetmettau_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v2.root";
   string vertex_weights_p8_4c_jetmet_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v2.root"; 
   string vertex_weights_p8_4c_jet_v2 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v2.root";
   string vertex_weights_p8_4c_jetmettau_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v3.root";
   string vertex_weights_p8_4c_jetmet_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v3.root"; 
   string vertex_weights_p8_4c_jet_v3 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v3.root";
   string vertex_weights_p8_4c_jetmettau_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v4.root";
   string vertex_weights_p8_4c_jetmet_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v4.root";
   string vertex_weights_p8_4c_jet_v4 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v4.root";
   string vertex_weights_p8_4c_jetmettau_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMETTau2010A_v5.root";
   string vertex_weights_p8_4c_jetmet_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_JetMET2010A_v5.root";
   string vertex_weights_p8_4c_jet_v5 = vertex_weights_dir + "vertex_Pythia8_Tune4C_Jet2010B_v5.root";
   string norm_p8_4c = pileup_norm_dir + "norm_Pythia8_Tune4C.root";
   string norm_p8_4c_jetmettau = pileup_norm_dir + "norm_Pythia8_Tune4C_JetMETTau_2010A.root";
   string norm_p8_4c_jetmet = pileup_norm_dir + "norm_Pythia8_Tune4C_JetMET_2010A.root";
   string norm_p8_4c_jet = pileup_norm_dir + "norm_Pythia8_Tune4C_Jet_2010B.root";
   string out_p8_4c_gen_nopileup = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup.root";
   string out_p8_4c_gen_allvertex = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex.root";
   string out_p8_4c_gen_nopileup_jetmettau = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup_JetMETTau_2010A.root";
   string out_p8_4c_gen_allvertex_jetmettau = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex_JetMETTau_2010A.root";
   string out_p8_4c_gen_nopileup_jetmet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup_JetMET_2010A.root";
   string out_p8_4c_gen_allvertex_jetmet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex_JetMET_2010A.root";
   string out_p8_4c_gen_nopileup_jet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_nopileup_Jet_2010B.root";
   string out_p8_4c_gen_allvertex_jet = mc_gen_dir + "xsec_Pythia8_Tune4C_gen_allvertex_Jet_2010B.root";
   string out_p8_4c_gen_norm = mc_norm_dir + "xsec_p8_4c_gen.root";
   string out_p8_4c_gen_jetmettau_norm = mc_norm_dir + "xsec_p8_4c_JetMETTau_2010A_gen.root";
   string out_p8_4c_gen_jetmet_norm = mc_norm_dir + "xsec_p8_4c_JetMET_2010A_gen.root";
   string out_p8_4c_gen_jet_norm = mc_norm_dir + "xsec_p8_4c_Jet_2010B_gen.root";
   string out_p8_4c_det_1vertex = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex.root";
   string out_p8_4c_det_allvertex = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex.root";
   string out_p8_4c_det_1vertex_jetmettau = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex_JetMETTau_2010A.root";
   string out_p8_4c_det_allvertex_jetmettau = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMETTau_2010A.root";
   string out_p8_4c_det_1vertex_jetmet = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex_JetMET_2010A.root";
   string out_p8_4c_det_allvertex_jetmet = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_JetMET_2010A.root";
   string out_p8_4c_det_1vertex_jet = mc_det_dir + "xsec_Pythia8_Tune4C_det_1vertex_Jet_2010B.root";
   string out_p8_4c_det_allvertex_jet = mc_det_dir + "xsec_Pythia8_Tune4C_det_allvertex_Jet_2010B.root";
   string out_p8_4c_det_norm = mc_norm_dir + "xsec_p8_4c_det.root";
   string out_p8_4c_det_jetmettau_norm = mc_norm_dir + "xsec_p8_4c_JetMETTau_2010A_det.root";
   string out_p8_4c_det_jetmet_norm = mc_norm_dir + "xsec_p8_4c_JetMET_2010A_det.root";
   string out_p8_4c_det_jet_norm = mc_norm_dir + "xsec_p8_4c_Jet_2010B_det.root";
   string corr_pileup_p8_4c = corrections_dir + "correction_pileup_p8_4c.root";
   string corr_pileup_p8_4c_jetmettau = corrections_dir + "correction_pileup_p8_4c_JetMETTau_2010A.root";
   string corr_pileup_p8_4c_jetmet = corrections_dir + "correction_pileup_p8_4c_JetMET_2010A.root";
   string corr_pileup_p8_4c_jet = corrections_dir + "correction_pileup_p8_4c_Jet_2010B.root";
   string corr_detector_p8_4c = corrections_dir + "correction_detector_p8_4c.root";
   string corr_detector_p8_4c_jetmettau = corrections_dir + "correction_detector_p8_4c_JetMETTau_2010A.root";
   string corr_detector_p8_4c_jetmet = corrections_dir + "correction_detector_p8_4c_JetMET_2010A.root";
   string corr_detector_p8_4c_jet = corrections_dir + "correction_detector_p8_4c_Jet_2010B.root";
   string corr_final_p8_4c = corrections_dir + "correction_final_p8_4c.root";
   string corr_final_p8_4c_jetmettau = corrections_dir + "correction_final_p8_4c_JetMETTau_2010A.root";
   string corr_final_p8_4c_jetmet = corrections_dir + "correction_final_p8_4c_JetMET_2010A.root";
   string corr_final_p8_4c_jet = corrections_dir + "correction_final_p8_4c_Jet_2010B.root";
   string out_p8_4c_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_jetmettau_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_jetmet_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_jet_unf_allvertex = response_dir + "unfresp_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_data = unfolded_dir + "Data_unfolded_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_p6_z2 = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_p8_4c = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia8_Tune4C_allvertex.root";
   string out_p8_4c_unfolded_data_jetmettau = unfolded_dir + "Data_unfolded_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_p6_z2_jetmettau = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_p8_4c_jetmettau = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia8_Tune4C_JetMETTau_2010A_allvertex.root";
   string out_p8_4c_unfolded_data_jetmet = unfolded_dir + "Data_unfolded_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_p6_z2_jetmet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_p8_4c_jetmet = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia8_Tune4C_JetMET_2010A_allvertex.root";
   string out_p8_4c_unfolded_data_jet = unfolded_dir + "Data_unfolded_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_p6_z2_jet = unfolded_dir + "Pythia6_TuneZ2star_unfolded_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   string out_p8_4c_unfolded_p8_4c_jet = unfolded_dir + "Pythia8_Tune4C_unfolded_Pythia8_Tune4C_Jet_2010B_allvertex.root";
   int n_files_p8_4c = 4;
   

//Data Ntuples locations and luminosities

//JetMETTauMonitor_2010A
   string jetmettaumon[1];
   double lumi_jetmettaumon[2];
   //jetmettaumon[0] = ntuple_trigeff_path + "JetTreeJETMETTAUMONITOR.root";
   jetmettaumon[0] = "../../../../../../../../work/c/cipriano/public/JetTreeJETMETTAUMONITOR.root";
   string systematics_trigger_jetmettaumon = trigger_syst_dir + "systematic_trigger_JetMETTauMonitor_2010A.root";
   string out_jetmettaumon_HLT_L1Jet6U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet6U.root";
   string out_jetmettaumon_HLT_L1Jet10U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet10U.root";
   string out_jetmettaumon_HLT_Jet15U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_Jet15U.root";
   string out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string out_jetmettaumon_HLT_Jet30U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_Jet30U.root";
   string out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_Jet15U_AND_HLT_Jet30U.root";
   string out_jetmettaumon_HLT_Jet50U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_Jet50U.root";
   string out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_Jet15U_AND_HLT_Jet50U.root";
   string out_jetmettaumon_HLT_DiJetAve15U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_DiJetAve15U.root";
   string out_jetmettaumon_HLT_DiJetAve30U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_DiJetAve30U.root";
   string out_jetmettaumon_HLT_FwdJet20U = trigger_dir + "events_JetMETTauMonitor2010A_HLT_FwdJet20U.root";
   string out_jetmettaumon_HLT_DoubleJet15U_ForwardBackward = trigger_dir + "events_JetMETTauMonitor2010A_HLT_DoubleJet15U_ForwardBackward.root";
   string out_jetmettaumon_all = trigger_dir + "events_JetMETTauMonitor2010A_all.root";
   string out_jetmettaumon_none = trigger_dir + "events_JetMETTauMonitor2010A_none.root";
   string eff_jetmettaumon_HLT_L1Jet6U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet6U.root";
   string eff_jetmettaumon_HLT_L1Jet10U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet10U.root";
   string eff_jetmettaumon_HLT_Jet15U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_Jet15U.root";
   string eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string eff_jetmettaumon_HLT_Jet30U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_Jet30U.root";
   string eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string eff_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_Jet15U_AND_HLT_Jet30U.root";
   string eff_jetmettaumon_HLT_Jet50U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_Jet50U.root";
   string eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string eff_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_Jet15U_AND_HLT_Jet50U.root";
   string eff_jetmettaumon_HLT_DiJetAve15U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_DiJetAve15U.root";
   string eff_jetmettaumon_HLT_DiJetAve30U = eff_dir + "efficiency_JetMETTauMonitor2010A_HLT_DiJetAve30U.root";
   string eff_jetmettaumon_HLT_FwdJet20U = eff_dir + "efficiency__JetMETTauMonitor2010A_HLT_FwdJet20U.root";
   string eff_jetmettaumon_HLT_DoubleJet15U_ForwardBackward = eff_dir + "efficiency__JetMETTauMonitor2010A_HLT_DoubleJet15U_ForwardBackward.root";
   string eff_jetmettaumon_all = eff_dir + "efficiency_JetMETTauMonitor2010A_all.root";
   string out_jetmettaumon_emu = trigger_dir + "events_JetMETTauMonitor2010A_emulated.root";
   string correction_jetmettaumon_emu_v1 = trigger_corr_dir + "correction_JetMETTauMonitor2010A_emulated_v1.txt";
   string correction_sj_jetmettaumon_emu_v1 = trigger_corr_dir + "correction_sj_JetMETTauMonitor2010A_emulated_v1.txt";
   string check_jetmettaumon_emu_v1 = trigger_dir + "check_JetMETTauMonitor2010A_emulated_v1.root";
   string check_sj_jetmettaumon_emu_v1 = trigger_dir + "check_sj_JetMETTauMonitor2010A_emulated_v1.root";
   string correction_jetmettaumon_emu_v2 = trigger_corr_dir + "correction_JetMETTauMonitor2010A_emulated_v2.txt";
   string correction_sj_jetmettaumon_emu_v2 = trigger_corr_dir + "correction_sj_JetMETTauMonitor2010A_emulated_v2.txt";
   string check_jetmettaumon_emu_v2 = trigger_dir + "check_JetMETTauMonitor2010A_emulated_v2.root";
   string check_sj_jetmettaumon_emu_v2 = trigger_dir + "check_sj_JetMETTauMonitor2010A_emulated_v2.root";
   int n_files_jetmettaumon = 1;
   lumi_jetmettaumon[0] = 1.0; // /pb cipriano

//JetMETTau_2010A
   string jetmettau[3];
   double lumi_jetmettau[3];
   jetmettau[0] = ntuple_data_path + "JetMETTau_Run2010A_Apr21ReReco_v1_AOD/JetMETTau_Run2010A_Apr21ReReco_v1_AOD_ProcessedTree_data.root"; //new
   jetmettau[1] = jetmettau[0];
   jetmettau[2] = jetmettau[0];
   string systematics_trigger_jetmettau = trigger_syst_dir + "systematic_trigger_JetMETTau_2010A.root";
   string vertex_jetmettau_allvertex = vertex_dir + "vertex_JetMETTau2010A_allvertex.root";
   string out_jetmettau_1vertex = raw_data_dir + "xsec_JetMETTau2010A_1vertex.root";
   string out_jetmettau_allvertex = raw_data_dir + "xsec_JetMETTau2010A_allvertex.root";
   string out_jetmettau_up = raw_data_dir + "xsec_JetMETTau2010A_up.root";
   string out_jetmettau_down = raw_data_dir + "xsec_JetMETTau2010A_down.root";
   string out_jetmettau_HLT_L1Jet6U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet6U.root";
   string out_jetmettau_HLT_L1Jet10U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet10U.root";
   string out_jetmettau_HLT_Jet15U = trigger_dir + "events_JetMETTau2010A_HLT_Jet15U.root";
   string out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string out_jetmettau_HLT_Jet30U = trigger_dir + "events_JetMETTau2010A_HLT_Jet30U.root";
   string out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string out_jetmettau_HLT_Jet15U_AND_HLT_Jet30U = trigger_dir + "events_JetMETTau2010A_HLT_Jet15U_AND_HLT_Jet30U.root";
   string out_jetmettau_HLT_Jet50U = trigger_dir + "events_JetMETTau2010A_HLT_Jet50U.root";
   string out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U = trigger_dir + "events_JetMETTau2010A_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string out_jetmettau_HLT_Jet15U_AND_HLT_Jet50U = trigger_dir + "events_JetMETTau2010A_HLT_Jet15U_AND_HLT_Jet50U.root";
   string out_jetmettau_HLT_DiJetAve15U = trigger_dir + "events_JetMETTau2010A_HLT_DiJetAve15U.root";
   string out_jetmettau_HLT_DiJetAve30U = trigger_dir + "events_JetMETTau2010A_HLT_DiJetAve30U.root";
   string out_jetmettau_HLT_FwdJet20U = trigger_dir + "events_JetMETTau2010A_HLT_FwdJet20U.root";
   string out_jetmettau_HLT_DoubleJet15U_ForwardBackward = trigger_dir + "events_JetMETTau2010A_HLT_DoubleJet15U_ForwardBackward.root";
   string out_jetmettau_all = trigger_dir + "events_JetMETTau2010A_all.root";
   string out_jetmettau_none = trigger_dir + "events_JetMETTau2010A_none.root";
   string eff_jetmettau_HLT_L1Jet6U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet6U.root";
   string eff_jetmettau_HLT_L1Jet10U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet10U.root";
   string eff_jetmettau_HLT_Jet15U = eff_dir + "efficiency_JetMETTau2010A_HLT_Jet15U.root";
   string eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string eff_jetmettau_HLT_Jet30U = eff_dir + "efficiency_JetMETTau2010A_HLT_Jet30U.root";
   string eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string eff_jetmettau_HLT_Jet15U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMETTau2010A_HLT_Jet15U_AND_HLT_Jet30U.root";
   string eff_jetmettau_HLT_Jet50U = eff_dir + "efficiency_JetMETTau2010A_HLT_Jet50U.root";
   string eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMETTau2010A_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string eff_jetmettau_HLT_Jet15U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMETTau2010A_HLT_Jet15U_AND_HLT_Jet50U.root";
   string eff_jetmettau_HLT_DiJetAve15U = eff_dir + "efficiency_JetMETTau2010A_HLT_DiJetAve15U.root";
   string eff_jetmettau_HLT_DiJetAve30U = eff_dir + "efficiency_JetMETTau2010A_HLT_DiJetAve30U.root";
   string eff_jetmettau_HLT_FwdJet20U = eff_dir + "efficiency__JetMETTau2010A_HLT_FwdJet20U.root";
   string eff_jetmettau_HLT_DoubleJet15U_ForwardBackward = eff_dir + "efficiency__JetMETTau2010A_HLT_DoubleJet15U_ForwardBackward.root";
   string eff_jetmettau_all = eff_dir + "efficiency_JetMETTau2010A_all.root";
   string out_jetmettau_emu = trigger_dir + "events_JetMETTau2010A_emulated.root";
   string correction_jetmettau_emu_v1 = trigger_corr_dir + "correction_JetMETTau2010A_emulated_v1.txt";
   string correction_sj_jetmettau_emu_v1 = trigger_corr_dir + "correction_sj_JetMETTau2010A_emulated_v1.txt";
   string check_jetmettau_emu_v1 = trigger_dir + "check_JetMETTau2010A_emulated_v1.root";
   string check_sj_jetmettau_emu_v1 = trigger_dir + "check_sj_JetMETTau2010A_emulated_v1.root";
   string correction_jetmettau_emu_v2 = trigger_corr_dir + "correction_JetMETTau2010A_emulated_v2.txt";
   string correction_sj_jetmettau_emu_v2 = trigger_corr_dir + "correction_sj_JetMETTau2010A_emulated_v2.txt";
   string check_jetmettau_emu_v2 = trigger_dir + "check_JetMETTau2010A_emulated_v2.root";
   string check_sj_jetmettau_emu_v2 = trigger_dir + "check_sj_JetMETTau2010A_emulated_v2.root";
   string jes_up_unc_jetmettau = jes_unc_dir + "jes_up_unc_JetMETTau_2010A.root";
   string jes_down_unc_jetmettau = jes_unc_dir + "jes_down_unc_JetMETTau_2010A.root";
   int n_files_jetmettau = 3;
   //lumi_jetmettau[0] = 0.284161; // /pb cipriano
   lumi_jetmettau[0] = 0.009645;
   lumi_jetmettau[1] = 0.192895;
   lumi_jetmettau[2] = 2.869;
   
   //JetMET_2010A
   string jetmet[3];
   double lumi_jetmet[3];
   jetmet[0] = ntuple_data_path + "JetMET_Run2010A_Apr21ReReco_v1_AOD/JETMET_AOD_21Apr_MERGE_ALL_ROOT_ProcessedTree_data.root";
   jetmet[1] = jetmet[0];
   jetmet[2] = jetmet[0];
   string systematics_trigger_jetmet = trigger_syst_dir + "systematic_trigger_JetMET_2010A.root";
   string vertex_jetmet_allvertex = vertex_dir + "vertex_JetMET2010A_allvertex.root";
   string out_jetmet_1vertex = raw_data_dir + "xsec_JetMET2010A_1vertex.root";
   string out_jetmet_allvertex = raw_data_dir + "xsec_JetMET2010A_allvertex.root";
   string out_jetmet_up = raw_data_dir + "xsec_JetMET2010A_up.root";
   string out_jetmet_down = raw_data_dir + "xsec_JetMET2010A_down.root";
   string out_jetmet_HLT_L1Jet6U = trigger_dir + "events_JetMET2010A_HLT_L1Jet6U.root";
   string out_jetmet_HLT_L1Jet10U = trigger_dir + "events_JetMET2010A_HLT_L1Jet10U.root";
   string out_jetmet_HLT_Jet15U = trigger_dir + "events_JetMET2010A_HLT_Jet15U.root";
   string out_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U = trigger_dir + "events_JetMET2010A_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string out_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U = trigger_dir + "events_JetMET2010A_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string out_jetmet_HLT_Jet30U = trigger_dir + "events_JetMET2010A_HLT_Jet30U.root";
   string out_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U = trigger_dir + "events_JetMET2010A_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string out_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U = trigger_dir + "events_JetMET2010A_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string out_jetmet_HLT_Jet15U_AND_HLT_Jet30U = trigger_dir + "events_JetMET2010A_HLT_Jet15U_AND_HLT_Jet30U.root";
   string out_jetmet_HLT_Jet50U = trigger_dir + "events_JetMET2010A_HLT_Jet50U.root";
   string out_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U = trigger_dir + "events_JetMET2010A_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string out_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U = trigger_dir + "events_JetMET2010A_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string out_jetmet_HLT_Jet15U_AND_HLT_Jet50U = trigger_dir + "events_JetMET2010A_HLT_Jet15U_AND_HLT_Jet50U.root";
   string out_jetmet_HLT_DiJetAve15U = trigger_dir + "events_JetMET2010A_HLT_DiJetAve15U.root";
   string out_jetmet_HLT_DiJetAve30U = trigger_dir + "events_JetMET2010A_HLT_DiJetAve30U.root";
   string out_jetmet_HLT_FwdJet20U = trigger_dir + "events_JetMET2010A_HLT_FwdJet20U.root";
   string out_jetmet_HLT_DoubleJet15U_ForwardBackward = trigger_dir + "events_JetMET2010A_HLT_DoubleJet15U_ForwardBackward.root";
   string out_jetmet_all = trigger_dir + "events_JetMET2010A_all.root";
   string out_jetmet_none = trigger_dir + "events_JetMET2010A_none.root";
   string eff_jetmet_HLT_L1Jet6U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet6U.root";
   string eff_jetmet_HLT_L1Jet10U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet10U.root";
   string eff_jetmet_HLT_Jet15U = eff_dir + "efficiency_JetMET2010A_HLT_Jet15U.root";
   string eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string eff_jetmet_HLT_Jet30U = eff_dir + "efficiency_JetMET2010A_HLT_Jet30U.root";
   string eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string eff_jetmet_HLT_Jet15U_AND_HLT_Jet30U = eff_dir + "efficiency_JetMET2010A_HLT_Jet15U_AND_HLT_Jet30U.root";
   string eff_jetmet_HLT_Jet50U = eff_dir + "efficiency_JetMET2010A_HLT_Jet50U.root";
   string eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMET2010A_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string eff_jetmet_HLT_Jet15U_AND_HLT_Jet50U = eff_dir + "efficiency_JetMET2010A_HLT_Jet15U_AND_HLT_Jet50U.root";
   string eff_jetmet_HLT_DiJetAve15U = eff_dir + "efficiency_JetMET2010A_HLT_DiJetAve15U.root";
   string eff_jetmet_HLT_DiJetAve30U = eff_dir + "efficiency_JetMET2010A_HLT_DiJetAve30U.root";
   string eff_jetmet_HLT_FwdJet20U = eff_dir + "efficiency_JetMET2010A_HLT_FwdJet20U.root";
   string eff_jetmet_HLT_DoubleJet15U_ForwardBackward = eff_dir + "efficiency_JetMET2010A_HLT_DoubleJet15U_ForwardBackward.root";
   string eff_jetmet_all = eff_dir + "efficiency_JetMET2010A_all.root";
   string out_jetmet_emu = trigger_dir + "events_JetMET2010A_emulated.root";
   string correction_jetmet_emu_v1 = trigger_corr_dir + "correction_JetMET2010A_emulated_v1.txt";
   string correction_sj_jetmet_emu_v1 = trigger_corr_dir + "correction_sj_JetMET2010A_emulated_v1.txt";
   string check_jetmet_emu_v1 = trigger_dir + "check_JetMET2010A_emulated_v1.root";
   string check_sj_jetmet_emu_v1 = trigger_dir + "check_sj_JetMET2010A_emulated_v1.root";
   string correction_jetmet_emu_v2 = trigger_corr_dir + "correction_JetMET2010A_emulated_v2.txt";
   string correction_sj_jetmet_emu_v2 = trigger_corr_dir + "correction_sj_JetMET2010A_emulated_v2.txt";
   string check_jetmet_emu_v2 = trigger_dir + "check_JetMET2010A_emulated_v2.root";
   string check_sj_jetmet_emu_v2 = trigger_dir + "check_sj_JetMET2010A_emulated_v2.root";
   string jes_up_unc_jetmet = jes_unc_dir + "jes_up_unc_JetMET_2010A.root";
   string jes_down_unc_jetmet = jes_unc_dir + "jes_down_unc_JetMET_2010A.root";
   int n_files_jetmet = 3;
   //lumi_jetmet[0] = 2.896; // /pb cipriano
   lumi_jetmet[0] = 0.013799;
   lumi_jetmet[1] = 0.117223;
   lumi_jetmet[2] = 0.278789;
   
   //Jet_2010B
   string jet[3];
   double lumi_jet[3];
   string jet[0] = ntuple_data_path + "Jet_Run2010B_Apr21ReReco_v1_AOD_New/JET_2010_AOD_21APR_MERGE_ALL_ROOT_ProcessedTree_data.root";
   jet[1] = jet[0];
   jet[2] = jet[0];   
   string systematics_trigger_jet = trigger_syst_dir + "systematic_trigger_Jet2010B.root";
   string vertex_jet_allvertex = vertex_dir + "vertex_Jet2010B_allvertex.root";
   string out_jet_1vertex = raw_data_dir + "xsec_Jet2010B_1vertex.root";
   string out_jet_allvertex = raw_data_dir + "xsec_Jet2010B_allvertex.root";
   string out_jet_up = raw_data_dir + "xsec_Jet2010B_up.root";
   string out_jet_down = raw_data_dir + "xsec_Jet2010B_down.root";
   string out_jet_HLT_L1Jet6U = trigger_dir + "events_Jet2010B_HLT_L1Jet6U.root";
   string out_jet_HLT_L1Jet10U = trigger_dir + "events_Jet2010B_HLT_L1Jet10U.root";
   string out_jet_HLT_Jet15U = trigger_dir + "events_Jet2010B_HLT_Jet15U.root";
   string out_jet_HLT_L1Jet6U_AND_HLT_Jet15U = trigger_dir + "events_Jet2010B_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string out_jet_HLT_L1Jet10U_AND_HLT_Jet15U = trigger_dir + "events_Jet2010B_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string out_jet_HLT_Jet30U = trigger_dir + "events_Jet2010B_HLT_Jet30U.root";
   string out_jet_HLT_L1Jet6U_AND_HLT_Jet30U = trigger_dir + "events_Jet2010B_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string out_jet_HLT_L1Jet10U_AND_HLT_Jet30U = trigger_dir + "events_Jet2010B_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string out_jet_HLT_Jet15U_AND_HLT_Jet30U = trigger_dir + "events_Jet2010B_HLT_Jet15U_AND_HLT_Jet30U.root";
   string out_jet_HLT_Jet50U = trigger_dir + "events_Jet2010B_HLT_Jet50U.root";
   string out_jet_HLT_L1Jet6U_AND_HLT_Jet50U = trigger_dir + "events_Jet2010B_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string out_jet_HLT_L1Jet10U_AND_HLT_Jet50U = trigger_dir + "events_Jet2010B_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string out_jet_HLT_Jet15U_AND_HLT_Jet50U = trigger_dir + "events_Jet2010B_HLT_Jet15U_AND_HLT_Jet50U.root";
   string out_jet_HLT_DiJetAve15U = trigger_dir + "events_Jet2010B_HLT_DiJetAve15U.root";
   string out_jet_HLT_DiJetAve30U = trigger_dir + "events_Jet2010B_HLT_DiJetAve30U.root";
   string out_jet_HLT_FwdJet20U = trigger_dir + "events_Jet2010B_HLT_FwdJet20U.root";
   string out_jet_HLT_DoubleJet15U_ForwardBackward = trigger_dir + "events_Jet2010B_HLT_DoubleJet15U_ForwardBackward.root";
   string out_jet_all = trigger_dir + "events_Jet2010B_all.root";
   string out_jet_none = trigger_dir + "events_Jet2010B_none.root";
   string eff_jet_HLT_L1Jet6U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet6U.root";
   string eff_jet_HLT_L1Jet10U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet10U.root";
   string eff_jet_HLT_Jet15U = eff_dir + "efficiency_Jet2010B_HLT_Jet15U.root";
   string eff_jet_HLT_L1Jet6U_AND_HLT_Jet15U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet6U_AND_HLT_Jet15U.root";
   string eff_jet_HLT_L1Jet10U_AND_HLT_Jet15U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet10U_AND_HLT_Jet15U.root";
   string eff_jet_HLT_Jet30U = eff_dir + "efficiency_Jet2010B_HLT_Jet30U.root";
   string eff_jet_HLT_L1Jet6U_AND_HLT_Jet30U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet6U_AND_HLT_Jet30U.root";
   string eff_jet_HLT_L1Jet10U_AND_HLT_Jet30U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet10U_AND_HLT_Jet30U.root";
   string eff_jet_HLT_Jet15U_AND_HLT_Jet30U = eff_dir + "efficiency_Jet2010B_HLT_Jet15U_AND_HLT_Jet30U.root";
   string eff_jet_HLT_Jet50U = eff_dir + "efficiency_Jet2010B_HLT_Jet50U.root";
   string eff_jet_HLT_L1Jet6U_AND_HLT_Jet50U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet6U_AND_HLT_Jet50U.root";
   string eff_jet_HLT_L1Jet10U_AND_HLT_Jet50U = eff_dir + "efficiency_Jet2010B_HLT_L1Jet10U_AND_HLT_Jet50U.root";
   string eff_jet_HLT_Jet15U_AND_HLT_Jet50U = eff_dir + "efficiency_Jet2010B_HLT_Jet15U_AND_HLT_Jet50U.root";
   string eff_jet_HLT_DiJetAve15U = eff_dir + "efficiency_Jet2010B_HLT_DiJetAve15U.root";
   string eff_jet_HLT_DiJetAve30U = eff_dir + "efficiency_Jet2010B_HLT_DiJetAve30U.root";
   string eff_jet_HLT_FwdJet20U = eff_dir + "efficiency_Jet2010B_HLT_FwdJet20U.root";
   string eff_jet_HLT_DoubleJet15U_ForwardBackward = eff_dir + "efficiency_Jet2010B_HLT_DoubleJet15U_ForwardBackward.root";
   string eff_jet_all = eff_dir + "efficiency_Jet2010B_all.root";
   string out_jet_emu = trigger_dir + "events_Jet2010B_emulated.root";
   string correction_jet_emu_v1 = trigger_corr_dir + "correction_Jet2010B_emulated_v1.txt";
   string correction_sj_jet_emu_v1 = trigger_corr_dir + "correction_sj_Jet2010B_emulated_v1.txt";
   string check_jet_emu_v1 = trigger_dir + "check_Jet2010B_emulated_v1.root";
   string check_sj_jet_emu_v1 = trigger_dir + "check_sj_Jet2010B_emulated_v1.root";
   string correction_jet_emu_v2 = trigger_corr_dir + "correction_Jet2010B_emulated_v2.txt";
   string correction_sj_jet_emu_v2 = trigger_corr_dir + "correction_sj_Jet2010B_emulated_v2.txt";
   string check_jet_emu_v2 = trigger_dir + "check_Jet2010B_emulated_v2.root";
   string check_sj_jet_emu_v2 = trigger_dir + "check_sj_Jet2010B_emulated_v2.root";
   string jes_up_unc_jet = jes_unc_dir + "jes_up_unc_Jet_2010B.root";
   string jes_down_unc_jet = jes_unc_dir + "jes_down_unc_Jet_2010B.root";
   int n_files_jet = 3;
   //lumi_jet[0] = 32.072; // /pb cipriano
   lumi_jet[0] = 0.001855;
   lumi_jet[1] = 0.026783;
   lumi_jet[2] = 0.239874;

//histogram locations
   string corrected_jetmettau = corrected_dir + "corrected_xsec_data_JetMETTau_2010A.root";
   string corrected_jetmet = corrected_dir + "corrected_xsec_data_JetMET_2010A.root";
   string corrected_jet = corrected_dir + "corrected_xsec_data_Jet_2010B.root";

   string merged_data = histograms_dir + "xsec_data_2010.root";
   string ratio_merged_data = ratio_dir + "xsec_data_2010.root";
   string data_unc = histograms_dir + "xsec_data_2010_unc.root";
   string ratio_data_unc = ratio_dir + "xsec_data_2010_unc.root";
   string merged_trig_comb_syst = trigger_syst_dir + "trig_comb_syst_data_2010.root";

   string model_unc_jetmettau = model_unc_dir + "model_uncertainty_JetMETTau_2010A.root";
   string model_unc_jetmet = model_unc_dir + "model_uncertainty_JetMET_2010A.root";
   string model_unc_jet = model_unc_dir + "model_uncertainty_Jet_2010B.root";
   string model_unc_merged = model_unc_dir + "model_uncertainty.root";
   string jes_unc_up_merged = jes_unc_dir + "jes_uncertainty_up.root";
   string jes_unc_down_merged = jes_unc_dir + "jes_uncertainty_down.root";
   string total_unc_up = total_unc_dir + "total_uncertainty_up.root";
   string total_unc_down = total_unc_dir + "total_uncertainty_down.root";

//monte carlo predictions locations
   string albert_public_dir = "../../../../../../../../work/k/knutsson/public/mcforpedro/";

   string mc_pred_herwig6 = albert_public_dir + "1M/herwig6pt30.root";
   string mc_pred_herwigpp = albert_public_dir + "1M/herwigpppt30.root";
   string mc_pred_py6_ambt1 = albert_public_dir + "1M/ambt1pt30.root";
   string mc_pred_py6_p11 = albert_public_dir + "1M/py6p11pt30.root";
   string mc_pred_py6_z2_nompi = albert_public_dir + "1M/z2nompipt30.root";
   string mc_pred_py6_z2 = albert_public_dir + "1M/z2starpt30.root";
   string mc_pred_py8_4c = albert_public_dir + "1M/py84cpt30.root";

   string mc_ratio_herwig6 = ratio_dir + "herwig6_pt30.root";
   string mc_ratio_herwigpp = ratio_dir + "herwigpp_pt30.root";
   string mc_ratio_py6_ambt1 = ratio_dir + "pythia6_ambt1_pt30.root";
   string mc_ratio_py6_p11 = ratio_dir + "pythia6_p11_pt30.root";
   string mc_ratio_py6_z2_nompi = ratio_dir + "pythia6_z2starnompi_pt30.root";
   string mc_ratio_py6_z2 = ratio_dir + "pythia6_z2star_pt30.root";
   string mc_ratio_py8_4c = ratio_dir + "pythia8_4c_pt30.root";

//predictions list
   string final_mc_list[10];
   string final_label_list[10];
   string final_prefix_list[10];

   final_mc_list[0] = out_p6_z2_gen_allvertex_jetmet;
   final_mc_list[1] = out_p8_4c_gen_allvertex_jetmet;
   final_mc_list[2] = mc_pred_herwig6;
   final_mc_list[3] = mc_pred_herwigpp;
   final_mc_list[4] = mc_pred_py6_p11;
   final_mc_list[5] = mc_pred_py6_z2_nompi;
   final_mc_list[6] = mc_pred_py6_z2;
   final_mc_list[7] = mc_pred_py8_4c;
   final_mc_list[8] = "";
   final_mc_list[9] = "";

   final_label_list[0] = "Pythia 6 - Z2* Tune (Excess)";
   final_label_list[1] = "Pythia 8 - 4C Tune (Excess)";
   final_label_list[2] = "Herwig 6";
   final_label_list[3] = "Herwig pp";
   final_label_list[4] = "Pythia 6 - P11 Tune";
   final_label_list[5] = "Pythia 6 - Z2* Tune (No MPI)";
   final_label_list[6] = "Pythia 6 - Z2* Tune";
   final_label_list[7] = "Pythia 8 - 4C";
   final_label_list[8] = "";
   final_label_list[9] = "";

   final_prefix_list[0] = "ak5Gen_";
   final_prefix_list[1] = "ak5Gen_";
   final_prefix_list[2] = "";
   final_prefix_list[3] = "";
   final_prefix_list[4] = "";
   final_prefix_list[5] = "";
   final_prefix_list[6] = "";
   final_prefix_list[7] = "";
   final_prefix_list[8] = "";
   final_prefix_list[9] = "";


//rations list
   string ratio_mc_list[10];
   string ratio_label_list[10];
   string ratio_prefix_list[10];

   ratio_mc_list[0] = out_p6_z2_gen_allvertex_jetmet;
   ratio_mc_list[1] = out_p8_4c_gen_allvertex_jetmet;
   ratio_mc_list[2] = mc_ratio_herwig6;
   ratio_mc_list[3] = mc_ratio_herwigpp;
   ratio_mc_list[4] = mc_pred_py6_p11;
   ratio_mc_list[5] = mc_pred_py6_z2_nompi;
   ratio_mc_list[6] = mc_pred_py6_z2;
   ratio_mc_list[7] = mc_pred_py8_4c;
   ratio_mc_list[8] = "";
   ratio_mc_list[9] = "";

   ratio_label_list[0] = "Pythia 6 - Z2* Tune (Excess)";
   ratio_label_list[1] = "Pythia 8 - 4C Tune (Excess)";
   ratio_label_list[2] = "Herwig 6";
   ratio_label_list[3] = "Herwig pp";
   ratio_label_list[4] = "Pythia 6 - P11 Tune";
   ratio_label_list[5] = "Pythia 6 - Z2* Tune (No MPI)";
   ratio_label_list[6] = "Pythia 6 - Z2* Tune";
   ratio_label_list[7] = "Pythia 8 - 4C";
   ratio_label_list[8] = "";
   ratio_label_list[9] = "";

   ratio_prefix_list[0] = "ak5Gen_";
   ratio_prefix_list[1] = "ak5Gen_";
   ratio_prefix_list[2] = "";
   ratio_prefix_list[3] = "";
   ratio_prefix_list[4] = "";
   ratio_prefix_list[5] = "";
   ratio_prefix_list[6] = "";
   ratio_prefix_list[7] = "";
   ratio_prefix_list[8] = "";
   ratio_prefix_list[9] = "";

//compile the different routines
if (show_steps) { cout << "Compiling code..."<<endl; }
  gROOT -> ProcessLine(".L create_directories.C++");
  gROOT -> ProcessLine(".L common_methods.h++");
  if (estimate_combination_systematic) { gROOT -> ProcessLine(".L estimate_combination_systematic.C++"); }
  if (get_trigger_efficiencies) { gROOT -> ProcessLine(".L get_triggered_events.C++"); }
  if (compute_trigger_efficiency) { gROOT -> ProcessLine(".L compute_trigger_efficiency.C++"); }
  if (compute_trigger_turn_on || check_trigger_correction_v1 || check_trigger_correction_v2 || apply_trigger_correction || plot_final_trigger_efficiency)
  {
  gROOT -> ProcessLine(".L compute_trigger_turn_on.C++");
  }
  if (counting_combination_statistics) { gROOT -> ProcessLine(".L counting_combination_statistics.C++"); }
  if (compute_trigger_correction_v1 || compute_trigger_correction_v2) { gROOT -> ProcessLine(".L compute_trigger_correction.C++"); }
  if (get_vertex_distribution_v0 || get_vertex_distribution_v1 || get_vertex_distribution_v2 || get_vertex_distribution_v3 || get_vertex_distribution_v4)
  {
  gROOT -> ProcessLine(".L get_vertex_distribution.C++");
  }
  if (compute_vertex_weights_v1 || compute_vertex_weights_v2 || compute_vertex_weights_v3 || compute_vertex_weights_v4 || compute_vertex_weights_v5)
  {
  gROOT -> ProcessLine(".L compute_vertex_weights.C++");
  gROOT -> ProcessLine(".L plot_status_vertex_reweight.C++");
  }
  if (plot_status_vertex_reweight) { gROOT -> ProcessLine(".L plot_evolution_vertex_reweight.C++"); }
  if (get_pileup_normalization) { gROOT -> ProcessLine(".L get_pileup_normalization.C++"); }
  if (read_mc_ntuples_gen || read_mc_ntuples_det || read_data_ntuples) { gROOT -> ProcessLine(".L read_ntuple.C++"); }
  if (merge_data) { gROOT -> ProcessLine(".L merge_data.C++"); }
  if (normalize_mc) { gROOT -> ProcessLine(".L normalize_mc.C++"); }
  
  if (plot_control_dist)
  {
  gROOT -> ProcessLine(".L plot_mc_gen_level.C++");
  gROOT -> ProcessLine(".L control_plots.C++");
  }
  if (compute_corrections)   { gROOT -> ProcessLine(".L compute_correction.C++"); }
  if (apply_corrections) { gROOT -> ProcessLine(".L apply_correction.C++"); }
  if (compute_model_uncertainty) { gROOT -> ProcessLine(".L estimate_model_uncertainty.C++"); }
  if (compute_jes_uncertainty) { gROOT -> ProcessLine(".L estimate_jes_uncertainty.C++"); }
  if (compute_model_uncertainty || compute_jes_uncertainty) { gROOT -> ProcessLine(".L merge_uncertainties.C++"); }
  if (compute_total_uncertainty) { gROOT -> ProcessLine(".L total_uncertainty.C++"); }
  if (apply_uncertainty) { gROOT -> ProcessLine(".L apply_uncertainty.C++"); }
  if (convert_predictions) { gROOT -> ProcessLine(".L convert_prediction.C++"); }
  if (create_ratio) { gROOT -> ProcessLine(".L create_ratios.C++"); }
  if (do_final_plots) { gROOT -> ProcessLine(".L final_plots.C++"); }
  if (create_unfolding_response) { gROOT -> ProcessLine(".L create_unfolding_response.C++"); }
  if (check_response_matrix) { gROOT -> ProcessLine(".L check_response_matrix.C++"); }
  if (unfold) { gROOT -> ProcessLine(".L unfolding.C++"); }


//get the trigger efficiencies
if (get_trigger_efficiencies)
{
if (show_steps) { cout << "Get Trigger Efficiencies..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//get_triggered_events(jetmettaumon, out_jetmettaumon_none, n_files_jetmettaumon, "NONE", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U, n_files_jetmettaumon, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U, n_files_jetmettaumon, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet15U, n_files_jetmettaumon, "HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jetmettaumon, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jetmettaumon, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet30U, n_files_jetmettaumon, "HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jetmettaumon, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jetmettaumon, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U, n_files_jetmettaumon, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet50U, n_files_jetmettaumon, "HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jetmettaumon, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jetmettaumon, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U, n_files_jetmettaumon, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_DiJetAve15U, n_files_jetmettaumon, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_DiJetAve30U, n_files_jetmettaumon, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_FwdJet20U, n_files_jetmettaumon, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_HLT_DoubleJet15U_ForwardBackward, n_files_jetmettaumon, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jetmettaumon, out_jetmettaumon_all, n_files_jetmettaumon, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "JetMETTau"<<endl; }
//get_triggered_events(jetmettau, out_jetmettau_none, n_files_jetmettau, "NONE", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U, n_files_jetmettau, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U, n_files_jetmettau, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet15U, n_files_jetmettau, "HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jetmettau, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jetmettau, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet30U, n_files_jetmettau, "HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jetmettau, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jetmettau, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet15U_AND_HLT_Jet30U, n_files_jetmettau, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet50U, n_files_jetmettau, "HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jetmettau, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jetmettau, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_Jet15U_AND_HLT_Jet50U, n_files_jetmettau, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_DiJetAve15U, n_files_jetmettau, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_DiJetAve30U, n_files_jetmettau, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_FwdJet20U, n_files_jetmettau, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_HLT_DoubleJet15U_ForwardBackward, n_files_jetmettau, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jetmettau, out_jetmettau_all, n_files_jetmettau, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
//get_triggered_events(jetmet, out_jetmet_none, n_files_jetmet, "NONE", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U, n_files_jetmet, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U, n_files_jetmet, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_Jet15U, n_files_jetmet, "HLT_Jet15U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jetmet, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jetmet, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_Jet30U, n_files_jetmet, "HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jetmet, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jetmet, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_Jet15U_AND_HLT_Jet30U, n_files_jetmet, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_Jet50U, n_files_jetmet, "HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jetmet, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jetmet, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jetmet, out_jetmet_HLT_Jet15U_AND_HLT_Jet50U, n_files_jetmet, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_DiJetAve15U, n_files_jetmet, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_DiJetAve30U, n_files_jetmet, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_FwdJet20U, n_files_jetmet, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_HLT_DoubleJet15U_ForwardBackward, n_files_jetmet, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jetmet, out_jetmet_all, n_files_jetmet, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "Jet"<<endl; }
//get_triggered_events(jet, out_jet_none, n_files_jet, "NONE", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet6U, n_files_jet, "HLT_L1Jet6U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet10U, n_files_jet, "HLT_L1Jet10U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet15U, n_files_jet, "HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet6U_AND_HLT_Jet15U, n_files_jet, "HLT_L1Jet6U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet10U_AND_HLT_Jet15U, n_files_jet, "HLT_L1Jet10U_AND_HLT_Jet15U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet30U, n_files_jet, "HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_L1Jet6U_AND_HLT_Jet30U, n_files_jet, "HLT_L1Jet6U_AND_HLT_Jet30U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_L1Jet10U_AND_HLT_Jet30U, n_files_jet, "HLT_L1Jet10U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet15U_AND_HLT_Jet30U, n_files_jet, "HLT_Jet15U_AND_HLT_Jet30U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_Jet50U, n_files_jet, "HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_L1Jet6U_AND_HLT_Jet50U, n_files_jet, "HLT_L1Jet6U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_L1Jet10U_AND_HLT_Jet50U, n_files_jet, "HLT_L1Jet10U_AND_HLT_Jet50U", "allvertex", detail, test);
///get_triggered_events(jet, out_jet_HLT_Jet15U_AND_HLT_Jet50U, n_files_jet, "HLT_Jet15U_AND_HLT_Jet50U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_DiJetAve15U, n_files_jet, "HLT_DiJetAve15U_8E29", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_DiJetAve30U, n_files_jet, "HLT_DiJetAve30U_8E29", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_FwdJet20U, n_files_jet, "HLT_FwdJet20U", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_HLT_DoubleJet15U_ForwardBackward, n_files_jet, "HLT_DoubleJet15U_ForwardBackward", "allvertex", detail, test);
//get_triggered_events(jet, out_jet_all, n_files_jet, "ALL", "allvertex", detail, test);
if (show_steps) { cout << "All triggered were sucessfully extracted!"<<endl; }
}


//compute the trigger efficiencies
if (compute_trigger_efficiency)
{
if (show_steps) { cout << "Compute Trigger Efficiencies..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, eff_dir, "root");
create_directories(output_dir, trigger_plots, "plots");

if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_L1Jet6U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_L1Jet10U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet15U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_Jet15U, "JetMETTauMonitor_2010A_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jetmettaumon_HLT_L1Jet6U, false, eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet15U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jetmettaumon_HLT_L1Jet10U, false, eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet15U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet30U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jetmettaumon_HLT_L1Jet6U, false, eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jetmettaumon_HLT_L1Jet10U, false, eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U, true, out_jetmettaumon_HLT_Jet15U, false, eff_jetmettaumon_HLT_Jet15U_AND_HLT_Jet30U, "JetMETTauMonitor_2010A_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet50U, true, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jetmettaumon_HLT_L1Jet6U, false, eff_jetmettaumon_HLT_L1Jet6U_AND_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jetmettaumon_HLT_L1Jet10U, false, eff_jetmettaumon_HLT_L1Jet10U_AND_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U, true, out_jetmettaumon_HLT_Jet15U, false, eff_jetmettaumon_HLT_Jet15U_AND_HLT_Jet50U, "JetMETTauMonitor_2010A_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_DiJetAve15U, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_DiJetAve15U, "JetMETTauMonitor_2010A_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_DiJetAve30U, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_DiJetAve30U, "JetMETTauMonitor_2010A_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_FwdJet20U, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_FwdJet20U, "JetMETTauMonitor_2010A_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_HLT_DoubleJet15U_ForwardBackward, false, out_jetmettaumon_none, false, eff_jetmettaumon_HLT_DoubleJet15U_ForwardBackward, "JetMETTauMonitor_2010A_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettaumon_all, false, out_jetmettaumon_none, false, eff_jetmettaumon_all, "JetMETTauMonitor_2010A_all_", trigger_plots, detail);

if (show_steps) { cout << "JetMETTau"<<endl; }
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U, false, out_jetmettau_none, false, eff_jetmettau_HLT_L1Jet6U, "JetMETTau_2010A_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U, false, out_jetmettau_none, false, eff_jetmettau_HLT_L1Jet10U, "JetMETTau_2010A_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet15U, false, out_jetmettau_none, false, eff_jetmettau_HLT_Jet15U, "JetMETTau_2010A_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jetmettau_HLT_L1Jet6U, false, eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet15U, "JetMETTau_2010A_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jetmettau_HLT_L1Jet10U, false, eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet15U, "JetMETTau_2010A_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet30U, true, out_jetmettau_none, false, eff_jetmettau_HLT_Jet30U, "JetMETTau_2010A_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jetmettau_HLT_L1Jet6U, false, eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet30U, "JetMETTau_2010A_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jetmettau_HLT_L1Jet10U, false, eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet30U, "JetMETTau_2010A_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet15U_AND_HLT_Jet30U, true, out_jetmettau_HLT_Jet15U, false, eff_jetmettau_HLT_Jet15U_AND_HLT_Jet30U, "JetMETTau_2010A_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet50U, true, out_jetmettau_none, false, eff_jetmettau_HLT_Jet50U, "JetMETTau_2010A_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jetmettau_HLT_L1Jet6U, false, eff_jetmettau_HLT_L1Jet6U_AND_HLT_Jet50U, "JetMETTau_2010A_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jetmettau_HLT_L1Jet10U, false, eff_jetmettau_HLT_L1Jet10U_AND_HLT_Jet50U, "JetMETTau_2010A_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_Jet15U_AND_HLT_Jet50U, true, out_jetmettau_HLT_Jet15U, false, eff_jetmettau_HLT_Jet15U_AND_HLT_Jet50U, "JetMETTau_2010A_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_DiJetAve15U, false, out_jetmettau_none, false, eff_jetmettau_HLT_DiJetAve15U, "JetMETTau_2010A_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_DiJetAve30U, false, out_jetmettau_none, false, eff_jetmettau_HLT_DiJetAve30U, "JetMETTau_2010A_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_FwdJet20U, false, out_jetmettau_none, false, eff_jetmettau_HLT_FwdJet20U, "JetMETTau_2010A_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_HLT_DoubleJet15U_ForwardBackward, false, out_jetmettau_none, false, eff_jetmettau_HLT_DoubleJet15U_ForwardBackward, "JetMETTau_2010A_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmettau_all, false, out_jetmettau_none, false, eff_jetmettau_all, "JetMETTau_2010A_all_", trigger_plots, detail);

if (show_steps) { cout << "JetMET"<<endl; }
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U, false, out_jetmet_none, false, eff_jetmet_HLT_L1Jet6U, "JetMET_2010A_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U, false, out_jetmet_none, false, eff_jetmet_HLT_L1Jet10U, "JetMET_2010A_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet15U, false, out_jetmet_none, false, eff_jetmet_HLT_Jet15U, "JetMET_2010A_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jetmet_HLT_L1Jet6U, false, eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet15U, "JetMET_2010A_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jetmet_HLT_L1Jet10U, false, eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet15U, "JetMET_2010A_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet30U, true, out_jetmet_none, false, eff_jetmet_HLT_Jet30U, "JetMET_2010A_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jetmet_HLT_L1Jet6U, false, eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet30U, "JetMET_2010A_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jetmet_HLT_L1Jet10U, false, eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet30U, "JetMET_2010A_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet15U_AND_HLT_Jet30U, true, out_jetmet_HLT_Jet15U, false, eff_jetmet_HLT_Jet15U_AND_HLT_Jet30U, "JetMET_2010A_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet50U, true, out_jetmet_none, false, eff_jetmet_HLT_Jet50U, "JetMET_2010A_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jetmet_HLT_L1Jet6U, false, eff_jetmet_HLT_L1Jet6U_AND_HLT_Jet50U, "JetMET_2010A_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jetmet_HLT_L1Jet10U, false, eff_jetmet_HLT_L1Jet10U_AND_HLT_Jet50U, "JetMET_2010A_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_Jet15U_AND_HLT_Jet50U, true, out_jetmet_HLT_Jet15U, false, eff_jetmet_HLT_Jet15U_AND_HLT_Jet50U, "JetMET_2010A_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_DiJetAve15U, false, out_jetmet_none, false, eff_jetmet_HLT_DiJetAve15U, "JetMET_2010A_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_DiJetAve30U, false, out_jetmet_none, false, eff_jetmet_HLT_DiJetAve30U, "JetMET_2010A_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_FwdJet20U, false, out_jetmet_none, false, eff_jetmet_HLT_FwdJet20U, "JetMET_2010A_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_HLT_DoubleJet15U_ForwardBackward, false, out_jetmet_none, false, eff_jetmet_HLT_DoubleJet15U_ForwardBackward, "JetMET_2010A_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jetmet_all, false, out_jetmet_none, false, eff_jetmet_all, "JetMET_2010A_all_", trigger_plots, detail);

if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_efficiency(out_jet_HLT_L1Jet6U, false, out_jet_none, false, eff_jet_HLT_L1Jet6U, "Jet_2010B_HLT_L1Jet6U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U, false, out_jet_none, false, eff_jet_HLT_L1Jet10U, "Jet_2010B_HLT_L1Jet10U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet15U, false, out_jet_none, false, eff_jet_HLT_Jet15U, "Jet_2010B_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet6U_AND_HLT_Jet15U, true, out_jet_HLT_L1Jet6U, false, eff_jet_HLT_L1Jet6U_AND_HLT_Jet15U, "Jet_2010B_HLT_L1Jet6U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U_AND_HLT_Jet15U, true, out_jet_HLT_L1Jet10U, false, eff_jet_HLT_L1Jet10U_AND_HLT_Jet15U, "Jet_2010B_HLT_L1Jet10U_AND_HLT_Jet15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet30U, true, out_jet_none, false, eff_jet_HLT_Jet30U, "Jet_2010B_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet6U_AND_HLT_Jet30U, true, out_jet_HLT_L1Jet6U, false, eff_jet_HLT_L1Jet6U_AND_HLT_Jet30U, "Jet_2010B_HLT_L1Jet6U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U_AND_HLT_Jet30U, true, out_jet_HLT_L1Jet10U, false, eff_jet_HLT_L1Jet10U_AND_HLT_Jet30U, "Jet_2010B_HLT_L1Jet10U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet15U_AND_HLT_Jet30U, true, out_jet_HLT_Jet15U, false, eff_jet_HLT_Jet15U_AND_HLT_Jet30U, "Jet_2010B_HLT_Jet15U_AND_HLT_Jet30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet50U, true, out_jet_none, false, eff_jet_HLT_Jet50U, "Jet_2010B_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet6U_AND_HLT_Jet50U, true, out_jet_HLT_L1Jet6U, false, eff_jet_HLT_L1Jet6U_AND_HLT_Jet50U, "Jet_2010B_HLT_L1Jet6U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_L1Jet10U_AND_HLT_Jet50U, true, out_jet_HLT_L1Jet10U, false, eff_jet_HLT_L1Jet10U_AND_HLT_Jet50U, "Jet_2010B_HLT_L1Jet10U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_Jet15U_AND_HLT_Jet50U, true, out_jet_HLT_Jet15U, false, eff_jet_HLT_Jet15U_AND_HLT_Jet50U, "Jet_2010B_HLT_Jet15U_AND_HLT_Jet50U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_DiJetAve15U, false, out_jet_none, false, eff_jet_HLT_DiJetAve15U, "Jet_2010B_HLT_DiJetAve15U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_DiJetAve30U, false, out_jet_none, false, eff_jet_HLT_DiJetAve30U, "Jet_2010B_HLT_DiJetAve30U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_FwdJet20U, false, out_jet_none, false, eff_jet_HLT_FwdJet20U, "Jet_2010B_HLT_FwdJet20U_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_HLT_DoubleJet15U_ForwardBackward, false, out_jet_none, false, eff_jet_HLT_DoubleJet15U_ForwardBackward, "Jet_2010B_HLT_DoubleJet15U_ForwardBackward_", trigger_plots, detail);
compute_trigger_efficiency(out_jet_all, false, out_jet_none, false, eff_jet_HLT_DiJetAve30U, "Jet_2010B_all_", trigger_plots, detail);
if (show_steps) { cout << "All trigger efficiencies were sucessfully computed!"<<endl; }
}


//get the trigger turn on
if (compute_trigger_turn_on)
{
if (show_steps) { cout << "Get Trigger turn on..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
///compute_trigger_turn_on(jetmettaumon, out_jetmettaumon_emu, n_files_jetmettaumon, "", "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
///compute_trigger_turn_on(jetmettau, out_jetmettau_emu, n_files_jetmettau, "", "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
///compute_trigger_turn_on(jetmet, out_jetmet_emu, n_files_jetmet, "", "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_turn_on(jet, out_jet_emu, n_files_jet, "", "allvertex", detail, test);
}


//get the trigger correction v1
if (compute_trigger_correction_v1)
{
if (show_steps) { cout << "Compute Trigger Correction v1..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_corr_dir, "root");
create_directories(output_dir, trigger_corr_plots, "plots");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
/*
prefix = "JetMETTauMonitor_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jetmettaumon_emu, correction_jetmettaumon_emu_v1, trigger_corr_plots, prefix, "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jetmettaumon_emu, correction_sj_jetmettaumon_emu_v1, trigger_corr_plots, prefix, "_v1", true, detail, test);
*/
if (show_steps) { cout << "JetMETTau"<<endl; }
/*
prefix = "JetMETTau_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jetmettau_emu, correction_jetmettau_emu_v1, trigger_corr_plots, prefix, "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jetmettau_emu, correction_sj_jetmettau_emu_v1, trigger_corr_plots, prefix, "_v1", true, detail, test);
*/
if (show_steps) { cout << "JetMET"<<endl; }
/*
prefix = "JetMET_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jetmet_emu, correction_jetmet_emu_v1, trigger_corr_plots, prefix, "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jetmet_emu, correction_sj_jetmet_emu_v1, trigger_corr_plots, prefix, "_v1", true, detail, test);
*/
if (show_steps) { cout << "Jet"<<endl; }
/*
prefix = "Jet_2010B";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(out_jet_emu, correction_jet_emu_v1, trigger_corr_plots, "Jet_2010B", "_v1", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(out_jet_emu, correction_sj_jet_emu_v1, trigger_corr_plots, "Jet_2010B", "_v1", true, detail, test);
*/
}

//check the trigger correction v1
if (check_trigger_correction_v1)
{
if (show_steps) { cout << "Check Trigger Correction v1..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
compute_trigger_turn_on(jetmettaumon, check_jetmettaumon_emu_v1, n_files_jetmettaumon, correction_jetmettaumon_emu_v1, "allvertex", detail, test);
compute_trigger_turn_on(jetmettaumon, check_sj_jetmettaumon_emu_v1, n_files_jetmettaumon, correction_sj_jetmettaumon_emu_v1, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//compute_trigger_turn_on(jetmettau, check_jetmettau_emu_v1, n_files_jetmettau, correction_jetmettau_emu_v1, "allvertex", detail, test);
//compute_trigger_turn_on(jetmettau, check_sj_jetmettau_emu_v1, n_files_jetmettau, correction_sj_jetmettau_emu_v1, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//compute_trigger_turn_on(jetmet, check_jetmet_emu_v1, n_files_jetmet, correction_jetmet_emu_v1, "allvertex", detail, test);
//compute_trigger_turn_on(jetmet, check_sj_jetmet_emu_v1, n_files_jetmet, correction_sj_jetmet_emu_v1, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
//compute_trigger_turn_on(jet, check_jet_emu_v1, n_files_jet, correction_jet_emu_v1, "allvertex", detail, test);
//compute_trigger_turn_on(jet, check_sj_jet_emu_v1, n_files_jet, correction_sj_jet_emu_v1, "allvertex", detail, test);
}


//get the trigger correction v2
if (compute_trigger_correction_v2)
{
if (show_steps) { cout << "Compute Trigger Correction v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_corr_dir, "root");
create_directories(output_dir, trigger_corr_plots, "plots");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }

prefix = "JetMETTauMonitor_2010A";
create_directories(trigger_corr_plots + "png/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", prefix, "root");
create_directories(trigger_corr_plots + "c/", prefix, "root");
create_directories(trigger_corr_plots + "eps/", prefix, "root");
compute_trigger_correction(check_jetmettaumon_emu_v1, correction_jetmettaumon_emu_v2, trigger_corr_plots, prefix, "_v2", false, detail, test);
create_directories(trigger_corr_plots + "png/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "control_sj_" + prefix, "root");
create_directories(trigger_corr_plots + "png/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "c/", "sj_" + prefix, "root");
create_directories(trigger_corr_plots + "eps/", "sj_" + prefix, "root");
compute_trigger_correction(check_sj_jetmettaumon_emu_v1, correction_sj_jetmettaumon_emu_v2, trigger_corr_plots, prefix, "_v2", true, detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//compute_trigger_correction(check_jetmettau_emu_v1, correction_jetmettau_emu_v2, trigger_corr_plots, "JetMETTau_2010A", "_v2", false, detail, test);
//compute_trigger_correction(check_jetmettau_emu_v1, correction_sj_jetmettau_emu_v2, trigger_corr_plots, "JetMETTau_2010A", "_v2", true, detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
//compute_trigger_correction(check_jetmet_emu_v1, correction_jetmet_emu_v2, trigger_corr_plots, "JetMET_2010A", "_v2", false, detail, test);
//compute_trigger_correction(check_jetmet_emu_v1, correction_sj_jetmet_emu_v2, trigger_corr_plots, "JetMET_2010A", "_v2", true, detail, test);
if (show_steps) { cout << "Jet"<<endl; }
//compute_trigger_correction(check_jet_emu_v1, correction_jet_emu_v2, trigger_corr_plots, "Jet_2010B", "_v2", false, detail, test);
//compute_trigger_correction(check_jet_emu_v1, correction_sj_jet_emu_v2, trigger_corr_plots, "Jet_2010B", "_v2", true, detail, test);
}


//check the trigger correction v2
if (check_trigger_correction_v2)
{
if (show_steps) { cout << "Check Trigger Correction v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
///compute_trigger_turn_on(jetmettaumon, check_jetmettaumon_emu_v2, n_files_jetmettaumon, correction_jetmettaumon_emu_v2, "allvertex", detail, test);
///compute_trigger_turn_on(jetmettaumon, check_sj_jetmettaumon_emu_v2, n_files_jetmettaumon, correction_sj_jetmettaumon_emu_v2, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//compute_trigger_turn_on(jetmettau, check_jetmettau_emu_v2, n_files_jetmettau, correction_jetmettau_emu_v2, "allvertex", detail, test);
//compute_trigger_turn_on(jetmettau, check_sj_jetmettau_emu_v2, n_files_jetmettau, correction_sj_jetmettau_emu_v2, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//compute_trigger_turn_on(jetmet, check_jetmet_emu_v2, n_files_jetmet, correction_jetmet_emu_v2, "allvertex", detail, test);
//compute_trigger_turn_on(jetmet, check_sj_jetmet_emu_v2, n_files_jetmet, correction_sj_jetmet_emu_v2, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
//compute_trigger_turn_on(jet, check_jet_emu_v2, n_files_jet, correction_jet_emu_v2, "allvertex", detail, test);
//compute_trigger_turn_on(jet, check_sj_jet_emu_v2, n_files_jet, correction_sj_jet_emu_v2, "allvertex", detail, test);
}


//apply the trigger correction
if (apply_trigger_correction)
{
if (show_steps) { cout << "Apply Trigger Correction..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_dir, "root");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
///compute_trigger_turn_on(jetmettaumon, check_jetmettaumon_emu_v1, n_files_jetmettaumon, out_jetmettaumon_emu, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
///compute_trigger_turn_on(jetmettau, check_jetmettau_emu_v1, n_files_jetmettau, out_jetmettau_emu, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
///compute_trigger_turn_on(jetmet, check_jetmet_emu_v1, n_files_jetmet, out_jetmet_emu, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
compute_trigger_turn_on(jet, check_jet_emu_v1, n_files_jet, out_jet_emu, "allvertex", detail, test);
}


// plot final trigger efficiency
if (plot_final_trigger_efficiency)
{
if (show_steps) { cout << "Plot Final Trigger Efficiencies..."<<endl; }
create_directories(output_dir, trigger_plots, "plots");

plot_final_trigger_efficiency(out_jetmettau_emu, out_jetmet_emu, out_jet_emu, trigger_plots, detail, test);

if (show_steps) { cout << "Done!"<<endl; }
}


//counting combination statistics
if (counting_combination_statistics)
{
if (show_steps) { cout << "Counting Combination Statistics..."<<endl; }
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//counting_combination_statistics(jetmettaumon, "", n_files_jetmettaumon, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//counting_combination_statistics(jetmettau, "", n_files_jetmettau, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//counting_combination_statistics(jetmet, "", n_files_jetmet, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
counting_combination_statistics(jet, "", n_files_jet, "allvertex", detail, test);
}

//get the vertex distribution
if (get_vertex_distribution_v0)
{
if (show_steps) { cout << "Get the vertex distribution..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia6 4CTune"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "DATA"<<endl; }
if (show_steps) { cout << "JetMETTAU 2010A"<<endl; }
///get_vertex_distribution(jetmettau, vertex_jetmettau_allvertex, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
if (show_steps) { cout << "JetMET 2010A"<<endl; }
///get_vertex_distribution(jetmet, vertex_jetmet_allvertex, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Jet 2010B"<<endl; }
///get_vertex_distribution(jet, vertex_jet_allvertex, lumi_jet, n_files_jet, "DATA", "allvertex", "", "", detail, test);
if (show_steps) { cout << "The vertex distribution was sucessfully computed!"<<endl; }
}


// compute the vertex weights v1
if (compute_vertex_weights_v1)
{
if (show_steps) { cout << "Compute the vertex weights..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jetmettau_v1_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jetmet_v1_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B",vertex_weights_p6_z2_jet_v1, "_v1", "", "", vertex_weights_plots, "p6_z2_jet_v1_", false, detail);
plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v1, vertex_weights_p6_z2_jetmet_v1, vertex_weights_p6_z2_jet_v1, "_v1", vertex_weights_plots, "p6_z2_v1_", detail);
if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jetmettau_v1_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jetmet_v1_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_v1, "_v1", "", "", vertex_weights_plots, "p8_4c_jet_v1_", false, detail);
plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v1, vertex_weights_p8_4c_jetmet_v1, vertex_weights_p8_4c_jet_v1, "_v1", vertex_weights_plots, "p8_4c_v1_", detail);
if (show_steps) { cout << "All the vertex weights were sucessfully read!"<<endl; }
}


//get the vertex distribution v1
if (get_vertex_distribution_v1)
{
if (show_steps) { cout << "Get the vertex distribution v1..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v1, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v1, "_v1", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v1, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v1, "_v1", detail, test);
if (show_steps) { cout << "The vertex distribution v1 was sucessfully computed!"<<endl; }
}

// compute the vertex weights v2
if (compute_vertex_weights_v2)
{
if (show_steps) { cout << "Compute the vertex weights v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v1, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_v2, "_v2", vertex_weights_p6_z2_jetmettau_v1,  "_v1",vertex_weights_plots, "p6_z2_jetmettau_v2_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v1, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_v2, "_v2", vertex_weights_p6_z2_jetmet_v1, "_v1", vertex_weights_plots, "p6_z2_jetmet_v2_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jet_v1, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p6_z2_jet_v2, "_v2", vertex_weights_p6_z2_jet_v1, "_v1", vertex_weights_plots, "p6_z2_jet_v2_", false, detail);
plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v2, vertex_weights_p6_z2_jetmet_v2, vertex_weights_p6_z2_jet_v2, "_v2", vertex_weights_plots, "p6_z2_v2_", detail);
if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v1, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_v2, "_v2", vertex_weights_p8_4c_jetmettau_v1, "_v1", vertex_weights_plots, "p8_4c_jetmettau_v2_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v1, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_v2, "_v2", vertex_weights_p8_4c_jetmet_v1, "_v1", vertex_weights_plots, "p8_4c_jetmet_v2_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jet_v1, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_p8_4c_jet_v1, "_v1", vertex_weights_plots, "p8_4c_jet_v2_", false, detail);
plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v2, vertex_weights_p8_4c_jetmet_v2, vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_plots, "p8_4c_v2_", detail);
if (show_steps) { cout << "All the vertex weights v2 were sucessfully read!"<<endl; }
}


//get the vertex distribution v2
if (get_vertex_distribution_v2)
{
if (show_steps) { cout << "Get the vertex distribution v2..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v2, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v2, "_v2", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v2, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v2, "_v2", detail, test);
if (show_steps) { cout << "The vertex distribution v2 was sucessfully computed!"<<endl; }
}


// compute the vertex weights
if (compute_vertex_weights_v3)
{
if (show_steps) { cout << "Compute the vertex weights v3..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v2, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_v3, "_v3", vertex_weights_p6_z2_jetmettau_v2,  "_v2",vertex_weights_plots, "p6_z2_jetmettau_v3_", true, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v2, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_v3, "_v3", vertex_weights_p6_z2_jetmet_v2, "_v2", vertex_weights_plots, "p6_z2_jetmet_v3_", true, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jet_v2, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p6_z2_jet_v3, "_v3", vertex_weights_p6_z2_jet_v2, "_v2", vertex_weights_plots, "p6_z2_jet_v3_", true, detail);
plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v3, vertex_weights_p6_z2_jetmet_v3, vertex_weights_p6_z2_jet_v3, "_v3", vertex_weights_plots, "p6_z2_v3_", detail);
if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v2, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_v3, "_v3", vertex_weights_p8_4c_jetmettau_v2, "_v2", vertex_weights_plots, "p8_4c_jetmettau_v3_", true, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v2, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_v3, "_v3", vertex_weights_p8_4c_jetmet_v2, "_v2", vertex_weights_plots, "p8_4c_jetmet_v3_", true, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jet_v2, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_p8_4c_jet_v2, "_v2", vertex_weights_plots, "p8_4c_jet_v3_", true, detail);
plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v3, vertex_weights_p8_4c_jetmet_v3, vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_plots, "p8_4c_v3_", detail);
if (show_steps) { cout << "All the new vertex weights v3 were sucessfully read!"<<endl; }
}

//get the vertex distribution v3
if (get_vertex_distribution_v3)
{
if (show_steps) { cout << "Get the vertex distribution v3..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v3, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v3, "_v3", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
///get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v3, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v3, "_v3", detail, test);
if (show_steps) { cout << "The vertex distribution v3 was sucessfully computed!"<<endl; }
}


// compute the vertex weights v4
if (compute_vertex_weights_v4)
{
if (show_steps) { cout << "Compute the new vertex weights v4..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v3, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_v4, "_v4", vertex_weights_p6_z2_jetmettau_v3,  "_v3",vertex_weights_plots, "p6_z2_jetmettau_v4_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v3, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_v4, "_v4", vertex_weights_p6_z2_jetmet_v3, "_v3", vertex_weights_plots, "p6_z2_jetmet_v4_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jet_v3, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p6_z2_jet_v4, "_v4", vertex_weights_p6_z2_jet_v3, "_v3", vertex_weights_plots, "p6_z2_jet_v4_", false, detail);
plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v4, vertex_weights_p6_z2_jetmet_v4, vertex_weights_p6_z2_jet_v4, "_v4", vertex_weights_plots, "p6_z2_v4_", detail);
if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v3, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_v4, "_v4", vertex_weights_p8_4c_jetmettau_v3, "_v3", vertex_weights_plots, "p8_4c_jetmettau_v4_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v3, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_v4, "_v4", vertex_weights_p8_4c_jetmet_v3, "_v3", vertex_weights_plots, "p8_4c_jetmet_v4_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jet_v3, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_p8_4c_jet_v3, "_v3", vertex_weights_plots, "p8_4c_jet_v4_", false, detail);
plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v4, vertex_weights_p8_4c_jetmet_v4, vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_plots, "p8_4c_v4_", detail);
if (show_steps) { cout << "All the new vertex weights v4 were sucessfully read!"<<endl; }
}


//get the vertex distribution v4
if (get_vertex_distribution_v4)
{
if (show_steps) { cout << "Get the vertex distribution v4..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_dir, "root");
if (show_steps) { cout << "MC on Detector Level Reweighted to match data"<<endl; }
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmettau_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jetmet_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
get_vertex_distribution(mc_p6_z2, vertex_p6_z2_allvertex_jet_v4, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmettau_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jetmet_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v4, "_v4", detail, test);
if (show_steps) { cout << "Pythia6 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//get_vertex_distribution(mc_p8_4c, vertex_p8_4c_allvertex_jet_v4, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v4, "_v4", detail, test);
if (show_steps) { cout << "The vertex distribution v4 was sucessfully computed!"<<endl; }
}

// compute the vertex weights v5
if (compute_vertex_weights_v5)
{
if (show_steps) { cout << "Compute the new vertex weights v5..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p6_z2_allvertex_jetmettau_v4, "Pythia6 Z2starTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p6_z2_jetmettau_v5, "_v5", vertex_weights_p6_z2_jetmettau_v4,  "_v4",vertex_weights_plots, "p6_z2_jetmettau_v5_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jetmet_v4, "Pythia6 Z2starTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p6_z2_jetmet_v5, "_v5", vertex_weights_p6_z2_jetmet_v4, "_v4", vertex_weights_plots, "p6_z2_jetmet_v5_", false, detail);
compute_vertex_weights(vertex_p6_z2_allvertex_jet_v4, "Pythia6 Z2starTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p6_z2_jet_v5, "_v5", vertex_weights_p6_z2_jet_v4, "_v4", vertex_weights_plots, "p6_z2_jet_v5_", false, detail);
plot_status_vertex_reweight(vertex_weights_p6_z2_jetmettau_v5, vertex_weights_p6_z2_jetmet_v5, vertex_weights_p6_z2_jet_v5, "_v5", vertex_weights_plots, "p6_z2_v5_", detail);
if (show_steps) { cout << "Pythia8 4CTune reweighted to data"<<endl; }
compute_vertex_weights(vertex_p8_4c_allvertex_jetmettau_v4, "Pythia8 4CTune", vertex_jetmettau_allvertex, "JetMETTau 2010A", vertex_weights_p8_4c_jetmettau_v5, "_v5", vertex_weights_p8_4c_jetmettau_v4, "_v4", vertex_weights_plots, "p8_4c_jetmettau_v5_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jetmet_v4, "Pythia8 4CTune", vertex_jetmet_allvertex, "JetMET 2010A", vertex_weights_p8_4c_jetmet_v5, "_v5", vertex_weights_p8_4c_jetmet_v4, "_v4", vertex_weights_plots, "p8_4c_jetmet_v5_", false, detail);
compute_vertex_weights(vertex_p8_4c_allvertex_jet_v4, "Pythia8 4CTune", vertex_jet_allvertex, "Jet 2010B", vertex_weights_p8_4c_jet_v5, "_v5", vertex_weights_p8_4c_jet_v4, "_v4", vertex_weights_plots, "p8_4c_jet_v5_", false, detail);
plot_status_vertex_reweight(vertex_weights_p8_4c_jetmettau_v5, vertex_weights_p8_4c_jetmet_v5, vertex_weights_p8_4c_jet_v5, "_v5", vertex_weights_plots, "p8_4c_v5_", detail);
if (show_steps) { cout << "All the new vertex weights v5 were sucessfully read!"<<endl; }
}


//plots the status of the vertex reweighting process
if (plot_status_vertex_reweight)
{
if (show_steps) { cout << "Plot the Status of the Vertex Reweight..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, vertex_weights_dir, "root");
create_directories(output_dir, vertex_weights_plots, "plots");
if (show_steps) { cout << "Pythia6 Z2starTune reweighting evolution..."<<endl; }
plot_evolution_vertex_reweight(vertex_weights_p6_z2_jetmettau_v1, vertex_weights_p6_z2_jetmettau_v2, vertex_weights_p6_z2_jetmettau_v3, vertex_weights_p6_z2_jetmettau_v4, vertex_weights_p6_z2_jetmettau_v5, vertex_weights_plots, "p6_z2_jetmettau_", detail);
plot_evolution_vertex_reweight(vertex_weights_p6_z2_jetmet_v1, vertex_weights_p6_z2_jetmet_v2, vertex_weights_p6_z2_jetmet_v3, vertex_weights_p6_z2_jetmet_v4, vertex_weights_p6_z2_jetmet_v5, vertex_weights_plots, "p6_z2_jetmet_", detail);
plot_evolution_vertex_reweight(vertex_weights_p6_z2_jet_v1, vertex_weights_p6_z2_jet_v2, vertex_weights_p6_z2_jet_v3, vertex_weights_p6_z2_jet_v4, vertex_weights_p6_z2_jet_v5, vertex_weights_plots, "p6_z2_jet_", detail);
if (show_steps) { cout << "Pythia8 4CTune reweighting evolution..."<<endl; }
plot_evolution_vertex_reweight(vertex_weights_p8_4c_jetmettau_v1, vertex_weights_p8_4c_jetmettau_v2, vertex_weights_p8_4c_jetmettau_v3, vertex_weights_p8_4c_jetmettau_v4, vertex_weights_p8_4c_jetmettau_v5, vertex_weights_plots, "p8_4c_jetmettau_", detail);
plot_evolution_vertex_reweight(vertex_weights_p8_4c_jetmet_v1, vertex_weights_p8_4c_jetmet_v2, vertex_weights_p8_4c_jetmet_v3, vertex_weights_p8_4c_jetmet_v4, vertex_weights_p8_4c_jetmet_v5, vertex_weights_plots, "p8_4c_jetmet_",detail);
plot_evolution_vertex_reweight(vertex_weights_p8_4c_jet_v1, vertex_weights_p8_4c_jet_v2, vertex_weights_p8_4c_jet_v3, vertex_weights_p8_4c_jet_v4, vertex_weights_p8_4c_jet_v5, vertex_weights_plots, "p8_4c_jet_", detail);
if (show_steps) { cout << "The status was sucessfully plotted!"<<endl; }
}


if (get_pileup_normalization)
{
if (show_steps) { cout << "Getting the Pileup Normalization..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, pileup_norm_dir, "root");
if (show_steps) { cout << "For Pythia 6 - Tune Z2Star..."<<endl; }
//get_pileup_normalization(mc_p6_z2, norm_p6_z2, lumi_p6_z2, n_files_p6_z2, "", "", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_jetmettau, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p6_z2, norm_p6_z2_jetmet, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
get_pileup_normalization(mc_p6_z2, norm_p6_z2_jet, lumi_p6_z2, n_files_p6_z2, vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "For Pythia 8 - Tune 4C..."<<endl; }
//get_pileup_normalization(mc_p8_4c, norm_p8_4c, lumi_p8_4c, n_files_p8_4c, "", "", detail, test);
///get_pileup_normalization(mc_p8_4c, norm_p8_4c_jetmettau, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p8_4c, norm_p8_4c_jetmet, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
//get_pileup_normalization(mc_p8_4c, norm_p8_4c_jet, lumi_p8_4c, n_files_p8_4c, vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "The Pileup Normalization was been sucessfully obtained!"<<endl; }
if (show_steps) { cout << "Showing Pileup Normalization ..."<<endl; }
show_pileup_normalization(norm_p6_z2, norm_p6_z2_jetmettau, norm_p6_z2_jetmet, norm_p6_z2_jet, norm_p8_4c, norm_p8_4c_jetmettau, norm_p8_4c_jetmet, norm_p8_4c_jet, detail);
if (show_steps) { cout << "The Pileup Normalization was been sucessfully outputed!"<<endl; }
}

//read the MC Ntuples on the Generator Level
if (read_mc_ntuples_gen)
{
if (show_steps) { cout << "Reading MC Ntuples on Generator Level..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_gen_dir, "root");
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", "", "", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_gen_nopileup_jet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "nopileup", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_gen_allvertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_GEN", "allvertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", "", "", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test); //bus error?
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_gen_nopileup_jet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "nopileup", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_gen_allvertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_GEN", "allvertex", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "All the MC Ntuples were sucessfully read!"<<endl; }
}


//read the MC Ntuples on the Detector Level
if (read_mc_ntuples_det)
{
if (show_steps) { cout << "Reading MC Ntuples on Detector Level..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_det_dir, "root");
if (show_steps) { cout << "Pythia6 Z2starTune"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", "", "", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmettau, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmettau_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jetmet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jetmet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia6 Z2starTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p6_z2, out_p6_z2_det_1vertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "1vertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p6_z2, out_p6_z2_det_allvertex_jet, lumi_p6_z2, n_files_p6_z2, "MC_DET", "allvertex", vertex_weights_p6_z2_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex", "", "", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", "", "", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMETTau_2010A"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex",  vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmettau, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmettau_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match JetMET_2010A"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex",  vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jetmet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jetmet_v5, "_v5", detail, test);
if (show_steps) { cout << "Pythia8 4CTune with Vertex Reweight to match Jet_2010B"<<endl; }
//read_ntuple(mc_p8_4c, out_p8_4c_det_1vertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "1vertex",  vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
//read_ntuple(mc_p8_4c, out_p8_4c_det_allvertex_jet, lumi_p8_4c, n_files_p8_4c, "MC_DET", "allvertex", vertex_weights_p8_4c_jet_v5, "_v5", detail, test);
if (show_steps) { cout << "All the MC Ntuples were sucessfully read!"<<endl; }
}


//read the data Ntuples
if (read_data_ntuples)
{
if (show_steps) { cout << "Reading Data Ntuples..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, raw_data_dir, "root");
if (show_steps) { cout << "JetMETTau"<<endl; }
//read_ntuple(jetmettau, out_jetmettau_1vertex, lumi_jetmettau, n_files_jetmettau, "DATA", "1vertex", "", "", detail, test);
//read_ntuple(jetmettau, out_jetmettau_allvertex, lumi_jetmettau, n_files_jetmettau, "DATA", "allvertex", "", "", detail, test);
//read_ntuple(jetmettau, out_jetmettau_up, lumi_jetmettau, n_files_jetmettau, "DATA", "up", "", "", detail, test);
//read_ntuple(jetmettau, out_jetmettau_down, lumi_jetmettau, n_files_jetmettau, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "JetMET"<<endl; }
//read_ntuple(jetmet, out_jetmet_1vertex, lumi_jetmet, n_files_jetmet, "DATA", "1vertex", "", "", detail, test);
//read_ntuple(jetmet, out_jetmet_allvertex, lumi_jetmet, n_files_jetmet, "DATA", "allvertex", "", "", detail, test);
//read_ntuple(jetmet, out_jetmet_up, lumi_jetmet, n_files_jetmet, "DATA", "up", "", "", detail, test);
read_ntuple(jetmet, out_jetmet_down, lumi_jetmet, n_files_jetmet, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "Jet"<<endl; }
//read_ntuple(jet, out_jet_1vertex, lumi_jet, n_files_jet, "DATA", "1vertex", "", "", detail, test);
//read_ntuple(jet, out_jet_allvertex, lumi_jet, n_files_jet, "DATA", "allvertex", "", "", detail, test);
//read_ntuple(jet, out_jet_up, lumi_jet, n_files_jet, "DATA", "up", "", "", detail, test);
//read_ntuple(jet, out_jet_down, lumi_jet, n_files_jet, "DATA", "down", "", "", detail, test);
if (show_steps) { cout << "All the Data Ntuples were sucessfully read!"<<endl; }
}


//normalize the different mc dataset
if (normalize_mc)
{
if (show_steps) { cout << "Normalize MC..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_norm_dir, "root");
create_directories(output_dir, mc_norm_plots, "plots");
normalize_mc(out_p6_z2_gen_allvertex, out_p6_z2_gen_nopileup, mc_norm_plots, out_p6_z2_gen_norm , "pythia6_z2_", "ak5Gen_", detail,show_results);
normalize_mc(out_p8_4c_gen_allvertex, out_p8_4c_gen_nopileup, mc_norm_plots, out_p8_4c_gen_norm , "pythia8_4c_", "ak5Gen_", detail,show_results);
normalize_mc(out_p6_z2_gen_allvertex_jetmettau, out_p6_z2_gen_nopileup_jetmettau, mc_norm_plots, out_p6_z2_gen_jetmettau_norm , "pythia6_z2_JetMETTau_2010A_", "ak5Gen_", detail,show_results);
normalize_mc(out_p8_4c_gen_allvertex_jetmettau, out_p8_4c_gen_nopileup_jetmettau, mc_norm_plots, out_p8_4c_gen_jetmettau_norm , "pythia8_4c_JetMETTau_2010A", "ak5Gen_", detail,show_results);
normalize_mc(out_p6_z2_gen_allvertex_jetmet, out_p6_z2_gen_nopileup_jetmet, mc_norm_plots, out_p6_z2_gen_jetmet_norm , "pythia6_z2_JetMET_2010A_", "ak5Gen_", detail,show_results);
normalize_mc(out_p8_4c_gen_allvertex_jetmet, out_p8_4c_gen_nopileup_jetmet, mc_norm_plots, out_p8_4c_gen_jetmet_norm , "pythia8_4c_JetMET_2010A", "ak5Gen_", detail,show_results);
normalize_mc(out_p6_z2_gen_allvertex_jet, out_p6_z2_gen_nopileup_jetmet, mc_norm_plots, out_p6_z2_gen_jet_norm , "pythia6_z2_Jet_2010B_", "ak5Gen_", detail,show_results);
normalize_mc(out_p8_4c_gen_allvertex_jet, out_p8_4c_gen_nopileup_jetmet, mc_norm_plots, out_p8_4c_gen_jet_norm , "pythia8_4c_Jet_2010B_", "ak5Gen_", detail,show_results);
normalize_mc(out_p6_z2_det_allvertex, out_p6_z2_det_1vertex, mc_norm_plots, out_p6_z2_det_norm , "pythia6_z2_", "ak5PF_", detail,show_results);
normalize_mc(out_p8_4c_det_allvertex, out_p8_4c_det_1vertex, mc_norm_plots, out_p8_4c_det_norm , "pythia8_4c_", "ak5PF_", detail,show_results);
normalize_mc(out_p6_z2_det_allvertex_jetmettau, out_p6_z2_det_1vertex_jetmettau, mc_norm_plots, out_p6_z2_det_jetmettau_norm , "pythia6_z2_JetMETTau_2010A_", "ak5PF_", detail,show_results);
normalize_mc(out_p8_4c_det_allvertex_jetmettau, out_p8_4c_det_1vertex_jetmettau, mc_norm_plots, out_p8_4c_det_jetmettau_norm , "pythia8_4c_JetMETTau_2010A_", "ak5PF_", detail,show_results);
normalize_mc(out_p6_z2_det_allvertex_jetmet, out_p6_z2_det_1vertex_jetmet, mc_norm_plots, out_p6_z2_det_jetmet_norm , "pythia6_z2_JetMET_2010A_", "ak5PF_", detail,show_results);
normalize_mc(out_p8_4c_det_allvertex_jetmet, out_p8_4c_det_1vertex_jetmet, mc_norm_plots, out_p8_4c_det_jetmet_norm , "pythia8_4c_JetMET_2010A_", "ak5PF_", detail,show_results);
normalize_mc(out_p6_z2_det_allvertex_jet, out_p6_z2_det_1vertex_jet, mc_norm_plots, out_p6_z2_det_jet_norm , "pythia6_z2_Jet_2010B_", "ak5PF_", detail,show_results);
normalize_mc(out_p8_4c_det_allvertex_jet, out_p8_4c_det_1vertex_jet, mc_norm_plots, out_p8_4c_det_jet_norm , "pythia8_4c_Jet_2010B_", "ak5PF_", detail,show_results);
if (show_steps) { cout << "The MC datasets were sucessfully normalized!"<<endl; }
}

//compute correction
if (compute_corrections)
{
if (show_steps) { cout << "Computing Corrections..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, corrections_dir, "root");
create_directories(output_dir, corrections_plots, "plots");

if (show_steps) { cout << "Computing the Pileup Correction..."<<endl; }
compute_correction(out_p6_z2_gen_nopileup, "ak5Gen_", "Pythia 6 Tune Z2* - Generator Level with no pileup", out_p6_z2_gen_allvertex, "ak5Gen_", "Pythia 6 Tune Z2* - Generator Level with all vertex", corrections_plots, corr_pileup_p6_z2, "pileup_pythia6_z2_", "Pileup Correction with Pythia 6 - Tune Z2*", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup, "ak5Gen_", "Pythia 8 Tune 4C - Generator Level with no pileup", out_p8_4c_gen_allvertex, "ak5Gen_", "Pythia 8 Tune 4C - Generator Level with all vertex", corrections_plots, corr_pileup_p8_4c, "pileup_pythia8_4c_", "Pileup Correction with Pythia 8 - Tune 4C", false, detail, show_results);
compute_correction(out_p6_z2_gen_nopileup_jetmettau, "ak5Gen_", "Pythia 6 Tune Z2* -  Gen. Level no pileup - JetMETTau_2010A", out_p6_z2_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level all vertex - JetMETTau_2010A", corrections_plots, corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pileup Correction with Pythia 6 - Tune Z2* - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup_jetmettau, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level no pileup - JetMETTau_2010A", out_p8_4c_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level all vertex - JetMETTau_2010A", corrections_plots, corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pileup Correction with Pythia 8 - Tune 4C - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_p6_z2_gen_nopileup_jetmet, "ak5Gen_", "Pythia 6 Tune Z2* -  Gen. Level no pileup JetMET_2010A", out_p6_z2_gen_allvertex_jetmet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level all vertex - JetMET_2010A", corrections_plots, corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pileup Correction with Pythia 6 - Tune Z2* - JetMET_2010A", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup_jetmet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level no pileup JetMET_2010A", out_p8_4c_gen_allvertex_jetmet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level all vertex - JetMET_2010A", corrections_plots, corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pileup Correction with Pythia 8 - Tune 4C - JetMET_2010A", false, detail, show_results);
compute_correction(out_p6_z2_gen_nopileup_jet, "ak5Gen_", "Pythia 6 Tune Z2* -  Gen. Level no pileup Jet_2010B", out_p6_z2_gen_allvertex_jet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level all vertex - Jet_2010B", corrections_plots, corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pileup Correction with Pythia 6 - Tune Z2* - Jet_2010B", false, detail, show_results);
compute_correction(out_p8_4c_gen_nopileup_jet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level no pileup Jet_2010B", out_p8_4c_gen_allvertex_jet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level all vertex - Jet_2010B", corrections_plots, corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pileup Correction with Pythia 8 - Tune 4C - Jet_2010B", false, detail, show_results);

if (show_steps) { cout << "Computing the Detector Correction..."<<endl; }
compute_correction(out_p6_z2_det_allvertex, "ak5PF_", "Pythia 6 Tune Z2* - Detector Level", out_p6_z2_gen_allvertex, "ak5Gen_", "Pythia 6 Tune Z2* - Generator Level", corrections_plots, corr_detector_p6_z2, "detector_pythia6_z2_", "Detector Correction with Pythia 6 - Tune Z2*", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex, "ak5PF_", "Pythia 8  - Tune 4C - Detector Level", out_p8_4c_gen_allvertex, "ak5Gen_", "Pythia 8 Tune 4C - Generator Level", corrections_plots, corr_detector_p8_4c, "detector_pythia8_4c_", "Detector Correction with Pythia 8 - Tune 4C", false, detail, show_results);
compute_correction(out_p6_z2_det_allvertex_jetmettau, "ak5PF_", "Pythia 6 Tune Z2* - Det. Level JetMETTau_2010A", out_p6_z2_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level JetMETTau_2010A", corrections_plots, corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Detector Correction with Pythia 6 - Tune Z2* - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex_jetmettau, "ak5PF_", "Pythia 8  - Tune 4C - Det. Level JetMETTau_2010A", out_p8_4c_gen_allvertex_jetmettau, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level JetMETTau_2010A", corrections_plots, corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Detector Correction with Pythia 8 - Tune 4C - JetMETTau_2010A", false, detail, show_results);
compute_correction(out_p6_z2_det_allvertex_jetmet, "ak5PF_", "Pythia 6 Tune Z2* - Det. Level JetMET_2010A", out_p6_z2_gen_allvertex_jetmet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level JetMET_2010A", corrections_plots, corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Detector Correction with Pythia 6 - Tune Z2* - JetMET_2010A", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex_jetmet, "ak5PF_", "Pythia 8  - Tune 4C - Det. Level JetMET_2010A", out_p8_4c_gen_allvertex_jetmet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level JetMET_2010A", corrections_plots, corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Detector Correction with Pythia 8 - Tune 4C - JetMET_2010A", false, detail, show_results);
compute_correction(out_p6_z2_det_allvertex_jet, "ak5PF_", "Pythia 6 Tune Z2* - Det. Level Jet_2010B", out_p6_z2_gen_allvertex_jet, "ak5Gen_", "Pythia 6 Tune Z2* - Gen. Level Jet_2010B", corrections_plots, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Detector Correction with Pythia 6 - Tune Z2* - Jet_2010B", false, detail, show_results);
compute_correction(out_p8_4c_det_allvertex_jet, "ak5PF_", "Pythia 8  - Tune 4C - Det. Level Jet_2010B", out_p8_4c_gen_allvertex_jet, "ak5Gen_", "Pythia 8 Tune 4C - Gen. Level Jet_2010B", corrections_plots, corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Detector Correction with Pythia 8 - Tune 4C - Jet_2010B", false, detail, show_results);

if (show_steps) { cout << "Computing the Overall Correction..."<<endl; }
compute_correction(corr_detector_p6_z2, "detector_pythia6_z2_", "Pythia 6 Tune Z2* - Detector Correction", corr_pileup_p6_z2, "pileup_pythia6_z2_", "Pythia 6 Tune Z2* - Pileup Correction", corrections_plots, corr_final_p6_z2, "final_pythia6_z2_", "Final Correction with Pythia 6 - Tune Z2*", true, detail, show_results);
compute_correction(corr_detector_p8_4c, "detector_pythia8_4c_", "Pythia 8  - Tune 4C - Detector Correction", corr_pileup_p8_4c, "pileup_pythia8_4c_", "Pythia 8 Tune 4C - Pileup Correction", corrections_plots, corr_final_p8_4c, "final_pythia8_4c_", "Final Correction with Pythia 8 - Tune 4C", true, detail, show_results);
compute_correction(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Pythia 6 Tune Z2* - Detector Correction - JetMETTau_2010A", corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pythia 6 Tune Z2* - Pileup Correction - JetMETTau_2010A", corrections_plots, corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "Final Correction with Pythia 6 - Tune Z2* - JetMETTau_2010A", true, detail, show_results);
compute_correction(corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Pythia 8 Tune 4C - Detector Correction - JetMETTau_2010A", corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pythia 8 Tune 4C - Pileup Correction - JetMETTau_2010A", corrections_plots, corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "Final Correction with Pythia 8 - Tune 4C - JetMETTau_2010A", true, detail, show_results);
compute_correction(corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Pythia 6 Tune Z2* - Detector Correction - JetMET_2010A", corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pythia 6 Tune Z2* - Pileup Correction - JetMET_2010A", corrections_plots, corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "Final Correction with Pythia 6 - Tune Z2* - JetMET_2010A", true, detail, show_results);
compute_correction(corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Pythia 8 Tune 4C - Detector Correction - JetMET_2010A", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pythia 8 Tune 4C - Pileup Correction - JetMET_2010A", corrections_plots, corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "Final Correction with Pythia 8 - Tune 4C - JetMET_2010A", true, detail, show_results);
compute_correction(corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Pythia 6 Tune Z2* - Detector Correction - Jet_2010B", corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pythia 6 Tune Z2* - Pileup Correction - Jet_2010B", corrections_plots, corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Final Correction with Pythia 6 - Tune Z2* - Jet_2010B", true, detail, show_results);
compute_correction(corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Pythia 8 Tune 4C - Detector Correction - Jet_2010B", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pythia 8 Tune 4C - Pileup Correction - Jet_2010B", corrections_plots, corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Final Correction with Pythia 8 - Tune 4C - Jet_2010B", true, detail, show_results);

if (show_steps) { cout << "Ploting Final Results..."<<endl; }
display_final_corrections(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Detector Correction", corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pileup Correction", corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "Final Correction", corrections_plots, "corrections_pythia6_z2_JetMETTau_2010A_", detail);
display_final_corrections(corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Detector Correction", corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pileup Correction", corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "Final Correction", corrections_plots, "corrections_pythia6_z2_JetMET_2010A_", detail);
display_final_corrections(corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Detector Correction", corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pileup Correction", corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Final Correction", corrections_plots, "corrections_pythia6_z2_Jet_2010B_", detail);

display_final_corrections(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "JetMETTau_2010A", corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "JetMET_2010A", corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia6_z2_detector_", detail);
display_final_corrections(corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "JetMETTau_2010A", corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "JetMET_2010A", corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia6_z2_pileup_", detail);
display_final_corrections(corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "JetMETTau_2010A", corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "JetMET_2010A", corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia6_z2_final_", detail);

display_final_corrections(corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Detector Correction", corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pileup Correction", corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "Final Correction", corrections_plots, "corrections_JetMETTau_2010A_", detail);
display_final_corrections(corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Detector Correction", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pileup Correction", corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "Final Correction", corrections_plots, "corrections_pythia8_4c_JetMET_2010A_", detail);
display_final_corrections(corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Detector Correction", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pileup Correction", corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Final Correction", corrections_plots, "corrections_pythia8_4c_Jet_2010B_", detail);

display_final_corrections(corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "JetMETTau_2010A", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "JetMET_2010A", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia8_4c_detector_", detail);
display_final_corrections(corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMETTau_2010A_", "JetMETTau_2010A", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "JetMET_2010A", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia8_4c_pileup_", detail);
display_final_corrections(corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "JetMETTau_2010A", corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "JetMET_2010A", corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Jet_2010B", corrections_plots, "corrections_pythia8_4c_final_", detail);

display_final_corrections(corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", "Pythia6 - Z2*", corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_detector_JetMETTau_2010A_", detail);
display_final_corrections(corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", "Pythia6 - Z2*", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_detector_JetMET_2010_", detail);
display_final_corrections(corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", "Pythia6 - Z2*", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_detector_Jet_2010B_", detail);

display_final_corrections(corr_pileup_p6_z2_jetmettau, "pileup_pythia6_z2_JetMETTau_2010A_", "Pythia6 - Z2*", corr_pileup_p8_4c_jetmettau, "pileup_pythia8_4c_JetMETTau_2010A_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_pileup_JetMETTau_2010A_", detail);
display_final_corrections(corr_pileup_p6_z2_jetmet, "pileup_pythia6_z2_JetMET_2010A_", "Pythia6 - Z2*", corr_pileup_p8_4c_jetmet, "pileup_pythia8_4c_JetMET_2010A_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_pileup_JetMET_2010_", detail);
display_final_corrections(corr_pileup_p6_z2_jet, "pileup_pythia6_z2_Jet_2010B_", "Pythia6 - Z2*", corr_pileup_p8_4c_jet, "pileup_pythia8_4c_Jet_2010B_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_pileup_Jet_2010B_", detail);

display_final_corrections(corr_final_p6_z2_jetmettau, "final_pythia6_z2_JetMETTau_2010A_", "Pythia6 - Z2*", corr_final_p8_4c_jetmettau, "final_pythia8_4c_JetMETTau_2010A_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_final_JetMETTau_2010A_", detail);
display_final_corrections(corr_final_p6_z2_jetmet, "final_pythia6_z2_JetMET_2010A_", "Pythia6 - Z2*", corr_final_p8_4c_jetmet, "final_pythia8_4c_JetMET_2010A_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_final_JetMET_2010_", detail);
display_final_corrections(corr_final_p6_z2_jet, "final_pythia6_z2_Jet_2010B_", "Pythia6 - Z2*", corr_final_p8_4c_jet, "final_pythia8_4c_Jet_2010B_", "Pythia8 - 4C", "", "", "", corrections_plots, "corrections_final_Jet_2010B_", detail);
if (show_steps) { cout << "The correction factors were sucessfully computed!"<<endl; }
}


//apply correction
if (apply_corrections)
{Anton Stark E que tal o facto de ser um órgão visualmente horrível e com todo o aspecto de uma ferida cicatrizada? O pénis e os testículos podem ser feios, aye, mas são-no de uma maneira tola. A vagina é um órgão feio como tudo.
if (show_steps) { cout << "Apply Corrections..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, corrected_dir, "root");
create_directories(output_dir, corrected_plots, "plots");

apply_correction(out_jetmettau_allvertex, corr_detector_p6_z2_jetmettau, "detector_pythia6_z2_JetMETTau_2010A_", corr_detector_p8_4c_jetmettau, "detector_pythia8_4c_JetMETTau_2010A_", corrected_jetmettau, corrected_plots, "corrected_jetmettau_", detail, test);

apply_correction(out_jetmet_allvertex, corr_detector_p6_z2_jetmet, "detector_pythia6_z2_JetMET_2010A_", corr_detector_p8_4c_jetmet, "detector_pythia8_4c_JetMET_2010A_", corrected_jetmet, corrected_plots, "corrected_jetmet_", detail, test);

apply_correction(out_jet_allvertex, corr_detector_p6_z2_jet, "detector_pythia6_z2_Jet_2010B_", corr_detector_p8_4c_jet, "detector_pythia8_4c_Jet_2010B_", corrected_jet, corrected_plots, "corrected_jet_", detail, test);
if (show_steps) { cout << "The correction factors were sucessfully applied!"<<endl; }
}


//plot uncorrected control distributions
if (plot_control_dist)
{
if (show_steps) { cout << "Ploting Control Distributions..."<<endl; }
create_directories(output_dir, control_plots, "plots");
if (show_steps) { cout << "Ploting MC on Generator Level..."<<endl; }
//plot_mc_gen_level(output_dir+"histograms/xsec_mc_gen/",output_dir+"control_dist/",detail,show_results);
control_plots(out_p6_z2_gen_allvertex, out_p8_4c_gen_allvertex, out_p6_z2_det_allvertex, out_p8_4c_det_allvertex, merged_data, control_plots, "control_", detail, show_results);

control_plots(out_p6_z2_gen_allvertex_jetmettau, out_p8_4c_gen_allvertex_jetmettau, out_p6_z2_det_allvertex_jetmettau, out_p8_4c_det_allvertex_jetmettau, out_jetmettau_allvertex, control_plots, "control_JetMETTau_2010A_", detail, show_results);

control_plots(out_p6_z2_gen_allvertex_jetmet, out_p8_4c_gen_allvertex_jetmet, out_p6_z2_det_allvertex_jetmet, out_p8_4c_det_allvertex_jetmet, out_jetmet_allvertex, control_plots, "control_JetMET_2010A_", detail, show_results);

control_plots(out_p6_z2_gen_allvertex_jet, out_p8_4c_gen_allvertex_jet, out_p6_z2_det_allvertex_jet, out_p8_4c_det_allvertex_jet, out_jet_allvertex, control_plots, "control_Jet_2010B_", detail, show_results);
if (show_steps) { cout << "The Control Plots have been sucessfully generated!"<<endl; }
}


//merge the different datasets
if (merge_data)
{
if (show_steps) { cout << "Merging Data..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, merged_data_plots, "plots");
merge_data(corrected_jetmettau, corrected_jetmet, corrected_jet, merged_data_plots, merged_data,detail, show_results, test);
if (show_steps) { cout << "The data was sucessfully merged!"<<endl; }
}


//estimate model uncertainty
if (compute_model_uncertainty)
{
if (show_steps) { cout << "Calculating Uncertainties..."<<endl; }
if (show_steps) { cout << "Estimating Model Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, model_unc_dir, "root");
create_directories(output_dir, model_unc_plots, "plots");
//estimate_model_uncertainty(corr_final_p6_z2, corr_final_p8_4c, out_jetmettau_allvertex, model_uncertainty, model_unc_plots, detail, show_results, test); //needed merged - maybe not needed

estimate_model_uncertainty(corr_detector_p6_z2_jetmettau, corr_detector_p8_4c_jetmettau, out_jetmettau_allvertex, model_unc_jetmettau, "JetMETTau_2010A", model_unc_plots, detail, show_results, test);
estimate_model_uncertainty(corr_detector_p6_z2_jetmet, corr_detector_p8_4c_jetmet, out_jetmet_allvertex, model_unc_jetmet, "JetMET_2010A", model_unc_plots, detail, show_results, test);
estimate_model_uncertainty(corr_detector_p6_z2_jet, corr_detector_p8_4c_jet, out_jet_allvertex, model_unc_jet, "Jet_2010B", model_unc_plots, detail, show_results, test);

//merging
merge_uncertainties(model_unc_jetmettau, model_unc_jetmet, model_unc_jet, model_unc_merged, "model_unc_", "merged_model_unc_", model_unc_plots, detail, show_results, test);
if (show_steps) { cout << "The model uncertainty was sucessfully calculated!"<<endl; }
}

//estimate jes uncertainty
if (compute_jes_uncertainty)
{
if (show_steps) { cout << "Calculating Uncertainties..."<<endl; }
if (show_steps) { cout << "Estimating JES Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, jes_unc_dir, "root");
create_directories(output_dir, jes_unc_plots, "plots");
//estimate_model_uncertainty(corr_final_p6_z2, corr_final_p8_4c, out_jetmettau_allvertex, model_uncertainty, model_unc_plots, detail, show_results, test); //needed merged - maybe not needed

estimate_jes_uncertainty(out_jetmettau_allvertex, out_jetmettau_up, jes_up_unc_jetmettau, "JetMETTau_2010A_up_", "jes_unc_up_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jetmettau_allvertex, out_jetmettau_down, jes_down_unc_jetmettau, "JetMETTau_2010A_down_", "jes_unc_down_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jetmet_allvertex, out_jetmet_up, jes_up_unc_jetmet, "JetMET_2010A_up", "jes_unc_up_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jetmet_allvertex, out_jetmet_down, jes_down_unc_jetmet, "JetMET_2010A_down", "jes_unc_down_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jet_allvertex, out_jet_up, jes_up_unc_jet, "Jet_2010A_up", "jes_unc_up_", jes_unc_plots, detail, show_results, test);
estimate_jes_uncertainty(out_jet_allvertex, out_jet_down, jes_down_unc_jet, "Jet_2010A_down", "jes_unc_down_", jes_unc_plots, detail, show_results, test);

//merging
merge_uncertainties(jes_up_unc_jetmettau, jes_up_unc_jetmet, jes_up_unc_jet, jes_unc_up_merged, "jes_unc_up_", "merged_jes_unc_up_", jes_unc_plots, detail,  show_results, test);
merge_uncertainties(jes_down_unc_jetmettau, jes_down_unc_jetmet, jes_down_unc_jet, jes_unc_down_merged, "jes_unc_down_", "merged_jes_unc_down_", jes_unc_plots, detail, show_results, test);
if (show_steps) { cout << "The JES uncertainty was sucessfully calculated!"<<endl; }
}


//estimate total uncertainty
if (compute_total_uncertainty)
{
if (show_steps) { cout << "Computing Total Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, total_unc_dir, "root");
create_directories(output_dir, total_unc_plots, "plots");

total_uncertainty(jes_unc_up_merged, model_unc_merged, merged_data, total_unc_up, "merged_jes_unc_up_", "total_unc_up_", total_unc_plots, detail, show_results, test);
total_uncertainty(jes_unc_down_merged, model_unc_merged, merged_data, total_unc_down, "merged_jes_unc_down_", "total_unc_down_", total_unc_plots, detail, show_results, test);
}

//estimate_combination_systematic
if (estimate_combination_systematic)
{
if (show_steps) { cout << "Estimate Combination Systematics..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, trigger_syst_dir, "root");
create_directories(output_dir, trigger_syst_plots, "plots");
if (show_steps) { cout << "JetMETTauMonitor"<<endl; }
//estimate_combination_systematic(jetmettaumon, systematics_trigger_jetmettaumon, n_files_jetmettaumon, "allvertex", detail, test);

if (show_steps) { cout << "JetMETTau"<<endl; }
//estimate_combination_systematic(jetmettau, systematics_trigger_jetmettau, n_files_jetmettau, "allvertex", detail, test);

if (show_steps) { cout << "JetMET"<<endl; }
//estimate_combination_systematic(jetmet, systematics_trigger_jetmet, n_files_jetmet, "allvertex", detail, test);

if (show_steps) { cout << "Jet"<<endl; }
estimate_combination_systematic(jet, systematics_trigger_jet, n_files_jet, "allvertex", detail, test);

if (show_steps) { cout << "Merging and ploting"<<endl; }
merge_plot_combination_systematic(systematics_trigger_jetmettau, systematics_trigger_jetmet, systematics_trigger_jet, merged_trig_comb_syst, trigger_syst_plots, detail, test);
if (show_steps) { cout << "The Estimation was sucessfully done!"<<endl; }
}


//apply uncertainty
if (apply_uncertainty)
{
if (show_steps) { cout << "Apply Uncertainty..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, data_unc_plots, "plots");

apply_uncertainty(merged_data, total_unc_up, total_unc_down, data_unc, data_unc_plots, detail, test);

if (show_steps) { cout << "The uncertainty was sucessfully applied!"<<endl; }
}


//convert predictions
if (convert_predictions)
{
if (show_steps) { cout << "Converting MC Predictions..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, mc_pred_dir, "root");

convert_prediction(source_herwig6, mc_pred_herwig6, detail);
//convert_prediction(source_herwigpp, mc_pred_herwigpp, detail);
convert_prediction(source_py6_ambt1, mc_pred_py6_ambt1, detail); //broken for 1M
//convert_prediction(source_py6_p11, mc_pred_py6_p11, detail); //completely broken
convert_prediction(source_py6_z2_nompi, mc_pred_py6_z2_nompi, detail); //broken for 1M
//convert_prediction(source_py6_z2, mc_pred_py6_z2, detail); //completely broken
//convert_prediction(source_py6_4c, mc_pred_py6_4c, detail); //completely broken!

if (show_steps) { cout << "The final plots were sucessfully generated!"<<endl; }
}


//ratio plots
if (create_ratio)
{
if (show_steps) { cout << "Ratio Plots..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, ratio_dir, "root");
create_directories(output_dir, ratio_plots, "plots");

create_ratios(merged_data, merged_data, "ak5PF_", ratio_merged_data, ratio_plots, "data_data_", detail, test);
create_ratios(data_unc, merged_data, "ak5PF_", ratio_data_unc, ratio_plots, "dataunc_data_", detail, test);
create_ratios(mc_pred_herwig6, merged_data, "", mc_ratio_herwig6, ratio_plots, "herwig6_data_", detail, test);
create_ratios(mc_pred_herwigpp, merged_data, "", mc_ratio_herwigpp, ratio_plots, "herwigpp_data_", detail, test);

if (show_steps) { cout << "The ratios were sucessfully computed!"<<endl; }
}

//final plots
if (do_final_plots)
{
if (show_steps) { cout << "Final Plots..."<<endl; }
create_directories(output_dir, final_plots, "plots");

final_plots(merged_data, data_unc, final_mc_list, final_label_list, final_prefix_list, final_plots, "xsec", detail, test);
final_plots(ratio_merged_data, ratio_data_unc, ratio_mc_list, ratio_label_list, ratio_prefix_list, final_plots, "ratio", detail, test);

if (show_steps) { cout << "The final plots were sucessfully generated!"<<endl; }
}

//creates the unfolding response matrix
if (create_unfolding_response)
{
if (show_steps) { cout << "Creating the RooUnfold Response MAtrixes..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, response_dir, "root");
if (show_steps) { cout << "MC without vertex reweighting..."<<endl; }
///create_unfolding_response(mc_p6_z2, out_p6_z2_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", "",  "", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", "",  "", detail, test);
if (show_steps) { cout << "MC with vertex reweighting as JetMETTau_2010A..."<<endl; }
///create_unfolding_response(mc_p6_z2, out_p6_z2_jetmettau_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmettau_v5,  "_v5", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_jetmettau_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmettau_v5,  "_v5", detail, test);
if (show_steps) { cout << "MC with vertex reweighting as JetMET_2010A..."<<endl; }
///create_unfolding_response(mc_p6_z2, out_p6_z2_jetmet_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jetmet_v5,  "_v5", detail, test);
///create_unfolding_response(mc_p8_4c, out_p8_4c_jetmet_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jetmet_v5,  "_v5", detail, test);
if (show_steps) { cout << "MC with vertex reweighting as Jet_2010B..."<<endl; }
///create_unfolding_response(mc_p6_z2, out_p6_z2_jet_unf_allvertex, lumi_p6_z2, n_files_p6_z2, "allvertex", vertex_weights_p6_z2_jet_v5,  "_v5", detail, test);
create_unfolding_response(mc_p8_4c, out_p8_4c_jet_unf_allvertex, lumi_p8_4c, n_files_p8_4c, "allvertex", vertex_weights_p8_4c_jet_v5,  "_v5", detail, test);
if (show_steps) { cout << "The RooUnfold Response Matrixes were been sucessfully created!"<<endl; }
}


//checks the response matrix
if (check_response_matrix)
{
if (show_steps) { cout << "Check the Response Matrixes..."<<endl; }
create_directories(output_dir, check_response_dir, "plots");
if (show_steps) { cout << "Pythia 6 - Tune Z2*..."<<endl; }
check_response_matrix(out_p6_z2_unf_allvertex, check_response_dir, "pythia6_z2_", detail, test);
check_response_matrix(out_p6_z2_jetmettau_unf_allvertex, check_response_dir, "pythia6_z2_jetmettau_", detail, test);
check_response_matrix(out_p6_z2_jetmet_unf_allvertex, check_response_dir, "pythia6_z2_jetmet_", detail, test);
check_response_matrix(out_p6_z2_jet_unf_allvertex, check_response_dir, "pythia6_z2_jet_", detail, test);
if (show_steps) { cout << "Pythia 8 - Tune 4C..."<<endl; }
check_response_matrix(out_p8_4c_unf_allvertex, check_response_dir, "pythia8_4c_", detail, test);
check_response_matrix(out_p8_4c_jetmettau_unf_allvertex, check_response_dir, "pythia8_4c_jetmettau_", detail, test);
check_response_matrix(out_p8_4c_jetmet_unf_allvertex, check_response_dir, "pythia8_4c_jetmet_", detail, test);
check_response_matrix(out_p8_4c_jet_unf_allvertex, check_response_dir, "pythia8_4c_jet_", detail, test);
if (show_steps) { cout << "The Response Matrixes were been sucessfully checked!"<<endl; }
}


//unfolds the distribution
if (unfold)
{
if (show_steps) { cout << "Unfolding..."<<endl; }
create_directories(output_dir, histograms_dir, "root");
create_directories(output_dir, unfolded_dir, "root");
create_directories(output_dir, unfolding_plots, "plots");
double tau[24];
for (int i=0; i < 24; i++)
{ tau[i] = 0.0; }
if (show_steps) { cout << "Unfolding without Vertex Reweight..."<<endl; }
unfolding(out_p6_z2_unf_allvertex, merged_data, out_p6_z2_gen_norm, out_p6_z2_unfolded_data, tau[0], "data_unfolded_pythia6_z2_", detail, test);
unfolding(out_p8_4c_unf_allvertex, merged_data, out_p8_4c_gen_norm, out_p8_4c_unfolded_data, tau[1], "data_unfolded_pythia8_4c_", detail, test);
unfolding(out_p6_z2_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_norm, out_p6_z2_unfolded_p6_z2, tau[2], "pythia6_z2_unfolded_pythia6_z2_", detail, test);
unfolding(out_p8_4c_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_norm, out_p8_4c_unfolded_p8_4c, tau[3], "pythia8_4c_unfolded_pythia8_4c_", detail, test);
unfolding(out_p8_4c_unf_allvertex, out_p6_z2_det_allvertex, out_p6_z2_gen_norm, out_p8_4c_unfolded_p6_z2, tau[4], "pythia6_z2_unfolded_pythia8_4c_", detail, test);
unfolding(out_p6_z2_unf_allvertex, out_p8_4c_det_allvertex, out_p8_4c_gen_norm, out_p6_z2_unfolded_p8_4c, tau[5], "pythia8_4c_unfolded_pythia6_z2_", detail, test);
if (show_steps) { cout << "Unfolding for JetMETTau_2010A..."<<endl; }
unfolding(out_p6_z2_jetmettau_unf_allvertex, out_jetmettau_allvertex, out_p6_z2_gen_jetmettau_norm, out_p6_z2_unfolded_data_jetmettau, tau[6], "data_unfolded_pythia6_z2_jetmettau_", detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_jetmettau_allvertex, out_p8_4c_gen_jetmettau_norm, out_p8_4c_unfolded_data_jetmettau, tau[7], "data_unfolded_pythia8_4c_jetmettau_", detail, test);
unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_jetmettau_norm, out_p6_z2_unfolded_p6_z2_jetmettau, tau[8], "pythia6_z2_unfolded_pythia6_z2_jetmettau_", detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_jetmettau_norm, out_p8_4c_unfolded_p8_4c_jetmettau, tau[9], "pythia8_4c_unfolded_pythia8_4c_jetmettau_", detail, test);
unfolding(out_p8_4c_jetmettau_unf_allvertex, out_p6_z2_det_allvertex_jetmettau, out_p6_z2_gen_jetmettau_norm, out_p8_4c_unfolded_p6_z2_jetmettau, tau[10], "pythia6_z2_unfolded_pythia8_4c_jetmettau_", detail, test);
unfolding(out_p6_z2_jetmettau_unf_allvertex, out_p8_4c_det_allvertex_jetmettau, out_p8_4c_gen_jetmettau_norm, out_p6_z2_unfolded_p8_4c_jetmettau, tau[11], "pythia8_4c_unfolded_pythia6_z2_jetmettau_", detail, test);
if (show_steps) { cout << "Unfolding for JetMET_2010A..."<<endl; }
unfolding(out_p6_z2_jetmet_unf_allvertex, out_jetmet_allvertex, out_p6_z2_gen_jetmet_norm, out_p6_z2_unfolded_data_jetmet, tau[12], "data_unfolded_pythia6_z2_jetmet_", detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_jetmet_allvertex, out_p8_4c_gen_jetmet_norm, out_p8_4c_unfolded_data_jetmet, tau[13], "data_unfolded_pythia8_4c_jetmet_", detail, test);
unfolding(out_p6_z2_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_jetmet_norm, out_p6_z2_unfolded_p6_z2_jetmet, tau[14], "pythia6_z2_unfolded_pythia6_z2_jetmet_", detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_jetmet_norm, out_p8_4c_unfolded_p8_4c_jetmet, tau[15], "pythia8_4c_unfolded_pythia8_4c_jetmet_", detail, test);
unfolding(out_p8_4c_jetmet_unf_allvertex, out_p6_z2_det_allvertex_jetmet, out_p6_z2_gen_jetmet_norm, out_p8_4c_unfolded_p6_z2_jetmet, tau[16], "pythia6_z2_unfolded_pythia8_4c_jetmet_", detail, test);
unfolding(out_p6_z2_jetmet_unf_allvertex, out_p8_4c_det_allvertex_jetmet, out_p8_4c_gen_jetmet_norm, out_p6_z2_unfolded_p8_4c_jetmet, tau[17], "pythia8_4c_unfolded_pythia6_z2_jetmet_", detail, test);
if (show_steps) { cout << "Unfolding for Jet_2010B..."<<endl; }
unfolding(out_p6_z2_jet_unf_allvertex, out_jet_allvertex, out_p6_z2_gen_jet_norm, out_p6_z2_unfolded_data_jet, tau[18], "data_unfolded_pythia6_z2_jet_", detail, test);
unfolding(out_p8_4c_jet_unf_allvertex, out_jet_allvertex, out_p8_4c_gen_jet_norm, out_p8_4c_unfolded_data_jet, tau[19], "data_unfolded_pythia8_4c_jet_", detail, test);
unfolding(out_p6_z2_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p6_z2_unfolded_p6_z2_jet, tau[20], "pythia6_z2_unfolded_pythia6_z2_jet_", detail, test);
unfolding(out_p8_4c_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p8_4c_unfolded_p8_4c_jet, tau[21], "pythia8_4c_unfolded_pythia8_4c_jet_", detail, test);
unfolding(out_p8_4c_jet_unf_allvertex, out_p6_z2_det_allvertex_jet, out_p6_z2_gen_jet_norm, out_p8_4c_unfolded_p6_z2_jet, tau[22], "pythia6_z2_unfolded_pythia8_4c_jet_", detail, test);
unfolding(out_p6_z2_jet_unf_allvertex, out_p8_4c_det_allvertex_jet, out_p8_4c_gen_jet_norm, out_p6_z2_unfolded_p8_4c_jet, tau[23], "pythia8_4c_unfolded_pythia6_z2_jet_", detail, test);
if (show_steps) { cout << "Show Tau Values..."<<endl; }
cout << "Tau Values No Reweight:     " << tau[0] << " " << tau[1] << " " << tau[2] << " " << tau[3] << " " << tau[4] << " " << tau[5] << endl;
cout << "Tau Values JetMETTau_2010A: " << tau[6] << " " << tau[7] << " " << tau[8] << " " << tau[9] << " " << tau[10] << " " << tau[11] << endl;
cout << "Tau Values JetMET_2010A:    " << tau[12] << " " << tau[13] << " " << tau[14] << " " << tau[15] << " " << tau[16] << " " << tau[17] << endl;
cout << "Tau Values Jet_2010B:       " << tau[18] << " " << tau[19] << " " << tau[20] << " " << tau[21] << " " << tau[22] << " " << tau[23] << endl;
if (show_steps) { cout << "The Unfolding has sucessful!"<<endl; }
}

}
