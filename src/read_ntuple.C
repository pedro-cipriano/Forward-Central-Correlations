// Pedro Cipriano, Dec 2011
// DESY, CMS
// Last Update: 13 Dez 2012
//
//read_ntuple(string *data_in, string data_out, double *data_lumi, int n_files, string data_type = "DATA", string sel_mode = "1vertex", bool detail = false)
// reads the Ntuples and created a root file with cross-section histograms

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <string>

#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "common_methods.h"

using namespace std;

void read_ntuple(string *data_in, string data_out, double *data_lumi, int n_files, string data_type = "DATA", string sel_mode = "1vertex", string vertex_weights = "", TString vertex_sufix = "_v0", bool detail = false, bool test = false)
{ 

//output the configuration
   if (detail) { cout<<"Read NTuple Configuration"<<endl; }
   if (detail) { cout<<"Data Type :            "<<data_type<<endl; }
   if (detail) { cout<<"Selection mode :       "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files :      "<<n_files<<endl; }
   if (detail) { cout<<"Vertex Weights File :  "<<vertex_weights<<endl; }
   if (detail) { cout<<"Vertex Weights Sufix : "<<vertex_sufix<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :          "<<data_out<<endl; }
   if (detail) { cout<<"Detail level :         "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :            "<<test<<endl; }

//opening the vertex weight file
    bool vertex_weighted = false;
    TH1D *vertex_hist = 0;
    if ((vertex_weights != "") and (data_type == "MC_DET" or data_type == "MC_GEN"))
    {
    if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
    TFile *vertex_file = new TFile( vertex_weights.c_str() );
    vertex_file->GetObject("vertex_weights"+vertex_sufix,vertex_hist);
    if (vertex_hist != 0) { vertex_weighted = true; }
    }

//setting the prefix
   TString prefix;
   if (data_type == "DATA" or data_type == "MC_DET") { prefix = "ak5PF_"; }
   if (data_type == "MC_GEN") { prefix = "ak5Gen_"; }
   if (detail) { cout<<"Prefix :          "<<prefix<<endl; }

//main cuts
   double pt_min = 35.0;
   double gap_req = 20.0;
 
   double pt_central, eta_central, eta_central2, phi_central, chm_central, elm_central;
   double pt_forward, eta_forward, eta_forward2, phi_forward, chm_forward, elm_forward;
   double pt2_central, eta2_central, eta2_central2, phi2_central, pt_central_check;
   double pt2_forward, eta2_forward, eta2_forward2, phi2_forward;
   double pt = 0.0, eta = 0.0, phi = 0.0, gen_pt = 0.0, gen_eta = 0.0, gen_phi = 0.0, chm = 0.0, elm = 0.0;
   double unc, pt_up, pt_down;
   double delta_eta, delta_phi, delta_phi2, delta_pt;
   bool pv_pass, pass_gap, pass_tight, pass_check;
   double pt_total_gap, pt_leading_gap, eta_gap, eta_star_inside, phi_gap, chm_gap, elm_gap;
   double pt_total_outside, pt_leading_outside, eta_outside, phi_outside, chm_outside;
   int multiplicity_inside, multiplicity_outside;
   double elm_outside, deta_out1, deta_out2;
   double prescale, prescale_prov = 0.0;
   unsigned int number_of_jets = 0.0;
   double leading_pt;
   double vertex_factor = 1.0;

   int selected_gap = 0;
   int selected_nogap = 0;

//declare the main binning
   int deta_nbins = 4;
   int dphi_nbins = 7;

   double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

   int all_nbins = 11;
   int all2_nbins = 13;
   int in_nbins = 9;
   int out_nbins = 9;
   int dpt_nbins = 10;

   double all_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
   double all2_bins[14] = {-10, 5, 10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
   double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
   double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
   double dpt_bins[11] = {0, 2.5, 5, 10, 15, 20, 25, 30, 40, 60, 100};

   int cent_nbins = 7;
   int forw_nbins = 7;

   double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
   double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};

   int etac_nbins = 6;
   double etac_bins[7] = {-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8};

   int etaf_nbins = 7;
   double etaf_bins[8] = {-4.7,-4.2,-3.7,-3.2,3.2,3.7,4.2,4.7};

   int eta_nbins = 14;
   double eta_bins[15] = {-4.7,-4.2,-3.7,-3.2,-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8,3.2,3.7,4.2,4.7};

   int etastar_nbins = 12;
   double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

   int deta_out_nbins = 6;
   double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//declare the histograms
   TH1D *hist_delta_phi;
   TH1D *hist_delta_phi_fine;
   TH1D *hist_delta_phi_gap;
   TH1D *hist_delta_phi_nogap;
   TH1D *hist_delta_phi_deta1;
   TH1D *hist_delta_phi_deta2;
   TH1D *hist_delta_phi_deta3;
   TH1D *hist_delta_phi_deta4;
   TH1D *hist_delta_phi_deta1_gap;
   TH1D *hist_delta_phi_deta2_gap;
   TH1D *hist_delta_phi_deta3_gap;
   TH1D *hist_delta_phi_deta4_gap;
   TH1D *hist_delta_phi_deta1_nogap;
   TH1D *hist_delta_phi_deta2_nogap;
   TH1D *hist_delta_phi_deta3_nogap;
   TH1D *hist_delta_phi_deta4_nogap;

   TH1D *hist_delta_eta;
   TH1D *hist_delta_eta_gap;
   TH1D *hist_delta_eta_nogap;

   TH1D *hist_total_pt_inside_gap;
   TH1D *hist_leading_pt_inside_gap;
   TH1D *hist_leading_eta_inside_gap;
   TH1D *hist_leading_phi_inside_gap;
   TH1D *hist_leading_chm_inside_gap;
   TH1D *hist_leading_elm_inside_gap;
   TH1D *hist_multiplicity_inside_gap;
   TH1D *hist_total_pt_outside_gap;
   TH1D *hist_leading_pt_outside_gap;
   TH1D *hist_delta_eta_outside_gap;
   TH1D *hist_leading_eta_outside_gap;
   TH1D *hist_leading_eta_star_inside_gap;
   TH1D *hist_leading_phi_outside_gap;
   TH1D *hist_leading_chm_outside_gap;
   TH1D *hist_leading_elm_outside_gap;
   TH1D *hist_multiplicity_outside_gap;

   TH1D *hist_leading_pt;
   TH1D *hist_leading_pt_fine;
   TH1D *hist_leading_central_pt_fine;
   TH1D *hist_leading_forward_pt_fine;

   TH1D *hist_leading_central_pt;
   TH1D *hist_leading_central_pt_gap;
   TH1D *hist_leading_central_pt_nogap;
   TH1D *hist_leading_central_eta;
   TH1D *hist_leading_central_phi;
   TH1D *hist_leading_central_chm;
   TH1D *hist_leading_central_elm;
   TH1D *hist_leading_forward_pt;
   TH1D *hist_leading_forward_pt_gap;
   TH1D *hist_leading_forward_pt_nogap;
   TH1D *hist_leading_forward_eta;
   TH1D *hist_leading_forward_phi;
   TH1D *hist_leading_forward_chm;
   TH1D *hist_leading_forward_elm;
   TH1D *hist_leading_eta;

   TH1D *hist_vertex_selected;
   TH1D *hist_pu_selected;
   TH1D *hist_pvz_selected;
   TH1D *hist_multiplicity;
   TH1D *hist_chm_all;
   TH1D *hist_elm_all;
   TH1D *hist_pt_all;
   TH1D *hist_eta_all;
   TH1D *hist_phi_all;

   TH1D *hist_inclusive_leading_central_pt;
   TH1D *hist_inclusive_leading_forward_pt;
   TH1D *hist_inclusive_central_pt;
   TH1D *hist_inclusive_forward_pt;
   TH1D *hist_leading_central_pt_check;
   
   TH2D *hist_central_pt_rel;
   TH2D *hist_forward_pt_rel;
   TH1D *hist_delta_pt_central_rel;
   TH1D *hist_delta_pt_forward_rel;
   TH1D *hist_delta_phi_central_rel;
   TH1D *hist_delta_phi_forward_rel;
   TH1D *hist_delta_phi_central_rel_small;
   TH1D *hist_delta_phi_forward_rel_small;
   TH1D *hist_delta_phi_central_rel_medium;
   TH1D *hist_delta_phi_forward_rel_medium;
   TH1D *hist_delta_phi_central_rel_large;
   TH1D *hist_delta_phi_forward_rel_large;

   TH1D *hist_trigger_passed;
   TH1D *hist_events;

  hist_delta_phi =  new TH1D(prefix+"delta_phi","#Delta#phi;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_fine =  new TH1D(prefix+"delta_phi_fine","#Delta#phi;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", 14, 0, 3.15);
  hist_delta_phi_gap =  new TH1D(prefix+"delta_phi_gap","#Delta#phi when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_nogap =  new TH1D(prefix+"delta_phi_nogap","#Delta#phi when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1 =  new TH1D(prefix+"delta_phi_deta1","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2 =  new TH1D(prefix+"delta_phi_deta2","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3 =  new TH1D(prefix+"delta_phi_deta3","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4 =  new TH1D(prefix+"delta_phi_deta4","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_gap =  new TH1D(prefix+"delta_phi_deta1_gap","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_gap =  new TH1D(prefix+"delta_phi_deta2_gap","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_gap =  new TH1D(prefix+"delta_phi_deta3_gap","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_gap =  new TH1D(prefix+"delta_phi_deta4_gap","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_nogap =  new TH1D(prefix+"delta_phi_deta1_nogap","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]i", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_nogap =  new TH1D(prefix+"delta_phi_deta2_nogap","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_nogap =  new TH1D(prefix+"delta_phi_deta3_nogap","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_nogap =  new TH1D(prefix+"delta_phi_deta4_nogap","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_eta =  new TH1D(prefix+"delta_eta","#Delta#eta;|#Delta#eta|;#frac{d#sigma}{d#Delta#eta} [pb]", deta_nbins, deta_bins);
  hist_delta_eta_gap =  new TH1D(prefix+"delta_eta_gap","#Delta#eta when requiring a gap;|#Delta#eta|;#frac{d#sigma}{d#Delta#eta} [pb]", deta_nbins, deta_bins);
  hist_delta_eta_nogap =  new TH1D(prefix+"delta_eta_nogap","#Delta#eta when vetoing a gap;|#Delta#eta|;#frac{d#sigma}{d#Delta#eta} [pb]", deta_nbins, deta_bins);
  hist_leading_pt_inside_gap =  new TH1D(prefix+"leading_pt_inside_gap","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  hist_total_pt_inside_gap =  new TH1D(prefix+"total_pt_inside_gap","Total Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  hist_leading_pt_outside_gap =  new TH1D(prefix+"leading_pt_outside_gap","Leading Jet p_{T} outside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", out_nbins, out_bins);
  hist_total_pt_outside_gap =  new TH1D(prefix+"total_pt_outside_gap","Total Jet p_{T} outside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", out_nbins, out_bins);
  hist_delta_eta_outside_gap =  new TH1D(prefix+"delta_eta_outside_gap","#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
  hist_leading_eta_inside_gap =  new TH1D(prefix+"leading_eta_inside_gap","Leading Jet #eta inside the gap;#eta;#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
  hist_leading_eta_outside_gap =  new TH1D(prefix+"leading_eta_outside_gap","Leading Jet #eta outside the gap;eta;#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
  hist_leading_phi_inside_gap =  new TH1D(prefix+"leading_phi_inside_gap","Leading Jet #phi inside the gap;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
  hist_leading_phi_outside_gap =  new TH1D(prefix+"leading_phi_outside_gap","Leading Jet #phi ouside the gap;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
  hist_leading_chm_inside_gap =  new TH1D(prefix+"leading_chm_inside_gap","charged hadron multiplicity central;Charged Hadron Multiplicity for inside Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_chm_outside_gap =  new TH1D(prefix+"leading_chm_outside_gap","charged hadron multiplicity forward;Charged Hadron Multiplicity for outside Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_elm_inside_gap =  new TH1D(prefix+"leading_elm_inside_gap","electron multiplicity inside;Electron Multiplicity for inside Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_elm_outside_gap =  new TH1D(prefix+"leading_elm_outside_gap","electron multiplicity outside;Electron Multiplicity for outside Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_eta_star_inside_gap =  new TH1D(prefix+"leading_eta_star_inside_gap","Leading Jet #eta^{*} inside the gap;#eta^{*};#frac{d#sigma}{d#eta^{*}} [pb]", etastar_nbins, etastar_bins);
  hist_leading_pt =  new TH1D(prefix+"leading_pt","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", all_nbins, all_bins);
  hist_multiplicity_inside_gap =  new TH1D(prefix+"multiplicity_inside_gap","Jet Multiplicity inside the gap;N_{jets};#frac{d#sigma}{d#N} [pb]", 10, 0, 10);
  hist_multiplicity_outside_gap =  new TH1D(prefix+"multiplicity_outside_gap","Jet Multiplicity outside the gap;N_{jets};#frac{d#sigma}{d#N} [pb]", 10, 0, 10);
  hist_leading_pt_fine =  new TH1D(prefix+"leading_pt_fine","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", 85, 30, 200);
  hist_leading_central_pt =  new TH1D(prefix+"leading_central_pt","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
  hist_leading_central_pt_fine =  new TH1D(prefix+"leading_central_pt_fine","Leading Central Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", 85, 30, 200);
  hist_leading_central_pt_gap =  new TH1D(prefix+"leading_central_pt_gap","Leading Central Jet p_{T} when requiring a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
  hist_leading_central_pt_nogap =  new TH1D(prefix+"leading_central_pt_nogap","Leading Central Jet p_{T} when vetoing a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
  hist_leading_forward_pt =  new TH1D(prefix+"leading_forward_pt","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
  hist_leading_forward_pt_fine =  new TH1D(prefix+"leading_forward_pt_fine","Leading Forward Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", 85, 30, 200);
  hist_leading_forward_pt_gap =  new TH1D(prefix+"leading_forward_pt_gap","Leading Forward Jet p_{T} when requiring a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
  hist_leading_forward_pt_nogap =  new TH1D(prefix+"leading_forward_pt_nogap","Leading Forward Jet p_{T} when vetoing a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
  hist_leading_central_phi =  new TH1D(prefix+"leading_central_phi","Leading Central Jet #phi;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
  hist_leading_forward_phi =  new TH1D(prefix+"leading_forward_phi","Leading Forward Jet #phi;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
  hist_leading_central_eta =  new TH1D(prefix+"leading_central_eta","Leading Central Jet #eta;|#eta|;#frac{d#sigma}{d#eta} [pb]", etac_nbins, etac_bins);
  hist_leading_forward_eta =  new TH1D(prefix+"leading_forward_eta","Leading Forward Jet #eta;|#eta|;#frac{d#sigma}{d#eta} [pb]", etaf_nbins, etaf_bins);
  hist_leading_central_chm =  new TH1D(prefix+"leading_central_chm","charged hadron multiplicity central;Charged Hadron Multiplicity for Central Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_forward_chm =  new TH1D(prefix+"leading_forward_chm","charged hadron multiplicity forward;Charged Hadron Multiplicity for Forward Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_central_elm =  new TH1D(prefix+"leading_central_elm","electron multiplicity forward;Electron Multiplicity for Central Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_forward_elm =  new TH1D(prefix+"leading_forward_elm","electron multiplicity forward;Electron Multiplicity for Forward Jet;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_leading_eta =  new TH1D(prefix+"leading_eta","Leading Jet #eta;|#eta|;#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
  hist_multiplicity =  new TH1D(prefix+"multiplicity","Jet multiplicity;Multiplicity;#frac{d#sigma}{dN} [pb]", 30, 0, 30);
  hist_vertex_selected =  new TH1D(prefix+"vertex_selected","Vertex multiplicity in selected events;Multiplicity;#frac{d#sigma}{dN} [pb]", 15, 0, 15);
  hist_pu_selected =  new TH1D(prefix+"pu_selected","Pileup in selected events;Pileup;#frac{d#sigma}{dPileup} [pb]", 15, 0, 15);
  hist_pvz_selected =  new TH1D(prefix+"pvz_selected","z of the Primary Vertex in selected events;z-position;#frac{d#sigma}{dz} [pb]", 200, -100, 100);
  hist_chm_all =  new TH1D(prefix+"chm_all","charged hadron multiplicity;Charged Hadron Multiplicity;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_elm_all =  new TH1D(prefix+"elm_all","electron multiplicity;Electron Multiplicity;#frac{d#sigma}{dN} [pb]", 50, 0, 50);
  hist_pt_all =  new TH1D(prefix+"pt_all","p_{T} all;p_{T}^{all};#frac{d#sigma}{dp_{T}} [#frac{pb}{GeV}]", all_nbins, all_bins);
  hist_eta_all =  new TH1D(prefix+"eta_all","#eta;#eta^{all};#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
  hist_phi_all =  new TH1D(prefix+"phi_all","#phi;#phi^{all};#frac{d#sigma}{d#phi} [#frac{pb}{rad}]", 14, -3.15, 3.15);
  hist_inclusive_leading_central_pt =  new TH1D(prefix+"inclusive_leading_central_pt","Inclusive Leading Central Jet pT;p_T^{central};#frac{d#sigma}{d#phi} [pb]", cent_nbins, cent_bins);
  hist_inclusive_leading_forward_pt =  new TH1D(prefix+"inclusive_leading_forward_pt","Inclusive Leading Forward Jet pT;p_T^{forward};#frac{d#sigma}{d#phi} [pb]", forw_nbins, forw_bins);
  hist_inclusive_central_pt =  new TH1D(prefix+"inclusive_central_pt","Inclusive Central Jet pT;p_T^{central};#frac{d#sigma}{d#phi} [pb]", cent_nbins, cent_bins);
  hist_inclusive_forward_pt =  new TH1D(prefix+"inclusive_forward_pt","Inclusive Forward Jet pT;p_T^{forward};#frac{d#sigma}{d#phi} [pb]", forw_nbins, forw_bins);
  hist_leading_central_pt_check =  new TH1D(prefix+"leading_central_pt_check","Leading Central Jet pT with eta < 2.5;p_T^{forward};#frac{d#sigma}{d#phi} [pb]", 250, 0, 500);
  
  hist_central_pt_rel =  new TH2D(prefix+"central_pt_rel","Central Jet pT relations;Leading 1 p_T^{central};Leading 2 p_T^{central}", all_nbins, all_bins, all2_nbins, all2_bins);
  hist_forward_pt_rel =  new TH2D(prefix+"forward_pt_rel","Forward Jet pT relations;Leading 1 p_T^{forward};Leading 2 p_T^{forward}", all_nbins, all_bins, all2_nbins, all2_bins);
  hist_delta_phi_central_rel =  new TH1D(prefix+"delta_phi_central_rel","#Delta#phi central relations;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_forward_rel =  new TH1D(prefix+"delta_phi_forward_rel","#Delta#phi forward relations;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_pt_central_rel =  new TH1D(prefix+"delta_pt_central_rel","#Deltap_{T} central relations;|#Deltap_{T}|;#frac{d#sigma}{d#Deltap_{T}} [pb]", dpt_nbins, dpt_bins);
  hist_delta_pt_forward_rel =  new TH1D(prefix+"delta_pt_forward_rel","#Deltap_{T} forward relations;|#Deltap_{T}|;#frac{d#sigma}{d#Deltap_{T}} [pb]", dpt_nbins, dpt_bins);
  hist_delta_phi_central_rel_small =  new TH1D(prefix+"delta_phi_central_rel_small","#Delta#phi central relations small;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_forward_rel_small =  new TH1D(prefix+"delta_phi_forward_rel_small","#Delta#phi forward relations small;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_central_rel_medium =  new TH1D(prefix+"delta_phi_central_rel_medium","#Delta#phi central relations medium;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_forward_rel_medium =  new TH1D(prefix+"delta_phi_forward_rel_medium","#Delta#phi forward relations medium;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_central_rel_large =  new TH1D(prefix+"delta_phi_central_rel_large","#Delta#phi central relations large;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_forward_rel_large =  new TH1D(prefix+"delta_phi_forward_rel_large","#Delta#phi forward relations large;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);

  hist_trigger_passed =  new TH1D(prefix+"trigger_passed","Events from each trigger;Trigger;Events", 3, 0 ,3);
  hist_events =  new TH1D(prefix+"events","Selection Chain;Type;# Events", 8, 0, 8);

     hist_delta_phi->Sumw2();
     hist_delta_phi_fine->Sumw2();
     hist_delta_phi_gap->Sumw2();
     hist_delta_phi_nogap->Sumw2();
     hist_delta_phi_deta1->Sumw2();
     hist_delta_phi_deta2->Sumw2();
     hist_delta_phi_deta3->Sumw2();
     hist_delta_phi_deta4->Sumw2();
     hist_delta_phi_deta1_gap->Sumw2();
     hist_delta_phi_deta2_gap->Sumw2();
     hist_delta_phi_deta3_gap->Sumw2();
     hist_delta_phi_deta4_gap->Sumw2();
     hist_delta_phi_deta1_nogap->Sumw2();
     hist_delta_phi_deta2_nogap->Sumw2();
     hist_delta_phi_deta3_nogap->Sumw2();
     hist_delta_phi_deta4_nogap->Sumw2();
     hist_delta_eta->Sumw2();
     hist_delta_eta_gap->Sumw2();
     hist_delta_eta_nogap->Sumw2();
     hist_total_pt_inside_gap->Sumw2();
     hist_leading_pt_inside_gap->Sumw2();
     hist_leading_eta_inside_gap->Sumw2();
     hist_leading_phi_inside_gap->Sumw2();
     hist_leading_chm_inside_gap->Sumw2();
     hist_leading_elm_inside_gap->Sumw2();
     hist_multiplicity_inside_gap->Sumw2();
     hist_total_pt_outside_gap->Sumw2();
     hist_delta_eta_outside_gap->Sumw2();
     hist_leading_pt_outside_gap->Sumw2();
     hist_leading_eta_outside_gap->Sumw2();
     hist_leading_phi_outside_gap->Sumw2();
     hist_leading_chm_outside_gap->Sumw2();
     hist_leading_elm_outside_gap->Sumw2();
     hist_multiplicity_outside_gap->Sumw2();
     hist_leading_eta_star_inside_gap->Sumw2();
     hist_leading_pt->Sumw2();
     hist_leading_pt_fine->Sumw2();
     hist_leading_central_pt_fine->Sumw2();
     hist_leading_forward_pt_fine->Sumw2();
     hist_leading_central_pt->Sumw2();
     hist_leading_central_pt_gap->Sumw2();
     hist_leading_central_pt_nogap->Sumw2();
     hist_leading_central_eta->Sumw2();
     hist_leading_central_phi->Sumw2();
     hist_leading_central_chm->Sumw2();
     hist_leading_central_elm->Sumw2();
     hist_leading_forward_pt->Sumw2();
     hist_leading_forward_pt_gap->Sumw2();
     hist_leading_forward_pt_nogap->Sumw2();
     hist_leading_forward_eta->Sumw2();
     hist_leading_forward_phi->Sumw2();
     hist_leading_forward_chm->Sumw2();
     hist_leading_forward_elm->Sumw2();
     hist_leading_eta->Sumw2();
     hist_vertex_selected->Sumw2();
     hist_pvz_selected->Sumw2();
     hist_pu_selected->Sumw2();
     hist_multiplicity->Sumw2();
     hist_chm_all->Sumw2();
     hist_elm_all->Sumw2();
     hist_pt_all->Sumw2();
     hist_eta_all->Sumw2();
     hist_phi_all->Sumw2();
     hist_inclusive_leading_central_pt->Sumw2();
     hist_inclusive_leading_forward_pt->Sumw2();
     hist_inclusive_central_pt->Sumw2();
     hist_inclusive_forward_pt->Sumw2();
     hist_leading_central_pt_check->Sumw2(); 
     hist_central_pt_rel->Sumw2();
     hist_forward_pt_rel->Sumw2();
     hist_delta_phi_central_rel->Sumw2();
     hist_delta_phi_forward_rel->Sumw2();
     hist_delta_pt_central_rel->Sumw2();
     hist_delta_pt_forward_rel->Sumw2();
     hist_delta_phi_central_rel_small->Sumw2();
     hist_delta_phi_forward_rel_small->Sumw2();
     hist_delta_phi_central_rel_medium->Sumw2();
     hist_delta_phi_forward_rel_medium->Sumw2();
     hist_delta_phi_central_rel_large->Sumw2();
     hist_delta_phi_forward_rel_large->Sumw2();

     hist_trigger_passed->Sumw2();
     hist_trigger_passed->Fill("HLT_Jet15U",0);
     hist_trigger_passed->Fill("HLT_Jet30U",0);
     hist_trigger_passed->Fill("HLT_Jet50U",0);

     hist_events->Sumw2();
     hist_events->Fill("Total Events",0);
     hist_events->Fill("PV Selection",0);
     hist_events->Fill("Trigger Selection",0);
     hist_events->Fill("Selected",0);
     hist_events->Fill("Selected with prescales",0);
     hist_events->Fill("Events Scenario anti-jet veto",0);
     hist_events->Fill("Events Scenario jet-veto",0);
     hist_events->Fill("Number of Jets",0);

  int counter_hlt = 0, counter_pv = 0, counter_jet = 0, counter_selected = 0, counter_entries = 0;
  int counter_jet15u = 0, counter_jet30u = 0, counter_jet50u = 0;
  double counter_weigth = 0.0, counter_check = 0.0;
  int nentries = 0;
  double pu_scale = 0.0;

// multiple triggers
  int number_triggers = 7;
  string trigger_list[number_triggers];
  int hlt_list[number_triggers];
//  double probs[number_triggers];
//  double prescales[number_triggers];
  bool hltPass(false);
//  bool trig_pass[number_triggers];
  trigger_list[0] = "HLT_Jet15U";
  trigger_list[1] = "HLT_Jet30U";
  trigger_list[2] = "HLT_Jet50U";
  trigger_list[3] = "HLT_DiJetAve15U_8E29";
  trigger_list[4] = "HLT_DiJetAve30U_8E29";
  trigger_list[5] = "HLT_DoubleJet15U_ForwardBackward";
  trigger_list[6] = "HLT_FwdJet20U";

//loop over the files
for (int z = 0; z < n_files; z++)
{
//open the file
  string file = data_in[z];
  cout<<z+1<<"/"<<n_files<<" -> "<<file<<endl;
  TFile *inf = TFile::Open( file.c_str() );

//define the tree branch
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);

//get the trigger information
  if (data_type == "DATA")
  {
  TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
  cout<<"Finding trigger mapping: "<<endl;
  //----------- loop over the X-axis labels -----------------
  for (int trg = 0; trg < number_triggers; trg++)
	{
	for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
		TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
		if (ss == trigger_list[trg]) {
			hlt_list[trg] = ibin;
			continue;
			}
		}
  	if (hlt_list[trg] == -1) {
		cout<<"The requested trigger ("<<trigger_list[trg]<<") is not found." << endl;
    		//break;
  		}
	else {
		cout<<trigger_list[trg]<<" --> "<<hlt_list[trg]<<endl;
		}
	}
  }

  prescale = 0.0;
  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  if (data_type == "MC_DET" or data_type == "MC_GEN")
  {
  prescale = data_lumi[z];
  cout<<"Cross-section per event : "<<prescale<<endl;
  }
  if (data_type == "DATA")
  {
  prescale_prov =  1.0 / data_lumi[z];
  cout<<"Luminosity for " << trigger_list[z] << " = " << data_lumi[z] << endl;
  cout<<"Cross-section per event : "<<prescale_prov<<endl;
  }

  int decade = 0;

//  if (test) { nentries = 100; } //reduced number of read entries, usefull for testing

   for (Int_t i=0;i<nentries;i++) {
     counter_entries++; 
     double progress = 10.0*i/(1.0*nentries);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
      cout<<k*10<<" % "<<endl;
     decade = k; 
     
     tr->GetEntry(i);

      //-------- check if the primary vertex is good ----  && Event->evtHdr().nVtxGood() == 1
      pv_pass = false;
      if (data_type == "MC_GEN")
	{
	if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)	{ pv_pass = true; }
	if (sel_mode == "allvertex") 					{ pv_pass = true; }
	}
      if (data_type == "MC_DET")
	{
	if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)						{ pv_pass = true; }
	if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) 					{ pv_pass = true; }
	if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)	{ pv_pass = true; }
	}
      if (data_type == "DATA")
	{
	if ((sel_mode == "up" || sel_mode == "down" || sel_mode == "allvertex") && Event->evtHdr().isPVgood() == 1)	{ pv_pass = true; }
	if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)    	  	{ pv_pass = true; }
	}

      if (pv_pass == true)
      {
        counter_pv++;
        if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1);} 

     if (data_type == "DATA")
     {
     // multiple triggers
     hltPass = false;
     prescale = 0.0;
	/*
	for (int trg = 0; trg < 3; trg++)
		{
		prescales[trg] = 1.0;
		trig_pass[trg] = false;
		if (hlt_list[trg] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[trg]) > 0) {
				hltPass = true;
				trig_pass[trg] = true;
				prescales[trg] = Event->preL1(hlt_list[trg]) * Event->preHLT(hlt_list[trg]);
				if (trg == 0) { counter_jet15u++; }
				if (trg == 1) { counter_jet30u++; }
				if (trg == 2) { counter_jet50u++; }
				}
			}
		}
	*/
	if (Event->fired(hlt_list[z]) > 0) { hltPass = true; }
     }

     if (hltPass and data_type == "DATA")
     {
     
     leading_pt = 0.0;

     for(unsigned int j=0; j<Event->nPFJets(); j++)
        {
	pt = Event->pfjet(j).ptCor();
	if (leading_pt < pt && Event->pfjet(j).tightID()) { leading_pt = pt; }
	}

	if (z == 0 and leading_pt > 35.0 and leading_pt < 70.0) { prescale = prescale_prov; }
	if (z == 1 and leading_pt > 70.0 and leading_pt < 110.0) { prescale = prescale_prov; }
	if (z == 2 and leading_pt > 110.0) { prescale = prescale_prov; }

     //if (leading_pt > 30.0 and leading_pt < 70.0 and trig_pass[0]) { prescale = prescales[0]; hist_trigger_passed->Fill(0.5); }
     //if (leading_pt > 70.0 and leading_pt < 110.0 and trig_pass[1]) { prescale = prescales[1]; hist_trigger_passed->Fill(1.5); }
     //if (leading_pt > 110.0 and trig_pass[2]) { prescale = prescales[2]; hist_trigger_passed->Fill(2.5); }
     
     //if (prescales[0] > 0.0 and prescales[1] == 0.0 and prescales[2] == 0.0) { prescale = prescales[0]; hist_trigger_passed->Fill(0.5); }
     //if (prescales[0] == 0.0 and prescales[1] > 0.0 and prescales[2] == 0.0) { prescale = prescales[1]; hist_trigger_passed->Fill(1.5); }
     //if (prescales[0] == 0.0 and prescales[1] == 0.0 and prescales[2] > 0.0) { prescale = prescales[2]; hist_trigger_passed->Fill(2.5); }
     //cout << "Event : " << i << endl;
     //cout << "Prescales : " << prescales[0] << " " << prescales[1] << " " << prescales[2] << endl;
     //cout << "Leading Jet pT : " << leading_pt << " prescale : " << prescale << endl;

// old combination

     //for(int l=0; l<number_triggers; l++)
//	{
//	probs[l] = 1.0;
//	if (prescales[l] > 0 and trig_pass[l]) { probs[l] = (1 - 1/prescales[l]); }
	//if (prescales[l] > 0 and trig_pass[l]) { cout<<"l = "<<l<<" | prob = "<<probs[l]<<" | prescale = "<<prescales[l]<<endl; }
//	}
//     prescale = 1.0 / ( 1 - (probs[0] * probs[1] * probs[2]));
     //prescale = 1.0 / ( 1 - (probs[0] * probs[1] * probs[2] * probs[3] * probs[4] * probs[5] * probs[6]));
     //if (prescale > 1) { cout<<"probs = ["<<probs[0]<<" ; "<<probs[1]<<" ; "<<probs[2]<<" ; "<<probs[3]<<" ; "<<probs[4]<<"]"<<endl;   
     //if (prescale > 1) {cout<<"Prescale = "<<prescale<<endl; }
     
     hist_leading_pt->Fill(leading_pt,prescale * vertex_factor);
     hist_leading_pt_fine->Fill(leading_pt,prescale * vertex_factor);

     }

     pt_forward = 0.0;
     eta_forward = 0.0;
     eta_forward2 = 0.0;
     phi_forward = 0.0;
     pt_central = 0.0;
     pt_central_check = 0.0;
     eta_central = 0.0;
     eta_central2 = 0.0;
     phi_central = 0.0;
     pt2_forward = 0.0;
     eta2_forward = 0.0;
     eta2_forward2 = 0.0;
     phi2_forward = 0.0;
     pt2_central = 0.0;
     eta2_central = 0.0;
     eta2_central2 = 0.0;
     phi2_central = 0.0;
     pass_gap = false;
     pt_total_gap = 0.0;
     eta_gap = -10.0;
     phi_gap = -10.0;
     chm_gap = -10.0;
     elm_gap = -10.0;
     multiplicity_inside = 0.0;
     pt_total_outside = 0.0;
     deta_out1 = -1.0;
     deta_out2 = -1.0;
     eta_outside = -10.0;
     phi_outside = -10.0;
     chm_outside = -10.0;
     elm_outside = -10.0;
     multiplicity_outside = 0.0;
     pt_leading_gap = 0.0;
     pt_leading_outside = 0.0;
     eta_star_inside = 0;
     chm_central = 0;
     chm_forward = 0;
     elm_central = 0;
     elm_forward = 0;
     pass_tight = false;
     pass_check = false;

     if (prescale > 0.0) { //or prescale > 0 hltPass_monitor
      counter_hlt++;
      //cout<<"Prescale = "<<prescale<<endl;

     if (data_type == "DATA" or data_type == "MC_DET") { number_of_jets = Event->nPFJets(); }
     if (data_type == "MC_GEN") { number_of_jets = Event->nGenJets(); }

    for(unsigned int j=0; j<number_of_jets; j++) {
	if (data_type == "DATA" or data_type == "MC_DET")
	{
        pass_tight = Event->pfjet(j).tightID();
	pt = Event->pfjet(j).ptCor();
	eta = Event->pfjet(j).eta();
	phi = Event->pfjet(j).phi();
	}
	if (data_type == "MC_GEN")
	{
	pass_tight = true;
	pt = Event->genjet(j).pt();
	eta = Event->genjet(j).eta();
	phi = Event->genjet(j).phi();
	}
	if (data_type == "MC_DET")
		{
		gen_pt = Event->genjet(j).pt();
		gen_eta = Event->genjet(j).eta();
		gen_phi = Event->genjet(j).phi();
		pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,test);
		}
	if (data_type == "DATA")
	{
	unc = Event->pfjet(j).unc();
	pt_up = pt * (1+unc);
	pt_down = pt * (1-unc);
	if (sel_mode == "up") { pt = pt_up;}
	if (sel_mode == "down") { pt = pt_down;}
	chm = Event->pfjet(j).chm();
	elm = Event->pfjet(j).elm();
	hist_chm_all->Fill(chm,prescale * vertex_factor);
	hist_elm_all->Fill(elm,prescale * vertex_factor);
	}
     hist_pt_all->Fill(pt,prescale * vertex_factor);
     hist_eta_all->Fill(eta,prescale * vertex_factor);
     hist_phi_all->Fill(phi,prescale * vertex_factor);
     if (pt >= pt_min && pass_tight) {
     counter_jet++;
     if (eta <= 4.7 && eta >= -4.7 && pt > 30) { pass_check = true; }
     if (eta <= 2.5 && eta >= -2.5 && pt > pt_central_check) {  pt_central_check = pt; }
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central)
     { pt_central = pt; eta_central = eta; phi_central = phi; chm_central = chm; elm_central = elm; }
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_min)
     { hist_inclusive_central_pt->Fill(pt,prescale * vertex_factor); }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward)
     { pt_forward = pt; eta_forward = eta; phi_forward = phi; chm_forward = chm; elm_forward = elm; }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_min)
     { hist_inclusive_forward_pt->Fill(pt,prescale * vertex_factor); }
     }
     if (eta <= 2.8 && eta >= -2.8 && pt > pt2_central && pt < pt_central)
     { pt2_central = pt; eta2_central = eta; phi2_central = phi;}
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt < pt_forward && pt > pt2_forward)
     { pt2_forward = pt; eta2_forward = eta; phi2_forward = phi;}
     }
     
     if (pt_central_check > pt_min) { hist_leading_central_pt_check->Fill(pt_central_check,prescale * vertex_factor); }
     if (pt_central > pt_min) { hist_inclusive_leading_central_pt->Fill(pt_central,prescale * vertex_factor); }
     if (pt_forward > pt_min) { hist_inclusive_leading_forward_pt->Fill(pt_forward,prescale * vertex_factor); }
     if (pass_check) { counter_check = counter_check + prescale * vertex_factor; }     

     if (pt_forward > pt_min && pt_central > pt_min)
     {
     
     hist_central_pt_rel->Fill(pt_central,pt2_central,prescale * vertex_factor);
     hist_forward_pt_rel->Fill(pt_forward,pt2_forward,prescale * vertex_factor);
     if (pt_central > 20.0 and pt2_central > 0.0)
     {
     delta_phi2 = calc_delta_phi(phi_central, phi2_central);
     hist_delta_phi_central_rel->Fill(delta_phi2,prescale * vertex_factor);
     delta_pt = pt_central - pt2_central;
     hist_delta_pt_central_rel->Fill(delta_pt,prescale * vertex_factor);
     if (delta_pt < 5.0) { hist_delta_phi_central_rel_small->Fill(delta_phi2,prescale * vertex_factor); }
     if (delta_pt > 5.0 && delta_pt < 15.0) { hist_delta_phi_central_rel_medium->Fill(delta_phi2,prescale * vertex_factor); }
     if (delta_pt > 15.0) { hist_delta_phi_central_rel_large->Fill(delta_phi2,prescale * vertex_factor); }
     }
     if (pt_forward > 20.0 and pt2_forward > 0.0)
     {
     delta_phi2 = calc_delta_phi(phi_forward, phi2_forward);
     hist_delta_phi_forward_rel->Fill(delta_phi2,prescale * vertex_factor);
     delta_pt = pt_forward - pt2_forward;
     hist_delta_pt_forward_rel->Fill(delta_pt,prescale * vertex_factor);
     if (delta_pt < 5.0) { hist_delta_phi_forward_rel_small->Fill(delta_phi2,prescale * vertex_factor); }
     if (delta_pt > 5.0 && delta_pt < 15.0) { hist_delta_phi_forward_rel_medium->Fill(delta_phi2,prescale * vertex_factor); }
     if (delta_pt > 15.0) { hist_delta_phi_forward_rel_large->Fill(delta_phi2,prescale * vertex_factor); }
     }
     
     hist_leading_central_pt->Fill(pt_central,prescale * vertex_factor);
     hist_leading_forward_pt->Fill(pt_forward,prescale * vertex_factor);
     hist_leading_central_pt_fine->Fill(pt_central,prescale * vertex_factor);
     hist_leading_forward_pt_fine->Fill(pt_forward,prescale * vertex_factor);
     hist_leading_central_phi->Fill(phi_central,prescale * vertex_factor);
     hist_leading_forward_phi->Fill(phi_forward,prescale * vertex_factor);
     eta_central2 = eta_central;     
     if (eta_central2 < 0) { eta_central2 = -eta_central2; }
     eta_forward2 = eta_forward;     
     if (eta_forward2 < 0) { eta_forward2 = -eta_forward2; }
     hist_leading_central_eta->Fill(eta_central,prescale * vertex_factor);
     hist_leading_forward_eta->Fill(eta_forward,prescale * vertex_factor);
     hist_leading_eta->Fill(eta_central,prescale * vertex_factor);
     hist_leading_eta->Fill(eta_forward,prescale * vertex_factor);
     hist_multiplicity->Fill(Event->nPFJets(),prescale * vertex_factor);
     hist_leading_central_chm->Fill(chm_central,prescale * vertex_factor);
     hist_leading_forward_chm->Fill(chm_forward,prescale * vertex_factor);
     hist_leading_central_elm->Fill(elm_central,prescale * vertex_factor);
     hist_leading_forward_elm->Fill(elm_forward,prescale * vertex_factor);
     hist_vertex_selected->Fill(Event->evtHdr().nVtxGood(),prescale * vertex_factor);
     hist_pvz_selected->Fill(Event->evtHdr().PVz(),prescale * vertex_factor);
     if (data_type == "MC_GEN" or data_type == "MC_DET") { hist_pu_selected->Fill(Event->evtHdr().pu(),prescale * vertex_factor); }
     counter_selected++;
     counter_weigth = counter_weigth + prescale;
     delta_eta = eta_forward - eta_central;
     if (delta_eta < 0) { delta_eta = -delta_eta; }
     delta_phi = calc_delta_phi(phi_forward, phi_central);
     hist_delta_phi->Fill(delta_phi,prescale * vertex_factor);
     hist_delta_phi_fine->Fill(delta_phi,prescale * vertex_factor);
     hist_delta_eta->Fill(delta_eta,prescale * vertex_factor);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { hist_delta_phi_deta1->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { hist_delta_phi_deta2->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { hist_delta_phi_deta3->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { hist_delta_phi_deta4->Fill(delta_phi,prescale * vertex_factor); }    

     for(unsigned int j=0;j<number_of_jets;j++) {
     if (data_type == "DATA" or data_type == "MC_DET")
	{
	eta = Event->pfjet(j).eta();
	phi = Event->pfjet(j).phi();
	pt = Event->pfjet(j).ptCor();
	}
     if (data_type == "MC_DET")
	{
	gen_eta = Event->genjet(j).eta();
	gen_phi = Event->genjet(j).phi();
	gen_pt = Event->genjet(j).pt();
	pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,test);
	}
     if (data_type == "MC_GEN") { pt = Event->genjet(j).pt(); }
	if (data_type == "DATA")
	{
	unc = Event->pfjet(j).unc();
	pt_up = pt * (1+unc);
	pt_down = pt * (1-unc);
	if (sel_mode == "up") { pt = pt_up;}
	if (sel_mode == "down") { pt = pt_down;}
	}
     if (data_type == "DATA" or data_type == "MC_DET") { pass_tight = Event->pfjet(j).tightID(); }
     if (data_type == "MC_GEN") { pass_tight = true; }

     if (pt >= 5 && pass_tight) {
     if (data_type == "MC_GEN")
     {
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
     } 
	if (data_type == "DATA")
	{
	chm = Event->pfjet(j).chm();
	elm = Event->pfjet(j).elm();
	}
     if (eta_central > eta_forward && eta > eta_forward && eta < eta_central ) {
        if (pt > pt_leading_gap) {pt_leading_gap = pt; eta_gap = eta; phi_gap = phi; chm_gap = chm; elm_gap = elm; }
        pt_total_gap = pt_total_gap + pt;
	multiplicity_inside = multiplicity_inside + 1;
        }
     if (eta_central < eta_forward && eta < eta_forward && eta > eta_central ) {
        if (pt > pt_leading_gap) {pt_leading_gap = pt; eta_gap = eta; phi_gap = phi; chm_gap = chm; elm_gap = elm; }
        pt_total_gap = pt_total_gap + pt;
	multiplicity_inside = multiplicity_inside + 1;
        }
     if (eta_central > eta_forward && (eta < eta_forward || eta > eta_central) ) {
        if (pt > pt_leading_outside) {pt_leading_outside = pt; eta_outside = eta; phi_outside = phi; chm_outside = chm; elm_outside = elm; }
        pt_total_outside = pt_total_outside + pt;
	multiplicity_outside = multiplicity_outside + 1;
        }
     if (eta_central < eta_forward && (eta > eta_forward || eta < eta_central) ) {
        if (pt > pt_leading_outside) {pt_leading_outside = pt; eta_outside = eta; phi_outside = phi; chm_outside = chm; elm_outside = elm; }
        pt_total_outside = pt_total_outside + pt;
	multiplicity_outside = multiplicity_outside + 1;
        }
     }
     }
     if (pt_leading_gap > gap_req )
     {
     pass_gap = true;
   //  cout<<i<<" [1]pt_gap = "<<pt_leading_gap<<endl;
     eta_star_inside = eta_gap - (eta_forward + eta_central)/2;
     hist_leading_pt_inside_gap->Fill(pt_leading_gap,prescale * vertex_factor);
     hist_total_pt_inside_gap->Fill(pt_total_gap,prescale * vertex_factor);
     hist_leading_eta_inside_gap->Fill(eta_gap,prescale * vertex_factor);
     hist_leading_phi_inside_gap->Fill(phi_gap,prescale * vertex_factor);
     hist_leading_chm_inside_gap->Fill(chm_gap,prescale * vertex_factor);
     hist_leading_elm_inside_gap->Fill(elm_gap,prescale * vertex_factor);
     hist_multiplicity_inside_gap->Fill(multiplicity_inside,prescale * vertex_factor);
     hist_leading_eta_star_inside_gap->Fill(eta_star_inside,prescale * vertex_factor);
     }
     if (pt_leading_outside > gap_req )
     {
     hist_leading_pt_outside_gap->Fill(pt_leading_outside,prescale * vertex_factor);
     hist_total_pt_outside_gap->Fill(pt_total_outside,prescale * vertex_factor);
     hist_leading_eta_outside_gap->Fill(eta_outside,prescale * vertex_factor);
     hist_leading_phi_outside_gap->Fill(phi_outside,prescale * vertex_factor);
     hist_leading_chm_outside_gap->Fill(chm_outside,prescale * vertex_factor);
     hist_leading_elm_outside_gap->Fill(elm_outside,prescale * vertex_factor);
     hist_multiplicity_outside_gap->Fill(multiplicity_outside,prescale * vertex_factor);
     deta_out1 = eta_central - eta_outside;
     if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     deta_out2 = eta_forward - eta_outside;
     if (deta_out2 < 0) { deta_out2 = -deta_out2; }
     if (deta_out1 < deta_out2) { hist_delta_eta_outside_gap->Fill(deta_out1,prescale * vertex_factor); }
     if (deta_out2 < deta_out1) { hist_delta_eta_outside_gap->Fill(deta_out2,prescale * vertex_factor); }
    // if ((deta_out1 < deta_out2 && deta_out1 > 7.5) || (deta_out2 < deta_out1 && deta_out2 > 7.5)) { cout<<deta_out1<<" - "<<deta_out2<<endl; }
     }
     if (pass_gap == false)
     {
     selected_gap = selected_gap + 1;
     hist_delta_phi_gap->Fill(delta_phi,prescale * vertex_factor);
     hist_delta_eta_gap->Fill(delta_eta,prescale * vertex_factor);
     hist_leading_central_pt_gap->Fill(pt_central,prescale * vertex_factor);
     hist_leading_forward_pt_gap->Fill(pt_forward,prescale * vertex_factor);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { hist_delta_phi_deta1_gap->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { hist_delta_phi_deta2_gap->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { hist_delta_phi_deta3_gap->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { hist_delta_phi_deta4_gap->Fill(delta_phi,prescale * vertex_factor); }
     }
     if (pass_gap == true)
     {
     selected_nogap = selected_nogap + 1;
   //  cout<<selected_nogap<<" [2] pt_gap = "<<pt_leading_gap<<endl;
     hist_delta_phi_nogap->Fill(delta_phi,prescale * vertex_factor);
     hist_delta_eta_nogap->Fill(delta_eta,prescale * vertex_factor);
     hist_leading_central_pt_nogap->Fill(pt_central,prescale * vertex_factor);
     hist_leading_forward_pt_nogap->Fill(pt_forward,prescale * vertex_factor);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { hist_delta_phi_deta1_nogap->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { hist_delta_phi_deta2_nogap->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { hist_delta_phi_deta3_nogap->Fill(delta_phi,prescale * vertex_factor); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { hist_delta_phi_deta4_nogap->Fill(delta_phi,prescale * vertex_factor); }
     }
     
     }
  } 


  }  } 
  
    cout<<"Events read:                      "<<nentries<<endl; 
 } 
  
  cout<<endl<<endl;
  cout<<"Events read:                      "<<counter_entries<<endl;
  cout<<"Events after the PV cut:          "<<counter_pv<<endl;
  if (data_type == "DATA")
	{
	cout<<"Events after the trigger cut:     "<<counter_hlt<<endl;
	cout<<"Events which fired HLT_Jet15U:    "<<counter_jet15u<<endl;
	cout<<"Events which fired HLT_Jet30U:    "<<counter_jet30u<<endl;
	cout<<"Events which fired HLT_Jet50U:    "<<counter_jet50u<<endl;
	}
  cout<<"Events Main Selection:            "<<counter_selected<<endl;
  cout<<"Events Scenario anti-jet veto:    "<<selected_nogap<<endl;
  cout<<"Events Scenario jet-veto:         "<<selected_gap<<endl;
  cout<<"Events weight count:              "<<counter_weigth<<endl;
  cout<<"Total number of Jets:             "<<counter_jet<<endl;
  cout<<"Cross-section Check:              "<<counter_check<<endl;

     pu_scale = (double) nentries/ (double) counter_pv;
     cout<<endl;
     cout<<"Pile up correction =  "<<pu_scale<<endl;
     cout<<endl;

//fill the trigger_passed histogram
     hist_trigger_passed->SetBinContent(1,counter_jet15u);
     hist_trigger_passed->SetBinContent(2,counter_jet30u);
     hist_trigger_passed->SetBinContent(3,counter_jet50u);

//fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_hlt);
     hist_events->SetBinContent(4,counter_selected);
     hist_events->SetBinContent(5,counter_weigth);
     hist_events->SetBinContent(6,selected_nogap);
     hist_events->SetBinContent(7,selected_gap);
     hist_events->SetBinContent(8,counter_jet);

//normalize the histograms
     normalize_histogram(hist_delta_phi, "Delta Phi");
     normalize_histogram(hist_delta_phi_fine, "Delta Phi Fine");
     normalize_histogram(hist_delta_phi_gap, "Delta Phi Gap");
     normalize_histogram(hist_delta_phi_nogap, "Delta Phi noGap");
     normalize_histogram(hist_delta_phi_deta1, "Delta Phi Deta1");
     normalize_histogram(hist_delta_phi_deta2, "Delta Phi Deta2");
     normalize_histogram(hist_delta_phi_deta3, "Delta Phi Deta3");
     normalize_histogram(hist_delta_phi_deta4, "Delta Phi Deta4");
     normalize_histogram(hist_delta_phi_deta1_gap, "Delta Phi Deta1 Gap");
     normalize_histogram(hist_delta_phi_deta2_gap, "Delta Phi Deta2 Gap");
     normalize_histogram(hist_delta_phi_deta3_gap, "Delta Phi Deta3 Gap");
     normalize_histogram(hist_delta_phi_deta4_gap, "Delta Phi Deta4 Gap");
     normalize_histogram(hist_delta_phi_deta1_nogap, "Delta Phi Deta1 noGap");
     normalize_histogram(hist_delta_phi_deta2_nogap, "Delta Phi Deta2 noGap");
     normalize_histogram(hist_delta_phi_deta3_nogap, "Delta Phi Deta3 noGap");
     normalize_histogram(hist_delta_phi_deta4_nogap, "Delta Phi Deta4 noGap");
     normalize_histogram(hist_delta_eta, "Delta Eta");
     normalize_histogram(hist_delta_eta_gap, "Delta Eta Gap");
     normalize_histogram(hist_delta_eta_nogap, "Delta Eta noGap");
     normalize_histogram(hist_total_pt_inside_gap, "Total pT Inside Gap");
     normalize_histogram(hist_leading_pt_inside_gap, "Leading pT Inside Gap");
     normalize_histogram(hist_leading_eta_inside_gap, "Leading Eta Inside Gap");
     normalize_histogram(hist_leading_phi_inside_gap, "Leading Phi Inside Gap");
     normalize_histogram(hist_leading_chm_inside_gap, "Leading Charged Hadron Multiplicity Inside Gap");
     normalize_histogram(hist_leading_elm_inside_gap, "Leading Electron Multiplicity Inside Gap");
     normalize_histogram(hist_multiplicity_inside_gap, "Jet Multiplicity Inside Gap");
     normalize_histogram(hist_total_pt_outside_gap, "Total pT Outside Gap");
     normalize_histogram(hist_delta_eta_outside_gap, "Delta Eta Outside Gap");
     normalize_histogram(hist_leading_pt_outside_gap, "Leading pT Outside Gap");
     normalize_histogram(hist_leading_eta_outside_gap, "Leading Eta Outside Gap");
     normalize_histogram(hist_leading_phi_outside_gap, "Leading Phi Outside Gap");
     normalize_histogram(hist_leading_chm_outside_gap, "Leading Charged Hadron Multiplicity Outside Gap");
     normalize_histogram(hist_leading_elm_outside_gap, "Leading Electron Multiplicity Outside Gap");
     normalize_histogram(hist_multiplicity_outside_gap, "Jet Multiplicity Outside Gap");
     normalize_histogram(hist_leading_eta_star_inside_gap, "Leading Eta* Outside Gap");
     normalize_histogram(hist_leading_pt, "Leading pT");
     normalize_histogram(hist_leading_pt_fine, "Leading pT Fine");
     normalize_histogram(hist_leading_central_pt_fine, "Leading Central pT Fine");
     normalize_histogram(hist_leading_forward_pt_fine, "Leading Forward pT Fine");
     normalize_histogram(hist_leading_central_pt, "Leading Central pT");
     normalize_histogram(hist_leading_central_pt_gap, "Leading Central pT Gap");
     normalize_histogram(hist_leading_central_pt_nogap, "Leading Central pT noGap");
     normalize_histogram(hist_leading_central_eta, "Leading Central Eta");
     normalize_histogram(hist_leading_central_phi, "Leading Central Phi");
     normalize_histogram(hist_leading_central_chm, "Leading Central Charged Hadron Multiplicity");
     normalize_histogram(hist_leading_central_elm, "Leading Central Electron Multiplicity");
     normalize_histogram(hist_leading_forward_pt, "Leading Forward pT");
     normalize_histogram(hist_leading_forward_pt_gap, "Leading Forward pT Gap");
     normalize_histogram(hist_leading_forward_pt_nogap, "Leading Forward pT noGap");
     normalize_histogram(hist_leading_forward_eta, "Leading Forward Eta");
     normalize_histogram(hist_leading_forward_phi, "Leading Forward Phi");
     normalize_histogram(hist_leading_forward_chm, "Leading Forward Charged Hadron Multiplicity");
     normalize_histogram(hist_leading_forward_elm, "Leading Forward Electron Multiplicity");
     normalize_histogram(hist_leading_eta, "Leading Eta");
     normalize_histogram(hist_vertex_selected, "Vertex Selected");
     normalize_histogram(hist_pvz_selected, "z-position Selected");
     normalize_histogram(hist_pu_selected, "Pileup Selected");
     normalize_histogram(hist_multiplicity, "Jet Multiplicity");
     normalize_histogram(hist_chm_all, "Charged Hadron Multiplicity All");
     normalize_histogram(hist_elm_all, "Electron Multiplicity All");
     normalize_histogram(hist_pt_all, "pT All");
     normalize_histogram(hist_eta_all, "Eta All");
     normalize_histogram(hist_phi_all, "Phi All");
     normalize_histogram(hist_inclusive_leading_central_pt, "Inclusive Leading Central pT");
     normalize_histogram(hist_inclusive_leading_forward_pt, "Inclusive Leading Forward pT");
     normalize_histogram(hist_inclusive_central_pt, "Inclusive Central pT");
     normalize_histogram(hist_inclusive_forward_pt, "Inclusive Forward pT");
     normalize_histogram(hist_delta_phi_central_rel, "Delta Phi Central Relative");
     normalize_histogram(hist_delta_phi_forward_rel, "Delta Phi Forward Relative");
     normalize_histogram(hist_delta_pt_central_rel, "Delta pT Central Relative");
     normalize_histogram(hist_delta_pt_forward_rel, "Delta pT Forward Relative");
     normalize_histogram(hist_delta_phi_central_rel_small, "Delta Phi Central Relative Small");
     normalize_histogram(hist_delta_phi_forward_rel_small, "Delta Phi Forward Relative Small");
     normalize_histogram(hist_delta_phi_central_rel_medium, "Delta Phi Central Relative Medium");
     normalize_histogram(hist_delta_phi_forward_rel_medium, "Delta Phi Forward Relative Medium");
     normalize_histogram(hist_delta_phi_central_rel_large, "Delta Phi Central Relative Large");
     normalize_histogram(hist_delta_phi_forward_rel_large, "Delta Phi Forward Relative Large");
    
     //Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     //save histograms on the file
     hist_delta_phi->Write();
     hist_delta_phi_fine->Write();
     hist_delta_phi_gap->Write();
     hist_delta_phi_nogap->Write();
     hist_delta_phi_deta1->Write();
     hist_delta_phi_deta2->Write();
     hist_delta_phi_deta3->Write();
     hist_delta_phi_deta4->Write();
     hist_delta_phi_deta1_gap->Write();
     hist_delta_phi_deta2_gap->Write();
     hist_delta_phi_deta3_gap->Write();
     hist_delta_phi_deta4_gap->Write();
     hist_delta_phi_deta1_nogap->Write();
     hist_delta_phi_deta2_nogap->Write();
     hist_delta_phi_deta3_nogap->Write();
     hist_delta_phi_deta4_nogap->Write();
     hist_delta_eta->Write();
     hist_delta_eta_gap->Write();
     hist_delta_eta_nogap->Write();
     hist_total_pt_inside_gap->Write();
     hist_leading_pt_inside_gap->Write();
     hist_leading_eta_inside_gap->Write();
     hist_leading_phi_inside_gap->Write();
     hist_multiplicity_inside_gap->Write();
     hist_total_pt_outside_gap->Write();
     hist_delta_eta_outside_gap->Write();
     hist_leading_eta_outside_gap->Write();
     hist_leading_phi_outside_gap->Write();
     hist_leading_pt_outside_gap->Write();
     hist_multiplicity_outside_gap->Write();
     hist_leading_eta_star_inside_gap->Write();
     hist_leading_central_pt_fine->Write();
     hist_leading_forward_pt_fine->Write();
     hist_leading_central_pt->Write();
     hist_leading_central_pt_gap->Write();
     hist_leading_central_pt_nogap->Write();
     hist_leading_central_eta->Write();
     hist_leading_central_phi->Write();
     hist_leading_forward_pt->Write();
     hist_leading_forward_pt_gap->Write();
     hist_leading_forward_pt_nogap->Write();
     hist_leading_forward_eta->Write();
     hist_leading_forward_phi->Write();
     hist_leading_eta->Write();
     hist_vertex_selected->Write();
     hist_pvz_selected->Write();
     hist_multiplicity->Write();
     hist_pt_all->Write();
     hist_eta_all->Write();
     hist_phi_all->Write();
     hist_inclusive_leading_central_pt->Write();
     hist_inclusive_leading_forward_pt->Write();
     hist_inclusive_central_pt->Write();
     hist_inclusive_forward_pt->Write();
     hist_leading_central_pt_check->Write();
     hist_central_pt_rel->Write();
     hist_forward_pt_rel->Write();
     hist_delta_phi_central_rel->Write();
     hist_delta_phi_forward_rel->Write();
     hist_delta_pt_central_rel->Write();
     hist_delta_pt_forward_rel->Write();
     hist_delta_phi_central_rel_small->Write();
     hist_delta_phi_forward_rel_small->Write();
     hist_delta_phi_central_rel_medium->Write();
     hist_delta_phi_forward_rel_medium->Write();
     hist_delta_phi_central_rel_large->Write();
     hist_delta_phi_forward_rel_large->Write();
     hist_events->Write();

	if (data_type == "DATA")
	{
        hist_leading_pt->Write();
        hist_leading_pt_fine->Write();
	hist_leading_chm_inside_gap->Write();
	hist_leading_elm_inside_gap->Write();
	hist_leading_elm_outside_gap->Write();
	hist_leading_chm_outside_gap->Write();
	hist_leading_central_chm->Write();
	hist_leading_central_elm->Write();
	hist_leading_forward_chm->Write();
	hist_leading_forward_elm->Write();
	hist_chm_all->Write();
	hist_elm_all->Write();
        hist_trigger_passed->Write();
	}

	if (data_type == "MC_GEN" or data_type == "MC_DET")
	{
	hist_pu_selected->Write();
	}

     //close the output file
     data_output->Close();

     //delete the histograms to avoid memory leak
     delete(hist_delta_phi);
     delete(hist_delta_phi_fine);
     delete(hist_delta_phi_gap);
     delete(hist_delta_phi_nogap);
     delete(hist_delta_phi_deta1);
     delete(hist_delta_phi_deta2);
     delete(hist_delta_phi_deta3);
     delete(hist_delta_phi_deta4);
     delete(hist_delta_phi_deta1_gap);
     delete(hist_delta_phi_deta2_gap);
     delete(hist_delta_phi_deta3_gap);
     delete(hist_delta_phi_deta4_gap);
     delete(hist_delta_phi_deta1_nogap);
     delete(hist_delta_phi_deta2_nogap);
     delete(hist_delta_phi_deta3_nogap);
     delete(hist_delta_phi_deta4_nogap);
     delete(hist_delta_eta);
     delete(hist_delta_eta_gap);
     delete(hist_delta_eta_nogap);
     delete(hist_total_pt_inside_gap);
     delete(hist_leading_pt_inside_gap);
     delete(hist_leading_eta_inside_gap);
     delete(hist_leading_phi_inside_gap);
     delete(hist_leading_chm_inside_gap);
     delete(hist_leading_elm_inside_gap);
     delete(hist_multiplicity_inside_gap);
     delete(hist_total_pt_outside_gap);
     delete(hist_delta_eta_outside_gap);
     delete(hist_leading_pt_outside_gap);
     delete(hist_leading_eta_outside_gap);
     delete(hist_leading_phi_outside_gap);
     delete(hist_leading_chm_outside_gap);
     delete(hist_leading_elm_outside_gap);
     delete(hist_multiplicity_outside_gap);
     delete(hist_leading_eta_star_inside_gap);
     delete(hist_leading_pt);
     delete(hist_leading_pt_fine);
     delete(hist_leading_central_pt_fine);
     delete(hist_leading_forward_pt_fine);
     delete(hist_leading_central_pt);
     delete(hist_leading_central_pt_gap);
     delete(hist_leading_central_pt_nogap);
     delete(hist_leading_central_eta);
     delete(hist_leading_central_phi);
     delete(hist_leading_central_chm);
     delete(hist_leading_central_elm);
     delete(hist_leading_forward_pt);
     delete(hist_leading_forward_pt_gap);
     delete(hist_leading_forward_pt_nogap);
     delete(hist_leading_forward_eta);
     delete(hist_leading_forward_phi);
     delete(hist_leading_forward_chm);
     delete(hist_leading_forward_elm);
     delete(hist_leading_eta);
     delete(hist_vertex_selected);
     delete(hist_pvz_selected);
     delete(hist_pu_selected);
     delete(hist_multiplicity);
     delete(hist_chm_all);
     delete(hist_elm_all);
     delete(hist_pt_all);
     delete(hist_eta_all);
     delete(hist_phi_all);
     delete(hist_inclusive_leading_central_pt);
     delete(hist_inclusive_leading_forward_pt);
     delete(hist_inclusive_central_pt);
     delete(hist_inclusive_forward_pt);
     delete(hist_leading_central_pt_check);
     delete(hist_central_pt_rel);
     delete(hist_forward_pt_rel);
     delete(hist_delta_phi_central_rel);
     delete(hist_delta_phi_forward_rel);
     delete(hist_delta_pt_central_rel);
     delete(hist_delta_pt_forward_rel);
     delete(hist_delta_phi_central_rel_small);
     delete(hist_delta_phi_forward_rel_small);
     delete(hist_delta_phi_central_rel_medium);
     delete(hist_delta_phi_forward_rel_medium);
     delete(hist_delta_phi_central_rel_large);
     delete(hist_delta_phi_forward_rel_large);
     delete(hist_trigger_passed);
     delete(hist_events);
}
