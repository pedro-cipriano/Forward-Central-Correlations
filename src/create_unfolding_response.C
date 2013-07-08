// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 11 Fev 2012
//

#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <string>

#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldResponse.h"

using namespace std;

#include "common_methods.h"

void create_unfolding_response(string *data_in, string data_out, double *data_lumi, int n_files,  string sel_mode = "allvertex", string vertex_weights = "", TString vertex_sufix = "_v0", bool detail = false, bool test = false)
{ 

//output the configuration
   if (detail) { cout<<"Create Unfolding Response Configuration"<<endl; }
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
    if (vertex_weights != "")
    {
    if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
    TFile *vertex_file = new TFile( vertex_weights.c_str() );
    vertex_file->GetObject("vertex_weights"+vertex_sufix,vertex_hist);
    if (vertex_hist != 0) { vertex_weighted = true; }
    }

   //main cuts
double pt_min = 35.0;
double gap_req = 20.0;

double pt, eta, phi, gen_pt, gen_eta, gen_phi;
double prescale, vertex_factor = 1.0; //val;

double pt_central_gen, eta_central_gen, eta_central2_gen, phi_central_gen;
double pt_forward_gen, eta_forward_gen, eta_forward2_gen, phi_forward_gen;
double pt2_central_gen, eta2_central_gen, eta2_central2_gen, phi2_central_gen;
double pt2_forward_gen, eta2_forward_gen, eta2_forward2_gen, phi2_forward_gen;
double delta_eta_gen, delta_phi_gen, delta_phi2_gen, delta_pt_gen, deta_out1_gen, deta_out2_gen, deta_out_gen;
double pt_total_gap_gen, pt_leading_gap_gen, eta_gap_gen, eta_star_inside_gen, phi_gap_gen;
double pt_total_outside_gen, pt_leading_outside_gen, eta_outside_gen, phi_outside_gen;

double pt_central_det, eta_central_det, eta_central2_det, phi_central_det;
double pt_forward_det, eta_forward_det, eta_forward2_det, phi_forward_det;
double pt2_central_det, eta2_central_det, eta2_central2_det, phi2_central_det;
double pt2_forward_det, eta2_forward_det, eta2_forward2_det, phi2_forward_det;
double delta_eta_det, delta_phi_det, delta_phi2_det, delta_pt_det, deta_out1_det, deta_out2_det, deta_out_det;
double pt_total_gap_det, pt_leading_gap_det, eta_gap_det, eta_star_inside_det, phi_gap_det;
double pt_total_outside_det, pt_leading_outside_det, eta_outside_det, phi_outside_det;

bool pass_gen, pass_pv_gen, pass_det, pass_pv_det;
bool pass_gap_gen, pass_nogap_gen, pass_gap_det, pass_nogap_det;
bool pass_out_gen, pass_out_det;
bool pass_deta1_gen, pass_deta2_gen, pass_deta3_gen, pass_deta4_gen;
bool pass_deta1_det, pass_deta2_det, pass_deta3_det, pass_deta4_det;

int counter_pv = 0, counter_selected_gen = 0, counter_selected_det = 0, counter_entries = 0;
double counter_weigth_gen = 0.0, counter_weigth_det = 0.0;
int counter_response = 0, counter_miss = 0, counter_truth = 0, counter_truth_missed = 0;
int counter_truth_stayed = 0, counter_truth_fakes = 0, counter_truth_fakes_slightly = 0;
int counter_gap_gen = 0, counter_nogap_gen = 0, counter_gap_det = 0, counter_nogap_det = 0;
int nentries = 0;

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};
int dphi_nbinsb = 14;
double dphi_binsb[15] = {0.0, 0.225, 0.45, 0.675, 0.9, 1.125, 1.35, 1.575, 1.8, 2.025, 2.25, 2.475, 2.7, 2.925, 3.15};
int dphi_nbins_true = 8;
double dphi_bins_true[9] = {-1.0, 0.0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};
int dphi_nbins_meas = 15;
double dphi_bins_meas[16] = {-1.0, 0.0, 0.225, 0.45, 0.675, 0.9, 1.125, 1.35, 1.575, 1.8, 2.025, 2.25, 2.475, 2.7, 2.925, 3.15};


int in_nbins = 9;
double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
int in_nbinsb = 18;
double in_binsb[19] = {20, 24, 27, 31, 35, 40, 45, 51, 57, 65, 72, 80, 90, 105, 120, 130, 150, 170, 200};
int in_nbins_true = 10;
double in_bins_true[11] = {0, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
int in_nbins_meas = 19;
double in_bins_meas[20] = {0, 20, 24, 27, 31, 35, 40, 45, 51, 57, 65, 72, 80, 90, 105, 120, 130, 150, 170, 200};

int out_nbins = 9;
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
int out_nbinsb = 18;
double out_binsb[19] = {20, 24, 27, 31, 35, 40, 45, 51, 57, 65, 72, 80, 90, 105, 120, 130, 150, 170, 200};
int out_nbins_true = 10;
double out_bins_true[11] = {0, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
int out_nbins_meas = 19;
double out_bins_meas[20] = {0, 20, 24, 27, 31, 35, 40, 45, 51, 57, 65, 72, 80, 90, 105, 120, 130, 150, 170, 200};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};
int etastar_nbinsb = 24;
double etastar_binsb[25] = {-3.6, -3.0, -2.5, -2.3, -2.0, -1.7, -1.5, -1.3, -1.0, -0.8, -0.5, -0.3, 0.0, 0.3, 0.5, 0.8, 1.0, 1.3, 1.5, 1.8, 2.0, 2.3, 2.5, 3.0, 3.6};
int etastar_nbins_true = 13;
double etastar_bins_true[14] = {-10,-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};
int etastar_nbins_meas = 25;
double etastar_bins_meas[26] = {-10,-3.6, -3.0, -2.5, -2.3, -2.0, -1.7, -1.5, -1.3, -1.0, -0.8, -0.5, -0.3, 0.0, 0.3, 0.5, 0.8, 1.0, 1.3, 1.5, 1.8, 2.0, 2.3, 2.5, 3.0, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};
int deta_out_nbinsb = 12;
double deta_out_binsb[13] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7.5};
int deta_out_nbins_true = 7;
double deta_out_bins_true[8] = {-1.0, 0, 1, 2, 3, 4, 5, 7.5};
int deta_out_nbins_meas = 13;
double deta_out_bins_meas[14] = {-1.0, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7.5};

TH1D *dummy_delta_phi;
TH1D *dummy_delta_phib;
TH1D *dummy_delta_phi_true;
TH1D *dummy_delta_phi_meas;
dummy_delta_phi =  new TH1D("dummy_delta_phi","#Delta#phi;#Delta#phi [rad];Events", dphi_nbins, dphi_bins);
dummy_delta_phib =  new TH1D("dummy_delta_phib","#Delta#phi;#Delta#phi [rad];Events", dphi_nbinsb, dphi_binsb);
dummy_delta_phi_true =  new TH1D("dummy_delta_phi_true","#Delta#phi;#Delta#phi;Events", dphi_nbins_true, dphi_bins_true);
dummy_delta_phi_meas =  new TH1D("dummy_delta_phi_meas","#Delta#phi;#Delta#phi [rad];Events", dphi_nbins_meas, dphi_bins_meas);

TH1D *dummy_pt_in;
TH1D *dummy_pt_inb;
TH1D *dummy_pt_in_true;
TH1D *dummy_pt_in_meas;
dummy_pt_in =  new TH1D("dummy_pt_in","p_T^{in};p_T^{in} [GeV];Events", in_nbins, in_bins);
dummy_pt_inb =  new TH1D("dummy_pt_inb","p_T^{in};p_T^{in} [GeV];Events", in_nbinsb, in_binsb);
dummy_pt_in_true =  new TH1D("dummy_pt_in_true","p_T^{in};p_T^{in} [GeV];Events", in_nbins_true, in_bins_true);
dummy_pt_in_meas =  new TH1D("dummy_pt_in_meas","p_T^{in};p_T^{in} [GeV];Events", in_nbins_meas, in_bins_meas);

TH1D *dummy_pt_out;
TH1D *dummy_pt_outb;
TH1D *dummy_pt_out_true;
TH1D *dummy_pt_out_meas;
dummy_pt_out =  new TH1D("dummy_pt_out","p_T^{out};p_T^{out} [GeV];Events", out_nbins, out_bins);
dummy_pt_outb =  new TH1D("dummy_pt_outb","p_T^{out};p_T^{out} [GeV];Events", out_nbinsb, out_binsb);
dummy_pt_out_true =  new TH1D("dummy_pt_out_true","p_T^{out};p_T^{out} [GeV];Events", out_nbins_true, out_bins_true);
dummy_pt_out_meas =  new TH1D("dummy_pt_out_meas","p_T^{out};p_T^{out} [GeV];Events", out_nbins_meas, out_bins_meas);

TH1D *dummy_etastar;
TH1D *dummy_etastarb;
TH1D *dummy_etastar_true;
TH1D *dummy_etastar_meas;
dummy_etastar =  new TH1D("dummy_etastar","#eta*;#eta*;Events", etastar_nbins, etastar_bins);
dummy_etastarb =  new TH1D("dummy_etastarb","#eta*;#eta*;Events", etastar_nbinsb, etastar_binsb);
dummy_etastar_true =  new TH1D("dummy_etastar_true","#eta*;#eta*;Events", etastar_nbins_true, etastar_bins_true);
dummy_etastar_meas =  new TH1D("dummy_etastar_meas","#eta*;#eta*;Events", etastar_nbins_meas, etastar_bins_meas);

TH1D *dummy_deta_out;
TH1D *dummy_deta_outb;
TH1D *dummy_deta_out_true;
TH1D *dummy_deta_out_meas;
dummy_deta_out =  new TH1D("dummy_deta_out","#Delta#eta^{out};#Delta#eta^{out};Events", deta_out_nbins, deta_out_bins);
dummy_deta_outb =  new TH1D("dummy_deta_outb","#Delta#eta^{out};#Delta#eta^{out};Events", deta_out_nbinsb, deta_out_binsb);
dummy_deta_out_true =  new TH1D("dummy_deta_out_true","#Delta#eta^{out};#Delta#eta^{out};Events", deta_out_nbins_true, deta_out_bins_true);
dummy_deta_out_meas =  new TH1D("dummy_deta_out_meas","#Delta#eta^{out};#Delta#eta^{out} [GeV];Events", deta_out_nbins_meas, deta_out_bins_meas);

  TH1D *ak5Gen_delta_phi_truth;
  TH1D *ak5PF_delta_phi_measured;
  TH1D *ak5Gen_delta_phi_truth_missed;
  TH1D *ak5Gen_delta_phi_truth_stayed;
  TH1D *ak5Gen_delta_phi_truth_fakes;
  TH1D *ak5PF_delta_phi_miss;
  TH1D *ak5PF_delta_phi_notmissed;
  TH1D *hist_events;

  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_25;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_26;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_27;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_28;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_29;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_30;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_31;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_32;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_33;
  TH1D *ak5Gen_delta_phi_truth_fakes_slightly_34;

  ak5Gen_delta_phi_truth =  new TH1D("ak5Gen_delta_phi_truth","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5PF_delta_phi_measured =  new TH1D("ak5PF_delta_phi_measured","#Delta#phi;#Delta#phi [rad];Events", 14, 0, 3.15);
  ak5Gen_delta_phi_truth_missed =  new TH1D("ak5Gen_delta_phi_truth_missed","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_stayed =  new TH1D("ak5Gen_delta_phi_truth_stayed","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes =  new TH1D("ak5Gen_delta_phi_truth_fakes","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5PF_delta_phi_miss =  new TH1D("ak5PF_delta_phi_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5PF_delta_phi_notmissed =  new TH1D("ak5PF_delta_phi_notmissed","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  hist_events =  new TH1D("Events","Selection Chain;Type;# Events", 13, 0, 13);

  ak5Gen_delta_phi_truth_fakes_slightly_25 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_25","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_26 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_26","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_27 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_27","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_28 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_28","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_29 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_29","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_30 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_30","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_31 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_31","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_32 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_32","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_33 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_33","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_truth_fakes_slightly_34 =  new TH1D("ak5Gen_delta_phi_truth_fakes_slightly_34","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);

  ak5Gen_delta_phi_truth->Sumw2();
  ak5PF_delta_phi_measured->Sumw2();
  ak5Gen_delta_phi_truth_missed->Sumw2();
  ak5Gen_delta_phi_truth_stayed->Sumw2();
  ak5Gen_delta_phi_truth_fakes->Sumw2();
  ak5PF_delta_phi_miss->Sumw2();
  ak5PF_delta_phi_notmissed->Sumw2();

  ak5Gen_delta_phi_truth_fakes_slightly_25->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_26->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_27->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_28->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_29->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_30->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_31->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_32->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_33->Sumw2();
  ak5Gen_delta_phi_truth_fakes_slightly_34->Sumw2();

  hist_events->Sumw2();
  hist_events->Fill("Total Events",0);
  hist_events->Fill("PV Selection",0);
  hist_events->Fill("Hadron Level Selected",0);
  hist_events->Fill("Hadron Level Cross-section",0);
  hist_events->Fill("Detector Level Selected",0);
  hist_events->Fill("Detector Level Cross-section",0);
  hist_events->Fill("Events in the Response Matrix",0);
  hist_events->Fill("Events in the Miss Histogram",0);
  hist_events->Fill("Events in the Truth Histogram",0);
  hist_events->Fill("Events in the Truth Missed Histogram",0);
  hist_events->Fill("Events in the Truth Stayed Histogram",0);
  hist_events->Fill("Events in the Truth Fakes Histogram",0);
  hist_events->Fill("Events in the Truth Fakes Slightly Histogram",0);


  TH1D *ak5Gen_delta_phi_miss;
  TH1D *ak5Gen_delta_phi_deta1_miss;
  TH1D *ak5Gen_delta_phi_deta2_miss;
  TH1D *ak5Gen_delta_phi_deta3_miss;
  TH1D *ak5Gen_delta_phi_deta4_miss;
  TH1D *ak5Gen_delta_phi_gap_miss;
  TH1D *ak5Gen_delta_phi_deta1_gap_miss;
  TH1D *ak5Gen_delta_phi_deta2_gap_miss;
  TH1D *ak5Gen_delta_phi_deta3_gap_miss;
  TH1D *ak5Gen_delta_phi_deta4_gap_miss;
  TH1D *ak5Gen_delta_phi_nogap_miss;
  TH1D *ak5Gen_delta_phi_deta1_nogap_miss;
  TH1D *ak5Gen_delta_phi_deta2_nogap_miss;
  TH1D *ak5Gen_delta_phi_deta3_nogap_miss;
  TH1D *ak5Gen_delta_phi_deta4_nogap_miss;
  TH1D *ak5Gen_leading_pt_inside_gap_miss;
  TH1D *ak5Gen_leading_eta_star_inside_gap_miss;
  TH1D *ak5Gen_leading_pt_outside_gap_miss;
  TH1D *ak5Gen_delta_eta_outside_gap_miss;

  ak5Gen_delta_phi_miss =  new TH1D("ak5Gen_delta_phi_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta1_miss =  new TH1D("ak5Gen_delta_phi_deta1_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta2_miss =  new TH1D("ak5Gen_delta_phi_deta2_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta3_miss =  new TH1D("ak5Gen_delta_phi_deta3_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta4_miss =  new TH1D("ak5Gen_delta_phi_deta4_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_gap_miss =  new TH1D("ak5Gen_delta_phi_gap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta1_gap_miss =  new TH1D("ak5Gen_delta_phi_deta1_gap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta2_gap_miss =  new TH1D("ak5Gen_delta_phi_deta2_gap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta3_gap_miss =  new TH1D("ak5Gen_delta_phi_deta3_gap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta4_gap_miss =  new TH1D("ak5Gen_delta_phi_deta4_gap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_nogap_miss =  new TH1D("ak5Gen_delta_phi_nogap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta1_nogap_miss =  new TH1D("ak5Gen_delta_phi_deta1_nogap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta2_nogap_miss =  new TH1D("ak5Gen_delta_phi_deta2_nogap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta3_nogap_miss =  new TH1D("ak5Gen_delta_phi_deta3_nogap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_deta4_nogap_miss =  new TH1D("ak5Gen_delta_phi_deta4_nogap_miss","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_leading_pt_inside_gap_miss =  new TH1D("ak5Gen_leading_pt_inside_gap_miss","p_{T};p_{T} [GeV];Events", in_nbins, in_bins);
  ak5Gen_leading_eta_star_inside_gap_miss =  new TH1D("ak5Gen_leading_eta_star_inside_gap_miss","#eta*;#eta*;Events", etastar_nbins, etastar_bins);
  ak5Gen_leading_pt_outside_gap_miss =  new TH1D("ak5Gen_leading_pt_outside_gap_miss","p_{T};p_{T} [GeV];Events", out_nbins, out_bins);
  ak5Gen_delta_eta_outside_gap_miss =  new TH1D("ak5Gen_delta_eta_outside_gap_miss","#Delta#eta^{out};#Delta#eta^{out};Events", deta_out_nbins, deta_out_bins);


  ak5Gen_delta_phi_miss->Sumw2();
  ak5Gen_delta_phi_deta1_miss->Sumw2();
  ak5Gen_delta_phi_deta2_miss->Sumw2();
  ak5Gen_delta_phi_deta3_miss->Sumw2();
  ak5Gen_delta_phi_deta4_miss->Sumw2();
  ak5Gen_delta_phi_gap_miss->Sumw2();
  ak5Gen_delta_phi_deta1_gap_miss->Sumw2();
  ak5Gen_delta_phi_deta2_gap_miss->Sumw2();
  ak5Gen_delta_phi_deta3_gap_miss->Sumw2();
  ak5Gen_delta_phi_deta4_gap_miss->Sumw2();
  ak5Gen_delta_phi_nogap_miss->Sumw2();
  ak5Gen_delta_phi_deta1_nogap_miss->Sumw2();
  ak5Gen_delta_phi_deta2_nogap_miss->Sumw2();
  ak5Gen_delta_phi_deta3_nogap_miss->Sumw2();
  ak5Gen_delta_phi_deta4_nogap_miss->Sumw2();
  ak5Gen_leading_pt_inside_gap_miss->Sumw2();
  ak5Gen_leading_eta_star_inside_gap_miss->Sumw2();
  ak5Gen_leading_pt_outside_gap_miss->Sumw2();
  ak5Gen_delta_eta_outside_gap_miss->Sumw2();

  RooUnfoldResponse resp_delta_phi_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi", "delta_phi");

  RooUnfoldResponse resp_delta_phi_deta1_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta1_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta1(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta1", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta2_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta2_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta2(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta2", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta3_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta3_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta3(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta3", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta4_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta4_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta4(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta4", "delta_phi");

  RooUnfoldResponse resp_delta_phi_gap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_gap_all", "delta_phi_gap");
  RooUnfoldResponse resp_delta_phi_gap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_gap", "delta_phi");

  RooUnfoldResponse resp_delta_phi_deta1_gap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta1_gap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta1_gap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta1_gap", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta2_gap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta2_gap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta2_gap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta2_gap", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta3_gap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta3_gap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta3_gap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta3_gap", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta4_gap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta4_gap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta4_gap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta4_gap", "delta_phi");

  RooUnfoldResponse resp_delta_phi_nogap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_nogap_all", "delta_phi_nogap");
  RooUnfoldResponse resp_delta_phi_nogap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_nogap", "delta_phi");

  RooUnfoldResponse resp_delta_phi_deta1_nogap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta1_nogap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta1_nogap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta1_nogap", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta2_nogap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta2_nogap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta2_nogap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta2_nogap", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta3_nogap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta3_nogap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta3_nogap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta3_nogap", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta4_nogap_all(dummy_delta_phi_meas, dummy_delta_phi_true, "resp_delta_phi_deta4_nogap_all", "delta_phi");
  RooUnfoldResponse resp_delta_phi_deta4_nogap(dummy_delta_phib, dummy_delta_phi, "resp_delta_phi_deta4_nogap", "delta_phi");

  RooUnfoldResponse resp_leading_pt_inside_gap_all(dummy_pt_in_meas, dummy_pt_in_true, "resp_leading_pt_inside_gap_all", "leading_pt_inside");
  RooUnfoldResponse resp_leading_pt_inside_gap(dummy_pt_inb, dummy_pt_in, "resp_leading_pt_inside_gap", "leading_pt_inside");

  RooUnfoldResponse resp_leading_eta_star_inside_gap_all(dummy_etastar_meas, dummy_etastar_true, "resp_leading_eta_star_inside_gap_all", "leading_pt_inside");
  RooUnfoldResponse resp_leading_eta_star_inside_gap(dummy_etastarb, dummy_etastar, "resp_leading_eta_star_inside_gap", "leading_pt_inside");

  RooUnfoldResponse resp_leading_pt_outside_gap_all(dummy_pt_out_meas, dummy_pt_out_true, "resp_leading_pt_outside_gap_all", "leading_pt_outside");
  RooUnfoldResponse resp_leading_pt_outside_gap(dummy_pt_outb, dummy_pt_out, "resp_leading_pt_outside_gap", "leading_pt_outside");

  RooUnfoldResponse resp_delta_eta_outside_gap_all(dummy_deta_out_meas, dummy_deta_out_true, "resp_delta_eta_outside_gap_all", "delta_eta_outside");
  RooUnfoldResponse resp_delta_eta_outside_gap(dummy_deta_outb, dummy_deta_out, "resp_delta_eta_outside_gap", "delta_eta_outside");

for (int z = 0; z < n_files; z++)
{
  string file = data_in[z];
  cout<<z+1<<"/"<<n_files<<" -> "<<file<<endl;
  TFile *inf = TFile::Open( file.c_str() );
  //cout<<"opened"<<endl;
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  //----------- define the tree branch -
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);
  //----------- settings ---------------

  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  prescale = data_lumi[z];
  cout<<"Cross-section per event : "<<prescale<<endl;
  int decade = 0;

  if (test) { nentries = 1000; } //reduced number of read entries, usefull for testing

   for (Int_t i=0;i<nentries;i++) {
     counter_entries++; 
     double progress = 10.0*i/(1.0*nentries);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
      cout<<10*k<<" % "<<endl;
     decade = k; 
   
     tr->GetEntry(i);

      //-------- check if the primary vertex is good ----
        pass_pv_det = false;
        pass_pv_gen = false;
	if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)	{ pass_pv_gen = true; }
	if (sel_mode == "allvertex") 					{ pass_pv_gen = true; }
	if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)				{ pass_pv_det = true; }
	if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) 			{ pass_pv_det = true; }
	if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1) { pass_pv_det = true; }
      
	vertex_factor = 1.0;
	if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1); }
     
     delta_phi_gen = -0.9;
     pt_forward_gen = 0.0;
     eta_forward_gen = 0.0;
     eta_forward2_gen = 0.0;
     phi_forward_gen = 0.0;
     pt_central_gen = 0.0;
     eta_central_gen = 0.0;
     eta_central2_gen = 0.0;
     phi_central_gen = 0.0;
     pt2_forward_gen = 0.0;
     eta2_forward_gen = 0.0;
     eta2_forward2_gen = 0.0;
     phi2_forward_gen = 0.0;
     pt2_central_gen = 0.0;
     eta2_central_gen = 0.0;
     eta2_central2_gen = 0.0;
     phi2_central_gen = 0.0;
     pass_gen = false;
     pass_gap_gen = false;
     pass_nogap_gen = false;
     pass_deta1_gen = false;
     pass_deta2_gen = false;
     pass_deta3_gen = false;
     pass_deta4_gen = false;
     pass_out_gen = false;
     pt_total_gap_gen = 0.0;
     eta_gap_gen = -10.0;
     phi_gap_gen = -10.0;
     pt_total_outside_gen = 0.0;
     deta_out1_gen = -1.0;
     deta_out2_gen = -1.0;
     deta_out_gen = -1.0;
     eta_outside_gen = -10.0;
     phi_outside_gen = -10.0;
     pt_leading_gap_gen = 0.0;
     pt_leading_outside_gen = 0.0;
     eta_star_inside_gen = 0;


     if (pass_pv_gen) {
     
    for(unsigned int j=0;j<Event->nGenJets();j++) {
     pt = Event->genjet(j).pt();
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
     if (pt >= 25.0) {
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central_gen)
     { pt_central_gen = pt; eta_central_gen = eta; phi_central_gen = phi; }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward_gen)
     { pt_forward_gen = pt; eta_forward_gen = eta; phi_forward_gen = phi; }
     }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt < pt_forward_gen && pt > pt2_forward_gen)
     { pt2_forward_gen = pt; eta2_forward_gen = eta; phi2_forward_gen = phi;}
     if (eta <= 2.8 && eta >= -2.8 && pt > pt2_central_gen && pt < pt_central_gen)
     { pt2_central_gen = pt; eta2_central_gen = eta; phi2_central_gen = phi;}
     }

     if (pt_forward_gen > 25.0 && pt_central_gen > 25.0)
     {
     if (pt_central_gen > 20.0 and pt2_central_gen > 0.0)
     {
     delta_phi2_gen = calc_delta_phi(phi_central_gen, phi2_central_gen);
     delta_pt_gen = pt_central_gen - pt2_central_gen;
     }
     if (pt_forward_gen > 20.0 and pt2_forward_gen > 0.0)
     {
     delta_phi2_gen = calc_delta_phi(phi_forward_gen, phi2_forward_gen);
     delta_pt_gen = pt_forward_gen - pt2_forward_gen;
     }
     
     eta_central2_gen = eta_central_gen;     
     if (eta_central2_gen < 0) { eta_central2_gen = -eta_central2_gen; }
     eta_forward2_gen = eta_forward_gen;     
     if (eta_forward2_gen < 0) { eta_forward2_gen = -eta_forward2_gen; }
     delta_eta_gen = eta_forward_gen - eta_central_gen;
     if (delta_eta_gen < 0) { delta_eta_gen = -delta_eta_gen; }
     delta_phi_gen = calc_delta_phi(phi_forward_gen, phi_central_gen);

	if (pt_forward_gen > pt_min && pt_central_gen > pt_min)
     	{
     	counter_selected_gen++;
	pass_gen = true;
     	counter_weigth_gen = counter_weigth_gen + prescale * vertex_factor;
	if (delta_eta_gen > 0.4 and delta_eta_gen <= 2.5) {pass_deta1_gen = true; }
	if (delta_eta_gen > 2.5 and delta_eta_gen <= 3.5) {pass_deta2_gen = true; }
	if (delta_eta_gen > 3.5 and delta_eta_gen <= 4.5) {pass_deta3_gen = true; }
	if (delta_eta_gen > 4.5 and delta_eta_gen <= 7.5) {pass_deta4_gen = true; }
	}

     for(unsigned int j=0;j<Event->nGenJets();j++) {
     pt = Event->genjet(j).pt();
     if (pt >= 5) {
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
     if (eta_central_gen > eta_forward_gen && eta > eta_forward_gen && eta < eta_central_gen ) {
        if (pt > pt_leading_gap_gen) {pt_leading_gap_gen = pt; eta_gap_gen = eta; phi_gap_gen = phi; }
        pt_total_gap_gen = pt_total_gap_gen + pt;
        }
     if (eta_central_gen < eta_forward_gen && eta < eta_forward_gen && eta > eta_central_gen ) {
        if (pt > pt_leading_gap_gen) {pt_leading_gap_gen = pt; eta_gap_gen = eta; phi_gap_gen = phi; }
        pt_total_gap_gen = pt_total_gap_gen + pt;
        }
     if (eta_central_gen > eta_forward_gen && (eta < eta_forward_gen || eta > eta_central_gen) ) {
        if (pt > pt_leading_outside_gen) {pt_leading_outside_gen = pt; eta_outside_gen = eta; phi_outside_gen = phi; }
        pt_total_outside_gen = pt_total_outside_gen + pt;
        }
     if (eta_central_gen < eta_forward_gen && (eta > eta_forward_gen || eta < eta_central_gen) ) {
        if (pt > pt_leading_outside_gen) {pt_leading_outside_gen = pt; eta_outside_gen = eta; phi_outside_gen = phi; }
        pt_total_outside_gen = pt_total_outside_gen + pt;
        }
  //   cout<<"here6 "<<gap_req<<endl;
  //   if (pt > gap_req && eta_central > eta_forward && eta > eta_forward && eta > eta_central )
  //   {      cout<<i<<" [0]pt_gap = "<<pt<<endl; pass_gap = 1; }
  //   cout<<"here1"<<endl;
  //   if (pt > gap_req && eta_central < eta_forward && eta < eta_forward && eta < eta_central )
  //   {      cout<<i<<" [0]pt_gap = "<<pt<<endl; pass_gap = 1; }
  //   cout<<"here2"<<endl;
     }
     }
     if (pt_leading_gap_gen > gap_req )
     	{
    	pass_nogap_gen = true;
	eta_star_inside_gen = eta_gap_gen - (eta_forward_gen + eta_central_gen)/2;
     	}
     else
     	{
     	pass_gap_gen = true;
     	}

     if (pt_leading_outside_gen > gap_req )
     	{
     	pass_out_gen = true;
     	deta_out1_gen = eta_central_gen - eta_outside_gen;
     	if (deta_out1_gen < 0) { deta_out1_gen = -deta_out1_gen; }
     	deta_out2_gen = eta_forward_gen - eta_outside_gen;
     	if (deta_out2_gen < 0) { deta_out2_gen = -deta_out2_gen; }
	if (deta_out1_gen < deta_out2_gen) { deta_out_gen = deta_out1_gen; }
	if (deta_out2_gen < deta_out1_gen) { deta_out_gen = deta_out2_gen; }
     	}

     if (pass_gap_gen == true)
     	{
     	counter_gap_gen = counter_gap_gen + 1;
     	}
     if (pass_nogap_gen == true)
     	{
     	counter_nogap_gen = counter_nogap_gen + 1;
     	}
     }
     
     delta_phi_det = -.9;
     pt_forward_det = 0.0;
     eta_forward_det = 0.0;
     eta_forward2_det = 0.0;
     phi_forward_det = 0.0;
     pt_central_det = 0.0;
     eta_central_det = 0.0;
     eta_central2_det = 0.0;
     phi_central_det = 0.0;
     pt2_forward_det = 0.0;
     eta2_forward_det = 0.0;
     eta2_forward2_det = 0.0;
     phi2_forward_det = 0.0;
     pt2_central_det = 0.0;
     eta2_central_det = 0.0;
     eta2_central2_det = 0.0;
     phi2_central_det = 0.0;
     pass_det = false;
     pass_gap_det = false;
     pass_nogap_det = false;
     pass_deta1_det = false;
     pass_deta2_det = false;
     pass_deta3_det = false;
     pass_deta4_det = false;
     pass_out_det = false;
     pt_total_gap_det = 0.0;
     eta_gap_det = -10.0;
     phi_gap_det = -10.0;
     pt_total_outside_det = 0.0;
     deta_out1_det = -1.0;
     deta_out2_det = -1.0;
     deta_out_det = -1.0;
     eta_outside_det = -10.0;
     phi_outside_det = -10.0;
     pt_leading_gap_det = 0.0;
     pt_leading_outside_det = 0.0;
     eta_star_inside_det = 0;
     
     if (pass_pv_det) {

     for(unsigned int j=0;j<Event->nPFJets();j++) {
     pt = Event->pfjet(j).ptCor();
     eta = Event->pfjet(j).eta();
     phi = Event->pfjet(j).phi();
     gen_pt = Event->genjet(j).pt();
     gen_eta = Event->genjet(j).eta();
     gen_phi = Event->genjet(j).phi();
     pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,false);
     if (pt >= 25.0 && Event->pfjet(j).tightID()) {
     //cout<<"id = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" correction = "<<cor<<endl;
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central_det)
     { pt_central_det = pt; eta_central_det = eta; phi_central_det = phi; } // chm_central = chm; elm_central = elm; }
     // if (eta <= 2.8 && eta >= -2.8 && pt > pt_min)
     // { ak5PF_inclusive_central_pt->Fill(pt,prescale); }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward_det)
     { pt_forward_det = pt; eta_forward_det = eta; phi_forward_det = phi; } // chm_forward = chm; elm_forward = elm; }
     // if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_min)
     // { ak5PF_inclusive_forward_pt->Fill(pt,prescale); }
     }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt < pt_forward_det && pt > pt2_forward_det)
     { pt2_forward_det = pt; eta2_forward_det = eta; phi2_forward_det = phi;}
     if (eta <= 2.8 && eta >= -2.8 && pt > pt2_central_det && pt < pt_central_det)
     { pt2_central_det = pt; eta2_central_det = eta; phi2_central_det = phi;}
     }
     
     // if (pt_central > pt_min) { ak5PF_inclusive_leading_central_pt->Fill(pt_central,prescale); }
     // if (pt_forward > pt_min) { ak5PF_inclusive_leading_forward_pt->Fill(pt_forward,prescale); }

     if (pt_forward_det > 25.0 && pt_central_det > 25.0)
     {
     
     // ak5PF_central_pt_rel->Fill(pt_central,pt2_central,prescale);
     // ak5PF_forward_pt_rel->Fill(pt_forward,pt2_forward,prescale);
     if (pt_central_det > 20.0 and pt2_central_det > 0.0)
     {
     delta_phi2_det = calc_delta_phi(phi_central_gen, phi2_central_gen);
     // ak5PF_delta_phi_central_rel->Fill(delta_phi2,prescale);
     delta_pt_det = pt_central_det - pt2_central_det;
     // ak5PF_delta_pt_central_rel->Fill(delta_pt,prescale);
     // if (delta_pt < 5.0) { ak5PF_delta_phi_central_rel_small->Fill(delta_phi2,prescale); }
     // if (delta_pt > 5.0 && delta_pt < 15.0) { ak5Gen_delta_phi_central_rel_medium->Fill(delta_phi2,prescale); }
     // if (delta_pt > 15.0) { ak5Gen_delta_phi_central_rel_large->Fill(delta_phi2,prescale); }
     }
     if (pt_forward_det > 20.0 and pt2_forward_det > 0.0)
     {
     delta_phi2_det = calc_delta_phi(phi_forward_gen, phi2_forward_gen);
     delta_pt_det = pt_forward_det - pt2_forward_det;
     }
     
     //cout<<"id = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" correction = "<<cor<<endl;
     eta_central2_det = eta_central_det;     
     if (eta_central2_det < 0) { eta_central2_det = -eta_central2_det; }
     eta_forward2_det = eta_forward_det;     
     if (eta_forward2_det < 0) { eta_forward2_det = -eta_forward2_det; }
     delta_eta_det = eta_forward_det - eta_central_det;
     if (delta_eta_det < 0) { delta_eta_det = -delta_eta_det; }
     delta_phi_det = calc_delta_phi(phi_forward_det, phi_central_det);

     if (pt_forward_det > pt_min and pt_central_det > pt_min)
	{	
	counter_selected_det++;
	pass_det = true;
	counter_weigth_det = counter_weigth_det + prescale * vertex_factor;
	if (delta_eta_det > 0.4 and delta_eta_det <= 2.5) {pass_deta1_det = true; }
	if (delta_eta_det > 2.5 and delta_eta_det <= 3.5) {pass_deta2_det = true; }
	if (delta_eta_det > 3.5 and delta_eta_det <= 4.5) {pass_deta3_det = true; }
	if (delta_eta_det > 4.5 and delta_eta_det <= 7.5) {pass_deta4_det = true; }
	}

     
  //   cout<<"here4 "<<gap_req<<endl;     

     for(unsigned int j=0;j<Event->nGenJets();j++) {
     pt = Event->genjet(j).pt();
     if (pt >= 5) {
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
     gen_pt = Event->genjet(j).pt();
     gen_eta = Event->genjet(j).eta();
     gen_phi = Event->genjet(j).phi();
     pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,false);
     if (eta_central_det > eta_forward_det && eta > eta_forward_det && eta < eta_central_det ) {
        if (pt > pt_leading_gap_det) {pt_leading_gap_det = pt; eta_gap_det = eta; phi_gap_det = phi; }
        pt_total_gap_det = pt_total_gap_det + pt;
        }
     if (eta_central_det < eta_forward_det && eta < eta_forward_det && eta > eta_central_det ) {
        if (pt > pt_leading_gap_det) {pt_leading_gap_det = pt; eta_gap_det = eta; phi_gap_det = phi; }
        pt_total_gap_det = pt_total_gap_det + pt;
        }
     if (eta_central_det > eta_forward_det && (eta < eta_forward_det || eta > eta_central_det) ) {
        if (pt > pt_leading_outside_det) {pt_leading_outside_det = pt; eta_outside_det = eta; phi_outside_det = phi; }
        pt_total_outside_det = pt_total_outside_det + pt;
        }
     if (eta_central_det < eta_forward_det && (eta > eta_forward_det || eta < eta_central_det) ) {
        if (pt > pt_leading_outside_det) {pt_leading_outside_det = pt; eta_outside_det = eta; phi_outside_det = phi; }
        pt_total_outside_det = pt_total_outside_det + pt;
        }
  //   cout<<"here6 "<<gap_req<<endl;
  //   if (pt > gap_req && eta_central > eta_forward && eta > eta_forward && eta > eta_central )
  //   {      cout<<i<<" [0]pt_gap = "<<pt<<endl; pass_gap = 1; }
  //   cout<<"here1"<<endl;
  //   if (pt > gap_req && eta_central < eta_forward && eta < eta_forward && eta < eta_central )
  //   {      cout<<i<<" [0]pt_gap = "<<pt<<endl; pass_gap = 1; }
  //   cout<<"here2"<<endl;
     }
     }
     if (pt_leading_gap_det > gap_req )
     	{
	pass_nogap_det = true;
	//     cout<<i<<" [1]pt_gap = "<<pt_leading_gap<<endl;
     	eta_star_inside_det = eta_gap_det - (eta_forward_det + eta_central_det)/2;
     	}
     else
     	{
     	pass_gap_det = true;
     	}

     if (pt_leading_outside_det > gap_req )
     	{
     	pass_out_det = true;
     	deta_out1_det = eta_central_det - eta_outside_det;
     	if (deta_out1_det < 0) { deta_out1_det = -deta_out1_det; }
     	deta_out2_det = eta_forward_det - eta_outside_det;
     	if (deta_out2_det < 0) { deta_out2_det = -deta_out2_det; }
	if (deta_out1_det < deta_out2_det) { deta_out_det = deta_out1_det; }
	if (deta_out2_det < deta_out1_det) { deta_out_det = deta_out2_det; }
     	}

     if (pass_gap_det == true)
     	{
     	counter_gap_det = counter_gap_det + 1;
     	}
     if (pass_nogap_det == true)
     	{
     	counter_nogap_det = counter_nogap_det + 1;
     	}
     }
     
     if ((pt_central_det > pt_min and pt_forward_det > pt_min) or (pt_central_gen > pt_min and pt_forward_gen > pt_min))
     {
     if (test)
	{
	cout<< "pT central det = " << pt_central_det << endl;
	cout<< "pT forward det = " << pt_forward_det << endl;
	cout<< "delta phi det =  " << delta_phi_det << endl;
	cout<< "pT central gen = " << pt_central_gen << endl;
	cout<< "pT forward gen = " << pt_forward_gen << endl;
	cout<< "delta phi gen =  " << delta_phi_gen << endl;
	}
     //filling response matrix
     if (pass_det and pass_gen)
	{
	if (test) { cout<< "Response filled" << endl; }
	counter_response++;
	resp_delta_phi.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
	resp_delta_phi_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
	ak5PF_delta_phi_notmissed->Fill(delta_phi_det,prescale * vertex_factor);

	if (pass_deta1_det and pass_deta1_gen)
		{
		resp_delta_phi_deta1.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_deta2_det and pass_deta2_gen)
		{
		resp_delta_phi_deta2.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_deta3_det and pass_deta3_gen)
		{
		resp_delta_phi_deta3.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_deta4_det and pass_deta4_gen)
		{
		resp_delta_phi_deta4.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		}

	if (pass_gap_det and pass_gap_gen)
		{
		resp_delta_phi_gap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_gap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		if (pass_deta1_det and pass_deta1_gen)
			{
			resp_delta_phi_deta1_gap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta1_gap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		if (pass_deta2_det and pass_deta2_gen)
			{
			resp_delta_phi_deta2_gap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta2_gap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		if (pass_deta3_det and pass_deta3_gen)
			{
			resp_delta_phi_deta3_gap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta3_gap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		if (pass_deta4_det and pass_deta4_gen)
			{
			resp_delta_phi_deta4_gap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta4_gap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		}
	if (pass_nogap_det and pass_nogap_det)
		{
		resp_delta_phi_nogap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_nogap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
		resp_leading_pt_inside_gap.Fill(pt_leading_gap_det,pt_leading_gap_gen,prescale * vertex_factor);
		resp_leading_pt_inside_gap_all.Fill(pt_leading_gap_det,pt_leading_gap_gen,prescale * vertex_factor);
		resp_leading_eta_star_inside_gap.Fill(eta_star_inside_det,eta_star_inside_gen,prescale * vertex_factor);
		resp_leading_eta_star_inside_gap_all.Fill(eta_star_inside_det,eta_star_inside_gen,prescale * vertex_factor);
		if (pass_deta1_det and pass_deta1_gen)
			{
			resp_delta_phi_deta1_nogap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta1_nogap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		if (pass_deta2_det and pass_deta2_gen)
			{
			resp_delta_phi_deta2_nogap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta2_nogap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		if (pass_deta3_det and pass_deta3_gen)
			{
			resp_delta_phi_deta3_nogap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta3_nogap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		if (pass_deta4_det and pass_deta4_gen)
			{
			resp_delta_phi_deta4_nogap.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			resp_delta_phi_deta4_nogap_all.Fill(delta_phi_det,delta_phi_gen,prescale * vertex_factor);
			}
		}
	if (pass_out_det and pass_out_gen)
		{
		resp_leading_pt_outside_gap.Fill(pt_leading_outside_det,pt_leading_outside_gen,prescale * vertex_factor);
		resp_leading_pt_outside_gap_all.Fill(pt_leading_outside_det,pt_leading_outside_gen,prescale * vertex_factor);
		resp_delta_eta_outside_gap.Fill(deta_out_det,deta_out_gen,prescale * vertex_factor);
		resp_delta_eta_outside_gap_all.Fill(deta_out_det,deta_out_gen,prescale * vertex_factor);
		}

	}


        if (!pass_gen and pass_det)
		{
		if (test) { cout<< "Fake filled" << endl; }
		//resp_delta_phi.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		resp_delta_phi.Fake(delta_phi_det,prescale * vertex_factor);
		}

	if (pass_deta1_det and !pass_deta1_gen)
		{
		resp_delta_phi_deta1.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta1_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta2_det and !pass_deta2_gen)
		{
		resp_delta_phi_deta2.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta2_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta3_det and !pass_deta3_gen)
		{
		resp_delta_phi_deta3.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta3_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta4_det and !pass_deta4_gen)
		{
		resp_delta_phi_deta4.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta4_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}

	if (!pass_gap_gen and pass_gap_det)
		{
		resp_delta_phi_gap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		resp_delta_phi_gap.Fake(delta_phi_det,prescale * vertex_factor);
		}

	if (pass_deta1_det and pass_gap_det and (!pass_deta1_gen or !pass_gap_gen))
		{
		resp_delta_phi_deta1_gap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta1_gap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta2_det and pass_gap_det and (!pass_deta2_gen or !pass_gap_gen))
		{
		resp_delta_phi_deta2_gap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta2_gap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta3_det and pass_gap_det and (!pass_deta3_gen or !pass_gap_gen))
		{
		resp_delta_phi_deta3_gap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta3_gap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta4_det and pass_gap_det and (!pass_deta4_gen or !pass_gap_gen))
		{
		resp_delta_phi_deta4_gap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta4_gap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}

	if (!pass_nogap_gen and pass_nogap_det)
		{
		resp_delta_phi_nogap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		resp_delta_phi_nogap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_leading_pt_inside_gap.Fake(pt_leading_gap_det,prescale * vertex_factor);
		resp_leading_pt_inside_gap_all.Fill(pt_leading_gap_det,10.0,prescale * vertex_factor);
		resp_leading_eta_star_inside_gap.Fake(eta_star_inside_det,prescale * vertex_factor);
		resp_leading_eta_star_inside_gap_all.Fill(eta_star_inside_det,-9.0,prescale * vertex_factor);
		}

	if (pass_deta1_det and pass_nogap_det and (!pass_deta1_gen or !pass_nogap_gen))
		{
		resp_delta_phi_deta1_nogap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta1_nogap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta2_det and pass_nogap_det and (!pass_deta2_gen or !pass_nogap_gen))
		{
		resp_delta_phi_deta2_nogap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta2_nogap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta3_det and pass_nogap_det and (!pass_deta3_gen or !pass_nogap_gen))
		{
		resp_delta_phi_deta3_nogap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta3_nogap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}
	if (pass_deta4_det and pass_nogap_det and (!pass_deta4_gen or !pass_nogap_gen))
		{
		resp_delta_phi_deta4_nogap.Fake(delta_phi_det,prescale * vertex_factor);
		resp_delta_phi_deta4_nogap_all.Fill(delta_phi_det,-0.5,prescale * vertex_factor);
		}

        if (!pass_det and pass_gen)
		{
		if (test) { cout<< "Miss filled" << endl; }
		counter_miss++;
		resp_delta_phi.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		ak5Gen_delta_phi_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		ak5PF_delta_phi_miss->Fill(delta_phi_det,prescale * vertex_factor);
		}

	if (!pass_deta1_det and pass_deta1_gen)
		{
		ak5Gen_delta_phi_deta1_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (!pass_deta2_det and pass_deta2_gen)
		{
		ak5Gen_delta_phi_deta2_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (!pass_deta3_det and pass_deta3_gen)
		{
		ak5Gen_delta_phi_deta3_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (!pass_deta4_det and pass_deta4_gen)
		{
		ak5Gen_delta_phi_deta4_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}

	if (!pass_gap_det and pass_gap_gen)
		{
		ak5Gen_delta_phi_gap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_gap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_gap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}

	if (pass_gap_gen and pass_deta1_gen and (!pass_deta1_det or !pass_gap_det))
		{
		ak5Gen_delta_phi_deta1_gap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1_gap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1_gap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_gap_gen and pass_deta2_gen and (!pass_deta2_det or !pass_gap_det))
		{
		ak5Gen_delta_phi_deta2_gap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2_gap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2_gap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_gap_gen and pass_deta3_gen and (!pass_deta3_det or !pass_gap_det))
		{
		ak5Gen_delta_phi_deta3_gap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3_gap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3_gap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_gap_gen and pass_deta4_gen and (!pass_deta4_det or !pass_gap_det))
		{
		ak5Gen_delta_phi_deta4_gap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4_gap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4_gap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}

	if (!pass_nogap_det and pass_nogap_gen)
		{
		ak5Gen_delta_phi_nogap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_nogap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_nogap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		ak5Gen_leading_pt_inside_gap_miss->Fill(pt_leading_gap_gen,prescale * vertex_factor);
		resp_leading_pt_inside_gap.Miss(pt_leading_gap_gen,prescale * vertex_factor);
		resp_leading_pt_inside_gap_all.Fill(10.0,pt_leading_gap_gen,prescale * vertex_factor);
		ak5Gen_leading_eta_star_inside_gap_miss->Fill(eta_star_inside_gen,prescale * vertex_factor);
		resp_leading_eta_star_inside_gap.Miss(eta_star_inside_gen,prescale * vertex_factor);
		resp_leading_eta_star_inside_gap_all.Fill(-9.0,eta_star_inside_gen,prescale * vertex_factor);
		}

	if (pass_nogap_gen and pass_deta1_gen and (!pass_deta1_det or !pass_nogap_det))
		{
		ak5Gen_delta_phi_deta1_nogap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1_nogap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta1_nogap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_nogap_gen and pass_deta2_gen and (!pass_deta2_det or !pass_nogap_det))
		{
		ak5Gen_delta_phi_deta2_nogap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2_nogap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta2_nogap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_nogap_gen and pass_deta3_gen and (!pass_deta3_det or !pass_nogap_det))
		{
		ak5Gen_delta_phi_deta3_nogap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3_nogap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta3_nogap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}
	if (pass_nogap_gen and pass_deta4_gen and (!pass_deta4_det or !pass_nogap_det))
		{
		ak5Gen_delta_phi_deta4_nogap_miss->Fill(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4_nogap.Miss(delta_phi_gen,prescale * vertex_factor);
		resp_delta_phi_deta4_nogap_all.Fill(-0.5,delta_phi_gen,prescale * vertex_factor);
		}

	if (!pass_out_det and pass_out_gen)
		{
		ak5Gen_leading_pt_outside_gap_miss->Fill(pt_leading_outside_gen,prescale * vertex_factor);
		resp_leading_pt_outside_gap.Miss(pt_leading_outside_gen,prescale * vertex_factor);
		resp_leading_pt_outside_gap_all.Fill(10.0,pt_leading_outside_gen,prescale * vertex_factor);
		ak5Gen_delta_eta_outside_gap_miss->Fill(deta_out_gen,prescale * vertex_factor);
		resp_delta_eta_outside_gap.Miss(deta_out_gen,prescale * vertex_factor);
		resp_delta_eta_outside_gap_all.Fill(-0.5,deta_out_gen,prescale * vertex_factor);
		}


        if (pt_central_det > pt_min and pt_forward_det > pt_min)
     		{
     		ak5PF_delta_phi_measured->Fill(delta_phi_det,prescale * vertex_factor);
     		}
     //fill additional histograms
     if (pt_central_gen > pt_min and pt_forward_gen > pt_min)
     {
     if (test) { cout<< "Truth filled" << endl; }
     counter_truth++;
     ak5Gen_delta_phi_truth->Fill(delta_phi_gen,prescale * vertex_factor);
     if (pt_central_det < pt_min or pt_forward_det < pt_min)
	{
	if (test) { cout<< "Truth Miss filled" << endl; }
	counter_truth_missed++;
	ak5Gen_delta_phi_truth_missed->Fill(delta_phi_gen,prescale * vertex_factor);
	}
     if (pt_central_det > pt_min and pt_forward_det > pt_min)
	{
	if (test) { cout<< "Truth Stayed filled" << endl; }
	counter_truth_stayed++;
	ak5Gen_delta_phi_truth_stayed->Fill(delta_phi_gen,prescale * vertex_factor);
	}
     }
     else
     {
     if (pt_central_det > pt_min and pt_forward_det > pt_min)
	{
	if (test) { cout<< "Truth Fakes filled" << endl; }
	counter_truth_fakes++;
	ak5Gen_delta_phi_truth_fakes->Fill(delta_phi_gen,prescale * vertex_factor);
	if (pt_central_gen > 25.0 and pt_forward_gen > 25.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 25 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_25->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 26.0 and pt_forward_gen > 26.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 26 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_26->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 27.0 and pt_forward_gen > 27.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 27 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_27->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 28.0 and pt_forward_gen > 28.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 28 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_28->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 29.0 and pt_forward_gen > 29.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 29 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_29->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 30.0 and pt_forward_gen > 30.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled" << endl; }
		counter_truth_fakes_slightly++;
		ak5Gen_delta_phi_truth_fakes_slightly_30->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 31.0 and pt_forward_gen > 31.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 31 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_31->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 32.0 and pt_forward_gen > 32.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 32 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_32->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 33.0 and pt_forward_gen > 33.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 33 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_33->Fill(delta_phi_gen,prescale * vertex_factor);
		}
	if (pt_central_gen > 34.0 and pt_forward_gen > 34.0)
		{
		if (test) { cout<< "Truth Fakes Slighty filled over 34 GeV" << endl; }
		ak5Gen_delta_phi_truth_fakes_slightly_34->Fill(delta_phi_gen,prescale * vertex_factor);
		}
     	}
	}
     }
     
     
  } 


  }  }
  
}

//output a summary
  cout<<endl<<endl;
  cout<<"Events read:                                  "<<counter_entries<<endl;
  cout<<"Events after the PV cut:                      "<<counter_pv<<endl;
  cout<<"Events Main Selection Hadron:                 "<<counter_selected_gen<<endl;
  cout<<"Events Weight Hadron:                         "<<counter_weigth_gen<<endl;
  cout<<"Events Main Selection Detector:               "<<counter_selected_det<<endl;
  cout<<"Events Weight Detector:                       "<<counter_weigth_det<<endl;
  cout<<"Events in the Response Matrix:                "<<counter_response<<endl;
  cout<<"Events in the Miss Histogram:                 "<<counter_miss<<endl;
  cout<<"Events in the Truth Histogram:                "<<counter_truth<<endl;
  cout<<"Events in the Truth Missed Histogram:         "<<counter_truth_missed<<endl;
  cout<<"Events in the Truth Stayed Histogram:         "<<counter_truth_stayed<<endl;
  cout<<"Events in the Truth Fakes Histogram:          "<<counter_truth_fakes<<endl;
  cout<<"Events in the Truth Fakes Slightly Histogram: "<<counter_truth_fakes_slightly<<endl;

  cout<<endl;
  cout<<"Total Cross-section Hadron   = "<<counter_weigth_gen<<" pb"<<endl;
  cout<<"Total Cross-section Detector = "<<counter_weigth_det<<" pb"<<endl;
  cout<<endl;

//fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_selected_gen);
     hist_events->SetBinContent(4,counter_weigth_gen);
     hist_events->SetBinContent(5,counter_selected_det);
     hist_events->SetBinContent(6,counter_weigth_det);
     hist_events->SetBinContent(7,counter_response);
     hist_events->SetBinContent(8,counter_miss);
     hist_events->SetBinContent(9,counter_truth);
     hist_events->SetBinContent(10,counter_truth_missed);
     hist_events->SetBinContent(11,counter_truth_stayed);
     hist_events->SetBinContent(12,counter_truth_fakes);
     hist_events->SetBinContent(13,counter_truth_fakes_slightly);

     //recreating the output file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     // saving response matrixes to file
     resp_delta_phi.Write();
     resp_delta_phi_all.Write();

     resp_delta_phi_deta1.Write();
     resp_delta_phi_deta1_all.Write();
     resp_delta_phi_deta2.Write();
     resp_delta_phi_deta2_all.Write();
     resp_delta_phi_deta3.Write();
     resp_delta_phi_deta3_all.Write();
     resp_delta_phi_deta4.Write();
     resp_delta_phi_deta4_all.Write();

     resp_delta_phi_gap.Write();
     resp_delta_phi_gap_all.Write();

     resp_delta_phi_deta1_gap.Write();
     resp_delta_phi_deta1_gap_all.Write();
     resp_delta_phi_deta2_gap.Write();
     resp_delta_phi_deta2_gap_all.Write();
     resp_delta_phi_deta3_gap.Write();
     resp_delta_phi_deta3_gap_all.Write();
     resp_delta_phi_deta4_gap.Write();
     resp_delta_phi_deta4_gap_all.Write();

     resp_delta_phi_nogap.Write();
     resp_delta_phi_nogap_all.Write();

     resp_delta_phi_deta1_nogap.Write();
     resp_delta_phi_deta1_nogap_all.Write();
     resp_delta_phi_deta2_nogap.Write();
     resp_delta_phi_deta2_nogap_all.Write();
     resp_delta_phi_deta3_nogap.Write();
     resp_delta_phi_deta3_nogap_all.Write();
     resp_delta_phi_deta4_nogap.Write();
     resp_delta_phi_deta4_nogap_all.Write();

     resp_leading_pt_inside_gap.Write();
     resp_leading_pt_inside_gap_all.Write();

     resp_leading_eta_star_inside_gap.Write();
     resp_leading_eta_star_inside_gap_all.Write();

     resp_leading_pt_outside_gap.Write();
     resp_leading_pt_outside_gap_all.Write();

     resp_delta_eta_outside_gap.Write();
     resp_delta_eta_outside_gap_all.Write();

     ak5Gen_delta_phi_truth->Write();
     ak5PF_delta_phi_measured->Write();
     ak5Gen_delta_phi_truth_missed->Write();
     ak5Gen_delta_phi_truth_stayed->Write();
     ak5Gen_delta_phi_truth_fakes->Write();
     ak5PF_delta_phi_miss->Write();
     ak5PF_delta_phi_notmissed->Write();
     hist_events->Write();

     ak5Gen_delta_phi_miss->Write();
     ak5Gen_delta_phi_deta1_miss->Write();
     ak5Gen_delta_phi_deta2_miss->Write();
     ak5Gen_delta_phi_deta3_miss->Write();
     ak5Gen_delta_phi_deta4_miss->Write();
     ak5Gen_delta_phi_gap_miss->Write();
     ak5Gen_delta_phi_deta1_gap_miss->Write();
     ak5Gen_delta_phi_deta2_gap_miss->Write();
     ak5Gen_delta_phi_deta3_gap_miss->Write();
     ak5Gen_delta_phi_deta4_gap_miss->Write();
     ak5Gen_delta_phi_nogap_miss->Write();
     ak5Gen_delta_phi_deta1_nogap_miss->Write();
     ak5Gen_delta_phi_deta2_nogap_miss->Write();
     ak5Gen_delta_phi_deta3_nogap_miss->Write();
     ak5Gen_delta_phi_deta4_nogap_miss->Write();
     ak5Gen_leading_pt_inside_gap_miss->Write();
     ak5Gen_leading_eta_star_inside_gap_miss->Write();
     ak5Gen_leading_pt_outside_gap_miss->Write();
     ak5Gen_delta_eta_outside_gap_miss->Write();

     ak5Gen_delta_phi_truth_fakes_slightly_25->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_26->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_27->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_28->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_29->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_30->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_31->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_32->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_33->Write();
     ak5Gen_delta_phi_truth_fakes_slightly_34->Write();

     //closing the file
     data_output->Close();
     
     //delete variable to avoid memory leak
     //delete(resp_delta_phi);
     //delete(resp_delta_phi_fine);
     delete(ak5Gen_delta_phi_truth);
     delete(ak5PF_delta_phi_measured);
     delete(ak5Gen_delta_phi_truth_missed);
     delete(ak5Gen_delta_phi_truth_stayed);
     delete(ak5Gen_delta_phi_truth_fakes);
     delete(ak5PF_delta_phi_miss);
     delete(ak5PF_delta_phi_notmissed);
     delete(hist_events);

     delete(ak5Gen_delta_phi_miss);
     delete(ak5Gen_delta_phi_deta1_miss);
     delete(ak5Gen_delta_phi_deta2_miss);
     delete(ak5Gen_delta_phi_deta3_miss);
     delete(ak5Gen_delta_phi_deta4_miss);
     delete(ak5Gen_delta_phi_gap_miss);
     delete(ak5Gen_delta_phi_deta1_gap_miss);
     delete(ak5Gen_delta_phi_deta2_gap_miss);
     delete(ak5Gen_delta_phi_deta3_gap_miss);
     delete(ak5Gen_delta_phi_deta4_gap_miss);
     delete(ak5Gen_delta_phi_nogap_miss);
     delete(ak5Gen_delta_phi_deta1_nogap_miss);
     delete(ak5Gen_delta_phi_deta2_nogap_miss);
     delete(ak5Gen_delta_phi_deta3_nogap_miss);
     delete(ak5Gen_delta_phi_deta4_nogap_miss);
     delete(ak5Gen_leading_pt_inside_gap_miss);
     delete(ak5Gen_leading_eta_star_inside_gap_miss);
     delete(ak5Gen_leading_pt_outside_gap_miss);
     delete(ak5Gen_delta_eta_outside_gap_miss);

     delete(ak5Gen_delta_phi_truth_fakes_slightly_25);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_26);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_27);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_28);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_29);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_30);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_31);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_32);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_33);
     delete(ak5Gen_delta_phi_truth_fakes_slightly_34);
}
