#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <TLegend.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TF1.h>
#include <TMath.h>
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>
#include <vector>
#include <string>

#include "KKousour/QCDAnalysis/interface/QCDEvent.h"

using namespace std;

void cross_sec_mc_match(string *data_in, string data_out, double *data_lumi, int n_files)
{ 
// total lumi 2010: 35.243pb

   //main cuts
double pt_min = 35.0;
double gap_req = 10.0;

double delta_eta_match = 0.5;

//double cleaning_cut = 10.0;
 
     double pt_central, eta_central, eta_central2, phi_central; // chm_central, elm_central;
     //double et_central, energy_central, mass_central, y_central, area_central; 
     double pt_forward, eta_forward, eta_forward2, phi_forward; // chm_forward, elm_forward;
     //double et_forward, energy_forward, mass_forward, y_forward, area_forward;
     double pt2_central, eta2_central, eta2_central2, phi2_central;
     double pt2_forward, eta2_forward, eta2_forward2, phi2_forward;
     double pt, eta, phi; // chm, elm; // multiplicity;
     double pt_forward_gen, pt_forward_det, pt_central_gen, pt_central_det;
     //double unc, pt_up, pt_down;
     //double et, energy, mass, y, area;
     double eta_for_gen, eta_cen_gen, eta_for_det, eta_cen_det;
     double delta_eta, delta_phi_gen, delta_phi_det, delta_phi_aux, delta_phi2, delta_pt;
     int pass_gap;
     double pt_total_gap, pt_leading_gap, eta_gap, eta_star_inside, phi_gap; // chm_gap, elm_gap;
     double pt_total_outside, pt_leading_outside, eta_outside, phi_outside; // chm_outside;
     //double elm_outside
     double deta_out1, deta_out2;
     double integral, val;        

   //int selected = 0;
   int selected_gap = 0;
   int selected_nogap = 0;
 
int deta_nbins = 4;
int dphi_nbins = 7;
int dphi_nbins_match = 8;
int dphi_nbins_det = 15;
int dphi_nbins_det_cut = 16;
int dphi_nbins_gen_cut = 9;

double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};
double dphi_bins_match[9] = {-1,0.00, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};
double dphi_bins_det[16] = {-1, 0.00, 0.225, 0.45, 0.675, 0.9, 1.125, 1.35, 1.575, 1.8, 2.025, 2.25, 2.475, 2.7, 2.925, 3.15};
double dphi_bins_gen_cut[10] = {0.00, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15, 3.55, 3.95};
double dphi_bins_det_cut[17] = {0.00, 0.225, 0.45, 0.675, 0.9, 1.125, 1.35, 1.575, 1.8, 2.025, 2.25, 2.475, 2.7, 2.925, 3.15, 3.55, 3.95};

int all_nbins = 11;
int all2_nbins = 13;
int in_nbins = 11;
int out_nbins = 11;
int dpt_nbins = 10;

double all_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double all2_bins[14] = {-10, 5, 10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double in_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
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

  TH2D *ak5Match_delta_phi;
  TH2D *ak5Match_delta_phi_matched;
  TH2D *ak5Match_delta_phi_norm;
  TH2D *ak5Match_delta_phi_matched_norm;
  TH2D *ak5Match_delta_phi_unfold;
  TH2D *ak5Match_delta_phi_unfold_norm;
  TH2D *ak5Match_delta_phi_unfold_cut;
  TH2D *ak5Match_delta_phi_unfold_cut_norm;
  
  TH1D *ak5Gen_delta_phi;
  TH1D *ak5Gen_delta_phi_matched;
  TH1D *ak5Gen_delta_phi_cut;

  ak5Match_delta_phi =  new TH2D("ak5Match_delta_phi","#Delta#phi;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events", dphi_nbins_match, dphi_bins_match, dphi_nbins_match, dphi_bins_match);
  ak5Match_delta_phi_matched =  new TH2D("ak5Match_delta_phi_matched","#Delta#phi matched;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events matched", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
  ak5Match_delta_phi_norm =  new TH2D("ak5Match_delta_phi_norm","#Delta#phi;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events", dphi_nbins_match, dphi_bins_match, dphi_nbins_match, dphi_bins_match);
  ak5Match_delta_phi_matched_norm =  new TH2D("ak5Match_delta_phi_matched_norm","#Delta#phi matched;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events matched", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
  ak5Match_delta_phi_unfold =  new TH2D("ak5Match_delta_phi_unfold","#Delta#phi;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events", dphi_nbins_match, dphi_bins_match, dphi_nbins_det, dphi_bins_det);
  ak5Match_delta_phi_unfold_norm =  new TH2D("ak5Match_delta_phi_unfold_norm","#Delta#phi matched;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events matched", dphi_nbins_match, dphi_bins_match, dphi_nbins_det, dphi_bins_det);
  ak5Match_delta_phi_unfold_cut =  new TH2D("ak5Match_delta_phi_unfold_cut","#Delta#phi;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events", dphi_nbins_gen_cut, dphi_bins_gen_cut, dphi_nbins_det_cut, dphi_bins_det_cut);
  ak5Match_delta_phi_unfold_cut_norm =  new TH2D("ak5Match_delta_phi_unfold_cut_norm","#Delta#phi matched;|#Delta#phi| Hadron;|#Delta#phi| Detector;Events matched", dphi_nbins_gen_cut, dphi_bins_gen_cut, dphi_nbins_det_cut, dphi_bins_det_cut);
    
  ak5Gen_delta_phi =  new TH1D("ak5Gen_delta_phi","#Delta#phi;|#Delta#phi|;Events", dphi_nbins_match, dphi_bins_match);
  ak5Gen_delta_phi_matched =  new TH1D("ak5Gen_delta_phi_matched","#Delta#phi;|#Delta#phi|;Events", dphi_nbins, dphi_bins);
  ak5Gen_delta_phi_cut =  new TH1D("ak5Gen_delta_phi_cut","#Delta#phi;|#Delta#phi|;Events", dphi_nbins_gen_cut, dphi_bins_gen_cut);

  ak5Match_delta_phi->Sumw2();
  ak5Match_delta_phi_matched->Sumw2();
  ak5Match_delta_phi_norm->Sumw2();
  ak5Match_delta_phi_matched_norm->Sumw2();
  ak5Match_delta_phi_unfold->Sumw2();
  ak5Match_delta_phi_unfold_norm->Sumw2();
  
  ak5Gen_delta_phi->Sumw2();
  ak5Gen_delta_phi_matched->Sumw2();
  ak5Gen_delta_phi_cut->Sumw2();

  int counter_hlt(0),counter_pv(0),counter_jet(0),counter_selected_gen(0),counter_selected_det(0),counter_entries(0);
  double counter_weigth = 0;
  int nentries = 0;
  double xsec_event = 0.0;

//single trigger
//TString HLT("No_Trigger"); //general not trigger
//TString HLT("HLT_Jet15U"); //2010

// multiple triggers
TString HLT1("HLT_Jet15U"); //2010
TString HLT2("HLT_Jet30U"); //2010
TString HLT3("HLT_Jet50U"); //2010
TString HLT4("HLT_DiJetAve15U_8E29"); //2010
TString HLT5("HLT_DiJetAve30U_8E29"); //2010
//TString HLT6("HLT_FwdJet20U"); //2010

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

  TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
  int ihlt1(-1),ihlt2(-1),ihlt3(-1), ihlt4(-1),ihlt5(-1); //ihlt6(-1);
  cout<<"Finding trigger mapping: "<<endl;
  //----------- loop over the X-axis labels -----------------
  for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT1) {
      ihlt1 = ibin;
      continue;
    }
  }
  if (ihlt1 == -1) {
    cout<<"The requested trigger ("<<HLT1<<") is not found ";
    //break;
  }
  else {
    cout<<HLT1<<" --> "<<ihlt1<<endl;
  }

for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT2) {
      ihlt2 = ibin;
      continue;
    }
  }
  if (ihlt2 == -1) {
    cout<<"The requested trigger ("<<HLT2<<") is not found ";
    //break;
  }
  else {
    cout<<HLT2<<" --> "<<ihlt2<<endl;
  }

for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT3) {
      ihlt3 = ibin;
      continue;
    }
  }
  if (ihlt3 == -1) {
    cout<<"The requested trigger ("<<HLT3<<") is not found ";
    //break;
  }
  else {
    cout<<HLT3<<" --> "<<ihlt3<<endl;
  }

for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT4) {
      ihlt4 = ibin;
      continue;
    }
  }
  if (ihlt4 == -1) {
    cout<<"The requested trigger ("<<HLT4<<") is not found ";
    //break;
  }
  else {
    cout<<HLT4<<" --> "<<ihlt4<<endl;
  }

for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT5) {
      ihlt5 = ibin;
      continue;
    }
  }
  if (ihlt5 == -1) {
    cout<<"The requested trigger ("<<HLT5<<") is not found ";
    //break;
  }
  else {
    cout<<HLT5<<" --> "<<ihlt5<<endl;
  }

/*
for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT6) {
      ihlt6 = ibin;
      continue;
    }
  }
  if (ihlt6 == -1) {
    cout<<"The requested trigger ("<<HLT6<<") is not found ";
    //break;
  }
  else {
    cout<<HLT6<<" --> "<<ihlt6<<endl;
  } */

//trigger studies
/*TString HLT_monitor("HLT_L1Jet6U"); //2010
int ihlt_monitor = -1;

  for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT_monitor) {
      ihlt_monitor = ibin;
      continue;
    }
  }
  if (ihlt_monitor == -1) {
    cout<<"The requested trigger ("<<HLT_monitor<<") is not found ";
    //break;
  }
  else {
    cout<<HLT_monitor<<" --> "<<ihlt_monitor<<endl;
  } */

  /*
TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
  int ihlt(-1);
  cout<<"Finding trigger mapping: "<<endl;
  //----------- loop over the X-axis labels -----------------
  for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == HLT) {
      ihlt = ibin;
      continue;
    }
  }
  if (ihlt == -1) {
    cout<<"The requested trigger ("<<HLT<<") is not found ";
    //break;
  }
  else {
    cout<<HLT<<" --> "<<ihlt<<endl;
  } */
  

  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  xsec_event = data_lumi[z];
  cout<<"Cross-section per event : "<<xsec_event<<endl;
  int decade = 0;

   // nentries = 10000;

   for (Int_t i=0;i<nentries;i++) {
     counter_entries++; 
     double progress = 20.0*i/(1.0*nentries);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
      cout<<5*k<<" % "<<endl;
     decade = k; 
   
     tr->GetEntry(i);
     //single trigger
     //bool hltPass(false);
     //int prescale(1);
     /* if (ihlt == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt) > 0) {
        hltPass = true;
        prescale = Event->preL1(ihlt) * Event->preHLT(ihlt);
      }
     } */

      //-------- check if the primary vertex is good ----  && Event->evtHdr().nVtxGood() == 1
      if (Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1) {
        counter_pv++; 


     // multiple triggers
     double prescales[5];
     double probs[5];
     double prescale = 0.0;
     bool hltPass(false);
     prescales[0] = 0.0;
     if (ihlt1 == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt1) > 0) {
        hltPass = true;
        prescales[0] = Event->preL1(ihlt1) * Event->preHLT(ihlt1);
      }
     }

    prescales[1] = 0.0;
     if (ihlt2 == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt2) > 0) {
        hltPass = true;
        prescales[1] = Event->preL1(ihlt2) * Event->preHLT(ihlt2);
      }
     }

     prescales[2] = 0.0;
     if (ihlt3 == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt3) > 0) {
        hltPass = true;
        prescales[2] = Event->preL1(ihlt3) * Event->preHLT(ihlt3);
      }
     }

     prescales[3] = 0.0;
     if (ihlt4 == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt4) > 0) {
        hltPass = true;
        prescales[3] = Event->preL1(ihlt4) * Event->preHLT(ihlt4);
      }
     }

     prescales[4] = 0.0;
     if (ihlt5 == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt5) > 0) {
        hltPass = true;
        prescales[4] = Event->preL1(ihlt5) * Event->preHLT(ihlt5);
      }
     }

     /* prescales[5] = 0;
     if (ihlt6 == -1)
      hltPass = true; // no trigger set
     else {
      if (Event->fired(ihlt6) > 0) {
        hltPass = true;
        prescales[5] = Event->preL1(ihlt6) * Event->preHLT(ihlt6);
      }
     } */

     //trigger studies
     /*bool hltPass_monitor(false);
     int prescale_monitor = 0;
     if (ihlt_monitor == -1)
      hltPass_monitor = true; // no trigger set
     else {
      if (Event->fired(ihlt_monitor) > 0) {
        hltPass_monitor = true;
        prescale_monitor = Event->preL1(ihlt_monitor) * Event->preHLT(ihlt_monitor);
      }
     }

     prescale = 1.0;*/

     if (hltPass)
     {
     for(unsigned int l=0;l<5;l++) {
     probs[l] = 1;
     if (prescales[l] > 0) { probs[l] = (1 - 1/prescales[l]); }
     //if (prescales[k] > 0) { cout<<"k = "<<k<<" | prob = "<<probs[k]<<" | prescale = "<<prescales[k]<<endl; }
     }
     prescale = 1.0 / ( 1 - (probs[0] * probs[1] * probs[2] * probs[3] * probs[4]));
     //if (prescale > 1) { cout<<"probs = ["<<probs[0]<<" ; "<<probs[1]<<" ; "<<probs[2]<<" ; "<<probs[3]<<" ; "<<probs[4]<<"]"<<endl;   
     //if (prescale > 1) {cout<<"Prescale = "<<prescale<<endl; }
     }
     /*if (hltPass && hltPass_monitor)
     {
     prescale = 1.0 / (1 - (1 - 1/prescale_monitor) * (1 - 1/prescale));  
     }
     if (!hltPass && hltPass_monitor)
     else
     {
     prescale = prescale_monitor;
     }
     //if (prescale  > 0) { cout<<"Prescale = "<<prescale<<endl; }*/
     
     eta_for_gen = 0.0;
     eta_cen_gen = 0.0;
     eta_for_det = 0.0;
     eta_cen_det = 0.0;
     
     delta_phi_gen = -0.9;
     pt_forward_gen = 0.0;
     pt_forward = 0.0;
     eta_forward = 0.0;
     eta_forward2 = 0.0;
     phi_forward = 0.0;
     pt_central = 0.0;
     pt_central_gen = 0.0;
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
     pass_gap = 0;
     pt_total_gap = 0.0;
     eta_gap = -10.0;
     phi_gap = -10.0;
//     chm_gap = -10.0;
//     elm_gap = -10.0;
     pt_total_outside = 0.0;
     deta_out1 = -1.0;
     deta_out2 = -1.0;
     eta_outside = -10.0;
     phi_outside = -10.0;
//     chm_outside = -10.0;
//     elm_outside = -10.0;
     pt_leading_gap = 0.0;
     pt_leading_outside = 0.0;
     eta_star_inside = 0;
//     chm_central = 0;
//     chm_forward = 0;
//     elm_central = 0;
//     elm_forward = 0;

    // cout<<prescale<<endl;

     if (prescale > 0) { //or prescale > 0 hltPass_monitor
      counter_hlt++;
      //cout<<"Prescale = "<<prescale<<endl;
     
    for(unsigned int j=0;j<Event->nGenJets();j++) {
//     beta = Event->genjet(j).beta();
//     chm = Event->genjet(j).chm();
//     elm = Event->genjet(j).elm();
     pt = Event->genjet(j).pt();
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
//     unc = Event->genjet(j).unc();
//     pt_up = pt * (1+unc);
//     pt_down = pt * (1-unc);
     //pt = pt_up;
//     cout<<"Event = "<<i<<" id jet = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" pt up = "<<pt_up<<" pt down = "<<pt_down<<endl;
//     ak5Gen_chm_all->Fill(chm,xsec_event);
//     ak5Gen_elm_all->Fill(elm,xsec_event);
//     ak5Gen_pt_all->Fill(pt,xsec_event);
//     ak5Gen_eta_all->Fill(eta,xsec_event);
//     ak5Gen_phi_all->Fill(phi,xsec_event);
     if (pt >= pt_min) {
     //cout<<"id = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" correction = "<<cor<<endl;
     counter_jet++;
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central)
     { pt_central = pt; eta_central = eta; phi_central = phi; } // chm_central = chm; elm_central = elm; }
     // if (eta <= 2.8 && eta >= -2.8 && pt > pt_min)
     // { ak5Gen_inclusive_central_pt->Fill(pt,prescale); }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward)
     { pt_forward = pt; eta_forward = eta; phi_forward = phi; } // chm_forward = chm; elm_forward = elm; }
     // if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_min)
     // { ak5Gen_inclusive_forward_pt->Fill(pt,prescale); }
     }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt < pt_forward && pt > pt2_forward)
     { pt2_forward = pt; eta2_forward = eta; phi2_forward = phi;}
     if (eta <= 2.8 && eta >= -2.8 && pt > pt2_central && pt < pt_central)
     { pt2_central = pt; eta2_central = eta; phi2_central = phi;}
     }
     pt_forward_gen = pt_forward;
     pt_central_gen = pt_central;
     // if (pt_central > pt_min) { ak5Gen_inclusive_leading_central_pt->Fill(pt_central,prescale); }
     // if (pt_forward > pt_min) { ak5Gen_inclusive_leading_forward_pt->Fill(pt_forward,prescale); }

     if (pt_forward > pt_min && pt_central > pt_min)
     {
     
     // ak5Gen_central_pt_rel->Fill(pt_central,pt2_central,prescale);
     // ak5Gen_forward_pt_rel->Fill(pt_forward,pt2_forward,prescale);
     if (pt_central > 20.0 and pt2_central > 0.0)
     {
     delta_phi2 = phi_central - phi2_central;
     if (delta_phi2 < 0) { delta_phi2 = -delta_phi2; }
     if (delta_phi2 > 3.14159265) {delta_phi_aux = delta_phi2 - 3.14159265; delta_phi2 = 3.14159265 - delta_phi_aux; }
     // ak5Gen_delta_phi_central_rel->Fill(delta_phi2,prescale);
     delta_pt = pt_central - pt2_central;
     // ak5Gen_delta_pt_central_rel->Fill(delta_pt,prescale);
     // if (delta_pt < 5.0) { ak5Gen_delta_phi_central_rel_small->Fill(delta_phi2,prescale); }
     // if (delta_pt > 5.0 && delta_pt < 15.0) { ak5Gen_delta_phi_central_rel_medium->Fill(delta_phi2,prescale); }
     // if (delta_pt > 15.0) { ak5Gen_delta_phi_central_rel_large->Fill(delta_phi2,prescale); }
     }
     if (pt_forward > 20.0 and pt2_forward > 0.0)
     {
     delta_phi2 = phi_forward - phi2_forward;
     if (delta_phi2 < 0) { delta_phi2 = -delta_phi2; }
     if (delta_phi2 > 3.14159265) {delta_phi_aux = delta_phi2 - 3.14159265; delta_phi2 = 3.14159265 - delta_phi_aux; }
     // ak5Gen_delta_phi_forward_rel->Fill(delta_phi2,prescale);
     delta_pt = pt_forward - pt2_forward;
     // ak5Gen_delta_pt_forward_rel->Fill(delta_pt,prescale);
     // if (delta_pt < 5.0) { ak5Gen_delta_phi_forward_rel_small->Fill(delta_phi2,prescale); }
     // if (delta_pt > 5.0 && delta_pt < 15.0) { ak5Gen_delta_phi_forward_rel_medium->Fill(delta_phi2,prescale); }
     // if (delta_pt > 15.0) { ak5Gen_delta_phi_forward_rel_large->Fill(delta_phi2,prescale); }
     }
     
     //cout<<"id = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" correction = "<<cor<<endl;
     // ak5Gen_leading_central_pt->Fill(pt_central,xsec_event);
     // ak5Gen_leading_forward_pt->Fill(pt_forward,xsec_event);
     // ak5Gen_leading_central_phi->Fill(phi_central,xsec_event);
     // ak5Gen_leading_forward_phi->Fill(phi_forward,xsec_event);
     eta_central2 = eta_central;     
     if (eta_central2 < 0) { eta_central2 = -eta_central2; }
     eta_forward2 = eta_forward;     
     if (eta_forward2 < 0) { eta_forward2 = -eta_forward2; }
     // ak5Gen_leading_central_eta->Fill(eta_central,xsec_event);
     // ak5Gen_leading_forward_eta->Fill(eta_forward,xsec_event);
     // ak5Gen_leading_eta->Fill(eta_central,xsec_event);
     // ak5Gen_leading_eta->Fill(eta_forward,data_lumi[z]);
     // ak5Gen_multiplicity->Fill(Event->nGenJets(),xsec_event);
//     ak5Gen_leading_central_beta->Fill(beta_central,xsec_event);
//     ak5Gen_leading_forward_beta->Fill(beta_forward,xsec_event);
//     ak5Gen_leading_central_chm->Fill(chm_central,xsec_event);
//     ak5Gen_leading_forward_chm->Fill(chm_forward,xsec_event);
//     ak5Gen_leading_central_elm->Fill(elm_central,xsec_event);
//     ak5Gen_leading_forward_elm->Fill(elm_forward,xsec_event);
//     vertex_selected->Fill(Event->evtHdr().nVtxGood(),xsec_event);
     counter_selected_gen++;
     counter_weigth = counter_weigth + xsec_event;
     delta_eta = eta_forward - eta_central;
     if (delta_eta < 0) { delta_eta = -delta_eta; }
     delta_phi_gen = phi_forward - phi_central;
     if (delta_phi_gen < 0) { delta_phi_gen = -delta_phi_gen; }
     if (delta_phi_gen > 3.14159265) {delta_phi_aux = delta_phi_gen - 3.14159265; delta_phi_gen = 3.14159265 - delta_phi_aux; }
     
     eta_cen_gen = eta_central;
     eta_for_gen = eta_forward;
     // ak5Gen_delta_eta->Fill(delta_eta,xsec_event);
     // if (delta_eta >= 0.4 && delta_eta < 2.5) { ak5Gen_delta_phi_deta1->Fill(delta_phi,xsec_event); }
     // if (delta_eta >= 2.5 && delta_eta < 3.5) { ak5Gen_delta_phi_deta2->Fill(delta_phi,xsec_event); }
     // if (delta_eta >= 3.5 && delta_eta < 4.5) { ak5Gen_delta_phi_deta3->Fill(delta_phi,xsec_event); }
     // if (delta_eta >= 4.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta4->Fill(delta_phi,xsec_event); }
    // if (delta_eta >= 4.5 && delta_eta < 5.5) { ak5Gen_delta_phi_deta5->Fill(delta_phi,xsec_event); }
    // if (delta_eta >= 5.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta6->Fill(delta_phi,xsec_event); }
     
  //   cout<<"here4 "<<gap_req<<endl;     

  /*   for(unsigned int j=0;j<Event->nGenJets();j++) {
     pt = Event->genjet(j).pt();
     //unc = Event->genjet(j).unc();
     //pt_up = pt * (1+unc);
     //pt_down = pt * (1-unc);
     //pt = pt_up;
     if (pt >= 5) {
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
  //   chm = Event->genjet(j).chm();
  //   elm = Event->genjet(j).elm();
     if (eta_central > eta_forward && eta > eta_forward && eta < eta_central ) {
        if (pt > pt_leading_gap) {pt_leading_gap = pt; eta_gap = eta; phi_gap = phi; } // chm_gap = chm; elm_gap = elm; }
        pt_total_gap = pt_total_gap + pt;
        }
     if (eta_central < eta_forward && eta < eta_forward && eta > eta_central ) {
        if (pt > pt_leading_gap) {pt_leading_gap = pt; eta_gap = eta; phi_gap = phi; } // chm_gap = chm; elm_gap = elm; }
        pt_total_gap = pt_total_gap + pt;
        }
     if (eta_central > eta_forward && (eta < eta_forward || eta > eta_central) ) {
        if (pt > pt_leading_outside) {pt_leading_outside = pt; eta_outside = eta; phi_outside = phi; } // chm_outside = chm; elm_outside = elm; }
        pt_total_outside = pt_total_outside + pt;
        }
     if (eta_central < eta_forward && (eta > eta_forward || eta < eta_central) ) {
        if (pt > pt_leading_outside) {pt_leading_outside = pt; eta_outside = eta; phi_outside = phi; } // chm_outside = chm; elm_outside = elm; }
        pt_total_outside = pt_total_outside + pt;
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
     if (pt_leading_gap > gap_req )
     {
     pass_gap = 1;
//     cout<<i<<" [1]pt_gap = "<<pt_leading_gap<<endl;
     eta_star_inside = eta_gap - (eta_forward + eta_central)/2;
     ak5Gen_leading_pt_inside_gap->Fill(pt_leading_gap,xsec_event);
     ak5Gen_total_pt_inside_gap->Fill(pt_total_gap,xsec_event);
     ak5Gen_leading_eta_inside_gap->Fill(eta_gap,xsec_event);
     ak5Gen_leading_phi_inside_gap->Fill(phi_gap,xsec_event);
//     ak5Gen_leading_chm_inside_gap->Fill(chm_gap,xsec_event);
//     ak5Gen_leading_elm_inside_gap->Fill(elm_gap,xsec_event);
     ak5Gen_leading_eta_star_inside_gap->Fill(eta_star_inside,xsec_event);
     }
     if (pt_leading_outside > gap_req )
     {
     ak5Gen_leading_pt_outside_gap->Fill(pt_leading_outside,xsec_event);
     ak5Gen_total_pt_outside_gap->Fill(pt_total_outside,xsec_event);
     ak5Gen_leading_eta_outside_gap->Fill(eta_outside,xsec_event);
     ak5Gen_leading_phi_outside_gap->Fill(phi_outside,xsec_event);
//     ak5Gen_leading_chm_outside_gap->Fill(chm_outside,xsec_event);
//     ak5Gen_leading_elm_outside_gap->Fill(elm_outside,xsec_event);
     deta_out1 = eta_central - eta_outside;
     if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     deta_out2 = eta_forward - eta_outside;
     if (deta_out2 < 0) { deta_out2 = -deta_out2; }
     if (deta_out1 < deta_out2) { ak5Gen_delta_eta_outside_gap->Fill(deta_out1,xsec_event); }
     if (deta_out2 < deta_out1) { ak5Gen_delta_eta_outside_gap->Fill(deta_out2,xsec_event); }
    // if ((deta_out1 < deta_out2 && deta_out1 > 7.5) || (deta_out2 < deta_out1 && deta_out2 > 7.5)) { cout<<deta_out1<<" - "<<deta_out2<<endl; }
     }
     if (pass_gap == 0)
     {
     selected_gap = selected_gap + 1;
     ak5Gen_delta_phi_gap->Fill(delta_phi,xsec_event);
     ak5Gen_delta_eta_gap->Fill(delta_eta,xsec_event);
     ak5Gen_leading_central_pt_gap->Fill(pt_central,xsec_event);
     ak5Gen_leading_forward_pt_gap->Fill(pt_forward,xsec_event);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { ak5Gen_delta_phi_deta1_gap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { ak5Gen_delta_phi_deta2_gap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { ak5Gen_delta_phi_deta3_gap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta4_gap->Fill(delta_phi,xsec_event); }
   //  if (delta_eta >= 4.5 && delta_eta < 5.5) { ak5Gen_delta_phi_deta5_gap->Fill(delta_phi,xsec_event); }
   //  if (delta_eta >= 5.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta6_gap->Fill(delta_phi,xsec_event); }
     }
     if (pass_gap == 1)
     {
     selected_nogap = selected_nogap + 1;
   //  cout<<selected_nogap<<" [2] pt_gap = "<<pt_leading_gap<<endl;
     ak5Gen_delta_phi_nogap->Fill(delta_phi,xsec_event);
     ak5Gen_delta_eta_nogap->Fill(delta_eta,xsec_event);
     ak5Gen_leading_central_pt_nogap->Fill(pt_central,xsec_event);
     ak5Gen_leading_forward_pt_nogap->Fill(pt_forward,xsec_event);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { ak5Gen_delta_phi_deta1_nogap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { ak5Gen_delta_phi_deta2_nogap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { ak5Gen_delta_phi_deta3_nogap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta4_nogap->Fill(delta_phi,xsec_event); }
  //   if (delta_eta >= 4.5 && delta_eta < 5.5) { ak5Gen_delta_phi_deta5_nogap->Fill(delta_phi,xsec_event); }
  //   if (delta_eta >= 5.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta6_nogap->Fill(delta_phi,xsec_event); }
     }
     */
     }
     
     delta_phi_det = -.9;
     pt_forward = 0.0;
     pt_forward_det = 0.0;
     eta_forward = 0.0;
     eta_forward2 = 0.0;
     phi_forward = 0.0;
     pt_central = 0.0;
     pt_central_det = 0.0;
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
     pass_gap = 0;
     pt_total_gap = 0.0;
     eta_gap = -10.0;
     phi_gap = -10.0;
//     chm_gap = -10.0;
//     elm_gap = -10.0;
     pt_total_outside = 0.0;
     deta_out1 = -1.0;
     deta_out2 = -1.0;
     eta_outside = -10.0;
     phi_outside = -10.0;
//     chm_outside = -10.0;
//     elm_outside = -10.0;
     pt_leading_gap = 0.0;
     pt_leading_outside = 0.0;
     eta_star_inside = 0;
//     chm_central = 0;
//     chm_forward = 0;
//     elm_central = 0;
//     elm_forward = 0;
     
       for(unsigned int j=0;j<Event->nPFJets();j++) {
//     beta = Event->pfjet(j).beta();
//     chm = Event->pfjet(j).chm();
//     elm = Event->pfjet(j).elm();
     pt = Event->pfjet(j).ptCor();
     eta = Event->pfjet(j).eta();
     phi = Event->pfjet(j).phi();
//     unc = Event->pfjet(j).unc();
//     pt_up = pt * (1+unc);
//     pt_down = pt * (1-unc);
     //pt = pt_up;
//     cout<<"Event = "<<i<<" id jet = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" pt up = "<<pt_up<<" pt down = "<<pt_down<<endl;
//     ak5PF_chm_all->Fill(chm,xsec_event);
//     ak5PF_elm_all->Fill(elm,xsec_event);
//     ak5PF_pt_all->Fill(pt,xsec_event);
//     ak5PF_eta_all->Fill(eta,xsec_event);
//     ak5PF_phi_all->Fill(phi,xsec_event);
     if (pt >= pt_min) {
     //cout<<"id = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" correction = "<<cor<<endl;
     counter_jet++;
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central)
     { pt_central = pt; eta_central = eta; phi_central = phi; } // chm_central = chm; elm_central = elm; }
     // if (eta <= 2.8 && eta >= -2.8 && pt > pt_min)
     // { ak5PF_inclusive_central_pt->Fill(pt,prescale); }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward)
     { pt_forward = pt; eta_forward = eta; phi_forward = phi; } // chm_forward = chm; elm_forward = elm; }
     // if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_min)
     // { ak5PF_inclusive_forward_pt->Fill(pt,prescale); }
     }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt < pt_forward && pt > pt2_forward)
     { pt2_forward = pt; eta2_forward = eta; phi2_forward = phi;}
     if (eta <= 2.8 && eta >= -2.8 && pt > pt2_central && pt < pt_central)
     { pt2_central = pt; eta2_central = eta; phi2_central = phi;}
     }
     
     pt_forward_det = pt_forward;
     pt_central_det = pt_central;
     // if (pt_central > pt_min) { ak5PF_inclusive_leading_central_pt->Fill(pt_central,prescale); }
     // if (pt_forward > pt_min) { ak5PF_inclusive_leading_forward_pt->Fill(pt_forward,prescale); }

     if (pt_forward > pt_min && pt_central > pt_min)
     {
     
     // ak5PF_central_pt_rel->Fill(pt_central,pt2_central,prescale);
     // ak5PF_forward_pt_rel->Fill(pt_forward,pt2_forward,prescale);
     if (pt_central > 20.0 and pt2_central > 0.0)
     {
     delta_phi2 = phi_central - phi2_central;
     if (delta_phi2 < 0) { delta_phi2 = -delta_phi2; }
     if (delta_phi2 > 3.14159265) {delta_phi_aux = delta_phi2 - 3.14159265; delta_phi2 = 3.14159265 - delta_phi_aux; }
     // ak5PF_delta_phi_central_rel->Fill(delta_phi2,prescale);
     delta_pt = pt_central - pt2_central;
     // ak5PF_delta_pt_central_rel->Fill(delta_pt,prescale);
     // if (delta_pt < 5.0) { ak5PF_delta_phi_central_rel_small->Fill(delta_phi2,prescale); }
     // if (delta_pt > 5.0 && delta_pt < 15.0) { ak5Gen_delta_phi_central_rel_medium->Fill(delta_phi2,prescale); }
     // if (delta_pt > 15.0) { ak5Gen_delta_phi_central_rel_large->Fill(delta_phi2,prescale); }
     }
     if (pt_forward > 20.0 and pt2_forward > 0.0)
     {
     delta_phi2 = phi_forward - phi2_forward;
     if (delta_phi2 < 0) { delta_phi2 = -delta_phi2; }
     if (delta_phi2 > 3.14159265) {delta_phi_aux = delta_phi2 - 3.14159265; delta_phi2 = 3.14159265 - delta_phi_aux; }
     // ak5Gen_delta_phi_forward_rel->Fill(delta_phi2,prescale);
     delta_pt = pt_forward - pt2_forward;
     // ak5Gen_delta_pt_forward_rel->Fill(delta_pt,prescale);
     // if (delta_pt < 5.0) { ak5Gen_delta_phi_forward_rel_small->Fill(delta_phi2,prescale); }
     // if (delta_pt > 5.0 && delta_pt < 15.0) { ak5Gen_delta_phi_forward_rel_medium->Fill(delta_phi2,prescale); }
     // if (delta_pt > 15.0) { ak5Gen_delta_phi_forward_rel_large->Fill(delta_phi2,prescale); }
     }
     
     //cout<<"id = "<<j<<" pt = "<<pt<<" eta = "<<eta<<" uncertainty = "<<unc<<" correction = "<<cor<<endl;
     // ak5Gen_leading_central_pt->Fill(pt_central,xsec_event);
     // ak5Gen_leading_forward_pt->Fill(pt_forward,xsec_event);
     // ak5Gen_leading_central_phi->Fill(phi_central,xsec_event);
     // ak5Gen_leading_forward_phi->Fill(phi_forward,xsec_event);
     eta_central2 = eta_central;     
     if (eta_central2 < 0) { eta_central2 = -eta_central2; }
     eta_forward2 = eta_forward;     
     if (eta_forward2 < 0) { eta_forward2 = -eta_forward2; }
     // ak5Gen_leading_central_eta->Fill(eta_central,xsec_event);
     // ak5Gen_leading_forward_eta->Fill(eta_forward,xsec_event);
     // ak5Gen_leading_eta->Fill(eta_central,xsec_event);
     // ak5Gen_leading_eta->Fill(eta_forward,data_lumi[z]);
     // ak5Gen_multiplicity->Fill(Event->nGenJets(),xsec_event);
//     ak5Gen_leading_central_beta->Fill(beta_central,xsec_event);
//     ak5Gen_leading_forward_beta->Fill(beta_forward,xsec_event);
//     ak5Gen_leading_central_chm->Fill(chm_central,xsec_event);
//     ak5Gen_leading_forward_chm->Fill(chm_forward,xsec_event);
//     ak5Gen_leading_central_elm->Fill(elm_central,xsec_event);
//     ak5Gen_leading_forward_elm->Fill(elm_forward,xsec_event);
//     vertex_selected->Fill(Event->evtHdr().nVtxGood(),xsec_event);
     counter_selected_det++;
     counter_weigth = counter_weigth + xsec_event;
     delta_eta = eta_forward - eta_central;
     if (delta_eta < 0) { delta_eta = -delta_eta; }
     delta_phi_det = phi_forward - phi_central;
     if (delta_phi_det < 0) { delta_phi_det = -delta_phi_det; }
     if (delta_phi_det > 3.14159265) {delta_phi_aux = delta_phi_det - 3.14159265; delta_phi_det = 3.14159265 - delta_phi_aux; }
     
     eta_cen_det = eta_central;
     eta_for_det = eta_forward;
     // ak5PF_delta_phi->Fill(delta_phi_det,xsec_event);
     // ak5Gen_delta_eta->Fill(delta_eta,xsec_event);
     // if (delta_eta >= 0.4 && delta_eta < 2.5) { ak5Gen_delta_phi_deta1->Fill(delta_phi,xsec_event); }
     // if (delta_eta >= 2.5 && delta_eta < 3.5) { ak5Gen_delta_phi_deta2->Fill(delta_phi,xsec_event); }
     // if (delta_eta >= 3.5 && delta_eta < 4.5) { ak5Gen_delta_phi_deta3->Fill(delta_phi,xsec_event); }
     // if (delta_eta >= 4.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta4->Fill(delta_phi,xsec_event); }
    // if (delta_eta >= 4.5 && delta_eta < 5.5) { ak5Gen_delta_phi_deta5->Fill(delta_phi,xsec_event); }
    // if (delta_eta >= 5.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta6->Fill(delta_phi,xsec_event); }
     
  //   cout<<"here4 "<<gap_req<<endl;     

  /*   for(unsigned int j=0;j<Event->nGenJets();j++) {
     pt = Event->genjet(j).pt();
     //unc = Event->genjet(j).unc();
     //pt_up = pt * (1+unc);
     //pt_down = pt * (1-unc);
     //pt = pt_up;
     if (pt >= 5) {
     eta = Event->genjet(j).eta();
     phi = Event->genjet(j).phi();
  //   chm = Event->genjet(j).chm();
  //   elm = Event->genjet(j).elm();
     if (eta_central > eta_forward && eta > eta_forward && eta < eta_central ) {
        if (pt > pt_leading_gap) {pt_leading_gap = pt; eta_gap = eta; phi_gap = phi; } // chm_gap = chm; elm_gap = elm; }
        pt_total_gap = pt_total_gap + pt;
        }
     if (eta_central < eta_forward && eta < eta_forward && eta > eta_central ) {
        if (pt > pt_leading_gap) {pt_leading_gap = pt; eta_gap = eta; phi_gap = phi; } // chm_gap = chm; elm_gap = elm; }
        pt_total_gap = pt_total_gap + pt;
        }
     if (eta_central > eta_forward && (eta < eta_forward || eta > eta_central) ) {
        if (pt > pt_leading_outside) {pt_leading_outside = pt; eta_outside = eta; phi_outside = phi; } // chm_outside = chm; elm_outside = elm; }
        pt_total_outside = pt_total_outside + pt;
        }
     if (eta_central < eta_forward && (eta > eta_forward || eta < eta_central) ) {
        if (pt > pt_leading_outside) {pt_leading_outside = pt; eta_outside = eta; phi_outside = phi; } // chm_outside = chm; elm_outside = elm; }
        pt_total_outside = pt_total_outside + pt;
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
     if (pt_leading_gap > gap_req )
     {
     pass_gap = 1;
//     cout<<i<<" [1]pt_gap = "<<pt_leading_gap<<endl;
     eta_star_inside = eta_gap - (eta_forward + eta_central)/2;
     ak5Gen_leading_pt_inside_gap->Fill(pt_leading_gap,xsec_event);
     ak5Gen_total_pt_inside_gap->Fill(pt_total_gap,xsec_event);
     ak5Gen_leading_eta_inside_gap->Fill(eta_gap,xsec_event);
     ak5Gen_leading_phi_inside_gap->Fill(phi_gap,xsec_event);
//     ak5Gen_leading_chm_inside_gap->Fill(chm_gap,xsec_event);
//     ak5Gen_leading_elm_inside_gap->Fill(elm_gap,xsec_event);
     ak5Gen_leading_eta_star_inside_gap->Fill(eta_star_inside,xsec_event);
     }
     if (pt_leading_outside > gap_req )
     {
     ak5Gen_leading_pt_outside_gap->Fill(pt_leading_outside,xsec_event);
     ak5Gen_total_pt_outside_gap->Fill(pt_total_outside,xsec_event);
     ak5Gen_leading_eta_outside_gap->Fill(eta_outside,xsec_event);
     ak5Gen_leading_phi_outside_gap->Fill(phi_outside,xsec_event);
//     ak5Gen_leading_chm_outside_gap->Fill(chm_outside,xsec_event);
//     ak5Gen_leading_elm_outside_gap->Fill(elm_outside,xsec_event);
     deta_out1 = eta_central - eta_outside;
     if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     deta_out2 = eta_forward - eta_outside;
     if (deta_out2 < 0) { deta_out2 = -deta_out2; }
     if (deta_out1 < deta_out2) { ak5Gen_delta_eta_outside_gap->Fill(deta_out1,xsec_event); }
     if (deta_out2 < deta_out1) { ak5Gen_delta_eta_outside_gap->Fill(deta_out2,xsec_event); }
    // if ((deta_out1 < deta_out2 && deta_out1 > 7.5) || (deta_out2 < deta_out1 && deta_out2 > 7.5)) { cout<<deta_out1<<" - "<<deta_out2<<endl; }
     }
     if (pass_gap == 0)
     {
     selected_gap = selected_gap + 1;
     ak5Gen_delta_phi_gap->Fill(delta_phi,xsec_event);
     ak5Gen_delta_eta_gap->Fill(delta_eta,xsec_event);
     ak5Gen_leading_central_pt_gap->Fill(pt_central,xsec_event);
     ak5Gen_leading_forward_pt_gap->Fill(pt_forward,xsec_event);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { ak5Gen_delta_phi_deta1_gap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { ak5Gen_delta_phi_deta2_gap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { ak5Gen_delta_phi_deta3_gap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta4_gap->Fill(delta_phi,xsec_event); }
   //  if (delta_eta >= 4.5 && delta_eta < 5.5) { ak5Gen_delta_phi_deta5_gap->Fill(delta_phi,xsec_event); }
   //  if (delta_eta >= 5.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta6_gap->Fill(delta_phi,xsec_event); }
     }
     if (pass_gap == 1)
     {
     selected_nogap = selected_nogap + 1;
   //  cout<<selected_nogap<<" [2] pt_gap = "<<pt_leading_gap<<endl;
     ak5Gen_delta_phi_nogap->Fill(delta_phi,xsec_event);
     ak5Gen_delta_eta_nogap->Fill(delta_eta,xsec_event);
     ak5Gen_leading_central_pt_nogap->Fill(pt_central,xsec_event);
     ak5Gen_leading_forward_pt_nogap->Fill(pt_forward,xsec_event);
     if (delta_eta >= 0.4 && delta_eta < 2.5) { ak5Gen_delta_phi_deta1_nogap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 2.5 && delta_eta < 3.5) { ak5Gen_delta_phi_deta2_nogap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 3.5 && delta_eta < 4.5) { ak5Gen_delta_phi_deta3_nogap->Fill(delta_phi,xsec_event); }
     if (delta_eta >= 4.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta4_nogap->Fill(delta_phi,xsec_event); }
  //   if (delta_eta >= 4.5 && delta_eta < 5.5) { ak5Gen_delta_phi_deta5_nogap->Fill(delta_phi,xsec_event); }
  //   if (delta_eta >= 5.5 && delta_eta < 7.5) { ak5Gen_delta_phi_deta6_nogap->Fill(delta_phi,xsec_event); }
     }
     */
     }
     
     if ((pt_central_det > pt_min and pt_forward_det > pt_min) or (pt_central_gen > pt_min and pt_forward_gen > pt_min))
     { 
     if (delta_phi_gen == 0.0) { delta_phi_gen = -0.5; }
     if (delta_phi_det == 0.0) { delta_phi_det = -0.5; }
     // cout << "event in : " << delta_phi_gen << ";" << delta_phi_det << endl;
     ak5Gen_delta_phi->Fill(delta_phi_gen,xsec_event); 
     ak5Match_delta_phi->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     ak5Match_delta_phi_norm->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     ak5Match_delta_phi_unfold->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     ak5Match_delta_phi_unfold_norm->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     	if (abs(eta_for_gen - eta_for_det) < delta_eta_match and abs(eta_cen_gen - eta_cen_det) < delta_eta_match)
     	{
     	ak5Gen_delta_phi_matched->Fill(delta_phi_gen,xsec_event);
     	ak5Match_delta_phi_matched->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     	ak5Match_delta_phi_matched_norm->Fill(delta_phi_gen,delta_phi_det,xsec_event); 
     	}
     if (pt_forward_gen < pt_min and delta_phi_gen < 0.0) { delta_phi_gen = 3.4; }
     if (pt_central_gen < pt_min and delta_phi_gen < 0.0) { delta_phi_gen = 3.7; }
     if (pt_forward_det < pt_min and delta_phi_det < 0.0) { delta_phi_det = 3.4; }
     if (pt_central_det < pt_min and delta_phi_det < 0.0) { delta_phi_det = 3.7; }
     ak5Gen_delta_phi_cut->Fill(delta_phi_gen,xsec_event); 
     ak5Match_delta_phi_unfold_cut->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     ak5Match_delta_phi_unfold_cut_norm->Fill(delta_phi_gen,delta_phi_det,xsec_event);
     }
     
     
  } 


  }  }
  
}

  cout<<endl<<endl;
  cout<<"Events read:                      "<<counter_entries<<endl;
  cout<<"Events after the PV cut:          "<<counter_pv<<endl;
  cout<<"Events after the trigger cut:     "<<counter_hlt<<endl;
  cout<<"Events Main Selection Hadron:     "<<counter_selected_gen<<endl;
  cout<<"Events Main Selection Detector:   "<<counter_selected_det<<endl;
  cout<<"Events weight count               "<<counter_weigth<<endl;

//double cross_section = (double) counter_weigth / data_lumi[z];
cout<<endl;
//cout<<"Dataset luminosity =  "<<data_lumi<<" pb^-1"<<endl;
cout<<"Total Cross-section = "<<counter_weigth<<" pb"<<endl;
cout<<endl;


for (int i = 0; i <= ak5Gen_delta_phi->GetNbinsX(); i++)
	{
	cout << i <<  " : " << ak5Gen_delta_phi->GetBinContent(i) << endl;
	if (ak5Gen_delta_phi->GetBinContent(i) > 0)
		{
		for (int j = 0; j <= ak5Match_delta_phi_norm->GetNbinsY(); j++)
			{
			val = ak5Match_delta_phi_norm->GetBinContent(i,j) / ak5Gen_delta_phi->GetBinContent(i);
			cout << i << "-" << j <<  " : " << val << endl;
			ak5Match_delta_phi_norm->SetBinContent(i,j,val);
			}
		}
	}

for (int i = 0; i <= ak5Gen_delta_phi->GetNbinsX(); i++)
	{
	cout << i <<  " : " << ak5Gen_delta_phi->GetBinContent(i) << endl;
	if (ak5Gen_delta_phi->GetBinContent(i) > 0)
		{
		for (int j = 0; j <= ak5Match_delta_phi_unfold_norm->GetNbinsY(); j++)
			{
			val = ak5Match_delta_phi_unfold_norm->GetBinContent(i,j) / ak5Gen_delta_phi->GetBinContent(i);
			cout << i << "-" << j <<  " : " << val << endl;
			ak5Match_delta_phi_unfold_norm->SetBinContent(i,j,val);
			}
		}
	}
	
	
for (int i = 0; i <= ak5Gen_delta_phi_cut->GetNbinsX(); i++)
	{
	cout << i <<  " : " << ak5Gen_delta_phi_cut->GetBinContent(i) << endl;
	if (ak5Gen_delta_phi_cut->GetBinContent(i) > 0)
		{
		for (int j = 0; j <= ak5Match_delta_phi_unfold_cut_norm->GetNbinsY(); j++)
			{
			val = ak5Match_delta_phi_unfold_cut_norm->GetBinContent(i,j) / ak5Gen_delta_phi_cut->GetBinContent(i);
			cout << i << "-" << j <<  " : " << val << endl;
			ak5Match_delta_phi_unfold_cut_norm->SetBinContent(i,j,val);
			}
		}
	}
	
for (int i = 0; i <= ak5Gen_delta_phi_matched->GetNbinsX(); i++)
	{
	cout << i <<  " : " << ak5Gen_delta_phi_matched->GetBinContent(i) << endl;
	if (ak5Gen_delta_phi_matched->GetBinContent(i) > 0)
		{
		for (int j = 0; j <= ak5Match_delta_phi_matched_norm->GetNbinsY(); j++)
			{
			val = ak5Match_delta_phi_matched_norm->GetBinContent(i,j) / ak5Gen_delta_phi_matched->GetBinContent(i);
			cout << i << "-" << j <<  " : " << val << endl;
			ak5Match_delta_phi_matched_norm->SetBinContent(i,j,val);
			}
		}
	}



     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     ak5Match_delta_phi->Write();
     ak5Match_delta_phi_matched->Write();
     ak5Match_delta_phi_norm->Write();
     ak5Match_delta_phi_matched_norm->Write();
     ak5Match_delta_phi_unfold->Write();
     ak5Match_delta_phi_unfold_norm->Write();
     ak5Match_delta_phi_unfold_cut->Write();
     ak5Match_delta_phi_unfold_cut_norm->Write();

     data_output->Close();

}



