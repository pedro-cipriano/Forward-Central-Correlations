// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 30 Nov 2012
//
//get_triggered_events(string *data_in, string data_out, int n_files, string sel_mode = "1vertex", bool detail = false)
//reads the Ntuples and created a root file with cross-section histograms

#include "TROOT.h"
#include "TSystem.h"
#include <TFile.h>
#include <TH1.h>
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

void get_triggered_events(string *data_in, string data_out, int n_files, string trigger_label = "HLT_Jet15U", string sel_mode = "allvertex", bool detail = false, bool test = false)
{

   double pt_min = 35.0;

//output the configuration
   if (detail) { cout<<"Get Triggered Events Configuration"<<endl; }
   if (detail) { cout<<"Trigger :         "<<trigger_label<<endl; }
   if (detail) { cout<<"Selection mode :  "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files : "<<n_files<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :     "<<data_out<<endl; }
   if (detail) { cout<<"Pt Min :          "<<pt_min<<endl; }
   if (detail) { cout<<"Detail level :    "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :       "<<test<<endl; }

//binning
   int all_nbins = 7;
   double all_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};

   int cent_nbins = 7;
   int forw_nbins = 7;

   double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
   double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
   
   int deta_nbins = 4;
   int dphi_nbins = 7;

   double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

   int eta_nbins = 14;
   double eta_bins[15] = {-4.7,-4.2,-3.7,-3.2,-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8,3.2,3.7,4.2,4.7};

//declaring histograms
     TH1D *hist_leading_pt;
     TH1D *hist_leading_pt_fine;
     TH1D *hist_leading_pt_nopre;
     TH1D *hist_leading_pt_fine_nopre;
     TH1D *hist_leading_eta;
     TH1D *hist_leading_eta_nopre;
     TH1D *hist_leading_phi;
     TH1D *hist_leading_phi_nopre;
     TH1D *hist_leading_central_pt;
     TH1D *hist_leading_forward_pt;
     TH1D *hist_leading_central_phi;
     TH1D *hist_leading_forward_phi;
     TH1D *hist_leading_central_eta;
     TH1D *hist_leading_forward_eta;
     TH1D *hist_delta_phi;
     TH1D *hist_delta_eta;
     TH1D *hist_leading_central_pt_nopre;
     TH1D *hist_leading_forward_pt_nopre;
     TH1D *hist_leading_central_phi_nopre;
     TH1D *hist_leading_forward_phi_nopre;
     TH1D *hist_leading_central_eta_nopre;
     TH1D *hist_leading_forward_eta_nopre;
     TH1D *hist_delta_phi_nopre;
     TH1D *hist_delta_eta_nopre;
     
     TH1D *hist_sj_leading_pt;
     TH1D *hist_sj_leading_pt_fine;
     TH1D *hist_sj_leading_pt_nopre;
     TH1D *hist_sj_leading_pt_fine_nopre;
     TH1D *hist_sj_leading_eta;
     TH1D *hist_sj_leading_eta_nopre;
     TH1D *hist_sj_leading_phi;
     TH1D *hist_sj_leading_phi_nopre;
     TH1D *hist_sj_leading_central_pt;
     TH1D *hist_sj_leading_forward_pt;
     TH1D *hist_sj_leading_central_phi;
     TH1D *hist_sj_leading_forward_phi;
     TH1D *hist_sj_leading_central_eta;
     TH1D *hist_sj_leading_forward_eta;
     TH1D *hist_sj_leading_central_pt_nopre;
     TH1D *hist_sj_leading_forward_pt_nopre;
     TH1D *hist_sj_leading_central_phi_nopre;
     TH1D *hist_sj_leading_forward_phi_nopre;
     TH1D *hist_sj_leading_central_eta_nopre;
     TH1D *hist_sj_leading_forward_eta_nopre;
     
     TH1D *hist_events;

     hist_leading_pt =  new TH1D("ak5PF_leading_pt","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];# events * prescale", all_nbins, all_bins);
     hist_leading_pt_fine =  new TH1D("ak5PF_leading_pt_fine","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];# events * prescale", 140, 20, 300);
     hist_leading_pt_nopre =  new TH1D("ak5PF_leading_pt_nopre","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_fine_nopre =  new TH1D("ak5PF_leading_pt_fine_nopre","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];# events", 140, 20, 300);
     hist_leading_eta =  new TH1D("ak5PF_leading_eta","Leading Jet #eta;#eta;# events * prescale", eta_nbins, eta_bins);
     hist_leading_eta_nopre =  new TH1D("ak5PF_leading_eta_nopre","Leading Jet #eta with fine binning;#eta;# events", eta_nbins, eta_bins);
     hist_leading_phi =  new TH1D("ak5PF_leading_phi","Leading Jet #phi;#phi;# events * prescale", 14, -3.15, 3.15);
     hist_leading_phi_nopre =  new TH1D("ak5PF_leading_phi_nopre","Leading Jet #phi with fine binning;#phi;# events", 14, -3.15, 3.15);
     hist_leading_forward_pt =  new TH1D("ak5PF_leading_forward_pt","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];# events * prescale", forw_nbins, forw_bins);
     hist_leading_central_pt =  new TH1D("ak5PF_leading_central_pt","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];# events * prescale", cent_nbins, cent_bins);
     hist_leading_central_phi =  new TH1D("ak5PF_leading_central_phi","Leading Central Jet #phi;#phi;# events * prescale", 14, -3.15, 3.15);
     hist_leading_forward_phi =  new TH1D("ak5PF_leading_forward_phi","Leading Forward Jet #phi;#phi;# events * prescale", 14, -3.15, 3.15);
     hist_leading_central_eta =  new TH1D("ak5PF_leading_central_eta","Leading Central Jet #eta;|#eta|;# events * prescale", eta_nbins, eta_bins);
     hist_leading_forward_eta =  new TH1D("ak5PF_leading_forward_eta","Leading Forward Jet #eta;|#eta|;# events * prescale", eta_nbins, eta_bins);
     hist_delta_phi =  new TH1D("ak5PF_delta_phi","#Delta#phi;|#Delta#phi|;# events * prescale", dphi_nbins, dphi_bins);
     hist_delta_eta =  new TH1D("ak5PF_delta_eta","#Delta#eta;|#Delta#eta|;# events * prescale", deta_nbins, deta_bins);
     hist_leading_forward_pt_nopre =  new TH1D("ak5PF_leading_forward_pt_nopre","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];# events", forw_nbins, forw_bins);
     hist_leading_central_pt_nopre =  new TH1D("ak5PF_leading_central_pt_nopre","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];# events", cent_nbins, cent_bins);
     hist_leading_central_phi_nopre =  new TH1D("ak5PF_leading_central_phi_nopre","Leading Central Jet #phi;#phi;# event", 14, -3.15, 3.15);
     hist_leading_forward_phi_nopre =  new TH1D("ak5PF_leading_forward_phi_nopre","Leading Forward Jet #phi;#phi;# events", 14, -3.15, 3.15);
     hist_leading_central_eta_nopre =  new TH1D("ak5PF_leading_central_eta_nopre","Leading Central Jet #eta;|#eta|;# events", eta_nbins, eta_bins);
     hist_leading_forward_eta_nopre =  new TH1D("ak5PF_leading_forward_eta_nopre","Leading Forward Jet #eta;|#eta|;# events", eta_nbins, eta_bins);
     hist_delta_phi_nopre =  new TH1D("ak5PF_delta_phi_nopre","#Delta#phi;|#Delta#phi|;# events", dphi_nbins, dphi_bins);
     hist_delta_eta_nopre =  new TH1D("ak5PF_delta_eta_nopre","#Delta#eta;|#Delta#eta|;# events", deta_nbins, deta_bins);
     
     hist_sj_leading_pt =  new TH1D("ak5PF_sj_leading_pt","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];# events * prescale", all_nbins, all_bins);
     hist_sj_leading_pt_fine =  new TH1D("ak5PF_sj_leading_pt_fine","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];# events * prescale", 140, 20, 300);
     hist_sj_leading_pt_nopre =  new TH1D("ak5PF_sj_leading_pt_nopre","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_fine_nopre =  new TH1D("ak5PF_sj_leading_pt_fine_nopre","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];# events", 140, 20, 300);
     hist_sj_leading_eta =  new TH1D("ak5PF_sj_leading_eta","Leading Jet #eta;#eta;# events * prescale", eta_nbins, eta_bins);
     hist_sj_leading_eta_nopre =  new TH1D("ak5PF_sj_leading_eta_nopre","Leading Jet #eta with fine binning;#eta;# events", eta_nbins, eta_bins);
     hist_sj_leading_phi =  new TH1D("ak5PF_sj_leading_phi","Leading Jet #phi;#phi;# events * prescale", 14, -3.15, 3.15);
     hist_sj_leading_phi_nopre =  new TH1D("ak5PF_sj_leading_phi_nopre","Leading Jet #phi with fine binning;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_forward_pt =  new TH1D("ak5PF_sj_leading_forward_pt","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];# events * prescale", forw_nbins, forw_bins);
     hist_sj_leading_central_pt =  new TH1D("ak5PF_sj_leading_central_pt","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];# events * prescale", cent_nbins, cent_bins);
     hist_sj_leading_central_phi =  new TH1D("ak5PF_sj_leading_central_phi","Leading Central Jet #phi;#phi;# events * prescale", 14, -3.15, 3.15);
     hist_sj_leading_forward_phi =  new TH1D("ak5PF_sj_leading_forward_phi","Leading Forward Jet #phi;#phi;# events * prescale", 14, -3.15, 3.15);
     hist_sj_leading_central_eta =  new TH1D("ak5PF_sj_leading_central_eta","Leading Central Jet #eta;|#eta|;# events * prescale", eta_nbins, eta_bins);
     hist_sj_leading_forward_eta =  new TH1D("ak5PF_sj_leading_forward_eta","Leading Forward Jet #eta;|#eta|;# events * prescale", eta_nbins, eta_bins);
     hist_sj_leading_forward_pt_nopre =  new TH1D("ak5PF_sj_leading_forward_pt_nopre","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];# events", forw_nbins, forw_bins);
     hist_sj_leading_central_pt_nopre =  new TH1D("ak5PF_sj_leading_central_pt_nopre","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];# events", cent_nbins, cent_bins);
     hist_sj_leading_central_phi_nopre =  new TH1D("ak5PF_sj_leading_central_phi_nopre","Leading Central Jet #phi;#phi;# event", 14, -3.15, 3.15);
     hist_sj_leading_forward_phi_nopre =  new TH1D("ak5PF_sj_leading_forward_phi_nopre","Leading Forward Jet #phi;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_central_eta_nopre =  new TH1D("ak5PF_sj_leading_central_eta_nopre","Leading Central Jet #eta;|#eta|;# events", eta_nbins, eta_bins);
     hist_sj_leading_forward_eta_nopre =  new TH1D("ak5PF_sj_leading_forward_eta_nopre","Leading Forward Jet #eta;|#eta|;# events", eta_nbins, eta_bins);
     
     hist_events =  new TH1D("Events","Selection Chain;Type;# Events", 9, 0, 9);

     hist_leading_pt->Sumw2();
     hist_leading_pt_fine->Sumw2();
     hist_leading_pt_nopre->Sumw2();
     hist_leading_pt_fine_nopre->Sumw2();
     hist_leading_eta->Sumw2();
     hist_leading_eta_nopre->Sumw2();
     hist_leading_phi->Sumw2();
     hist_leading_phi_nopre->Sumw2();
     hist_leading_central_pt->Sumw2();
     hist_leading_forward_pt->Sumw2();
     hist_leading_central_phi->Sumw2();
     hist_leading_forward_phi->Sumw2();
     hist_leading_central_eta->Sumw2();
     hist_leading_forward_eta->Sumw2();
     hist_delta_eta->Sumw2();
     hist_delta_phi->Sumw2();
     hist_leading_central_pt_nopre->Sumw2();
     hist_leading_forward_pt_nopre->Sumw2();
     hist_leading_central_phi_nopre->Sumw2();
     hist_leading_forward_phi_nopre->Sumw2();
     hist_leading_central_eta_nopre->Sumw2();
     hist_leading_forward_eta_nopre->Sumw2();
     hist_delta_eta_nopre->Sumw2();
     hist_delta_phi_nopre->Sumw2();
     
     hist_sj_leading_pt->Sumw2();
     hist_sj_leading_pt_fine->Sumw2();
     hist_sj_leading_pt_nopre->Sumw2();
     hist_sj_leading_pt_fine_nopre->Sumw2();
     hist_sj_leading_eta->Sumw2();
     hist_sj_leading_eta_nopre->Sumw2();
     hist_sj_leading_phi->Sumw2();
     hist_sj_leading_phi_nopre->Sumw2();
     hist_sj_leading_central_pt->Sumw2();
     hist_sj_leading_forward_pt->Sumw2();
     hist_sj_leading_central_phi->Sumw2();
     hist_sj_leading_forward_phi->Sumw2();
     hist_sj_leading_central_eta->Sumw2();
     hist_sj_leading_forward_eta->Sumw2();
     hist_sj_leading_central_pt_nopre->Sumw2();
     hist_sj_leading_forward_pt_nopre->Sumw2();
     hist_sj_leading_central_phi_nopre->Sumw2();
     hist_sj_leading_forward_phi_nopre->Sumw2();
     hist_sj_leading_central_eta_nopre->Sumw2();
     hist_sj_leading_forward_eta_nopre->Sumw2();
     
     hist_events->Sumw2();
     hist_events->Fill("Total Events",0);
     hist_events->Fill("PV Selection",0);
     hist_events->Fill("Trigger Selection",0);
     hist_events->Fill("Single Jet Selected",0);
     hist_events->Fill("Single Jet Selected with prescales",0);
     hist_events->Fill("Selected",0);
     hist_events->Fill("Selected with prescales",0);
     hist_events->Fill("Number of Jets",0);
     hist_events->Fill("Leading in the H-HF gap",0);
     
//declaring variables
  int counter_hlt(0),counter_pv(0), counter_entries(0), counter_sj(0), counter_selected(0), counter_jet(0), counter_gap(0);
  double counter_weigth = 0.0, counter_wsj = 0.0;
  int nentries = 0;
  double pu_scale = 0.0;
  double eff = 0.0;
  double prescale = 0.0; // pre1 = 0.0, pre2 = 0.0;
  bool pv_pass = false;
  double leading_pt = 0.0, leading_eta = 0.0, leading_phi = 0.0, pt = 0.0, eta = 0.0, phi = 0.0;
  double pt_central = 0.0, pt_forward = 0.0, eta_central = 0.0, eta_forward = 0.0, phi_central = 0.0, phi_forward = 0.0;
  int ihlt(-1);

//setup for the special case with all triggers
  int number_triggers = 9;
  string trigger_list[number_triggers];
  int hlt_list[number_triggers];
  double probs[number_triggers];
  double prescales[number_triggers];
  bool hltPass(false);
  bool trig_pass[number_triggers];
  trigger_list[0] = "HLT_Jet15U";
  trigger_list[1] = "HLT_Jet30U";
  trigger_list[2] = "HLT_Jet50U";
  trigger_list[3] = "HLT_DiJetAve15U_8E29";
  trigger_list[4] = "HLT_DiJetAve30U_8E29";
  trigger_list[5] = "HLT_FwdJet20U";
  trigger_list[6] = "HLT_DoubleJet15U_ForwardBackward";
  trigger_list[7] = "HLT_L1Jet10U";
  trigger_list[8] = "HLT_L1Jet6U";

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
  TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
  cout<<"Finding trigger mapping: "<<endl;
  
  //----------- loop over the X-axis labels -----------------
  if (trigger_label == "ALL" or trigger_label == "HLT_Jet15U_AND_HLT_Jet30U" or trigger_label == "HLT_Jet15U_AND_HLT_Jet50U" or trigger_label == "HLT_L1Jet6U_AND_HLT_Jet15U" or trigger_label == "HLT_L1Jet10U_AND_HLT_Jet15U" or trigger_label == "HLT_L1Jet6U_AND_HLT_Jet30U" or trigger_label == "HLT_L1Jet10U_AND_HLT_Jet30U" or trigger_label == "HLT_L1Jet6U_AND_HLT_Jet50U" or trigger_label == "HLT_L1Jet10U_AND_HLT_Jet50U")
  {
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
  else
  {
  for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
    TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
    if (ss == trigger_label) {
      ihlt = ibin;
      continue;
    }
  }
  if (ihlt == -1) {
    cout<<"The requested trigger ("<<trigger_label<<") is not found." << endl;
    //break;
  }
  else {
    cout<<trigger_label<<" --> "<<ihlt<<endl;
  }
  }

  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  int decade = 0;

  if (test) { nentries = 10000; } //reduced number of read entries, usefull for testing

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
      if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) { pv_pass = true; }
      if (sel_mode == "1vertex"  && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)    	 { pv_pass = true; }
      
      if (pv_pass == true)
      {
        counter_pv++; 

	hltPass = false;
        prescale = 0.0;
        if (trigger_label == "ALL")
		{
	        for (int trg = 0; trg < 3; trg++)
			{
			prescales[trg] = 1.0;
			trig_pass[trg] = false;
			if (hlt_list[trg] == -1)
				hltPass = true; // no trigger set
        		else
        			{
				if (Event->fired(hlt_list[trg]) > 0)
					{
					hltPass = true;
					trig_pass[trg] = true;
					prescales[trg] = Event->preL1(hlt_list[trg]) * Event->preHLT(hlt_list[trg]);
					}
				}
			}
		if (hltPass)
			{
			for(int l = 0; l < 3; l++)
				{
				probs[l] = 1.0;
				if (trig_pass[l]) { probs[l] = (1 - 1/prescales[l]); }
				}
     			prescale = 1.0 / ( 1 - (probs[0] * probs[1] * probs[2]));
     			}
     		else	{ prescale = 1.0; }
		}
	if (trigger_label == "HLT_L1Jet6U_AND_HLT_Jet15U")
	{
        	if (hlt_list[0] == -1 or hlt_list[8] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[0]) > 0 and Event->fired(hlt_list[8]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[0]) * Event->preHLT(hlt_list[0]);
				// pre1 = Event->preL1(hlt_list[0]) * Event->preHLT(hlt_list[0]);
				// pre2 = Event->preL1(hlt_list[8]) * Event->preHLT(hlt_list[8]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
	if (trigger_label == "HLT_L1Jet10U_AND_HLT_Jet15U")
	{
        	if (hlt_list[0] == -1 or hlt_list[7] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[0]) > 0 and Event->fired(hlt_list[7]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[0]) * Event->preHLT(hlt_list[0]);
				// pre1 = Event->preL1(hlt_list[0]) * Event->preHLT(hlt_list[0]);
				// pre2 = Event->preL1(hlt_list[7]) * Event->preHLT(hlt_list[7]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
	if (trigger_label == "HLT_L1Jet6U_AND_HLT_Jet30U")
	{
        	if (hlt_list[1] == -1 or hlt_list[8] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[1]) > 0 and Event->fired(hlt_list[8]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[1]) * Event->preHLT(hlt_list[1]);
				// pre1 = Event->preL1(hlt_list[1]) * Event->preHLT(hlt_list[1]);
				// pre2 = Event->preL1(hlt_list[8]) * Event->preHLT(hlt_list[8]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
	if (trigger_label == "HLT_L1Jet10U_AND_HLT_Jet30U")
	{
        	if (hlt_list[1] == -1 or hlt_list[7] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[1]) > 0 and Event->fired(hlt_list[7]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[1]) * Event->preHLT(hlt_list[1]);
				// pre1 = Event->preL1(hlt_list[1]) * Event->preHLT(hlt_list[1]);
				// pre2 = Event->preL1(hlt_list[7]) * Event->preHLT(hlt_list[7]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
        if (trigger_label == "HLT_Jet15U_AND_HLT_Jet30U")
	{
        	if (hlt_list[0] == -1 or hlt_list[1] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[0]) > 0 and Event->fired(hlt_list[1]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[1]) * Event->preHLT(hlt_list[1]);
				// pre1 = Event->preL1(hlt_list[1]) * Event->preHLT(hlt_list[1]);
				// pre2 = Event->preL1(hlt_list[0]) * Event->preHLT(hlt_list[0]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
	if (trigger_label == "HLT_L1Jet6U_AND_HLT_Jet50U")
	{
        	if (hlt_list[2] == -1 or hlt_list[8] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[2]) > 0 and Event->fired(hlt_list[8]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[2]) * Event->preHLT(hlt_list[2]);
				// pre1 = Event->preL1(hlt_list[2]) * Event->preHLT(hlt_list[2]);
				// pre2 = Event->preL1(hlt_list[8]) * Event->preHLT(hlt_list[8]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
	if (trigger_label == "HLT_L1Jet10U_AND_HLT_Jet50U")
	{
        	if (hlt_list[2] == -1 or hlt_list[7] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[2]) > 0 and Event->fired(hlt_list[7]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[2]) * Event->preHLT(hlt_list[2]);
				// pre1 = Event->preL1(hlt_list[2]) * Event->preHLT(hlt_list[2]);
				// pre2 = Event->preL1(hlt_list[7]) * Event->preHLT(hlt_list[7]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
        if (trigger_label == "HLT_Jet15U_AND_HLT_Jet50U")
	{
        	if (hlt_list[0] == -1 or hlt_list[2] == -1)
			hltPass = true; // no trigger set
        	else {
			if (Event->fired(hlt_list[0]) > 0 and Event->fired(hlt_list[2]) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(hlt_list[2]) * Event->preHLT(hlt_list[2]);
				// pre1 = Event->preL1(hlt_list[2]) * Event->preHLT(hlt_list[2]);
				// pre2 = Event->preL1(hlt_list[0]) * Event->preHLT(hlt_list[0]);
				// prescale = 1.0 / (1.0 - (1.0 - 1.0/pre1) * (1.0 - 1.0/pre2));
			}
        		
		}
	}
	if (trigger_label != "HLT_Jet15U_AND_HLT_Jet30U" and trigger_label != "HLT_Jet15U_AND_HLT_Jet50U" and trigger_label != "HLT_L1Jet10U_AND_HLT_Jet15U" and trigger_label != "HLT_L1Jet6U_AND_HLT_Jet15U" and trigger_label != "HLT_L1Jet10U_AND_HLT_Jet30U" and trigger_label != "HLT_L1Jet6U_AND_HLT_Jet30U" and trigger_label != "HLT_L1Jet10U_AND_HLT_Jet50U" and trigger_label != "HLT_L1Jet6U_AND_HLT_Jet50U" and trigger_label != "ALL" and trigger_label != "NONE")
	{
        	if (ihlt == -1)
			hltPass = true; // no trigger set
        	else
		{
			if (Event->fired(ihlt) > 0)
			{
				hltPass = true;
				prescale = Event->preL1(ihlt) * Event->preHLT(ihlt);
			}
        		
		}
	}
        if (trigger_label == "NONE") { prescale = 1.0; hltPass = true; }
     

     if (hltPass)
     {
     counter_hlt++;
     leading_pt = -10.0;
     leading_eta = -10.0;
     leading_phi = -10.0;
     pt_forward = -10.0;
     pt_central = -10.0;
     eta_forward = -10.0;
     eta_central = -10.0;
     phi_forward = -10.0;
     eta_central = -10.0;

    for(unsigned int j=0; j<Event->nPFJets(); j++) {
	pt = Event->pfjet(j).ptCor();
	eta = Event->pfjet(j).eta();
	phi = Event->pfjet(j).phi();
     if (pt >= pt_min && Event->pfjet(j).tightID() ) {
     counter_jet++;
     if (leading_pt < pt) { leading_pt = pt; leading_eta = eta; leading_phi = phi; }
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central)
     { pt_central = pt; eta_central = eta; phi_central = phi; }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward)
     { pt_forward = pt; eta_forward = eta; phi_forward = phi; }
     }
     }
     
     if (leading_pt > pt_min)
     	{
    	if ((leading_eta >= 2.8 && leading_eta <= 3.2) || (leading_eta <= -2.8 && leading_eta >= -3.2))
     		{
     		counter_gap++;
     		}
     		else
     		{
     		counter_sj++;
     		counter_wsj = counter_wsj + prescale;
     		hist_sj_leading_forward_pt->Fill(pt_forward,prescale);
     		hist_sj_leading_central_pt->Fill(pt_central,prescale);
     		hist_sj_leading_forward_phi->Fill(phi_forward,prescale);
    		hist_sj_leading_central_phi->Fill(phi_central,prescale);
     		hist_sj_leading_forward_eta->Fill(eta_forward,prescale);
     		hist_sj_leading_central_eta->Fill(eta_central,prescale);
     		hist_sj_leading_forward_pt_nopre->Fill(pt_forward);
     		hist_sj_leading_central_pt_nopre->Fill(pt_central);
     		hist_sj_leading_forward_phi_nopre->Fill(phi_forward);
     		hist_sj_leading_central_phi_nopre->Fill(phi_central);
     		hist_sj_leading_forward_eta_nopre->Fill(eta_forward);
     		hist_sj_leading_central_eta_nopre->Fill(eta_central);
     		hist_sj_leading_pt->Fill(leading_pt,prescale);
     		hist_sj_leading_pt_fine->Fill(leading_pt,prescale);
     		hist_sj_leading_pt_nopre->Fill(leading_pt);
     		hist_sj_leading_pt_fine_nopre->Fill(leading_pt);
     		hist_sj_leading_eta->Fill(leading_eta,prescale);
     		hist_sj_leading_eta_nopre->Fill(leading_eta);
     		hist_sj_leading_phi->Fill(leading_phi,prescale);
     		hist_sj_leading_phi_nopre->Fill(leading_phi);
     		}
     	}
     
     if (pt_forward > pt_min && pt_central > pt_min)
     {
     counter_selected++;
     hist_leading_forward_pt->Fill(pt_forward,prescale);
     hist_leading_central_pt->Fill(pt_central,prescale);
     hist_leading_forward_phi->Fill(phi_forward,prescale);
     hist_leading_central_phi->Fill(phi_central,prescale);
     hist_leading_forward_eta->Fill(eta_forward,prescale);
     hist_leading_central_eta->Fill(eta_central,prescale);
     hist_delta_eta->Fill(abs(eta_central-eta_forward),prescale);
     hist_delta_phi->Fill(calc_delta_phi(phi_central,phi_forward),prescale);
     hist_leading_forward_pt_nopre->Fill(pt_forward);
     hist_leading_central_pt_nopre->Fill(pt_central);
     hist_leading_forward_phi_nopre->Fill(phi_forward);
     hist_leading_central_phi_nopre->Fill(phi_central);
     hist_leading_forward_eta_nopre->Fill(eta_forward);
     hist_leading_central_eta_nopre->Fill(eta_central);
     hist_delta_eta_nopre->Fill(abs(eta_central-eta_forward));
     hist_delta_phi_nopre->Fill(calc_delta_phi(phi_central,phi_forward));
     hist_leading_pt->Fill(leading_pt,prescale);
     hist_leading_pt_fine->Fill(leading_pt,prescale);
     hist_leading_pt_nopre->Fill(leading_pt);
     hist_leading_pt_fine_nopre->Fill(leading_pt);
     hist_leading_eta->Fill(leading_eta,prescale);
     hist_leading_eta_nopre->Fill(leading_eta);
     hist_leading_phi->Fill(leading_phi,prescale);
     hist_leading_phi_nopre->Fill(leading_phi);
     counter_weigth = counter_weigth + prescale;
     }


     }
     
     } //closing the pv if

     } //closing the event loop
     
     } // closing the file loop

//computing auxiliary outputs
     pu_scale = (double)counter_entries/(double)counter_pv;
     eff = (double)counter_hlt/(double)counter_pv;

//output a summary
     if (detail) { cout<<"Total Number of Events :         "<<counter_entries<<endl; }
     if (detail) { cout<<"Primary Vertex Filter :          "<<counter_pv<<endl; }
     if (detail) { cout<<"Pileup Scale :                   "<<pu_scale<<endl; }
     if (detail) { cout<<"Triggered Events :               "<<counter_hlt<<endl; }
     if (detail) { cout<<"Trigger Efficiency :             "<<eff<<endl; }
     if (detail) { cout<<"Single Jet Selected :            "<<counter_sj<<endl; }
     if (detail) { cout<<"Single Jet Weights :             "<<counter_wsj<<endl; }
     if (detail) { cout<<"Selected :                       "<<counter_selected<<endl; }
     if (detail) { cout<<"Total Weights :                  "<<counter_weigth<<endl; }
     if (detail) { cout<<"Total Jets :                     "<<counter_jet<<endl; }
     if (detail) { cout<<"Leading Jets in H-HF gap:        "<<counter_gap<<endl; }
     
//fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_hlt);
     hist_events->SetBinContent(4,counter_sj);
     hist_events->SetBinContent(5,counter_wsj);
     hist_events->SetBinContent(6,counter_selected);
     hist_events->SetBinContent(7,counter_weigth);
     hist_events->SetBinContent(8,counter_jet);
     hist_events->SetBinContent(9,counter_gap);
     
//Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

//write histograms on file
     hist_leading_pt->Write();
     hist_leading_pt_fine->Write();
     hist_leading_pt_nopre->Write();
     hist_leading_pt_fine_nopre->Write();
     hist_leading_central_pt->Write();
     hist_leading_forward_pt->Write();
     hist_leading_central_phi->Write();
     hist_leading_forward_phi->Write();
     hist_leading_central_eta->Write();
     hist_leading_forward_eta->Write();
     hist_delta_eta->Write();
     hist_delta_phi->Write();
     hist_leading_central_pt_nopre->Write();
     hist_leading_forward_pt_nopre->Write();
     hist_leading_central_phi_nopre->Write();
     hist_leading_forward_phi_nopre->Write();
     hist_leading_central_eta_nopre->Write();
     hist_leading_forward_eta_nopre->Write();
     hist_delta_eta_nopre->Write();
     hist_delta_phi_nopre->Write();
     hist_leading_eta->Write();
     hist_leading_eta_nopre->Write();
     hist_leading_phi->Write();
     hist_leading_phi_nopre->Write();
     
     hist_sj_leading_pt->Write();
     hist_sj_leading_pt_fine->Write();
     hist_sj_leading_pt_nopre->Write();
     hist_sj_leading_pt_fine_nopre->Write();
     hist_sj_leading_central_pt->Write();
     hist_sj_leading_forward_pt->Write();
     hist_sj_leading_central_phi->Write();
     hist_sj_leading_forward_phi->Write();
     hist_sj_leading_central_eta->Write();
     hist_sj_leading_forward_eta->Write();
     hist_sj_leading_central_pt_nopre->Write();
     hist_sj_leading_forward_pt_nopre->Write();
     hist_sj_leading_central_phi_nopre->Write();
     hist_sj_leading_forward_phi_nopre->Write();
     hist_sj_leading_central_eta_nopre->Write();
     hist_sj_leading_forward_eta_nopre->Write();
     hist_sj_leading_eta->Write();
     hist_sj_leading_eta_nopre->Write();
     hist_sj_leading_phi->Write();
     hist_sj_leading_phi_nopre->Write();
     
     hist_events->Write();

//close the output file
     data_output->Close();

//delete the histograms to avoid memory leak
     delete(hist_leading_pt);
     delete(hist_leading_pt_fine);
     delete(hist_leading_pt_nopre);
     delete(hist_leading_pt_fine_nopre);
     delete(hist_leading_eta);
     delete(hist_leading_eta_nopre);
     delete(hist_leading_phi);
     delete(hist_leading_phi_nopre);
     delete(hist_leading_central_pt);
     delete(hist_leading_forward_pt);
     delete(hist_leading_central_phi);
     delete(hist_leading_forward_phi);
     delete(hist_leading_central_eta);
     delete(hist_leading_forward_eta);
     delete(hist_delta_eta);
     delete(hist_delta_phi);
     delete(hist_leading_central_pt_nopre);
     delete(hist_leading_forward_pt_nopre);
     delete(hist_leading_central_phi_nopre);
     delete(hist_leading_forward_phi_nopre);
     delete(hist_leading_central_eta_nopre);
     delete(hist_leading_forward_eta_nopre);
     delete(hist_delta_eta_nopre);
     delete(hist_delta_phi_nopre);

     delete(hist_sj_leading_pt);
     delete(hist_sj_leading_pt_fine);
     delete(hist_sj_leading_pt_nopre);
     delete(hist_sj_leading_pt_fine_nopre);
     delete(hist_sj_leading_eta);
     delete(hist_sj_leading_eta_nopre);
     delete(hist_sj_leading_phi);
     delete(hist_sj_leading_phi_nopre);
     delete(hist_sj_leading_central_pt);
     delete(hist_sj_leading_forward_pt);
     delete(hist_sj_leading_central_phi);
     delete(hist_sj_leading_forward_phi);
     delete(hist_sj_leading_central_eta);
     delete(hist_sj_leading_forward_eta);
     delete(hist_sj_leading_central_pt_nopre);
     delete(hist_sj_leading_forward_pt_nopre);
     delete(hist_sj_leading_central_phi_nopre);
     delete(hist_sj_leading_forward_phi_nopre);
     delete(hist_sj_leading_central_eta_nopre);
     delete(hist_sj_leading_forward_eta_nopre);
     
     delete(hist_events);
}
