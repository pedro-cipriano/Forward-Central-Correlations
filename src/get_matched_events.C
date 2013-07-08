// Pedro Cipriano, Apr 2013
// DESY, CMS
// Last Update: 12 Apr 2012
//
//get_matched_events(string *data_in, string data_out, double *data_lumi, int n_files, string sel_mode = "1vertex", bool detail = false)
// reads the MC Ntuples and and matches the events

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

void get_matched_events(string *data_in, string data_out, double *data_lumi, int n_files, string sel_mode = "1vertex", string vertex_weights = "", TString vertex_sufix = "_v0", bool detail = false, bool test = false)
{

   //main cuts
   double pt_min = 35.0;
   double gap_req = 20.0;
   double match_radius = 0.3;

   //output the configuration
   if (detail) { cout<<"Get Matched Events Configuration"<<endl; }
   if (detail) { cout<<"Selection mode :       "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files :      "<<n_files<<endl; }
   if (detail) { cout<<"Vertex Weights File :  "<<vertex_weights<<endl; }
   if (detail) { cout<<"Vertex Weights Sufix : "<<vertex_sufix<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :          "<<data_out<<endl; }
   if (detail) { cout<<"Detail level :         "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :            "<<test<<endl; }
   if (detail) { cout<<"pT Min :               "<<pt_min<<endl; }
   if (detail) { cout<<"Match Radius :         "<<match_radius<<endl; }
   if (detail) { cout<<"Gap Requirement :      "<<gap_req<<endl; }

   double prescale, vertex_factor;
   double deta_out1, deta_out2;
   double pt, eta, phi, gen_pt, gen_eta, gen_phi;
   double pt_central_det, pt_forward_det, eta_central_det, eta_forward_det, phi_central_det, phi_forward_det;
   double pt_inside_det, eta_inside_det, phi_inside_det, pt_outside_det, eta_outside_det, phi_outside_det;
   double delta_phi_det, delta_eta_det, eta_star_inside_det, deta_out_det;
   double pt_central_gen, pt_forward_gen, eta_central_gen, eta_forward_gen, phi_central_gen, phi_forward_gen;
   double pt_inside_gen, eta_inside_gen, phi_inside_gen, pt_outside_gen, eta_outside_gen, phi_outside_gen;
   double delta_phi_gen, delta_eta_gen, eta_star_inside_gen, deta_out_gen;
   double deta_forward, deta_central, dphi_forward, dphi_central, forward_radius, central_radius;
   double deta_inside, deta_outside, dphi_inside, dphi_outside, inside_radius, outside_radius;

   bool pv_pass_det, pv_pass_gen, pass_tight, pass_gen, pass_det;
   bool pass_match, pass_inside_match, pass_outside_match;
   bool pass_deta1_gen, pass_deta2_gen, pass_deta3_gen, pass_deta4_gen, pass_gap_gen, pass_nogap_gen;
   bool pass_deta1_det, pass_deta2_det, pass_deta3_det, pass_deta4_det, pass_gap_det, pass_nogap_det;
   bool pass_out_gen, pass_out_det;

   int nentries = 0;
   int counter_pv_gen = 0, counter_pv_det = 0, counter_selected_gen = 0, counter_selected_det = 0;
   int counter_selected = 0, counter_matched = 0, counter_entries = 0;

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


   // binning

   int dphi_nbins = 7;
   int dphi_nbins_nomatch = 8;
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};
   double dphi_bins_nomatch[9] = {-1.0, 0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

   int in_nbins = 9;
   int in_nbins_nomatch = 10;
   double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
   double in_bins_nomatch[11] = {0, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

   int out_nbins = 9;
   int out_nbins_nomatch = 10;
   double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
   double out_bins_nomatch[11] = {0, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

   int etastar_nbins = 12;
   int etastar_nbins_nomatch = 13;
   double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};
   double etastar_bins_nomatch[14] = {-10,-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

   int deta_out_nbins = 6;
   int deta_out_nbins_nomatch = 7;
   double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};
   double deta_out_bins_nomatch[8] = {-1.0,0, 1, 2, 3, 4, 5, 7.5};

   //histograms

   TH2D *hist_delta_phi;
   TH2D *hist_delta_phi_deta1;
   TH2D *hist_delta_phi_deta2;
   TH2D *hist_delta_phi_deta3;
   TH2D *hist_delta_phi_deta4;
   TH2D *hist_delta_phi_gap;
   TH2D *hist_delta_phi_deta1_gap;
   TH2D *hist_delta_phi_deta2_gap;
   TH2D *hist_delta_phi_deta3_gap;
   TH2D *hist_delta_phi_deta4_gap;
   TH2D *hist_delta_phi_nogap;
   TH2D *hist_delta_phi_deta1_nogap;
   TH2D *hist_delta_phi_deta2_nogap;
   TH2D *hist_delta_phi_deta3_nogap;
   TH2D *hist_delta_phi_deta4_nogap;
   TH2D *hist_leading_pt_inside_gap;
   TH2D *hist_leading_eta_star_inside_gap;
   TH2D *hist_leading_pt_outside_gap;
   TH2D *hist_delta_eta_outside_gap;

   TH2D *hist_delta_phi_nomatch;
   TH2D *hist_delta_phi_deta1_nomatch;
   TH2D *hist_delta_phi_deta2_nomatch;
   TH2D *hist_delta_phi_deta3_nomatch;
   TH2D *hist_delta_phi_deta4_nomatch;
   TH2D *hist_delta_phi_gap_nomatch;
   TH2D *hist_delta_phi_deta1_gap_nomatch;
   TH2D *hist_delta_phi_deta2_gap_nomatch;
   TH2D *hist_delta_phi_deta3_gap_nomatch;
   TH2D *hist_delta_phi_deta4_gap_nomatch;
   TH2D *hist_delta_phi_nogap_nomatch;
   TH2D *hist_delta_phi_deta1_nogap_nomatch;
   TH2D *hist_delta_phi_deta2_nogap_nomatch;
   TH2D *hist_delta_phi_deta3_nogap_nomatch;
   TH2D *hist_delta_phi_deta4_nogap_nomatch;
   TH2D *hist_leading_pt_inside_gap_nomatch;
   TH2D *hist_leading_eta_star_inside_gap_nomatch;
   TH2D *hist_leading_pt_outside_gap_nomatch;
   TH2D *hist_delta_eta_outside_gap_nomatch;

   hist_delta_phi =  new TH2D("matched_delta_phi","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta1 =  new TH2D("matched_delta_phi_deta1","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta2 =  new TH2D("matched_delta_phi_deta2","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta3 =  new TH2D("matched_delta_phi_deta3","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta4 =  new TH2D("matched_delta_phi_deta4","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_gap =  new TH2D("matched_delta_phi_gap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta1_gap =  new TH2D("matched_delta_phi_deta1_gap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta2_gap =  new TH2D("matched_delta_phi_deta2_gap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta3_gap =  new TH2D("matched_delta_phi_deta3_gap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta4_gap =  new TH2D("matched_delta_phi_deta4_gap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_nogap =  new TH2D("matched_delta_phi_nogap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta1_nogap =  new TH2D("matched_delta_phi_deta1_nogap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta2_nogap =  new TH2D("matched_delta_phi_deta2_nogap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta3_nogap =  new TH2D("matched_delta_phi_deta3_nogap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_delta_phi_deta4_nogap =  new TH2D("matched_delta_phi_deta4_nogap","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins, dphi_nbins, dphi_bins);
   hist_leading_pt_inside_gap =  new TH2D("matched_leading_pt_inside_gap","p_{T};p_{T}^{gen} [GeV];p_{T}^{det} [GeV];#frac{d#sigma}{dp_{T}} [pb]", in_nbins, in_bins, in_nbins, in_bins);
   hist_leading_eta_star_inside_gap =  new TH2D("matched_leading_eta_star_inside_gap","#eta*;#eta_{gen};#eta_{det};#frac{d#sigma}{d#eta*}} [pb]", etastar_nbins, etastar_bins, etastar_nbins, etastar_bins);
   hist_leading_pt_outside_gap =  new TH2D("matched_leading_pt_outside_gap","p_{T};p_{T}^{gen} [GeV];p_{T}^{det} [GeV];#frac{d#sigma}{dp_{T}} [pb]", out_nbins, out_bins, out_nbins, out_bins);
   hist_delta_eta_outside_gap =  new TH2D("matched_delta_eta_outside_gap","#Delta#eta^{out};#Delta#eta^{out}_{gen};#Delta#eta^{out}_{det};#frac{d#sigma}{d#Delta#eta^{out}} [pb]", deta_out_nbins, deta_out_bins, deta_out_nbins, deta_out_bins);

   hist_delta_phi_nomatch =  new TH2D("matched_delta_phi_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta1_nomatch =  new TH2D("matched_delta_phi_deta1_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta2_nomatch =  new TH2D("matched_delta_phi_deta2_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta3_nomatch =  new TH2D("matched_delta_phi_deta3_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta4_nomatch =  new TH2D("matched_delta_phi_deta4_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_gap_nomatch =  new TH2D("matched_delta_phi_gap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta1_gap_nomatch =  new TH2D("matched_delta_phi_deta1_gap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta2_gap_nomatch =  new TH2D("matched_delta_phi_deta2_gap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta3_gap_nomatch =  new TH2D("matched_delta_phi_deta3_gap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta4_gap_nomatch =  new TH2D("matched_delta_phi_deta4_gap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_nogap_nomatch =  new TH2D("matched_delta_phi_nogap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta1_nogap_nomatch =  new TH2D("matched_delta_phi_deta1_nogap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta2_nogap_nomatch =  new TH2D("matched_delta_phi_deta2_nogap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta3_nogap_nomatch =  new TH2D("matched_delta_phi_deta3_nogap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_delta_phi_deta4_nogap_nomatch =  new TH2D("matched_delta_phi_deta4_nogap_nomatch","#Delta#phi;#Delta#phi_{gen} [rad];#Delta#phi_{det} [rad];#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins_nomatch, dphi_bins_nomatch, dphi_nbins_nomatch, dphi_bins_nomatch);
   hist_leading_pt_inside_gap_nomatch =  new TH2D("matched_leading_pt_inside_gap_nomatch","p_{T};p_{T}^{gen} [GeV];p_{T}^{det} [GeV];#frac{d#sigma}{dp_{T}} [pb]", in_nbins_nomatch, in_bins_nomatch, in_nbins_nomatch, in_bins_nomatch);
   hist_leading_eta_star_inside_gap_nomatch =  new TH2D("matched_leading_eta_star_inside_gap_nomatch","#eta*;#eta_{gen};#eta_{det};#frac{d#sigma}{d#eta*}} [pb]", etastar_nbins_nomatch, etastar_bins_nomatch, etastar_nbins_nomatch, etastar_bins_nomatch);
   hist_leading_pt_outside_gap_nomatch =  new TH2D("matched_leading_pt_outside_gap_nomatch","p_{T};p_{T}^{gen} [GeV];p_{T}^{det} [GeV];#frac{d#sigma}{dp_{T}} [pb]", out_nbins_nomatch, out_bins_nomatch, out_nbins_nomatch, out_bins_nomatch);
   hist_delta_eta_outside_gap_nomatch =  new TH2D("matched_delta_eta_outside_gap_nomatch","#Delta#eta^{out};#Delta#eta^{out}_{gen};#Delta#eta^{out}_{det};#frac{d#sigma}{d#Delta#eta^{out}} [pb]", deta_out_nbins_nomatch, deta_out_bins_nomatch, deta_out_nbins_nomatch, deta_out_bins_nomatch);

   hist_delta_phi->Sumw2();
   hist_delta_phi_deta1->Sumw2();
   hist_delta_phi_deta2->Sumw2();
   hist_delta_phi_deta3->Sumw2();
   hist_delta_phi_deta4->Sumw2();
   hist_delta_phi_gap->Sumw2();
   hist_delta_phi_deta1_gap->Sumw2();
   hist_delta_phi_deta2_gap->Sumw2();
   hist_delta_phi_deta3_gap->Sumw2();
   hist_delta_phi_deta4_gap->Sumw2();
   hist_delta_phi_nogap->Sumw2();
   hist_delta_phi_deta1_nogap->Sumw2();
   hist_delta_phi_deta2_nogap->Sumw2();
   hist_delta_phi_deta3_nogap->Sumw2();
   hist_delta_phi_deta4_nogap->Sumw2();
   hist_leading_pt_inside_gap->Sumw2();
   hist_leading_eta_star_inside_gap->Sumw2();
   hist_leading_pt_outside_gap->Sumw2();
   hist_delta_eta_outside_gap->Sumw2();

   hist_delta_phi_nomatch->Sumw2();
   hist_delta_phi_deta1_nomatch->Sumw2();
   hist_delta_phi_deta2_nomatch->Sumw2();
   hist_delta_phi_deta3_nomatch->Sumw2();
   hist_delta_phi_deta4_nomatch->Sumw2();
   hist_delta_phi_gap_nomatch->Sumw2();
   hist_delta_phi_deta1_gap_nomatch->Sumw2();
   hist_delta_phi_deta2_gap_nomatch->Sumw2();
   hist_delta_phi_deta3_gap_nomatch->Sumw2();
   hist_delta_phi_deta4_gap_nomatch->Sumw2();
   hist_delta_phi_nogap_nomatch->Sumw2();
   hist_delta_phi_deta1_nogap_nomatch->Sumw2();
   hist_delta_phi_deta2_nogap_nomatch->Sumw2();
   hist_delta_phi_deta3_nogap_nomatch->Sumw2();
   hist_delta_phi_deta4_nogap_nomatch->Sumw2();
   hist_leading_pt_inside_gap_nomatch->Sumw2();
   hist_leading_eta_star_inside_gap_nomatch->Sumw2();
   hist_leading_pt_outside_gap_nomatch->Sumw2();
   hist_delta_eta_outside_gap_nomatch->Sumw2();

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

  	prescale = 0.0;
  	nentries = tr->GetEntries();
  	cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  	prescale = data_lumi[z];
  	cout<<"Cross-section per event : "<<prescale<<endl;

  	int decade = 0;

  	if (test) { nentries = 100; }

   	for (Int_t i=0;i<nentries;i++)
		{
     		counter_entries++; 
     		double progress = 10.0*i/(1.0*nentries);
     		int k = TMath::FloorNint(progress); 
     		if (k > decade) 
      		cout<<k*10<<" % "<<endl;
     		decade = k; 
     
     		tr->GetEntry(i);

        	pv_pass_det = false;
        	pv_pass_gen = false;
		if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)	{ pv_pass_gen = true; }
		if (sel_mode == "allvertex") 					{ pv_pass_gen = true; }
		if (sel_mode == "nopileup" && Event->evtHdr().pu() == 0)				{ pv_pass_det = true; }
		if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) 			{ pv_pass_det = true; }
		if (sel_mode == "1vertex" && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1) { pv_pass_det = true; }

		vertex_factor = 1.0;
		if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1); }
		
		pt_central_det = 0.0;
		pt_forward_det = 0.0;
		pt_inside_det = 0.0;
		pt_outside_det = 0.0;
		eta_central_det = -9.0;
		eta_forward_det = 0.0;
		eta_inside_det = -9.0;
		eta_star_inside_det = -9.0;
		eta_outside_det = -9.0;
		phi_central_det = -9.0;
		phi_forward_det = -9.0;
		phi_inside_det = -9.0;
		phi_outside_det = -9.0;
		delta_phi_det = -0.5;
		delta_eta_det = -0.5;
		deta_out_det = -0.5;
		pass_det = false;
		pass_deta1_det = false;
		pass_deta2_det = false;
		pass_deta3_det = false;
		pass_deta4_det = false;
		pass_nogap_det = false;
		pass_gap_det = false;
		pass_out_det = false;

		pt_central_gen = 0.0;
		pt_forward_gen = 0.0;
		pt_inside_gen = 0.0;
		pt_outside_gen = 0.0;
		eta_central_gen = -9.0;
		eta_forward_gen = 0.0;
		eta_inside_gen = -9.0;
		eta_star_inside_gen = -9.0;
		eta_outside_gen = -9.0;
		phi_central_gen = -9.0;
		phi_forward_gen = -9.0;
		phi_inside_gen = -9.0;
		phi_outside_gen = -9.0;
		delta_phi_gen = -0.5;
		delta_eta_gen = -0.5;
		deta_out_gen = -0.5;
		pass_gen = false;
		pass_deta1_gen = false;
		pass_deta2_gen = false;
		pass_deta3_gen = false;
		pass_deta4_gen = false;
		pass_nogap_gen = false;
		pass_gap_gen = false;
		pass_out_gen = false;

		deta_forward = 5.0;
		deta_central = 5.0;
		deta_inside = 5.0;
		deta_outside = 5.0;
		dphi_forward = 5.0;
		dphi_central = 5.0;
		dphi_inside = 5.0;
		dphi_outside = 5.0;
		forward_radius = 5.0;
		central_radius = 5.0;
		inside_radius = 5.0;
		outside_radius = 5.0;
		pass_match = false;
		pass_inside_match = false;
		pass_outside_match = false;

		if (pv_pass_gen)
			{
			counter_pv_gen++;

    			for(unsigned int j=0; j<Event->nGenJets(); j++)
				{
				pt = Event->genjet(j).pt();
				eta = Event->genjet(j).eta();
				phi = Event->genjet(j).phi();
     				if (eta <= 2.8 && eta >= -2.8 && pt > pt_central_gen)
     				{ pt_central_gen = pt; eta_central_gen = eta; phi_central_gen = phi; }
     				if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward_gen)
     				{ pt_forward_gen = pt; eta_forward_gen = eta; phi_forward_gen = phi; }
				}
			if (pt_central_gen > pt_min && pt_forward_gen > pt_min)
				{
				pass_gen = true;
				counter_selected_gen++;
				delta_eta_gen = eta_forward_gen - eta_central_gen;
     				if (delta_eta_gen < 0) { delta_eta_gen = -delta_eta_gen; }
				if (delta_eta_gen > 0.4 and delta_eta_gen < 2.5) { pass_deta1_gen = true; }
				if (delta_eta_gen > 2.5 and delta_eta_gen < 3.5) { pass_deta2_gen = true; }
				if (delta_eta_gen > 3.5 and delta_eta_gen < 4.5) { pass_deta3_gen = true; }
				if (delta_eta_gen > 4.5 and delta_eta_gen < 7.5) { pass_deta4_gen = true; }
     				delta_phi_gen = calc_delta_phi(phi_forward_gen, phi_central_gen);
    				for(unsigned int j=0; j<Event->nGenJets(); j++)
					{
					pt = Event->genjet(j).pt();
					eta = Event->genjet(j).eta();
					phi = Event->genjet(j).phi();
     					if (eta_central_gen > eta_forward_gen && eta > eta_forward_gen && eta < eta_central_gen )
						{
        					if (pt > pt_inside_gen)
							{pt_inside_gen = pt; eta_inside_gen = eta; phi_inside_gen = phi; }
        					}
     					if (eta_central_gen < eta_forward_gen && eta < eta_forward_gen && eta > eta_central_gen )
						{
        					if (pt > pt_inside_gen)
							{pt_inside_gen = pt; eta_inside_gen = eta; phi_inside_gen = phi; }
        					}
     					if (eta_central_gen > eta_forward_gen && (eta < eta_forward_gen || eta > eta_central_gen) )
						{
        					if (pt > pt_outside_gen)
							{pt_outside_gen = pt; eta_outside_gen = eta; phi_outside_gen = phi; }
        					}
     					if (eta_central_gen < eta_forward_gen && (eta > eta_forward_gen || eta < eta_central_gen) )
						{
        					if (pt > pt_outside_gen)
							{pt_outside_gen = pt; eta_outside_gen = eta; phi_outside_gen = phi; }
						}
					}
				if (pt_inside_gen < gap_req) { pass_gap_gen = true; }
				if (pt_inside_gen > gap_req)
					{
					pass_nogap_gen = true;
					eta_star_inside_gen = eta_inside_gen - (eta_forward_gen + eta_central_gen)/2;
					}
				if (pt_outside_gen > gap_req)
					{
					pass_out_gen = true;
					deta_out1 = eta_central_gen - eta_outside_gen;
     					if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     					deta_out2 = eta_forward_gen - eta_outside_gen;
     					if (deta_out2 < 0) { deta_out2 = -deta_out2; }
					if (deta_out1 < deta_out2) { deta_out_gen = deta_out1; }
					if (deta_out1 > deta_out2) { deta_out_gen = deta_out2; }
					}
				}
			}


		if (pv_pass_det)
			{
			counter_pv_det++;
    			for(unsigned int j=0; j<Event->nPFJets(); j++)
				{
        			pass_tight = Event->pfjet(j).tightID();
				pt = Event->pfjet(j).ptCor();
				eta = Event->pfjet(j).eta();
				phi = Event->pfjet(j).phi();
				gen_pt = Event->genjet(j).pt();
				gen_eta = Event->genjet(j).eta();
				gen_phi = Event->genjet(j).phi();
				pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,false);
     				if (pass_tight)
					{
     					if (eta <= 2.8 && eta >= -2.8 && pt > pt_central_det)
     					{ pt_central_det = pt; eta_central_det = eta; phi_central_det = phi; }
     					if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward_det)
     					{ pt_forward_det = pt; eta_forward_det = eta; phi_forward_det = phi; }
					}
				}
			if (pt_central_det > pt_min && pt_forward_det > pt_min)
				{
				pass_det = true;
				counter_selected_det++;
				delta_eta_det = eta_forward_det - eta_central_det;
     				if (delta_eta_det < 0) { delta_eta_det = -delta_eta_det; }
				if (delta_eta_det > 0.4 and delta_eta_det < 2.5) { pass_deta1_det = true; }
				if (delta_eta_det > 2.5 and delta_eta_det < 3.5) { pass_deta2_det = true; }
				if (delta_eta_det > 3.5 and delta_eta_det < 4.5) { pass_deta3_det = true; }
				if (delta_eta_det > 4.5 and delta_eta_det < 7.5) { pass_deta4_det = true; }
     				delta_phi_det = calc_delta_phi(phi_forward_det, phi_central_det);
    				for(unsigned int j=0; j<Event->nPFJets(); j++)
					{
        				pass_tight = Event->pfjet(j).tightID();
					pt = Event->pfjet(j).ptCor();
					eta = Event->pfjet(j).eta();
					phi = Event->pfjet(j).phi();
					gen_pt = Event->genjet(j).pt();
					gen_eta = Event->genjet(j).eta();
					gen_phi = Event->genjet(j).phi();
					pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,false);
     					if (eta_central_det > eta_forward_det && eta > eta_forward_det && eta < eta_central_det )
						{
        					if (pt > pt_inside_det)
							{pt_inside_det = pt; eta_inside_det = eta; phi_inside_det = phi; }
        					}
     					if (eta_central_det < eta_forward_det && eta < eta_forward_det && eta > eta_central_det )
						{
        					if (pt > pt_inside_det)
							{pt_inside_det = pt; eta_inside_det = eta; phi_inside_det = phi; }
        					}
     					if (eta_central_det > eta_forward_det && (eta < eta_forward_det || eta > eta_central_det) )
						{
        					if (pt > pt_outside_det)
							{pt_outside_det = pt; eta_outside_det = eta; phi_outside_det = phi; }
        					}
     					if (eta_central_det < eta_forward_det && (eta > eta_forward_det || eta < eta_central_det) )
						{
        					if (pt > pt_outside_det)
							{pt_outside_det = pt; eta_outside_det = eta; phi_outside_det = phi; }
						}
					}
				if (pt_inside_det < gap_req) { pass_gap_det = true; }
				if (pt_inside_det > gap_req)
					{
					pass_nogap_det = true;
					eta_star_inside_det = eta_inside_det - (eta_forward_det + eta_central_det)/2;
					}
				if (pt_outside_det > gap_req)
					{
					pass_out_det = true;
     					deta_out1 = eta_central_det - eta_outside_det;
     					if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     					deta_out2 = eta_forward_det - eta_outside_det;
     					if (deta_out2 < 0) { deta_out2 = -deta_out2; }
					if (deta_out1 < deta_out2) { deta_out_det = deta_out1; }
					if (deta_out1 > deta_out2) { deta_out_det = deta_out2; }
					}
				}
			}

		if (pass_gen && pass_det)
			{
			counter_selected++;
			deta_forward = eta_forward_gen - eta_forward_det;
			if (deta_forward < 0.0) { deta_forward = -deta_forward; }
			deta_central = eta_central_gen - eta_central_det;
			if (deta_central < 0.0) { deta_central = -deta_central; }
			dphi_forward = phi_forward_gen - phi_forward_det;
			if (dphi_forward < 0.0) { dphi_forward = -dphi_forward; }
			dphi_central = phi_central_gen - phi_central_det;
			if (dphi_central < 0.0) { dphi_central = -dphi_central; }
			forward_radius = deta_forward;
			central_radius = deta_central;
			if (central_radius < match_radius && forward_radius < match_radius) { pass_match = true; }
			}
		if (pass_nogap_gen && pass_nogap_det)
			{
			deta_inside = eta_inside_gen - eta_inside_det;
			if (deta_inside < 0.0) { deta_inside = -deta_inside; }
			dphi_inside = phi_inside_gen - phi_inside_det;
			if (dphi_inside < 0.0) { dphi_inside = -dphi_inside; }
			inside_radius = deta_inside + dphi_inside;
			if (inside_radius < match_radius) { pass_inside_match = true; }
			}
		if (pass_out_gen && pass_out_det)
			{
			deta_outside = eta_outside_gen - eta_outside_det;
			if (deta_outside < 0.0) { deta_outside = -deta_outside; }
			dphi_outside = phi_outside_gen - phi_outside_det;
			if (dphi_outside < 0.0) { dphi_outside = -dphi_outside; }
			outside_radius = deta_outside + dphi_outside;
			if (outside_radius < match_radius) { pass_outside_match = true; }
			}

		if (pass_match)
			{
			counter_matched++;
			hist_delta_phi->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor);

			if (pass_deta1_gen && pass_deta1_det) { hist_delta_phi_deta1->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta2_gen && pass_deta2_det) { hist_delta_phi_deta2->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta3_gen && pass_deta3_det) { hist_delta_phi_deta3->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta4_gen && pass_deta4_det) { hist_delta_phi_deta4->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }

			if (pass_gap_gen && pass_gap_det) { hist_delta_phi_gap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta1_gen && pass_deta1_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta1_gap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta2_gen && pass_deta2_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta2_gap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta3_gen && pass_deta3_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta3_gap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta4_gen && pass_deta4_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta4_gap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }

			if (pass_inside_match)
				{
				hist_delta_phi_nogap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor);
				hist_leading_pt_inside_gap->Fill(pt_inside_gen,pt_inside_det,prescale*vertex_factor);
				hist_leading_eta_star_inside_gap->Fill(eta_star_inside_gen,eta_star_inside_det,prescale*vertex_factor);
				if (pass_deta1_gen && pass_deta1_det)
					{ hist_delta_phi_deta1_nogap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
				if (pass_deta2_gen && pass_deta2_det)
					{ hist_delta_phi_deta2_nogap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
				if (pass_deta3_gen && pass_deta3_det)
					{ hist_delta_phi_deta3_nogap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
				if (pass_deta4_gen && pass_deta4_det)
					{ hist_delta_phi_deta4_nogap->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
				}
			if (pass_outside_match)
				{
				hist_leading_pt_outside_gap->Fill(pt_outside_gen,pt_outside_det,prescale*vertex_factor);
				hist_delta_eta_outside_gap->Fill(deta_out_gen,deta_out_det,prescale*vertex_factor);
				}
			}
		if (pass_gen || pass_det)
			{
			hist_delta_phi_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor);

			if (pass_deta1_gen && pass_deta1_det)
				{ hist_delta_phi_deta1_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta1_gen && !pass_deta1_det) { hist_delta_phi_deta1_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (!pass_deta1_gen && pass_deta1_det) { hist_delta_phi_deta1_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); } 
			if (pass_deta2_gen && pass_deta2_det)
				{ hist_delta_phi_deta2_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta2_gen && !pass_deta2_det) { hist_delta_phi_deta2_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (!pass_deta2_gen && pass_deta2_det) { hist_delta_phi_deta2_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); } 
			if (pass_deta3_gen && pass_deta3_det)
				{ hist_delta_phi_deta3_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta3_gen && !pass_deta3_det) { hist_delta_phi_deta3_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (!pass_deta3_gen && pass_deta3_det) { hist_delta_phi_deta3_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); } 
			if (pass_deta4_gen && pass_deta4_det)
				{ hist_delta_phi_deta4_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_deta4_gen && !pass_deta4_det) { hist_delta_phi_deta4_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (!pass_deta4_gen && pass_deta4_det) { hist_delta_phi_deta4_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_nogap_gen && pass_nogap_det)
				{ hist_delta_phi_nogap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_nogap_gen && !pass_nogap_det) { hist_delta_phi_nogap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (!pass_nogap_gen && pass_det && pass_nogap_det) { hist_delta_phi_nogap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }


			if (pass_nogap_gen && pass_nogap_det)
				{ hist_leading_pt_inside_gap_nomatch->Fill(pt_inside_gen,pt_inside_det,prescale*vertex_factor); }
			if (pass_gen && pass_nogap_gen && !pass_nogap_det) { hist_leading_pt_inside_gap_nomatch->Fill(pt_inside_gen,0.0,prescale*vertex_factor); }
			if (!pass_nogap_gen && pass_det && pass_nogap_det) { hist_leading_pt_inside_gap_nomatch->Fill(0.0,pt_inside_det,prescale*vertex_factor); }


			if (pass_nogap_gen && pass_nogap_det && eta_star_inside_gen > -3.6 && eta_star_inside_det > -3.6)
				{ hist_leading_eta_star_inside_gap_nomatch->Fill(eta_star_inside_gen,eta_star_inside_det,prescale*vertex_factor); }
			if (pass_gen && pass_nogap_gen && !pass_nogap_det && eta_star_inside_gen > -3.6)
				{ hist_leading_eta_star_inside_gap_nomatch->Fill(eta_star_inside_gen,-9.0,prescale*vertex_factor); }
			if (!pass_nogap_gen && pass_det && pass_nogap_det && eta_star_inside_det > -3.6)
				{ hist_leading_eta_star_inside_gap_nomatch->Fill(-9.0,eta_star_inside_det,prescale*vertex_factor); }


			if (pass_deta1_gen && pass_deta1_det && pass_nogap_gen && pass_nogap_det)
				{ hist_delta_phi_deta1_nogap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta1_gen && pass_nogap_gen && (!pass_nogap_det || !pass_deta1_det))
				{ hist_delta_phi_deta1_nogap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta1_det && pass_nogap_det && (!pass_nogap_gen || !pass_deta1_gen))
				{ hist_delta_phi_deta1_nogap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta2_gen && pass_deta2_det && pass_nogap_gen && pass_nogap_det)
				{ hist_delta_phi_deta2_nogap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta2_gen && pass_nogap_gen && (!pass_nogap_det || !pass_deta2_det))
				{ hist_delta_phi_deta2_nogap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta2_det && pass_nogap_det && (!pass_nogap_gen || !pass_deta2_gen))
				{ hist_delta_phi_deta2_nogap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta3_gen && pass_deta3_det && pass_nogap_gen && pass_nogap_det)
				{ hist_delta_phi_deta3_nogap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta3_gen && pass_nogap_gen && (!pass_nogap_det || !pass_deta3_det))
				{ hist_delta_phi_deta3_nogap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta3_det && pass_nogap_det && (!pass_nogap_gen || !pass_deta3_gen))
				{ hist_delta_phi_deta3_nogap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta4_gen && pass_deta4_det && pass_nogap_gen && pass_nogap_det)
				{ hist_delta_phi_deta4_nogap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta4_gen && pass_nogap_gen && (!pass_nogap_det || !pass_deta4_det))
				{ hist_delta_phi_deta4_nogap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta4_det && pass_nogap_det && (!pass_nogap_gen || !pass_deta4_gen))
				{ hist_delta_phi_deta4_nogap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }


			if (pass_gap_gen && pass_gap_det) { hist_delta_phi_gap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_gap_gen && !pass_gap_det) { hist_delta_phi_gap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && !pass_gap_gen && pass_gap_det) { hist_delta_phi_gap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta1_gen && pass_deta1_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta1_gap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta1_gen && !pass_nogap_gen && (!pass_gap_det || !pass_deta1_det))
				{ hist_delta_phi_deta1_gap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta1_det && !pass_nogap_det && (!pass_gap_gen || !pass_deta1_gen))
				{ hist_delta_phi_deta1_gap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta2_gen && pass_deta2_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta2_gap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta2_gen && !pass_nogap_gen && (!pass_gap_det || !pass_deta2_det))
				{ hist_delta_phi_deta2_gap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta2_det && !pass_nogap_det && (!pass_gap_gen || !pass_deta2_gen))
				{ hist_delta_phi_deta2_gap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta3_gen && pass_deta3_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta3_gap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta3_gen && !pass_nogap_gen && (!pass_gap_det || !pass_deta3_det))
				{ hist_delta_phi_deta3_gap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta3_det && !pass_nogap_det && (!pass_gap_gen || !pass_deta3_gen))
				{ hist_delta_phi_deta3_gap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }

			if (pass_deta4_gen && pass_deta4_det && pass_gap_gen && pass_gap_det)
				{ hist_delta_phi_deta4_gap_nomatch->Fill(delta_phi_gen,delta_phi_det,prescale*vertex_factor); }
			if (pass_gen && pass_deta4_gen && !pass_nogap_gen && (!pass_gap_det || !pass_deta4_det))
				{ hist_delta_phi_deta4_gap_nomatch->Fill(delta_phi_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && pass_deta4_det && !pass_nogap_det && (!pass_gap_gen || !pass_deta4_gen))
				{ hist_delta_phi_deta4_gap_nomatch->Fill(-0.5,delta_phi_det,prescale*vertex_factor); }


			if (pass_out_gen && pass_out_det)
				{ hist_leading_pt_outside_gap_nomatch->Fill(pt_outside_gen,pt_outside_det,prescale*vertex_factor); }
			if (pass_gen && pass_out_gen && !pass_out_det) { hist_leading_pt_outside_gap_nomatch->Fill(pt_outside_gen,0.0,prescale*vertex_factor); }
			if (pass_det && !pass_out_gen && pass_out_det) { hist_leading_pt_outside_gap_nomatch->Fill(0.0,pt_outside_det,prescale*vertex_factor); }


			if (pass_out_gen && pass_out_det)
				{ hist_delta_eta_outside_gap_nomatch->Fill(deta_out_gen,deta_out_det,prescale*vertex_factor); }
			if (pass_gen && pass_out_gen && !pass_out_det) { hist_delta_eta_outside_gap_nomatch->Fill(deta_out_gen,-0.5,prescale*vertex_factor); }
			if (pass_det && !pass_out_gen && pass_out_det) { hist_delta_eta_outside_gap_nomatch->Fill(-0.5,deta_out_det,prescale*vertex_factor); }

			}
		}

	}

cout << "End of Job Report" << endl;
cout << "Entries read :  " << counter_entries << endl;
cout << "PV pass gen :   " << counter_pv_gen << endl;
cout << "Selected gen :  " << counter_selected_gen << endl;
cout << "PV pass det :   " << counter_pv_det << endl;
cout << "Selected det :  " << counter_selected_det << endl;
cout << "Selected both : " << counter_selected << endl;
cout << "Matched :       " << counter_matched << endl;

     //Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     //save histograms on the file
     hist_delta_phi->Write();
     hist_delta_phi_deta1->Write();
     hist_delta_phi_deta2->Write();
     hist_delta_phi_deta3->Write();
     hist_delta_phi_deta4->Write();
     hist_delta_phi_gap->Write();
     hist_delta_phi_deta1_gap->Write();
     hist_delta_phi_deta2_gap->Write();
     hist_delta_phi_deta3_gap->Write();
     hist_delta_phi_deta4_gap->Write();
     hist_delta_phi_nogap->Write();
     hist_delta_phi_deta1_nogap->Write();
     hist_delta_phi_deta2_nogap->Write();
     hist_delta_phi_deta3_nogap->Write();
     hist_delta_phi_deta4_nogap->Write();
     hist_leading_pt_inside_gap->Write();
     hist_leading_eta_star_inside_gap->Write();
     hist_leading_pt_outside_gap->Write();
     hist_delta_eta_outside_gap->Write();

     hist_delta_phi_nomatch->Write();
     hist_delta_phi_deta1_nomatch->Write();
     hist_delta_phi_deta2_nomatch->Write();
     hist_delta_phi_deta3_nomatch->Write();
     hist_delta_phi_deta4_nomatch->Write();
     hist_delta_phi_gap_nomatch->Write();
     hist_delta_phi_deta1_gap_nomatch->Write();
     hist_delta_phi_deta2_gap_nomatch->Write();
     hist_delta_phi_deta3_gap_nomatch->Write();
     hist_delta_phi_deta4_gap_nomatch->Write();
     hist_delta_phi_nogap_nomatch->Write();
     hist_delta_phi_deta1_nogap_nomatch->Write();
     hist_delta_phi_deta2_nogap_nomatch->Write();
     hist_delta_phi_deta3_nogap_nomatch->Write();
     hist_delta_phi_deta4_nogap_nomatch->Write();
     hist_leading_pt_inside_gap_nomatch->Write();
     hist_leading_eta_star_inside_gap_nomatch->Write();
     hist_leading_pt_outside_gap_nomatch->Write();
     hist_delta_eta_outside_gap_nomatch->Write();

     //close the output file
     data_output->Close();

     delete(hist_delta_phi);
     delete(hist_delta_phi_deta1);
     delete(hist_delta_phi_deta2);
     delete(hist_delta_phi_deta3);
     delete(hist_delta_phi_deta4);
     delete(hist_delta_phi_gap);
     delete(hist_delta_phi_deta1_gap);
     delete(hist_delta_phi_deta2_gap);
     delete(hist_delta_phi_deta3_gap);
     delete(hist_delta_phi_deta4_gap);
     delete(hist_delta_phi_nogap);
     delete(hist_delta_phi_deta1_nogap);
     delete(hist_delta_phi_deta2_nogap);
     delete(hist_delta_phi_deta3_nogap);
     delete(hist_delta_phi_deta4_nogap);
     delete(hist_leading_pt_inside_gap);
     delete(hist_leading_eta_star_inside_gap);
     delete(hist_leading_pt_outside_gap);
     delete(hist_delta_eta_outside_gap);

     delete(hist_delta_phi_nomatch);
     delete(hist_delta_phi_deta1_nomatch);
     delete(hist_delta_phi_deta2_nomatch);
     delete(hist_delta_phi_deta3_nomatch);
     delete(hist_delta_phi_deta4_nomatch);
     delete(hist_delta_phi_gap_nomatch);
     delete(hist_delta_phi_deta1_gap_nomatch);
     delete(hist_delta_phi_deta2_gap_nomatch);
     delete(hist_delta_phi_deta3_gap_nomatch);
     delete(hist_delta_phi_deta4_gap_nomatch);
     delete(hist_delta_phi_nogap_nomatch);
     delete(hist_delta_phi_deta1_nogap_nomatch);
     delete(hist_delta_phi_deta2_nogap_nomatch);
     delete(hist_delta_phi_deta3_nogap_nomatch);
     delete(hist_delta_phi_deta4_nogap_nomatch);
     delete(hist_leading_pt_inside_gap_nomatch);
     delete(hist_leading_eta_star_inside_gap_nomatch);
     delete(hist_leading_pt_outside_gap_nomatch);
     delete(hist_delta_eta_outside_gap_nomatch);
}
