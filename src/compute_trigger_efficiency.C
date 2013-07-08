// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 30 Nov 2012
//
//compute_trigger_efficiencies(string data_triggered, string data_all, string root_out, string prefix, bool detail)
//computes the trigger efficiency

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "common_methods.h"


void compute_eff(TH1D *eff, TH1D *triggered, TH1D *all, double scale = 1.0, bool detail = false)
{

if (detail) { cout << "Scale = " << scale << endl; }
all->Scale(scale);
eff->Divide(triggered, all, 1.0, 1.0, "B");
if (detail) { cout << "Efficiency = " << triggered->Integral()/all->Integral() << endl; }
}


void compute_trigger_efficiency(string data_triggered = "../output/histograms/trigger_eff/events_JetMETTau2010A_all.root", bool triggered_prescaled = false, string data_all = "../output/histograms/trigger_eff/events_JetMETTau2010A_HLT_Jet15U.root", bool all_prescaled = false, string root_out = "../output/histograms/trigger_eff/efficiency_JetMETTau2010A_HLT_Jet15U.root", string prefix = "JetMETTau_2010A_HLT_Jet15U_", string output_path = "../output/trigger_eff/", bool detail = false)
{
//computes the trigger efficiency

//check the input parameters
   if (detail) { cout<<"Compute Trigger Efficiency Configuration"<<endl; }
   if (detail) { cout<<"Triggered Data :       "<<data_triggered<<endl; }
   if (detail) { cout<<"Triggered Prescaled :  "<<triggered_prescaled<<endl; }
   if (detail) { cout<<"All Data :             "<<data_all<<endl; }
   if (detail) { cout<<"All Prescaled :        "<<all_prescaled<<endl; }
   if (detail) { cout<<"Output File :          "<<root_out<<endl; }
   if (detail) { cout<<"Prefix :               "<<prefix<<endl; }
   if (detail) { cout<<"Detail Level :         "<<detail<<endl; }

//setting the histogram sufixes
   string triggered_sufix = "";
   if (!triggered_prescaled) { triggered_sufix = "_nopre"; }
   string all_sufix = "";
   if (!all_prescaled) { all_sufix = "_nopre"; }
   if (detail) { cout<<"Triggered Sufix :      "<<triggered_sufix<<endl; }
   if (detail) { cout<<"All Sufix :            "<<all_sufix<<endl; }

//opening the input data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *trigered = new TFile( data_triggered.c_str() );
    TFile *all = new TFile( data_all.c_str() );
    if (detail) { cout<<"All files opened sucessfully!"<<endl; }

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

//declaring the histogram strings
   TString trig_str, all_str;

//compute the efficiency

//get the downscale factor for the triggered sample
    //TH1D *events = (TH1D*) trigered->Get("Events");
    //double selected = events->GetBinContent(3);
    //double normalization = events->GetBinContent(4);
    //double up_scale = normalization/selected;
    double up_scale = 1.0;


//single jet leading pt
    trig_str = "ak5PF_sj_leading_pt" + triggered_sufix;
    all_str = "ak5PF_sj_leading_pt" + all_sufix;
    TH1D *sj_leading_pt_triggered = 0;
    TH1D *sj_leading_pt_all = 0;
    TH1D *sj_leading_pt;
    trigered->GetObject(trig_str,sj_leading_pt_triggered);
    all->GetObject(all_str,sj_leading_pt_all);
    if (sj_leading_pt_triggered != 0 and sj_leading_pt_all != 0)
    	{
    	sj_leading_pt =  new TH1D("eff_sj_leading_pt","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];Trigger Efficiency", all_nbins, all_bins);
    	sj_leading_pt->Sumw2();
    
    	compute_eff(sj_leading_pt, sj_leading_pt_triggered, sj_leading_pt_all, up_scale, detail);
    	plot_2histograms(sj_leading_pt_triggered, "Triggered and Selected", sj_leading_pt_all, "Selected", output_path, "control_" + prefix + "sj_leading_pt", "top_right", true, detail);
    	plot_efficiency(sj_leading_pt, prefix, "sj_leading_pt", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading pT histogram not found!"<<endl; return; }


//single jet leading pt fine
    trig_str = "ak5PF_sj_leading_pt_fine" + triggered_sufix;
    all_str = "ak5PF_sj_leading_pt_fine" + all_sufix;
    TH1D *sj_leading_pt_fine_triggered = 0;
    TH1D *sj_leading_pt_fine_all = 0;
    TH1D *sj_leading_pt_fine;
    trigered->GetObject(trig_str,sj_leading_pt_fine_triggered);
    all->GetObject(all_str,sj_leading_pt_fine_all);
    if (sj_leading_pt_fine_triggered != 0 and sj_leading_pt_fine_all != 0)
    	{    
    	sj_leading_pt_fine =  new TH1D("eff_sj_leading_pt_fine","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];Trigger Efficiency", 140, 20, 300);
    	sj_leading_pt_fine->Sumw2();
    
    	compute_eff(sj_leading_pt_fine, sj_leading_pt_fine_triggered, sj_leading_pt_fine_all, up_scale, detail);
    	plot_2histograms(sj_leading_pt_fine_triggered, "Triggered and Selected", sj_leading_pt_fine_all, "Selected", output_path, "control_" + prefix + "sj_leading_pt_fine", "top_right", true, detail);
    	plot_efficiency(sj_leading_pt_fine, prefix, "sj_leading_pt_fine", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading pT fine histogram not found!"<<endl; return; }


//single jet leading central pt
    trig_str = "ak5PF_sj_leading_central_pt" + triggered_sufix;
    all_str = "ak5PF_sj_leading_central_pt" + all_sufix;
    TH1D *sj_leading_central_pt_triggered = 0;
    TH1D *sj_leading_central_pt_all = 0;
    TH1D *sj_leading_central_pt;
    trigered->GetObject(trig_str,sj_leading_central_pt_triggered);
    all->GetObject(all_str,sj_leading_central_pt_all);
    if (sj_leading_central_pt_triggered != 0 and sj_leading_central_pt_all != 0)
    	{
    	sj_leading_central_pt =  new TH1D("eff_sj_leading_central_pt","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];Trigger Efficiency", cent_nbins, cent_bins);
    	sj_leading_central_pt->Sumw2();
    
    	compute_eff(sj_leading_central_pt, sj_leading_central_pt_triggered, sj_leading_central_pt_all, up_scale, detail);
    	plot_2histograms(sj_leading_central_pt_triggered, "Triggered and Selected", sj_leading_central_pt_all, "Selected", output_path, "control_" + prefix + "sj_leading_central_pt", "top_right", true, detail);
    	plot_efficiency(sj_leading_central_pt, prefix, "sj_leading_central_pt", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading central pT histogram not found!"<<endl; return; }


//single jet leading forward pt
    trig_str = "ak5PF_sj_leading_forward_pt" + triggered_sufix;
    all_str = "ak5PF_sj_leading_forward_pt" + all_sufix;
    TH1D *sj_leading_forward_pt_triggered = 0;
    TH1D *sj_leading_forward_pt_all = 0;
    TH1D *sj_leading_forward_pt;
    trigered->GetObject(trig_str,sj_leading_forward_pt_triggered);
    all->GetObject(all_str,sj_leading_forward_pt_all);
    if (sj_leading_forward_pt_triggered != 0 and sj_leading_forward_pt_all != 0)
    	{
    sj_leading_forward_pt =  new TH1D("eff_sj_leading_forward_pt","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];Trigger Efficiency", forw_nbins, forw_bins);
    sj_leading_forward_pt->Sumw2();
    
    compute_eff(sj_leading_forward_pt, sj_leading_forward_pt_triggered, sj_leading_forward_pt_all, up_scale, detail);
    plot_2histograms(sj_leading_forward_pt_triggered, "Triggered and Selected", sj_leading_forward_pt_all, "Selected", output_path, "control_" + prefix + "sj_leading_forward_pt", "top_right", true, detail);
    plot_efficiency(sj_leading_forward_pt, prefix, "sj_leading_forward_pt", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading forward pT histogram not found!"<<endl; return; }


//single jet leading eta
    trig_str = "ak5PF_sj_leading_eta" + triggered_sufix;
    all_str = "ak5PF_sj_leading_eta" + all_sufix;
    TH1D *sj_leading_eta_triggered = 0;
    TH1D *sj_leading_eta_all = 0;
    TH1D *sj_leading_eta;
    trigered->GetObject(trig_str,sj_leading_eta_triggered);
    all->GetObject(all_str,sj_leading_eta_all);
    if (sj_leading_eta_triggered != 0 and sj_leading_eta_all != 0)
    	{
    	sj_leading_eta =  new TH1D("eff_sj:leading_eta","Leading Jet #eta;#eta;Trigger Efficiency", eta_nbins, eta_bins);
    	sj_leading_eta->Sumw2();
    
    	compute_eff(sj_leading_eta, sj_leading_eta_triggered, sj_leading_eta_all, up_scale, detail);
    	plot_2histograms(sj_leading_eta_triggered, "Triggered and Selected", sj_leading_eta_all, "Selected", output_path, "control_" + prefix + "sj_leading_eta", "bottom_middle", true, detail);
    	plot_efficiency(sj_leading_eta, prefix, "sj_leading_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading eta histogram not found!"<<endl; return; }


//single jet leading central eta
    trig_str = "ak5PF_sj_leading_central_eta" + triggered_sufix;
    all_str = "ak5PF_sj_leading_central_eta" + all_sufix;
    TH1D *sj_leading_central_eta_triggered = 0;
    TH1D *sj_leading_central_eta_all = 0;
    TH1D *sj_leading_central_eta;
    trigered->GetObject(trig_str,sj_leading_central_eta_triggered);
    all->GetObject(all_str,sj_leading_central_eta_all);
    if (sj_leading_central_eta_triggered != 0 and sj_leading_central_eta_all != 0)
    	{
	sj_leading_central_eta =  new TH1D("eff_sj_leading_central_eta","Leading Central Jet #eta;#eta_{central};Trigger Efficiency", eta_nbins, eta_bins);
	sj_leading_central_eta->Sumw2();
    
	compute_eff(sj_leading_central_eta, sj_leading_central_eta_triggered, sj_leading_central_eta_all, up_scale, detail);
	plot_2histograms(sj_leading_central_eta_triggered, "Triggered and Selected", sj_leading_central_eta_all, "Selected", output_path, "control_" + prefix + "sj_leading_central_eta", "bottom_middle", true, detail);
	plot_efficiency(sj_leading_central_eta, prefix, "sj_leading_central_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading central eta histogram not found!"<<endl; return; }


//single jet leading forward eta
    trig_str = "ak5PF_sj_leading_forward_eta" + triggered_sufix;
    all_str = "ak5PF_sj_leading_forward_eta" + all_sufix;
    TH1D *sj_leading_forward_eta_triggered = 0;
    TH1D *sj_leading_forward_eta_all = 0;
    TH1D *sj_leading_forward_eta;
    trigered->GetObject(trig_str,sj_leading_forward_eta_triggered);
    all->GetObject(all_str,sj_leading_forward_eta_all);
    if (sj_leading_forward_eta_triggered != 0 and sj_leading_forward_eta_all != 0)
    	{
	sj_leading_forward_eta =  new TH1D("eff_sj_leading_forward_eta","Leading Forward Jet #eta;#eta_{forward};Trigger Efficiency", eta_nbins, eta_bins);
	sj_leading_forward_eta->Sumw2();

	compute_eff(sj_leading_forward_eta, sj_leading_forward_eta_triggered, sj_leading_forward_eta_all, up_scale, detail);
	plot_2histograms(sj_leading_forward_eta_triggered, "Triggered and Selected", sj_leading_forward_eta_all, "Selected", output_path, "control_" + prefix + "sj_leading_forward_eta", "bottom_middle", true, detail);
	plot_efficiency(sj_leading_forward_eta, prefix, "sj_leading_forward_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading forward eta histogram not found!"<<endl; return; }


//single jet leading phi
    trig_str = "ak5PF_sj_leading_phi" + triggered_sufix;
    all_str = "ak5PF_sj_leading_phi" + all_sufix;
    TH1D *sj_leading_phi_triggered = 0;
    TH1D *sj_leading_phi_all = 0;
    TH1D *sj_leading_phi;
    trigered->GetObject(trig_str,sj_leading_phi_triggered);
    all->GetObject(all_str,sj_leading_phi_all);
    if (sj_leading_phi_triggered != 0 and sj_leading_phi_all != 0)
    	{
    	sj_leading_phi =  new TH1D("eff_sj_leading_phi","Leading Jet #phi;#phi;Trigger Efficiency", 14, -3.15, 3.15);
    	sj_leading_phi->Sumw2();
    
    	compute_eff(sj_leading_phi, sj_leading_phi_triggered, sj_leading_phi_all, up_scale, detail);
    	plot_2histograms(sj_leading_phi_triggered, "Triggered and Selected", sj_leading_phi_all, "Selected", output_path, "control_" + prefix + "sj_leading_phi", "top_right", true, detail);
    	plot_efficiency(sj_leading_phi, prefix, "sj_leading_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading phi histogram not found!"<<endl; return; }


//single jet leading central phi
    trig_str = "ak5PF_sj_leading_central_phi" + triggered_sufix;
    all_str = "ak5PF_sj_leading_central_phi" + all_sufix;
    TH1D *sj_leading_central_phi_triggered = 0;
    TH1D *sj_leading_central_phi_all = 0;
    TH1D *sj_leading_central_phi;
    trigered->GetObject(trig_str,sj_leading_central_phi_triggered);
    all->GetObject(all_str,sj_leading_central_phi_all);
    if (sj_leading_central_phi_triggered != 0 and sj_leading_central_phi_all != 0)
    	{
    	sj_leading_central_phi =  new TH1D("eff_sj_leading_central_phi","Leading Central Jet #phi;#phi_{central};Trigger Efficiency", 14, -3.15, 3.15);
    	sj_leading_central_phi->Sumw2();
    
    	compute_eff(sj_leading_central_phi, sj_leading_central_phi_triggered, sj_leading_central_phi_all, up_scale, detail);
    	plot_2histograms(sj_leading_central_phi_triggered, "Triggered and Selected", sj_leading_central_phi_all, "Selected", output_path, "control_" + prefix + "sj_leading_central_phi", "top_right", true, detail);
    	plot_efficiency(sj_leading_central_phi, prefix, "sj_leading_central_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading central phi histogram not found!"<<endl; return; }


//single jet leading forward phi
    trig_str = "ak5PF_sj_leading_forward_phi" + triggered_sufix;
    all_str = "ak5PF_sj_leading_forward_phi" + all_sufix;
    TH1D *sj_leading_forward_phi_triggered = 0;
    TH1D *sj_leading_forward_phi_all = 0;
    TH1D *sj_leading_forward_phi;
    trigered->GetObject(trig_str,sj_leading_forward_phi_triggered);
    all->GetObject(all_str,sj_leading_forward_phi_all);
    if (sj_leading_forward_phi_triggered != 0 and sj_leading_forward_phi_all != 0)
    	{
    	sj_leading_forward_phi =  new TH1D("eff_sj_leading_forward_phi","Leading Forward Jet #pho;#phi_{forward};Trigger Efficiency", 14, -3.15, 3.15);
    	sj_leading_forward_phi->Sumw2();
    
    	compute_eff(sj_leading_forward_phi, sj_leading_forward_phi_triggered, sj_leading_forward_phi_all, up_scale, detail);
    	plot_2histograms(sj_leading_forward_phi_triggered, "Triggered and Selected", sj_leading_forward_phi_all, "Selected", output_path, "control_" + prefix + "sj_leading_forward_phi", "top_right", true, detail);
    	plot_efficiency(sj_leading_forward_phi, prefix, "sj_leading_forward_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"SJ Leading forward phi histogram not found!"<<endl; return; }


//leading pt
    trig_str = "ak5PF_leading_pt" + triggered_sufix;
    all_str = "ak5PF_leading_pt" + all_sufix;
    TH1D *leading_pt_triggered = 0;
    TH1D *leading_pt_all = 0;
    TH1D *leading_pt;
    trigered->GetObject(trig_str,leading_pt_triggered);
    all->GetObject(all_str,leading_pt_all);
    if (leading_pt_triggered != 0 and leading_pt_all != 0)
    	{
    	leading_pt =  new TH1D("eff_leading_pt","Leading Jet p_{T};p_{T} [#frac{GeV}{c}];Trigger Efficiency", all_nbins, all_bins);
    	leading_pt->Sumw2();
    
    	compute_eff(leading_pt, leading_pt_triggered, leading_pt_all, up_scale, detail);
    	plot_2histograms(leading_pt_triggered, "Triggered and Selected", leading_pt_all, "Selected", output_path, "control_" + prefix + "leading_pt", "top_right", true, detail);
    	plot_efficiency(leading_pt, prefix, "leading_pt", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading pT histogram not found!"<<endl; return; }


//leading pt fine
    trig_str = "ak5PF_leading_pt_fine" + triggered_sufix;
    all_str = "ak5PF_leading_pt_fine" + all_sufix;
    TH1D *leading_pt_fine_triggered = 0;
    TH1D *leading_pt_fine_all = 0;
    TH1D *leading_pt_fine;
    trigered->GetObject(trig_str,leading_pt_fine_triggered);
    all->GetObject(all_str,leading_pt_fine_all);
    if (leading_pt_fine_triggered != 0 and leading_pt_fine_all != 0)
    	{
    	leading_pt_fine =  new TH1D("eff_leading_pt_fine","Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];Trigger Efficiency", 140, 20, 300);
    	leading_pt_fine->Sumw2();
    
    	compute_eff(leading_pt_fine, leading_pt_fine_triggered, leading_pt_fine_all, up_scale, detail);
    	plot_2histograms(leading_pt_fine_triggered, "Triggered and Selected", leading_pt_fine_all, "Selected", output_path, "control_" + prefix + "leading_pt_fine", "top_right", true, detail);
    	plot_efficiency(leading_pt_fine, prefix, "leading_pt_fine", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading pT fine histogram not found!"<<endl; return; }


//leading central pt
    trig_str = "ak5PF_leading_central_pt" + triggered_sufix;
    all_str = "ak5PF_leading_central_pt" + all_sufix;
    TH1D *leading_central_pt_triggered = 0;
    TH1D *leading_central_pt_all = 0;
    TH1D *leading_central_pt;
    trigered->GetObject(trig_str,leading_central_pt_triggered);
    all->GetObject(all_str,leading_central_pt_all);
    if (leading_central_pt_triggered != 0 and leading_central_pt_all != 0)
    	{
    	leading_central_pt =  new TH1D("eff_leading_central_pt","Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];Trigger Efficiency", cent_nbins, cent_bins);
    	leading_central_pt->Sumw2();
    
    	compute_eff(leading_central_pt, leading_central_pt_triggered, leading_central_pt_all, up_scale, detail);
    	plot_2histograms(leading_central_pt_triggered, "Triggered and Selected", leading_central_pt_all, "Selected", output_path, "control_" + prefix + "leading_central_pt", "top_right", true, detail);
    	plot_efficiency(leading_central_pt, prefix, "leading_central_pt", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading central pT histogram not found!"<<endl; return; }


//leading forward pt
    trig_str = "ak5PF_leading_forward_pt" + triggered_sufix;
    all_str = "ak5PF_leading_forward_pt" + all_sufix;
    TH1D *leading_forward_pt_triggered = 0;
    TH1D *leading_forward_pt_all = 0;
    TH1D *leading_forward_pt;
    trigered->GetObject(trig_str,leading_forward_pt_triggered);
    all->GetObject(all_str,leading_forward_pt_all);
    if (leading_forward_pt_triggered != 0 and leading_forward_pt_all != 0)
    	{
    	leading_forward_pt =  new TH1D("eff_leading_forward_pt","Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];Trigger Efficiency", forw_nbins, forw_bins);
    	leading_forward_pt->Sumw2();
    
    	compute_eff(leading_forward_pt, leading_forward_pt_triggered, leading_forward_pt_all, up_scale, detail);
    	plot_2histograms(leading_forward_pt_triggered, "Triggered and Selected", leading_forward_pt_all, "Selected", output_path, "control_" + prefix + "leading_forward_pt", "top_right", true, detail);
    	plot_efficiency(leading_forward_pt, prefix, "leading_forward_pt", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading forward pT histogram not found!"<<endl; return; }


//leading eta
    trig_str = "ak5PF_leading_eta" + triggered_sufix;
    all_str = "ak5PF_leading_eta" + all_sufix;
    TH1D *leading_eta_triggered = 0;
    TH1D *leading_eta_all = 0;
    TH1D *leading_eta;
    trigered->GetObject(trig_str,leading_eta_triggered);
    all->GetObject(all_str,leading_eta_all);
    if (leading_eta_triggered != 0 and leading_eta_all != 0)
    	{
    	leading_eta =  new TH1D("eff_leading_eta","Leading Jet #eta;#eta;Trigger Efficiency", eta_nbins, eta_bins);
    	leading_eta->Sumw2();
    
    	compute_eff(leading_eta, leading_eta_triggered, leading_eta_all, up_scale, detail);
    	plot_2histograms(leading_eta_triggered, "Triggered and Selected", leading_eta_all, "Selected", output_path, "control_" + prefix + "leading_eta", "bottom_middle", true, detail);
    	plot_efficiency(leading_eta, prefix, "leading_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading eta histogram not found!"<<endl; return; }


//leading central eta
    trig_str = "ak5PF_leading_central_eta" + triggered_sufix;
    all_str = "ak5PF_leading_central_eta" + all_sufix;
    TH1D *leading_central_eta_triggered = 0;
    TH1D *leading_central_eta_all = 0;
    TH1D *leading_central_eta;
    trigered->GetObject(trig_str,leading_central_eta_triggered);
    all->GetObject(all_str,leading_central_eta_all);
    if (leading_central_eta_triggered != 0 and leading_central_eta_all != 0)
    	{
    	leading_central_eta =  new TH1D("eff_leading_central_eta","Leading Central Jet #eta;#eta_{central};Trigger Efficiency", eta_nbins, eta_bins);
    	leading_central_eta->Sumw2();
    
    	compute_eff(leading_central_eta, leading_central_eta_triggered, leading_central_eta_all, up_scale, detail);
    	plot_2histograms(leading_central_eta_triggered, "Triggered and Selected", leading_central_eta_all, "Selected", output_path, "control_" + prefix + "leading_central_eta", "bottom_middle", true, detail);
    	plot_efficiency(leading_central_eta, prefix, "leading_central_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading central eta histogram not found!"<<endl; return; }


//leading forward eta
    trig_str = "ak5PF_leading_forward_eta" + triggered_sufix;
    all_str = "ak5PF_leading_forward_eta" + all_sufix;
    TH1D *leading_forward_eta_triggered = 0;
    TH1D *leading_forward_eta_all = 0;
    TH1D *leading_forward_eta;
    trigered->GetObject(trig_str,leading_forward_eta_triggered);
    all->GetObject(all_str,leading_forward_eta_all);
    if (leading_forward_eta_triggered != 0 and leading_forward_eta_all != 0)
    	{
    	leading_forward_eta =  new TH1D("eff_leading_forward_eta","Leading Forward Jet #eta;#eta_{forward};Trigger Efficiency", eta_nbins, eta_bins);
    	leading_forward_eta->Sumw2();
    
    	compute_eff(leading_forward_eta, leading_forward_eta_triggered, leading_forward_eta_all, up_scale, detail);
    	plot_2histograms(leading_forward_eta_triggered, "Triggered and Selected", leading_forward_eta_all, "Selected", output_path, "control_" + prefix + "leading_forward_eta", "bottom_middle", true, detail);
    	plot_efficiency(leading_forward_eta, prefix, "leading_forward_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading forward eta histogram not found!"<<endl; return; }


//leading phi
    trig_str = "ak5PF_leading_phi" + triggered_sufix;
    all_str = "ak5PF_leading_phi" + all_sufix;
    TH1D *leading_phi_triggered = 0;
    TH1D *leading_phi_all = 0;
    TH1D *leading_phi;
    trigered->GetObject(trig_str,leading_phi_triggered);
    all->GetObject(all_str,leading_phi_all);
    if (leading_phi_triggered != 0 and leading_phi_all != 0)
    	{
    	leading_phi =  new TH1D("eff_leading_phi","Leading Jet #phi;#phi;Trigger Efficiency", 14, -3.15, 3.15);
    	leading_phi->Sumw2();
    
    	compute_eff(leading_phi, leading_phi_triggered, leading_phi_all, up_scale, detail);
    	plot_2histograms(leading_phi_triggered, "Triggered and Selected", leading_phi_all, "Selected", output_path, "control_" + prefix + "leading_phi", "top_right", true, detail);
    	plot_efficiency(leading_phi, prefix, "leading_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading phi histogram not found!"<<endl; return; }


//leading central phi
    trig_str = "ak5PF_leading_central_phi" + triggered_sufix;
    all_str = "ak5PF_leading_central_phi" + all_sufix;
    TH1D *leading_central_phi_triggered = 0;
    TH1D *leading_central_phi_all = 0;
    TH1D *leading_central_phi;
    trigered->GetObject(trig_str,leading_central_phi_triggered);
    all->GetObject(all_str,leading_central_phi_all);
    if (leading_central_phi_triggered != 0 and leading_central_phi_all != 0)
    	{
    	leading_central_phi =  new TH1D("eff_leading_central_phi","Leading Central Jet #phi;#phi_{central};Trigger Efficiency", 14, -3.15, 3.15);
    	leading_central_phi->Sumw2();
    
    	compute_eff(leading_central_phi, leading_central_phi_triggered, leading_central_phi_all, up_scale, detail);
    	plot_2histograms(leading_central_phi_triggered, "Triggered and Selected", leading_central_phi_all, "Selected", output_path, "control_" + prefix + "leading_central_phi", "top_right", true, detail);
    	plot_efficiency(leading_central_phi, prefix, "leading_central_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading central phi histogram not found!"<<endl; return; }


//leading forward phi
    trig_str = "ak5PF_leading_forward_phi" + triggered_sufix;
    all_str = "ak5PF_leading_forward_phi" + all_sufix;
    TH1D *leading_forward_phi_triggered = 0;
    TH1D *leading_forward_phi_all = 0;
    TH1D *leading_forward_phi;
    trigered->GetObject(trig_str,leading_forward_phi_triggered);
    all->GetObject(all_str,leading_forward_phi_all);
    if (leading_forward_phi_triggered != 0 and leading_forward_phi_all != 0)
    	{
    	leading_forward_phi =  new TH1D("eff_leading_forward_phi","Leading Forward Jet #pho;#phi_{forward};Trigger Efficiency", 14, -3.15, 3.15);
    	leading_forward_phi->Sumw2();
    
    	compute_eff(leading_forward_phi, leading_forward_phi_triggered, leading_forward_phi_all, up_scale, detail);
    	plot_2histograms(leading_forward_phi_triggered, "Triggered and Selected", leading_forward_phi_all, "Selected", output_path, "control_" + prefix + "leading_forward_phi", "top_right", true, detail);
    	plot_efficiency(leading_forward_phi, prefix, "leading_forward_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading forward phi histogram not found!"<<endl; return; }


//delta eta
    trig_str = "ak5PF_delta_eta" + triggered_sufix;
    all_str = "ak5PF_delta_eta" + all_sufix;
    TH1D *delta_eta_triggered = 0;
    TH1D *delta_eta_all = 0;
    TH1D *delta_eta;
    trigered->GetObject(trig_str,delta_eta_triggered);
    all->GetObject(all_str,delta_eta_all);
    if (delta_eta_triggered != 0 and delta_eta_all != 0)
    	{
    	delta_eta =  new TH1D("eff_delta_eta","#Delta#eta;#Delta#eta;Trigger Efficiency", deta_nbins, deta_bins);
    	delta_eta->Sumw2();
    
    	compute_eff(delta_eta, delta_eta_triggered, delta_eta_all, up_scale, detail);
    	plot_2histograms(delta_eta_triggered, "Triggered and Selected", delta_eta_all, "Selected", output_path, "control_" + prefix + "delta_eta", "top_right", true, detail);
    	plot_efficiency(delta_eta, prefix, "delta_eta", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading delta eta histogram not found!"<<endl; return; }


//delta phi
    trig_str = "ak5PF_delta_phi" + triggered_sufix;
    all_str = "ak5PF_delta_phi" + all_sufix;
    TH1D *delta_phi_triggered = 0;
    TH1D *delta_phi_all = 0;
    TH1D *delta_phi;
    trigered->GetObject(trig_str,delta_phi_triggered);
    all->GetObject(all_str,delta_phi_all);
    if (delta_phi_triggered != 0 and delta_phi_all != 0)
    	{
    	delta_phi =  new TH1D("eff_delta_phi","#Delta#phi;#Delta#phi;Trigger Efficiency", dphi_nbins, dphi_bins);
    	delta_phi->Sumw2();
    
    	compute_eff(delta_phi, delta_phi_triggered, delta_phi_all, up_scale, detail);
    	plot_2histograms(delta_phi_triggered, "Triggered and Selected", delta_phi_all, "Selected", output_path, "control_" + prefix + "delta_phi", "top_left", true, detail);
    	plot_efficiency(delta_phi, prefix, "delta_phi", output_path, "top_left", detail);
	}
    else
	{ cout<<"Leading delta phi histogram not found!"<<endl; return; }


//creating the output file
    TFile data_output( root_out.c_str() , "RECREATE");
    
//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file "<<root_out<<" ..."<<endl; }
    leading_pt->Write();
    leading_pt_fine->Write();
    leading_central_pt->Write();
    leading_forward_pt->Write();
    leading_eta->Write();
    leading_central_eta->Write();
    leading_forward_eta->Write();
    leading_phi->Write();
    leading_central_phi->Write();
    leading_forward_phi->Write();
    delta_eta->Write();
    delta_phi->Write();
    sj_leading_pt->Write();
    sj_leading_pt_fine->Write();
    sj_leading_central_pt->Write();
    sj_leading_forward_pt->Write();
    sj_leading_eta->Write();
    sj_leading_central_eta->Write();
    sj_leading_forward_eta->Write();
    sj_leading_phi->Write();
    sj_leading_central_phi->Write();
    sj_leading_forward_phi->Write();
    if (detail) { cout<<"Writing on "<<root_out<<" was sucessfull!"<<endl; }
   
//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    trigered->Close();
    all->Close();
    data_output.Close();
    if (detail) { cout<<"Done!"<<endl; }
}
