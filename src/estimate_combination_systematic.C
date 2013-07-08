// Pedro Cipriano, Jan 2012
// DESY, CMS
// Last Update: 05 Mar 2013
//
// compute_trigger_turn_on(string *data_in, string data_out, int n_files, string sel_mode = "allvertex", bool detail = false, bool test = false)
// Computes the trigger turn on curves using trigger elements

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
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

void estimate_trigger_systematic(TH1D *sample, TH1D *all, TH1D *out)
{
double sample_int = sample->Integral();
double all_int = all->Integral();
double value = 0, error = 0;

cout << "Integral = " << sample_int << " " << all_int << endl;

out->Divide(sample,all,all_int,sample_int,"");

double out_int = out->Integral();
cout << "Integral = " << out_int << endl;

for (int i = 1; i <= out->GetNbinsX(); i++)
	{
	value = out->GetBinContent(i);
	//error = sample->GetBinError(i) / sample->GetBinContent(i);
	error = out->GetBinError(i);
	if (value > 2.0) { value = 2.0; }
	value = value - 1.0;
	cout << " i " << i << " val = " << value << " error = " << error << endl;
	out->SetBinContent(i,value);
	//out->SetBinError(i,error);
	}

out->SetEntries(sample->GetEntries());

}

void merge_and_plot_3histograms(TH1D *jetmettau, TH1D *jetmet, TH1D *jet, TH1D *merged, string path, string fileout, string legend_position, bool detail = false)
{
//plots the different datasets and the merged result
//declares the variables needed for the calculation
    double entries[3] = {0.,0.,0.};
    double total_entries = 0.0;
    double weigth[3] = {1.0, 1.0, 1.0};

//calculate the scale for the first dataset
    entries[0] = jetmettau->GetEntries();
    if (entries[0] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[0] << endl; } entries[0] = 1.0; }
    if (detail) { cout<<"JetMETTau_2010A "<<setw(8)<<entries[0]<<endl; }

//calculate the scale for the second dataset
    entries[1] = jetmet->GetEntries();
    if (entries[1] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[1] << endl; } entries[1] = 1.0; }
    if (detail) { cout<<"JetMET_2010A    "<<setw(8)<<entries[1]<<endl; }

//calculate the scale for the third dataset
    entries[2] = jet->GetEntries();
    if (entries[2] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[2] << endl; } entries[2] = 1.0; }
    if (detail) { cout<<"Jet_2010B       "<<setw(8)<<entries[2]<<endl; }

//calculate the total number of entries
    total_entries = entries[0] + entries[1] + entries[2];
    // total_entries = entries[0] + entries[1]; //old version when we used only 2 datasets
    if (total_entries <= 0) { if (detail) { cout<<"Bad number of total entries = " << total_entries << endl; } total_entries = 1.0; }
    weigth[0] = entries[0]/total_entries;
    weigth[1] = entries[1]/total_entries;
    weigth[2] = entries[2]/total_entries;

    //merging
    merged->Add(jetmettau,weigth[0]);
    merged->Add(jetmet,weigth[1]);
    merged->Add(jet,weigth[2]);

for (int i = 1; i <= merged->GetNbinsX(); i++)
	{
	if (merged->GetBinContent(i) < 0.0) { merged->SetBinContent(i,-merged->GetBinContent(i)); }
	}

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    //gPad->SetLogy();
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    
    double max = 1.0;
    double min = -1.0;
    
//plooting
    merged->SetMaximum(max);
    merged->SetMinimum(min);
    merged->SetFillColor(5);
    merged->SetLineWidth(4);
    merged->Draw("e2");
    merged->Draw("same");
    jetmettau->SetLineWidth(4);
    jetmettau->SetLineColor(2);
    jetmettau->SetLineStyle(1);
    jetmettau->Draw("e1 same");
    jetmet->SetLineWidth(4);
    jetmet->SetLineColor(3);
    jetmet->SetLineStyle(2);
    jetmet->Draw("e1same");
    jet->SetLineWidth(4);
    jet->SetLineColor(4);
    jet->SetLineStyle(4);
    jet->Draw("e1same");
    
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 4, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(jetmettau,"JetMETTau_2010A","l");
    leg01->AddEntry(jetmet,"JetMET_2010A ","l");
    leg01->AddEntry(jet,"Jet_2010B","l");
    leg01->AddEntry(merged,"Merged Data","f");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, prefix + fileout);
    
    merged->SetMaximum(0.3);
    merged->SetMinimum(0.0);

}

void merge_plot_combination_systematic(string data_jetmettau, string data_jetmet, string data_jet, string data_out, string output_path, bool detail = false, bool test = false)
{

   if (detail) { cout<<"Merge and Plot Combination Systematic Configuration"<<endl; }
   if (detail) { cout<<"Source for JetMETTau Systematic : "<<data_jetmettau<<endl; }
   if (detail) { cout<<"Source for JetMET Systematic :    "<<data_jetmet<<endl; }
   if (detail) { cout<<"Source for Jet Systematic :       "<<data_jet<<endl; }
   if (detail) { cout<<"Output File :                     "<<data_out<<endl; }
   if (detail) { cout<<"Output Path :                     "<<output_path<<endl; }
   if (detail) { cout<<"Detail level :                    "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :                       "<<test<<endl; }

   bool not_loaded = false;

   //binning
   int dphi_nbins = 7;
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

   int in_nbins = 9;
   double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

   int etastar_nbins = 12;
   double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

   int deta_out_nbins = 6;
   double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

   //Opening the source files
   TFile *jetmettau_file = new TFile( data_jetmettau.c_str() );
   TFile *jetmet_file = new TFile( data_jetmet.c_str() );
   TFile *jet_file = new TFile( data_jet.c_str() );

   //declaring empty pointers
   TH1D *jetmettau_delta_phi = 0;
   TH1D *jetmettau_delta_phi_deta1 = 0;
   TH1D *jetmettau_delta_phi_deta2 = 0;
   TH1D *jetmettau_delta_phi_deta3 = 0;
   TH1D *jetmettau_delta_phi_deta4 = 0;
   TH1D *jetmettau_delta_phi_gap = 0;
   TH1D *jetmettau_delta_phi_deta1_gap = 0;
   TH1D *jetmettau_delta_phi_deta2_gap = 0;
   TH1D *jetmettau_delta_phi_deta3_gap = 0;
   TH1D *jetmettau_delta_phi_deta4_gap = 0;
   TH1D *jetmettau_delta_phi_nogap = 0;
   TH1D *jetmettau_delta_phi_deta1_nogap = 0;
   TH1D *jetmettau_delta_phi_deta2_nogap = 0;
   TH1D *jetmettau_delta_phi_deta3_nogap = 0;
   TH1D *jetmettau_delta_phi_deta4_nogap = 0;
   TH1D *jetmettau_leading_pt_inside_gap = 0;
   TH1D *jetmettau_leading_eta_star_inside_gap = 0;
   TH1D *jetmettau_delta_eta_outside_gap = 0;

   TH1D *jetmet_delta_phi = 0;
   TH1D *jetmet_delta_phi_deta1 = 0;
   TH1D *jetmet_delta_phi_deta2 = 0;
   TH1D *jetmet_delta_phi_deta3 = 0;
   TH1D *jetmet_delta_phi_deta4 = 0;
   TH1D *jetmet_delta_phi_gap = 0;
   TH1D *jetmet_delta_phi_deta1_gap = 0;
   TH1D *jetmet_delta_phi_deta2_gap = 0;
   TH1D *jetmet_delta_phi_deta3_gap = 0;
   TH1D *jetmet_delta_phi_deta4_gap = 0;
   TH1D *jetmet_delta_phi_nogap = 0;
   TH1D *jetmet_delta_phi_deta1_nogap = 0;
   TH1D *jetmet_delta_phi_deta2_nogap = 0;
   TH1D *jetmet_delta_phi_deta3_nogap = 0;
   TH1D *jetmet_delta_phi_deta4_nogap = 0;
   TH1D *jetmet_leading_pt_inside_gap = 0;
   TH1D *jetmet_leading_eta_star_inside_gap = 0;
   TH1D *jetmet_delta_eta_outside_gap = 0;

   TH1D *jet_delta_phi = 0;
   TH1D *jet_delta_phi_deta1 = 0;
   TH1D *jet_delta_phi_deta2 = 0;
   TH1D *jet_delta_phi_deta3 = 0;
   TH1D *jet_delta_phi_deta4 = 0;
   TH1D *jet_delta_phi_gap = 0;
   TH1D *jet_delta_phi_deta1_gap = 0;
   TH1D *jet_delta_phi_deta2_gap = 0;
   TH1D *jet_delta_phi_deta3_gap = 0;
   TH1D *jet_delta_phi_deta4_gap = 0;
   TH1D *jet_delta_phi_nogap = 0;
   TH1D *jet_delta_phi_deta1_nogap = 0;
   TH1D *jet_delta_phi_deta2_nogap = 0;
   TH1D *jet_delta_phi_deta3_nogap = 0;
   TH1D *jet_delta_phi_deta4_nogap = 0;
   TH1D *jet_leading_pt_inside_gap = 0;
   TH1D *jet_leading_eta_star_inside_gap = 0;
   TH1D *jet_delta_eta_outside_gap = 0;

   //loading histograms
   jetmettau_file->GetObject("trigger_syst_delta_phi_inclusive",jetmettau_delta_phi);
   if (jetmettau_delta_phi == 0) { cout << "trigger_syst_delta_phi_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta1_inclusive",jetmettau_delta_phi_deta1);
   if (jetmettau_delta_phi_deta1 == 0) { cout << "trigger_syst_delta_phi_deta1_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta2_inclusive",jetmettau_delta_phi_deta2);
   if (jetmettau_delta_phi_deta2 == 0) { cout << "trigger_syst_delta_phi_deta2_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta3_inclusive",jetmettau_delta_phi_deta3);
   if (jetmettau_delta_phi_deta3 == 0) { cout << "trigger_syst_delta_phi_deta3_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta4_inclusive",jetmettau_delta_phi_deta4);
   if (jetmettau_delta_phi_deta4 == 0) { cout << "trigger_syst_delta_phi_deta4_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_gap_inclusive",jetmettau_delta_phi_gap);
   if (jetmettau_delta_phi_gap == 0) { cout << "trigger_syst_delta_phi_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta1_gap_inclusive",jetmettau_delta_phi_deta1_gap);
   if (jetmettau_delta_phi_deta1_gap == 0) { cout << "trigger_syst_delta_phi_gap_deta1_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta2_gap_inclusive",jetmettau_delta_phi_deta2_gap);
   if (jetmettau_delta_phi_deta2_gap == 0) { cout << "trigger_syst_delta_phi_deta2_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta3_gap_inclusive",jetmettau_delta_phi_deta3_gap);
   if (jetmettau_delta_phi_deta3_gap == 0) { cout << "trigger_syst_delta_phi_deta3_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta4_gap_inclusive",jetmettau_delta_phi_deta4_gap);
   if (jetmettau_delta_phi_deta4_gap == 0) { cout << "trigger_syst_delta_phi_deta4_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_nogap_inclusive",jetmettau_delta_phi_nogap);
   if (jetmettau_delta_phi_nogap == 0) { cout << "trigger_syst_delta_phi_nogap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta1_nogap_inclusive",jetmettau_delta_phi_deta1_nogap);
   if (jetmettau_delta_phi_deta1_nogap == 0) { cout << "trigger_syst_delta_phi_deta1_nogap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta2_nogap_inclusive",jetmettau_delta_phi_deta2_nogap);
   if (jetmettau_delta_phi_deta2_nogap == 0) { cout << "trigger_syst_delta_phi_deta2_nogap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta3_nogap_inclusive",jetmettau_delta_phi_deta3_nogap);
   if (jetmettau_delta_phi_deta3_nogap == 0) { cout << "trigger_syst_delta_phi_nogap_deta3_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_phi_deta4_nogap_inclusive",jetmettau_delta_phi_deta4_nogap);
   if (jetmettau_delta_phi_deta4_nogap == 0) { cout << "trigger_syst_delta_phi_deta4_nogap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_leading_pt_inside_gap_inclusive",jetmettau_leading_pt_inside_gap);   
   if (jetmettau_leading_pt_inside_gap == 0) { cout << "trigger_syst_leading_pt_inside_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_leading_eta_star_inside_gap_inclusive",jetmettau_leading_eta_star_inside_gap);   
   if (jetmettau_leading_eta_star_inside_gap == 0) { cout << "trigger_syst_leading_eta_star_inside_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("trigger_syst_delta_eta_outside_gap_inclusive",jetmettau_delta_eta_outside_gap);   
   if (jetmettau_delta_eta_outside_gap == 0) { cout << "trigger_syst_delta_eta_outside_gap_inclusive on JetMETTau_2010A not found!" << endl; not_loaded = true; }


   jetmet_file->GetObject("trigger_syst_delta_phi_inclusive",jetmet_delta_phi);
   if (jetmet_delta_phi == 0) { cout << "trigger_syst_delta_phi_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta1_inclusive",jetmet_delta_phi_deta1);
   if (jetmet_delta_phi_deta1 == 0) { cout << "trigger_syst_delta_phi_deta1_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta2_inclusive",jetmet_delta_phi_deta2);
   if (jetmet_delta_phi_deta2 == 0) { cout << "trigger_syst_delta_phi_deta2_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta3_inclusive",jetmet_delta_phi_deta3);
   if (jetmet_delta_phi_deta3 == 0) { cout << "trigger_syst_delta_phi_deta3_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta4_inclusive",jetmet_delta_phi_deta4);
   if (jetmet_delta_phi_deta4 == 0) { cout << "trigger_syst_delta_phi_deta4_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_gap_inclusive",jetmet_delta_phi_gap);
   if (jetmet_delta_phi_gap == 0) { cout << "trigger_syst_delta_phi_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta1_gap_inclusive",jetmet_delta_phi_deta1_gap);
   if (jetmet_delta_phi_deta1_gap == 0) { cout << "trigger_syst_delta_phi_deta1_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta2_gap_inclusive",jetmet_delta_phi_deta2_gap);
   if (jetmet_delta_phi_deta2_gap == 0) { cout << "trigger_syst_delta_phi_deta2_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta3_gap_inclusive",jetmet_delta_phi_deta3_gap);
   if (jetmet_delta_phi_deta3_gap == 0) { cout << "trigger_syst_delta_phi_deta3_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta4_gap_inclusive",jetmet_delta_phi_deta4_gap);
   if (jetmet_delta_phi_deta4_gap == 0) { cout << "trigger_syst_delta_phi_deta4_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_nogap_inclusive",jetmet_delta_phi_nogap);
   if (jetmet_delta_phi_nogap == 0) { cout << "trigger_syst_delta_phi_nogap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta1_nogap_inclusive",jetmet_delta_phi_deta1_nogap);
   if (jetmet_delta_phi_deta1_nogap == 0) { cout << "trigger_syst_delta_phi_deta1_nogap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta2_nogap_inclusive",jetmet_delta_phi_deta2_nogap);
   if (jetmet_delta_phi_deta2_nogap == 0) { cout << "trigger_syst_delta_phi_deta2_nogap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta3_nogap_inclusive",jetmet_delta_phi_deta3_nogap);
   if (jetmet_delta_phi_deta3_nogap == 0) { cout << "trigger_syst_delta_phi_deta3_nogap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_phi_deta4_nogap_inclusive",jetmet_delta_phi_deta4_nogap);
   if (jetmet_delta_phi_deta4_nogap == 0) { cout << "trigger_syst_delta_phi_deta4_nogap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_leading_pt_inside_gap_inclusive",jetmet_leading_pt_inside_gap);   
   if (jetmet_leading_pt_inside_gap == 0) { cout << "trigger_syst_leading_pt_inside_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_leading_eta_star_inside_gap_inclusive",jetmet_leading_eta_star_inside_gap);   
   if (jetmet_leading_eta_star_inside_gap == 0) { cout << "trigger_syst_leading_eta_star_inside_gap_exclusive on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("trigger_syst_delta_eta_outside_gap_inclusive",jetmet_delta_eta_outside_gap);   
   if (jetmet_delta_eta_outside_gap == 0) { cout << "trigger_syst_delta_eta_outside_gap_inclusive on JetMET_2010A not found!" << endl; not_loaded = true; }


   jet_file->GetObject("trigger_syst_delta_phi_inclusive",jet_delta_phi);
   if (jet_delta_phi == 0) { cout << "trigger_syst_delta_phi_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta1_inclusive",jet_delta_phi_deta1);
   if (jet_delta_phi_deta1 == 0) { cout << "trigger_syst_delta_phi_deta1_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta2_inclusive",jet_delta_phi_deta2);
   if (jet_delta_phi_deta2 == 0) { cout << "trigger_syst_delta_phi_deta2_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta3_inclusive",jet_delta_phi_deta3);
   if (jet_delta_phi_deta3 == 0) { cout << "trigger_syst_delta_phi_deta3_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta4_inclusive",jet_delta_phi_deta4);
   if (jet_delta_phi_deta4 == 0) { cout << "trigger_syst_delta_phi_deta4_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_gap_inclusive",jet_delta_phi_gap);
   if (jet_delta_phi_gap == 0) { cout << "trigger_syst_delta_phi_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta1_gap_inclusive",jet_delta_phi_deta1_gap);
   if (jet_delta_phi_deta1_gap == 0) { cout << "trigger_syst_delta_phi_deta1_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta2_gap_inclusive",jet_delta_phi_deta2_gap);
   if (jet_delta_phi_deta2_gap == 0) { cout << "trigger_syst_delta_phi_deta2_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta3_gap_inclusive",jet_delta_phi_deta3_gap);
   if (jet_delta_phi_deta3_gap == 0) { cout << "trigger_syst_delta_phi_deta3_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta4_gap_inclusive",jet_delta_phi_deta4_gap);
   if (jet_delta_phi_deta4_gap == 0) { cout << "trigger_syst_delta_phi_deta4_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_nogap_inclusive",jet_delta_phi_nogap);
   if (jet_delta_phi_nogap == 0) { cout << "trigger_syst_delta_phi_nogap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta1_nogap_inclusive",jet_delta_phi_deta1_nogap);
   if (jet_delta_phi_deta1_nogap == 0) { cout << "trigger_syst_delta_phi_deta1_nogap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta2_nogap_inclusive",jet_delta_phi_deta2_nogap);
   if (jet_delta_phi_deta2_nogap == 0) { cout << "trigger_syst_delta_phi_deta2_nogap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta3_nogap_inclusive",jet_delta_phi_deta3_nogap);
   if (jet_delta_phi_deta3_nogap == 0) { cout << "trigger_syst_delta_phi_deta3_nogap_exclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_phi_deta4_nogap_inclusive",jet_delta_phi_deta4_nogap);
   if (jet_delta_phi_deta4_nogap == 0) { cout << "trigger_syst_delta_phi_deta4_nogap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_leading_pt_inside_gap_inclusive",jet_leading_pt_inside_gap);   
   if (jet_leading_pt_inside_gap == 0) { cout << "trigger_syst_leading_pt_inside_gap_exclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_leading_eta_star_inside_gap_inclusive",jet_leading_eta_star_inside_gap);   
   if (jet_leading_eta_star_inside_gap == 0) { cout << "trigger_syst_leading_eta_star_inside_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("trigger_syst_delta_eta_outside_gap_inclusive",jet_delta_eta_outside_gap);   
   if (jet_delta_eta_outside_gap == 0) { cout << "trigger_syst_delta_eta_outside_gap_inclusive on Jet_2010B not found!" << endl; not_loaded = true; }

    //Summary regarding the loading of the histograms
    if (detail && !not_loaded) { cout<<"All histograms loaded sucessfully!"<<endl; }
    if (not_loaded)
	{
	cout<<"Some histograms were not loaded! Please check the input rootfile. The routine will terminate now!" << endl;
	return;
	}

   //declaring the merging histograms
   TH1D *merged_delta_phi;
   TH1D *merged_delta_phi_gap;
   TH1D *merged_delta_phi_nogap;
   TH1D *merged_delta_phi_deta1;
   TH1D *merged_delta_phi_deta2;
   TH1D *merged_delta_phi_deta3;
   TH1D *merged_delta_phi_deta4;
   TH1D *merged_delta_phi_deta1_gap;
   TH1D *merged_delta_phi_deta2_gap;
   TH1D *merged_delta_phi_deta3_gap;
   TH1D *merged_delta_phi_deta4_gap;
   TH1D *merged_delta_phi_deta1_nogap;
   TH1D *merged_delta_phi_deta2_nogap;
   TH1D *merged_delta_phi_deta3_nogap;
   TH1D *merged_delta_phi_deta4_nogap;
   TH1D *merged_leading_pt_inside_gap;
   TH1D *merged_leading_eta_star_inside_gap;
   TH1D *merged_delta_eta_outside_gap;

   //setup the merged histogram binning
   merged_delta_phi =  new TH1D("ak5PF_delta_phi","#Delta#phi;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","#Delta#phi when requiring a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","#Delta#phi when vetoing a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta1 =  new TH1D("ak5PF_delta_phi_deta1","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta2 =  new TH1D("ak5PF_delta_phi_deta2","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta3 =  new TH1D("ak5PF_delta_phi_deta3","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta4 =  new TH1D("ak5PF_delta_phi_deta4","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta1_gap =  new TH1D("ak5PF_delta_phi_deta1_gap","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta2_gap =  new TH1D("ak5PF_delta_phi_deta2_gap","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta3_gap =  new TH1D("ak5PF_delta_phi_deta3_gap","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta4_gap =  new TH1D("ak5PF_delta_phi_deta4_gap","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta1_nogap =  new TH1D("ak5PF_delta_phi_deta1_nogap","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta2_nogap =  new TH1D("ak5PF_delta_phi_deta2_nogap","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta3_nogap =  new TH1D("ak5PF_delta_phi_deta3_nogap","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_delta_phi_deta4_nogap =  new TH1D("ak5PF_delta_phi_deta4_nogap","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;Combination Uncertainty", dphi_nbins, dphi_bins);
   merged_leading_pt_inside_gap =  new TH1D("ak5PF_leading_pt_inside_gap","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];Combination Uncertainty", in_nbins, in_bins);
   merged_delta_eta_outside_gap =  new TH1D("ak5PF_delta_eta_outside_gap","#Delta#eta outside the gap;#Delta#eta;Combination Uncertainty", deta_out_nbins, deta_out_bins);
   merged_leading_eta_star_inside_gap =  new TH1D("ak5PF_leading_eta_star_inside_gap","Leading Jet #eta^{*} inside the gap;#eta^{*};Combination Uncertainty", etastar_nbins, etastar_bins);

     //initializing statistics
     merged_delta_phi->Sumw2();
     merged_delta_phi_gap->Sumw2();
     merged_delta_phi_nogap->Sumw2();
     merged_delta_phi_deta1->Sumw2();
     merged_delta_phi_deta2->Sumw2();
     merged_delta_phi_deta3->Sumw2();
     merged_delta_phi_deta4->Sumw2();
     merged_delta_phi_deta1_gap->Sumw2();
     merged_delta_phi_deta2_gap->Sumw2();
     merged_delta_phi_deta3_gap->Sumw2();
     merged_delta_phi_deta4_gap->Sumw2();
     merged_delta_phi_deta1_nogap->Sumw2();
     merged_delta_phi_deta2_nogap->Sumw2();
     merged_delta_phi_deta3_nogap->Sumw2();
     merged_delta_phi_deta4_nogap->Sumw2();
     merged_leading_pt_inside_gap->Sumw2();
     merged_leading_eta_star_inside_gap->Sumw2();
     merged_delta_eta_outside_gap->Sumw2();

      //combining and ploting
      merge_and_plot_3histograms(jetmettau_delta_phi, jetmet_delta_phi, jet_delta_phi, merged_delta_phi, output_path, "delta_phi", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_gap, jetmet_delta_phi_gap, jet_delta_phi_gap, merged_delta_phi_gap, output_path, "delta_phi_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_nogap, jetmet_delta_phi_nogap, jet_delta_phi_nogap, merged_delta_phi_nogap, output_path, "delta_phi_nogap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta1, jetmet_delta_phi_deta1, jet_delta_phi_deta1, merged_delta_phi_deta1, output_path, "delta_phi_deta1", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta2, jetmet_delta_phi_deta2, jet_delta_phi_deta2, merged_delta_phi_deta2, output_path, "delta_phi_deta2", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta3, jetmet_delta_phi_deta3, jet_delta_phi_deta3, merged_delta_phi_deta3, output_path, "delta_phi_deta3", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta4, jetmet_delta_phi_deta4, jet_delta_phi_deta4, merged_delta_phi_deta4, output_path, "delta_phi_deta4", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta1_gap, jetmet_delta_phi_deta1_gap, jet_delta_phi_deta1_gap, merged_delta_phi_deta1_gap, output_path, "delta_phi_deta1_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta2_gap, jetmet_delta_phi_deta2_gap, jet_delta_phi_deta2_gap, merged_delta_phi_deta2_gap, output_path, "delta_phi_deta2_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta3_gap, jetmet_delta_phi_deta3_gap, jet_delta_phi_deta3_gap, merged_delta_phi_deta3_gap, output_path, "delta_phi_deta3_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta4_gap, jetmet_delta_phi_deta4_gap, jet_delta_phi_deta4_gap, merged_delta_phi_deta4_gap, output_path, "delta_phi_deta4_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta1_nogap, jetmet_delta_phi_deta1_nogap, jet_delta_phi_deta1_nogap, merged_delta_phi_deta1_nogap, output_path, "delta_phi_deta1_nogap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta2_nogap, jetmet_delta_phi_deta2_nogap, jet_delta_phi_deta2_nogap, merged_delta_phi_deta2_nogap, output_path, "delta_phi_deta2_nogap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta3_nogap, jetmet_delta_phi_deta3_nogap, jet_delta_phi_deta3_nogap, merged_delta_phi_deta3_nogap, output_path, "delta_phi_deta3_nogap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_phi_deta4_nogap, jetmet_delta_phi_deta4_nogap, jet_delta_phi_deta4_nogap, merged_delta_phi_deta4_nogap, output_path, "delta_phi_deta4_nogap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_leading_pt_inside_gap, jetmet_leading_pt_inside_gap, jet_leading_pt_inside_gap, merged_leading_pt_inside_gap, output_path, "leading_pt_inside_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_leading_eta_star_inside_gap, jetmet_leading_eta_star_inside_gap, jet_leading_eta_star_inside_gap, merged_leading_eta_star_inside_gap, output_path, "delta_leading_eta_star_inside_gap", "top_right", detail);
      merge_and_plot_3histograms(jetmettau_delta_eta_outside_gap, jetmet_delta_eta_outside_gap, jet_delta_eta_outside_gap, merged_delta_eta_outside_gap, output_path, "delta_delta_eta_outside_gap", "top_right", detail);


     //Open the output root file
     TFile *data_output = TFile::Open( data_out.c_str() , "RECREATE");

     if (detail) { cout<<"Writing histograms on " << data_out << endl; }

     //writing histograms on file
     merged_delta_phi->Write();
     merged_delta_phi_gap->Write();
     merged_delta_phi_nogap->Write();
     merged_delta_phi_deta1->Write();
     merged_delta_phi_deta2->Write();
     merged_delta_phi_deta3->Write();
     merged_delta_phi_deta4->Write();
     merged_delta_phi_deta1_gap->Write();
     merged_delta_phi_deta2_gap->Write();
     merged_delta_phi_deta3_gap->Write();
     merged_delta_phi_deta4_gap->Write();
     merged_delta_phi_deta1_nogap->Write();
     merged_delta_phi_deta2_nogap->Write();
     merged_delta_phi_deta3_nogap->Write();
     merged_delta_phi_deta4_nogap->Write();
     merged_leading_pt_inside_gap->Write();
     merged_leading_eta_star_inside_gap->Write();
     merged_delta_eta_outside_gap->Write();

     if (detail) { cout<<"Histograms written sucessfully!"<<endl; }
     
     //close the output file
     data_output->Close();

     //deleting the histograms
     delete(merged_delta_phi);
     delete(merged_delta_phi_gap);
     delete(merged_delta_phi_nogap);
     delete(merged_delta_phi_deta1);
     delete(merged_delta_phi_deta2);
     delete(merged_delta_phi_deta3);
     delete(merged_delta_phi_deta4);
     delete(merged_delta_phi_deta1_gap);
     delete(merged_delta_phi_deta2_gap);
     delete(merged_delta_phi_deta3_gap);
     delete(merged_delta_phi_deta4_gap);
     delete(merged_delta_phi_deta1_nogap);
     delete(merged_delta_phi_deta2_nogap);
     delete(merged_delta_phi_deta3_nogap);
     delete(merged_delta_phi_deta4_nogap);
     delete(merged_leading_pt_inside_gap);
     delete(merged_leading_eta_star_inside_gap);
     delete(merged_delta_eta_outside_gap);

}

void estimate_combination_systematic(string *data_in, string data_out, int n_files, string sel_mode = "allvertex", bool detail = false, bool test = false)
{

   double pt_min = 35.0;
   double pt_min_gap = 20.0;
   int HLTJetPtN[3] = {15,30,50};
   int nJetTrig = 4;

//output the configuration
   if (detail) { cout<<"Estimate Combination Systematic Configuration"<<endl; }
   if (detail) { cout<<"Selection mode :  "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files : "<<n_files<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :     "<<data_out<<endl; }
   if (detail) { cout<<"Pt Min :          "<<pt_min<<endl; }
   if (detail) { cout<<"Detail level :    "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :       "<<test<<endl; }

//declaring the variables
   TString HLTJet[nJetTrig];

   int ihltj[nJetTrig];
   int prescalej[nJetTrig];
   int nentries = 0;
   int counter_entries = 0, counter_pv = 0, counter_hlt = 0, counter_selected = 0, counter_jet = 0;
   int triggered[7], selected[7];
   int events_selected[7][15];
   int index, index1, index2;

   char trigtitle[200];
   
   bool pv_pass = false;
   bool hltPass = false;
   bool emulation_pass = false;
   bool hltPassj[nJetTrig];
   bool hard_emission;
   
   double pt, eta, phi;
   double leading_pt, leading_eta, leading_phi;
   double forward_pt, forward_eta, forward_phi;
   double central_pt, central_eta, central_phi;
   double gap_leading_pt, gap_eta, gap_phi, gap_eta_star;
   double outside_leading_pt, outside_eta, outside_phi, outside_delta_eta;
   double deta_out1, deta_out2;
   double delta_eta, delta_phi;
   double events_loss[2][15];
   double prescale, prescale_old, prescale_new, pu_scale, eff, trigger_loss1, trigger_loss2;
   double probs[nJetTrig];

//declare the main binning
   int dphi_nbins = 7;
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

   int in_nbins = 9;
   double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

   int etastar_nbins = 12;
   double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

   int deta_out_nbins = 6;
   double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//declaring histograms
     TH1D *hist_events;
     TH1D *hist_events_HLT_Jet15U;
     TH1D *hist_events_HLT_Jet30U;
     TH1D *hist_events_HLT_Jet50U;
     TH1D *hist_events_all;
     TH1D *hist_events_inclusive_zones;
     TH1D *hist_events_exclusive_zones;

     TH1D *hist_delta_phi;
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
     TH1D *hist_leading_pt_inside_gap;
     TH1D *hist_leading_eta_star_inside_gap;
     TH1D *hist_delta_eta_outside_gap;

     TH1D *hist_delta_phi_inclusive;
     TH1D *hist_delta_phi_gap_inclusive;
     TH1D *hist_delta_phi_nogap_inclusive;
     TH1D *hist_delta_phi_deta1_inclusive;
     TH1D *hist_delta_phi_deta2_inclusive;
     TH1D *hist_delta_phi_deta3_inclusive;
     TH1D *hist_delta_phi_deta4_inclusive;
     TH1D *hist_delta_phi_deta1_gap_inclusive;
     TH1D *hist_delta_phi_deta2_gap_inclusive;
     TH1D *hist_delta_phi_deta3_gap_inclusive;
     TH1D *hist_delta_phi_deta4_gap_inclusive;
     TH1D *hist_delta_phi_deta1_nogap_inclusive;
     TH1D *hist_delta_phi_deta2_nogap_inclusive;
     TH1D *hist_delta_phi_deta3_nogap_inclusive;
     TH1D *hist_delta_phi_deta4_nogap_inclusive;
     TH1D *hist_leading_pt_inside_gap_inclusive;
     TH1D *hist_leading_eta_star_inside_gap_inclusive;
     TH1D *hist_delta_eta_outside_gap_inclusive;

     TH1D *hist_delta_phi_exclusive;
     TH1D *hist_delta_phi_gap_exclusive;
     TH1D *hist_delta_phi_nogap_exclusive;
     TH1D *hist_delta_phi_deta1_exclusive;
     TH1D *hist_delta_phi_deta2_exclusive;
     TH1D *hist_delta_phi_deta3_exclusive;
     TH1D *hist_delta_phi_deta4_exclusive;
     TH1D *hist_delta_phi_deta1_gap_exclusive;
     TH1D *hist_delta_phi_deta2_gap_exclusive;
     TH1D *hist_delta_phi_deta3_gap_exclusive;
     TH1D *hist_delta_phi_deta4_gap_exclusive;
     TH1D *hist_delta_phi_deta1_nogap_exclusive;
     TH1D *hist_delta_phi_deta2_nogap_exclusive;
     TH1D *hist_delta_phi_deta3_nogap_exclusive;
     TH1D *hist_delta_phi_deta4_nogap_exclusive;
     TH1D *hist_leading_pt_inside_gap_exclusive;
     TH1D *hist_leading_eta_star_inside_gap_exclusive;
     TH1D *hist_delta_eta_outside_gap_exclusive;

     TH1D *syst_delta_phi_inclusive;
     TH1D *syst_delta_phi_gap_inclusive;
     TH1D *syst_delta_phi_nogap_inclusive;
     TH1D *syst_delta_phi_deta1_inclusive;
     TH1D *syst_delta_phi_deta2_inclusive;
     TH1D *syst_delta_phi_deta3_inclusive;
     TH1D *syst_delta_phi_deta4_inclusive;
     TH1D *syst_delta_phi_deta1_gap_inclusive;
     TH1D *syst_delta_phi_deta2_gap_inclusive;
     TH1D *syst_delta_phi_deta3_gap_inclusive;
     TH1D *syst_delta_phi_deta4_gap_inclusive;
     TH1D *syst_delta_phi_deta1_nogap_inclusive;
     TH1D *syst_delta_phi_deta2_nogap_inclusive;
     TH1D *syst_delta_phi_deta3_nogap_inclusive;
     TH1D *syst_delta_phi_deta4_nogap_inclusive;
     TH1D *syst_leading_pt_inside_gap_inclusive;
     TH1D *syst_leading_eta_star_inside_gap_inclusive;
     TH1D *syst_delta_eta_outside_gap_inclusive;

     TH1D *syst_delta_phi_exclusive;
     TH1D *syst_delta_phi_gap_exclusive;
     TH1D *syst_delta_phi_nogap_exclusive;
     TH1D *syst_delta_phi_deta1_exclusive;
     TH1D *syst_delta_phi_deta2_exclusive;
     TH1D *syst_delta_phi_deta3_exclusive;
     TH1D *syst_delta_phi_deta4_exclusive;
     TH1D *syst_delta_phi_deta1_gap_exclusive;
     TH1D *syst_delta_phi_deta2_gap_exclusive;
     TH1D *syst_delta_phi_deta3_gap_exclusive;
     TH1D *syst_delta_phi_deta4_gap_exclusive;
     TH1D *syst_delta_phi_deta1_nogap_exclusive;
     TH1D *syst_delta_phi_deta2_nogap_exclusive;
     TH1D *syst_delta_phi_deta3_nogap_exclusive;
     TH1D *syst_delta_phi_deta4_nogap_exclusive;
     TH1D *syst_leading_pt_inside_gap_exclusive;
     TH1D *syst_leading_eta_star_inside_gap_exclusive;
     TH1D *syst_delta_eta_outside_gap_exclusive;

//declaring histogram binning
     hist_events =  new TH1D("Events","Selection Chain;Type;# Events", 7, 0, 7);
     hist_events_HLT_Jet15U =  new TH1D("Events_HLT_Jet15U","Selection Chain for HLT_Jet15U;Type;# Events", 2, 0, 2);
     hist_events_HLT_Jet30U =  new TH1D("Events_HLT_Jet30U","Selection Chain for HLT_Jet30U;Type;# Events", 2, 0, 2);
     hist_events_HLT_Jet50U =  new TH1D("Events_HLT_Jet50U","Selection Chain for HLT_Jet50U;Type;# Events", 2, 0, 2);
     hist_events_all =  new TH1D("Events_Combined","Selection Chain for combination;Type;# Events", 2, 0, 2);
     hist_events_inclusive_zones =  new TH1D("Events_Inclusive_Zones","Selection Chain for combination Inclusive Zones;Type;# Events", 3, 0, 3);
     hist_events_exclusive_zones =  new TH1D("Events_Exclusive_Zones","Selection Chain for combination Exclusive Zones;Type;# Events", 3, 0, 3);

  hist_delta_phi =  new TH1D("ak5PF_delta_phi","#Delta#phi;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","#Delta#phi when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","#Delta#phi when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1 =  new TH1D("ak5PF_delta_phi_deta1","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2 =  new TH1D("ak5PF_delta_phi_deta2","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3 =  new TH1D("ak5PF_delta_phi_deta3","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4 =  new TH1D("ak5PF_delta_phi_deta4","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_gap =  new TH1D("ak5PF_delta_phi_deta1_gap","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_gap =  new TH1D("ak5PF_delta_phi_deta2_gap","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_gap =  new TH1D("ak5PF_delta_phi_deta3_gap","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_gap =  new TH1D("ak5PF_delta_phi_deta4_gap","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_nogap =  new TH1D("ak5PF_delta_phi_deta1_nogap","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]i", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_nogap =  new TH1D("ak5PF_delta_phi_deta2_nogap","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_nogap =  new TH1D("ak5PF_delta_phi_deta3_nogap","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_nogap =  new TH1D("ak5PF_delta_phi_deta4_nogap","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_leading_pt_inside_gap =  new TH1D("ak5PF_leading_pt_inside_gap","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  hist_delta_eta_outside_gap =  new TH1D("ak5PF_delta_eta_outside_gap","#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
  hist_leading_eta_star_inside_gap =  new TH1D("ak5PF_leading_eta_star_inside_gap","Leading Jet #eta^{*} inside the gap;#eta^{*};#frac{d#sigma}{d#eta^{*}} [pb]", etastar_nbins, etastar_bins);

  hist_delta_phi_inclusive =  new TH1D("ak5PF_delta_phi_inclusive","#Delta#phi;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_gap_inclusive =  new TH1D("ak5PF_delta_phi_gap_inclusive","#Delta#phi when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_nogap_inclusive =  new TH1D("ak5PF_delta_phi_nogap_inclusive","#Delta#phi when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_inclusive =  new TH1D("ak5PF_delta_phi_deta1_inclusive","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_inclusive =  new TH1D("ak5PF_delta_phi_deta2_inclusive","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_inclusive =  new TH1D("ak5PF_delta_phi_deta3_inclusive","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_inclusive =  new TH1D("ak5PF_delta_phi_deta4_inclusive","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_gap_inclusive =  new TH1D("ak5PF_delta_phi_deta1_gap_inclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_gap_inclusive =  new TH1D("ak5PF_delta_phi_deta2_gap_inclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_gap_inclusive =  new TH1D("ak5PF_delta_phi_deta3_gap_inclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_gap_inclusive =  new TH1D("ak5PF_delta_phi_deta4_gap_inclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_nogap_inclusive =  new TH1D("ak5PF_delta_phi_deta1_nogap_inclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]i", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_nogap_inclusive =  new TH1D("ak5PF_delta_phi_deta2_nogap_inclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_nogap_inclusive =  new TH1D("ak5PF_delta_phi_deta3_nogap_inclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_nogap_inclusive =  new TH1D("ak5PF_delta_phi_deta4_nogap_inclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_leading_pt_inside_gap_inclusive =  new TH1D("ak5PF_leading_pt_inside_gap_inclusive","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  hist_delta_eta_outside_gap_inclusive =  new TH1D("ak5PF_delta_eta_outside_gap_inclusive","#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
  hist_leading_eta_star_inside_gap_inclusive =  new TH1D("ak5PF_leading_eta_star_inside_gap_inclusive","Leading Jet #eta^{*} inside the gap;#eta^{*};#frac{d#sigma}{d#eta^{*}} [pb]", etastar_nbins, etastar_bins);

  hist_delta_phi_exclusive =  new TH1D("ak5PF_delta_phi_exclusive","#Delta#phi;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_gap_exclusive =  new TH1D("ak5PF_delta_phi_gap_exclusive","#Delta#phi when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_nogap_exclusive =  new TH1D("ak5PF_delta_phi_nogap_exclusive","#Delta#phi when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_exclusive =  new TH1D("ak5PF_delta_phi_deta1_exclusive","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_exclusive =  new TH1D("ak5PF_delta_phi_deta2_exclusive","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_exclusive =  new TH1D("ak5PF_delta_phi_deta3_exclusive","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_exclusive =  new TH1D("ak5PF_delta_phi_deta4_exclusive","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_gap_exclusive =  new TH1D("ak5PF_delta_phi_deta1_gap_exclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_gap_exclusive =  new TH1D("ak5PF_delta_phi_deta2_gap_exclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_gap_exclusive =  new TH1D("ak5PF_delta_phi_deta3_gap_exclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_gap_exclusive =  new TH1D("ak5PF_delta_phi_deta4_gap_exclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta1_nogap_exclusive =  new TH1D("ak5PF_delta_phi_deta1_nogap_exclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]i", dphi_nbins, dphi_bins);
  hist_delta_phi_deta2_nogap_exclusive =  new TH1D("ak5PF_delta_phi_deta2_nogap_exclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta3_nogap_exclusive =  new TH1D("ak5PF_delta_phi_deta3_nogap_exclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_delta_phi_deta4_nogap_exclusive =  new TH1D("ak5PF_delta_phi_deta4_nogap_exclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;#frac{d#sigma}{d#Delta#phid#Delta#eta} [pb]", dphi_nbins, dphi_bins);
  hist_leading_pt_inside_gap_exclusive =  new TH1D("ak5PF_leading_pt_inside_gap_exclusive","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  hist_delta_eta_outside_gap_exclusive =  new TH1D("ak5PF_delta_eta_outside_gap_exclusive","#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
  hist_leading_eta_star_inside_gap_exclusive =  new TH1D("ak5PF_leading_eta_star_exside_gap_inclusive","Leading Jet #eta^{*} inside the gap;#eta^{*};#frac{d#sigma}{d#eta^{*}} [pb]", etastar_nbins, etastar_bins);

  syst_delta_phi_inclusive =  new TH1D("trigger_syst_delta_phi_inclusive","#Delta#phi;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_gap_inclusive =  new TH1D("trigger_syst_delta_phi_gap_inclusive","#Delta#phi when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_nogap_inclusive =  new TH1D("trigger_syst_delta_phi_nogap_inclusive","#Delta#phi when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta1_inclusive =  new TH1D("trigger_syst_delta_phi_deta1_inclusive","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta2_inclusive =  new TH1D("trigger_syst_delta_phi_deta2_inclusive","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta3_inclusive =  new TH1D("trigger_syst_delta_phi_deta3_inclusive","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta4_inclusive =  new TH1D("trigger_syst_delta_phi_deta4_inclusive","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta1_gap_inclusive =  new TH1D("trigger_syst_delta_phi_deta1_gap_inclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta2_gap_inclusive =  new TH1D("trigger_syst_delta_phi_deta2_gap_inclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta3_gap_inclusive =  new TH1D("trigger_syst_delta_phi_deta3_gap_inclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta4_gap_inclusive =  new TH1D("trigger_syst_delta_phi_deta4_gap_inclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta1_nogap_inclusive =  new TH1D("trigger_syst_delta_phi_deta1_nogap_inclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta2_nogap_inclusive =  new TH1D("trigger_syst_delta_phi_deta2_nogap_inclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta3_nogap_inclusive =  new TH1D("trigger_syst_delta_phi_deta3_nogap_inclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta4_nogap_inclusive =  new TH1D("trigger_syst_delta_phi_deta4_nogap_inclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_leading_pt_inside_gap_inclusive =  new TH1D("trigger_syst_leading_pt_inside_gap_inclusive","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  syst_delta_eta_outside_gap_inclusive =  new TH1D("trigger_syst_delta_eta_outside_gap_inclusive","#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
  syst_leading_eta_star_inside_gap_inclusive =  new TH1D("trigger_syst_leading_eta_star_inside_gap_inclusive","Leading Jet #eta^{*} inside the gap;#eta^{*};#frac{d#sigma}{d#eta^{*}} [pb]", etastar_nbins, etastar_bins);


  syst_delta_phi_exclusive =  new TH1D("trigger_syst_delta_phi_exclusive","#Delta#phi;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_gap_exclusive =  new TH1D("trigger_syst_delta_phi_gap_exclusive","#Delta#phi when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_nogap_exclusive =  new TH1D("trigger_syst_delta_phi_nogap_exclusive","#Delta#phi when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta1_exclusive =  new TH1D("trigger_syst_delta_phi_deta1_exclusive","#Delta#phi for 0.4 >= Delta Eta > 2.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta2_exclusive =  new TH1D("trigger_syst_delta_phi_deta2_exclusive","#Delta#phi for 2.5 >= Delta Eta > 3.5;|#Delta#phi|;Trigger Systematic]", dphi_nbins, dphi_bins);
  syst_delta_phi_deta3_exclusive =  new TH1D("trigger_syst_delta_phi_deta3_exclusive","#Delta#phi for 3.5 >= Delta Eta > 4.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta4_exclusive =  new TH1D("trigger_syst_delta_phi_deta4_exclusive","#Delta#phi for 4.5 >= Delta Eta > 7.5;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta1_gap_exclusive =  new TH1D("trigger_syst_delta_phi_deta1_gap_exclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta2_gap_exclusive =  new TH1D("trigger_syst_delta_phi_deta2_gap_exclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta3_gap_exclusive =  new TH1D("trigger_syst_delta_phi_deta3_gap_exclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta4_gap_exclusive =  new TH1D("trigger_syst_delta_phi_deta4_gap_exclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when requiring a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta1_nogap_exclusive =  new TH1D("trigger_syst_delta_phi_deta1_nogap_exclusive","#Delta#phi for 0.4 >= #Delta#eta > 2.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta2_nogap_exclusive =  new TH1D("trigger_syst_delta_phi_deta2_nogap_exclusive","#Delta#phi for 2.5 >= #Delta#eta > 3.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta3_nogap_exclusive =  new TH1D("trigger_syst_delta_phi_deta3_nogap_exclusive","#Delta#phi for 3.5 >= #Delta#eta > 4.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_delta_phi_deta4_nogap_exclusive =  new TH1D("trigger_syst_delta_phi_deta4_nogap_exclusive","#Delta#phi for 4.5 >= #Delta#eta > 7.5 when vetoing a gap;|#Delta#phi|;Trigger Systematic", dphi_nbins, dphi_bins);
  syst_leading_pt_inside_gap_exclusive =  new TH1D("trigger_syst_leading_pt_inside_gap_exclusive","Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
  syst_delta_eta_outside_gap_exclusive =  new TH1D("trigger_syst_delta_eta_outside_gap_exclusive","#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
  syst_leading_eta_star_inside_gap_exclusive =  new TH1D("trigger_syst_leading_eta_star_inside_gap_exclusive","Leading Jet #eta^{*} inside the gap;#eta^{*};#frac{d#sigma}{d#eta^{*}} [pb]", etastar_nbins, etastar_bins);

//inicializing the histogram statistics
     hist_events->Sumw2();
     hist_events->Fill("Total Events",0);
     hist_events->Fill("PV Selection",0);
     hist_events->Fill("Trigger Selection",0);
     hist_events->Fill("Selected",0);
     hist_events->Fill("Number of Jets",0);
     hist_events->Fill("Pileup Scale",0);
     hist_events->Fill("Efficiency",0);

     hist_events_HLT_Jet15U->Sumw2();
     hist_events_HLT_Jet15U->Fill("Triggered Events",0);
     hist_events_HLT_Jet15U->Fill("Selected",0);

     hist_events_HLT_Jet30U->Sumw2();
     hist_events_HLT_Jet30U->Fill("Triggered Events",0);
     hist_events_HLT_Jet30U->Fill("Selected",0);

     hist_events_HLT_Jet50U->Sumw2();
     hist_events_HLT_Jet50U->Fill("Triggered Events",0);
     hist_events_HLT_Jet50U->Fill("Selected",0);

     hist_events_all->Sumw2();
     hist_events_all->Fill("Triggered Events",0);
     hist_events_all->Fill("Selected",0);

     hist_events_inclusive_zones->Sumw2();
     hist_events_inclusive_zones->Fill("Triggered Events",0);
     hist_events_inclusive_zones->Fill("Selected",0);
     hist_events_inclusive_zones->Fill("Proportion",0);

     hist_events_exclusive_zones->Sumw2();
     hist_events_exclusive_zones->Fill("Triggered Events",0);
     hist_events_exclusive_zones->Fill("Selected",0);
     hist_events_exclusive_zones->Fill("Proportion",0);

     hist_delta_phi->Sumw2();
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
     hist_leading_pt_inside_gap->Sumw2();
     hist_leading_eta_star_inside_gap->Sumw2();
     hist_delta_eta_outside_gap->Sumw2();

     hist_delta_phi_inclusive->Sumw2();
     hist_delta_phi_gap_inclusive->Sumw2();
     hist_delta_phi_nogap_inclusive->Sumw2();
     hist_delta_phi_deta1_inclusive->Sumw2();
     hist_delta_phi_deta2_inclusive->Sumw2();
     hist_delta_phi_deta3_inclusive->Sumw2();
     hist_delta_phi_deta4_inclusive->Sumw2();
     hist_delta_phi_deta1_gap_inclusive->Sumw2();
     hist_delta_phi_deta2_gap_inclusive->Sumw2();
     hist_delta_phi_deta3_gap_inclusive->Sumw2();
     hist_delta_phi_deta4_gap_inclusive->Sumw2();
     hist_delta_phi_deta1_nogap_inclusive->Sumw2();
     hist_delta_phi_deta2_nogap_inclusive->Sumw2();
     hist_delta_phi_deta3_nogap_inclusive->Sumw2();
     hist_delta_phi_deta4_nogap_inclusive->Sumw2();
     hist_leading_pt_inside_gap_inclusive->Sumw2();
     hist_leading_eta_star_inside_gap_inclusive->Sumw2();
     hist_delta_eta_outside_gap_inclusive->Sumw2();

     hist_delta_phi_exclusive->Sumw2();
     hist_delta_phi_gap_exclusive->Sumw2();
     hist_delta_phi_nogap_exclusive->Sumw2();
     hist_delta_phi_deta1_exclusive->Sumw2();
     hist_delta_phi_deta2_exclusive->Sumw2();
     hist_delta_phi_deta3_exclusive->Sumw2();
     hist_delta_phi_deta4_exclusive->Sumw2();
     hist_delta_phi_deta1_gap_exclusive->Sumw2();
     hist_delta_phi_deta2_gap_exclusive->Sumw2();
     hist_delta_phi_deta3_gap_exclusive->Sumw2();
     hist_delta_phi_deta4_gap_exclusive->Sumw2();
     hist_delta_phi_deta1_nogap_exclusive->Sumw2();
     hist_delta_phi_deta2_nogap_exclusive->Sumw2();
     hist_delta_phi_deta3_nogap_exclusive->Sumw2();
     hist_delta_phi_deta4_nogap_exclusive->Sumw2();
     hist_leading_pt_inside_gap_exclusive->Sumw2();
     hist_leading_eta_star_inside_gap_exclusive->Sumw2();
     hist_delta_eta_outside_gap_exclusive->Sumw2();

     syst_delta_phi_inclusive->Sumw2();
     syst_delta_phi_gap_inclusive->Sumw2();
     syst_delta_phi_nogap_inclusive->Sumw2();
     syst_delta_phi_deta1_inclusive->Sumw2();
     syst_delta_phi_deta2_inclusive->Sumw2();
     syst_delta_phi_deta3_inclusive->Sumw2();
     syst_delta_phi_deta4_inclusive->Sumw2();
     syst_delta_phi_deta1_gap_inclusive->Sumw2();
     syst_delta_phi_deta2_gap_inclusive->Sumw2();
     syst_delta_phi_deta3_gap_inclusive->Sumw2();
     syst_delta_phi_deta4_gap_inclusive->Sumw2();
     syst_delta_phi_deta1_nogap_inclusive->Sumw2();
     syst_delta_phi_deta2_nogap_inclusive->Sumw2();
     syst_delta_phi_deta3_nogap_inclusive->Sumw2();
     syst_delta_phi_deta4_nogap_inclusive->Sumw2();
     syst_leading_pt_inside_gap_inclusive->Sumw2();
     syst_leading_eta_star_inside_gap_inclusive->Sumw2();
     syst_delta_eta_outside_gap_inclusive->Sumw2();

     syst_delta_phi_exclusive->Sumw2();
     syst_delta_phi_gap_exclusive->Sumw2();
     syst_delta_phi_nogap_exclusive->Sumw2();
     syst_delta_phi_deta1_exclusive->Sumw2();
     syst_delta_phi_deta2_exclusive->Sumw2();
     syst_delta_phi_deta3_exclusive->Sumw2();
     syst_delta_phi_deta4_exclusive->Sumw2();
     syst_delta_phi_deta1_gap_exclusive->Sumw2();
     syst_delta_phi_deta2_gap_exclusive->Sumw2();
     syst_delta_phi_deta3_gap_exclusive->Sumw2();
     syst_delta_phi_deta4_gap_exclusive->Sumw2();
     syst_delta_phi_deta1_nogap_exclusive->Sumw2();
     syst_delta_phi_deta2_nogap_exclusive->Sumw2();
     syst_delta_phi_deta3_nogap_exclusive->Sumw2();
     syst_delta_phi_deta4_nogap_exclusive->Sumw2();
     syst_leading_pt_inside_gap_exclusive->Sumw2();
     syst_leading_eta_star_inside_gap_exclusive->Sumw2();
     syst_delta_eta_outside_gap_exclusive->Sumw2();


for (int m = 0; m < 15; m++)
	{
	for (int y = 0; y < 6; y++)
		{
		events_selected[y][m] = 0;
		}
	events_loss[0][m] = 0.0;
	events_loss[1][m] = 0.0;
	}

//loop over the files
for (int z = 0; z < n_files; z++)
{
//open the file
  string file = data_in[z];
  cout<<z+1<<"/"<<n_files<<" -> "<<file<<endl;
  TFile *inf = 0;
  inf = TFile::Open( file.c_str() );
  if (inf == 0) { cout << "Ntuple not loaded!" << endl; return; }

//define the tree branch
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);

   //------------------------- Initializing The Trigger Variables -------------------------- //
      cout << "Trigger list : " << endl;
      for (int i = 0; i<=nJetTrig+3; i++)
         {
                triggered[i] = 0;
		selected[i] = 0;
         }

      for (int i = 0; i<nJetTrig-1; i++)
         {
                sprintf(trigtitle, "HLT_Jet%iU",HLTJetPtN[i]);             
		HLTJet[i] = trigtitle;
		//cout << "Testing " << trigtitle << endl;
         }
	HLTJet[3] = "HLT_L1Jet6U";
      for (int i = 0; i<nJetTrig; i++)
         {
                  ihltj[i] = -1 ;
         }

   TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
   // ------------ Assigning An Integer To Each Trigger Path---------------//
    for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++)
	{
       	TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));

       	for (int ii=0; ii<nJetTrig; ii++)
          {
          if (ss == HLTJet[ii])
		{
                ihltj[ii] = ibin;
                continue;
          	}
          } // for (int ii=0; ii<nJetTrig; ii++)
        } // for(int ibin=0;ibin<55;ibin++)

  // ----------------------- Checking For The Trigger Assignment ------------------- // 
    for (int ij=0; ij<nJetTrig; ij++)
        {
              if (ihltj[ij] == -1)
		{
              	cout<<"The requested trigger ("<<HLTJet[ij]<<") is not found "<<endl;
   //         	break;
               	}
              else
		{
              	cout<<HLTJet[ij]<<" --> "<<ihltj[ij]<<endl;
             	}
      	}
      
  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  int decade = 0;

  if (test) { nentries = 1000000; } //reduced number of read entries, usefull for testing

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
      
  //-------------------------- Initializing the Boolians and the Prescale Values----------------------//
        for (int j=0; j<=nJetTrig+2; j++)
        {
                   hltPassj[j] = false ;
        } // for (int i=0; i<nJetTrig; i++)

  //----------------------- Computing The Prescale Values For a given event ----------------------- //
    for (int j=0; j<nJetTrig; j++)
         {
         if (ihltj[j] == -1)
         hltPassj[j] = true; // no trigger set
         else
		{
                if (Event->fired(ihltj[j]) > 0)
			{
                 	hltPassj[j] = true;
		 	triggered[j] = triggered[j] + 1;
                 	hltPass = true;
                 	prescalej[j] = Event->preL1(ihltj[j]) * Event->preHLT(ihltj[j]);
			probs[j] = 1.0;
			if (prescalej[j] > 0.0) { probs[j] = (1.0 - 1.0/prescalej[j]); }
			//if (test) { cout << "Event: " << i << " trigger : " << j << " prescale = " << prescalej[j] << endl; }
			}
		}
      	}

	prescale_old = 1.0 / ( 1.0 - (probs[0] * probs[1] * probs[2]));

     if (hltPass)
	{
	counter_hlt++;
	if (prescale_old > 1000.0)
		{
		cout << "Event: " << i << " Prescale = " << prescale_old << endl;
		cout << "Prescales = " << prescalej[0] << " " << prescalej[1] << " " << prescalej[2] << endl;
		}
	triggered[4] = triggered[4] + 1;

     	leading_pt = -10.0;
     	forward_pt = -10.0;
     	central_pt = -10.0;
     	gap_leading_pt = -10.0;
	outside_leading_pt = -10.0;
     	forward_eta = -10.0;
     	central_eta = -10.0;
     	leading_eta = -10.0;
     	gap_eta = -10.0;
	gap_eta_star = -100.0;
	outside_eta = -10.0;
     	forward_phi = -10.0;
     	central_phi = -10.0;
     	leading_phi = -10.0;
     	gap_phi = -10.0;
	outside_phi = -10.0;
	hard_emission = false;
	deta_out1 = -10.0;
	deta_out2 = -10.0;
	delta_eta = -1.0;
	delta_phi = -1.0;
	outside_delta_eta = -1.0;
	index = -1;
	index1 = -1;
	index2 = -1;
	prescale = 1.0;
	prescale_new = 1.0;
	emulation_pass = false;

    for(unsigned int l=0; l<Event->nPFJets(); l++)
	{
	pt = Event->pfjet(l).ptCor();
	eta = Event->pfjet(l).eta();
	phi = Event->pfjet(l).phi();
     	if (pt >= pt_min && Event->pfjet(l).tightID() )
		{
     		counter_jet++;
     		if (leading_pt < pt) { leading_pt = pt; leading_eta = eta; leading_phi = phi; }
     		if (eta <= 2.8 && eta >= -2.8 && pt > central_pt) { central_pt = pt; central_eta = eta; central_phi = phi; }
     		if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > forward_pt)
     			{
			forward_pt = pt; forward_eta = eta; forward_phi = phi;
			}
     		}
     	}

	if (leading_pt > 35.0 and leading_pt < 70.0 and hltPassj[0])
		{
		hltPassj[5] = true;
		triggered[5] = triggered[5] + 1;
		prescale = prescalej[0];
		}
	if (leading_pt > 70.0 and leading_pt < 110.0 and hltPassj[1])
		{
		hltPassj[5] = true;
		triggered[5] = triggered[5] + 1;
		prescale = prescalej[1];
		}

	if (leading_pt > 110.0 and hltPassj[2])
		{
		hltPassj[5] = true;
		triggered[5] = triggered[5] + 1;
		prescale = prescalej[2];
		}

	if (leading_pt > 35.0 and leading_pt < 70.0 and hltPassj[4])
		{
		emulation_pass = true;
		prescale_new = prescalej[4];
		}
	if (leading_pt > 70.0 and leading_pt < 110.0 and hltPassj[0])
		{
		emulation_pass = true;
		prescale_new = prescalej[1];
		}

	if (leading_pt > 110.0 and hltPassj[1])
		{
		emulation_pass = true;
		prescale_new = prescalej[1];
		}

	if (leading_pt > 35.0 and leading_pt < 70.0 and hltPassj[0] and !hltPassj[1] and !hltPassj[2])
		{
		hltPassj[6] = true;
		triggered[6] = triggered[6] + 1;
		prescale = prescalej[0];
		}
	if (leading_pt > 70.0 and leading_pt < 110.0 and !hltPassj[0] and hltPassj[1] and !hltPassj[2])
		{
		hltPassj[6] = true;
		triggered[6] = triggered[6] + 1;
		prescale = prescalej[1];
		}

	if (leading_pt > 110.0 and !hltPassj[0] and !hltPassj[1] and hltPassj[2])
		{
		hltPassj[6] = true;
		triggered[6] = triggered[6] + 1;
		prescale = prescalej[2];
		}
     
          if (forward_pt > pt_min && central_pt > pt_min)
     		{
     		counter_selected++;
		delta_eta = abs(forward_eta - central_eta);
		delta_phi = calc_delta_phi(central_phi, forward_phi);
    		for(unsigned int n=0; n<Event->nPFJets(); n++)
			{
			pt = Event->pfjet(n).ptCor();
			eta = Event->pfjet(n).eta();
			phi = Event->pfjet(n).phi();
     			if (Event->pfjet(n).tightID() )
				{
     				if (central_eta > forward_eta && pt > gap_leading_pt && eta > forward_eta && eta < central_eta )
					{
        				gap_leading_pt = pt; gap_eta = eta; gap_phi = phi;
					}
     				if (central_eta < forward_eta && pt > gap_leading_pt && eta < forward_eta && eta > central_eta )
					{
       					gap_leading_pt = pt; gap_eta = eta; gap_phi = phi;
					}
     				if (pt > outside_leading_pt && central_eta > forward_eta && (eta < forward_eta || eta > central_eta))
					{
        				outside_leading_pt = pt; outside_eta = eta; outside_phi = phi;
        				}
     				if (pt > outside_leading_pt && central_eta < forward_eta && (eta > forward_eta || eta < central_eta)) 
					{
        				outside_leading_pt = pt; outside_eta = eta; outside_phi = phi;
					}
        			}
     			}

			gap_eta_star = gap_eta - (forward_eta + central_eta)/2;
     			deta_out1 = central_eta - outside_eta;
     			if (deta_out1 < 0) { deta_out1 = -deta_out1; }
     			deta_out2 = forward_eta - outside_eta;
     			if (deta_out2 < 0) { deta_out2 = -deta_out2; }
     			if (deta_out1 < deta_out2) { outside_delta_eta = deta_out1; }
     			if (deta_out2 < deta_out1) { outside_delta_eta = deta_out2; }

			if (gap_leading_pt > pt_min_gap) { index1 = 1; }
			else { index1 = 2; }
			if (delta_eta >= 0.4 and delta_eta < 2.5) { index2 = 0; }
			if (delta_eta >= 2.5 and delta_eta < 3.5) { index2 = 1; }
			if (delta_eta >= 3.5 and delta_eta < 4.5) { index2 = 2; }
			if (delta_eta >= 4.4 and delta_eta < 7.5) { index2 = 3; }

			index = index1 * 5 + index2;

    		for (int l=0; l<nJetTrig; l++)
         		{
			if (hltPassj[l])
				{
				selected[l] = selected[l] + 1;
				events_selected[l][index] = events_selected[l][index] + 1;
				events_selected[l][index2] = events_selected[l][index2] + 1;
				events_selected[l][index1 * 5 + 4] = events_selected[l][index1 * 5 + 4] + 1;
				events_selected[l][4] = events_selected[l][4] + 1;
				}
			}
		selected[4] = selected[4] + 1;
		events_selected[4][index] = events_selected[4][index] + 1;
		events_selected[4][index2] = events_selected[4][index2] + 1;
		events_selected[4][index1 * 5 + 4] = events_selected[4][index1 * 5 + 4] + 1;
		events_selected[4][4] = events_selected[4][4] + 1;
		if (emulation_pass)
			{
			hist_delta_phi->Fill(delta_phi, prescale_new);
     			if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1->Fill(delta_phi, prescale_new); }
     			if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2->Fill(delta_phi, prescale_new); }
     			if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3->Fill(delta_phi, prescale_new); }
     			if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4->Fill(delta_phi, prescale_new); }
     			if (gap_leading_pt < pt_min_gap)
				{
				hist_delta_phi_gap->Fill(delta_phi,prescale_new);
     				if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_gap->Fill(delta_phi, prescale_new); }
     				if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_gap->Fill(delta_phi, prescale_new); }
     				if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_gap->Fill(delta_phi, prescale_new); }
     				if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_gap->Fill(delta_phi, prescale_new); }
				}
     			if (gap_leading_pt > pt_min_gap)
				{
				hist_delta_phi_nogap->Fill(delta_phi, prescale_new);
     				if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_nogap->Fill(delta_phi, prescale_new); }
     				if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_nogap->Fill(delta_phi, prescale_new); }
     				if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_nogap->Fill(delta_phi, prescale_new); }
     				if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_nogap->Fill(delta_phi, prescale_new); }
				hist_leading_pt_inside_gap->Fill(gap_leading_pt, prescale_new);
				hist_leading_eta_star_inside_gap->Fill(gap_eta_star, prescale_new);
				}
			if (outside_leading_pt > pt_min_gap ) { hist_delta_eta_outside_gap->Fill(outside_delta_eta, prescale_new); }
			}

		if (hltPassj[5])
			{
			selected[5] = selected[5] + 1;
			events_selected[5][index] = events_selected[5][index] + 1;
			events_selected[5][index2] = events_selected[5][index2] + 1;
			events_selected[5][index1 * 5 + 4] = events_selected[5][index1 * 5 + 4] + 1;
			events_selected[5][4] = events_selected[5][4] + 1;
			hist_delta_phi_inclusive->Fill(delta_phi, prescale);
     			if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_inclusive->Fill(delta_phi, prescale); }
     			if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_inclusive->Fill(delta_phi, prescale); }
     			if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_inclusive->Fill(delta_phi, prescale); }
     			if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_inclusive->Fill(delta_phi, prescale); }
     			if (gap_leading_pt < pt_min_gap)
				{
				hist_delta_phi_gap_inclusive->Fill(delta_phi, prescale);
     				if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_gap_inclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_gap_inclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_gap_inclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_gap_inclusive->Fill(delta_phi, prescale); }
				}
     			if (gap_leading_pt > pt_min_gap)
				{
				hist_delta_phi_nogap_inclusive->Fill(delta_phi, prescale);
     				if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_nogap_inclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_nogap_inclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_nogap_inclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_nogap_inclusive->Fill(delta_phi, prescale); }
				hist_leading_pt_inside_gap_inclusive->Fill(gap_leading_pt, prescale);
				hist_leading_eta_star_inside_gap_inclusive->Fill(gap_eta_star, prescale);
				}
			if (outside_leading_pt > pt_min_gap ) { hist_delta_eta_outside_gap_inclusive->Fill(outside_delta_eta, prescale); }
			}

		if (hltPassj[6])
			{
			selected[6] = selected[6] + 1;
			events_selected[6][index] = events_selected[6][index] + 1;
			events_selected[6][index2] = events_selected[6][index2] + 1;
			events_selected[6][index1 * 5 + 4] = events_selected[6][index1 * 5 + 4] + 1;
			events_selected[6][4] = events_selected[6][4] + 1;
			hist_delta_phi_exclusive->Fill(delta_phi, prescale);
     			if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_exclusive->Fill(delta_phi, prescale); }
     			if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_exclusive->Fill(delta_phi, prescale); }
     			if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_exclusive->Fill(delta_phi, prescale); }
     			if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_exclusive->Fill(delta_phi, prescale); }
     			if (gap_leading_pt < pt_min_gap)
				{
				hist_delta_phi_gap_exclusive->Fill(delta_phi, prescale);
     				if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_gap_exclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_gap_exclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_gap_exclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_gap_exclusive->Fill(delta_phi, prescale); }
				}
     			if (gap_leading_pt > pt_min_gap)
				{
				hist_delta_phi_nogap_exclusive->Fill(delta_phi, prescale);
     				if (delta_eta >= 0.4 and delta_eta < 2.5) { hist_delta_phi_deta1_nogap_exclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 2.5 and delta_eta < 3.5) { hist_delta_phi_deta2_nogap_exclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 3.5 and delta_eta < 4.5) { hist_delta_phi_deta3_nogap_exclusive->Fill(delta_phi, prescale); }
     				if (delta_eta >= 4.4 and delta_eta < 7.5) { hist_delta_phi_deta4_nogap_exclusive->Fill(delta_phi, prescale); }
				hist_leading_pt_inside_gap_exclusive->Fill(gap_leading_pt, prescale);
				hist_leading_eta_star_inside_gap_exclusive->Fill(gap_eta_star, prescale);
				}
			if (outside_leading_pt > pt_min_gap ) { hist_delta_eta_outside_gap_exclusive->Fill(outside_delta_eta, prescale); }
			}

		}
	}

}

} 
}

//computing auxiliary outputs
     pu_scale = (double)counter_entries/(double)counter_pv;
     eff = (double)counter_hlt/(double)counter_pv;

//output a summary
     if (detail) { cout<<"Total Number of Events :         "<<counter_entries<<endl; }
     if (detail) { cout<<"Primary Vertex Filter :          "<<counter_pv<<endl; }
     if (detail) { cout<<"Pileup Scale :                   "<<pu_scale<<endl; }
     if (detail) { cout<<"Triggered Events :               "<<counter_hlt<<endl; }
     if (detail) { cout<<"Trigger Efficiency :             "<<eff<<endl; }
     if (detail) { cout<<"Selected :                       "<<counter_selected<<endl; }
     if (detail) { cout<<"Total Jets :                     "<<counter_jet<<endl; }
     
     //fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_hlt);
     hist_events->SetBinContent(4,counter_selected);
     hist_events->SetBinContent(7,counter_jet);
     hist_events->SetBinContent(8,pu_scale);
     hist_events->SetBinContent(9,eff);

     //output a summary for HLT_Jet15U
     if (detail) { cout<<"Results for HLT_Jet15U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet15U: "<<triggered[0]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[0]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[0][0]<< "        "<<setw(7)<<events_selected[0][5]<<"       "<<setw(7)<<events_selected[0][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[0][1]<< "        "<<setw(7)<<events_selected[0][6]<<"       "<<setw(7)<<events_selected[0][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[0][2]<< "        "<<setw(7)<<events_selected[0][7]<<"       "<<setw(7)<<events_selected[0][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[0][3]<< "        "<<setw(7)<<events_selected[0][8]<<"       "<<setw(7)<<events_selected[0][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[0][4]<< "        "<<setw(7)<<events_selected[0][9]<<"       "<<setw(7)<<events_selected[0][14]<<endl; }

     //fill the events in HLT_Jet15U histogram
     hist_events_HLT_Jet15U->SetBinContent(1,triggered[0]);
     hist_events_HLT_Jet15U->SetBinContent(2,selected[0]);

     //output a summary for HLT_Jet30U
     if (detail) { cout<<"Results for HLT_Jet30U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet30U:  "<<triggered[1]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[1]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[1][0]<< "        "<<setw(7)<<events_selected[1][5]<<"       "<<setw(7)<<events_selected[1][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[1][1]<< "        "<<setw(7)<<events_selected[1][6]<<"       "<<setw(7)<<events_selected[1][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[1][2]<< "        "<<setw(7)<<events_selected[1][7]<<"       "<<setw(7)<<events_selected[1][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[1][3]<< "        "<<setw(7)<<events_selected[1][8]<<"       "<<setw(7)<<events_selected[1][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[1][4]<< "        "<<setw(7)<<events_selected[1][9]<<"       "<<setw(7)<<events_selected[1][14]<<endl; }

     //fill the events in HLT_Jet30U histogram
     hist_events_HLT_Jet30U->SetBinContent(1,triggered[1]);
     hist_events_HLT_Jet30U->SetBinContent(2,selected[1]);

     //output a summary for HLT_Jet50U
     if (detail) { cout<<"Results for HLT_Jet50U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet50U:  "<<triggered[2]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[2]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[2][0]<< "        "<<setw(7)<<events_selected[2][5]<<"       "<<setw(7)<<events_selected[2][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[2][1]<< "        "<<setw(7)<<events_selected[2][6]<<"       "<<setw(7)<<events_selected[2][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[2][2]<< "        "<<setw(7)<<events_selected[2][7]<<"       "<<setw(7)<<events_selected[2][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[2][3]<< "        "<<setw(7)<<events_selected[2][8]<<"       "<<setw(7)<<events_selected[2][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[2][4]<< "        "<<setw(7)<<events_selected[2][9]<<"       "<<setw(7)<<events_selected[2][14]<<endl; }

     //fill the events in HLT_Jet50U histogram
     hist_events_HLT_Jet50U->SetBinContent(1,triggered[2]);
     hist_events_HLT_Jet50U->SetBinContent(2,selected[2]);

     //output a summary for combination
     if (detail) { cout<<"Results for combination" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[3]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[3]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[3][0]<< "        "<<setw(7)<<events_selected[3][5]<<"       "<<setw(7)<<events_selected[3][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[3][1]<< "        "<<setw(7)<<events_selected[3][6]<<"       "<<setw(7)<<events_selected[3][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[3][2]<< "        "<<setw(7)<<events_selected[3][7]<<"       "<<setw(7)<<events_selected[3][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[3][3]<< "        "<<setw(7)<<events_selected[3][8]<<"       "<<setw(7)<<events_selected[3][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[3][4]<< "        "<<setw(7)<<events_selected[3][9]<<"       "<<setw(7)<<events_selected[3][14]<<endl; }

     //fill the events in combination histogram
     hist_events_all->SetBinContent(1,triggered[3]);
     hist_events_all->SetBinContent(2,selected[3]);

     trigger_loss1 = (double) selected[4] / (double) selected[3];
     trigger_loss2 = (double) selected[5] / (double) selected[3];

for (int m = 0; m < 15; m++)
	{
		events_loss[0][m] = (double) events_selected[4][m] / (double) events_selected[3][m];
		events_loss[1][m] = (double) events_selected[5][m] / (double) events_selected[3][m];
	}

     //output a summary for combination
     if (detail) { cout<<"Results for combination with inclusive zones" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[4]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[4]<<endl; }
     if (detail) { cout<<"Proportion :                           "<<trigger_loss1<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[4][0]<< "        "<<setw(7)<<events_selected[4][5]<<"       "<<setw(7)<<events_selected[4][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[4][1]<< "        "<<setw(7)<<events_selected[4][6]<<"       "<<setw(7)<<events_selected[4][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[4][2]<< "        "<<setw(7)<<events_selected[4][7]<<"       "<<setw(7)<<events_selected[4][11]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[4][3]<< "        "<<setw(7)<<events_selected[4][8]<<"       "<<setw(7)<<events_selected[4][12]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[4][4]<< "        "<<setw(7)<<events_selected[4][9]<<"       "<<setw(7)<<events_selected[4][14]<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_loss[0][0]<< "        "<<setw(7)<<events_loss[0][5]<<"       "<<setw(7)<<events_loss[0][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_loss[0][1]<< "        "<<setw(7)<<events_loss[0][6]<<"       "<<setw(7)<<events_loss[0][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_loss[0][2]<< "        "<<setw(7)<<events_loss[0][7]<<"       "<<setw(7)<<events_loss[0][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_loss[0][3]<< "        "<<setw(7)<<events_loss[0][8]<<"       "<<setw(7)<<events_loss[0][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_loss[0][4]<< "        "<<setw(7)<<events_loss[0][9]<<"       "<<setw(7)<<events_loss[0][14]<<endl; }

     //fill the events in combination histogram
     hist_events_inclusive_zones->SetBinContent(1,triggered[4]);
     hist_events_inclusive_zones->SetBinContent(2,selected[4]);
     hist_events_inclusive_zones->SetBinContent(3,trigger_loss1);

     //output a summary for combination
     if (detail) { cout<<"Results for combination with exclusive zones" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[5]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[5]<<endl; }
     if (detail) { cout<<"Proportion :                           "<<trigger_loss2<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[5][0]<< "        "<<setw(7)<<events_selected[5][5]<<"       "<<setw(7)<<events_selected[5][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[5][1]<< "        "<<setw(7)<<events_selected[5][6]<<"       "<<setw(7)<<events_selected[5][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[5][2]<< "        "<<setw(7)<<events_selected[5][7]<<"       "<<setw(7)<<events_selected[5][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[5][3]<< "        "<<setw(7)<<events_selected[5][8]<<"       "<<setw(7)<<events_selected[5][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[5][4]<< "        "<<setw(7)<<events_selected[5][9]<<"       "<<setw(7)<<events_selected[5][14]<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_loss[1][0]<< "        "<<setw(7)<<events_loss[1][5]<<"       "<<setw(7)<<events_loss[1][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_loss[1][1]<< "        "<<setw(7)<<events_loss[1][6]<<"       "<<setw(7)<<events_loss[1][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_loss[1][2]<< "        "<<setw(7)<<events_loss[1][7]<<"       "<<setw(7)<<events_loss[1][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_loss[1][3]<< "        "<<setw(7)<<events_loss[1][8]<<"       "<<setw(7)<<events_loss[1][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_loss[1][4]<< "        "<<setw(7)<<events_loss[1][9]<<"       "<<setw(7)<<events_loss[1][14]<<endl; }


     //fill the events in combination histogram
     hist_events_exclusive_zones->SetBinContent(1,triggered[5]);
     hist_events_exclusive_zones->SetBinContent(2,selected[5]);
     hist_events_exclusive_zones->SetBinContent(3,trigger_loss2);


	estimate_trigger_systematic(hist_delta_phi_inclusive, hist_delta_phi, syst_delta_phi_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta1_inclusive, hist_delta_phi_deta1, syst_delta_phi_deta1_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta2_inclusive, hist_delta_phi_deta2, syst_delta_phi_deta2_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta3_inclusive, hist_delta_phi_deta3, syst_delta_phi_deta3_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta4_inclusive, hist_delta_phi_deta4, syst_delta_phi_deta4_inclusive);
	estimate_trigger_systematic(hist_delta_phi_gap_inclusive, hist_delta_phi_gap, syst_delta_phi_gap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta1_gap_inclusive, hist_delta_phi_deta1_gap, syst_delta_phi_deta1_gap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta2_gap_inclusive, hist_delta_phi_deta2_gap, syst_delta_phi_deta2_gap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta3_gap_inclusive, hist_delta_phi_deta3_gap, syst_delta_phi_deta3_gap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta4_gap_inclusive, hist_delta_phi_deta4_gap, syst_delta_phi_deta4_gap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_nogap_inclusive, hist_delta_phi_nogap, syst_delta_phi_nogap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta1_nogap_inclusive, hist_delta_phi_deta1_nogap, syst_delta_phi_deta1_nogap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta2_nogap_inclusive, hist_delta_phi_deta2_nogap, syst_delta_phi_deta2_nogap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta3_nogap_inclusive, hist_delta_phi_deta3_nogap, syst_delta_phi_deta3_nogap_inclusive);
	estimate_trigger_systematic(hist_delta_phi_deta4_nogap_inclusive, hist_delta_phi_deta4_nogap, syst_delta_phi_deta4_nogap_inclusive);
	estimate_trigger_systematic(hist_leading_pt_inside_gap_inclusive, hist_leading_pt_inside_gap, syst_leading_pt_inside_gap_inclusive);
	estimate_trigger_systematic(hist_leading_eta_star_inside_gap_inclusive, hist_leading_eta_star_inside_gap, syst_leading_eta_star_inside_gap_inclusive);
	estimate_trigger_systematic(hist_delta_eta_outside_gap_inclusive, hist_delta_eta_outside_gap, syst_delta_eta_outside_gap_inclusive);


	estimate_trigger_systematic(hist_delta_phi_exclusive, hist_delta_phi, syst_delta_phi_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta1_exclusive, hist_delta_phi_deta1, syst_delta_phi_deta1_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta2_exclusive, hist_delta_phi_deta2, syst_delta_phi_deta2_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta3_exclusive, hist_delta_phi_deta3, syst_delta_phi_deta3_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta4_exclusive, hist_delta_phi_deta4, syst_delta_phi_deta4_exclusive);
	estimate_trigger_systematic(hist_delta_phi_gap_exclusive, hist_delta_phi_gap, syst_delta_phi_gap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta1_gap_exclusive, hist_delta_phi_deta1_gap, syst_delta_phi_deta1_gap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta2_gap_exclusive, hist_delta_phi_deta2_gap, syst_delta_phi_deta2_gap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta3_gap_exclusive, hist_delta_phi_deta3_gap, syst_delta_phi_deta3_gap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta4_gap_exclusive, hist_delta_phi_deta4_gap, syst_delta_phi_deta4_gap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_nogap_exclusive, hist_delta_phi_nogap, syst_delta_phi_nogap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta1_nogap_exclusive, hist_delta_phi_deta1_nogap, syst_delta_phi_deta1_nogap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta2_nogap_exclusive, hist_delta_phi_deta2_nogap, syst_delta_phi_deta2_nogap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta3_nogap_exclusive, hist_delta_phi_deta3_nogap, syst_delta_phi_deta3_nogap_exclusive);
	estimate_trigger_systematic(hist_delta_phi_deta4_nogap_exclusive, hist_delta_phi_deta4_nogap, syst_delta_phi_deta4_nogap_exclusive);
	estimate_trigger_systematic(hist_leading_pt_inside_gap_exclusive, hist_leading_pt_inside_gap, syst_leading_pt_inside_gap_exclusive);
	estimate_trigger_systematic(hist_leading_eta_star_inside_gap_exclusive, hist_leading_eta_star_inside_gap, syst_leading_eta_star_inside_gap_exclusive);
	estimate_trigger_systematic(hist_delta_eta_outside_gap_exclusive, hist_delta_eta_outside_gap, syst_delta_eta_outside_gap_exclusive);


//Open the output root file
     if (detail) { cout<<"Opening "<<data_out<<endl; }
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

//write histograms on file
     if (detail) { cout<<"Writing histograms on file..."<<endl; }

//general control histograms
     hist_events->Write();
     hist_events_HLT_Jet15U->Write();
     hist_events_HLT_Jet30U->Write();
     hist_events_HLT_Jet50U->Write();
     hist_events_all->Write();
     hist_events_inclusive_zones->Write();
     hist_events_exclusive_zones->Write();

     hist_delta_phi->Write();
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
     hist_leading_pt_inside_gap->Write();
     hist_leading_eta_star_inside_gap->Write();
     hist_delta_eta_outside_gap->Write();

     hist_delta_phi_inclusive->Write();
     hist_delta_phi_gap_inclusive->Write();
     hist_delta_phi_nogap_inclusive->Write();
     hist_delta_phi_deta1_inclusive->Write();
     hist_delta_phi_deta2_inclusive->Write();
     hist_delta_phi_deta3_inclusive->Write();
     hist_delta_phi_deta4_inclusive->Write();
     hist_delta_phi_deta1_gap_inclusive->Write();
     hist_delta_phi_deta2_gap_inclusive->Write();
     hist_delta_phi_deta3_gap_inclusive->Write();
     hist_delta_phi_deta4_gap_inclusive->Write();
     hist_delta_phi_deta1_nogap_inclusive->Write();
     hist_delta_phi_deta2_nogap_inclusive->Write();
     hist_delta_phi_deta3_nogap_inclusive->Write();
     hist_delta_phi_deta4_nogap_inclusive->Write();
     hist_leading_pt_inside_gap_inclusive->Write();
     hist_leading_eta_star_inside_gap_inclusive->Write();
     hist_delta_eta_outside_gap_inclusive->Write();

     hist_delta_phi_exclusive->Write();
     hist_delta_phi_gap_exclusive->Write();
     hist_delta_phi_nogap_exclusive->Write();
     hist_delta_phi_deta1_exclusive->Write();
     hist_delta_phi_deta2_exclusive->Write();
     hist_delta_phi_deta3_exclusive->Write();
     hist_delta_phi_deta4_exclusive->Write();
     hist_delta_phi_deta1_gap_exclusive->Write();
     hist_delta_phi_deta2_gap_exclusive->Write();
     hist_delta_phi_deta3_gap_exclusive->Write();
     hist_delta_phi_deta4_gap_exclusive->Write();
     hist_delta_phi_deta1_nogap_exclusive->Write();
     hist_delta_phi_deta2_nogap_exclusive->Write();
     hist_delta_phi_deta3_nogap_exclusive->Write();
     hist_delta_phi_deta4_nogap_exclusive->Write();
     hist_leading_pt_inside_gap_exclusive->Write();
     hist_leading_eta_star_inside_gap_exclusive->Write();
     hist_delta_eta_outside_gap_exclusive->Write();

     syst_delta_phi_inclusive->Write();
     syst_delta_phi_gap_inclusive->Write();
     syst_delta_phi_nogap_inclusive->Write();
     syst_delta_phi_deta1_inclusive->Write();
     syst_delta_phi_deta2_inclusive->Write();
     syst_delta_phi_deta3_inclusive->Write();
     syst_delta_phi_deta4_inclusive->Write();
     syst_delta_phi_deta1_gap_inclusive->Write();
     syst_delta_phi_deta2_gap_inclusive->Write();
     syst_delta_phi_deta3_gap_inclusive->Write();
     syst_delta_phi_deta4_gap_inclusive->Write();
     syst_delta_phi_deta1_nogap_inclusive->Write();
     syst_delta_phi_deta2_nogap_inclusive->Write();
     syst_delta_phi_deta3_nogap_inclusive->Write();
     syst_delta_phi_deta4_nogap_inclusive->Write();
     syst_leading_pt_inside_gap_inclusive->Write();
     syst_leading_eta_star_inside_gap_inclusive->Write();
     syst_delta_eta_outside_gap_inclusive->Write();

     syst_delta_phi_exclusive->Write();
     syst_delta_phi_gap_exclusive->Write();
     syst_delta_phi_nogap_exclusive->Write();
     syst_delta_phi_deta1_exclusive->Write();
     syst_delta_phi_deta2_exclusive->Write();
     syst_delta_phi_deta3_exclusive->Write();
     syst_delta_phi_deta4_exclusive->Write();
     syst_delta_phi_deta1_gap_exclusive->Write();
     syst_delta_phi_deta2_gap_exclusive->Write();
     syst_delta_phi_deta3_gap_exclusive->Write();
     syst_delta_phi_deta4_gap_exclusive->Write();
     syst_delta_phi_deta1_nogap_exclusive->Write();
     syst_delta_phi_deta2_nogap_exclusive->Write();
     syst_delta_phi_deta3_nogap_exclusive->Write();
     syst_delta_phi_deta4_nogap_exclusive->Write();
     syst_leading_pt_inside_gap_exclusive->Write();
     syst_leading_eta_star_inside_gap_exclusive->Write();
     syst_delta_eta_outside_gap_exclusive->Write();

     if (detail) { cout<<"Histograms written sucessfully!"<<endl; }
     
     //close the output file
     data_output->Close();
     
//delete the histograms to avoid memory leak
//
     delete(hist_events);
     delete(hist_events_HLT_Jet15U);
     delete(hist_events_HLT_Jet30U);
     delete(hist_events_HLT_Jet50U);
     delete(hist_events_all);
     delete(hist_events_inclusive_zones);
     delete(hist_events_exclusive_zones);


     delete(hist_delta_phi);
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
     delete(hist_leading_pt_inside_gap);
     delete(hist_leading_eta_star_inside_gap);
     delete(hist_delta_eta_outside_gap);

     delete(hist_delta_phi_inclusive);
     delete(hist_delta_phi_gap_inclusive);
     delete(hist_delta_phi_nogap_inclusive);
     delete(hist_delta_phi_deta1_inclusive);
     delete(hist_delta_phi_deta2_inclusive);
     delete(hist_delta_phi_deta3_inclusive);
     delete(hist_delta_phi_deta4_inclusive);
     delete(hist_delta_phi_deta1_gap_inclusive);
     delete(hist_delta_phi_deta2_gap_inclusive);
     delete(hist_delta_phi_deta3_gap_inclusive);
     delete(hist_delta_phi_deta4_gap_inclusive);
     delete(hist_delta_phi_deta1_nogap_inclusive);
     delete(hist_delta_phi_deta2_nogap_inclusive);
     delete(hist_delta_phi_deta3_nogap_inclusive);
     delete(hist_delta_phi_deta4_nogap_inclusive);
     delete(hist_leading_pt_inside_gap_inclusive);
     delete(hist_leading_eta_star_inside_gap_inclusive);
     delete(hist_delta_eta_outside_gap_inclusive);

     delete(hist_delta_phi_exclusive);
     delete(hist_delta_phi_gap_exclusive);
     delete(hist_delta_phi_nogap_exclusive);
     delete(hist_delta_phi_deta1_exclusive);
     delete(hist_delta_phi_deta2_exclusive);
     delete(hist_delta_phi_deta3_exclusive);
     delete(hist_delta_phi_deta4_exclusive);
     delete(hist_delta_phi_deta1_gap_exclusive);
     delete(hist_delta_phi_deta2_gap_exclusive);
     delete(hist_delta_phi_deta3_gap_exclusive);
     delete(hist_delta_phi_deta4_gap_exclusive);
     delete(hist_delta_phi_deta1_nogap_exclusive);
     delete(hist_delta_phi_deta2_nogap_exclusive);
     delete(hist_delta_phi_deta3_nogap_exclusive);
     delete(hist_delta_phi_deta4_nogap_exclusive);
     delete(hist_leading_pt_inside_gap_exclusive);
     delete(hist_leading_eta_star_inside_gap_exclusive);
     delete(hist_delta_eta_outside_gap_exclusive);

     delete(syst_delta_phi_inclusive);
     delete(syst_delta_phi_gap_inclusive);
     delete(syst_delta_phi_nogap_inclusive);
     delete(syst_delta_phi_deta1_inclusive);
     delete(syst_delta_phi_deta2_inclusive);
     delete(syst_delta_phi_deta3_inclusive);
     delete(syst_delta_phi_deta4_inclusive);
     delete(syst_delta_phi_deta1_gap_inclusive);
     delete(syst_delta_phi_deta2_gap_inclusive);
     delete(syst_delta_phi_deta3_gap_inclusive);
     delete(syst_delta_phi_deta4_gap_inclusive);
     delete(syst_delta_phi_deta1_nogap_inclusive);
     delete(syst_delta_phi_deta2_nogap_inclusive);
     delete(syst_delta_phi_deta3_nogap_inclusive);
     delete(syst_delta_phi_deta4_nogap_inclusive);
     delete(syst_leading_pt_inside_gap_inclusive);
     delete(syst_leading_eta_star_inside_gap_inclusive);
     delete(syst_delta_eta_outside_gap_inclusive);

     delete(syst_delta_phi_exclusive);
     delete(syst_delta_phi_gap_exclusive);
     delete(syst_delta_phi_nogap_exclusive);
     delete(syst_delta_phi_deta1_exclusive);
     delete(syst_delta_phi_deta2_exclusive);
     delete(syst_delta_phi_deta3_exclusive);
     delete(syst_delta_phi_deta4_exclusive);
     delete(syst_delta_phi_deta1_gap_exclusive);
     delete(syst_delta_phi_deta2_gap_exclusive);
     delete(syst_delta_phi_deta3_gap_exclusive);
     delete(syst_delta_phi_deta4_gap_exclusive);
     delete(syst_delta_phi_deta1_nogap_exclusive);
     delete(syst_delta_phi_deta2_nogap_exclusive);
     delete(syst_delta_phi_deta3_nogap_exclusive);
     delete(syst_delta_phi_deta4_nogap_exclusive);
     delete(syst_leading_pt_inside_gap_exclusive);
     delete(syst_leading_eta_star_inside_gap_exclusive);
     delete(syst_delta_eta_outside_gap_exclusive);

//Sucess confirmation
     if (detail) { cout<<"Done!"<<endl; }
}
