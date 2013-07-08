// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 27 Nov 2012
//
// void compute_status_vertex_reweight(string v1, string v2, string v3, string v4, string v5, string path_plots = "../output/", string plot_prefix = "test_", bool detail = false)
// plots the diference between each iteration of the vertex reweight process


#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "common_methods.h"

void plot_status(TH1D *jetmettau, TH1D *jetmet, TH1D *jet, string path = "../output/", string fileout = "test", string legend_position = "top_left", bool log_scale = false, bool detail = false)
{

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    if (log_scale) { gPad->SetLogy(); }
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 1.0;
    double max = jetmettau->GetMaximum();
    if (jetmettau->GetMinimum() <= 0.0)
    {
    min = get_non0_minimum(jetmettau,detail);
    }
    else
    {
    min = jetmettau->GetMinimum();
    } 

    set_histogram_min_max(jetmet, min, max, detail);
    set_histogram_min_max(jet, min, max, detail);
    
    max = 1.5 * max;
    min = 0.7 * min;

    cout << "max = " << max << " min = " << min << endl;
//format and ploting the histogram
    if (detail) { cout << "Drawning on the canvas..." << endl; }
    jetmettau->SetMaximum(max);
    jetmettau->SetMinimum(min);
    if (!log_scale) { jetmettau->SetMinimum(0.0); }

    jetmettau->SetLineWidth(4);
    jetmettau->SetLineColor(2);
    jetmettau->SetLineStyle(1);
    jetmettau->Draw("e1");
    jetmet->SetMaximum(max);
    jetmet->SetMinimum(min);
    jetmet->SetLineWidth(4);
    jetmet->SetLineColor(1);
    jetmet->SetLineStyle(2);
    jetmet->Draw("e1 same");
    jet->SetMaximum(max);
    jet->SetMinimum(min);
    jet->SetLineWidth(4);
    jet->SetLineColor(4);
    jet->SetLineStyle(3);
    //jet->Draw("e1 same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 2, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(jetmettau,"JetMETTau_2010A","l");
    leg01->AddEntry(jetmet,"JetMET_2010A","l");
    //leg01->AddEntry(jet,"Jet_2010B","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}

void plot_status_vertex_reweight(string jetmettau, string jetmet, string jet, TString hist_sufix, string path_plots = "../output/", string plot_prefix = "test_", bool detail = false)
{

//output the configuration
   if (detail) { cout<<"Compute Status Vertex Reweighting configuration"<<endl; }
   if (detail) { cout<<"Input Dataset 1 :       "<<jetmettau<<endl; }
   if (detail) { cout<<"Input Dataset 2 :       "<<jetmet<<endl; }
   if (detail) { cout<<"Input Dataset 3 :       "<<jet<<endl; }
   if (detail) { cout<<"Input Histogram Sufix : "<<hist_sufix<<endl; }
   if (detail) { cout<<"Path Plots :            "<<path_plots<<endl; }
   if (detail) { cout<<"Plot Prefix :           "<<plot_prefix<<endl; }
   if (detail) { cout<<"Detail level :          "<<detail<<endl; }

//opening the files
   if (detail) { cout << "Opening files..." << endl; }
   TFile *wjetmettau = new TFile( jetmettau.c_str() );
   TFile *wjetmet = new TFile( jetmet.c_str() );
   TFile *wjet = new TFile( jet.c_str() );
   if (detail) { cout << "All files opened sucessfully!" << endl; }

//opening histograms
   if (detail) { cout << "Opening histograms..." << endl; }
   TH1D *weights_jetmettau = (TH1D*) wjetmettau->Get("vertex_weights"+hist_sufix);
   TH1D *weights_jetmet = (TH1D*) wjetmet->Get("vertex_weights"+hist_sufix);
   TH1D *weights_jet = (TH1D*) wjet->Get("vertex_weights"+hist_sufix);
   TH1D *factors_jetmettau = (TH1D*) wjetmettau->Get("vertex_factor"+hist_sufix);
   TH1D *factors_jetmet = (TH1D*) wjetmet->Get("vertex_factor"+hist_sufix);
   TH1D *factors_jet = (TH1D*) wjet->Get("vertex_factor"+hist_sufix);
   TH1D *variation_jetmettau = (TH1D*) wjetmettau->Get("weights_variation"+hist_sufix);
   TH1D *variation_jetmet = (TH1D*) wjetmet->Get("weights_variation"+hist_sufix);
   TH1D *variation_jet = (TH1D*) wjet->Get("weights_variation"+hist_sufix);
   if (detail) { cout << "Histograms opened sucessfully..." << endl; }

   if (detail) { cout << "Plotting Status..." << endl; }
   plot_status(weights_jetmettau, weights_jetmet, weights_jet, path_plots, plot_prefix+"weights_evolution", "bottom_left", true, detail);
   plot_status(factors_jetmettau, factors_jetmet, factors_jet,  path_plots, plot_prefix+"factors_evolution", "top_right", false, detail);
   plot_status(variation_jetmettau, variation_jetmet, variation_jet, path_plots, plot_prefix+"variation_evolution", "bottom_left", false, detail);
}
