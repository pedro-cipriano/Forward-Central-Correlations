// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 27 Nov 2012
//
// void compute_status_vertex_reweighting(string v1, string v2, string v3, string v4, string v5, string path_plots = "../output/", string plot_prefix = "test_", bool detail = false)
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

void plot_evolution(TH1D *v1, TH1D *v2, TH1D *v3, TH1D *v4, TH1D *v5, string path = "../output/", string fileout = "test", TString label = "test" , string legend_position = "top_left", bool log_scale = false, bool detail = false)
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
    double max = v1->GetMaximum();
    if (v1->GetMinimum() <= 0.0)
    {
    min = get_non0_minimum(v1,detail);
    }
    else
    {
    min = v1->GetMinimum();
    } 

    set_histogram_min_max(v2, min, max, detail);
    set_histogram_min_max(v3, min, max, detail);
    set_histogram_min_max(v4, min, max, detail);
    set_histogram_min_max(v5, min, max, detail);
    
    max = 1.5 * max;
    min = 0.7 * min;

    cout << "max = " << max << " min = " << min << endl;
//format and ploting the histogram
    if (detail) { cout << "Drawning on the canvas..." << endl; }
    v1->SetMaximum(max);
    v1->SetMinimum(min);
    if (!log_scale) { v1->SetMinimum(0.0); }

    v1->SetLineWidth(4);
    v1->SetLineColor(2);
    v1->SetLineStyle(1);
    v1->Draw("e1");
    v2->SetMaximum(max);
    v2->SetMinimum(min);
    v2->SetLineWidth(4);
    v2->SetLineColor(3);
    v2->SetLineStyle(2);
    v2->Draw("e1 same");
    v3->SetMaximum(max);
    v3->SetMinimum(min);
    v3->SetLineWidth(4);
    v3->SetLineColor(4);
    v3->SetLineStyle(3);
    v3->Draw("e1 same");
    v4->SetMaximum(max);
    v4->SetMinimum(min);
    v4->SetLineWidth(4);
    v4->SetLineColor(5);
    v4->SetLineStyle(4);
    v4->Draw("e1 same");
    v5->SetMaximum(max);
    v5->SetMinimum(min);
    v5->SetLineWidth(4);
    v5->SetLineColor(6);
    v5->SetLineStyle(5);
    v5->Draw("e1 same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 5, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(v1,label+" - first iteration","l");
    leg01->AddEntry(v2,label+" - second iteration","l");
    leg01->AddEntry(v3,label+" - third iteration","l");
    leg01->AddEntry(v4,label+" - forth iteration","l");
    leg01->AddEntry(v5,label+" - fifth iteration","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}

void plot_evolution_vertex_reweight(string v1, string v2, string v3, string v4, string v5, string path_plots = "../output/", string plot_prefix = "test_", bool detail = false)
{

//output the configuration
   if (detail) { cout<<"Compute Status Vertex Reweighting configuration"<<endl; }
   if (detail) { cout<<"Input Weights v1 : "<<v1<<endl; }
   if (detail) { cout<<"Input Weights v2 : "<<v2<<endl; }
   if (detail) { cout<<"Input Weights v3 : "<<v3<<endl; }
   if (detail) { cout<<"Input Weights v4 : "<<v4<<endl; }
   if (detail) { cout<<"Input Weights v5 : "<<v5<<endl; }
   if (detail) { cout<<"Path Plots :       "<<path_plots<<endl; }
   if (detail) { cout<<"Plot Prefix :      "<<plot_prefix<<endl; }
   if (detail) { cout<<"Detail level :     "<<detail<<endl; }

//opening the files
   if (detail) { cout << "Opening files..." << endl; }
   TFile *w_v1 = new TFile( v1.c_str() );
   TFile *w_v2 = new TFile( v2.c_str() );
   TFile *w_v3 = new TFile( v3.c_str() );
   TFile *w_v4 = new TFile( v4.c_str() );
   TFile *w_v5 = new TFile( v5.c_str() );
   if (detail) { cout << "All files opened sucessfully!" << endl; }

//opening histograms
   if (detail) { cout << "Opening histograms..." << endl; }
   TH1D *weights_v1 = (TH1D*) w_v1->Get("vertex_weights_v1");
   TH1D *weights_v2 = (TH1D*) w_v2->Get("vertex_weights_v2");
   TH1D *weights_v3 = (TH1D*) w_v3->Get("vertex_weights_v3");
   TH1D *weights_v4 = (TH1D*) w_v4->Get("vertex_weights_v4");
   TH1D *weights_v5 = (TH1D*) w_v5->Get("vertex_weights_v5");
   TH1D *factors_v1 = (TH1D*) w_v1->Get("vertex_factor_v1");
   TH1D *factors_v2 = (TH1D*) w_v2->Get("vertex_factor_v2");
   TH1D *factors_v3 = (TH1D*) w_v3->Get("vertex_factor_v3");
   TH1D *factors_v4 = (TH1D*) w_v4->Get("vertex_factor_v4");
   TH1D *factors_v5 = (TH1D*) w_v5->Get("vertex_factor_v5");
   TH1D *variation_v1 = (TH1D*) w_v1->Get("weights_variation_v1");
   TH1D *variation_v2 = (TH1D*) w_v2->Get("weights_variation_v2");
   TH1D *variation_v3 = (TH1D*) w_v3->Get("weights_variation_v3");
   TH1D *variation_v4 = (TH1D*) w_v4->Get("weights_variation_v4");
   TH1D *variation_v5 = (TH1D*) w_v5->Get("weights_variation_v5");
   if (detail) { cout << "Histograms opened sucessfully..." << endl; }

   if (detail) { cout << "Plotting Evolution..." << endl; }
   plot_evolution(weights_v1, weights_v2, weights_v3, weights_v4, weights_v5, path_plots, plot_prefix+"weights_evolution", "Weight", "bottom_left", true, detail);
   plot_evolution(factors_v1, factors_v2, factors_v3, factors_v4, factors_v5, path_plots, plot_prefix+"factors_evolution", "Factor", "top_right", false, detail);
   plot_evolution(variation_v1, variation_v2, variation_v3, variation_v4, variation_v5, path_plots, plot_prefix+"variation_evolution", "Variation", "bottom_left", false, detail);
}
