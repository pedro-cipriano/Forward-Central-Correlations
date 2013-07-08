// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 21 Nov 2012
//
// compute_vertex_weights(string input_mc, string label_mc, string input_data, string label_data, string output_root, string output_plots, string prefix, bool detail = false)
// computes the vertex weights


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
#include <iomanip>
using namespace std;

#include "common_methods.h"


void plot_ratio(TH1D *ratio, string path, string fileout, string legend_position = "top_left", bool log_scale = false, bool min_is_0 = false, bool detail = false)
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

    if (!log_scale and min_is_0) { ratio->SetMinimum(0.0); }
    //ratio->SetMaximum(1.0);
    ratio->SetLineWidth(4);
    ratio->SetLineColor(2);
    ratio->SetLineStyle(1);
    ratio->Draw("e1");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(ratio,"Vertex Ratio","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

//setting the output files
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
   
//save the file and close the canvas
    if (detail) { cout<<"Saving to "<<out_png<<endl; }
    c01->Print( out_png.c_str() );
    if (detail) { cout<<"Saving to "<<out_c<<endl; }
    c01->Print( out_c.c_str() );
    if (detail) { cout<<"Saving to "<<out_eps<<endl; }
    c01->Print( out_eps.c_str() );
    c01->Close();

}


void plot_test(TH1D *data, TString data_label, TH1D *mc, TString mc_label, TH1D *ratio, string path, string fileout, string legend_position = "top_left", bool detail = false)
{

    if (detail) { cout<<"Testing "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    for (int i = 2; i <= 10; i++)
    {
    double ratio_val = ratio->GetBinContent(i-1);
    double mc_val = mc->GetBinContent(i);
    //double ratio_error = ratio->GetBinError(i);
    double mc_error = mc->GetBinError(i);
    double new_val = ratio_val * mc_val;
    double new_error = ratio_val * mc_error;
    mc->SetBinContent(i,new_val);
    mc->SetBinError(i,new_error);
    //cout << i << " -> Ratio = " << ratio << endl; 
    }

    for (int i = 11; i <= mc->GetNbinsX(); i++)
    {
    mc->SetBinContent(i,0);
    mc->SetBinError(i,0);
    }

//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 0.0;
    double max = data->GetMaximum();
    if (data->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(data,detail);
    }
    else
    {
    min = data->GetMinimum();
    } 

    set_histogram_min_max(mc, min, max, detail);
    
    max = 2.0 * max;
    min = 0.7 * min;

    cout << "max = " << max << " min = " << min << endl;
//format and ploting the histogram
    if (detail) { cout << "Drawning on the canvas..." << endl; }
    data->SetMaximum(max);
    data->SetMinimum(min);

    data->SetLineWidth(4);
    data->SetLineColor(2);
    data->SetLineStyle(1);
    data->Draw("e1");
    mc->SetLineWidth(4);
    mc->SetLineColor(4);
    mc->SetLineStyle(2);
    mc->Draw("e1 same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 2, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(data,data_label,"l");
    leg01->AddEntry(mc,mc_label+" Reweighted","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    cout << "Integrals: data = " << data->Integral() << " mc = " << mc->Integral() << endl; 

    print_plots(c01, path, fileout);

}

void compute_vertex_weights(string input_mc, string mc_label, string input_data, string data_label, string original_mc, string output_root, TString output_sufix = "_v0", string input_root = "", TString input_sufix = "_v1", string output_plots = "../output/vertex_weights/", string prefix = "test_", bool keep_norm = false, bool detail = false)
{

//confingure the vertex range
   int v_begin = 2;
   int v_end = 10;

//output the configuration
   if (detail) { cout<<"Compute Vertex Weights configuration"<<endl; }
   if (detail) { cout<<"Input MC :           "<<input_mc<<endl; }
   if (detail) { cout<<"MC Label :           "<<mc_label<<endl; }
   if (detail) { cout<<"Input Data :         "<<input_data<<endl; }
   if (detail) { cout<<"Data Label :         "<<data_label<<endl; }
   if (detail) { cout<<"Original MC :        "<<original_mc<<endl; }
   if (detail) { cout<<"Output Root :        "<<output_root<<endl; }
   if (detail) { cout<<"Output Sufix :       "<<output_sufix<<endl; }
   if (detail) { cout<<"Input Root :         "<<input_root<<endl; }
   if (detail) { cout<<"Input Sufix :        "<<input_sufix<<endl; }
   if (detail) { cout<<"Output Plots :       "<<output_plots<<endl; }
   if (detail) { cout<<"Plot Prefix :        "<<prefix<<endl; }
   if (detail) { cout<<"Detail level :       "<<detail<<endl; }
   if (detail) { cout<<"Keep Normalization : "<<keep_norm<<endl; }
   if (detail) { cout<<"First Bin :          "<<v_begin<<endl; }
   if (detail) { cout<<"Last Bin :           "<<v_end<<endl; }

//opening the files
   if (detail) { cout << "Opening files..." << endl; }
   TFile *mc = new TFile( input_mc.c_str() );
   TFile *original = new TFile( original_mc.c_str() );
   TFile *data = new TFile( input_data.c_str() );
   bool vertex_weighted = false;
   TH1D *vertex_hist = 0;
   if (input_root != "")
   {
   if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
   TFile *vertex_file = new TFile( input_root.c_str() );
   vertex_file->GetObject("vertex_weights"+input_sufix,vertex_hist);
   if (vertex_hist != 0) { vertex_weighted = true; }
   }
   if (detail) { cout << "All files opened sucessfully!" << endl; }

//compute the vertex weights
    if (detail) { cout << "Computing the vertex weights..." << endl; }
    TH1D *mc_vertex = (TH1D*) mc->Get("ak5PF_vertex_selected");
    TH1D *original_vertex = (TH1D*) original->Get("ak5PF_vertex_selected");
    TH1D *data_vertex = (TH1D*) data->Get("ak5PF_vertex_selected");

    TH1D *vertex_weights;
    vertex_weights = new TH1D("vertex_weights"+output_sufix,"Vertex weights;Vertex;weights", 9, 1, 10);
    vertex_weights->Sumw2();
    
    TH1D *vertex_factor;
    vertex_factor = new TH1D("vertex_factor"+output_sufix,"Vertex factor;Vertex;factor", 9, 1, 10);
    vertex_factor->Sumw2();
    
    TH1D *weights_variation;
    weights_variation = new TH1D("weights_variation"+output_sufix,"Weights Variation;Vertex;variation", 9, 1, 10);
    weights_variation->Sumw2();
    
    if (keep_norm)
	{
	mc_vertex->Scale(1./mc_vertex->Integral());
	data_vertex->Scale(1./data_vertex->Integral());
	}
    if (detail) { cout << "mc_val = " << mc_vertex->Integral() << " data_val = " << data_vertex->Integral() << endl; }
    for (int i = v_begin; i <= v_end; i++)
    {
    double mc_val = mc_vertex->GetBinContent(i);
    double data_val = data_vertex->GetBinContent(i);
    double prev_value = 1.0;
    if (vertex_weighted) { prev_value = vertex_hist->GetBinContent(1 - v_begin + i); }
    double factor = data_val / mc_val;
    double ratio = prev_value * factor;
    //if (detail) { cout << "ratio = " << ratio  << " prve_value = "<< prev_value << " factor = " << factor << data_val << mc_val << endl; }
    vertex_weights->SetBinContent(i-1,ratio);
    vertex_factor->SetBinContent(i-1,factor);
    double variation = (ratio - prev_value) / ratio;
    if (detail) { cout << i << " -> Ratio = " << ratio << " previous value = " << prev_value << " change = " << 100*variation << "%" << endl; }
    weights_variation->SetBinContent(i-1,variation);
    }

    plot_2histograms(mc_vertex, mc_label, data_vertex, data_label, output_plots, prefix+"control_vertex_dist", "bottom_left", true, detail);
    plot_ratio(vertex_weights, output_plots, prefix+"vertex_weights_logscale", "bottom_left", true, false, detail);
    plot_ratio(vertex_weights, output_plots, prefix+"vertex_weights", "bottom_left", false, true, detail);
    plot_ratio(vertex_factor, output_plots, prefix+"vertex_factor_logscale", "bottom_left", true, false, detail);
    plot_ratio(vertex_factor, output_plots, prefix+"vertex_factor", "bottom_left", false, true, detail);
    plot_ratio(weights_variation, output_plots, prefix+"weights_variation", "bottom_left", false, false, detail);
    plot_test(data_vertex, data_label, mc_vertex, mc_label, vertex_weights, output_plots, prefix+"test_vertex_weights", "top_right", detail);

//creating the output file
    TFile data_output( output_root.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file "<<output_root<<" ..."<<endl; }
    vertex_weights->Write();
    vertex_factor->Write();
    weights_variation->Write();
    if (detail) { cout<<"Writing on "<<output_root<<" was sucessfull!"<<endl; }
    
//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    mc->Close();
    data->Close();
    data_output.Close();
    if (detail) { cout<<"Done!"<<endl; }
       

}
