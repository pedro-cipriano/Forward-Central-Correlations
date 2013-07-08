// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 13 Dez 2012
//
//compute_correction(string input_file_uncorrected, string hist_prefix_unc, TString label_unc, string input_file_true,  string hist_prefix_true, TString label_true, string output_path_plots, string output_rootfile, string plot_prefix, TString legend_label, bool detail, bool show_stats)
//computes the correction factors for detector effects or pile up using the uncorrected and true distributions

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
#include <iomanip>
using namespace std;

#include "common_methods.h"

void new_correction(TH1D *hist_out, TH1D *hist1, TH1D *hist2, bool merge = false, bool detail = false)
{

if (merge)
{
//merges the correction factors
if (detail) { cout << "Merging correction factors..." << endl; }
hist_out->Multiply(hist1, hist2, 1.0, 1.0);
}
else
{
//computes the correction factor between uncorrected and corrected data
if (detail) { cout << "Computing correction factor..." << endl; }
hist_out->Divide(hist2, hist1, 1.0, 1.0);
}

}

void plot_correction(TH1D *hratio, string prefix, string fileout, string path= "../output/corrections/", TString legend_label = "Correction", string legend_position = "top_left", bool detail = false)
{
//plots the correction factor

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

//    hratio->SetMinimum(0.0);
//    hratio->SetMaximum(2.0);
    hratio->SetLineWidth(8);
    hratio->SetLineColor(4);
    hratio->SetLineStyle(2);
    hratio->Draw("e1");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hratio,legend_label,"l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, prefix + fileout);

}

void display_final_corrections(string input_file1, string prefix1, string label1, string input_file2, string prefix2, string label2, string input_file3, string prefix3, string label3, string output_path, string output_prefix, bool detail = false)
{
//computes the correction factors for detector effects or pile up using the uncorrected and true distributions

   if (detail) { cout<<"Input file1 :             "<<input_file1<<endl; }
   if (detail) { cout<<"Prefix for file1 :        "<<prefix1<<endl; }
   if (detail) { cout<<"Label for file1 :         "<<label1<<endl; }
   if (detail) { cout<<"Input file2 :             "<<input_file2<<endl; }
   if (detail) { cout<<"Prefix for file2 :        "<<prefix2<<endl; }
   if (detail) { cout<<"Label for file2 :         "<<label2<<endl; }
   if (detail) { cout<<"Input file3 :             "<<input_file3<<endl; }
   if (detail) { cout<<"Prefix for file3 :        "<<prefix3<<endl; }
   if (detail) { cout<<"Label for file3 :         "<<label3<<endl; }
   if (detail) { cout<<"Output Folder for plots : "<<output_path<<endl; }
   if (detail) { cout<<"Output Prefix :           "<<output_prefix<<endl; }
   if (detail) { cout<<"Detail :                  "<<detail<<endl; }

   bool have3files = true;
   if (input_file3 == "") { have3files = false; cout<<"Only two files have been given as input!"<<endl; }

//opening the input and output data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *file1 = new TFile( input_file1.c_str() );
    TFile *file2 = new TFile( input_file2.c_str() );
    TFile *file3 = 0;
    if (have3files) { file3 = new TFile( input_file3.c_str() ); }


//plot the correction factor for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_file1 = 0;
    TH1D *delta_phi_file2 = 0;
    TH1D *delta_phi_file3 = 0;
    TString delta_phi_file1_name = prefix1 + "delta_phi";
    TString delta_phi_file2_name = prefix2 + "delta_phi";
    TString delta_phi_file3_name = prefix3 + "delta_phi";

    file1->GetObject(delta_phi_file1_name,delta_phi_file1);
    if (delta_phi_file1 == 0) { cout << delta_phi_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_file2_name,delta_phi_file2);
    if (delta_phi_file2 == 0) { cout << delta_phi_file2_name << " not found!" << endl; return; }
    if (have3files)
	{ 
    	file3->GetObject(delta_phi_file3_name,delta_phi_file3);
    	if (delta_phi_file3 == 0) { cout << delta_phi_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_file1, label1, delta_phi_file2, label2, delta_phi_file3, label3, output_path, output_prefix + "delta_phi", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_file1, label1, delta_phi_file2, label2, output_path, output_prefix + "delta_phi", "top_right", false, detail); }


//plot the correction factor for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_file1 = 0;
    TH1D *delta_phi_deta1_file2 = 0;
    TH1D *delta_phi_deta1_file3 = 0;
    TString delta_phi_deta1_file1_name = prefix1 + "delta_phi_deta1";
    TString delta_phi_deta1_file2_name = prefix2 + "delta_phi_deta1";
    TString delta_phi_deta1_file3_name = prefix3 + "delta_phi_deta1";

    file1->GetObject(delta_phi_deta1_file1_name,delta_phi_deta1_file1);
    if (delta_phi_deta1_file1 == 0) { cout << delta_phi_deta1_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta1_file2_name,delta_phi_deta1_file2);
    if (delta_phi_deta1_file2 == 0) { cout << delta_phi_deta1_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta1_file3_name,delta_phi_deta1_file3);
    	if (delta_phi_deta1_file3 == 0) { cout << delta_phi_deta1_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta1_file1, label1, delta_phi_deta1_file2, label2, delta_phi_deta1_file3, label3, output_path, output_prefix + "delta_phi_deta1", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta1_file1, label1, delta_phi_deta1_file2, label2, output_path, output_prefix + "delta_phi_deta1", "top_right", false, detail); }


//plot the correction factor for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_file1 = 0;
    TH1D *delta_phi_deta2_file2 = 0;
    TH1D *delta_phi_deta2_file3 = 0;
    TString delta_phi_deta2_file1_name = prefix1 + "delta_phi_deta2";
    TString delta_phi_deta2_file2_name = prefix2 + "delta_phi_deta2";
    TString delta_phi_deta2_file3_name = prefix3 + "delta_phi_deta2";

    file1->GetObject(delta_phi_deta2_file1_name,delta_phi_deta2_file1);
    if (delta_phi_deta2_file1 == 0) { cout << delta_phi_deta2_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta2_file2_name,delta_phi_deta2_file2);
    if (delta_phi_deta2_file2 == 0) { cout << delta_phi_deta2_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta2_file3_name,delta_phi_deta2_file3);
    	if (delta_phi_deta2_file3 == 0) { cout << delta_phi_deta2_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta2_file1, label1, delta_phi_deta2_file2, label2, delta_phi_deta2_file3, label3, output_path, output_prefix + "delta_phi_deta2", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta2_file1, label1, delta_phi_deta2_file2, label2, output_path, output_prefix + "delta_phi_deta2", "top_right", false, detail); }


//plot the correction factor for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_file1 = 0;
    TH1D *delta_phi_deta3_file2 = 0;
    TH1D *delta_phi_deta3_file3 = 0;
    TString delta_phi_deta3_file1_name = prefix1 + "delta_phi_deta3";
    TString delta_phi_deta3_file2_name = prefix2 + "delta_phi_deta3";
    TString delta_phi_deta3_file3_name = prefix3 + "delta_phi_deta3";

    file1->GetObject(delta_phi_deta3_file1_name,delta_phi_deta3_file1);
    if (delta_phi_deta3_file1 == 0) { cout << delta_phi_deta3_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta3_file2_name,delta_phi_deta3_file2);
    if (delta_phi_deta3_file2 == 0) { cout << delta_phi_deta3_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta3_file3_name,delta_phi_deta3_file3);
    	if (delta_phi_deta3_file3 == 0) { cout << delta_phi_deta3_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta3_file1, label1, delta_phi_deta3_file2, label2, delta_phi_deta3_file3, label3, output_path, output_prefix + "delta_phi_deta3", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta3_file1, label1, delta_phi_deta3_file2, label2, output_path, output_prefix + "delta_phi_deta3", "top_right", false, detail); }


//plot the correction factor for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_file1 = 0;
    TH1D *delta_phi_deta4_file2 = 0;
    TH1D *delta_phi_deta4_file3 = 0;
    TString delta_phi_deta4_file1_name = prefix1 + "delta_phi_deta4";
    TString delta_phi_deta4_file2_name = prefix2 + "delta_phi_deta4";
    TString delta_phi_deta4_file3_name = prefix3 + "delta_phi_deta4";

    file1->GetObject(delta_phi_deta4_file1_name,delta_phi_deta4_file1);
    if (delta_phi_deta4_file1 == 0) { cout << delta_phi_deta4_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta4_file2_name,delta_phi_deta4_file2);
    if (delta_phi_deta4_file2 == 0) { cout << delta_phi_deta4_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta4_file3_name,delta_phi_deta4_file3);
    	if (delta_phi_deta4_file3 == 0) { cout << delta_phi_deta4_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta4_file1, label1, delta_phi_deta4_file2, label2, delta_phi_deta4_file3, label3, output_path, output_prefix + "delta_phi_deta4", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta4_file1, label1, delta_phi_deta4_file2, label2, output_path, output_prefix + "delta_phi_deta4", "top_right", false, detail); }


//plot the correction factor for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_file1 = 0;
    TH1D *delta_phi_gap_file2 = 0;
    TH1D *delta_phi_gap_file3 = 0;
    TString delta_phi_gap_file1_name = prefix1 + "delta_phi_gap";
    TString delta_phi_gap_file2_name = prefix2 + "delta_phi_gap";
    TString delta_phi_gap_file3_name = prefix3 + "delta_phi_gap";

    file1->GetObject(delta_phi_gap_file1_name,delta_phi_gap_file1);
    if (delta_phi_gap_file1 == 0) { cout << delta_phi_gap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_gap_file2_name,delta_phi_gap_file2);
    if (delta_phi_gap_file2 == 0) { cout << delta_phi_gap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_gap_file3_name,delta_phi_gap_file3);
    	if (delta_phi_gap_file3 == 0) { cout << delta_phi_gap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_gap_file1, label1, delta_phi_gap_file2, label2, delta_phi_gap_file3, label3, output_path, output_prefix + "delta_phi_gap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_gap_file1, label1, delta_phi_gap_file2, label2, output_path, output_prefix + "delta_phi_gap", "top_right", false, detail); }


//plot the correction factor for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_file1 = 0;
    TH1D *delta_phi_deta1_gap_file2 = 0;
    TH1D *delta_phi_deta1_gap_file3 = 0;
    TString delta_phi_deta1_gap_file1_name = prefix1 + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_file2_name = prefix2 + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_file3_name = prefix3 + "delta_phi_deta1_gap";

    file1->GetObject(delta_phi_deta1_gap_file1_name,delta_phi_deta1_gap_file1);
    if (delta_phi_deta1_gap_file1 == 0) { cout << delta_phi_deta1_gap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta1_gap_file2_name,delta_phi_deta1_gap_file2);
    if (delta_phi_deta1_gap_file2 == 0) { cout << delta_phi_deta1_gap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta1_gap_file3_name,delta_phi_deta1_gap_file3);
    	if (delta_phi_deta1_gap_file3 == 0) { cout << delta_phi_deta1_gap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta1_gap_file1, label1, delta_phi_deta1_gap_file2, label2, delta_phi_deta1_gap_file3, label3, output_path, output_prefix + "delta_phi_deta1_gap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta1_gap_file1, label1, delta_phi_deta1_gap_file2, label2, output_path, output_prefix + "delta_phi_deta1_gap", "top_right", false, detail); }


//plot the correction factor for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_file1 = 0;
    TH1D *delta_phi_deta2_gap_file2 = 0;
    TH1D *delta_phi_deta2_gap_file3 = 0;
    TString delta_phi_deta2_gap_file1_name = prefix1 + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_file2_name = prefix2 + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_file3_name = prefix3 + "delta_phi_deta2_gap";

    file1->GetObject(delta_phi_deta2_gap_file1_name,delta_phi_deta2_gap_file1);
    if (delta_phi_deta2_gap_file1 == 0) { cout << delta_phi_deta2_gap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta2_gap_file2_name,delta_phi_deta2_gap_file2);
    if (delta_phi_deta2_gap_file2 == 0) { cout << delta_phi_deta2_gap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta2_gap_file3_name,delta_phi_deta2_gap_file3);
    	if (delta_phi_deta2_gap_file3 == 0) { cout << delta_phi_deta2_gap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta2_gap_file1, label1, delta_phi_deta2_gap_file2, label2, delta_phi_deta2_gap_file3, label3, output_path, output_prefix + "delta_phi_deta2_gap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta2_gap_file1, label1, delta_phi_deta2_gap_file2, label2, output_path, output_prefix + "delta_phi_deta2_gap", "top_right", false, detail); }


//plot the correction factor for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_file1 = 0;
    TH1D *delta_phi_deta3_gap_file2 = 0;
    TH1D *delta_phi_deta3_gap_file3 = 0;
    TString delta_phi_deta3_gap_file1_name = prefix1 + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_file2_name = prefix2 + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_file3_name = prefix3 + "delta_phi_deta3_gap";

    file1->GetObject(delta_phi_deta3_gap_file1_name,delta_phi_deta3_gap_file1);
    if (delta_phi_deta3_gap_file1 == 0) { cout << delta_phi_deta3_gap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta3_gap_file2_name,delta_phi_deta3_gap_file2);
    if (delta_phi_deta3_gap_file2 == 0) { cout << delta_phi_deta3_gap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta3_gap_file3_name,delta_phi_deta3_gap_file3);
    	if (delta_phi_deta3_gap_file3 == 0) { cout << delta_phi_deta3_gap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta3_gap_file1, label1, delta_phi_deta3_gap_file2, label2, delta_phi_deta3_gap_file3, label3, output_path, output_prefix + "delta_phi_deta3_gap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta3_gap_file1, label1, delta_phi_deta3_gap_file2, label2, output_path, output_prefix + "delta_phi_deta3_gap", "top_right", false, detail); }


//plot the correction factor for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_file1 = 0;
    TH1D *delta_phi_deta4_gap_file2 = 0;
    TH1D *delta_phi_deta4_gap_file3 = 0;
    TString delta_phi_deta4_gap_file1_name = prefix1 + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_file2_name = prefix2 + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_file3_name = prefix3 + "delta_phi_deta4_gap";

    file1->GetObject(delta_phi_deta4_gap_file1_name,delta_phi_deta4_gap_file1);
    if (delta_phi_deta4_gap_file1 == 0) { cout << delta_phi_deta4_gap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta4_gap_file2_name,delta_phi_deta4_gap_file2);
    if (delta_phi_deta4_gap_file2 == 0) { cout << delta_phi_deta4_gap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta4_gap_file3_name,delta_phi_deta4_gap_file3);
    	if (delta_phi_deta4_gap_file3 == 0) { cout << delta_phi_deta4_gap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta4_gap_file1, label1, delta_phi_deta4_gap_file2, label2, delta_phi_deta4_gap_file3, label3, output_path, output_prefix + "delta_phi_deta4_gap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta4_gap_file1, label1, delta_phi_deta4_gap_file2, label2, output_path, output_prefix + "delta_phi_deta4_gap", "top_right", false, detail); }


//plot the correction factor for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_file1 = 0;
    TH1D *delta_phi_nogap_file2 = 0;
    TH1D *delta_phi_nogap_file3 = 0;
    TString delta_phi_nogap_file1_name = prefix1 + "delta_phi_nogap";
    TString delta_phi_nogap_file2_name = prefix2 + "delta_phi_nogap";
    TString delta_phi_nogap_file3_name = prefix3 + "delta_phi_nogap";

    file1->GetObject(delta_phi_nogap_file1_name,delta_phi_nogap_file1);
    if (delta_phi_nogap_file1 == 0) { cout << delta_phi_nogap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_nogap_file2_name,delta_phi_nogap_file2);
    if (delta_phi_nogap_file2 == 0) { cout << delta_phi_nogap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_nogap_file3_name,delta_phi_nogap_file3);
    	if (delta_phi_nogap_file3 == 0) { cout << delta_phi_nogap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_nogap_file1, label1, delta_phi_nogap_file2, label2, delta_phi_nogap_file3, label3, output_path, output_prefix + "delta_phi_nogap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_nogap_file1, label1, delta_phi_nogap_file2, label2, output_path, output_prefix + "delta_phi_nogap", "top_right", false, detail); }


//plot the correction factor for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_file1 = 0;
    TH1D *delta_phi_deta1_nogap_file2 = 0;
    TH1D *delta_phi_deta1_nogap_file3 = 0;
    TString delta_phi_deta1_nogap_file1_name = prefix1 + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_file2_name = prefix2 + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_file3_name = prefix3 + "delta_phi_deta1_nogap";

    file1->GetObject(delta_phi_deta1_nogap_file1_name,delta_phi_deta1_nogap_file1);
    if (delta_phi_deta1_nogap_file1 == 0) { cout << delta_phi_deta1_nogap_file1_name << " not found!" << endl; return; }
    file2->GetObject(delta_phi_deta1_nogap_file2_name,delta_phi_deta1_nogap_file2);
    if (delta_phi_deta1_nogap_file2 == 0) { cout << delta_phi_deta1_nogap_file2_name << " not found!" << endl; return; }
    if (have3files)
	{
    	file3->GetObject(delta_phi_deta1_nogap_file3_name,delta_phi_deta1_nogap_file3);
    	if (delta_phi_deta1_nogap_file3 == 0) { cout << delta_phi_deta1_nogap_file3_name << " not found!" << endl; return; }
	}

    if (have3files) { plot_3histograms(delta_phi_deta1_nogap_file1, label1, delta_phi_deta1_nogap_file2, label2, delta_phi_deta1_nogap_file3, label3, output_path, output_prefix + "delta_phi_deta1_nogap", "top_right", false, detail); }
    else { plot_2histograms(delta_phi_deta1_nogap_file1, label1, delta_phi_deta1_nogap_file2, label2, output_path, output_prefix + "delta_phi_deta1_nogap", "top_right", false, detail); }


}


void compute_correction(string input_file_uncorrected = "../output/xsec_mc_gen/xsec_Pythia6_TuneZ2star_gen_allvertex.root", string hist_prefix_unc = "ak5Gen_", TString label_unc = "Pythia6 TuneZ2* - All vertex on Generator Level", string input_file_true = "../output/normalized_mc/xsec_p6_z2_gen.root",  string hist_prefix_true = "ak5Gen_", TString label_true = "Pythia6 TuneZ2* - 1 vertex on Generator Level", string output_path_plots = "../output/corrections/", string output_rootfile = "../output/histograms/corrections/pileup_p6_z2.root", string plot_prefix = "pileup_pythia6_z2_", TString legend_label = "Pileup Correction with Pythia 6 - Tune Z2*", bool merge = false, bool detail = false, bool show_stats = true)
{
//computes the correction factors for detector effects or pile up using the uncorrected and true distributions

   if (detail) { cout<<"Input File with uncorrected data : "<<input_file_uncorrected<<endl; }
   if (detail) { cout<<"Prefix for uncorrected data :      "<<hist_prefix_unc<<endl; }
   if (detail) { cout<<"Label for uncorrected data :       "<<label_unc<<endl; }
   if (detail) { cout<<"Input File with true data :        "<<input_file_true<<endl; }
   if (detail) { cout<<"Prefix for true data :             "<<hist_prefix_true<<endl; }
   if (detail) { cout<<"Label for true data :              "<<label_true<<endl; }
   if (detail) { cout<<"Output Folder for plots :          "<<output_path_plots<<endl; }
   if (detail) { cout<<"Output Rootfile :                  "<<output_rootfile<<endl; }
   if (detail) { cout<<"Prefix :                           "<<plot_prefix<<endl; }
   if (detail) { cout<<"Legend Label :                     "<<legend_label<<endl; }
   if (detail) { cout<<"Detail :                           "<<detail<<endl; }
   if (detail) { cout<<"Show Statistics :                  "<<show_stats<<endl; }

//opening the input and output data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *mc_uncorrected = new TFile( input_file_uncorrected.c_str() );
    TFile *mc_true = new TFile( input_file_true.c_str() );

//declare histogram bins
int cent_nbins = 7;
int forw_nbins = 7;

double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};

int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int deta_nbins = 4;
int dphi_nbins = 7;

double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

   string  hist_name;
   TString hist_name_unc;
   TString hist_name_true;
   TString hist_name_out;

//correction factor for the delta phi distribution
    hist_name = "delta_phi";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi;
    delta_phi = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi->Sumw2();
    
    new_correction(delta_phi, delta_phi_uncorrected, delta_phi_true, merge, detail);
    plot_correction(delta_phi, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_uncorrected, label_unc, delta_phi_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);
    
//correction factor for the delta phi deta1 distribution
    hist_name = "delta_phi_deta1";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta1_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta1_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta1;
    delta_phi_deta1 = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta1->Sumw2();
    
    new_correction(delta_phi_deta1, delta_phi_deta1_uncorrected, delta_phi_deta1_true, merge, detail);
    plot_correction(delta_phi_deta1, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta1_uncorrected, label_unc, delta_phi_deta1_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);
    
//correction factor for the delta phi deta2 distribution
    hist_name = "delta_phi_deta2";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta2_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta2_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta2;
    delta_phi_deta2 = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta2->Sumw2();
    
    new_correction(delta_phi_deta2, delta_phi_deta2_uncorrected, delta_phi_deta2_true, merge, detail);
    plot_correction(delta_phi_deta2, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta2_uncorrected, label_unc, delta_phi_deta2_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);
    
//correction factor for the delta phi deta3 distribution
    hist_name = "delta_phi_deta3";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta3_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta3_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta3;
    delta_phi_deta3 = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta3->Sumw2();
    
    new_correction(delta_phi_deta3, delta_phi_deta3_uncorrected, delta_phi_deta3_true, merge, detail);
    plot_correction(delta_phi_deta3, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta3_uncorrected, label_unc, delta_phi_deta3_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name,  "top_left", true, detail);
    
//correction factor for the delta phi deta4 distribution
    hist_name = "delta_phi_deta4";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta4_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta4_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta4;
    delta_phi_deta4 = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta4->Sumw2();
    
    new_correction(delta_phi_deta4, delta_phi_deta4_uncorrected, delta_phi_deta4_true, merge, detail);
    plot_correction(delta_phi_deta4, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta4_uncorrected, label_unc, delta_phi_deta4_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi gap distribution
    hist_name = "delta_phi_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_gap;
    delta_phi_gap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_gap->Sumw2();
    
    new_correction(delta_phi_gap, delta_phi_gap_uncorrected, delta_phi_gap_true, merge, detail);
    plot_correction(delta_phi_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_gap_uncorrected, label_unc, delta_phi_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);
    
//correction factor for the delta phi deta1 gap distribution
    hist_name = "delta_phi_deta1_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta1_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta1_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap->Sumw2();
    
    new_correction(delta_phi_deta1_gap, delta_phi_deta1_gap_uncorrected, delta_phi_deta1_gap_true, merge, detail);
    plot_correction(delta_phi_deta1_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta1_gap_uncorrected, label_unc, delta_phi_deta1_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta2 gap distribution
    hist_name = "delta_phi_deta2_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta2_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta2_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap->Sumw2();
    
    new_correction(delta_phi_deta2_gap, delta_phi_deta2_gap_uncorrected, delta_phi_deta2_gap_true, merge, detail);
    plot_correction(delta_phi_deta2_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta2_gap_uncorrected, label_unc, delta_phi_deta2_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta3 gap distribution
    hist_name = "delta_phi_deta3_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta3_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta3_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap->Sumw2();
    
    new_correction(delta_phi_deta3_gap, delta_phi_deta3_gap_uncorrected, delta_phi_deta3_gap_true, merge, detail);
    plot_correction(delta_phi_deta3_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta3_gap_uncorrected, label_unc, delta_phi_deta3_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta4 gap distribution
    hist_name = "delta_phi_deta4_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta4_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta4_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap->Sumw2();
    
    new_correction(delta_phi_deta4_gap, delta_phi_deta4_gap_uncorrected, delta_phi_deta4_gap_true, merge, detail);
    plot_correction(delta_phi_deta4_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta4_gap_uncorrected, label_unc, delta_phi_deta4_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi nogap distribution
    hist_name = "delta_phi_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_nogap;
    delta_phi_nogap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_nogap->Sumw2();
    
    new_correction(delta_phi_nogap, delta_phi_nogap_uncorrected, delta_phi_nogap_true, merge, detail);
    plot_correction(delta_phi_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_nogap_uncorrected, label_unc, delta_phi_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta1 nogap distribution
    hist_name = "delta_phi_deta1_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta1_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta1_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap->Sumw2();
    
    new_correction(delta_phi_deta1_nogap, delta_phi_deta1_nogap_uncorrected, delta_phi_deta1_nogap_true, merge, detail);
    plot_correction(delta_phi_deta1_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta1_nogap_uncorrected, label_unc, delta_phi_deta1_nogap_true, label_true, output_path_plots,  "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta2 nogap distribution
    hist_name = "delta_phi_deta2_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta2_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta2_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap->Sumw2();
    
    new_correction(delta_phi_deta2_nogap, delta_phi_deta2_nogap_uncorrected, delta_phi_deta2_nogap_true, merge, detail);
    plot_correction(delta_phi_deta2_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta2_nogap_uncorrected, label_unc, delta_phi_deta2_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta3 nogap distribution
    hist_name = "delta_phi_deta3_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta3_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta3_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap->Sumw2();
    
    new_correction(delta_phi_deta3_nogap, delta_phi_deta3_nogap_uncorrected, delta_phi_deta3_nogap_true, merge, detail);
    plot_correction(delta_phi_deta3_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta3_nogap_uncorrected, label_unc, delta_phi_deta3_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta phi deta4 nogap distribution
    hist_name = "delta_phi_deta4_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_phi_deta4_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_phi_deta4_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap = new TH1D(hist_name_out,"#Delta#phi;|#Delta#phi| [rad];Correction Factor", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap->Sumw2();
    
    new_correction(delta_phi_deta4_nogap, delta_phi_deta4_nogap_uncorrected, delta_phi_deta4_nogap_true, merge, detail);
    plot_correction(delta_phi_deta4_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_right", detail);
    plot_2histograms(delta_phi_deta4_nogap_uncorrected, label_unc, delta_phi_deta4_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "top_left", true, detail);

//correction factor for the delta eta distribution
    hist_name = "delta_eta";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_eta_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_eta_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_eta;
    delta_eta = new TH1D(hist_name_out,"#Delta#eta;|#Delta#eta| [rad];Correction Factor", deta_nbins, deta_bins);
    delta_eta->Sumw2();
    
    new_correction(delta_eta, delta_eta_uncorrected, delta_eta_true, merge, detail);
    plot_correction(delta_eta, plot_prefix, hist_name, output_path_plots, legend_label, "bottom_left", detail);
    plot_2histograms(delta_eta_uncorrected, label_unc, delta_eta_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the delta eta gap distribution
    hist_name = "delta_eta_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_eta_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_eta_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_eta_gap;
    delta_eta_gap = new TH1D(hist_name_out,"#Delta#eta when a gap is present;|#Delta#eta| [rad];Correction Factor", deta_nbins, deta_bins);
    delta_eta_gap->Sumw2();
    
    new_correction(delta_eta_gap, delta_eta_gap_uncorrected, delta_eta_gap_true, merge, detail);
    plot_correction(delta_eta_gap, plot_prefix, hist_name, output_path_plots, legend_label, "bottom_left", detail);
    plot_2histograms(delta_eta_gap_uncorrected, label_unc, delta_eta_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the delta eta nogap distribution
    hist_name = "delta_eta_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_eta_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_eta_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_eta_nogap;
    delta_eta_nogap = new TH1D(hist_name_out,"#Delta#eta when there is no gap;|#Delta#eta| [rad];Correction Factor", deta_nbins, deta_bins);
    delta_eta_nogap->Sumw2();
    
    new_correction(delta_eta_nogap, delta_eta_nogap_uncorrected, delta_eta_nogap_true, merge, detail);
    plot_correction(delta_eta_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "bottom_left", detail);
    plot_2histograms(delta_eta_nogap_uncorrected, label_unc, delta_eta_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the total pt inside gap distribution
    hist_name = "total_pt_inside_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *total_pt_inside_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *total_pt_inside_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *total_pt_inside_gap;
    total_pt_inside_gap = new TH1D(hist_name_out,"Total Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];Correction Factor", in_nbins, in_bins);
    total_pt_inside_gap->Sumw2();
    
    new_correction(total_pt_inside_gap, total_pt_inside_gap_uncorrected, total_pt_inside_gap_true, merge, detail);
    plot_correction(total_pt_inside_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(total_pt_inside_gap_uncorrected, label_unc, total_pt_inside_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading pt inside gap distribution
    hist_name = "leading_pt_inside_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_pt_inside_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_pt_inside_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap = new TH1D(hist_name_out,"Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];Correction Factor", in_nbins, in_bins);
    leading_pt_inside_gap->Sumw2();
    
    new_correction(leading_pt_inside_gap, leading_pt_inside_gap_uncorrected, leading_pt_inside_gap_true, merge, detail);
    plot_correction(leading_pt_inside_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_pt_inside_gap_uncorrected, label_unc, leading_pt_inside_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the total pt outside gap distribution
    hist_name = "total_pt_outside_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *total_pt_outside_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *total_pt_outside_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *total_pt_outside_gap;
    total_pt_outside_gap = new TH1D(hist_name_out,"Total Jet p_{T} outside the gap;p_{T} [#frac{GeV}{c}];Correction Factor", out_nbins, out_bins);
    total_pt_outside_gap->Sumw2();
    
    new_correction(total_pt_outside_gap, total_pt_outside_gap_uncorrected, total_pt_outside_gap_true, merge, detail);
    plot_correction(total_pt_outside_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(total_pt_outside_gap_uncorrected, label_unc, total_pt_outside_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading pt outside gap distribution
    hist_name = "leading_pt_outside_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_pt_outside_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_pt_outside_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap = new TH1D(hist_name_out,"Leading Jet p_{T} outside the gap;p_{T} [#frac{GeV}{c}];Correction Factor", out_nbins, out_bins);
    leading_pt_outside_gap->Sumw2();
    
    new_correction(leading_pt_outside_gap, leading_pt_outside_gap_uncorrected, leading_pt_outside_gap_true, merge, detail);
    plot_correction(leading_pt_outside_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_pt_outside_gap_uncorrected, label_unc, leading_pt_outside_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the delta eta outside gap distribution
    hist_name = "delta_eta_outside_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *delta_eta_outside_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *delta_eta_outside_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap = new TH1D(hist_name_out,"#Delta#eta outside the gap;#Delta#eta;Correction Factor", deta_out_nbins, deta_out_bins);
    delta_eta_outside_gap->Sumw2();
    
    new_correction(delta_eta_outside_gap, delta_eta_outside_gap_uncorrected, delta_eta_outside_gap_true, merge, detail);
    plot_correction(delta_eta_outside_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(delta_eta_outside_gap_uncorrected, label_unc, delta_eta_outside_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading eta star inside gap distribution
    hist_name = "leading_eta_star_inside_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_eta_star_inside_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_eta_star_inside_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap = new TH1D(hist_name_out,"Leading Jet #eta* inside the gap;#eta*;Correction Factor", etastar_nbins, etastar_bins);
    leading_eta_star_inside_gap->Sumw2();
    
    new_correction(leading_eta_star_inside_gap, leading_eta_star_inside_gap_uncorrected, leading_eta_star_inside_gap_true, merge, detail);
    plot_correction(leading_eta_star_inside_gap, plot_prefix, hist_name, output_path_plots, legend_label, "bottom_right", detail);
    plot_2histograms(leading_eta_star_inside_gap_uncorrected, label_unc, leading_eta_star_inside_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_middle", true, detail);

//correction factor for the leading central pt distribution
    hist_name = "leading_central_pt";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_central_pt_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_central_pt_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_central_pt;
    leading_central_pt = new TH1D(hist_name_out,"Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];Correction Factor", cent_nbins, cent_bins);
    leading_central_pt->Sumw2();
    
    new_correction(leading_central_pt, leading_central_pt_uncorrected, leading_central_pt_true, merge, detail);
    plot_correction(leading_central_pt, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_central_pt_uncorrected, label_unc, leading_central_pt_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading forward pt distribution
    hist_name = "leading_forward_pt";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_forward_pt_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_forward_pt_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_forward_pt;
    leading_forward_pt = new TH1D(hist_name_out,"Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];Correction Factor", forw_nbins, forw_bins);
    leading_forward_pt->Sumw2();
    
    new_correction(leading_forward_pt, leading_forward_pt_uncorrected, leading_forward_pt_true, merge, detail);
    plot_correction(leading_forward_pt, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_forward_pt_uncorrected, label_unc, leading_forward_pt_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading central pt gap distribution
    hist_name = "leading_central_pt_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_central_pt_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_central_pt_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_central_pt_gap;
    leading_central_pt_gap = new TH1D(hist_name_out,"Leading Central Jet p_{T} when a gap is present;p_{T} [#frac{GeV}{c}];Correction Factor", cent_nbins, cent_bins);
    leading_central_pt_gap->Sumw2();
    
    new_correction(leading_central_pt_gap, leading_central_pt_gap_uncorrected, leading_central_pt_gap_true, merge, detail);
    plot_correction(leading_central_pt_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_central_pt_gap_uncorrected, label_unc, leading_central_pt_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading forward pt gap distribution
    hist_name = "leading_forward_pt_gap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_forward_pt_gap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_forward_pt_gap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_forward_pt_gap;
    leading_forward_pt_gap = new TH1D(hist_name_out,"Leading Forward Jet p_{T} when a gap is present;p_{T} [#frac{GeV}{c}];Correction Factor", forw_nbins, forw_bins);
    leading_forward_pt_gap->Sumw2();
    
    new_correction(leading_forward_pt_gap, leading_forward_pt_gap_uncorrected, leading_forward_pt_gap_true, merge, detail);
    plot_correction(leading_forward_pt_gap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_forward_pt_gap_uncorrected, label_unc, leading_forward_pt_gap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

//correction factor for the leading central pt nogap distribution
    hist_name = "leading_central_pt_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_central_pt_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_central_pt_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_central_pt_nogap;
    leading_central_pt_nogap = new TH1D(hist_name_out,"Leading Central Jet p_{T} when a gap is not present;p_{T} [#frac{GeV}{c}];Correction Factor", cent_nbins, cent_bins);
    leading_central_pt_nogap->Sumw2();
    
    new_correction(leading_central_pt_nogap, leading_central_pt_nogap_uncorrected, leading_central_pt_nogap_true, merge, detail);
    plot_correction(leading_central_pt_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_central_pt_nogap_uncorrected, label_unc, leading_central_pt_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name,  "bottom_left", true, detail);

//correction factor for the leading forward pt nogap distribution
    hist_name = "leading_forward_pt_nogap";
    hist_name_unc = hist_prefix_unc + hist_name;
    hist_name_true = hist_prefix_true + hist_name;
    hist_name_out = plot_prefix + hist_name;
    TH1D *leading_forward_pt_nogap_uncorrected = (TH1D*) mc_uncorrected->Get(hist_name_unc);
    TH1D *leading_forward_pt_nogap_true = (TH1D*) mc_true->Get(hist_name_true);

    TH1D *leading_forward_pt_nogap;
    leading_forward_pt_nogap = new TH1D(hist_name_out,"Leading Forward Jet p_{T} when a gap is not present;p_{T} [#frac{GeV}{c}];Correction Factor", forw_nbins, forw_bins);
    leading_forward_pt_nogap->Sumw2();
    
    new_correction(leading_forward_pt_nogap, leading_forward_pt_nogap_uncorrected, leading_forward_pt_nogap_true, merge, detail);
    plot_correction(leading_forward_pt_nogap, plot_prefix, hist_name, output_path_plots, legend_label, "top_left", detail);
    plot_2histograms(leading_forward_pt_nogap_uncorrected, label_unc, leading_forward_pt_nogap_true, label_true, output_path_plots, "control_" + plot_prefix + hist_name, "bottom_left", true, detail);

    //creating the output file
    TFile data_output( output_rootfile.c_str() , "RECREATE");
    
    //save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file "<<output_rootfile<<" ..."<<endl; }
    delta_phi->Write();
    delta_phi_deta1->Write();
    delta_phi_deta2->Write();
    delta_phi_deta3->Write();
    delta_phi_deta4->Write();
    delta_phi_gap->Write();
    delta_phi_deta1_gap->Write();
    delta_phi_deta2_gap->Write();
    delta_phi_deta3_gap->Write();
    delta_phi_deta4_gap->Write();
    delta_phi_nogap->Write();
    delta_phi_deta1_nogap->Write();
    delta_phi_deta2_nogap->Write();
    delta_phi_deta3_nogap->Write();
    delta_phi_deta4_nogap->Write();
    delta_eta->Write();
    delta_eta_gap->Write();
    delta_eta_nogap->Write();
    total_pt_inside_gap->Write();
    leading_pt_inside_gap->Write();
    total_pt_outside_gap->Write();
    leading_pt_outside_gap->Write();
    delta_eta_outside_gap->Write();
    leading_eta_star_inside_gap->Write();
    leading_central_pt->Write();
    leading_central_pt_gap->Write();
    leading_central_pt_nogap->Write();
    leading_forward_pt->Write();
    leading_forward_pt_gap->Write();
    leading_forward_pt_nogap->Write();
    if (detail) { cout<<"Writing on "<<output_rootfile<<" was sucessfull!"<<endl; }
    
    //close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    mc_uncorrected->Close();
    mc_true->Close();
    data_output.Close();
    if (detail) { cout<<"Done!"<<endl; }

}
