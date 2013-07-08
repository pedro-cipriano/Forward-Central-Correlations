// DESY, CMS
// Last Update: 24 Mar 2013
//
// apply correction()

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TString.h>
#include <TF1.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "common_methods.h"

using namespace std;


void correct_histogram(TH1D *corrected, TH1D *data, TH1D *p6_z2, TH1D *p8_4c, bool detail = false)
{

double factor_p6_z2 = 0.0;
double factor_p8_4c = 0.0;
double correction_factor = 0.0;
double xsec = 0.0;
double error = 0.0;
double xsec_final = 0.0;
double error_final = 0.0;

for (int i = 1; i <= corrected->GetNbinsX();i++)
	{
	factor_p6_z2 = p6_z2->GetBinContent(i);
	factor_p8_4c = p8_4c->GetBinContent(i);
	correction_factor = (factor_p6_z2 + factor_p8_4c)/2.0;
	xsec = data->GetBinContent(i);
	error = data->GetBinError(i);
	xsec_final = xsec * correction_factor;
	error_final = error * correction_factor;
	corrected->SetBinContent(i,xsec_final);
	corrected->SetBinError(i,error_final);
	if (detail) { cout << "Uncorrected: " << xsec << "+-" << error << " Corrections: " << factor_p6_z2 << "/" << factor_p8_4c << "=" << correction_factor << " Corrected:" << xsec_final << "+-" << error_final << endl; }
	}

corrected->SetEntries(data->GetEntries());

}


void apply_correction(string path_data, string path_correction_p6_z2, string label_p6_z2, string path_correction_p8_4c, string label_p8_4c, string path_corrected, string output_path_plots, string plot_label, bool detail = false, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Apply Correction Configuration"<<endl; }
    if (detail) { cout << "Input path for data:                  " << path_data << endl; }
    if (detail) { cout << "Input path for correction Pythia6 Z2: " << path_correction_p6_z2 << endl; }
    if (detail) { cout << "Label for Pythia6 Z2:                 " << label_p6_z2 << endl; }
    if (detail) { cout << "Input path for correction Pythia8 4C: " << path_correction_p8_4c << endl; }
    if (detail) { cout << "Label for Pythia8 4C:                 " << label_p8_4c << endl; }
    if (detail) { cout << "Output path:                          " << path_corrected << endl; }
    if (detail) { cout << "Output Path Plots:                    " << output_path_plots << endl; }
    if (detail) { cout << "Plot Label:                           " << plot_label << endl; }
    if (detail) { cout << "Detail Level:                         " << detail << endl; }
    if (detail) { cout << "Test Mode:                            " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data = new TFile( path_data.c_str() );
    TFile *p6_z2 = new TFile( path_correction_p6_z2.c_str() );
    TFile *p8_4c = new TFile( path_correction_p8_4c.c_str() );


//histogram bins
int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};


//compute corrected data for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_data = 0;
    TH1D *delta_phi_p6_z2 = 0;
    TH1D *delta_phi_p8_4c = 0;
    TString delta_phi_name_p6_z2 = label_p6_z2 + "delta_phi";
    TString delta_phi_name_p8_4c = label_p8_4c + "delta_phi";

    data->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_name_p6_z2,delta_phi_p6_z2);
    if (delta_phi_p6_z2 == 0) { cout << delta_phi_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_name_p8_4c,delta_phi_p8_4c);
    if (delta_phi_p8_4c == 0) { cout << delta_phi_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi;
    delta_phi =  new TH1D("ak5PF_delta_phi","Corrected Data;#Delta#phi;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi, delta_phi_data, delta_phi_p6_z2, delta_phi_p8_4c, detail);
    plot_2histograms(delta_phi, "Data Corrected", delta_phi_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi", "top_left", detail);


//compute corrected data for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_data = 0;
    TH1D *delta_phi_deta1_p6_z2 = 0;
    TH1D *delta_phi_deta1_p8_4c = 0;
    TString delta_phi_deta1_name_p6_z2 = label_p6_z2 + "delta_phi_deta1";
    TString delta_phi_deta1_name_p8_4c = label_p8_4c + "delta_phi_deta1";

    data->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta1_name_p6_z2,delta_phi_deta1_p6_z2);
    if (delta_phi_deta1_p6_z2 == 0) { cout << delta_phi_deta1_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta1_name_p8_4c,delta_phi_deta1_p8_4c);
    if (delta_phi_deta1_p8_4c == 0) { cout << delta_phi_deta1_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D("ak5PF_delta_phi_deta1","Corrected Data;#Delta#phi deta1;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta1, delta_phi_deta1_data, delta_phi_deta1_p6_z2, delta_phi_deta1_p8_4c, detail);
    plot_2histograms(delta_phi_deta1, "Data Corrected", delta_phi_deta1_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta1", "top_left", detail);


//compute corrected data for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_data = 0;
    TH1D *delta_phi_deta2_p6_z2 = 0;
    TH1D *delta_phi_deta2_p8_4c = 0;
    TString delta_phi_deta2_name_p6_z2 = label_p6_z2 + "delta_phi_deta2";
    TString delta_phi_deta2_name_p8_4c = label_p8_4c + "delta_phi_deta2";

    data->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta2_name_p6_z2,delta_phi_deta2_p6_z2);
    if (delta_phi_deta2_p6_z2 == 0) { cout << delta_phi_deta2_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta2_name_p8_4c,delta_phi_deta2_p8_4c);
    if (delta_phi_deta2_p8_4c == 0) { cout << delta_phi_deta2_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D("ak5PF_delta_phi_deta2","Corrected Data;#Delta#phi deta2;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta2, delta_phi_deta2_data, delta_phi_deta2_p6_z2, delta_phi_deta2_p8_4c, detail);
    plot_2histograms(delta_phi_deta2, "Data Corrected", delta_phi_deta2_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta2", "top_left", detail);


//compute corrected data for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_data = 0;
    TH1D *delta_phi_deta3_p6_z2 = 0;
    TH1D *delta_phi_deta3_p8_4c = 0;
    TString delta_phi_deta3_name_p6_z2 = label_p6_z2 + "delta_phi_deta3";
    TString delta_phi_deta3_name_p8_4c = label_p8_4c + "delta_phi_deta3";

    data->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta3_name_p6_z2,delta_phi_deta3_p6_z2);
    if (delta_phi_deta3_p6_z2 == 0) { cout << delta_phi_deta3_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta3_name_p8_4c,delta_phi_deta3_p8_4c);
    if (delta_phi_deta3_p8_4c == 0) { cout << delta_phi_deta3_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D("ak5PF_delta_phi_deta3","Corrected Data;#Delta#phi deta3;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta3, delta_phi_deta3_data, delta_phi_deta3_p6_z2, delta_phi_deta3_p8_4c, detail);
    plot_2histograms(delta_phi_deta3, "Data Corrected", delta_phi_deta3_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta3", "top_left", detail);


//compute corrected data for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_data = 0;
    TH1D *delta_phi_deta4_p6_z2 = 0;
    TH1D *delta_phi_deta4_p8_4c = 0;
    TString delta_phi_deta4_name_p6_z2 = label_p6_z2 + "delta_phi_deta4";
    TString delta_phi_deta4_name_p8_4c = label_p8_4c + "delta_phi_deta4";

    data->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta4_name_p6_z2,delta_phi_deta4_p6_z2);
    if (delta_phi_deta4_p6_z2 == 0) { cout << delta_phi_deta4_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta4_name_p8_4c,delta_phi_deta4_p8_4c);
    if (delta_phi_deta4_p8_4c == 0) { cout << delta_phi_deta4_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D("ak5PF_delta_phi_deta4","Corrected Data;#Delta#phi deta4;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta4, delta_phi_deta4_data, delta_phi_deta4_p6_z2, delta_phi_deta4_p8_4c, detail);
    plot_2histograms(delta_phi_deta4, "Data Corrected", delta_phi_deta4_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta4", "top_left", detail);


//compute corrected data for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_data = 0;
    TH1D *delta_phi_gap_p6_z2 = 0;
    TH1D *delta_phi_gap_p8_4c = 0;
    TString delta_phi_gap_name_p6_z2 = label_p6_z2 + "delta_phi_gap";
    TString delta_phi_gap_name_p8_4c = label_p8_4c + "delta_phi_gap";

    data->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_gap_name_p6_z2,delta_phi_gap_p6_z2);
    if (delta_phi_gap_p6_z2 == 0) { cout << delta_phi_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_gap_name_p8_4c,delta_phi_gap_p8_4c);
    if (delta_phi_gap_p8_4c == 0) { cout << delta_phi_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","Corrected Data;#Delta#phi gap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_gap, delta_phi_gap_data, delta_phi_gap_p6_z2, delta_phi_gap_p8_4c, detail);
    plot_2histograms(delta_phi_gap, "Data Corrected", delta_phi_gap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_gap", "top_left", detail);


//compute corrected data for delta phi deta1 gap distribution
   if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_data = 0;
    TH1D *delta_phi_deta1_gap_p6_z2 = 0;
    TH1D *delta_phi_deta1_gap_p8_4c = 0;
    TString delta_phi_deta1_gap_name_p6_z2 = label_p6_z2 + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_name_p8_4c = label_p8_4c + "delta_phi_deta1_gap";

    data->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta1_gap_name_p6_z2,delta_phi_deta1_gap_p6_z2);
    if (delta_phi_deta1_gap_p6_z2 == 0) { cout << delta_phi_deta1_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta1_gap_name_p8_4c,delta_phi_deta1_gap_p8_4c);
    if (delta_phi_deta1_gap_p8_4c == 0) { cout << delta_phi_deta1_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D("ak5PF_delta_phi_deta1_gap","Corrected Data;#Delta#phi deta1 gap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta1_gap, delta_phi_deta1_gap_data, delta_phi_deta1_gap_p6_z2, delta_phi_deta1_gap_p8_4c, detail);
    plot_2histograms(delta_phi_deta1_gap, "Data Corrected", delta_phi_deta1_gap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta1_gap", "top_left", detail);


//compute corrected data for delta phi deta2 gap distribution
   if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_data = 0;
    TH1D *delta_phi_deta2_gap_p6_z2 = 0;
    TH1D *delta_phi_deta2_gap_p8_4c = 0;
    TString delta_phi_deta2_gap_name_p6_z2 = label_p6_z2 + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_name_p8_4c = label_p8_4c + "delta_phi_deta2_gap";

    data->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta2_gap_name_p6_z2,delta_phi_deta2_gap_p6_z2);
    if (delta_phi_deta2_gap_p6_z2 == 0) { cout << delta_phi_deta2_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta2_gap_name_p8_4c,delta_phi_deta2_gap_p8_4c);
    if (delta_phi_deta2_gap_p8_4c == 0) { cout << delta_phi_deta2_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D("ak5PF_delta_phi_deta2_gap","Corrected Data;#Delta#phi deta2 gap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta2_gap, delta_phi_deta2_gap_data, delta_phi_deta2_gap_p6_z2, delta_phi_deta2_gap_p8_4c, detail);
    plot_2histograms(delta_phi_deta2_gap, "Data Corrected", delta_phi_deta2_gap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta2_gap", "top_left", detail);


//compute corrected data for delta phi deta3 gap distribution
   if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_data = 0;
    TH1D *delta_phi_deta3_gap_p6_z2 = 0;
    TH1D *delta_phi_deta3_gap_p8_4c = 0;
    TString delta_phi_deta3_gap_name_p6_z2 = label_p6_z2 + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_name_p8_4c = label_p8_4c + "delta_phi_deta3_gap";

    data->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta3_gap_name_p6_z2,delta_phi_deta3_gap_p6_z2);
    if (delta_phi_deta3_gap_p6_z2 == 0) { cout << delta_phi_deta3_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta3_gap_name_p8_4c,delta_phi_deta3_gap_p8_4c);
    if (delta_phi_deta3_gap_p8_4c == 0) { cout << delta_phi_deta3_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D("ak5PF_delta_phi_deta3_gap","Corrected Data;#Delta#phi deta3 gap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta3_gap, delta_phi_deta3_gap_data, delta_phi_deta3_gap_p6_z2, delta_phi_deta3_gap_p8_4c, detail);
    plot_2histograms(delta_phi_deta3_gap, "Data Corrected", delta_phi_deta3_gap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta3_gap", "top_left", detail);


//compute corrected data for delta phi deta4 gap distribution
   if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_data = 0;
    TH1D *delta_phi_deta4_gap_p6_z2 = 0;
    TH1D *delta_phi_deta4_gap_p8_4c = 0;
    TString delta_phi_deta4_gap_name_p6_z2 = label_p6_z2 + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_name_p8_4c = label_p8_4c + "delta_phi_deta4_gap";

    data->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta4_gap_name_p6_z2,delta_phi_deta4_gap_p6_z2);
    if (delta_phi_deta4_gap_p6_z2 == 0) { cout << delta_phi_deta4_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta4_gap_name_p8_4c,delta_phi_deta4_gap_p8_4c);
    if (delta_phi_deta4_gap_p8_4c == 0) { cout << delta_phi_deta4_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D("ak5PF_delta_phi_deta4_gap","Corrected Data;#Delta#phi deta4 gap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta4_gap, delta_phi_deta4_gap_data, delta_phi_deta4_gap_p6_z2, delta_phi_deta4_gap_p8_4c, detail);
    plot_2histograms(delta_phi_deta4_gap, "Data Corrected", delta_phi_deta4_gap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta4_gap", "top_left", detail);


//compute corrected data for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_data = 0;
    TH1D *delta_phi_nogap_p6_z2 = 0;
    TH1D *delta_phi_nogap_p8_4c = 0;
    TString delta_phi_nogap_name_p6_z2 = label_p6_z2 + "delta_phi_nogap";
    TString delta_phi_nogap_name_p8_4c = label_p8_4c + "delta_phi_nogap";

    data->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_nogap_name_p6_z2,delta_phi_nogap_p6_z2);
    if (delta_phi_nogap_p6_z2 == 0) { cout << delta_phi_nogap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_nogap_name_p8_4c,delta_phi_nogap_p8_4c);
    if (delta_phi_nogap_p8_4c == 0) { cout << delta_phi_nogap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","Corrected Data;#Delta#phi nogap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_nogap, delta_phi_nogap_data, delta_phi_nogap_p6_z2, delta_phi_nogap_p8_4c, detail);
    plot_2histograms(delta_phi_nogap, "Data Corrected", delta_phi_nogap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_nogap", "top_left", detail);


//compute corrected data for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_data = 0;
    TH1D *delta_phi_deta1_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta1_nogap_p8_4c = 0;
    TString delta_phi_deta1_nogap_name_p6_z2 = label_p6_z2 + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_name_p8_4c = label_p8_4c + "delta_phi_deta1_nogap";

    data->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta1_nogap_name_p6_z2,delta_phi_deta1_nogap_p6_z2);
    if (delta_phi_deta1_nogap_p6_z2 == 0) { cout << delta_phi_deta1_nogap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta1_nogap_name_p8_4c,delta_phi_deta1_nogap_p8_4c);
    if (delta_phi_deta1_nogap_p8_4c == 0) { cout << delta_phi_deta1_nogap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D("ak5PF_delta_phi_deta1_nogap","Corrected Data;#Delta#phi deta1 nogap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta1_nogap, delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_p6_z2, delta_phi_deta1_nogap_p8_4c, detail);
    plot_2histograms(delta_phi_deta1_nogap, "Data Corrected", delta_phi_deta1_nogap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta1_nogap", "top_left", detail);


//compute corrected data for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_data = 0;
    TH1D *delta_phi_deta2_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta2_nogap_p8_4c = 0;
    TString delta_phi_deta2_nogap_name_p6_z2 = label_p6_z2 + "delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_name_p8_4c = label_p8_4c + "delta_phi_deta2_nogap";

    data->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta2_nogap_name_p6_z2,delta_phi_deta2_nogap_p6_z2);
    if (delta_phi_deta2_nogap_p6_z2 == 0) { cout << delta_phi_deta2_nogap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta2_nogap_name_p8_4c,delta_phi_deta2_nogap_p8_4c);
    if (delta_phi_deta2_nogap_p8_4c == 0) { cout << delta_phi_deta2_nogap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D("ak5PF_delta_phi_deta2_nogap","Corrected Data;#Delta#phi deta2 nogap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta2_nogap, delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_p6_z2, delta_phi_deta2_nogap_p8_4c, detail);
    plot_2histograms(delta_phi_deta2_nogap, "Data Corrected", delta_phi_deta2_nogap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta2_nogap", "top_left", detail);


//compute corrected data for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_data = 0;
    TH1D *delta_phi_deta3_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta3_nogap_p8_4c = 0;
    TString delta_phi_deta3_nogap_name_p6_z2 = label_p6_z2 + "delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_name_p8_4c = label_p8_4c + "delta_phi_deta3_nogap";

    data->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta3_nogap_name_p6_z2,delta_phi_deta3_nogap_p6_z2);
    if (delta_phi_deta3_nogap_p6_z2 == 0) { cout << delta_phi_deta3_nogap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta3_nogap_name_p8_4c,delta_phi_deta3_nogap_p8_4c);
    if (delta_phi_deta3_nogap_p8_4c == 0) { cout << delta_phi_deta3_nogap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D("ak5PF_delta_phi_deta3_nogap","Corrected Data;#Delta#phi deta3 nogap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta3_nogap, delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_p6_z2, delta_phi_deta3_nogap_p8_4c, detail);
    plot_2histograms(delta_phi_deta3_nogap, "Data Corrected", delta_phi_deta3_nogap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta3_nogap", "top_left", detail);


//compute corrected data for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_data = 0;
    TH1D *delta_phi_deta4_nogap_p6_z2 = 0;
    TH1D *delta_phi_deta4_nogap_p8_4c = 0;
    TString delta_phi_deta4_nogap_name_p6_z2 = label_p6_z2 + "delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_name_p8_4c = label_p8_4c + "delta_phi_deta4_nogap";

    data->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    p6_z2->GetObject(delta_phi_deta4_nogap_name_p6_z2,delta_phi_deta4_nogap_p6_z2);
    if (delta_phi_deta4_nogap_p6_z2 == 0) { cout << delta_phi_deta4_nogap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_phi_deta4_nogap_name_p8_4c,delta_phi_deta4_nogap_p8_4c);
    if (delta_phi_deta4_nogap_p8_4c == 0) { cout << delta_phi_deta4_nogap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D("ak5PF_delta_phi_deta4_nogap","Corrected Data;#Delta#phi deta4 nogap;#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);

    correct_histogram(delta_phi_deta4_nogap, delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_p6_z2, delta_phi_deta4_nogap_p8_4c, detail);
    plot_2histograms(delta_phi_deta4_nogap, "Data Corrected", delta_phi_deta4_nogap_data, "Data Uncorrected", output_path_plots, plot_label + "delta_phi_deta4_nogap", "top_left", detail);


//compute corrected data for leading pt inside gap distribution
    if (detail) { cout<<"Leading pt inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_data = 0;
    TH1D *leading_pt_inside_gap_p6_z2 = 0;
    TH1D *leading_pt_inside_gap_p8_4c = 0;
    TString leading_pt_inside_gap_name_p6_z2 = label_p6_z2 + "leading_pt_inside_gap";
    TString leading_pt_inside_gap_name_p8_4c = label_p8_4c + "leading_pt_inside_gap";

    data->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    p6_z2->GetObject(leading_pt_inside_gap_name_p6_z2,leading_pt_inside_gap_p6_z2);
    if (leading_pt_inside_gap_p6_z2 == 0) { cout << leading_pt_inside_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(leading_pt_inside_gap_name_p8_4c,leading_pt_inside_gap_p8_4c);
    if (leading_pt_inside_gap_p8_4c == 0) { cout << leading_pt_inside_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D("ak5PF_leading_pt_inside_gap","Corrected Data;p_{T} inside gap;#frac{d#sigma}{d#p_{T}}", in_nbins, in_bins);

    correct_histogram(leading_pt_inside_gap, leading_pt_inside_gap_data, leading_pt_inside_gap_p6_z2, leading_pt_inside_gap_p8_4c, detail);
    plot_2histograms(leading_pt_inside_gap, "Data Corrected", leading_pt_inside_gap_data, "Data Uncorrected", output_path_plots, plot_label + "leading_pt_inside_gap", "top_left", detail);


//compute corrected data for leading eta star inside gap distribution
    if (detail) { cout<<"Leading eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_data = 0;
    TH1D *leading_eta_star_inside_gap_p6_z2 = 0;
    TH1D *leading_eta_star_inside_gap_p8_4c = 0;
    TString leading_eta_star_inside_gap_name_p6_z2 = label_p6_z2 + "leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_name_p8_4c = label_p8_4c + "leading_eta_star_inside_gap";

    data->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    p6_z2->GetObject(leading_eta_star_inside_gap_name_p6_z2,leading_eta_star_inside_gap_p6_z2);
    if (leading_eta_star_inside_gap_p6_z2 == 0) { cout << leading_eta_star_inside_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(leading_eta_star_inside_gap_name_p8_4c,leading_eta_star_inside_gap_p8_4c);
    if (leading_eta_star_inside_gap_p8_4c == 0) { cout << leading_eta_star_inside_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D("ak5PF_leading_eta_star_inside_gap","Corrected Data;#eta* inside gap;#frac{d#sigma}{d#eta*}", etastar_nbins, etastar_bins);

    correct_histogram(leading_eta_star_inside_gap, leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_p6_z2, leading_eta_star_inside_gap_p8_4c, detail);
    plot_2histograms(leading_eta_star_inside_gap, "Data Corrected", leading_eta_star_inside_gap_data, "Data Uncorrected", output_path_plots, plot_label + "leading_eta_star_inside_gap", "top_left", detail);


//compute corrected data for leading pt outside gap distribution
    if (detail) { cout<<"Leading pt outside gap"<<endl; }

    TH1D *delta_eta_outside_gap_data = 0;
    TH1D *delta_eta_outside_gap_p6_z2 = 0;
    TH1D *delta_eta_outside_gap_p8_4c = 0;
    TString delta_eta_outside_gap_name_p6_z2 = label_p6_z2 + "delta_eta_outside_gap";
    TString delta_eta_outside_gap_name_p8_4c = label_p8_4c + "delta_eta_outside_gap";

    data->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    p6_z2->GetObject(delta_eta_outside_gap_name_p6_z2,delta_eta_outside_gap_p6_z2);
    if (delta_eta_outside_gap_p6_z2 == 0) { cout << delta_eta_outside_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(delta_eta_outside_gap_name_p8_4c,delta_eta_outside_gap_p8_4c);
    if (delta_eta_outside_gap_p8_4c == 0) { cout << delta_eta_outside_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D("ak5PF_delta_eta_outside_gap","Corrected Data;#Delta#eta outside gap;#frac{d#sigma}{d#Delta#eta}", deta_out_nbins, deta_out_bins);

    correct_histogram(delta_eta_outside_gap, delta_eta_outside_gap_data, delta_eta_outside_gap_p6_z2, delta_eta_outside_gap_p8_4c, detail);
    plot_2histograms(delta_eta_outside_gap, "Data Corrected", delta_eta_outside_gap, "Data Uncorrected", output_path_plots, plot_label + "delta_eta_outside_gap", "top_left", detail);


//compute corrected data for leading pt outside gap distribution
    if (detail) { cout<<"Leading pt outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_data = 0;
    TH1D *leading_pt_outside_gap_p6_z2 = 0;
    TH1D *leading_pt_outside_gap_p8_4c = 0;
    TString leading_pt_outside_gap_name_p6_z2 = label_p6_z2 + "leading_pt_outside_gap";
    TString leading_pt_outside_gap_name_p8_4c = label_p8_4c + "leading_pt_outside_gap";

    data->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    p6_z2->GetObject(leading_pt_outside_gap_name_p6_z2,leading_pt_outside_gap_p6_z2);
    if (leading_pt_outside_gap_p6_z2 == 0) { cout << leading_pt_outside_gap_name_p6_z2 << " not found!" << endl; return; }
    p8_4c->GetObject(leading_pt_outside_gap_name_p8_4c,leading_pt_outside_gap_p8_4c);
    if (leading_pt_outside_gap_p8_4c == 0) { cout << leading_pt_outside_gap_name_p8_4c << " not found!" << endl; return; }

    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D("ak5PF_leading_pt_outside_gap","Corrected Data;p_{T} outside gap;#frac{d#sigma}{d#p_{T}}", out_nbins, out_bins);

    correct_histogram(leading_pt_outside_gap, leading_pt_outside_gap_data, leading_pt_outside_gap_p6_z2, leading_pt_outside_gap_p8_4c, detail);
    plot_2histograms(leading_pt_outside_gap, "Data Corrected", leading_pt_outside_gap_data, "Data Uncorrected", output_path_plots, plot_label + "leading_pt_outside_gap", "top_left", detail);


//Opening the output root file
    if (detail) { cout<<"Creating " << path_corrected << "..."<<endl; }
    TFile *corrected = TFile::Open( path_corrected.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file..."<<endl; }
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
    leading_pt_inside_gap->Write();
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();
    leading_pt_outside_gap->Write();
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data->Close();
    p6_z2->Close();
    p8_4c->Close();
    corrected->Close();

//deleting all histograms to avoid memory leak - causes a memory leak!
/*    delete(delta_phi);
    delete(delta_phi_deta1);
    delete(delta_phi_deta2);
    delete(delta_phi_deta3);
    delete(delta_phi_deta4);
    delete(delta_phi_gap);
    delete(delta_phi_deta1_gap);
    delete(delta_phi_deta2_gap);
    delete(delta_phi_deta3_gap);
    delete(delta_phi_deta4_gap);
    delete(delta_phi_nogap);
    delete(delta_phi_deta1_nogap);
    delete(delta_phi_deta2_nogap);
    delete(delta_phi_deta3_nogap);
    delete(delta_phi_deta4_nogap);
    delete(leading_pt_inside_gap);
    delete(leading_pt_outside_gap);
    delete(delta_eta_outside_gap);
    delete(leading_eta_star_inside_gap); */

    if (detail) { cout<<"Done!"<<endl; }

}
