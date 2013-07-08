// Pedro Cipriano, Mar 2013
// DESY, CMS
// Last Update: 13 May 2013
//
// create_ratios()

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

void create_ratios(string path_source, string path_reference, string prefix_source, string prefix_reference, string path_output, string output_path_plots, string prefix_output, bool detail = false, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Create Ratios Configuration"<<endl; }
    if (detail) { cout << "Input path for source:    " << path_source << endl; }
    if (detail) { cout << "Input path for reference: " << path_reference << endl; }
    if (detail) { cout << "Prefix for source:        " << prefix_source << endl; }
    if (detail) { cout << "Prefix for reference:     " << prefix_reference << endl; }
    if (detail) { cout << "Output path:              " << path_output << endl; }
    if (detail) { cout << "Output Path Plots:        " << output_path_plots << endl; }
    if (detail) { cout << "Prefix Output:            " << prefix_output << endl; }
    if (detail) { cout << "Detail Level:             " << detail << endl; }
    if (detail) { cout << "Test Mode:                " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *source = new TFile( path_source.c_str() );
    TFile *reference = new TFile( path_reference.c_str() );

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
double deta_out_bins[7] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5};

//compute ratio for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1 *delta_phi_source = 0;
    TH1 *delta_phi_reference = 0;
    TString delta_phi_name_source;
    TString delta_phi_name_reference;

    if (prefix_source != "") { delta_phi_name_source = prefix_source + "delta_phi"; }
    else { delta_phi_name_source = "d01-x01-y01"; }
    delta_phi_name_reference = prefix_reference + "delta_phi"; 

    source->GetObject(delta_phi_name_source,delta_phi_source);
    if (delta_phi_source == 0) { cout << delta_phi_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_name_reference,delta_phi_reference);
    if (delta_phi_reference == 0) { cout << delta_phi_name_reference + " reference not found!" << endl; return; }

    TH1D *delta_phi_ratio;
    delta_phi_ratio =  new TH1D(delta_phi_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_ratio->Divide(delta_phi_source, delta_phi_reference, 1., 1., "");
    plot_2histograms(delta_phi_source, "Source", delta_phi_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi", "top_left", true, detail);
    plot_histogram(delta_phi_ratio, output_path_plots, prefix_output + "delta_phi", "ratio", "top_left");


//compute ratio for delta phi nogap distribution
    if (detail) { cout<<"Delta phi noGap"<<endl; }

    TH1 *delta_phi_nogap_source = 0;
    TH1 *delta_phi_nogap_reference = 0;
    TString delta_phi_nogap_name_source;
    TString delta_phi_nogap_name_reference;

    if (prefix_source != "") { delta_phi_nogap_name_source = prefix_source + "delta_phi_nogap"; }
    else { delta_phi_nogap_name_source = "d02-x01-y01"; }
    delta_phi_nogap_name_reference = prefix_reference + "delta_phi_nogap";

    source->GetObject(delta_phi_nogap_name_source,delta_phi_nogap_source);
    if (delta_phi_nogap_source == 0) { cout << delta_phi_nogap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_nogap_name_reference,delta_phi_nogap_reference);
    if (delta_phi_nogap_reference == 0) { cout << delta_phi_nogap_name_reference + " reference not found!" << endl; return; }

    TH1D *delta_phi_nogap_ratio;
    delta_phi_nogap_ratio =  new TH1D(delta_phi_nogap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_nogap_ratio->Divide(delta_phi_nogap_source, delta_phi_nogap_reference, 1., 1., "");
    plot_2histograms(delta_phi_nogap_source, "Source", delta_phi_nogap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_nogap_ratio, output_path_plots, prefix_output + "delta_phi_nogap", "ratio", "top_left");


//compute ratio for delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1 *delta_phi_gap_source = 0;
    TH1 *delta_phi_gap_reference = 0;
    TString delta_phi_gap_name_source;
    TString delta_phi_gap_name_reference;

    if (prefix_source != "") { delta_phi_gap_name_source = prefix_source + "delta_phi_gap"; }
    else { delta_phi_gap_name_source = "d03-x01-y01"; }
    delta_phi_gap_name_reference = prefix_reference + "delta_phi_gap";

    source->GetObject(delta_phi_gap_name_source,delta_phi_gap_source);
    if (delta_phi_gap_source == 0) { cout << delta_phi_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_gap_name_reference,delta_phi_gap_reference);
    if (delta_phi_gap_reference == 0) { cout << delta_phi_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_gap_ratio;
    delta_phi_gap_ratio =  new TH1D(delta_phi_gap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_gap_ratio->Divide(delta_phi_gap_source, delta_phi_gap_reference, 1., 1., "");
    plot_2histograms(delta_phi_gap_source, "Source", delta_phi_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_gap", "top_left", true, detail);
    plot_histogram(delta_phi_gap_ratio, output_path_plots, prefix_output + "delta_phi_gap", "ratio", "top_left");


//compute ratio for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1 *delta_phi_deta1_source = 0;
    TH1 *delta_phi_deta1_reference = 0;
    TString delta_phi_deta1_name_source;
    TString delta_phi_deta1_name_reference;

    if (prefix_source != "") { delta_phi_deta1_name_source = prefix_source + "delta_phi_deta1"; }
    else { delta_phi_deta1_name_source = "d01-x02-y01"; }
    delta_phi_deta1_name_reference = prefix_reference + "delta_phi_deta1";

    source->GetObject(delta_phi_deta1_name_source,delta_phi_deta1_source);
    if (delta_phi_deta1_source == 0) { cout << delta_phi_deta1_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta1_name_reference,delta_phi_deta1_reference);
    if (delta_phi_deta1_reference == 0) { cout << "ak5PF_delta_phi_deta1 reference not found!" << endl; return; }

    TH1D *delta_phi_deta1_ratio;
    delta_phi_deta1_ratio =  new TH1D(delta_phi_deta1_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta1_ratio->Divide(delta_phi_deta1_source, delta_phi_deta1_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta1_source, "Source", delta_phi_deta1_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta1", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_ratio, output_path_plots, prefix_output + "delta_phi_deta1", "ratio", "top_left");


//compute ratio for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1 *delta_phi_deta2_source = 0;
    TH1 *delta_phi_deta2_reference = 0;
    TString delta_phi_deta2_name_source;
    TString delta_phi_deta2_name_reference;

    if (prefix_source != "") { delta_phi_deta2_name_source = prefix_source + "delta_phi_deta2"; }
    else { delta_phi_deta2_name_source = "d01-x02-y02"; }
    delta_phi_deta2_name_reference = prefix_reference + "delta_phi_deta2";

    source->GetObject(delta_phi_deta2_name_source,delta_phi_deta2_source);
    if (delta_phi_deta2_source == 0) { cout << delta_phi_deta2_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta2_name_reference,delta_phi_deta2_reference);
    if (delta_phi_deta2_reference == 0) { cout << delta_phi_deta2_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta2_ratio;
    delta_phi_deta2_ratio =  new TH1D(delta_phi_deta2_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta2_ratio->Divide(delta_phi_deta2_source, delta_phi_deta2_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta2_source, "Source", delta_phi_deta2_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta2", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_ratio, output_path_plots, prefix_output + "delta_phi_deta2", "ratio","top_left");

//compute ratio for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1 *delta_phi_deta3_source = 0;
    TH1 *delta_phi_deta3_reference = 0;
    TString delta_phi_deta3_name_source;
    TString delta_phi_deta3_name_reference;

    if (prefix_source != "") { delta_phi_deta3_name_source = prefix_source + "delta_phi_deta3"; }
    else { delta_phi_deta3_name_source = "d01-x02-y03"; }
    delta_phi_deta3_name_reference = prefix_reference + "delta_phi_deta3";

    source->GetObject(delta_phi_deta3_name_source,delta_phi_deta3_source);
    if (delta_phi_deta3_source == 0) { cout << delta_phi_deta3_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta3_name_reference,delta_phi_deta3_reference);
    if (delta_phi_deta3_reference == 0) { cout << delta_phi_deta3_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta3_ratio;
    delta_phi_deta3_ratio =  new TH1D(delta_phi_deta3_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta3_ratio->Divide(delta_phi_deta3_source, delta_phi_deta3_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta3_source, "Source", delta_phi_deta3_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta3", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_ratio, output_path_plots, prefix_output + "delta_phi_deta3", "ratio", "top_left");


//compute ratio for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1 *delta_phi_deta4_source = 0;
    TH1 *delta_phi_deta4_reference = 0;
    TString delta_phi_deta4_name_source;
    TString delta_phi_deta4_name_reference;

    if (prefix_source != "") { delta_phi_deta4_name_source = prefix_source + "delta_phi_deta4"; }
    else { delta_phi_deta4_name_source = "d01-x02-y04"; }
    delta_phi_deta4_name_reference = prefix_reference + "delta_phi_deta4";

    source->GetObject(delta_phi_deta4_name_source,delta_phi_deta4_source);
    if (delta_phi_deta4_source == 0) { cout << delta_phi_deta4_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta4_name_reference,delta_phi_deta4_reference);
    if (delta_phi_deta4_reference == 0) { cout << delta_phi_deta4_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta4_ratio;
    delta_phi_deta4_ratio =  new TH1D(delta_phi_deta4_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta4_ratio->Divide(delta_phi_deta4_source, delta_phi_deta4_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta4_source, "Source", delta_phi_deta4_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta4", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_ratio, output_path_plots, prefix_output + "delta_phi_deta4", "ratio", "top_left");


//compute ratio for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1 *delta_phi_deta1_nogap_source = 0;
    TH1 *delta_phi_deta1_nogap_reference = 0;
    TString delta_phi_deta1_nogap_name_source;
    TString delta_phi_deta1_nogap_name_reference;

    if (prefix_source != "") { delta_phi_deta1_nogap_name_source = prefix_source + "delta_phi_deta1_nogap"; }
    else { delta_phi_deta1_nogap_name_source = "d02-x02-y01"; }
    delta_phi_deta1_nogap_name_reference = prefix_reference + "delta_phi_deta1_nogap";

    source->GetObject(delta_phi_deta1_nogap_name_source,delta_phi_deta1_nogap_source);
    if (delta_phi_deta1_nogap_source == 0) { cout << delta_phi_deta1_nogap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta1_nogap_name_reference,delta_phi_deta1_nogap_reference);
    if (delta_phi_deta1_nogap_reference == 0) { cout << delta_phi_deta1_nogap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta1_nogap_ratio;
    delta_phi_deta1_nogap_ratio =  new TH1D(delta_phi_deta1_nogap_name_source,"Ratio plot;#Delta#phi [rad];MC/ DATA", dphi_nbins, dphi_bins);

    delta_phi_deta1_nogap_ratio->Divide(delta_phi_deta1_nogap_source, delta_phi_deta1_nogap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta1_nogap_source, "Source", delta_phi_deta1_nogap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta1_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_nogap_ratio, output_path_plots, prefix_output + "delta_phi_deta1_nogap", "ratio", "top_left");


//compute ratio for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1 *delta_phi_deta2_nogap_source = 0;
    TH1 *delta_phi_deta2_nogap_reference = 0;
    TString delta_phi_deta2_nogap_name_source;
    TString delta_phi_deta2_nogap_name_reference;

    if (prefix_source != "") { delta_phi_deta2_nogap_name_source = prefix_source + "delta_phi_deta2_nogap"; }
    else { delta_phi_deta2_nogap_name_source = "d02-x02-y02"; }
    delta_phi_deta2_nogap_name_reference = prefix_reference + "delta_phi_deta2_nogap";

    source->GetObject(delta_phi_deta2_nogap_name_source,delta_phi_deta2_nogap_source);
    if (delta_phi_deta2_nogap_source == 0) { cout << delta_phi_deta2_nogap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta2_nogap_name_reference,delta_phi_deta2_nogap_reference);
    if (delta_phi_deta2_nogap_reference == 0) { cout << delta_phi_deta2_nogap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta2_nogap_ratio;
    delta_phi_deta2_nogap_ratio =  new TH1D(delta_phi_deta2_nogap_name_source,"Ratio plot;#Delta#phi [rad];MC/ DATA", dphi_nbins, dphi_bins);

    delta_phi_deta2_nogap_ratio->Divide(delta_phi_deta2_nogap_source, delta_phi_deta2_nogap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta2_nogap_source, "Source", delta_phi_deta2_nogap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta2_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_nogap_ratio, output_path_plots, prefix_output + "delta_phi_deta2_nogap", "ratio", "top_left");


//compute ratio for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1 *delta_phi_deta3_nogap_source = 0;
    TH1 *delta_phi_deta3_nogap_reference = 0;
    TString delta_phi_deta3_nogap_name_source;
    TString delta_phi_deta3_nogap_name_reference;

    if (prefix_source != "") { delta_phi_deta3_nogap_name_source = prefix_source + "delta_phi_deta3_nogap"; }
    else { delta_phi_deta3_nogap_name_source = "d02-x02-y03"; }
    delta_phi_deta3_nogap_name_reference = prefix_reference + "delta_phi_deta3_nogap";

    source->GetObject(delta_phi_deta3_nogap_name_source,delta_phi_deta3_nogap_source);
    if (delta_phi_deta3_nogap_source == 0) { cout << delta_phi_deta3_nogap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta3_nogap_name_reference,delta_phi_deta3_nogap_reference);
    if (delta_phi_deta3_nogap_reference == 0) { cout << delta_phi_deta3_nogap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta3_nogap_ratio;
    delta_phi_deta3_nogap_ratio =  new TH1D(delta_phi_deta3_nogap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta3_nogap_ratio->Divide(delta_phi_deta3_nogap_source, delta_phi_deta3_nogap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta3_nogap_source, "Source", delta_phi_deta3_nogap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta3_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_nogap_ratio, output_path_plots, prefix_output + "delta_phi_deta3_nogap", "ratio", "top_left");


//compute ratio for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1 *delta_phi_deta4_nogap_source = 0;
    TH1 *delta_phi_deta4_nogap_reference = 0;
    TString delta_phi_deta4_nogap_name_source;
    TString delta_phi_deta4_nogap_name_reference;

    if (prefix_source != "") { delta_phi_deta4_nogap_name_source = prefix_source + "delta_phi_deta4_nogap"; }
    else { delta_phi_deta4_nogap_name_source = "d02-x02-y04"; }
    delta_phi_deta4_nogap_name_reference = prefix_reference + "delta_phi_deta4_nogap";

    source->GetObject(delta_phi_deta4_nogap_name_source,delta_phi_deta4_nogap_source);
    if (delta_phi_deta4_nogap_source == 0) { cout << delta_phi_deta4_nogap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta4_nogap_name_reference,delta_phi_deta4_nogap_reference);
    if (delta_phi_deta4_nogap_reference == 0) { cout << delta_phi_deta4_nogap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta4_nogap_ratio;
    delta_phi_deta4_nogap_ratio =  new TH1D(delta_phi_deta4_nogap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta4_nogap_ratio->Divide(delta_phi_deta4_nogap_source, delta_phi_deta4_nogap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta4_nogap_source, "Source", delta_phi_deta4_nogap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta4_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_nogap_ratio, output_path_plots, prefix_output + "delta_phi_deta4_nogap", "ratio", "top_left");


//compute ratio for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1 *delta_phi_deta1_gap_source = 0;
    TH1 *delta_phi_deta1_gap_reference = 0;
    TString delta_phi_deta1_gap_name_source;
    TString delta_phi_deta1_gap_name_reference;

    if (prefix_source != "") { delta_phi_deta1_gap_name_source = prefix_source + "delta_phi_deta1_gap"; }
    else { delta_phi_deta1_gap_name_source = "d03-x02-y01"; }
    delta_phi_deta1_gap_name_reference = prefix_reference + "delta_phi_deta1_gap";

    source->GetObject(delta_phi_deta1_gap_name_source,delta_phi_deta1_gap_source);
    if (delta_phi_deta1_gap_source == 0) { cout << delta_phi_deta1_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta1_gap_name_reference,delta_phi_deta1_gap_reference);
    if (delta_phi_deta1_gap_reference == 0) { cout << delta_phi_deta1_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta1_gap_ratio;
    delta_phi_deta1_gap_ratio =  new TH1D(delta_phi_deta1_gap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta1_gap_ratio->Divide(delta_phi_deta1_gap_source, delta_phi_deta1_gap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta1_gap_source, "Source", delta_phi_deta1_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta1_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_gap_ratio, output_path_plots, prefix_output + "delta_phi_deta1_gap", "ratio", "top_left");


//compute ratio for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1 *delta_phi_deta2_gap_source = 0;
    TH1 *delta_phi_deta2_gap_reference = 0;
    TString delta_phi_deta2_gap_name_source;
    TString delta_phi_deta2_gap_name_reference;

    if (prefix_source != "") { delta_phi_deta2_gap_name_source = prefix_source + "delta_phi_deta2_gap"; }
    else { delta_phi_deta2_gap_name_source = "d03-x02-y02"; }
    delta_phi_deta2_gap_name_reference = prefix_reference + "delta_phi_deta2_gap";

    source->GetObject(delta_phi_deta2_gap_name_source,delta_phi_deta2_gap_source);
    if (delta_phi_deta2_gap_source == 0) { cout << delta_phi_deta2_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta2_gap_name_reference,delta_phi_deta2_gap_reference);
    if (delta_phi_deta2_gap_reference == 0) { cout << delta_phi_deta2_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta2_gap_ratio;
    delta_phi_deta2_gap_ratio =  new TH1D(delta_phi_deta2_gap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta2_gap_ratio->Divide(delta_phi_deta2_gap_source, delta_phi_deta2_gap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta2_gap_source, "Source", delta_phi_deta2_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta2_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_gap_ratio, output_path_plots, prefix_output + "delta_phi_deta2_gap", "ratio", "top_left");


//compute ratio for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1 *delta_phi_deta3_gap_source = 0;
    TH1 *delta_phi_deta3_gap_reference = 0;
    TString delta_phi_deta3_gap_name_source;
    TString delta_phi_deta3_gap_name_reference;

    if (prefix_source != "") { delta_phi_deta3_gap_name_source = prefix_source + "delta_phi_deta3_gap"; }
    else { delta_phi_deta3_gap_name_source = "d03-x02-y03"; }
    delta_phi_deta3_gap_name_reference = prefix_reference + "delta_phi_deta3_gap";

    source->GetObject(delta_phi_deta3_gap_name_source,delta_phi_deta3_gap_source);
    if (delta_phi_deta3_gap_source == 0) { cout << delta_phi_deta3_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta3_gap_name_reference,delta_phi_deta3_gap_reference);
    if (delta_phi_deta3_gap_reference == 0) { cout << delta_phi_deta3_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta3_gap_ratio;
    delta_phi_deta3_gap_ratio =  new TH1D(delta_phi_deta3_gap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta3_gap_ratio->Divide(delta_phi_deta3_gap_source, delta_phi_deta3_gap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta3_gap_source, "Source", delta_phi_deta3_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta3_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_gap_ratio, output_path_plots, prefix_output + "delta_phi_deta3_gap", "ratio", "top_left");


//compute ratio for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1 *delta_phi_deta4_gap_source = 0;
    TH1 *delta_phi_deta4_gap_reference = 0;
    TString delta_phi_deta4_gap_name_source;
    TString delta_phi_deta4_gap_name_reference;

    if (prefix_source != "") { delta_phi_deta4_gap_name_source = prefix_source + "delta_phi_deta4_gap"; }
    else { delta_phi_deta4_gap_name_source = "d03-x02-y04"; }
    delta_phi_deta4_gap_name_reference = prefix_reference + "delta_phi_deta4_gap";

    source->GetObject(delta_phi_deta4_gap_name_source,delta_phi_deta4_gap_source);
    if (delta_phi_deta4_gap_source == 0) { cout << delta_phi_deta4_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_phi_deta4_gap_name_reference,delta_phi_deta4_gap_reference);
    if (delta_phi_deta4_gap_reference == 0) { cout << delta_phi_deta4_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_phi_deta4_gap_ratio;
    delta_phi_deta4_gap_ratio =  new TH1D(delta_phi_deta4_gap_name_source,"Ratio plot;#Delta#phi [rad];MC/DATA", dphi_nbins, dphi_bins);

    delta_phi_deta4_gap_ratio->Divide(delta_phi_deta4_gap_source, delta_phi_deta4_gap_reference, 1., 1., "");
    plot_2histograms(delta_phi_deta4_gap_source, "Source", delta_phi_deta4_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_phi_deta4_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_gap_ratio, output_path_plots, prefix_output + "delta_phi_deta4_gap", "ratio", "top_left");


//compute ratio for leading pt inside gap distribution
    if (detail) { cout<<"Leading pT inside Gap"<<endl; }

    TH1 *leading_pt_inside_gap_source = 0;
    TH1 *leading_pt_inside_gap_reference = 0;
    TString leading_pt_inside_gap_name_source;
    TString leading_pt_inside_gap_name_reference;

    if (prefix_source != "") { leading_pt_inside_gap_name_source = prefix_source + "leading_pt_inside_gap"; }
    else { leading_pt_inside_gap_name_source = "d04-x01-y01"; }
    leading_pt_inside_gap_name_reference = prefix_reference + "leading_pt_inside_gap";

    source->GetObject(leading_pt_inside_gap_name_source,leading_pt_inside_gap_source);
    if (leading_pt_inside_gap_source == 0) { cout << leading_pt_inside_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(leading_pt_inside_gap_name_reference,leading_pt_inside_gap_reference);
    if (leading_pt_inside_gap_reference == 0) { cout << leading_pt_inside_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *leading_pt_inside_gap_ratio;
    leading_pt_inside_gap_ratio =  new TH1D(leading_pt_inside_gap_name_source,"Ratio plot;p_{T} [GeV];MC/DATA", in_nbins, in_bins);

    leading_pt_inside_gap_ratio->Divide(leading_pt_inside_gap_source, leading_pt_inside_gap_reference, 1., 1., "");
    plot_2histograms(leading_pt_inside_gap_source, "Source", leading_pt_inside_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "leading_pt_inside_gap", "top_left", true, detail);
    plot_histogram(leading_pt_inside_gap_ratio, output_path_plots, prefix_output + "leading_pt_inside_gap", "ratio", "top_right");


//compute ratio for leading eta* inside gap distribution
    if (detail) { cout<<"Leading Eta* inside Gap"<<endl; }

    TH1 *leading_eta_star_inside_gap_source = 0;
    TH1 *leading_eta_star_inside_gap_reference = 0;
    TString leading_eta_star_inside_gap_name_source;
    TString leading_eta_star_inside_gap_name_reference;

    if (prefix_source != "") { leading_eta_star_inside_gap_name_source = prefix_source + "leading_eta_star_inside_gap"; }
    else { leading_eta_star_inside_gap_name_source = "d06-x01-y01"; }
    leading_eta_star_inside_gap_name_reference = prefix_reference + "leading_eta_star_inside_gap";

    source->GetObject(leading_eta_star_inside_gap_name_source,leading_eta_star_inside_gap_source);
    if (leading_eta_star_inside_gap_source == 0) { cout << leading_eta_star_inside_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(leading_eta_star_inside_gap_name_reference,leading_eta_star_inside_gap_reference);
    if (leading_eta_star_inside_gap_reference == 0) { cout << leading_eta_star_inside_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *leading_eta_star_inside_gap_ratio;
    leading_eta_star_inside_gap_ratio =  new TH1D(leading_eta_star_inside_gap_name_source,"Ratio plot;#eta*;MC/DATA", etastar_nbins, etastar_bins);

    leading_eta_star_inside_gap_ratio->Divide(leading_eta_star_inside_gap_source, leading_eta_star_inside_gap_reference, 1., 1., "");
    plot_2histograms(leading_eta_star_inside_gap_source, "Source", leading_eta_star_inside_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "leading_eta_star_inside_gap", "top_left", true, detail);
    plot_histogram(leading_eta_star_inside_gap_ratio, output_path_plots, prefix_output + "leading_eta_star_inside_gap", "ratio","top_right");


//compute ratio for leading pt outside gap distribution
    if (detail) { cout<<"Leading pT outside Gap"<<endl; }

    TH1 *leading_pt_outside_gap_source = 0;
    TH1 *leading_pt_outside_gap_reference = 0;
    TString leading_pt_outside_gap_name_source;
    TString leading_pt_outside_gap_name_reference;

    if (prefix_source != "") { leading_pt_outside_gap_name_source = prefix_source + "leading_pt_outside_gap"; }
    else { leading_pt_outside_gap_name_source = "d05-x01-y01"; }
    leading_pt_outside_gap_name_reference = prefix_reference + "leading_pt_outside_gap";

    source->GetObject(leading_pt_outside_gap_name_source,leading_pt_outside_gap_source);
    if (leading_pt_outside_gap_source == 0) { cout << leading_pt_outside_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(leading_pt_outside_gap_name_reference,leading_pt_outside_gap_reference);
    if (leading_pt_outside_gap_reference == 0) { cout << leading_pt_outside_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *leading_pt_outside_gap_ratio;
    leading_pt_outside_gap_ratio =  new TH1D(leading_pt_outside_gap_name_source,"Ratio plot;p_{T} [GeV];MC/DATA", out_nbins, out_bins);

    leading_pt_outside_gap_ratio->Divide(leading_pt_outside_gap_source, leading_pt_outside_gap_reference, 1., 1., "");
    plot_2histograms(leading_pt_outside_gap_source, "Source", leading_pt_outside_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "leading_pt_outside_gap", "top_left", true, detail);
    plot_histogram(leading_pt_outside_gap_ratio, output_path_plots, prefix_output + "leading_pt_outside_gap", "ratio","top_right");


//compute ratio for delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta outside Gap"<<endl; }

    TH1 *delta_eta_outside_gap_source = 0;
    TH1 *delta_eta_outside_gap_reference = 0;
    TString delta_eta_outside_gap_name_source;
    TString delta_eta_outside_gap_name_reference;

    if (prefix_source != "") { delta_eta_outside_gap_name_source = prefix_source + "delta_eta_outside_gap"; }
    else { delta_eta_outside_gap_name_source = "d07-x01-y01"; }
    delta_eta_outside_gap_name_reference = prefix_reference + "delta_eta_outside_gap"; 

    source->GetObject(delta_eta_outside_gap_name_source,delta_eta_outside_gap_source);
    if (delta_eta_outside_gap_source == 0) { cout << delta_eta_outside_gap_name_source << " source not found!" << endl; return; }
    reference->GetObject(delta_eta_outside_gap_name_reference,delta_eta_outside_gap_reference);
    if (delta_eta_outside_gap_reference == 0) { cout << delta_eta_outside_gap_name_reference << " reference not found!" << endl; return; }

    TH1D *delta_eta_outside_gap_ratio;
    delta_eta_outside_gap_ratio =  new TH1D(delta_eta_outside_gap_name_source,"Ratio plot;#Delta#eta^{out};MC/DATA", deta_out_nbins, deta_out_bins);

    delta_eta_outside_gap_ratio->Divide(delta_eta_outside_gap_source, delta_eta_outside_gap_reference, 1., 1., "");
    plot_2histograms(delta_eta_outside_gap_source, "Source", delta_eta_outside_gap_reference, "Reference", output_path_plots, "control_" + prefix_output + "delta_eta_outside_gap", "top_left", true, detail);
    plot_histogram(delta_eta_outside_gap_ratio, output_path_plots, prefix_output + "delta_eta_outside_gap", "ratio","top_right");


//Opening the output root file
    if (detail) { cout<<"Creating " << path_output << "..."<<endl; }
    TFile *output = TFile::Open( path_output.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file..."<<endl; }
    delta_phi_ratio->Write();
    delta_phi_deta1_ratio->Write();
    delta_phi_deta2_ratio->Write();
    delta_phi_deta3_ratio->Write();
    delta_phi_deta4_ratio->Write();
    delta_phi_gap_ratio->Write();
    delta_phi_deta1_gap_ratio->Write();
    delta_phi_deta2_gap_ratio->Write();
    delta_phi_deta3_gap_ratio->Write();
    delta_phi_deta4_gap_ratio->Write();
    delta_phi_nogap_ratio->Write();
    delta_phi_deta1_nogap_ratio->Write();
    delta_phi_deta2_nogap_ratio->Write();
    delta_phi_deta3_nogap_ratio->Write();
    delta_phi_deta4_nogap_ratio->Write();
    leading_pt_inside_gap_ratio->Write();
    leading_eta_star_inside_gap_ratio->Write();
    delta_eta_outside_gap_ratio->Write();
    leading_pt_outside_gap_ratio->Write();
    if (detail) { cout<<"Writing was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    source->Close();
    reference->Close();
    output->Close();

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
