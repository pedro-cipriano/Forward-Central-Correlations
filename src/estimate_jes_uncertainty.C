// Pedro Cipriano, Mar 2013
// DESY, CMS
// Last Update: 22 Mar 2013
//
// estimate_jes_uncertainty()
// calculate the jes uncertainty for pythia 8 - tune 1, pyhtia 6 - tune z2 and herwig 6

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

void calc_jes_unc(TH1D *histogram, TH1D *hist_data, TH1D *hist_data_jes, double *final_jes_unc, int index, bool detail)
{
// calculates the jes uncertainty for a given distribution taking the difference between p6_z2 and h6

//declaring variables
double min = 0.0;
double max = 0.0;
double tot = 0.0;
double ave = 0.0;
double sist = 0.0;
double sist_error = 0.0;
double cont = 0.0;

histogram->Divide(hist_data_jes,hist_data,1.0,1.0,"");

    for(Int_t i=1;i<=histogram->GetNbinsX();i++)
    {
    cont = histogram->GetBinContent(i);
    sist_error = histogram->GetBinError(i);
    sist = 1 - cont;
    if (sist < 0.0) { sist = -sist; }
    histogram->SetBinContent(i,sist);

//control output of the calculation
//    if (detail) { cout << "cont data = " << hist_data->GetBinContent(i) << ", cont data jes = " << hist_data_jes->GetBinContent(i) << ", division = " << cont << ", jes = " <<sist<< " jes error = " << sist_error << endl; }

//checks for minimum, maximum and sums the uncertainties
    if (sist > max) {max = sist;}
    tot = tot + (1-cont);
    if (sist < min || i == 1) { min = sist;}
    }

    histogram->SetEntries(hist_data->GetEntries());

    ave = tot/hist_data->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    final_jes_unc[index*3+0] = ave*100;
    final_jes_unc[index*3+1] = min*100;
    final_jes_unc[index*3+2] = max*100;

}


void show_jes_uncertainties(double *final_jes_unc)
{
//shows the computed jes uncertainties
    cout<<" "<<endl;
    cout<<"jes Uncertainty"<<endl;
    cout<<"Observable                  Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   "<<final_jes_unc[0]<<"  "<<final_jes_unc[1]<<" "<<final_jes_unc[2]<<endl;
    cout<<"Delta phi deta1             "<<final_jes_unc[3]<<"  "<<final_jes_unc[4]<<"  "<<final_jes_unc[5]<<endl;
    cout<<"Delta phi deta2             "<<final_jes_unc[6]<<"   "<<final_jes_unc[7]<<" "<<final_jes_unc[8]<<endl;
    cout<<"Delta phi deta3             "<<final_jes_unc[9]<<"  "<<final_jes_unc[10]<<" "<<final_jes_unc[11]<<endl;
    cout<<"Delta phi deta4             "<<final_jes_unc[12]<<"   "<<final_jes_unc[13]<<" "<<final_jes_unc[14]<<endl;
    cout<<"Delta phi gap               "<<final_jes_unc[15]<<"  "<<final_jes_unc[16]<<" "<<final_jes_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         "<<final_jes_unc[18]<<"  "<<final_jes_unc[19]<<"  "<<final_jes_unc[20]<<endl;
    cout<<"Delta phi deta2 gap         "<<final_jes_unc[21]<<"  "<<final_jes_unc[22]<<" "<<final_jes_unc[23]<<endl;
    cout<<"Delta phi deta3 gap         "<<final_jes_unc[24]<<"  "<<final_jes_unc[25]<<" "<<final_jes_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         "<<final_jes_unc[27]<<"  "<<final_jes_unc[28]<<" "<<final_jes_unc[29]<<endl;
    cout<<"Delta phi nogap             "<<final_jes_unc[30]<<"  "<<final_jes_unc[31]<<" "<<final_jes_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap       "<<final_jes_unc[33]<<"  "<<final_jes_unc[34]<<"  "<<final_jes_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap       "<<final_jes_unc[36]<<"   "<<final_jes_unc[37]<<"  "<<final_jes_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap       "<<final_jes_unc[39]<<"  "<<final_jes_unc[40]<<"  "<<final_jes_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap       "<<final_jes_unc[42]<<"  "<<final_jes_unc[43]<<"  "<<final_jes_unc[44]<<endl;
    cout<<"Leading pT inside gap       "<<final_jes_unc[45]<<"  "<<final_jes_unc[46]<<" "<<final_jes_unc[47]<<endl;
    cout<<"Leading eta star inside gap "<<final_jes_unc[48]<<"  "<<final_jes_unc[49]<<"  "<<final_jes_unc[50]<<endl;
    cout<<"Delta eta outside gap       "<<final_jes_unc[51]<<"  "<<final_jes_unc[52]<<" "<<final_jes_unc[53]<<endl;
    cout<<"Leading pT outside gap      "<<final_jes_unc[54]<<"  "<<final_jes_unc[55]<<" "<<final_jes_unc[56]<<endl;
}


void estimate_jes_uncertainty(string path_data, string path_data_jes, string out_jes_uncertainty, string label, string prefix, string output_path_plots = "../output/jes_uncertainty/", bool detail = false, bool disp_uncertainty = true, bool test = false)
{
//main estimate jes uncertainty routine
//output detail: true-output all the major steps; false-run quietly
//disp_uncertainty: true-display the jes uncartainty estimation values; false-dont show anything

//outputs the configuration
    if (detail) { cout << "Estimate JES Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for Data:          " << path_data << endl; }
    if (detail) { cout << "Input path for Data with JES: " << path_data_jes << endl; }
    if (detail) { cout << "Histogram label:              " << label << endl; }
    if (detail) { cout << "Output path:                  " << out_jes_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:            " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:                 " << detail << endl; }
    if (detail) { cout << "Display Results:              " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:                    " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_file = new TFile( path_data.c_str() );
    TFile *data_jes_file = new TFile( path_data_jes.c_str() );


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

//declaring and initializing the uncertainties array
    double final_jes_unc[19*3];

    for (int i=0; i<= 19*3-1;i++)
    {
    final_jes_unc[i] = 0.0;
    }

//starts the estimation of the jes uncertainty

//estimate the jes uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_data = 0;
    TH1D *delta_phi_data_jes = 0;
    TString delta_phi_name = prefix + "delta_phi";

    data_file->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi",delta_phi_data_jes);
    if (delta_phi_data_jes == 0) { cout << "ak5PF_delta_phi for JES not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D(delta_phi_name,"JES Uncertainty;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi, delta_phi_data, delta_phi_data_jes, final_jes_unc, 0, detail);
    plot_histogram(delta_phi, output_path_plots, "jes_uncertainty_" + label + "_delta_phi", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_data = 0;
    TH1D *delta_phi_deta1_data_jes = 0;
    TString delta_phi_deta1_name = prefix + "delta_phi_deta1";

    data_file->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data_jes);
    if (delta_phi_deta1_data_jes == 0) { cout << "ak5PF_delta_phi_deta1 for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D(delta_phi_deta1_name,"JES Uncertainty;#Delta#phi for 0.4 < #Delta#eta < 2.5;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta1, delta_phi_deta1_data, delta_phi_deta1_data_jes, final_jes_unc, 1, detail);
    plot_histogram(delta_phi_deta1, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta1", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_data = 0;
    TH1D *delta_phi_deta2_data_jes = 0;
    TString delta_phi_deta2_name = prefix + "delta_phi_deta2";

    data_file->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data_jes);
    if (delta_phi_deta2_data_jes == 0) { cout << "ak5PF_delta_phi_deta2 for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D(delta_phi_deta2_name,"JES Uncertainty;#Delta#phi for 2.5 < #Delta#eta < 3.5;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta2, delta_phi_deta2_data, delta_phi_deta2_data_jes, final_jes_unc, 2, detail);
    plot_histogram(delta_phi_deta2, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta2", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_data = 0;
    TH1D *delta_phi_deta3_data_jes = 0;
    TString delta_phi_deta3_name = prefix + "delta_phi_deta3";

    data_file->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data_jes);
    if (delta_phi_deta3_data_jes == 0) { cout << "ak5PF_delta_phi_deta3 for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D(delta_phi_deta3_name,"JES Uncertainty;#Delta#phi for 3.5 < #Delta#eta < 4.5;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta3, delta_phi_deta3_data, delta_phi_deta3_data_jes, final_jes_unc, 3, detail);
    plot_histogram(delta_phi_deta3, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta3", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_data = 0;
    TH1D *delta_phi_deta4_data_jes = 0;
    TString delta_phi_deta4_name = prefix + "delta_phi_deta4";

    data_file->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data_jes);
    if (delta_phi_deta4_data_jes == 0) { cout << "ak5PF_delta_phi_deta4 for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D(delta_phi_deta4_name,"JES Uncertainty;#Delta#phi for 4.5 < #Delta#eta < 7.5;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta4, delta_phi_deta4_data, delta_phi_deta4_data_jes, final_jes_unc, 4, detail);
    plot_histogram(delta_phi_deta4, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta4", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_data = 0;
    TH1D *delta_phi_gap_data_jes = 0;
    TString delta_phi_gap_name = prefix + "delta_phi_gap";

    data_file->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data_jes);
    if (delta_phi_gap_data_jes == 0) { cout << "ak5PF_delta_phi_gap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D(delta_phi_gap_name,"JES Uncertainty;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_gap, delta_phi_gap_data, delta_phi_gap_data_jes, final_jes_unc, 5, detail);
    plot_histogram(delta_phi_gap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_data = 0;
    TH1D *delta_phi_deta1_gap_data_jes = 0;
    TString delta_phi_deta1_gap_name = prefix + "delta_phi_deta1_gap";

    data_file->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data_jes);
    if (delta_phi_deta1_gap_data_jes == 0) { cout << "ak5PF_delta_phi_deta1_gap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D(delta_phi_deta1_gap_name,"JES Uncertainty for 0.4 < #Delta#eta < 2.5;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta1_gap, delta_phi_deta1_gap_data, delta_phi_deta1_gap_data_jes, final_jes_unc, 6, detail);
    plot_histogram(delta_phi_deta1_gap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta1_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_data = 0;
    TH1D *delta_phi_deta2_gap_data_jes = 0;
    TString delta_phi_deta2_gap_name = prefix + "delta_phi_deta2_gap";

    data_file->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data_jes);
    if (delta_phi_deta2_gap_data_jes == 0) { cout << "ak5PF_delta_phi_deta2_gap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D(delta_phi_deta2_gap_name,"JES Uncertainty for 2.5 < #Delta#eta < 3.5;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta2_gap, delta_phi_deta2_gap_data, delta_phi_deta2_gap_data_jes, final_jes_unc, 7, detail);
    plot_histogram(delta_phi_deta2_gap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta2_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_data = 0;
    TH1D *delta_phi_deta3_gap_data_jes = 0;
    TString delta_phi_deta3_gap_name = prefix + "delta_phi_deta3_gap";

    data_file->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data_jes);
    if (delta_phi_deta3_gap_data_jes == 0) { cout << "ak5PF_delta_phi_deta3_gap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D(delta_phi_deta3_gap_name,"JES Uncertainty for 3.5 < #Delta#eta < 4.5;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta3_gap, delta_phi_deta3_gap_data, delta_phi_deta3_gap_data_jes, final_jes_unc, 8, detail);
    plot_histogram(delta_phi_deta3_gap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta3_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_data = 0;
    TH1D *delta_phi_deta4_gap_data_jes = 0;
    TString delta_phi_deta4_gap_name = prefix + "delta_phi_deta4_gap";

    data_file->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data_jes);
    if (delta_phi_deta4_gap_data_jes == 0) { cout << "ak5PF_delta_phi_deta4_gap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D(delta_phi_deta4_gap_name,"JES Uncertainty for 4.5 < #Delta#eta < 7.5;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta4_gap, delta_phi_deta4_gap_data, delta_phi_deta4_gap_data_jes, final_jes_unc, 9, detail);
    plot_histogram(delta_phi_deta4_gap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta4_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_data = 0;
    TH1D *delta_phi_nogap_data_jes = 0;
    TString delta_phi_nogap_name = prefix + "delta_phi_nogap";

    data_file->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data_jes);
    if (delta_phi_nogap_data_jes == 0) { cout << "ak5PF_delta_phi_nogap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D(delta_phi_nogap_name,"JES Uncertainty nogap;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_nogap, delta_phi_nogap_data, delta_phi_nogap_data_jes, final_jes_unc, 10, detail);
    plot_histogram(delta_phi_nogap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_nogap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_data = 0;
    TH1D *delta_phi_deta1_nogap_data_jes = 0;
    TString delta_phi_deta1_nogap_name = prefix + "delta_phi_deta1_nogap";

    data_file->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data_jes);
    if (delta_phi_deta1_nogap_data_jes == 0) { cout << "ak5PF_delta_phi_deta1_nogap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D(delta_phi_deta1_nogap_name,"JES Uncertainty for 0.4 < #Delta#eta < 2.5 nogap;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta1_nogap, delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_data_jes, final_jes_unc, 11, detail);
    plot_histogram(delta_phi_deta1_nogap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta1_nogap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_data = 0;
    TH1D *delta_phi_deta2_nogap_data_jes = 0;
    TString delta_phi_deta2_nogap_name = prefix + "delta_phi_deta2_nogap";

    data_file->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data_jes);
    if (delta_phi_deta2_nogap_data_jes == 0) { cout << "ak5PF_delta_phi_deta2_nogap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D(delta_phi_deta2_nogap_name,"JES Uncertainty for 2.5 < #Delta#eta < 3.5 nogap;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta2_nogap, delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_data_jes, final_jes_unc, 12, detail);
    plot_histogram(delta_phi_deta2_nogap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta2_nogap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_data = 0;
    TH1D *delta_phi_deta3_nogap_data_jes = 0;
    TString delta_phi_deta3_nogap_name = prefix + "delta_phi_deta3_nogap";

    data_file->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data_jes);
    if (delta_phi_deta3_nogap_data_jes == 0) { cout << "ak5PF_delta_phi_deta3_nogap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D(delta_phi_deta3_nogap_name,"JES Uncertainty for 3.5 < #Delta#eta < 4.5 nogap;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta3_nogap, delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_data_jes, final_jes_unc, 13, detail);
    plot_histogram(delta_phi_deta3_nogap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta3_nogap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_data = 0;
    TH1D *delta_phi_deta4_nogap_data_jes = 0;
    TString delta_phi_deta4_nogap_name = prefix + "delta_phi_deta4_nogap";

    data_file->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data_jes);
    if (delta_phi_deta4_nogap_data_jes == 0) { cout << "ak5PF_delta_phi_deta4_nogap for JES not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D(delta_phi_deta4_nogap_name,"JES Uncertainty for 4.5 < #Delta#eta < 7.5 nogap;#Delta#phi;JES Uncertainty", dphi_nbins, dphi_bins);

    calc_jes_unc(delta_phi_deta4_nogap, delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_data_jes, final_jes_unc, 14, detail);
    plot_histogram(delta_phi_deta4_nogap, output_path_plots, "jes_uncertainty_" + label + "_delta_phi_deta4_nogap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pt inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_data = 0;
    TH1D *leading_pt_inside_gap_data_jes = 0;
    TString leading_pt_inside_gap_name = prefix + "leading_pt_inside_gap";

    data_file->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data_jes);
    if (leading_pt_inside_gap_data_jes == 0) { cout << "ak5PF_leading_pt_inside_gap for JES not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D(leading_pt_inside_gap_name,"JES Uncertainty for leading pt inside gap;p_{T} inside;JES Uncertainty", in_nbins, in_bins);

    calc_jes_unc(leading_pt_inside_gap, leading_pt_inside_gap_data, leading_pt_inside_gap_data_jes, final_jes_unc, 15, detail);
    plot_histogram(leading_pt_inside_gap, output_path_plots, "jes_uncertainty_" + label + "_leading_pt_inside_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_data = 0;
    TH1D *leading_eta_star_inside_gap_data_jes = 0;
    TString leading_eta_star_inside_gap_name = prefix + "leading_eta_star_inside_gap";

    data_file->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data_jes);
    if (leading_eta_star_inside_gap_data_jes == 0) { cout << "ak5PF_leading_eta_star_inside_gap for JES not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D(leading_eta_star_inside_gap_name,"JES Uncertainty for leading eta* inside gap;#eta* inside;JES Uncertainty", etastar_nbins, etastar_bins);

    calc_jes_unc(leading_eta_star_inside_gap, leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_data_jes, final_jes_unc, 16, detail);
    plot_histogram(leading_eta_star_inside_gap, output_path_plots, "jes_uncertainty_" + label + "_leading_eta_star_inside_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta eta outside gap"<<endl; }

    TH1D *delta_eta_outside_gap_data = 0;
    TH1D *delta_eta_outside_gap_data_jes = 0;
    TString delta_eta_outside_gap_name = prefix + "delta_eta_outside_gap";

    data_file->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data_jes);
    if (delta_eta_outside_gap_data_jes == 0) { cout << "ak5PF_delta_eta_outside_gap for JES not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D(delta_eta_outside_gap_name,"JES Uncertainty for leading pt inside gap;p_{T} inside;JES Uncertainty", deta_out_nbins, deta_out_bins);

    calc_jes_unc(delta_eta_outside_gap, delta_eta_outside_gap_data, delta_eta_outside_gap_data_jes, final_jes_unc, 17, detail);
    plot_histogram(delta_eta_outside_gap, output_path_plots, "jes_uncertainty_" + label + "_delta_eta_outside_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//estimate the jes uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pt outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_data = 0;
    TH1D *leading_pt_outside_gap_data_jes = 0;
    TString leading_pt_outside_gap_name = prefix + "leading_pt_outside_gap";

    data_file->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    data_jes_file->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data_jes);
    if (leading_pt_outside_gap_data_jes == 0) { cout << "ak5PF_leading_pt_outside_gap for JES not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D(leading_pt_outside_gap_name,"JES Uncertainty for leading pt outside gap;p_{T} outside;JES Uncertainty", out_nbins, out_bins);

    calc_jes_unc(leading_pt_outside_gap, leading_pt_outside_gap_data, leading_pt_outside_gap_data_jes, final_jes_unc, 18, detail);
    plot_histogram(leading_pt_outside_gap, output_path_plots, "jes_uncertainty_" + label + "_leading_pt_outside_gap", "Jet Energy Scale Uncertainty " + label, "top_left", true);


//output the error variation
    if (detail) { cout<<"Display the final jes uncertainties..."<<endl; }
    if (disp_uncertainty) { show_jes_uncertainties(final_jes_unc); }

//Opening the output root file
    if (detail) { cout<<"Creating " << out_jes_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( out_jes_uncertainty.c_str() , "RECREATE");

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
    data_file->Close();
    data_jes_file->Close();
    data_output->Close();

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
