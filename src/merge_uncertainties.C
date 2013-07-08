// Pedro Cipriano, Mar 2011
// DESY, CMS
// Last Update: 22 Mar 2013
//
// merge_uncertaintes()

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


void merge_uncertainty(TH1D *hist, TH1D *hist_jetmettau, TH1D *hist_jetmet, TH1D *hist_jet, double *merged_unc, int index, bool detail)
{

double min = 0.0;
double max = 0.0;
double tot = 0.0;
double ave = 0.0;
double syst = 0.0;

int entries_jetmettau = hist_jetmettau->GetEntries();
int entries_jetmet = hist_jetmet->GetEntries();
//int entries_jet = hist_jet->GetEntries();
int entries_total = entries_jetmettau + entries_jetmet; // + entries_jet;

double frac_jetmettau = (double) entries_jetmettau / (double) entries_total;
double frac_jetmet = (double) entries_jetmet / (double) entries_total;
//double frac_jet = (double) entries_jet / (double) entries_total;

//if (detail) { cout << "Fractions : " << frac_jetmettau << " " << frac_jetmet << " " << frac_jet << endl; }

hist->Add(hist_jetmettau,frac_jetmettau);
hist->Add(hist_jetmet,frac_jetmet);
//hist->Add(hist_jet,frac_jet);

//loops over the bins, gets and saves the deviations
    for(Int_t i=1;i<=hist->GetNbinsX();i++)
    {
    syst = hist->GetBinContent(i);
    if (syst > max) {max = syst;}
    tot = tot + syst;
    if (syst < min || i == 1) { min = syst;}
    }

//calculates the average model uncertainty
    ave = tot/hist->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    merged_unc[index*3+0] = ave*100;
    merged_unc[index*3+1] = min*100;
    merged_unc[index*3+2] = max*100;

}


void show_merged_uncertainties(double *merged_unc)
{
//shows the computed jes uncertainties
    cout<<" "<<endl;
    cout<<"Merged Uncertainty"<<endl;
    cout<<"Observable                  Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   "<<merged_unc[0]<<"  "<<merged_unc[1]<<" "<<merged_unc[2]<<endl;
    cout<<"Delta phi deta1             "<<merged_unc[3]<<"  "<<merged_unc[4]<<"  "<<merged_unc[5]<<endl;
    cout<<"Delta phi deta2             "<<merged_unc[6]<<"   "<<merged_unc[7]<<" "<<merged_unc[8]<<endl;
    cout<<"Delta phi deta3             "<<merged_unc[9]<<"  "<<merged_unc[10]<<" "<<merged_unc[11]<<endl;
    cout<<"Delta phi deta4             "<<merged_unc[12]<<"   "<<merged_unc[13]<<" "<<merged_unc[14]<<endl;
    cout<<"Delta phi gap               "<<merged_unc[15]<<"  "<<merged_unc[16]<<" "<<merged_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         "<<merged_unc[18]<<"  "<<merged_unc[19]<<"  "<<merged_unc[20]<<endl;
    cout<<"Delta phi deta2 gap         "<<merged_unc[21]<<"  "<<merged_unc[22]<<" "<<merged_unc[23]<<endl;
    cout<<"Delta phi deta3 gap         "<<merged_unc[24]<<"  "<<merged_unc[25]<<" "<<merged_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         "<<merged_unc[27]<<"  "<<merged_unc[28]<<" "<<merged_unc[29]<<endl;
    cout<<"Delta phi nogap             "<<merged_unc[30]<<"  "<<merged_unc[31]<<" "<<merged_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap       "<<merged_unc[33]<<"  "<<merged_unc[34]<<"  "<<merged_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap       "<<merged_unc[36]<<"   "<<merged_unc[37]<<"  "<<merged_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap       "<<merged_unc[39]<<"  "<<merged_unc[40]<<"  "<<merged_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap       "<<merged_unc[42]<<"  "<<merged_unc[43]<<"  "<<merged_unc[44]<<endl;
    cout<<"Leading pT inside gap       "<<merged_unc[45]<<"  "<<merged_unc[46]<<" "<<merged_unc[47]<<endl;
    cout<<"Leading eta star inside gap "<<merged_unc[48]<<"  "<<merged_unc[49]<<"  "<<merged_unc[50]<<endl;
    cout<<"Delta eta outside gap       "<<merged_unc[51]<<"  "<<merged_unc[52]<<" "<<merged_unc[53]<<endl;
    cout<<"Leading pT outside gap      "<<merged_unc[54]<<"  "<<merged_unc[55]<<" "<<merged_unc[56]<<endl;
}


void merge_uncertainties(string path_data_jetmettau, string path_data_jetmet, string path_data_jet, string out_uncertainty, string label_in, string label_out, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Merge Uncertainties Configuration"<<endl; }
    if (detail) { cout << "Input path for JetMETTau file: " << path_data_jetmettau << endl; }
    if (detail) { cout << "Input path for JetMET file:    " << path_data_jetmet << endl; }
    if (detail) { cout << "Input path for Jet file:       " << path_data_jet << endl; }
    if (detail) { cout << "Label In:                      " << label_in << endl; }
    if (detail) { cout << "Label Out:                     " << label_out << endl; }
    if (detail) { cout << "Output path:                   " << out_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:             " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:                  " << detail << endl; }
    if (detail) { cout << "Display Results:               " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:                     " << test << endl; }


//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_jetmettau = new TFile( path_data_jetmettau.c_str() );
    TFile *data_jetmet = new TFile( path_data_jetmet.c_str() );
    TFile *data_jet = new TFile( path_data_jet.c_str() );


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
    double merged_unc[19*3];

    for (int i=0; i<= 19*3-1;i++)
    {
    merged_unc[i] = 0.0;
    }

//starts merging uncertainties

//merge the uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_jetmettau = 0;
    TH1D *delta_phi_jetmet = 0;
    TH1D *delta_phi_jet = 0;
    TString delta_phi_name_in = label_in + "delta_phi";
    TString delta_phi_name_out = label_out + "delta_phi";

    data_jetmettau->GetObject(delta_phi_name_in,delta_phi_jetmettau);
    if (delta_phi_jetmettau == 0) { cout << delta_phi_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_name_in,delta_phi_jetmet);
    if (delta_phi_jetmet == 0) { cout << delta_phi_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_name_in,delta_phi_jet);
    if (delta_phi_jet == 0) { cout << delta_phi_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D(delta_phi_name_out,"Uncertainty;#Delta#phi;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi, delta_phi_jetmettau, delta_phi_jetmet, delta_phi_jet, merged_unc, 0, detail);
    plot_histogram(delta_phi, output_path_plots, label_out + "delta_phi",  label_out + "Uncertainty ", "top_left", true);

//merge the uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_jetmettau = 0;
    TH1D *delta_phi_deta1_jetmet = 0;
    TH1D *delta_phi_deta1_jet = 0;
    TString delta_phi_deta1_name_in = label_in + "delta_phi_deta1";
    TString delta_phi_deta1_name_out = label_out + "delta_phi_deta1";

    data_jetmettau->GetObject(delta_phi_deta1_name_in,delta_phi_deta1_jetmettau);
    if (delta_phi_deta1_jetmettau == 0) { cout << delta_phi_deta1_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta1_name_in,delta_phi_deta1_jetmet);
    if (delta_phi_deta1_jetmet == 0) { cout << delta_phi_deta1_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta1_name_in,delta_phi_deta1_jet);
    if (delta_phi_deta1_jet == 0) { cout << delta_phi_deta1_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D(delta_phi_deta1_name_out,"Uncertainty;#Delta#phi deta1;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta1, delta_phi_deta1_jetmettau, delta_phi_deta1_jetmet, delta_phi_deta1_jet, merged_unc, 1, detail);
    plot_histogram(delta_phi_deta1, output_path_plots, label_out + "delta_phi_deta1",  label_out + "Uncertainty", "top_left", true);

//merge the uncertainty for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_jetmettau = 0;
    TH1D *delta_phi_deta2_jetmet = 0;
    TH1D *delta_phi_deta2_jet = 0;
    TString delta_phi_deta2_name_in = label_in + "delta_phi_deta2";
    TString delta_phi_deta2_name_out = label_out + "delta_phi_deta2";

    data_jetmettau->GetObject(delta_phi_deta2_name_in,delta_phi_deta2_jetmettau);
    if (delta_phi_deta2_jetmettau == 0) { cout << delta_phi_deta2_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta2_name_in,delta_phi_deta2_jetmet);
    if (delta_phi_deta2_jetmet == 0) { cout << delta_phi_deta2_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta2_name_in,delta_phi_deta2_jet);
    if (delta_phi_deta2_jet == 0) { cout << delta_phi_deta2_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D(delta_phi_deta2_name_out,"Uncertainty;#Delta#phi deta2;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta2, delta_phi_deta2_jetmettau, delta_phi_deta2_jetmet, delta_phi_deta2_jet, merged_unc, 2, detail);
    plot_histogram(delta_phi_deta2, output_path_plots, label_out + "delta_phi_deta2",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_jetmettau = 0;
    TH1D *delta_phi_deta3_jetmet = 0;
    TH1D *delta_phi_deta3_jet = 0;
    TString delta_phi_deta3_name_in = label_in + "delta_phi_deta3";
    TString delta_phi_deta3_name_out = label_out + "delta_phi_deta3";

    data_jetmettau->GetObject(delta_phi_deta3_name_in,delta_phi_deta3_jetmettau);
    if (delta_phi_deta3_jetmettau == 0) { cout << delta_phi_deta3_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta3_name_in,delta_phi_deta3_jetmet);
    if (delta_phi_deta3_jetmet == 0) { cout << delta_phi_deta3_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta3_name_in,delta_phi_deta3_jet);
    if (delta_phi_deta3_jet == 0) { cout << delta_phi_deta3_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D(delta_phi_deta3_name_out,"Uncertainty;#Delta#phi deta3;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta3, delta_phi_deta3_jetmettau, delta_phi_deta3_jetmet, delta_phi_deta3_jet, merged_unc, 3, detail);
    plot_histogram(delta_phi_deta3, output_path_plots, label_out + "delta_phi_deta3",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_jetmettau = 0;
    TH1D *delta_phi_deta4_jetmet = 0;
    TH1D *delta_phi_deta4_jet = 0;
    TString delta_phi_deta4_name_in = label_in + "delta_phi_deta4";
    TString delta_phi_deta4_name_out = label_out + "delta_phi_deta4";

    data_jetmettau->GetObject(delta_phi_deta4_name_in,delta_phi_deta4_jetmettau);
    if (delta_phi_deta4_jetmettau == 0) { cout << delta_phi_deta4_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta4_name_in,delta_phi_deta4_jetmet);
    if (delta_phi_deta4_jetmet == 0) { cout << delta_phi_deta4_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta4_name_in,delta_phi_deta4_jet);
    if (delta_phi_deta4_jet == 0) { cout << delta_phi_deta4_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D(delta_phi_deta4_name_out,"Uncertainty;#Delta#phi deta4;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta4, delta_phi_deta4_jetmettau, delta_phi_deta4_jetmet, delta_phi_deta4_jet, merged_unc, 4, detail);
    plot_histogram(delta_phi_deta4, output_path_plots, label_out + "delta_phi_deta4",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_jetmettau = 0;
    TH1D *delta_phi_gap_jetmet = 0;
    TH1D *delta_phi_gap_jet = 0;
    TString delta_phi_gap_name_in = label_in + "delta_phi_gap";
    TString delta_phi_gap_name_out = label_out + "delta_phi_gap";

    data_jetmettau->GetObject(delta_phi_gap_name_in,delta_phi_gap_jetmettau);
    if (delta_phi_gap_jetmettau == 0) { cout << delta_phi_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_gap_name_in,delta_phi_gap_jetmet);
    if (delta_phi_gap_jetmet == 0) { cout << delta_phi_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_gap_name_in,delta_phi_gap_jet);
    if (delta_phi_gap_jet == 0) { cout << delta_phi_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D(delta_phi_gap_name_out,"Uncertainty;#Delta#phi gap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_gap, delta_phi_gap_jetmettau, delta_phi_gap_jetmet, delta_phi_gap_jet, merged_unc, 5, detail);
    plot_histogram(delta_phi_gap, output_path_plots, label_out + "delta_phi_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_jetmettau = 0;
    TH1D *delta_phi_deta1_gap_jetmet = 0;
    TH1D *delta_phi_deta1_gap_jet = 0;
    TString delta_phi_deta1_gap_name_in = label_in + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_name_out = label_out + "delta_phi_deta1_gap";

    data_jetmettau->GetObject(delta_phi_deta1_gap_name_in,delta_phi_deta1_gap_jetmettau);
    if (delta_phi_deta1_gap_jetmettau == 0) { cout << delta_phi_deta1_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta1_gap_name_in,delta_phi_deta1_gap_jetmet);
    if (delta_phi_deta1_gap_jetmet == 0) { cout << delta_phi_deta1_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta1_gap_name_in,delta_phi_deta1_gap_jet);
    if (delta_phi_deta1_gap_jet == 0) { cout << delta_phi_deta1_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D(delta_phi_deta1_gap_name_out,"Uncertainty;#Delta#phi deta1 gap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta1_gap, delta_phi_deta1_gap_jetmettau, delta_phi_deta1_gap_jetmet, delta_phi_deta1_gap_jet, merged_unc, 6, detail);
    plot_histogram(delta_phi_deta1_gap, output_path_plots, label_out + "delta_phi_deta1_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_jetmettau = 0;
    TH1D *delta_phi_deta2_gap_jetmet = 0;
    TH1D *delta_phi_deta2_gap_jet = 0;
    TString delta_phi_deta2_gap_name_in = label_in + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_name_out = label_out + "delta_phi_deta2_gap";

    data_jetmettau->GetObject(delta_phi_deta2_gap_name_in,delta_phi_deta2_gap_jetmettau);
    if (delta_phi_deta2_gap_jetmettau == 0) { cout << delta_phi_deta2_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta2_gap_name_in,delta_phi_deta2_gap_jetmet);
    if (delta_phi_deta2_gap_jetmet == 0) { cout << delta_phi_deta2_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta2_gap_name_in,delta_phi_deta2_gap_jet);
    if (delta_phi_deta2_gap_jet == 0) { cout << delta_phi_deta2_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D(delta_phi_deta2_gap_name_out,"Uncertainty;#Delta#phi deta2 gap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta2_gap, delta_phi_deta2_gap_jetmettau, delta_phi_deta2_gap_jetmet, delta_phi_deta2_gap_jet, merged_unc, 7, detail);
    plot_histogram(delta_phi_deta2_gap, output_path_plots, label_out + "delta_phi_deta2_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_jetmettau = 0;
    TH1D *delta_phi_deta3_gap_jetmet = 0;
    TH1D *delta_phi_deta3_gap_jet = 0;
    TString delta_phi_deta3_gap_name_in = label_in + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_name_out = label_out + "delta_phi_deta3_gap";

    data_jetmettau->GetObject(delta_phi_deta3_gap_name_in,delta_phi_deta3_gap_jetmettau);
    if (delta_phi_deta3_gap_jetmettau == 0) { cout << delta_phi_deta3_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta3_gap_name_in,delta_phi_deta3_gap_jetmet);
    if (delta_phi_deta3_gap_jetmet == 0) { cout << delta_phi_deta3_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta3_gap_name_in,delta_phi_deta3_gap_jet);
    if (delta_phi_deta3_gap_jet == 0) { cout << delta_phi_deta3_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D(delta_phi_deta3_gap_name_out,"Uncertainty;#Delta#phi deta3 gap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta3_gap, delta_phi_deta3_gap_jetmettau, delta_phi_deta3_gap_jetmet, delta_phi_deta3_gap_jet, merged_unc, 8, detail);
    plot_histogram(delta_phi_deta3_gap, output_path_plots, label_out + "delta_phi_deta3_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_jetmettau = 0;
    TH1D *delta_phi_deta4_gap_jetmet = 0;
    TH1D *delta_phi_deta4_gap_jet = 0;
    TString delta_phi_deta4_gap_name_in = label_in + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_name_out = label_out + "delta_phi_deta4_gap";

    data_jetmettau->GetObject(delta_phi_deta4_gap_name_in,delta_phi_deta4_gap_jetmettau);
    if (delta_phi_deta4_gap_jetmettau == 0) { cout << delta_phi_deta4_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta4_gap_name_in,delta_phi_deta4_gap_jetmet);
    if (delta_phi_deta4_gap_jetmet == 0) { cout << delta_phi_deta4_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta4_gap_name_in,delta_phi_deta4_gap_jet);
    if (delta_phi_deta4_gap_jet == 0) { cout << delta_phi_deta4_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D(delta_phi_deta4_gap_name_out,"Uncertainty;#Delta#phi deta4 gap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta4_gap, delta_phi_deta4_gap_jetmettau, delta_phi_deta4_gap_jetmet, delta_phi_deta4_gap_jet, merged_unc, 9, detail);
    plot_histogram(delta_phi_deta4_gap, output_path_plots, label_out + "delta_phi_deta4_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_jetmettau = 0;
    TH1D *delta_phi_nogap_jetmet = 0;
    TH1D *delta_phi_nogap_jet = 0;
    TString delta_phi_nogap_name_in = label_in + "delta_phi_nogap";
    TString delta_phi_nogap_name_out = label_out + "delta_phi_nogap";

    data_jetmettau->GetObject(delta_phi_nogap_name_in,delta_phi_nogap_jetmettau);
    if (delta_phi_nogap_jetmettau == 0) { cout << delta_phi_nogap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_nogap_name_in,delta_phi_nogap_jetmet);
    if (delta_phi_nogap_jetmet == 0) { cout << delta_phi_nogap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_nogap_name_in,delta_phi_nogap_jet);
    if (delta_phi_nogap_jet == 0) { cout << delta_phi_nogap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D(delta_phi_nogap_name_out,"Uncertainty;#Delta#phi nogap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_nogap, delta_phi_nogap_jetmettau, delta_phi_nogap_jetmet, delta_phi_nogap_jet, merged_unc, 10, detail);
    plot_histogram(delta_phi_nogap, output_path_plots, label_out + "delta_phi_nogap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_jetmettau = 0;
    TH1D *delta_phi_deta1_nogap_jetmet = 0;
    TH1D *delta_phi_deta1_nogap_jet = 0;
    TString delta_phi_deta1_nogap_name_in = label_in + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_name_out = label_out + "delta_phi_deta1_nogap";

    data_jetmettau->GetObject(delta_phi_deta1_nogap_name_in,delta_phi_deta1_nogap_jetmettau);
    if (delta_phi_deta1_nogap_jetmettau == 0) { cout << delta_phi_deta1_nogap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta1_nogap_name_in,delta_phi_deta1_nogap_jetmet);
    if (delta_phi_deta1_nogap_jetmet == 0) { cout << delta_phi_deta1_nogap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta1_nogap_name_in,delta_phi_deta1_nogap_jet);
    if (delta_phi_deta1_nogap_jet == 0) { cout << delta_phi_deta1_nogap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D(delta_phi_deta1_nogap_name_out,"Uncertainty;#Delta#phi deta1 nogap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jetmettau, delta_phi_deta1_nogap_jetmet, delta_phi_deta1_nogap_jet, merged_unc, 11, detail);
    plot_histogram(delta_phi_deta1_nogap, output_path_plots, label_out + "delta_phi_deta1_nogap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_jetmettau = 0;
    TH1D *delta_phi_deta2_nogap_jetmet = 0;
    TH1D *delta_phi_deta2_nogap_jet = 0;
    TString delta_phi_deta2_nogap_name_in = label_in + "delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_name_out = label_out + "delta_phi_deta2_nogap";

    data_jetmettau->GetObject(delta_phi_deta2_nogap_name_in,delta_phi_deta2_nogap_jetmettau);
    if (delta_phi_deta2_nogap_jetmettau == 0) { cout << delta_phi_deta2_nogap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta2_nogap_name_in,delta_phi_deta2_nogap_jetmet);
    if (delta_phi_deta2_nogap_jetmet == 0) { cout << delta_phi_deta2_nogap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta2_nogap_name_in,delta_phi_deta2_nogap_jet);
    if (delta_phi_deta2_nogap_jet == 0) { cout << delta_phi_deta2_nogap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D(delta_phi_deta2_nogap_name_out,"Uncertainty;#Delta#phi deta2 nogap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jetmettau, delta_phi_deta2_nogap_jetmet, delta_phi_deta2_nogap_jet, merged_unc, 12, detail);
    plot_histogram(delta_phi_deta2_nogap, output_path_plots, label_out + "delta_phi_deta2_nogap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_jetmettau = 0;
    TH1D *delta_phi_deta3_nogap_jetmet = 0;
    TH1D *delta_phi_deta3_nogap_jet = 0;
    TString delta_phi_deta3_nogap_name_in = label_in + "delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_name_out = label_out + "delta_phi_deta3_nogap";

    data_jetmettau->GetObject(delta_phi_deta3_nogap_name_in,delta_phi_deta3_nogap_jetmettau);
    if (delta_phi_deta3_nogap_jetmettau == 0) { cout << delta_phi_deta3_nogap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta3_nogap_name_in,delta_phi_deta3_nogap_jetmet);
    if (delta_phi_deta3_nogap_jetmet == 0) { cout << delta_phi_deta3_nogap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta3_nogap_name_in,delta_phi_deta3_nogap_jet);
    if (delta_phi_deta3_nogap_jet == 0) { cout << delta_phi_deta3_nogap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D(delta_phi_deta3_nogap_name_out,"Uncertainty;#Delta#phi deta3 nogap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jetmettau, delta_phi_deta3_nogap_jetmet, delta_phi_deta3_nogap_jet, merged_unc, 13, detail);
    plot_histogram(delta_phi_deta3_nogap, output_path_plots, label_out + "delta_phi_deta3_nogap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_jetmettau = 0;
    TH1D *delta_phi_deta4_nogap_jetmet = 0;
    TH1D *delta_phi_deta4_nogap_jet = 0;
    TString delta_phi_deta4_nogap_name_in = label_in + "delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_name_out = label_out + "delta_phi_deta4_nogap";

    data_jetmettau->GetObject(delta_phi_deta4_nogap_name_in,delta_phi_deta4_nogap_jetmettau);
    if (delta_phi_deta4_nogap_jetmettau == 0) { cout << delta_phi_deta4_nogap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_phi_deta4_nogap_name_in,delta_phi_deta4_nogap_jetmet);
    if (delta_phi_deta4_nogap_jetmet == 0) { cout << delta_phi_deta4_nogap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_phi_deta4_nogap_name_in,delta_phi_deta4_nogap_jet);
    if (delta_phi_deta4_nogap_jet == 0) { cout << delta_phi_deta4_nogap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D(delta_phi_deta4_nogap_name_out,"Uncertainty;#Delta#phi deta4 nogap;Uncertainty", dphi_nbins, dphi_bins);

    merge_uncertainty(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jetmettau, delta_phi_deta4_nogap_jetmet, delta_phi_deta4_nogap_jet, merged_unc, 14, detail);
    plot_histogram(delta_phi_deta4_nogap, output_path_plots, label_out + "delta_phi_deta4_nogap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pt inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_jetmettau = 0;
    TH1D *leading_pt_inside_gap_jetmet = 0;
    TH1D *leading_pt_inside_gap_jet = 0;
    TString leading_pt_inside_gap_name_in = label_in + "leading_pt_inside_gap";
    TString leading_pt_inside_gap_name_out = label_out + "leading_pt_inside_gap";

    data_jetmettau->GetObject(leading_pt_inside_gap_name_in,leading_pt_inside_gap_jetmettau);
    if (leading_pt_inside_gap_jetmettau == 0) { cout << leading_pt_inside_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(leading_pt_inside_gap_name_in,leading_pt_inside_gap_jetmet);
    if (leading_pt_inside_gap_jetmet == 0) { cout << leading_pt_inside_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(leading_pt_inside_gap_name_in,leading_pt_inside_gap_jet);
    if (leading_pt_inside_gap_jet == 0) { cout << leading_pt_inside_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D(leading_pt_inside_gap_name_out,"Uncertainty;Leading pt inside gap;Uncertainty", in_nbins, in_bins);

    merge_uncertainty(leading_pt_inside_gap, leading_pt_inside_gap_jetmettau, leading_pt_inside_gap_jetmet, leading_pt_inside_gap_jet, merged_unc, 15, detail);
    plot_histogram(leading_pt_inside_gap, output_path_plots, label_out + "leading_pt_inside_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_jetmettau = 0;
    TH1D *leading_eta_star_inside_gap_jetmet = 0;
    TH1D *leading_eta_star_inside_gap_jet = 0;
    TString leading_eta_star_inside_gap_name_in = label_in + "leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_name_out = label_out + "leading_eta_star_inside_gap";

    data_jetmettau->GetObject(leading_eta_star_inside_gap_name_in,leading_eta_star_inside_gap_jetmettau);
    if (leading_eta_star_inside_gap_jetmettau == 0) { cout << leading_eta_star_inside_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(leading_eta_star_inside_gap_name_in,leading_eta_star_inside_gap_jetmet);
    if (leading_eta_star_inside_gap_jetmet == 0) { cout << leading_eta_star_inside_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(leading_eta_star_inside_gap_name_in,leading_eta_star_inside_gap_jet);
    if (leading_eta_star_inside_gap_jet == 0) { cout << leading_eta_star_inside_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D(leading_eta_star_inside_gap_name_out,"Uncertainty;Leading eta* inside gap;Uncertainty", etastar_nbins, etastar_bins);

    merge_uncertainty(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jetmettau, leading_eta_star_inside_gap_jetmet, leading_eta_star_inside_gap_jet, merged_unc, 16, detail);
    plot_histogram(leading_eta_star_inside_gap, output_path_plots, label_out + "leading_eta_star_inside_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta eta outside gap"<<endl; }

    TH1D *delta_eta_outside_gap_jetmettau = 0;
    TH1D *delta_eta_outside_gap_jetmet = 0;
    TH1D *delta_eta_outside_gap_jet = 0;
    TString delta_eta_outside_gap_name_in = label_in + "delta_eta_outside_gap";
    TString delta_eta_outside_gap_name_out = label_out + "delta_eta_outside_gap";

    data_jetmettau->GetObject(delta_eta_outside_gap_name_in,delta_eta_outside_gap_jetmettau);
    if (delta_eta_outside_gap_jetmettau == 0) { cout << delta_eta_outside_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(delta_eta_outside_gap_name_in,delta_eta_outside_gap_jetmet);
    if (delta_eta_outside_gap_jetmet == 0) { cout << delta_eta_outside_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(delta_eta_outside_gap_name_in,delta_eta_outside_gap_jet);
    if (delta_eta_outside_gap_jet == 0) { cout << delta_eta_outside_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D(delta_eta_outside_gap_name_out,"Uncertainty;#Delta#eta outside gap;Uncertainty", deta_out_nbins, deta_out_bins);

    merge_uncertainty(delta_eta_outside_gap, delta_eta_outside_gap_jetmettau, delta_eta_outside_gap_jetmet, delta_eta_outside_gap_jet, merged_unc, 17, detail);
    plot_histogram(delta_eta_outside_gap, output_path_plots, label_out + "delta_eta_outside_gap",  label_out + "Uncertainty", "top_left", true);


//merge the uncertainty for leading pt outside gap distribution
    if (detail) { cout<<"Leading pt outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_jetmettau = 0;
    TH1D *leading_pt_outside_gap_jetmet = 0;
    TH1D *leading_pt_outside_gap_jet = 0;
    TString leading_pt_outside_gap_name_in = label_in + "leading_pt_outside_gap";
    TString leading_pt_outside_gap_name_out = label_out + "leading_pt_outside_gap";

    data_jetmettau->GetObject(leading_pt_outside_gap_name_in,leading_pt_outside_gap_jetmettau);
    if (leading_pt_outside_gap_jetmettau == 0) { cout << leading_pt_outside_gap_name_in << " on jetmettau not found!" << endl; return; }
    data_jetmet->GetObject(leading_pt_outside_gap_name_in,leading_pt_outside_gap_jetmet);
    if (leading_pt_outside_gap_jetmet == 0) { cout << leading_pt_outside_gap_name_in << " on jetmet not found!" << endl; return; }
    data_jet->GetObject(leading_pt_outside_gap_name_in,leading_pt_outside_gap_jet);
    if (leading_pt_outside_gap_jet == 0) { cout << leading_pt_outside_gap_name_in << " on jet not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D(leading_pt_outside_gap_name_out,"Uncertainty;Leading pt outside gap;Uncertainty", out_nbins, out_bins);

    merge_uncertainty(leading_pt_outside_gap, leading_pt_outside_gap_jetmettau, leading_pt_outside_gap_jetmet, leading_pt_outside_gap_jet, merged_unc, 18, detail);
    plot_histogram(leading_pt_outside_gap, output_path_plots, label_out + "leading_pt_outside_gap",  label_out + "Uncertainty", "top_left", true);


//output the error variation
    if (detail) { cout<<"Display the merged uncertainties..."<<endl; }
    if (disp_uncertainty) { show_merged_uncertainties(merged_unc); }

//Opening the output root file
    if (detail) { cout<<"Creating " << out_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( out_uncertainty.c_str() , "RECREATE");

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
    data_jetmettau->Close();
    data_jetmet->Close();
    data_jet->Close();
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
