// Pedro Cipriano, Mar 2011
// DESY, CMS
// Last Update: 22 Mar 2013
//
// total_uncertainty()

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

void calc_total_uncertainty(TH1D *total, TH1D *jes, TH1D *model, TH1D *data, double *total_unc, int index, bool detail)
{

double max = 0.0;
double min = 0.0;
double tot = 0.0;
double ave = 0.0;

double val_lumi = 0.04;
double val_jes = 0.0;
double val_model = 0.0;
double val_stat = 0.0;
double val_total = 0.0;

for (int i = 1; i <= data->GetNbinsX();i++)
	{
	val_jes = jes->GetBinContent(i);
	val_model = model->GetBinContent(i);
	val_stat = data->GetBinError(i)/data->GetBinContent(i);
	val_total = sqrt(val_jes * val_jes + val_model * val_model + val_stat * val_stat + val_lumi * val_lumi);
	total->SetBinContent(i,val_total);
        if (val_total > max) {max = val_total;}
        tot = tot + val_total;
        if (val_total < min || i == 1) { min = val_total;}
        if (detail) { cout<< "Bin " << i << " Jes = " << val_jes << ", Model = " << val_model << ", Stat = " << val_stat << " and Total = " << val_total << endl; }
	}

//calculates the average model uncertainty
    ave = tot/total->GetNbinsX();

//displays the result
    if (detail) { cout<<"Result: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl; }

//saves the results in the results array
    total_unc[index*3+0] = ave*100;
    total_unc[index*3+1] = min*100;
    total_unc[index*3+2] = max*100;

}


void show_total_uncertainties(double *total_unc)
{
//shows the computed jes uncertainties
    cout<<" "<<endl;
    cout<<"Total Uncertainty"<<endl;
    cout<<"Observable                  Average  Minimum  Maximum"<<endl;
    cout<<"Delta phi                   "<<total_unc[0]<<"  "<<total_unc[1]<<" "<<total_unc[2]<<endl;
    cout<<"Delta phi deta1             "<<total_unc[3]<<"  "<<total_unc[4]<<"  "<<total_unc[5]<<endl;
    cout<<"Delta phi deta2             "<<total_unc[6]<<"   "<<total_unc[7]<<" "<<total_unc[8]<<endl;
    cout<<"Delta phi deta3             "<<total_unc[9]<<"  "<<total_unc[10]<<" "<<total_unc[11]<<endl;
    cout<<"Delta phi deta4             "<<total_unc[12]<<"   "<<total_unc[13]<<" "<<total_unc[14]<<endl;
    cout<<"Delta phi gap               "<<total_unc[15]<<"  "<<total_unc[16]<<" "<<total_unc[17]<<endl;
    cout<<"Delta phi deta1 gap         "<<total_unc[18]<<"  "<<total_unc[19]<<"  "<<total_unc[20]<<endl;
    cout<<"Delta phi deta2 gap         "<<total_unc[21]<<"  "<<total_unc[22]<<" "<<total_unc[23]<<endl;
    cout<<"Delta phi deta3 gap         "<<total_unc[24]<<"  "<<total_unc[25]<<" "<<total_unc[26]<<endl;
    cout<<"Delta phi deta4 gap         "<<total_unc[27]<<"  "<<total_unc[28]<<" "<<total_unc[29]<<endl;
    cout<<"Delta phi nogap             "<<total_unc[30]<<"  "<<total_unc[31]<<" "<<total_unc[32]<<endl;
    cout<<"Delta phi deta1 nogap       "<<total_unc[33]<<"  "<<total_unc[34]<<"  "<<total_unc[35]<<endl;
    cout<<"Delta phi deta2 nogap       "<<total_unc[36]<<"   "<<total_unc[37]<<"  "<<total_unc[38]<<endl;
    cout<<"Delta phi deta3 nogap       "<<total_unc[39]<<"  "<<total_unc[40]<<"  "<<total_unc[41]<<endl;
    cout<<"Delta phi deta4 nogap       "<<total_unc[42]<<"  "<<total_unc[43]<<"  "<<total_unc[44]<<endl;
    cout<<"Leading pT inside gap       "<<total_unc[45]<<"  "<<total_unc[46]<<" "<<total_unc[47]<<endl;
    cout<<"Leading eta star inside gap "<<total_unc[48]<<"  "<<total_unc[49]<<"  "<<total_unc[50]<<endl;
    cout<<"Delta eta outside gap       "<<total_unc[51]<<"  "<<total_unc[52]<<" "<<total_unc[53]<<endl;
    cout<<"Leading pT outside gap      "<<total_unc[54]<<"  "<<total_unc[55]<<" "<<total_unc[56]<<endl;
}

void total_uncertainty(string path_jes, string path_model, string path_data, string total_uncertainty, string label_in, string label_out, string output_path_plots, bool detail = false, bool disp_uncertainty = true, bool test = false)
{

//outputs the configuration
    if (detail) { cout << "Total Uncertainty Configuration"<<endl; }
    if (detail) { cout << "Input path for JES:   " << path_jes << endl; }
    if (detail) { cout << "Input path for Model: " << path_model << endl; }
    if (detail) { cout << "Input path for Stat:  " << path_data << endl; }
    if (detail) { cout << "Label In:             " << label_in << endl; }
    if (detail) { cout << "Label Out:            " << label_out << endl; }
    if (detail) { cout << "Output path:          " << total_uncertainty << endl; }
    if (detail) { cout << "Output Path Plots:    " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:         " << detail << endl; }
    if (detail) { cout << "Display Results:      " << disp_uncertainty << endl; }
    if (detail) { cout << "Test Mode:            " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data_jes = new TFile( path_jes.c_str() );
    TFile *data_model = new TFile( path_model.c_str() );
    TFile *data = new TFile( path_data.c_str() );


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
    double total_unc[19*3];

    for (int i=0; i<= 19*3-1;i++)
    {
    total_unc[i] = 0.0;
    }


//compute total uncertainty for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_jes = 0;
    TH1D *delta_phi_model = 0;
    TH1D *delta_phi_data = 0;
    TString delta_phi_name_in = label_in + "delta_phi";
    TString delta_phi_name_out = label_out + "delta_phi";

    data_jes->GetObject(delta_phi_name_in,delta_phi_jes);
    if (delta_phi_jes == 0) { cout << delta_phi_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi",delta_phi_model);
    if (delta_phi_model == 0) { cout << "merged_model_unc_delta_phi not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    
    TH1D *delta_phi;
    delta_phi =  new TH1D(delta_phi_name_out,"Total Uncertainty;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi, delta_phi_jes, delta_phi_model, delta_phi_data, total_unc, 0, detail);
    plot_histogram(delta_phi, output_path_plots, label_out + "delta_phi",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_jes = 0;
    TH1D *delta_phi_deta1_model = 0;
    TH1D *delta_phi_deta1_data = 0;
    TString delta_phi_deta1_name_in = label_in + "delta_phi_deta1";
    TString delta_phi_deta1_name_out = label_out + "delta_phi_deta1";

    data_jes->GetObject(delta_phi_deta1_name_in,delta_phi_deta1_jes);
    if (delta_phi_deta1_jes == 0) { cout << delta_phi_deta1_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1",delta_phi_deta1_model);
    if (delta_phi_deta1_model == 0) { cout << "merged_model_unc_delta_phi_deta1 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 =  new TH1D(delta_phi_deta1_name_out,"Total Uncertainty deta1;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta1, delta_phi_deta1_jes, delta_phi_deta1_model, delta_phi_deta1_data, total_unc, 1, detail);
    plot_histogram(delta_phi_deta1, output_path_plots, label_out + "delta_phi_deta1",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_jes = 0;
    TH1D *delta_phi_deta2_model = 0;
    TH1D *delta_phi_deta2_data = 0;
    TString delta_phi_deta2_name_in = label_in + "delta_phi_deta2";
    TString delta_phi_deta2_name_out = label_out + "delta_phi_deta2";

    data_jes->GetObject(delta_phi_deta2_name_in,delta_phi_deta2_jes);
    if (delta_phi_deta2_jes == 0) { cout << delta_phi_deta2_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2",delta_phi_deta2_model);
    if (delta_phi_deta2_model == 0) { cout << "merged_model_unc_delta_phi_deta2 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 =  new TH1D(delta_phi_deta2_name_out,"Total Uncertainty deta2;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta2, delta_phi_deta2_jes, delta_phi_deta2_model, delta_phi_deta2_data, total_unc, 2, detail);
    plot_histogram(delta_phi_deta2, output_path_plots, label_out + "delta_phi_deta2",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_jes = 0;
    TH1D *delta_phi_deta3_model = 0;
    TH1D *delta_phi_deta3_data = 0;
    TString delta_phi_deta3_name_in = label_in + "delta_phi_deta3";
    TString delta_phi_deta3_name_out = label_out + "delta_phi_deta3";

    data_jes->GetObject(delta_phi_deta3_name_in,delta_phi_deta3_jes);
    if (delta_phi_deta3_jes == 0) { cout << delta_phi_deta3_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3",delta_phi_deta3_model);
    if (delta_phi_deta3_model == 0) { cout << "merged_model_unc_delta_phi_deta3 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 =  new TH1D(delta_phi_deta3_name_out,"Total Uncertainty deta3;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta3, delta_phi_deta3_jes, delta_phi_deta3_model, delta_phi_deta3_data, total_unc, 3, detail);
    plot_histogram(delta_phi_deta3, output_path_plots, label_out + "delta_phi_deta3",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_jes = 0;
    TH1D *delta_phi_deta4_model = 0;
    TH1D *delta_phi_deta4_data = 0;
    TString delta_phi_deta4_name_in = label_in + "delta_phi_deta4";
    TString delta_phi_deta4_name_out = label_out + "delta_phi_deta4";

    data_jes->GetObject(delta_phi_deta4_name_in,delta_phi_deta4_jes);
    if (delta_phi_deta4_jes == 0) { cout << delta_phi_deta4_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4",delta_phi_deta4_model);
    if (delta_phi_deta4_model == 0) { cout << "merged_model_unc_delta_phi_deta4 not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 =  new TH1D(delta_phi_deta4_name_out,"Total Uncertainty deta4;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta4, delta_phi_deta4_jes, delta_phi_deta4_model, delta_phi_deta4_data, total_unc, 4, detail);
    plot_histogram(delta_phi_deta4, output_path_plots, label_out + "delta_phi_deta4",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi gap distribution
    if (detail) { cout<<"Delta phi gap"<<endl; }

    TH1D *delta_phi_gap_jes = 0;
    TH1D *delta_phi_gap_model = 0;
    TH1D *delta_phi_gap_data = 0;
    TString delta_phi_gap_name_in = label_in + "delta_phi_gap";
    TString delta_phi_gap_name_out = label_out + "delta_phi_gap";

    data_jes->GetObject(delta_phi_gap_name_in,delta_phi_gap_jes);
    if (delta_phi_gap_jes == 0) { cout << delta_phi_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_gap",delta_phi_gap_model);
    if (delta_phi_gap_model == 0) { cout << "merged_model_unc_delta_phi_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D(delta_phi_gap_name_out,"Total Uncertainty gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_gap, delta_phi_gap_jes, delta_phi_gap_model, delta_phi_gap_data, total_unc, 5, detail);
    plot_histogram(delta_phi_gap, output_path_plots, label_out + "delta_phi_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_jes = 0;
    TH1D *delta_phi_deta1_gap_model = 0;
    TH1D *delta_phi_deta1_gap_data = 0;
    TString delta_phi_deta1_gap_name_in = label_in + "delta_phi_deta1_gap";
    TString delta_phi_deta1_gap_name_out = label_out + "delta_phi_deta1_gap";

    data_jes->GetObject(delta_phi_deta1_gap_name_in,delta_phi_deta1_gap_jes);
    if (delta_phi_deta1_gap_jes == 0) { cout << delta_phi_deta1_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_gap",delta_phi_deta1_gap_model);
    if (delta_phi_deta1_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta1_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap =  new TH1D(delta_phi_deta1_gap_name_out,"Total Uncertainty deta1 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta1_gap, delta_phi_deta1_gap_jes, delta_phi_deta1_gap_model, delta_phi_deta1_gap_data, total_unc, 6, detail);
    plot_histogram(delta_phi_deta1_gap, output_path_plots, label_out + "delta_phi_deta1_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_jes = 0;
    TH1D *delta_phi_deta2_gap_model = 0;
    TH1D *delta_phi_deta2_gap_data = 0;
    TString delta_phi_deta2_gap_name_in = label_in + "delta_phi_deta2_gap";
    TString delta_phi_deta2_gap_name_out = label_out + "delta_phi_deta2_gap";

    data_jes->GetObject(delta_phi_deta2_gap_name_in,delta_phi_deta2_gap_jes);
    if (delta_phi_deta2_gap_jes == 0) { cout << delta_phi_deta2_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_gap",delta_phi_deta2_gap_model);
    if (delta_phi_deta2_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta2_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap =  new TH1D(delta_phi_deta2_gap_name_out,"Total Uncertainty deta2 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta2_gap, delta_phi_deta2_gap_jes, delta_phi_deta2_gap_model, delta_phi_deta2_gap_data, total_unc, 7, detail);
    plot_histogram(delta_phi_deta2_gap, output_path_plots, label_out + "delta_phi_deta2_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_jes = 0;
    TH1D *delta_phi_deta3_gap_model = 0;
    TH1D *delta_phi_deta3_gap_data = 0;
    TString delta_phi_deta3_gap_name_in = label_in + "delta_phi_deta3_gap";
    TString delta_phi_deta3_gap_name_out = label_out + "delta_phi_deta3_gap";

    data_jes->GetObject(delta_phi_deta3_gap_name_in,delta_phi_deta3_gap_jes);
    if (delta_phi_deta3_gap_jes == 0) { cout << delta_phi_deta3_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_gap",delta_phi_deta3_gap_model);
    if (delta_phi_deta3_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta3_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap =  new TH1D(delta_phi_deta3_gap_name_out,"Total Uncertainty deta3 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta3_gap, delta_phi_deta3_gap_jes, delta_phi_deta3_gap_model, delta_phi_deta3_gap_data, total_unc, 8, detail);
    plot_histogram(delta_phi_deta3_gap, output_path_plots, label_out + "delta_phi_deta3_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_jes = 0;
    TH1D *delta_phi_deta4_gap_model = 0;
    TH1D *delta_phi_deta4_gap_data = 0;
    TString delta_phi_deta4_gap_name_in = label_in + "delta_phi_deta4_gap";
    TString delta_phi_deta4_gap_name_out = label_out + "delta_phi_deta4_gap";

    data_jes->GetObject(delta_phi_deta4_gap_name_in,delta_phi_deta4_gap_jes);
    if (delta_phi_deta4_gap_jes == 0) { cout << delta_phi_deta4_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_gap",delta_phi_deta4_gap_model);
    if (delta_phi_deta4_gap_model == 0) { cout << "merged_model_unc_delta_phi_deta4_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap =  new TH1D(delta_phi_deta4_gap_name_out,"Total Uncertainty deta4 gap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta4_gap, delta_phi_deta4_gap_jes, delta_phi_deta4_gap_model, delta_phi_deta4_gap_data, total_unc, 9, detail);
    plot_histogram(delta_phi_deta4_gap, output_path_plots, label_out + "delta_phi_deta4_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi nogap distribution
    if (detail) { cout<<"Delta phi nogap"<<endl; }

    TH1D *delta_phi_nogap_jes = 0;
    TH1D *delta_phi_nogap_model = 0;
    TH1D *delta_phi_nogap_data = 0;
    TString delta_phi_nogap_name_in = label_in + "delta_phi_nogap";
    TString delta_phi_nogap_name_out = label_out + "delta_phi_nogap";

    data_jes->GetObject(delta_phi_nogap_name_in,delta_phi_nogap_jes);
    if (delta_phi_nogap_jes == 0) { cout << delta_phi_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_nogap",delta_phi_nogap_model);
    if (delta_phi_nogap_model == 0) { cout << "merged_model_unc_delta_phi_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D(delta_phi_nogap_name_out,"Total Uncertainty nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_nogap, delta_phi_nogap_jes, delta_phi_nogap_model, delta_phi_nogap_data, total_unc, 10, detail);
    plot_histogram(delta_phi_nogap, output_path_plots, label_out + "delta_phi_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_jes = 0;
    TH1D *delta_phi_deta1_nogap_model = 0;
    TH1D *delta_phi_deta1_nogap_data = 0;
    TString delta_phi_deta1_nogap_name_in = label_in + "delta_phi_deta1_nogap";
    TString delta_phi_deta1_nogap_name_out = label_out + "delta_phi_deta1_nogap";

    data_jes->GetObject(delta_phi_deta1_nogap_name_in,delta_phi_deta1_nogap_jes);
    if (delta_phi_deta1_nogap_jes == 0) { cout << delta_phi_deta1_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta1_nogap",delta_phi_deta1_nogap_model);
    if (delta_phi_deta1_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta1_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap =  new TH1D(delta_phi_deta1_nogap_name_out,"Total Uncertainty deta1 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jes, delta_phi_deta1_nogap_model, delta_phi_deta1_nogap_data, total_unc, 11, detail);
    plot_histogram(delta_phi_deta1_nogap, output_path_plots, label_out + "delta_phi_deta1_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_jes = 0;
    TH1D *delta_phi_deta2_nogap_model = 0;
    TH1D *delta_phi_deta2_nogap_data = 0;
    TString delta_phi_deta2_nogap_name_in = label_in + "delta_phi_deta2_nogap";
    TString delta_phi_deta2_nogap_name_out = label_out + "delta_phi_deta2_nogap";

    data_jes->GetObject(delta_phi_deta2_nogap_name_in,delta_phi_deta2_nogap_jes);
    if (delta_phi_deta2_nogap_jes == 0) { cout << delta_phi_deta2_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta2_nogap",delta_phi_deta2_nogap_model);
    if (delta_phi_deta2_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta2_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap =  new TH1D(delta_phi_deta2_nogap_name_out,"Total Uncertainty deta2 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jes, delta_phi_deta2_nogap_model, delta_phi_deta2_nogap_data, total_unc, 12, detail);
    plot_histogram(delta_phi_deta2_nogap, output_path_plots, label_out + "delta_phi_deta2_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_jes = 0;
    TH1D *delta_phi_deta3_nogap_model = 0;
    TH1D *delta_phi_deta3_nogap_data = 0;
    TString delta_phi_deta3_nogap_name_in = label_in + "delta_phi_deta3_nogap";
    TString delta_phi_deta3_nogap_name_out = label_out + "delta_phi_deta3_nogap";

    data_jes->GetObject(delta_phi_deta3_nogap_name_in,delta_phi_deta3_nogap_jes);
    if (delta_phi_deta3_nogap_jes == 0) { cout << delta_phi_deta3_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta3_nogap",delta_phi_deta3_nogap_model);
    if (delta_phi_deta3_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta3_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap =  new TH1D(delta_phi_deta3_nogap_name_out,"Total Uncertainty deta3 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jes, delta_phi_deta3_nogap_model, delta_phi_deta3_nogap_data, total_unc, 13, detail);
    plot_histogram(delta_phi_deta3_nogap, output_path_plots, label_out + "delta_phi_deta3_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_jes = 0;
    TH1D *delta_phi_deta4_nogap_model = 0;
    TH1D *delta_phi_deta4_nogap_data = 0;
    TString delta_phi_deta4_nogap_name_in = label_in + "delta_phi_deta4_nogap";
    TString delta_phi_deta4_nogap_name_out = label_out + "delta_phi_deta4_nogap";

    data_jes->GetObject(delta_phi_deta4_nogap_name_in,delta_phi_deta4_nogap_jes);
    if (delta_phi_deta4_nogap_jes == 0) { cout << delta_phi_deta4_nogap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_phi_deta4_nogap",delta_phi_deta4_nogap_model);
    if (delta_phi_deta4_nogap_model == 0) { cout << "merged_model_unc_delta_phi_deta4_nogap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap =  new TH1D(delta_phi_deta4_nogap_name_out,"Total Uncertainty deta4 nogap;#Delta#phi;Total Uncertainty", dphi_nbins, dphi_bins);

    calc_total_uncertainty(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jes, delta_phi_deta4_nogap_model, delta_phi_deta4_nogap_data, total_unc, 14, detail);
    plot_histogram(delta_phi_deta4_nogap, output_path_plots, label_out + "delta_phi_deta4_nogap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for leading pt inside gap distribution
    if (detail) { cout<<"Leading pt inside gap"<<endl; }

    TH1D *leading_pt_inside_gap_jes = 0;
    TH1D *leading_pt_inside_gap_model = 0;
    TH1D *leading_pt_inside_gap_data = 0;
    TString leading_pt_inside_gap_name_in = label_in + "leading_pt_inside_gap";
    TString leading_pt_inside_gap_name_out = label_out + "leading_pt_inside_gap";

    data_jes->GetObject(leading_pt_inside_gap_name_in,leading_pt_inside_gap_jes);
    if (leading_pt_inside_gap_jes == 0) { cout << leading_pt_inside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_inside_gap",leading_pt_inside_gap_model);
    if (leading_pt_inside_gap_model == 0) { cout << "merged_model_unc_leading_pt_inside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap =  new TH1D(leading_pt_inside_gap_name_out,"Total Uncertainty pt leading inside gap;p_{T};Total Uncertainty", in_nbins, in_bins);

    calc_total_uncertainty(leading_pt_inside_gap, leading_pt_inside_gap_jes, leading_pt_inside_gap_model, leading_pt_inside_gap_data, total_unc, 15, detail);
    plot_histogram(leading_pt_inside_gap, output_path_plots, label_out + "leading_pt_inside_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for leading eta star inside gap distribution
    if (detail) { cout<<"Leading eta star inside gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_jes = 0;
    TH1D *leading_eta_star_inside_gap_model = 0;
    TH1D *leading_eta_star_inside_gap_data = 0;
    TString leading_eta_star_inside_gap_name_in = label_in + "leading_eta_star_inside_gap";
    TString leading_eta_star_inside_gap_name_out = label_out + "leading_eta_star_inside_gap";

    data_jes->GetObject(leading_eta_star_inside_gap_name_in,leading_eta_star_inside_gap_jes);
    if (leading_eta_star_inside_gap_jes == 0) { cout << leading_eta_star_inside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_eta_star_inside_gap",leading_eta_star_inside_gap_model);
    if (leading_eta_star_inside_gap_model == 0) { cout << "merged_model_unc_leading_eta_star_inside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap =  new TH1D(leading_eta_star_inside_gap_name_out,"Total Uncertainty eta* leading inside gap;#eta*;Total Uncertainty", etastar_nbins, etastar_bins);

    calc_total_uncertainty(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jes, leading_eta_star_inside_gap_model, leading_eta_star_inside_gap_data, total_unc, 16, detail);
    plot_histogram(leading_eta_star_inside_gap, output_path_plots, label_out + "leading_eta_star_inside_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for delta eta outside gap distribution
    if (detail) { cout<<"Delta eta outside gap"<<endl; }

    TH1D *delta_eta_outside_gap_jes = 0;
    TH1D *delta_eta_outside_gap_model = 0;
    TH1D *delta_eta_outside_gap_data = 0;
    TString delta_eta_outside_gap_name_in = label_in + "delta_eta_outside_gap";
    TString delta_eta_outside_gap_name_out = label_out + "delta_eta_outside_gap";

    data_jes->GetObject(delta_eta_outside_gap_name_in,delta_eta_outside_gap_jes);
    if (delta_eta_outside_gap_jes == 0) { cout << delta_eta_outside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_delta_eta_outside_gap",delta_eta_outside_gap_model);
    if (delta_eta_outside_gap_model == 0) { cout << "merged_model_unc_delta_eta_outside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap =  new TH1D(delta_eta_outside_gap_name_out,"Total Uncertainty #Delta#eta leading outside gap;#Delta#eta;Total Uncertainty", deta_out_nbins, deta_out_bins);

    calc_total_uncertainty(delta_eta_outside_gap, delta_eta_outside_gap_jes, delta_eta_outside_gap_model, delta_eta_outside_gap_data, total_unc, 17, detail);
    plot_histogram(delta_eta_outside_gap, output_path_plots, label_out + "delta_eta_outside_gap",  label_out + "Uncertainty ", "top_right", true);


//compute total uncertainty for leading pt outside gap distribution
    if (detail) { cout<<"Leading pt outside gap"<<endl; }

    TH1D *leading_pt_outside_gap_jes = 0;
    TH1D *leading_pt_outside_gap_model = 0;
    TH1D *leading_pt_outside_gap_data = 0;
    TString leading_pt_outside_gap_name_in = label_in + "leading_pt_outside_gap";
    TString leading_pt_outside_gap_name_out = label_out + "leading_pt_outside_gap";

    data_jes->GetObject(leading_pt_outside_gap_name_in,leading_pt_outside_gap_jes);
    if (leading_pt_outside_gap_jes == 0) { cout << leading_pt_outside_gap_name_in << " not found!" << endl; return; }
    data_model->GetObject("merged_model_unc_leading_pt_outside_gap",leading_pt_outside_gap_model);
    if (leading_pt_outside_gap_model == 0) { cout << "merged_model_unc_leading_pt_outside_gap not found!" << endl; return; }
    data->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap =  new TH1D(leading_pt_outside_gap_name_out,"Total Uncertainty pt leading outside gap;p_{T};Total Uncertainty", out_nbins, out_bins);

    calc_total_uncertainty(leading_pt_outside_gap, leading_pt_outside_gap_jes, leading_pt_outside_gap_model, leading_pt_outside_gap_data, total_unc, 18, detail);
    plot_histogram(leading_pt_outside_gap, output_path_plots, label_out + "leading_pt_outside_gap",  label_out + "Uncertainty ", "top_right", true);


//output the error variation
    if (detail) { cout<<"Display the total uncertainties..."<<endl; }
    if (disp_uncertainty) { show_total_uncertainties(total_unc); }


//Opening the output root file
    if (detail) { cout<<"Creating " << total_uncertainty << "..."<<endl; }
    TFile *data_output = TFile::Open( total_uncertainty.c_str() , "RECREATE");

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
    data_jes->Close();
    data_model->Close();
    data->Close();
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
