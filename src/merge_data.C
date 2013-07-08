// Pedro Cipriano, Mar 2012
// DESY, CMS
// Last Update: 20 Nov 2012
//
// merge_data(output_dir,detail,disp_errors)
// merge the different data datasets in one root file and plots the control distributions

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "common_methods.h"

void show_error_variation(double *final_errors)
{
//print the error variation on the screen

    cout<<" "<<endl;
    cout<<setiosflags(ios::fixed)<<setprecision(3)<<setfill('0'); //setprecision(3);
    cout<<"Merging data errors"<<endl;
    cout<<"Observable                Merged  ->  JetMETTau  JetMET  Jet"<<endl;
    cout<<"Delta Phi                 "<<setw(6)<<final_errors[0]<<"  ->  "<<setw(6)<<final_errors[1]<<"     "<<setw(6)<<final_errors[2]<<"  "<<setw(6)<<final_errors[3]<<endl;
    cout<<"Delta Phi Gap             "<<setw(6)<<final_errors[4]<<"  ->  "<<setw(6)<<final_errors[5]<<"     "<<setw(6)<<final_errors[6]<<"  "<<setw(6)<<final_errors[7]<<endl; 
    cout<<"Delta Phi Nogap           "<<setw(6)<<final_errors[8]<<"  ->  "<<setw(6)<<final_errors[9]<<"     "<<setw(6)<<final_errors[10]<<"  "<<setw(6)<<final_errors[11]<<endl;
    cout<<"Leading pT inside gap     "<<setw(6)<<final_errors[12]<<"  ->  "<<setw(6)<<final_errors[13]<<"     "<<setw(6)<<final_errors[14]<<"  "<<setw(6)<<final_errors[15]<<endl;
    cout<<"Leading pT outside gap    "<<setw(6)<<final_errors[16]<<"  ->  "<<setw(6)<<final_errors[17]<<"     "<<setw(6)<<final_errors[18]<<"  "<<setw(6)<<final_errors[19]<<endl;
    cout<<"Delta Phi Deta1           "<<setw(6)<<final_errors[20]<<"  ->  "<<setw(6)<<final_errors[21]<<"     "<<setw(6)<<final_errors[22]<<"  "<<setw(6)<<final_errors[23]<<endl;
    cout<<"Delta Phi Deta2           "<<setw(6)<<final_errors[24]<<"  ->  "<<setw(6)<<final_errors[25]<<"     "<<setw(6)<<final_errors[26]<<"  "<<setw(6)<<final_errors[27]<<endl;
    cout<<"Delta Phi Deta3           "<<setw(6)<<final_errors[28]<<"  ->  "<<setw(6)<<final_errors[29]<<"     "<<setw(6)<<final_errors[30]<<"  "<<setw(6)<<final_errors[31]<<endl;
    cout<<"Delta Phi Deta4           "<<setw(6)<<final_errors[32]<<"  ->  "<<setw(6)<<final_errors[33]<<"     "<<setw(6)<<final_errors[34]<<"  "<<setw(6)<<final_errors[35]<<endl;
    cout<<"Delta Phi Gap Deta1       "<<setw(6)<<final_errors[36]<<"  ->  "<<setw(6)<<final_errors[37]<<"     "<<setw(6)<<final_errors[38]<<"  "<<setw(6)<<final_errors[39]<<endl;
    cout<<"Delta Phi Gap Deta2       "<<setw(6)<<final_errors[40]<<"  ->  "<<setw(6)<<final_errors[41]<<"     "<<setw(6)<<final_errors[42]<<"  "<<setw(6)<<final_errors[43]<<endl;
    cout<<"Delta Phi Gap Deta3       "<<setw(6)<<final_errors[44]<<"  ->  "<<setw(6)<<final_errors[45]<<"     "<<setw(6)<<final_errors[46]<<"  "<<setw(6)<<final_errors[47]<<endl;
    cout<<"Delta Phi Gap Deta4       "<<setw(6)<<final_errors[48]<<"  ->  "<<setw(6)<<final_errors[49]<<"     "<<setw(6)<<final_errors[50]<<"  "<<setw(6)<<final_errors[51]<<endl;
    cout<<"Delta Phi noGap Deta1     "<<setw(6)<<final_errors[52]<<"  ->  "<<setw(6)<<final_errors[53]<<"     "<<setw(6)<<final_errors[54]<<"  "<<setw(6)<<final_errors[55]<<endl;
    cout<<"Delta Phi noGap Deta2     "<<setw(6)<<final_errors[56]<<"  ->  "<<setw(6)<<final_errors[57]<<"     "<<setw(6)<<final_errors[58]<<"  "<<setw(6)<<final_errors[59]<<endl;
    cout<<"Delta Phi noGap Deta3     "<<setw(6)<<final_errors[60]<<"  ->  "<<setw(6)<<final_errors[61]<<"     "<<setw(6)<<final_errors[62]<<"  "<<setw(6)<<final_errors[63]<<endl;
    cout<<"Delta Phi noGap Deta4     "<<setw(6)<<final_errors[64]<<"  ->  "<<setw(6)<<final_errors[65]<<"     "<<setw(6)<<final_errors[66]<<"  "<<setw(6)<<final_errors[67]<<endl;
    cout<<"Leading Forward pT        "<<setw(6)<<final_errors[68]<<"  ->  "<<setw(6)<<final_errors[69]<<"     "<<setw(6)<<final_errors[70]<<"  "<<setw(6)<<final_errors[71]<<endl;
    cout<<"Leading Central pT        "<<setw(6)<<final_errors[72]<<"  ->  "<<setw(6)<<final_errors[73]<<"     "<<setw(6)<<final_errors[74]<<"  "<<setw(6)<<final_errors[75]<<endl;
    cout<<"Leading eta* inside gap   "<<setw(6)<<final_errors[76]<<"  ->  "<<setw(6)<<final_errors[77]<<"     "<<setw(6)<<final_errors[78]<<"  "<<setw(6)<<final_errors[79]<<endl;
    cout<<"Delta eta outside gap     "<<setw(6)<<final_errors[80]<<"  ->  "<<setw(6)<<final_errors[81]<<"     "<<setw(6)<<final_errors[82]<<"  "<<setw(6)<<final_errors[83]<<endl;
    cout<<"Leading Forward Eta       "<<setw(6)<<final_errors[84]<<"  ->  "<<setw(6)<<final_errors[85]<<"     "<<setw(6)<<final_errors[86]<<"  "<<setw(6)<<final_errors[87]<<endl;
    cout<<"Leading Central Eta       "<<setw(6)<<final_errors[88]<<"  ->  "<<setw(6)<<final_errors[89]<<"     "<<setw(6)<<final_errors[90]<<"  "<<setw(6)<<final_errors[91]<<endl;
    cout<<"Delta eta                 "<<setw(6)<<final_errors[92]<<"  ->  "<<setw(6)<<final_errors[93]<<"     "<<setw(6)<<final_errors[94]<<"  "<<setw(6)<<final_errors[95]<<endl;
    cout<<"Delta eta Gap             "<<setw(6)<<final_errors[96]<<"  ->  "<<setw(6)<<final_errors[97]<<"     "<<setw(6)<<final_errors[98]<<"  "<<setw(6)<<final_errors[99]<<endl;
    cout<<"Delta eta noGap           "<<setw(6)<<final_errors[100]<<"  ->  "<<setw(6)<<final_errors[101]<<"     "<<setw(6)<<final_errors[102]<<"  "<<setw(6)<<final_errors[103]<<endl;
    cout<<"Vertex Selected           "<<setw(6)<<final_errors[104]<<"  ->  "<<setw(6)<<final_errors[105]<<"     "<<setw(6)<<final_errors[106]<<"  "<<setw(6)<<final_errors[107]<<endl;
    cout<<"Delta Phi Fine            "<<setw(6)<<final_errors[108]<<"  ->  "<<setw(6)<<final_errors[109]<<"     "<<setw(6)<<final_errors[110]<<"  "<<setw(6)<<final_errors[111]<<endl;
    cout<<"Leading pT                "<<setw(6)<<final_errors[112]<<"  ->  "<<setw(6)<<final_errors[113]<<"     "<<setw(6)<<final_errors[114]<<"  "<<setw(6)<<final_errors[115]<<endl;
    cout<<"Leading pT Fine           "<<setw(6)<<final_errors[116]<<"  ->  "<<setw(6)<<final_errors[117]<<"     "<<setw(6)<<final_errors[118]<<"  "<<setw(6)<<final_errors[119]<<endl;
    cout<<"Leading Central pT Fine   "<<setw(6)<<final_errors[120]<<"  ->  "<<setw(6)<<final_errors[121]<<"     "<<setw(6)<<final_errors[122]<<"  "<<setw(6)<<final_errors[123]<<endl;
    cout<<"Leading Forward pT Fine   "<<setw(6)<<final_errors[124]<<"  ->  "<<setw(6)<<final_errors[125]<<"     "<<setw(6)<<final_errors[126]<<"  "<<setw(6)<<final_errors[127]<<endl;
    cout<<"Primary Vertex Z-position "<<setw(6)<<final_errors[128]<<"  ->  "<<setw(6)<<final_errors[129]<<"     "<<setw(6)<<final_errors[130]<<"  "<<setw(6)<<final_errors[131]<<endl;

}

void plot_histogram(TH1D *merged, TH1D *jetmettau, TH1D *jetmet, TH1D *jet, string path, string fileout, string legend_position = "top_left", bool logscale = true, bool detail = false)
{
//plots the different datasets and the merged result

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    if (logscale) { gPad->SetLogy(); }
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 0.0;
    double max = merged->GetMaximum();
    if (merged->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(merged,detail);
    }
    else
    {
    min = merged->GetMinimum();
    }
    
    set_histogram_min_max(jetmettau, min, max, detail);
    set_histogram_min_max(jetmet, min, max, detail);
    //set_histogram_min_max(jet, min, max, detail);
    
    if (logscale)
	{
	max = 1.3 * max;
    	min = 0.7 * min;
    	}
    else
	{
	max = 0.7 * max;
    	min = 0;
	}
    
//plooting
    merged->SetMaximum(max);
    merged->SetMinimum(min);
    merged->SetFillColor(5);
    merged->Draw("e2");
    jetmettau->SetLineWidth(4);
    jetmettau->SetLineColor(2);
    jetmettau->SetLineStyle(1);
    jetmettau->Draw("e1same");
    jetmet->SetLineWidth(4);
    jetmet->SetLineColor(3);
    jetmet->SetLineStyle(2);
    jetmet->Draw("e1same");
    jet->SetLineWidth(4);
    jet->SetLineColor(4);
    jet->SetLineStyle(4);
    //jet->Draw("e1same");
    
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(jetmettau,"JetMETTau_2010A - L = 284 nb^{-1}","l");
    leg01->AddEntry(jetmet,"JetMET_2010A - L = 2.90 pb^{-1}","l");
  //  leg01->AddEntry(jet,"Jet_2010B - L = 32.1 pb^{-1}","l");
    leg01->AddEntry(merged,"Merged Data - L = 3.2 pb^{-1}","f");
  //  leg01->AddEntry(merged,"Merged Data - L = 35.2 pb^{-1}","f");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}

void calc_nomalization(TH1D *jetmettau, TH1D *jetmet, TH1D *jet, double *weigth, bool detail)
{
//computes the merging factors

//declares the variables needed for the calculation
    double entries[3] = {0.,0.,0.};
    double total_entries = 0.0;
    double eff_entries[3] = {0.,0.,0.};
    double eff_total_entries = 0.;
    double integral[3] = {0.,0.,0.};
    
//reset the weigth and scale arrays
    weigth[0] = 1.0;
    weigth[1] = 1.0;
    weigth[2] = 1.0;

//calculate the scale for the first dataset
    entries[0] = jetmettau->GetEntries();
    if (entries[0] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[0] << endl; } entries[0] = 1.0; }
    eff_entries[0] = jetmettau->GetEffectiveEntries();
    integral[0] = jetmettau->Integral();
    if (integral[0] <= 0) { if (detail) { cout<<"Bad integral = " << integral[0] << endl; } integral[0] = 1.0; }
    if (detail) { cout<<"JetMETTau_2010A "<<setw(8)<<entries[0]<<endl; }

//calculate the scale for the second dataset
    entries[1] = jetmet->GetEntries();
    if (entries[1] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[1] << endl; } entries[1] = 1.0; }
    eff_entries[1] = jetmet->GetEffectiveEntries();
    integral[1] = jetmet->Integral();
    if (integral[1] <= 0) { if (detail) { cout<<"Bad integral = " << integral[1] << endl; } integral[1] = 1.0; }
    if (detail) { cout<<"JetMET_2010A    "<<setw(8)<<entries[1]<<endl; }

//calculate the scale for the third dataset
    entries[2] = jet->GetEntries();
    if (entries[2] <= 0) { if (detail) { cout<<"Bad number of entries = " << entries[2] << endl; } entries[2] = 1.0; }
    eff_entries[2] = jet->GetEffectiveEntries();
    integral[2] = jet->Integral();
    if (integral[2] <= 0) { if (detail) { cout<<"Bad integral = " << integral[2] << endl; } integral[2] = 1.0; }
    if (detail) { cout<<"Jet_2010B       "<<setw(8)<<entries[2]<<endl; }

//calculate the total number of entries
    //total_entries = entries[0] + entries[1] + entries[2];
    total_entries = entries[0] + entries[1]; //old version when we used only 2 datasets
    //eff_total_entries = eff_entries[0] + eff_entries[1] + eff_entries[2];
    eff_total_entries = eff_entries[0] + eff_entries[1]; //old version when we used only 2 datasets
    if (total_entries <= 0) { if (detail) { cout<<"Bad number of total entries = " << total_entries << endl; } total_entries = 1.0; }

//a few control outputs
    if (detail)
	{
	cout<<"Integral   "<<setw(6)<<integral[0]<<"    "<<setw(6)<<integral[1]<<" "<<setw(6)<<integral[2]<<endl;
	cout<<"Entries    "<<setw(6)<<entries[0]<<"     "<<setw(6)<<entries[1]<<"  "<<setw(6)<<entries[2]<<"  "<<setw(6)<<total_entries<<endl;
	}

//calculate the relative weight    
    weigth[0] = entries[0]/total_entries;
    weigth[1] = entries[1]/total_entries;
  //  weigth[2] = entries[2]/total_entries;
    if (detail) { cout<<"Weigths "<<weigth[0]<<" "<<weigth[1]<<" "<<weigth[2]<<endl; }

}

void calc_error(TH1D *merged, TH1D *jetmettau, TH1D *jetmet, TH1D *jet, double *final_errors, int index, bool detail)
{
//evaluate the statistical errors on the merged samples and the individual samples

if (detail) { cout<<"Evaluating statistical errors ..."<<endl; }
//variable declaration
double error_merged = 0.0;
double error[3] = {0.,0.,0.};

//setting the number of bins
int tot_merged = merged->GetNbinsX();
int tot_jetmettau = jetmettau->GetNbinsX();
int tot_jetmet = jetmet->GetNbinsX();
//int tot_jet = jet->GetNbinsX();

//remove the empty bins from the bin count
for (int i=1; i<= merged->GetNbinsX();i++)
    {
	if (merged->GetBinContent(i) == 0) { tot_merged = tot_merged - 1; }
	if (jetmettau->GetBinContent(i) == 0) { tot_jetmettau = tot_jetmettau - 1; }
	if (jetmet->GetBinContent(i) == 0) { tot_jetmet = tot_jetmet - 1; }
	// if (jet->GetBinContent(i) == 0) { tot_jet = tot_jet - 1; }
    }

//loop over the bins to get and average the statistical error
    for (int i=1; i<= merged->GetNbinsX();i++)
    {
	if (merged->GetBinContent(i) > 0) { error_merged = error_merged + (merged->GetBinError(i)/merged->GetBinContent(i))/tot_merged; }
	if (jetmettau->GetBinContent(i) > 0) { error[0] = error[0] + (jetmettau->GetBinError(i)/jetmettau->GetBinContent(i))/tot_jetmettau; }
	if (jetmet->GetBinContent(i) > 0) { error[1] = error[1] + (jetmet->GetBinError(i)/jetmet->GetBinContent(i))/tot_jetmet; }
	// if (jet->GetBinContent(i) > 0) { error[2] = error[2] + (jet->GetBinError(i)/jet->GetBinContent(i))/tot_jet; }
    }
    if (detail) { cout<<"Errors          "<<error[0]<<" "<<error[1]<<" "<<error[2]<<" -> "<<error_merged<<endl; }

//saving the errors on final_errors array
    final_errors[index*4+0] = error_merged * 100;
    final_errors[index*4+1] = error[0] * 100;
    final_errors[index*4+2] = error[1] * 100;
  //  final_errors[index*4+3] = error[2] * 100;
}


void merge(TH1D *merged, TH1D *jetmettau, TH1D *jetmet, TH1D *jet, double weigth[3])
{
//merges the different datasets

    merged->Add(jetmettau,weigth[0]);
    merged->Add(jetmet,weigth[1]);
  //  merged->Add(jet,weigth[2]);
}

void merge_data(string path_jetmettau, string path_jetmet, string path_jet, string output_path_plots = "../output/merge_data/", string output_rootfile = "../output/histograms/xsec_data_2010.root", bool detail = false, bool disp_errors = true, bool test = false)
{
//main merging routine
//output detail: true-output all the major steps; false-run quietly
//display errors: true-display merging final errors; false-dont show anything

   if (detail) { cout<<"Merge Data Configuration"<<endl; }
   if (detail) { cout<<"Input Path :       "<<path_jetmettau<<endl; }
   if (detail) { cout<<"Input Path :       "<<path_jetmet<<endl; }
   if (detail) { cout<<"Input Path :       "<<path_jet<<endl; }
   if (detail) { cout<<"Plot Output Path : "<<output_path_plots<<endl; }
   if (detail) { cout<<"Root Output Path : "<<output_rootfile<<endl; }
   if (detail) { cout<<"Detail :           "<<detail<<endl; }
   if (detail) { cout<<"Disp Errors :      "<<disp_errors<<endl; }
   if (detail) { cout<<"Test Mode :        "<<test<<endl; }

//opening the input and output data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *data_jetmettau = new TFile( path_jetmettau.c_str() );
    TFile *data_jetmet = new TFile( path_jetmet.c_str() );
    TFile *data_jet = new TFile( path_jet.c_str() );

//setting the labels for the ploting
    TString label1 = "JetMETTau2010A - L = 284 nb^{-1}";
    TString label2 = "JetMET2010A - L = 2.90 pb^{-1}";
    TString label3 = "Jet2010B - L = 32.1 pb^{-1}";
    TString label4 = "Merged Data - L = 35.2 pb^{-1}";

//declare histogram bins
//int cent_nbins = 7;
//int forw_nbins = 7;
//int all_nbins = 11;

//double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
//double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
//double all_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int in_nbins = 9;
int out_nbins = 9;

double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

//int deta_nbins = 4;
int dphi_nbins = 7;

//double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

//int etac_nbins = 6;
//double etac_bins[7] = {-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8};

//int etaf_nbins = 7;
//double etaf_bins[8] = {-4.7,-4.2,-3.7,-3.2,3.2,3.7,4.2,4.7};

//int eta_nbins = 14;
//double eta_bins[15] = {-4.7,-4.2,-3.7,-3.2,-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8,3.2,3.7,4.2,4.7};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//declaring and setting variables
double weigth[3] = {0.,0.,0.};
double final_errors[33*4];

    for (int i=0; i<= 33*4-1;i++)
    {
    final_errors[i] = 0.0;
    }

//start the merging of the data

//merging the delta phi distribution
    TH1D *delta_phi_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi");
    TH1D *delta_phi_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi");
    TH1D *delta_phi_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi");
    
    TH1D *delta_phi;
    delta_phi = new TH1D("ak5PF_delta_phi","#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi"<<endl; }
    calc_nomalization(delta_phi_jetmettau, delta_phi_jetmet, delta_phi_jet, weigth, detail);
    merge(delta_phi, delta_phi_jetmettau, delta_phi_jetmet, delta_phi_jet, weigth);    
    plot_histogram(delta_phi, delta_phi_jetmettau, delta_phi_jetmet, delta_phi_jet, output_path_plots, "merging_delta_phi", "top_left", true, detail);
    plot_histogram(delta_phi, delta_phi_jetmettau, delta_phi_jetmet, delta_phi_jet, output_path_plots, "merging_delta_phi_linear", "top_left", false, detail);
    calc_error(delta_phi, delta_phi_jetmettau, delta_phi_jetmet, delta_phi_jet, final_errors, 0, detail);    
    
//merging the delta phi with gap distribution
    TH1D *delta_phi_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_gap");
    TH1D *delta_phi_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_gap");
    TH1D *delta_phi_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_gap");
    
    TH1D *delta_phi_gap;
    delta_phi_gap =  new TH1D("ak5PF_delta_phi_gap","#Delta#phi^{gap};|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_gap->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi gap"<<endl; }
    calc_nomalization(delta_phi_gap_jetmettau, delta_phi_gap_jetmet, delta_phi_gap_jet, weigth, detail);
    merge(delta_phi_gap, delta_phi_gap_jetmettau, delta_phi_gap_jetmet, delta_phi_gap_jet, weigth);    
    plot_histogram(delta_phi_gap, delta_phi_gap_jetmettau, delta_phi_gap_jetmet, delta_phi_gap_jet, output_path_plots, "merging_delta_phi_gap", "top_left", true, detail);
    plot_histogram(delta_phi_gap, delta_phi_gap_jetmettau, delta_phi_gap_jetmet, delta_phi_gap_jet, output_path_plots, "merging_delta_phi_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_gap, delta_phi_gap_jetmettau, delta_phi_gap_jetmet, delta_phi_gap_jet, final_errors, 1, detail); 
    
//merging the delta phi without gap distribution 
    TH1D *delta_phi_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_nogap");
    TH1D *delta_phi_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_nogap");
    TH1D *delta_phi_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_nogap");
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap =  new TH1D("ak5PF_delta_phi_nogap","#Delta#phi^{no gap};|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Phi nogap"<<endl; }
    calc_nomalization(delta_phi_nogap_jetmettau, delta_phi_nogap_jetmet, delta_phi_nogap_jet, weigth, detail);
    merge(delta_phi_nogap, delta_phi_nogap_jetmettau, delta_phi_nogap_jetmet, delta_phi_nogap_jet, weigth);    
    plot_histogram(delta_phi_nogap, delta_phi_nogap_jetmettau, delta_phi_nogap_jetmet, delta_phi_nogap_jet, output_path_plots, "merging_delta_phi_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_nogap, delta_phi_nogap_jetmettau, delta_phi_nogap_jetmet, delta_phi_nogap_jet, output_path_plots, "merging_delta_phi_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_nogap, delta_phi_nogap_jetmettau, delta_phi_nogap_jetmet, delta_phi_nogap_jet, final_errors, 2, detail); 
    
//merging the leading pt inside the gap distribution 
    TH1D *leading_pt_inside_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_jet = (TH1D*) data_jet->Get("ak5PF_leading_pt_inside_gap");

    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap = new TH1D("ak5PF_leading_pt_inside_gap","Leading_pt_inside_gap;p_{T} [GeV];#frac{d#sigma}{dp_{T}}", in_nbins, in_bins);
    leading_pt_inside_gap->Sumw2();
    
    if (detail) { cout<<"Leading pt inside gap"<<endl; }
    calc_nomalization(leading_pt_inside_gap_jetmettau, leading_pt_inside_gap_jetmet, leading_pt_inside_gap_jet, weigth, detail);
    merge(leading_pt_inside_gap, leading_pt_inside_gap_jetmettau, leading_pt_inside_gap_jetmet, leading_pt_inside_gap_jet, weigth);
    plot_histogram(leading_pt_inside_gap, leading_pt_inside_gap_jetmettau, leading_pt_inside_gap_jetmet, leading_pt_inside_gap_jet, output_path_plots, "merging_leading_pt_inside_gap", "top_right", true, detail);
    plot_histogram(leading_pt_inside_gap, leading_pt_inside_gap_jetmettau, leading_pt_inside_gap_jetmet, leading_pt_inside_gap_jet, output_path_plots, "merging_leading_pt_inside_gap_linear", "top_right", false, detail);
    calc_error(leading_pt_inside_gap, leading_pt_inside_gap_jetmettau, leading_pt_inside_gap_jetmet, leading_pt_inside_gap_jet, final_errors, 3, detail);
    
//merging the leading pt outside the gap distribution    
    TH1D *leading_pt_outside_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_jet = (TH1D*) data_jet->Get("ak5PF_leading_pt_outside_gap");
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap = new TH1D("ak5PF_leading_pt_outside_gap","Leading_pt_outside_gap;p_{T} [GeV];#frac{d#sigma}{dp_{T}}", out_nbins, out_bins);
    leading_pt_outside_gap->Sumw2();
    
    if (detail) { cout<<"Leading pt outside gap "<<endl; }
    calc_nomalization(leading_pt_outside_gap_jetmettau, leading_pt_outside_gap_jetmet, leading_pt_outside_gap_jet, weigth, detail);
    merge(leading_pt_outside_gap, leading_pt_outside_gap_jetmettau, leading_pt_outside_gap_jetmet, leading_pt_outside_gap_jet, weigth);
    plot_histogram(leading_pt_outside_gap, leading_pt_outside_gap_jetmettau, leading_pt_outside_gap_jetmet, leading_pt_outside_gap_jet, output_path_plots, "merging_leading_pt_outside_gap", "top_right", true, detail);
    plot_histogram(leading_pt_outside_gap, leading_pt_outside_gap_jetmettau, leading_pt_outside_gap_jetmet, leading_pt_outside_gap_jet, output_path_plots, "merging_leading_pt_outside_gap_linear", "top_right", false, detail);
    calc_error(leading_pt_outside_gap, leading_pt_outside_gap_jetmettau, leading_pt_outside_gap_jetmet, leading_pt_outside_gap_jet, final_errors, 4, detail);
    
//merging the delta phi deta1 distributions
    TH1D *delta_phi_deta1_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta1");
    TH1D *delta_phi_deta1_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta1");
    TH1D *delta_phi_deta1_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta1");
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 = new TH1D("ak5PF_delta_phi_deta1","Delta_phi_deta1;|#Delta#phi| [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta1 "<<endl; }
    calc_nomalization(delta_phi_deta1_jetmettau, delta_phi_deta1_jetmet, delta_phi_deta1_jet, weigth, detail);
    merge(delta_phi_deta1, delta_phi_deta1_jetmettau, delta_phi_deta1_jetmet, delta_phi_deta1_jet, weigth);
    plot_histogram(delta_phi_deta1, delta_phi_deta1_jetmettau, delta_phi_deta1_jetmet, delta_phi_deta1_jet, output_path_plots, "merging_delta_phi_deta1", "top_left", true, detail);
    plot_histogram(delta_phi_deta1, delta_phi_deta1_jetmettau, delta_phi_deta1_jetmet, delta_phi_deta1_jet, output_path_plots, "merging_delta_phi_deta1_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1, delta_phi_deta1_jetmettau, delta_phi_deta1_jetmet, delta_phi_deta1_jet, final_errors, 5, detail);

//merging the delta phi deta2 distributions    
    TH1D *delta_phi_deta2_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta2");
    TH1D *delta_phi_deta2_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta2");
    TH1D *delta_phi_deta2_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta2");
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 = new TH1D("ak5PF_delta_phi_deta2","Delta_phi_deta2;|#Delta#phi| [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta2"<<endl; }
    calc_nomalization(delta_phi_deta2_jetmettau, delta_phi_deta2_jetmet, delta_phi_deta2_jet, weigth, detail);
    merge(delta_phi_deta2, delta_phi_deta2_jetmettau, delta_phi_deta2_jetmet, delta_phi_deta2_jet, weigth);
    plot_histogram(delta_phi_deta2, delta_phi_deta2_jetmettau, delta_phi_deta2_jetmet, delta_phi_deta2_jet, output_path_plots, "merging_delta_phi_deta2", "top_left", true, detail);
    plot_histogram(delta_phi_deta2, delta_phi_deta2_jetmettau, delta_phi_deta2_jetmet, delta_phi_deta2_jet, output_path_plots, "merging_delta_phi_deta2_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2, delta_phi_deta2_jetmettau, delta_phi_deta2_jetmet, delta_phi_deta2_jet, final_errors, 6, detail);
    

//merging the delta phi deta3 distributions    
    TH1D *delta_phi_deta3_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta3");
    TH1D *delta_phi_deta3_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta3");
    TH1D *delta_phi_deta3_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta3");

    TH1D *delta_phi_deta3;
    delta_phi_deta3 = new TH1D("ak5PF_delta_phi_deta3","Delta_phi_deta3;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta3 "<<endl; }
    calc_nomalization(delta_phi_deta3_jetmettau, delta_phi_deta3_jetmet, delta_phi_deta3_jet, weigth, detail);
    merge(delta_phi_deta3, delta_phi_deta3_jetmettau, delta_phi_deta3_jetmet, delta_phi_deta3_jet, weigth);
    plot_histogram(delta_phi_deta3, delta_phi_deta3_jetmettau, delta_phi_deta3_jetmet, delta_phi_deta3_jet, output_path_plots, "merging_delta_phi_deta3", "top_left", true, detail);
    plot_histogram(delta_phi_deta3, delta_phi_deta3_jetmettau, delta_phi_deta3_jetmet, delta_phi_deta3_jet, output_path_plots, "merging_delta_phi_deta3_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3, delta_phi_deta3_jetmettau, delta_phi_deta3_jetmet, delta_phi_deta3_jet, final_errors, 7, detail);

//merging the delta phi deta4 distributions         
    TH1D *delta_phi_deta4_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta4");
    TH1D *delta_phi_deta4_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta4");
    TH1D *delta_phi_deta4_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta4");

    TH1D *delta_phi_deta4;
    delta_phi_deta4 = new TH1D("ak5PF_delta_phi_deta4","Delta_phi_deta4;|#Delta#phi| [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta4 "<<endl; }
    calc_nomalization(delta_phi_deta4_jetmettau, delta_phi_deta4_jetmet, delta_phi_deta4_jet, weigth, detail);
    merge(delta_phi_deta4, delta_phi_deta4_jetmettau, delta_phi_deta4_jetmet, delta_phi_deta4_jet, weigth);
    plot_histogram(delta_phi_deta4, delta_phi_deta4_jetmettau, delta_phi_deta4_jetmet, delta_phi_deta4_jet, output_path_plots, "merging_delta_phi_deta4", "top_left", true, detail);
    plot_histogram(delta_phi_deta4, delta_phi_deta4_jetmettau, delta_phi_deta4_jetmet, delta_phi_deta4_jet, output_path_plots, "merging_delta_phi_deta4_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4, delta_phi_deta4_jetmettau, delta_phi_deta4_jetmet, delta_phi_deta4_jet, final_errors, 8, detail);

//merging the delta phi with gap deta1 distributions
    TH1D *delta_phi_deta1_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta1_gap");
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap = new TH1D("ak5PF_delta_phi_deta1_gap","Delta_phi_deta1_gap;|#Delta#phi| [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi deta1 gap"<<endl; }
    calc_nomalization(delta_phi_deta1_gap_jetmettau, delta_phi_deta1_gap_jetmet, delta_phi_deta1_gap_jet, weigth, detail);
    merge(delta_phi_deta1_gap, delta_phi_deta1_gap_jetmettau, delta_phi_deta1_gap_jetmet, delta_phi_deta1_gap_jet, weigth);
    plot_histogram(delta_phi_deta1_gap, delta_phi_deta1_gap_jetmettau, delta_phi_deta1_gap_jetmet, delta_phi_deta1_gap_jet, output_path_plots, "merging_delta_phi_deta1_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_gap, delta_phi_deta1_gap_jetmettau, delta_phi_deta1_gap_jetmet, delta_phi_deta1_gap_jet, output_path_plots, "merging_delta_phi_deta1_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_gap, delta_phi_deta1_gap_jetmettau, delta_phi_deta1_gap_jetmet, delta_phi_deta1_gap_jet, final_errors, 9, detail);

//merging the delta phi with gap deta2 distributions    
    TH1D *delta_phi_deta2_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta2_gap");

    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap = new TH1D("ak5PF_delta_phi_deta2_gap","lala;|#Delta#phi| [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta2"<<endl; }
    calc_nomalization(delta_phi_deta2_gap_jetmettau, delta_phi_deta2_gap_jetmet, delta_phi_deta2_gap_jet, weigth, detail);
    merge(delta_phi_deta2_gap, delta_phi_deta2_gap_jetmettau, delta_phi_deta2_gap_jetmet, delta_phi_deta2_gap_jet, weigth);
    plot_histogram(delta_phi_deta2_gap, delta_phi_deta2_gap_jetmettau, delta_phi_deta2_gap_jetmet, delta_phi_deta2_gap_jet, output_path_plots, "merging_delta_phi_deta2_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_gap, delta_phi_deta2_gap_jetmettau, delta_phi_deta2_gap_jetmet, delta_phi_deta2_gap_jet, output_path_plots, "merging_delta_phi_deta2_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_gap, delta_phi_deta2_gap_jetmettau, delta_phi_deta2_gap_jetmet, delta_phi_deta2_gap_jet, final_errors, 10, detail);

//merging the delta phi with gap deta3 distributions    
    TH1D *delta_phi_deta3_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta3_gap");
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap = new TH1D("ak5PF_delta_phi_deta3_gap","Delta_phi_deta3_gap;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta3"<<endl; }
    calc_nomalization(delta_phi_deta3_gap_jetmettau, delta_phi_deta3_gap_jetmet, delta_phi_deta3_gap_jet, weigth, detail);
    merge(delta_phi_deta3_gap, delta_phi_deta3_gap_jetmettau, delta_phi_deta3_gap_jetmet, delta_phi_deta3_gap_jet, weigth);
    plot_histogram(delta_phi_deta3_gap, delta_phi_deta3_gap_jetmettau, delta_phi_deta3_gap_jetmet, delta_phi_deta3_gap_jet, output_path_plots, "merging_delta_phi_deta3_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_gap, delta_phi_deta3_gap_jetmettau, delta_phi_deta3_gap_jetmet, delta_phi_deta3_gap_jet, output_path_plots, "merging_delta_phi_deta3_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_gap, delta_phi_deta3_gap_jetmettau, delta_phi_deta3_gap_jetmet, delta_phi_deta3_gap_jet, final_errors, 11, detail);
    
//merging the delta phi with gap deta4 distributions
    TH1D *delta_phi_deta4_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta4_gap");

    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap = new TH1D("ak5PF_delta_phi_deta4_gap","Delta_phi_deta4_gap;|#Delta#phi| [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi gap deta4"<<endl; }
    calc_nomalization(delta_phi_deta4_gap_jetmettau, delta_phi_deta4_gap_jetmet, delta_phi_deta4_gap_jet, weigth, detail);
    merge(delta_phi_deta4_gap, delta_phi_deta4_gap_jetmettau, delta_phi_deta4_gap_jetmet, delta_phi_deta4_gap_jet, weigth);
    plot_histogram(delta_phi_deta4_gap, delta_phi_deta4_gap_jetmettau, delta_phi_deta4_gap_jetmet, delta_phi_deta4_gap_jet, output_path_plots, "merging_delta_phi_deta4_gap", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_gap, delta_phi_deta4_gap_jetmettau, delta_phi_deta4_gap_jetmet, delta_phi_deta4_gap_jet, output_path_plots, "merging_delta_phi_deta4_gap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_gap, delta_phi_deta4_gap_jetmettau, delta_phi_deta4_gap_jetmet, delta_phi_deta4_gap_jet, final_errors, 12, detail);

//merging the delta phi without gap deta1 distributions    
    TH1D *delta_phi_deta1_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta1_nogap");

    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap = new TH1D("ak5PF_delta_phi_deta1_nogap","lala;|#Delta#phi| [rad] 0.4 < |#Delta#eta| < 2.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta1"<<endl; }
    calc_nomalization(delta_phi_deta1_nogap_jetmettau, delta_phi_deta1_nogap_jetmet, delta_phi_deta1_nogap_jet, weigth, detail);
    merge(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jetmettau, delta_phi_deta1_nogap_jetmet, delta_phi_deta1_nogap_jet, weigth);
    plot_histogram(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jetmettau, delta_phi_deta1_nogap_jetmet, delta_phi_deta1_nogap_jet, output_path_plots, "merging_delta_phi_deta1_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jetmettau, delta_phi_deta1_nogap_jetmet, delta_phi_deta1_nogap_jet, output_path_plots, "merging_delta_phi_deta1_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta1_nogap, delta_phi_deta1_nogap_jetmettau, delta_phi_deta1_nogap_jetmet, delta_phi_deta1_nogap_jet, final_errors, 13, detail);

//merging the delta phi without gap deta2 distributions
    TH1D *delta_phi_deta2_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta2_nogap");

    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap = new TH1D("ak5PF_delta_phi_deta2_nogap","lala;|#Delta#phi| [rad] 2.5 < |#Delta#eta| < 3.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta2"<<endl; }
    calc_nomalization(delta_phi_deta2_nogap_jetmettau, delta_phi_deta2_nogap_jetmet, delta_phi_deta2_nogap_jet, weigth, detail);
    merge(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jetmettau, delta_phi_deta2_nogap_jetmet, delta_phi_deta2_nogap_jet, weigth);
    plot_histogram(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jetmettau, delta_phi_deta2_nogap_jetmet, delta_phi_deta2_nogap_jet, output_path_plots, "merging_delta_phi_deta2_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jetmettau, delta_phi_deta2_nogap_jetmet, delta_phi_deta2_nogap_jet, output_path_plots, "merging_delta_phi_deta2_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta2_nogap, delta_phi_deta2_nogap_jetmettau, delta_phi_deta2_nogap_jetmet, delta_phi_deta2_nogap_jet, final_errors, 14, detail);

//merging the delta phi without gap deta3 distributions
    TH1D *delta_phi_deta3_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta3_nogap");

    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap = new TH1D("ak5PF_delta_phi_deta3_nogap","lala;|#Delta#phi| [rad] 3.5 < |#Delta#eta| < 4.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta3"<<endl; }
    calc_nomalization(delta_phi_deta3_nogap_jetmettau, delta_phi_deta3_nogap_jetmet, delta_phi_deta3_nogap_jet, weigth, detail);
    merge(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jetmettau, delta_phi_deta3_nogap_jetmet, delta_phi_deta3_nogap_jet, weigth);
    plot_histogram(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jetmettau, delta_phi_deta3_nogap_jetmet, delta_phi_deta3_nogap_jet, output_path_plots, "merging_delta_phi_deta3_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jetmettau, delta_phi_deta3_nogap_jetmet, delta_phi_deta3_nogap_jet, output_path_plots, "merging_delta_phi_deta3_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta3_nogap, delta_phi_deta3_nogap_jetmettau, delta_phi_deta3_nogap_jetmet, delta_phi_deta3_nogap_jet, final_errors, 15, detail);

//merging the delta phi without gap deta4 distributions
    TH1D *delta_phi_deta4_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_deta4_nogap");

    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap = new TH1D("ak5PF_delta_phi_deta4_nogap","lala;|#Delta#phi| [rad] 4.5 < |#Delta#eta| < 7.5;#frac{d#sigma}{d#Delta#phi#Delta#eta}", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta phi nogap deta4"<<endl; }
    calc_nomalization(delta_phi_deta4_nogap_jetmettau, delta_phi_deta4_nogap_jetmet, delta_phi_deta4_nogap_jet, weigth, detail);
    merge(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jetmettau, delta_phi_deta4_nogap_jetmet, delta_phi_deta4_nogap_jet, weigth);
    plot_histogram(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jetmettau, delta_phi_deta4_nogap_jetmet, delta_phi_deta4_nogap_jet, output_path_plots, "merging_delta_phi_deta4_nogap", "top_left", true, detail);
    plot_histogram(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jetmettau, delta_phi_deta4_nogap_jetmet, delta_phi_deta4_nogap_jet, output_path_plots, "merging_delta_phi_deta4_nogap_linear", "top_left", false, detail);
    calc_error(delta_phi_deta4_nogap, delta_phi_deta4_nogap_jetmettau, delta_phi_deta4_nogap_jetmet, delta_phi_deta4_nogap_jet, final_errors, 16, detail);

/*
//merging the leading forward pt distributions
    TH1D *leading_forward_pt_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_forward_pt");
    TH1D *leading_forward_pt_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_forward_pt");
    TH1D *leading_forward_pt_jet = (TH1D*) data_jet->Get("ak5PF_leading_forward_pt");

    TH1D *leading_forward_pt;
    leading_forward_pt = new TH1D("ak5PF_leading_forward_pt","Forward_pt;p_{T}^{forward} [GeV];#frac{d#sigma}{dp_{T}}", forw_nbins, forw_bins);
    leading_forward_pt->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading_forward pT"<<endl; }
    calc_nomalization(leading_forward_pt_jetmettau, leading_forward_pt_jetmet, leading_forward_pt_jet, weigth, detail);
    merge(leading_forward_pt, leading_forward_pt_jetmettau, leading_forward_pt_jetmet, leading_forward_pt_jet, weigth);    
    plot_histogram(leading_forward_pt, leading_forward_pt_jetmettau, leading_forward_pt_jetmet, leading_forward_pt_jet, output_path_plots, "merging_leading_forward_pt", "top_right", detail);
    calc_error(leading_forward_pt, leading_forward_pt_jetmettau, leading_forward_pt_jetmet, leading_forward_pt_jet, final_errors, 17, detail);

//merging the leading central pt distributions
    TH1D *leading_central_pt_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_central_pt");
    TH1D *leading_central_pt_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_central_pt");
    TH1D *leading_central_pt_jet = (TH1D*) data_jet->Get("ak5PF_leading_central_pt");

    TH1D *leading_central_pt;
    leading_central_pt = new TH1D("ak5PF_leading_central_pt","Central_pt;p_{T}^{central} [GeV];#frac{d#sigma}{dp_{T}}", cent_nbins, cent_bins);
    leading_central_pt->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading_central pT"<<endl; }
    calc_nomalization(leading_central_pt_jetmettau, leading_central_pt_jetmet, leading_central_pt_jet, weigth, detail);
    merge(leading_central_pt, leading_central_pt_jetmettau, leading_central_pt_jetmet, leading_central_pt_jet, weigth);    
    plot_histogram(leading_central_pt, leading_central_pt_jetmettau, leading_central_pt_jetmet, leading_central_pt_jet, output_path_plots, "merging_leading_central_pt", "top_right", detail);
    calc_error(leading_central_pt, leading_central_pt_jetmettau, leading_central_pt_jetmet, leading_central_pt_jet, final_errors, 18, detail);
*/


//merging the leading eta star inside gap   
    TH1D *leading_eta_star_inside_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_jet = (TH1D*) data_jet->Get("ak5PF_leading_eta_star_inside_gap");
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap = new TH1D("ak5PF_leading_eta_star_inside_gap","Leading_eta_star_inside_gap;#eta* ;#frac{d#sigma}{d#eta*}", etastar_nbins, etastar_bins);
    leading_eta_star_inside_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading Eta Star Gap"<<endl; }
    calc_nomalization(leading_eta_star_inside_gap_jetmettau, leading_eta_star_inside_gap_jetmet, leading_eta_star_inside_gap_jet, weigth, detail);
    merge(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jetmettau, leading_eta_star_inside_gap_jetmet, leading_eta_star_inside_gap_jet, weigth);
    plot_histogram(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jetmettau, leading_eta_star_inside_gap_jetmet, leading_eta_star_inside_gap_jet, output_path_plots, "merging_leading_eta_star_inside_gap", "bottom_middle", true, detail);
    plot_histogram(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jetmettau, leading_eta_star_inside_gap_jetmet, leading_eta_star_inside_gap_jet, output_path_plots, "merging_leading_eta_star_inside_gap_linear", "bottom_left", false, detail);
    calc_error(leading_eta_star_inside_gap, leading_eta_star_inside_gap_jetmettau, leading_eta_star_inside_gap_jetmet, leading_eta_star_inside_gap_jet, final_errors, 19, detail);

//merging the delta eta outside gap distribution
    TH1D *delta_eta_outside_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta_outside_gap");

    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap = new TH1D("ak5PF_delta_eta_outside_gap","delta_eta_outside_gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta}", deta_out_nbins, deta_out_bins);
    delta_eta_outside_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta eta outside gap"<<endl; }
    calc_nomalization(delta_eta_outside_gap_jetmettau, delta_eta_outside_gap_jetmet, delta_eta_outside_gap_jet, weigth, detail);
    merge(delta_eta_outside_gap, delta_eta_outside_gap_jetmettau, delta_eta_outside_gap_jetmet, delta_eta_outside_gap_jet, weigth);
    plot_histogram(delta_eta_outside_gap, delta_eta_outside_gap_jetmettau, delta_eta_outside_gap_jetmet, delta_eta_outside_gap_jet, output_path_plots, "merging_delta_eta_outside_gap", "bottom_left", true, detail);
    plot_histogram(delta_eta_outside_gap, delta_eta_outside_gap_jetmettau, delta_eta_outside_gap_jetmet, delta_eta_outside_gap_jet, output_path_plots, "merging_delta_eta_outside_gap_linear", "top_right", false, detail);
    calc_error(delta_eta_outside_gap, delta_eta_outside_gap_jetmettau, delta_eta_outside_gap_jetmet, delta_eta_outside_gap_jet, final_errors, 20, detail);


/*
//merging the leading forward eta distributions    
    TH1D *leading_forward_eta_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_forward_eta");
    TH1D *leading_forward_eta_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_forward_eta");
    TH1D *leading_forward_eta_jet = (TH1D*) data_jet->Get("ak5PF_leading_forward_eta");

    TH1D *leading_forward_eta;
    leading_forward_eta = new TH1D("ak5PF_leading_forward_eta","Forward_eta;#eta^{forward};#frac{d#sigma}{d#eta}", etaf_nbins, etaf_bins);
    leading_forward_eta->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading_forward Eta"<<endl; }
    calc_nomalization(leading_forward_eta_jetmettau, leading_forward_eta_jetmet, leading_forward_eta_jet, weigth, detail);
    merge(leading_forward_eta, leading_forward_eta_jetmettau, leading_forward_eta_jetmet, leading_forward_eta_jet, weigth);
    plot_histogram(leading_forward_eta, leading_forward_eta_jetmettau, leading_forward_eta_jetmet, leading_forward_eta_jet, output_path_plots, "merging_leading_forward_eta", "bottom_middle", detail);
    calc_error(leading_forward_eta, leading_forward_eta_jetmettau, leading_forward_eta_jetmet, leading_forward_eta_jet, final_errors, 21, detail);

//merging the leading central eta distributions
    TH1D *leading_central_eta_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_central_eta");
    TH1D *leading_central_eta_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_central_eta");
    TH1D *leading_central_eta_jet = (TH1D*) data_jet->Get("ak5PF_leading_central_eta");

    TH1D *leading_central_eta;
    leading_central_eta = new TH1D("ak5PF_leading_central_eta","Central_eta;#eta^{central};#frac{d#sigma}{d#eta}", etac_nbins, etac_bins);
    leading_central_eta->Sumw2();
    
    if (detail) { cout<<"ak5PF_leading_central Eta"<<endl; }
    calc_nomalization(leading_central_eta_jetmettau, leading_central_eta_jetmet, leading_central_eta_jet, weigth, detail);
    merge(leading_central_eta, leading_central_eta_jetmettau, leading_central_eta_jetmet, leading_central_eta_jet, weigth);
    plot_histogram(leading_central_eta, leading_central_eta_jetmettau, leading_central_eta_jetmet, leading_central_eta_jet, output_path_plots, "merging_leading_central_eta", "bottom_middle", detail);
    calc_error(leading_central_eta, leading_central_eta_jetmettau, leading_central_eta_jetmet, leading_central_eta_jet, final_errors, 22, detail);
    
//ploting the eta distributions
//declaring the canvas
    TCanvas *canvas = new TCanvas("canvas","Canvas",0,29,1450,870);
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

//draw the plot
    leading_forward_eta->SetTitle("Eta;#eta;#frac{d#sigma}{d#eta}");
    leading_forward_eta->Draw("e2");
    leading_forward_eta_jetmettau->Draw("e1same");
    leading_forward_eta_jetmet->Draw("e1same");
    leading_forward_eta_jet->Draw("e1same");
    leading_central_eta->Draw("e2same");
    leading_central_eta_jetmettau->Draw("e1same");
    leading_central_eta_jetmet->Draw("e1same");
    leading_central_eta_jet->Draw("e1same");
    
//draw the legend
    TLegend *legend = new TLegend(0.28,0.13,0.78,0.38);
    legend->AddEntry(leading_central_eta_jetmettau,"JetMETTau_2010A - L = 284 nb^{-1}","l");
    legend->AddEntry(leading_central_eta_jetmet,"JetMET_2010A - L = 2.90 pb^{-1}","l");
    legend->AddEntry(leading_central_eta_jet,"Jet_2010B - L = 32.1 pb^{-1}","l");
    legend->AddEntry(leading_central_eta,"Merged Data - L = 35.2 pb^{-1}","f");
    legend->SetFillColor(0);
    legend->SetLineWidth(1);
    legend->SetLineColor(0);
    legend->SetFillStyle(1001);
    legend->Draw();
    
//setting save paths
    string out_png = output_path_plots+"png/merging_leading_eta.png";
    string out_c = output_path_plots+"c/merging_leading_eta.C";
    string out_eps = output_path_plots+"eps/merging_leading_eta.eps";

//save the control plots
    canvas->Print( out_png.c_str() );
    canvas->Print( out_c.c_str() );
    canvas->Print( out_eps.c_str() );

//close the camvas
    canvas->Close();
    
//merging the delta eta distributions
    TH1D *delta_eta_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta");
    TH1D *delta_eta_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta");
    TH1D *delta_eta_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta");

    TH1D *delta_eta;
    delta_eta = new TH1D("ak5PF_delta_eta","Delta_eta;#Delta#eta;#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Eta"<<endl; }
    calc_nomalization(delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, weigth, detail);
    merge(delta_eta, delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, weigth);
    plot_histogram(delta_eta, delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, output_path_plots, "merging_delta_eta", "bottom_left");
    calc_error(delta_eta, delta_eta_jetmettau, delta_eta_jetmet, delta_eta_jet, final_errors, 23, detail);

//merging the delta eta with gap distributions
    TH1D *delta_eta_gap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta_gap");
    TH1D *delta_eta_gap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta_gap");
    TH1D *delta_eta_gap_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta_gap");

    TH1D *delta_eta_gap;
    delta_eta_gap = new TH1D("ak5PF_delta_eta_gap","Delta_eta_gap;#Delta#eta^{jet-veto};#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta_gap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Eta Gap"<<endl; }
    calc_nomalization(delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, weigth, detail);
    merge(delta_eta_gap, delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, weigth);
    plot_histogram(delta_eta_gap, delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, output_path_plots, "merging_delta_eta_gap", "bottom_left");
    calc_error(delta_eta_gap, delta_eta_gap_jetmettau, delta_eta_gap_jetmet, delta_eta_gap_jet, final_errors, 24, detail);

//merging the delta eta without gap distributions
    TH1D *delta_eta_nogap_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_eta_nogap");
    TH1D *delta_eta_nogap_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_eta_nogap");
    TH1D *delta_eta_nogap_jet = (TH1D*) data_jet->Get("ak5PF_delta_eta_nogap");

    TH1D *delta_eta_nogap;
    delta_eta_nogap = new TH1D("ak5PF_delta_eta_nogap","Delta_eta_nogap;#Delta#eta^{requiring a jet};#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta_nogap->Sumw2();
    
    if (detail) { cout<<"ak5PF_delta Eta nogap"<<endl; }
    calc_nomalization(delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, weigth, detail);
    merge(delta_eta_nogap, delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, weigth);
    plot_histogram(delta_eta_nogap, delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, output_path_plots, "merging_delta_eta_nogap", "bottom_left");
    calc_error(delta_eta_nogap, delta_eta_nogap_jetmettau, delta_eta_nogap_jetmet, delta_eta_nogap_jet, final_errors, 25, detail);

//merging the delta eta without gap distributions
    TH1D *vertex_selected_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_vertex_selected");
    TH1D *vertex_selected_jetmet = (TH1D*) data_jetmet->Get("ak5PF_vertex_selected");
    TH1D *vertex_selected_jet = (TH1D*) data_jet->Get("ak5PF_vertex_selected");

    TH1D *vertex_selected;
    vertex_selected = new TH1D("ak5PF_vertex_selected","Vertex multiplicity in selected events;Multiplicity;#frac{d#sigma}{dN} [pb]", 15, 0, 15);
    vertex_selected->Sumw2();
    
    if (detail) { cout<<"ak5PF_vertex_selected"<<endl; }
    calc_nomalization(vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, weigth, detail);
    merge(vertex_selected, vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, weigth);
    plot_histogram(vertex_selected, vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, output_path_plots, "merging_vertex_selected", "bottom_left", detail);
    calc_error(vertex_selected, vertex_selected_jetmettau, vertex_selected_jetmet, vertex_selected_jet, final_errors, 26, detail);

//merging the delta phi fine distribution
    TH1D *delta_phi_fine_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_delta_phi_fine");
    TH1D *delta_phi_fine_jetmet = (TH1D*) data_jetmet->Get("ak5PF_delta_phi_fine");
    TH1D *delta_phi_fine_jet = (TH1D*) data_jet->Get("ak5PF_delta_phi_fine");
    
    TH1D *delta_phi_fine;
    delta_phi_fine = new TH1D("ak5PF_delta_phi_fine","#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", 14, 0, 3.15);
    delta_phi_fine->Sumw2();

    if (detail) { cout<<"ak5PF_delta Phi"<<endl; }
    calc_nomalization(delta_phi_fine_jetmettau, delta_phi_fine_jetmet, delta_phi_fine_jet, weigth, detail);
    merge(delta_phi_fine, delta_phi_fine_jetmettau, delta_phi_fine_jetmet, delta_phi_fine_jet, weigth);    
    plot_histogram(delta_phi_fine, delta_phi_fine_jetmettau, delta_phi_fine_jetmet, delta_phi_fine_jet, output_path_plots, "merging_delta_phi_fine", "top_left", detail);
    calc_error(delta_phi_fine, delta_phi_fine_jetmettau, delta_phi_fine_jetmet, delta_phi_fine_jet, final_errors, 27, detail);

//merging the leading pt distribution
    TH1D *leading_pt_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_pt");
    TH1D *leading_pt_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_pt");
    TH1D *leading_pt_jet = (TH1D*) data_jet->Get("ak5PF_leading_pt");
    
    TH1D *leading_pt;
    leading_pt = new TH1D("ak5PF_leading_pt","Leading p_{T}^{leading};p_{T}^{leading} [GeV];#frac{d#sigma}{dp_{T}}", all_nbins, all_bins);
    leading_pt->Sumw2();

    if (detail) { cout<<"ak5PF_leading_pt"<<endl; }
    calc_nomalization(leading_pt_jetmettau, leading_pt_jetmet, leading_pt_jet, weigth, detail);
    merge(leading_pt, leading_pt_jetmettau, leading_pt_jetmet, leading_pt_jet, weigth);    
    plot_histogram(leading_pt, leading_pt_jetmettau, leading_pt_jetmet, leading_pt_jet, output_path_plots, "merging_leading_pt", "top_right", detail);
    calc_error(leading_pt, leading_pt_jetmettau, leading_pt_jetmet, leading_pt_jet, final_errors, 28, detail); 

//merging the leading pt distribution
    TH1D *leading_pt_fine_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_pt_fine");
    TH1D *leading_pt_fine_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_pt_fine");
    TH1D *leading_pt_fine_jet = (TH1D*) data_jet->Get("ak5PF_leading_pt_fine");
    
    TH1D *leading_pt_fine;
    leading_pt_fine = new TH1D("ak5PF_leading_pt_fine","Leading p_{T}^{leading} with fine binning;p_{T}^{leading} [GeV];#frac{d#sigma}{dp_{T}}", 85, 30, 200);
    leading_pt_fine->Sumw2();

    if (detail) { cout<<"ak5PF_leading_pt_fine"<<endl; }
    calc_nomalization(leading_pt_fine_jetmettau, leading_pt_fine_jetmet, leading_pt_fine_jet, weigth, detail);
    merge(leading_pt_fine, leading_pt_fine_jetmettau, leading_pt_fine_jetmet, leading_pt_fine_jet, weigth);    
    plot_histogram(leading_pt_fine, leading_pt_fine_jetmettau, leading_pt_fine_jetmet, leading_pt_fine_jet, output_path_plots, "merging_leading_pt_fine", "top_right", detail);
    calc_error(leading_pt_fine, leading_pt_fine_jetmettau, leading_pt_fine_jetmet, leading_pt_fine_jet, final_errors, 29, detail); 

//merging the leading central pt distribution
    TH1D *leading_central_pt_fine_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_central_pt_fine");
    TH1D *leading_central_pt_fine_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_central_pt_fine");
    TH1D *leading_central_pt_fine_jet = (TH1D*) data_jet->Get("ak5PF_leading_central_pt_fine");
    
    TH1D *leading_central_pt_fine;
    leading_central_pt_fine = new TH1D("ak5PF_leading_central_pt_fine","Leading p_{T}^{central-leading} with fine binning;p_{T}^{central-leading} [GeV];#frac{d#sigma}{dp_{T}}", 85, 30, 200);
    leading_central_pt_fine->Sumw2();

    if (detail) { cout<<"ak5PF_leading_central_pt_fine"<<endl; }
    calc_nomalization(leading_central_pt_fine_jetmettau, leading_central_pt_fine_jetmet, leading_central_pt_fine_jet, weigth, detail);
    merge(leading_central_pt_fine, leading_central_pt_fine_jetmettau, leading_central_pt_fine_jetmet, leading_central_pt_fine_jet, weigth);    
    plot_histogram(leading_central_pt_fine, leading_central_pt_fine_jetmettau, leading_central_pt_fine_jetmet, leading_central_pt_fine_jet, output_path_plots, "merging_leading_central_pt_fine", "top_right", detail);
    calc_error(leading_central_pt_fine, leading_central_pt_fine_jetmettau, leading_central_pt_fine_jetmet, leading_central_pt_fine_jet, final_errors, 30, detail); 

//merging the leading forward pt distribution
    TH1D *leading_forward_pt_fine_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_leading_forward_pt_fine");
    TH1D *leading_forward_pt_fine_jetmet = (TH1D*) data_jetmet->Get("ak5PF_leading_forward_pt_fine");
    TH1D *leading_forward_pt_fine_jet = (TH1D*) data_jet->Get("ak5PF_leading_forward_pt_fine");
    
    TH1D *leading_forward_pt_fine;
    leading_forward_pt_fine = new TH1D("ak5PF_leading_forward_pt_fine","Leading p_{T}^{forward-leading} with fine binning;p_{T}^{forward-leading} [GeV];#frac{d#sigma}{dp_{T}}", 85, 30, 200);
    leading_forward_pt_fine->Sumw2();

    if (detail) { cout<<"ak5PF_leading_forward_pt_fine"<<endl; }
    calc_nomalization(leading_forward_pt_fine_jetmettau, leading_forward_pt_fine_jetmet, leading_forward_pt_fine_jet, weigth, detail);
    merge(leading_forward_pt_fine, leading_forward_pt_fine_jetmettau, leading_forward_pt_fine_jetmet, leading_forward_pt_fine_jet, weigth);    
    plot_histogram(leading_forward_pt_fine, leading_forward_pt_fine_jetmettau, leading_forward_pt_fine_jetmet, leading_forward_pt_fine_jet, output_path_plots, "merging_leading_forward_pt_fine", "top_right", detail);
    calc_error(leading_forward_pt_fine, leading_forward_pt_fine_jetmettau, leading_forward_pt_fine_jetmet, leading_forward_pt_fine_jet, final_errors, 31, detail); 
    
    //merging the pvz selected distribution
    TH1D *pvz_selected_jetmettau = (TH1D*) data_jetmettau->Get("ak5PF_pvz_selected");
    TH1D *pvz_selected_jetmet = (TH1D*) data_jetmet->Get("ak5PF_pvz_selected");
    TH1D *pvz_selected_jet = (TH1D*) data_jet->Get("ak5PF_pvz_selected");
    
    TH1D *pvz_selected;
    pvz_selected = new TH1D("ak5PF_pvz_selected","z of the Primary Vertex in selected events;z-position;#frac{d#sigma}{dz} [pb]", 200, -100, 100);
    pvz_selected->Sumw2();

    if (detail) { cout<<"ak5PF_pvz_selected"<<endl; }
    calc_nomalization(pvz_selected_jetmettau, pvz_selected_jetmet, pvz_selected_jet, weigth, detail);
    merge(pvz_selected, pvz_selected_jetmettau, pvz_selected_jetmet, pvz_selected_jet, weigth);    
    plot_histogram(pvz_selected, pvz_selected_jetmettau, pvz_selected_jetmet, pvz_selected_jet, output_path_plots, "merging_pvz_selected", "bottom_right", detail);
    calc_error(pvz_selected, pvz_selected_jetmettau, pvz_selected_jetmet, pvz_selected_jet, final_errors, 32, detail);
*/


//output the error variation
    if (detail) { cout<<"Display the final errors..."<<endl; }
    if (disp_errors) { show_error_variation(final_errors); }

//creating the output file
    TFile data_output( output_rootfile.c_str() , "RECREATE");

//save the histograms in a root file
    if (detail) { cout<<"Writing histograms on file "<<output_rootfile<<" ..."<<endl; }
    delta_phi->Write();
    //delta_phi_fine->Write();
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
    leading_pt_outside_gap->Write();
    //leading_central_pt->Write();
    //leading_forward_pt->Write();
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();
    //leading_central_eta->Write();
    //leading_forward_eta->Write();
    //delta_eta->Write();
    //delta_eta_gap->Write();
    //delta_eta_nogap->Write();
    //vertex_selected->Write();
    //leading_pt->Write();
    //leading_pt_fine->Write();
    //leading_central_pt_fine->Write();
    //leading_forward_pt_fine->Write();
    //pvz_selected->Write();
    if (detail) { cout<<"Writing on "<<output_rootfile<<" was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    data_jetmettau->Close();
    data_jetmet->Close();
    data_jet->Close();
    data_output.Close();
    if (detail) { cout<<"Done!"<<endl; }
       
}
