// Pedro Cipriano, Mar 2012
// DESY, CMS
// Last Update: 20 Nov 2012
//
// normalize_mc(string input_file_allvertex, string input_file_1vertex, string output_path_plots, string output_file_rootfile, string plot_prefix, string hist_prefix, bool detail, bool disp_errors)
// Normalize the MC to get the correct cross-section 

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

double normalize_histogram(TH1D *histogram, TH1D *hist_allvertex, TH1D *hist_1vertex, string path, string plot_prefix, string hist_name, string legend_position = "top_left", bool detail = false)
{
//normalize hist_1vertex using the integral of hist_allvertex and saves the result in histogram

//declaring variables
   double integral_allvertex, integral_1vertex, scale;

   if (detail) { cout<<"Processing "<<hist_name<<" ..."<<endl; }

//clone the hist_1vertex
    for(Int_t i=1;i<=hist_1vertex->GetNbinsX();i++)
    {
        histogram->SetBinContent(i,hist_1vertex->GetBinContent(i));
        histogram->SetBinError(i,hist_1vertex->GetBinError(i));
        //if (detail) { cout << i << " -> " << hist_1vertex->GetBinContent(i) << " +- " << hist_1vertex->GetBinError(i) << endl; }
    }
    histogram->SetEntries(hist_1vertex->GetEntries());

//compute the scaling factor
    integral_1vertex = hist_1vertex->Integral();
    integral_allvertex = hist_allvertex->Integral();
    scale = integral_allvertex/integral_1vertex;
    if (detail) { cout<<"Scale = "<<scale<<endl; }

//apply the scale factor
    histogram->Scale(scale);

//declare and configure the canvas
    string fileout = plot_prefix + hist_name;
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

//calculate the plooting range
    double min = 0.0;
    double max = histogram->GetMaximum();
    if (histogram->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(histogram,detail);
    }
    else
    {
    min = histogram->GetMinimum();
    }    
    
    set_histogram_min_max(hist_1vertex, min, max, detail);
    set_histogram_min_max(hist_allvertex, min, max, detail);

    max = 1.3 * max;
    min = 0.7 * min;
    
//plooting
    histogram->SetMaximum(max);
    histogram->SetMinimum(min);
    histogram->SetLineWidth(4);
    histogram->SetLineColor(2);
    histogram->SetLineStyle(1);
    histogram->Draw("e1");
    hist_allvertex->SetLineWidth(4);
    hist_allvertex->SetLineColor(3);
    hist_allvertex->SetLineStyle(2);
    hist_allvertex->Draw("e1same");
    hist_1vertex->SetLineWidth(4);
    hist_1vertex->SetLineColor(4);
    hist_1vertex->SetLineStyle(4);
    hist_1vertex->Draw("e1same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(histogram,"Final cross-section","l");
    leg01->AddEntry(hist_allvertex,"Normalization cross-section","l");
    leg01->AddEntry(hist_1vertex,"Shape cross-section","l");
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
    c01->Print( out_png.c_str() );
    c01->Print( out_c.c_str() );
    c01->Print( out_eps.c_str() );
    c01->Close();

//return the normalization scale
    return scale;
}

void show_scaling_factors(double *scales)
{
//shows the scaling factors applied to the histograms

cout << "Normalization scales" << endl;
cout << "Observable                         Scaling factor" << endl;
cout << "Delta Phi                          " << scales[0] << endl;
cout << "Delta Phi Deta1                    " << scales[1] << endl;
cout << "Delta Phi Deta2                    " << scales[2] << endl;
cout << "Delta Phi Deta3                    " << scales[3] << endl;
cout << "Delta Phi Deta4                    " << scales[4] << endl;
cout << "Delta Phi Gap                      " << scales[5] << endl;
cout << "Delta Phi Deta1 Gap                " << scales[6] << endl;
cout << "Delta Phi Deta2 Gap                " << scales[7] << endl;
cout << "Delta Phi Deta3 Gap                " << scales[8] << endl;
cout << "Delta Phi Deta4 Gap                " << scales[9] << endl;
cout << "Delta Phi noGap                    " << scales[10] << endl;
cout << "Delta Phi Deta1 noGap              " << scales[11] << endl;
cout << "Delta Phi Deta2 noGap              " << scales[12] << endl;
cout << "Delta Phi Deta3 noGap              " << scales[13] << endl;
cout << "Delta Phi Deta4 noGap              " << scales[14] << endl;
cout << "Delta Eta                          " << scales[15] << endl;
cout << "Delta Eta Gap                      " << scales[16] << endl;
cout << "Delta Eta noGap                    " << scales[17] << endl;
cout << "Total pT Inside Gap                " << scales[18] << endl;
cout << "Leading pT Inside Gap              " << scales[19] << endl;
cout << "Leading Eta Inside Gap             " << scales[20] << endl;
cout << "Leading Phi Inside Gap             " << scales[21] << endl;
cout << "Total pT Outside Gap               " << scales[22] << endl;
cout << "Leading pT Outside Gap             " << scales[23] << endl;
cout << "Leading Eta Outside Gap            " << scales[24] << endl;
cout << "Leading Phi Outside Gap            " << scales[25] << endl;
cout << "Leading Eta* Inside Gap            " << scales[26] << endl;
cout << "Delta Eta Outside Gap              " << scales[27] << endl;
cout << "Leading Central pT                 " << scales[28] << endl;
cout << "Leading Central pT Gap             " << scales[29] << endl;
cout << "Leading Central pT noGap           " << scales[30] << endl;
cout << "Leading Central Eta                " << scales[31] << endl;
cout << "Leading Central Phi                " << scales[32] << endl;
cout << "Leading Forward pT                 " << scales[33] << endl;
cout << "Leading Forward pT Gap             " << scales[34] << endl;
cout << "Leading Forward pT noGap           " << scales[35] << endl;
cout << "Leading Forward Eta                " << scales[36] << endl;
cout << "Leading forward Phi                " << scales[37] << endl;
cout << "Leading Eta                        " << scales[38] << endl;
cout << "Vertex Selected                    " << scales[39] << endl;
cout << "Jet Multiplicity                   " << scales[40] << endl;
cout << "pT all                             " << scales[41] << endl;
cout << "Eta all                            " << scales[42] << endl;
cout << "Phi all                            " << scales[43] << endl;
cout << "Inclusive Leading Forward pT       " << scales[44] << endl;
cout << "Inclusive Leading Central pT       " << scales[45] << endl;
cout << "Inclusive Forward pT               " << scales[46] << endl;
cout << "Inclusive Central pT               " << scales[47] << endl;
cout << "Delta Phi Central Relative         " << scales[48] << endl;
cout << "Delta Phi Forward Relative         " << scales[49] << endl;
cout << "Delta Phi Central Relative Small   " << scales[50] << endl;
cout << "Delta Phi Forward Relative Small   " << scales[51] << endl;
cout << "Delta Phi Central Relative Medium  " << scales[52] << endl;
cout << "Delta Phi Forward Relative Medium  " << scales[53] << endl;
cout << "Delta Phi Central Relative Large   " << scales[54] << endl;
cout << "Delta Phi Forward Relative Large   " << scales[55] << endl;
cout << "Delta pT Central Relative          " << scales[56] << endl;
cout << "Delta pT Forward Relative          " << scales[57] << endl;
cout << "Delta Phi Fine                     " << scales[58] << endl;
//cout << "Leading pT                         " << scales[59] << endl;
//cout << "Leading pT Fine                    " << scales[60] << endl;
cout << "Leading Central pT Fine            " << scales[61] << endl;
cout << "Leading Forward pT Fine            " << scales[62] << endl;
cout << "Primary Vertex Z-position          " << scales[63] << endl;
}

void normalize_mc(string input_file_allvertex = "../output/xsec_mc_gen/xsec_Pythia6_TuneZ2star_gen_allvertex.root", string input_file_1vertex = "../output/xsec_mc_gen/xsec_Pythia6_TuneZ2star_gen_1vertex.root", string output_path_plots = "../output/normalize_mc/", string output_file_rootfile = "../output/histograms/xsec_p6_z2_gen.root", string plot_prefix = "pythia6_z2_", string hist_prefix = "ak5Gen_", bool detail = false, bool show_scales = true)
{

   if (detail) { cout<<"Input File with normalization : "<<input_file_allvertex<<endl; }
   if (detail) { cout<<"Input File with shape         : "<<input_file_1vertex<<endl; }
   if (detail) { cout<<"Plot Output Path              : "<<output_path_plots<<endl; }
   if (detail) { cout<<"Root Output File              : "<<output_file_rootfile<<endl; }
   if (detail) { cout<<"Histograms Prefix             : "<<hist_prefix<<endl; }
   if (detail) { cout<<"Plots Prefix                  : "<<plot_prefix<<endl; }
   if (detail) { cout<<"Show Details                  : "<<detail<<endl; }
   if (detail) { cout<<"Show Scaling Factors          : "<<show_scales<<endl; }

//opening the input data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *mc_allvertex = 0;
    mc_allvertex = new TFile( input_file_allvertex.c_str() );
    if (mc_allvertex == 0) { cout << "File " << input_file_allvertex << " not found!" << endl; return; }
    TFile *mc_1vertex = 0;
    mc_1vertex = new TFile( input_file_1vertex.c_str() );
    if (mc_1vertex == 0) { cout << "File " << input_file_1vertex << " not found!" << endl; return; }
    if (detail) { cout<<"All files opened sucessfully!"<<endl; }

//declare histogram bins
int cent_nbins = 7;
int forw_nbins = 7;

double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};

int all_nbins = 11;
int in_nbins = 9;
int out_nbins = 9;
int dpt_nbins = 10;

double all_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double dpt_bins[11] = {0, 2.5, 5, 10, 15, 20, 25, 30, 40, 60, 100};

int deta_nbins = 4;
int dphi_nbins = 7;

double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etac_nbins = 6;
double etac_bins[7] = {-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8};

int etaf_nbins = 7;
double etaf_bins[8] = {-4.7,-4.2,-3.7,-3.2,3.2,3.7,4.2,4.7};

int eta_nbins = 14;
double eta_bins[15] = {-4.7,-4.2,-3.7,-3.2,-2.8,-2.0,-1.0,0.0,1.0,2.0,2.8,3.2,3.7,4.2,4.7};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//declaring variables
   TString hist_name_in;
   string  hist_name_out;
   double scales[64];

//normalizing the delta phi distribution
    hist_name_in = hist_prefix+"delta_phi";
    hist_name_out = hist_prefix+"delta_phi";
    TH1D *delta_phi_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi;
    delta_phi = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi->Sumw2();

    scales[0] = normalize_histogram(delta_phi, delta_phi_allvertex, delta_phi_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta1 distribution
    hist_name_in = hist_prefix+"delta_phi_deta1";
    hist_name_out = hist_prefix+"delta_phi_deta1";
    TH1D *delta_phi_deta1_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta1_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta1;
    delta_phi_deta1 = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta1->Sumw2();

    scales[1] = normalize_histogram(delta_phi_deta1, delta_phi_deta1_allvertex, delta_phi_deta1_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta2 distribution
    hist_name_in = hist_prefix+"delta_phi_deta2";
    hist_name_out = hist_prefix+"delta_phi_deta2";
    TH1D *delta_phi_deta2_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta2_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta2;
    delta_phi_deta2 = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta2->Sumw2();

    scales[2] = normalize_histogram(delta_phi_deta2, delta_phi_deta2_allvertex, delta_phi_deta2_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta3 distribution
    hist_name_in = hist_prefix+"delta_phi_deta3";
    hist_name_out = hist_prefix+"delta_phi_deta3";
    TH1D *delta_phi_deta3_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta3_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta3;
    delta_phi_deta3 = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta3->Sumw2();

    scales[3] = normalize_histogram(delta_phi_deta3, delta_phi_deta3_allvertex, delta_phi_deta3_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta4 distribution
    hist_name_in = hist_prefix+"delta_phi_deta4";
    hist_name_out = hist_prefix+"delta_phi_deta4";
    TH1D *delta_phi_deta4_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta4_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta4;
    delta_phi_deta4 = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta4->Sumw2();

    scales[4] = normalize_histogram(delta_phi_deta4, delta_phi_deta4_allvertex, delta_phi_deta4_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi gap distribution
    hist_name_in = hist_prefix+"delta_phi_gap";
    hist_name_out = hist_prefix+"delta_phi_gap";
    TH1D *delta_phi_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_gap;
    delta_phi_gap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_gap->Sumw2();

    scales[5] = normalize_histogram(delta_phi_gap, delta_phi_gap_allvertex, delta_phi_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);


//normalizing the delta phi deta1 gap distribution
    hist_name_in = hist_prefix+"delta_phi_deta1_gap";
    hist_name_out = hist_prefix+"delta_phi_deta1_gap";
    TH1D *delta_phi_deta1_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta1_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta1_gap;
    delta_phi_deta1_gap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap->Sumw2();

    scales[6] = normalize_histogram(delta_phi_deta1_gap, delta_phi_deta1_gap_allvertex, delta_phi_deta1_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta2 gap distribution
    hist_name_in = hist_prefix+"delta_phi_deta2_gap";
    hist_name_out = hist_prefix+"delta_phi_deta2_gap";
    TH1D *delta_phi_deta2_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta2_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta2_gap;
    delta_phi_deta2_gap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap->Sumw2();

    scales[7] = normalize_histogram(delta_phi_deta2_gap, delta_phi_deta2_gap_allvertex, delta_phi_deta2_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta3 gap distribution
    hist_name_in = hist_prefix+"delta_phi_deta3_gap";
    hist_name_out = hist_prefix+"delta_phi_deta3_gap";
    TH1D *delta_phi_deta3_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta3_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta3_gap;
    delta_phi_deta3_gap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap->Sumw2();

    scales[8] = normalize_histogram(delta_phi_deta3_gap, delta_phi_deta3_gap_allvertex, delta_phi_deta3_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta4 gap distribution
    hist_name_in = hist_prefix+"delta_phi_deta4_gap";
    hist_name_out = hist_prefix+"delta_phi_deta4_gap";
    TH1D *delta_phi_deta4_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta4_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta4_gap;
    delta_phi_deta4_gap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap->Sumw2();

    scales[9] = normalize_histogram(delta_phi_deta4_gap, delta_phi_deta4_gap_allvertex, delta_phi_deta4_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi nogap distribution
    hist_name_in = hist_prefix+"delta_phi_nogap";
    hist_name_out = hist_prefix+"delta_phi_nogap";
    TH1D *delta_phi_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_nogap;
    delta_phi_nogap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_nogap->Sumw2();

    scales[10] = normalize_histogram(delta_phi_nogap, delta_phi_nogap_allvertex, delta_phi_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta1 nogap distribution
    hist_name_in = hist_prefix+"delta_phi_deta1_nogap";
    hist_name_out = hist_prefix+"delta_phi_deta1_nogap";
    TH1D *delta_phi_deta1_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta1_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta1_nogap;
    delta_phi_deta1_nogap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap->Sumw2();

    scales[11] = normalize_histogram(delta_phi_deta1_nogap, delta_phi_deta1_nogap_allvertex, delta_phi_deta1_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta2 nogap distribution
    hist_name_in = hist_prefix+"delta_phi_deta2_nogap";
    hist_name_out = hist_prefix+"delta_phi_deta2_nogap";
    TH1D *delta_phi_deta2_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta2_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta2_nogap;
    delta_phi_deta2_nogap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap->Sumw2();

    scales[12] = normalize_histogram(delta_phi_deta2_nogap, delta_phi_deta2_nogap_allvertex, delta_phi_deta2_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta3 nogap distribution
    hist_name_in = hist_prefix+"delta_phi_deta3_nogap";
    hist_name_out = hist_prefix+"delta_phi_deta3_nogap";
    TH1D *delta_phi_deta3_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta3_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta3_nogap;
    delta_phi_deta3_nogap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap->Sumw2();

    scales[13] = normalize_histogram(delta_phi_deta3_nogap, delta_phi_deta3_nogap_allvertex, delta_phi_deta3_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi deta4 nogap distribution
    hist_name_in = hist_prefix+"delta_phi_deta4_nogap";
    hist_name_out = hist_prefix+"delta_phi_deta4_nogap";
    TH1D *delta_phi_deta4_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_deta4_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_deta4_nogap;
    delta_phi_deta4_nogap = new TH1D(hist_name_in,"#Delta#phi;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap->Sumw2();

    scales[14] = normalize_histogram(delta_phi_deta4_nogap, delta_phi_deta4_nogap_allvertex, delta_phi_deta4_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta eta distribution
    hist_name_in = hist_prefix+"delta_eta";
    hist_name_out = hist_prefix+"delta_eta";
    TH1D *delta_eta_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_eta_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_eta;
    delta_eta = new TH1D(hist_name_in,"#Delta#eta;|#Delta#eta| [rad];#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta->Sumw2();

    scales[15] = normalize_histogram(delta_eta, delta_eta_allvertex, delta_eta_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the delta eta gap distribution
    hist_name_in = hist_prefix+"delta_eta_gap";
    hist_name_out = hist_prefix+"delta_eta_gap";
    TH1D *delta_eta_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_eta_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_eta_gap;
    delta_eta_gap = new TH1D(hist_name_in,"#Delta#eta;|#Delta#eta| [rad];#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta_gap->Sumw2();

    scales[16] = normalize_histogram(delta_eta_gap, delta_eta_gap_allvertex, delta_eta_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the delta eta nogap distribution
    hist_name_in = hist_prefix+"delta_eta_nogap";
    hist_name_out = hist_prefix+"delta_eta_nogap";
    TH1D *delta_eta_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_eta_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_eta_nogap;
    delta_eta_nogap = new TH1D(hist_name_in,"#Delta#eta;|#Delta#eta| [rad];#frac{d#sigma}{d#Delta#eta}", deta_nbins, deta_bins);
    delta_eta_nogap->Sumw2();

    scales[17] = normalize_histogram(delta_eta_nogap, delta_eta_nogap_allvertex, delta_eta_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the total pt inside gap distribution
    hist_name_in = hist_prefix+"total_pt_inside_gap";
    hist_name_out = hist_prefix+"total_pt_inside_gap";
    TH1D *total_pt_inside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *total_pt_inside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *total_pt_inside_gap;
    total_pt_inside_gap = new TH1D(hist_name_in,"Total Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
    total_pt_inside_gap->Sumw2();

    scales[18] = normalize_histogram(total_pt_inside_gap, total_pt_inside_gap_allvertex, total_pt_inside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading pt inside gap distribution
    hist_name_in = hist_prefix+"leading_pt_inside_gap";
    hist_name_out = hist_prefix+"leading_pt_inside_gap";
    TH1D *leading_pt_inside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_pt_inside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_pt_inside_gap;
    leading_pt_inside_gap = new TH1D(hist_name_in,"Leading Jet p_{T} inside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", in_nbins, in_bins);
    leading_pt_inside_gap->Sumw2();

    scales[19] = normalize_histogram(leading_pt_inside_gap, leading_pt_inside_gap_allvertex, leading_pt_inside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading eta inside gap distribution
    hist_name_in = hist_prefix+"leading_eta_inside_gap";
    hist_name_out = hist_prefix+"leading_eta_inside_gap";
    TH1D *leading_eta_inside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_eta_inside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_eta_inside_gap;
    leading_eta_inside_gap = new TH1D(hist_name_in,"Leading Jet #eta inside the gap;#eta;#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
    leading_eta_inside_gap->Sumw2();

    scales[20] = normalize_histogram(leading_eta_inside_gap, leading_eta_inside_gap_allvertex, leading_eta_inside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading phi inside gap distribution
    hist_name_in = hist_prefix+"leading_phi_inside_gap";
    hist_name_out = hist_prefix+"leading_phi_inside_gap";
    TH1D *leading_phi_inside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_phi_inside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_phi_inside_gap;
    leading_phi_inside_gap = new TH1D(hist_name_in,"Leading Jet #phi inside the gap;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
    leading_phi_inside_gap->Sumw2();

    scales[21] = normalize_histogram(leading_phi_inside_gap, leading_phi_inside_gap_allvertex, leading_phi_inside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the total pt outside gap distribution
    hist_name_in = hist_prefix+"total_pt_outside_gap";
    hist_name_out = hist_prefix+"total_pt_outside_gap";
    TH1D *total_pt_outside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *total_pt_outside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *total_pt_outside_gap;
    total_pt_outside_gap = new TH1D(hist_name_in,"Total Jet p_{T} outside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", out_nbins, out_bins);
    total_pt_outside_gap->Sumw2();

    scales[22] = normalize_histogram(total_pt_outside_gap, total_pt_outside_gap_allvertex, total_pt_outside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading pt outside gap distribution
    hist_name_in = hist_prefix+"leading_pt_outside_gap";
    hist_name_out = hist_prefix+"leading_pt_outside_gap";
    TH1D *leading_pt_outside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_pt_outside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_pt_outside_gap;
    leading_pt_outside_gap = new TH1D(hist_name_in,"Leading Jet p_{T} outside the gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", out_nbins, out_bins);
    leading_pt_outside_gap->Sumw2();

    scales[23] = normalize_histogram(leading_pt_outside_gap, leading_pt_outside_gap_allvertex, leading_pt_outside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading eta outside gap distribution
    hist_name_in = hist_prefix+"leading_eta_outside_gap";
    hist_name_out = hist_prefix+"leading_eta_outside_gap";
    TH1D *leading_eta_outside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_eta_outside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_eta_outside_gap;
    leading_eta_outside_gap = new TH1D(hist_name_in,"Leading Jet #eta outside the gap;#eta;#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
    leading_eta_outside_gap->Sumw2();

    scales[24] = normalize_histogram(leading_eta_outside_gap, leading_eta_outside_gap_allvertex, leading_eta_outside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading phi outside gap distribution
    hist_name_in = hist_prefix+"leading_phi_outside_gap";
    hist_name_out = hist_prefix+"leading_phi_outside_gap";
    TH1D *leading_phi_outside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_phi_outside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_phi_outside_gap;
    leading_phi_outside_gap = new TH1D(hist_name_in,"Leading Jet #phi outside the gap;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
    leading_phi_outside_gap->Sumw2();

    scales[25] = normalize_histogram(leading_phi_outside_gap, leading_phi_outside_gap_allvertex, leading_phi_outside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading eta star inside gap distribution
    hist_name_in = hist_prefix+"leading_eta_star_inside_gap";
    hist_name_out = hist_prefix+"leading_eta_star_inside_gap";
    TH1D *leading_eta_star_inside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_eta_star_inside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_eta_star_inside_gap;
    leading_eta_star_inside_gap = new TH1D(hist_name_in,"Leading Jet #eta* inside the gap;#eta*;#frac{d#sigma}{d#eta*} [pb]", etastar_nbins, etastar_bins);
    leading_eta_star_inside_gap->Sumw2();

    scales[26] = normalize_histogram(leading_eta_star_inside_gap, leading_eta_star_inside_gap_allvertex, leading_eta_star_inside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the delta eta outside gap distribution
    hist_name_in = hist_prefix+"delta_eta_outside_gap";
    hist_name_out = hist_prefix+"delta_eta_outside_gap";
    TH1D *delta_eta_outside_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_eta_outside_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_eta_outside_gap;
    delta_eta_outside_gap = new TH1D(hist_name_in,"#Delta#eta outside the gap;#Delta#eta;#frac{d#sigma}{d#Delta#eta} [pb]", deta_out_nbins, deta_out_bins);
    delta_eta_outside_gap->Sumw2();

    scales[27] = normalize_histogram(delta_eta_outside_gap, delta_eta_outside_gap_allvertex, delta_eta_outside_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading central pt distribution
    hist_name_in = hist_prefix+"leading_central_pt";
    hist_name_out = hist_prefix+"leading_central_pt";
    TH1D *leading_central_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_central_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_central_pt;
    leading_central_pt = new TH1D(hist_name_in,"Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
    leading_central_pt->Sumw2();

    scales[28] = normalize_histogram(leading_central_pt, leading_central_pt_allvertex, leading_central_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading central pt gap distribution
    hist_name_in = hist_prefix+"leading_central_pt_gap";
    hist_name_out = hist_prefix+"leading_central_pt_gap";
    TH1D *leading_central_pt_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_central_pt_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_central_pt_gap;
    leading_central_pt_gap = new TH1D(hist_name_in,"Leading Central Jet p_{T} when requiring a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
    leading_central_pt_gap->Sumw2();

    scales[29] = normalize_histogram(leading_central_pt_gap, leading_central_pt_gap_allvertex, leading_central_pt_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading central pt nogap distribution
    hist_name_in = hist_prefix+"leading_central_pt_nogap";
    hist_name_out = hist_prefix+"leading_central_pt_nogap";
    TH1D *leading_central_pt_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_central_pt_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_central_pt_nogap;
    leading_central_pt_nogap = new TH1D(hist_name_in,"Leading Central Jet p_{T} when vetoing a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
    leading_central_pt_nogap->Sumw2();

    scales[30] = normalize_histogram(leading_central_pt_nogap, leading_central_pt_nogap_allvertex, leading_central_pt_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading central eta distribution
    hist_name_in = hist_prefix+"leading_central_eta";
    hist_name_out = hist_prefix+"leading_central_eta";
    TH1D *leading_central_eta_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_central_eta_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_central_eta;
    leading_central_eta = new TH1D(hist_name_in,"Leading Central Jet #eta;|#eta|;#frac{d#sigma}{d#eta} [pb]", etac_nbins, etac_bins);
    leading_central_eta->Sumw2();

    scales[31] = normalize_histogram(leading_central_eta, leading_central_eta_allvertex, leading_central_eta_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading central phi distribution
    hist_name_in = hist_prefix+"leading_central_phi";
    hist_name_out = hist_prefix+"leading_central_phi";
    TH1D *leading_central_phi_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_central_phi_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_central_phi;
    leading_central_phi = new TH1D(hist_name_in,"Leading Central Jet #phi;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
    leading_central_phi->Sumw2();

    scales[32] = normalize_histogram(leading_central_phi, leading_central_phi_allvertex, leading_central_phi_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading forward pt distribution
    hist_name_in = hist_prefix+"leading_forward_pt";
    hist_name_out = hist_prefix+"leading_forward_pt";
    TH1D *leading_forward_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_forward_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_forward_pt;
    leading_forward_pt = new TH1D(hist_name_in,"Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
    leading_forward_pt->Sumw2();

    scales[33] = normalize_histogram(leading_forward_pt, leading_forward_pt_allvertex, leading_forward_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);
    
//normalizing the leading forward pt gap distribution
    hist_name_in = hist_prefix+"leading_forward_pt_gap";
    hist_name_out = hist_prefix+"leading_forward_pt_gap";
    TH1D *leading_forward_pt_gap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_forward_pt_gap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_forward_pt_gap;
    leading_forward_pt_gap = new TH1D(hist_name_in,"Leading Forward Jet p_{T} when requiring a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
    leading_forward_pt_gap->Sumw2();

    scales[34] = normalize_histogram(leading_forward_pt_gap, leading_forward_pt_gap_allvertex, leading_forward_pt_gap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading forward pt nogap distribution
    hist_name_in = hist_prefix+"leading_forward_pt_nogap";
    hist_name_out = hist_prefix+"leading_forward_pt_nogap";
    TH1D *leading_forward_pt_nogap_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_forward_pt_nogap_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_forward_pt_nogap;
    leading_forward_pt_nogap = new TH1D(hist_name_in,"Leading Forward Jet p_{T} when vetoing a gap;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
    leading_forward_pt_nogap->Sumw2();

    scales[35] = normalize_histogram(leading_forward_pt_nogap, leading_forward_pt_nogap_allvertex, leading_forward_pt_nogap_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the leading forward eta distribution
    hist_name_in = hist_prefix+"leading_forward_eta";
    hist_name_out = hist_prefix+"leading_forward_eta";
    TH1D *leading_forward_eta_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_forward_eta_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_forward_eta;
    leading_forward_eta = new TH1D(hist_name_in,"Leading Forward Jet #eta;|#eta|;#frac{d#sigma}{d#eta} [pb]", etaf_nbins, etaf_bins);
    leading_forward_eta->Sumw2();

    scales[36] = normalize_histogram(leading_forward_eta, leading_forward_eta_allvertex, leading_forward_eta_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading forward phi distribution
    hist_name_in = hist_prefix+"leading_forward_phi";
    hist_name_out = hist_prefix+"leading_forward_phi";
    TH1D *leading_forward_phi_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_forward_phi_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_forward_phi;
    leading_forward_phi = new TH1D(hist_name_in,"Leading Forward Jet #phi;#phi;#frac{d#sigma}{d#phi} [pb]", 14, -3.15, 3.15);
    leading_forward_phi->Sumw2();

    scales[37] = normalize_histogram(leading_forward_phi, leading_forward_phi_allvertex, leading_forward_phi_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the leading eta distribution
    hist_name_in = hist_prefix+"leading_eta";
    hist_name_out = hist_prefix+"leading_eta";
    TH1D *leading_eta_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_eta_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_eta;
    leading_eta = new TH1D(hist_name_in,"Leading Jet #eta;|#eta|;#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
    leading_eta->Sumw2();

    scales[38] = normalize_histogram(leading_eta, leading_eta_allvertex, leading_eta_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_middle", detail);

//normalizing the vertex selected distribution
    hist_name_in = hist_prefix+"vertex_selected";
    hist_name_out = hist_prefix+"vertex_selected";
    TH1D *vertex_selected_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *vertex_selected_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *vertex_selected;
    vertex_selected = new TH1D(hist_name_in,"Vertex multiplicity in selected events;Multiplicity;#frac{d#sigma}{dN} [pb]", 15, 0, 15);
    vertex_selected->Sumw2();

    scales[39] = normalize_histogram(vertex_selected, vertex_selected_allvertex, vertex_selected_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the jet multiplicity distribution
    hist_name_in = hist_prefix+"multiplicity";
    hist_name_out = hist_prefix+"multiplicity";
    TH1D *multiplicity_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *multiplicity_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *multiplicity;
    multiplicity = new TH1D(hist_name_in,"Jet multiplicity;Multiplicity;#frac{d#sigma}{dN} [pb]", 30, 0, 30);
    multiplicity->Sumw2();

    scales[40] = normalize_histogram(multiplicity, multiplicity_allvertex, multiplicity_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the pt all distribution
    hist_name_in = hist_prefix+"pt_all";
    hist_name_out = hist_prefix+"pt_all";
    TH1D *pt_all_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *pt_all_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *pt_all;
    pt_all = new TH1D(hist_name_in,"p_{T} all;p_{T}^{all};#frac{d#sigma}{dp_{T}} [#frac{pb}{GeV}]", all_nbins, all_bins);
    pt_all->Sumw2();

    scales[41] = normalize_histogram(pt_all, pt_all_allvertex, pt_all_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_right", detail);

//normalizing the eta all distribution
    hist_name_in = hist_prefix+"eta_all";
    hist_name_out = hist_prefix+"eta_all";
    TH1D *eta_all_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *eta_all_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *eta_all;
    eta_all = new TH1D(hist_name_in,"#eta;#eta^{all};#frac{d#sigma}{d#eta} [pb]", eta_nbins, eta_bins);
    eta_all->Sumw2();

    scales[42] = normalize_histogram(eta_all, eta_all_allvertex, eta_all_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the phi all distribution
    hist_name_in = hist_prefix+"phi_all";
    hist_name_out = hist_prefix+"phi_all";
    TH1D *phi_all_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *phi_all_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *phi_all;
    phi_all = new TH1D(hist_name_in,"#phi;#phi^{all};#frac{d#sigma}{d#phi} [#frac{pb}{rad}]", 14, -3.15, 3.15);
    phi_all->Sumw2();

    scales[43] = normalize_histogram(phi_all, phi_all_allvertex, phi_all_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the inclusive leading forward pt distribution
    hist_name_in = hist_prefix+"inclusive_leading_forward_pt";
    hist_name_out = hist_prefix+"inclusive_leading_forward_pt";
    TH1D *inclusive_leading_forward_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *inclusive_leading_forward_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *inclusive_leading_forward_pt;
    inclusive_leading_forward_pt = new TH1D(hist_name_in,"Inclusive Leading Forward Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
    inclusive_leading_forward_pt->Sumw2();

    scales[44] = normalize_histogram(inclusive_leading_forward_pt, inclusive_leading_forward_pt_allvertex, inclusive_leading_forward_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the inclusive leading central pt distribution
    hist_name_in = hist_prefix+"inclusive_leading_central_pt";
    hist_name_out = hist_prefix+"inclusive_leading_central_pt";
    TH1D *inclusive_leading_central_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *inclusive_leading_central_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *inclusive_leading_central_pt;
    inclusive_leading_central_pt = new TH1D(hist_name_in,"Inclusive Leading Central Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
    inclusive_leading_central_pt->Sumw2();

    scales[45] = normalize_histogram(inclusive_leading_central_pt, inclusive_leading_central_pt_allvertex, inclusive_leading_central_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the inclusive forward pt distribution
    hist_name_in = hist_prefix+"inclusive_forward_pt";
    hist_name_out = hist_prefix+"inclusive_forward_pt";
    TH1D *inclusive_forward_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *inclusive_forward_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *inclusive_forward_pt;
    inclusive_forward_pt = new TH1D(hist_name_in,"Inclusive Forward Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", forw_nbins, forw_bins);
    inclusive_forward_pt->Sumw2();

    scales[46] = normalize_histogram(inclusive_forward_pt, inclusive_forward_pt_allvertex, inclusive_forward_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the inclusive central pt distribution
    hist_name_in = hist_prefix+"inclusive_central_pt";
    hist_name_out = hist_prefix+"inclusive_central_pt";
    TH1D *inclusive_central_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *inclusive_central_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *inclusive_central_pt;
    inclusive_central_pt = new TH1D(hist_name_in,"Inclusive Central Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", cent_nbins, cent_bins);
    inclusive_central_pt->Sumw2();

    scales[47] = normalize_histogram(inclusive_central_pt, inclusive_central_pt_allvertex, inclusive_central_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the delta phi central relative distribution
    hist_name_in = hist_prefix+"delta_phi_central_rel";
    hist_name_out = hist_prefix+"delta_phi_central_rel";
    TH1D *delta_phi_central_rel_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_central_rel_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_central_rel;
    delta_phi_central_rel = new TH1D(hist_name_in,"#Delta#phi central relations;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_central_rel->Sumw2();

    scales[48] = normalize_histogram(delta_phi_central_rel, delta_phi_central_rel_allvertex, delta_phi_central_rel_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi forward relative distribution
    hist_name_in = hist_prefix+"delta_phi_forward_rel";
    hist_name_out = hist_prefix+"delta_phi_forward_rel";
    TH1D *delta_phi_forward_rel_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_forward_rel_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_forward_rel;
    delta_phi_forward_rel = new TH1D(hist_name_in,"#Delta#phi forward relations;|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_forward_rel->Sumw2();

    scales[49] = normalize_histogram(delta_phi_forward_rel, delta_phi_forward_rel_allvertex, delta_phi_forward_rel_1vertex, output_path_plots, plot_prefix, hist_name_out, "middle_left", detail);

//normalizing the delta phi central relative for small dpt distribution
    hist_name_in = hist_prefix+"delta_phi_central_rel_small";
    hist_name_out = hist_prefix+"delta_phi_central_rel_small";
    TH1D *delta_phi_central_rel_small_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_central_rel_small_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_central_rel_small;
    delta_phi_central_rel_small = new TH1D(hist_name_in,"#Delta#phi central relations for small #Delta#p_{T};|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_central_rel_small->Sumw2();

    scales[50] = normalize_histogram(delta_phi_central_rel_small, delta_phi_central_rel_small_allvertex, delta_phi_central_rel_small_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi forward relative for small dpt distribution
    hist_name_in = hist_prefix+"delta_phi_forward_rel_small";
    hist_name_out = hist_prefix+"delta_phi_forward_rel_small";
    TH1D *delta_phi_forward_rel_small_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_forward_rel_small_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_forward_rel_small;
    delta_phi_forward_rel_small = new TH1D(hist_name_in,"#Delta#phi forward relations for small #Delta#p_{T};|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_forward_rel_small->Sumw2();

    scales[51] = normalize_histogram(delta_phi_forward_rel_small, delta_phi_forward_rel_small_allvertex, delta_phi_forward_rel_small_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the delta phi central relative for medium dpt distribution
    hist_name_in = hist_prefix+"delta_phi_central_rel_medium";
    hist_name_out = hist_prefix+"delta_phi_central_rel_medium";
    TH1D *delta_phi_central_rel_medium_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_central_rel_medium_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_central_rel_medium;
    delta_phi_central_rel_medium = new TH1D(hist_name_in,"#Delta#phi central relations for medium #Delta#p_{T};|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_central_rel_medium->Sumw2();

    scales[52] = normalize_histogram(delta_phi_central_rel_medium, delta_phi_central_rel_medium_allvertex, delta_phi_central_rel_medium_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi forward relative for medium dpt distribution
    hist_name_in = hist_prefix+"delta_phi_forward_rel_medium";
    hist_name_out = hist_prefix+"delta_phi_forward_rel_medium";
    TH1D *delta_phi_forward_rel_medium_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_forward_rel_medium_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_forward_rel_medium;
    delta_phi_forward_rel_medium = new TH1D(hist_name_in,"#Delta#phi forward relations for medium #Delta#p_{T};|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_forward_rel_medium->Sumw2();

    scales[53] = normalize_histogram(delta_phi_forward_rel_medium, delta_phi_forward_rel_medium_allvertex, delta_phi_forward_rel_medium_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_middle", detail);

//normalizing the delta phi central relative for large dpt distribution
    hist_name_in = hist_prefix+"delta_phi_central_rel_large";
    hist_name_out = hist_prefix+"delta_phi_central_rel_large";
    TH1D *delta_phi_central_rel_large_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_central_rel_large_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_central_rel_large;
    delta_phi_central_rel_large = new TH1D(hist_name_in,"#Delta#phi central relations for large #Delta#p_{T};|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_central_rel_large->Sumw2();

    scales[54] = normalize_histogram(delta_phi_central_rel_large, delta_phi_central_rel_large_allvertex, delta_phi_central_rel_large_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

//normalizing the delta phi forward relative for large dpt distribution
    hist_name_in = hist_prefix+"delta_phi_forward_rel_large";
    hist_name_out = hist_prefix+"delta_phi_forward_rel_large";
    TH1D *delta_phi_forward_rel_large_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_forward_rel_large_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_forward_rel_large;
    delta_phi_forward_rel_large = new TH1D(hist_name_in,"#Delta#phi forward relations for large #Delta#p_{T};|#Delta#phi|;#frac{d#sigma}{d#Delta#phi} [pb]", dphi_nbins, dphi_bins);
    delta_phi_forward_rel_large->Sumw2();

    scales[55] = normalize_histogram(delta_phi_forward_rel_large, delta_phi_forward_rel_large_allvertex, delta_phi_forward_rel_large_1vertex, output_path_plots, plot_prefix, hist_name_out, "middle_left", detail);

//normalizing the delta pT central relative distribution
    hist_name_in = hist_prefix+"delta_pt_central_rel";
    hist_name_out = hist_prefix+"delta_pt_central_rel";
    TH1D *delta_pt_central_rel_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_pt_central_rel_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_pt_central_rel;
    delta_pt_central_rel = new TH1D(hist_name_in,"#Deltap_{T} central relations;|#Deltap_{T}|;#frac{d#sigma}{d#Deltap_{T}} [pb]", dpt_nbins, dpt_bins);
    delta_pt_central_rel->Sumw2();

    scales[56] = normalize_histogram(delta_pt_central_rel, delta_pt_central_rel_allvertex, delta_pt_central_rel_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//normalizing the delta pt forward relative distribution
    hist_name_in = hist_prefix+"delta_pt_forward_rel";
    hist_name_out = hist_prefix+"delta_pt_forward_rel";
    TH1D *delta_pt_forward_rel_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_pt_forward_rel_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_pt_forward_rel;
    delta_pt_forward_rel = new TH1D(hist_name_in,"#Deltap_{T} central relations;|#Deltap_{T}|;#frac{d#sigma}{d#Deltap_{T}} [pb]", dpt_nbins, dpt_bins);
    delta_pt_forward_rel->Sumw2();

    scales[57] = normalize_histogram(delta_pt_forward_rel, delta_pt_forward_rel_allvertex, delta_pt_forward_rel_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_right", detail);

//normalizing the delta phi fine distribution
    hist_name_in = hist_prefix+"delta_phi_fine";
    hist_name_out = hist_prefix+"delta_phi_fine";
    TH1D *delta_phi_fine_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *delta_phi_fine_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *delta_phi_fine;
    delta_phi_fine = new TH1D(hist_name_in,"#Delta#phi fine;|#Delta#phi| [rad];#frac{d#sigma}{d#Delta#phi}", 14, 0, 3.15);
    delta_phi_fine->Sumw2();

    scales[58] = normalize_histogram(delta_phi_fine, delta_phi_fine_allvertex, delta_phi_fine_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_left", detail);

/* //normalizing the leading pt distribution
    hist_name_in = hist_prefix+"leading_pt";
    hist_name_out = hist_prefix+"leading_pt";
    TH1D *leading_pt_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_pt_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_pt;
    leading_pt = new TH1D(hist_name_in,"Leading Jet p_{T};p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", all_nbins, all_bins);
    leading_pt->Sumw2();

    scales[59] = normalize_histogram(leading_pt, leading_pt_allvertex, leading_pt_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_right", detail);

//normalizing the leading pt fine distribution
    hist_name_in = hist_prefix+"leading_pt_fine";
    hist_name_out = hist_prefix+"leading_pt_fine";
    TH1D *leading_pt_fine_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_pt_fine_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_pt_fine;
    leading_pt_fine = new TH1D(hist_name_in,"Leading Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", 85, 30, 200);
    leading_pt_fine->Sumw2();

    scales[60] = normalize_histogram(leading_pt_fine, leading_pt_fine_allvertex, leading_pt_fine_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_right", detail); */

//normalizing the leading central pt fine distribution
    hist_name_in = hist_prefix+"leading_central_pt_fine";
    hist_name_out = hist_prefix+"leading_central_pt_fine";
    TH1D *leading_central_pt_fine_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_central_pt_fine_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_central_pt_fine;
    leading_central_pt_fine = new TH1D(hist_name_in,"Leading Central Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", 85, 30, 200);
    leading_central_pt_fine->Sumw2();

    scales[61] = normalize_histogram(leading_central_pt_fine, leading_central_pt_fine_allvertex, leading_central_pt_fine_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_right", detail);

//normalizing the leading forward pt fine distribution
    hist_name_in = hist_prefix+"leading_forward_pt_fine";
    hist_name_out = hist_prefix+"leading_forward_pt_fine";
    TH1D *leading_forward_pt_fine_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *leading_forward_pt_fine_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *leading_forward_pt_fine;
    leading_forward_pt_fine = new TH1D(hist_name_in,"Leading Forward Jet p_{T} with fine binning;p_{T} [#frac{GeV}{c}];#frac{d#sigma}{dp_{T}} [#frac{pb . c}{GeV}]", 85, 30, 200);
    leading_forward_pt_fine->Sumw2();

    scales[62] = normalize_histogram(leading_forward_pt_fine, leading_forward_pt_fine_allvertex, leading_forward_pt_fine_1vertex, output_path_plots, plot_prefix, hist_name_out, "top_right", detail);

//normalizing the pvz selected distribution
    hist_name_in = hist_prefix+"pvz_selected";
    hist_name_out = hist_prefix+"pvz_selected";
    TH1D *pvz_selected_allvertex = (TH1D*) mc_allvertex->Get(hist_name_in);
    TH1D *pvz_selected_1vertex = (TH1D*) mc_1vertex->Get(hist_name_in);
    
    TH1D *pvz_selected;
    pvz_selected = new TH1D(hist_name_in,"z of the Primary Vertex in selected events;z-position;#frac{d#sigma}{dz} [pb]", 200, -100, 100);
    pvz_selected->Sumw2();

    scales[63] = normalize_histogram(pvz_selected, pvz_selected_allvertex, pvz_selected_1vertex, output_path_plots, plot_prefix, hist_name_out, "bottom_left", detail);

//show the scale factors
    if (show_scales) { show_scaling_factors(scales); }

//recreating the output file
    TFile data_output( output_file_rootfile.c_str() , "RECREATE");

//write histograms on file
    if (detail) { cout<<"Writing histograms on file "<<output_file_rootfile<<" ..."<<endl; }
    delta_phi->Write();
    delta_phi_fine->Write();
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
    leading_eta_inside_gap->Write();
    leading_phi_inside_gap->Write();
    total_pt_outside_gap->Write();
    leading_pt_outside_gap->Write();
    leading_eta_outside_gap->Write();
    leading_phi_outside_gap->Write();
    leading_eta_star_inside_gap->Write();
    delta_eta_outside_gap->Write();

    //leading_pt->Write();
    //leading_pt_fine->Write();
    leading_central_pt_fine->Write();
    leading_forward_pt_fine->Write();
    leading_central_pt->Write();
    leading_central_pt_gap->Write();
    leading_central_pt_nogap->Write();
    leading_central_eta->Write();
    leading_central_phi->Write();
    leading_forward_pt->Write();
    leading_forward_pt_gap->Write();
    leading_forward_pt_nogap->Write();
    leading_forward_eta->Write();
    leading_forward_phi->Write();
    leading_eta->Write();
    
    vertex_selected->Write();
    pvz_selected->Write();
    multiplicity->Write();
    pt_all->Write();
    eta_all->Write();
    phi_all->Write();

    inclusive_leading_central_pt->Write();
    inclusive_leading_forward_pt->Write();
    inclusive_central_pt->Write();
    inclusive_forward_pt->Write();

    delta_phi_central_rel->Write();
    delta_phi_forward_rel->Write();
    delta_phi_central_rel_small->Write();
    delta_phi_forward_rel_small->Write();
    delta_phi_central_rel_medium->Write();
    delta_phi_forward_rel_medium->Write();
    delta_phi_central_rel_large->Write();
    delta_phi_forward_rel_large->Write();
    delta_pt_central_rel->Write();
    delta_pt_forward_rel->Write();
    // central_pt_rel->Write(); //2d histogram not normalized or passed forward
    // forward_pt_rel->Write(); //2d histogram not normalized or passed forward
    if (detail) { cout<<"Writing on "<<output_file_rootfile<<" was sucessfull!"<<endl; }

//close all TFiles
    if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    mc_allvertex->Close();
    mc_1vertex->Close();
    data_output.Close();

//delete the variables to prevent memory leak
/*    delete(delta_phi);
    delete(delta_phi_fine);
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

    delete(delta_eta);
    delete(delta_eta_gap);
    delete(delta_eta_nogap);

    delete(total_pt_inside_gap);
    delete(leading_pt_inside_gap);
    delete(leading_eta_inside_gap);
    delete(leading_phi_inside_gap);
    delete(total_pt_outside_gap);
    delete(leading_pt_outside_gap);
    delete(leading_eta_outside_gap);
    delete(leading_phi_outside_gap);
    delete(leading_eta_star_inside_gap);
    delete(delta_eta_outside_gap);

    //delete(leading_pt);
    //delete(leading_pt_fine);
    delete(leading_central_pt_fine);
    delete(leading_forward_pt_fine);
    delete(leading_central_pt);
    delete(leading_central_pt_gap);
    delete(leading_central_pt_nogap);
    delete(leading_central_eta);
    delete(leading_central_phi);
    delete(leading_forward_pt);
    delete(leading_forward_pt_gap);
    delete(leading_forward_pt_nogap);
    delete(leading_forward_eta);
    delete(leading_forward_phi);
    delete(leading_eta);

    delete(vertex_selected);
    delete(multiplicity);
    delete(pt_all);
    delete(eta_all);
    delete(phi_all);

    delete(inclusive_leading_central_pt);
    delete(inclusive_leading_forward_pt);
    delete(inclusive_central_pt);
    delete(inclusive_forward_pt);

    delete(delta_phi_central_rel);
    delete(delta_phi_forward_rel);
    delete(delta_phi_central_rel_small);
    delete(delta_phi_forward_rel_small);
    delete(delta_phi_central_rel_medium);
    delete(delta_phi_forward_rel_medium);
    delete(delta_phi_central_rel_large);
    delete(delta_phi_forward_rel_large); */

    if (detail) { cout<<"Done!"<<endl; }
}
