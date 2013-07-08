// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 20 Nov 2012
//
// control_plots(string input_path_MC_GEN, string input_path_MC_DET, string input_path_data, string output_path_plots, bool detail, bool disp_output)
// makes control plots of the uncorrected data

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TPaveText.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "common_methods.h"

void plot_histogram(TH1D *p6_z2_gen, TH1D *p8_4c_gen, TH1D *p6_z2_det, TH1D *p8_4c_det, TH1D *data, string path, string fileout, string legend_position = "top_left", string label_position = "top_right", TString scenario = "INCLUSIVE", bool detail = false)
{
//plots the model uncertainty control plots

//declaring the canvas
    if (detail) { cout << "Ploting " << fileout << endl; }
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,800,600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.07);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogy();

    p6_z2_det->GetXaxis()->SetLabelSize(0.04);
    p6_z2_det->GetYaxis()->SetLabelSize(0.05);
    p6_z2_det->GetXaxis()->SetTitleSize(0.04);
    p6_z2_det->GetYaxis()->SetTitleSize(0.05);

//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 0.0;
    double max = p6_z2_det->GetMaximum();
    if (p6_z2_det->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(p6_z2_det,detail);
    }
    else
    {
    min = p6_z2_det->GetMinimum();
    } 

    //set_histogram_min_max(p8_4c_gen, min, max, detail);
    //set_histogram_min_max(p6_z2_det, min, max, detail);
    set_histogram_min_max(p8_4c_det, min, max, detail);
    set_histogram_min_max(data, min, max, detail);
    
    max = 1.3 * max;
    min = 0.7 * min;

    cout << "max = " << max << " min = " << min << endl;
//format and ploting the histogram
    if (detail) { cout << "Drawning on the canvas..." << endl; }
    //p6_z2_gen->SetLineColor(2);
    //p6_z2_gen->SetLineStyle(1);
    //p6_z2_gen->SetLineWidth(3);
    //p6_z2_gen->Draw("e1");
    //p8_4c_gen->SetLineColor(1);
    //p8_4c_gen->SetLineStyle(2);
    //p8_4c_gen->SetLineWidth(3);
    //p8_4c_gen->Draw("e1 same");
    p6_z2_det->SetMaximum(max);
    p6_z2_det->SetMinimum(min);
    p6_z2_det->SetLineColor(2);
    p6_z2_det->SetLineStyle(2);
    p6_z2_det->SetLineWidth(4);
    p6_z2_det->Draw("e1");
    p8_4c_det->SetLineColor(4);
    p8_4c_det->SetLineStyle(3);
    p8_4c_det->SetLineWidth(4);
    p8_4c_det->Draw("e1 same");
    data->SetLineColor(1);
    data->SetLineStyle(1);
    data->SetLineWidth(4);
    data->SetMarkerColor(1);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(2);
    data->Draw("e1 same");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1,y1-0.06,x2,y2-0.06);
    //leg00->AddEntry(p6_z2_gen,"Pythia6 Z2* Tune - Generator Level","l");
    //leg00->AddEntry(p8_4c_gen,"Pythia8 4C Tune - Generator Level","l");
    leg00->AddEntry(p6_z2_det,"Pythia6 Z2* Tune","l");
    leg00->AddEntry(p8_4c_det,"Pythia8 4C Tune","l");
    leg00->AddEntry(data,"Data 2010","lep");
    leg00->SetFillColor(0);
    leg00->SetLineStyle(1);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();

//    if (mode == "xsec") { x1 = 0.65; y1 = 0.13; x2 = 0.97; y2 = 0.23; }
    if (label_position == "top_right")    { x1 = 0.60; y1 = 0.77; x2 = 0.98; y2 = 0.92; }
    if (label_position == "top_left")     { x1 = 0.13; y1 = 0.77; x2 = 0.49; y2 = 0.92; }
    if (label_position == "bottom_left")  { x1 = 0.13; y1 = 0.13; x2 = 0.49; y2 = 0.28; }
    if (label_position == "bottom_right") { x1 = 0.60; y1 = 0.13; x2 = 0.98; y2 = 0.28; }
    if (label_position == "middle")       { x1 = 0.50; y1 = 0.38; x2 = 0.60; y2 = 0.53; }
    if (label_position == "top_middle")   { x1 = 0.50; y1 = 0.77; x2 = 0.60; y2 = 0.92; }

    TPaveText *cms_topleft = new TPaveText(0.00,0.99,0.44,0.94,"NDC"); // NDC sets coords
    cms_topleft->SetTextSize(0.04);
    cms_topleft->SetBorderSize(0); 
    cms_topleft->SetTextAlign(12);
    cms_topleft->SetLineWidth(1);
    cms_topleft->SetLineColor(0);
    cms_topleft->SetFillColor(0);
    cms_topleft->SetFillStyle(1001);
    cms_topleft->AddText("CMS Preliminary, pp #rightarrow 2 jets + X");
    cms_topleft->Draw();


    TPaveText *cms_topmiddle = new TPaveText(0.44,0.99,0.75,0.94,"NDC"); // NDC sets coords
    cms_topmiddle->SetTextSize(0.04);
    cms_topmiddle->SetBorderSize(0); 
    cms_topmiddle->SetTextAlign(12);
    cms_topmiddle->SetLineWidth(1);
    cms_topmiddle->SetLineColor(0);
    cms_topmiddle->SetFillColor(0);
    cms_topmiddle->SetFillStyle(1001);
    cms_topmiddle->AddText("["+scenario+"]");
    cms_topmiddle->Draw();

    TPaveText *cms_topright = new TPaveText(0.85,0.99,1.00,0.94,"NDC"); // NDC sets coords
    cms_topright->SetTextSize(0.04);
    cms_topright->SetBorderSize(0); 
    cms_topright->SetTextAlign(12);
    cms_topright->SetLineWidth(1);
    cms_topright->SetLineColor(0);
    cms_topright->SetFillColor(0);
    cms_topright->SetFillStyle(1001);
    cms_topright->AddText("#sqrt{s} = 7 TeV");
    cms_topright->Draw();

    TPaveText *cms = new TPaveText(x1,y1,x2,y2,"NDC"); // NDC sets coords
    cms->SetTextSize(0.04);
    cms->SetBorderSize(0); 
    cms->SetTextAlign(22);
    cms->SetLineWidth(1);
    cms->SetLineColor(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(1001);
//    cms->AddText("L_{int} = 34.5 pb^{-1}, Anti-k_{T} (R = 0.5)");
    cms->AddText("L_{int} = 3.2 pb^{-1}, Anti-k_{T} (R = 0.5)");
    cms->AddText("p_{T}^{jet} > 35 GeV and |#eta| < 2.8");
    cms->AddText("p_{T}^{jet} > 35 GeV and 3.2 < |#eta| < 4.7");
    cms->Draw();

    print_plots(c1, path, fileout);

}

void control_plots(string path_p6_z2_gen = " ../output/histograms/normalized_mc/xsec_p6_z2_gen.root", string path_p8_4c_gen = " ../output/histograms/normalized_mc/xsec_p8_4c_gen.root", string path_p6_z2_det = "../output/histograms/xsec_mc_det/xsec_Pythia6_TuneZ2star_det_allvertex.root", string path_p8_4c_det = "../output/histograms/xsec_mc_det/xsec_Pythia8_Tune4C_det_allvertex.root", string input_data = "../output/histograms/xsec_data_2010.root", string output_path_plots = "../output/control_dist/", string plot_prefix = "control_", bool detail = false, bool disp_output = true)
{
//plots the mc generator level cross-section

//outputs the configuration
    if (detail) { cout << "Plot Control Distributions Configuration"<<endl; }
    if (detail) { cout << "Input path for Pythia 6 - Tune Z2* on Generator Level: " << path_p6_z2_gen << endl; }
    if (detail) { cout << "Input path for Pythia 8 - Tune 4C on Generator Level:  " << path_p8_4c_gen << endl; }
    if (detail) { cout << "Input path for Pythia 6 - Tune Z2* on Detector Level:  " << path_p6_z2_det << endl; }
    if (detail) { cout << "Input path for Pythia 8 - Tune 4C on Detector Level:   " << path_p8_4c_det << endl; }
    if (detail) { cout << "Input path for Data:                                   " << input_data << endl; }
    if (detail) { cout << "Output Path Plots:                                     " << output_path_plots << endl; }
    if (detail) { cout << "Plot Prefix:                                           " << plot_prefix << endl; }
    if (detail) { cout << "Detail Level:                                          " << detail << endl; }
    if (detail) { cout << "Display Results:                                       " << disp_output << endl; }


    if (detail) { cout << "Opening files..." << endl; }
//opens the MC files
    if (detail) { cout << "Opening MC files" << endl; }
    TFile *p6_z2_gen = new TFile( path_p6_z2_gen.c_str() );
    TFile *p8_4c_gen = new TFile( path_p8_4c_gen.c_str() );
    TFile *p6_z2_det = new TFile( path_p6_z2_det.c_str() );
    TFile *p8_4c_det = new TFile( path_p8_4c_det.c_str() );

//opens the data files
    if (detail) { cout << "Opening data files" << endl; }
    TFile *data_file = new TFile( input_data.c_str() );

    if (detail) { cout << "All files opened sucessfully!" << endl; }

//plot delta phi distribution
    if (detail) { cout << "Ploting Delta Phi..." << endl; }
    TH1D *delta_phi_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi");
    TH1D *delta_phi_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi");
    TH1D *delta_phi_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi");
    TH1D *delta_phi_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi");
    TH1D *delta_phi_data = (TH1D*) data_file->Get("ak5PF_delta_phi");

    delta_phi_det_p6_z2->SetTitle("#Delta#phi;#Delta#phi [rad];dN/d#Delta#phi [rad^{-1}]");
    plot_histogram(delta_phi_gen_p6_z2, delta_phi_gen_p8_4c, delta_phi_det_p6_z2, delta_phi_det_p8_4c, delta_phi_data, output_path_plots, plot_prefix + "xsec_all_delta_phi", "top_left", "bottom_right", "INCLUSIVE", detail);


//plot delta phi fine distribution
    if (detail) { cout << "Ploting Delta Phi Fine..." << endl; }
    TH1D *delta_phi_fine_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_fine");
    TH1D *delta_phi_fine_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_fine");
    TH1D *delta_phi_fine_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_fine");
    TH1D *delta_phi_fine_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_fine");
    TH1D *delta_phi_fine_data = (TH1D*) data_file->Get("ak5PF_delta_phi_fine");

    plot_histogram(delta_phi_fine_gen_p6_z2, delta_phi_fine_gen_p8_4c, delta_phi_fine_det_p6_z2, delta_phi_fine_det_p8_4c, delta_phi_fine_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_fine", "top_left", "bottom_right", "INCLUSIVE", detail);


//plot delta phi deta1 distribution
    if (detail) { cout << "Ploting Delta Phi Deta1..." << endl; }
    TH1D *delta_phi_deta1_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta1");
    TH1D *delta_phi_deta1_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta1");
    TH1D *delta_phi_deta1_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta1");
    TH1D *delta_phi_deta1_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta1");
    TH1D *delta_phi_deta1_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta1");

    plot_histogram(delta_phi_deta1_gen_p6_z2, delta_phi_deta1_gen_p8_4c, delta_phi_deta1_det_p6_z2, delta_phi_deta1_det_p8_4c, delta_phi_deta1_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta1", "top_left", "bottom_right", "INCLUSIVE", detail);

//plot delta phi deta2 distribution
    if (detail) { cout << "Ploting Delta Phi Deta1..." << endl; }
    TH1D *delta_phi_deta2_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta2");
    TH1D *delta_phi_deta2_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta2");
    TH1D *delta_phi_deta2_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta2");
    TH1D *delta_phi_deta2_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta2");
    TH1D *delta_phi_deta2_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta2");

    plot_histogram(delta_phi_deta2_gen_p6_z2, delta_phi_deta2_gen_p8_4c, delta_phi_deta2_det_p6_z2, delta_phi_deta2_det_p8_4c, delta_phi_deta2_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta2", "top_left", "bottom_right", "INCLUSIVE", detail);

//plot delta phi deta3 distribution
    if (detail) { cout << "Ploting Delta Phi Deta3..." << endl; }
    TH1D *delta_phi_deta3_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta3");
    TH1D *delta_phi_deta3_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta3");
    TH1D *delta_phi_deta3_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta3");
    TH1D *delta_phi_deta3_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta3");
    TH1D *delta_phi_deta3_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta3");

    plot_histogram(delta_phi_deta3_gen_p6_z2, delta_phi_deta3_gen_p8_4c, delta_phi_deta3_det_p6_z2, delta_phi_deta3_det_p8_4c, delta_phi_deta3_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta3", "top_left", "bottom_right", "INCLUSIVE", detail);

//plot delta phi deta4 distribution
    if (detail) { cout << "Ploting Delta Phi Deta4..." << endl; }
    TH1D *delta_phi_deta4_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta4");
    TH1D *delta_phi_deta4_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta4");
    TH1D *delta_phi_deta4_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta4");
    TH1D *delta_phi_deta4_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta4");
    TH1D *delta_phi_deta4_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta4");

    plot_histogram(delta_phi_deta4_gen_p6_z2, delta_phi_deta4_gen_p8_4c, delta_phi_deta4_det_p6_z2, delta_phi_deta4_det_p8_4c, delta_phi_deta4_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta4", "top_left", "bottom_right", "INCLUSIVE", detail);

//plot delta phi gap distribution
    if (detail) { cout << "Ploting Delta Phi Gap..." << endl; }
    TH1D *delta_phi_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_gap");
    TH1D *delta_phi_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_gap");
    TH1D *delta_phi_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_gap");
    TH1D *delta_phi_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_gap");
    TH1D *delta_phi_gap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_gap");

    delta_phi_gap_det_p6_z2->SetTitle("#Delta#phi;#Delta#phi [rad];dN/d#Delta#phi [rad^{-1}]");
    plot_histogram(delta_phi_gap_gen_p6_z2, delta_phi_gap_gen_p8_4c, delta_phi_gap_det_p6_z2, delta_phi_gap_det_p8_4c, delta_phi_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_gap", "top_left", "bottom_right", "INSIDE-JET VETO", detail);

//plot delta phi deta1 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta1 Gap..." << endl; }
    TH1D *delta_phi_deta1_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta1_gap");

    plot_histogram(delta_phi_deta1_gap_gen_p6_z2, delta_phi_deta1_gap_gen_p8_4c, delta_phi_deta1_gap_det_p6_z2, delta_phi_deta1_gap_det_p8_4c, delta_phi_deta1_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta1_gap", "top_left", "bottom_right", "INSIDE-JET VETO", detail);

//plot delta phi deta2 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta2 Gap..." << endl; }
    TH1D *delta_phi_deta2_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta2_gap");

    plot_histogram(delta_phi_deta2_gap_gen_p6_z2, delta_phi_deta2_gap_gen_p8_4c, delta_phi_deta2_gap_det_p6_z2, delta_phi_deta2_gap_det_p8_4c, delta_phi_deta2_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta2_gap", "top_left", "bottom_right", "INSIDE-JET VETO", detail);

//plot delta phi deta3 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta3 Gap..." << endl; }
    TH1D *delta_phi_deta3_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta3_gap");

    plot_histogram(delta_phi_deta3_gap_gen_p6_z2, delta_phi_deta3_gap_gen_p8_4c, delta_phi_deta3_gap_det_p6_z2, delta_phi_deta3_gap_det_p8_4c, delta_phi_deta3_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta3_gap", "top_left", "bottom_right", "INSIDE-JET VETO", detail);

//plot delta phi deta4 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta4 Gap..." << endl; }
    TH1D *delta_phi_deta4_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta4_gap");

    plot_histogram(delta_phi_deta4_gap_gen_p6_z2, delta_phi_deta4_gap_gen_p8_4c, delta_phi_deta4_gap_det_p6_z2, delta_phi_deta4_gap_det_p8_4c, delta_phi_deta4_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta4_gap", "top_left", "bottom_right", "INSIDE-JET VETO", detail);


//plot delta eta distribution
    if (detail) { cout << "Ploting Delta Eta..." << endl; }
    TH1D *delta_eta_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_eta");
    TH1D *delta_eta_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_eta");
    TH1D *delta_eta_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_eta");
    TH1D *delta_eta_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_eta");
    TH1D *delta_eta_data = (TH1D*) data_file->Get("ak5PF_delta_eta");

    plot_histogram(delta_eta_gen_p6_z2, delta_eta_gen_p8_4c, delta_eta_det_p6_z2, delta_eta_det_p8_4c, delta_eta_data, output_path_plots, plot_prefix + "xsec_all_delta_eta", "bottom_left", "top_right", "INCLUSIVE", detail);


//plot delta eta gap distribution
    if (detail) { cout << "Ploting Delta Eta Gap..." << endl; }
    TH1D *delta_eta_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_eta_gap");
    TH1D *delta_eta_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_eta_gap");
    TH1D *delta_eta_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_eta_gap");
    TH1D *delta_eta_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_eta_gap");
    TH1D *delta_eta_gap_data = (TH1D*) data_file->Get("ak5PF_delta_eta_gap");

    plot_histogram(delta_eta_gap_gen_p6_z2, delta_eta_gap_gen_p8_4c, delta_eta_gap_det_p6_z2, delta_eta_gap_det_p8_4c, delta_eta_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_eta_gap", "bottom_left", "top_right", "INSIDE-JET VETO", detail);


//plot delta eta nogap distribution
    if (detail) { cout << "Ploting Delta Eta noGap..." << endl; }
    TH1D *delta_eta_nogap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_eta_nogap");
    TH1D *delta_eta_nogap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_eta_nogap");
    TH1D *delta_eta_nogap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_eta_nogap");
    TH1D *delta_eta_nogap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_eta_nogap");
    TH1D *delta_eta_nogap_data = (TH1D*) data_file->Get("ak5PF_delta_eta_nogap");

    plot_histogram(delta_eta_nogap_gen_p6_z2, delta_eta_nogap_gen_p8_4c, delta_eta_nogap_det_p6_z2, delta_eta_nogap_det_p8_4c, delta_eta_nogap_data, output_path_plots, plot_prefix + "xsec_all_delta_eta_nogap", "bottom_left", "top_right", "INSIDE-JET TAG", detail);


//plot delta phi nogap distribution
    if (detail) { cout << "Ploting Delta Phi noGap..." << endl; }
    TH1D *delta_phi_nogap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_nogap");
    TH1D *delta_phi_nogap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_nogap");
    TH1D *delta_phi_nogap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_nogap");
    TH1D *delta_phi_nogap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_nogap");
    TH1D *delta_phi_nogap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_nogap");

    delta_phi_nogap_det_p6_z2->SetTitle("#Delta#phi;#Delta#phi [rad];dN/d#Delta#phi [rad^{-1}]");
    plot_histogram(delta_phi_nogap_gen_p6_z2, delta_phi_nogap_gen_p8_4c, delta_phi_nogap_det_p6_z2, delta_phi_nogap_det_p8_4c, delta_phi_nogap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_nogap", "top_left", "bottom_right", "INSIDE-JET TAG", detail);

//plot delta phi deta1 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta1 noGap..." << endl; }
    TH1D *delta_phi_deta1_nogap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta1_nogap");

    plot_histogram(delta_phi_deta1_nogap_gen_p6_z2, delta_phi_deta1_nogap_gen_p8_4c, delta_phi_deta1_nogap_det_p6_z2, delta_phi_deta1_nogap_det_p8_4c, delta_phi_deta1_nogap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta1_nogap", "top_left", "bottom_right", "INSIDE-JET TAG", detail);

//plot delta phi deta2 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta2 noGap..." << endl; }
    TH1D *delta_phi_deta2_nogap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta2_nogap");

    plot_histogram(delta_phi_deta2_nogap_gen_p6_z2, delta_phi_deta2_nogap_gen_p8_4c, delta_phi_deta2_nogap_det_p6_z2, delta_phi_deta2_nogap_det_p8_4c, delta_phi_deta2_nogap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta2_nogap", "top_left", "bottom_right", "INSIDE-JET TAG", detail);

//plot delta phi deta3 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta3 noGap..." << endl; }
    TH1D *delta_phi_deta3_nogap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta3_nogap");

    plot_histogram(delta_phi_deta3_nogap_gen_p6_z2, delta_phi_deta3_nogap_gen_p8_4c, delta_phi_deta3_nogap_det_p6_z2, delta_phi_deta3_nogap_det_p8_4c, delta_phi_deta3_nogap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta3_nogap", "top_left", "bottom_right", "INSIDE-JET TAG", detail);

//plot delta phi deta4 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta4 noGap..." << endl; }
    TH1D *delta_phi_deta4_nogap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_data = (TH1D*) data_file->Get("ak5PF_delta_phi_deta4_nogap");

    plot_histogram(delta_phi_deta4_nogap_gen_p6_z2, delta_phi_deta4_nogap_gen_p8_4c, delta_phi_deta4_nogap_det_p6_z2, delta_phi_deta4_nogap_det_p8_4c, delta_phi_deta4_nogap_data, output_path_plots, plot_prefix + "xsec_all_delta_phi_deta4_nogap", "top_left", "bottom_right", "INSIDE-JET TAG", detail);


//plot leading forward pT distribution
    if (detail) { cout << "Ploting Leading Forward pT..." << endl; }
    TH1D *leading_forward_pt_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_forward_pt");
    TH1D *leading_forward_pt_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_forward_pt");
    TH1D *leading_forward_pt_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_forward_pt");
    TH1D *leading_forward_pt_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_forward_pt");
    TH1D *leading_forward_pt_data = (TH1D*) data_file->Get("ak5PF_leading_forward_pt");

    plot_histogram(leading_forward_pt_gen_p6_z2, leading_forward_pt_gen_p8_4c, leading_forward_pt_det_p6_z2, leading_forward_pt_det_p8_4c, leading_forward_pt_data, output_path_plots, plot_prefix + "xsec_all_leading_forward_pt", "bottom_left", "top_right", "INCLUSIVE", detail);
    
    
//plot leading forward pT fine distribution
    if (detail) { cout << "Ploting Leading Forward pT fine..." << endl; }
    TH1D *leading_forward_pt_fine_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_forward_pt_fine");
    TH1D *leading_forward_pt_fine_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_forward_pt_fine");
    TH1D *leading_forward_pt_fine_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_forward_pt_fine");
    TH1D *leading_forward_pt_fine_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_forward_pt_fine");
    TH1D *leading_forward_pt_fine_data = (TH1D*) data_file->Get("ak5PF_leading_forward_pt_fine");

    plot_histogram(leading_forward_pt_fine_gen_p6_z2, leading_forward_pt_fine_gen_p8_4c, leading_forward_pt_fine_det_p6_z2, leading_forward_pt_fine_det_p8_4c, leading_forward_pt_fine_data, output_path_plots, plot_prefix + "xsec_all_leading_forward_pt_fine", "bottom_left", "top_right", "INCLUSIVE",  detail);


//plot leading forward eta distribution
    if (detail) { cout << "Ploting Leading Forward eta..." << endl; }
    TH1D *leading_forward_eta_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_forward_eta");
    TH1D *leading_forward_eta_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_forward_eta");
    TH1D *leading_forward_eta_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_forward_eta");
    TH1D *leading_forward_eta_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_forward_eta");
    TH1D *leading_forward_eta_data = (TH1D*) data_file->Get("ak5PF_leading_forward_eta");

    plot_histogram(leading_forward_eta_gen_p6_z2, leading_forward_eta_gen_p8_4c, leading_forward_eta_det_p6_z2, leading_forward_eta_det_p8_4c, leading_forward_eta_data, output_path_plots, plot_prefix + "xsec_all_leading_forward_eta", "bottom_middle", "middle", "INCLUSIVE", detail);


//plot leading forward phi distribution
    if (detail) { cout << "Ploting Leading Forward phi..." << endl; }
    TH1D *leading_forward_phi_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_forward_phi");
    TH1D *leading_forward_phi_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_forward_phi");
    TH1D *leading_forward_phi_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_forward_phi");
    TH1D *leading_forward_phi_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_forward_phi");
    TH1D *leading_forward_phi_data = (TH1D*) data_file->Get("ak5PF_leading_forward_phi");

    plot_histogram(leading_forward_phi_gen_p6_z2, leading_forward_phi_gen_p8_4c, leading_forward_phi_det_p6_z2, leading_forward_phi_det_p8_4c, leading_forward_phi_data, output_path_plots, plot_prefix + "xsec_all_leading_forward_phi", "bottom_middle", "top_middle", "INCLUSIVE", detail);


//plot leading central pT fine distribution
    if (detail) { cout << "Ploting Leading Central pT..." << endl; }
    TH1D *leading_central_pt_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_central_pt");
    TH1D *leading_central_pt_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_central_pt");
    TH1D *leading_central_pt_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_central_pt");
    TH1D *leading_central_pt_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_central_pt");
    TH1D *leading_central_pt_data = (TH1D*) data_file->Get("ak5PF_leading_central_pt");

    plot_histogram(leading_central_pt_gen_p6_z2, leading_central_pt_gen_p8_4c, leading_central_pt_det_p6_z2, leading_central_pt_det_p8_4c, leading_central_pt_data, output_path_plots, plot_prefix + "xsec_all_leading_central_pt", "bottom_left", "top_right", "INCLUSIVE", detail);
    
    //plot leading central pT distribution
    if (detail) { cout << "Ploting Leading Central pT fine ..." << endl; }
    TH1D *leading_central_pt_fine_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_central_pt_fine");
    TH1D *leading_central_pt_fine_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_central_pt_fine");
    TH1D *leading_central_pt_fine_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_central_pt_fine");
    TH1D *leading_central_pt_fine_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_central_pt_fine");
    TH1D *leading_central_pt_fine_data = (TH1D*) data_file->Get("ak5PF_leading_central_pt_fine");

    plot_histogram(leading_central_pt_fine_gen_p6_z2, leading_central_pt_fine_gen_p8_4c, leading_central_pt_fine_det_p6_z2, leading_central_pt_fine_det_p8_4c, leading_central_pt_fine_data, output_path_plots, plot_prefix + "xsec_all_leading_central_pt_fine", "bottom_left", "top_right", "INCLUSIVE", detail);

//plot leading central eta distribution
    if (detail) { cout << "Ploting Leading Central eta..." << endl; }
    TH1D *leading_central_eta_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_central_eta");
    TH1D *leading_central_eta_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_central_eta");
    TH1D *leading_central_eta_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_central_eta");
    TH1D *leading_central_eta_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_central_eta");
    TH1D *leading_central_eta_data = (TH1D*) data_file->Get("ak5PF_leading_central_eta");

    plot_histogram(leading_central_eta_gen_p6_z2, leading_central_eta_gen_p8_4c, leading_central_eta_det_p6_z2, leading_central_eta_det_p8_4c, leading_central_eta_data, output_path_plots, plot_prefix + "xsec_all_leading_central_eta", "bottom_middle", "middle", "INCLUSIVE", detail);


//plot leading central phi distribution
    if (detail) { cout << "Ploting Leading Central phi..." << endl; }
    TH1D *leading_central_phi_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_central_phi");
    TH1D *leading_central_phi_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_central_phi");
    TH1D *leading_central_phi_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_central_phi");
    TH1D *leading_central_phi_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_central_phi");
    TH1D *leading_central_phi_data = (TH1D*) data_file->Get("ak5PF_leading_central_phi");

    plot_histogram(leading_central_phi_gen_p6_z2, leading_central_phi_gen_p8_4c, leading_central_phi_det_p6_z2, leading_central_phi_det_p8_4c, leading_central_phi_data, output_path_plots, plot_prefix + "xsec_all_leading_central_phi", "bottom_middle", "top_middle", "INCLUSIVE", detail);


//multiplicity selected distribution
    if (detail) { cout << "Ploting Vertex Selected..." << endl; }
    TH1D *multiplicity_gen_p6_z2 = 0;
    TH1D *multiplicity_gen_p8_4c = 0;
    TH1D *multiplicity_det_p6_z2 = 0;
    TH1D *multiplicity_det_p8_4c = 0;
    TH1D *multiplicity_data = 0;
    p6_z2_gen->GetObject("ak5Gen_multiplicity",multiplicity_gen_p6_z2);
    if (multiplicity_gen_p6_z2 == 0) { cout << "ak5Gen_multiplicity in Pythia 6 Z2* generator not found!" << endl; return; }
    p8_4c_gen->GetObject("ak5Gen_multiplicity",multiplicity_gen_p8_4c);
    if (multiplicity_gen_p8_4c == 0) { cout << "ak5Gen_multiplicity in Pythia 8 4C generator not found!" << endl; return; }
    p6_z2_det->GetObject("ak5PF_multiplicity",multiplicity_det_p6_z2);
    if (multiplicity_det_p6_z2 == 0) { cout << "ak5PF_multiplicity in Pythia 6 Z2* detector not found!" << endl; return; }
    p8_4c_det->GetObject("ak5PF_multiplicity",multiplicity_det_p8_4c);
    if (multiplicity_det_p8_4c == 0) { cout << "ak5PF_multiplicity in Pythia 8 4C detector not found!" << endl; return; }
    data_file->GetObject("ak5PF_multiplicity",multiplicity_data);
    if (multiplicity_data == 0) { cout << "ak5PF_multiplicity in data detector not found!" << endl; return; }

    plot_histogram(multiplicity_gen_p6_z2, multiplicity_gen_p8_4c, multiplicity_det_p6_z2, multiplicity_det_p8_4c, multiplicity_data, output_path_plots, plot_prefix + "xsec_all_multiplicity", "bottom_left", "bottom_right", "INCLUSIVE", detail);

//plot vertex selected distribution
    if (detail) { cout << "Ploting Vertex Selected..." << endl; }
    TH1D *vertex_selected_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_vertex_selected");
    TH1D *vertex_selected_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_vertex_selected");
    TH1D *vertex_selected_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_vertex_selected");
    TH1D *vertex_selected_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_vertex_selected");
    TH1D *vertex_selected_data = (TH1D*) data_file->Get("ak5PF_vertex_selected");

    plot_histogram(vertex_selected_gen_p6_z2, vertex_selected_gen_p8_4c, vertex_selected_det_p6_z2, vertex_selected_det_p8_4c, vertex_selected_data, output_path_plots, plot_prefix + "xsec_all_vertex_selected", "bottom_left", "top_right", "INCLUSIVE", detail);
    
//plot pvz selected distribution
    if (detail) { cout << "Ploting pvz Selected..." << endl; }
    TH1D *pvz_selected_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_pvz_selected");
    TH1D *pvz_selected_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_pvz_selected");
    TH1D *pvz_selected_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_pvz_selected");
    TH1D *pvz_selected_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_pvz_selected");
    TH1D *pvz_selected_data = (TH1D*) data_file->Get("ak5PF_pvz_selected");

    plot_histogram(pvz_selected_gen_p6_z2, pvz_selected_gen_p8_4c, pvz_selected_det_p6_z2, pvz_selected_det_p8_4c, pvz_selected_data, output_path_plots, plot_prefix + "xsec_all_pvz_selected", "bottom_left", "bottom_right", "INCLUSIVE", detail);


//plot leading_pt_inside_gap distribution
    if (detail) { cout << "Ploting leading_pt_inside_gap..." << endl; }
    TH1D *leading_pt_inside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_pt_inside_gap");

    leading_pt_inside_gap_det_p6_z2->SetTitle("Inside-jet Leading p_{T};p_{T}^{inside} [GeV];dN/dp_{T}^{inside} [GeV^{-1}]");
    plot_histogram(leading_pt_inside_gap_gen_p6_z2, leading_pt_inside_gap_gen_p8_4c, leading_pt_inside_gap_det_p6_z2, leading_pt_inside_gap_det_p8_4c, leading_pt_inside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_pt_inside_gap", "bottom_left", "top_right", "INSIDE-JET TAG", detail);


//plot total_pt_inside_gap distribution
    if (detail) { cout << "Ploting total_pt_inside_gap..." << endl; }
    TH1D *total_pt_inside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_data = (TH1D*) data_file->Get("ak5PF_total_pt_inside_gap");

    total_pt_inside_gap_det_p6_z2->SetTitle("Inside-jet Total p_{T};p_{T}^{inside} [GeV];dN/dp_{T}^{inside} [GeV^{-1}]");
    plot_histogram(total_pt_inside_gap_gen_p6_z2, total_pt_inside_gap_gen_p8_4c, total_pt_inside_gap_det_p6_z2, total_pt_inside_gap_det_p8_4c, total_pt_inside_gap_data, output_path_plots, plot_prefix + "xsec_all_total_pt_inside_gap", "bottom_left", "top_right", "INSIDE-JET TAG", detail);


//plot leading_eta_inside_gap distribution
    if (detail) { cout << "Ploting leading_eta_inside_gap..." << endl; }
    TH1D *leading_eta_inside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_eta_inside_gap");

    leading_eta_inside_gap_det_p6_z2->SetTitle("Inside-jet Leading #eta;#eta^{inside};dN/d#eta^{inside}");
    plot_histogram(leading_eta_inside_gap_gen_p6_z2, leading_eta_inside_gap_gen_p8_4c, leading_eta_inside_gap_det_p6_z2, leading_eta_inside_gap_det_p8_4c, leading_eta_inside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_eta_inside_gap", "bottom_middle", "middle", "INSIDE-JET TAG", detail);


//plot leading_phi_inside_gap distribution
    if (detail) { cout << "Ploting leading_phi_inside_gap..." << endl; }
    TH1D *leading_phi_inside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_phi_inside_gap");

    leading_phi_inside_gap_det_p6_z2->SetTitle("Inside-jet Leading #phi;#phi^{inside} [rad];dN/d#phi^{inside}");
    plot_histogram(leading_phi_inside_gap_gen_p6_z2, leading_phi_inside_gap_gen_p8_4c, leading_phi_inside_gap_det_p6_z2, leading_phi_inside_gap_det_p8_4c, leading_phi_inside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_phi_inside_gap", "bottom_left", "bottom_right", "INSIDE-JET TAG", detail);


//plot multiplicity_inside_gap distribution
    if (detail) { cout << "Ploting multiplicity_inside_gap..." << endl; }
    TH1D *multiplicity_inside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_multiplicity_inside_gap");
    TH1D *multiplicity_inside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_multiplicity_inside_gap");
    TH1D *multiplicity_inside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_multiplicity_inside_gap");
    TH1D *multiplicity_inside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_multiplicity_inside_gap");
    TH1D *multiplicity_inside_gap_data = (TH1D*) data_file->Get("ak5PF_multiplicity_inside_gap");

    multiplicity_inside_gap_det_p6_z2->SetTitle("Multiplicity Inside;Multiplicity^{inside};dN/dMultiplicity^{inside}");
    plot_histogram(multiplicity_inside_gap_gen_p6_z2, multiplicity_inside_gap_gen_p8_4c, multiplicity_inside_gap_det_p6_z2, multiplicity_inside_gap_det_p8_4c, multiplicity_inside_gap_data, output_path_plots, plot_prefix + "xsec_all_multiplicity_inside_gap", "bottom_left", "top_right", "INSIDE-JET TAG", detail);


//plot leading_eta_star_inside_gap distribution
    if (detail) { cout << "Ploting leading_eta_star_inside_gap..." << endl; }
    TH1D *leading_eta_star_inside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_eta_star_inside_gap");

    leading_eta_star_inside_gap_det_p6_z2->SetTitle("Inside-jet Leading #eta*;#eta*;dN/d#eta*");
    plot_histogram(leading_eta_star_inside_gap_gen_p6_z2, leading_eta_star_inside_gap_gen_p8_4c, leading_eta_star_inside_gap_det_p6_z2, leading_eta_star_inside_gap_det_p8_4c, leading_eta_star_inside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_eta_star_inside_gap", "bottom_middle", "middle", "INSIDE-JET TAG", detail);


//plot leading_pt_outside_gap distribution
    if (detail) { cout << "Ploting leading_pt_outside_gap..." << endl; }
    TH1D *leading_pt_outside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_pt_outside_gap");

    leading_pt_outside_gap_det_p6_z2->SetTitle("Outside-jet Leading p_{T};p_{T}^{outside} [GeV];dN/dp_{T}^{outside} [GeV^{-1}]");
    plot_histogram(leading_pt_outside_gap_gen_p6_z2, leading_pt_outside_gap_gen_p8_4c, leading_pt_outside_gap_det_p6_z2, leading_pt_outside_gap_det_p8_4c, leading_pt_outside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_pt_outside_gap", "bottom_left", "top_right", "OUTSIDE-JET TAG", detail);


//plot total_pt_outside_gap distribution
    if (detail) { cout << "Ploting total_pt_outside_gap..." << endl; }
    TH1D *total_pt_outside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_data = (TH1D*) data_file->Get("ak5PF_total_pt_outside_gap");

    total_pt_outside_gap_det_p6_z2->SetTitle("Outside-jet Total p_{T};p_{T}^{outside} [GeV];dN/dp_{T}^{outside} [GeV^{-1}]");
    plot_histogram(total_pt_outside_gap_gen_p6_z2, total_pt_outside_gap_gen_p8_4c, total_pt_outside_gap_det_p6_z2, total_pt_outside_gap_det_p8_4c, total_pt_outside_gap_data, output_path_plots, plot_prefix + "xsec_all_total_pt_outside_gap", "bottom_left", "top_right", "OUTSIDE-JET TAG", detail);


//plot leading_eta_outside_gap distribution
    if (detail) { cout << "Ploting leading_eta_outside_gap..." << endl; }
    TH1D *leading_eta_outside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_eta_outside_gap");

    leading_eta_outside_gap_det_p6_z2->SetTitle("Outside-jet Leading #eta;#eta^{outside};dN/d#eta^{outside}");
    plot_histogram(leading_eta_outside_gap_gen_p6_z2, leading_eta_outside_gap_gen_p8_4c, leading_eta_outside_gap_det_p6_z2, leading_eta_outside_gap_det_p8_4c, leading_eta_outside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_eta_outside_gap", "bottom_middle", "middle", "OUTSIDE-JET TAG", detail);


//plot leading_phi_outside_gap distribution
    if (detail) { cout << "Ploting leading_phi_outside_gap..." << endl; }
    TH1D *leading_phi_outside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_data = (TH1D*) data_file->Get("ak5PF_leading_phi_outside_gap");

    leading_phi_outside_gap_det_p6_z2->SetTitle("Outside-jet Leading #phi;#phi^{outside} [rad];dN/d#phi^{outside}");
    plot_histogram(leading_phi_outside_gap_gen_p6_z2, leading_phi_outside_gap_gen_p8_4c, leading_phi_outside_gap_det_p6_z2, leading_phi_outside_gap_det_p8_4c, leading_phi_outside_gap_data, output_path_plots, plot_prefix + "xsec_all_leading_phi_outside_gap", "bottom_left", "bottom_right", "OUTSIDE-JET TAG", detail);


//plot multiplicity_outside_gap distribution
    if (detail) { cout << "Ploting multiplicity_outside_gap..." << endl; }
    TH1D *multiplicity_outside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_multiplicity_outside_gap");
    TH1D *multiplicity_outside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_multiplicity_outside_gap");
    TH1D *multiplicity_outside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_multiplicity_outside_gap");
    TH1D *multiplicity_outside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_multiplicity_outside_gap");
    TH1D *multiplicity_outside_gap_data = (TH1D*) data_file->Get("ak5PF_multiplicity_outside_gap");

    multiplicity_outside_gap_det_p6_z2->SetTitle("Multiplicity Outside;Multiplicity^{outside};dN/dMultiplicity^{outside}");
    plot_histogram(multiplicity_outside_gap_gen_p6_z2, multiplicity_outside_gap_gen_p8_4c, multiplicity_outside_gap_det_p6_z2, multiplicity_outside_gap_det_p8_4c, multiplicity_outside_gap_data, output_path_plots, plot_prefix + "xsec_all_multiplicity_outside_gap", "bottom_left", "top_right", "OUTSIDE-JET TAG", detail);

//plot leading_pt_outside_gap distribution
    if (detail) { cout << "Ploting delta_eta_outside_gap..." << endl; }
    TH1D *delta_eta_outside_gap_gen_p6_z2 = (TH1D*) p6_z2_gen->Get("ak5Gen_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_gen_p8_4c = (TH1D*) p8_4c_gen->Get("ak5Gen_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_det_p6_z2 = (TH1D*) p6_z2_det->Get("ak5PF_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_det_p8_4c = (TH1D*) p8_4c_det->Get("ak5PF_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_data = (TH1D*) data_file->Get("ak5PF_delta_eta_outside_gap");

    delta_eta_outside_gap_det_p6_z2->SetTitle("#Delta#eta';#Delta#eta^{out};dN/d#Delta#eta^{out}");
    plot_histogram(delta_eta_outside_gap_gen_p6_z2, delta_eta_outside_gap_gen_p8_4c, delta_eta_outside_gap_det_p6_z2, delta_eta_outside_gap_det_p8_4c, delta_eta_outside_gap_data, output_path_plots, plot_prefix + "xsec_all_delta_eta_outside_gap", "bottom_left", "top_right", "OUTSIDE-JET TAG", detail);

}
