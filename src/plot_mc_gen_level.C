// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 24 Oct 2012
//
// plot_mc_gen_level()
// compiles the analysis routines and runs them to change the analysis

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "common_methods.h"

void scale_histogram(TH1D *shape, TH1D *normalization, bool detail)
{
//scale a histogram to the total cross section

    double integral = 0.0;
    double norm_integral = 0.0;
    double scale = 0.0;

    integral = shape->Integral();
    norm_integral = normalization->Integral();
    scale = norm_integral/integral;

    if (detail) { cout << "Scaling Factor = " << scale << endl; }
    shape->Scale(scale);

}

void plot_histogram(TH1D *mc1, TH1D *mc2, TH1D *mc3, TH1D *mc4, string path, string fileout, string legend_position = "top_left", bool detail = false)
{
//plots the model uncertainty control plots

//declaring the canvas
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogy();

//calculate the plooting range
    double min = 0.0, max = 0.0;
    min = mc1->GetMinimum();
    max = mc1->GetMaximum();
    if (min > mc2->GetMinimum()) { min = mc2->GetMinimum(); }
    if (max < mc2->GetMaximum()) { max = mc2->GetMaximum(); }
    if (min > mc4->GetMinimum()) { min = mc4->GetMinimum(); }
    if (max < mc4->GetMaximum()) { max = mc4->GetMaximum(); }
    
    max = 1.3 * max;
    min = 0.7 * min;

//format and ploting the histogram
    mc1->SetMaximum(max);
    mc1->SetMinimum(min);
    mc1->SetLineColor(2);
    mc1->SetLineStyle(1);
    mc1->SetLineWidth(3);
    mc1->Draw("e1");
    mc2->SetLineColor(3);
    mc2->SetLineStyle(2);
    mc2->SetLineWidth(3);
    mc2->Draw("e1 same");
    mc3->SetLineColor(4);
    mc3->SetLineStyle(3);
    mc3->SetLineWidth(3);
    mc3->Draw("e1 same");
    mc4->SetLineColor(5);
    mc4->SetLineStyle(4);
    mc4->SetLineWidth(3);
    mc4->Draw("e1 same");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 4, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1,y1,x2,y2);
    leg00->AddEntry(mc1,"Out of the box","l");
    leg00->AddEntry(mc2,"JetMETTau","l");
    leg00->AddEntry(mc3,"JetMET","l");
    leg00->AddEntry(mc4,"Jet","l");
    leg00->SetFillColor(0);
    leg00->SetLineStyle(1);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();

    print_plots(c1, path, fileout);

}


void plot_mc_gen_level(string path_mc1, string path_mc2, string path_mc3, string path_mc4, string output_path_plots = "../output/control_dist/", bool detail = false, bool disp_output = true)
{
//plots the mc generator level cross-section

//output configurations
    if (detail) { cout << "Plot MC Gen Configuration" << endl; }
    if (detail) { cout << "Input path mc1:  " << path_mc1 << endl; }
    if (detail) { cout << "Input path mc2:  " << path_mc2 << endl; }
    if (detail) { cout << "Input path mc3:  " << path_mc3 << endl; }
    if (detail) { cout << "Input path mc4:  " << path_mc4 << endl; }
    if (detail) { cout << "Output path:     " << output_path_plots << endl; }
    if (detail) { cout << "Detail Level:    " << detail << endl; }
    if (detail) { cout << "Display Results: " << disp_output << endl; }

    if (detail) { cout << "Opening files..." << endl; }

    TFile *file_mc1 = new TFile( path_mc1.c_str() );
    TFile *file_mc2 = new TFile( path_mc2.c_str() );
    TFile *file_mc3 = new TFile( path_mc3.c_str() );
    TFile *file_mc4= new TFile( path_mc4.c_str() );

    if (detail) { cout << "Files Opened sucessfully!" << endl; }

//plot delta phi distribution
    if (detail) { cout << "Ploting Delta Phi..." << endl; }
    TH1D *delta_phi_mc1 = (TH1D*) file_mc1->Get("ak5Gen_delta_phi");
    TH1D *delta_phi_mc2 = (TH1D*) file_mc2->Get("ak5Gen_delta_phi");
    TH1D *delta_phi_mc3 = (TH1D*) file_mc3->Get("ak5Gen_delta_phi");
    TH1D *delta_phi_mc4 = (TH1D*) file_mc4->Get("ak5Gen_delta_phi");

    plot_histogram(delta_phi_mc1, delta_phi_mc2, delta_phi_mc3, delta_phi_mc4, output_path_plots, "xsec_gen_delta_phi","top_left",detail);

/*

//plot delta phi deta1 distribution
    if (detail) { cout << "Ploting Delta Phi Deta1..." << endl; }
    TH1D *delta_phi_deta1_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta1");
    TH1D *delta_phi_deta1_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta1");
    TH1D *delta_phi_deta1_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta1");
    TH1D *delta_phi_deta1_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta1");

    plot_histogram(delta_phi_deta1_p6_z2_1vertex, delta_phi_deta1_p6_z2_allvertex, delta_phi_deta1_p8_4c_1vertex, delta_phi_deta1_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta1","top_left",detail);

//plot delta phi deta2 distribution
    if (detail) { cout << "Ploting Delta Phi Deta2..." << endl; }
    TH1D *delta_phi_deta2_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta2");
    TH1D *delta_phi_deta2_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta2");
    TH1D *delta_phi_deta2_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta2");
    TH1D *delta_phi_deta2_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta2");

    plot_histogram(delta_phi_deta2_p6_z2_1vertex, delta_phi_deta2_p6_z2_allvertex, delta_phi_deta2_p8_4c_1vertex, delta_phi_deta2_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta2","top_left",detail);

//plot delta phi deta3 distribution
    if (detail) { cout << "Ploting Delta Phi Deta3..." << endl; }
    TH1D *delta_phi_deta3_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta3");
    TH1D *delta_phi_deta3_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta3");
    TH1D *delta_phi_deta3_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta3");
    TH1D *delta_phi_deta3_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta3");

    plot_histogram(delta_phi_deta3_p6_z2_1vertex, delta_phi_deta3_p6_z2_allvertex, delta_phi_deta3_p8_4c_1vertex, delta_phi_deta3_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta3","top_left",detail);

//plot delta phi deta4 distribution
    if (detail) { cout << "Ploting Delta Phi Deta3..." << endl; }
    TH1D *delta_phi_deta4_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta4");
    TH1D *delta_phi_deta4_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta4");
    TH1D *delta_phi_deta4_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta4");
    TH1D *delta_phi_deta4_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta4");

    plot_histogram(delta_phi_deta4_p6_z2_1vertex, delta_phi_deta4_p6_z2_allvertex, delta_phi_deta4_p8_4c_1vertex, delta_phi_deta4_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta4","top_left",detail);

//plot delta phi gap distribution
    if (detail) { cout << "Ploting Delta Phi Gap..." << endl; }
    TH1D *delta_phi_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_gap");
    TH1D *delta_phi_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_gap");
    TH1D *delta_phi_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_gap");
    TH1D *delta_phi_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_gap");

    plot_histogram(delta_phi_gap_p6_z2_1vertex, delta_phi_gap_p6_z2_allvertex, delta_phi_gap_p8_4c_1vertex, delta_phi_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_gap","top_left",detail);

//plot delta phi deta1 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta1 Gap..." << endl; }
    TH1D *delta_phi_deta1_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta1_gap");
    TH1D *delta_phi_deta1_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta1_gap");

    plot_histogram(delta_phi_deta1_gap_p6_z2_1vertex, delta_phi_deta1_gap_p6_z2_allvertex, delta_phi_deta1_gap_p8_4c_1vertex, delta_phi_deta1_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta1_gap","top_left",detail);

//plot delta phi deta2 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta2 Gap..." << endl; }
    TH1D *delta_phi_deta2_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta2_gap");
    TH1D *delta_phi_deta2_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta2_gap");

    plot_histogram(delta_phi_deta2_gap_p6_z2_1vertex, delta_phi_deta2_gap_p6_z2_allvertex, delta_phi_deta2_gap_p8_4c_1vertex, delta_phi_deta2_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta2_gap","top_left",detail);

//plot delta phi deta3 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta3 Gap..." << endl; }
    TH1D *delta_phi_deta3_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta3_gap");
    TH1D *delta_phi_deta3_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta3_gap");

    plot_histogram(delta_phi_deta3_gap_p6_z2_1vertex, delta_phi_deta3_gap_p6_z2_allvertex, delta_phi_deta3_gap_p8_4c_1vertex, delta_phi_deta3_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta3_gap","top_left",detail);

//plot delta phi deta4 gap distribution
    if (detail) { cout << "Ploting Delta Phi Deta4 Gap..." << endl; }
    TH1D *delta_phi_deta4_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta4_gap");
    TH1D *delta_phi_deta4_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta4_gap");

    plot_histogram(delta_phi_deta4_gap_p6_z2_1vertex, delta_phi_deta4_gap_p6_z2_allvertex, delta_phi_deta4_gap_p8_4c_1vertex, delta_phi_deta4_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta4_gap","top_left",detail);

//plot delta phi nogap distribution
    if (detail) { cout << "Ploting Delta Phi noGap..." << endl; }
    TH1D *delta_phi_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_nogap");
    TH1D *delta_phi_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_nogap");
    TH1D *delta_phi_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_nogap");
    TH1D *delta_phi_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_nogap");

    plot_histogram(delta_phi_nogap_p6_z2_1vertex, delta_phi_nogap_p6_z2_allvertex, delta_phi_nogap_p8_4c_1vertex, delta_phi_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_nogap","top_left",detail);

//plot delta phi deta1 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta1 noGap..." << endl; }
    TH1D *delta_phi_deta1_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta1_nogap");
    TH1D *delta_phi_deta1_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta1_nogap");

    plot_histogram(delta_phi_deta1_nogap_p6_z2_1vertex, delta_phi_deta1_nogap_p6_z2_allvertex, delta_phi_deta1_nogap_p8_4c_1vertex, delta_phi_deta1_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta1_nogap","top_left",detail);

//plot delta phi deta2 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta2 noGap..." << endl; }
    TH1D *delta_phi_deta2_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta2_nogap");
    TH1D *delta_phi_deta2_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta2_nogap");

    plot_histogram(delta_phi_deta2_nogap_p6_z2_1vertex, delta_phi_deta2_nogap_p6_z2_allvertex, delta_phi_deta2_nogap_p8_4c_1vertex, delta_phi_deta2_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta2_nogap","top_left",detail);

//plot delta phi deta3 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta3 noGap..." << endl; }
    TH1D *delta_phi_deta3_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta3_nogap");
    TH1D *delta_phi_deta3_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta3_nogap");

    plot_histogram(delta_phi_deta3_nogap_p6_z2_1vertex, delta_phi_deta3_nogap_p6_z2_allvertex, delta_phi_deta3_nogap_p8_4c_1vertex, delta_phi_deta3_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta3_nogap","top_left",detail);

//plot delta phi deta4 nogap distribution
    if (detail) { cout << "Ploting Delta Phi Deta4 noGap..." << endl; }
    TH1D *delta_phi_deta4_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_phi_deta4_nogap");
    TH1D *delta_phi_deta4_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_phi_deta4_nogap");

    plot_histogram(delta_phi_deta4_nogap_p6_z2_1vertex, delta_phi_deta4_nogap_p6_z2_allvertex, delta_phi_deta4_nogap_p8_4c_1vertex, deltusing namespace std;

#include "common_methods.h"a_phi_deta4_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_phi_deta4_nogap","top_left",detail);

//plot delta eta distribution
    if (detail) { cout << "Ploting Delta Eta..." << endl; }
    TH1D *delta_eta_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_eta");
    TH1D *delta_eta_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_eta");
    TH1D *delta_eta_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_eta");
    TH1D *delta_eta_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_eta");

    plot_histogram(delta_eta_p6_z2_1vertex, delta_eta_p6_z2_allvertex, delta_eta_p8_4c_1vertex, delta_eta_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_eta","bottom_left",detail);

//plot delta eta gap distribution
    if (detail) { cout << "Ploting Delta Eta Gap..." << endl; }
    TH1D *delta_eta_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_eta_gap");
    TH1D *delta_eta_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_eta_gap");
    TH1D *delta_eta_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_eta_gap");
    TH1D *delta_eta_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_eta_gap");

    plot_histogram(delta_eta_gap_p6_z2_1vertex, delta_eta_gap_p6_z2_allvertex, delta_eta_gap_p8_4c_1vertex, delta_eta_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_eta_gap","bottom_left",detail);

//plot delta eta nogap distribution
    if (detail) { cout << "Ploting Delta Eta noGap..." << endl; }
    TH1D *delta_eta_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_eta_nogap");
    TH1D *delta_eta_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_eta_nogap");
    TH1D *delta_eta_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_eta_nogap");
    TH1D *delta_eta_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_eta_nogap");

    plot_histogram(delta_eta_nogap_p6_z2_1vertex, delta_eta_nogap_p6_z2_allvertex, delta_eta_nogap_p8_4c_1vertex, delta_eta_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_eta_nogap","bottom_left",detail);

//plot total pT inside gap distribution
    if (detail) { cout << "Ploting Total pT Inside Gap..." << endl; }
    TH1D *total_pt_inside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_total_pt_inside_gap");
    TH1D *total_pt_inside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_total_pt_inside_gap");

    plot_histogram(total_pt_inside_gap_p6_z2_1vertex, total_pt_inside_gap_p6_z2_allvertex, total_pt_inside_gap_p8_4c_1vertex, total_pt_inside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_total_pt_inside_gap","bottom_left",detail);

//plot leading pT inside gap distribution
    if (detail) { cout << "Ploting Leading pT Inside Gap..." << endl; }
    TH1D *leading_pt_inside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_pt_inside_gap");
    TH1D *leading_pt_inside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_pt_inside_gap");

    plot_histogram(leading_pt_inside_gap_p6_z2_1vertex, leading_pt_inside_gap_p6_z2_allvertex, leading_pt_inside_gap_p8_4c_1vertex, leading_pt_inside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_pt_inside_gap","bottom_left",detail);

//plot leading eta inside gap distribution
    if (detail) { cout << "Ploting Leading Eta Inside Gap..." << endl; }
    TH1D *leading_eta_inside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_eta_inside_gap");
    TH1D *leading_eta_inside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_eta_inside_gap");

    plot_histogram(leading_eta_inside_gap_p6_z2_1vertex, leading_eta_inside_gap_p6_z2_allvertex, leading_eta_inside_gap_p8_4c_1vertex, leading_eta_inside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_eta_inside_gap","bottom_middle",detail);

//plot leading phi inside gap distribution
    if (detail) { cout << "Ploting Leading Phi Inside Gap..." << endl; }
    TH1D *leading_phi_inside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_phi_inside_gap");
    TH1D *leading_phi_inside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_phi_inside_gap");

    plot_histogram(leading_phi_inside_gap_p6_z2_1vertex, leading_phi_inside_gap_p6_z2_allvertex, leading_phi_inside_gap_p8_4c_1vertex, leading_phi_inside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_phi_inside_gap","bottom_middle",detail);

//plot total pT outside gap distribution
    if (detail) { cout << "Ploting Total pT Outside Gap..." << endl; }
    TH1D *total_pt_outside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_total_pt_outside_gap");
    TH1D *total_pt_outside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_total_pt_outside_gap");

    plot_histogram(total_pt_outside_gap_p6_z2_1vertex, total_pt_outside_gap_p6_z2_allvertex, total_pt_outside_gap_p8_4c_1vertex, total_pt_outside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_total_pt_outside_gap","bottom_left",detail);

//plot leading pT outside gap distribution
    if (detail) { cout << "Ploting Leading pT Outside Gap..." << endl; }
    TH1D *leading_pt_outside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_pt_outside_gap");
    TH1D *leading_pt_outside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_pt_outside_gap");

    plot_histogram(leading_pt_outside_gap_p6_z2_1vertex, leading_pt_outside_gap_p6_z2_allvertex, leading_pt_outside_gap_p8_4c_1vertex, leading_pt_outside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_pt_outside_gap","bottom_left",detail);

//plot delta eta outside gap distribution
    if (detail) { cout << "Ploting Delta Eta Outside Gap..." << endl; }
    TH1D *delta_eta_outside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_delta_eta_outside_gap");
    TH1D *delta_eta_outside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_delta_eta_outside_gap");

    plot_histogram(delta_eta_outside_gap_p6_z2_1vertex, delta_eta_outside_gap_p6_z2_allvertex, delta_eta_outside_gap_p8_4c_1vertex, delta_eta_outside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_delta_eta_outside_gap","bottom_left",detail);

//plot leading eta outside gap distribution
    if (detail) { cout << "Ploting Leading Eta Outside Gap..." << endl; }
    TH1D *leading_eta_outside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_eta_outside_gap");
    TH1D *leading_eta_outside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_eta_outside_gap");

    plot_histogram(leading_eta_outside_gap_p6_z2_1vertex, leading_eta_outside_gap_p6_z2_allvertex, leading_eta_outside_gap_p8_4c_1vertex, leading_eta_outside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_eta_outside_gap","bottom_middle",detail);

//plot leading phi outside gap distribution
    if (detail) { cout << "Ploting Leading Phi Outside Gap..." << endl; }
    TH1D *leading_phi_outside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_phi_outside_gap");
    TH1D *leading_phi_outside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_phi_outside_gap");

    plot_histogram(leading_phi_outside_gap_p6_z2_1vertex, leading_phi_outside_gap_p6_z2_allvertex, leading_phi_outside_gap_p8_4c_1vertex, leading_phi_outside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_phi_outside_gap","bottom_middle",detail);

//plot leading eta* inside gap distribution
    if (detail) { cout << "Ploting Leading Eta* Inside Gap..." << endl; }
    TH1D *leading_eta_star_inside_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_eta_star_inside_gap");
    TH1D *leading_eta_star_inside_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_eta_star_inside_gap");

    plot_histogram(leading_eta_star_inside_gap_p6_z2_1vertex, leading_eta_star_inside_gap_p6_z2_allvertex, leading_eta_star_inside_gap_p8_4c_1vertex, leading_eta_star_inside_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_eta_star_inside_gap","bottom_middle",detail);

//plot leading central pT distribution
    if (detail) { cout << "Ploting Leading Central pT..." << endl; }
    TH1D *leading_central_pt_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_central_pt");
    TH1D *leading_central_pt_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_central_pt");
    TH1D *leading_central_pt_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_central_pt");
    TH1D *leading_central_pt_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_central_pt");

    plot_histogram(leading_central_pt_p6_z2_1vertex, leading_central_pt_p6_z2_allvertex, leading_central_pt_p8_4c_1vertex, leading_central_pt_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_central_pt","bottom_left",detail);

//plot leading central pT gap distribution
    if (detail) { cout << "Ploting Leading Central pT Gap..." << endl; }
    TH1D *leading_central_pt_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_central_pt_gap");
    TH1D *leading_central_pt_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_central_pt_gap");
    TH1D *leading_central_pt_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_central_pt_gap");
    TH1D *leading_central_pt_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_central_pt_gap");

    plot_histogram(leading_central_pt_gap_p6_z2_1vertex, leading_central_pt_gap_p6_z2_allvertex, leading_central_pt_gap_p8_4c_1vertex, leading_central_pt_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_central_pt_gap","bottom_left",detail);

//plot leading central pT nogap distribution
    if (detail) { cout << "Ploting Leading Central pT noGap..." << endl; }
    TH1D *leading_central_pt_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_central_pt_nogap");
    TH1D *leading_central_pt_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_central_pt_nogap");
    TH1D *leading_central_pt_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_central_pt_nogap");
    TH1D *leading_central_pt_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_central_pt_nogap");

    plot_histogram(leading_central_pt_nogap_p6_z2_1vertex, leading_central_pt_nogap_p6_z2_allvertex, leading_central_pt_nogap_p8_4c_1vertex, leading_central_pt_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_central_pt_nogap","bottom_left",detail);

//plot leading central eta distribution
    if (detail) { cout << "Ploting Leading Central Eta..." << endl; }
    TH1D *leading_central_eta_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_central_eta");
    TH1D *leading_central_eta_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_central_eta");
    TH1D *leading_central_eta_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_central_eta");
    TH1D *leading_central_eta_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_central_eta");

    plot_histogram(leading_central_eta_p6_z2_1vertex, leading_central_eta_p6_z2_allvertex, leading_central_eta_p8_4c_1vertex, leading_central_eta_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_central_eta","bottom_middle",detail);

//plot leading central phi distribution
    if (detail) { cout << "Ploting Leading Central Phi..." << endl; }
    TH1D *leading_central_phi_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_central_phi");
    TH1D *leading_central_phi_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_central_phi");
    TH1D *leading_central_phi_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_central_phi");
    TH1D *leading_central_phi_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_central_phi");

    plot_histogram(leading_central_phi_p6_z2_1vertex, leading_central_phi_p6_z2_allvertex, leading_central_phi_p8_4c_1vertex, leading_central_phi_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_central_phi","bottom_middle",detail);

//plot leading forward pT distribution
    if (detail) { cout << "Ploting Leading Forward pT..." << endl; }
    TH1D *leading_forward_pt_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_forward_pt");
    TH1D *leading_forward_pt_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_forward_pt");
    TH1D *leading_forward_pt_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_forward_pt");
    TH1D *leading_forward_pt_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_forward_pt");

    plot_histogram(leading_forward_pt_p6_z2_1vertex, leading_forward_pt_p6_z2_allvertex, leading_forward_pt_p8_4c_1vertex, leading_forward_pt_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_forward_pt","bottom_left",detail);

//plot leading forward pT gap distribution
    if (detail) { cout << "Ploting Leading Forward pT Gap..." << endl; }
    TH1D *leading_forward_pt_gap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_forward_pt_gap");
    TH1D *leading_forward_pt_gap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_forward_pt_gap");
    TH1D *leading_forward_pt_gap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_forward_pt_gap");
    TH1D *leading_forward_pt_gap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_forward_pt_gap");

    plot_histogram(leading_forward_pt_gap_p6_z2_1vertex, leading_forward_pt_gap_p6_z2_allvertex, leading_forward_pt_gap_p8_4c_1vertex, leading_forward_pt_gap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_forward_pt_gap","bottom_left",detail);

//plot leading forward pT nogap distribution
    if (detail) { cout << "Ploting Leading Forward pT noGap..." << endl; }
    TH1D *leading_forward_pt_nogap_p6_z2_1vertex = (TH1D*) p6_z2_file_1vertex->Get("ak5Gen_leading_forward_pt_nogap");
    TH1D *leading_forward_pt_nogap_p6_z2_allvertex = (TH1D*) p6_z2_file_allvertex->Get("ak5Gen_leading_forward_pt_nogap");
    TH1D *leading_forward_pt_nogap_p8_4c_1vertex = (TH1D*) p8_4c_file_1vertex->Get("ak5Gen_leading_forward_pt_nogap");
    TH1D *leading_forward_pt_nogap_p8_4c_allvertex = (TH1D*) p8_4c_file_allvertex->Get("ak5Gen_leading_forward_pt_nogap");

    plot_histogram(leading_forward_pt_nogap_p6_z2_1vertex, leading_forward_pt_nogap_p6_z2_allvertex, leading_forward_pt_nogap_p8_4c_1vertex, leading_forward_pt_nogap_p8_4c_allvertex, output_path_plots, "xsec_gen_leading_forward_pt_nogap","bottom_left",detail);
*/

}
