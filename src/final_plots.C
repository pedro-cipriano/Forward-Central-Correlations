// DESY, CMS
// Last Update: 24 Mar 2013
//
// final_plots()

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TPad.h>
#include <TString.h>
#include <TPaveText.h>

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "common_methods.h"

using namespace std;

void create_final_plot(TH1D *data, TH1D *data_unc, TH1 *mc1, TString label1, TH1 *mc2, TString label2, TH1 *mc3, TString label3, TH1 *mc4, TString label4, TH1 *mc5, TString label5, TString scenario, string output_path, string mode, string fileout, string legend_position, string label_position, TString extra_label = "", string extra_label_position = "", bool detail = false)
{
//creates the final plots

if (detail) { cout << "Positions : " << legend_position << " and " << label_position << endl; }

int n_items = 0;
if (mc1 != 0) { n_items = 1; }
if (mc2 != 0 and n_items == 1) { n_items = 2; }
if (mc3 != 0 and n_items == 2) { n_items = 3; }
if (mc4 != 0 and n_items == 3) { n_items = 4; }
if (mc5 != 0 and n_items == 4) { n_items = 5; }

//declaring the canvas

    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,800,600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    if (mode == "xsec") { gPad->SetLogy(1); }
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.07);
    gPad->SetFrameBorderMode(0);

    data_unc->GetXaxis()->SetLabelSize(0.05);
    data_unc->GetYaxis()->SetLabelSize(0.05);
    data_unc->GetXaxis()->SetTitleSize(0.05);
    data_unc->GetYaxis()->SetTitleSize(0.05);

//calculate the plooting range
    	double min = 0.0;
    	double max = data->GetMaximum();
    if (mode == "xsec")
	{
    	if (data->GetMinimum() == 0.0)
    	{
    	min = get_non0_minimum(data,detail);
   	}
    	else
    	{
    	min = data->GetMinimum();
    	}    
    
    	set_histogram_min_max(data_unc, min, max, detail);
    	if (n_items >= 1) { set_histogram_min_max(mc1, min, max, detail); }
    	if (n_items >= 2) { set_histogram_min_max(mc2, min, max, detail); }
    	if (n_items >= 3) { set_histogram_min_max(mc3, min, max, detail); }
    	if (n_items >= 4) { set_histogram_min_max(mc4, min, max, detail); }
    	if (n_items == 5) { set_histogram_min_max(mc5, min, max, detail); }
    
    	max = 1.3 * max;
    	min = 0.7 * min;
	if (detail) { cout << "Min = " << min << " Max = " << max << endl; }
	}
    if (mode == "ratio") { max = 2.8; } 

//format and ploting the histogram
    data_unc->SetMinimum(min);
    data_unc->SetMaximum(max);
    data_unc->SetFillColor(5);
    data_unc->Draw("e2");
    data->SetMinimum(min);
    data->SetMaximum(max);
    data->SetLineStyle(1);
    data->SetLineWidth(4);
    data->SetLineColor(1);
    data->SetMarkerColor(1);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(2);
    data->DrawCopy("e1 same");
	if (n_items >= 1)
	{
    	mc1->SetLineStyle(1);
    	mc1->SetLineWidth(4);
    	mc1->SetLineColor(2);
    	//mc1->SetMarkerColor(2);
    	//mc1->SetMarkerStyle(20);
    	//mc1->SetMarkerSize(3);
    	mc1->DrawCopy("hsame");
    	mc1->DrawCopy("esame");
	}
	if (n_items >= 2)
	{
    	mc2->SetLineStyle(2);
    	mc2->SetLineWidth(4);
    	mc2->SetLineColor(4);
    	//mc2->SetMarkerColor(4);
    	//mc2->SetMarkerStyle(21);
    	//mc2->SetMarkerSize(3);
    	mc2->DrawCopy("hsame");
    	mc2->DrawCopy("esame");
	}
	if (n_items >= 3)
	{
    	mc3->SetLineStyle(10);
    	mc3->SetLineWidth(4);
    	mc3->SetLineColor(6);
    	//mc3->SetMarkerColor(6);
    	//mc3->SetMarkerStyle(22);
    	//mc3->SetMarkerSize(3);
    	mc3->DrawCopy("hsame");
    	mc3->DrawCopy("esame");
	}
	if (n_items >= 4)
	{
    	mc4->SetLineStyle(3);
    	mc4->SetLineWidth(4);
    	mc4->SetLineColor(9);
    	//mc4->SetMarkerColor(7);
    	//mc4->SetMarkerStyle(23);
    	//mc4->SetMarkerSize(3);
    	mc4->DrawCopy("hsame");
    	mc4->DrawCopy("esame");
	}
	if (n_items >= 5)
	{
    	mc5->SetLineStyle(9);
    	mc5->SetLineWidth(4);
    	mc5->SetLineColor(11);
    	//mc5->SetMarkerColor(8);
    	//mc5->SetMarkerStyle(24);
    	//mc5->SetMarkerSize(3);
    	mc5->DrawCopy("hsame");
    	mc5->DrawCopy("esame");
	}
 
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, n_items + 1, x1, y1, x2, y2);

//sets and draw the legend
    TLegend *leg00 = new TLegend(x1, y1-0.06, x2, y2-0.06);
    leg00->AddEntry(data,"Data","lep");
    if (n_items >= 1) { leg00->AddEntry(mc1,label1,"l"); }
    if (n_items >= 2) {	leg00->AddEntry(mc2,label2,"l"); }
    if (n_items >= 3) {	leg00->AddEntry(mc3,label3,"l"); }
    if (n_items >= 4) { leg00->AddEntry(mc4,label4,"l"); }
    if (n_items == 5) { leg00->AddEntry(mc5,label5,"l"); }
    leg00->SetFillColor(0);
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
  //  cms->AddText("L_{int} = 34.5 pb^{-1}, Anti-k_{T} (R = 0.5)");
    cms->AddText("L_{int} = 3.2 pb^{-1}, Anti-k_{T} (R = 0.5)");
    cms->AddText("p_{T}^{jet} > 35 GeV and |#eta| < 2.8");
    cms->AddText("p_{T}^{jet} > 35 GeV and 3.2 < |#eta| < 4.7");
    cms->Draw();

    if (extra_label_position == "top_right")    { x1 = 0.79; y1 = 0.77; x2 = 0.98; y2 = 0.92; }
    if (extra_label_position == "bottom_right") { x1 = 0.75; y1 = 0.13; x2 = 0.98; y2 = 0.18; }
    if (extra_label_position == "top_left")     { x1 = 0.13; y1 = 0.77; x2 = 0.49; y2 = 0.92; }
    if (extra_label_position == "middle_left")  { x1 = 0.13; y1 = 0.50; x2 = 0.30; y2 = 0.55; }
    if (extra_label_position == "bottom_left")  { x1 = 0.13; y1 = 0.13; x2 = 0.47; y2 = 0.18; }
    if (extra_label_position == "middle")       { x1 = 0.50; y1 = 0.48; x2 = 0.60; y2 = 0.53; }

    TPaveText *extra = new TPaveText(x1,y1,x2,y2,"NDC"); // NDC sets coords
    extra->SetTextSize(0.05);
    extra->SetBorderSize(0); 
    extra->SetTextAlign(12);
    extra->SetLineWidth(1);
    extra->SetLineColor(0);
    extra->SetFillColor(0);
    extra->SetFillStyle(1001);
    extra->AddText(extra_label);
    if (extra_label != "") { extra->Draw(); }

 print_plots(c1, output_path, mode + fileout);

}


void create_scenario_plot(TH1D *data_A, TString label_A, TH1D *data_B, TString label_B, TH1D *data_C, TString label_C, TH1D *data_unc_A, TH1D *data_unc_B, TH1D *data_unc_C, TH1 *mc1_A, TH1 *mc1_B, TH1 *mc1_C, TString label1, TH1 *mc2_A, TH1 *mc2_B, TH1 *mc2_C, TString label2, TH1 *mc3_A, TH1 *mc3_B, TH1 *mc3_C, TString label3, TH1 *mc4_A, TH1 *mc4_B, TH1 *mc4_C, TString label4, TH1 *mc5_A, TH1 *mc5_B, TH1 *mc5_C, TString label5, string output_path, string fileout, string legend1_position, string legend2_position, string label_position, bool detail = false)
{
//plots the model uncertainty control plots

if (detail) { cout << "Positions : " << legend1_position << " and " << label_position << endl; }

int n_items = 0;
if (mc1_A != 0) { n_items = 1; }
if (mc2_A != 0 and n_items == 1) { n_items = 2; }
if (mc3_A != 0 and n_items == 2) { n_items = 3; }
if (mc4_A != 0 and n_items == 3) { n_items = 4; }
if (mc5_A != 0 and n_items == 4) { n_items = 5; }

//declaring the canvas
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,800,600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetLogy();
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.07);
    gPad->SetFrameBorderMode(0);

//scaling factors
	double factor_B = 100;
	double factor_C = 10000;

//scaling histograms
	data_B->Scale(factor_B);
	data_unc_B->Scale(factor_B);
	if (n_items >= 1) { mc1_B->Scale(factor_B); }
	if (n_items >= 2) { mc2_B->Scale(factor_B); }
	if (n_items >= 3) { mc3_B->Scale(factor_B); }
	if (n_items >= 4) { mc4_B->Scale(factor_B); }
	if (n_items == 5) { mc5_B->Scale(factor_B); }

	data_C->Scale(factor_C);
	data_unc_C->Scale(factor_C);
	if (n_items >= 1) { mc1_C->Scale(factor_C); }
	if (n_items >= 2) { mc2_C->Scale(factor_C); }
	if (n_items >= 3) { mc3_C->Scale(factor_C); }
	if (n_items >= 4) { mc4_C->Scale(factor_C); }
	if (n_items == 5) { mc5_C->Scale(factor_C); }

//calculate the plooting range
    	double min = 0.0;
    	double max = data_A->GetMaximum();
    	if (data_A->GetMinimum() == 0.0)
    	{
    	min = get_non0_minimum(data_A,detail);
   	}
    	else
    	{
    	min = data_A->GetMinimum();
    	}    
    
    	set_histogram_min_max(data_unc_A, min, max, detail);
    	if (n_items >= 1) { set_histogram_min_max(mc1_A, min, max, detail); }
    	if (n_items >= 2) { set_histogram_min_max(mc2_A, min, max, detail); }
    	if (n_items >= 3) { set_histogram_min_max(mc3_A, min, max, detail); }
    	if (n_items >= 4) { set_histogram_min_max(mc4_A, min, max, detail); }
    	if (n_items == 5) { set_histogram_min_max(mc5_A, min, max, detail); }
    
    	set_histogram_min_max(data_unc_C, min, max, detail);
    	if (n_items >= 1) { set_histogram_min_max(mc1_C, min, max, detail); }
    	if (n_items >= 2) { set_histogram_min_max(mc2_C, min, max, detail); }
    	if (n_items >= 3) { set_histogram_min_max(mc3_C, min, max, detail); }
    	if (n_items >= 4) { set_histogram_min_max(mc4_C, min, max, detail); }
    	if (n_items == 5) { set_histogram_min_max(mc5_C, min, max, detail); }

    	max = 1000 * max;
    	min = 0.7 * min;

//format and ploting the histogram
    data_unc_A->SetMinimum(min);
    data_unc_A->SetMaximum(max);
    data_unc_A->SetFillColor(5);
    data_unc_A->DrawCopy("e2");
    data_unc_B->SetMinimum(min);
    data_unc_B->SetMaximum(max);
    data_unc_B->SetFillColor(5);
    data_unc_B->DrawCopy("e2 same");
    data_unc_C->SetFillColor(5);
    data_unc_C->DrawCopy("e2 same");
    data_A->SetLineStyle(1);
    data_A->SetLineWidth(3);
    data_A->SetLineColor(1);
    data_A->SetMarkerColor(1);
    data_A->SetMarkerStyle(20);
    data_A->SetMarkerSize(2);
    data_A->DrawCopy("e1 same");
    data_B->SetLineStyle(1);
    data_B->SetLineWidth(3);
    data_B->SetLineColor(1);
    data_B->SetMarkerColor(2);
    data_B->SetMarkerStyle(21);
    data_B->SetMarkerSize(2);
    data_B->DrawCopy("e1 same");
    data_C->SetLineStyle(1);
    data_C->SetLineWidth(3);
    data_C->SetLineColor(1);
    data_C->SetMarkerColor(3);
    data_C->SetMarkerStyle(22);
    data_C->SetMarkerSize(2);
    data_C->DrawCopy("e1 same");
	if (n_items >= 1)
	{
    	mc1_A->SetLineStyle(1);
    	mc1_A->SetLineWidth(3);
    	mc1_A->SetLineColor(2);
    	//mc1_A->SetMarkerColor(2);
    	//mc1_A->SetMarkerStyle(20);
    	//mc1_A->SetMarkerSize(3);
    	mc1_A->DrawCopy("hsame");
	//mc1_A->DrawCopy("esame");
    	mc1_B->SetLineStyle(1);
    	mc1_B->SetLineWidth(3);
    	mc1_B->SetLineColor(2);
    	//mc1_B->SetMarkerColor(2);
    	//mc1_B->SetMarkerStyle(20);
    	//mc1_B->SetMarkerSize(3);
    	mc1_B->DrawCopy("hsame");
	//mc1_B->DrawCopy("esame");
    	mc1_C->SetLineStyle(1);
    	mc1_C->SetLineWidth(3);
    	mc1_C->SetLineColor(2);
    	//mc1_C->SetMarkerColor(2);
    	//mc1_C->SetMarkerStyle(20);
    	//mc1_C->SetMarkerSize(3);
    	mc1_C->DrawCopy("hsame");
	//mc1_C->DrawCopy("esame");
	}
	if (n_items >= 2)
	{
    	mc2_A->SetLineStyle(2);
    	mc2_A->SetLineWidth(3);
    	mc2_A->SetLineColor(4);
    	//mc2_A->SetMarkerColor(4);
    	//mc2_A->SetMarkerStyle(21);
    	//mc2_A->SetMarkerSize(3);
    	mc2_A->DrawCopy("hsame");
	mc2_A->DrawCopy("esame");
    	mc2_B->SetLineStyle(2);
    	mc2_B->SetLineWidth(3);
    	mc2_B->SetLineColor(4);
    	//mc2_B->SetMarkerColor(4);
    	//mc2_B->SetMarkerStyle(21);
    	//mc2_B->SetMarkerSize(3);
    	mc2_B->DrawCopy("hsame");
	mc2_B->DrawCopy("esame");
    	mc2_C->SetLineStyle(2);
    	mc2_C->SetLineWidth(3);
    	mc2_C->SetLineColor(4);
    	//mc2_C->SetMarkerColor(4);
    	//mc2_C->SetMarkerStyle(21);
    	//mc2_C->SetMarkerSize(3);
    	mc2_C->DrawCopy("hsame");
	mc2_C->DrawCopy("esame");
	}
	if (n_items >= 3)
	{
    	mc3_A->SetLineStyle(10);
    	mc3_A->SetLineWidth(3);
    	mc3_A->SetLineColor(6);
    	//mc3_A->SetMarkerColor(6);
    	//mc3_A->SetMarkerStyle(22);
    	//mc3_A->SetMarkerSize(3);
    	mc3_A->DrawCopy("hsame");
	mc3_A->DrawCopy("esame");
    	mc3_B->SetLineStyle(10);
    	mc3_B->SetLineWidth(3);
    	mc3_B->SetLineColor(6);
    	//mc3_B->SetMarkerColor(6);
    	//mc3_B->SetMarkerStyle(22);
    	//mc3_B->SetMarkerSize(3);
    	mc3_B->DrawCopy("hsame");
	mc3_B->DrawCopy("esame");
    	mc3_C->SetLineStyle(10);
    	mc3_C->SetLineWidth(3);
    	mc3_C->SetLineColor(6);
    	//mc3_C->SetMarkerColor(6);
    	//mc3_C->SetMarkerStyle(22);
    	//mc3_C->SetMarkerSize(3);
    	mc3_C->DrawCopy("hsame");
	mc3_C->DrawCopy("esame");
	}
	if (n_items >= 4)
	{
    	mc4_A->SetLineStyle(3);
    	mc4_A->SetLineWidth(3);
    	mc4_A->SetLineColor(9);
    	//mc4_A->SetMarkerColor(7);
    	//mc4_A->SetMarkerStyle(23);
    	//mc4_A->SetMarkerSize(3);
    	mc4_A->DrawCopy("hsame");
	mc4_A->DrawCopy("esame");
    	mc4_B->SetLineStyle(3);
    	mc4_B->SetLineWidth(3);
    	mc4_B->SetLineColor(9);
    	//mc4_B->SetMarkerColor(7);
    	//mc4_B->SetMarkerStyle(23);
    	//mc4_B->SetMarkerSize(3);
    	mc4_B->DrawCopy("hsame");
	mc4_B->DrawCopy("esame");
    	mc4_C->SetLineStyle(3);
    	mc4_C->SetLineWidth(3);
    	mc4_C->SetLineColor(9);
    	//mc4_C->SetMarkerColor(7);
    	//mc4_C->SetMarkerStyle(23);
    	//mc4_C->SetMarkerSize(3);
    	mc4_C->DrawCopy("hsame");
	mc4_C->DrawCopy("esame");
	}
	if (n_items >= 5)
	{
    	mc5_A->SetLineStyle(10);
    	mc5_A->SetLineWidth(3);
    	mc5_A->SetLineColor(8);
    	//mc5_A->SetMarkerColor(8);
    	//mc5_A->SetMarkerStyle(24);
    	//mc5_A->SetMarkerSize(3);
    	mc5_A->DrawCopy("hsame");
	mc5_A->DrawCopy("esame");
    	mc5_B->SetLineStyle(10);
    	mc5_B->SetLineWidth(3);
    	mc5_B->SetLineColor(8);
    	//mc5_B->SetMarkerColor(8);
    	//mc5_B->SetMarkerStyle(24);
    	//mc5_B->SetMarkerSize(3);
    	mc5_B->DrawCopy("hsame");
	mc5_B->DrawCopy("esame");
    	mc5_C->SetLineStyle(10);
    	mc5_C->SetLineWidth(3);
    	mc5_C->SetLineColor(8);
    	//mc5_C->SetMarkerColor(8);
    	//mc5_C->SetMarkerStyle(24);
    	//mc5_C->SetMarkerSize(3);
    	mc5_C->DrawCopy("hsame");
	mc5_C->DrawCopy("esame");
	}
    
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

//sets and draw the legend
    set_legend_position(legend1_position, n_items, x1, y1, x2, y2);
    TLegend *leg00 = new TLegend(x1, y1-0.06, x2-0.10, y2-0.06);
    if (n_items >= 1) { leg00->AddEntry(mc1_A,label1,"lp"); }
    if (n_items >= 2) {	leg00->AddEntry(mc2_A,label2,"lp"); }
    if (n_items >= 3) {	leg00->AddEntry(mc3_A,label3,"lp"); }
    if (n_items >= 4) { leg00->AddEntry(mc4_A,label4,"lp"); }
    if (n_items == 5) { leg00->AddEntry(mc5_A,label5,"lp"); }
    leg00->SetFillColor(0);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();

    set_legend_position(legend2_position, 3, x1, y1, x2, y2);
    TLegend *leg01 = new TLegend(x1, y1-0.06, x2, y2-0.06);
    leg01->AddEntry(data_C,label_C+"(x10000)","lep");
    leg01->AddEntry(data_B,label_B+"(x100)","lep");
    leg01->AddEntry(data_A,label_A+"(x1)","lep");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();
    
    if (label_position == "top_right")    { x1 = 0.63; y1 = 0.82; x2 = 0.98; y2 = 0.92; }
    if (label_position == "top_left")     { x1 = 0.13; y1 = 0.82; x2 = 0.47; y2 = 0.92; }
    if (label_position == "bottom_left")  { x1 = 0.13; y1 = 0.13; x2 = 0.47; y2 = 0.23; }
    if (label_position == "bottom_right") { x1 = 0.65; y1 = 0.13; x2 = 0.98; y2 = 0.23; }
    if (label_position == "middle")       { x1 = 0.50; y1 = 0.40; x2 = 0.60; y2 = 0.50; }

    TPaveText *cms_topleft = new TPaveText(0.01,0.99,0.35,0.94,"NDC"); // NDC sets coords
    cms_topleft->SetTextSize(0.03);
    cms_topleft->SetBorderSize(0); 
    cms_topleft->SetTextAlign(22);
    cms_topleft->SetLineWidth(1);
    cms_topleft->SetLineColor(0);
    cms_topleft->SetFillColor(0);
    cms_topleft->SetFillStyle(1001);
    cms_topleft->AddText("CMS Preliminary, #sqrt{s} = 7 TeV");
    cms_topleft->Draw();

    TPaveText *cms_toplmiddle = new TPaveText(0.35,0.99,0.49,0.94,"NDC"); // NDC sets coords
    cms_toplmiddle->SetTextSize(0.03);
    cms_toplmiddle->SetBorderSize(0); 
    cms_toplmiddle->SetTextAlign(22);
    cms_toplmiddle->SetLineWidth(1);
    cms_toplmiddle->SetLineColor(0);
    cms_toplmiddle->SetFillColor(0);
    cms_toplmiddle->SetFillStyle(1001);
  //  cms_toplmiddle->AddText("L = 34.5 pb^{-1}");
    cms_toplmiddle->AddText("L = 3.2 pb^{-1}");
    cms_toplmiddle->Draw();

    TPaveText *cms_toprmiddle = new TPaveText(0.51,0.99,0.74,0.94,"NDC"); // NDC sets coords
    cms_toprmiddle->SetTextSize(0.03);
    cms_toprmiddle->SetBorderSize(0); 
    cms_toprmiddle->SetTextAlign(22);
    cms_toprmiddle->SetLineWidth(1);
    cms_toprmiddle->SetLineColor(0);
    cms_toprmiddle->SetFillColor(0);
    cms_toprmiddle->SetFillStyle(1001);
    cms_toprmiddle->AddText("pp #rightarrow 2 jets + X");
    cms_toprmiddle->Draw();

    TPaveText *cms_topright = new TPaveText(0.76,0.99,0.99,0.94,"NDC"); // NDC sets coords
    cms_topright->SetTextSize(0.03);
    cms_topright->SetBorderSize(0); 
    cms_topright->SetTextAlign(22);
    cms_topright->SetLineWidth(1);
    cms_topright->SetLineColor(0);
    cms_topright->SetFillColor(0);
    cms_topright->SetFillStyle(1001);
    cms_topright->AddText("Anti-k_{T} (R = 0.5)");
    cms_topright->Draw();

    TPaveText *cms = new TPaveText(x1,y1,x2,y2,"NDC"); // NDC sets coords
    cms->SetTextSize(0.03);
    cms->SetBorderSize(0); 
    cms->SetTextAlign(22);
    cms->SetLineWidth(1);
    cms->SetLineColor(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(1001);
    cms->AddText("p_{T}^{jet} > 35 GeV and |#eta| < 2.8");
    cms->AddText("p_{T}^{jet} > 35 GeV and 3.2 < |#eta| < 4.7");
    cms->Draw();

 print_plots(c1, output_path, "xsec" + fileout);


//scaling back histograms
        data_unc_A->SetMaximum(max/1000.);

	data_B->Scale(1./factor_B);
	data_unc_B->Scale(1./factor_B);
	if (n_items >= 1) { mc1_B->Scale(1./factor_B); }
	if (n_items >= 2) { mc2_B->Scale(1./factor_B); }
	if (n_items >= 3) { mc3_B->Scale(1./factor_B); }
	if (n_items >= 4) { mc4_B->Scale(1./factor_B); }
	if (n_items == 5) { mc5_B->Scale(1./factor_B); }

	data_C->Scale(1./factor_C);
	data_unc_C->Scale(1./factor_C);
	if (n_items >= 1) { mc1_C->Scale(1./factor_C); }
	if (n_items >= 2) { mc2_C->Scale(1./factor_C); }
	if (n_items >= 3) { mc3_C->Scale(1./factor_C); }
	if (n_items >= 4) { mc4_C->Scale(1./factor_C); }
	if (n_items == 5) { mc5_C->Scale(1./factor_C); }

}


void create_detamerged_plot(TH1D *data_deta1, TString label_deta1, TH1D *data_deta2, TString label_deta2, TH1D *data_deta3, TString label_deta3, TH1D *data_deta4, TString label_deta4, TH1D *data_unc_deta1, TH1D *data_unc_deta2, TH1D *data_unc_deta3, TH1D *data_unc_deta4, TH1 *mc1_deta1, TH1 *mc1_deta2, TH1 *mc1_deta3, TH1 *mc1_deta4, TH1 *mc2_deta1, TH1 *mc2_deta2, TH1 *mc2_deta3, TH1 *mc2_deta4, TH1 *mc3_deta1, TH1 *mc3_deta2, TH1 *mc3_deta3, TH1 *mc3_deta4, TH1 *mc4_deta1, TH1 *mc4_deta2, TH1 *mc4_deta3, TH1 *mc4_deta4, TH1 *mc5_deta1, TH1 *mc5_deta2, TH1 *mc5_deta3, TH1 *mc5_deta4, TString *list_mc, TString scenario, string output_path, string fileout, string legend1_position, string legend2_position, string label_position, bool detail = false)
{
//plots the model uncertainty control plots

if (detail) { cout << "Positions : " << legend1_position << " and " << label_position << endl; }

int n_items = 0;
if (mc1_deta1 != 0) { n_items = 1; }
if (mc2_deta1 != 0 and n_items == 1) { n_items = 2; }
if (mc3_deta1 != 0 and n_items == 2) { n_items = 3; }
if (mc4_deta1 != 0 and n_items == 3) { n_items = 4; }
if (mc5_deta1 != 0 and n_items == 4) { n_items = 5; }

//declaring the canvas
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,800,600);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetLogy();
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.07);
    gPad->SetFrameBorderMode(0);

//scaling factors
	double factor_deta2 = 100;
	double factor_deta3 = 10000;
	double factor_deta4 = 1000000;

//scaling histograms
	data_deta2->Scale(factor_deta2);
	data_unc_deta2->Scale(factor_deta2);
	if (n_items >= 1) { mc1_deta2->Scale(factor_deta2); }
	if (n_items >= 2) { mc2_deta2->Scale(factor_deta2); }
	if (n_items >= 3) { mc3_deta2->Scale(factor_deta2); }
	if (n_items >= 4) { mc4_deta2->Scale(factor_deta2); }
	if (n_items == 5) { mc5_deta2->Scale(factor_deta2); }

	data_deta3->Scale(factor_deta3);
	data_unc_deta3->Scale(factor_deta3);
	if (n_items >= 1) { mc1_deta3->Scale(factor_deta3); }
	if (n_items >= 2) { mc2_deta3->Scale(factor_deta3); }
	if (n_items >= 3) { mc3_deta3->Scale(factor_deta3); }
	if (n_items >= 4) { mc4_deta3->Scale(factor_deta3); }
	if (n_items == 5) { mc5_deta3->Scale(factor_deta3); }

	data_deta4->Scale(factor_deta4);
	data_unc_deta4->Scale(factor_deta4);
	if (n_items >= 1) { mc1_deta4->Scale(factor_deta4); }
	if (n_items >= 2) { mc2_deta4->Scale(factor_deta4); }
	if (n_items >= 3) { mc3_deta4->Scale(factor_deta4); }
	if (n_items >= 4) { mc4_deta4->Scale(factor_deta4); }
	if (n_items == 5) { mc5_deta4->Scale(factor_deta4); }

//calculate the plooting range
    	double min = 0.0;
    	double max = data_deta1->GetMaximum();
    	if (data_deta1->GetMinimum() == 0.0)
    	{
    	min = get_non0_minimum(data_deta1,detail);
   	}
    	else
    	{
    	min = data_deta1->GetMinimum();
    	}    
    
    	set_histogram_min_max(data_unc_deta1, min, max, detail);
    	if (n_items >= 1) { set_histogram_min_max(mc1_deta1, min, max, detail); }
    	if (n_items >= 2) { set_histogram_min_max(mc2_deta1, min, max, detail); }
    	if (n_items >= 3) { set_histogram_min_max(mc3_deta1, min, max, detail); }
    	if (n_items >= 4) { set_histogram_min_max(mc4_deta1, min, max, detail); }
    	if (n_items == 5) { set_histogram_min_max(mc5_deta1, min, max, detail); }

    	max = 10e15;
    	min = 0.1 * min;

//format and ploting the histogram
    data_unc_deta1->SetMinimum(min);
    data_unc_deta1->SetMaximum(max);
    data_unc_deta1->SetFillColor(5);
    data_unc_deta1->DrawCopy("e2");
    data_unc_deta2->SetFillColor(5);
    data_unc_deta2->DrawCopy("e2 same");
    data_unc_deta3->SetFillColor(5);
    data_unc_deta3->DrawCopy("e2 same");
    data_unc_deta4->SetFillColor(5);
    data_unc_deta4->DrawCopy("e2 same");
    data_deta1->SetLineStyle(1);
    data_deta1->SetLineWidth(3);
    data_deta1->SetLineColor(1);
    data_deta1->SetMarkerColor(1);
    data_deta1->SetMarkerStyle(20);
    data_deta1->SetMarkerSize(2);
    data_deta1->DrawCopy("e1 same");
    data_deta2->SetLineStyle(1);
    data_deta2->SetLineWidth(3);
    data_deta2->SetLineColor(1);
    data_deta2->SetMarkerColor(2);
    data_deta2->SetMarkerStyle(21);
    data_deta2->SetMarkerSize(2);
    data_deta2->DrawCopy("e1 same");
    data_deta3->SetLineStyle(1);
    data_deta3->SetLineWidth(3);
    data_deta3->SetLineColor(1);
    data_deta3->SetMarkerColor(3);
    data_deta3->SetMarkerStyle(22);
    data_deta3->SetMarkerSize(2);
    data_deta3->DrawCopy("e1 same");
    data_deta4->SetLineStyle(1);
    data_deta4->SetLineWidth(3);
    data_deta4->SetLineColor(1);
    data_deta4->SetMarkerColor(4);
    data_deta4->SetMarkerStyle(23);
    data_deta4->SetMarkerSize(2);
    data_deta4->DrawCopy("e1 same");
	if (n_items >= 1)
	{
    	mc1_deta1->SetLineStyle(1);
    	mc1_deta1->SetLineWidth(3);
    	mc1_deta1->SetLineColor(2);
    	//mc1_deta1->SetMarkerColor(2);
    	//mc1_deta1->SetMarkerStyle(20);
    	//mc1_deta1->SetMarkerSize(3);
    	mc1_deta1->DrawCopy("hsame");
	//mc1_deta1->DrawCopy("esame");
    	mc1_deta2->SetLineStyle(1);
    	mc1_deta2->SetLineWidth(3);
    	mc1_deta2->SetLineColor(2);
    	//mc1_deta2->SetMarkerColor(2);
    	//mc1_deta2->SetMarkerStyle(20);
    	//mc1_deta2->SetMarkerSize(3);
    	mc1_deta2->DrawCopy("hsame");
	//mc1_deta2->DrawCopy("esame");
    	mc1_deta3->SetLineStyle(1);
    	mc1_deta3->SetLineWidth(3);
    	mc1_deta3->SetLineColor(2);
    	//mc1_deta3->SetMarkerColor(2);
    	//mc1_deta3->SetMarkerStyle(20);
    	//mc1_deta3->SetMarkerSize(3);
    	mc1_deta3->DrawCopy("hsame");
	//mc1_deta3->DrawCopy("esame");
    	mc1_deta4->SetLineStyle(1);
    	mc1_deta4->SetLineWidth(3);
    	mc1_deta4->SetLineColor(2);
    	//mc1_deta4->SetMarkerColor(2);
    	//mc1_deta4->SetMarkerStyle(20);
    	//mc1_deta4->SetMarkerSize(3);
    	mc1_deta4->DrawCopy("hsame");
	//mc1_deta4->DrawCopy("esame");
	}
	if (n_items >= 2)
	{
    	mc2_deta1->SetLineStyle(2);
    	mc2_deta1->SetLineWidth(3);
    	mc2_deta1->SetLineColor(4);
    	//mc2_deta1->SetMarkerColor(4);
    	//mc2_deta1->SetMarkerStyle(21);
    	//mc2_deta1->SetMarkerSize(3);
    	mc2_deta1->DrawCopy("hsame");
	mc2_deta1->DrawCopy("esame");
    	mc2_deta2->SetLineStyle(2);
    	mc2_deta2->SetLineWidth(3);
    	mc2_deta2->SetLineColor(4);
    	//mc2_deta2->SetMarkerColor(4);
    	//mc2_deta2->SetMarkerStyle(21);
    	//mc2_deta2->SetMarkerSize(3);
    	mc2_deta2->DrawCopy("hsame");
	mc2_deta2->DrawCopy("esame");
    	mc2_deta3->SetLineStyle(2);
    	mc2_deta3->SetLineWidth(3);
    	mc2_deta3->SetLineColor(4);
    	//mc2_deta3->SetMarkerColor(4);
    	//mc2_deta3->SetMarkerStyle(21);
    	//mc2_deta3->SetMarkerSize(3);
    	mc2_deta3->DrawCopy("hsame");
	mc2_deta3->DrawCopy("esame");
    	mc2_deta4->SetLineStyle(2);
    	mc2_deta4->SetLineWidth(3);
    	mc2_deta4->SetLineColor(4);
    	//mc2_deta4->SetMarkerColor(4);
    	//mc2_deta4->SetMarkerStyle(21);
    	//mc2_deta4->SetMarkerSize(3);
    	mc2_deta4->DrawCopy("hsame");
	mc2_deta4->DrawCopy("esame");
	}
	if (n_items >= 3)
	{
    	mc3_deta1->SetLineStyle(10);
    	mc3_deta1->SetLineWidth(3);
    	mc3_deta1->SetLineColor(6);
    	//mc3_deta1->SetMarkerColor(6);
    	//mc3_deta1->SetMarkerStyle(22);
    	//mc3_deta1->SetMarkerSize(3);
    	mc3_deta1->DrawCopy("hsame");
	mc3_deta1->DrawCopy("esame");
    	mc3_deta2->SetLineStyle(10);
    	mc3_deta2->SetLineWidth(3);
    	mc3_deta2->SetLineColor(6);
    	//mc3_deta2->SetMarkerColor(6);
    	//mc3_deta2->SetMarkerStyle(22);
    	//mc3_deta2->SetMarkerSize(3);
    	mc3_deta2->DrawCopy("hsame");
	mc3_deta2->DrawCopy("esame");
    	mc3_deta3->SetLineStyle(10);
    	mc3_deta3->SetLineWidth(3);
    	mc3_deta3->SetLineColor(6);
    	//mc3_deta3->SetMarkerColor(6);
    	//mc3_deta3->SetMarkerStyle(22);
    	//mc3_deta3->SetMarkerSize(3);
    	mc3_deta3->DrawCopy("hsame");
	mc3_deta3->DrawCopy("esame");
    	mc3_deta4->SetLineStyle(10);
    	mc3_deta4->SetLineWidth(3);
    	mc3_deta4->SetLineColor(6);
    	//mc3_deta4->SetMarkerColor(6);
    	//mc3_deta4->SetMarkerStyle(22);
    	//mc3_deta4->SetMarkerSize(3);
    	mc3_deta4->DrawCopy("hsame");
	mc3_deta4->DrawCopy("esame");
	}
	if (n_items >= 4)
	{
    	mc4_deta1->SetLineStyle(3);
    	mc4_deta1->SetLineWidth(3);
    	mc4_deta1->SetLineColor(9);
    	//mc4_deta1->SetMarkerColor(7);
    	//mc4_deta1->SetMarkerStyle(23);
    	//mc4_deta1->SetMarkerSize(3);
    	mc4_deta1->DrawCopy("hsame");
	mc4_deta1->DrawCopy("esame");
    	mc4_deta2->SetLineStyle(3);
    	mc4_deta2->SetLineWidth(3);
    	mc4_deta2->SetLineColor(9);
    	//mc4_deta2->SetMarkerColor(7);
    	//mc4_deta2->SetMarkerStyle(23);
    	//mc4_deta2->SetMarkerSize(3);
    	mc4_deta2->DrawCopy("hsame");
	mc4_deta2->DrawCopy("esame");
    	mc4_deta3->SetLineStyle(3);
    	mc4_deta3->SetLineWidth(3);
    	mc4_deta3->SetLineColor(9);
    	//mc4_deta3->SetMarkerColor(7);
    	//mc4_deta3->SetMarkerStyle(23);
    	//mc4_deta3->SetMarkerSize(3);
    	mc4_deta3->DrawCopy("hsame");
	mc4_deta3->DrawCopy("esame");
    	mc4_deta4->SetLineStyle(3);
    	mc4_deta4->SetLineWidth(3);
    	mc4_deta4->SetLineColor(9);
    	//mc4_deta4->SetMarkerColor(7);
    	//mc4_deta4->SetMarkerStyle(23);
    	//mc4_deta4->SetMarkerSize(3);
    	mc4_deta4->DrawCopy("hsame");
	mc4_deta4->DrawCopy("esame");
	}
	if (n_items >= 5)
	{
    	mc5_deta1->SetLineStyle(10);
    	mc5_deta1->SetLineWidth(3);
    	mc5_deta1->SetLineColor(8);
    	//mc5_deta1->SetMarkerColor(8);
    	//mc5_deta1->SetMarkerStyle(24);
    	//mc5_deta1->SetMarkerSize(3);
    	mc5_deta1->DrawCopy("hsame");
	mc5_deta1->DrawCopy("esame");
    	mc5_deta2->SetLineStyle(10);
    	mc5_deta2->SetLineWidth(3);
    	mc5_deta2->SetLineColor(8);
    	//mc5_deta2->SetMarkerColor(8);
    	//mc5_deta2->SetMarkerStyle(24);
    	//mc5_deta2->SetMarkerSize(3);
    	mc5_deta2->DrawCopy("hsame");
	mc5_deta2->DrawCopy("esame");
    	mc5_deta3->SetLineStyle(10);
    	mc5_deta3->SetLineWidth(3);
    	mc5_deta3->SetLineColor(8);
    	//mc5_deta3->SetMarkerColor(8);
    	//mc5_deta3->SetMarkerStyle(24);
    	//mc5_deta3->SetMarkerSize(3);
    	mc5_deta3->DrawCopy("hsame");
	mc5_deta3->DrawCopy("esame");
    	mc5_deta4->SetLineStyle(10);
    	mc5_deta4->SetLineWidth(3);
    	mc5_deta4->SetLineColor(8);
    	//mc5_deta4->SetMarkerColor(8);
    	//mc5_deta4->SetMarkerStyle(24);
    	//mc5_deta4->SetMarkerSize(3);
    	mc5_deta4->DrawCopy("hsame");
	mc5_deta4->DrawCopy("esame");
	}
    
//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

//sets and draw the legend
    set_legend_position(legend1_position, n_items, x1, y1, x2, y2);
    TLegend *leg00 = new TLegend(x1, y1-0.06, x2-0.10, y2-0.06);
    if (n_items >= 1) { leg00->AddEntry(mc1_deta1,list_mc[0],"lp"); }
    if (n_items >= 2) {	leg00->AddEntry(mc2_deta1,list_mc[1],"lp"); }
    if (n_items >= 3) {	leg00->AddEntry(mc3_deta1,list_mc[2],"lp"); }
    if (n_items >= 4) { leg00->AddEntry(mc4_deta1,list_mc[3],"lp"); }
    if (n_items == 5) { leg00->AddEntry(mc5_deta1,list_mc[4],"lp"); }
    leg00->SetFillColor(0);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();

    set_legend_position(legend2_position, 3, x1, y1, x2, y2);
    TLegend *leg01 = new TLegend(x1, y1-0.06, x2, y2-0.06);
    leg01->AddEntry(data_deta4,label_deta4+" (x10^{6})","lep");
    leg01->AddEntry(data_deta3,label_deta3+" (x10^{4})","lep");
    leg01->AddEntry(data_deta2,label_deta2+" (x10^{2})","lep");
    leg01->AddEntry(data_deta1,label_deta1+" (x1)","lep");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();
    
    if (label_position == "top_right")    { x1 = 0.63; y1 = 0.82; x2 = 0.98; y2 = 0.92; }
    if (label_position == "top_left")     { x1 = 0.13; y1 = 0.82; x2 = 0.47; y2 = 0.92; }
    if (label_position == "bottom_left")  { x1 = 0.13; y1 = 0.13; x2 = 0.47; y2 = 0.23; }
    if (label_position == "bottom_right") { x1 = 0.60; y1 = 0.13; x2 = 0.98; y2 = 0.23; }
    if (label_position == "middle")       { x1 = 0.50; y1 = 0.40; x2 = 0.60; y2 = 0.50; }

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
  //  cms->AddText("L_{int} = 34.5 pb^{-1}, Anti-k_{T} (R = 0.5)");
    cms->AddText("L_{int} = 3.2 pb^{-1}, Anti-k_{T} (R = 0.5)");
    cms->AddText("p_{T}^{jet} > 35 GeV and |#eta| < 2.8");
    cms->AddText("p_{T}^{jet} > 35 GeV and 3.2 < |#eta| < 4.7");
    cms->Draw();

 print_plots(c1, output_path, "xsec" + fileout);


//scaling back histograms
        data_unc_deta1->SetMaximum(max/1000.);

	data_deta2->Scale(1./factor_deta2);
	data_unc_deta2->Scale(1./factor_deta2);
	if (n_items >= 1) { mc1_deta2->Scale(1./factor_deta2); }
	if (n_items >= 2) { mc2_deta2->Scale(1./factor_deta2); }
	if (n_items >= 3) { mc3_deta2->Scale(1./factor_deta2); }
	if (n_items >= 4) { mc4_deta2->Scale(1./factor_deta2); }
	if (n_items == 5) { mc5_deta2->Scale(1./factor_deta2); }

	data_deta3->Scale(1./factor_deta3);
	data_unc_deta3->Scale(1./factor_deta3);
	if (n_items >= 1) { mc1_deta3->Scale(1./factor_deta3); }
	if (n_items >= 2) { mc2_deta3->Scale(1./factor_deta3); }
	if (n_items >= 3) { mc3_deta3->Scale(1./factor_deta3); }
	if (n_items >= 4) { mc4_deta3->Scale(1./factor_deta3); }
	if (n_items == 5) { mc5_deta3->Scale(1./factor_deta3); }

	data_deta4->Scale(1./factor_deta4);
	data_unc_deta4->Scale(1./factor_deta4);
	if (n_items >= 1) { mc1_deta4->Scale(1./factor_deta4); }
	if (n_items >= 2) { mc2_deta4->Scale(1./factor_deta4); }
	if (n_items >= 3) { mc3_deta4->Scale(1./factor_deta4); }
	if (n_items >= 4) { mc4_deta4->Scale(1./factor_deta4); }
	if (n_items == 5) { mc5_deta4->Scale(1./factor_deta4); }

}


void final_plots(string path_data, string path_data_unc, string *path_mc, string *labels, string *prefixes, string output_path_plots, string mode, bool detail = false, bool test = false)
{

    int i = 0;
    TString list_mc[5];

//outputs the configuration
    if (detail) { cout << "Final Plots Configuration"<<endl; }
    if (detail) { cout << "Input path for data:                    " << path_data << endl; }
    if (detail) { cout << "Input path for data with uncertainties: " << path_data_unc << endl; }
    while (path_mc[i] != "")
	{
	if (detail)
		{
		cout << "MC Input " << i+1 << " :                            " << path_mc[i] << endl;
		cout << "Label " << i+1 << " :                               " << labels[i] << " and prefix " << prefixes[i] << endl;
		}
	i = i + 1;
	}
    if (detail) { cout << "Output Path Plots:                      " << output_path_plots << endl; }
    if (detail) { cout << "Mode:                                   " << mode << endl; }
    if (detail) { cout << "Detail Level:                           " << detail << endl; }
    if (detail) { cout << "Test Mode:                              " << test << endl; }

//opening the files
    if (detail) { cout << "Opening the data files" << endl; }
    TFile *data = new TFile( path_data.c_str() );
    TFile *data_unc = new TFile( path_data_unc.c_str() );
    TFile *mc1 = 0;
    if (i >= 1)
	{
	mc1 = new TFile( path_mc[0].c_str() );
	if (mc1 == 0) { cout << "File " << path_mc[0] << "not found!" << endl; i = 0;}
	}
    TFile *mc2 = 0;
    if (i >= 2)
	{
	mc2 = new TFile( path_mc[1].c_str() );
	if (mc2 == 0) { cout << "File " << path_mc[1] << "not found!" << endl; i = 1;}
	}
    TFile *mc3 = 0;
    if (i >= 3)
	{
	mc3 = new TFile( path_mc[2].c_str() );
	if (mc3 == 0) { cout << "File " << path_mc[2] << "not found!" << endl; i = 2;}
	}
    TFile *mc4 = 0;
    if (i >= 4)
	{
	mc4 = new TFile( path_mc[3].c_str() );
	if (mc4 == 0) { cout << "File " << path_mc[3] << "not found!" << endl; i = 3;}
	}
    TFile *mc5 = 0;
    if (i >= 5)
	{
	mc5 = new TFile( path_mc[4].c_str() );
	if (mc5 == 0) { cout << "File " << path_mc[4] << "not found!" << endl; i = 4;}
	}
    TFile *mc6 = 0;
    if (i >= 6)
	{
	mc6 = new TFile( path_mc[5].c_str() );
	if (mc6 == 0) { cout << "File " << path_mc[5] << "not found!" << endl; i = 5;}
	}
    TFile *mc7 = 0;
    if (i >= 7)
	{
	mc7 = new TFile( path_mc[6].c_str() );
	if (mc7 == 0) { cout << "File " << path_mc[6] << "not found!" << endl; i = 6;}
	}
    TFile *mc8 = 0;
    if (i >= 8)
	{
	mc8 = new TFile( path_mc[7].c_str() );
	if (mc8 == 0) { cout << "File " << path_mc[7] << "not found!" << endl; i = 7;}
	}
    TFile *mc9 = 0;
    if (i >= 9)
	{
	mc9 = new TFile( path_mc[8].c_str() );
	if (mc9 == 0) { cout << "File " << path_mc[8] << "not found!" << endl; i = 8;}
	}
    TFile *mc10 = 0;
    if (i >= 10)
	{
	mc10 = new TFile( path_mc[9].c_str() );
	if (mc10 == 0) { cout << "File " << path_mc[9] << "not found!" << endl; i = 9;}
	}
    TFile *mc11 = 0;
    if (i >= 11)
	{
	mc11 = new TFile( path_mc[10].c_str() );
	if (mc11 == 0) { cout << "File " << path_mc[10] << " not found!" << endl; i = 10;}
	}
    TFile *mc12 = 0;
    if (i >= 12)
	{
	mc12 = new TFile( path_mc[11].c_str() );
	if (mc12 == 0) { cout << "File " << path_mc[11] << "not found!" << endl; i = 11;}
	}
    TFile *mc13 = 0;
    if (i >= 13)
	{
	mc13 = new TFile( path_mc[12].c_str() );
	if (mc13 == 0) { cout << "File " << path_mc[12] << "not found!" << endl; i = 12;}
	}
    TFile *mc14 = 0;
    if (i >= 14)
	{
	mc14 = new TFile( path_mc[13].c_str() );
	if (mc14 == 0) { cout << "File " << path_mc[13] << "not found!" << endl; i = 13;}
	}
    TFile *mc15 = 0;
    if (i >= 15)
	{
	mc15 = new TFile( path_mc[14].c_str() );
	if (mc15 == 0) { cout << "File " << path_mc[14] << "not found!" << endl; i = 14;}
	}
    TFile *mc16 = 0;
    if (i >= 16)
	{
	mc16 = new TFile( path_mc[15].c_str() );
	if (mc16 == 0) { cout << "File " << path_mc[15] << "not found!" << endl; i = 15;}
	}
    TFile *mc17 = 0;
    if (i >= 17)
	{
	mc17 = new TFile( path_mc[16].c_str() );
	if (mc17 == 0) { cout << "File " << path_mc[16] << "not found!" << endl; i = 16;}
	}
    TFile *mc18 = 0;
    if (i >= 18)
	{
	mc18 = new TFile( path_mc[17].c_str() );
	if (mc18 == 0) { cout << "File " << path_mc[17] << "not found!" << endl; i = 17;}
	}
    TFile *mc19 = 0;
    if (i >= 19)
	{
	mc19 = new TFile( path_mc[18].c_str() );
	if (mc19 == 0) { cout << "File " << path_mc[18] << "not found!" << endl; i = 18;}
	}
    TFile *mc20 = 0;
    if (i >= 20)
	{
	mc20 = new TFile( path_mc[19].c_str() );
	if (mc20 == 0) { cout << "File " << path_mc[19] << "not found!" << endl; i = 19;}
	}


//declaring temp_strings
    string legend_position, legend1_position, label_position, extra_label_position;

//final plot for delta phi distribution
    if (detail) { cout<<"Delta phi"<<endl; }

    TH1D *delta_phi_data = 0;
    TH1D *delta_phi_data_unc = 0;
    TH1D *delta_phi_mc1 = 0;
    TH1D *delta_phi_mc2 = 0;
    TH1 *delta_phi_mc3 = 0;
    TH1 *delta_phi_mc4 = 0;
    TH1 *delta_phi_mc5 = 0;
    TH1 *delta_phi_mc6 = 0;
    TH1 *delta_phi_mc7 = 0;
    TH1 *delta_phi_mc8 = 0;
    TH1 *delta_phi_mc9 = 0;
    TH1 *delta_phi_mc10 = 0;
    TH1 *delta_phi_mc11 = 0;
    TH1 *delta_phi_mc12 = 0;
    TH1 *delta_phi_mc13 = 0;
    TH1 *delta_phi_mc14 = 0;
    TH1 *delta_phi_mc15 = 0;
    TH1 *delta_phi_mc16 = 0;
    TH1 *delta_phi_mc17 = 0;
    TH1 *delta_phi_mc18 = 0;
    TH1 *delta_phi_mc19 = 0;
    TH1 *delta_phi_mc20 = 0;

    data->GetObject("ak5PF_delta_phi",delta_phi_data);
    if (delta_phi_data == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi",delta_phi_data_unc);
    if (delta_phi_data_unc == 0) { cout << "ak5PF_delta_phi unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi",delta_phi_mc1);
	if (delta_phi_mc1 == 0) { cout << "ak5Gen_delta_phi mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi",delta_phi_mc2);
	if (delta_phi_mc2 == 0) { cout << "ak5Gen_delta_phi mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d01-x01-y01",delta_phi_mc3);
	if (delta_phi_mc3 == 0) { cout << "ak5Gen_delta_phi mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d01-x01-y01",delta_phi_mc4);
	if (delta_phi_mc4 == 0) { cout << "ak5Gen_delta_phi mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d01-x01-y01",delta_phi_mc5);
	if (delta_phi_mc5 == 0) { cout << "ak5Gen_delta_phi mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d01-x01-y01",delta_phi_mc6);
	if (delta_phi_mc6 == 0) { cout << "ak5Gen_delta_phi mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d01-x01-y01",delta_phi_mc7);
	if (delta_phi_mc7 == 0) { cout << "ak5Gen_delta_phi mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d01-x01-y01",delta_phi_mc8);
	if (delta_phi_mc8 == 0) { cout << "ak5Gen_delta_phi mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d01-x01-y01",delta_phi_mc9);
	if (delta_phi_mc9 == 0) { cout << "ak5Gen_delta_phi mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d01-x01-y01",delta_phi_mc10);
	if (delta_phi_mc10 == 0) { cout << "ak5Gen_delta_phi mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d01-x01-y01",delta_phi_mc11);
	if (delta_phi_mc11 == 0) { cout << "ak5Gen_delta_phi mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d01-x01-y01",delta_phi_mc12);
	if (delta_phi_mc12 == 0) { cout << "ak5Gen_delta_phi mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d01-x01-y01",delta_phi_mc13);
	if (delta_phi_mc13 == 0) { cout << "ak5Gen_delta_phi mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d01-x01-y01",delta_phi_mc14);
	if (delta_phi_mc14 == 0) { cout << "ak5Gen_delta_phi mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d01-x01-y01",delta_phi_mc15);
	if (delta_phi_mc15 == 0) { cout << "ak5Gen_delta_phi mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d01-x01-y01",delta_phi_mc16);
	if (delta_phi_mc16 == 0) { cout << "ak5Gen_delta_phi mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d01-x01-y01",delta_phi_mc17);
	if (delta_phi_mc17 == 0) { cout << "ak5Gen_delta_phi mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d01-x01-y01",delta_phi_mc18);
	if (delta_phi_mc18 == 0) { cout << "ak5Gen_delta_phi mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d01-x01-y01",delta_phi_mc19);
	if (delta_phi_mc19 == 0) { cout << "ak5Gen_delta_phi mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d01-x01-y01",delta_phi_mc20);
	if (delta_phi_mc20 == 0) { cout << "ak5Gen_delta_phi mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	legend_position = "top_left";
	label_position = "bottom_right";
	}
    if (mode == "ratio")
	{
	legend_position = "top_right";
	label_position = "top_left";
	}

    create_final_plot(delta_phi_data, delta_phi_data_unc, delta_phi_mc3, labels[2], delta_phi_mc4, labels[3], delta_phi_mc9, labels[8], delta_phi_mc20, labels[19], delta_phi_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_all", legend_position, label_position, "", "", detail);
    create_final_plot(delta_phi_data, delta_phi_data_unc, delta_phi_mc5, labels[4], delta_phi_mc8, labels[7], delta_phi_mc7, labels[6], delta_phi_mc6, labels[5], delta_phi_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_pythia6", legend_position, label_position, "", "", detail);
    create_final_plot(delta_phi_data, delta_phi_data_unc, delta_phi_mc1, labels[0], delta_phi_mc2, labels[1], delta_phi_mc7, labels[6], delta_phi_mc9, labels[8], delta_phi_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_check", legend_position, label_position, "", "", detail);
    create_final_plot(delta_phi_data, delta_phi_data_unc, delta_phi_mc7, labels[6], delta_phi_mc6, labels[5], delta_phi_mc10, labels[9], delta_phi_mc11, labels[10], delta_phi_mc12, labels[11], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_hadronization", legend_position, label_position, "", "", detail);

//final plot for delta phi nogap distribution
    if (detail) { cout<<"Delta phi noGap"<<endl; }

    TH1D *delta_phi_nogap_data = 0;
    TH1D *delta_phi_nogap_data_unc = 0;
    TH1D *delta_phi_nogap_mc1 = 0;
    TH1D *delta_phi_nogap_mc2 = 0;
    TH1 *delta_phi_nogap_mc3 = 0;
    TH1 *delta_phi_nogap_mc4 = 0;
    TH1 *delta_phi_nogap_mc5 = 0;
    TH1 *delta_phi_nogap_mc6 = 0;
    TH1 *delta_phi_nogap_mc7 = 0;
    TH1 *delta_phi_nogap_mc8 = 0;
    TH1 *delta_phi_nogap_mc9 = 0;
    TH1 *delta_phi_nogap_mc10 = 0;
    TH1 *delta_phi_nogap_mc11 = 0;
    TH1 *delta_phi_nogap_mc12 = 0;
    TH1 *delta_phi_nogap_mc13 = 0;
    TH1 *delta_phi_nogap_mc14 = 0;
    TH1 *delta_phi_nogap_mc15 = 0;
    TH1 *delta_phi_nogap_mc16 = 0;
    TH1 *delta_phi_nogap_mc17 = 0;
    TH1 *delta_phi_nogap_mc18 = 0;
    TH1 *delta_phi_nogap_mc19 = 0;
    TH1 *delta_phi_nogap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data);
    if (delta_phi_nogap_data == 0) { cout << "ak5PF_delta_phi_nogap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_nogap",delta_phi_nogap_data_unc);
    if (delta_phi_nogap_data_unc == 0) { cout << "ak5PF_delta_phi_nogap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_nogap",delta_phi_nogap_mc1);
	if (delta_phi_nogap_mc1 == 0) { cout << "ak5Gen_delta_phi_nogap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_nogap",delta_phi_nogap_mc2);
	if (delta_phi_nogap_mc2 == 0) { cout << "ak5Gen_delta_phi_nogap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d02-x01-y01",delta_phi_nogap_mc3);
	if (delta_phi_nogap_mc3 == 0) { cout << "ak5Gen_delta_phi_nogap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d02-x01-y01",delta_phi_nogap_mc4);
	if (delta_phi_nogap_mc4 == 0) { cout << "ak5Gen_delta_phi_nogap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d02-x01-y01",delta_phi_nogap_mc5);
	if (delta_phi_nogap_mc5 == 0) { cout << "ak5Gen_delta_phi_nogap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d02-x01-y01",delta_phi_nogap_mc6);
	if (delta_phi_nogap_mc6 == 0) { cout << "ak5Gen_delta_phi_nogap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d02-x01-y01",delta_phi_nogap_mc7);
	if (delta_phi_nogap_mc7 == 0) { cout << "ak5Gen_delta_phi_nogap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d02-x01-y01",delta_phi_nogap_mc8);
	if (delta_phi_nogap_mc8 == 0) { cout << "ak5Gen_delta_phi_nogap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d02-x01-y01",delta_phi_nogap_mc9);
	if (delta_phi_nogap_mc9 == 0) { cout << "ak5Gen_delta_phi_nogap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d02-x01-y01",delta_phi_nogap_mc10);
	if (delta_phi_nogap_mc10 == 0) { cout << "ak5Gen_delta_phi_nogap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d02-x01-y01",delta_phi_nogap_mc11);
	if (delta_phi_nogap_mc11 == 0) { cout << "ak5Gen_delta_phi_nogap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d02-x01-y01",delta_phi_nogap_mc12);
	if (delta_phi_nogap_mc12 == 0) { cout << "ak5Gen_delta_phi_nogap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d02-x01-y01",delta_phi_nogap_mc13);
	if (delta_phi_nogap_mc13 == 0) { cout << "ak5Gen_delta_phi_nogap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d02-x01-y01",delta_phi_nogap_mc14);
	if (delta_phi_nogap_mc14 == 0) { cout << "ak5Gen_delta_phi_nogap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d02-x01-y01",delta_phi_nogap_mc15);
	if (delta_phi_nogap_mc15 == 0) { cout << "ak5Gen_delta_phi_nogap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d02-x01-y01",delta_phi_nogap_mc16);
	if (delta_phi_nogap_mc16 == 0) { cout << "ak5Gen_delta_phi_nogap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d02-x01-y01",delta_phi_nogap_mc17);
	if (delta_phi_nogap_mc17 == 0) { cout << "ak5Gen_delta_phi_nogap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d02-x01-y01",delta_phi_nogap_mc18);
	if (delta_phi_nogap_mc18 == 0) { cout << "ak5Gen_delta_phi_nogap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d02-x01-y01",delta_phi_nogap_mc19);
	if (delta_phi_nogap_mc19 == 0) { cout << "ak5Gen_delta_phi_nogap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d02-x01-y01",delta_phi_nogap_mc20);
	if (delta_phi_nogap_mc20 == 0) { cout << "ak5Gen_delta_phi_nogap mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_nogap_data, delta_phi_nogap_data_unc, delta_phi_nogap_mc3, labels[2], delta_phi_nogap_mc4, labels[3], delta_phi_nogap_mc9, labels[8], delta_phi_nogap_mc20, labels[19], delta_phi_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_nogap_all", legend_position, label_position);
    create_final_plot(delta_phi_nogap_data, delta_phi_nogap_data_unc, delta_phi_nogap_mc5, labels[4], delta_phi_nogap_mc8, labels[7], delta_phi_nogap_mc7, labels[6], delta_phi_nogap_mc6, labels[5], delta_phi_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_nogap_pythia6", legend_position, label_position);
    create_final_plot(delta_phi_nogap_data, delta_phi_nogap_data_unc, delta_phi_nogap_mc1, labels[0], delta_phi_nogap_mc2, labels[1], delta_phi_nogap_mc7, labels[6], delta_phi_nogap_mc9, labels[8], delta_phi_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_nogap_check", legend_position, label_position);
    create_final_plot(delta_phi_nogap_data, delta_phi_nogap_data_unc, delta_phi_nogap_mc7, labels[6], delta_phi_nogap_mc6, labels[5], delta_phi_nogap_mc10, labels[9], delta_phi_nogap_mc11, labels[10], delta_phi_nogap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_nogap_hadronization", legend_position, label_position);

//final plot for delta phi gap distribution
    if (detail) { cout<<"Delta phi Gap"<<endl; }

    TH1D *delta_phi_gap_data = 0;
    TH1D *delta_phi_gap_data_unc = 0;
    TH1D *delta_phi_gap_mc1 = 0;
    TH1D *delta_phi_gap_mc2 = 0;
    TH1 *delta_phi_gap_mc3 = 0;
    TH1 *delta_phi_gap_mc4 = 0;
    TH1 *delta_phi_gap_mc5 = 0;
    TH1 *delta_phi_gap_mc6 = 0;
    TH1 *delta_phi_gap_mc7 = 0;
    TH1 *delta_phi_gap_mc8 = 0;
    TH1 *delta_phi_gap_mc9 = 0;
    TH1 *delta_phi_gap_mc10 = 0;
    TH1 *delta_phi_gap_mc11 = 0;
    TH1 *delta_phi_gap_mc12 = 0;
    TH1 *delta_phi_gap_mc13 = 0;
    TH1 *delta_phi_gap_mc14 = 0;
    TH1 *delta_phi_gap_mc15 = 0;
    TH1 *delta_phi_gap_mc16 = 0;
    TH1 *delta_phi_gap_mc17 = 0;
    TH1 *delta_phi_gap_mc18 = 0;
    TH1 *delta_phi_gap_mc19 = 0;
    TH1 *delta_phi_gap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data);
    if (delta_phi_gap_data == 0) { cout << "ak5PF_delta_phi_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_gap",delta_phi_gap_data_unc);
    if (delta_phi_gap_data_unc == 0) { cout << "ak5PF_delta_phi_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_gap",delta_phi_gap_mc1);
	if (delta_phi_gap_mc1 == 0) { cout << "ak5Gen_delta_phi_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_gap",delta_phi_gap_mc2);
	if (delta_phi_gap_mc2 == 0) { cout << "ak5Gen_delta_phi_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d03-x01-y01",delta_phi_gap_mc3);
	if (delta_phi_gap_mc3 == 0) { cout << "ak5Gen_delta_phi_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d03-x01-y01",delta_phi_gap_mc4);
	if (delta_phi_gap_mc4 == 0) { cout << "ak5Gen_delta_phi_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d03-x01-y01",delta_phi_gap_mc5);
	if (delta_phi_gap_mc5 == 0) { cout << "ak5Gen_delta_phi_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d03-x01-y01",delta_phi_gap_mc6);
	if (delta_phi_gap_mc6 == 0) { cout << "ak5Gen_delta_phi_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d03-x01-y01",delta_phi_gap_mc7);
	if (delta_phi_gap_mc7 == 0) { cout << "ak5Gen_delta_phi_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d03-x01-y01",delta_phi_gap_mc8);
	if (delta_phi_gap_mc8 == 0) { cout << "ak5Gen_delta_phi_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d03-x01-y01",delta_phi_gap_mc9);
	if (delta_phi_gap_mc9 == 0) { cout << "ak5Gen_delta_phi_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d03-x01-y01",delta_phi_gap_mc10);
	if (delta_phi_gap_mc10 == 0) { cout << "ak5Gen_delta_phi_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d03-x01-y01",delta_phi_gap_mc11);
	if (delta_phi_gap_mc11 == 0) { cout << "ak5Gen_delta_phi_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d03-x01-y01",delta_phi_gap_mc12);
	if (delta_phi_gap_mc12 == 0) { cout << "ak5Gen_delta_phi_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d03-x01-y01",delta_phi_gap_mc13);
	if (delta_phi_gap_mc13 == 0) { cout << "ak5Gen_delta_phi_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d03-x01-y01",delta_phi_gap_mc14);
	if (delta_phi_gap_mc14 == 0) { cout << "ak5Gen_delta_phi_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d03-x01-y01",delta_phi_gap_mc15);
	if (delta_phi_gap_mc15 == 0) { cout << "ak5Gen_delta_phi_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d03-x01-y01",delta_phi_gap_mc16);
	if (delta_phi_gap_mc16 == 0) { cout << "ak5Gen_delta_phi_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d03-x01-y01",delta_phi_gap_mc17);
	if (delta_phi_gap_mc17 == 0) { cout << "ak5Gen_delta_phi_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d03-x01-y01",delta_phi_gap_mc18);
	if (delta_phi_gap_mc18 == 0) { cout << "ak5Gen_delta_phi_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d03-x01-y01",delta_phi_gap_mc19);
	if (delta_phi_gap_mc19 == 0) { cout << "ak5Gen_delta_phi_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d03-x01-y01",delta_phi_gap_mc20);
	if (delta_phi_gap_mc20 == 0) { cout << "ak5Gen_delta_phi_gap mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_gap_data, delta_phi_gap_data_unc, delta_phi_gap_mc3, labels[2], delta_phi_gap_mc4, labels[3], delta_phi_gap_mc9, labels[8], delta_phi_gap_mc20, labels[19], delta_phi_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_gap_all", legend_position, label_position);
    create_final_plot(delta_phi_gap_data, delta_phi_gap_data_unc, delta_phi_gap_mc5, labels[4], delta_phi_gap_mc8, labels[7], delta_phi_gap_mc7, labels[6], delta_phi_gap_mc6, labels[5], delta_phi_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_gap_pythia6", legend_position, label_position);
    create_final_plot(delta_phi_gap_data, delta_phi_gap_data_unc, delta_phi_gap_mc1, labels[0], delta_phi_gap_mc2, labels[1], delta_phi_gap_mc7, labels[6], delta_phi_gap_mc9, labels[8], delta_phi_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_gap_check", legend_position, label_position);
    create_final_plot(delta_phi_gap_data, delta_phi_gap_data_unc, delta_phi_gap_mc7, labels[6], delta_phi_gap_mc6, labels[5], delta_phi_gap_mc10, labels[9], delta_phi_gap_mc11, labels[10], delta_phi_gap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_gap_hadronization", legend_position, label_position);


//final plot for delta phi deta1 distribution
    if (detail) { cout<<"Delta phi deta1"<<endl; }

    TH1D *delta_phi_deta1_data = 0;
    TH1D *delta_phi_deta1_data_unc = 0;
    TH1D *delta_phi_deta1_mc1 = 0;
    TH1D *delta_phi_deta1_mc2 = 0;
    TH1 *delta_phi_deta1_mc3 = 0;
    TH1 *delta_phi_deta1_mc4 = 0;
    TH1 *delta_phi_deta1_mc5 = 0;
    TH1 *delta_phi_deta1_mc6 = 0;
    TH1 *delta_phi_deta1_mc7 = 0;
    TH1 *delta_phi_deta1_mc8 = 0;
    TH1 *delta_phi_deta1_mc9 = 0;
    TH1 *delta_phi_deta1_mc10 = 0;
    TH1 *delta_phi_deta1_mc11 = 0;
    TH1 *delta_phi_deta1_mc12 = 0;
    TH1 *delta_phi_deta1_mc13 = 0;
    TH1 *delta_phi_deta1_mc14 = 0;
    TH1 *delta_phi_deta1_mc15 = 0;
    TH1 *delta_phi_deta1_mc16 = 0;
    TH1 *delta_phi_deta1_mc17 = 0;
    TH1 *delta_phi_deta1_mc18 = 0;
    TH1 *delta_phi_deta1_mc19 = 0;
    TH1 *delta_phi_deta1_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data);
    if (delta_phi_deta1_data == 0) { cout << "ak5PF_delta_phi_deta1 not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta1",delta_phi_deta1_data_unc);
    if (delta_phi_deta1_data_unc == 0) { cout << "ak5PF_delta_phi_deta1 unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta1",delta_phi_deta1_mc1);
	if (delta_phi_deta1_mc1 == 0) { cout << "ak5Gen_delta_phi_deta1 mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta1",delta_phi_deta1_mc2);
	if (delta_phi_deta1_mc2 == 0) { cout << "ak5Gen_delta_phi_deta1 mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d01-x02-y01",delta_phi_deta1_mc3);
	if (delta_phi_deta1_mc3 == 0) { cout << "ak5Gen_delta_phi_deta1 mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d01-x02-y01",delta_phi_deta1_mc4);
	if (delta_phi_deta1_mc4 == 0) { cout << "ak5Gen_delta_phi_deta1 mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d01-x02-y01",delta_phi_deta1_mc5);
	if (delta_phi_deta1_mc5 == 0) { cout << "ak5Gen_delta_phi_deta1 mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d01-x02-y01",delta_phi_deta1_mc6);
	if (delta_phi_deta1_mc6 == 0) { cout << "ak5Gen_delta_phi_deta1 mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d01-x02-y01",delta_phi_deta1_mc7);
	if (delta_phi_deta1_mc7 == 0) { cout << "ak5Gen_delta_phi_deta1 mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d01-x02-y01",delta_phi_deta1_mc8);
	if (delta_phi_deta1_mc8 == 0) { cout << "ak5Gen_delta_phi_deta1 mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d01-x02-y01",delta_phi_deta1_mc9);
	if (delta_phi_deta1_mc9 == 0) { cout << "ak5Gen_delta_phi_deta1 mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d01-x02-y01",delta_phi_deta1_mc10);
	if (delta_phi_deta1_mc10 == 0) { cout << "ak5Gen_delta_phi_deta1 mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d01-x02-y01",delta_phi_deta1_mc11);
	if (delta_phi_deta1_mc11 == 0) { cout << "ak5Gen_delta_phi_deta1 mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d01-x02-y01",delta_phi_deta1_mc12);
	if (delta_phi_deta1_mc12 == 0) { cout << "ak5Gen_delta_phi_deta1 mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d01-x02-y01",delta_phi_deta1_mc13);
	if (delta_phi_deta1_mc13 == 0) { cout << "ak5Gen_delta_phi_deta1 mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d01-x02-y01",delta_phi_deta1_mc14);
	if (delta_phi_deta1_mc14 == 0) { cout << "ak5Gen_delta_phi_deta1 mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d01-x02-y01",delta_phi_deta1_mc15);
	if (delta_phi_deta1_mc15 == 0) { cout << "ak5Gen_delta_phi_deta1 mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d01-x02-y01",delta_phi_deta1_mc16);
	if (delta_phi_deta1_mc16 == 0) { cout << "ak5Gen_delta_phi_deta1 mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d01-x02-y01",delta_phi_deta1_mc17);
	if (delta_phi_deta1_mc17 == 0) { cout << "ak5Gen_delta_phi_deta1 mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d01-x02-y01",delta_phi_deta1_mc18);
	if (delta_phi_deta1_mc18 == 0) { cout << "ak5Gen_delta_phi_deta1 mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d01-x02-y01",delta_phi_deta1_mc19);
	if (delta_phi_deta1_mc19 == 0) { cout << "ak5Gen_delta_phi_deta1 mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d01-x02-y01",delta_phi_deta1_mc20);
	if (delta_phi_deta1_mc20 == 0) { cout << "ak5Gen_delta_phi_deta1 mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	extra_label_position = "middle_left";
	}
    if (mode == "ratio")
	{
	extra_label_position = "bottom_right";
	}

    create_final_plot(delta_phi_deta1_data, delta_phi_deta1_data_unc, delta_phi_deta1_mc3, labels[2], delta_phi_deta1_mc4, labels[3], delta_phi_deta1_mc9, labels[8], delta_phi_deta1_mc20, labels[19], delta_phi_deta1_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta1_all", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_data, delta_phi_deta1_data_unc, delta_phi_deta1_mc5, labels[4], delta_phi_deta1_mc8, labels[7], delta_phi_deta1_mc7, labels[6], delta_phi_deta1_mc6, labels[5], delta_phi_deta1_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta1_pythia6", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_data, delta_phi_deta1_data_unc, delta_phi_deta1_mc1, labels[0], delta_phi_deta1_mc2, labels[1], delta_phi_deta1_mc7, labels[6], delta_phi_deta1_mc9, labels[8], delta_phi_deta1_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta1_check", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_data, delta_phi_deta1_data_unc, delta_phi_deta1_mc7, labels[6], delta_phi_deta1_mc6, labels[5], delta_phi_deta1_mc10, labels[9], delta_phi_deta1_mc11, labels[10], delta_phi_deta1_mc12, labels[11], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta1_hadronization", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
   
//final plot for delta phi deta2 distribution
    if (detail) { cout<<"Delta phi deta2"<<endl; }

    TH1D *delta_phi_deta2_data = 0;
    TH1D *delta_phi_deta2_data_unc = 0;
    TH1D *delta_phi_deta2_mc1 = 0;
    TH1D *delta_phi_deta2_mc2 = 0;
    TH1 *delta_phi_deta2_mc3 = 0;
    TH1 *delta_phi_deta2_mc4 = 0;
    TH1 *delta_phi_deta2_mc5 = 0;
    TH1 *delta_phi_deta2_mc6 = 0;
    TH1 *delta_phi_deta2_mc7 = 0;
    TH1 *delta_phi_deta2_mc8 = 0;
    TH1 *delta_phi_deta2_mc9 = 0;
    TH1 *delta_phi_deta2_mc10 = 0;
    TH1 *delta_phi_deta2_mc11 = 0;
    TH1 *delta_phi_deta2_mc12 = 0;
    TH1 *delta_phi_deta2_mc13 = 0;
    TH1 *delta_phi_deta2_mc14 = 0;
    TH1 *delta_phi_deta2_mc15 = 0;
    TH1 *delta_phi_deta2_mc16 = 0;
    TH1 *delta_phi_deta2_mc17 = 0;
    TH1 *delta_phi_deta2_mc18 = 0;
    TH1 *delta_phi_deta2_mc19 = 0;
    TH1 *delta_phi_deta2_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data);
    if (delta_phi_deta2_data == 0) { cout << "ak5PF_delta_phi_deta2 not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta2",delta_phi_deta2_data_unc);
    if (delta_phi_deta2_data_unc == 0) { cout << "ak5PF_delta_phi_deta2 unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta2",delta_phi_deta2_mc1);
	if (delta_phi_deta2_mc1 == 0) { cout << "ak5Gen_delta_phi_deta2 mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta2",delta_phi_deta2_mc2);
	if (delta_phi_deta2_mc2 == 0) { cout << "ak5Gen_delta_phi_deta2 mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d01-x02-y02",delta_phi_deta2_mc3);
	if (delta_phi_deta2_mc3 == 0) { cout << "ak5Gen_delta_phi_deta2 mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d01-x02-y02",delta_phi_deta2_mc4);
	if (delta_phi_deta2_mc4 == 0) { cout << "ak5Gen_delta_phi_deta2 mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d01-x02-y02",delta_phi_deta2_mc5);
	if (delta_phi_deta2_mc5 == 0) { cout << "ak5Gen_delta_phi_deta2 mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d01-x02-y02",delta_phi_deta2_mc6);
	if (delta_phi_deta2_mc6 == 0) { cout << "ak5Gen_delta_phi_deta2 mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d01-x02-y02",delta_phi_deta2_mc7);
	if (delta_phi_deta2_mc7 == 0) { cout << "ak5Gen_delta_phi_deta2 mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d01-x02-y02",delta_phi_deta2_mc8);
	if (delta_phi_deta2_mc8 == 0) { cout << "ak5Gen_delta_phi_deta2 mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d01-x02-y02",delta_phi_deta2_mc9);
	if (delta_phi_deta2_mc9 == 0) { cout << "ak5Gen_delta_phi_deta2 mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d01-x02-y02",delta_phi_deta2_mc10);
	if (delta_phi_deta2_mc10 == 0) { cout << "ak5Gen_delta_phi_deta2 mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d01-x02-y02",delta_phi_deta2_mc11);
	if (delta_phi_deta2_mc11 == 0) { cout << "ak5Gen_delta_phi_deta2 mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d01-x02-y02",delta_phi_deta2_mc12);
	if (delta_phi_deta2_mc12 == 0) { cout << "ak5Gen_delta_phi_deta2 mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d01-x02-y02",delta_phi_deta2_mc13);
	if (delta_phi_deta2_mc13 == 0) { cout << "ak5Gen_delta_phi_deta2 mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d01-x02-y02",delta_phi_deta2_mc14);
	if (delta_phi_deta2_mc14 == 0) { cout << "ak5Gen_delta_phi_deta2 mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d01-x02-y02",delta_phi_deta2_mc15);
	if (delta_phi_deta2_mc15 == 0) { cout << "ak5Gen_delta_phi_deta2 mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d01-x02-y02",delta_phi_deta2_mc16);
	if (delta_phi_deta2_mc16 == 0) { cout << "ak5Gen_delta_phi_deta2 mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d01-x02-y02",delta_phi_deta2_mc17);
	if (delta_phi_deta2_mc17 == 0) { cout << "ak5Gen_delta_phi_deta2 mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d01-x02-y02",delta_phi_deta2_mc18);
	if (delta_phi_deta2_mc18 == 0) { cout << "ak5Gen_delta_phi_deta2 mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d01-x02-y02",delta_phi_deta2_mc19);
	if (delta_phi_deta2_mc19 == 0) { cout << "ak5Gen_delta_phi_deta2 mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d01-x02-y02",delta_phi_deta2_mc20);
	if (delta_phi_deta2_mc20 == 0) { cout << "ak5Gen_delta_phi_deta2 mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_deta2_data, delta_phi_deta2_data_unc, delta_phi_deta2_mc3, labels[2], delta_phi_deta2_mc4, labels[3], delta_phi_deta2_mc9, labels[8], delta_phi_deta2_mc20, labels[19], delta_phi_deta2_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta2_all", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_data, delta_phi_deta2_data_unc, delta_phi_deta2_mc5, labels[4], delta_phi_deta2_mc8, labels[7], delta_phi_deta2_mc7, labels[6], delta_phi_deta2_mc6, labels[5], delta_phi_deta2_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta2_pythia6", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_data, delta_phi_deta2_data_unc, delta_phi_deta2_mc1, labels[0], delta_phi_deta2_mc2, labels[1], delta_phi_deta2_mc7, labels[6], delta_phi_deta2_mc9, labels[8], delta_phi_deta2_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta2_check", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_data, delta_phi_deta2_data_unc, delta_phi_deta2_mc7, labels[6], delta_phi_deta2_mc6, labels[5], delta_phi_deta2_mc10, labels[9], delta_phi_deta2_mc11, labels[10], delta_phi_deta2_mc12, labels[11], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta2_hadronization", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);


//final plot for delta phi deta3 distribution
    if (detail) { cout<<"Delta phi deta3"<<endl; }

    TH1D *delta_phi_deta3_data = 0;
    TH1D *delta_phi_deta3_data_unc = 0;
    TH1D *delta_phi_deta3_mc1 = 0;
    TH1D *delta_phi_deta3_mc2 = 0;
    TH1 *delta_phi_deta3_mc3 = 0;
    TH1 *delta_phi_deta3_mc4 = 0;
    TH1 *delta_phi_deta3_mc5 = 0;
    TH1 *delta_phi_deta3_mc6 = 0;
    TH1 *delta_phi_deta3_mc7 = 0;
    TH1 *delta_phi_deta3_mc8 = 0;
    TH1 *delta_phi_deta3_mc9 = 0;
    TH1 *delta_phi_deta3_mc10 = 0;
    TH1 *delta_phi_deta3_mc11 = 0;
    TH1 *delta_phi_deta3_mc12 = 0;
    TH1 *delta_phi_deta3_mc13 = 0;
    TH1 *delta_phi_deta3_mc14 = 0;
    TH1 *delta_phi_deta3_mc15 = 0;
    TH1 *delta_phi_deta3_mc16 = 0;
    TH1 *delta_phi_deta3_mc17 = 0;
    TH1 *delta_phi_deta3_mc18 = 0;
    TH1 *delta_phi_deta3_mc19 = 0;
    TH1 *delta_phi_deta3_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data);
    if (delta_phi_deta3_data == 0) { cout << "ak5PF_delta_phi_deta3 not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta3",delta_phi_deta3_data_unc);
    if (delta_phi_deta3_data_unc == 0) { cout << "ak5PF_delta_phi_deta3 unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta3",delta_phi_deta3_mc1);
	if (delta_phi_deta3_mc1 == 0) { cout << "ak5Gen_delta_phi_deta3 mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta3",delta_phi_deta3_mc2);
	if (delta_phi_deta3_mc2 == 0) { cout << "ak5Gen_delta_phi_deta3 mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d01-x02-y03",delta_phi_deta3_mc3);
	if (delta_phi_deta3_mc3 == 0) { cout << "ak5Gen_delta_phi_deta3 mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d01-x02-y03",delta_phi_deta3_mc4);
	if (delta_phi_deta3_mc4 == 0) { cout << "ak5Gen_delta_phi_deta3 mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d01-x02-y03",delta_phi_deta3_mc5);
	if (delta_phi_deta3_mc5 == 0) { cout << "ak5Gen_delta_phi_deta3 mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d01-x02-y03",delta_phi_deta3_mc6);
	if (delta_phi_deta3_mc6 == 0) { cout << "ak5Gen_delta_phi_deta3 mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d01-x02-y03",delta_phi_deta3_mc7);
	if (delta_phi_deta3_mc7 == 0) { cout << "ak5Gen_delta_phi_deta3 mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d01-x02-y03",delta_phi_deta3_mc8);
	if (delta_phi_deta3_mc8 == 0) { cout << "ak5Gen_delta_phi_deta3 mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d01-x02-y03",delta_phi_deta3_mc9);
	if (delta_phi_deta3_mc9 == 0) { cout << "ak5Gen_delta_phi_deta3 mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d01-x02-y03",delta_phi_deta3_mc10);
	if (delta_phi_deta3_mc10 == 0) { cout << "ak5Gen_delta_phi_deta3 mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d01-x02-y03",delta_phi_deta3_mc11);
	if (delta_phi_deta3_mc11 == 0) { cout << "ak5Gen_delta_phi_deta3 mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d01-x02-y03",delta_phi_deta3_mc12);
	if (delta_phi_deta3_mc12 == 0) { cout << "ak5Gen_delta_phi_deta3 mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d01-x02-y03",delta_phi_deta3_mc13);
	if (delta_phi_deta3_mc13 == 0) { cout << "ak5Gen_delta_phi_deta3 mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d01-x02-y03",delta_phi_deta3_mc14);
	if (delta_phi_deta3_mc14 == 0) { cout << "ak5Gen_delta_phi_deta3 mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d01-x02-y03",delta_phi_deta3_mc15);
	if (delta_phi_deta3_mc15 == 0) { cout << "ak5Gen_delta_phi_deta3 mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d01-x02-y03",delta_phi_deta3_mc16);
	if (delta_phi_deta3_mc16 == 0) { cout << "ak5Gen_delta_phi_deta3 mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d01-x02-y03",delta_phi_deta3_mc17);
	if (delta_phi_deta3_mc17 == 0) { cout << "ak5Gen_delta_phi_deta3 mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d01-x02-y03",delta_phi_deta3_mc18);
	if (delta_phi_deta3_mc18 == 0) { cout << "ak5Gen_delta_phi_deta3 mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d01-x02-y03",delta_phi_deta3_mc19);
	if (delta_phi_deta3_mc19 == 0) { cout << "ak5Gen_delta_phi_deta3 mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d01-x02-y03",delta_phi_deta3_mc20);
	if (delta_phi_deta3_mc20 == 0) { cout << "ak5Gen_delta_phi_deta3 mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_deta3_data, delta_phi_deta3_data_unc, delta_phi_deta3_mc3, labels[2], delta_phi_deta3_mc4, labels[3], delta_phi_deta3_mc9, labels[8], delta_phi_deta3_mc20, labels[19], delta_phi_deta3_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta3_all", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_data, delta_phi_deta3_data_unc, delta_phi_deta3_mc5, labels[4], delta_phi_deta3_mc8, labels[7], delta_phi_deta3_mc7, labels[6], delta_phi_deta3_mc6, labels[5], delta_phi_deta3_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta3_pythia6", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_data, delta_phi_deta3_data_unc, delta_phi_deta3_mc1, labels[0], delta_phi_deta3_mc2, labels[1], delta_phi_deta3_mc7, labels[6], delta_phi_deta3_mc9, labels[8], delta_phi_deta3_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta3_check", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_data, delta_phi_deta3_data_unc, delta_phi_deta3_mc7, labels[6], delta_phi_deta3_mc6, labels[5], delta_phi_deta3_mc10, labels[9], delta_phi_deta3_mc11, labels[10], delta_phi_deta3_mc12, labels[11], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta3_hadronization", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);


//final plot for delta phi deta4 distribution
    if (detail) { cout<<"Delta phi deta4"<<endl; }

    TH1D *delta_phi_deta4_data = 0;
    TH1D *delta_phi_deta4_data_unc = 0;
    TH1D *delta_phi_deta4_mc1 = 0;
    TH1D *delta_phi_deta4_mc2 = 0;
    TH1 *delta_phi_deta4_mc3 = 0;
    TH1 *delta_phi_deta4_mc4 = 0;
    TH1 *delta_phi_deta4_mc5 = 0;
    TH1 *delta_phi_deta4_mc6 = 0;
    TH1 *delta_phi_deta4_mc7 = 0;
    TH1 *delta_phi_deta4_mc8 = 0;
    TH1 *delta_phi_deta4_mc9 = 0;
    TH1 *delta_phi_deta4_mc10 = 0;
    TH1 *delta_phi_deta4_mc11 = 0;
    TH1 *delta_phi_deta4_mc12 = 0;
    TH1 *delta_phi_deta4_mc13 = 0;
    TH1 *delta_phi_deta4_mc14 = 0;
    TH1 *delta_phi_deta4_mc15 = 0;
    TH1 *delta_phi_deta4_mc16 = 0;
    TH1 *delta_phi_deta4_mc17 = 0;
    TH1 *delta_phi_deta4_mc18 = 0;
    TH1 *delta_phi_deta4_mc19 = 0;
    TH1 *delta_phi_deta4_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data);
    if (delta_phi_deta4_data == 0) { cout << "ak5PF_delta_phi_deta4 not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta4",delta_phi_deta4_data_unc);
    if (delta_phi_deta4_data_unc == 0) { cout << "ak5PF_delta_phi_deta4 unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta4",delta_phi_deta4_mc1);
	if (delta_phi_deta4_mc1 == 0) { cout << "ak5Gen_delta_phi_deta4 mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta4",delta_phi_deta4_mc2);
	if (delta_phi_deta4_mc2 == 0) { cout << "ak5Gen_delta_phi_deta4 mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d01-x02-y04",delta_phi_deta4_mc3);
	if (delta_phi_deta4_mc3 == 0) { cout << "ak5Gen_delta_phi_deta4 mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d01-x02-y04",delta_phi_deta4_mc4);
	if (delta_phi_deta4_mc4 == 0) { cout << "ak5Gen_delta_phi_deta4 mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d01-x02-y04",delta_phi_deta4_mc5);
	if (delta_phi_deta4_mc5 == 0) { cout << "ak5Gen_delta_phi_deta4 mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d01-x02-y04",delta_phi_deta4_mc6);
	if (delta_phi_deta4_mc6 == 0) { cout << "ak5Gen_delta_phi_deta4 mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d01-x02-y04",delta_phi_deta4_mc7);
	if (delta_phi_deta4_mc7 == 0) { cout << "ak5Gen_delta_phi_deta4 mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d01-x02-y04",delta_phi_deta4_mc8);
	if (delta_phi_deta4_mc8 == 0) { cout << "ak5Gen_delta_phi_deta4 mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d01-x02-y04",delta_phi_deta4_mc9);
	if (delta_phi_deta4_mc9 == 0) { cout << "ak5Gen_delta_phi_deta4 mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d01-x02-y04",delta_phi_deta4_mc10);
	if (delta_phi_deta4_mc10 == 0) { cout << "ak5Gen_delta_phi_deta4 mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d01-x02-y04",delta_phi_deta4_mc11);
	if (delta_phi_deta4_mc11 == 0) { cout << "ak5Gen_delta_phi_deta4 mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d01-x02-y04",delta_phi_deta4_mc12);
	if (delta_phi_deta4_mc12 == 0) { cout << "ak5Gen_delta_phi_deta4 mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d01-x02-y04",delta_phi_deta4_mc13);
	if (delta_phi_deta4_mc13 == 0) { cout << "ak5Gen_delta_phi_deta4 mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d01-x02-y04",delta_phi_deta4_mc14);
	if (delta_phi_deta4_mc14 == 0) { cout << "ak5Gen_delta_phi_deta4 mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d01-x02-y04",delta_phi_deta4_mc15);
	if (delta_phi_deta4_mc15 == 0) { cout << "ak5Gen_delta_phi_deta4 mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d01-x02-y04",delta_phi_deta4_mc16);
	if (delta_phi_deta4_mc16 == 0) { cout << "ak5Gen_delta_phi_deta4 mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d01-x02-y04",delta_phi_deta4_mc17);
	if (delta_phi_deta4_mc17 == 0) { cout << "ak5Gen_delta_phi_deta4 mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d01-x02-y04",delta_phi_deta4_mc18);
	if (delta_phi_deta4_mc18 == 0) { cout << "ak5Gen_delta_phi_deta4 mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d01-x02-y04",delta_phi_deta4_mc19);
	if (delta_phi_deta4_mc19 == 0) { cout << "ak5Gen_delta_phi_deta4 mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d01-x02-y04",delta_phi_deta4_mc20);
	if (delta_phi_deta4_mc20 == 0) { cout << "ak5Gen_delta_phi_deta4 mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_deta4_data, delta_phi_deta4_data_unc, delta_phi_deta4_mc3, labels[2], delta_phi_deta4_mc4, labels[3], delta_phi_deta4_mc9, labels[8], delta_phi_deta4_mc20, labels[19], delta_phi_deta4_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta4_all", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_data, delta_phi_deta4_data_unc, delta_phi_deta4_mc5, labels[4], delta_phi_deta4_mc8, labels[7], delta_phi_deta4_mc7, labels[6], delta_phi_deta4_mc6, labels[5], delta_phi_deta4_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta4_pythia6", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_data, delta_phi_deta4_data_unc, delta_phi_deta4_mc1, labels[0], delta_phi_deta4_mc2, labels[1], delta_phi_deta4_mc7, labels[6], delta_phi_deta4_mc9, labels[8], delta_phi_deta4_mc20, labels[19], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta4_check", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_data, delta_phi_deta4_data_unc, delta_phi_deta4_mc7, labels[6], delta_phi_deta4_mc6, labels[5], delta_phi_deta4_mc10, labels[9], delta_phi_deta4_mc11, labels[10], delta_phi_deta4_mc12, labels[11], "INCLUSIVE", output_path_plots, mode, "_final_delta_phi_deta4_hadronization", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);

//final plot for delta phi deta1 nogap distribution
    if (detail) { cout<<"Delta phi deta1 nogap"<<endl; }

    TH1D *delta_phi_deta1_nogap_data = 0;
    TH1D *delta_phi_deta1_nogap_data_unc = 0;
    TH1D *delta_phi_deta1_nogap_mc1 = 0;
    TH1D *delta_phi_deta1_nogap_mc2 = 0;
    TH1 *delta_phi_deta1_nogap_mc3 = 0;
    TH1 *delta_phi_deta1_nogap_mc4 = 0;
    TH1 *delta_phi_deta1_nogap_mc5 = 0;
    TH1 *delta_phi_deta1_nogap_mc6 = 0;
    TH1 *delta_phi_deta1_nogap_mc7 = 0;
    TH1 *delta_phi_deta1_nogap_mc8 = 0;
    TH1 *delta_phi_deta1_nogap_mc9 = 0;
    TH1 *delta_phi_deta1_nogap_mc10 = 0;
    TH1 *delta_phi_deta1_nogap_mc11 = 0;
    TH1 *delta_phi_deta1_nogap_mc12 = 0;
    TH1 *delta_phi_deta1_nogap_mc13 = 0;
    TH1 *delta_phi_deta1_nogap_mc14 = 0;
    TH1 *delta_phi_deta1_nogap_mc15 = 0;
    TH1 *delta_phi_deta1_nogap_mc16 = 0;
    TH1 *delta_phi_deta1_nogap_mc17 = 0;
    TH1 *delta_phi_deta1_nogap_mc18 = 0;
    TH1 *delta_phi_deta1_nogap_mc19 = 0;
    TH1 *delta_phi_deta1_nogap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data);
    if (delta_phi_deta1_nogap_data == 0) { cout << "ak5PF_delta_phi_deta1_nogap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta1_nogap",delta_phi_deta1_nogap_data_unc);
    if (delta_phi_deta1_nogap_data_unc == 0) { cout << "ak5PF_delta_phi_deta1_nogap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta1_nogap",delta_phi_deta1_nogap_mc1);
	if (delta_phi_deta1_nogap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta1_nogap",delta_phi_deta1_nogap_mc2);
	if (delta_phi_deta1_nogap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc3);
	if (delta_phi_deta1_nogap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc4);
	if (delta_phi_deta1_nogap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc5);
	if (delta_phi_deta1_nogap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc6);
	if (delta_phi_deta1_nogap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc7);
	if (delta_phi_deta1_nogap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc8);
	if (delta_phi_deta1_nogap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc9);
	if (delta_phi_deta1_nogap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc10);
	if (delta_phi_deta1_nogap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc11);
	if (delta_phi_deta1_nogap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc12);
	if (delta_phi_deta1_nogap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc13);
	if (delta_phi_deta1_nogap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc14);
	if (delta_phi_deta1_nogap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc15);
	if (delta_phi_deta1_nogap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc16);
	if (delta_phi_deta1_nogap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc17);
	if (delta_phi_deta1_nogap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc18);
	if (delta_phi_deta1_nogap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc19);
	if (delta_phi_deta1_nogap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d02-x02-y01",delta_phi_deta1_nogap_mc20);
	if (delta_phi_deta1_nogap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta1_nogap mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_data_unc, delta_phi_deta1_nogap_mc3, labels[2], delta_phi_deta1_nogap_mc4, labels[3], delta_phi_deta1_nogap_mc9, labels[8], delta_phi_deta1_nogap_mc20, labels[19], delta_phi_deta1_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta1_nogap_all", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_data_unc, delta_phi_deta1_nogap_mc5, labels[4], delta_phi_deta1_nogap_mc8, labels[7], delta_phi_deta1_nogap_mc7, labels[6], delta_phi_deta1_nogap_mc6, labels[5], delta_phi_deta1_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta1_nogap_pythia6", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_data_unc, delta_phi_deta1_nogap_mc1, labels[0], delta_phi_deta1_nogap_mc2, labels[1], delta_phi_deta1_nogap_mc7, labels[6], delta_phi_deta1_nogap_mc9, labels[8], delta_phi_deta1_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta1_nogap_check", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_nogap_data, delta_phi_deta1_nogap_data_unc, delta_phi_deta1_nogap_mc7, labels[6], delta_phi_deta1_nogap_mc6, labels[5], delta_phi_deta1_nogap_mc10, labels[9], delta_phi_deta1_nogap_mc11, labels[10], delta_phi_deta1_nogap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta1_nogap_hadronization", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);

//final plot for delta phi deta2 nogap distribution
    if (detail) { cout<<"Delta phi deta2 nogap"<<endl; }

    TH1D *delta_phi_deta2_nogap_data = 0;
    TH1D *delta_phi_deta2_nogap_data_unc = 0;
    TH1D *delta_phi_deta2_nogap_mc1 = 0;
    TH1D *delta_phi_deta2_nogap_mc2 = 0;
    TH1 *delta_phi_deta2_nogap_mc3 = 0;
    TH1 *delta_phi_deta2_nogap_mc4 = 0;
    TH1 *delta_phi_deta2_nogap_mc5 = 0;
    TH1 *delta_phi_deta2_nogap_mc6 = 0;
    TH1 *delta_phi_deta2_nogap_mc7 = 0;
    TH1 *delta_phi_deta2_nogap_mc8 = 0;
    TH1 *delta_phi_deta2_nogap_mc9 = 0;
    TH1 *delta_phi_deta2_nogap_mc10 = 0;
    TH1 *delta_phi_deta2_nogap_mc11 = 0;
    TH1 *delta_phi_deta2_nogap_mc12 = 0;
    TH1 *delta_phi_deta2_nogap_mc13 = 0;
    TH1 *delta_phi_deta2_nogap_mc14 = 0;
    TH1 *delta_phi_deta2_nogap_mc15 = 0;
    TH1 *delta_phi_deta2_nogap_mc16 = 0;
    TH1 *delta_phi_deta2_nogap_mc17 = 0;
    TH1 *delta_phi_deta2_nogap_mc18 = 0;
    TH1 *delta_phi_deta2_nogap_mc19 = 0;
    TH1 *delta_phi_deta2_nogap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data);
    if (delta_phi_deta2_nogap_data == 0) { cout << "ak5PF_delta_phi_deta2_nogap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta2_nogap",delta_phi_deta2_nogap_data_unc);
    if (delta_phi_deta2_nogap_data_unc == 0) { cout << "ak5PF_delta_phi_deta2_nogap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta2_nogap",delta_phi_deta2_nogap_mc1);
	if (delta_phi_deta2_nogap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta2_nogap",delta_phi_deta2_nogap_mc2);
	if (delta_phi_deta2_nogap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc3);
	if (delta_phi_deta2_nogap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc4);
	if (delta_phi_deta2_nogap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc5);
	if (delta_phi_deta2_nogap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc6);
	if (delta_phi_deta2_nogap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc7);
	if (delta_phi_deta2_nogap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc8);
	if (delta_phi_deta2_nogap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc9);
	if (delta_phi_deta2_nogap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc10);
	if (delta_phi_deta2_nogap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc11);
	if (delta_phi_deta2_nogap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc12);
	if (delta_phi_deta2_nogap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc13);
	if (delta_phi_deta2_nogap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc14);
	if (delta_phi_deta2_nogap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc15);
	if (delta_phi_deta2_nogap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc16);
	if (delta_phi_deta2_nogap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc17);
	if (delta_phi_deta2_nogap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc18);
	if (delta_phi_deta2_nogap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc19);
	if (delta_phi_deta2_nogap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d02-x02-y02",delta_phi_deta2_nogap_mc20);
	if (delta_phi_deta2_nogap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta2_nogap mc20 not found!" << endl; return; }
	}


    if (mode == "xsec")
	{
	extra_label_position = "bottom_left";
	}

    create_final_plot(delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_data_unc, delta_phi_deta2_nogap_mc3, labels[2], delta_phi_deta2_nogap_mc4, labels[3], delta_phi_deta2_nogap_mc9, labels[8], delta_phi_deta2_nogap_mc20, labels[19], delta_phi_deta2_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta2_nogap_all", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_data_unc, delta_phi_deta2_nogap_mc5, labels[4], delta_phi_deta2_nogap_mc8, labels[7], delta_phi_deta2_nogap_mc7, labels[6], delta_phi_deta2_nogap_mc6, labels[5], delta_phi_deta2_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta2_nogap_pythia6", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_data_unc, delta_phi_deta2_nogap_mc1, labels[0], delta_phi_deta2_nogap_mc2, labels[1], delta_phi_deta2_nogap_mc7, labels[6], delta_phi_deta2_nogap_mc9, labels[8], delta_phi_deta2_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta2_nogap_check", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_nogap_data, delta_phi_deta2_nogap_data_unc, delta_phi_deta2_nogap_mc7, labels[6], delta_phi_deta2_nogap_mc6, labels[5], delta_phi_deta2_nogap_mc10, labels[9], delta_phi_deta2_nogap_mc11, labels[10], delta_phi_deta2_nogap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta2_nogap_hadronization", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);

//final plot for delta phi deta3 nogap distribution
    if (detail) { cout<<"Delta phi deta3 nogap"<<endl; }

    TH1D *delta_phi_deta3_nogap_data = 0;
    TH1D *delta_phi_deta3_nogap_data_unc = 0;
    TH1D *delta_phi_deta3_nogap_mc1 = 0;
    TH1D *delta_phi_deta3_nogap_mc2 = 0;
    TH1 *delta_phi_deta3_nogap_mc3 = 0;
    TH1 *delta_phi_deta3_nogap_mc4 = 0;
    TH1 *delta_phi_deta3_nogap_mc5 = 0;
    TH1 *delta_phi_deta3_nogap_mc6 = 0;
    TH1 *delta_phi_deta3_nogap_mc7 = 0;
    TH1 *delta_phi_deta3_nogap_mc8 = 0;
    TH1 *delta_phi_deta3_nogap_mc9 = 0;
    TH1 *delta_phi_deta3_nogap_mc10 = 0;
    TH1 *delta_phi_deta3_nogap_mc11 = 0;
    TH1 *delta_phi_deta3_nogap_mc12 = 0;
    TH1 *delta_phi_deta3_nogap_mc13 = 0;
    TH1 *delta_phi_deta3_nogap_mc14 = 0;
    TH1 *delta_phi_deta3_nogap_mc15 = 0;
    TH1 *delta_phi_deta3_nogap_mc16 = 0;
    TH1 *delta_phi_deta3_nogap_mc17 = 0;
    TH1 *delta_phi_deta3_nogap_mc18 = 0;
    TH1 *delta_phi_deta3_nogap_mc19 = 0;
    TH1 *delta_phi_deta3_nogap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data);
    if (delta_phi_deta3_nogap_data == 0) { cout << "ak5PF_delta_phi_deta3_nogap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta3_nogap",delta_phi_deta3_nogap_data_unc);
    if (delta_phi_deta3_nogap_data_unc == 0) { cout << "ak5PF_delta_phi_deta3_nogap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta3_nogap",delta_phi_deta3_nogap_mc1);
	if (delta_phi_deta3_nogap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta3_nogap",delta_phi_deta3_nogap_mc2);
	if (delta_phi_deta3_nogap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc3);
	if (delta_phi_deta3_nogap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc4);
	if (delta_phi_deta3_nogap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc5);
	if (delta_phi_deta3_nogap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc6);
	if (delta_phi_deta3_nogap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc7);
	if (delta_phi_deta3_nogap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc8);
	if (delta_phi_deta3_nogap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc9);
	if (delta_phi_deta3_nogap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc10);
	if (delta_phi_deta3_nogap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc11);
	if (delta_phi_deta3_nogap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc12);
	if (delta_phi_deta3_nogap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc13);
	if (delta_phi_deta3_nogap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc14);
	if (delta_phi_deta3_nogap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc15);
	if (delta_phi_deta3_nogap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc16);
	if (delta_phi_deta3_nogap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc17);
	if (delta_phi_deta3_nogap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc18);
	if (delta_phi_deta3_nogap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc19);
	if (delta_phi_deta3_nogap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d02-x02-y03",delta_phi_deta3_nogap_mc20);
	if (delta_phi_deta3_nogap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta3_nogap mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_data_unc, delta_phi_deta3_nogap_mc3, labels[2], delta_phi_deta3_nogap_mc4, labels[3], delta_phi_deta3_nogap_mc9, labels[8], delta_phi_deta3_nogap_mc20, labels[19], delta_phi_deta3_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta3_nogap_all", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_data_unc, delta_phi_deta3_nogap_mc5, labels[4], delta_phi_deta3_nogap_mc8, labels[7], delta_phi_deta3_nogap_mc7, labels[6], delta_phi_deta3_nogap_mc6, labels[5], delta_phi_deta3_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta3_nogap_pythia6", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_data_unc, delta_phi_deta3_nogap_mc1, labels[0], delta_phi_deta3_nogap_mc2, labels[1], delta_phi_deta3_nogap_mc7, labels[6], delta_phi_deta3_nogap_mc9, labels[8], delta_phi_deta3_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta3_nogap_check", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_nogap_data, delta_phi_deta3_nogap_data_unc, delta_phi_deta3_nogap_mc7, labels[6], delta_phi_deta3_nogap_mc6, labels[5], delta_phi_deta3_nogap_mc10, labels[9], delta_phi_deta3_nogap_mc11, labels[10], delta_phi_deta3_nogap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta3_nogap_hadronization", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);

//final plot for delta phi deta4 nogap distribution
    if (detail) { cout<<"Delta phi deta4 nogap"<<endl; }

    TH1D *delta_phi_deta4_nogap_data = 0;
    TH1D *delta_phi_deta4_nogap_data_unc = 0;
    TH1D *delta_phi_deta4_nogap_mc1 = 0;
    TH1D *delta_phi_deta4_nogap_mc2 = 0;
    TH1 *delta_phi_deta4_nogap_mc3 = 0;
    TH1 *delta_phi_deta4_nogap_mc4 = 0;
    TH1 *delta_phi_deta4_nogap_mc5 = 0;
    TH1 *delta_phi_deta4_nogap_mc6 = 0;
    TH1 *delta_phi_deta4_nogap_mc7 = 0;
    TH1 *delta_phi_deta4_nogap_mc8 = 0;
    TH1 *delta_phi_deta4_nogap_mc9 = 0;
    TH1 *delta_phi_deta4_nogap_mc10 = 0;
    TH1 *delta_phi_deta4_nogap_mc11 = 0;
    TH1 *delta_phi_deta4_nogap_mc12 = 0;
    TH1 *delta_phi_deta4_nogap_mc13 = 0;
    TH1 *delta_phi_deta4_nogap_mc14 = 0;
    TH1 *delta_phi_deta4_nogap_mc15 = 0;
    TH1 *delta_phi_deta4_nogap_mc16 = 0;
    TH1 *delta_phi_deta4_nogap_mc17 = 0;
    TH1 *delta_phi_deta4_nogap_mc18 = 0;
    TH1 *delta_phi_deta4_nogap_mc19 = 0;
    TH1 *delta_phi_deta4_nogap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data);
    if (delta_phi_deta4_nogap_data == 0) { cout << "ak5PF_delta_phi_deta4_nogap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta4_nogap",delta_phi_deta4_nogap_data_unc);
    if (delta_phi_deta4_nogap_data_unc == 0) { cout << "ak5PF_delta_phi_deta4_nogap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta4_nogap",delta_phi_deta4_nogap_mc1);
	if (delta_phi_deta4_nogap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta4_nogap",delta_phi_deta4_nogap_mc2);
	if (delta_phi_deta4_nogap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc3);
	if (delta_phi_deta4_nogap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc4);
	if (delta_phi_deta4_nogap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc5);
	if (delta_phi_deta4_nogap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc6);
	if (delta_phi_deta4_nogap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc7);
	if (delta_phi_deta4_nogap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc8);
	if (delta_phi_deta4_nogap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc9);
	if (delta_phi_deta4_nogap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc10);
	if (delta_phi_deta4_nogap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc11);
	if (delta_phi_deta4_nogap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc12);
	if (delta_phi_deta4_nogap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc13);
	if (delta_phi_deta4_nogap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc14);
	if (delta_phi_deta4_nogap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc15);
	if (delta_phi_deta4_nogap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc16);
	if (delta_phi_deta4_nogap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc17);
	if (delta_phi_deta4_nogap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc18);
	if (delta_phi_deta4_nogap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc19);
	if (delta_phi_deta4_nogap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d02-x02-y04",delta_phi_deta4_nogap_mc20);
	if (delta_phi_deta4_nogap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta4_nogap mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	extra_label_position = "middle_left";
	}

    create_final_plot(delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_data_unc, delta_phi_deta4_nogap_mc3, labels[2], delta_phi_deta4_nogap_mc4, labels[3], delta_phi_deta4_nogap_mc9, labels[8], delta_phi_deta4_nogap_mc20, labels[19], delta_phi_deta4_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta4_nogap_all", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_data_unc, delta_phi_deta4_nogap_mc5, labels[4], delta_phi_deta4_nogap_mc8, labels[7], delta_phi_deta4_nogap_mc7, labels[6], delta_phi_deta4_nogap_mc6, labels[5], delta_phi_deta4_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta4_nogap_pythia6", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_data_unc, delta_phi_deta4_nogap_mc1, labels[0], delta_phi_deta4_nogap_mc2, labels[1], delta_phi_deta4_nogap_mc7, labels[6], delta_phi_deta4_nogap_mc9, labels[8], delta_phi_deta4_nogap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta4_nogap_check", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_nogap_data, delta_phi_deta4_nogap_data_unc, delta_phi_deta4_nogap_mc7, labels[6], delta_phi_deta4_nogap_mc6, labels[5], delta_phi_deta4_nogap_mc10, labels[9], delta_phi_deta4_nogap_mc11, labels[10], delta_phi_deta4_nogap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_delta_phi_deta4_nogap_hadronization", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);

//final plot for delta phi deta1 gap distribution
    if (detail) { cout<<"Delta phi deta1 gap"<<endl; }

    TH1D *delta_phi_deta1_gap_data = 0;
    TH1D *delta_phi_deta1_gap_data_unc = 0;
    TH1D *delta_phi_deta1_gap_mc1 = 0;
    TH1D *delta_phi_deta1_gap_mc2 = 0;
    TH1 *delta_phi_deta1_gap_mc3 = 0;
    TH1 *delta_phi_deta1_gap_mc4 = 0;
    TH1 *delta_phi_deta1_gap_mc5 = 0;
    TH1 *delta_phi_deta1_gap_mc6 = 0;
    TH1 *delta_phi_deta1_gap_mc7 = 0;
    TH1 *delta_phi_deta1_gap_mc8 = 0;
    TH1 *delta_phi_deta1_gap_mc9 = 0;
    TH1 *delta_phi_deta1_gap_mc10 = 0;
    TH1 *delta_phi_deta1_gap_mc11 = 0;
    TH1 *delta_phi_deta1_gap_mc12 = 0;
    TH1 *delta_phi_deta1_gap_mc13 = 0;
    TH1 *delta_phi_deta1_gap_mc14 = 0;
    TH1 *delta_phi_deta1_gap_mc15 = 0;
    TH1 *delta_phi_deta1_gap_mc16 = 0;
    TH1 *delta_phi_deta1_gap_mc17 = 0;
    TH1 *delta_phi_deta1_gap_mc18 = 0;
    TH1 *delta_phi_deta1_gap_mc19 = 0;
    TH1 *delta_phi_deta1_gap_mc20 = 0;


    data->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data);
    if (delta_phi_deta1_gap_data == 0) { cout << "ak5PF_delta_phi_deta1_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta1_gap",delta_phi_deta1_gap_data_unc);
    if (delta_phi_deta1_gap_data_unc == 0) { cout << "ak5PF_delta_phi_deta1_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta1_gap",delta_phi_deta1_gap_mc1);
	if (delta_phi_deta1_gap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta1_gap",delta_phi_deta1_gap_mc2);
	if (delta_phi_deta1_gap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc3);
	if (delta_phi_deta1_gap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc4);
	if (delta_phi_deta1_gap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc5);
	if (delta_phi_deta1_gap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc6);
	if (delta_phi_deta1_gap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc7);
	if (delta_phi_deta1_gap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc8);
	if (delta_phi_deta1_gap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc9);
	if (delta_phi_deta1_gap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc10);
	if (delta_phi_deta1_gap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc11);
	if (delta_phi_deta1_gap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc12);
	if (delta_phi_deta1_gap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc13);
	if (delta_phi_deta1_gap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc14);
	if (delta_phi_deta1_gap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc15);
	if (delta_phi_deta1_gap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc16);
	if (delta_phi_deta1_gap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc17);
	if (delta_phi_deta1_gap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc18);
	if (delta_phi_deta1_gap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc19);
	if (delta_phi_deta1_gap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d03-x02-y01",delta_phi_deta1_gap_mc20);
	if (delta_phi_deta1_gap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta1_gap mc20 not found!" << endl; return; }
	}


    create_final_plot(delta_phi_deta1_gap_data, delta_phi_deta1_gap_data_unc, delta_phi_deta1_gap_mc3, labels[2], delta_phi_deta1_gap_mc4, labels[3], delta_phi_deta1_gap_mc9, labels[8], delta_phi_deta1_gap_mc20, labels[19], delta_phi_deta1_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta1_gap_all", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_gap_data, delta_phi_deta1_gap_data_unc, delta_phi_deta1_gap_mc5, labels[4], delta_phi_deta1_gap_mc8, labels[7], delta_phi_deta1_gap_mc7, labels[6], delta_phi_deta1_gap_mc6, labels[5], delta_phi_deta1_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta1_gap_pythia6", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_gap_data, delta_phi_deta1_gap_data_unc, delta_phi_deta1_gap_mc1, labels[0], delta_phi_deta1_gap_mc2, labels[1], delta_phi_deta1_gap_mc7, labels[6], delta_phi_deta1_gap_mc9, labels[8], delta_phi_deta1_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta1_gap_check", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);
    create_final_plot(delta_phi_deta1_gap_data, delta_phi_deta1_gap_data_unc, delta_phi_deta1_gap_mc7, labels[6], delta_phi_deta1_gap_mc6, labels[5], delta_phi_deta1_gap_mc10, labels[9], delta_phi_deta1_gap_mc11, labels[10], delta_phi_deta1_gap_mc12, labels[11], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta1_gap_hadronization", legend_position, label_position, "0.4 < #Delta#eta < 2.5", extra_label_position);

//final plot for delta phi deta2 gap distribution
    if (detail) { cout<<"Delta phi deta2 gap"<<endl; }

    TH1D *delta_phi_deta2_gap_data = 0;
    TH1D *delta_phi_deta2_gap_data_unc = 0;
    TH1D *delta_phi_deta2_gap_mc1 = 0;
    TH1D *delta_phi_deta2_gap_mc2 = 0;
    TH1 *delta_phi_deta2_gap_mc3 = 0;
    TH1 *delta_phi_deta2_gap_mc4 = 0;
    TH1 *delta_phi_deta2_gap_mc5 = 0;
    TH1 *delta_phi_deta2_gap_mc6 = 0;
    TH1 *delta_phi_deta2_gap_mc7 = 0;
    TH1 *delta_phi_deta2_gap_mc8 = 0;
    TH1 *delta_phi_deta2_gap_mc9 = 0;
    TH1 *delta_phi_deta2_gap_mc10 = 0;
    TH1 *delta_phi_deta2_gap_mc11 = 0;
    TH1 *delta_phi_deta2_gap_mc12 = 0;
    TH1 *delta_phi_deta2_gap_mc13 = 0;
    TH1 *delta_phi_deta2_gap_mc14 = 0;
    TH1 *delta_phi_deta2_gap_mc15 = 0;
    TH1 *delta_phi_deta2_gap_mc16 = 0;
    TH1 *delta_phi_deta2_gap_mc17 = 0;
    TH1 *delta_phi_deta2_gap_mc18 = 0;
    TH1 *delta_phi_deta2_gap_mc19 = 0;
    TH1 *delta_phi_deta2_gap_mc20 = 0;


    data->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data);
    if (delta_phi_deta2_gap_data == 0) { cout << "ak5PF_delta_phi_deta2_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta2_gap",delta_phi_deta2_gap_data_unc);
    if (delta_phi_deta2_gap_data_unc == 0) { cout << "ak5PF_delta_phi_deta2_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta2_gap",delta_phi_deta2_gap_mc1);
	if (delta_phi_deta2_gap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta2_gap",delta_phi_deta2_gap_mc2);
	if (delta_phi_deta2_gap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc3);
	if (delta_phi_deta2_gap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc4);
	if (delta_phi_deta2_gap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc5);
	if (delta_phi_deta2_gap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc6);
	if (delta_phi_deta2_gap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc7);
	if (delta_phi_deta2_gap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc8);
	if (delta_phi_deta2_gap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc9);
	if (delta_phi_deta2_gap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc10);
	if (delta_phi_deta2_gap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc11);
	if (delta_phi_deta2_gap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc12);
	if (delta_phi_deta2_gap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc13);
	if (delta_phi_deta2_gap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc14);
	if (delta_phi_deta2_gap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc15);
	if (delta_phi_deta2_gap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc16);
	if (delta_phi_deta2_gap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc17);
	if (delta_phi_deta2_gap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc18);
	if (delta_phi_deta2_gap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc19);
	if (delta_phi_deta2_gap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d03-x02-y02",delta_phi_deta2_gap_mc20);
	if (delta_phi_deta2_gap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta2_gap mc20 not found!" << endl; return; }
	}

    create_final_plot(delta_phi_deta2_gap_data, delta_phi_deta2_gap_data_unc, delta_phi_deta2_gap_mc3, labels[2], delta_phi_deta2_gap_mc4, labels[3], delta_phi_deta2_gap_mc9, labels[8], delta_phi_deta2_gap_mc20, labels[19], delta_phi_deta2_gap_mc10, labels[9], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta2_gap_all", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_gap_data, delta_phi_deta2_gap_data_unc, delta_phi_deta2_gap_mc5, labels[4], delta_phi_deta2_gap_mc8, labels[7], delta_phi_deta2_gap_mc7, labels[6], delta_phi_deta2_gap_mc6, labels[5], delta_phi_deta2_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta2_gap_pythia6", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_gap_data, delta_phi_deta2_gap_data_unc, delta_phi_deta2_gap_mc1, labels[0], delta_phi_deta2_gap_mc2, labels[1], delta_phi_deta2_gap_mc7, labels[6], delta_phi_deta2_gap_mc9, labels[8], delta_phi_deta2_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta2_gap_check", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);
    create_final_plot(delta_phi_deta2_gap_data, delta_phi_deta2_gap_data_unc, delta_phi_deta2_gap_mc7, labels[6], delta_phi_deta2_gap_mc6, labels[5], delta_phi_deta2_gap_mc10, labels[9], delta_phi_deta2_gap_mc11, labels[10], delta_phi_deta2_gap_mc12, labels[11], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta2_gap_hadronization", legend_position, label_position, "2.5 < #Delta#eta < 3.5", extra_label_position);

//final plot for delta phi deta3 gap distribution
    if (detail) { cout<<"Delta phi deta3 gap"<<endl; }

    TH1D *delta_phi_deta3_gap_data = 0;
    TH1D *delta_phi_deta3_gap_data_unc = 0;
    TH1D *delta_phi_deta3_gap_mc1 = 0;
    TH1D *delta_phi_deta3_gap_mc2 = 0;
    TH1 *delta_phi_deta3_gap_mc3 = 0;
    TH1 *delta_phi_deta3_gap_mc4 = 0;
    TH1 *delta_phi_deta3_gap_mc5 = 0;
    TH1 *delta_phi_deta3_gap_mc6 = 0;
    TH1 *delta_phi_deta3_gap_mc7 = 0;
    TH1 *delta_phi_deta3_gap_mc8 = 0;
    TH1 *delta_phi_deta3_gap_mc9 = 0;
    TH1 *delta_phi_deta3_gap_mc10 = 0;
    TH1 *delta_phi_deta3_gap_mc11 = 0;
    TH1 *delta_phi_deta3_gap_mc12 = 0;
    TH1 *delta_phi_deta3_gap_mc13 = 0;
    TH1 *delta_phi_deta3_gap_mc14 = 0;
    TH1 *delta_phi_deta3_gap_mc15 = 0;
    TH1 *delta_phi_deta3_gap_mc16 = 0;
    TH1 *delta_phi_deta3_gap_mc17 = 0;
    TH1 *delta_phi_deta3_gap_mc18 = 0;
    TH1 *delta_phi_deta3_gap_mc19 = 0;
    TH1 *delta_phi_deta3_gap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data);
    if (delta_phi_deta3_gap_data == 0) { cout << "ak5PF_delta_phi_deta3_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta3_gap",delta_phi_deta3_gap_data_unc);
    if (delta_phi_deta3_gap_data_unc == 0) { cout << "ak5PF_delta_phi_deta3_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta3_gap",delta_phi_deta3_gap_mc1);
	if (delta_phi_deta3_gap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta3_gap",delta_phi_deta3_gap_mc2);
	if (delta_phi_deta3_gap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc3);
	if (delta_phi_deta3_gap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc4);
	if (delta_phi_deta3_gap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc5);
	if (delta_phi_deta3_gap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc6);
	if (delta_phi_deta3_gap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc7);
	if (delta_phi_deta3_gap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc8);
	if (delta_phi_deta3_gap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc9);
	if (delta_phi_deta3_gap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc10);
	if (delta_phi_deta3_gap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc11);
	if (delta_phi_deta3_gap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc12);
	if (delta_phi_deta3_gap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc13);
	if (delta_phi_deta3_gap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc14);
	if (delta_phi_deta3_gap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc15);
	if (delta_phi_deta3_gap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc16);
	if (delta_phi_deta3_gap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc17);
	if (delta_phi_deta3_gap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc18);
	if (delta_phi_deta3_gap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc19);
	if (delta_phi_deta3_gap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d03-x02-y03",delta_phi_deta3_gap_mc20);
	if (delta_phi_deta3_gap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta3_gap mc20 not found!" << endl; return; }
	}


    create_final_plot(delta_phi_deta3_gap_data, delta_phi_deta3_gap_data_unc, delta_phi_deta3_gap_mc3, labels[2], delta_phi_deta3_gap_mc4, labels[3], delta_phi_deta3_gap_mc9, labels[8], delta_phi_deta3_gap_mc20, labels[19], delta_phi_deta3_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta3_gap_all", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_gap_data, delta_phi_deta3_gap_data_unc, delta_phi_deta3_gap_mc5, labels[4], delta_phi_deta3_gap_mc8, labels[7], delta_phi_deta3_gap_mc7, labels[6], delta_phi_deta3_gap_mc6, labels[5], delta_phi_deta3_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta3_gap_pythia6", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_gap_data, delta_phi_deta3_gap_data_unc, delta_phi_deta3_gap_mc1, labels[0], delta_phi_deta3_gap_mc2, labels[1], delta_phi_deta3_gap_mc7, labels[6], delta_phi_deta3_gap_mc9, labels[8], delta_phi_deta3_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta3_gap_check", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    create_final_plot(delta_phi_deta3_gap_data, delta_phi_deta3_gap_data_unc, delta_phi_deta3_gap_mc7, labels[6], delta_phi_deta3_gap_mc6, labels[5], delta_phi_deta3_gap_mc10, labels[9], delta_phi_deta3_gap_mc11, labels[10], delta_phi_deta3_gap_mc12, labels[11], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta3_gap_hadronization", legend_position, label_position, "3.5 < #Delta#eta < 4.5", extra_label_position);
    
//final plot for delta phi deta4 gap distribution
    if (detail) { cout<<"Delta phi deta4 gap"<<endl; }

    TH1D *delta_phi_deta4_gap_data = 0;
    TH1D *delta_phi_deta4_gap_data_unc = 0;
    TH1D *delta_phi_deta4_gap_mc1 = 0;
    TH1D *delta_phi_deta4_gap_mc2 = 0;
    TH1 *delta_phi_deta4_gap_mc3 = 0;
    TH1 *delta_phi_deta4_gap_mc4 = 0;
    TH1 *delta_phi_deta4_gap_mc5 = 0;
    TH1 *delta_phi_deta4_gap_mc6 = 0;
    TH1 *delta_phi_deta4_gap_mc7 = 0;
    TH1 *delta_phi_deta4_gap_mc8 = 0;
    TH1 *delta_phi_deta4_gap_mc9 = 0;
    TH1 *delta_phi_deta4_gap_mc10 = 0;
    TH1 *delta_phi_deta4_gap_mc11 = 0;
    TH1 *delta_phi_deta4_gap_mc12 = 0;
    TH1 *delta_phi_deta4_gap_mc13 = 0;
    TH1 *delta_phi_deta4_gap_mc14 = 0;
    TH1 *delta_phi_deta4_gap_mc15 = 0;
    TH1 *delta_phi_deta4_gap_mc16 = 0;
    TH1 *delta_phi_deta4_gap_mc17 = 0;
    TH1 *delta_phi_deta4_gap_mc18 = 0;
    TH1 *delta_phi_deta4_gap_mc19 = 0;
    TH1 *delta_phi_deta4_gap_mc20 = 0;

    data->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data);
    if (delta_phi_deta4_gap_data == 0) { cout << "ak5PF_delta_phi_deta4_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_phi_deta4_gap",delta_phi_deta4_gap_data_unc);
    if (delta_phi_deta4_gap_data_unc == 0) { cout << "ak5PF_delta_phi_deta4_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_phi_deta4_gap",delta_phi_deta4_gap_mc1);
	if (delta_phi_deta4_gap_mc1 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_phi_deta4_gap",delta_phi_deta4_gap_mc2);
	if (delta_phi_deta4_gap_mc2 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc3);
	if (delta_phi_deta4_gap_mc3 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc4);
	if (delta_phi_deta4_gap_mc4 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc5);
	if (delta_phi_deta4_gap_mc5 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc6);
	if (delta_phi_deta4_gap_mc6 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc7);
	if (delta_phi_deta4_gap_mc7 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc8);
	if (delta_phi_deta4_gap_mc8 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc9);
	if (delta_phi_deta4_gap_mc9 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc10);
	if (delta_phi_deta4_gap_mc10 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc11);
	if (delta_phi_deta4_gap_mc11 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc12);
	if (delta_phi_deta4_gap_mc12 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc13);
	if (delta_phi_deta4_gap_mc13 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc14);
	if (delta_phi_deta4_gap_mc14 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc15);
	if (delta_phi_deta4_gap_mc15 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc16);
	if (delta_phi_deta4_gap_mc16 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc17);
	if (delta_phi_deta4_gap_mc17 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc18);
	if (delta_phi_deta4_gap_mc18 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc19);
	if (delta_phi_deta4_gap_mc19 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d03-x02-y04",delta_phi_deta4_gap_mc20);
	if (delta_phi_deta4_gap_mc20 == 0) { cout << "ak5Gen_delta_phi_deta4_gap mc20 not found!" << endl; return; }
	}


    create_final_plot(delta_phi_deta4_gap_data, delta_phi_deta4_gap_data_unc, delta_phi_deta4_gap_mc3, labels[2], delta_phi_deta4_gap_mc4, labels[3], delta_phi_deta4_gap_mc9, labels[8], delta_phi_deta4_gap_mc20, labels[19], delta_phi_deta4_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta4_gap_all", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_gap_data, delta_phi_deta4_gap_data_unc, delta_phi_deta4_gap_mc5, labels[4], delta_phi_deta4_gap_mc8, labels[7], delta_phi_deta4_gap_mc7, labels[6], delta_phi_deta4_gap_mc6, labels[5], delta_phi_deta4_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta4_gap_pythia6", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_gap_data, delta_phi_deta4_gap_data_unc, delta_phi_deta4_gap_mc1, labels[0], delta_phi_deta4_gap_mc2, labels[1], delta_phi_deta4_gap_mc7, labels[6], delta_phi_deta4_gap_mc9, labels[8], delta_phi_deta4_gap_mc20, labels[19], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta4_gap_check", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);
    create_final_plot(delta_phi_deta4_gap_data, delta_phi_deta4_gap_data_unc, delta_phi_deta4_gap_mc7, labels[6], delta_phi_deta4_gap_mc6, labels[5], delta_phi_deta4_gap_mc10, labels[9], delta_phi_deta4_gap_mc11, labels[10], delta_phi_deta4_gap_mc12, labels[11], "INSIDE-JET VETO", output_path_plots, mode, "_final_delta_phi_deta4_gap_hadronization", legend_position, label_position, "4.5 < #Delta#eta < 7.5", extra_label_position);

//final plot for leading pt inside gap distribution
    if (detail) { cout<<"Leading pT Inside Gap"<<endl; }

    TH1D *leading_pt_inside_gap_data = 0;
    TH1D *leading_pt_inside_gap_data_unc = 0;
    TH1D *leading_pt_inside_gap_mc1 = 0;
    TH1D *leading_pt_inside_gap_mc2 = 0;
    TH1 *leading_pt_inside_gap_mc3 = 0;
    TH1 *leading_pt_inside_gap_mc4 = 0;
    TH1 *leading_pt_inside_gap_mc5 = 0;
    TH1 *leading_pt_inside_gap_mc6 = 0;
    TH1 *leading_pt_inside_gap_mc7 = 0;
    TH1 *leading_pt_inside_gap_mc8 = 0;
    TH1 *leading_pt_inside_gap_mc9 = 0;
    TH1 *leading_pt_inside_gap_mc10 = 0;
    TH1 *leading_pt_inside_gap_mc11 = 0;
    TH1 *leading_pt_inside_gap_mc12 = 0;
    TH1 *leading_pt_inside_gap_mc13 = 0;
    TH1 *leading_pt_inside_gap_mc14 = 0;
    TH1 *leading_pt_inside_gap_mc15 = 0;
    TH1 *leading_pt_inside_gap_mc16 = 0;
    TH1 *leading_pt_inside_gap_mc17 = 0;
    TH1 *leading_pt_inside_gap_mc18 = 0;
    TH1 *leading_pt_inside_gap_mc19 = 0;
    TH1 *leading_pt_inside_gap_mc20 = 0;

    data->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data);
    if (leading_pt_inside_gap_data == 0) { cout << "ak5PF_leading_pt_inside_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_leading_pt_inside_gap",leading_pt_inside_gap_data_unc);
    if (leading_pt_inside_gap_data_unc == 0) { cout << "ak5PF_leading_pt_inside_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_leading_pt_inside_gap",leading_pt_inside_gap_mc1);
	if (leading_pt_inside_gap_mc1 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_leading_pt_inside_gap",leading_pt_inside_gap_mc2);
	if (leading_pt_inside_gap_mc2 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d04-x01-y01",leading_pt_inside_gap_mc3);
	if (leading_pt_inside_gap_mc3 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d04-x01-y01",leading_pt_inside_gap_mc4);
	if (leading_pt_inside_gap_mc4 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d04-x01-y01",leading_pt_inside_gap_mc5);
	if (leading_pt_inside_gap_mc5 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d04-x01-y01",leading_pt_inside_gap_mc6);
	if (leading_pt_inside_gap_mc6 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d04-x01-y01",leading_pt_inside_gap_mc7);
	if (leading_pt_inside_gap_mc7 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d04-x01-y01",leading_pt_inside_gap_mc8);
	if (leading_pt_inside_gap_mc8 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d04-x01-y01",leading_pt_inside_gap_mc9);
	if (leading_pt_inside_gap_mc9 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d04-x01-y01",leading_pt_inside_gap_mc10);
	if (leading_pt_inside_gap_mc10 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d04-x01-y01",leading_pt_inside_gap_mc11);
	if (leading_pt_inside_gap_mc11 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d04-x01-y01",leading_pt_inside_gap_mc12);
	if (leading_pt_inside_gap_mc12 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d04-x01-y01",leading_pt_inside_gap_mc13);
	if (leading_pt_inside_gap_mc13 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d04-x01-y01",leading_pt_inside_gap_mc14);
	if (leading_pt_inside_gap_mc14 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d04-x01-y01",leading_pt_inside_gap_mc15);
	if (leading_pt_inside_gap_mc15 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d04-x01-y01",leading_pt_inside_gap_mc16);
	if (leading_pt_inside_gap_mc16 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d04-x01-y01",leading_pt_inside_gap_mc17);
	if (leading_pt_inside_gap_mc17 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d04-x01-y01",leading_pt_inside_gap_mc18);
	if (leading_pt_inside_gap_mc18 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d04-x01-y01",leading_pt_inside_gap_mc19);
	if (leading_pt_inside_gap_mc19 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d04-x01-y01",leading_pt_inside_gap_mc20);
	if (leading_pt_inside_gap_mc20 == 0) { cout << "ak5Gen_leading_pt_inside_gap mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	legend_position = "top_right";
	label_position = "bottom_left";
	}

    create_final_plot(leading_pt_inside_gap_data, leading_pt_inside_gap_data_unc, leading_pt_inside_gap_mc3, labels[2], leading_pt_inside_gap_mc4, labels[3], leading_pt_inside_gap_mc9, labels[8], leading_pt_inside_gap_mc20, labels[19], leading_pt_inside_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_inside_gap_all", legend_position, label_position);
    create_final_plot(leading_pt_inside_gap_data, leading_pt_inside_gap_data_unc, leading_pt_inside_gap_mc5, labels[4], leading_pt_inside_gap_mc8, labels[7], leading_pt_inside_gap_mc7, labels[6], leading_pt_inside_gap_mc6, labels[5], leading_pt_inside_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_inside_gap_pythia6", legend_position, label_position);
    create_final_plot(leading_pt_inside_gap_data, leading_pt_inside_gap_data_unc, leading_pt_inside_gap_mc1, labels[0], leading_pt_inside_gap_mc2, labels[1], leading_pt_inside_gap_mc7, labels[6], leading_pt_inside_gap_mc9, labels[8], leading_pt_inside_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_inside_gap_check", legend_position, label_position);
    create_final_plot(leading_pt_inside_gap_data, leading_pt_inside_gap_data_unc, leading_pt_inside_gap_mc7, labels[6], leading_pt_inside_gap_mc6, labels[5], leading_pt_inside_gap_mc10, labels[9], leading_pt_inside_gap_mc11, labels[10], leading_pt_inside_gap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_inside_gap_hadronization", legend_position, label_position);

//final plot for leading Eta* inside gap distribution
    if (detail) { cout<<"Leading Eta* Inside Gap"<<endl; }

    TH1D *leading_eta_star_inside_gap_data = 0;
    TH1D *leading_eta_star_inside_gap_data_unc = 0;
    TH1D *leading_eta_star_inside_gap_mc1 = 0;
    TH1D *leading_eta_star_inside_gap_mc2 = 0;
    TH1 *leading_eta_star_inside_gap_mc3 = 0;
    TH1 *leading_eta_star_inside_gap_mc4 = 0;
    TH1 *leading_eta_star_inside_gap_mc5 = 0;
    TH1 *leading_eta_star_inside_gap_mc6 = 0;
    TH1 *leading_eta_star_inside_gap_mc7 = 0;
    TH1 *leading_eta_star_inside_gap_mc8 = 0;
    TH1 *leading_eta_star_inside_gap_mc9 = 0;
    TH1 *leading_eta_star_inside_gap_mc10 = 0;
    TH1 *leading_eta_star_inside_gap_mc11 = 0;
    TH1 *leading_eta_star_inside_gap_mc12 = 0;
    TH1 *leading_eta_star_inside_gap_mc13 = 0;
    TH1 *leading_eta_star_inside_gap_mc14 = 0;
    TH1 *leading_eta_star_inside_gap_mc15 = 0;
    TH1 *leading_eta_star_inside_gap_mc16 = 0;
    TH1 *leading_eta_star_inside_gap_mc17 = 0;
    TH1 *leading_eta_star_inside_gap_mc18 = 0;
    TH1 *leading_eta_star_inside_gap_mc19 = 0;
    TH1 *leading_eta_star_inside_gap_mc20 = 0;

    data->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data);
    if (leading_eta_star_inside_gap_data == 0) { cout << "ak5PF_leading_eta_star_inside_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_leading_eta_star_inside_gap",leading_eta_star_inside_gap_data_unc);
    if (leading_eta_star_inside_gap_data_unc == 0) { cout << "ak5PF_leading_eta_star_inside_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_leading_eta_star_inside_gap",leading_eta_star_inside_gap_mc1);
	if (leading_eta_star_inside_gap_mc1 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_leading_eta_star_inside_gap",leading_eta_star_inside_gap_mc2);
	if (leading_eta_star_inside_gap_mc2 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc3);
	if (leading_eta_star_inside_gap_mc3 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc4);
	if (leading_eta_star_inside_gap_mc4 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc5);
	if (leading_eta_star_inside_gap_mc5 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc6);
	if (leading_eta_star_inside_gap_mc6 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc7);
	if (leading_eta_star_inside_gap_mc7 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc8);
	if (leading_eta_star_inside_gap_mc8 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc9);
	if (leading_eta_star_inside_gap_mc9 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc10);
	if (leading_eta_star_inside_gap_mc10 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc11);
	if (leading_eta_star_inside_gap_mc11 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc12);
	if (leading_eta_star_inside_gap_mc12 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc13);
	if (leading_eta_star_inside_gap_mc13 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc14);
	if (leading_eta_star_inside_gap_mc14 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc15);
	if (leading_eta_star_inside_gap_mc15 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc16);
	if (leading_eta_star_inside_gap_mc16 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc17);
	if (leading_eta_star_inside_gap_mc17 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc18);
	if (leading_eta_star_inside_gap_mc18 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc19);
	if (leading_eta_star_inside_gap_mc19 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d06-x01-y01",leading_eta_star_inside_gap_mc20);
	if (leading_eta_star_inside_gap_mc20 == 0) { cout << "ak5Gen_leading_eta_star_inside_gap mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	legend_position = "bottom_middle";
	label_position = "middle";
	}


    create_final_plot(leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_data_unc, leading_eta_star_inside_gap_mc3, labels[2], leading_eta_star_inside_gap_mc4, labels[3], leading_eta_star_inside_gap_mc9, labels[8], leading_eta_star_inside_gap_mc20, labels[19], leading_eta_star_inside_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_eta_star_inside_gap_all", legend_position, label_position);
    create_final_plot(leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_data_unc, leading_eta_star_inside_gap_mc5, labels[4], leading_eta_star_inside_gap_mc8, labels[7], leading_eta_star_inside_gap_mc7, labels[6], leading_eta_star_inside_gap_mc6, labels[5], leading_eta_star_inside_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_eta_star_inside_gap_pythia6", legend_position, label_position);
    create_final_plot(leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_data_unc, leading_eta_star_inside_gap_mc1, labels[0], leading_eta_star_inside_gap_mc2, labels[1], leading_eta_star_inside_gap_mc7, labels[6], leading_eta_star_inside_gap_mc9, labels[8], leading_eta_star_inside_gap_mc20, labels[19], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_eta_star_inside_gap_check", legend_position, label_position);
    create_final_plot(leading_eta_star_inside_gap_data, leading_eta_star_inside_gap_data_unc, leading_eta_star_inside_gap_mc7, labels[6], leading_eta_star_inside_gap_mc6, labels[5], leading_eta_star_inside_gap_mc10, labels[9], leading_eta_star_inside_gap_mc11, labels[10], leading_eta_star_inside_gap_mc12, labels[11], "INSIDE-JET TAG", output_path_plots, mode, "_final_leading_eta_star_inside_gap_hadronization", legend_position, label_position);

//final plot for leading pt outside gap distribution
    if (detail) { cout<<"Leading pT Outside Gap"<<endl; }

    TH1D *leading_pt_outside_gap_data = 0;
    TH1D *leading_pt_outside_gap_data_unc = 0;
    TH1D *leading_pt_outside_gap_mc1 = 0;
    TH1D *leading_pt_outside_gap_mc2 = 0;
    TH1 *leading_pt_outside_gap_mc3 = 0;
    TH1 *leading_pt_outside_gap_mc4 = 0;
    TH1 *leading_pt_outside_gap_mc5 = 0;
    TH1 *leading_pt_outside_gap_mc6 = 0;
    TH1 *leading_pt_outside_gap_mc7 = 0;
    TH1 *leading_pt_outside_gap_mc8 = 0;
    TH1 *leading_pt_outside_gap_mc9 = 0;
    TH1 *leading_pt_outside_gap_mc10 = 0;
    TH1 *leading_pt_outside_gap_mc11 = 0;
    TH1 *leading_pt_outside_gap_mc12 = 0;
    TH1 *leading_pt_outside_gap_mc13 = 0;
    TH1 *leading_pt_outside_gap_mc14 = 0;
    TH1 *leading_pt_outside_gap_mc15 = 0;
    TH1 *leading_pt_outside_gap_mc16 = 0;
    TH1 *leading_pt_outside_gap_mc17 = 0;
    TH1 *leading_pt_outside_gap_mc18 = 0;
    TH1 *leading_pt_outside_gap_mc19 = 0;
    TH1 *leading_pt_outside_gap_mc20 = 0;

    data->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data);
    if (leading_pt_outside_gap_data == 0) { cout << "ak5PF_leading_pt_outside_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_leading_pt_outside_gap",leading_pt_outside_gap_data_unc);
    if (leading_pt_outside_gap_data_unc == 0) { cout << "ak5PF_leading_pt_outside_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_leading_pt_outside_gap",leading_pt_outside_gap_mc1);
	if (leading_pt_outside_gap_mc1 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_leading_pt_outside_gap",leading_pt_outside_gap_mc2);
	if (leading_pt_outside_gap_mc2 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d05-x01-y01",leading_pt_outside_gap_mc3);
	if (leading_pt_outside_gap_mc3 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d05-x01-y01",leading_pt_outside_gap_mc4);
	if (leading_pt_outside_gap_mc4 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d05-x01-y01",leading_pt_outside_gap_mc5);
	if (leading_pt_outside_gap_mc5 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d05-x01-y01",leading_pt_outside_gap_mc6);
	if (leading_pt_outside_gap_mc6 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d05-x01-y01",leading_pt_outside_gap_mc7);
	if (leading_pt_outside_gap_mc7 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d05-x01-y01",leading_pt_outside_gap_mc8);
	if (leading_pt_outside_gap_mc8 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d05-x01-y01",leading_pt_outside_gap_mc9);
	if (leading_pt_outside_gap_mc9 == 0) { cout << "ak5Gen_leading_pt_outnside_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d05-x01-y01",leading_pt_outside_gap_mc10);
	if (leading_pt_outside_gap_mc10 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d05-x01-y01",leading_pt_outside_gap_mc11);
	if (leading_pt_inside_gap_mc11 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d05-x01-y01",leading_pt_outside_gap_mc12);
	if (leading_pt_outside_gap_mc12 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d05-x01-y01",leading_pt_outside_gap_mc13);
	if (leading_pt_outside_gap_mc13 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d05-x01-y01",leading_pt_outside_gap_mc14);
	if (leading_pt_outside_gap_mc14 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d05-x01-y01",leading_pt_outside_gap_mc15);
	if (leading_pt_outside_gap_mc15 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d05-x01-y01",leading_pt_outside_gap_mc16);
	if (leading_pt_outside_gap_mc16 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d05-x01-y01",leading_pt_outside_gap_mc17);
	if (leading_pt_outside_gap_mc17 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d05-x01-y01",leading_pt_outside_gap_mc18);
	if (leading_pt_outside_gap_mc18 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d05-x01-y01",leading_pt_outside_gap_mc19);
	if (leading_pt_outside_gap_mc19 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d05-x01-y01",leading_pt_outside_gap_mc20);
	if (leading_pt_outside_gap_mc20 == 0) { cout << "ak5Gen_leading_pt_outside_gap mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	legend_position = "top_right";
	label_position = "bottom_left";
	}

    create_final_plot(leading_pt_outside_gap_data, leading_pt_outside_gap_data_unc, leading_pt_outside_gap_mc3, labels[2], leading_pt_outside_gap_mc4, labels[3], leading_pt_outside_gap_mc9, labels[8], leading_pt_outside_gap_mc20, labels[19], leading_pt_outside_gap_mc20, labels[19], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_outside_gap_all", legend_position, label_position);
    create_final_plot(leading_pt_outside_gap_data, leading_pt_outside_gap_data_unc, leading_pt_outside_gap_mc5, labels[4], leading_pt_outside_gap_mc8, labels[7], leading_pt_outside_gap_mc7, labels[6], leading_pt_outside_gap_mc6, labels[5], leading_pt_outside_gap_mc20, labels[19], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_outside_gap_pythia6", legend_position, label_position);
    create_final_plot(leading_pt_outside_gap_data, leading_pt_outside_gap_data_unc, leading_pt_outside_gap_mc1, labels[0], leading_pt_outside_gap_mc2, labels[1], leading_pt_outside_gap_mc7, labels[6], leading_pt_outside_gap_mc9, labels[8], leading_pt_outside_gap_mc20, labels[19], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_outside_gap_check", legend_position, label_position);
    create_final_plot(leading_pt_outside_gap_data, leading_pt_outside_gap_data_unc, leading_pt_outside_gap_mc7, labels[6], leading_pt_outside_gap_mc6, labels[5], leading_pt_outside_gap_mc10, labels[9], leading_pt_outside_gap_mc11, labels[10], leading_pt_outside_gap_mc12, labels[11], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_leading_pt_outside_gap_hadronization", legend_position, label_position);


//final plot for delta eta outside gap distribution
    if (detail) { cout<<"Delta Eta Outside Gap"<<endl; }

    TH1D *delta_eta_outside_gap_data = 0;
    TH1D *delta_eta_outside_gap_data_unc = 0;
    TH1D *delta_eta_outside_gap_mc1 = 0;
    TH1D *delta_eta_outside_gap_mc2 = 0;
    TH1 *delta_eta_outside_gap_mc3 = 0;
    TH1 *delta_eta_outside_gap_mc4 = 0;
    TH1 *delta_eta_outside_gap_mc5 = 0;
    TH1 *delta_eta_outside_gap_mc6 = 0;
    TH1 *delta_eta_outside_gap_mc7 = 0;
    TH1 *delta_eta_outside_gap_mc8 = 0;
    TH1 *delta_eta_outside_gap_mc9 = 0;
    TH1 *delta_eta_outside_gap_mc10 = 0;
    TH1 *delta_eta_outside_gap_mc11 = 0;
    TH1 *delta_eta_outside_gap_mc12 = 0;
    TH1 *delta_eta_outside_gap_mc13 = 0;
    TH1 *delta_eta_outside_gap_mc14 = 0;
    TH1 *delta_eta_outside_gap_mc15 = 0;
    TH1 *delta_eta_outside_gap_mc16 = 0;
    TH1 *delta_eta_outside_gap_mc17 = 0;
    TH1 *delta_eta_outside_gap_mc18 = 0;
    TH1 *delta_eta_outside_gap_mc19 = 0;
    TH1 *delta_eta_outside_gap_mc20 = 0;


    data->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data);
    if (delta_eta_outside_gap_data == 0) { cout << "ak5PF_delta_eta_outside_gap not found!" << endl; return; }
    data_unc->GetObject("ak5PF_delta_eta_outside_gap",delta_eta_outside_gap_data_unc);
    if (delta_eta_outside_gap_data_unc == 0) { cout << "ak5PF_delta_eta_outside_gap unc not found!" << endl; return; }
    if (i >= 1)
	{
	mc1->GetObject("ak5Gen_delta_eta_outside_gap",delta_eta_outside_gap_mc1);
	if (delta_eta_outside_gap_mc1 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc1 not found!" << endl; return; }
	}
    if (i >= 2)
	{
	mc2->GetObject("ak5Gen_delta_eta_outside_gap",delta_eta_outside_gap_mc2);
	if (delta_eta_outside_gap_mc2 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc2 not found!" << endl; return; }
	}
    if (i >= 3)
	{
	mc3->GetObject("d07-x01-y01",delta_eta_outside_gap_mc3);
	if (delta_eta_outside_gap_mc3 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc3 not found!" << endl; return; }
	}
    if (i >= 4)
	{
	mc4->GetObject("d07-x01-y01",delta_eta_outside_gap_mc4);
	if (delta_eta_outside_gap_mc4 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc4 not found!" << endl; return; }
	}
    if (i >= 5)
	{
	mc5->GetObject("d07-x01-y01",delta_eta_outside_gap_mc5);
	if (delta_eta_outside_gap_mc5 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc5 not found!" << endl; return; }
	}
    if (i >= 6)
	{
	mc6->GetObject("d07-x01-y01",delta_eta_outside_gap_mc6);
	if (delta_eta_outside_gap_mc6 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc6 not found!" << endl; return; }
	}
    if (i >= 7)
	{
	mc7->GetObject("d07-x01-y01",delta_eta_outside_gap_mc7);
	if (delta_eta_outside_gap_mc7 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc7 not found!" << endl; return; }
	}
    if (i >= 8)
	{
	mc8->GetObject("d07-x01-y01",delta_eta_outside_gap_mc8);
	if (delta_eta_outside_gap_mc8 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc8 not found!" << endl; return; }
	}
    if (i >= 9)
	{
	mc9->GetObject("d07-x01-y01",delta_eta_outside_gap_mc9);
	if (delta_eta_outside_gap_mc9 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc9 not found!" << endl; return; }
	}
    if (i >= 10)
	{
	mc10->GetObject("d07-x01-y01",delta_eta_outside_gap_mc10);
	if (delta_eta_outside_gap_mc10 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc10 not found!" << endl; return; }
	}
    if (i >= 11)
	{
	mc11->GetObject("d07-x01-y01",delta_eta_outside_gap_mc11);
	if (delta_eta_outside_gap_mc11 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc11 not found!" << endl; return; }
	}
    if (i >= 12)
	{
	mc12->GetObject("d07-x01-y01",delta_eta_outside_gap_mc12);
	if (delta_eta_outside_gap_mc12 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc12 not found!" << endl; return; }
	}
    if (i >= 13)
	{
	mc13->GetObject("d07-x01-y01",delta_eta_outside_gap_mc13);
	if (delta_eta_outside_gap_mc13 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc13 not found!" << endl; return; }
	}
    if (i >= 14)
	{
	mc14->GetObject("d07-x01-y01",delta_eta_outside_gap_mc14);
	if (delta_eta_outside_gap_mc14 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc14 not found!" << endl; return; }
	}
    if (i >= 15)
	{
	mc15->GetObject("d07-x01-y01",delta_eta_outside_gap_mc15);
	if (delta_eta_outside_gap_mc15 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc15 not found!" << endl; return; }
	}
    if (i >= 16)
	{
	mc16->GetObject("d07-x01-y01",delta_eta_outside_gap_mc16);
	if (delta_eta_outside_gap_mc16 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc16 not found!" << endl; return; }
	}
    if (i >= 17)
	{
	mc17->GetObject("d07-x01-y01",delta_eta_outside_gap_mc17);
	if (delta_eta_outside_gap_mc17 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc17 not found!" << endl; return; }
	}
    if (i >= 18)
	{
	mc18->GetObject("d07-x01-y01",delta_eta_outside_gap_mc18);
	if (delta_eta_outside_gap_mc18 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc18 not found!" << endl; return; }
	}
    if (i >= 19)
	{
	mc19->GetObject("d07-x01-y01",delta_eta_outside_gap_mc19);
	if (delta_eta_outside_gap_mc19 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc19 not found!" << endl; return; }
	}
    if (i >= 20)
	{
	mc20->GetObject("d07-x01-y01",delta_eta_outside_gap_mc20);
	if (delta_eta_outside_gap_mc20 == 0) { cout << "ak5Gen_delta_eta_outside_gap mc20 not found!" << endl; return; }
	}

    if (mode == "xsec")
	{
	legend_position = "bottom_left";
	label_position = "top_right";
	}

    create_final_plot(delta_eta_outside_gap_data, delta_eta_outside_gap_data_unc, delta_eta_outside_gap_mc3, labels[2], delta_eta_outside_gap_mc4, labels[3], delta_eta_outside_gap_mc9, labels[8], delta_eta_outside_gap_mc20, labels[19], delta_eta_outside_gap_mc20, labels[19], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_delta_eta_outside_gap_all", legend_position, label_position);
    create_final_plot(delta_eta_outside_gap_data, delta_eta_outside_gap_data_unc, delta_eta_outside_gap_mc5, labels[4], delta_eta_outside_gap_mc8, labels[7], delta_eta_outside_gap_mc7, labels[6], delta_eta_outside_gap_mc6, labels[5], delta_eta_outside_gap_mc20, labels[19], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_delta_eta_outside_gap_pythia6", legend_position, label_position);
    create_final_plot(delta_eta_outside_gap_data, delta_eta_outside_gap_data_unc, delta_eta_outside_gap_mc1, labels[0], delta_eta_outside_gap_mc2, labels[1], delta_eta_outside_gap_mc7, labels[6], delta_eta_outside_gap_mc9, labels[8], delta_eta_outside_gap_mc20, labels[19], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_delta_eta_outside_gap_check", legend_position, label_position);
    create_final_plot(delta_eta_outside_gap_data, delta_eta_outside_gap_data_unc, delta_eta_outside_gap_mc7, labels[6], delta_eta_outside_gap_mc6, labels[5], delta_eta_outside_gap_mc10, labels[9], delta_eta_outside_gap_mc11, labels[10], delta_eta_outside_gap_mc12, labels[11], "OUTSIDE-JET TAG", output_path_plots, mode, "_final_delta_eta_outside_gap_hadronization", legend_position, label_position);

//some extra plots
    if (mode == "xsec")
    {
    legend_position = "top_left";
    legend1_position = "top_right";
    label_position = "bottom_right";

    if (detail) { cout<<"Delta phi Scenarios"<<endl; }
    create_scenario_plot(delta_phi_data, "Inclusive scenario", delta_phi_gap_data, "Inside-jet veto scenario", delta_phi_nogap_data, "Inside-jet tag scenario", delta_phi_data_unc, delta_phi_gap_data_unc, delta_phi_nogap_data_unc, delta_phi_mc3, delta_phi_gap_mc3, delta_phi_nogap_mc3, labels[2], delta_phi_mc4, delta_phi_gap_mc4, delta_phi_nogap_mc4, labels[3], delta_phi_mc8, delta_phi_gap_mc8, delta_phi_nogap_mc9, labels[8], delta_phi_mc20,  delta_phi_gap_mc20, delta_phi_nogap_mc20,labels[19], delta_phi_mc20, delta_phi_gap_mc20, delta_phi_nogap_mc20, labels[19], output_path_plots, "_scenario_delta_phi_all", legend_position, legend1_position, label_position);
    create_scenario_plot(delta_phi_data, "Inclusive scenario", delta_phi_gap_data, "Inside-jet veto scenario", delta_phi_nogap_data, "Inside-jet tag scenario", delta_phi_data_unc, delta_phi_gap_data_unc, delta_phi_nogap_data_unc, delta_phi_mc5, delta_phi_gap_mc5, delta_phi_nogap_mc5, labels[4], delta_phi_mc8, delta_phi_gap_mc8, delta_phi_nogap_mc8, labels[7], delta_phi_mc7, delta_phi_gap_mc7, delta_phi_nogap_mc7, labels[6], delta_phi_mc6, delta_phi_gap_mc6, delta_phi_nogap_mc6, labels[5], delta_phi_mc20, delta_phi_gap_mc20, delta_phi_nogap_mc20, labels[19], output_path_plots, "_scenario_delta_phi_pythia6", legend_position, legend1_position, label_position);
    create_scenario_plot(delta_phi_data, "Inclusive scenario", delta_phi_gap_data, "Inside-jet veto scenario", delta_phi_nogap_data, "Inside-jet tag scenario", delta_phi_data_unc, delta_phi_gap_data_unc, delta_phi_nogap_data_unc, delta_phi_mc1, delta_phi_gap_mc1, delta_phi_nogap_mc1,labels[0], delta_phi_mc2, delta_phi_gap_mc2, delta_phi_nogap_mc2, labels[1], delta_phi_mc7, delta_phi_gap_mc7, delta_phi_nogap_mc7, labels[6], delta_phi_mc8, delta_phi_gap_mc8, delta_phi_nogap_mc8, labels[7], delta_phi_mc20, delta_phi_gap_mc20, delta_phi_nogap_mc20, labels[19], output_path_plots, "_scenario_delta_phi_check", legend_position, legend1_position, label_position);
    
    if (detail) { cout<<"Delta phi Deta"<<endl; }
    list_mc[0] = labels[2];
    list_mc[1] = labels[3];
    list_mc[2] = labels[8];
    list_mc[3] = labels[19];
    list_mc[4] = labels[19];        

    create_detamerged_plot(delta_phi_deta1_data, "0.4 < #Delta#eta < 2.5", delta_phi_deta2_data, "2.5 < #Delta#eta < 3.5", delta_phi_deta3_data, "3.5 < #Delta#eta < 4.5", delta_phi_deta4_data, "4.5 < #Delta#eta < 7.5", delta_phi_deta1_data_unc, delta_phi_deta2_data_unc, delta_phi_deta3_data_unc, delta_phi_deta4_data_unc, delta_phi_deta1_mc3, delta_phi_deta2_mc3, delta_phi_deta3_mc3, delta_phi_deta4_mc3, delta_phi_deta1_mc4, delta_phi_deta2_mc4, delta_phi_deta3_mc4, delta_phi_deta4_mc4, delta_phi_deta1_mc8, delta_phi_deta2_mc8, delta_phi_deta3_mc9, delta_phi_deta4_mc9, delta_phi_deta1_mc20,  delta_phi_deta2_mc20, delta_phi_deta3_mc20, delta_phi_deta4_mc20, delta_phi_deta1_mc20, delta_phi_deta2_mc20, delta_phi_deta3_mc20, delta_phi_deta4_mc20, list_mc, "INCLUSIVE", output_path_plots, "_detamerged_delta_phi_all", legend_position, legend1_position, label_position, detail);

    list_mc[0] = labels[4];
    list_mc[1] = labels[7];
    list_mc[2] = labels[6];
    list_mc[3] = labels[5];
    list_mc[4] = labels[19];        

    create_detamerged_plot(delta_phi_deta1_data, "0.4 < #Delta#eta < 2.5", delta_phi_deta2_data, "2.5 < #Delta#eta < 3.5", delta_phi_deta3_data, "3.5 < #Delta#eta < 4.5", delta_phi_deta4_data, "4.5 < #Delta#eta < 7.5", delta_phi_deta1_data_unc, delta_phi_deta2_data_unc, delta_phi_deta3_data_unc, delta_phi_deta4_data_unc, delta_phi_deta1_mc5, delta_phi_deta2_mc5, delta_phi_deta3_mc5, delta_phi_deta4_mc5, delta_phi_deta1_mc8, delta_phi_deta2_mc8, delta_phi_deta3_mc8, delta_phi_deta4_mc8, delta_phi_deta1_mc7, delta_phi_deta2_mc7, delta_phi_deta3_mc7, delta_phi_deta4_mc7, delta_phi_deta1_mc6,  delta_phi_deta2_mc6, delta_phi_deta3_mc6, delta_phi_deta4_mc6, delta_phi_deta1_mc20, delta_phi_deta2_mc20, delta_phi_deta3_mc20, delta_phi_deta4_mc20, list_mc, "INCLUSIVE", output_path_plots, "_detamerged_delta_phi_pythia6", legend_position, legend1_position, label_position, detail);


    if (detail) { cout<<"Delta phi Deta Gap"<<endl; }

    list_mc[0] = labels[2];
    list_mc[1] = labels[3];
    list_mc[2] = labels[8];
    list_mc[3] = labels[19];
    list_mc[4] = labels[19];        

    create_detamerged_plot(delta_phi_deta1_gap_data, "0.4 < #Delta#eta < 2.5", delta_phi_deta2_gap_data, "2.5 < #Delta#eta < 3.5", delta_phi_deta3_gap_data, "3.5 < #Delta#eta < 4.5", delta_phi_deta4_gap_data, "4.5 < #Delta#eta < 7.5", delta_phi_deta1_gap_data_unc, delta_phi_deta2_gap_data_unc, delta_phi_deta3_gap_data_unc, delta_phi_deta4_gap_data_unc, delta_phi_deta1_gap_mc3, delta_phi_deta2_gap_mc3, delta_phi_deta3_gap_mc3, delta_phi_deta4_gap_mc3, delta_phi_deta1_gap_mc4, delta_phi_deta2_gap_mc4, delta_phi_deta3_gap_mc4, delta_phi_deta4_gap_mc4, delta_phi_deta1_gap_mc8, delta_phi_deta2_gap_mc8, delta_phi_deta3_gap_mc9, delta_phi_deta4_gap_mc9, delta_phi_deta1_gap_mc20, delta_phi_deta2_gap_mc20, delta_phi_deta3_gap_mc20, delta_phi_deta4_gap_mc20, delta_phi_deta1_gap_mc20, delta_phi_deta2_gap_mc20, delta_phi_deta3_gap_mc20, delta_phi_deta4_gap_mc20, list_mc, "INSIDE-JET VETO", output_path_plots, "_detamerged_delta_phi_gap_all", legend_position, legend1_position, label_position, detail);

    list_mc[0] = labels[4];
    list_mc[1] = labels[7];
    list_mc[2] = labels[6];
    list_mc[3] = labels[5];
    list_mc[4] = labels[19];        

    create_detamerged_plot(delta_phi_deta1_gap_data, "0.4 < #Delta#eta < 2.5", delta_phi_deta2_gap_data, "2.5 < #Delta#eta < 3.5", delta_phi_deta3_gap_data, "3.5 < #Delta#eta < 4.5", delta_phi_deta4_gap_data, "4.5 < #Delta#eta < 7.5", delta_phi_deta1_gap_data_unc, delta_phi_deta2_gap_data_unc, delta_phi_deta3_gap_data_unc, delta_phi_deta4_gap_data_unc, delta_phi_deta1_gap_mc5, delta_phi_deta2_gap_mc5, delta_phi_deta3_gap_mc5, delta_phi_deta4_gap_mc5, delta_phi_deta1_gap_mc8, delta_phi_deta2_gap_mc8, delta_phi_deta3_gap_mc8, delta_phi_deta4_gap_mc8, delta_phi_deta1_gap_mc7, delta_phi_deta2_gap_mc7, delta_phi_deta3_gap_mc7, delta_phi_deta4_gap_mc7, delta_phi_deta1_gap_mc6, delta_phi_deta2_gap_mc6, delta_phi_deta3_gap_mc6, delta_phi_deta4_gap_mc6, delta_phi_deta1_gap_mc20, delta_phi_deta2_gap_mc20, delta_phi_deta3_gap_mc20, delta_phi_deta4_gap_mc20, list_mc, "INSIDE-JET VETO", output_path_plots, "_detamerged_delta_phi_gap_pythia6", legend_position, legend1_position, label_position, detail);


    if (detail) { cout<<"Delta phi Deta NoGap"<<endl; }
    list_mc[0] = labels[2];
    list_mc[1] = labels[3];
    list_mc[2] = labels[8];
    list_mc[3] = labels[19];
    list_mc[4] = labels[19];        

    create_detamerged_plot(delta_phi_deta1_nogap_data, "0.4 < #Delta#eta < 2.5", delta_phi_deta2_nogap_data, "2.5 < #Delta#eta < 3.5", delta_phi_deta3_nogap_data, "3.5 < #Delta#eta < 4.5", delta_phi_deta4_nogap_data, "4.5 < #Delta#eta < 7.5", delta_phi_deta1_nogap_data_unc, delta_phi_deta2_nogap_data_unc, delta_phi_deta3_nogap_data_unc, delta_phi_deta4_nogap_data_unc, delta_phi_deta1_nogap_mc3, delta_phi_deta2_nogap_mc3, delta_phi_deta3_nogap_mc3, delta_phi_deta4_nogap_mc3, delta_phi_deta1_nogap_mc4, delta_phi_deta2_nogap_mc4, delta_phi_deta3_nogap_mc4, delta_phi_deta4_nogap_mc4, delta_phi_deta1_nogap_mc8, delta_phi_deta2_nogap_mc8, delta_phi_deta3_nogap_mc9, delta_phi_deta4_nogap_mc9, delta_phi_deta1_nogap_mc20, delta_phi_deta2_nogap_mc20, delta_phi_deta3_nogap_mc20, delta_phi_deta4_nogap_mc20, delta_phi_deta1_nogap_mc20, delta_phi_deta2_nogap_mc20, delta_phi_deta3_nogap_mc20, delta_phi_deta4_nogap_mc20, list_mc, "INSIDE-JET TAG", output_path_plots, "_detamerged_delta_phi_nogap_all", legend_position, legend1_position, label_position, detail);

    list_mc[0] = labels[4];
    list_mc[1] = labels[7];
    list_mc[2] = labels[6];
    list_mc[3] = labels[5];
    list_mc[4] = labels[19];        

    create_detamerged_plot(delta_phi_deta1_nogap_data, "0.4 < #Delta#eta < 2.5", delta_phi_deta2_nogap_data, "2.5 < #Delta#eta < 3.5", delta_phi_deta3_nogap_data, "3.5 < #Delta#eta < 4.5", delta_phi_deta4_nogap_data, "4.5 < #Delta#eta < 7.5", delta_phi_deta1_nogap_data_unc, delta_phi_deta2_nogap_data_unc, delta_phi_deta3_nogap_data_unc, delta_phi_deta4_nogap_data_unc, delta_phi_deta1_nogap_mc5, delta_phi_deta2_nogap_mc5, delta_phi_deta3_nogap_mc5, delta_phi_deta4_nogap_mc5, delta_phi_deta1_nogap_mc8, delta_phi_deta2_nogap_mc8, delta_phi_deta3_nogap_mc8, delta_phi_deta4_nogap_mc8, delta_phi_deta1_nogap_mc7, delta_phi_deta2_nogap_mc7, delta_phi_deta3_nogap_mc7, delta_phi_deta4_nogap_mc7, delta_phi_deta1_nogap_mc6, delta_phi_deta2_nogap_mc6, delta_phi_deta3_nogap_mc6, delta_phi_deta4_nogap_mc6, delta_phi_deta1_nogap_mc20, delta_phi_deta2_nogap_mc20, delta_phi_deta3_nogap_mc20, delta_phi_deta4_nogap_mc20, list_mc, "INSIDE-JET TAG", output_path_plots, "_detamerged_delta_phi_nogap_pythia6", legend_position, legend1_position, label_position, detail);


    }
}
