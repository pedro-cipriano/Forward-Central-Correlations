// Pedro Cipriano, May 2013
// DESY, CMS
// Last Update: 06 May 2012
//
//

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveText.h>
#include <TString.h>

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

#include "common_methods.h"

void do_psba(TH2* matched, TH2* nonmatched, TH1* check_correction, int nbins, double *bins, TString label, string path, string prefix, string name, bool detail = true)
{

    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogz();

    matched->Draw("colz");
    matched->Draw("text same");

//setting the output files
   print_plots(c1, path, prefix+name+"_matched");
   c1->Close();

    TCanvas *c2 = new TCanvas("c2","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogz();

    nonmatched->Draw("colz");
    nonmatched->Draw("text same");

//setting the output files
   print_plots(c2, path, prefix+name+"_nonmatched");
   c2->Close();

//declaring variables

   int a;
   double n_matched, n_all;
   double wdt, cent;
   double accept, backg, purit, stab, corr;
   double acc, bac, pur, sta;

//declaring histograms
   TH1D *acceptance;
   TH1D *background;
   TH1D *purity;
   TH1D *stability;
   TH1D *correction;

   TH1D *slice_nonmatched;
   TH1D *slice_matched;
   TH1D *slice;

   acceptance =  new TH1D("acceptance","Acceptance;"+label+";Acceptance", nbins, bins);
   background =  new TH1D("background","Background;"+label+";Background", nbins, bins);
   purity =  new TH1D("purity","Purity;"+label+";Purity", nbins, bins);
   stability =  new TH1D("stability","Stability;"+label+";Stability", nbins, bins);
   correction =  new TH1D("correction","Correction Factor;"+label+";Correction Factor", nbins, bins);

   acceptance->Sumw2();
   background->Sumw2();
   purity->Sumw2();
   stability->Sumw2();
   correction->Sumw2();

   if (detail) { cout<<"Acceptance"<<endl; }
   for (a = 1; a <= nbins;a++)
	{
	slice_nonmatched = nonmatched->ProjectionY("slice_nonmatched",a+1,a+1);
	slice_matched = matched->ProjectionY("slice_matched",a,a);
	//n_matched = slice_nonmatched->Integral(2,nbins+1);
	//n_matched = slice_matched->Integral(1,nbins);
	n_matched = slice_matched->Integral();
	//n_all = slice_nonmatched->Integral(1,nbins+1);
	n_all = slice_nonmatched->Integral();
	wdt = slice_nonmatched->GetBinWidth(a+1);
	cent = slice_nonmatched->GetBinCenter(a+1);
	accept = 0.0;
	if (n_all > 0) {accept = n_matched/n_all; }
	acceptance->SetBinContent(a,accept);
	acceptance->SetBinError(a,0);
	//if (detail)
	//{ cout<<"bin -> "<<a<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Acceptance = "<<accept<<endl; }
	}
   acceptance->SetMinimum(0.0);
   acceptance->SetMaximum(1.2);

   if (detail) { cout<<"Background"<<endl; }
   for (a = 1; a <= nbins;a++)
	{
	slice_nonmatched = nonmatched->ProjectionX("slice_nonmatched",a+1,a+1);
	slice_matched = matched->ProjectionX("slice_matched",a,a);
	//n_matched = slice_nonmatched->Integral(2,nbins+1);
	//n_matched = slice_matched->Integral(1,nbins);
	n_matched = slice_matched->Integral();
	//n_all = slice_nonmatched->Integral(1,nbins+1);
	n_all = slice_nonmatched->Integral();
	wdt = slice_nonmatched->GetBinWidth(a+1);
	cent = slice_nonmatched->GetBinCenter(a+1);
	backg = 0.0;
	if (n_all > 0) { backg = 1.-(n_matched/n_all); }
	background->SetBinContent(a,backg);
	background->SetBinError(a,0);
	//if (detail)
	//{ cout<<"bin -> "<<a<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Background = "<<backg<<endl; }
	}
   background->SetMinimum(0.0);
   background->SetMaximum(1.2);

   if (detail) { cout<<"purity"<<endl; }
   for (a = 1; a <= nbins;a++)
   //for (a = 2; a <= nbins+1;a++)
	{
	slice = matched->ProjectionX("slice",a,a);
	//TH1D *slice = nonmatched->ProjectionX("slice",a,a);
	n_matched = slice->GetBinContent(a);
	//n_all = slice->Integral(2,nbins+1);
	//n_all = slice->Integral(1,nbins);
	n_all = slice->Integral();
	wdt = slice->GetBinWidth(a);
	cent = slice->GetBinCenter(a);
	purit = 0.0;
	if (n_all > 0) { purit = n_matched/n_all; }
	purity->SetBinContent(a,purit);
	purity->SetBinError(a,0);
	//purity->SetBinContent(a-1,purit);
	//purity->SetBinError(a-1,0);
	//if (detail) 
	//{ cout<<"bin -> "<<a<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Purity = "<<purit<<endl; }
	//{ cout<<"bin -> "<<a-1<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Purity = "<<purit<<endl; }
	}
   purity->SetMinimum(0.0);
   purity->SetMaximum(1.2);

   if (detail) { cout<<"Stability"<<endl; }
   for (a = 1; a <= nbins;a++)
   //for (a = 2; a <= nbins+1;a++)
	{
	slice = matched->ProjectionY("slice",a,a);
	//TH1D *slice = nonmatched->ProjectionY("slice",a,a);
	n_matched = slice->GetBinContent(a);
	//n_all = slice->Integral(2,nbins+1);
	n_all = slice->Integral(1,nbins);
	n_all = slice->Integral();
	wdt = slice->GetBinWidth(a);
	cent = slice->GetBinCenter(a);
	stab = 0.0;
	if (n_all > 0) { stab = n_matched/n_all; }
	stability->SetBinContent(a,stab);
	stability->SetBinError(a,0);
	//stability->SetBinContent(a-1,stab);
	//stability->SetBinError(a-1,0);
	//if (detail)
	//{ cout<<"bin -> "<<a<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Stability = "<<stab<<endl; }
	//{ cout<<"bin -> "<<a-1<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Stability = "<<stab<<endl; }
	}
   stability->SetMinimum(0.0);
   stability->SetMaximum(1.2);

	plot_histogram(acceptance, path, prefix+name+"_acceptance", "Acceptance", "top_left", false);
	plot_histogram(background, path, prefix+name+"_background", "Background", "top_left", false);
	plot_histogram(purity, path, prefix+name+"_purity", "Purity", "top_left", false);
	plot_histogram(stability, path, prefix+name+"_stability", "Stability", "top_left", false);

   if (detail) { cout<<"Correction"<<endl; }
   for (a = 1; a <= nbins;a++)
	{
	acc = acceptance->GetBinContent(a);
	bac = background->GetBinContent(a);
	pur = purity->GetBinContent(a);
	sta = stability->GetBinContent(a);
	corr = ((1.-bac)/acc) * (pur/sta);
	correction->SetBinContent(a,corr);
	correction->SetBinError(a,0);
	if (detail)
		{
		//cout<<"bin -> "<<a<<" correction =  "<<corr<< " check = " << check_correction->GetBinContent(a) << " ratio = " << corr/check_correction->GetBinContent(a) << endl;
		//cout << "1-bac / acc = " << (1.-bac)/acc << " put/sta = " << pur/sta << endl;
		}
	}

	plot_2histograms(correction, "PSBA Correction", check_correction, "Had/Det Correction", path, prefix+name+"_correction", "top_left", false, false);

	delete(slice_nonmatched);
	delete(slice_matched);
	delete(slice);
	delete(acceptance);
	delete(background);
	delete(purity);
	delete(stability);
	delete(correction);
}

void do_psba_unmatched(TH2* nonmatched, TH1* check_correction, int nbins, double *bins, TString label, string path, string prefix, string name, bool detail = true)
{

    TCanvas *c2 = new TCanvas("c2","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogz();

    nonmatched->Draw("colz");
    nonmatched->Draw("text same");

//setting the output files
   print_plots(c2, path, prefix+name+"_nonmatched");
   c2->Close();

//declaring variables

   int a;
   double n_matched, n_all;
   double wdt, cent;
   double accept, backg, purit, stab, corr;
   double acc, bac, pur, sta;

//declaring histograms
   TH1D *acceptance;
   TH1D *background;
   TH1D *purity;
   TH1D *stability;
   TH1D *correction;

   TH1D *slice1;
   TH1D *slice2;
   TH1D *slice3;
   TH1D *slice4;

   acceptance =  new TH1D("acceptance","Acceptance;"+label+";Acceptance", nbins, bins);
   background =  new TH1D("background","Background;"+label+";Background", nbins, bins);
   purity =  new TH1D("purity","Purity;"+label+";Purity", nbins, bins);
   stability =  new TH1D("stability","Stability;"+label+";Stability", nbins, bins);
   correction =  new TH1D("correction","Correction Factor;"+label+";Correction Factor", nbins, bins);

   acceptance->Sumw2();
   background->Sumw2();
   purity->Sumw2();
   stability->Sumw2();
   correction->Sumw2();

   if (detail) { cout<<"Acceptance"<<endl; }
   for (a = 2; a <= nbins+1;a++)
	{
	slice1 = nonmatched->ProjectionY("slice1",a,a);
	n_matched = slice1->Integral(2,nbins+1);
	n_all = slice1->Integral(1,nbins+1);
	wdt = slice1->GetBinWidth(a);
	cent = slice1->GetBinCenter(a);
	accept = 0.0;
	if (n_all > 0) {accept = n_matched/n_all; }
	acceptance->SetBinContent(a-1,accept);
	acceptance->SetBinError(a-1,0);
	if (detail)
	{ cout<<"bin -> "<<a-1<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Acceptance = "<<accept<<endl; }
	}
   acceptance->SetMinimum(0.0);
   acceptance->SetMaximum(1.2);

   if (detail) { cout<<"Background"<<endl; }
   for (a = 2; a <= nbins+1;a++)
	{
	slice2 = nonmatched->ProjectionX("slice2",a,a);
	n_matched = slice2->Integral(2,nbins+1);
	n_all = slice2->Integral(1,nbins+1);
	wdt = slice2->GetBinWidth(a);
	cent = slice2->GetBinCenter(a);
	backg = 0.0;
	if (n_all > 0) { backg = 1.-(n_matched/n_all); }
	background->SetBinContent(a-1,backg);
	background->SetBinError(a-1,0);
	if (detail)
	{ cout<<"bin -> "<<a-1<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Background = "<<backg<<endl; }
	}
   background->SetMinimum(0.0);
   background->SetMaximum(1.2);

   if (detail) { cout<<"purity"<<endl; }
   for (a = 2; a <= nbins+1;a++)
	{
	slice3 = nonmatched->ProjectionX("slice3",a,a);
	n_matched = slice3->GetBinContent(a);
	n_all = slice3->Integral(2,nbins+1);
	wdt = slice3->GetBinWidth(a);
	cent = slice3->GetBinCenter(a);
	purit = 0.0;
	if (n_all > 0) { purit = n_matched/n_all; }
	purity->SetBinContent(a-1,purit);
	purity->SetBinError(a-1,0);
	if (detail) 
	{ cout<<"bin -> "<<a-1<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Purity = "<<purit<<endl; }
	}
   purity->SetMinimum(0.0);
   purity->SetMaximum(1.2);

   if (detail) { cout<<"Stability"<<endl; }
   for (a = 2; a <= nbins+1;a++)
	{
	slice4 = nonmatched->ProjectionY("slice4",a,a);
	n_matched = slice4->GetBinContent(a);
	n_all = slice4->Integral(2,nbins+1);
	wdt = slice4->GetBinWidth(a);
	cent = slice4->GetBinCenter(a);
	stab = 0.0;
	if (n_all > 0) { stab = n_matched/n_all; }
	stability->SetBinContent(a-1,stab);
	stability->SetBinError(a-1,0);
	if (detail)
	{ cout<<"bin -> "<<a-1<<" n_matched =  "<<n_matched<<" n_all = "<<n_all<<" Stability = "<<stab<<endl; }
	}
   stability->SetMinimum(0.0);
   stability->SetMaximum(1.2);

	plot_histogram(acceptance, path, prefix+name+"_acceptance", "Acceptance", "top_left", false);
	plot_histogram(background, path, prefix+name+"_background", "Background", "top_left", false);
	plot_histogram(purity, path, prefix+name+"_purity", "Purity", "top_left", false);
	plot_histogram(stability, path, prefix+name+"_stability", "Stability", "top_left", false);

   if (detail) { cout<<"Correction"<<endl; }
   for (a = 1; a <= nbins;a++)
	{
	acc = acceptance->GetBinContent(a);
	bac = background->GetBinContent(a);
	pur = purity->GetBinContent(a);
	sta = stability->GetBinContent(a);
	corr = ((1.-bac)/acc) * (pur/sta);
	correction->SetBinContent(a,corr);
	correction->SetBinError(a,0);
	if (detail)
		{
		cout<<"bin -> "<<a<<" correction =  "<<corr<< " check = " << check_correction->GetBinContent(a) << " ratio = " << corr/check_correction->GetBinContent(a) << endl;
		//cout << "1-bac / acc = " << (1.-bac)/acc << " put/sta = " << pur/sta << endl;
		}
	}

	plot_2histograms(correction, "PSBA Correction", check_correction, "Had/Det Correction", path, prefix+name+"_correction", "top_left", false, false);

	delete(slice1);
	delete(slice2);
	delete(slice3);
	delete(slice4);
	delete(acceptance);
	delete(background);
	delete(purity);
	delete(stability);
	delete(correction);
}

void plot_matched_matrix(TH2* matched, string path, string prefix, string name, bool detail = true)
{

    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2g");
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.10);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    gPad->SetLogz();

    matched->Draw("colz");
    matched->Draw("text same");

//setting the output files
   print_plots(c1, path, prefix+name+"_matched");
   c1->Close();
}


void compute_psba(string root_in, string correction, string corr_prefix, string path_plots, string prefix, bool detail = true, bool test = false)
{

//output the configuration
   if (detail) { cout<<"Compute PSBA Configuration"<<endl; }
   if (detail) { cout<<"Input :             "<<root_in<<endl; }
   if (detail) { cout<<"Corrections :       "<<correction<<endl; }
   if (detail) { cout<<"Correction Prefix : "<<corr_prefix<<endl; }
   if (detail) { cout<<"Path Plots :        "<<path_plots<<endl; }
   if (detail) { cout<<"Prefix :            "<<prefix<<endl; }
   if (detail) { cout<<"Detail level :      "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :         "<<test<<endl; }

    if (detail) { cout << "Opening files..." << endl; }
//opens the MC files
    TFile *mc= new TFile( root_in.c_str() );
    TFile *corr= new TFile( correction.c_str() );

//binning
   int deta_nbins = 4;
   double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};

   int dphi_nbins = 7;
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

   int in_nbins = 9;
   double in_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

   int etastar_nbins = 12;
   double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

   int out_nbins = 9;
   double out_bins[10] = {20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

   int deta_out_nbins = 6;
   double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

//compute delta phi psba
    if (detail) { cout << "Computing Delta Phi PSBA..." << endl; }

    TH2D *matched_delta_phi_matched = 0;
    mc->GetObject("matched_delta_phi",matched_delta_phi_matched);
    if (matched_delta_phi_matched == 0) { cout << "matched_delta_phi not found!" << endl; return; }
    TH2D *matched_delta_phi_nonmatched = 0;
    mc->GetObject("matched_delta_phi_nomatch",matched_delta_phi_nonmatched);
    if (matched_delta_phi_nonmatched == 0) { cout << "matched_delta_phi_nomatch not found!" << endl; return; }
    TH1D *corr_delta_phi = 0;
    TString corr_delta_phi_name = corr_prefix + "delta_phi";
    corr->GetObject(corr_delta_phi_name,corr_delta_phi);
    if (corr_delta_phi == 0) { cout << corr_delta_phi_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_matched, matched_delta_phi_nonmatched, corr_delta_phi, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi", detail);
    do_psba_unmatched(matched_delta_phi_nonmatched, corr_delta_phi, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi", detail);
    plot_matched_matrix(matched_delta_phi_matched, path_plots, prefix, "delta_phi", detail);


//compute delta phi deta1 psba
    if (detail) { cout << "Computing Delta Phi Deta1 PSBA..." << endl; }

    TH2D *matched_delta_phi_deta1_matched = 0;
    mc->GetObject("matched_delta_phi_deta1",matched_delta_phi_deta1_matched);
    if (matched_delta_phi_deta1_matched == 0) { cout << "matched_delta_phi_deta1_match source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta1_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta1_nomatch",matched_delta_phi_deta1_nonmatched);
    if (matched_delta_phi_deta1_nonmatched == 0) { cout << "matched_delta_phi_deta1_nomatch source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta1 = 0;
    TString corr_delta_phi_deta1_name = corr_prefix + "delta_phi_deta1";
    corr->GetObject(corr_delta_phi_deta1_name,corr_delta_phi_deta1);
    if (corr_delta_phi_deta1 == 0) { cout << corr_delta_phi_deta1_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta1_matched, matched_delta_phi_deta1_nonmatched, corr_delta_phi_deta1, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta1", detail);
    do_psba_unmatched(matched_delta_phi_deta1_nonmatched, corr_delta_phi_deta1, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta1", detail);
    plot_matched_matrix(matched_delta_phi_deta1_matched, path_plots, prefix, "delta_phi_deta1", detail);

//compute delta phi deta2 psba
    if (detail) { cout << "Computing Delta Phi Deta2 PSBA..." << endl; }

    TH2D *matched_delta_phi_deta2_matched = 0;
    mc->GetObject("matched_delta_phi_deta2",matched_delta_phi_deta2_matched);
    if (matched_delta_phi_deta2_matched == 0) { cout << "matched_delta_phi_deta2 source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta2_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta2_nomatch",matched_delta_phi_deta2_nonmatched);
    if (matched_delta_phi_deta2_nonmatched == 0) { cout << "matched_delta_phi_deta2_nomatch source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta2 = 0;
    TString corr_delta_phi_deta2_name = corr_prefix + "delta_phi_deta2";
    corr->GetObject(corr_delta_phi_deta2_name,corr_delta_phi_deta2);
    if (corr_delta_phi_deta2 == 0) { cout << corr_delta_phi_deta2_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta2_matched, matched_delta_phi_deta2_nonmatched, corr_delta_phi_deta2, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta2", detail);
    do_psba_unmatched(matched_delta_phi_deta2_nonmatched, corr_delta_phi_deta2, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta2", detail);
    plot_matched_matrix(matched_delta_phi_deta2_matched, path_plots, prefix, "delta_phi_deta2", detail);

//compute delta phi deta3 psba
    if (detail) { cout << "Computing Delta Phi Deta3 PSBA..." << endl; }

    TH2D *matched_delta_phi_deta3_matched = 0;
    mc->GetObject("matched_delta_phi_deta3",matched_delta_phi_deta3_matched);
    if (matched_delta_phi_deta3_matched == 0) { cout << "matched_delta_phi_deta3 source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta3_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta3_nomatch",matched_delta_phi_deta3_nonmatched);
    if (matched_delta_phi_deta3_nonmatched == 0) { cout << "matched_delta_phi_deta3_nomatch source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta3 = 0;
    TString corr_delta_phi_deta3_name = corr_prefix + "delta_phi_deta3";
    corr->GetObject(corr_delta_phi_deta3_name,corr_delta_phi_deta3);
    if (corr_delta_phi_deta3 == 0) { cout << corr_delta_phi_deta3_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta3_matched, matched_delta_phi_deta3_nonmatched, corr_delta_phi_deta3, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta3", detail);
    do_psba_unmatched(matched_delta_phi_deta3_nonmatched, corr_delta_phi_deta3, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta3", detail);
    plot_matched_matrix(matched_delta_phi_deta3_matched, path_plots, prefix, "delta_phi_deta3", detail);


//compute delta phi deta4 psba
    if (detail) { cout << "Computing Delta Phi Deta4 PSBA..." << endl; }

    TH2D *matched_delta_phi_deta4_matched = 0;
    mc->GetObject("matched_delta_phi_deta4",matched_delta_phi_deta4_matched);
    if (matched_delta_phi_deta4_matched == 0) { cout << "matched_delta_phi_deta4 source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta4_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta4_nomatch",matched_delta_phi_deta4_nonmatched);
    if (matched_delta_phi_deta4_nonmatched == 0) { cout << "matched_delta_phi_deta4_nomatch source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta4 = 0;
    TString corr_delta_phi_deta4_name = corr_prefix + "delta_phi_deta4";
    corr->GetObject(corr_delta_phi_deta4_name,corr_delta_phi_deta4);
    if (corr_delta_phi_deta4 == 0) { cout << corr_delta_phi_deta4_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta4_matched, matched_delta_phi_deta4_nonmatched, corr_delta_phi_deta4, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta4", detail);
    do_psba_unmatched(matched_delta_phi_deta4_nonmatched, corr_delta_phi_deta4, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta4", detail);
    plot_matched_matrix(matched_delta_phi_deta4_matched, path_plots, prefix, "delta_phi_deta4", detail);


//compute delta phi gap psba
    if (detail) { cout << "Computing Delta Phi Gap PSBA..." << endl; }

    TH2D *matched_delta_phi_gap_matched = 0;
    mc->GetObject("matched_delta_phi_gap",matched_delta_phi_gap_matched);
    if (matched_delta_phi_gap_matched == 0) { cout << "matched_delta_phi_gap source not found!" << endl; return; }
    TH2D *matched_delta_phi_gap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_gap_nomatch",matched_delta_phi_gap_nonmatched);
    if (matched_delta_phi_gap_nonmatched == 0) { cout << "matched_delta_phi_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_gap = 0;
    TString corr_delta_phi_gap_name = corr_prefix + "delta_phi_gap";
    corr->GetObject(corr_delta_phi_gap_name,corr_delta_phi_gap);
    if (corr_delta_phi_gap == 0) { cout << corr_delta_phi_gap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_gap_matched, matched_delta_phi_gap_nonmatched, corr_delta_phi_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_gap", detail);
    do_psba_unmatched(matched_delta_phi_gap_nonmatched, corr_delta_phi_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_gap", detail);
    plot_matched_matrix(matched_delta_phi_gap_matched, path_plots, prefix, "delta_phi_gap", detail);


//compute delta phi deta1 gap psba
    if (detail) { cout << "Computing Delta Phi Deta1 Gap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta1_gap_matched = 0;
    mc->GetObject("matched_delta_phi_deta1_gap",matched_delta_phi_deta1_gap_matched);
    if (matched_delta_phi_deta1_gap_matched == 0) { cout << "matched_delta_phi_deta1_gap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta1_gap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta1_gap_nomatch",matched_delta_phi_deta1_gap_nonmatched);
    if (matched_delta_phi_deta1_gap_nonmatched == 0) { cout << "matched_delta_phi_deta1_gap_nonmatche source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta1_gap = 0;
    TString corr_delta_phi_deta1_gap_name = corr_prefix + "delta_phi_deta1_gap";
    corr->GetObject(corr_delta_phi_deta1_gap_name,corr_delta_phi_deta1_gap);
    if (corr_delta_phi_deta1_gap == 0) { cout << corr_delta_phi_deta1_gap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta1_gap_matched, matched_delta_phi_deta1_gap_nonmatched, corr_delta_phi_deta1_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta1_gap", detail);
    do_psba_unmatched(matched_delta_phi_deta1_gap_nonmatched, corr_delta_phi_deta1_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta1_gap", detail);
    plot_matched_matrix(matched_delta_phi_deta1_gap_matched, path_plots, prefix, "delta_phi_deta1_gap", detail);


//compute delta phi deta2 gap psba
    if (detail) { cout << "Computing Delta Phi Deta2 Gap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta2_gap_matched = 0;
    mc->GetObject("matched_delta_phi_deta2_gap",matched_delta_phi_deta2_gap_matched);
    if (matched_delta_phi_deta2_gap_matched == 0) { cout << "matched_delta_phi_deta2_gap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta2_gap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta2_gap_nomatch",matched_delta_phi_deta2_gap_nonmatched);
    if (matched_delta_phi_deta2_gap_nonmatched == 0) { cout << "matched_delta_phi_deta2_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta2_gap = 0;
    TString corr_delta_phi_deta2_gap_name = corr_prefix + "delta_phi_deta2_gap";
    corr->GetObject(corr_delta_phi_deta2_gap_name,corr_delta_phi_deta2_gap);
    if (corr_delta_phi_deta2_gap == 0) { cout << corr_delta_phi_deta2_gap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta2_gap_matched, matched_delta_phi_deta2_gap_nonmatched, corr_delta_phi_deta2_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta2_gap", detail);
    do_psba_unmatched(matched_delta_phi_deta2_gap_nonmatched, corr_delta_phi_deta2_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta2_gap", detail);
    plot_matched_matrix(matched_delta_phi_deta2_gap_matched, path_plots, prefix, "delta_phi_deta2_gap", detail);


//compute delta phi deta3 gap psba
    if (detail) { cout << "Computing Delta Phi Deta3 Gap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta3_gap_matched = 0;
    mc->GetObject("matched_delta_phi_deta3_gap",matched_delta_phi_deta3_gap_matched);
    if (matched_delta_phi_deta3_gap_matched == 0) { cout << "matched_delta_phi_deta3_gap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta3_gap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta3_gap_nomatch",matched_delta_phi_deta3_gap_nonmatched);
    if (matched_delta_phi_deta3_gap_nonmatched == 0) { cout << "matched_delta_phi_deta3_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta3_gap = 0;
    TString corr_delta_phi_deta3_gap_name = corr_prefix + "delta_phi_deta3_gap";
    corr->GetObject(corr_delta_phi_deta3_gap_name,corr_delta_phi_deta3_gap);
    if (corr_delta_phi_deta3_gap == 0) { cout << corr_delta_phi_deta3_gap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta3_gap_matched, matched_delta_phi_deta3_gap_nonmatched, corr_delta_phi_deta3_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta3_gap", detail);
    do_psba_unmatched(matched_delta_phi_deta3_gap_nonmatched, corr_delta_phi_deta3_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta3_gap", detail);
    plot_matched_matrix(matched_delta_phi_deta3_gap_matched, path_plots, prefix, "delta_phi_deta3_gap", detail);


//compute delta phi deta4 gap psba
    if (detail) { cout << "Computing Delta Phi Deta4 Gap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta4_gap_matched = 0;
    mc->GetObject("matched_delta_phi_deta4_gap",matched_delta_phi_deta4_gap_matched);
    if (matched_delta_phi_deta4_gap_matched == 0) { cout << "matched_delta_phi_deta4_gap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta4_gap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta4_gap_nomatch",matched_delta_phi_deta4_gap_nonmatched);
    if (matched_delta_phi_deta4_gap_nonmatched == 0) { cout << "matched_delta_phi_deta4_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta4_gap = 0;
    TString corr_delta_phi_deta4_gap_name = corr_prefix + "delta_phi_deta4_gap";
    corr->GetObject(corr_delta_phi_deta4_gap_name,corr_delta_phi_deta4_gap);
    if (corr_delta_phi_deta4_gap == 0) { cout << corr_delta_phi_deta4_gap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta4_gap_matched, matched_delta_phi_deta4_gap_nonmatched, corr_delta_phi_deta4_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta4_gap", detail);
    do_psba_unmatched(matched_delta_phi_deta4_gap_nonmatched, corr_delta_phi_deta4_gap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta4_gap", detail);
    plot_matched_matrix(matched_delta_phi_deta4_gap_matched, path_plots, prefix, "delta_phi_deta4_gap", detail);


//compute delta phi nogap psba
    if (detail) { cout << "Computing Delta Phi NoGap PSBA..." << endl; }

    TH2D *matched_delta_phi_nogap_matched = 0;
    mc->GetObject("matched_delta_phi_nogap",matched_delta_phi_nogap_matched);
    if (matched_delta_phi_nogap_matched == 0) { cout << "matched_delta_phi_nogap source not found!" << endl; return; }
    TH2D *matched_delta_phi_nogap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_nogap_nomatch",matched_delta_phi_nogap_nonmatched);
    if (matched_delta_phi_nogap_nonmatched == 0) { cout << "matched_delta_phi_nogap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_nogap = 0;
    TString corr_delta_phi_nogap_name = corr_prefix + "delta_phi_nogap";
    corr->GetObject(corr_delta_phi_nogap_name,corr_delta_phi_nogap);
    if (corr_delta_phi_nogap == 0) { cout << corr_delta_phi_nogap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_nogap_matched, matched_delta_phi_nogap_nonmatched, corr_delta_phi_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_nogap", detail);
    do_psba_unmatched(matched_delta_phi_nogap_nonmatched, corr_delta_phi_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_nogap", detail);
    plot_matched_matrix(matched_delta_phi_nogap_matched, path_plots, prefix, "delta_phi_nogap", detail);


//compute delta phi deta1 nogap psba
    if (detail) { cout << "Computing Delta Phi Deta1 NoGap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta1_nogap_matched = 0;
    mc->GetObject("matched_delta_phi_deta1_nogap",matched_delta_phi_deta1_nogap_matched);
    if (matched_delta_phi_deta1_nogap_matched == 0) { cout << "matched_delta_phi_deta1_nogap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta1_nogap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta1_nogap_nomatch",matched_delta_phi_deta1_nogap_nonmatched);
    if (matched_delta_phi_deta1_nogap_nonmatched == 0) { cout << "matched_delta_phi_deta1_nogap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta1_nogap = 0;
    TString corr_delta_phi_deta1_nogap_name = corr_prefix + "delta_phi_deta1_nogap";
    corr->GetObject(corr_delta_phi_deta1_nogap_name,corr_delta_phi_deta1_nogap);
    if (corr_delta_phi_deta1_nogap == 0) { cout << corr_delta_phi_deta1_nogap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta1_nogap_matched, matched_delta_phi_deta1_nogap_nonmatched, corr_delta_phi_deta1_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta1_nogap", detail);
    do_psba_unmatched(matched_delta_phi_deta1_nogap_nonmatched, corr_delta_phi_deta1_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta1_nogap", detail);
    plot_matched_matrix(matched_delta_phi_deta1_nogap_matched, path_plots, prefix, "delta_phi_deta1_nogap", detail);


//compute delta phi deta2 nogap psba
    if (detail) { cout << "Computing Delta Phi Deta2 NoGap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta2_nogap_matched = 0;
    mc->GetObject("matched_delta_phi_deta2_nogap",matched_delta_phi_deta2_nogap_matched);
    if (matched_delta_phi_deta2_nogap_matched == 0) { cout << "matched_delta_phi_deta2_nogap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta2_nogap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta2_nogap_nomatch",matched_delta_phi_deta2_nogap_nonmatched);
    if (matched_delta_phi_deta2_nogap_nonmatched == 0) { cout << "matched_delta_phi_deta2_nogap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta2_nogap = 0;
    TString corr_delta_phi_deta2_nogap_name = corr_prefix + "delta_phi_deta2_nogap";
    corr->GetObject(corr_delta_phi_deta2_nogap_name,corr_delta_phi_deta2_nogap);
    if (corr_delta_phi_deta2_nogap == 0) { cout << corr_delta_phi_deta2_nogap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta2_nogap_matched, matched_delta_phi_deta2_nogap_nonmatched, corr_delta_phi_deta2_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta2_nogap", detail);
    do_psba_unmatched(matched_delta_phi_deta2_nogap_nonmatched, corr_delta_phi_deta2_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta2_nogap", detail);
    plot_matched_matrix(matched_delta_phi_deta2_nogap_matched, path_plots, prefix, "delta_phi_deta2_nogap", detail);


//compute delta phi deta3 nogap psba
    if (detail) { cout << "Computing Delta Phi Deta3 NoGap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta3_nogap_matched = 0;
    mc->GetObject("matched_delta_phi_deta3_nogap",matched_delta_phi_deta3_nogap_matched);
    if (matched_delta_phi_deta3_nogap_matched == 0) { cout << "matched_delta_phi_deta3_nogap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta3_nogap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta3_nogap_nomatch",matched_delta_phi_deta3_nogap_nonmatched);
    if (matched_delta_phi_deta3_nogap_nonmatched == 0) { cout << "matched_delta_phi_deta3_nogap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta3_nogap = 0;
    TString corr_delta_phi_deta3_nogap_name = corr_prefix + "delta_phi_deta3_nogap";
    corr->GetObject(corr_delta_phi_deta3_nogap_name,corr_delta_phi_deta3_nogap);
    if (corr_delta_phi_deta3_nogap == 0) { cout << corr_delta_phi_deta3_nogap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta3_nogap_matched, matched_delta_phi_deta3_nogap_nonmatched, corr_delta_phi_deta3_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta3_nogap", detail);
    do_psba_unmatched(matched_delta_phi_deta3_nogap_nonmatched, corr_delta_phi_deta3_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta3_nogap", detail);
    plot_matched_matrix(matched_delta_phi_deta3_nogap_matched, path_plots, prefix, "delta_phi_deta3_nogap", detail);


//compute delta phi deta4 nogap psba
    if (detail) { cout << "Computing Delta Phi Deta4 NoGap PSBA..." << endl; }

    TH2D *matched_delta_phi_deta4_nogap_matched = 0;
    mc->GetObject("matched_delta_phi_deta4_nogap",matched_delta_phi_deta4_nogap_matched);
    if (matched_delta_phi_deta4_nogap_matched == 0) { cout << "matched_delta_phi_deta4_nogap source not found!" << endl; return; }
    TH2D *matched_delta_phi_deta4_nogap_nonmatched = 0;
    mc->GetObject("matched_delta_phi_deta4_nogap_nomatch",matched_delta_phi_deta4_nogap_nonmatched);
    if (matched_delta_phi_deta4_nogap_nonmatched == 0) { cout << "matched_delta_phi_deta4_nogap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_phi_deta4_nogap = 0;
    TString corr_delta_phi_deta4_nogap_name = corr_prefix + "delta_phi_deta4_nogap";
    corr->GetObject(corr_delta_phi_deta4_nogap_name,corr_delta_phi_deta4_nogap);
    if (corr_delta_phi_deta4_nogap == 0) { cout << corr_delta_phi_deta4_nogap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_phi_deta4_nogap_matched, matched_delta_phi_deta4_nogap_nonmatched, corr_delta_phi_deta4_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta4_nogap", detail);
    do_psba_unmatched(matched_delta_phi_deta4_nogap_nonmatched, corr_delta_phi_deta4_nogap, dphi_nbins, dphi_bins, "#Delta#phi [rad]", path_plots, prefix, "delta_phi_deta4_nogap", detail);
    plot_matched_matrix(matched_delta_phi_deta4_nogap_matched, path_plots, prefix, "delta_phi_deta4_nogap", detail);


//compute leading pt inside gap psba
    if (detail) { cout << "Computing Leading pT Inside Gap PSBA..." << endl; }

    TH2D *matched_leading_pt_inside_gap_matched = 0;
    mc->GetObject("matched_leading_pt_inside_gap",matched_leading_pt_inside_gap_matched);
    if (matched_leading_pt_inside_gap_matched == 0) { cout << "matched_leading_pt_inside_gap source not found!" << endl; return; }
    TH2D *matched_leading_pt_inside_gap_nonmatched = 0;
    mc->GetObject("matched_leading_pt_inside_gap_nomatch",matched_leading_pt_inside_gap_nonmatched);
    if (matched_leading_pt_inside_gap_nonmatched == 0) { cout << "matched_leading_pt_inside_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_leading_pt_inside_gap = 0;
    TString corr_leading_pt_inside_gap_name = corr_prefix + "leading_pt_inside_gap";
    corr->GetObject(corr_leading_pt_inside_gap_name,corr_leading_pt_inside_gap);
    if (corr_leading_pt_inside_gap == 0) { cout << corr_leading_pt_inside_gap_name + " not found!" << endl; return; }

    //do_psba(matched_leading_pt_inside_gap_matched, matched_leading_pt_inside_gap_nonmatched, corr_leading_pt_inside_gap, in_nbins, in_bins, "p_{T} [GeV]", path_plots, prefix, "leading_pt_inside_gap", detail);
    do_psba_unmatched(matched_leading_pt_inside_gap_nonmatched, corr_leading_pt_inside_gap, in_nbins, in_bins, "p_{T} [GeV]", path_plots, prefix, "leading_pt_inside_gap", detail);
    plot_matched_matrix(matched_leading_pt_inside_gap_matched, path_plots, prefix, "leading_pt_inside_gap", detail);

//compute leading eta star inside gap psba
    if (detail) { cout << "Computing Leading Eta* Inside Gap PSBA..." << endl; }

    TH2D *matched_leading_eta_star_inside_gap_matched = 0;
    mc->GetObject("matched_leading_eta_star_inside_gap",matched_leading_eta_star_inside_gap_matched);
    if (matched_leading_eta_star_inside_gap_matched == 0) { cout << "matched_leading_eta_star_inside_gap source not found!" << endl; return; }
    TH2D *matched_leading_eta_star_inside_gap_nonmatched = 0;
    mc->GetObject("matched_leading_eta_star_inside_gap_nomatch",matched_leading_eta_star_inside_gap_nonmatched);
    if (matched_leading_eta_star_inside_gap_nonmatched == 0) { cout << "matched_leading_eta_star_inside_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_leading_eta_star_inside_gap = 0;
    TString corr_leading_eta_star_inside_gap_name = corr_prefix + "leading_eta_star_inside_gap";
    corr->GetObject(corr_leading_eta_star_inside_gap_name,corr_leading_eta_star_inside_gap);
    if (corr_leading_eta_star_inside_gap == 0) { cout << corr_leading_eta_star_inside_gap_name + " not found!" << endl; return; }

    //do_psba(matched_leading_eta_star_inside_gap_matched, matched_leading_eta_star_inside_gap_nonmatched, corr_leading_eta_star_inside_gap, etastar_nbins, etastar_bins, "#Eta*", path_plots, prefix, "leading_eta_star_inside_gap", detail);
    do_psba_unmatched(matched_leading_eta_star_inside_gap_nonmatched, corr_leading_eta_star_inside_gap, etastar_nbins, etastar_bins, "#Eta*", path_plots, prefix, "leading_eta_star_inside_gap", detail);
    plot_matched_matrix(matched_leading_eta_star_inside_gap_matched, path_plots, prefix, "leading_eta_star_inside_gap", detail);

//compute leading pt outside gap psba
    if (detail) { cout << "Computing Leading pT Outside Gap PSBA..." << endl; }

    TH2D *matched_leading_pt_outside_gap_matched = 0;
    mc->GetObject("matched_leading_pt_outside_gap",matched_leading_pt_outside_gap_matched);
    if (matched_leading_pt_outside_gap_matched == 0) { cout << "matched_leading_pt_outside_gap source not found!" << endl; return; }
    TH2D *matched_leading_pt_outside_gap_nonmatched = 0;
    mc->GetObject("matched_leading_pt_outside_gap_nomatch",matched_leading_pt_outside_gap_nonmatched);
    if (matched_leading_pt_outside_gap_nonmatched == 0) { cout << "matched_leading_pt_outside_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_leading_pt_outside_gap = 0;
    TString corr_leading_pt_outside_gap_name = corr_prefix + "leading_pt_outside_gap";
    corr->GetObject(corr_leading_pt_outside_gap_name,corr_leading_pt_outside_gap);
    if (corr_leading_pt_outside_gap == 0) { cout << corr_leading_pt_outside_gap_name + " not found!" << endl; return; }

    //do_psba(matched_leading_pt_outside_gap_matched, matched_leading_pt_outside_gap_nonmatched, corr_leading_pt_outside_gap, out_nbins, out_bins, "p_{T} [GeV]", path_plots, prefix, "leading_pt_outside_gap", detail);
    do_psba_unmatched(matched_leading_pt_outside_gap_nonmatched, corr_leading_pt_outside_gap, out_nbins, out_bins, "p_{T} [GeV]", path_plots, prefix, "leading_pt_outside_gap", detail);
    plot_matched_matrix(matched_leading_pt_outside_gap_matched, path_plots, prefix, "leading_pt_outside_gap", detail);


//compute delta eta outside gap psba
    if (detail) { cout << "Computing Delta Eta Outside Gap PSBA..." << endl; }

    TH2D *matched_delta_eta_outside_gap_matched = 0;
    mc->GetObject("matched_delta_eta_outside_gap",matched_delta_eta_outside_gap_matched);
    if (matched_delta_eta_outside_gap_matched == 0) { cout << "matched_delta_eta_outside_gap source not found!" << endl; return; }
    TH2D *matched_delta_eta_outside_gap_nonmatched = 0;
    mc->GetObject("matched_delta_eta_outside_gap_nomatch",matched_delta_eta_outside_gap_nonmatched);
    if (matched_delta_eta_outside_gap_nonmatched == 0) { cout << "matched_delta_eta_outside_gap_nonmatched source not found!" << endl; return; }
    TH1D *corr_delta_eta_outside_gap = 0;
    TString corr_delta_eta_outside_gap_name = corr_prefix + "delta_eta_outside_gap";
    corr->GetObject(corr_delta_eta_outside_gap_name,corr_delta_eta_outside_gap);
    if (corr_delta_eta_outside_gap == 0) { cout << corr_delta_eta_outside_gap_name + " not found!" << endl; return; }

    //do_psba(matched_delta_eta_outside_gap_matched, matched_delta_eta_outside_gap_nonmatched, corr_delta_eta_outside_gap, deta_out_nbins, deta_out_bins, "#delta#Eta^{outside}", path_plots, prefix, "delta_eta_outside_gap", detail);
    do_psba_unmatched(matched_delta_eta_outside_gap_nonmatched, corr_delta_eta_outside_gap, deta_out_nbins, deta_out_bins, "#delta#Eta^{outside}", path_plots, prefix, "delta_eta_outside_gap", detail);
    plot_matched_matrix(matched_delta_eta_outside_gap_matched, path_plots, prefix, "delta_eta_outside_gap", detail);
}
