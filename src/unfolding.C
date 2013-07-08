// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 17 Jun 2013
//
// unfolding(string response_file, string path_data, string mc_file, string root_out, bool detail = false)
// unfolds the delta_phi distribution present in the path_data with TUnfold method using the response matrix suplied by response_file and saves the result in root_out

#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TPad.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TSpline.h>
#include "TLegend.h"
#include "TDecompSVD.h"

#include "../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldTUnfold.h" 

#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include "common_methods.h"

// Global definitions
const int n_loop = 10000;
const int seed = 5;
//TUnfold::ERegMode regmode=TUnfold::kRegModeNone;
TUnfold::ERegMode regmode=TUnfold::kRegModeSize;
//TUnfold::ERegMode regmode=TUnfold::kRegModeDerivative;
//TUnfold::ERegMode regmode=TUnfold::kRegModeCurvature;

void plot_covariance(TMatrixD cov, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{
//declaring the canvas
    if (detail) { cout << "Ploting " << name << endl; }
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

    cov.Draw("colz");
    cov.Draw("text same");


//setting the output files
   string fileout = prefix + name;
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c1->Print( out_png.c_str() );
    c1->Print( out_c.c_str() );
    c1->Print( out_eps.c_str() );
    c1->Close();

}


void plot_lcurve(TGraph *lcurve, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{
//declaring the canvas
    if (detail) { cout << "Ploting " << name << endl; }
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
    //gPad->SetLogy();

    lcurve->SetTitle("L-curve;log(#chi^2);log(reg_cond)");
    lcurve->Draw("AC*");

//setting the output files
   string fileout = prefix + name;
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c1->Print( out_png.c_str() );
    c1->Print( out_c.c_str() );
    c1->Print( out_eps.c_str() );
    c1->Close();

}


void plot_taucurves(TSpline *logTauX, TSpline *logTauY, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{
    if (detail) { cout << "Ploting " << name << endl; }
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
    //gPad->SetLogy();
    c1->Divide(1,2);

    c1->cd(1);
    logTauX->SetTitle("#tau as function of #chi^2;log(#tau);log(#chi^2)");
    logTauX->Draw("");
    c1->cd(2);
    logTauY->SetTitle("#tau as function of regularization condition;log(#tau);log(#reg_cond)");
    logTauY->Draw("");


//setting the output files
   string fileout = prefix + name;
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c1->Print( out_png.c_str() );
    c1->Print( out_c.c_str() );
    c1->Print( out_eps.c_str() );
    c1->Close();
}

void plot_unfolded_result(TH1D *unfolded, TH1D *true_dist, TH1D *measured_dist, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
{

//declaring the canvas
    if (detail) { cout << "Ploting " << name << endl; }
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

for (int j = 1; j <= true_dist->GetNbinsX(); j++)
	{
	//cout << j << " : " << true_dist->GetBinContent(j) << " +- " << true_dist->GetBinError(j) << endl;
	}

for (int j = 1; j <= measured_dist->GetNbinsX(); j++)
	{
	//cout << j << " : " << measured_dist->GetBinContent(j) << " +- " << measured_dist->GetBinError(j) << endl;
	}

//calculate the plooting range
    if (detail) { cout << "Getting the minimum and maximum for the plot..." << endl; }
    double min = 0.0;
    double max = unfolded->GetMaximum();
    if (unfolded->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(unfolded,detail);
    }
    else
    {
    min = unfolded->GetMinimum();
    } 

    set_histogram_min_max(true_dist, min, max, detail);
    set_histogram_min_max(measured_dist, min, max, detail);
//    set_histogram_min_max(measured_plot, min, max, detail);

    max = 1.3 * max;
    min = 0.7 * min;

//format and ploting the histogram
    if (detail) { cout << "Drawning on the canvas..." << endl; }
    unfolded->SetMaximum(max);
    unfolded->SetMinimum(min);
    unfolded->SetLineColor(2);
    unfolded->SetLineStyle(1);
    unfolded->SetLineWidth(3);
    unfolded->Draw("e1");
    true_dist->SetLineColor(3);
    true_dist->SetLineStyle(2);
    true_dist->SetLineWidth(3);
    true_dist->Draw("e1 same");
    measured_dist->SetLineColor(4);
    measured_dist->SetLineStyle(3);
    measured_dist->SetLineWidth(3);
    measured_dist->Draw("e1 same");
//    measured_plot->SetLineColor(6);
//    measured_plot->SetLineStyle(4);
//    measured_plot->SetLineWidth(3);
//    measured_plot->Draw("e1 same");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1,y1,x2,y2);
    leg00->AddEntry(unfolded,"Unfolded result","l");
    leg00->AddEntry(true_dist,"True distribution","l");
    leg00->AddEntry(measured_dist,"Measured distribution","l");
//    leg00->AddEntry(measured_plot,"Measured distribution","l");
    leg00->SetFillColor(0);
    leg00->SetLineStyle(1);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();

//setting the output files
   print_plots(c1, path, prefix + name);

}

void unfolding(string response_file, string data_file, string mc_file, string root_out, string method, double scale, double& tau, string prefix, bool true_is_corr_data = false, bool detail = false, bool test = false)
{

if (detail)
{
cout << "Unfolding Configuration" << endl;
cout << "Response File :  " << response_file << endl;
cout << "Data File :      " << data_file << endl;
cout << "MC File :        " << mc_file << endl;
cout << "Root output :    " << root_out << endl;
cout << "Method :         " << method << endl;
cout << "Scaling Factor : " << scale << endl;
cout << "Prefix :         " << prefix << endl;
cout << "Detail :         " << detail << endl;
cout << "Test Mode :      " << test << endl;
}

//opening the files
TFile *respfile = new TFile( response_file.c_str() );
TFile *datafile = new TFile( data_file.c_str() );
TFile *mcfile = new TFile( mc_file.c_str() );

//declare the objects names
TString hname = "resp_delta_phi";
TString hname_all = "resp_delta_phi_all";
TString h_measured = "ak5PF_delta_phi_fine";
TString h_measured_plot = "ak5PF_delta_phi";
TString h_true;
if (true_is_corr_data)
	{
	h_true = "ak5PF_delta_phi";
	}
else
	{
	h_true = "ak5Gen_delta_phi";
	}

//get the objects on a safe way
RooUnfoldResponse *response = 0;
respfile->GetObject(hname,response);
if (response == 0) { cout<<"Response not found!"<<endl; return; }

RooUnfoldResponse *response_all = 0;
respfile->GetObject(hname_all,response_all);
if (response_all == 0) { cout<<"Response_all not found!"<<endl; return; }

TH1D *detected = 0;
datafile->GetObject(h_measured,detected);
//respfile->GetObject("ak5PF_delta_phi_measured",detected);
if (detected == 0) { cout<<"Measured data " << h_measured << " not found!"<<endl; return; }

TH1D *detected_plot = 0;
datafile->GetObject(h_measured_plot,detected_plot);
if (detected_plot == 0) { cout<<"Measured data not found!"<<endl; return; }

TH1D *generated = 0;
mcfile->GetObject(h_true,generated);
//respfile->GetObject("ak5Gen_delta_phi_truth",generated);
if (generated == 0) { cout<< h_true << " true data not found!"<<endl; return; }

//setting the binning for the unfolding inputs
int dphi_nbins_meas = 15;
double dphi_bins_meas[16] = {-1.0, 0.0, 0.225, 0.45, 0.675, 0.9, 1.125, 1.35, 1.575, 1.8, 2.025, 2.25, 2.475, 2.7, 2.925, 3.15};

  TH1D *input_detected;
  input_detected =  new TH1D("input_detected_delta_phi","#Delta#phi;|#Delta#phi|;Events", dphi_nbins_meas, dphi_bins_meas);

	input_detected->SetBinContent(1,0.0);
	input_detected->SetBinError(1,0.0);

for (int a = 1; a <= detected->GetNbinsX();a++)
	{
	input_detected->SetBinContent(a+1,detected->GetBinContent(a));
	input_detected->SetBinError(a+1,detected->GetBinError(a));
	cout << "Bin : " << a << " -> " << detected->GetBinContent(a) << " +- " << detected->GetBinError(a) << endl;
	}


//initialize the output objects
TH1D *hReco = 0;
TMatrixD cov(7,7);
TMatrixD invertida(7,7);
TMatrixD chk(7,7);
TGraph *lcurve = 0;
TSpline *logTauX = 0, *logTauY = 0;

//  Creating the RooUnfoldTUnfold object
if (method == "TUnfold")
	{
	RooUnfoldTUnfold unfold(response, detected, regmode);
	hReco = (TH1D*) unfold.Hreco();
	tau = unfold.GetTau();
	cout<<"Tau : "<<tau<<endl;
	lcurve = static_cast<TGraph*>(unfold.GetLCurve()->Clone());
	logTauX = static_cast<TSpline*>(unfold.GetLogTauX()->Clone());
	logTauY = static_cast<TSpline*>(unfold.GetLogTauY()->Clone());
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (method == "BinByBin")
	{
	RooUnfoldBinByBin unfold(response_all, input_detected);
	hReco = (TH1D*) unfold.Hreco();
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (method == "SVD")
	{
	RooUnfoldSvd unfold(response_all, input_detected, 5);
	hReco = (TH1D*) unfold.Hreco();
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (method == "Bayes")
	{
	RooUnfoldBayes unfold(response, detected, 1);
	hReco = (TH1D*) unfold.Hreco();
	cov = static_cast<TMatrixD>(unfold.Ereco(RooUnfold::kCovariance));
	}

if (hReco == 0) { cout << "Invalid unfolding method! Routine termination!!!" << endl; return; }

// unfolding


//because we input the double of the bins the result is the double of normal, we need to scale the histogram
hReco->Scale(scale);

//matrix investion and cross-check multiplication
invertida = cov;
invertida.Invert();
chk = cov * invertida;

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

  TH1D *output_true;
  output_true =  new TH1D("output_true_delta_phi","#Delta#phi;|#Delta#phi| [rad];Events", dphi_nbins, dphi_bins);
  int fac = 0;

for (int a = 1; a <= output_true->GetNbinsX();a++)
	{
	if (method == "TUnfold" or method == "Bayes") { fac = a;} else { fac = a+1; }
	output_true->SetBinContent(a,hReco->GetBinContent(fac));
	output_true->SetBinError(a,hReco->GetBinError(fac));
	cout << "Bin : " << a << " -> " << output_true->GetBinContent(a) << " +- " << output_true->GetBinError(a) << endl;
	}


//ploting the result
if (method == "TUnfold" or method == "Bayes")
	{
	plot_unfolded_result(hReco, generated, detected, "../output/unfolding/", prefix, "delta_phi", "top_left", detail);

	}
else
	{
	plot_unfolded_result(output_true, generated, detected, "../output/unfolding/", prefix, "delta_phi", "top_left", detail);
	}

if (method == "TUnfold")
	{
	plot_lcurve(lcurve, "../output/unfolding/", prefix, "lcurve", "top_left", detail);
	plot_taucurves(logTauX, logTauY, "../output/unfolding/", prefix, "taucurves", "top_left", detail);
	}

plot_covariance(cov, "../output/unfolding/", prefix, "covariance", "top_left", detail);
plot_covariance(invertida, "../output/unfolding/", prefix, "inv_covariance", "top_left", detail);
plot_covariance(chk, "../output/unfolding/", prefix, "inv_covariance_mult_convariance", "top_left", detail);


// recreating the output file
TFile *f = TFile::Open(root_out.c_str(), "RECREATE");

// save the unfolded result in the output file
output_true->Write();
cout<<"Unfolding result saved to : "<<root_out<<endl;

//delete the variables to avoid memory leak
delete f;
delete hReco;
//delete lcurve;
//delete logTauX;
//delete logTauY;

}
