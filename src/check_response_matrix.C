// Pedro Cipriano, Dez 2012
// DESY, CMS
// Last Update: 11 Fev 2012
//
// check_response_matrix(string response_file, string check_response_dir = "../output/check_response/", string prefix = "test_", bool detail = false, bool test = false)
// checks the behaviour of the response martrix supplied

#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TPad.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TSpline.h>
#include <TLegend.h>

#include "../RooUnfold-1.1.1/src/RooUnfold.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldResponse.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

#include "common_methods.h"

void plot_response(TH2 *response, string path, string prefix, string name, string legend_position = "top_left", bool normalized = false, bool detail = false, bool test = false)
{

//to avoid memory leaks
TH1D* proj = 0;

//normalizing the matrix
if (normalized)
	{
	if (detail) { cout << "Normalizing response matrix" << endl; }
	for (int i = 1; i <= response->GetNbinsY(); i++)
		{
		proj = response->ProjectionX("proj",i,i);
		double integral = proj->Integral();
		if (test) { cout << "i = " << i << " -- Integral = " << integral << endl; }
		for (int j = 1; j <= response->GetNbinsX(); j++)
			{
			double val = response->GetBinContent(j,i) / integral;
			if (test) { cout << j << "-" << i <<  " : " << val << endl; }
			response->SetBinContent(j,i,val);
			}
		}
	}


for (int i = 1; i <= response->GetNbinsY(); i++)
	{
	for (int j = 1; j <= response->GetNbinsX(); j++)
		{
		double val = response->GetBinContent(j,i);
		//cout << j << "-" << i <<  " : " << val << endl;
		}
	}


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
    if (!normalized) { gPad->SetLogz(); }

    response->Draw("colz");
    response->Draw("text same");

    if (detail) { cout << "No legend will be used on this plot. If needed the " << legend_position << " position will be used!" << endl; }

//setting the output files
   string fileout = prefix + name;
   print_plots(c1, path, fileout);

//to avoid memory leaks
   delete(proj);

}


void plot_dist(TH1 *dist, string path, string prefix, string name, string legend_position = "top_left", TString dist_label = "distribution", bool detail = false)
{


for (int j = 1; j <= dist->GetNbinsX(); j++)
	{
	//cout << j << " : " << dist->GetBinContent(j) << " +- " << dist->GetBinError(j) << endl;
	}

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
    gPad->SetLogy();

    dist->Draw("e1");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1,y1,x2,y2);
    leg00->AddEntry(dist,dist_label,"l");
    leg00->SetFillColor(0);
    leg00->SetLineStyle(1);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();


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

void check_response_matrix(string response_file, string check_response_dir = "../output/check_response/", string prefix = "test_", bool detail = false, bool test = false)
{

if (detail)
{
cout << "Check Response Matrix Configuration" << endl;
cout << "Response File : " << response_file << endl;
cout << "Plots output :  " << check_response_dir << endl;
cout << "Prefix :        " << prefix << endl;
cout << "Detail :        " << detail << endl;
cout << "Test Mode :     " << test << endl;
}

//opening the files
if (detail) { cout << "Opening Response Root File..." << endl; }
TFile *respfile = new TFile( response_file.c_str() );

//declare the objects names
if (detail) { cout << "Leading histograms..." << endl; }
TString hname = "resp_delta_phi";
TString hname_all = "resp_delta_phi_all";
TString delta_phi_miss_name = "ak5Gen_delta_phi_miss";
TString miss_det_name = "ak5PF_delta_phi_miss";
TString notmissed_det_name = "ak5PF_delta_phi_notmissed";
TString truth_name = "ak5Gen_delta_phi_truth";
TString measured_name = "ak5PF_delta_phi_measured";
TString truth_missed_name = "ak5Gen_delta_phi_truth_missed";
TString truth_stayed_name = "ak5Gen_delta_phi_truth_stayed";
TString truth_fakes_name = "ak5Gen_delta_phi_truth_fakes";
TString truth_fakes_slightly_25_name = "ak5Gen_delta_phi_truth_fakes_slightly_25";
TString truth_fakes_slightly_26_name = "ak5Gen_delta_phi_truth_fakes_slightly_26";
TString truth_fakes_slightly_27_name = "ak5Gen_delta_phi_truth_fakes_slightly_27";
TString truth_fakes_slightly_28_name = "ak5Gen_delta_phi_truth_fakes_slightly_28";
TString truth_fakes_slightly_29_name = "ak5Gen_delta_phi_truth_fakes_slightly_29";
TString truth_fakes_slightly_30_name = "ak5Gen_delta_phi_truth_fakes_slightly_30";
TString truth_fakes_slightly_31_name = "ak5Gen_delta_phi_truth_fakes_slightly_31";
TString truth_fakes_slightly_32_name = "ak5Gen_delta_phi_truth_fakes_slightly_32";
TString truth_fakes_slightly_33_name = "ak5Gen_delta_phi_truth_fakes_slightly_33";
TString truth_fakes_slightly_34_name = "ak5Gen_delta_phi_truth_fakes_slightly_34";

//get the objects on a safe way
RooUnfoldResponse *delta_phi_response = 0;
respfile->GetObject(hname,delta_phi_response);
if (delta_phi_response == 0) { cout<<"Response delta_phi not found!"<<endl; return; }

RooUnfoldResponse *delta_phi_response_all = 0;
respfile->GetObject(hname_all,delta_phi_response_all);
if (delta_phi_response_all == 0) { cout<<"Response_all delta_phi not found!"<<endl; return; }

TH1D *delta_phi_miss = 0;
respfile->GetObject(delta_phi_miss_name,delta_phi_miss);
if (delta_phi_miss == 0) { cout<<"Miss delta_phi histogram not found!"<<endl; return; }

TH1D *miss_det = 0;
respfile->GetObject(miss_det_name,miss_det);
if (miss_det == 0) { cout<<"Miss on detector level histogram not found!"<<endl; return; }

TH1D *notmissed_det = 0;
respfile->GetObject(notmissed_det_name,notmissed_det);
if (notmissed_det == 0) { cout<<"NotMissed on detector level histogram not found!"<<endl; return; }

TH1D *truth = 0;
respfile->GetObject(truth_name,truth);
if (truth == 0) { cout<<"Truth histogram not found!"<<endl; return; }

TH1D *measured = 0;
respfile->GetObject(measured_name,measured);
if (measured == 0) { cout<<"Measured histogram not found!"<<endl; return; }

TH1D *truth_missed = 0;
respfile->GetObject(truth_missed_name,truth_missed);
if (truth_missed == 0) { cout<<"Truth Missed histogram not found!"<<endl; return; }

TH1D *truth_stayed = 0;
respfile->GetObject(truth_stayed_name,truth_stayed);
if (truth_stayed == 0) { cout<<"Truth Stayed histogram not found!"<<endl; return; }

TH1D *truth_fakes = 0;
respfile->GetObject(truth_fakes_name,truth_fakes);
if (truth_fakes == 0) { cout<<"Truth Fakes histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_25 = 0;
respfile->GetObject(truth_fakes_slightly_25_name,truth_fakes_slightly_25);
if (truth_fakes_slightly_25 == 0) { cout<<"Truth Fakes Slightly 25 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_26 = 0;
respfile->GetObject(truth_fakes_slightly_26_name,truth_fakes_slightly_26);
if (truth_fakes_slightly_26 == 0) { cout<<"Truth Fakes Slightly 26 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_27 = 0;
respfile->GetObject(truth_fakes_slightly_27_name,truth_fakes_slightly_27);
if (truth_fakes_slightly_27 == 0) { cout<<"Truth Fakes Slightly 27 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_28 = 0;
respfile->GetObject(truth_fakes_slightly_28_name,truth_fakes_slightly_28);
if (truth_fakes_slightly_28 == 0) { cout<<"Truth Fakes Slightly 28 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_29 = 0;
respfile->GetObject(truth_fakes_slightly_29_name,truth_fakes_slightly_29);
if (truth_fakes_slightly_29 == 0) { cout<<"Truth Fakes Slightly 29 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_30 = 0;
respfile->GetObject(truth_fakes_slightly_30_name,truth_fakes_slightly_30);
if (truth_fakes_slightly_30 == 0) { cout<<"Truth Fakes Slightly 30 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_31 = 0;
respfile->GetObject(truth_fakes_slightly_31_name,truth_fakes_slightly_31);
if (truth_fakes_slightly_31 == 0) { cout<<"Truth Fakes Slightly 31 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_32 = 0;
respfile->GetObject(truth_fakes_slightly_32_name,truth_fakes_slightly_32);
if (truth_fakes_slightly_32 == 0) { cout<<"Truth Fakes Slightly 32 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_33 = 0;
respfile->GetObject(truth_fakes_slightly_33_name,truth_fakes_slightly_33);
if (truth_fakes_slightly_33 == 0) { cout<<"Truth Fakes Slightly 33 GeV histogram not found!"<<endl; return; }

TH1D *truth_fakes_slightly_34 = 0;
respfile->GetObject(truth_fakes_slightly_34_name,truth_fakes_slightly_34);
if (truth_fakes_slightly_34 == 0) { cout<<"Truth Fakes Slightly 34 GeV histogram not found!"<<endl; return; }


TH2 *delta_phi_hResponse = (TH2*) delta_phi_response->Hresponse();
TH2 *delta_phi_hResponse_all = (TH2*) delta_phi_response_all->Hresponse();
TH2 *delta_phi_hResponse_trim = (TH2*) delta_phi_response->HresponseNoOverflow();
TH1D *delta_phi_hFakes = (TH1D*) delta_phi_response->Hfakes();
TH1D *delta_phi_hTruth = (TH1D*) delta_phi_response->Htruth();

//plotting histograms
if (detail) { cout << "Ploting histograms..." << endl; }
plot_response(delta_phi_hResponse, check_response_dir, prefix, "delta_phi_response", "top_left", false, detail, test);
plot_response(delta_phi_hResponse_all, check_response_dir, prefix, "delta_phi_response_all", "top_left", false, detail, test);
plot_response(delta_phi_hResponse_trim, check_response_dir, prefix, "detal_phi_response_nooverflow", "top_left", false, detail, test);
plot_response(delta_phi_hResponse, check_response_dir, prefix, "delta_phi_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_hResponse_all, check_response_dir, prefix, "delta_phi_response_all_norm", "top_left", true, detail, test);
plot_response(delta_phi_hResponse_trim, check_response_dir, prefix, "delta_phi_response_nooverflow_norm", "top_left", true, detail, test);
plot_dist(delta_phi_hFakes, check_response_dir, prefix, "delta_phi_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_miss, check_response_dir, prefix, "delta_phi_miss", "top_left", "Miss True", detail);
plot_dist(miss_det, check_response_dir, prefix, "delta_phi_miss_det", "top_left", "Miss Detector", detail);
plot_dist(notmissed_det, check_response_dir, prefix, "delta_phi_notmissed_det", "top_left", "Miss Detector", detail);
plot_dist(truth, check_response_dir, prefix, "delta_phi_truth", "top_left", "Truth", detail);
plot_dist(measured, check_response_dir, prefix, "delta_phi_measured", "top_left", "Measured", detail);
plot_dist(truth_missed, check_response_dir, prefix, "delta_phi_truth_missed", "top_left", "Truth of the missed events", detail);
plot_dist(truth_stayed, check_response_dir, prefix, "delta_phi_truth_stayed", "top_left", "Truth of the events in response", detail);
plot_dist(truth_fakes, check_response_dir, prefix, "delta_phi_truth_fakes", "top_left", "Truth of the fakes", detail);
plot_dist(delta_phi_hTruth, check_response_dir, prefix, "delta_phi_truth_from_response", "top_left", "Truth from response matrix", detail);

plot_dist(truth_fakes_slightly_25, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_25", "top_left", "Truth of the fakes missed slightly 25 GeV", detail);
plot_dist(truth_fakes_slightly_26, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_26", "top_left", "Truth of the fakes missed slightly 26 GeV", detail);
plot_dist(truth_fakes_slightly_27, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_27", "top_left", "Truth of the fakes missed slightly 27 GeV", detail);
plot_dist(truth_fakes_slightly_28, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_28", "top_left", "Truth of the fakes missed slightly 28 GeV", detail);
plot_dist(truth_fakes_slightly_29, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_29", "top_left", "Truth of the fakes missed slightly 29 GeV", detail);
plot_dist(truth_fakes_slightly_30, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_30", "top_left", "Truth of the fakes missed slightly 30 GeV", detail);
plot_dist(truth_fakes_slightly_31, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_31", "top_left", "Truth of the fakes missed slightly 31 GeV", detail);
plot_dist(truth_fakes_slightly_32, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_32", "top_left", "Truth of the fakes missed slightly 32 GeV", detail);
plot_dist(truth_fakes_slightly_33, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_33", "top_left", "Truth of the fakes missed slightly 33 GeV", detail);
plot_dist(truth_fakes_slightly_34, check_response_dir, prefix, "delta_phi_truth_fakes_slightly_34", "top_left", "Truth of the fakes missed slightly 34 GeV", detail);

plot_six_dist(truth_fakes, "truth_fakes", truth_fakes_slightly_25, "truth_fakes_slightly 25 GeV", truth_fakes_slightly_26, "truth_fakes_slightly 26 GeV", truth_fakes_slightly_27, "truth_fakes_slightly 27 GeV", truth_fakes_slightly_28, "truth_fakes_slightly 28 GeV", truth_fakes_slightly_29, "truth_fakes_slightly 29 GeV", check_response_dir, prefix, "delta_phi_fakes_slightly_all_far", "top_left", detail);

plot_six_dist(truth_fakes, "truth_fakes", truth_fakes_slightly_30, "truth_fakes_slightly 30 GeV", truth_fakes_slightly_31, "truth_fakes_slightly 31 GeV", truth_fakes_slightly_32, "truth_fakes_slightly 32 GeV", truth_fakes_slightly_33, "truth_fakes_slightly 33 GeV", truth_fakes_slightly_34, "truth_fakes_slightly 34 GeV", check_response_dir, prefix, "delta_phi_fakes_slightly_all_close", "top_left", detail);

ratio_2histograms(truth_fakes_slightly_25, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_25_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_26, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_26_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_27, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_27_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_28, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_28_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_29, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_29_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_30, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_30_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_31, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_31_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_32, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_32_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_33, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_33_over_fakes", "top_left", detail);
ratio_2histograms(truth_fakes_slightly_34, truth_fakes, "Fakes Slightly / Fakes", check_response_dir, prefix +  "delta_phi_ratio_fakes_slightly_34_over_fakes", "top_left", detail);

ratio_2histograms(truth_missed, truth_stayed, "Truth Missed/ Truth Stayed", check_response_dir, prefix + "delta_phi_ratio_truth_missed_over_truth_stayed", "bottom_left", detail);
ratio_2histograms(miss_det, notmissed_det, "Measured Fakes/ Measured Without Fakes", check_response_dir, prefix + "delta_phi_ratio_measured_fakes_over_measured_without_fakes", "bottom_left", detail);
ratio_2histograms(miss_det, truth_missed, "Miss Detector/ Truth Missed", check_response_dir, prefix + "delta_phi_ratio_miss_over_truth_missed", "top_left", detail);
plot_six_dist(miss_det, "miss_detector", truth, "truth", truth_missed, "truth_missed", truth_stayed, "truth_stayed", truth_fakes, "truth_fakes", truth_fakes_slightly_30, "truth_fakes_slightly 30 GeV", check_response_dir, prefix, "delta_phi_all", "top_left", detail);


//get the delta_phi_deta1 response on a safe way
TString delta_phi_deta1_hname = "resp_delta_phi_deta1";
RooUnfoldResponse *delta_phi_deta1_response = 0;
respfile->GetObject(delta_phi_deta1_hname,delta_phi_deta1_response);
if (delta_phi_deta1_response == 0) { cout<<"Response delta_phi_deta1 not found!"<<endl; return; }

TString delta_phi_deta1_hname_all = "resp_delta_phi_deta1_all";
RooUnfoldResponse *delta_phi_deta1_response_all = 0;
respfile->GetObject(delta_phi_deta1_hname_all,delta_phi_deta1_response_all);
if (delta_phi_deta1_response_all == 0) { cout<<"Response_all delta_phi_deta1 not found!"<<endl; return; }

TString delta_phi_deta1_miss_name = "ak5Gen_delta_phi_deta1_miss";
TH1D *delta_phi_deta1_miss = 0;
respfile->GetObject(delta_phi_deta1_miss_name,delta_phi_deta1_miss);
if (delta_phi_deta1_miss == 0) { cout<<"Miss delta_phi_deta1 histogram not found!"<<endl; return; }

TH2 *delta_phi_deta1_hResponse = (TH2*) delta_phi_deta1_response->Hresponse();
TH2 *delta_phi_deta1_hResponse_all = (TH2*) delta_phi_deta1_response_all->Hresponse();
TH1 *delta_phi_deta1_hFakes = (TH1D*) delta_phi_deta1_response->Hfakes();
TH1 *delta_phi_deta1_hTruth = (TH1D*) delta_phi_deta1_response->Htruth();

plot_response(delta_phi_deta1_hResponse, check_response_dir, prefix, "delta_phi_deta1_response", "top_left", false, detail, test);
plot_response(delta_phi_deta1_hResponse_all, check_response_dir, prefix, "delta_phi_deta1_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta1_hResponse, check_response_dir, prefix, "delta_phi_deta1_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta1_hResponse_all, check_response_dir, prefix, "delta_phi_deta1_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta1_hFakes, check_response_dir, prefix, "delta_phi_deta1_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta1_hTruth, check_response_dir, prefix, "delta_phi_deta1_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta1_miss, check_response_dir, prefix, "delta_phi_deta1_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta2 response on a safe way
TString delta_phi_deta2_hname = "resp_delta_phi_deta2";
RooUnfoldResponse *delta_phi_deta2_response = 0;
respfile->GetObject(delta_phi_deta2_hname,delta_phi_deta2_response);
if (delta_phi_deta2_response == 0) { cout<<"Response delta_phi_deta2 not found!"<<endl; return; }

TString delta_phi_deta2_hname_all = "resp_delta_phi_deta2_all";
RooUnfoldResponse *delta_phi_deta2_response_all = 0;
respfile->GetObject(delta_phi_deta2_hname_all,delta_phi_deta2_response_all);
if (delta_phi_deta2_response_all == 0) { cout<<"Response_all delta_phi_deta2 not found!"<<endl; return; }

TString delta_phi_deta2_miss_name = "ak5Gen_delta_phi_deta2_miss";
TH1D *delta_phi_deta2_miss = 0;
respfile->GetObject(delta_phi_deta2_miss_name,delta_phi_deta2_miss);
if (delta_phi_deta2_miss == 0) { cout<<"Miss delta_phi_deta2 histogram not found!"<<endl; return; }

TH2 *delta_phi_deta2_hResponse = (TH2*) delta_phi_deta2_response->Hresponse();
TH2 *delta_phi_deta2_hResponse_all = (TH2*) delta_phi_deta2_response_all->Hresponse();
TH1 *delta_phi_deta2_hFakes = (TH1D*) delta_phi_deta2_response->Hfakes();
TH1 *delta_phi_deta2_hTruth = (TH1D*) delta_phi_deta2_response->Htruth();

plot_response(delta_phi_deta2_hResponse, check_response_dir, prefix, "delta_phi_deta2_response", "top_left", false, detail, test);
plot_response(delta_phi_deta2_hResponse_all, check_response_dir, prefix, "delta_phi_deta2_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta2_hResponse, check_response_dir, prefix, "delta_phi_deta2_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta2_hResponse_all, check_response_dir, prefix, "delta_phi_deta2_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta2_hFakes, check_response_dir, prefix, "delta_phi_deta2_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta2_hTruth, check_response_dir, prefix, "delta_phi_deta2_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta2_miss, check_response_dir, prefix, "delta_phi_deta2_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta3 response on a safe way
TString delta_phi_deta3_hname = "resp_delta_phi_deta3";
RooUnfoldResponse *delta_phi_deta3_response = 0;
respfile->GetObject(delta_phi_deta3_hname,delta_phi_deta3_response);
if (delta_phi_deta3_response == 0) { cout<<"Response delta_phi_deta3 not found!"<<endl; return; }

TString delta_phi_deta3_hname_all = "resp_delta_phi_deta3_all";
RooUnfoldResponse *delta_phi_deta3_response_all = 0;
respfile->GetObject(delta_phi_deta3_hname_all,delta_phi_deta3_response_all);
if (delta_phi_deta3_response_all == 0) { cout<<"Response_all delta_phi_deta3 not found!"<<endl; return; }

TString delta_phi_deta3_miss_name = "ak5Gen_delta_phi_deta3_miss";
TH1D *delta_phi_deta3_miss = 0;
respfile->GetObject(delta_phi_deta3_miss_name,delta_phi_deta3_miss);
if (delta_phi_deta3_miss == 0) { cout<<"Miss delta_phi_deta3 histogram not found!"<<endl; return; }

TH2 *delta_phi_deta3_hResponse = (TH2*) delta_phi_deta3_response->Hresponse();
TH2 *delta_phi_deta3_hResponse_all = (TH2*) delta_phi_deta3_response_all->Hresponse();
TH1 *delta_phi_deta3_hFakes = (TH1D*) delta_phi_deta3_response->Hfakes();
TH1 *delta_phi_deta3_hTruth = (TH1D*) delta_phi_deta3_response->Htruth();

plot_response(delta_phi_deta3_hResponse, check_response_dir, prefix, "delta_phi_deta3_response", "top_left", false, detail, test);
plot_response(delta_phi_deta3_hResponse_all, check_response_dir, prefix, "delta_phi_deta3_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta3_hResponse, check_response_dir, prefix, "delta_phi_deta3_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta3_hResponse_all, check_response_dir, prefix, "delta_phi_deta3_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta3_hFakes, check_response_dir, prefix, "delta_phi_deta3_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta3_hTruth, check_response_dir, prefix, "delta_phi_deta3_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta3_miss, check_response_dir, prefix, "delta_phi_deta3_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta4 response on a safe way
TString delta_phi_deta4_hname = "resp_delta_phi_deta4";
RooUnfoldResponse *delta_phi_deta4_response = 0;
respfile->GetObject(delta_phi_deta4_hname,delta_phi_deta4_response);
if (delta_phi_deta4_response == 0) { cout<<"Response delta_phi_deta4 not found!"<<endl; return; }

TString delta_phi_deta4_hname_all = "resp_delta_phi_deta4_all";
RooUnfoldResponse *delta_phi_deta4_response_all = 0;
respfile->GetObject(delta_phi_deta4_hname_all,delta_phi_deta4_response_all);
if (delta_phi_deta4_response_all == 0) { cout<<"Response_all delta_phi_deta4 not found!"<<endl; return; }

TString delta_phi_deta4_miss_name = "ak5Gen_delta_phi_deta4_miss";
TH1D *delta_phi_deta4_miss = 0;
respfile->GetObject(delta_phi_deta4_miss_name,delta_phi_deta4_miss);
if (delta_phi_deta4_miss == 0) { cout<<"Miss delta_phi_deta4 histogram not found!"<<endl; return; }

TH2 *delta_phi_deta4_hResponse = (TH2*) delta_phi_deta4_response->Hresponse();
TH2 *delta_phi_deta4_hResponse_all = (TH2*) delta_phi_deta4_response_all->Hresponse();
TH1 *delta_phi_deta4_hFakes = (TH1D*) delta_phi_deta4_response->Hfakes();
TH1 *delta_phi_deta4_hTruth = (TH1D*) delta_phi_deta4_response->Htruth();

plot_response(delta_phi_deta4_hResponse, check_response_dir, prefix, "delta_phi_deta4_response", "top_left", false, detail, test);
plot_response(delta_phi_deta4_hResponse_all, check_response_dir, prefix, "delta_phi_deta4_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta4_hResponse, check_response_dir, prefix, "delta_phi_deta4_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta4_hResponse_all, check_response_dir, prefix, "delta_phi_deta4_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta4_hFakes, check_response_dir, prefix, "delta_phi_deta4_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta4_hTruth, check_response_dir, prefix, "delta_phi_deta4_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta4_miss, check_response_dir, prefix, "delta_phi_deta4_miss", "top_left", "Miss True", detail);

//get the delta_phi_gap response on a safe way
TString delta_phi_gap_hname = "resp_delta_phi_gap";
RooUnfoldResponse *delta_phi_gap_response = 0;
respfile->GetObject(delta_phi_gap_hname,delta_phi_gap_response);
if (delta_phi_gap_response == 0) { cout<<"Response delta_phi_gap not found!"<<endl; return; }

TString delta_phi_gap_hname_all = "resp_delta_phi_gap_all";
RooUnfoldResponse *delta_phi_gap_response_all = 0;
respfile->GetObject(delta_phi_gap_hname_all,delta_phi_gap_response_all);
if (delta_phi_gap_response_all == 0) { cout<<"Response_all delta_phi_gap not found!"<<endl; return; }

TString delta_phi_gap_miss_name = "ak5Gen_delta_phi_gap_miss";
TH1D *delta_phi_gap_miss = 0;
respfile->GetObject(delta_phi_gap_miss_name,delta_phi_gap_miss);
if (delta_phi_gap_miss == 0) { cout<<"Miss delta_phi_gap histogram not found!"<<endl; return; }

TH2 *delta_phi_gap_hResponse = (TH2*) delta_phi_gap_response->Hresponse();
TH2 *delta_phi_gap_hResponse_all = (TH2*) delta_phi_gap_response_all->Hresponse();
TH1 *delta_phi_gap_hFakes = (TH1D*) delta_phi_gap_response->Hfakes();
TH1 *delta_phi_gap_hTruth = (TH1D*) delta_phi_gap_response->Htruth();

plot_response(delta_phi_gap_hResponse, check_response_dir, prefix, "delta_phi_gap_response", "top_left", false, detail, test);
plot_response(delta_phi_gap_hResponse_all, check_response_dir, prefix, "delta_phi_gap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_gap_hResponse, check_response_dir, prefix, "delta_phi_gap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_gap_hResponse_all, check_response_dir, prefix, "delta_phi_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_gap_hFakes, check_response_dir, prefix, "delta_phi_gap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_gap_hTruth, check_response_dir, prefix, "delta_phi_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_gap_miss, check_response_dir, prefix, "delta_phi_gap_miss", "top_left", "Miss True", detail);

//get the delta_phi_deta1_gap response on a safe way
TString delta_phi_deta1_gap_hname = "resp_delta_phi_deta1_gap";
RooUnfoldResponse *delta_phi_deta1_gap_response = 0;
respfile->GetObject(delta_phi_deta1_gap_hname,delta_phi_deta1_gap_response);
if (delta_phi_deta1_gap_response == 0) { cout<<"Response delta_phi_deta1_gap not found!"<<endl; return; }

TString delta_phi_deta1_gap_hname_all = "resp_delta_phi_deta1_gap_all";
RooUnfoldResponse *delta_phi_deta1_gap_response_all = 0;
respfile->GetObject(delta_phi_deta1_gap_hname_all,delta_phi_deta1_gap_response_all);
if (delta_phi_deta1_gap_response_all == 0) { cout<<"Response_all delta_phi_deta1_gap not found!"<<endl; return; }

TString delta_phi_deta1_gap_miss_name = "ak5Gen_delta_phi_deta1_gap_miss";
TH1D *delta_phi_deta1_gap_miss = 0;
respfile->GetObject(delta_phi_deta1_gap_miss_name,delta_phi_deta1_gap_miss);
if (delta_phi_deta1_gap_miss == 0) { cout<<"Miss delta_phi_deta1_gap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta1_gap_hResponse = (TH2*) delta_phi_deta1_gap_response->Hresponse();
TH2 *delta_phi_deta1_gap_hResponse_all = (TH2*) delta_phi_deta1_gap_response_all->Hresponse();
TH1 *delta_phi_deta1_gap_hFakes = (TH1D*) delta_phi_deta1_gap_response->Hfakes();
TH1 *delta_phi_deta1_gap_hTruth = (TH1D*) delta_phi_deta1_gap_response->Htruth();

plot_response(delta_phi_deta1_gap_hResponse, check_response_dir, prefix, "delta_phi_deta1_gap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta1_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta1_gap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta1_gap_hResponse, check_response_dir, prefix, "delta_phi_deta1_gap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta1_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta1_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta1_gap_hFakes, check_response_dir, prefix, "delta_phi_deta1_gap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta1_gap_hTruth, check_response_dir, prefix, "delta_phi_deta1_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta1_gap_miss, check_response_dir, prefix, "delta_phi_deta1_gap_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta2_gap response on a safe way
TString delta_phi_deta2_gap_hname = "resp_delta_phi_deta2_gap";
RooUnfoldResponse *delta_phi_deta2_gap_response = 0;
respfile->GetObject(delta_phi_deta2_gap_hname,delta_phi_deta2_gap_response);
if (delta_phi_deta2_gap_response == 0) { cout<<"Response delta_phi_deta2_gap not found!"<<endl; return; }

TString delta_phi_deta2_gap_hname_all = "resp_delta_phi_deta2_gap_all";
RooUnfoldResponse *delta_phi_deta2_gap_response_all = 0;
respfile->GetObject(delta_phi_deta2_gap_hname_all,delta_phi_deta2_gap_response_all);
if (delta_phi_deta2_gap_response_all == 0) { cout<<"Response_all delta_phi_deta2_gap not found!"<<endl; return; }

TString delta_phi_deta2_gap_miss_name = "ak5Gen_delta_phi_deta2_gap_miss";
TH1D *delta_phi_deta2_gap_miss = 0;
respfile->GetObject(delta_phi_deta2_gap_miss_name,delta_phi_deta2_gap_miss);
if (delta_phi_deta2_gap_miss == 0) { cout<<"Miss delta_phi_deta2_gap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta2_gap_hResponse = (TH2*) delta_phi_deta2_gap_response->Hresponse();
TH2 *delta_phi_deta2_gap_hResponse_all = (TH2*) delta_phi_deta2_gap_response_all->Hresponse();
TH1 *delta_phi_deta2_gap_hFakes = (TH1D*) delta_phi_deta2_gap_response->Hfakes();
TH1 *delta_phi_deta2_gap_hTruth = (TH1D*) delta_phi_deta2_gap_response->Htruth();

plot_response(delta_phi_deta2_gap_hResponse, check_response_dir, prefix, "delta_phi_deta2_gap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta2_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta2_gap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta2_gap_hResponse, check_response_dir, prefix, "delta_phi_deta2_gap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta2_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta2_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta2_gap_hFakes, check_response_dir, prefix, "delta_phi_deta2_gap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta2_gap_hTruth, check_response_dir, prefix, "delta_phi_deta2_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta2_gap_miss, check_response_dir, prefix, "delta_phi_deta2_gap_miss", "top_left", "Miss True", detail);

//get the delta_phi_deta3_gap response on a safe way
TString delta_phi_deta3_gap_hname = "resp_delta_phi_deta3_gap";
RooUnfoldResponse *delta_phi_deta3_gap_response = 0;
respfile->GetObject(delta_phi_deta3_gap_hname,delta_phi_deta3_gap_response);
if (delta_phi_deta3_gap_response == 0) { cout<<"Response delta_phi_deta3_gap not found!"<<endl; return; }

TString delta_phi_deta3_gap_hname_all = "resp_delta_phi_deta3_gap_all";
RooUnfoldResponse *delta_phi_deta3_gap_response_all = 0;
respfile->GetObject(delta_phi_deta3_gap_hname_all,delta_phi_deta3_gap_response_all);
if (delta_phi_deta3_gap_response_all == 0) { cout<<"Response_all delta_phi_deta3_gap not found!"<<endl; return; }

TString delta_phi_deta3_gap_miss_name = "ak5Gen_delta_phi_deta3_gap_miss";
TH1D *delta_phi_deta3_gap_miss = 0;
respfile->GetObject(delta_phi_deta3_gap_miss_name,delta_phi_deta3_gap_miss);
if (delta_phi_deta3_gap_miss == 0) { cout<<"Miss delta_phi_deta3_gap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta3_gap_hResponse = (TH2*) delta_phi_deta3_gap_response->Hresponse();
TH2 *delta_phi_deta3_gap_hResponse_all = (TH2*) delta_phi_deta3_gap_response_all->Hresponse();
TH1 *delta_phi_deta3_gap_hFakes = (TH1D*) delta_phi_deta3_gap_response->Hfakes();
TH1 *delta_phi_deta3_gap_hTruth = (TH1D*) delta_phi_deta3_gap_response->Htruth();

plot_response(delta_phi_deta3_gap_hResponse, check_response_dir, prefix, "delta_phi_deta3_gap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta3_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta3_gap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta3_gap_hResponse, check_response_dir, prefix, "delta_phi_deta3_gap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta3_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta3_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta3_gap_hFakes, check_response_dir, prefix, "delta_phi_deta3_gap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta3_gap_hTruth, check_response_dir, prefix, "delta_phi_deta3_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta3_gap_miss, check_response_dir, prefix, "delta_phi_deta3_gap_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta4_gap response on a safe way
TString delta_phi_deta4_gap_hname = "resp_delta_phi_deta4_gap";
RooUnfoldResponse *delta_phi_deta4_gap_response = 0;
respfile->GetObject(delta_phi_deta4_gap_hname,delta_phi_deta4_gap_response);
if (delta_phi_deta4_gap_response == 0) { cout<<"Response delta_phi_deta4_gap not found!"<<endl; return; }

TString delta_phi_deta4_gap_hname_all = "resp_delta_phi_deta4_gap_all";
RooUnfoldResponse *delta_phi_deta4_gap_response_all = 0;
respfile->GetObject(delta_phi_deta4_gap_hname_all,delta_phi_deta4_gap_response_all);
if (delta_phi_deta4_gap_response_all == 0) { cout<<"Response_all delta_phi_deta4_gap not found!"<<endl; return; }

TString delta_phi_deta4_gap_miss_name = "ak5Gen_delta_phi_deta4_gap_miss";
TH1D *delta_phi_deta4_gap_miss = 0;
respfile->GetObject(delta_phi_deta4_gap_miss_name,delta_phi_deta4_gap_miss);
if (delta_phi_deta4_gap_miss == 0) { cout<<"Miss delta_phi_deta4_gap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta4_gap_hResponse = (TH2*) delta_phi_deta4_gap_response->Hresponse();
TH2 *delta_phi_deta4_gap_hResponse_all = (TH2*) delta_phi_deta4_gap_response_all->Hresponse();
TH1 *delta_phi_deta4_gap_hFakes = (TH1D*) delta_phi_deta4_gap_response->Hfakes();
TH1 *delta_phi_deta4_gap_hTruth = (TH1D*) delta_phi_deta4_gap_response->Htruth();

plot_response(delta_phi_deta4_gap_hResponse, check_response_dir, prefix, "delta_phi_deta4_gap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta4_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta4_gap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta4_gap_hResponse, check_response_dir, prefix, "delta_phi_deta4_gap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta4_gap_hResponse_all, check_response_dir, prefix, "delta_phi_deta4_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta4_gap_hFakes, check_response_dir, prefix, "delta_phi_deta4_gap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta4_gap_hTruth, check_response_dir, prefix, "delta_phi_deta4_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta4_gap_miss, check_response_dir, prefix, "delta_phi_deta4_gap_miss", "top_left", "Miss True", detail);


//get the delta_phi_nogap response on a safe way
TString delta_phi_nogap_hname = "resp_delta_phi_nogap";
RooUnfoldResponse *delta_phi_nogap_response = 0;
respfile->GetObject(delta_phi_nogap_hname,delta_phi_nogap_response);
if (delta_phi_nogap_response == 0) { cout<<"Response delta_phi_nogap not found!"<<endl; return; }

TString delta_phi_nogap_hname_all = "resp_delta_phi_nogap_all";
RooUnfoldResponse *delta_phi_nogap_response_all = 0;
respfile->GetObject(delta_phi_nogap_hname_all,delta_phi_nogap_response_all);
if (delta_phi_nogap_response_all == 0) { cout<<"Response_all delta_phi_nogap not found!"<<endl; return; }

TString delta_phi_nogap_miss_name = "ak5Gen_delta_phi_nogap_miss";
TH1D *delta_phi_nogap_miss = 0;
respfile->GetObject(delta_phi_nogap_miss_name,delta_phi_nogap_miss);
if (delta_phi_nogap_miss == 0) { cout<<"Miss delta_phi_nogap histogram not found!"<<endl; return; }

TH2 *delta_phi_nogap_hResponse = (TH2*) delta_phi_nogap_response->Hresponse();
TH2 *delta_phi_nogap_hResponse_all = (TH2*) delta_phi_nogap_response_all->Hresponse();
TH1 *delta_phi_nogap_hFakes = (TH1D*) delta_phi_nogap_response->Hfakes();
TH1 *delta_phi_nogap_hTruth = (TH1D*) delta_phi_nogap_response->Htruth();

plot_response(delta_phi_nogap_hResponse, check_response_dir, prefix, "delta_phi_nogap_response", "top_left", false, detail, test);
plot_response(delta_phi_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_nogap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_nogap_hResponse, check_response_dir, prefix, "delta_phi_nogap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_nogap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_nogap_hFakes, check_response_dir, prefix, "delta_phi_nogap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_nogap_hTruth, check_response_dir, prefix, "delta_phi_nogap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_nogap_miss, check_response_dir, prefix, "delta_phi_nogap_miss", "top_left", "Miss True", detail);

//get the delta_phi_deta1_nogap response on a safe way
TString delta_phi_deta1_nogap_hname = "resp_delta_phi_deta1_nogap";
RooUnfoldResponse *delta_phi_deta1_nogap_response = 0;
respfile->GetObject(delta_phi_deta1_nogap_hname,delta_phi_deta1_nogap_response);
if (delta_phi_deta1_nogap_response == 0) { cout<<"Response delta_phi_deta1_nogap not found!"<<endl; return; }

TString delta_phi_deta1_nogap_hname_all = "resp_delta_phi_deta1_nogap_all";
RooUnfoldResponse *delta_phi_deta1_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta1_nogap_hname_all,delta_phi_deta1_nogap_response_all);
if (delta_phi_deta1_nogap_response_all == 0) { cout<<"Response_all delta_phi_deta1_nogap not found!"<<endl; return; }

TString delta_phi_deta1_nogap_miss_name = "ak5Gen_delta_phi_deta1_nogap_miss";
TH1D *delta_phi_deta1_nogap_miss = 0;
respfile->GetObject(delta_phi_deta1_nogap_miss_name,delta_phi_deta1_nogap_miss);
if (delta_phi_deta1_nogap_miss == 0) { cout<<"Miss delta_phi_deta1_nogap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta1_nogap_hResponse = (TH2*) delta_phi_deta1_nogap_response->Hresponse();
TH2 *delta_phi_deta1_nogap_hResponse_all = (TH2*) delta_phi_deta1_nogap_response_all->Hresponse();
TH1 *delta_phi_deta1_nogap_hFakes = (TH1D*) delta_phi_deta1_nogap_response->Hfakes();
TH1 *delta_phi_deta1_nogap_hTruth = (TH1D*) delta_phi_deta1_nogap_response->Htruth();

plot_response(delta_phi_deta1_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta1_nogap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta1_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta1_nogap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta1_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta1_nogap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta1_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta1_nogap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta1_nogap_hFakes, check_response_dir, prefix, "delta_phi_deta1_nogap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta1_nogap_hTruth, check_response_dir, prefix, "delta_phi_deta1_nogap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta1_nogap_miss, check_response_dir, prefix, "delta_phi_deta1_nogap_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta2_nogap response on a safe way
TString delta_phi_deta2_nogap_hname = "resp_delta_phi_deta2_nogap";
RooUnfoldResponse *delta_phi_deta2_nogap_response = 0;
respfile->GetObject(delta_phi_deta2_nogap_hname,delta_phi_deta2_nogap_response);
if (delta_phi_deta2_nogap_response == 0) { cout<<"Response delta_phi_deta2_nogap not found!"<<endl; return; }

TString delta_phi_deta2_nogap_hname_all = "resp_delta_phi_deta2_nogap_all";
RooUnfoldResponse *delta_phi_deta2_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta2_nogap_hname_all,delta_phi_deta2_nogap_response_all);
if (delta_phi_deta2_nogap_response_all == 0) { cout<<"Response_all delta_phi_deta2_nogap not found!"<<endl; return; }

TString delta_phi_deta2_nogap_miss_name = "ak5Gen_delta_phi_deta2_nogap_miss";
TH1D *delta_phi_deta2_nogap_miss = 0;
respfile->GetObject(delta_phi_deta2_nogap_miss_name,delta_phi_deta2_nogap_miss);
if (delta_phi_deta2_nogap_miss == 0) { cout<<"Miss delta_phi_deta2_nogap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta2_nogap_hResponse = (TH2*) delta_phi_deta2_nogap_response->Hresponse();
TH2 *delta_phi_deta2_nogap_hResponse_all = (TH2*) delta_phi_deta2_nogap_response_all->Hresponse();
TH1 *delta_phi_deta2_nogap_hFakes = (TH1D*) delta_phi_deta2_nogap_response->Hfakes();
TH1 *delta_phi_deta2_nogap_hTruth = (TH1D*) delta_phi_deta2_nogap_response->Htruth();

plot_response(delta_phi_deta2_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta2_nogap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta2_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta2_nogap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta2_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta2_nogap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta2_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta2_nogap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta2_nogap_hFakes, check_response_dir, prefix, "delta_phi_deta2_nogap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta2_nogap_hTruth, check_response_dir, prefix, "delta_phi_deta2_nogap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta2_nogap_miss, check_response_dir, prefix, "delta_phi_deta2_nogap_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta3_nogap response on a safe way
TString delta_phi_deta3_nogap_hname = "resp_delta_phi_deta3_nogap";
RooUnfoldResponse *delta_phi_deta3_nogap_response = 0;
respfile->GetObject(delta_phi_deta3_nogap_hname,delta_phi_deta3_nogap_response);
if (delta_phi_deta3_nogap_response == 0) { cout<<"Response delta_phi_deta3_nogap not found!"<<endl; return; }

TString delta_phi_deta3_nogap_hname_all = "resp_delta_phi_deta3_nogap_all";
RooUnfoldResponse *delta_phi_deta3_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta3_nogap_hname_all,delta_phi_deta3_nogap_response_all);
if (delta_phi_deta3_nogap_response_all == 0) { cout<<"Response_all delta_phi_deta3_nogap not found!"<<endl; return; }

TString delta_phi_deta3_nogap_miss_name = "ak5Gen_delta_phi_deta3_nogap_miss";
TH1D *delta_phi_deta3_nogap_miss = 0;
respfile->GetObject(delta_phi_deta3_nogap_miss_name,delta_phi_deta3_nogap_miss);
if (delta_phi_deta3_nogap_miss == 0) { cout<<"Miss delta_phi_deta3_nogap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta3_nogap_hResponse = (TH2*) delta_phi_deta3_nogap_response->Hresponse();
TH2 *delta_phi_deta3_nogap_hResponse_all = (TH2*) delta_phi_deta3_nogap_response_all->Hresponse();
TH1 *delta_phi_deta3_nogap_hFakes = (TH1D*) delta_phi_deta3_nogap_response->Hfakes();
TH1 *delta_phi_deta3_nogap_hTruth = (TH1D*) delta_phi_deta3_nogap_response->Htruth();

plot_response(delta_phi_deta3_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta3_nogap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta3_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta3_nogap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta3_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta3_nogap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta3_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta3_nogap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta3_nogap_hFakes, check_response_dir, prefix, "delta_phi_deta3_nogap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta3_nogap_hTruth, check_response_dir, prefix, "delta_phi_deta3_nogap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta3_nogap_miss, check_response_dir, prefix, "delta_phi_deta3_nogap_miss", "top_left", "Miss True", detail);


//get the delta_phi_deta4_nogap response on a safe way
TString delta_phi_deta4_nogap_hname = "resp_delta_phi_deta4_nogap";
RooUnfoldResponse *delta_phi_deta4_nogap_response = 0;
respfile->GetObject(delta_phi_deta4_nogap_hname,delta_phi_deta4_nogap_response);
if (delta_phi_deta4_nogap_response == 0) { cout<<"Response delta_phi_deta4_nogap not found!"<<endl; return; }

TString delta_phi_deta4_nogap_hname_all = "resp_delta_phi_deta4_nogap_all";
RooUnfoldResponse *delta_phi_deta4_nogap_response_all = 0;
respfile->GetObject(delta_phi_deta4_nogap_hname_all,delta_phi_deta4_nogap_response_all);
if (delta_phi_deta4_nogap_response_all == 0) { cout<<"Response_all delta_phi_deta4_nogap not found!"<<endl; return; }

TString delta_phi_deta4_nogap_miss_name = "ak5Gen_delta_phi_deta4_nogap_miss";
TH1D *delta_phi_deta4_nogap_miss = 0;
respfile->GetObject(delta_phi_deta4_nogap_miss_name,delta_phi_deta4_nogap_miss);
if (delta_phi_deta4_nogap_miss == 0) { cout<<"Miss delta_phi_deta4_nogap histogram not found!"<<endl; return; }

TH2 *delta_phi_deta4_nogap_hResponse = (TH2*) delta_phi_deta4_nogap_response->Hresponse();
TH2 *delta_phi_deta4_nogap_hResponse_all = (TH2*) delta_phi_deta4_nogap_response_all->Hresponse();
TH1 *delta_phi_deta4_nogap_hFakes = (TH1D*) delta_phi_deta4_nogap_response->Hfakes();
TH1 *delta_phi_deta4_nogap_hTruth = (TH1D*) delta_phi_deta4_nogap_response->Htruth();

plot_response(delta_phi_deta4_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta4_nogap_response", "top_left", false, detail, test);
plot_response(delta_phi_deta4_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta4_nogap_response_all", "top_left", false, detail, test);
plot_response(delta_phi_deta4_nogap_hResponse, check_response_dir, prefix, "delta_phi_deta4_nogap_response_norm", "top_left", true, detail, test);
plot_response(delta_phi_deta4_nogap_hResponse_all, check_response_dir, prefix, "delta_phi_deta4_nogap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_phi_deta4_nogap_hFakes, check_response_dir, prefix, "delta_phi_deta4_nogap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_phi_deta4_nogap_hTruth, check_response_dir, prefix, "delta_phi_deta4_nogap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_phi_deta4_nogap_miss, check_response_dir, prefix, "delta_phi_deta4_nogap_miss", "top_left", "Miss True", detail);


//get the leading_pt_inside_gap response on a safe way
TString leading_pt_inside_gap_hname = "resp_leading_pt_inside_gap";
RooUnfoldResponse *leading_pt_inside_gap_response = 0;
respfile->GetObject(leading_pt_inside_gap_hname,leading_pt_inside_gap_response);
if (leading_pt_inside_gap_response == 0) { cout<<"Response leading_pt_inside_gap not found!"<<endl; return; }

TString leading_pt_inside_gap_hname_all = "resp_leading_pt_inside_gap_all";
RooUnfoldResponse *leading_pt_inside_gap_response_all = 0;
respfile->GetObject(leading_pt_inside_gap_hname_all,leading_pt_inside_gap_response_all);
if (leading_pt_inside_gap_response_all == 0) { cout<<"Response_all leading_pt_inside_gap not found!"<<endl; return; }

TString leading_pt_inside_gap_miss_name = "ak5Gen_leading_pt_inside_gap_miss";
TH1D *leading_pt_inside_gap_miss = 0;
respfile->GetObject(leading_pt_inside_gap_miss_name,leading_pt_inside_gap_miss);
if (leading_pt_inside_gap_miss == 0) { cout<<"Miss leading_pt_inside_gap histogram not found!"<<endl; return; }

TH2 *leading_pt_inside_gap_hResponse = (TH2*) leading_pt_inside_gap_response->Hresponse();
TH2 *leading_pt_inside_gap_hResponse_all = (TH2*) leading_pt_inside_gap_response_all->Hresponse();
TH1 *leading_pt_inside_gap_hFakes = (TH1D*) leading_pt_inside_gap_response->Hfakes();
TH1 *leading_pt_inside_gap_hTruth = (TH1D*) leading_pt_inside_gap_response->Htruth();

plot_response(leading_pt_inside_gap_hResponse, check_response_dir, prefix, "leading_pt_inside_gap_response", "top_left", false, detail, test);
plot_response(leading_pt_inside_gap_hResponse_all, check_response_dir, prefix, "leading_pt_inside_gap_response_all", "top_left", false, detail, test);
plot_response(leading_pt_inside_gap_hResponse, check_response_dir, prefix, "leading_pt_inside_gap_response_norm", "top_left", true, detail, test);
plot_response(leading_pt_inside_gap_hResponse_all, check_response_dir, prefix, "leading_pt_inside_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(leading_pt_inside_gap_hFakes, check_response_dir, prefix, "leading_pt_inside_gap_fakes", "top_left", "Fakes", detail);
plot_dist(leading_pt_inside_gap_hTruth, check_response_dir, prefix, "leading_pt_inside_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(leading_pt_inside_gap_miss, check_response_dir, prefix, "leading_pt_inside_gap_miss", "top_right", "Miss True", detail);


//get the leading_eta_star_inside_gap response on a safe way
TString leading_eta_star_inside_gap_hname = "resp_leading_eta_star_inside_gap";
RooUnfoldResponse *leading_eta_star_inside_gap_response = 0;
respfile->GetObject(leading_eta_star_inside_gap_hname,leading_eta_star_inside_gap_response);
if (leading_eta_star_inside_gap_response == 0) { cout<<"Response leading_eta_star_inside_gap not found!"<<endl; return; }

TString leading_eta_star_inside_gap_hname_all = "resp_leading_eta_star_inside_gap_all";
RooUnfoldResponse *leading_eta_star_inside_gap_response_all = 0;
respfile->GetObject(leading_eta_star_inside_gap_hname_all,leading_eta_star_inside_gap_response_all);
if (leading_eta_star_inside_gap_response_all == 0) { cout<<"Response_all leading_eta_star_inside_gap not found!"<<endl; return; }

TString leading_eta_star_inside_gap_miss_name = "ak5Gen_leading_eta_star_inside_gap_miss";
TH1D *leading_eta_star_inside_gap_miss = 0;
respfile->GetObject(leading_eta_star_inside_gap_miss_name,leading_eta_star_inside_gap_miss);
if (leading_eta_star_inside_gap_miss == 0) { cout<<"Miss leading_eta_star_inside_gap histogram not found!"<<endl; return; }

TH2 *leading_eta_star_inside_gap_hResponse = (TH2*) leading_eta_star_inside_gap_response->Hresponse();
TH2 *leading_eta_star_inside_gap_hResponse_all = (TH2*) leading_eta_star_inside_gap_response_all->Hresponse();
TH1 *leading_eta_star_inside_gap_hFakes = (TH1D*) leading_eta_star_inside_gap_response->Hfakes();
TH1 *leading_eta_star_inside_gap_hTruth = (TH1D*) leading_eta_star_inside_gap_response->Htruth();

plot_response(leading_eta_star_inside_gap_hResponse, check_response_dir, prefix, "leading_eta_star_inside_gap_response", "top_left", false, detail, test);
plot_response(leading_eta_star_inside_gap_hResponse_all, check_response_dir, prefix, "leading_eta_star_inside_gap_response_all", "top_left", false, detail, test);
plot_response(leading_eta_star_inside_gap_hResponse, check_response_dir, prefix, "leading_eta_star_inside_gap_response_norm", "top_left", true, detail, test);
plot_response(leading_eta_star_inside_gap_hResponse_all, check_response_dir, prefix, "leading_eta_star_inside_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(leading_eta_star_inside_gap_hFakes, check_response_dir, prefix, "leading_eta_star_inside_gap_fakes", "top_left", "Fakes", detail);
plot_dist(leading_eta_star_inside_gap_hTruth, check_response_dir, prefix, "leading_eta_star_inside_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(leading_eta_star_inside_gap_miss, check_response_dir, prefix, "leading_eta_star_inside_gap_miss", "top_right", "Miss True", detail);


//get the leading_pt_outside_gap response on a safe way
TString leading_pt_outside_gap_hname = "resp_leading_pt_outside_gap";
RooUnfoldResponse *leading_pt_outside_gap_response = 0;
respfile->GetObject(leading_pt_outside_gap_hname,leading_pt_outside_gap_response);
if (leading_pt_outside_gap_response == 0) { cout<<"Response leading_pt_outside_gap not found!"<<endl; return; }

TString leading_pt_outside_gap_hname_all = "resp_leading_pt_outside_gap_all";
RooUnfoldResponse *leading_pt_outside_gap_response_all = 0;
respfile->GetObject(leading_pt_outside_gap_hname_all,leading_pt_outside_gap_response_all);
if (leading_pt_outside_gap_response_all == 0) { cout<<"Response_all leading_pt_outside_gap not found!"<<endl; return; }

TString leading_pt_outside_gap_miss_name = "ak5Gen_leading_pt_outside_gap_miss";
TH1D *leading_pt_outside_gap_miss = 0;
respfile->GetObject(leading_pt_outside_gap_miss_name,leading_pt_outside_gap_miss);
if (leading_pt_outside_gap_miss == 0) { cout<<"Miss leading_pt_outside_gap histogram not found!"<<endl; return; }

TH2 *leading_pt_outside_gap_hResponse = (TH2*) leading_pt_outside_gap_response->Hresponse();
TH2 *leading_pt_outside_gap_hResponse_all = (TH2*) leading_pt_outside_gap_response_all->Hresponse();
TH1 *leading_pt_outside_gap_hFakes = (TH1D*) leading_pt_outside_gap_response->Hfakes();
TH1 *leading_pt_outside_gap_hTruth = (TH1D*) leading_pt_outside_gap_response->Htruth();

plot_response(leading_pt_outside_gap_hResponse, check_response_dir, prefix, "leading_pt_outside_gap_response", "top_left", false, detail, test);
plot_response(leading_pt_outside_gap_hResponse_all, check_response_dir, prefix, "leading_pt_outside_gap_response_all", "top_left", false, detail, test);
plot_response(leading_pt_outside_gap_hResponse, check_response_dir, prefix, "leading_pt_outside_gap_response_norm", "top_left", true, detail, test);
plot_response(leading_pt_outside_gap_hResponse_all, check_response_dir, prefix, "leading_pt_outside_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(leading_pt_outside_gap_hFakes, check_response_dir, prefix, "leading_pt_outside_gap_fakes", "top_left", "Fakes", detail);
plot_dist(leading_pt_outside_gap_hTruth, check_response_dir, prefix, "leading_pt_outside_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(leading_pt_outside_gap_miss, check_response_dir, prefix, "leading_pt_outside_gap_miss", "top_right", "Miss True", detail);


//get the delta_eta_outside_gap response on a safe way
TString delta_eta_outside_gap_hname = "resp_delta_eta_outside_gap";
RooUnfoldResponse *delta_eta_outside_gap_response = 0;
respfile->GetObject(delta_eta_outside_gap_hname,delta_eta_outside_gap_response);
if (delta_eta_outside_gap_response == 0) { cout<<"Response delta_eta_outside_gap not found!"<<endl; return; }

TString delta_eta_outside_gap_hname_all = "resp_delta_eta_outside_gap_all";
RooUnfoldResponse *delta_eta_outside_gap_response_all = 0;
respfile->GetObject(delta_eta_outside_gap_hname_all,delta_eta_outside_gap_response_all);
if (delta_eta_outside_gap_response_all == 0) { cout<<"Response_all delta_eta_outside_gap not found!"<<endl; return; }

TString delta_eta_outside_gap_miss_name = "ak5Gen_delta_eta_outside_gap_miss";
TH1D *delta_eta_outside_gap_miss = 0;
respfile->GetObject(delta_eta_outside_gap_miss_name,delta_eta_outside_gap_miss);
if (delta_eta_outside_gap_miss == 0) { cout<<"Miss delta_eta_outside_gap histogram not found!"<<endl; return; }

TH2 *delta_eta_outside_gap_hResponse = (TH2*) delta_eta_outside_gap_response->Hresponse();
TH2 *delta_eta_outside_gap_hResponse_all = (TH2*) delta_eta_outside_gap_response_all->Hresponse();
TH1 *delta_eta_outside_gap_hFakes = (TH1D*) delta_eta_outside_gap_response->Hfakes();
TH1 *delta_eta_outside_gap_hTruth = (TH1D*) delta_eta_outside_gap_response->Htruth();

plot_response(delta_eta_outside_gap_hResponse, check_response_dir, prefix, "delta_eta_outside_gap_response", "top_left", false, detail, test);
plot_response(delta_eta_outside_gap_hResponse_all, check_response_dir, prefix, "delta_eta_outside_gap_response_all", "top_left", false, detail, test);
plot_response(delta_eta_outside_gap_hResponse, check_response_dir, prefix, "delta_eta_outside_gap_response_norm", "top_left", true, detail, test);
plot_response(delta_eta_outside_gap_hResponse_all, check_response_dir, prefix, "delta_eta_outside_gap_response_all_norm", "top_left", true, detail, test);
plot_dist(delta_eta_outside_gap_hFakes, check_response_dir, prefix, "delta_eta_outside_gap_fakes", "top_left", "Fakes", detail);
plot_dist(delta_eta_outside_gap_hTruth, check_response_dir, prefix, "delta_eta_outside_gap_truth_from_response", "top_left", "Truth from response matrix", detail);
plot_dist(delta_eta_outside_gap_miss, check_response_dir, prefix, "delta_eta_outside_gap_miss", "top_right", "Miss True", detail);

//delete the variables to avoid memory leak

if (detail) { cout << "Done!" << endl; }
}
