// Pedro Cipriano, Dec 2011
// DESY, CMS
// Last Update: 24 Oct 2012
//
// estimate_pileup_uncertainty()
// calculate the pileup uncertianty

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
#include <TF1.h>

#include <TChain.h>
#include <TChainElement.h>

#include <iostream>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void alisar(TH1D *histo_1, TH1D *histo_2, double max_dif)
{

int max_bins = histo_1->GetNbinsX();

    for(Int_t i=1; i<=max_bins-1;i++)
    {
    Float_t val1 = histo_1->GetBinContent(i);
    Float_t val2 = histo_1->GetBinContent(i+1);
    Float_t diff = val1 - val2;
    if (diff < 0) {diff = - diff; }
    if (diff > max_dif)
	{
	Float_t ave = (val1 + val2)/2;
	histo_2->SetBinContent(i,ave);
	histo_2->SetBinContent(i+1,ave);
	i = i + 1;
	}
	else
	{
	histo_2->SetBinContent(i,val1);
	}
    cout<<i<<" ="<<val1<<" . "<<val2<<" . "<<diff<<endl;
    }

if (histo_2->GetBinContent(max_bins) == 0)
{ 
histo_2->SetBinContent(max_bins,histo_1->GetBinContent(max_bins));
}

}


void number_of_vertices_ratio()
{

gROOT->Reset();
gROOT->SetStyle("Plain");

    TFile *data_file1a = new TFile("output/histograms/raw_xsec_data/xsec_JetMETTau2010A_allvertex.root");
    TString labela = "JetMETTau 2010A"; 

    TFile *data_file1b = new TFile("output/histograms/raw_xsec_data/xsec_JetMET2010A_allvertex.root");
    TString labelb = "JetMET 2010A";

    TFile *data_file1c = new TFile("output/histograms/raw_xsec_data/xsec_Jet2010B_allvertex.root");
    TString labelc = "Jet 2010B";

    TFile *data_file2a = new TFile("output/histograms/xsec_official_mc/xsec_Herwig6_allvertex.root");
    TString label1 = "Herwig6"; 

    TFile *data_file2b = new TFile("output/histograms/xsec_official_mc/xsec_Pythia6_D6T_allvertex.root");
    TString label2 = "Pythia6 D6T";

    TFile *data_file2c = new TFile("output/histograms/xsec_official_mc/xsec_Pythia6_Z2_allvertex.root");
    TString label3 = "Pythia6 Z2";
    
    TFile *data_file2d = new TFile("output/histograms/xsec_official_mc/xsec_Pythia8_1_allvertex.root");
    TString label4 = "Pythia8 1";

    TH1D *x1a = (TH1D*) data_file1a->Get("vertex_selected");
    TH1D *x1b = (TH1D*) data_file1b->Get("vertex_selected");
    TH1D *x1c = (TH1D*) data_file1c->Get("vertex_selected");
    TH1D *x2a = (TH1D*) data_file2a->Get("vertex_selected");
    TH1D *x2b = (TH1D*) data_file2b->Get("vertex_selected");
    TH1D *x2c = (TH1D*) data_file2c->Get("vertex_selected");
    TH1D *x2d = (TH1D*) data_file2d->Get("vertex_selected");
    
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetLogy();
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    TH1D *ratio1a;
    ratio1a =  new TH1D("Number of vertices ratio_a_a","JetMETTau_2010A / Herwig6;#Vertices;Ratio", 10,0,10);
    ratio1a->Divide(x1a,x2a,1.,1.,"");
    ratio1a->SetLineColor(2);
    ratio1a->SetLineStyle(1);
    ratio1a->SetLineWidth(3);
    ratio1a->SetMaximum(1e6);
    ratio1a->Draw("hist");
    
    TH1D *ratio2a;
    ratio2a =  new TH1D("Number of vertices ratio_a_b","JetMETTau_2010A / Pythia6 D6T;#Vertices;Ratio", 10,0,10);
    ratio2a->Divide(x1a,x2b,1.,1.,"");
    ratio2a->SetLineColor(3);
    ratio2a->SetLineStyle(2);
    ratio2a->SetLineWidth(3);
    ratio2a->Draw("hist same");
    
    TH1D *ratio3a;
    ratio3a =  new TH1D("Number of vertices ratio_a_c","JetMETTau_2010A / Pythia6 Z2;#Vertices;Ratio", 10,0,10);
    ratio3a->Divide(x1a,x2c,1.,1.,"");
    ratio3a->SetLineColor(4);
    ratio3a->SetLineStyle(3);
    ratio3a->SetLineWidth(3);
    ratio3a->Draw("hist same");
    
    TH1D *ratio4a;
    ratio4a =  new TH1D("Number of vertices ratio_a_d","JetMETTau_2010A / Pythia8 1;#Vertices;Ratio", 10,0,10);
    ratio4a->Divide(x1a,x2d,1.,1.,"");
    ratio4a->SetLineColor(6);
    ratio4a->SetLineStyle(4);
    ratio4a->SetLineWidth(3);
    ratio4a->Draw("hist same");
    
    TLegend *leg1 = new TLegend(0.6,0.6,0.97,0.97);
    leg1->AddEntry(x1a,labela,"");
    leg1->AddEntry(ratio1a,label1,"l");
    leg1->AddEntry(ratio2a,label2,"l");
    leg1->AddEntry(ratio3a,label3,"l");
    leg1->AddEntry(ratio4a,label4,"l");
    leg1->SetFillColor(0);
    leg1->Draw();
    
    c1->Print("output/pile_up_studies/vertex_ratio_JetMETTAU_2010A_log.png");
    c1->Close();
    
    TCanvas *c1b = new TCanvas("c1b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    ratio1a->Draw("hist");
    ratio2a->Draw("hist same");
    ratio3a->Draw("hist same");
    ratio4a->Draw("hist same");
    
    leg1->Draw();
    
    c1b->Print("output/pile_up_studies/vertex_ratio_JetMETTAU_2010A.png");
    c1b->Close();
    
    TCanvas *c2 = new TCanvas("c2","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetLogy();
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    TH1D *ratio1b;
    ratio1b =  new TH1D("Number of vertices ratio_b_a","JetMETTau_2010A / Herwig6;#Vertices;Ratio", 10,0,10);
    ratio1b->Divide(x1b,x2a,1.,1.,"");
    ratio1b->SetLineColor(2);
    ratio1b->SetLineStyle(1);
    ratio1b->SetLineWidth(3);
    ratio1b->SetMaximum(1e6);
    ratio1b->Draw("hist");
    
    TH1D *ratio2b;
    ratio2b =  new TH1D("Number of vertices ratio_b_b","JetMETTau_2010A / Pythia6 D6T;#Vertices;Ratio", 10,0,10);
    ratio2b->Divide(x1b,x2b,1.,1.,"");
    ratio2b->SetLineColor(3);
    ratio2b->SetLineStyle(2);
    ratio2b->SetLineWidth(3);
    ratio2b->Draw("hist same");
    
    TH1D *ratio3b;
    ratio3b =  new TH1D("Number of vertices ratio_b_c","JetMETTau_2010A / Pythia6 Z2;#Vertices;Ratio", 10,0,10);
    ratio3b->Divide(x1b,x2c,1.,1.,"");
    ratio3b->SetLineColor(4);
    ratio3b->SetLineStyle(3);
    ratio3b->SetLineWidth(3);
    ratio3b->Draw("hist same");
    
    TH1D *ratio4b;
    ratio4b =  new TH1D("Number of vertices ratio_b_d","JetMETTau_2010A / Pythia8 1;#Vertices;Ratio", 10,0,10);
    ratio4b->Divide(x1b,x2d,1.,1.,"");
    ratio4b->SetLineColor(6);
    ratio4b->SetLineStyle(4);
    ratio4b->SetLineWidth(3);
    ratio4b->Draw("hist same");
    
    TLegend *leg2 = new TLegend(0.6,0.6,0.97,0.97);
    leg2->AddEntry(x1b,labelb,"");
    leg2->AddEntry(ratio1b,label1,"l");
    leg2->AddEntry(ratio2b,label2,"l");
    leg2->AddEntry(ratio3b,label3,"l");
    leg2->AddEntry(ratio4b,label4,"l");
    leg2->SetFillColor(0);
    leg2->Draw();
    
    c2->Print("output/pile_up_studies/vertex_ratio_JetMET_2010A_log.png");
    c2->Close();
    
    TCanvas *c2b = new TCanvas("c2b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    ratio1b->Draw("hist");
    ratio2b->Draw("hist same");
    ratio3b->Draw("hist same");
    ratio4b->Draw("hist same");
    
    leg2->Draw();
    
    c2b->Print("output/pile_up_studies/vertex_ratio_JetMET_2010A.png");
    c2b->Close();
    
    TCanvas *c3 = new TCanvas("c3","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetLogy();
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    TH1D *ratio1c;
    ratio1c =  new TH1D("Number of vertices ratio_c_a","JetMETTau_2010A / Herwig6;#Vertices;Ratio", 10,0,10);
    ratio1c->Divide(x1c,x2a,1.,1.,"");
    ratio1c->SetLineColor(2);
    ratio1c->SetLineStyle(1);
    ratio1c->SetLineWidth(3);
    ratio1c->SetMaximum(1e6);
    ratio1c->Draw("hist");
    
    TH1D *ratio2c;
    ratio2c =  new TH1D("Number of vertices ratio_c_b","JetMETTau_2010A / Pythia6 D6T;#Vertices;Ratio", 10,0,10);
    ratio2c->Divide(x1c,x2b,1.,1.,"");
    ratio2c->SetLineColor(3);
    ratio2c->SetLineStyle(2);
    ratio2c->SetLineWidth(3);
    ratio2c->Draw("hist same");
    
    TH1D *ratio3c;
    ratio3c =  new TH1D("Number of vertices ratio_c_c","JetMETTau_2010A / Pythia6 Z2;#Vertices;Ratio", 10,0,10);
    ratio3c->Divide(x1c,x2c,1.,1.,"");
    ratio3c->SetLineColor(4);
    ratio3c->SetLineStyle(3);
    ratio3c->SetLineWidth(3);
    ratio3c->Draw("hist same");
    
    TH1D *ratio4c;
    ratio4c =  new TH1D("Number of vertices ratio_c_d","JetMETTau_2010A / Pythia8 1;#Vertices;Ratio", 10,0,10);
    ratio4c->Divide(x1c,x2d,1.,1.,"");
    ratio4c->SetLineColor(6);
    ratio4c->SetLineStyle(4);
    ratio4c->SetLineWidth(3);
    ratio4c->Draw("hist same");
    
    TLegend *leg3 = new TLegend(0.6,0.6,0.97,0.97);
    leg3->AddEntry(x1c,labelc,"");
    leg3->AddEntry(ratio1c,label1,"l");
    leg3->AddEntry(ratio2c,label2,"l");
    leg3->AddEntry(ratio3c,label3,"l");
    leg3->AddEntry(ratio4c,label4,"l");
    leg3->SetFillColor(0);
    leg3->Draw();
    
    c3->Print("output/pile_up_studies/vertex_ratio_Jet_2010B_log.png");
    c3->Close();
    
    TCanvas *c3b = new TCanvas("c3b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    ratio1c->Draw("hist");
    ratio2c->Draw("hist same");
    ratio3c->Draw("hist same");
    ratio4c->Draw("hist same");
    
    leg3->Draw();
    
    c3b->Print("output/pile_up_studies/vertex_ratio_Jet_2010B.png");
    c3b->Close();

}


void estimate_pileup_uncertainty(string output_path = "../output/pileup_uncertainty/", bool detail = false, bool disp_uncertainty = true)
{
//main estimate pileup uncertainty routine


//setting the root style
gROOT->Reset();
gROOT->SetStyle("Plain");

//histogram bins
int cent_nbins = 7;
int forw_nbins = 7;

double cent_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};
double forw_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};

int in_nbins = 11;
int out_nbins = 11;

double in_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};
double out_bins[12] = {10, 15, 20, 27, 35, 45, 57, 72, 90, 120, 150, 200};

int dphi_nbins = 7;
double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

int etastar_nbins = 12;
double etastar_bins[13] = {-3.6, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.6};

int deta_out_nbins = 6;
double deta_out_bins[7] = {0, 1, 2, 3, 4, 5, 7.5};

    TFile *data_file1a = new TFile("../output/histograms/raw_xsec_data/xsec_JetMETTau2010A_allvertex.root");
    TString label1a = "Data JetMETTau 2010A - All Vertex";  
    TFile *data_file2a = new TFile("../output/histograms/raw_xsec_data/xsec_JetMETTau2010A_1vertex.root");
    TString label2a = "Data JetMETTau 2010A - 1 Vertex";
    TString labela = "Data JetMETTau 2010A"; 

    TFile *data_file1b = new TFile("../output/histograms/raw_xsec_data/xsec_JetMET2010A_allvertex.root");
    TString label1b = "Data JetMET 2010A - All Vertex";
    TFile *data_file2b = new TFile("../output/histograms/raw_xsec_data/xsec_JetMET2010A_1vertex.root");
    TString label2b = "Data JetMET 2010A - 1 Vertex";
    TString labelb = "Data JetMET 2010A";

    TFile *data_file1c = new TFile("../output/histograms/raw_xsec_data/xsec_Jet2010B_allvertex.root");
    TString label1c = "Data Jet 2010B - All Vertex";
    TFile *data_file2c = new TFile("../output/histograms/raw_xsec_data/xsec_Jet2010B_1vertex.root");
    TString label2c = "Data Jet 2010B - 1 Vertex";
    TString labelc = "Data Jet 2010B";
    
    TFile data_output("../output/histograms/pile_up_uncertainty_2010.root", "RECREATE"); 

    double int1a, int1b, int2a, int2b, int1c, int2c;
    double factora, factorb, factorc;
    double enta, entb, entc;
    double min, max, tot, ave, abs_cont, tot_old, ave_old;

    double min_row[10], max_row[10], ave_row[10], old_row[10];

    double step = 0.05;

    factora = 1.0;
    factorb = 1.0;
    factorc = 1.0;
    
    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x1a = (TH1D*) data_file1a->Get("ak5PF_delta_phi");
     TH1D *x2a = (TH1D*) data_file2a->Get("ak5PF_delta_phi");
     TH1D *x1b = (TH1D*) data_file1b->Get("ak5PF_delta_phi");
     TH1D *x2b = (TH1D*) data_file2b->Get("ak5PF_delta_phi");
     TH1D *x1c = (TH1D*) data_file1c->Get("ak5PF_delta_phi");
     TH1D *x2c = (TH1D*) data_file2c->Get("ak5PF_delta_phi");
     
     enta = x1a->GetEntries()-x1a->GetNbinsX();
     entb = x1b->GetEntries()-x1b->GetNbinsX();
     entc = x1c->GetEntries()-x1c->GetNbinsX();

     int1a = x1a->Integral();
     int2a = x2a->Integral();
     factora = int2a/int1a;

     int1b = x1b->Integral();
     int2b = x2b->Integral();
     factorb = int2b/int1b;

     int1c = x1c->Integral();
     int2c = x2c->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi, *delta_phi2, *delta_phi_old, *delta_phi_a, *delta_phi_b, *delta_phi_c;
    delta_phi =  new TH1D("Delta_phi","delta_phi;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi2 =  new TH1D("Delta_phi2","delta_phi2;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_old =  new TH1D("Delta_phi_old","delta_phi_old;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_a =  new TH1D("Delta_phi_a","delta_phi_a;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_b =  new TH1D("Delta_phi_b","delta_phi_b;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_c =  new TH1D("Delta_phi_c","delta_phi_c;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_a->Divide(x1a,x2a,1.,1./factora,"");
    delta_phi_b->Divide(x1b,x2b,1.,1./factorb,"");
    delta_phi_c->Divide(x1c,x2c,1.,1./factorc,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_phi_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    //cout<<i<<" ="<<conta<<" . "<<contb<<" . "<<contc<<endl;
    delta_phi_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi2, delta_phi, step);

    delta_phi->SetMinimum(0.0);
    delta_phi->SetMaximum(0.2);
    delta_phi->SetLineColor(2);
    delta_phi->SetLineStyle(1);
    delta_phi->SetLineWidth(3);
    delta_phi->Draw("hist");

    TLegend *leg1 = new TLegend(0.5,0.9,0.97,0.97);
    leg1->AddEntry(delta_phi,"Average Pile-up uncertainty for 2010 Data","l");
    leg1->SetFillColor(0);
    leg1->Draw();
    
    ave = tot/delta_phi->GetNbinsX();
    ave_old = tot_old/delta_phi->GetNbinsX();
    cout<<"Main selection: min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[0] = min*100;
    max_row[0] = max*100;
    ave_row[0] = ave*100;
    old_row[0] = ave_old*100;
    
    c1->Print("output/uncertainities/png/pile_up_delta_phi.png");
    c1->Print("output/uncertainities/c/pile_up_delta_phi.C");
    c1->Print("output/uncertainities/eps/pile_up_delta_phi.eps");
    c1->Close();


    TCanvas *c1o = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    delta_phi->SetMinimum(0.0);
    delta_phi->SetMaximum(0.2);
    delta_phi->Draw("hist");
    delta_phi2->SetLineColor(2);
    delta_phi2->SetLineStyle(3);
    delta_phi2->SetLineWidth(3);
    delta_phi2->Draw("hist same");
    delta_phi_old->SetLineColor(2);
    delta_phi_old->SetLineStyle(2);
    delta_phi_old->SetLineWidth(3);
    delta_phi_old->Draw("hist same");

    TLegend *leg1o = new TLegend(0.5,0.84,0.97,0.97);
    leg1o->AddEntry(delta_phi,"Average Pile-up uncertainty for 2010 Data","l");
    leg1o->AddEntry(delta_phi2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg1o->AddEntry(delta_phi_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg1o->SetFillColor(0);
    leg1o->Draw();

    c1o->Print("output/uncertainities/png/pile_up_delta_phi_a.png");
    c1o->Print("output/uncertainities/c/pile_up_delta_phi_a.C");
    c1o->Print("output/uncertainities/eps/pile_up_delta_phi_a.eps");
    c1o->Close();

    TCanvas *c1b = new TCanvas("c1b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    delta_phi_a->SetMaximum(1.6);
    delta_phi_a->SetMinimum(0.6);
    delta_phi_a->SetLineColor(2);
    delta_phi_a->SetLineStyle(1);
    delta_phi_a->SetLineWidth(3);
    delta_phi_a->Draw("e1");
    delta_phi_b->SetLineColor(3);
    delta_phi_b->SetLineStyle(1);
    delta_phi_b->SetLineWidth(3);
    delta_phi_b->Draw("e1 same");
    delta_phi_c->SetLineColor(4);
    delta_phi_c->SetLineStyle(1);
    delta_phi_c->SetLineWidth(3);
    delta_phi_c->Draw("e1 same");

    TLegend *leg1b = new TLegend(0.5,0.84,0.97,0.97);
    leg1b->AddEntry(delta_phi_a,labela,"l");
    leg1b->AddEntry(delta_phi_b,labelb,"l");
    leg1b->AddEntry(delta_phi_c,labelc,"l");
    leg1b->SetFillColor(0);
    leg1b->Draw();
    
    c1b->Print("output/uncertainities/png/pile_up_delta_phi_b.png");
    c1b->Print("output/uncertainities/c/pile_up_delta_phi_b.C");
    c1b->Print("output/uncertainities/eps/pile_up_delta_phi_b.eps");
    c1b->Close();
    
    
    TCanvas *c2 = new TCanvas("c2","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x1ag = (TH1D*) data_file1a->Get("ak5PF_delta_phi_gap");
     TH1D *x2ag = (TH1D*) data_file2a->Get("ak5PF_delta_phi_gap");
     TH1D *x1bg = (TH1D*) data_file1b->Get("ak5PF_delta_phi_gap");
     TH1D *x2bg = (TH1D*) data_file2b->Get("ak5PF_delta_phi_gap");
     TH1D *x1cg = (TH1D*) data_file1c->Get("ak5PF_delta_phi_gap");
     TH1D *x2cg = (TH1D*) data_file2c->Get("ak5PF_delta_phi_gap");
     
     enta = x1ag->GetEntries()-x1ag->GetNbinsX();
     entb = x1bg->GetEntries()-x1bg->GetNbinsX();
     entc = x1cg->GetEntries()-x1cg->GetNbinsX();

     int1a = x1ag->Integral();
     int2a = x2ag->Integral();
     factora = int2a/int1a;

     int1b = x1bg->Integral();
     int2b = x2bg->Integral();
     factorb = int2b/int1b;

     int1c = x1cg->Integral();
     int2c = x2cg->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_gap, *delta_phi_gap2, *delta_phi_gap_old, *delta_phi_gap_a, *delta_phi_gap_b, *delta_phi_gap_c;
    delta_phi_gap =  new TH1D("Delta_phi_gap","delta_phi_gap;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_gap2 =  new TH1D("Delta_phi_gap2","delta_phi_gap2;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_gap_old =  new TH1D("Delta_phi_gap_old","delta_phi_gap_old;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_gap_a =  new TH1D("Delta_phi_gap_a","delta_phi_gap_a;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_gap_b =  new TH1D("Delta_phi_gap_b","delta_phi_gap_b;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_gap_c =  new TH1D("Delta_phi_gap_c","delta_phi_gap_c;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_gap_a->Divide(x1ag,x2ag,1.,1./factora,"");
    delta_phi_gap_b->Divide(x1bg,x2bg,1.,1./factorb,"");
    delta_phi_gap_c->Divide(x1cg,x2cg,1.,1./factorc,"");

    min = 0;
    max = 0;
    tot = 0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_phi_gap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_gap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_gap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    //cout<<i<<" ="<<conta<<" . "<<contb<<" . "<<contc<<endl;
    delta_phi_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_gap2, delta_phi_gap, step);

    delta_phi_gap->SetMinimum(0.0);
    delta_phi_gap->SetMaximum(0.25);
    delta_phi_gap->SetLineColor(2);
    delta_phi_gap->SetLineStyle(1);
    delta_phi_gap->SetLineWidth(3);
    delta_phi_gap->Draw("hist");
    
    TLegend *leg2 = new TLegend(0.5,0.9,0.97,0.97);
    leg2->AddEntry(delta_phi_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg2->SetFillColor(0);
    leg2->Draw();

    ave = tot/delta_phi_gap->GetNbinsX();
    ave_old = tot_old/delta_phi->GetNbinsX();
    cout<<"Gap: min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[2] = min*100;
    max_row[2] = max*100;
    ave_row[2] = ave*100;
    old_row[2] = ave_old*100;   

    c2->Print("output/uncertainities/png/pile_up_delta_phi_gap.png");
    c2->Print("output/uncertainities/c/pile_up_delta_phi_gap.C");
    c2->Print("output/uncertainities/eps/pile_up_delta_phi_gap.eps");
    c2->Close();
    
    TCanvas *c2o = new TCanvas("c2b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    delta_phi_gap->SetMinimum(0.0);
    delta_phi_gap->SetMaximum(0.25);
    delta_phi_gap->Draw("hist");
    delta_phi_gap2->SetLineColor(2);
    delta_phi_gap2->SetLineStyle(3);
    delta_phi_gap2->SetLineWidth(3);
    delta_phi_gap2->Draw("hist same");
    delta_phi_gap_old->SetLineColor(2);
    delta_phi_gap_old->SetLineStyle(2);
    delta_phi_gap_old->SetLineWidth(3);
    delta_phi_gap_old->Draw("hist same");

    TLegend *leg2o = new TLegend(0.5,0.84,0.97,0.97);
    leg2o->AddEntry(delta_phi_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg2o->AddEntry(delta_phi_gap2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg2o->AddEntry(delta_phi_gap_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg2o->SetFillColor(0);
    leg2o->Draw();

    c2o->Print("output/uncertainities/png/pile_up_delta_phi_gap_a.png");
    c2o->Print("output/uncertainities/c/pile_up_delta_phi_gap_a.C");
    c2o->Print("output/uncertainities/eps/pile_up_delta_phi_gap_a.eps");
    c2o->Close();

    TCanvas *c2b = new TCanvas("c2b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    delta_phi_gap_a->SetMaximum(1.6);
    delta_phi_gap_a->SetMinimum(0.4);
    delta_phi_gap_a->SetLineColor(2);
    delta_phi_gap_a->SetLineStyle(1);
    delta_phi_gap_a->SetLineWidth(3);
    delta_phi_gap_a->Draw("e1");
    delta_phi_gap_b->SetLineColor(3);
    delta_phi_gap_b->SetLineStyle(1);
    delta_phi_gap_b->SetLineWidth(3);
    delta_phi_gap_b->Draw("e1 same");
    delta_phi_gap_c->SetLineColor(4);
    delta_phi_gap_c->SetLineStyle(1);
    delta_phi_gap_c->SetLineWidth(3);
    delta_phi_gap_c->Draw("e1 same");

    TLegend *leg2b = new TLegend(0.5,0.84,0.97,0.97);
    leg2b->AddEntry(delta_phi_gap_a,labela,"l");
    leg2b->AddEntry(delta_phi_gap_b,labelb,"l");
    leg2b->AddEntry(delta_phi_gap_c,labelc,"l");
    leg2b->SetFillColor(0);
    leg2b->Draw();
    
    c2b->Print("output/uncertainities/png/pile_up_delta_phi_gap_b.png");
    c2b->Print("output/uncertainities/c/pile_up_delta_phi_gap_b.C");
    c2b->Print("output/uncertainities/eps/pile_up_delta_phi_gap_b.eps");
    c2b->Close();
    
    
    TCanvas *c3 = new TCanvas("c3","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x1an = (TH1D*) data_file1a->Get("ak5PF_delta_phi_nogap");
     TH1D *x2an = (TH1D*) data_file2a->Get("ak5PF_delta_phi_nogap");
     TH1D *x1bn = (TH1D*) data_file1b->Get("ak5PF_delta_phi_nogap");
     TH1D *x2bn = (TH1D*) data_file2b->Get("ak5PF_delta_phi_nogap");
     TH1D *x1cn = (TH1D*) data_file1c->Get("ak5PF_delta_phi_nogap");
     TH1D *x2cn = (TH1D*) data_file2c->Get("ak5PF_delta_phi_nogap");
     
     enta = x1an->GetEntries()-x1an->GetNbinsX();
     entb = x1bn->GetEntries()-x1bn->GetNbinsX();
     entc = x1cn->GetEntries()-x1cn->GetNbinsX();

     int1a = x1an->Integral();
     int2a = x2an->Integral();
     factora = int2a/int1a;

     int1b = x1bn->Integral();
     int2b = x2bn->Integral();
     factorb = int2b/int1b;

     int1c = x1cn->Integral();
     int2c = x2cn->Integral();
     factorc = int2c/int1c;
     
     cout<<"factor"<<factora<<endl;
    
    TH1D *delta_phi_nogap, *delta_phi_nogap2, *delta_phi_nogap_old, *delta_phi_nogap_a, *delta_phi_nogap_b, *delta_phi_nogap_c;
    delta_phi_nogap =  new TH1D("Delta_phi_nogap","delta_phi_nogap;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_nogap2 =  new TH1D("Delta_phi_nogap2","delta_phi_nogap2;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_nogap_old =  new TH1D("Delta_phi_nogap_old","delta_phi_nogap_old;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_nogap_a =  new TH1D("Delta_phi_nogap_a","delta_phi_nogap_a;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_nogap_b =  new TH1D("Delta_phi_nogap_b","delta_phi_nogap_b;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_nogap_c =  new TH1D("Delta_phi_nogap_c","delta_phi_nogap_c;#Delta#phi [rad];Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_nogap_a->Divide(x1an,x2an,1.,1./factora,"");
    delta_phi_nogap_b->Divide(x1bn,x2bn,1.,1./factorb,"");
    delta_phi_nogap_c->Divide(x1cn,x2cn,1.,1./factorc,"");

    min = 0;
    max = 0;
    tot = 0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_phi_nogap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_nogap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_nogap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_nogap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_nogap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_nogap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_nogap2, delta_phi_nogap, step);

    delta_phi_nogap->SetMinimum(0.0);
    delta_phi_nogap->SetMaximum(0.25);
    delta_phi_nogap->SetLineColor(2);
    delta_phi_nogap->SetLineStyle(1);
    delta_phi_nogap->SetLineWidth(3);
    delta_phi_nogap->Draw("hist");

    TLegend *leg3 = new TLegend(0.5,0.9,0.97,0.97);
    leg3->AddEntry(delta_phi_nogap,"Average Pile-up uncertainty for 2010 Data","l");
    leg3->SetFillColor(0);
    leg3->Draw();
    
    ave = tot/delta_phi_nogap->GetNbinsX();
    ave_old = tot_old/delta_phi_nogap->GetNbinsX();
    cout<<"No gap: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[4] = min*100;
    max_row[4] = max*100;
    ave_row[4] = ave*100;
    old_row[4] = ave_old*100;

    c3->Print("output/uncertainities/png/pile_up_delta_phi_nogap.png");
    c3->Print("output/uncertainities/c/pile_up_delta_phi_nogap.C");
    c3->Print("output/uncertainities/eps/pile_up_delta_phi_nogap.eps");
    c3->Close();
    
    TCanvas *c3o = new TCanvas("c3b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    delta_phi_nogap->SetMinimum(0.0);
    delta_phi_nogap->SetMaximum(0.25);
    delta_phi_nogap->Draw("hist");
    delta_phi_nogap2->SetLineColor(2);
    delta_phi_nogap2->SetLineStyle(3);
    delta_phi_nogap2->SetLineWidth(3);
    delta_phi_nogap2->Draw("hist same");
    delta_phi_nogap_old->SetLineColor(2);
    delta_phi_nogap_old->SetLineStyle(2);
    delta_phi_nogap_old->SetLineWidth(3);
    delta_phi_nogap_old->Draw("hist same");

    TLegend *leg3o = new TLegend(0.5,0.84,0.97,0.97);
    leg3o->AddEntry(delta_phi_nogap,"Average Pile-up uncertainty for 2010 Data","l");
    leg3o->AddEntry(delta_phi_nogap2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg3o->AddEntry(delta_phi_nogap_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg3o->SetFillColor(0);
    leg3o->Draw();

    c3o->Print("output/uncertainities/png/pile_up_delta_phi_nogap_a.png");
    c3o->Print("output/uncertainities/c/pile_up_delta_phi_nogap_a.C");
    c3o->Print("output/uncertainities/eps/pile_up_delta_phi_nogap_a.eps");
    c3o->Close();

    TCanvas *c3b = new TCanvas("c3b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    delta_phi_nogap_a->SetMaximum(1.5);
    delta_phi_nogap_a->SetMinimum(0.5);
    delta_phi_nogap_a->SetLineColor(2);
    delta_phi_nogap_a->SetLineStyle(1);
    delta_phi_nogap_a->SetLineWidth(3);
    delta_phi_nogap_a->Draw("e1");
    delta_phi_nogap_b->SetLineColor(3);
    delta_phi_nogap_b->SetLineStyle(1);
    delta_phi_nogap_b->SetLineWidth(3);
    delta_phi_nogap_b->Draw("e1 same");
    delta_phi_nogap_c->SetLineColor(4);
    delta_phi_nogap_c->SetLineStyle(1);
    delta_phi_nogap_c->SetLineWidth(3);
    delta_phi_nogap_c->Draw("e1 same");

    TLegend *leg3b = new TLegend(0.5,0.84,0.97,0.97);
    leg3b->AddEntry(delta_phi_nogap_a,labela,"l");
    leg3b->AddEntry(delta_phi_nogap_b,labelb,"l");
    leg3b->AddEntry(delta_phi_nogap_c,labelc,"l");
    leg3b->SetFillColor(0);
    leg3b->Draw();
    
    c3b->Print("output/uncertainities/png/pile_up_delta_phi_nogap_b.png");
    c3b->Print("output/uncertainities/c/pile_up_delta_phi_nogap_b.C");
    c3b->Print("output/uncertainities/eps/pile_up_delta_phi_nogap_b.eps");
    c3b->Close();
    
    TCanvas *c4 = new TCanvas("c4","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x3a = (TH1D*) data_file1a->Get("ak5PF_leading_pt_inside_gap");
     TH1D *x4a = (TH1D*) data_file2a->Get("ak5PF_leading_pt_inside_gap");
     TH1D *x3b = (TH1D*) data_file1b->Get("ak5PF_leading_pt_inside_gap");
     TH1D *x4b = (TH1D*) data_file2b->Get("ak5PF_leading_pt_inside_gap");
     TH1D *x3c = (TH1D*) data_file1c->Get("ak5PF_leading_pt_inside_gap");
     TH1D *x4c = (TH1D*) data_file2c->Get("ak5PF_leading_pt_inside_gap");
     
     enta = x3a->GetEntries()-x3a->GetNbinsX();
     entb = x3b->GetEntries()-x3b->GetNbinsX();
     entc = x3c->GetEntries()-x3c->GetNbinsX();

     int1a = x3a->Integral();
     int2a = x4a->Integral();
     factora = int2a/int1a;

     int1b = x3b->Integral();
     int2b = x4b->Integral();
     factorb = int2b/int1b;

     int1c = x3c->Integral();
     int2c = x4c->Integral();
     factorc = int2c/int1c;
     
     cout<<"factor"<<factora<<endl;
    
    TH1D *leading_pt_inside_gap, *leading_pt_inside_gap2, *leading_pt_inside_gap_old, *leading_pt_inside_gap_a, *leading_pt_inside_gap_b, *leading_pt_inside_gap_c;
    leading_pt_inside_gap =  new TH1D("Leading_pt_inside_gap","Leading_pt_inside_gap;p_{T} [GeV];Pile Up Uncertainty", in_nbins, in_bins);
    leading_pt_inside_gap2 =  new TH1D("Leading_pt_inside_gap2","Leading_pt_inside_gap2;p_{T} [GeV];Pile Up Uncertainty", in_nbins, in_bins);
    leading_pt_inside_gap_old =  new TH1D("Leading_pt_inside_gap_old","Leading_pt_inside_gap_old;p_{T} [GeV];Pile Up Uncertainty", in_nbins, in_bins);
    leading_pt_inside_gap_a =  new TH1D("Leading_pt_inside_gap_a","Leading_pt_inside_gap_a;p_{T} [GeV];Ratio to all vertexes selection", in_nbins, in_bins);
    leading_pt_inside_gap_b =  new TH1D("Leading_pt_inside_gap_b","Leading_pt_inside_gap_b;p_{T} [GeV];Ratio to all vertexes selection", in_nbins, in_bins);
    leading_pt_inside_gap_c =  new TH1D("Leading_pt_inside_gap_c","Leading_pt_inside_gap_c;p_{T} [GeV];Ratio to all vertexes selection", in_nbins, in_bins);
    leading_pt_inside_gap_a->Divide(x3a,x4a,1.,1./factora,"");
    leading_pt_inside_gap_b->Divide(x3b,x4b,1.,1./factorb,"");
    leading_pt_inside_gap_c->Divide(x3c,x4c,1.,1./factorc,"");
    
    min = 0;
    max = 0;
    tot = 0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=leading_pt_inside_gap_a->GetNbinsX();i++)
    {
    Float_t conta = leading_pt_inside_gap_a->GetBinContent(i) - 1;
    Float_t contb = leading_pt_inside_gap_b->GetBinContent(i) - 1;
    Float_t contc = leading_pt_inside_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    leading_pt_inside_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    leading_pt_inside_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(leading_pt_inside_gap2, leading_pt_inside_gap, step);

    leading_pt_inside_gap->SetMinimum(0.0);
    leading_pt_inside_gap->SetMaximum(0.25);
    leading_pt_inside_gap->SetLineColor(2);
    leading_pt_inside_gap->SetLineStyle(1);
    leading_pt_inside_gap->SetLineWidth(3);
    leading_pt_inside_gap->Draw("hist");
    
    TLegend *leg4 = new TLegend(0.13,0.9,0.55,0.97);
    leg4->AddEntry(leading_pt_inside_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg4->SetFillColor(0);
    leg4->Draw();

    ave = tot/(leading_pt_inside_gap->GetNbinsX());
    ave_old = tot_old/leading_pt_inside_gap->GetNbinsX();
    cout<<"leading_pt_inside_gap: min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[6] = min*100;
    max_row[6] = max*100;
    ave_row[6] = ave*100;
    old_row[6] = ave_old*100;   

    c4->Print("output/uncertainities/png/pile_up_leading_pt_inside_gap.png");
    c4->Print("output/uncertainities/c/pile_up_leading_pt_inside_gap.C");
    c4->Print("output/uncertainities/eps/pile_up_leading_pt_inside_gap.eps");
    c4->Close();
    
    TCanvas *c4o = new TCanvas("c4b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    leading_pt_inside_gap->SetMinimum(0.0);
    leading_pt_inside_gap->SetMaximum(0.25);
    leading_pt_inside_gap->Draw("hist");
    leading_pt_inside_gap2->SetLineColor(2);
    leading_pt_inside_gap2->SetLineStyle(3);
    leading_pt_inside_gap2->SetLineWidth(3);
    leading_pt_inside_gap2->Draw("hist same");
    leading_pt_inside_gap_old->SetLineColor(2);
    leading_pt_inside_gap_old->SetLineStyle(2);
    leading_pt_inside_gap_old->SetLineWidth(3);
    leading_pt_inside_gap_old->Draw("hist same");

    TLegend *leg4o = new TLegend(0.13,0.84,0.55,0.97);
    leg4o->AddEntry(leading_pt_inside_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg4o->AddEntry(leading_pt_inside_gap2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg4o->AddEntry(leading_pt_inside_gap_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg4o->SetFillColor(0);
    leg4o->Draw();

    c4o->Print("output/uncertainities/png/pile_up_leading_pt_inside_gap_a.png");
    c4o->Print("output/uncertainities/c/pile_up_leading_pt_inside_gap_a.C");
    c4o->Print("output/uncertainities/eps/pile_up_leading_pt_inside_gap_a.eps");
    c4o->Close();
    
    TCanvas *c4b = new TCanvas("c4b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    leading_pt_inside_gap_a->SetMaximum(1.5);
    leading_pt_inside_gap_a->SetMinimum(0.5);
    leading_pt_inside_gap_a->SetLineColor(2);
    leading_pt_inside_gap_a->SetLineStyle(1);
    leading_pt_inside_gap_a->SetLineWidth(3);
    leading_pt_inside_gap_a->Draw("e1");
    leading_pt_inside_gap_b->SetLineColor(3);
    leading_pt_inside_gap_b->SetLineStyle(1);
    leading_pt_inside_gap_b->SetLineWidth(3);
    leading_pt_inside_gap_b->Draw("e1 same");
    leading_pt_inside_gap_c->SetLineColor(4);
    leading_pt_inside_gap_c->SetLineStyle(1);
    leading_pt_inside_gap_c->SetLineWidth(3);
    leading_pt_inside_gap_c->Draw("e1 same");

    TLegend *leg4b = new TLegend(0.5,0.84,0.97,0.97);
    leg4b->AddEntry(leading_pt_inside_gap_a,labela,"l");
    leg4b->AddEntry(leading_pt_inside_gap_b,labelb,"l");
    leg4b->AddEntry(leading_pt_inside_gap_c,labelc,"l");
    leg4b->SetFillColor(0);
    leg4b->Draw();
    
    c4b->Print("output/uncertainities/png/pile_up_leading_pt_inside_gap_b.png");
    c4b->Print("output/uncertainities/c/pile_up_leading_pt_inside_gap_b.C");
    c4b->Print("output/uncertainities/eps/pile_up_leading_pt_inside_gap_b.eps");
    c4b->Close();


    TCanvas *c5 = new TCanvas("c5","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x3a2 = (TH1D*) data_file1a->Get("ak5PF_leading_pt_outside_gap");
     TH1D *x4a2 = (TH1D*) data_file2a->Get("ak5PF_leading_pt_outside_gap");
     TH1D *x3b2 = (TH1D*) data_file1b->Get("ak5PF_leading_pt_outside_gap");
     TH1D *x4b2 = (TH1D*) data_file2b->Get("ak5PF_leading_pt_outside_gap");
     TH1D *x3c2 = (TH1D*) data_file1c->Get("ak5PF_leading_pt_outside_gap");
     TH1D *x4c2 = (TH1D*) data_file2c->Get("ak5PF_leading_pt_outside_gap");
     
     enta = x3a2->GetEntries()-x3a2->GetNbinsX();
     entb = x3b2->GetEntries()-x3b2->GetNbinsX();
     entc = x3c2->GetEntries()-x3c2->GetNbinsX();

     int1a = x3a2->Integral();
     int2a = x4a2->Integral();
     factora = int2a/int1a;

     int1b = x3b2->Integral();
     int2b = x4b2->Integral();
     factorb = int2b/int1b;

     int1c = x3c2->Integral();
     int2c = x4c2->Integral();
     factorc = int2c/int1c;
     
     cout<<"factor"<<factora<<endl;
    
    TH1D *leading_pt_outside_gap, *leading_pt_outside_gap2, *leading_pt_outside_gap_old, *leading_pt_outside_gap_a, *leading_pt_outside_gap_b, *leading_pt_outside_gap_c;
    leading_pt_outside_gap =  new TH1D("Leading_pt_outside_gap","Leading_pt_outside_gap;p_{T} [GeV];Pile Up Uncertainty", in_nbins, in_bins);
    leading_pt_outside_gap2 =  new TH1D("Leading_pt_outside_gap2","Leading_pt_outside_gap2;p_{T} [GeV];Pile Up Uncertainty", in_nbins, in_bins);
    leading_pt_outside_gap_old =  new TH1D("Leading_pt_outside_gap_old","Leading_pt_outside_gap_old;p_{T} [GeV];Pile Up Uncertainty", in_nbins, in_bins);
    leading_pt_outside_gap_a =  new TH1D("Leading_pt_outside_gap_a","Leading_pt_outside_gap_a;p_{T} [GeV];Ratio to all vertexes selection", in_nbins, in_bins);
    leading_pt_outside_gap_b =  new TH1D("Leading_pt_outside_gap_b","Leading_pt_outside_gap_b;p_{T} [GeV];Ratio to all vertexes selection", out_nbins, out_bins);
    leading_pt_outside_gap_c =  new TH1D("Leading_pt_outside_gap_c","Leading_pt_outside_gap_c;p_{T} [GeV];Ratio to all vertexes selection", out_nbins, out_bins);
    leading_pt_outside_gap_a->Divide(x3a2,x4a2,1.,1./factora,"");
    leading_pt_outside_gap_b->Divide(x3b2,x4b2,1.,1./factorb,"");
    leading_pt_outside_gap_c->Divide(x3c2,x4c2,1.,1./factorc,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=leading_pt_outside_gap_a->GetNbinsX();i++)
    {
    Float_t conta = leading_pt_outside_gap_a->GetBinContent(i) - 1;
    Float_t contb = leading_pt_outside_gap_b->GetBinContent(i) - 1;
    Float_t contc = leading_pt_outside_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    leading_pt_outside_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    leading_pt_outside_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }

    alisar(leading_pt_outside_gap2, leading_pt_outside_gap, step);
    
    leading_pt_outside_gap->SetMinimum(0.0);
    leading_pt_outside_gap->SetMaximum(0.25);
    leading_pt_outside_gap->SetLineColor(2);
    leading_pt_outside_gap->SetLineStyle(1);
    leading_pt_outside_gap->SetLineWidth(3);
    leading_pt_outside_gap->Draw("hist");
    
    TLegend *leg5 = new TLegend(0.13,0.9,0.55,0.97);
    leg5->AddEntry(leading_pt_outside_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg5->SetFillColor(0);
    leg5->Draw();

    ave = tot/(leading_pt_outside_gap->GetNbinsX());
    ave_old = tot_old/leading_pt_outside_gap->GetNbinsX();
    cout<<"leading_pt_outside_gap: min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[7] = min*100;
    max_row[7] = max*100;
    ave_row[7] = ave*100;
    old_row[7] = ave_old*100;  

    c5->Print("output/uncertainities/png/pile_up_leading_pt_outside_gap.png");
    c5->Print("output/uncertainities/c/pile_up_leading_pt_outside_gap.C");
    c5->Print("output/uncertainities/eps/pile_up_leading_pt_outside_gap.eps");
    c5->Close();
    
    TCanvas *c5o = new TCanvas("c5b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    leading_pt_outside_gap->SetMinimum(0.0);
    leading_pt_outside_gap->SetMaximum(0.25);
    leading_pt_outside_gap->Draw("hist");
    leading_pt_outside_gap2->SetLineColor(2);
    leading_pt_outside_gap2->SetLineStyle(3);
    leading_pt_outside_gap2->SetLineWidth(3);
    leading_pt_outside_gap2->Draw("hist same");
    leading_pt_outside_gap_old->SetLineColor(2);
    leading_pt_outside_gap_old->SetLineStyle(2);
    leading_pt_outside_gap_old->SetLineWidth(3);
    leading_pt_outside_gap_old->Draw("hist same");

    TLegend *leg5o = new TLegend(0.13,0.84,0.55,0.97);
    leg5o->AddEntry(leading_pt_outside_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg5o->AddEntry(leading_pt_outside_gap2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg5o->AddEntry(leading_pt_outside_gap_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg5o->SetFillColor(0);
    leg5o->Draw();

    c5o->Print("output/uncertainities/png/pile_up_leading_pt_outside_gap_a.png");
    c5o->Print("output/uncertainities/c/pile_up_leading_pt_outside_gap_a.C");
    c5o->Print("output/uncertainities/eps/pile_up_leading_pt_outside_gap_a.eps");
    c5o->Close();
    
    TCanvas *c5b = new TCanvas("c5b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    leading_pt_outside_gap_a->SetMaximum(1.4);
    leading_pt_outside_gap_a->SetMinimum(0.4);
    leading_pt_outside_gap_a->SetLineColor(2);
    leading_pt_outside_gap_a->SetLineStyle(1);
    leading_pt_outside_gap_a->SetLineWidth(3);
    leading_pt_outside_gap_a->Draw("e1");
    leading_pt_outside_gap_b->SetLineColor(3);
    leading_pt_outside_gap_b->SetLineStyle(1);
    leading_pt_outside_gap_b->SetLineWidth(3);
    leading_pt_outside_gap_b->Draw("e1 same");
    leading_pt_outside_gap_c->SetLineColor(4);
    leading_pt_outside_gap_c->SetLineStyle(1);
    leading_pt_outside_gap_c->SetLineWidth(3);
    leading_pt_outside_gap_c->Draw("e1 same");

    TLegend *leg5b = new TLegend(0.5,0.84,0.97,0.97);
    leg5b->AddEntry(leading_pt_outside_gap_a,labela,"l");
    leg5b->AddEntry(leading_pt_outside_gap_b,labelb,"l");
    leg5b->AddEntry(leading_pt_outside_gap_c,labelc,"l");
    leg5b->SetFillColor(0);
    leg5b->Draw();
    
    c5b->Print("output/uncertainities/png/pile_up_leading_pt_outside_gap_b.png");
    c5b->Print("output/uncertainities/c/pile_up_leading_pt_outside_gap_b.C");
    c5b->Print("output/uncertainities/eps/pile_up_leading_pt_outside_gap_b.eps");
    c5b->Close();
    
    TCanvas *c6 = new TCanvas("c6","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c6->Divide(2,2);
    
    c6->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1a1 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta1");
     TH1D *x2a1 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta1");
     TH1D *x1b1 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta1");
     TH1D *x2b1 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta1");
     TH1D *x1c1 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta1");
     TH1D *x2c1 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta1");
     
     enta = x1a1->GetEntries()-x1a1->GetNbinsX();
     entb = x1b1->GetEntries()-x1b1->GetNbinsX();
     entc = x1c1->GetEntries()-x1c1->GetNbinsX();

     int1a = x1a1->Integral();
     int2a = x2a1->Integral();
     factora = int2a/int1a;

     int1b = x1b1->Integral();
     int2b = x2b1->Integral();
     factorb = int2b/int1b;

     int1c = x1c1->Integral();
     int2c = x2c1->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta1, *delta_phi_deta12, *delta_phi_deta1_old, *delta_phi_deta1_a, *delta_phi_deta1_b, *delta_phi_deta1_c;
    delta_phi_deta1 =  new TH1D("Delta_phi_deta1","delta_phi_deta1;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta12 =  new TH1D("Delta_phi_deta12","delta_phi_deta12;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_old =  new TH1D("Delta_phi_deta1_old","delta_phi_deta1_old;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_a =  new TH1D("Delta_phi_deta1_a","delta_phi_deta1_a;|#Delta#phi| [rad] for 0.3 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_b =  new TH1D("Delta_phi_deta1_b","delta_phi_deta1_b;|#Delta#phi| [rad] for 0.3 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_c =  new TH1D("Delta_phi_deta1_c","delta_phi_deta1_c;|#Delta#phi| [rad] for 0.3 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_a->Divide(x1a1,x2a1,1.,1./factora,"");
    delta_phi_deta1_b->Divide(x1b1,x2b1,1.,1./factorb,"");
    delta_phi_deta1_c->Divide(x1c1,x2c1,1.,1./factorc,"");

    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_phi_deta1_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta1_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta1_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta1_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta12->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta1_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }

    alisar(delta_phi_deta12, delta_phi_deta1, step);
    
    delta_phi_deta1->SetMinimum(0.0);
    delta_phi_deta1->SetMaximum(0.45);
    delta_phi_deta1->SetLineColor(2);
    delta_phi_deta1->SetLineStyle(1);
    delta_phi_deta1->SetLineWidth(3);
    delta_phi_deta1->Draw("hist");
    
    TLegend *leg6 = new TLegend(0.40,0.88,0.98,0.98);
    leg6->AddEntry(delta_phi_deta1,"Average Pile-up uncertainty for 2010 Data","l");
    leg6->SetFillColor(0);
    leg6->Draw();

    c6->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1a2 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta2");
     TH1D *x2a2 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta2");
     TH1D *x1b2 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta2");
     TH1D *x2b2 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta2");
     TH1D *x1c2 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta2");
     TH1D *x2c2 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta2");
     
     enta = x1a2->GetEntries()-x1a2->GetNbinsX();
     entb = x1b2->GetEntries()-x1b2->GetNbinsX();
     entc = x1c2->GetEntries()-x1c2->GetNbinsX();

     int1a = x1a2->Integral();
     int2a = x2a2->Integral();
     factora = int2a/int1a;

     int1b = x1b2->Integral();
     int2b = x2b2->Integral();
     factorb = int2b/int1b;

     int1c = x1c2->Integral();
     int2c = x2c2->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta2, *delta_phi_deta22, *delta_phi_deta2_old, *delta_phi_deta2_a, *delta_phi_deta2_b, *delta_phi_deta2_c;
    delta_phi_deta2 =  new TH1D("Delta_phi_deta2","delta_phi_deta2;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta22 =  new TH1D("Delta_phi_deta22","delta_phi_deta22;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_old =  new TH1D("Delta_phi_deta2_old","delta_phi_deta2_old;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_a =  new TH1D("Delta_phi_deta2_a","delta_phi_deta2_a;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_b =  new TH1D("Delta_phi_deta2_b","delta_phi_deta2_b;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_c =  new TH1D("Delta_phi_deta2_c","delta_phi_deta2_c;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_a->Divide(x1a2,x2a2,1.,1./factora,"");
    delta_phi_deta2_b->Divide(x1b2,x2b2,1.,1./factorb,"");
    delta_phi_deta2_c->Divide(x1c2,x2c2,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta2_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta2_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta2_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta2_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta22->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta2_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta22, delta_phi_deta2, step);

    delta_phi_deta2->SetMinimum(0.0);
    delta_phi_deta2->SetMaximum(0.45);
    delta_phi_deta2->SetLineColor(2);
    delta_phi_deta2->SetLineStyle(1);
    delta_phi_deta2->SetLineWidth(3);
    delta_phi_deta2->Draw("hist");
    
    c6->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1a3 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta3");
     TH1D *x2a3 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta3");
     TH1D *x1b3 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta3");
     TH1D *x2b3 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta3");
     TH1D *x1c3 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta3");
     TH1D *x2c3 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta3");
     
     enta = x1a3->GetEntries()-x1a3->GetNbinsX();
     entb = x1b3->GetEntries()-x1b3->GetNbinsX();
     entc = x1c3->GetEntries()-x1c3->GetNbinsX();

     int1a = x1a3->Integral();
     int2a = x2a3->Integral();
     factora = int2a/int1a;

     int1b = x1b3->Integral();
     int2b = x2b3->Integral();
     factorb = int2b/int1b;

     int1c = x1c3->Integral();
     int2c = x2c3->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta3, *delta_phi_deta32, *delta_phi_deta3_old, *delta_phi_deta3_a, *delta_phi_deta3_b, *delta_phi_deta3_c;
    delta_phi_deta3 =  new TH1D("Delta_phi_deta3","delta_phi_deta3;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta32 =  new TH1D("Delta_phi_deta32","delta_phi_deta32;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_old =  new TH1D("Delta_phi_deta3_old","delta_phi_deta3_old;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_a =  new TH1D("Delta_phi_deta3_a","delta_phi_deta3_a;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_b =  new TH1D("Delta_phi_deta3_b","delta_phi_deta3_b;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_c =  new TH1D("Delta_phi_deta3_c","delta_phi_deta3_c;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_a->Divide(x1a3,x2a3,1.,1./factora,"");
    delta_phi_deta3_b->Divide(x1b3,x2b3,1.,1./factorb,"");
    delta_phi_deta3_c->Divide(x1c3,x2c3,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta3_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta3_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta3_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta3_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta32->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta3_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta32, delta_phi_deta3, step);

    delta_phi_deta3->SetMinimum(0.0);
    delta_phi_deta3->SetMaximum(0.45);
    delta_phi_deta3->SetLineColor(2);
    delta_phi_deta3->SetLineStyle(1);
    delta_phi_deta3->SetLineWidth(3);
    delta_phi_deta3->Draw("hist");
    
    c6->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1a4 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta4");
     TH1D *x2a4 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta4");
     TH1D *x1b4 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta4");
     TH1D *x2b4 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta4");
     TH1D *x1c4 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta4");
     TH1D *x2c4 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta4");
     
     enta = x1a4->GetEntries()-x1a4->GetNbinsX();
     entb = x1b4->GetEntries()-x1b4->GetNbinsX();
     entc = x1c4->GetEntries()-x1c4->GetNbinsX();

     int1a = x1a4->Integral();
     int2a = x2a4->Integral();
     factora = int2a/int1a;

     int1b = x1b4->Integral();
     int2b = x2b4->Integral();
     factorb = int2b/int1b;

     int1c = x1c4->Integral();
     int2c = x2c4->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta4, *delta_phi_deta42, *delta_phi_deta4_old, *delta_phi_deta4_a, *delta_phi_deta4_b, *delta_phi_deta4_c;
    delta_phi_deta4 =  new TH1D("Delta_phi_deta4","delta_phi_deta4;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta42 =  new TH1D("Delta_phi_deta42","delta_phi_deta42;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_old =  new TH1D("Delta_phi_deta4_old","delta_phi_deta4_old;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_a =  new TH1D("Delta_phi_deta4_a","delta_phi_deta4_a;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_b =  new TH1D("Delta_phi_deta4_b","delta_phi_deta4_b;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_c =  new TH1D("Delta_phi_deta4_c","delta_phi_deta4_c;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_a->Divide(x1a4,x2a4,1.,1./factora,"");
    delta_phi_deta4_b->Divide(x1b4,x2b4,1.,1./factorb,"");
    delta_phi_deta4_c->Divide(x1c4,x2c4,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta4_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta4_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta4_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta4_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta42->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta4_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta42, delta_phi_deta4, step);

    delta_phi_deta4->SetMinimum(0.0);
    delta_phi_deta4->SetMaximum(0.45);
    delta_phi_deta4->SetLineColor(2);
    delta_phi_deta4->SetLineStyle(1);
    delta_phi_deta4->SetLineWidth(3);
    delta_phi_deta4->Draw("hist"); 
    
    
  /*      c5->cd(5);
     TH1D *x1a5 = (TH1D*) data_file1->Get("ak5PF_delta_phi_deta5");
     TH1D *x2a5 = (TH1D*) data_file2->Get("ak5PF_delta_phi_deta5");
     
     double int1 = x1a5->Integral();
     double int2 = x2a5->Integral();
     double factor = int2/int1;
     
     //cout<<"factor"<<factor<<endl;
    
    delta_phi_deta5 =  new TH1D("Delta_phi_deta5","delta_phi_deta5;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta5_b =  new TH1D("Delta_phi_deta5_b","delta_phi_deta5_b;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta5_b->Divide(x1a5,x2a5,1.,1./factor,"");
    
    for(Int_t i=1;i<=delta_phi_deta5_b->GetNbinsX();i++)
    {
    Float_t cont = delta_phi_deta5_b->GetBinContent(i) - 1;
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min) { min = cont;}
    delta_phi_deta5->SetBinContent(i,cont);
    }
    
    delta_phi_deta5->SetLineColor(2);
    delta_phi_deta5->SetLineStyle(1);
    delta_phi_deta5->SetLineWidth(2);
    delta_phi_deta5->Draw("hist");
    
            c5->cd(6);
     TH1D *x1a6 = (TH1D*) data_file1->Get("ak5PF_delta_phi_deta6");
     TH1D *x2a6 = (TH1D*) data_file2->Get("ak5PF_delta_phi_deta6");
     
     double int1 = x1a6->Integral();
     double int2 = x2a6->Integral();
     double factor = int2/int1;
     
     //cout<<"factor"<<factor<<endl;
    
    delta_phi_deta6 =  new TH1D("Delta_phi_deta6","delta_phi_deta6;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta6_b =  new TH1D("Delta_phi_deta6_b","delta_phi_deta6_b;#Delta#phi [rad];Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta6_b->Divide(x1a6,x2a6,1.,1./factor,"");
    
    for(Int_t i=1;i<=delta_phi_deta6_b->GetNbinsX();i++)
    {
    Float_t cont = delta_phi_deta6_b->GetBinContent(i) - 1;
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min) { min = cont;}
    delta_phi_deta6->SetBinContent(i,cont);
    }
    
    delta_phi_deta6->SetLineColor(2);
    delta_phi_deta6->SetLineStyle(1);
    delta_phi_deta6->SetLineWidth(2);
    delta_phi_deta6->Draw("hist"); */
    
    ave = tot/(delta_phi_deta1->GetNbinsX()*4);
    ave_old = tot_old/(delta_phi_deta1->GetNbinsX()*4);
    cout<<"Main selection deta: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[1] = min*100;
    max_row[1] = max*100;
    ave_row[1] = ave*100;
    old_row[1] = ave_old*100;      

    c6->Print("output/uncertainities/png/pile_up_delta_phi_deta.png");
    c6->Print("output/uncertainities/c/pile_up_delta_phi_deta.C");
    c6->Print("output/uncertainities/eps/pile_up_delta_phi_deta.eps");
    c6->Close();
    
    TCanvas *c6o = new TCanvas("c6o","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c6o->Divide(2,2);

    c6o->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta1->SetMinimum(0.0);
    delta_phi_deta1->Draw("hist");
    delta_phi_deta12->SetLineColor(2);
    delta_phi_deta12->SetLineStyle(3);
    delta_phi_deta12->SetLineWidth(3);
    delta_phi_deta12->Draw("hist same");
    delta_phi_deta1_old->SetLineColor(2);
    delta_phi_deta1_old->SetLineStyle(2);
    delta_phi_deta1_old->SetLineWidth(3);
    delta_phi_deta1_old->Draw("hist same");

    TLegend *leg6o = new TLegend(0.30,0.80,0.98,0.98);
    leg6o->AddEntry(delta_phi_deta1,"Average Pile-up uncertainty for 2010 Data","l");
    leg6o->AddEntry(delta_phi_deta12,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg6o->AddEntry(delta_phi_deta1_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg6o->SetFillColor(0);
    leg6o->Draw();

    c6o->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta2->SetMinimum(0.0);
    delta_phi_deta2->Draw("hist");
    delta_phi_deta22->SetLineColor(2);
    delta_phi_deta22->SetLineStyle(3);
    delta_phi_deta22->SetLineWidth(3);
    delta_phi_deta22->Draw("hist same");
    delta_phi_deta2_old->SetLineColor(2);
    delta_phi_deta2_old->SetLineStyle(2);
    delta_phi_deta2_old->SetLineWidth(3);
    delta_phi_deta2_old->Draw("hist same");

    c6o->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta3->SetMinimum(0.0);
    delta_phi_deta3->Draw("hist");
    delta_phi_deta32->SetLineColor(2);
    delta_phi_deta32->SetLineStyle(3);
    delta_phi_deta32->SetLineWidth(3);
    delta_phi_deta32->Draw("hist same");
    delta_phi_deta3_old->SetLineColor(2);
    delta_phi_deta3_old->SetLineStyle(2);
    delta_phi_deta3_old->SetLineWidth(3);
    delta_phi_deta3_old->Draw("hist same");

    c6o->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta4->SetMinimum(0.0);
    delta_phi_deta4->Draw("hist");
    delta_phi_deta42->SetLineColor(2);
    delta_phi_deta42->SetLineStyle(3);
    delta_phi_deta42->SetLineWidth(3);
    delta_phi_deta42->Draw("hist same");
    delta_phi_deta4_old->SetLineColor(2);
    delta_phi_deta4_old->SetLineStyle(2);
    delta_phi_deta4_old->SetLineWidth(3);
    delta_phi_deta4_old->Draw("hist same");

    c6o->Print("output/uncertainities/png/pile_up_delta_phi_deta_a.png");
    c6o->Print("output/uncertainities/c/pile_up_delta_phi_deta_a.C");
    c6o->Print("output/uncertainities/eps/pile_up_delta_phi_deta_a.eps");
    c6o->Close();

    TCanvas *c6b = new TCanvas("c6b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c6b->Divide(2,2);

    c6b->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta1_a->SetMaximum(1.8);
    delta_phi_deta1_a->SetMinimum(0.8);
    delta_phi_deta1_a->SetLineColor(2);
    delta_phi_deta1_a->SetLineStyle(1);
    delta_phi_deta1_a->SetLineWidth(3);
    delta_phi_deta1_a->Draw("e1");
    delta_phi_deta1_b->SetLineColor(3);
    delta_phi_deta1_b->SetLineStyle(1);
    delta_phi_deta1_b->SetLineWidth(3);
    delta_phi_deta1_b->Draw("e1 same");
    delta_phi_deta1_c->SetLineColor(4);
    delta_phi_deta1_c->SetLineStyle(1);
    delta_phi_deta1_c->SetLineWidth(3);
    delta_phi_deta1_c->Draw("e1 same");

    TLegend *leg6b = new TLegend(0.5,0.75,0.98,0.98);
    leg6b->AddEntry(delta_phi_deta1_a,labela,"l");
    leg6b->AddEntry(delta_phi_deta1_b,labelb,"l");
    leg6b->AddEntry(delta_phi_deta1_c,labelc,"l");
    leg6b->SetFillColor(0);
    leg6b->Draw();

    c6b->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta2_a->SetMaximum(1.8);
    delta_phi_deta2_a->SetMinimum(0.8);
    delta_phi_deta2_a->SetLineColor(2);
    delta_phi_deta2_a->SetLineStyle(1);
    delta_phi_deta2_a->SetLineWidth(3);
    delta_phi_deta2_a->Draw("e1");
    delta_phi_deta2_b->SetLineColor(3);
    delta_phi_deta2_b->SetLineStyle(1);
    delta_phi_deta2_b->SetLineWidth(3);
    delta_phi_deta2_b->Draw("e1 same");
    delta_phi_deta2_c->SetLineColor(4);
    delta_phi_deta2_c->SetLineStyle(1);
    delta_phi_deta2_c->SetLineWidth(3);
    delta_phi_deta2_c->Draw("e1 same");

    c6b->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta3_a->SetMaximum(1.8);
    delta_phi_deta3_a->SetMinimum(0.8);
    delta_phi_deta3_a->SetLineColor(2);
    delta_phi_deta3_a->SetLineStyle(1);
    delta_phi_deta3_a->SetLineWidth(3);
    delta_phi_deta3_a->Draw("e1");
    delta_phi_deta3_b->SetLineColor(3);
    delta_phi_deta3_b->SetLineStyle(1);
    delta_phi_deta3_b->SetLineWidth(3);
    delta_phi_deta3_b->Draw("e1 same");
    delta_phi_deta3_c->SetLineColor(4);
    delta_phi_deta3_c->SetLineStyle(1);
    delta_phi_deta3_c->SetLineWidth(3);
    delta_phi_deta3_c->Draw("e1 same");

    c6b->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta4_a->SetMaximum(1.8);
    delta_phi_deta4_a->SetMinimum(0.8);
    delta_phi_deta4_a->SetLineColor(2);
    delta_phi_deta4_a->SetLineStyle(1);
    delta_phi_deta4_a->SetLineWidth(3);
    delta_phi_deta4_a->Draw("e1");
    delta_phi_deta4_b->SetLineColor(3);
    delta_phi_deta4_b->SetLineStyle(1);
    delta_phi_deta4_b->SetLineWidth(3);
    delta_phi_deta4_b->Draw("e1 same");
    delta_phi_deta4_c->SetLineColor(4);
    delta_phi_deta4_c->SetLineStyle(1);
    delta_phi_deta4_c->SetLineWidth(3);
    delta_phi_deta4_c->Draw("e1 same");

    c6b->Print("output/uncertainities/png/pile_up_delta_phi_deta_b.png");
    c6b->Print("output/uncertainities/c/pile_up_delta_phi_deta_b.C");
    c6b->Print("output/uncertainities/eps/pile_up_delta_phi_deta_b.eps");
    c6b->Close();

    
    TCanvas *c7 = new TCanvas("c7","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c7->Divide(2,2);
    
     c7->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1ag1 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta1_gap");
     TH1D *x2ag1 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta1_gap");
     TH1D *x1bg1 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta1_gap");
     TH1D *x2bg1 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta1_gap");
     TH1D *x1cg1 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta1_gap");
     TH1D *x2cg1 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta1_gap");
     
     enta = x1ag1->GetEntries()-x1ag1->GetNbinsX();
     entb = x1bg1->GetEntries()-x1bg1->GetNbinsX();
     entc = x1cg1->GetEntries()-x1cg1->GetNbinsX();

     int1a = x1ag1->Integral();
     int2a = x2ag1->Integral();
     factora = int2a/int1a;

     int1b = x1bg1->Integral();
     int2b = x2bg1->Integral();
     factorb = int2b/int1b;

     int1c = x1cg1->Integral();
     int2c = x2cg1->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta1_gap, *delta_phi_deta1_gap2, *delta_phi_deta1_gap_old, *delta_phi_deta1_gap_a, *delta_phi_deta1_gap_b, *delta_phi_deta1_gap_c;
    delta_phi_deta1_gap =  new TH1D("Delta_phi_deta1_gap","delta_phi_deta1_gap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap2 =  new TH1D("Delta_phi_deta1_gap2","delta_phi_deta1_gap2;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap_old =  new TH1D("Delta_phi_deta1_gap_old","delta_phi_deta1_gap_old;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap_a =  new TH1D("Delta_phi_deta1_gap_a","delta_phi_deta1_gap_a;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap_b =  new TH1D("Delta_phi_deta1_gap_b","delta_phi_deta1_gap_b;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap_c =  new TH1D("Delta_phi_deta1_gap_c","delta_phi_deta1_gap_c;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_gap_a->Divide(x1ag1,x2ag1,1.,1./factora,"");
    delta_phi_deta1_gap_b->Divide(x1bg1,x2bg1,1.,1./factorb,"");
    delta_phi_deta1_gap_c->Divide(x1cg1,x2cg1,1.,1./factorc,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_phi_deta1_gap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta1_gap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta1_gap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta1_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta1_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta1_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }

    alisar(delta_phi_deta1_gap2, delta_phi_deta1_gap, step);
    
    delta_phi_deta1_gap->SetMinimum(0.0);
    delta_phi_deta1_gap->SetMaximum(1.2);
    delta_phi_deta1_gap->SetLineColor(2);
    delta_phi_deta1_gap->SetLineStyle(1);
    delta_phi_deta1_gap->SetLineWidth(3);
    delta_phi_deta1_gap->Draw("hist");
    
    TLegend *leg7 = new TLegend(0.40,0.88,0.98,0.98);
    leg7->AddEntry(delta_phi_deta1_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg7->SetFillColor(0);
    leg7->Draw();

     c7->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1ag2 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta2_gap");
     TH1D *x2ag2 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta2_gap");
     TH1D *x1bg2 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta2_gap");
     TH1D *x2bg2 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta2_gap");
     TH1D *x1cg2 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta2_gap");
     TH1D *x2cg2 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta2_gap");
     
     enta = x1ag2->GetEntries()-x1ag2->GetNbinsX();
     entb = x1bg2->GetEntries()-x1bg2->GetNbinsX();
     entc = x1cg2->GetEntries()-x1cg2->GetNbinsX();

     int1a = x1ag2->Integral();
     int2a = x2ag2->Integral();
     factora = int2a/int1a;

     int1b = x1bg2->Integral();
     int2b = x2bg2->Integral();
     factorb = int2b/int1b;

     int1c = x1cg2->Integral();
     int2c = x2cg2->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;

    TH1D *delta_phi_deta2_gap, *delta_phi_deta2_gap2, *delta_phi_deta2_gap_old, *delta_phi_deta2_gap_a, *delta_phi_deta2_gap_b, *delta_phi_deta2_gap_c;
    delta_phi_deta2_gap =  new TH1D("Delta_phi_deta2_gap","delta_phi_deta2_gap;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap2 =  new TH1D("Delta_phi_deta2_gap2","delta_phi_deta2_gap2;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap_old =  new TH1D("Delta_phi_deta2_gap_old","delta_phi_deta2_gap_old;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap_a =  new TH1D("Delta_phi_deta2_gap_a","delta_phi_deta2_gap_a;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap_b =  new TH1D("Delta_phi_deta2_gap_b","delta_phi_deta2_gap_b;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap_c =  new TH1D("Delta_phi_deta2_gap_c","delta_phi_deta2_gap_c;|#Delta#phi| [rad] for 2.5 > |#Delta#eta| > 3.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_gap_a->Divide(x1ag2,x2ag2,1.,1./factora,"");
    delta_phi_deta2_gap_b->Divide(x1bg2,x2bg2,1.,1./factorb,"");
    delta_phi_deta2_gap_c->Divide(x1cg2,x2cg2,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta2_gap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta2_gap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta2_gap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta2_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta2_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta2_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta2_gap2, delta_phi_deta2_gap, step);

    delta_phi_deta2_gap->SetMinimum(0.0);
    delta_phi_deta2_gap->SetMaximum(1.2);
    delta_phi_deta2_gap->SetLineColor(2);
    delta_phi_deta2_gap->SetLineStyle(1);
    delta_phi_deta2_gap->SetLineWidth(3);
    delta_phi_deta2_gap->Draw("hist");
    
    c7->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1ag3 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta3_gap");
     TH1D *x2ag3 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta3_gap");
     TH1D *x1bg3 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta3_gap");
     TH1D *x2bg3 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta3_gap");
     TH1D *x1cg3 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta3_gap");
     TH1D *x2cg3 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta3_gap");
     
     enta = x1ag3->GetEntries()-x1ag3->GetNbinsX();
     entb = x1bg3->GetEntries()-x1bg3->GetNbinsX();
     entc = x1cg3->GetEntries()-x1cg3->GetNbinsX();

     int1a = x1ag3->Integral();
     int2a = x2ag3->Integral();
     factora = int2a/int1a;

     int1b = x1bg3->Integral();
     int2b = x2bg3->Integral();
     factorb = int2b/int1b;

     int1c = x1cg3->Integral();
     int2c = x2cg3->Integral();
     factorc = int2c/int1c;
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta3_gap, *delta_phi_deta3_gap2, *delta_phi_deta3_gap_old, *delta_phi_deta3_gap_a, *delta_phi_deta3_gap_b, *delta_phi_deta3_gap_c;
    delta_phi_deta3_gap =  new TH1D("Delta_phi_deta3_gap","delta_phi_deta3_gap;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap2 =  new TH1D("Delta_phi_deta3_gap2","delta_phi_deta3_gap2;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap_old =  new TH1D("Delta_phi_deta3_gap_old","delta_phi_deta3_gap_old;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap_a =  new TH1D("Delta_phi_deta3_gap_a","delta_phi_deta3_gap_a;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap_b =  new TH1D("Delta_phi_deta3_gap_b","delta_phi_deta3_gap_b;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap_c =  new TH1D("Delta_phi_deta3_gap_c","delta_phi_deta3_gap_c;|#Delta#phi| [rad] for 3.5 > |#Delta#eta| > 4.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_gap_a->Divide(x1ag3,x2ag3,1.,1./factora,"");
    delta_phi_deta3_gap_b->Divide(x1bg3,x2bg3,1.,1./factorb,"");
    delta_phi_deta3_gap_c->Divide(x1cg3,x2cg3,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta3_gap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta3_gap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta3_gap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta3_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta3_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta3_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta3_gap2, delta_phi_deta3_gap, step);

    delta_phi_deta3_gap->SetMinimum(0.0);
    delta_phi_deta3_gap->SetMaximum(1.2);
    delta_phi_deta3_gap->SetLineColor(2);
    delta_phi_deta3_gap->SetLineStyle(1);
    delta_phi_deta3_gap->SetLineWidth(3);
    delta_phi_deta3_gap->Draw("hist");
    
     c7->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1ag4 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta4_gap");
     TH1D *x2ag4 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta4_gap");
     TH1D *x1bg4 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta4_gap");
     TH1D *x2bg4 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta4_gap");
     TH1D *x1cg4 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta4_gap");
     TH1D *x2cg4 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta4_gap");
     
     enta = x1ag4->GetEntries()-x1ag4->GetNbinsX();
     entb = x1bg4->GetEntries()-x1bg4->GetNbinsX();
     entc = x1cg4->GetEntries()-x1cg4->GetNbinsX();

     int1a = x1ag4->Integral();
     int2a = x2ag4->Integral();
     factora = int2a/int1a;

     int1b = x1bg4->Integral();
     int2b = x2bg4->Integral();
     factorb = int2b/int1b;

     int1c = x1cg4->Integral();
     int2c = x2cg4->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta4_gap, *delta_phi_deta4_gap2, *delta_phi_deta4_gap_old, *delta_phi_deta4_gap_a, *delta_phi_deta4_gap_b, *delta_phi_deta4_gap_c;
    delta_phi_deta4_gap =  new TH1D("Delta_phi_deta4_gap","delta_phi_deta4_gap;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap2 =  new TH1D("Delta_phi_deta4_gap2","delta_phi_deta4_gap2;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap_old =  new TH1D("Delta_phi_deta4_gap_old","delta_phi_deta4_gap_old;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap_a =  new TH1D("Delta_phi_deta4_gap_a","delta_phi_deta4_gap_a;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap_b =  new TH1D("Delta_phi_deta4_gap_b","delta_phi_deta4_gap_b;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap_c =  new TH1D("Delta_phi_deta4_gap_c","delta_phi_deta4_gap_c;|#Delta#phi| [rad] for 4.5 > |#Delta#eta| > 7.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_gap_a->Divide(x1ag4,x2ag4,1.,1./factora,"");
    delta_phi_deta4_gap_b->Divide(x1bg4,x2bg4,1.,1./factorb,"");
    delta_phi_deta4_gap_c->Divide(x1cg4,x2cg4,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta4_gap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta4_gap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta4_gap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta4_gap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta4_gap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta4_gap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta4_gap2, delta_phi_deta4_gap, step);

    delta_phi_deta4_gap->SetMinimum(0.0);
    delta_phi_deta4_gap->SetMaximum(1.2);
    delta_phi_deta4_gap->SetLineColor(2);
    delta_phi_deta4_gap->SetLineStyle(1);
    delta_phi_deta4_gap->SetLineWidth(3);
    delta_phi_deta4_gap->Draw("hist"); 
    
    
  /*   c6->cd(5);
     TH1D *x1ag5 = (TH1D*) data_file1->Get("ak5PF_delta_phi_deta5_gap");
     TH1D *x2ag5 = (TH1D*) data_file2->Get("ak5PF_delta_phi_deta5_gap");
     
     double int1 = x1ag5->Integral();
     double int2 = x2ag5->Integral();
     double factor = int2/int1;
     
     //cout<<"factor"<<factor<<endl;
    
    delta_phi_deta5_gap =  new TH1D("Delta_phi_deta5_gap","delta_phi_deta5_gap;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta5_gap_b =  new TH1D("Delta_phi_deta5_gap_b","delta_phi_deta5_gap_b;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta5_gap_b->Divide(x1ag5,x2ag5,1.,1./factor,"");
    
    for(Int_t i=1;i<=delta_phi_deta5_gap_b->GetNbinsX();i++)
    {
    Float_t cont = delta_phi_deta5_gap_b->GetBinContent(i) - 1;
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min) { min = cont;}
    delta_phi_deta5_gap->SetBinContent(i,cont);
    }
    
    delta_phi_deta5_gap->SetLineColor(2);
    delta_phi_deta5_gap->SetLineStyle(1);
    delta_phi_deta5_gap->SetLineWidth(3);
    delta_phi_deta5_gap->Draw("hist");
    
     c6->cd(6);
     TH1D *x1ag6 = (TH1D*) data_file1->Get("ak5PF_delta_phi_deta6_gap");
     TH1D *x2ag6 = (TH1D*) data_file2->Get("ak5PF_delta_phi_deta6_gap");
     
     double int1 = x1ag6->Integral();
     double int2 = x2ag6->Integral();
     double factor = int2/int1;
     
     //cout<<"factor"<<factor<<endl;
    
    delta_phi_deta6_gap =  new TH1D("Delta_phi_deta6_gap","delta_phi_deta6_gap;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta6_gap_b =  new TH1D("Delta_phi_deta6_gap_b","delta_phi_deta6_gap_b;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta6_gap_b->Divide(x1ag6,x2ag6,1.,1./factor,"");
    
    for(Int_t i=1;i<=delta_phi_deta6_gap_b->GetNbinsX();i++)
    {
    Float_t cont = delta_phi_deta6_gap_b->GetBinContent(i) - 1;
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min) { min = cont;}
    delta_phi_deta6_gap->SetBinContent(i,cont);
    }
    
    delta_phi_deta6_gap->SetLineColor(2);
    delta_phi_deta6_gap->SetLineStyle(1);
    delta_phi_deta6_gap->SetLineWidth(3);
    delta_phi_deta6_gap->Draw("hist"); */
    
    ave = tot/(delta_phi_deta1_gap->GetNbinsX()*4);
    ave_old = tot_old/(delta_phi_deta1_gap->GetNbinsX()*4);
    cout<<"Gap deta: min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[3] = min*100;
    max_row[3] = max*100;
    ave_row[3] = ave*100;
    old_row[3] = ave_old*100; 

    c7->Print("output/uncertainities/png/pile_up_delta_phi_deta_gap.png");
    c7->Print("output/uncertainities/c/pile_up_delta_phi_deta_gap.C");
    c7->Print("output/uncertainities/eps/pile_up_delta_phi_deta_gap.eps");
    c7->Close();

    TCanvas *c7o = new TCanvas("c7o","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c7o->Divide(2,2);

    c7o->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta1_gap->SetMinimum(0.0);
    delta_phi_deta1_gap->Draw("hist");
    delta_phi_deta1_gap2->SetLineColor(2);
    delta_phi_deta1_gap2->SetLineStyle(3);
    delta_phi_deta1_gap2->SetLineWidth(3);
    delta_phi_deta1_gap2->Draw("hist same");
    delta_phi_deta1_gap_old->SetLineColor(2);
    delta_phi_deta1_gap_old->SetLineStyle(2);
    delta_phi_deta1_gap_old->SetLineWidth(3);
    delta_phi_deta1_gap_old->Draw("hist same");

    TLegend *leg7o = new TLegend(0.30,0.80,0.98,0.98);
    leg7o->AddEntry(delta_phi_deta1_gap,"Average Pile-up uncertainty for 2010 Data","l");
    leg7o->AddEntry(delta_phi_deta1_gap2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg7o->AddEntry(delta_phi_deta1_gap_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg7o->SetFillColor(0);
    leg7o->Draw();

    c7o->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta2_gap->SetMinimum(0.0);
    delta_phi_deta2_gap->Draw("hist");
    delta_phi_deta2_gap2->SetLineColor(2);
    delta_phi_deta2_gap2->SetLineStyle(3);
    delta_phi_deta2_gap2->SetLineWidth(3);
    delta_phi_deta2_gap2->Draw("hist same");
    delta_phi_deta2_gap_old->SetLineColor(2);
    delta_phi_deta2_gap_old->SetLineStyle(2);
    delta_phi_deta2_gap_old->SetLineWidth(3);
    delta_phi_deta2_gap_old->Draw("hist same");

    c7o->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta3_gap->SetMinimum(0.0);
    delta_phi_deta3_gap->Draw("hist");
    delta_phi_deta3_gap2->SetLineColor(2);
    delta_phi_deta3_gap2->SetLineStyle(3);
    delta_phi_deta3_gap2->SetLineWidth(3);
    delta_phi_deta3_gap2->Draw("hist same");
    delta_phi_deta3_gap_old->SetLineColor(2);
    delta_phi_deta3_gap_old->SetLineStyle(2);
    delta_phi_deta3_gap_old->SetLineWidth(3);
    delta_phi_deta3_gap_old->Draw("hist same");

    c7o->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta4_gap->SetMinimum(0.0);
    delta_phi_deta4_gap->Draw("hist");
    delta_phi_deta4_gap2->SetLineColor(2);
    delta_phi_deta4_gap2->SetLineStyle(3);
    delta_phi_deta4_gap2->SetLineWidth(3);
    delta_phi_deta4_gap2->Draw("hist same");
    delta_phi_deta4_gap_old->SetLineColor(2);
    delta_phi_deta4_gap_old->SetLineStyle(2);
    delta_phi_deta4_gap_old->SetLineWidth(3);
    delta_phi_deta4_gap_old->Draw("hist same");

    c7o->Print("output/uncertainities/png/pile_up_delta_phi_deta_gap_a.png");
    c7o->Print("output/uncertainities/c/pile_up_delta_phi_deta_gap_a.C");
    c7o->Print("output/uncertainities/eps/pile_up_delta_phi_deta_gap_a.eps");
    c7o->Close();

    TCanvas *c7b = new TCanvas("c7b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c7b->Divide(2,2);

    c7b->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta1_gap_a->SetMaximum(2.0);
    delta_phi_deta1_gap_a->SetMinimum(0.0);
    delta_phi_deta1_gap_a->SetLineColor(2);
    delta_phi_deta1_gap_a->SetLineStyle(1);
    delta_phi_deta1_gap_a->SetLineWidth(3);
    delta_phi_deta1_gap_a->Draw("e1");
    delta_phi_deta1_gap_b->SetLineColor(3);
    delta_phi_deta1_gap_b->SetLineStyle(1);
    delta_phi_deta1_gap_b->SetLineWidth(3);
    delta_phi_deta1_gap_b->Draw("e1 same");
    delta_phi_deta1_gap_c->SetLineColor(4);
    delta_phi_deta1_gap_c->SetLineStyle(1);
    delta_phi_deta1_gap_c->SetLineWidth(3);
    delta_phi_deta1_gap_c->Draw("e1 same");

    TLegend *leg7b = new TLegend(0.5,0.75,0.98,0.98);
    leg7b->AddEntry(delta_phi_deta1_gap_a,labela,"l");
    leg7b->AddEntry(delta_phi_deta1_gap_b,labelb,"l");
    leg7b->AddEntry(delta_phi_deta1_gap_c,labelc,"l");
    leg7b->SetFillColor(0);
    leg7b->Draw();

    c7b->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta2_gap_a->SetMaximum(2.0);
    delta_phi_deta2_gap_a->SetMinimum(0.0);
    delta_phi_deta2_gap_a->SetLineColor(2);
    delta_phi_deta2_gap_a->SetLineStyle(1);
    delta_phi_deta2_gap_a->SetLineWidth(3);
    delta_phi_deta2_gap_a->Draw("e1");
    delta_phi_deta2_gap_b->SetLineColor(3);
    delta_phi_deta2_gap_b->SetLineStyle(1);
    delta_phi_deta2_gap_b->SetLineWidth(3);
    delta_phi_deta2_gap_b->Draw("e1 same");
    delta_phi_deta2_gap_c->SetLineColor(4);
    delta_phi_deta2_gap_c->SetLineStyle(1);
    delta_phi_deta2_gap_c->SetLineWidth(3);
    delta_phi_deta2_gap_c->Draw("e1 same");

    c7b->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta3_gap_a->SetMaximum(2.0);
    delta_phi_deta3_gap_a->SetMinimum(0.0);
    delta_phi_deta3_gap_a->SetLineColor(2);
    delta_phi_deta3_gap_a->SetLineStyle(1);
    delta_phi_deta3_gap_a->SetLineWidth(3);
    delta_phi_deta3_gap_a->Draw("e1");
    delta_phi_deta3_gap_b->SetLineColor(3);
    delta_phi_deta3_gap_b->SetLineStyle(1);
    delta_phi_deta3_gap_b->SetLineWidth(3);
    delta_phi_deta3_gap_b->Draw("e1 same");
    delta_phi_deta3_gap_c->SetLineColor(4);
    delta_phi_deta3_gap_c->SetLineStyle(1);
    delta_phi_deta3_gap_c->SetLineWidth(3);
    delta_phi_deta3_gap_c->Draw("e1 same");

    c7b->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta4_gap_a->SetMaximum(2.0);
    delta_phi_deta4_gap_a->SetMinimum(0.0);
    delta_phi_deta4_gap_a->SetLineColor(2);
    delta_phi_deta4_gap_a->SetLineStyle(1);
    delta_phi_deta4_gap_a->SetLineWidth(3);
    delta_phi_deta4_gap_a->Draw("e1");
    delta_phi_deta4_gap_b->SetLineColor(3);
    delta_phi_deta4_gap_b->SetLineStyle(1);
    delta_phi_deta4_gap_b->SetLineWidth(3);
    delta_phi_deta4_gap_b->Draw("e1 same");
    delta_phi_deta4_gap_c->SetLineColor(4);
    delta_phi_deta4_gap_c->SetLineStyle(1);
    delta_phi_deta4_gap_c->SetLineWidth(3);
    delta_phi_deta4_gap_c->Draw("e1 same");

    c7b->Print("output/uncertainities/png/pile_up_delta_phi_deta_gap_b.png");
    c7b->Print("output/uncertainities/c/pile_up_delta_phi_deta_gap_b.C");
    c7b->Print("output/uncertainities/eps/pile_up_delta_phi_deta_gap_b.eps");
    c7b->Close();

    TCanvas *c8 = new TCanvas("c8","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c8->Divide(2,2);
    
     c8->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1an1 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta1_nogap");
     TH1D *x2an1 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta1_nogap");
     TH1D *x1bn1 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta1_nogap");
     TH1D *x2bn1 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta1_nogap");
     TH1D *x1cn1 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta1_nogap");
     TH1D *x2cn1 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta1_nogap");
     
     enta = x1an1->GetEntries()-x1an1->GetNbinsX();
     entb = x1bn1->GetEntries()-x1bn1->GetNbinsX();
     entc = x1cn1->GetEntries()-x1cn1->GetNbinsX();

     int1a = x1an1->Integral();
     int2a = x2an1->Integral();
     factora = int2a/int1a;

     int1b = x1bn1->Integral();
     int2b = x2bn1->Integral();
     factorb = int2b/int1b;

     int1c = x1cn1->Integral();
     int2c = x2cn1->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta1_nogap, *delta_phi_deta1_nogap2, *delta_phi_deta1_nogap_old, *delta_phi_deta1_nogap_a, *delta_phi_deta1_nogap_b, *delta_phi_deta1_nogap_c;
    delta_phi_deta1_nogap =  new TH1D("Delta_phi_deta1_nogap","delta_phi_deta1_nogap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap2 =  new TH1D("Delta_phi_deta1_nogap2","delta_phi_deta1_nogap2;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap_old =  new TH1D("Delta_phi_deta1_nogap_old","delta_phi_deta1_nogap_old;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap_a =  new TH1D("Delta_phi_deta1_nogap_a","delta_phi_deta1_nogap_a;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap_b =  new TH1D("Delta_phi_deta1_nogap_b","delta_phi_deta1_nogap_b;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap_c =  new TH1D("Delta_phi_deta1_nogap_c","delta_phi_deta1_nogap_c;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta1_nogap_a->Divide(x1an1,x2an1,1.,1./factora,"");
    delta_phi_deta1_nogap_b->Divide(x1bn1,x2bn1,1.,1./factorb,"");
    delta_phi_deta1_nogap_c->Divide(x1cn1,x2cn1,1.,1./factorc,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_phi_deta1_nogap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta1_nogap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta1_nogap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta1_nogap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta1_nogap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta1_nogap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta1_nogap2, delta_phi_deta1_nogap, step);

    delta_phi_deta1_nogap->SetMinimum(0.0);
    delta_phi_deta1_nogap->SetMaximum(0.35);
    delta_phi_deta1_nogap->SetLineColor(2);
    delta_phi_deta1_nogap->SetLineStyle(1);
    delta_phi_deta1_nogap->SetLineWidth(3);
    delta_phi_deta1_nogap->Draw("hist");
    
    TLegend *leg8 = new TLegend(0.40,0.88,0.98,0.98);
    leg8->AddEntry(delta_phi_deta1_nogap,"Average Pile-up uncertainty for 2010 Data","l");
    leg8->SetFillColor(0);
    leg8->Draw();

     c8->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1an2 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta2_nogap");
     TH1D *x2an2 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta2_nogap");
     TH1D *x1bn2 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta2_nogap");
     TH1D *x2bn2 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta2_nogap");
     TH1D *x1cn2 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta2_nogap");
     TH1D *x2cn2 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta2_nogap");
     
     enta = x1an2->GetEntries()-x1an2->GetNbinsX();
     entb = x1bn2->GetEntries()-x1bn2->GetNbinsX();
     entc = x1cn2->GetEntries()-x1cn2->GetNbinsX();

     int1a = x1an2->Integral();
     int2a = x2an2->Integral();
     factora = int2a/int1a;

     int1b = x1bn2->Integral();
     int2b = x2bn2->Integral();
     factorb = int2b/int1b;

     int1c = x1cn2->Integral();
     int2c = x2cn2->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta2_nogap, *delta_phi_deta2_nogap2, *delta_phi_deta2_nogap_old, *delta_phi_deta2_nogap_a, *delta_phi_deta2_nogap_b, *delta_phi_deta2_nogap_c;
    delta_phi_deta2_nogap =  new TH1D("Delta_phi_deta2_nogap","delta_phi_deta2_nogap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap2 =  new TH1D("Delta_phi_deta2_nogap2","delta_phi_deta2_nogap2;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap_old =  new TH1D("Delta_phi_deta2_nogap_old","delta_phi_deta2_nogap_old;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap_a =  new TH1D("Delta_phi_deta2_nogap_a","delta_phi_deta2_nogap_a;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap_b =  new TH1D("Delta_phi_deta2_nogap_b","delta_phi_deta2_nogap_b;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap_c =  new TH1D("Delta_phi_deta2_nogap_c","delta_phi_deta2_nogap_c;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta2_nogap_a->Divide(x1an2,x2an2,1.,1./factora,"");
    delta_phi_deta2_nogap_b->Divide(x1bn2,x2bn2,1.,1./factorb,"");
    delta_phi_deta2_nogap_c->Divide(x1cn2,x2cn2,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta2_nogap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta2_nogap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta2_nogap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta2_nogap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta2_nogap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta2_nogap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta2_nogap2, delta_phi_deta2_nogap, step);

    delta_phi_deta2_nogap->SetMinimum(0.0);
    delta_phi_deta2_nogap->SetMaximum(0.35);
    delta_phi_deta2_nogap->SetLineColor(2);
    delta_phi_deta2_nogap->SetLineStyle(1);
    delta_phi_deta2_nogap->SetLineWidth(3);
    delta_phi_deta2_nogap->Draw("hist");
    
    c8->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1an3 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta3_nogap");
     TH1D *x2an3 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta3_nogap");
     TH1D *x1bn3 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta3_nogap");
     TH1D *x2bn3 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta3_nogap");
     TH1D *x1cn3 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta3_nogap");
     TH1D *x2cn3 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta3_nogap");
     
     enta = x1an3->GetEntries()-x1an3->GetNbinsX();
     entb = x1bn3->GetEntries()-x1bn3->GetNbinsX();
     entc = x1cn3->GetEntries()-x1cn3->GetNbinsX();

     int1a = x1an3->Integral();
     int2a = x2an3->Integral();
     factora = int2a/int1a;

     int1b = x1bn3->Integral();
     int2b = x2bn3->Integral();
     factorb = int2b/int1b;

     int1c = x1cn3->Integral();
     int2c = x2cn3->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta3_nogap, *delta_phi_deta3_nogap2, *delta_phi_deta3_nogap_old, *delta_phi_deta3_nogap_a, *delta_phi_deta3_nogap_b, *delta_phi_deta3_nogap_c;
    delta_phi_deta3_nogap =  new TH1D("Delta_phi_deta3_nogap","delta_phi_deta3_nogap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap2 =  new TH1D("Delta_phi_deta3_nogap2","delta_phi_deta3_nogap2;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap_old =  new TH1D("Delta_phi_deta3_nogap_old","delta_phi_deta3_nogap_old;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap_a =  new TH1D("Delta_phi_deta3_nogap_a","delta_phi_deta3_nogap_a;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap_b =  new TH1D("Delta_phi_deta3_nogap_b","delta_phi_deta3_nogap_b;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap_c =  new TH1D("Delta_phi_deta3_nogap_c","delta_phi_deta3_nogap_c;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta3_nogap_a->Divide(x1an3,x2an3,1.,1./factora,"");
    delta_phi_deta3_nogap_b->Divide(x1bn3,x2bn3,1.,1./factorb,"");
    delta_phi_deta3_nogap_c->Divide(x1cn3,x2cn3,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta3_nogap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta3_nogap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta3_nogap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta3_nogap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta3_nogap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta3_nogap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta3_nogap2, delta_phi_deta3_nogap, step);

    delta_phi_deta3_nogap->SetMinimum(0.0);
    delta_phi_deta3_nogap->SetMaximum(0.35);
    delta_phi_deta3_nogap->SetLineColor(2);
    delta_phi_deta3_nogap->SetLineStyle(1);
    delta_phi_deta3_nogap->SetLineWidth(3);
    delta_phi_deta3_nogap->Draw("hist");
    
     c8->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

     TH1D *x1an4 = (TH1D*) data_file1a->Get("ak5PF_delta_phi_deta4_nogap");
     TH1D *x2an4 = (TH1D*) data_file2a->Get("ak5PF_delta_phi_deta4_nogap");
     TH1D *x1bn4 = (TH1D*) data_file1b->Get("ak5PF_delta_phi_deta4_nogap");
     TH1D *x2bn4 = (TH1D*) data_file2b->Get("ak5PF_delta_phi_deta4_nogap");
     TH1D *x1cn4 = (TH1D*) data_file1c->Get("ak5PF_delta_phi_deta4_nogap");
     TH1D *x2cn4 = (TH1D*) data_file2c->Get("ak5PF_delta_phi_deta4_nogap");
     
     enta = x1an4->GetEntries()-x1an4->GetNbinsX();
     entb = x1bn4->GetEntries()-x1bn4->GetNbinsX();
     entc = x1cn4->GetEntries()-x1cn4->GetNbinsX();

     int1a = x1an4->Integral();
     int2a = x2an4->Integral();
     factora = int2a/int1a;

     int1b = x1bn4->Integral();
     int2b = x2bn4->Integral();
     factorb = int2b/int1b;

     int1c = x1cn4->Integral();
     int2c = x2cn4->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_phi_deta4_nogap, *delta_phi_deta4_nogap2, *delta_phi_deta4_nogap_old, *delta_phi_deta4_nogap_a, *delta_phi_deta4_nogap_b, *delta_phi_deta4_nogap_c;
    delta_phi_deta4_nogap =  new TH1D("Delta_phi_deta4_nogap","delta_phi_deta4_nogap;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap2 =  new TH1D("Delta_phi_deta4_nogap2","delta_phi_deta4_nogap2;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap_old =  new TH1D("Delta_phi_deta4_nogap_old","delta_phi_deta4_nogap_old;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap_a =  new TH1D("Delta_phi_deta4_nogap_a","delta_phi_deta4_nogap_a;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap_b =  new TH1D("Delta_phi_deta4_nogap_b","delta_phi_deta4_nogap_b;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap_c =  new TH1D("Delta_phi_deta4_nogap_c","delta_phi_deta4_nogap_c;|#Delta#phi| [rad] for 0.4 > |#Delta#eta| > 2.5;Ratio to all vertexes selection", dphi_nbins, dphi_bins);
    delta_phi_deta4_nogap_a->Divide(x1an4,x2an4,1.,1./factora,"");
    delta_phi_deta4_nogap_b->Divide(x1bn4,x2bn4,1.,1./factorb,"");
    delta_phi_deta4_nogap_c->Divide(x1cn4,x2cn4,1.,1./factorc,"");
    
    for(Int_t i=1;i<=delta_phi_deta4_nogap_a->GetNbinsX();i++)
    {
    Float_t conta = delta_phi_deta4_nogap_a->GetBinContent(i) - 1;
    Float_t contb = delta_phi_deta4_nogap_b->GetBinContent(i) - 1;
    Float_t contc = delta_phi_deta4_nogap_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_phi_deta4_nogap2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    delta_phi_deta4_nogap_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_phi_deta4_nogap2, delta_phi_deta4_nogap, step);

    delta_phi_deta4_nogap->SetMinimum(0.0);
    delta_phi_deta4_nogap->SetMaximum(0.35);
    delta_phi_deta4_nogap->SetLineColor(2);
    delta_phi_deta4_nogap->SetLineStyle(1);
    delta_phi_deta4_nogap->SetLineWidth(3);
    delta_phi_deta4_nogap->Draw("hist"); 
    
    
   /*  c7->cd(5);
     TH1D *x1an5 = (TH1D*) data_file1->Get("ak5PF_delta_phi_deta5_nogap");
     TH1D *x2an5 = (TH1D*) data_file2->Get("ak5PF_delta_phi_deta5_nogap");
     
     double int1 = x1an5->Integral();
     double int2 = x2an5->Integral();
     double factor = int2/int1;
     
     //cout<<"factor"<<factor<<endl;
    
    delta_phi_deta5_nogap =  new TH1D("Delta_phi_deta5_nogap","delta_phi_deta5_nogap;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta5_nogap_b =  new TH1D("Delta_phi_deta5_nogap_b","delta_phi_deta5_nogap_b;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta5_nogap_b->Divide(x1an5,x2an5,1.,1./factor,"");
    
    for(Int_t i=1;i<=delta_phi_deta5_nogap_b->GetNbinsX();i++)
    {
    Float_t cont = delta_phi_deta5_nogap_b->GetBinContent(i) - 1;
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min) { min = cont;}
    delta_phi_deta5_nogap->SetBinContent(i,cont);
    }
    
    delta_phi_deta5_nogap->SetLineColor(2);
    delta_phi_deta5_nogap->SetLineStyle(1);
    delta_phi_deta5_nogap->SetLineWidth(3);
    delta_phi_deta5_nogap->Draw("hist");
    
     c7->cd(6);
     TH1D *x1an6 = (TH1D*) data_file1->Get("ak5PF_delta_phi_deta6_nogap");
     TH1D *x2an6 = (TH1D*) data_file2->Get("ak5PF_delta_phi_deta6_nogap");
     
     double int1 = x1an6->Integral();
     double int2 = x2an6->Integral();
     double factor = int2/int1;
     
     //cout<<"factor"<<factor<<endl;
    
    delta_phi_deta6_nogap =  new TH1D("Delta_phi_deta6_nogap","delta_phi_deta6_nogap;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta6_nogap_b =  new TH1D("Delta_phi_deta6_nogap_b","delta_phi_deta6_nogap_b;#Delta#phi;Pile Up Uncertainty", dphi_nbins, dphi_bins);
    delta_phi_deta6_nogap_b->Divide(x1an6,x2an6,1.,1./factor,"");
    
    for(Int_t i=1;i<=delta_phi_deta6_nogap_b->GetNbinsX();i++)
    {
    Float_t cont = delta_phi_deta6_nogap_b->GetBinContent(i) - 1;
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min) { min = cont;}
    delta_phi_deta6_nogap->SetBinContent(i,cont);
    }
    
    delta_phi_deta6_nogap->SetLineColor(2);
    delta_phi_deta6_nogap->SetLineStyle(1);
    delta_phi_deta6_nogap->SetLineWidth(3);
    delta_phi_deta6_nogap->Draw("hist"); */
    
    ave = tot/(delta_phi_deta1_nogap->GetNbinsX()*4);
    ave_old = tot_old/(delta_phi_deta1_nogap->GetNbinsX()*4);
    cout<<"No gap deta: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[5] = min*100;
    max_row[5] = max*100;
    ave_row[5] = ave*100;
    old_row[5] = ave_old*100;
    
    c8->Print("output/uncertainities/png/pile_up_delta_phi_deta_nogap.png");
    c8->Print("output/uncertainities/c/pile_up_delta_phi_deta_nogap.C");
    c8->Print("output/uncertainities/eps/pile_up_delta_phi_deta_nogap.eps");
    c8->Close();

    TCanvas *c8o = new TCanvas("c8o","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c8o->Divide(2,2);

    c8o->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta1_nogap->SetMinimum(0.0);
    delta_phi_deta1_nogap->Draw("hist");
    delta_phi_deta1_nogap2->SetLineColor(2);
    delta_phi_deta1_nogap2->SetLineStyle(3);
    delta_phi_deta1_nogap2->SetLineWidth(3);
    delta_phi_deta1_nogap2->Draw("hist same");
    delta_phi_deta1_nogap_old->SetLineColor(2);
    delta_phi_deta1_nogap_old->SetLineStyle(2);
    delta_phi_deta1_nogap_old->SetLineWidth(3);
    delta_phi_deta1_nogap_old->Draw("hist same");

    TLegend *leg8o = new TLegend(0.30,0.80,0.98,0.98);
    leg8o->AddEntry(delta_phi_deta1_nogap,"Average Pile-up uncertainty for 2010 Data","l");
    leg8o->AddEntry(delta_phi_deta1_nogap2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg8o->AddEntry(delta_phi_deta1_nogap_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg8o->SetFillColor(0);
    leg8o->Draw();

    c8o->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta2_nogap->SetMinimum(0.0);
    delta_phi_deta2_nogap->Draw("hist");
    delta_phi_deta2_nogap2->SetLineColor(2);
    delta_phi_deta2_nogap2->SetLineStyle(3);
    delta_phi_deta2_nogap2->SetLineWidth(3);
    delta_phi_deta2_nogap2->Draw("hist same");
    delta_phi_deta2_nogap_old->SetLineColor(2);
    delta_phi_deta2_nogap_old->SetLineStyle(2);
    delta_phi_deta2_nogap_old->SetLineWidth(3);
    delta_phi_deta2_nogap_old->Draw("hist same");

    c8o->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta3_nogap->SetMinimum(0.0);
    delta_phi_deta3_nogap->Draw("hist");
    delta_phi_deta3_nogap2->SetLineColor(2);
    delta_phi_deta3_nogap2->SetLineStyle(3);
    delta_phi_deta3_nogap2->SetLineWidth(3);
    delta_phi_deta3_nogap2->Draw("hist same");
    delta_phi_deta3_nogap_old->SetLineColor(2);
    delta_phi_deta3_nogap_old->SetLineStyle(2);
    delta_phi_deta3_nogap_old->SetLineWidth(3);
    delta_phi_deta3_nogap_old->Draw("hist same");

    c8o->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta4_nogap->SetMinimum(0.0);
    delta_phi_deta4_nogap->Draw("hist");
    delta_phi_deta4_nogap2->SetLineColor(2);
    delta_phi_deta4_nogap2->SetLineStyle(3);
    delta_phi_deta4_nogap2->SetLineWidth(3);
    delta_phi_deta4_nogap2->Draw("hist same");
    delta_phi_deta4_nogap_old->SetLineColor(2);
    delta_phi_deta4_nogap_old->SetLineStyle(2);
    delta_phi_deta4_nogap_old->SetLineWidth(3);
    delta_phi_deta4_nogap_old->Draw("hist same");

    c8o->Print("output/uncertainities/png/pile_up_delta_phi_deta_nogap_a.png");
    c8o->Print("output/uncertainities/c/pile_up_delta_phi_deta_nogap_a.C");
    c8o->Print("output/uncertainities/eps/pile_up_delta_phi_deta_nogap_a.eps");
    c8o->Close();

    TCanvas *c8b = new TCanvas("c8b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    c8b->Divide(2,2);

    c8b->cd(1);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta1_nogap_a->SetMaximum(2.0);
    delta_phi_deta1_nogap_a->SetMinimum(0.0);
    delta_phi_deta1_nogap_a->SetLineColor(2);
    delta_phi_deta1_nogap_a->SetLineStyle(1);
    delta_phi_deta1_nogap_a->SetLineWidth(3);
    delta_phi_deta1_nogap_a->Draw("e1");
    delta_phi_deta1_nogap_b->SetLineColor(3);
    delta_phi_deta1_nogap_b->SetLineStyle(1);
    delta_phi_deta1_nogap_b->SetLineWidth(3);
    delta_phi_deta1_nogap_b->Draw("e1 same");
    delta_phi_deta1_nogap_c->SetLineColor(4);
    delta_phi_deta1_nogap_c->SetLineStyle(1);
    delta_phi_deta1_nogap_c->SetLineWidth(3);
    delta_phi_deta1_nogap_c->Draw("e1 same");

    TLegend *leg8b = new TLegend(0.5,0.75,0.98,0.98);
    leg8b->AddEntry(delta_phi_deta1_nogap_a,labela,"l");
    leg8b->AddEntry(delta_phi_deta1_nogap_b,labelb,"l");
    leg8b->AddEntry(delta_phi_deta1_nogap_c,labelc,"l");
    leg8b->SetFillColor(0);
    leg8b->Draw();

    c8b->cd(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta2_nogap_a->SetMaximum(2.0);
    delta_phi_deta2_nogap_a->SetMinimum(0.0);
    delta_phi_deta2_nogap_a->SetLineColor(2);
    delta_phi_deta2_nogap_a->SetLineStyle(1);
    delta_phi_deta2_nogap_a->SetLineWidth(3);
    delta_phi_deta2_nogap_a->Draw("e1");
    delta_phi_deta2_nogap_b->SetLineColor(3);
    delta_phi_deta2_nogap_b->SetLineStyle(1);
    delta_phi_deta2_nogap_b->SetLineWidth(3);
    delta_phi_deta2_nogap_b->Draw("e1 same");
    delta_phi_deta2_nogap_c->SetLineColor(4);
    delta_phi_deta2_nogap_c->SetLineStyle(1);
    delta_phi_deta2_nogap_c->SetLineWidth(3);
    delta_phi_deta2_nogap_c->Draw("e1 same");

    c8b->cd(3);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta3_nogap_a->SetMaximum(2.0);
    delta_phi_deta3_nogap_a->SetMinimum(0.0);
    delta_phi_deta3_nogap_a->SetLineColor(2);
    delta_phi_deta3_nogap_a->SetLineStyle(1);
    delta_phi_deta3_nogap_a->SetLineWidth(3);
    delta_phi_deta3_nogap_a->Draw("e1");
    delta_phi_deta3_nogap_b->SetLineColor(3);
    delta_phi_deta3_nogap_b->SetLineStyle(1);
    delta_phi_deta3_nogap_b->SetLineWidth(3);
    delta_phi_deta3_nogap_b->Draw("e1 same");
    delta_phi_deta3_nogap_c->SetLineColor(4);
    delta_phi_deta3_nogap_c->SetLineStyle(1);
    delta_phi_deta3_nogap_c->SetLineWidth(3);
    delta_phi_deta3_nogap_c->Draw("e1 same");

    c8b->cd(4);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);

    delta_phi_deta4_nogap_a->SetMaximum(2.0);
    delta_phi_deta4_nogap_a->SetMinimum(0.0);
    delta_phi_deta4_nogap_a->SetLineColor(2);
    delta_phi_deta4_nogap_a->SetLineStyle(1);
    delta_phi_deta4_nogap_a->SetLineWidth(3);
    delta_phi_deta4_nogap_a->Draw("e1");
    delta_phi_deta4_nogap_b->SetLineColor(3);
    delta_phi_deta4_nogap_b->SetLineStyle(1);
    delta_phi_deta4_nogap_b->SetLineWidth(3);
    delta_phi_deta4_nogap_b->Draw("e1 same");
    delta_phi_deta4_nogap_c->SetLineColor(4);
    delta_phi_deta4_nogap_c->SetLineStyle(1);
    delta_phi_deta4_nogap_c->SetLineWidth(3);
    delta_phi_deta4_nogap_c->Draw("e1 same");

    c8b->Print("output/uncertainities/png/pile_up_delta_phi_deta_nogap_b.png");
    c8b->Print("output/uncertainities/c/pile_up_delta_phi_deta_nogap_b.C");
    c8b->Print("output/uncertainities/eps/pile_up_delta_phi_deta_nogap_b.eps");
    c8b->Close();

    TCanvas *c9 = new TCanvas("c9","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x5c = (TH1D*) data_file1a->Get("ak5PF_leading_central_pt");
     TH1D *x6c = (TH1D*) data_file2a->Get("ak5PF_leading_central_pt");
     
     int1a = x5c->Integral();
     int2a = x6c->Integral();
     //factora = int2a/int1a;
     
     cout<<"factor"<<factora<<endl;
    
    TH1D *leading_central_pt, *leading_central_pt_b;
    leading_central_pt =  new TH1D("Leading_central_pt","Leading_central_pt;p_{T};Pile Up Uncertainty", cent_nbins, cent_bins);
    leading_central_pt->Divide(x5c,x6c,1.,1./factora,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;

    for(Int_t i=3;i<=leading_central_pt->GetNbinsX();i++)
    {
    Float_t cont = 1 - leading_central_pt->GetBinContent(i);
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min || i == 3) { min = cont;}
    }
    
    leading_central_pt->SetLineColor(2);
    leading_central_pt->SetLineStyle(1);
    leading_central_pt->SetLineWidth(3);
    leading_central_pt->Draw("hist");
    
    ave = tot/(leading_central_pt->GetNbinsX());
    cout<<"leading_pt_central: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl;
    
    c9->Print("output/uncertainities/png/pile_up_leading_central_pt.png");
    c9->Print("output/uncertainities/c/pile_up_leading_central_pt.C");
    c9->Print("output/uncertainities/eps/pile_up_leading_central_pt.eps");
    c9->Close();

    TCanvas *c10 = new TCanvas("c10","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x5f = (TH1D*) data_file1a->Get("ak5PF_leading_forward_pt");
     TH1D *x6f = (TH1D*) data_file2a->Get("ak5PF_leading_forward_pt");
     
     int1a = x5f->Integral();
     int2a = x6f->Integral();
     //factora = int2a/int1a;
     
     cout<<"factor"<<factora<<endl;
    
    TH1D *leading_forward_pt, *leading_forward_pt_b;
    leading_forward_pt =  new TH1D("Leading_forward_pt","Leading_forward_pt;p_{T};Pile Up Uncertainty", forw_nbins, forw_bins);
    leading_forward_pt->Divide(x5f,x6f,1.,1./factora,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=3;i<=leading_forward_pt->GetNbinsX();i++)
    {
    Float_t cont = 1 - leading_forward_pt->GetBinContent(i);
    if (cont < 0) { cont = -cont;}
    if (cont > max) {max = cont;}
    tot = tot + cont;
    if (cont < min || i == 3) { min = cont;}
    }
    
    leading_forward_pt->SetLineColor(2);
    leading_forward_pt->SetLineStyle(1);
    leading_forward_pt->SetLineWidth(3);
    leading_forward_pt->Draw("hist");
    
    ave = tot/(leading_forward_pt->GetNbinsX());
    cout<<"leading_forward_pt: min ="<<min*100<<", max = "<<max*100<<", ave= "<<ave*100<<endl;
    
    c10->Print("output/uncertainities/png/pile_up_leading_forward_pt.png");
    c10->Print("output/uncertainities/c/pile_up_leading_forward_pt.C");
    c10->Print("output/uncertainities/eps/pile_up_leading_forward_pt.eps");
    c10->Close();

    TCanvas *c11 = new TCanvas("c11","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x7a1 = (TH1D*) data_file1a->Get("ak5PF_leading_eta_star_inside_gap");
     TH1D *x8a1 = (TH1D*) data_file2a->Get("ak5PF_leading_eta_star_inside_gap");
     TH1D *x7b1 = (TH1D*) data_file1b->Get("ak5PF_leading_eta_star_inside_gap");
     TH1D *x8b1 = (TH1D*) data_file2b->Get("ak5PF_leading_eta_star_inside_gap");
     TH1D *x7c1 = (TH1D*) data_file1c->Get("ak5PF_leading_eta_star_inside_gap");
     TH1D *x8c1 = (TH1D*) data_file2c->Get("ak5PF_leading_eta_star_inside_gap");
     
     enta = x7a1->GetEntries()-x7a1->GetNbinsX();
     entb = x7b1->GetEntries()-x7b1->GetNbinsX();
     entc = x7c1->GetEntries()-x7c1->GetNbinsX();

     int1a = x7a1->Integral();
     int2a = x8a1->Integral();
     factora = int2a/int1a;

     int1b = x7b1->Integral();
     int2b = x8b1->Integral();
     factorb = int2b/int1b;

     int1c = x7c1->Integral();
     int2c = x8c1->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *eta_star, *eta_star2, *eta_star_old, *eta_star_a, *eta_star_b, *eta_star_c;
    eta_star =  new TH1D("Eta_star","eta_star;#eta*;Pile Up Uncertainty", etastar_nbins, etastar_bins);
    eta_star2 =  new TH1D("Eta_star2","eta_star2;#eta*;Pile Up Uncertainty", etastar_nbins, etastar_bins);
    eta_star_old =  new TH1D("Eta_star_old","eta_star_old;#eta*;Pile Up Uncertainty", etastar_nbins, etastar_bins);
    eta_star_a =  new TH1D("Eta_star_a","eta_star_a;#eta*;Ratio to all vertexes selection", etastar_nbins, etastar_bins);
    eta_star_b =  new TH1D("Eta_star_b","eta_star_b;#eta*;Ratio to all vertexes selection", etastar_nbins, etastar_bins);
    eta_star_c =  new TH1D("Eta_star_c","eta_star_c;#eta*;Ratio to all vertexes selection", etastar_nbins, etastar_bins);
    eta_star_a->Divide(x7a1,x8a1,1.,1./factora,"");
    eta_star_b->Divide(x7b1,x8b1,1.,1./factorb,"");
    eta_star_c->Divide(x7c1,x8c1,1.,1./factorc,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=eta_star_a->GetNbinsX();i++)
    {
    Float_t conta = eta_star_a->GetBinContent(i) - 1;
    Float_t contb = eta_star_b->GetBinContent(i) - 1;
    Float_t contc = eta_star_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    eta_star2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    //cout<<i<<" ="<<conta<<" . "<<contb<<" . "<<contc<<endl;
    eta_star_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(eta_star2, eta_star, step);

    eta_star->SetLineColor(2);
    eta_star->SetLineStyle(1);
    eta_star->SetLineWidth(3);
    eta_star->Draw("hist");

    TLegend *leg11 = new TLegend(0.13,0.9,0.55,0.97);
    leg11->AddEntry(eta_star,"Average Pile-up uncertainty for 2010 Data","l");
    leg11->SetFillColor(0);
    leg11->Draw();
    
    ave = tot/eta_star->GetNbinsX();
    ave_old = tot_old/eta_star->GetNbinsX();
    cout<<"Main selection: min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[8] = min*100;
    max_row[8] = max*100;
    ave_row[8] = ave*100;
    old_row[8] = ave_old*100;
    
    c11->Print("output/uncertainities/png/pile_up_eta_star.png");
    c11->Print("output/uncertainities/c/pile_up_eta_star.C");
    c11->Print("output/uncertainities/eps/pile_up_eta_star.eps");
    c11->Close();


    TCanvas *c11o = new TCanvas("c11","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    eta_star->SetMinimum(0.0);
    eta_star->SetMaximum(0.15);
    eta_star->Draw("hist");
    eta_star2->SetLineColor(2);
    eta_star2->SetLineStyle(3);
    eta_star2->SetLineWidth(3);
    eta_star2->Draw("hist same");
    eta_star_old->SetLineColor(2);
    eta_star_old->SetLineStyle(2);
    eta_star_old->SetLineWidth(3);
    eta_star_old->Draw("hist same");

    TLegend *leg11o = new TLegend(0.5,0.84,0.97,0.97);
    leg11o->AddEntry(eta_star,"Average Pile-up uncertainty for 2010 Data","l");
    leg11o->AddEntry(eta_star2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg11o->AddEntry(eta_star_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg11o->SetFillColor(0);
    leg11o->Draw();

    c11o->Print("output/uncertainities/png/pile_up_eta_star_a.png");
    c11o->Print("output/uncertainities/c/pile_up_eta_star_a.C");
    c11o->Print("output/uncertainities/eps/pile_up_eta_star_a.eps");
    c11o->Close();

    TCanvas *c11b = new TCanvas("c11b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    eta_star_a->SetMaximum(1.5);
    eta_star_a->SetMinimum(0.5);
    eta_star_a->SetLineColor(2);
    eta_star_a->SetLineStyle(1);
    eta_star_a->SetLineWidth(3);
    eta_star_a->Draw("e1");
    eta_star_b->SetLineColor(3);
    eta_star_b->SetLineStyle(1);
    eta_star_b->SetLineWidth(3);
    eta_star_b->Draw("e1 same");
    eta_star_c->SetLineColor(4);
    eta_star_c->SetLineStyle(1);
    eta_star_c->SetLineWidth(3);
    eta_star_c->Draw("e1 same");

    TLegend *leg11b = new TLegend(0.5,0.84,0.97,0.97);
    leg11b->AddEntry(eta_star_a,labela,"l");
    leg11b->AddEntry(eta_star_b,labelb,"l");
    leg11b->AddEntry(eta_star_c,labelc,"l");
    leg11b->SetFillColor(0);
    leg11b->Draw();
    
    c11b->Print("output/uncertainities/png/pile_up_eta_star_b.png");
    c11b->Print("output/uncertainities/c/pile_up_eta_star_b.C");
    c11b->Print("output/uncertainities/eps/pile_up_eta_star_b.eps");
    c11b->Close();


    TCanvas *c12 = new TCanvas("c12","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
     TH1D *x7a2 = (TH1D*) data_file1a->Get("ak5PF_delta_eta_outside_gap");
     TH1D *x8a2 = (TH1D*) data_file2a->Get("ak5PF_delta_eta_outside_gap");
     TH1D *x7b2 = (TH1D*) data_file1b->Get("ak5PF_delta_eta_outside_gap");
     TH1D *x8b2 = (TH1D*) data_file2b->Get("ak5PF_delta_eta_outside_gap");
     TH1D *x7c2 = (TH1D*) data_file1c->Get("ak5PF_delta_eta_outside_gap");
     TH1D *x8c2 = (TH1D*) data_file2c->Get("ak5PF_delta_eta_outside_gap");
     
     enta = x7a2->GetEntries()-x7a2->GetNbinsX();
     entb = x7b2->GetEntries()-x7b2->GetNbinsX();
     entc = x7c2->GetEntries()-x7c2->GetNbinsX();

     int1a = x7a2->Integral();
     int2a = x8a2->Integral();
     factora = int2a/int1a;

     int1b = x7b2->Integral();
     int2b = x8b2->Integral();
     factorb = int2b/int1b;

     int1c = x7c2->Integral();
     int2c = x8c2->Integral();
     factorc = int2c/int1c;
     
     //cout<<"factor"<<factor<<endl;
    
    TH1D *delta_eta_out, *delta_eta_out2, *delta_eta_out_old, *delta_eta_out_a, *delta_eta_out_b, *delta_eta_out_c;
    delta_eta_out =  new TH1D("Delta_eta_out","delta_eta_out;#Delta#eta;Pile Up Uncertainty", deta_out_nbins, deta_out_bins);
    delta_eta_out2 =  new TH1D("Delta_eta_out2","delta_eta_out2;#Delta#eta;Pile Up Uncertainty", deta_out_nbins, deta_out_bins);
    delta_eta_out_old =  new TH1D("Delta_eta_out_old","delta_eta_out_old;#Delta#eta;Pile Up Uncertainty", deta_out_nbins, deta_out_bins);
    delta_eta_out_a =  new TH1D("Delta_eta_out_a","delta_eta_out_a;#Delta#eta;Ratio to all vertexes selection", deta_out_nbins, deta_out_bins);
    delta_eta_out_b =  new TH1D("Delta_eta_out_b","delta_eta_out_b;#Delta#eta;Ratio to all vertexes selection", deta_out_nbins, deta_out_bins);
    delta_eta_out_c =  new TH1D("Delta_eta_out_c","delta_eta_out_c;#Delta#eta;Ratio to all vertexes selection", deta_out_nbins, deta_out_bins);
    delta_eta_out_a->Divide(x7a2,x8a2,1.,1./factora,"");
    delta_eta_out_b->Divide(x7b2,x8b2,1.,1./factorb,"");
    delta_eta_out_c->Divide(x7c2,x8c2,1.,1./factorc,"");
    
    min = 0.0;
    max = 0.0;
    tot = 0.0;
    tot_old = 0.0;
    
    for(Int_t i=1;i<=delta_eta_out_a->GetNbinsX();i++)
    {
    Float_t conta = delta_eta_out_a->GetBinContent(i) - 1;
    Float_t contb = delta_eta_out_b->GetBinContent(i) - 1;
    Float_t contc = delta_eta_out_c->GetBinContent(i) - 1;
    abs_cont = (enta*conta+entb*contb+entc*contc)/(enta+entb+entc);
    if (abs_cont < 0) { abs_cont = -abs_cont;}
    tot = tot + abs_cont;
    delta_eta_out2->SetBinContent(i,abs_cont);
    if (conta < 0) { conta = -conta;}
    if (conta > max) {max = conta;}
    if (conta < min || i == 1) { min = conta;}
    if (contb < 0) { contb = -contb;}
    if (contb > max) {max = contb;}
    if (contb < min) { min = contb;}
    if (contc < 0) { contc = -contc;}
    if (contc > max) {max = contc;}
    if (contc < min) { min = contc;}
    //cout<<i<<" ="<<conta<<" . "<<contb<<" . "<<contc<<endl;
    delta_eta_out_old->SetBinContent(i,conta);
    tot_old = tot_old + conta;
    }
    
    alisar(delta_eta_out2, delta_eta_out, step);

    delta_eta_out->SetLineColor(2);
    delta_eta_out->SetLineStyle(1);
    delta_eta_out->SetLineWidth(3);
    delta_eta_out->Draw("hist");

    TLegend *leg12 = new TLegend(0.13,0.9,0.55,0.97);
    leg12->AddEntry(delta_eta_out,"Average Pile-up uncertainty for 2010 Data","l");
    leg12->SetFillColor(0);
    leg12->Draw();
    
    ave = tot/delta_eta_out->GetNbinsX();
    ave_old = tot_old/delta_eta_out->GetNbinsX();
    cout<<"Delta eta outside gap : min ="<<min*100<<", max = "<<max*100<<", ave = "<<ave*100<<", ave old = "<<ave_old*100<<endl;
    min_row[9] = min*100;
    max_row[9] = max*100;
    ave_row[9] = ave*100;
    old_row[9] = ave_old*100;
    
    c12->Print("output/uncertainities/png/pile_up_delta_eta_out.png");
    c12->Print("output/uncertainities/c/pile_up_delta_eta_out.C");
    c12->Print("output/uncertainities/eps/pile_up_delta_eta_out.eps");
    c12->Close();


    TCanvas *c12o = new TCanvas("c12o","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    delta_eta_out->SetMinimum(0.0);
    delta_eta_out->SetMaximum(0.2);
    delta_eta_out->Draw("hist");
    delta_eta_out2->SetLineColor(2);
    delta_eta_out2->SetLineStyle(3);
    delta_eta_out2->SetLineWidth(3);
    delta_eta_out2->Draw("hist same");
    delta_eta_out_old->SetLineColor(2);
    delta_eta_out_old->SetLineStyle(2);
    delta_eta_out_old->SetLineWidth(3);
    delta_eta_out_old->Draw("hist same");

    TLegend *leg12o = new TLegend(0.5,0.84,0.97,0.97);
    leg12o->AddEntry(delta_eta_out,"Average Pile-up uncertainty for 2010 Data","l");
    leg12o->AddEntry(delta_eta_out2,"Average Pile-up uncertainty for 2010 Data (no smoothing)","l");
    leg12o->AddEntry(delta_eta_out_old,"Average Pile-up uncertainty for JetMETTau_2010A dataset","l");
    leg12o->SetFillColor(0);
    leg12o->Draw();

    c12o->Print("output/uncertainities/png/pile_up_delta_eta_out_a.png");
    c12o->Print("output/uncertainities/c/pile_up_delta_eta_out_a.C");
    c12o->Print("output/uncertainities/eps/pile_up_delta_eta_out_a.eps");
    c12o->Close();

    TCanvas *c12b = new TCanvas("c12b","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);
    
    delta_eta_out_a->SetMaximum(1.5);
    delta_eta_out_a->SetMinimum(0.5);
    delta_eta_out_a->SetLineColor(2);
    delta_eta_out_a->SetLineStyle(1);
    delta_eta_out_a->SetLineWidth(3);
    delta_eta_out_a->Draw("e1");
    delta_eta_out_b->SetLineColor(3);
    delta_eta_out_b->SetLineStyle(1);
    delta_eta_out_b->SetLineWidth(3);
    delta_eta_out_b->Draw("e1 same");
    delta_eta_out_c->SetLineColor(4);
    delta_eta_out_c->SetLineStyle(1);
    delta_eta_out_c->SetLineWidth(3);
    delta_eta_out_c->Draw("e1 same");

    TLegend *leg12b = new TLegend(0.13,0.84,0.45,0.97);
    leg12b->AddEntry(delta_eta_out_a,labela,"l");
    leg12b->AddEntry(delta_eta_out_b,labelb,"l");
    leg12b->AddEntry(delta_eta_out_c,labelc,"l");
    leg12b->SetFillColor(0);
    leg12b->Draw();
    
    c12b->Print("output/uncertainities/png/pile_up_delta_eta_out_b.png");
    c12b->Print("output/uncertainities/c/pile_up_delta_eta_out_b.C");
    c12b->Print("output/uncertainities/eps/pile_up_delta_eta_out_b.eps");
    c12b->Close();

    delta_phi->Write();
    delta_phi_deta1->Write();
    delta_phi_deta2->Write();
    delta_phi_deta3->Write();
    delta_phi_deta4->Write();
  //  delta_phi_deta5->Write();
  //  delta_phi_deta6->Write();
    delta_phi_gap->Write();
    delta_phi_deta1_gap->Write();
    delta_phi_deta2_gap->Write();
    delta_phi_deta3_gap->Write();
    delta_phi_deta4_gap->Write();
  //  delta_phi_deta5_gap->Write();
  //  delta_phi_deta6_gap->Write();
    delta_phi_nogap->Write();
    delta_phi_deta1_nogap->Write();
    delta_phi_deta2_nogap->Write();
    delta_phi_deta3_nogap->Write();
    delta_phi_deta4_nogap->Write();
  //  delta_phi_deta5_nogap->Write();
  //  delta_phi_deta6_nogap->Write();
    leading_pt_inside_gap->Write();
    leading_pt_outside_gap->Write();
    eta_star->Write();
    delta_eta_out->Write();


   cout<<"Delta Phi              : min = "<<min_row[0]<<",    max = "<<max_row[0]<<", ave = "<<ave_row[0]<<", ave old = "<<old_row[0]<<endl;
   cout<<"Delta Phi Deta         : min = "<<min_row[1]<<",   max = "<<max_row[1]<<", ave = "<<ave_row[1]<<", ave old = "<<old_row[1]<<endl;
   cout<<"Delta Phi Gap          : min = "<<min_row[2]<<",   max = "<<max_row[2]<<", ave = "<<ave_row[2]<<", ave old = "<<old_row[2]<<endl;
   cout<<"Delta Phi Gap Deta     : min = "<<min_row[3]<<",   max = "<<max_row[3]<<", ave = "<<ave_row[3]<<", ave old = "<<old_row[3]<<endl;
   cout<<"Delta Phi Nogap        : min = "<<min_row[4]<<",   max = "<<max_row[4]<<", ave = "<<ave_row[4]<<", ave old = "<<old_row[4]<<endl;
   cout<<"Delta Phi Nogap Deta   : min = "<<min_row[5]<<",   max = "<<max_row[5]<<",  ave = "<<ave_row[5]<<", ave old = "<<old_row[5]<<endl;
   cout<<"Leading pt inside gap  : min = "<<min_row[6]<<",  max = "<<max_row[6]<<", ave = "<<ave_row[6]<<", ave old = "<<old_row[6]<<endl;
   cout<<"Leading pt outside gap : min = "<<min_row[7]<<",   max = "<<max_row[7]<<", ave = "<<ave_row[7]<<", ave old = "<<old_row[7]<<endl;
   cout<<"Eta star inside gap    : min = "<<min_row[8]<<", max = "<<max_row[8]<<", ave = "<<ave_row[8]<<", ave old = "<<old_row[8]<<endl;
   cout<<"Delta eta outside gap  : min = "<<min_row[9]<<",   max = "<<max_row[9]<<", ave = "<<ave_row[9]<<", ave old = "<<old_row[9]<<endl;

}

void estimate_nvertices_pileup()
{
//estimates the relation between number of vertices and pileup
//this method has been used as a crosscheck
//not part of the main analysis, it was not refactored or adapted

gROOT->Reset();
gROOT->SetStyle("Plain");

   Int_t Vertex_multiplicity;
   std::vector<float> *num_PU_vertices;

    TH1D *vertex, *num_pu;
    TH2D *rel;

    vertex =  new TH1D("vertex","vertex;vertex;events", 10, 0, 10);
    num_pu =  new TH1D("num_PU","num_pu;num_pu;events", 10, 0, 10);
    rel =  new TH2D("rel","rel;offline;true;events", 10, 0, 10, 10, 0, 10);

TChain mc("demo/Event;1");
mc.Add("trees/private_mc/CFTree_v17_oldvertexreq.root");

TObjArray *fileElements = 0;
fileElements = mc.GetListOfFiles();
TIter next(fileElements);
TChainElement *chEl=0;
while (( chEl=(TChainElement*)next() )) {
   TFile f(chEl->GetTitle());

   TTree *data_Event = (TTree*)f.Get("demo/Event");

   data_Event->SetBranchAddress("Vertex_multiplicity",&Vertex_multiplicity);  
   data_Event->SetBranchAddress("num_PU_vertices",&num_PU_vertices); 

   cout<<"file:/"<<chEl->GetTitle()<<endl;
   
   Int_t nentries = (Int_t)data_Event->GetEntries();
   cout<<endl;
   //nentries = 10;

   for (Int_t i=0;i<nentries;i++) {
     data_Event->GetEntry(i);
     vertex->Fill(Vertex_multiplicity-1);
     num_pu->Fill(num_PU_vertices->at(0));
     rel->Fill(Vertex_multiplicity-1,num_PU_vertices->at(0));
}

}

    TCanvas *c1 = new TCanvas("c1","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);


    vertex->SetLineColor(3);
    vertex->SetLineStyle(2);
    vertex->SetLineWidth(3);
    vertex->Draw("hist");
    num_pu->SetLineColor(2);
    num_pu->SetLineStyle(1);
    num_pu->SetLineWidth(3);
    num_pu->Draw("hist same");

    TLegend *leg = new TLegend(0.6,0.8,0.96,0.95);
    leg->AddEntry(vertex,"offline num_PU_vertices","l");
    leg->AddEntry(num_pu,"True num_PU_vertices","l");
    leg->SetFillColor(0);
    leg->Draw();

    c1->Print("output/pile_up_studies/true_vs_offline_old.png");
    c1->Close();

    TCanvas *c2 = new TCanvas("c2","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    //gStyle->SetErrorX(0);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    rel->Draw("lego2z");

    c2->Print("output/pile_up_studies/true_vs_offline_2d_old.png");
    c2->Close();

}
