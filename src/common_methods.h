// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 22 Nov 2012
//
// provides functions to be used in the other parts of the analisys
//
// double get_non0_minimum(TH1D, bool detail = false)
// similar to TH1 GetMinimum but ignores the entries with value 0
//
// set_histogram_min_max(TH1D *histogram, double &min, double &max, bool detail = false)
// gets the max and minimum for and histogram
//
// double calc_delta_phi(double, double)
// computes delta phi
//
// plot_2histograms(TH1D *hist1, TString label1 = "label1", TH1D *hist2, TString label2 = "label2", string file_out = "test", string path= "../output/", string legend_position = "top_left", bool logscale = true, bool detail = false)
// default plotter for a pair of histograms

#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
#include <string>

using namespace std;

void set_legend_position(string legend_position, int size, double &x1, double &y1, double &x2, double &y2)
{
//sets the x1, y1, x2, y2 to draw the legend according to legend_positions and size

    if (legend_position == "top_left")
    {
    x1 = 0.13;
    y1 = 0.98;
    x2 = 0.61;
    if (size == 1) { y2 = 0.88; }
    if (size == 2) { y2 = 0.83; }
    if (size == 3) { y2 = 0.78; }
    if (size == 4) { y2 = 0.73; }
    if (size == 5) { y2 = 0.68; }
    if (size == 6) { y2 = 0.63; }
    }

    if (legend_position == "bottom_left")
    {
    x1 = 0.13;
    y1 = 0.20;
    x2 = 0.58;
    if (size == 1) { y2 = 0.30; }
    if (size == 2) { y2 = 0.35; }
    if (size == 3) { y2 = 0.40; }
    if (size == 4) { y2 = 0.45; }
    if (size == 5) { y2 = 0.50; }
    if (size == 6) { y2 = 0.55; }
    }
    
    if (legend_position == "bottom_middle")
    {
    x1 = 0.33;
    y1 = 0.18;
    x2 = 0.76;
    if (size == 1) { y2 = 0.28; }
    if (size == 2) { y2 = 0.33; }
    if (size == 3) { y2 = 0.38; }
    if (size == 4) { y2 = 0.43; }
    if (size == 5) { y2 = 0.43; }
    if (size == 6) { y2 = 0.43; }
    }
    
    if (legend_position == "top_right")
    {
    x1 = 0.50;
    y1 = 0.98;
    x2 = 0.98;
    if (size == 1) { y2 = 0.88; }
    if (size == 2) { y2 = 0.83; }
    if (size == 3) { y2 = 0.78; }
    if (size == 4) { y2 = 0.73; }
    if (size == 5) { y2 = 0.68; }
    if (size == 6) { y2 = 0.63; }
    }

    if (legend_position == "bottom_right")
    {
    x1 = 0.50;
    y1 = 0.18;
    x2 = 0.98;
    if (size == 1) { y2 = 0.28; }
    if (size == 2) { y2 = 0.33; }
    if (size == 3) { y2 = 0.38; }
    if (size == 4) { y2 = 0.43; }
    if (size == 5) { y2 = 0.48; }
    if (size == 6) { y2 = 0.53; }
    }
    
    if (legend_position == "middle_left")
    {
    x1 = 0.13;
    y1 = 0.48;
    x2 = 0.50;
    y2 = 0.68;
    }

    if (legend_position == "top_middle")
    {
    x1 = 0.33;
    y1 = 0.73;
    x2 = 0.76;
    y2 = 0.98;
    }
}


void print_plots(TCanvas *c, string path, string fileout)
{
//saves a given TCanvas c into path path with filename fileout

//setting the output files
   string out_png = path + "png/" + fileout + ".png";
   string out_c = path + "c/" + fileout + ".C";
   string out_eps = path + "eps/" + fileout + ".eps";
    
//save the file and close the canvas
    c->Print( out_png.c_str() );
    c->Print( out_c.c_str() );
    c->Print( out_eps.c_str() );
    c->Close();
}


void format_histogram(TH1 *hist, int color, int style)
{
//formats an histogram with color, style and line width

    hist->SetLineColor(color);
    hist->SetLineStyle(style);
    hist->SetLineWidth(4);
}


double CorrectFactorSmear(double eta)
{
  double corrfactor=0;
  if(fabs(eta)<1.1) corrfactor=1.066;
  if(fabs(eta)>=1.1 && fabs(eta)<1.7) corrfactor=1.191;
  if(fabs(eta)>=1.7 && fabs(eta)<2.3) corrfactor=1.096;
  if(fabs(eta)>=2.3 && fabs(eta)<5.0) corrfactor=1.166;

  return corrfactor;
}


double smearpt(double PFphi, double PFeta, double PFpt, double GenPhi, double GenEta, double GenPt, bool detail = false){
  
  double deltaPhiMatch1=-100;
  double deltaR1=100;
  double deltaRtest=100;
  double smear1=-100;
  double Pi=3.141592653;
  double newpt1=-100;
  
  deltaPhiMatch1=PFphi-GenPhi;
  if(deltaPhiMatch1<-Pi) deltaPhiMatch1=deltaPhiMatch1+2*Pi;
  if(deltaPhiMatch1>Pi) deltaPhiMatch1=deltaPhiMatch1-2*Pi;
  deltaPhiMatch1=fabs(deltaPhiMatch1);
  
  deltaRtest=sqrt(pow(deltaPhiMatch1,2)+pow(PFeta-GenEta,2));
      
  if(deltaRtest<deltaR1){
    newpt1=GenPt;
    deltaR1=deltaRtest;
  }
    
  if(deltaR1<0.3){
    smear1=newpt1+CorrectFactorSmear(PFeta)*(PFpt-newpt1);
  }   
  
  if(deltaR1>0.3) smear1=PFpt;
  
  if (detail) { cout << "Old Pt = " << GenPt << " New Pt = " << PFpt << endl; }
  return smear1;
}

double get_non0_minimum(TH1 *histogram, bool detail = false)
{
//similar to TH1 GetMinimum but ignores the entries with value 0
double min = 1e50;

    for(Int_t i=1;i<=histogram->GetNbinsX();i++)
    {
	if (min > histogram->GetBinContent(i) && histogram->GetBinContent(i) > 0.0) { min = histogram->GetBinContent(i); }
    }
    
    if (min == 0.0) { min = 1.0; }
    if (detail) { cout << "Non zero minimum required = " << min << endl; }    

return min;
}

void set_histogram_min_max(TH1 *histogram, double &min, double &max, bool detail = false)
{
//given a min and a max, this funtion will check if they fit on the histogram range, if not it will change them
    if (min > histogram->GetMinimum() && histogram->GetMinimum() > 0.0) { min = histogram->GetMinimum(); }
    if (min > get_non0_minimum(histogram) && histogram->GetMinimum() == 0.0) { min = get_non0_minimum(histogram,detail); }
    if (max < histogram->GetMaximum()) { max = histogram->GetMaximum(); }
}

double calc_delta_phi(double phi1, double phi2)
{
//computes delta phi
double delta_phi = 0.0;
double delta_phi_aux = 0.0;
const double pi = 3.14159265;

     delta_phi = phi1 - phi2;
     if (delta_phi < 0.0) { delta_phi = -delta_phi; }
     if (delta_phi > pi)
	{
	delta_phi_aux = delta_phi - pi;
	delta_phi = pi - delta_phi_aux;
	}

return delta_phi;
}


void normalize_histogram(TH1D *histogram, string name, bool detail = false)
{
//normalize a histogram

    double integral = 0.0;
    integral = histogram->Integral();
    if (detail) { cout<<name<<" = "<<integral<<" pb^-1"<<endl; }
    
    for(Int_t i=1;i<=histogram->GetNbinsX();i++)
    {
        Float_t cont = histogram->GetBinContent(i);
        Float_t width = histogram->GetBinWidth(i);
        Float_t error = histogram->GetBinError(i);
        cont /= width;
        error /= width;
        histogram->SetBinContent(i,cont);
        histogram->SetBinError(i,error);
    }
    
histogram->SetEntries(histogram->GetEntries() - histogram->GetNbinsX());

}


void plot_histogram(TH1D *histogram, string path, string fileout, TString label, string legend_position, bool isunc = false)
{
//plots the model uncertainty control plots

//declaring the canvas
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

//format and ploting the histogram
	if (isunc) { histogram->SetMinimum(0); histogram->SetMaximum(0.7); }
    format_histogram(histogram, 2, 1);
    histogram->Draw("e1");
    
//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg00 = new TLegend(x1, y1, x2, y2);
    leg00->AddEntry(histogram,label,"l");
    leg00->SetFillColor(0);
    leg00->SetLineWidth(1);
    leg00->SetLineColor(0);
    leg00->SetFillStyle(1001);
    leg00->Draw();
    
    print_plots(c1, path, fileout);
}


void plot_3histograms(TH1D *hist1 = 0, TString label1 = "label1", TH1D *hist2 = 0, TString label2 = "label2", TH1D *hist3 = 0, TString label3 = "label3", string path= "../output/", string fileout = "test", string legend_position = "top_left", bool logscale = true, bool detail = false)
{
// plots three distributions

// check if there are any histograms inputed
if (hist1 == 0) { cout << "Histogram1 is not provided!" << endl; return; } 
if (hist2 == 0) { cout << "Histogram2 is not provided!" << endl; return; } 
if (hist3 == 0) { cout << "Histogram3 is not provided!" << endl; return; } 

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
// declare and configure the canvas
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
    double min = 0.0;
    double max = hist1->GetMaximum();
    if (hist1->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(hist1,detail);
    }
    else
    {
    min = hist1->GetMinimum();
    }    
    
    set_histogram_min_max(hist2, min, max, detail);
    set_histogram_min_max(hist3, min, max, detail);
    
    max = 1.3 * max;
    min = 0.7 * min;

//plooting
    hist1->SetMaximum(max);
    hist1->SetMinimum(min);
    format_histogram(hist1, 1, 1);
    hist1->Draw("e1");
    format_histogram(hist2, 2, 2);
    hist2->Draw("e1 same");
    format_histogram(hist3, 4, 4);
    hist3->Draw("e1 same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 3, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hist1,label1,"l");
    leg01->AddEntry(hist2,label2,"l");
    leg01->AddEntry(hist3,label3,"l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}


void plot_4histograms(TH1D *hist1 = 0, TString label1 = "label1", TH1D *hist2 = 0, TString label2 = "label2", TH1D *hist3 = 0, TString label3 = "label3", TH1D *hist4 = 0, TString label4 = "label4", string path= "../output/", string fileout = "test", string legend_position = "top_left", bool logscale = true, bool detail = false)
{
// plots three distributions

// check if there are any histograms inputed
if (hist1 == 0) { cout << "Histogram1 is not provided!" << endl; return; } 
if (hist2 == 0) { cout << "Histogram2 is not provided!" << endl; return; } 
if (hist3 == 0) { cout << "Histogram3 is not provided!" << endl; return; }
if (hist4 == 0) { cout << "Histogram4 is not provided!" << endl; return; } 

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
// declare and configure the canvas
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
    double min = 0.0;
    double max = hist1->GetMaximum();
    if (hist1->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(hist1,detail);
    }
    else
    {
    min = hist1->GetMinimum();
    }    
    
    set_histogram_min_max(hist2, min, max, detail);
    set_histogram_min_max(hist3, min, max, detail);
    set_histogram_min_max(hist4, min, max, detail);
    
    max = 1.5 * max;
    min = 0.7 * min;

//plooting
    hist1->SetMaximum(max);
    hist1->SetMinimum(min);
    format_histogram(hist1, 1, 1);
    hist1->Draw("e1");
    format_histogram(hist2, 2, 2);
    hist2->Draw("e1 same");
    format_histogram(hist3, 4, 4);
    hist3->Draw("e1 same");
    format_histogram(hist4, 6, 6);
    hist4->Draw("e1 same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 4, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hist1,label1,"l");
    leg01->AddEntry(hist2,label2,"l");
    leg01->AddEntry(hist3,label3,"l");
    leg01->AddEntry(hist4,label4,"l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}


void plot_2histograms(TH1 *hist1 = 0, TString label1 = "label1", TH1 *hist2 = 0, TString label2 = "label2", string path= "../output/", string fileout = "test", string legend_position = "top_left", bool logscale = true, bool detail = false)
{
// plots the two distributions

// check if there are any histograms inputed
if (hist1 == 0) { cout << "Histogram1 is not provided!" << endl; return; } 
if (hist2 == 0) { cout << "Histogram2 is not provided!" << endl; return; } 

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
// declare and configure the canvas
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
    double min = 0.0;
    double max = hist1->GetMaximum();
    if (hist1->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(hist1,detail);
    }
    else
    {
    min = hist1->GetMinimum();
    }    
    
    set_histogram_min_max(hist2, min, max, detail);
    
    max = 1.3 * max;
    min = 0.7 * min;

//plooting
    hist1->SetMaximum(max);
    hist1->SetMinimum(min);
    format_histogram(hist1, 1, 1);
    hist1->Draw("e1");
    format_histogram(hist2, 2, 2);
    hist2->Draw("e1same");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 2, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hist1,label1,"l");
    leg01->AddEntry(hist2,label2,"l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}

void ratio_2histograms(TH1D *hist1 = 0, TH1D *hist2 = 0, TString label = "label", string path= "../output/", string fileout = "test", string legend_position = "top_left", bool detail = false)
{
// plots the two distributions

// check if there are any histograms inputed
if (hist1 == 0) { cout << "Histogram1 is not provided!" << endl; return; } 
if (hist2 == 0) { cout << "Histogram2 is not provided!" << endl; return; } 

TH1D *ratio;
ratio = static_cast<TH1D*>(hist1->Clone());
ratio->Divide(hist1,hist2,1.,1.,"B");

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
// declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

//calculate the plooting range
    //double min = 0.0;
    //double max = ratio->GetMaximum();
    //if (ratio->GetMinimum() == 0.0)
    //{
    //min = get_non0_minimum(ratio,detail);
    //}
    //else
    //{
    //min = ratio->GetMinimum();
    //}
    
    //max = 1.3 * max;
    //min = 0.7 * min;

//plooting
    //ratio->SetMaximum(max);
    //ratio->SetMinimum(min);
    format_histogram(ratio, 2, 1);
    ratio->Draw("e1");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(ratio,label,"l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, fileout);
}

void plot_efficiency(TH1D *hratio, string prefix, string fileout, string path= "../output/trigger_eff/", string legend_position = "top_left", bool detail = false)
{
//plots the correction factor

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    //gPad->SetLogy();
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    hratio->SetMinimum(0);
    hratio->SetMaximum(1.5);
    format_histogram(hratio, 4, 2);
    hratio->Draw("e1");

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hratio,"Trigger Efficiency","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, prefix + fileout);
}

void fit_and_plot(TH1D *hratio, double *fit, char *fit_function, double fit_begin, double fit_end, string prefix, string fileout, string path= "../output/trigger_eff/", string legend_position = "top_left", bool detail = false, bool test = false)
{
//plots the correction factor

    if (detail) { cout<<"Ploting "<<fileout<<endl; }
//declare and configure the canvas
    TCanvas *c01 = new TCanvas("c01","Canvas",0,29,1450,870);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(kFALSE);
    gPad->SetFillColor(0);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(2);
    //gPad->SetLogy();
    gPad->SetLeftMargin(0.10);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    gPad->SetFrameBorderMode(0);

    hratio->SetMinimum(0);
    hratio->SetMaximum(1.2);
    format_histogram(hratio, 4, 2);
    hratio->Draw("e1");
    hratio->Fit(fit_function, "", "", fit_begin, fit_end);
    TF1 *tfit = 0;
    tfit = hratio->GetFunction(fit_function);

	for (int i=0; i <=8; i++)
		{
   		fit[i] = tfit->GetParameter(i);
		fit[9+i] = tfit->GetParError(i);
		if (test) { cout << "Parameter = " << fit[i] << "+-" << fit[9+i] << endl; }
		}

//assign the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
    set_legend_position(legend_position, 1, x1, y1, x2, y2);

    TLegend *leg01 = new TLegend(x1, y1, x2, y2);
    leg01->AddEntry(hratio,"Trigger Efficiency","l");
    leg01->SetFillColor(0);
    leg01->SetLineWidth(1);
    leg01->SetLineColor(0);
    leg01->SetFillStyle(1001);
    leg01->Draw();

    print_plots(c01, path, prefix + fileout);
}


void plot_six_dist(TH1D *dist1, TString label1, TH1D *dist2, TString label2, TH1D *dist3, TString label3, TH1D *dist4, TString label4, TH1D *dist5, TString label5, TH1D *dist6, TString label6, string path, string prefix, string name, string legend_position = "top_left", bool detail = false)
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
    gPad->SetLogy();

//calculate the plooting range
    double min = 0.0;
    double max = dist1->GetMaximum();
    if (dist1->GetMinimum() == 0.0)
    {
    min = get_non0_minimum(dist1,detail);
    }
    else
    {
    min = dist1->GetMinimum();
    }    
    
    set_histogram_min_max(dist2, min, max, detail);
    set_histogram_min_max(dist3, min, max, detail);
    set_histogram_min_max(dist4, min, max, detail);
    set_histogram_min_max(dist5, min, max, detail);
    set_histogram_min_max(dist6, min, max, detail);
        
    max = 1.3 * max;
    min = 0.7 * min;

//plooting
    dist1->SetMaximum(max);
    dist1->SetMinimum(min);
    dist1->SetLineWidth(4);
    dist1->SetLineColor(2);
    dist1->SetLineStyle(1);
    dist1->Draw("e1");
    dist2->SetLineWidth(4);
    dist2->SetLineColor(3);
    dist2->SetLineStyle(2);
    dist2->Draw("e1same");
    dist3->SetLineWidth(4);
    dist3->SetLineColor(4);
    dist3->SetLineStyle(3);
    dist3->Draw("e1same");
    dist4->SetLineWidth(4);
    dist4->SetLineColor(8);
    dist4->SetLineStyle(4);
    dist4->Draw("e1same");
    dist5->SetLineWidth(4);
    dist5->SetLineColor(6);
    dist5->SetLineStyle(5);
    dist5->Draw("e1same");
    dist6->SetLineWidth(4);
    dist6->SetLineColor(7);
    dist6->SetLineStyle(6);
    dist6->Draw("e1same");

//sets and draw the legend
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;

    if (legend_position == "top_left")
    {
    x1 = 0.13;
    y1 = 0.98;
    x2 = 0.58;
    y2 = 0.68;
    }

    if (legend_position == "bottom_left")
    {
    x1 = 0.13;
    y1 = 0.13;
    x2 = 0.58;
    y2 = 0.23;
    }

    if (legend_position == "bottom_right")
    {
    x1 = 0.50;
    y1 = 0.13;
    x2 = 0.98;
    y2 = 0.23;
    }
    
    if (legend_position == "bottom_middle")
    {
    x1 = 0.33;
    y1 = 0.13;
    x2 = 0.76;
    y2 = 0.23;
    }

    TLegend *leg00 = new TLegend(x1,y1,x2,y2);
    leg00->AddEntry(dist1,label1,"l");
    leg00->AddEntry(dist2,label2,"l");
    leg00->AddEntry(dist3,label3,"l");
    leg00->AddEntry(dist4,label4,"l");
    leg00->AddEntry(dist5,label5,"l");
    leg00->AddEntry(dist6,label6,"l");
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
