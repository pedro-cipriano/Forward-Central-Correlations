// Pedro Cipriano, Jan 2012
// DESY, CMS
// Last Update: 05 Mar 2013
//
// compute_trigger_turn_on(string *data_in, string data_out, int n_files, string sel_mode = "allvertex", bool detail = false, bool test = false)
// Computes the trigger turn on curves using trigger elements

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TPad.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <string>

#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
#include "common_methods.h"

using namespace std;

double correction_factor(double *p, int index, double *range, double x, bool test = false)
{

double y = 0.0;
bool special_case = false;

if (x < range[index*2+0])
	{
	if (test) { cout << "X value is too low!" << endl; }
	y = 0.05;
	special_case = true;
	}
if (x > range[index*2+1])
	{
	if (test) { cout << "X value is too high!" << endl; }
	y = 1.0;
	special_case = true;
	}

if (!special_case) { y = p[index*9+8]*pow(x,8) + p[index*9+7]*pow(x,7) + p[index*9+6]*pow(x,6) + p[index*9+5]*pow(x,5) + p[index*9+4]*pow(x,4) + p[index*9+3]*pow(x,3) + p[index*9+2]*pow(x,2) + p[index*9+1]*x + p[index*9+0]; }

//if (test) { cout << "Correction factor = " << y << endl; }
return y;
}

void check_correction_function(double *p, int index, double *range, bool test = false)
{

double y = 0.0, x = 0.0, last_value = 0.0;
double npoints = 1000;
bool endloop = false;
double cutoff = 0.0;


if (test) { cout << "Index = " << index << endl; }

if (index == 0) { cutoff = 50.0; }
if (index == 1) { cutoff = 80.0; }
double new_min = range[index*2+0];
double new_max = range[index*2+1];

if (test) { cout << "Range: Minimum = " << range[index*2+0] << " Maximum = " << range[index*2+1] << endl; }

for (int i=0; i <= npoints; i++)
	{
	if (!endloop)
		{
		x = range[index*2+0] + i*(range[index*2+1] - range[index*2+0])/(double)npoints;
		y = correction_factor(p, index, range, x, test);
		//if (test) { cout << "x = " << x << " -> " << y << endl; }
		if (y < 0.02 or last_value < 0.0)
			{
			if (test) { cout << "New minimum found! x = " << x << " -> " << y << endl; }
			new_min = x;
			}
		if ((y > 1.01 or y < last_value) && x > cutoff)
			{
			if (test) { cout << "New maximum found! x = " << x << " -> " << y << endl; }
			new_max = x;
			endloop = true;
			}
		last_value = y;
		}
	}

range[index*2+0] = new_min;
range[index*2+1] = new_max;

cout << "New Range: Minimum = " << range[index*2+0] << " Maximum = " << range[index*2+1] << endl;
}

void correction_factor(TH1* hist1, TH1* hist2, double x, double *factor)
{
factor[0] = 1.0;
factor[1] = 1.0;
double bin_min, bin_max, bin_center, bin_width, bin_selected = 0, bin_content1, bin_content2;


for (int i=1; i <= hist1->GetNbinsX(); i++)
{
bin_center = hist1->GetBinCenter(i);
bin_width = hist1->GetBinWidth(i);
bin_min = bin_center - bin_width;
bin_max = bin_center + bin_width;
if (x > bin_min && x < bin_max) { bin_selected = i; }
}

bin_content1 = hist1->GetBinContent(bin_selected);
bin_content2 = hist2->GetBinContent(bin_selected);

if (bin_content1 > 0)
{
factor[0] = 1.0/bin_content1;
}

if (bin_content2 > 0)
{
factor[1] = 1.0/bin_content2;
}

}


void plot_final_trigger_efficiency(string data_jetmettau, string data_jetmet, string data_jet, string path, bool detail = false, bool test = false)
{

//output the configuration
   if (detail) { cout<<"Plot Final Trigger Efficiency Configuration"<<endl; }
   if (detail) { cout<<"Data JetMETTau :  "<<data_jetmettau<<endl; }
   if (detail) { cout<<"Data JetMET :     "<<data_jetmet<<endl; }
   if (detail) { cout<<"Data Jet :        "<<data_jet<<endl; }
   if (detail) { cout<<"Output path :     "<<path<<endl; }
   if (detail) { cout<<"Detail level :    "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :       "<<test<<endl; }

   //opening files
   TFile *jetmettau_file = new TFile( data_jetmettau.c_str() );
   TFile *jetmet_file = new TFile( data_jetmet.c_str() );
   TFile *jet_file = new TFile( data_jet.c_str() );

   //declaring empty pointers
   TH1D *jetmettau_15U = 0;
   TH1D *jetmettau_30U = 0;
   TH1D *jetmettau_50U = 0;
   TH1D *jetmet_15U = 0;
   TH1D *jetmet_30U = 0;
   TH1D *jetmet_50U = 0;
   TH1D *jet_15U = 0;
   TH1D *jet_30U = 0;
   TH1D *jet_50U = 0;

   //switch to check if all histograms were loaded sucessfully
   bool not_loaded = false;

   //loading histograms
   jetmettau_file->GetObject("ak5PF_leading_pt_HLT15U_eff_fine",jetmettau_15U);
   if (jetmettau_15U == 0) { cout << "ak5PF_leading_pt_HLT15U_eff_fine on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("ak5PF_leading_pt_HLT30U_eff_fine",jetmettau_30U);
   if (jetmettau_30U == 0) { cout << "ak5PF_leading_pt_HLT30U_eff_fine on JetMETTau_2010A not found!" << endl; not_loaded = true; }
   jetmettau_file->GetObject("ak5PF_leading_pt_HLT50U_eff_fine",jetmettau_50U);
   if (jetmettau_50U == 0) { cout << "ak5PF_leading_pt_HLT50U_eff_fine on JetMETTau_2010A not found!" << endl; not_loaded = true; }

   jetmet_file->GetObject("ak5PF_leading_pt_HLT15U_eff_fine",jetmet_15U);
   if (jetmet_15U == 0) { cout << "ak5PF_leading_pt_HLT15U_eff_fine on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("ak5PF_leading_pt_HLT30U_eff_fine",jetmet_30U);
   if (jetmet_30U == 0) { cout << "ak5PF_leading_pt_HLT30U_eff_fine on JetMET_2010A not found!" << endl; not_loaded = true; }
   jetmet_file->GetObject("ak5PF_leading_pt_HLT50U_eff_fine",jetmet_50U);
   if (jetmet_50U == 0) { cout << "ak5PF_leading_pt_HLT50U_eff_fine on JetMET_2010A not found!" << endl; not_loaded = true; }

   jet_file->GetObject("ak5PF_leading_pt_HLT15U_eff_fine",jet_15U);
   if (jet_15U == 0) { cout << "ak5PF_leading_pt_HLT15U_eff_fine on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("ak5PF_leading_pt_HLT30U_eff_fine",jet_30U);
   if (jet_30U == 0) { cout << "ak5PF_leading_pt_HLT30U_eff_fine on Jet_2010B not found!" << endl; not_loaded = true; }
   jet_file->GetObject("ak5PF_leading_pt_HLT50U_eff_fine",jet_50U);
   if (jet_50U == 0) { cout << "ak5PF_leading_pt_HLT50U_eff_fine on Jet_2010B not found!" << endl; not_loaded = true; }

   //Summary regarding the loading of the histograms
   if (detail && !not_loaded) { cout<<"All histograms loaded sucessfully!"<<endl; }
   if (not_loaded)
	{
	cout<<"Some histograms were not loaded! Please check the input rootfile. The routine will terminate now!" << endl;
	return;
	}

	//for (int i = 1; i <= jetmettau_15U->GetBinsX(); i++)
	//{

	//}

   //plotting
   plot_3histograms(jetmettau_15U, "HLT_Jet15U", jetmettau_30U, "HLT_Jet30U", jetmettau_50U, "HLT_Jet50U", path, "trigger_efficiency_final_jetmettau", "bottom_right", false, detail);
   plot_3histograms(jetmet_15U, "HLT_Jet15U", jetmet_30U, "HLT_Jet30U", jetmet_50U, "HLT_Jet50U", path, "trigger_efficiency_final_jetmet", "bottom_right", false, detail);
   plot_3histograms(jet_15U, "HLT_Jet15U", jet_30U, "HLT_Jet30U", jet_50U, "HLT_Jet50U", path, "trigger_efficiency_final_jet", "bottom_right", false, detail);

}

void compute_trigger_turn_on(string *data_in, string data_out, int n_files, string correction, string sel_mode = "allvertex", bool detail = false, bool test = false)
{

   double pt_min = 35.0;
   int nJetTrig = 4;
   int HLTJetPtN[3] = {15,30,50};
   int ATrig[3] = {6,20,30};
   bool corr = false;
   //double p[2*9];
   //double range[2*2];
   TH1D *trigger_hist1 = 0;
   TH1D *trigger_hist2 = 0;

//output the configuration
   if (detail) { cout<<"Get Triggered Events Configuration"<<endl; }
   if (detail) { cout<<"Selection mode :  "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files : "<<n_files<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :     "<<data_out<<endl; }
   if (detail) { cout<<"Correction File : "<<correction<<endl; }
   if (detail) { cout<<"Pt Min :          "<<pt_min<<endl; }
   if (detail) { cout<<"Detail level :    "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :       "<<test<<endl; }

//loading correction file
/*
   if (correction != "")
	{
	if (detail) { cout << "Loading correction file..." << endl; }
  	FILE *f = fopen( correction.c_str() , "r");
  	if(f==NULL)
		{
    		cout << "Can't find file: " << correction << endl;
    		return;
  		}
	if (detail) { cout << "Correction file loaded sucessfully!" << endl; }
	if (detail) { cout << "Loading fitting..." << endl; }
	int i = 0;
	float par30, par50, err30, err50;
	while(!feof(f))
		{
		fscanf(f,"%f %f %f %f", &par30, &err30, &par50, &err50);
		if (i < 9)
			{
			p[i] = par30;
			p[9+i] = par50;
			}
		i = i + 1;
		}
	range[0] = 36.0;
	range[1] = 80.0;
	range[2] = 36.0;
	range[3] = 120.0;
	if (detail) { cout << "HLT_Jet30U : y = " << p[8] << "*x^8 + " << p[7] << "*x^7 + " << p[6] << "*x^6 + " << p[5] << "*x^5 +" << endl; }
	if (detail) { cout << "                 " << p[4] << "*x^4 + " << p[3] << "*x^3 + " << p[2] << "*x^2 + " << p[1] << "*x + " << p[0] << endl; }
	if (detail) { cout << "HLT_Jet50U : y = " << p[17] << "*x^8 + " << p[16] << "*x^7 + " << p[15] << "*x^6 + " << p[14] << "*x^5 +" << endl; }
	if (detail) { cout << "                 " << p[13] << "*x^4 + " << p[12] << "*x^3 + " << p[11] << "*x^2 + " << p[10] << "*x + " << p[9] << endl; }
	if (detail) { cout << "Testing fitting..." << endl; }
	check_correction_function(p, 0, range, test);
	check_correction_function(p, 1, range, test);
	corr = true;
	}
*/

   if (correction != "")
	{
	TFile *trigger_file = new TFile( correction.c_str() );
   	trigger_file->GetObject("ak5PF_leading_central_pt_HLT30U_eff_fine",trigger_hist1);
	trigger_file->GetObject("ak5PF_leading_central_pt_HLT50U_eff_fine",trigger_hist2);
   	if (trigger_hist1 != 0 && trigger_hist2 != 0) { corr = true; }
	else { cout << "Correction histograms not found! The routine will terminate now!" << endl; return;}
   	if (detail) { cout << "Correction factors loaded sucessfully!" << endl; }
	}

//binning
   int all_nbins = 7;
   double all_bins[8] = {35, 45, 57, 72, 90, 120, 150, 200};

   int deta_nbins = 4;
   double deta_bins[5] = {0.4, 2.5, 3.5, 4.5, 7.5};

   int dphi_nbins = 7;
   double dphi_bins[8] = {0, 0.45, 0.9, 1.35, 1.8, 2.25, 2.7, 3.15};

//declaring the variables
   TString HLTJet[nJetTrig];

   int ihltj[nJetTrig];
   int prescalej[nJetTrig];
   int nentries = 0;
   int counter_entries = 0, counter_pv = 0, counter_hlt = 0, counter_selected = 0, counter_jet = 0, counter_l1obj = 0, counter_hltobj = 0;
   int triggered[4], l1pass[4], hltpass[4], condpass[4], selected[4], final[4];   

   char trigtitle[200];
   
   bool pv_pass = false;
   bool hltPass = false;
   bool l1cut[nJetTrig], hltcut[nJetTrig];
   bool hltPassj[nJetTrig];
   bool fill_all, fill_emu, emu_l1, emu_hlt;
   bool hard_emission;
   
   double pt, eta, phi;
   double leading_pt, leading_eta, leading_phi;
   double forward_pt, forward_eta, forward_phi;
   double central_pt, central_eta, central_phi;
   double hard_pt, hard_eta, hard_phi;
   double delta_eta, delta_phi;
   double pu_scale, eff;
   double trigger_factor, trigger_factors[2];

//declaring histograms
//general control event count
     TH1D *hist_events;
     TH1D *hist_events_HLT_Jet15U;
     TH1D *hist_events_HLT_Jet30U;
     TH1D *hist_events_HLT_Jet50U;
     TH1D *hist_events_all;

//monitor distributions with broad binning for dijets
     TH1D *hist_leading_pt_HLT_Jet15U_all;
     TH1D *hist_leading_pt_HLT_Jet30U_all;
     TH1D *hist_leading_pt_HLT_Jet50U_all;
     TH1D *hist_leading_pt_all;

     TH1D *hist_leading_central_pt_HLT_Jet15U_all;
     TH1D *hist_leading_central_pt_HLT_Jet30U_all;
     TH1D *hist_leading_central_pt_HLT_Jet50U_all;
     TH1D *hist_leading_central_pt_all;

     TH1D *hist_leading_eta_HLT_Jet15U_all;
     TH1D *hist_leading_eta_HLT_Jet30U_all;
     TH1D *hist_leading_eta_HLT_Jet50U_all;
     TH1D *hist_leading_eta_all;
     TH1D *hist_leading_eta_HLT_Jet30U_all_lowpt;
     TH1D *hist_leading_eta_HLT_Jet50U_all_lowpt;
     TH1D *hist_leading_eta_all_lowpt;
     TH1D *hist_leading_eta_HLT_Jet30U_all_highpt;
     TH1D *hist_leading_eta_HLT_Jet50U_all_highpt;
     TH1D *hist_leading_eta_all_highpt;

     TH1D *hist_leading_phi_HLT_Jet15U_all;
     TH1D *hist_leading_phi_HLT_Jet30U_all;
     TH1D *hist_leading_phi_HLT_Jet50U_all;
     TH1D *hist_leading_phi_all;

//monitor distributions with broad binning for leading jet
     TH1D *hist_sj_leading_pt_HLT_Jet15U_all;
     TH1D *hist_sj_leading_pt_HLT_Jet30U_all;
     TH1D *hist_sj_leading_pt_HLT_Jet50U_all;
     TH1D *hist_sj_leading_pt_all;

     TH1D *hist_sj_leading_eta_HLT_Jet15U_all;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_all;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_all;
     TH1D *hist_sj_leading_eta_all;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_all_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_all_lowpt;
     TH1D *hist_sj_leading_eta_all_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_all_highpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_all_highpt;
     TH1D *hist_sj_leading_eta_all_highpt;

     TH1D *hist_sj_leading_phi_HLT_Jet15U_all;
     TH1D *hist_sj_leading_phi_HLT_Jet30U_all;
     TH1D *hist_sj_leading_phi_HLT_Jet50U_all;
     TH1D *hist_sj_leading_phi_all;
     
//triggered distributions with broad binning for dijets
     TH1D *hist_leading_pt_HLT_Jet15U_emulated;
     TH1D *hist_leading_pt_HLT_Jet30U_emulated;
     TH1D *hist_leading_pt_HLT_Jet50U_emulated;
     TH1D *hist_leading_pt_emulated;

     TH1D *hist_leading_central_pt_HLT_Jet15U_emulated;
     TH1D *hist_leading_central_pt_HLT_Jet30U_emulated;
     TH1D *hist_leading_central_pt_HLT_Jet50U_emulated;
     TH1D *hist_leading_central_pt_emulated;

     TH1D *hist_leading_eta_HLT_Jet15U_emulated;
     TH1D *hist_leading_eta_HLT_Jet30U_emulated;
     TH1D *hist_leading_eta_HLT_Jet50U_emulated;
     TH1D *hist_leading_eta_emulated;
     TH1D *hist_leading_eta_HLT_Jet30U_emulated_lowpt;
     TH1D *hist_leading_eta_HLT_Jet50U_emulated_lowpt;
     TH1D *hist_leading_eta_emulated_lowpt;
     TH1D *hist_leading_eta_HLT_Jet30U_emulated_highpt;
     TH1D *hist_leading_eta_HLT_Jet50U_emulated_highpt;
     TH1D *hist_leading_eta_emulated_highpt;

     TH1D *hist_leading_phi_HLT_Jet15U_emulated;
     TH1D *hist_leading_phi_HLT_Jet30U_emulated;
     TH1D *hist_leading_phi_HLT_Jet50U_emulated;
     TH1D *hist_leading_phi_emulated;

//triggered distributions with broad binning for leading jet
     TH1D *hist_sj_leading_pt_HLT_Jet15U_emulated;
     TH1D *hist_sj_leading_pt_HLT_Jet30U_emulated;
     TH1D *hist_sj_leading_pt_HLT_Jet50U_emulated;
     TH1D *hist_sj_leading_pt_emulated;

     TH1D *hist_sj_leading_eta_HLT_Jet15U_emulated;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_emulated;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_emulated;
     TH1D *hist_sj_leading_eta_emulated;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt;
     TH1D *hist_sj_leading_eta_emulated_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_emulated_highpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_emulated_highpt;
     TH1D *hist_sj_leading_eta_emulated_highpt;

     TH1D *hist_sj_leading_phi_HLT_Jet15U_emulated;
     TH1D *hist_sj_leading_phi_HLT_Jet30U_emulated;
     TH1D *hist_sj_leading_phi_HLT_Jet50U_emulated;
     TH1D *hist_sj_leading_phi_emulated;

//efficiency distributions with broad binning for dijets
     TH1D *hist_leading_pt_HLT_Jet15U_eff;
     TH1D *hist_leading_pt_HLT_Jet30U_eff;
     TH1D *hist_leading_pt_HLT_Jet50U_eff;
     TH1D *hist_leading_pt_eff;

     TH1D *hist_leading_central_pt_HLT_Jet15U_eff;
     TH1D *hist_leading_central_pt_HLT_Jet30U_eff;
     TH1D *hist_leading_central_pt_HLT_Jet50U_eff;
     TH1D *hist_leading_central_pt_eff;

     TH1D *hist_leading_eta_HLT_Jet15U_eff;
     TH1D *hist_leading_eta_HLT_Jet30U_eff;
     TH1D *hist_leading_eta_HLT_Jet50U_eff;
     TH1D *hist_leading_eta_eff;
     TH1D *hist_leading_eta_HLT_Jet30U_eff_lowpt;
     TH1D *hist_leading_eta_HLT_Jet50U_eff_lowpt;
     TH1D *hist_leading_eta_eff_lowpt;
     TH1D *hist_leading_eta_HLT_Jet30U_eff_highpt;
     TH1D *hist_leading_eta_HLT_Jet50U_eff_highpt;
     TH1D *hist_leading_eta_eff_highpt;

     TH1D *hist_leading_phi_HLT_Jet15U_eff;
     TH1D *hist_leading_phi_HLT_Jet30U_eff;
     TH1D *hist_leading_phi_HLT_Jet50U_eff;
     TH1D *hist_leading_phi_eff;

//efficiency distributions with broad binning for leading jet
     TH1D *hist_sj_leading_pt_HLT_Jet15U_eff;
     TH1D *hist_sj_leading_pt_HLT_Jet30U_eff;
     TH1D *hist_sj_leading_pt_HLT_Jet50U_eff;
     TH1D *hist_sj_leading_pt_eff;

     TH1D *hist_sj_leading_eta_HLT_Jet15U_eff;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_eff;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_eff;
     TH1D *hist_sj_leading_eta_eff;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_eff_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_eff_lowpt;
     TH1D *hist_sj_leading_eta_eff_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_eff_highpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_eff_highpt;
     TH1D *hist_sj_leading_eta_eff_highpt;

     TH1D *hist_sj_leading_phi_HLT_Jet15U_eff;
     TH1D *hist_sj_leading_phi_HLT_Jet30U_eff;
     TH1D *hist_sj_leading_phi_HLT_Jet50U_eff;
     TH1D *hist_sj_leading_phi_eff;

//monitor distributions with fine binning for dijets
     TH1D *hist_leading_pt_HLT_Jet15U_all_fine;
     TH1D *hist_leading_pt_HLT_Jet30U_all_fine;
     TH1D *hist_leading_pt_HLT_Jet50U_all_fine;
     TH1D *hist_leading_pt_all_fine;

     TH1D *hist_leading_central_pt_HLT_Jet15U_all_fine;
     TH1D *hist_leading_central_pt_HLT_Jet30U_all_fine;
     TH1D *hist_leading_central_pt_HLT_Jet50U_all_fine;
     TH1D *hist_leading_central_pt_all_fine;

     TH1D *hist_leading_eta_HLT_Jet15U_all_fine;
     TH1D *hist_leading_eta_HLT_Jet30U_all_fine;
     TH1D *hist_leading_eta_HLT_Jet50U_all_fine;
     TH1D *hist_leading_eta_all_fine;
     TH1D *hist_leading_eta_HLT_Jet30U_all_fine_lowpt;
     TH1D *hist_leading_eta_HLT_Jet50U_all_fine_lowpt;
     TH1D *hist_leading_eta_all_fine_lowpt;
     TH1D *hist_leading_eta_HLT_Jet30U_all_fine_highpt;
     TH1D *hist_leading_eta_HLT_Jet50U_all_fine_highpt;
     TH1D *hist_leading_eta_all_fine_highpt;

     TH1D *hist_leading_phi_HLT_Jet15U_all_fine;
     TH1D *hist_leading_phi_HLT_Jet30U_all_fine;
     TH1D *hist_leading_phi_HLT_Jet50U_all_fine;
     TH1D *hist_leading_phi_all_fine;

//monitor distributions with fine binning for leading jet
     TH1D *hist_sj_leading_pt_HLT_Jet15U_all_fine;
     TH1D *hist_sj_leading_pt_HLT_Jet30U_all_fine;
     TH1D *hist_sj_leading_pt_HLT_Jet50U_all_fine;
     TH1D *hist_sj_leading_pt_all_fine;

     TH1D *hist_sj_leading_eta_HLT_Jet15U_all_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_all_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_all_fine;
     TH1D *hist_sj_leading_eta_all_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt;
     TH1D *hist_sj_leading_eta_all_fine_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt;
     TH1D *hist_sj_leading_eta_all_fine_highpt;

     TH1D *hist_sj_leading_phi_HLT_Jet15U_all_fine;
     TH1D *hist_sj_leading_phi_HLT_Jet30U_all_fine;
     TH1D *hist_sj_leading_phi_HLT_Jet50U_all_fine;
     TH1D *hist_sj_leading_phi_all_fine;

//triggered distributions with fine binning for dijets
     TH1D *hist_leading_pt_HLT_Jet15U_emulated_fine;
     TH1D *hist_leading_pt_HLT_Jet30U_emulated_fine;
     TH1D *hist_leading_pt_HLT_Jet50U_emulated_fine;
     TH1D *hist_leading_pt_emulated_fine;

     TH1D *hist_leading_central_pt_HLT_Jet15U_emulated_fine;
     TH1D *hist_leading_central_pt_HLT_Jet30U_emulated_fine;
     TH1D *hist_leading_central_pt_HLT_Jet50U_emulated_fine;
     TH1D *hist_leading_central_pt_emulated_fine;

     TH1D *hist_leading_eta_HLT_Jet15U_emulated_fine;
     TH1D *hist_leading_eta_HLT_Jet30U_emulated_fine;
     TH1D *hist_leading_eta_HLT_Jet50U_emulated_fine;
     TH1D *hist_leading_eta_emulated_fine;
     TH1D *hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt;
     TH1D *hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt;
     TH1D *hist_leading_eta_emulated_fine_lowpt;
     TH1D *hist_leading_eta_HLT_Jet30U_emulated_fine_highpt;
     TH1D *hist_leading_eta_HLT_Jet50U_emulated_fine_highpt;
     TH1D *hist_leading_eta_emulated_fine_highpt;

     TH1D *hist_leading_phi_HLT_Jet15U_emulated_fine;
     TH1D *hist_leading_phi_HLT_Jet30U_emulated_fine;
     TH1D *hist_leading_phi_HLT_Jet50U_emulated_fine;
     TH1D *hist_leading_phi_emulated_fine;

//triggered distributions with fine binning for leading jet
     TH1D *hist_sj_leading_pt_HLT_Jet15U_emulated_fine;
     TH1D *hist_sj_leading_pt_HLT_Jet30U_emulated_fine;
     TH1D *hist_sj_leading_pt_HLT_Jet50U_emulated_fine;
     TH1D *hist_sj_leading_pt_emulated_fine;

     TH1D *hist_sj_leading_eta_HLT_Jet15U_emulated_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_emulated_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_emulated_fine;
     TH1D *hist_sj_leading_eta_emulated_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt;
     TH1D *hist_sj_leading_eta_emulated_fine_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt;
     TH1D *hist_sj_leading_eta_emulated_fine_highpt;

     TH1D *hist_sj_leading_phi_HLT_Jet15U_emulated_fine;
     TH1D *hist_sj_leading_phi_HLT_Jet30U_emulated_fine;
     TH1D *hist_sj_leading_phi_HLT_Jet50U_emulated_fine;
     TH1D *hist_sj_leading_phi_emulated_fine;

//efficiency distributions with fine binning for leading jet
     TH1D *hist_sj_leading_pt_HLT_Jet15U_eff_fine;
     TH1D *hist_sj_leading_pt_HLT_Jet30U_eff_fine;
     TH1D *hist_sj_leading_pt_HLT_Jet50U_eff_fine;
     TH1D *hist_sj_leading_pt_eff_fine;

     TH1D *hist_sj_leading_eta_HLT_Jet15U_eff_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_eff_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_eff_fine;
     TH1D *hist_sj_leading_eta_eff_fine;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_eff_fine_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_eff_fine_lowpt;
     TH1D *hist_sj_leading_eta_eff_fine_lowpt;
     TH1D *hist_sj_leading_eta_HLT_Jet30U_eff_fine_highpt;
     TH1D *hist_sj_leading_eta_HLT_Jet50U_eff_fine_highpt;
     TH1D *hist_sj_leading_eta_eff_fine_highpt;

     TH1D *hist_sj_leading_phi_HLT_Jet15U_eff_fine;
     TH1D *hist_sj_leading_phi_HLT_Jet30U_eff_fine;
     TH1D *hist_sj_leading_phi_HLT_Jet50U_eff_fine;
     TH1D *hist_sj_leading_phi_eff_fine;

//distributions for the observables
     TH1D *hist_delta_phi_all;
     TH1D *hist_delta_eta_all;
     TH1D *hist_delta_phi_emulated;
     TH1D *hist_delta_eta_emulated;
     TH1D *hist_delta_phi_eff;
     TH1D *hist_delta_eta_eff;

     TH1D *hist_delta_phi_HLT_Jet15U_all;
     TH1D *hist_delta_eta_HLT_Jet15U_all;
     TH1D *hist_delta_phi_HLT_Jet15U_emulated;
     TH1D *hist_delta_eta_HLT_Jet15U_emulated;
     TH1D *hist_delta_phi_HLT_Jet15U_eff;
     TH1D *hist_delta_eta_HLT_Jet15U_eff;

     TH1D *hist_delta_phi_HLT_Jet30U_all;
     TH1D *hist_delta_eta_HLT_Jet30U_all;
     TH1D *hist_delta_phi_HLT_Jet30U_emulated;
     TH1D *hist_delta_eta_HLT_Jet30U_emulated;
     TH1D *hist_delta_phi_HLT_Jet30U_eff;
     TH1D *hist_delta_eta_HLT_Jet30U_eff;

     TH1D *hist_delta_phi_HLT_Jet50U_all;
     TH1D *hist_delta_eta_HLT_Jet50U_all;
     TH1D *hist_delta_phi_HLT_Jet50U_emulated;
     TH1D *hist_delta_eta_HLT_Jet50U_emulated;
     TH1D *hist_delta_phi_HLT_Jet50U_eff;
     TH1D *hist_delta_eta_HLT_Jet50U_eff;

//declaring histogram binning
//general control event count
     hist_events =  new TH1D("Events","Selection Chain;Type;# Events", 9, 0, 9);
     hist_events_HLT_Jet15U =  new TH1D("Events_HLT_Jet15U","Selection Chain for HLT_Jet15U;Type;# Events", 6, 0, 6);
     hist_events_HLT_Jet30U =  new TH1D("Events_HLT_Jet30U","Selection Chain for HLT_Jet30U;Type;# Events", 6, 0, 6);
     hist_events_HLT_Jet50U =  new TH1D("Events_HLT_Jet50U","Selection Chain for HLT_Jet50U;Type;# Events", 6, 0, 6);
     hist_events_all =  new TH1D("Events_Combined","Selection Chain for combination;Type;# Events", 6, 0, 6);
 
//monitor distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_all =  new TH1D("ak5PF_leading_pt_HLT15U_all","Leading Jet p_{T} for HLT_Jet15U all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_HLT_Jet30U_all =  new TH1D("ak5PF_leading_pt_HLT30U_all","Leading Jet p_{T} for HLT_Jet30U all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_HLT_Jet50U_all =  new TH1D("ak5PF_leading_pt_HLT50U_all","Leading Jet p_{T} for HLT_Jet50U all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_all =  new TH1D("ak5PF_leading_pt_all","Leading Jet p_{T} for triggers combination in all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);

     hist_leading_central_pt_HLT_Jet15U_all =  new TH1D("ak5PF_leading_central_pt_HLT15U_all","Leading Central Jet p_{T} for HLT_Jet15U all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_central_pt_HLT_Jet30U_all =  new TH1D("ak5PF_leading_central_pt_HLT30U_all","Leading Central Jet p_{T} for HLT_Jet30U all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_central_pt_HLT_Jet50U_all =  new TH1D("ak5PF_leading_central_pt_HLT50U_all","Leading Central Jet p_{T} for HLT_Jet50U all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_central_pt_all =  new TH1D("ak5PF_leading_central_pt_all","Leading Central Jet p_{T} for triggers combination in all events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);

     hist_leading_eta_HLT_Jet15U_all =  new TH1D("ak5PF_leading_eta_HLT15U_all","Leading Jet #eta for HLT_Jet15U all events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_all =  new TH1D("ak5PF_leading_eta_HLT30U_all","Leading Jet #eta for HLT_Jet30U all events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_all =  new TH1D("ak5PF_leading_eta_HLT50U_all","Leading Jet #eta for HLT_Jet50U all events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_all =  new TH1D("ak5PF_leading_eta_all","Leading Jet #eta for triggers combination in all events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_all_lowpt =  new TH1D("ak5PF_leading_eta_HLT30U_all_lowpt","Leading Jet #eta for HLT_Jet30U in all low pt events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_all_lowpt =  new TH1D("ak5PF_leading_eta_HLT50U_all_lowpt","Leading Jet #eta for HLT_Jet50U in all low pt events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_all_lowpt =  new TH1D("ak5PF_leading_eta_all_lowpt","Leading Jet #eta for triggers combination in all low pt events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_all_highpt =  new TH1D("ak5PF_leading_eta_HLT30U_all_highpt","Leading Jet #eta for HLT_Jet30U in all high pt events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_all_highpt =  new TH1D("ak5PF_leading_eta_HLT50U_all_highpt","Leading Jet #eta for HLT_Jet50U in all high pt events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_all_highpt =  new TH1D("ak5PF_leading_eta_all_highpt","Leading Jet #eta for triggers combination in all high pt events;#eta;# events", 20, -5.0, 5.0);

     hist_leading_phi_HLT_Jet15U_all =  new TH1D("ak5PF_leading_phi_HLT15U_all","Leading Jet #phi for HLT_Jet15U all events;#phi;# events", 14, -3.15, 3.15);
     hist_leading_phi_HLT_Jet30U_all =  new TH1D("ak5PF_leading_phi_HLT30U_all","Leading Jet #phi for HLT_Jet30U all events;#phi;# events", 14, -3.15, 3.15);
     hist_leading_phi_HLT_Jet50U_all =  new TH1D("ak5PF_leading_phi_HLT50U_all","Leading Jet #phi for HLT_Jet50U all events;#phi;# events", 14, -3.15, 3.15);
     hist_leading_phi_all =  new TH1D("ak5PF_leading_phi_all","Leading Jet #phi for triggers combination in all events;#phi;# events", 14, -3.15, 3.15);

//monitor distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_all =  new TH1D("ak5PF_sj_leading_pt_HLT15U_all","Leading Jet p_{T} for HLT_Jet15U all events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_HLT_Jet30U_all =  new TH1D("ak5PF_sj_leading_pt_HLT30U_all","Leading Jet p_{T} for HLT_Jet30U all events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_HLT_Jet50U_all =  new TH1D("ak5PF_sj_leading_pt_HLT50U_all","Leading Jet p_{T} for HLT_Jet50U all events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_all =  new TH1D("ak5PF_sj_leading_pt_all","Leading Jet p_{T} for triggers combination in all events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);

     hist_sj_leading_eta_HLT_Jet15U_all =  new TH1D("ak5PF_sj_leading_eta_HLT15U_all","Leading Jet #eta for HLT_Jet15U all events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_all =  new TH1D("ak5PF_sj_leading_eta_HLT30U_all","Leading Jet #eta for HLT_Jet30U all events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_all =  new TH1D("ak5PF_sj_leading_eta_HLT50U_all","Leading Jet #eta for HLT_Jet50U all events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_all =  new TH1D("ak5PF_sj_leading_eta_all","Leading Jet #eta for triggers combination in all events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_all_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_all_lowpt","Leading Jet #eta for HLT_Jet30U in all low pt events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_all_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_all_lowpt","Leading Jet #eta for HLT_Jet50U in all low pt events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_all_lowpt =  new TH1D("ak5PF_sj_leading_eta_all_lowpt","Leading Jet #eta for triggers combination in all low pt events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_all_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_all_highpt","Leading Jet #eta for HLT_Jet30U in all high pt events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_all_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_all_highpt","Leading Jet #eta for HLT_Jet50U in all high pt events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_all_highpt =  new TH1D("ak5PF_sj_leading_eta_all_highpt","Leading Jet #eta for triggers combination in all high pt events SJ Selection;#eta;# events", 20, -5.0, 5.0);

     hist_sj_leading_phi_HLT_Jet15U_all =  new TH1D("ak5PF_sj_leading_phi_HLT15U_all","Leading Jet #phi for HLT_Jet15U all events SJ Selection;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet30U_all =  new TH1D("ak5PF_sj_leading_phi_HLT30U_all","Leading Jet #phi for HLT_Jet30U all events SJ Selection;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet50U_all =  new TH1D("ak5PF_sj_leading_phi_HLT50U_all","Leading Jet #phi for HLT_Jet50U all events SJ Selection;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_phi_all =  new TH1D("ak5PF_sj_leading_phi_all","Leading Jet #phi for triggers combination in all events SJ Selection;#phi;# events", 14, -3.15, 3.15);

//triggered distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_emulated =  new TH1D("ak5PF_leading_pt_HLT15U_emulated","Leading Jet p_{T} for HLT_Jet15U emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_HLT_Jet30U_emulated =  new TH1D("ak5PF_leading_pt_HLT30U_emulated","Leading Jet p_{T} for HLT_Jet30U emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_HLT_Jet50U_emulated =  new TH1D("ak5PF_leading_pt_HLT50U_emulated","Leading Jet p_{T} for HLT_Jet50U emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_pt_emulated =  new TH1D("ak5PF_leading_pt_emulated","Leading Jet p_{T} for trigger combination in emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);

     hist_leading_central_pt_HLT_Jet15U_emulated =  new TH1D("ak5PF_leading_central_pt_HLT15U_emulated","Leading Central Jet p_{T} for HLT_Jet15U emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_central_pt_HLT_Jet30U_emulated =  new TH1D("ak5PF_leading_central_pt_HLT30U_emulated","Leading Central Jet p_{T} for HLT_Jet30U emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_central_pt_HLT_Jet50U_emulated =  new TH1D("ak5PF_leading_central_pt_HLT50U_emulated","Leading Central Jet p_{T} for HLT_Jet50U emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_leading_central_pt_emulated =  new TH1D("ak5PF_leading_central_pt_emulated","Leading Central Jet p_{T} for trigger combination in emulated events;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);

     hist_leading_eta_HLT_Jet15U_emulated =  new TH1D("ak5PF_leading_eta_HLT15U_emulated","Leading Jet #eta for HLT_Jet15U emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_emulated =  new TH1D("ak5PF_leading_eta_HLT30U_emulated","Leading Jet #eta for HLT_Jet30U emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_emulated =  new TH1D("ak5PF_leading_eta_HLT50U_emulated","Leading Jet #eta for HLT_Jet50U emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_emulated =  new TH1D("ak5PF_leading_eta_emulated","Leading Jet #eta for trigger combination in emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_emulated_lowpt =  new TH1D("ak5PF_leading_eta_HLT30U_emulated_lowpt","Leading Jet #eta for HLT_Jet30U in low pt emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_emulated_lowpt =  new TH1D("ak5PF_leading_eta_HLT50U_emulated_lowpt","Leading Jet #eta for HLT_Jet50U in low pt emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_emulated_lowpt =  new TH1D("ak5PF_leading_eta_emulated_lowpt","Leading Jet #eta for trigger combination in low pt emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_emulated_highpt =  new TH1D("ak5PF_leading_eta_HLT30U_emulated_highpt","Leading Jet #eta for HLT_Jet30U in high pt emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_emulated_highpt =  new TH1D("ak5PF_leading_eta_HLT50U_emulated_highpt","Leading Jet #eta for HLT_Jet50U in high pt emulated events;#eta;# events", 20, -5.0, 5.0);
     hist_leading_eta_emulated_highpt =  new TH1D("ak5PF_leading_eta_emulated_highpt","Leading Jet #eta for trigger combination in high pt emulated events;#eta;# events", 20, -5.0, 5.0);     

     hist_leading_phi_HLT_Jet15U_emulated =  new TH1D("ak5PF_leading_phi_HLT15U_emulated","Leading Jet #phi for HLT_Jet15U emulated events;#phi;# events", 14, -3.15, 3.15);
     hist_leading_phi_HLT_Jet30U_emulated =  new TH1D("ak5PF_leading_phi_HLT30U_emulated","Leading Jet #phi for HLT_Jet30U emulated events;#phi;# events", 14, -3.15, 3.15);
     hist_leading_phi_HLT_Jet50U_emulated =  new TH1D("ak5PF_leading_phi_HLT50U_emulated","Leading Jet #phi for HLT_Jet50U emulated events;#phi;# events", 14, -3.15, 3.15);
     hist_leading_phi_emulated =  new TH1D("ak5PF_leading_phi_emulated","Leading Jet #phi for trigger combination in emulated events;#phi;# events", 14, -3.15, 3.15);

//triggered distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_emulated =  new TH1D("ak5PF_sj_leading_pt_HLT15U_emulated","Leading Jet p_{T} for HLT_Jet15U emulated events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_HLT_Jet30U_emulated =  new TH1D("ak5PF_sj_leading_pt_HLT30U_emulated","Leading Jet p_{T} for HLT_Jet30U emulated events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_HLT_Jet50U_emulated =  new TH1D("ak5PF_sj_leading_pt_HLT50U_emulated","Leading Jet p_{T} for HLT_Jet50U emulated events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);
     hist_sj_leading_pt_emulated =  new TH1D("ak5PF_sj_leading_pt_emulated","Leading Jet p_{T} for trigger combination in emulated events SJ Selection;p_{T} [#frac{GeV}{c}];# events", all_nbins, all_bins);

     hist_sj_leading_eta_HLT_Jet15U_emulated =  new TH1D("ak5PF_sj_leading_eta_HLT15U_emulated","Leading Jet #eta for HLT_Jet15U emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_emulated =  new TH1D("ak5PF_sj_leading_eta_HLT30U_emulated","Leading Jet #eta for HLT_Jet30U emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_emulated =  new TH1D("ak5PF_sj_leading_eta_HLT50U_emulated","Leading Jet #eta for HLT_Jet50U emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_emulated =  new TH1D("ak5PF_sj_leading_eta_emulated","Leading Jet #eta for trigger combination in emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_emulated_lowpt","Leading Jet #eta for HLT_Jet30U in low pt emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_emulated_lowpt","Leading Jet #eta for HLT_Jet50U in low pt emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_emulated_lowpt =  new TH1D("ak5PF_sj_leading_eta_emulated_lowpt","Leading Jet #eta for trigger combination in low pt emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_emulated_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_emulated_highpt","Leading Jet #eta for HLT_Jet30U in high pt emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_emulated_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_emulated_highpt","Leading Jet #eta for HLT_Jet50U in high pt emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);
     hist_sj_leading_eta_emulated_highpt =  new TH1D("ak5PF_sj_leading_eta_emulated_highpt","Leading Jet #eta for trigger combination in high pt emulated events SJ Selection;#eta;# events", 20, -5.0, 5.0);     

     hist_sj_leading_phi_HLT_Jet15U_emulated =  new TH1D("ak5PF_sj_leading_phi_HLT15U_emulated","Leading Jet #phi for HLT_Jet15U emulated events SJ Selection;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet30U_emulated =  new TH1D("ak5PF_sj_leading_phi_HLT30U_emulated","Leading Jet #phi for HLT_Jet30U emulated events SJ Selection;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet50U_emulated =  new TH1D("ak5PF_sj_leading_phi_HLT50U_emulated","Leading Jet #phi for HLT_Jet50U emulated events SJ Selection;#phi;# events", 14, -3.15, 3.15);
     hist_sj_leading_phi_emulated =  new TH1D("ak5PF_sj_leading_phi_emulated","Leading Jet #phi for trigger combination in emulated events SJ Selection;#phi;# events", 14, -3.15, 3.15);

//efficiency distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_eff =  new TH1D("ak5PF_leading_pt_HLT15U_eff","Leading Jet p_{T} for HLT_Jet15U efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_leading_pt_HLT_Jet30U_eff =  new TH1D("ak5PF_leading_pt_HLT30U_eff","Leading Jet p_{T} for HLT_Jet30U efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_leading_pt_HLT_Jet50U_eff =  new TH1D("ak5PF_leading_pt_HLT50U_eff","Leading Jet p_{T} for HLT_Jet50U efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_leading_pt_eff =  new TH1D("ak5PF_leading_pt_eff","Leading Jet p_{T} for trigger combination efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);

     hist_leading_central_pt_HLT_Jet15U_eff =  new TH1D("ak5PF_leading_central_pt_HLT15U_eff","Leading Central Jet p_{T} for HLT_Jet15U efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_leading_central_pt_HLT_Jet30U_eff =  new TH1D("ak5PF_leading_central_pt_HLT30U_eff","Leading Central Jet p_{T} for HLT_Jet30U efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_leading_central_pt_HLT_Jet50U_eff =  new TH1D("ak5PF_leading_central_pt_HLT50U_eff","Leading Central Jet p_{T} for HLT_Jet50U efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_leading_central_pt_eff =  new TH1D("ak5PF_leading_central_pt_eff","Leading Central Jet p_{T} for trigger combination efficiency;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);

     hist_leading_eta_HLT_Jet15U_eff =  new TH1D("ak5PF_leading_eta_HLT15U_eff","Leading Jet #eta for HLT_Jet15U efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_eff =  new TH1D("ak5PF_leading_eta_HLT30U_eff","Leading Jet #eta for HLT_Jet30U efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_eff =  new TH1D("ak5PF_leading_eta_HLT50U_eff","Leading Jet #eta for HLT_Jet50U efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_eff =  new TH1D("ak5PF_leading_eta_eff","Leading Jet #eta for trigger combination efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_eff_lowpt =  new TH1D("ak5PF_leading_eta_HLT30U_eff_lowpt","Leading Jet #eta for HLT_Jet30U low pt efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_eff_lowpt =  new TH1D("ak5PF_leading_eta_HLT50U_eff_lowpt","Leading Jet #eta for HLT_Jet50U lowe pt efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_eff_lowpt =  new TH1D("ak5PF_leading_eta_eff_lowpt","Leading Jet #eta for trigger combination low pt efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_eff_highpt =  new TH1D("ak5PF_leading_eta_HLT30U_eff_highpt","Leading Jet #eta for HLT_Jet30U high pt efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_eff_highpt =  new TH1D("ak5PF_leading_eta_HLT50U_eff_highpt","Leading Jet #eta for HLT_Jet50U high pt efficiency;#eta;efficiency", 20, -5.0, 5.0);
     hist_leading_eta_eff_highpt =  new TH1D("ak5PF_leading_eta_eff_highpt","Leading Jet #eta for trigger combination high pt efficiency;#eta;efficiency", 20, -5.0, 5.0);

     hist_leading_phi_HLT_Jet15U_eff =  new TH1D("ak5PF_leading_phi_HLT15U_eff","Leading Jet #phi for HLT_Jet15U efficiency;#phi;efficiency",  14, -3.15, 3.15);
     hist_leading_phi_HLT_Jet30U_eff =  new TH1D("ak5PF_leading_phi_HLT30U_eff","Leading Jet #phi for HLT_Jet30U efficiency;#phi;efficiency",  14, -3.15, 3.15);
     hist_leading_phi_HLT_Jet50U_eff =  new TH1D("ak5PF_leading_phi_HLT50U_eff","Leading Jet #phi for HLT_Jet50U efficiency;#phi;efficiency",  14, -3.15, 3.15);
     hist_leading_phi_eff =  new TH1D("ak5PF_leading_phi_eff","Leading Jet #phi for trigger combination efficiency;#phi;efficiency",  14, -3.15, 3.15);

//efficiency distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff =  new TH1D("ak5PF_sj_leading_pt_HLT15U_eff","Leading Jet p_{T} for HLT_Jet15U efficiency SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_sj_leading_pt_HLT_Jet30U_eff =  new TH1D("ak5PF_sj_leading_pt_HLT30U_eff","Leading Jet p_{T} for HLT_Jet30U efficiency SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_sj_leading_pt_HLT_Jet50U_eff =  new TH1D("ak5PF_sj_leading_pt_HLT50U_eff","Leading Jet p_{T} for HLT_Jet50U efficiency SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);
     hist_sj_leading_pt_eff =  new TH1D("ak5PF_sj_leading_pt_eff","Leading Jet p_{T} for trigger combination efficiency SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", all_nbins, all_bins);

     hist_sj_leading_eta_HLT_Jet15U_eff =  new TH1D("ak5PF_sj_leading_eta_HLT15U_eff","Leading Jet #eta for HLT_Jet15U efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_eff =  new TH1D("ak5PF_sj_leading_eta_HLT30U_eff","Leading Jet #eta for HLT_Jet30U efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_eff =  new TH1D("ak5PF_sj_leading_eta_HLT50U_eff","Leading Jet #eta for HLT_Jet50U efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_eff =  new TH1D("ak5PF_sj_leading_eta_eff","Leading Jet #eta for trigger combination efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_eff_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_eff_lowpt","Leading Jet #eta for HLT_Jet30U low pt efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_eff_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_eff_lowpt","Leading Jet #eta for HLT_Jet50U lowe pt efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_eff_lowpt =  new TH1D("ak5PF_sj_leading_eta_eff_lowpt","Leading Jet #eta for trigger combination low pt efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_eff_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_eff_highpt","Leading Jet #eta for HLT_Jet30U high pt efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_eff_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_eff_highpt","Leading Jet #eta for HLT_Jet50U high pt efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);
     hist_sj_leading_eta_eff_highpt =  new TH1D("ak5PF_sj_leading_eta_eff_highpt","Leading Jet #eta for trigger combination high pt efficiency SJ Selection;#eta;efficiency", 20, -5.0, 5.0);

     hist_sj_leading_phi_HLT_Jet15U_eff =  new TH1D("ak5PF_sj_leading_phi_HLT15U_eff","Leading Jet #phi for HLT_Jet15U efficiency SJ Selection;#phi;efficiency",  14, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet30U_eff =  new TH1D("ak5PF_sj_leading_phi_HLT30U_eff","Leading Jet #phi for HLT_Jet30U efficiency SJ Selection;#phi;efficiency",  14, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet50U_eff =  new TH1D("ak5PF_sj_leading_phi_HLT50U_eff","Leading Jet #phi for HLT_Jet50U efficiency SJ Selection;#phi;efficiency",  14, -3.15, 3.15);
     hist_sj_leading_phi_eff =  new TH1D("ak5PF_sj_leading_phi_eff","Leading Jet #phi for trigger combination efficiency SJ Selection;#phi;efficiency",  14, -3.15, 3.15);

//monitor distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_all_fine =  new TH1D("ak5PF_leading_pt_HLT15U_all_fine","Leading Jet p_{T} for HLT_Jet15U all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_pt_HLT_Jet30U_all_fine =  new TH1D("ak5PF_leading_pt_HLT30U_all_fine","Leading Jet p_{T} for HLT_Jet30U all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_pt_HLT_Jet50U_all_fine =  new TH1D("ak5PF_leading_pt_HLT50U_all_fine","Leading Jet p_{T} for HLT_Jet50U all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_pt_all_fine =  new TH1D("ak5PF_leading_pt_all_fine","Leading Jet p_{T} for triggers combination in all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);

     hist_leading_central_pt_HLT_Jet15U_all_fine =  new TH1D("ak5PF_leading_central pt_HLT15U_all_fine","Leading Central Jet p_{T} for HLT_Jet15U all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_central_pt_HLT_Jet30U_all_fine =  new TH1D("ak5PF_leading_central_pt_HLT30U_all_fine","Leading Central Jet p_{T} for HLT_Jet30U all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_central_pt_HLT_Jet50U_all_fine =  new TH1D("ak5PF_leading_central_pt_HLT50U_all_fine","Leading Central Jet p_{T} for HLT_Jet50U all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_central_pt_all_fine =  new TH1D("ak5PF_leading_central_pt_all_fine","Leading Central Jet p_{T} for triggers combination in all events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);


     hist_leading_eta_HLT_Jet15U_all_fine =  new TH1D("ak5PF_leading_eta_HLT15U_all_fine","Leading Jet #eta for HLT_Jet15U all events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_all_fine =  new TH1D("ak5PF_leading_eta_HLT30U_all_fine","Leading Jet #eta for HLT_Jet30U all events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_all_fine =  new TH1D("ak5PF_leading_eta_HLT50U_all_fine","Leading Jet #eta for HLT_Jet50U all events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_all_fine =  new TH1D("ak5PF_leading_eta_all_fine","Leading Jet #eta for triggers combination in all events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_all_fine_lowpt =  new TH1D("ak5PF_leading_eta_HLT30U_all_fine_lowpt","Leading Jet #eta for HLT_Jet30U in all low pt events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_all_fine_lowpt =  new TH1D("ak5PF_leading_eta_HLT50U_all_fine_lowpt","Leading Jet #eta for HLT_Jet50U in all low pt events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_all_fine_lowpt =  new TH1D("ak5PF_leading_eta_all_fine_lowpt","Leading Jet #eta for triggers combination in all low pt events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_all_fine_highpt =  new TH1D("ak5PF_leading_eta_HLT30U_all_fine_highpt","Leading Jet #eta for HLT_Jet30U in all high pt events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_all_fine_highpt =  new TH1D("ak5PF_leading_eta_HLT50U_all_fine_highpt","Leading Jet #eta for HLT_Jet50U in all high pt events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_all_fine_highpt =  new TH1D("ak5PF_leading_eta_all_fine_highpt","Leading Jet #eta for triggers combination in all high pt events fine;#eta;# events", 50, -5.0, 5.0);

     hist_leading_phi_HLT_Jet15U_all_fine =  new TH1D("ak5PF_leading_phi_HLT15U_all_fine","Leading Jet #phi for HLT_Jet15U all events fine;#phi;# events", 64, -3.15, 3.15);
     hist_leading_phi_HLT_Jet30U_all_fine =  new TH1D("ak5PF_leading_phi_HLT30U_all_fine","Leading Jet #phi for HLT_Jet30U all events fine;#phi;# events", 64, -3.15, 3.15);
     hist_leading_phi_HLT_Jet50U_all_fine =  new TH1D("ak5PF_leading_phi_HLT50U_all_fine","Leading Jet #phi for HLT_Jet50U all events fine;#phi;# events", 64, -3.15, 3.15);
     hist_leading_phi_all_fine =  new TH1D("ak5PF_leading_phi_all_fine","Leading Jet #phi for triggers combination in all events fine;#phi;# events", 64, -3.15, 3.15);

//monitor distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_all_fine =  new TH1D("ak5PF_sj_leading_pt_HLT15U_all_fine","Leading Jet p_{T} for HLT_Jet15U all events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_sj_leading_pt_HLT_Jet30U_all_fine =  new TH1D("ak5PF_sj_leading_pt_HLT30U_all_fine","Leading Jet p_{T} for HLT_Jet30U all events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_sj_leading_pt_HLT_Jet50U_all_fine =  new TH1D("ak5PF_sj_leading_pt_HLT50U_all_fine","Leading Jet p_{T} for HLT_Jet50U all events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_sj_leading_pt_all_fine =  new TH1D("ak5PF_sj_leading_pt_all_fine","Leading Jet p_{T} for triggers combination in all events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);

     hist_sj_leading_eta_HLT_Jet15U_all_fine =  new TH1D("ak5PF_sj_leading_eta_HLT15U_all_fine","Leading Jet #eta for HLT_Jet15U all events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_all_fine =  new TH1D("ak5PF_sj_leading_eta_HLT30U_all_fine","Leading Jet #eta for HLT_Jet30U all events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_all_fine =  new TH1D("ak5PF_sj_leading_eta_HLT50U_all_fine","Leading Jet #eta for HLT_Jet50U all events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_all_fine =  new TH1D("ak5PF_sj_leading_eta_all_fine","Leading Jet #eta for triggers combination in all events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_all_fine_lowpt","Leading Jet #eta for HLT_Jet30U in all low pt events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_all_fine_lowpt","Leading Jet #eta for HLT_Jet50U in all low pt events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_all_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_all_fine_lowpt","Leading Jet #eta for triggers combination in all low pt events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_all_fine_highpt","Leading Jet #eta for HLT_Jet30U in all high pt events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_all_fine_highpt","Leading Jet #eta for HLT_Jet50U in all high pt events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_all_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_all_fine_highpt","Leading Jet #eta for triggers combination in all high pt events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);

     hist_sj_leading_phi_HLT_Jet15U_all_fine =  new TH1D("ak5PF_sj_leading_phi_HLT15U_all_fine","Leading Jet #phi for HLT_Jet15U all events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet30U_all_fine =  new TH1D("ak5PF_sj_leading_phi_HLT30U_all_fine","Leading Jet #phi for HLT_Jet30U all events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet50U_all_fine =  new TH1D("ak5PF_sj_leading_phi_HLT50U_all_fine","Leading Jet #phi for HLT_Jet50U all events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);
     hist_sj_leading_phi_all_fine =  new TH1D("ak5PF_sj_leading_phi_all_fine","Leading Jet #phi for triggers combination in all events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);

//triggered distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_leading_pt_HLT15U_emulated_fine","Leading Jet p_{T} for HLT_Jet15U emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_pt_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_leading_pt_HLT30U_emulated_fine","Leading Jet p_{T} for HLT_Jet30U emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_pt_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_leading_pt_HLT50U_emulated_fine","Leading Jet p_{T} for HLT_Jet50U emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_pt_emulated_fine =  new TH1D("ak5PF_leading_pt_emulated_fine","Leading Jet p_{T} for trigger combination in emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);

     hist_leading_central_pt_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_leading_central_pt_HLT15U_emulated_fine","Leading Central Jet p_{T} for HLT_Jet15U emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_central_pt_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_leading_central_pt_HLT30U_emulated_fine","Leading Central Jet p_{T} for HLT_Jet30U emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_central_pt_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_leading_central_pt_HLT50U_emulated_fine","Leading Central Jet p_{T} for HLT_Jet50U emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_leading_central_pt_emulated_fine =  new TH1D("ak5PF_leading_central_pt_emulated_fine","Leading Central Jet p_{T} for trigger combination in emulated events fine;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);

     hist_leading_eta_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_leading_eta_HLT15U_emulated_fine","Leading Jet #eta for HLT_Jet15U emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_leading_eta_HLT30U_emulated_fine","Leading Jet #eta for HLT_Jet30U emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_leading_eta_HLT50U_emulated_fine","Leading Jet #eta for HLT_Jet50U emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_emulated_fine =  new TH1D("ak5PF_leading_eta_emulated_fine","Leading Jet #eta for trigger combination in emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt =  new TH1D("ak5PF_leading_eta_HLT30U_emulated_fine_lowpt","Leading Jet #eta for HLT_Jet30U in low pt emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt =  new TH1D("ak5PF_leading_eta_HLT50U_emulated_fine_lowpt","Leading Jet #eta for HLT_Jet50U in low pt emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_emulated_fine_lowpt =  new TH1D("ak5PF_leading_eta_emulated_fine_lowpt","Leading Jet #eta for trigger combination in low pt emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_emulated_fine_highpt =  new TH1D("ak5PF_leading_eta_HLT30U_emulated_fine_highpt","Leading Jet #eta for HLT_Jet30U in high pt emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_emulated_fine_highpt =  new TH1D("ak5PF_leading_eta_HLT50U_emulated_fine_highpt","Leading Jet #eta for HLT_Jet50U in high pt emulated events fine;#eta;# events", 50, -5.0, 5.0);
     hist_leading_eta_emulated_fine_highpt =  new TH1D("ak5PF_leading_eta_emulated_fine_highpt","Leading Jet #eta for trigger combination in high pt emulated events fine;#eta;# events", 50, -5.0, 5.0);

     hist_leading_phi_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_leading_phi_HLT15U_emulated_fine","Leading Jet #phi for HLT_Jet15U emulated events fine;#phi;# events", 64, -3.15, 3.15);
     hist_leading_phi_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_leading_phi_HLT30U_emulated_fine","Leading Jet #phi for HLT_Jet30U emulated events fine;#phi;# events", 64, -3.15, 3.15);
     hist_leading_phi_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_leading_phi_HLT50U_emulated_fine","Leading Jet #phi for HLT_Jet50U emulated events fine;#phi;# events", 64, -3.15, 3.15);
     hist_leading_phi_emulated_fine =  new TH1D("ak5PF_leading_phi_emulated_fine","Leading Jet #phi for trigger combination in emulated events fine;#phi;# events", 64, -3.15, 3.15);

//triggered distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_sj_leading_pt_HLT15U_emulated_fine","Leading Jet p_{T} for HLT_Jet15U emulated events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_sj_leading_pt_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_sj_leading_pt_HLT30U_emulated_fine","Leading Jet p_{T} for HLT_Jet30U emulated events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_sj_leading_pt_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_sj_leading_pt_HLT50U_emulated_fine","Leading Jet p_{T} for HLT_Jet50U emulated events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);
     hist_sj_leading_pt_emulated_fine =  new TH1D("ak5PF_sj_leading_pt_emulated_fine","Leading Jet p_{T} for trigger combination in emulated events fine SJ Selection;p_{T} [#frac{GeV}{c}];# events", 165, 35, 200);

     hist_sj_leading_eta_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_sj_leading_eta_HLT15U_emulated_fine","Leading Jet #eta for HLT_Jet15U emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_sj_leading_eta_HLT30U_emulated_fine","Leading Jet #eta for HLT_Jet30U emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_sj_leading_eta_HLT50U_emulated_fine","Leading Jet #eta for HLT_Jet50U emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_emulated_fine =  new TH1D("ak5PF_sj_leading_eta_emulated_fine","Leading Jet #eta for trigger combination in emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_emulated_fine_lowpt","Leading Jet #eta for HLT_Jet30U in low pt emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_emulated_fine_lowpt","Leading Jet #eta for HLT_Jet50U in low pt emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_emulated_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_emulated_fine_lowpt","Leading Jet #eta for trigger combination in low pt emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_emulated_fine_highpt","Leading Jet #eta for HLT_Jet30U in high pt emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_emulated_fine_highpt","Leading Jet #eta for HLT_Jet50U in high pt emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);
     hist_sj_leading_eta_emulated_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_emulated_fine_highpt","Leading Jet #eta for trigger combination in high pt emulated events fine SJ Selection;#eta;# events", 50, -5.0, 5.0);

     hist_sj_leading_phi_HLT_Jet15U_emulated_fine =  new TH1D("ak5PF_sj_leading_phi_HLT15U_emulated_fine","Leading Jet #phi for HLT_Jet15U emulated events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet30U_emulated_fine =  new TH1D("ak5PF_sj_leading_phi_HLT30U_emulated_fine","Leading Jet #phi for HLT_Jet30U emulated events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet50U_emulated_fine =  new TH1D("ak5PF_sj_leading_phi_HLT50U_emulated_fine","Leading Jet #phi for HLT_Jet50U emulated events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);
     hist_sj_leading_phi_emulated_fine =  new TH1D("ak5PF_sj_leading_phi_emulated_fine","Leading Jet #phi for trigger combination in emulated events fine SJ Selection;#phi;# events", 64, -3.15, 3.15);

//efficiency distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_sj_leading_pt_HLT15U_eff_fine","Leading Jet p_{T} for HLT_Jet15U efficiency fine SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_sj_leading_pt_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_sj_leading_pt_HLT30U_eff_fine","Leading Jet p_{T} for HLT_Jet30U efficiency fine SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_sj_leading_pt_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_sj_leading_pt_HLT50U_eff_fine","Leading Jet p_{T} for HLT_Jet50U efficiency fine SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_sj_leading_pt_eff_fine =  new TH1D("ak5PF_sj_leading_pt_eff_fine","Leading Jet p_{T} for combination efficiency fine SJ Selection;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);

     hist_sj_leading_eta_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_sj_leading_eta_HLT15U_eff_fine","Leading Jet #eta for HLT_Jet15U efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_sj_leading_eta_HLT30U_eff_fine","Leading Jet #eta for HLT_Jet30U efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_sj_leading_eta_HLT50U_eff_fine","Leading Jet #eta for HLT_Jet50U efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_eff_fine =  new TH1D("ak5PF_sj_leading_eta_eff_fine","Leading Jet #eta for combination efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_eff_fine_lowpt","Leading Jet #eta for HLT_Jet30U low pt efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_eff_fine_lowpt","Leading Jet #eta for HLT_Jet50U low pt efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_eff_fine_lowpt =  new TH1D("ak5PF_sj_leading_eta_eff_fine_lowpt","Leading Jet #eta for combination low pt efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT30U_eff_fine_highpt","Leading Jet #eta for HLT_Jet30U high pt efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_HLT50U_eff_fine_highpt","Leading Jet #eta for HLT_Jet50U high pt efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);
     hist_sj_leading_eta_eff_fine_highpt =  new TH1D("ak5PF_sj_leading_eta_eff_fine_highpt","Leading Jet #eta for combination high pt efficiency fine SJ Selection;#eta;efficiency", 50, -5.0, 5.0);

     hist_sj_leading_phi_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_sj_leading_phi_HLT15U_eff_fine","Leading Jet #phi for HLT_Jet15U efficiency fine SJ Selection;#phi;efficiency", 64, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_sj_leading_phi_HLT30U_eff_fine","Leading Jet #phi for HLT_Jet30U efficiency fine SJ Selection;#phi;efficiency", 64, -3.15, 3.15);
     hist_sj_leading_phi_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_sj_leading_phi_HLT50U_eff_fine","Leading Jet #phi for HLT_Jet50U efficiency fine SJ Selection;#phi;efficiency", 64, -3.15, 3.15);
     hist_sj_leading_phi_eff_fine =  new TH1D("ak5PF_sj_leading_phi_eff_fine","Leading Jet #phi for combination efficiency fine SJ Selection;#phi;efficiency", 64, -3.15, 3.15);


//distributions for the observables
     hist_delta_phi_all =  new TH1D("ak5PF_delta_phi_all","#Delta#phi in all events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_all =  new TH1D("ak5PF_delta_eta_all","#Delta#eta in all events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_emulated =  new TH1D("ak5PF_delta_phi_emulated","#Delta#phi in emulated events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_emulated =  new TH1D("ak5PF_delta_eta_emulated","#Delta#eta in emulated events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_eff =  new TH1D("ak5PF_delta_phi_eff","#Delta#phi in efficiency;#Delta#phi [rad];efficiency", dphi_nbins, dphi_bins);
     hist_delta_eta_eff =  new TH1D("ak5PF_delta_eta_eff","#Delta#eta in efficiency;#Delta#eta;efficiency", deta_nbins, deta_bins);

     hist_delta_phi_HLT_Jet15U_all =  new TH1D("ak5PF_delta_phi_HLT_Jet15U_all","#Delta#phi for HLT_Jet15U in all events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet15U_all =  new TH1D("ak5PF_delta_eta_HLT_Jet15U_all","#Delta#eta for HLT_Jet15U in all events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_HLT_Jet15U_emulated =  new TH1D("ak5PF_delta_phi_HLT_Jet15U_emulated","#Delta#phi for HLT_Jet15U in emulated events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet15U_emulated =  new TH1D("ak5PF_delta_eta_HLT_Jet15U_emulated","#Delta#eta for HLT_Jet15U in emulated events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_HLT_Jet15U_eff =  new TH1D("ak5PF_delta_phi_HLT_Jet15U_eff","#Delta#phi for HLT_Jet15U in efficiency;#Delta#phi [rad];efficiency", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet15U_eff =  new TH1D("ak5PF_delta_eta_HLT_Jet15U_eff","#Delta#eta for HLT_Jet15U in efficiency;#Delta#eta;efficiency", deta_nbins, deta_bins);

     hist_delta_phi_HLT_Jet30U_all =  new TH1D("ak5PF_delta_phi_HLT_Jet30U_all","#Delta#phi for HLT_Jet30U in all events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet30U_all =  new TH1D("ak5PF_delta_eta_HLT_Jet30U_all","#Delta#eta for HLT_Jet30U in all events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_HLT_Jet30U_emulated =  new TH1D("ak5PF_delta_phi_HLT_Jet30U_emulated","#Delta#phi for HLT_Jet30U in emulated events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet30U_emulated =  new TH1D("ak5PF_delta_eta_HLT_Jet30U_emulated","#Delta#eta for HLT_Jet30U in emulated events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_HLT_Jet30U_eff =  new TH1D("ak5PF_delta_phi_HLT_Jet30U_eff","#Delta#phi for HLT_Jet30U in efficiency;#Delta#phi [rad];efficiency", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet30U_eff =  new TH1D("ak5PF_delta_eta_HLT_Jet30U_eff","#Delta#eta for HLT_Jet30U in efficiency;#Delta#eta [rad];efficiency", deta_nbins, deta_bins);

     hist_delta_phi_HLT_Jet50U_all =  new TH1D("ak5PF_delta_phi_HLT_Jet50U_all","#Delta#phi for HLT_Jet50U in all events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet50U_all =  new TH1D("ak5PF_delta_eta_HLT_Jet50U_all","#Delta#eta for HLT_Jet50U in all events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_HLT_Jet50U_emulated =  new TH1D("ak5PF_delta_phi_HLT_Jet50U_emulated","#Delta#phi for HLT_Jet50U in emulated events;#Delta#phi [rad];# events", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet50U_emulated =  new TH1D("ak5PF_delta_eta_HLT_Jet50U_emulated","#Delta#eta for HLT_Jet50U in emulated events;#Delta#eta;# events", deta_nbins, deta_bins);
     hist_delta_phi_HLT_Jet50U_eff =  new TH1D("ak5PF_delta_phi_HLT_Jet50U_eff","#Delta#phi for HLT_Jet50U in efficiency;#Delta#phi [rad];efficiency", dphi_nbins, dphi_bins);
     hist_delta_eta_HLT_Jet50U_eff =  new TH1D("ak5PF_delta_eta_HLT_Jet50U_eff","#Delta#eta for HLT_Jet50U in efficiency;#Delta#eta;efficiency", deta_nbins, deta_bins);

//inicializing the histogram statistics
//general control label setup
     hist_events->Sumw2();
     hist_events->Fill("Total Events",0);
     hist_events->Fill("PV Selection",0);
     hist_events->Fill("Trigger Selection",0);
     hist_events->Fill("Selected",0);
     hist_events->Fill("L1 Objects",0);
     hist_events->Fill("HLT Objects",0);
     hist_events->Fill("Number of Jets",0);
     hist_events->Fill("Pileup Scale",0);
     hist_events->Fill("Efficiency",0);

     hist_events_HLT_Jet15U->Sumw2();
     hist_events_HLT_Jet15U->Fill("Triggered Events",0);
     hist_events_HLT_Jet15U->Fill("Passed L1 Selection",0);
     hist_events_HLT_Jet15U->Fill("Passed HLT Selection",0);
     hist_events_HLT_Jet15U->Fill("Passed Condition",0);
     hist_events_HLT_Jet15U->Fill("Selected",0);
     hist_events_HLT_Jet15U->Fill("Selected & Emulation",0);

     hist_events_HLT_Jet30U->Sumw2();
     hist_events_HLT_Jet30U->Fill("Triggered Events",0);
     hist_events_HLT_Jet30U->Fill("Passed L1 Selection",0);
     hist_events_HLT_Jet30U->Fill("Passed HLT Selection",0);
     hist_events_HLT_Jet30U->Fill("Passed Condition",0);
     hist_events_HLT_Jet30U->Fill("Selected",0);
     hist_events_HLT_Jet30U->Fill("Selected & Emulation",0);

     hist_events_HLT_Jet50U->Sumw2();
     hist_events_HLT_Jet50U->Fill("Triggered Events",0);
     hist_events_HLT_Jet50U->Fill("Passed L1 Selection",0);
     hist_events_HLT_Jet50U->Fill("Passed HLT Selection",0);
     hist_events_HLT_Jet50U->Fill("Passed Condition",0);
     hist_events_HLT_Jet50U->Fill("Selected",0);
     hist_events_HLT_Jet50U->Fill("Selected & Emulation",0);

     hist_events_all->Sumw2();
     hist_events_all->Fill("Triggered Events",0);
     hist_events_all->Fill("Passed L1 Selection",0);
     hist_events_all->Fill("Passed HLT Selection",0);
     hist_events_all->Fill("Passed Condition",0);
     hist_events_all->Fill("Selected",0);
     hist_events_all->Fill("Selected & Emulation",0);

//monitor distributions with broad binning for dijet
     hist_leading_pt_HLT_Jet15U_all->Sumw2();
     hist_leading_pt_HLT_Jet30U_all->Sumw2();
     hist_leading_pt_HLT_Jet50U_all->Sumw2();
     hist_leading_pt_all->Sumw2();

     hist_leading_central_pt_HLT_Jet15U_all->Sumw2();
     hist_leading_central_pt_HLT_Jet30U_all->Sumw2();
     hist_leading_central_pt_HLT_Jet50U_all->Sumw2();
     hist_leading_central_pt_all->Sumw2();

     hist_leading_eta_HLT_Jet15U_all->Sumw2();
     hist_leading_eta_HLT_Jet30U_all->Sumw2();
     hist_leading_eta_HLT_Jet50U_all->Sumw2();
     hist_leading_eta_all->Sumw2();
     hist_leading_eta_HLT_Jet30U_all_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_all_lowpt->Sumw2();
     hist_leading_eta_all_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet30U_all_highpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_all_highpt->Sumw2();
     hist_leading_eta_all_highpt->Sumw2();

     hist_leading_phi_HLT_Jet15U_all->Sumw2();
     hist_leading_phi_HLT_Jet30U_all->Sumw2();
     hist_leading_phi_HLT_Jet50U_all->Sumw2();
     hist_leading_phi_all->Sumw2();

//monitor distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_all->Sumw2();
     hist_sj_leading_pt_HLT_Jet30U_all->Sumw2();
     hist_sj_leading_pt_HLT_Jet50U_all->Sumw2();
     hist_sj_leading_pt_all->Sumw2();

     hist_sj_leading_eta_HLT_Jet15U_all->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_all->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_all->Sumw2();
     hist_sj_leading_eta_all->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_all_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_all_lowpt->Sumw2();
     hist_sj_leading_eta_all_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_all_highpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_all_highpt->Sumw2();
     hist_sj_leading_eta_all_highpt->Sumw2();

     hist_sj_leading_phi_HLT_Jet15U_all->Sumw2();
     hist_sj_leading_phi_HLT_Jet30U_all->Sumw2();
     hist_sj_leading_phi_HLT_Jet50U_all->Sumw2();
     hist_sj_leading_phi_all->Sumw2();

//triggered distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_emulated->Sumw2();
     hist_leading_pt_HLT_Jet30U_emulated->Sumw2();
     hist_leading_pt_HLT_Jet50U_emulated->Sumw2();
     hist_leading_pt_emulated->Sumw2();

     hist_leading_central_pt_HLT_Jet15U_emulated->Sumw2();
     hist_leading_central_pt_HLT_Jet30U_emulated->Sumw2();
     hist_leading_central_pt_HLT_Jet50U_emulated->Sumw2();
     hist_leading_central_pt_emulated->Sumw2();

     hist_leading_eta_HLT_Jet15U_emulated->Sumw2();
     hist_leading_eta_HLT_Jet30U_emulated->Sumw2();
     hist_leading_eta_HLT_Jet50U_emulated->Sumw2();
     hist_leading_eta_emulated->Sumw2();
     hist_leading_eta_HLT_Jet30U_emulated_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_emulated_lowpt->Sumw2();
     hist_leading_eta_emulated_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet30U_emulated_highpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_emulated_highpt->Sumw2();
     hist_leading_eta_emulated_highpt->Sumw2();
     
     hist_leading_phi_HLT_Jet15U_emulated->Sumw2();
     hist_leading_phi_HLT_Jet30U_emulated->Sumw2();
     hist_leading_phi_HLT_Jet50U_emulated->Sumw2();
     hist_leading_phi_emulated->Sumw2();

//triggered distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_emulated->Sumw2();
     hist_sj_leading_pt_HLT_Jet30U_emulated->Sumw2();
     hist_sj_leading_pt_HLT_Jet50U_emulated->Sumw2();
     hist_sj_leading_pt_emulated->Sumw2();

     hist_sj_leading_eta_HLT_Jet15U_emulated->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_emulated->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_emulated->Sumw2();
     hist_sj_leading_eta_emulated->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt->Sumw2();
     hist_sj_leading_eta_emulated_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_emulated_highpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_emulated_highpt->Sumw2();
     hist_sj_leading_eta_emulated_highpt->Sumw2();
     
     hist_sj_leading_phi_HLT_Jet15U_emulated->Sumw2();
     hist_sj_leading_phi_HLT_Jet30U_emulated->Sumw2();
     hist_sj_leading_phi_HLT_Jet50U_emulated->Sumw2();
     hist_sj_leading_phi_emulated->Sumw2();

//efficiency distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_eff->Sumw2();
     hist_leading_pt_HLT_Jet30U_eff->Sumw2();
     hist_leading_pt_HLT_Jet50U_eff->Sumw2();
     hist_leading_pt_eff->Sumw2();

     hist_leading_central_pt_HLT_Jet15U_eff->Sumw2();
     hist_leading_central_pt_HLT_Jet30U_eff->Sumw2();
     hist_leading_central_pt_HLT_Jet50U_eff->Sumw2();
     hist_leading_central_pt_eff->Sumw2();

     hist_leading_eta_HLT_Jet15U_eff->Sumw2();
     hist_leading_eta_HLT_Jet30U_eff->Sumw2();
     hist_leading_eta_HLT_Jet50U_eff->Sumw2();
     hist_leading_eta_eff->Sumw2();
     hist_leading_eta_HLT_Jet30U_eff_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_eff_lowpt->Sumw2();
     hist_leading_eta_eff_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet30U_eff_highpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_eff_highpt->Sumw2();
     hist_leading_eta_eff_highpt->Sumw2();

     hist_leading_phi_HLT_Jet15U_eff->Sumw2();
     hist_leading_phi_HLT_Jet30U_eff->Sumw2();
     hist_leading_phi_HLT_Jet50U_eff->Sumw2();
     hist_leading_phi_eff->Sumw2();

//efficiency distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff->Sumw2();
     hist_sj_leading_pt_HLT_Jet30U_eff->Sumw2();
     hist_sj_leading_pt_HLT_Jet50U_eff->Sumw2();
     hist_sj_leading_pt_eff->Sumw2();

     hist_sj_leading_eta_HLT_Jet15U_eff->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_eff->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_eff->Sumw2();
     hist_sj_leading_eta_eff->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_eff_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_eff_lowpt->Sumw2();
     hist_sj_leading_eta_eff_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_eff_highpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_eff_highpt->Sumw2();
     hist_sj_leading_eta_eff_highpt->Sumw2();

     hist_sj_leading_phi_HLT_Jet15U_eff->Sumw2();
     hist_sj_leading_phi_HLT_Jet30U_eff->Sumw2();
     hist_sj_leading_phi_HLT_Jet50U_eff->Sumw2();
     hist_sj_leading_phi_eff->Sumw2();

//monitor distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_all_fine->Sumw2();
     hist_leading_pt_HLT_Jet30U_all_fine->Sumw2();
     hist_leading_pt_HLT_Jet50U_all_fine->Sumw2();
     hist_leading_pt_all_fine->Sumw2();

     hist_leading_central_pt_HLT_Jet15U_all_fine->Sumw2();
     hist_leading_central_pt_HLT_Jet30U_all_fine->Sumw2();
     hist_leading_central_pt_HLT_Jet50U_all_fine->Sumw2();
     hist_leading_central_pt_all_fine->Sumw2();

     hist_leading_eta_HLT_Jet15U_all_fine->Sumw2();
     hist_leading_eta_HLT_Jet30U_all_fine->Sumw2();
     hist_leading_eta_HLT_Jet50U_all_fine->Sumw2();
     hist_leading_eta_all_fine->Sumw2();
     hist_leading_eta_HLT_Jet30U_all_fine_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_all_fine_lowpt->Sumw2();
     hist_leading_eta_all_fine_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet30U_all_fine_highpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_all_fine_highpt->Sumw2();
     hist_leading_eta_all_fine_highpt->Sumw2();

     hist_leading_phi_HLT_Jet15U_all_fine->Sumw2();
     hist_leading_phi_HLT_Jet30U_all_fine->Sumw2();
     hist_leading_phi_HLT_Jet50U_all_fine->Sumw2();
     hist_leading_phi_all_fine->Sumw2();

//monitor distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_all_fine->Sumw2();
     hist_sj_leading_pt_HLT_Jet30U_all_fine->Sumw2();
     hist_sj_leading_pt_HLT_Jet50U_all_fine->Sumw2();
     hist_sj_leading_pt_all_fine->Sumw2();

     hist_sj_leading_eta_HLT_Jet15U_all_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_all_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_all_fine->Sumw2();
     hist_sj_leading_eta_all_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt->Sumw2();
     hist_sj_leading_eta_all_fine_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt->Sumw2();
     hist_sj_leading_eta_all_fine_highpt->Sumw2();

     hist_sj_leading_phi_HLT_Jet15U_all_fine->Sumw2();
     hist_sj_leading_phi_HLT_Jet30U_all_fine->Sumw2();
     hist_sj_leading_phi_HLT_Jet50U_all_fine->Sumw2();
     hist_sj_leading_phi_all_fine->Sumw2();

//triggered distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_emulated_fine->Sumw2();
     hist_leading_pt_HLT_Jet30U_emulated_fine->Sumw2();
     hist_leading_pt_HLT_Jet50U_emulated_fine->Sumw2();
     hist_leading_pt_emulated_fine->Sumw2();

     hist_leading_central_pt_HLT_Jet15U_emulated_fine->Sumw2();
     hist_leading_central_pt_HLT_Jet30U_emulated_fine->Sumw2();
     hist_leading_central_pt_HLT_Jet50U_emulated_fine->Sumw2();
     hist_leading_central_pt_emulated_fine->Sumw2();

     hist_leading_eta_HLT_Jet15U_emulated_fine->Sumw2();
     hist_leading_eta_HLT_Jet30U_emulated_fine->Sumw2();
     hist_leading_eta_HLT_Jet50U_emulated_fine->Sumw2();
     hist_leading_eta_emulated_fine->Sumw2();
     hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt->Sumw2();
     hist_leading_eta_emulated_fine_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet30U_emulated_fine_highpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_emulated_fine_highpt->Sumw2();
     hist_leading_eta_emulated_fine_highpt->Sumw2();

     hist_leading_phi_HLT_Jet15U_emulated_fine->Sumw2();
     hist_leading_phi_HLT_Jet30U_emulated_fine->Sumw2();
     hist_leading_phi_HLT_Jet50U_emulated_fine->Sumw2();
     hist_leading_phi_emulated_fine->Sumw2();

//triggered distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_emulated_fine->Sumw2();
     hist_sj_leading_pt_HLT_Jet30U_emulated_fine->Sumw2();
     hist_sj_leading_pt_HLT_Jet50U_emulated_fine->Sumw2();
     hist_sj_leading_pt_emulated_fine->Sumw2();

     hist_sj_leading_eta_HLT_Jet15U_emulated_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine->Sumw2();
     hist_sj_leading_eta_emulated_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt->Sumw2();
     hist_sj_leading_eta_emulated_fine_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt->Sumw2();
     hist_sj_leading_eta_emulated_fine_highpt->Sumw2();

     hist_sj_leading_phi_HLT_Jet15U_emulated_fine->Sumw2();
     hist_sj_leading_phi_HLT_Jet30U_emulated_fine->Sumw2();
     hist_sj_leading_phi_HLT_Jet50U_emulated_fine->Sumw2();
     hist_sj_leading_phi_emulated_fine->Sumw2();

//efficiency distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff_fine->Sumw2();
     hist_sj_leading_pt_HLT_Jet30U_eff_fine->Sumw2();
     hist_sj_leading_pt_HLT_Jet50U_eff_fine->Sumw2();
     hist_sj_leading_pt_eff_fine->Sumw2();

     hist_sj_leading_eta_HLT_Jet15U_eff_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_eff_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_eff_fine->Sumw2();
     hist_sj_leading_eta_eff_fine->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_lowpt->Sumw2();
     hist_sj_leading_eta_eff_fine_lowpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_highpt->Sumw2();
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_highpt->Sumw2();
     hist_sj_leading_eta_eff_fine_highpt->Sumw2();

     hist_sj_leading_phi_HLT_Jet15U_eff_fine->Sumw2();
     hist_sj_leading_phi_HLT_Jet30U_eff_fine->Sumw2();
     hist_sj_leading_phi_HLT_Jet50U_eff_fine->Sumw2();
     hist_sj_leading_phi_eff_fine->Sumw2();

//distributions for the observables
     hist_delta_phi_all->Sumw2();
     hist_delta_eta_all->Sumw2();
     hist_delta_phi_emulated->Sumw2();
     hist_delta_eta_emulated->Sumw2();
     hist_delta_phi_eff->Sumw2();
     hist_delta_eta_eff->Sumw2();

     hist_delta_phi_HLT_Jet15U_all->Sumw2();
     hist_delta_eta_HLT_Jet15U_all->Sumw2();
     hist_delta_phi_HLT_Jet15U_emulated->Sumw2();
     hist_delta_eta_HLT_Jet15U_emulated->Sumw2();
     hist_delta_phi_HLT_Jet15U_eff->Sumw2();
     hist_delta_eta_HLT_Jet15U_eff->Sumw2();

     hist_delta_phi_HLT_Jet30U_all->Sumw2();
     hist_delta_eta_HLT_Jet30U_all->Sumw2();
     hist_delta_phi_HLT_Jet30U_emulated->Sumw2();
     hist_delta_eta_HLT_Jet30U_emulated->Sumw2();
     hist_delta_phi_HLT_Jet30U_eff->Sumw2();
     hist_delta_eta_HLT_Jet30U_eff->Sumw2();

     hist_delta_phi_HLT_Jet50U_all->Sumw2();
     hist_delta_eta_HLT_Jet50U_all->Sumw2();
     hist_delta_phi_HLT_Jet50U_emulated->Sumw2();
     hist_delta_eta_HLT_Jet50U_emulated->Sumw2();
     hist_delta_phi_HLT_Jet50U_eff->Sumw2();
     hist_delta_eta_HLT_Jet50U_eff->Sumw2();

//loop over the files
for (int z = 0; z < n_files; z++)
{
//open the file
  string file = data_in[z];
  cout<<z+1<<"/"<<n_files<<" -> "<<file<<endl;
  TFile *inf = 0;
  inf = TFile::Open( file.c_str() );
  if (inf == 0) { cout << "Ntuple not loaded!" << endl; return; }

//define the tree branch
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);

   //------------------------- Initializing The Trigger Variables -------------------------- //
      cout << "Trigger list : " << endl;
      HLTJet[0] = "HLT_L1Jet6U";
      for (int i=0; i<=nJetTrig; i++)
         {
                triggered[i] = 0;
		l1pass[i] = 0;
		hltpass[i] = 0;
		condpass[i] = 0;
		selected[i] = 0;
		final[i] = 0;
         }

      for (int i=1; i<nJetTrig; i++)
         {
                sprintf(trigtitle, "HLT_Jet%iU",HLTJetPtN[i-1]);             
		HLTJet[i] = trigtitle;
		//cout << "Testing " << trigtitle << endl;
         }

      for (int i=0; i<=nJetTrig; i++)
         {
                  ihltj[i] = -1 ;
         }

   TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
   // ------------ Assigning An Integer To Each Trigger Path---------------//
    for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
       TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));

       for (int ii=0; ii<nJetTrig; ii++)
          {
                 if (ss == HLTJet[ii]) {
                 ihltj[ii] = ibin;
                 continue;
               }
         } // for (int ii=0; ii<nJetTrig; ii++)
      } // for(int ibin=0;ibin<55;ibin++)

  // ----------------------- Checking For The Trigger Assignment ------------------- // 
    for (int ij=0; ij<nJetTrig; ij++)
        {
              if (ihltj[ij] == -1) {
              cout<<"The requested trigger ("<<HLTJet[ij]<<") is not found "<<endl;
   //         break;
               }
              else {
              cout<<HLTJet[ij]<<" --> "<<ihltj[ij]<<endl;
             }
      }
      
  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  int decade = 0;

  if (test) { nentries = 1; } //reduced number of read entries, usefull for testing

   for (Int_t i=0;i<nentries;i++) {
     counter_entries++; 
     double progress = 10.0*i/(1.0*nentries);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
       cout<<k*10<<" % "<<endl;
     decade = k; 
     
     tr->GetEntry(i);
     
      //-------- check if the primary vertex is good ----  && Event->evtHdr().nVtxGood() == 1
      pv_pass = false;
      if (sel_mode == "allvertex" && Event->evtHdr().isPVgood() == 1) { pv_pass = true; }
      if (sel_mode == "1vertex"  && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)    	 { pv_pass = true; }
      
      if (pv_pass == true)
      {
        counter_pv++; 
        hltPass = false;
      
  //-------------------------- Initializing the Boolians and the Prescale Values----------------------//
        for (int j=0; j<nJetTrig; j++)
        {
                   hltPassj[j] = false ;
                   l1cut[j] = false;
                   hltcut[j] = false;
                   prescalej[j] = 1;
        } // for (int i=0; i<nJetTrig; i++)
	emu_l1 = false;
	emu_hlt = false;

  //----------------------- Computing The Prescale Values For a given event ----------------------- //
    for (int j=0; j<nJetTrig; j++)
         {
                 if (ihltj[j] == -1)
                 hltPassj[j] = true; // no trigger set
                 else {
                 if (Event->fired(ihltj[j]) > 0) {
                 hltPassj[j] = true;
		 triggered[j] = triggered[j] + 1;
                 hltPass = true;
                 prescalej[j] = Event->preL1(ihltj[j]) * Event->preHLT(ihltj[j]);
        //         if(i==0)
        //         cout<<HLTJet[i][j]<<" has prescale = "<<prescalej[i][j]<<endl;
                }
           }
      }
      
      
         //--------------------------- Filling Up The Histograms ------------------------------ //
      for (int j=0; j<nJetTrig-1; j++)
        {
            if (hltPassj[j])
             {
             if (test) { cout << "Event : " << i << " - Testing: " << HLTJet[j+1] << endl; }
               //----------------- L1 Theshold Checking --------------------------//
               for ( unsigned l1iobj=0; l1iobj<Event->nL1Obj(ihltj[j]); l1iobj++ )
                 {
                 counter_l1obj++;
                 if (test) { cout << "L1 object id: "<<l1iobj << " -> " << Event->l1obj(ihltj[j],l1iobj).pt() << " > " << ATrig[j] << endl; }
                   if (Event->l1obj(ihltj[j],l1iobj).pt() > ATrig[j])
                    {
                      l1cut[j] = true;
		      emu_l1 = true;
                    } // if (Event->l1obj(ihltj[i][j],l1iobj).pt() > ATrig[i+1])
                    
                   if(j==0)
                    {
                    if (Event->l1obj(ihltj[0],l1iobj).pt() > HLTJetPtN[0])
                     {
                     hltcut[0] = true;
                     emu_hlt = true;
                     }
                    }
                    
                  } // for ( unsigned l1iobj=0; l1iobj<Event->nL1Obj(ihltj[i][j]); l1iobj++ )

		if (l1cut[j] && test) { cout << "L1 passed!" << endl; }
		if (l1cut[j]) { l1pass[j] = l1pass[j] + 1; }

               //------------ HLT Threshold Checking For HLTObj ------------------ //
               for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[j]); hltiobj++ )
                 {
                 counter_hltobj++;
                 if (test) { cout << "HLT object id: " << hltiobj << " -> " << Event->hltobj(ihltj[j],hltiobj).pt() << " > " << HLTJetPtN[j] << endl; }
                  if (Event->hltobj(ihltj[j],hltiobj).pt() >  HLTJetPtN[j])
                   {
                     hltcut[j] = true;
		     emu_hlt = true;
                  } // if (Event->hltobj(ihltj[i][j],hltiobj).pt() >  HLTJetPtN[i+1])
                 
                } // for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[i][j]); hltiobj++ ) 
              
                if (hltcut[j] && test) { cout << "HLT passed!" << endl; }
		if (hltcut[j]) { hltpass[j] = hltpass[j] + 1; }
		if (l1cut[j] && hltcut[j]) { condpass[j] = condpass[j] + 1; }    

    	}
    
    }

	if (emu_l1) { l1pass[3] = l1pass[3] + 1; }
	if (emu_hlt) { hltpass[3] = hltpass[3] + 1; }
	if (emu_l1 && emu_hlt) { condpass[3] = condpass[3] + 1; }

     if (hltPass)
	{
	counter_hlt++;
	triggered[3] = triggered[3] + 1;

     	leading_pt = -10.0;
     	leading_eta = -10.0;
     	leading_phi = -10.0;
     	forward_pt = -10.0;
     	central_pt = -10.0;
     	forward_eta = -10.0;
     	central_eta = -10.0;
     	forward_phi = -10.0;
     	central_phi = -10.0;
	hard_pt = -10.0;
	hard_eta = -10.0;
	hard_phi = -10.0;
	trigger_factors[0] = 1.0;
	trigger_factors[1] = 1.0;
	trigger_factor = 1.0;
	delta_phi = -1.0;
	delta_eta = -1.0;
	fill_all = false;
	fill_emu = false;
	hard_emission = false;

    for(unsigned int j=0; j<Event->nPFJets(); j++) {
	pt = Event->pfjet(j).ptCor();
	eta = Event->pfjet(j).eta();
	phi = Event->pfjet(j).phi();
     if (pt >= pt_min && Event->pfjet(j).tightID() ) {
     counter_jet++;
     if (leading_pt < pt) { leading_pt = pt; leading_eta = eta; leading_phi = phi; }
     if (eta <= 2.8 && eta >= -2.8 && pt > central_pt)
     { central_pt = pt; central_eta = eta; central_phi = phi; }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > forward_pt)
     { forward_pt = pt; forward_eta = eta; forward_phi = phi; }
     }
     }


          if (leading_pt > pt_min)
     		{
     		//counter_selected++;
		//if (test) { cout << "selection in" << endl; }
		if (hltPassj[0])
			{
			hist_sj_leading_pt_HLT_Jet15U_all->Fill(leading_pt);
			hist_sj_leading_pt_HLT_Jet15U_all_fine->Fill(leading_pt);
			hist_sj_leading_phi_HLT_Jet15U_all->Fill(leading_phi);
			hist_sj_leading_phi_HLT_Jet15U_all_fine->Fill(leading_phi);
			hist_sj_leading_eta_HLT_Jet15U_all->Fill(leading_eta);
			hist_sj_leading_eta_HLT_Jet15U_all_fine->Fill(leading_eta);
			//selected[0] = selected[0] + 1;
			fill_all = true;
			//if (test) { cout << "HLT_Jet15 all: " << leading_pt << endl; }
			}
		if (hltPassj[1])
			{
			hist_sj_leading_pt_HLT_Jet30U_all->Fill(leading_pt);
			hist_sj_leading_pt_HLT_Jet30U_all_fine->Fill(leading_pt);
			hist_sj_leading_phi_HLT_Jet30U_all->Fill(leading_phi);
			hist_sj_leading_phi_HLT_Jet30U_all_fine->Fill(leading_phi);
			hist_sj_leading_eta_HLT_Jet30U_all->Fill(leading_eta);
			hist_sj_leading_eta_HLT_Jet30U_all_fine->Fill(leading_eta);
			if (leading_pt <= 70.0)
				{
				hist_sj_leading_eta_HLT_Jet30U_all_lowpt->Fill(leading_eta);
				hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt->Fill(leading_eta);
				}
				else
				{
				hist_sj_leading_eta_HLT_Jet30U_all_highpt->Fill(leading_eta);
				hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt->Fill(leading_eta);
				}
			//selected[1] = selected[1] + 1;
			fill_all = true;
			//if (test) { cout << "HLT_Jet30 all: " << leading_pt << endl; }
			}   
                if (hltPassj[2])
                	{
                	hist_sj_leading_pt_HLT_Jet50U_all->Fill(leading_pt);
			hist_sj_leading_pt_HLT_Jet50U_all_fine->Fill(leading_pt);
                	hist_sj_leading_phi_HLT_Jet50U_all->Fill(leading_phi);
			hist_sj_leading_phi_HLT_Jet50U_all_fine->Fill(leading_phi);
                	hist_sj_leading_eta_HLT_Jet50U_all->Fill(leading_eta);
			hist_sj_leading_eta_HLT_Jet50U_all_fine->Fill(leading_eta);
			if (leading_pt <= 110.0)
				{
                		hist_sj_leading_eta_HLT_Jet50U_all_lowpt->Fill(leading_eta);
				hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt->Fill(leading_eta);
				}
				else
				{
                		hist_sj_leading_eta_HLT_Jet50U_all_highpt->Fill(leading_eta);
				hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt->Fill(leading_eta);
				}
			//selected[2] = selected[2] + 1;
			fill_all = true;
                	//if (test) { cout << "HLT_Jet50 all: " << leading_pt << endl; }
                	}
                if (fill_all)
                	{
                	hist_sj_leading_pt_all->Fill(leading_pt);
			hist_sj_leading_pt_all_fine->Fill(leading_pt);
                	hist_sj_leading_phi_all->Fill(leading_phi);
			hist_sj_leading_phi_all_fine->Fill(leading_phi);
                	hist_sj_leading_eta_all->Fill(leading_eta);
			hist_sj_leading_eta_all_fine->Fill(leading_eta);
			if (leading_pt <= 110.0)
				{
                		hist_sj_leading_eta_all_lowpt->Fill(leading_eta);
				hist_sj_leading_eta_all_fine_lowpt->Fill(leading_eta);
				}
				else
				{
                		hist_sj_leading_eta_all_highpt->Fill(leading_eta);
				hist_sj_leading_eta_all_fine_highpt->Fill(leading_eta);
				}
			//selected[3] = selected[3] + 1;
                	//if (test) { cout << "Combination all: " << leading_pt << " and correction factor = " << trigger_factor << endl; }
                	} 
                if (hltPassj[0] && l1cut[0] && hltcut[0])
                	{
                	hist_sj_leading_pt_HLT_Jet15U_emulated->Fill(leading_pt);
			hist_sj_leading_pt_HLT_Jet15U_emulated_fine->Fill(leading_pt);
                	hist_sj_leading_phi_HLT_Jet15U_emulated->Fill(leading_phi);
			hist_sj_leading_phi_HLT_Jet15U_emulated_fine->Fill(leading_phi);
                	hist_sj_leading_eta_HLT_Jet15U_emulated->Fill(leading_eta);
			hist_sj_leading_eta_HLT_Jet15U_emulated_fine->Fill(leading_eta);
			//final[0] = final[0] + 1;
			fill_emu = true;
                	//if (test) { cout << "HLT_Jet15 emu" << endl; }
                	}   
                if (hltPassj[1] && l1cut[1] && hltcut[1])
                	{
			//trigger_factor = 1.0/correction_factor(p, 0, range, leading_pt, test);
                	hist_sj_leading_pt_HLT_Jet30U_emulated->Fill(leading_pt, trigger_factor);
			hist_sj_leading_pt_HLT_Jet30U_emulated_fine->Fill(leading_pt, trigger_factor);
                	hist_sj_leading_phi_HLT_Jet30U_emulated->Fill(leading_phi, trigger_factor);
			hist_sj_leading_phi_HLT_Jet30U_emulated_fine->Fill(leading_phi, trigger_factor);
                	hist_sj_leading_eta_HLT_Jet30U_emulated->Fill(leading_eta, trigger_factor);
			hist_sj_leading_eta_HLT_Jet30U_emulated_fine->Fill(leading_eta, trigger_factor);
			if (leading_pt <= 70.0)
				{
                		hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt->Fill(leading_eta, trigger_factor);
				hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt->Fill(leading_eta, trigger_factor);
				}
				else
				{
                		hist_sj_leading_eta_HLT_Jet30U_emulated_highpt->Fill(leading_eta, trigger_factor);
				hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt->Fill(leading_eta, trigger_factor);
				}
			//final[1] = final[1] + 1;
			fill_emu = true;
                	//if (test) { cout << "HLT_Jet30 emu" << endl; }
                	}
                if (hltPassj[2] && l1cut[2] && hltcut[2])
                	{
			//trigger_factor = 1.0/correction_factor(p, 1, range, leading_pt, test);
                	hist_sj_leading_pt_HLT_Jet50U_emulated->Fill(leading_pt, trigger_factor);
			hist_sj_leading_pt_HLT_Jet50U_emulated_fine->Fill(leading_pt, trigger_factor);
                	hist_sj_leading_phi_HLT_Jet50U_emulated->Fill(leading_phi, trigger_factor);
			hist_sj_leading_phi_HLT_Jet50U_emulated_fine->Fill(leading_phi, trigger_factor);
                	hist_sj_leading_eta_HLT_Jet50U_emulated->Fill(leading_eta, trigger_factor);
			hist_sj_leading_eta_HLT_Jet50U_emulated_fine->Fill(leading_eta, trigger_factor);
			if (leading_pt <= 110.0)
				{
                		hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt->Fill(leading_eta, trigger_factor);
				hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt->Fill(leading_eta, trigger_factor);
				}
				else
				{
                		hist_sj_leading_eta_HLT_Jet50U_emulated_highpt->Fill(leading_eta, trigger_factor);
				hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt->Fill(leading_eta, trigger_factor);
				}
			//final[2] = final[2] + 1;
			fill_emu = true;
                	//if (test) { cout << "HLT_Jet50 emu" << endl; }
                	}
                if (fill_emu)
                	{
			trigger_factor = 1.0;
			if (hltPassj[2] && l1cut[2] && hltcut[2])
                		{
				//trigger_factor = correction_factor(p, 0, range, leading_pt, test);
				}
			if (hltPassj[1] && l1cut[1] && hltcut[1])
                		{
				//trigger_factor = correction_factor(p, 1, range, leading_pt, test);
				}
			//trigger_factor = 1.0/trigger_factor;
                	hist_sj_leading_pt_emulated->Fill(leading_pt, trigger_factor);
			hist_sj_leading_pt_emulated_fine->Fill(leading_pt, trigger_factor);
                	hist_sj_leading_phi_emulated->Fill(leading_phi, trigger_factor);
			hist_sj_leading_phi_emulated_fine->Fill(leading_phi, trigger_factor);
                	hist_sj_leading_eta_emulated->Fill(leading_eta, trigger_factor);
			hist_sj_leading_eta_emulated_fine->Fill(leading_eta, trigger_factor);
			if (leading_pt <= 110.0)
				{
                		hist_sj_leading_eta_emulated_lowpt->Fill(leading_eta, trigger_factor);
				hist_sj_leading_eta_emulated_fine_lowpt->Fill(leading_eta, trigger_factor);
				}
				else
				{
                		hist_sj_leading_eta_emulated_highpt->Fill(leading_eta, trigger_factor);
				hist_sj_leading_eta_emulated_fine_highpt->Fill(leading_eta, trigger_factor);
				}
			//final[3] = final[3] + 1;
                	//if (test) { cout << "Combination emu" << endl; }
                	}

		}

		fill_all = false;
		fill_emu = false;
     
          if (forward_pt > pt_min && central_pt > pt_min)
     		{
     		counter_selected++;
		delta_phi = calc_delta_phi(forward_phi, central_phi);
		delta_eta = forward_eta - central_eta;
     		if (delta_eta < 0) { delta_eta = -delta_eta; }
		if (test) { cout << "selection in" << endl; }
		if (hltPassj[0])
			{
			hist_leading_pt_HLT_Jet15U_all->Fill(leading_pt);
			hist_leading_pt_HLT_Jet15U_all_fine->Fill(leading_pt);
			hist_leading_central_pt_HLT_Jet15U_all->Fill(central_pt);
			hist_leading_central_pt_HLT_Jet15U_all_fine->Fill(central_pt);
			hist_leading_phi_HLT_Jet15U_all->Fill(leading_phi);
			hist_leading_phi_HLT_Jet15U_all_fine->Fill(leading_phi);
			hist_leading_eta_HLT_Jet15U_all->Fill(leading_eta);
			hist_leading_eta_HLT_Jet15U_all_fine->Fill(leading_eta);
			hist_delta_phi_HLT_Jet15U_all->Fill(delta_phi);
			hist_delta_eta_HLT_Jet15U_all->Fill(delta_eta);
			selected[0] = selected[0] + 1;
			fill_all = true;
			if (test) { cout << "HLT_Jet15 all: " << leading_pt << endl; }
			}
		if (hltPassj[1])
			{
			hist_leading_pt_HLT_Jet30U_all->Fill(leading_pt);
			hist_leading_pt_HLT_Jet30U_all_fine->Fill(leading_pt);
			hist_leading_central_pt_HLT_Jet30U_all->Fill(central_pt);
			hist_leading_central_pt_HLT_Jet30U_all_fine->Fill(central_pt);
			hist_leading_phi_HLT_Jet30U_all->Fill(leading_phi);
			hist_leading_phi_HLT_Jet30U_all_fine->Fill(leading_phi);
			hist_leading_eta_HLT_Jet30U_all->Fill(leading_eta);
			hist_leading_eta_HLT_Jet30U_all_fine->Fill(leading_eta);
			hist_delta_phi_HLT_Jet30U_all->Fill(delta_phi);
			hist_delta_eta_HLT_Jet30U_all->Fill(delta_eta);
			if (leading_pt <= 70.0)
				{
				hist_leading_eta_HLT_Jet30U_all_lowpt->Fill(leading_eta);
				hist_leading_eta_HLT_Jet30U_all_fine_lowpt->Fill(leading_eta);
				}
				else
				{
				hist_leading_eta_HLT_Jet30U_all_highpt->Fill(leading_eta);
				hist_leading_eta_HLT_Jet30U_all_fine_highpt->Fill(leading_eta);
				}
			selected[1] = selected[1] + 1;
			fill_all = true;
			if (test) { cout << "HLT_Jet30 all: " << leading_pt << endl; }
			}   
                if (hltPassj[2])
                	{
                	hist_leading_pt_HLT_Jet50U_all->Fill(leading_pt);
			hist_leading_pt_HLT_Jet50U_all_fine->Fill(leading_pt);
                	hist_leading_central_pt_HLT_Jet50U_all->Fill(central_pt);
			hist_leading_central_pt_HLT_Jet50U_all_fine->Fill(central_pt);
                	hist_leading_phi_HLT_Jet50U_all->Fill(leading_phi);
			hist_leading_phi_HLT_Jet50U_all_fine->Fill(leading_phi);
                	hist_leading_eta_HLT_Jet50U_all->Fill(leading_eta);
			hist_leading_eta_HLT_Jet50U_all_fine->Fill(leading_eta);
			hist_delta_phi_HLT_Jet50U_all->Fill(delta_phi);
			hist_delta_eta_HLT_Jet50U_all->Fill(delta_eta);
			if (leading_pt <= 110.0)
				{
                		hist_leading_eta_HLT_Jet50U_all_lowpt->Fill(leading_eta);
				hist_leading_eta_HLT_Jet50U_all_fine_lowpt->Fill(leading_eta);
				}
				else
				{
                		hist_leading_eta_HLT_Jet50U_all_highpt->Fill(leading_eta);
				hist_leading_eta_HLT_Jet50U_all_fine_highpt->Fill(leading_eta);
				}
			selected[2] = selected[2] + 1;
			fill_all = true;
                	if (test) { cout << "HLT_Jet50 all: " << leading_pt << endl; }
                	}
                if (fill_all)
                	{
                	hist_leading_pt_all->Fill(leading_pt);
			hist_leading_pt_all_fine->Fill(leading_pt);
                	hist_leading_central_pt_all->Fill(central_pt);
			hist_leading_central_pt_all_fine->Fill(central_pt);
                	hist_leading_phi_all->Fill(leading_phi);
			hist_leading_phi_all_fine->Fill(leading_phi);
                	hist_leading_eta_all->Fill(leading_eta);
			hist_leading_eta_all_fine->Fill(leading_eta);
			hist_delta_phi_all->Fill(delta_phi);
			hist_delta_eta_all->Fill(delta_eta);
			if (leading_pt <= 110.0)
				{
                		hist_leading_eta_all_lowpt->Fill(leading_eta);
				hist_leading_eta_all_fine_lowpt->Fill(leading_eta);
				}
				else
				{
                		hist_leading_eta_all_highpt->Fill(leading_eta);
				hist_leading_eta_all_fine_highpt->Fill(leading_eta);
				}
			selected[3] = selected[3] + 1;
                	if (test) { cout << "Combination all: " << leading_pt << " and correction factor = " << trigger_factor << endl; }
                	} 
                if (hltPassj[0] && l1cut[0])
                	{
                	hist_leading_pt_HLT_Jet15U_emulated->Fill(leading_pt);
			hist_leading_pt_HLT_Jet15U_emulated_fine->Fill(leading_pt);
                	hist_leading_central_pt_HLT_Jet15U_emulated->Fill(central_pt);
			hist_leading_central_pt_HLT_Jet15U_emulated_fine->Fill(central_pt);
                	hist_leading_phi_HLT_Jet15U_emulated->Fill(leading_phi);
			hist_leading_phi_HLT_Jet15U_emulated_fine->Fill(leading_phi);
                	hist_leading_eta_HLT_Jet15U_emulated->Fill(leading_eta);
			hist_leading_eta_HLT_Jet15U_emulated_fine->Fill(leading_eta);
			hist_delta_phi_HLT_Jet15U_emulated->Fill(delta_phi);
			hist_delta_eta_HLT_Jet15U_emulated->Fill(delta_eta);
			final[0] = final[0] + 1;
			fill_emu = true;
                	if (test) { cout << "HLT_Jet15 emu" << endl; }
                	} 
		if (corr) { correction_factor(trigger_hist1, trigger_hist2, central_pt, trigger_factors); }  
                if (hltPassj[1] && l1cut[1] && hltcut[1])
                	{
			trigger_factor = trigger_factors[0];
                	hist_leading_pt_HLT_Jet30U_emulated->Fill(leading_pt, trigger_factor);
			hist_leading_pt_HLT_Jet30U_emulated_fine->Fill(leading_pt, trigger_factor);
               		hist_leading_central_pt_HLT_Jet30U_emulated->Fill(central_pt, trigger_factor);
			hist_leading_central_pt_HLT_Jet30U_emulated_fine->Fill(central_pt, trigger_factor);
                	hist_leading_phi_HLT_Jet30U_emulated->Fill(leading_phi, trigger_factor);
			hist_leading_phi_HLT_Jet30U_emulated_fine->Fill(leading_phi, trigger_factor);
                	hist_leading_eta_HLT_Jet30U_emulated->Fill(leading_eta, trigger_factor);
			hist_leading_eta_HLT_Jet30U_emulated->Fill(leading_eta, trigger_factor);
                	hist_delta_phi_HLT_Jet30U_emulated->Fill(delta_phi, trigger_factor);
			hist_delta_eta_HLT_Jet30U_emulated->Fill(delta_eta, trigger_factor);
			if (leading_pt <= 70.0)
				{
                		hist_leading_eta_HLT_Jet30U_emulated_lowpt->Fill(leading_eta, trigger_factor);
				hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt->Fill(leading_eta, trigger_factor);
				}
				else
				{
                		hist_leading_eta_HLT_Jet30U_emulated_highpt->Fill(leading_eta, trigger_factor);
				hist_leading_eta_HLT_Jet30U_emulated_fine_highpt->Fill(leading_eta, trigger_factor);
				}
			final[1] = final[1] + 1;
			fill_emu = true;
                	if (test) { cout << "HLT_Jet30 emu" << endl; }
                	}
                if (hltPassj[2] && l1cut[2] && hltcut[2])
                	{
			trigger_factor = trigger_factors[1];
                	hist_leading_pt_HLT_Jet50U_emulated->Fill(leading_pt, trigger_factor);
			hist_leading_pt_HLT_Jet50U_emulated_fine->Fill(leading_pt, trigger_factor);
                	hist_leading_central_pt_HLT_Jet50U_emulated->Fill(central_pt, trigger_factor);
			hist_leading_central_pt_HLT_Jet50U_emulated_fine->Fill(central_pt, trigger_factor);
                	hist_leading_phi_HLT_Jet50U_emulated->Fill(leading_phi, trigger_factor);
			hist_leading_phi_HLT_Jet50U_emulated_fine->Fill(leading_phi, trigger_factor);
                	hist_leading_eta_HLT_Jet50U_emulated->Fill(leading_eta, trigger_factor);
			hist_leading_eta_HLT_Jet50U_emulated_fine->Fill(leading_eta, trigger_factor);
                	hist_delta_phi_HLT_Jet50U_emulated->Fill(delta_phi, trigger_factor);
			hist_delta_eta_HLT_Jet50U_emulated->Fill(delta_eta, trigger_factor);
			if (leading_pt <= 110.0)
				{
                		hist_leading_eta_HLT_Jet50U_emulated_lowpt->Fill(leading_eta, trigger_factor);
				hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt->Fill(leading_eta, trigger_factor);
				}
				else
				{
                		hist_leading_eta_HLT_Jet50U_emulated_highpt->Fill(leading_eta, trigger_factor);
				hist_leading_eta_HLT_Jet50U_emulated_fine_highpt->Fill(leading_eta, trigger_factor);
				}
			final[2] = final[2] + 1;
			fill_emu = true;
                	if (test) { cout << "HLT_Jet50 emu" << endl; }
                	}
                if (fill_emu)
                	{
			trigger_factor = 1.0;
			if (hltPassj[1] && l1cut[1] && hltcut[1] && !hltPassj[2] && !l1cut[2] && !hltcut[2])
                		{
				trigger_factor = trigger_factors[0];
				}
			if (!hltPassj[1] && !l1cut[1] && !hltcut[1] && hltPassj[2] && l1cut[2] && hltcut[2])
                		{
				trigger_factor = trigger_factors[1];
				}
			if (hltPassj[1] && l1cut[1] && hltcut[1] && hltPassj[2] && l1cut[2] && hltcut[2])
                		{
				trigger_factor = trigger_factors[0]/2.0 + trigger_factors[1]/2.0;
				}
                	hist_leading_pt_emulated->Fill(leading_pt,trigger_factor);
			hist_leading_pt_emulated_fine->Fill(leading_pt,trigger_factor);
                	hist_leading_central_pt_emulated->Fill(central_pt,trigger_factor);
			hist_leading_central_pt_emulated_fine->Fill(central_pt,trigger_factor);
                	hist_leading_phi_emulated->Fill(leading_phi,trigger_factor);
			hist_leading_phi_emulated_fine->Fill(leading_phi,trigger_factor);
                	hist_leading_eta_emulated->Fill(leading_eta,trigger_factor);
			hist_leading_eta_emulated_fine->Fill(leading_eta,trigger_factor);
                	hist_delta_phi_emulated->Fill(delta_phi,trigger_factor);
			hist_delta_eta_emulated->Fill(delta_eta,trigger_factor);
			if (leading_pt <= 110.0)
				{
                		hist_leading_eta_emulated_lowpt->Fill(leading_eta,trigger_factor);
				hist_leading_eta_emulated_fine_lowpt->Fill(leading_eta,trigger_factor);
				}
				else
				{
                		hist_leading_eta_emulated_highpt->Fill(leading_eta,trigger_factor);
				hist_leading_eta_emulated_fine_highpt->Fill(leading_eta,trigger_factor);
				}
			final[3] = final[3] + 1;
                	if (test) { cout << "Combination emu" << endl; }
                	}

		}
	}

}

}

}

//cleaning the correction histogram
     if (corr)
	{
	delete(trigger_hist1);
	delete(trigger_hist2);
	}

//computing auxiliary outputs
     pu_scale = (double)counter_entries/(double)counter_pv;
     eff = (double)counter_hlt/(double)counter_pv;

//output a summary
     if (detail) { cout<<"Total Number of Events :         "<<counter_entries<<endl; }
     if (detail) { cout<<"Primary Vertex Filter :          "<<counter_pv<<endl; }
     if (detail) { cout<<"Pileup Scale :                   "<<pu_scale<<endl; }
     if (detail) { cout<<"Triggered Events :               "<<counter_hlt<<endl; }
     if (detail) { cout<<"Trigger Efficiency :             "<<eff<<endl; }
     if (detail) { cout<<"Selected :                       "<<counter_selected<<endl; }
     if (detail) { cout<<"Total L1 Objects :               "<<counter_l1obj<<endl;}
     if (detail) { cout<<"Total HLT Objects :              "<<counter_hltobj<<endl;}
     if (detail) { cout<<"Total Jets :                     "<<counter_jet<<endl; }
     
     //fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_hlt);
     hist_events->SetBinContent(4,counter_selected);
     hist_events->SetBinContent(5,counter_l1obj);
     hist_events->SetBinContent(6,counter_hltobj);
     hist_events->SetBinContent(7,counter_jet);
     hist_events->SetBinContent(8,pu_scale);
     hist_events->SetBinContent(9,eff);

     //output a summary for HLT_Jet15U
     if (detail) { cout<<"Results for HLT_Jet15U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_L1Jet6U: "<<triggered[0]<<endl; }
     if (detail) { cout<<"Passed the L1 Selection :              "<<l1pass[0]<<endl; }
     if (detail) { cout<<"Passed the HLT Selection :             "<<hltpass[0]<<endl; }
     if (detail) { cout<<"Passed the emulation condition :       "<<condpass[0]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[0]<<endl; }
     if (detail) { cout<<"Selected & Emulation:                  "<<final[0]<<endl; }

     //fill the events in HLT_Jet15U histogram
     hist_events_HLT_Jet15U->SetBinContent(1,triggered[0]);
     hist_events_HLT_Jet15U->SetBinContent(2,l1pass[0]);
     hist_events_HLT_Jet15U->SetBinContent(3,hltpass[0]);
     hist_events_HLT_Jet15U->SetBinContent(4,condpass[0]);
     hist_events_HLT_Jet15U->SetBinContent(5,selected[0]);
     hist_events_HLT_Jet15U->SetBinContent(6,final[0]);

     //output a summary for HLT_Jet30U
     if (detail) { cout<<"Results for HLT_Jet30U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet15U:  "<<triggered[1]<<endl; }
     if (detail) { cout<<"Passed the L1 Selection :              "<<l1pass[1]<<endl; }
     if (detail) { cout<<"Passed the HLT Selection :             "<<hltpass[1]<<endl; }
     if (detail) { cout<<"Passed the emulation condition :       "<<condpass[1]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[1]<<endl; }
     if (detail) { cout<<"Selected & Emulation:                  "<<final[1]<<endl; }

     //fill the events in HLT_Jet30U histogram
     hist_events_HLT_Jet30U->SetBinContent(1,triggered[1]);
     hist_events_HLT_Jet30U->SetBinContent(2,l1pass[1]);
     hist_events_HLT_Jet30U->SetBinContent(3,hltpass[1]);
     hist_events_HLT_Jet30U->SetBinContent(4,condpass[1]);
     hist_events_HLT_Jet30U->SetBinContent(5,selected[1]);
     hist_events_HLT_Jet30U->SetBinContent(6,final[1]);

     //output a summary for HLT_Jet50U
     if (detail) { cout<<"Results for HLT_Jet50U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet30U:  "<<triggered[2]<<endl; }
     if (detail) { cout<<"Passed the L1 Selection :              "<<l1pass[2]<<endl; }
     if (detail) { cout<<"Passed the HLT Selection :             "<<hltpass[2]<<endl; }
     if (detail) { cout<<"Passed the emulation condition :       "<<condpass[2]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[2]<<endl; }
     if (detail) { cout<<"Selected & Emulation:                  "<<final[2]<<endl; }

     //fill the events in HLT_Jet50U histogram
     hist_events_HLT_Jet50U->SetBinContent(1,triggered[2]);
     hist_events_HLT_Jet50U->SetBinContent(2,l1pass[2]);
     hist_events_HLT_Jet50U->SetBinContent(3,hltpass[2]);
     hist_events_HLT_Jet50U->SetBinContent(4,condpass[2]);
     hist_events_HLT_Jet50U->SetBinContent(5,selected[2]);
     hist_events_HLT_Jet50U->SetBinContent(6,final[2]);

     //output a summary for combination
     if (detail) { cout<<"Results for combination" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[3]<<endl; }
     if (detail) { cout<<"Passed the L1 Selection :              "<<l1pass[3]<<endl; }
     if (detail) { cout<<"Passed the HLT Selection :             "<<hltpass[3]<<endl; }
     if (detail) { cout<<"Passed the emulation condition :       "<<condpass[3]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[3]<<endl; }
     if (detail) { cout<<"Selected & Emulation:                  "<<final[3]<<endl; }

     //fill the events in combination histogram
     hist_events_all->SetBinContent(1,triggered[3]);
     hist_events_all->SetBinContent(2,l1pass[3]);
     hist_events_all->SetBinContent(3,hltpass[3]);
     hist_events_all->SetBinContent(4,condpass[3]);
     hist_events_all->SetBinContent(5,selected[3]);
     hist_events_all->SetBinContent(6,final[3]);

     // creating the efficiency histograms
     if (detail) { cout<<"Creating the efficiency histograms..."<<endl; }

//efficiency distributions with fine binning for dijets
     TH1D *hist_leading_pt_HLT_Jet15U_eff_fine;
     TH1D *hist_leading_pt_HLT_Jet30U_eff_fine;
     TH1D *hist_leading_pt_HLT_Jet50U_eff_fine;
     TH1D *hist_leading_pt_eff_fine;

     TH1D *hist_leading_central_pt_HLT_Jet15U_eff_fine;
     TH1D *hist_leading_central_pt_HLT_Jet30U_eff_fine;
     TH1D *hist_leading_central_pt_HLT_Jet50U_eff_fine;
     TH1D *hist_leading_central_pt_eff_fine;

     TH1D *hist_leading_eta_HLT_Jet15U_eff_fine;
     TH1D *hist_leading_eta_HLT_Jet30U_eff_fine;
     TH1D *hist_leading_eta_HLT_Jet50U_eff_fine;
     TH1D *hist_leading_eta_eff_fine;
     TH1D *hist_leading_eta_HLT_Jet30U_eff_fine_lowpt;
     TH1D *hist_leading_eta_HLT_Jet50U_eff_fine_lowpt;
     TH1D *hist_leading_eta_eff_fine_lowpt;
     TH1D *hist_leading_eta_HLT_Jet30U_eff_fine_highpt;
     TH1D *hist_leading_eta_HLT_Jet50U_eff_fine_highpt;
     TH1D *hist_leading_eta_eff_fine_highpt;

     TH1D *hist_leading_phi_HLT_Jet15U_eff_fine;
     TH1D *hist_leading_phi_HLT_Jet30U_eff_fine;
     TH1D *hist_leading_phi_HLT_Jet50U_eff_fine;
     TH1D *hist_leading_phi_eff_fine;

//efficiency distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_leading_pt_HLT15U_eff_fine","Leading Jet p_{T} for HLT_Jet15U efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_leading_pt_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_leading_pt_HLT30U_eff_fine","Leading Jet p_{T} for HLT_Jet30U efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_leading_pt_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_leading_pt_HLT50U_eff_fine","Leading Jet p_{T} for HLT_Jet50U efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_leading_pt_eff_fine =  new TH1D("ak5PF_leading_pt_eff_fine","Leading Jet p_{T} for combination efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);

     hist_leading_central_pt_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_leading_central_pt_HLT15U_eff_fine","Leading Central Jet p_{T} for HLT_Jet15U efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_leading_central_pt_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_leading_central_pt_HLT30U_eff_fine","Leading Central Jet p_{T} for HLT_Jet30U efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_leading_central_pt_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_leading_central_pt_HLT50U_eff_fine","Leading Central Jet p_{T} for HLT_Jet50U efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);
     hist_leading_central_pt_eff_fine =  new TH1D("ak5PF_leading_central_pt_eff_fine","Leading Central Jet p_{T} for combination efficiency fine;p_{T} [#frac{GeV}{c}];efficiency", 165, 35, 200);

     hist_leading_eta_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_leading_eta_HLT15U_eff_fine","Leading Jet #eta for HLT_Jet15U efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_leading_eta_HLT30U_eff_fine","Leading Jet #eta for HLT_Jet30U efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_leading_eta_HLT50U_eff_fine","Leading Jet #eta for HLT_Jet50U efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_eff_fine =  new TH1D("ak5PF_leading_eta_eff_fine","Leading Jet #eta for combination efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_eff_fine_lowpt =  new TH1D("ak5PF_leading_eta_HLT30U_eff_fine_lowpt","Leading Jet #eta for HLT_Jet30U low pt efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_eff_fine_lowpt =  new TH1D("ak5PF_leading_eta_HLT50U_eff_fine_lowpt","Leading Jet #eta for HLT_Jet50U low pt efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_eff_fine_lowpt =  new TH1D("ak5PF_leading_eta_eff_fine_lowpt","Leading Jet #eta for combination low pt efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet30U_eff_fine_highpt =  new TH1D("ak5PF_leading_eta_HLT30U_eff_fine_highpt","Leading Jet #eta for HLT_Jet30U high pt efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_HLT_Jet50U_eff_fine_highpt =  new TH1D("ak5PF_leading_eta_HLT50U_eff_fine_highpt","Leading Jet #eta for HLT_Jet50U high pt efficiency fine;#eta;efficiency", 50, -5.0, 5.0);
     hist_leading_eta_eff_fine_highpt =  new TH1D("ak5PF_leading_eta_eff_fine_highpt","Leading Jet #eta for combination high pt efficiency fine;#eta;efficiency", 50, -5.0, 5.0);

     hist_leading_phi_HLT_Jet15U_eff_fine =  new TH1D("ak5PF_leading_phi_HLT15U_eff_fine","Leading Jet #phi for HLT_Jet15U efficiency fine;#phi;efficiency", 64, -3.15, 3.15);
     hist_leading_phi_HLT_Jet30U_eff_fine =  new TH1D("ak5PF_leading_phi_HLT30U_eff_fine","Leading Jet #phi for HLT_Jet30U efficiency fine;#phi;efficiency", 64, -3.15, 3.15);
     hist_leading_phi_HLT_Jet50U_eff_fine =  new TH1D("ak5PF_leading_phi_HLT50U_eff_fine","Leading Jet #phi for HLT_Jet50U efficiency fine;#phi;efficiency", 64, -3.15, 3.15);
     hist_leading_phi_eff_fine =  new TH1D("ak5PF_leading_phi_eff_fine","Leading Jet #phi for combination efficiency fine;#phi;efficiency", 64, -3.15, 3.15);

//efficiency distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_eff_fine->Sumw2();
     hist_leading_pt_HLT_Jet30U_eff_fine->Sumw2();
     hist_leading_pt_HLT_Jet50U_eff_fine->Sumw2();
     hist_leading_pt_eff_fine->Sumw2();

     hist_leading_central_pt_HLT_Jet15U_eff_fine->Sumw2();
     hist_leading_central_pt_HLT_Jet30U_eff_fine->Sumw2();
     hist_leading_central_pt_HLT_Jet50U_eff_fine->Sumw2();
     hist_leading_central_pt_eff_fine->Sumw2();

     hist_leading_eta_HLT_Jet15U_eff_fine->Sumw2();
     hist_leading_eta_HLT_Jet30U_eff_fine->Sumw2();
     hist_leading_eta_HLT_Jet50U_eff_fine->Sumw2();
     hist_leading_eta_eff_fine->Sumw2();
     hist_leading_eta_HLT_Jet30U_eff_fine_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_eff_fine_lowpt->Sumw2();
     hist_leading_eta_eff_fine_lowpt->Sumw2();
     hist_leading_eta_HLT_Jet30U_eff_fine_highpt->Sumw2();
     hist_leading_eta_HLT_Jet50U_eff_fine_highpt->Sumw2();
     hist_leading_eta_eff_fine_highpt->Sumw2();

     hist_leading_phi_HLT_Jet15U_eff_fine->Sumw2();
     hist_leading_phi_HLT_Jet30U_eff_fine->Sumw2();
     hist_leading_phi_HLT_Jet50U_eff_fine->Sumw2();
     hist_leading_phi_eff_fine->Sumw2();


     // broad binning histograms for dijets
     hist_leading_pt_HLT_Jet15U_eff->Divide(hist_leading_pt_HLT_Jet15U_emulated,hist_leading_pt_HLT_Jet15U_all,1.,1.,"B");
     hist_leading_pt_HLT_Jet30U_eff->Divide(hist_leading_pt_HLT_Jet30U_emulated,hist_leading_pt_HLT_Jet30U_all,1.,1.,"B");
     hist_leading_pt_HLT_Jet50U_eff->Divide(hist_leading_pt_HLT_Jet50U_emulated,hist_leading_pt_HLT_Jet50U_all,1.,1.,"B");
     hist_leading_pt_eff->Divide(hist_leading_pt_emulated,hist_leading_pt_all,1.,1.,"B");

     hist_leading_central_pt_HLT_Jet15U_eff->Divide(hist_leading_central_pt_HLT_Jet15U_emulated,hist_leading_central_pt_HLT_Jet15U_all,1.,1.,"B");
     hist_leading_central_pt_HLT_Jet30U_eff->Divide(hist_leading_central_pt_HLT_Jet30U_emulated,hist_leading_central_pt_HLT_Jet30U_all,1.,1.,"B");
     hist_leading_central_pt_HLT_Jet50U_eff->Divide(hist_leading_central_pt_HLT_Jet50U_emulated,hist_leading_central_pt_HLT_Jet50U_all,1.,1.,"B");
     hist_leading_central_pt_eff->Divide(hist_leading_central_pt_emulated,hist_leading_central_pt_all,1.,1.,"B");

     hist_leading_eta_HLT_Jet15U_eff->Divide(hist_leading_eta_HLT_Jet15U_emulated,hist_leading_eta_HLT_Jet15U_all,1.,1.,"B");
     hist_leading_eta_HLT_Jet30U_eff->Divide(hist_leading_eta_HLT_Jet30U_emulated,hist_leading_eta_HLT_Jet30U_all,1.,1.,"B");
     hist_leading_eta_HLT_Jet50U_eff->Divide(hist_leading_eta_HLT_Jet50U_emulated,hist_leading_eta_HLT_Jet50U_all,1.,1.,"B");
     hist_leading_eta_eff->Divide(hist_leading_eta_emulated,hist_leading_eta_all,1.,1.,"B");
     hist_leading_eta_HLT_Jet30U_eff_lowpt->Divide(hist_leading_eta_HLT_Jet30U_emulated_lowpt,hist_leading_eta_HLT_Jet30U_all_lowpt,1.,1.,"B");
     hist_leading_eta_HLT_Jet50U_eff_lowpt->Divide(hist_leading_eta_HLT_Jet50U_emulated_lowpt,hist_leading_eta_HLT_Jet50U_all_lowpt,1.,1.,"B");
     hist_leading_eta_eff_lowpt->Divide(hist_leading_eta_emulated_lowpt,hist_leading_eta_all_lowpt,1.,1.,"B");
     hist_leading_eta_HLT_Jet30U_eff_highpt->Divide(hist_leading_eta_HLT_Jet30U_emulated_highpt,hist_leading_eta_HLT_Jet30U_all_highpt,1.,1.,"B");
     hist_leading_eta_HLT_Jet50U_eff_highpt->Divide(hist_leading_eta_HLT_Jet50U_emulated_highpt,hist_leading_eta_HLT_Jet50U_all_highpt,1.,1.,"B");
     hist_leading_eta_eff_highpt->Divide(hist_leading_eta_emulated_highpt,hist_leading_eta_all_highpt,1.,1.,"B");

     hist_leading_phi_HLT_Jet15U_eff->Divide(hist_leading_phi_HLT_Jet15U_emulated,hist_leading_phi_HLT_Jet15U_all,1.,1.,"B");
     hist_leading_phi_HLT_Jet30U_eff->Divide(hist_leading_phi_HLT_Jet30U_emulated,hist_leading_phi_HLT_Jet30U_all,1.,1.,"B");
     hist_leading_phi_HLT_Jet50U_eff->Divide(hist_leading_phi_HLT_Jet50U_emulated,hist_leading_phi_HLT_Jet50U_all,1.,1.,"B");
     hist_leading_phi_eff->Divide(hist_leading_phi_emulated,hist_leading_phi_all,1.,1.,"B");

     // broad binning histograms for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff->Divide(hist_sj_leading_pt_HLT_Jet15U_emulated,hist_sj_leading_pt_HLT_Jet15U_all,1.,1.,"B");
     hist_sj_leading_pt_HLT_Jet30U_eff->Divide(hist_sj_leading_pt_HLT_Jet30U_emulated,hist_sj_leading_pt_HLT_Jet30U_all,1.,1.,"B");
     hist_sj_leading_pt_HLT_Jet50U_eff->Divide(hist_sj_leading_pt_HLT_Jet50U_emulated,hist_sj_leading_pt_HLT_Jet50U_all,1.,1.,"B");
     hist_sj_leading_pt_eff->Divide(hist_sj_leading_pt_emulated,hist_sj_leading_pt_all,1.,1.,"B");

     hist_sj_leading_eta_HLT_Jet15U_eff->Divide(hist_sj_leading_eta_HLT_Jet15U_emulated,hist_sj_leading_eta_HLT_Jet15U_all,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet30U_eff->Divide(hist_sj_leading_eta_HLT_Jet30U_emulated,hist_sj_leading_eta_HLT_Jet30U_all,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet50U_eff->Divide(hist_sj_leading_eta_HLT_Jet50U_emulated,hist_sj_leading_eta_HLT_Jet50U_all,1.,1.,"B");
     hist_sj_leading_eta_eff->Divide(hist_sj_leading_eta_emulated,hist_sj_leading_eta_all,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet30U_eff_lowpt->Divide(hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt,hist_sj_leading_eta_HLT_Jet30U_all_lowpt,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet50U_eff_lowpt->Divide(hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt,hist_sj_leading_eta_HLT_Jet50U_all_lowpt,1.,1.,"B");
     hist_sj_leading_eta_eff_lowpt->Divide(hist_sj_leading_eta_emulated_lowpt,hist_sj_leading_eta_all_lowpt,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet30U_eff_highpt->Divide(hist_sj_leading_eta_HLT_Jet30U_emulated_highpt,hist_sj_leading_eta_HLT_Jet30U_all_highpt,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet50U_eff_highpt->Divide(hist_sj_leading_eta_HLT_Jet50U_emulated_highpt,hist_sj_leading_eta_HLT_Jet50U_all_highpt,1.,1.,"B");
     hist_sj_leading_eta_eff_highpt->Divide(hist_sj_leading_eta_emulated_highpt,hist_sj_leading_eta_all_highpt,1.,1.,"B");

     hist_sj_leading_phi_HLT_Jet15U_eff->Divide(hist_sj_leading_phi_HLT_Jet15U_emulated,hist_sj_leading_phi_HLT_Jet15U_all,1.,1.,"B");
     hist_sj_leading_phi_HLT_Jet30U_eff->Divide(hist_sj_leading_phi_HLT_Jet30U_emulated,hist_sj_leading_phi_HLT_Jet30U_all,1.,1.,"B");
     hist_sj_leading_phi_HLT_Jet50U_eff->Divide(hist_sj_leading_phi_HLT_Jet50U_emulated,hist_sj_leading_phi_HLT_Jet50U_all,1.,1.,"B");
     hist_sj_leading_phi_eff->Divide(hist_sj_leading_phi_emulated,hist_sj_leading_phi_all,1.,1.,"B");

     // fine binning histograms for dijets
     hist_leading_pt_HLT_Jet15U_eff_fine->Divide(hist_leading_pt_HLT_Jet15U_emulated_fine,hist_leading_pt_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_leading_pt_HLT_Jet30U_eff_fine->Divide(hist_leading_pt_HLT_Jet30U_emulated_fine,hist_leading_pt_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_leading_pt_HLT_Jet50U_eff_fine->Divide(hist_leading_pt_HLT_Jet50U_emulated_fine,hist_leading_pt_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_leading_pt_eff_fine->Divide(hist_leading_pt_emulated_fine,hist_leading_pt_all_fine,1.,1.,"B");

     hist_leading_central_pt_HLT_Jet15U_eff_fine->Divide(hist_leading_central_pt_HLT_Jet15U_emulated_fine,hist_leading_central_pt_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_leading_central_pt_HLT_Jet30U_eff_fine->Divide(hist_leading_central_pt_HLT_Jet30U_emulated_fine,hist_leading_central_pt_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_leading_central_pt_HLT_Jet50U_eff_fine->Divide(hist_leading_central_pt_HLT_Jet50U_emulated_fine,hist_leading_central_pt_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_leading_central_pt_eff_fine->Divide(hist_leading_central_pt_emulated_fine,hist_leading_central_pt_all_fine,1.,1.,"B");

     hist_leading_eta_HLT_Jet15U_eff_fine->Divide(hist_leading_eta_HLT_Jet15U_emulated_fine,hist_leading_eta_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_leading_eta_HLT_Jet30U_eff_fine->Divide(hist_leading_eta_HLT_Jet30U_emulated_fine,hist_leading_eta_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_leading_eta_HLT_Jet50U_eff_fine->Divide(hist_leading_eta_HLT_Jet50U_emulated_fine,hist_leading_eta_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_leading_eta_eff_fine->Divide(hist_leading_eta_emulated_fine,hist_leading_eta_all_fine,1.,1.,"B");
     hist_leading_eta_HLT_Jet30U_eff_fine_lowpt->Divide(hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt,hist_leading_eta_HLT_Jet30U_all_fine_lowpt,1.,1.,"B");
     hist_leading_eta_HLT_Jet50U_eff_fine_lowpt->Divide(hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt,hist_leading_eta_HLT_Jet50U_all_fine_lowpt,1.,1.,"B");
     hist_leading_eta_eff_fine_lowpt->Divide(hist_leading_eta_emulated_fine_lowpt,hist_leading_eta_all_fine_lowpt,1.,1.,"B");
     hist_leading_eta_HLT_Jet30U_eff_fine_highpt->Divide(hist_leading_eta_HLT_Jet30U_emulated_fine_highpt,hist_leading_eta_HLT_Jet30U_all_fine_highpt,1.,1.,"B");
     hist_leading_eta_HLT_Jet50U_eff_fine_highpt->Divide(hist_leading_eta_HLT_Jet50U_emulated_fine_highpt,hist_leading_eta_HLT_Jet50U_all_fine_highpt,1.,1.,"B");
     hist_leading_eta_eff_fine_highpt->Divide(hist_leading_eta_emulated_fine_highpt,hist_leading_eta_all_fine_highpt,1.,1.,"B");

     hist_leading_phi_HLT_Jet15U_eff_fine->Divide(hist_leading_phi_HLT_Jet15U_emulated_fine,hist_leading_phi_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_leading_phi_HLT_Jet30U_eff_fine->Divide(hist_leading_phi_HLT_Jet30U_emulated_fine,hist_leading_phi_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_leading_phi_HLT_Jet50U_eff_fine->Divide(hist_leading_phi_HLT_Jet50U_emulated_fine,hist_leading_phi_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_leading_phi_eff_fine->Divide(hist_leading_phi_emulated_fine,hist_leading_phi_all_fine,1.,1.,"B");

     // fine binning histograms for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff_fine->Divide(hist_sj_leading_pt_HLT_Jet15U_emulated_fine,hist_sj_leading_pt_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_sj_leading_pt_HLT_Jet30U_eff_fine->Divide(hist_sj_leading_pt_HLT_Jet30U_emulated_fine,hist_sj_leading_pt_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_sj_leading_pt_HLT_Jet50U_eff_fine->Divide(hist_sj_leading_pt_HLT_Jet50U_emulated_fine,hist_sj_leading_pt_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_sj_leading_pt_eff_fine->Divide(hist_sj_leading_pt_emulated_fine,hist_sj_leading_pt_all_fine,1.,1.,"B");

     hist_sj_leading_eta_HLT_Jet15U_eff_fine->Divide(hist_sj_leading_eta_HLT_Jet15U_emulated_fine,hist_sj_leading_eta_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet30U_eff_fine->Divide(hist_sj_leading_eta_HLT_Jet30U_emulated_fine,hist_sj_leading_eta_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet50U_eff_fine->Divide(hist_sj_leading_eta_HLT_Jet50U_emulated_fine,hist_sj_leading_eta_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_sj_leading_eta_eff_fine->Divide(hist_sj_leading_eta_emulated_fine,hist_sj_leading_eta_all_fine,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_lowpt->Divide(hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt,hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_lowpt->Divide(hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt,hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt,1.,1.,"B");
     hist_sj_leading_eta_eff_fine_lowpt->Divide(hist_sj_leading_eta_emulated_fine_lowpt,hist_sj_leading_eta_all_fine_lowpt,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_highpt->Divide(hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt,hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt,1.,1.,"B");
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_highpt->Divide(hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt,hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt,1.,1.,"B");
     hist_sj_leading_eta_eff_fine_highpt->Divide(hist_sj_leading_eta_emulated_fine_highpt,hist_sj_leading_eta_all_fine_highpt,1.,1.,"B");

     hist_sj_leading_phi_HLT_Jet15U_eff_fine->Divide(hist_sj_leading_phi_HLT_Jet15U_emulated_fine,hist_sj_leading_phi_HLT_Jet15U_all_fine,1.,1.,"B");
     hist_sj_leading_phi_HLT_Jet30U_eff_fine->Divide(hist_sj_leading_phi_HLT_Jet30U_emulated_fine,hist_sj_leading_phi_HLT_Jet30U_all_fine,1.,1.,"B");
     hist_sj_leading_phi_HLT_Jet50U_eff_fine->Divide(hist_sj_leading_phi_HLT_Jet50U_emulated_fine,hist_sj_leading_phi_HLT_Jet50U_all_fine,1.,1.,"B");
     hist_sj_leading_phi_eff_fine->Divide(hist_sj_leading_phi_emulated_fine,hist_sj_leading_phi_all_fine,1.,1.,"B");

     hist_delta_phi_HLT_Jet15U_eff->Divide(hist_delta_phi_HLT_Jet15U_emulated,hist_delta_phi_HLT_Jet15U_all,1.,1.,"B");
     hist_delta_phi_HLT_Jet30U_eff->Divide(hist_delta_phi_HLT_Jet30U_emulated,hist_delta_phi_HLT_Jet30U_all,1.,1.,"B");
     hist_delta_phi_HLT_Jet50U_eff->Divide(hist_delta_phi_HLT_Jet50U_emulated,hist_delta_phi_HLT_Jet50U_all,1.,1.,"B");
     hist_delta_phi_eff->Divide(hist_delta_phi_emulated,hist_delta_phi_all,1.,1.,"B");

     hist_delta_eta_HLT_Jet15U_eff->Divide(hist_delta_eta_HLT_Jet15U_emulated,hist_delta_eta_HLT_Jet15U_all,1.,1.,"B");
     hist_delta_eta_HLT_Jet30U_eff->Divide(hist_delta_eta_HLT_Jet30U_emulated,hist_delta_eta_HLT_Jet30U_all,1.,1.,"B");
     hist_delta_eta_HLT_Jet50U_eff->Divide(hist_delta_eta_HLT_Jet50U_emulated,hist_delta_eta_HLT_Jet50U_all,1.,1.,"B");
     hist_delta_eta_eff->Divide(hist_delta_eta_emulated,hist_delta_eta_all,1.,1.,"B");

     if (detail) { cout<<"Efficiency histograms created sucessfully!"<<endl; }

//Open the output root file
     if (detail) { cout<<"Opening "<<data_out<<endl; }
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

//write histograms on file
     if (detail) { cout<<"Writing histograms on file..."<<endl; }

//general control histograms
     hist_events->Write();
     hist_events_HLT_Jet15U->Write();
     hist_events_HLT_Jet30U->Write();
     hist_events_HLT_Jet50U->Write();
     hist_events_all->Write();
     
//monitor distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_all->Write();
     hist_leading_pt_HLT_Jet30U_all->Write();
     hist_leading_pt_HLT_Jet50U_all->Write();
     hist_leading_pt_all->Write();

     hist_leading_central_pt_HLT_Jet15U_all->Write();
     hist_leading_central_pt_HLT_Jet30U_all->Write();
     hist_leading_central_pt_HLT_Jet50U_all->Write();
     hist_leading_central_pt_all->Write();

     hist_leading_eta_HLT_Jet15U_all->Write();
     hist_leading_eta_HLT_Jet30U_all->Write();
     hist_leading_eta_HLT_Jet50U_all->Write();
     hist_leading_eta_all->Write();
     hist_leading_eta_HLT_Jet30U_all_lowpt->Write();
     hist_leading_eta_HLT_Jet50U_all_lowpt->Write();
     hist_leading_eta_all_lowpt->Write();
     hist_leading_eta_HLT_Jet30U_all_highpt->Write();
     hist_leading_eta_HLT_Jet50U_all_highpt->Write();
     hist_leading_eta_all_highpt->Write();

     hist_leading_phi_HLT_Jet15U_all->Write();
     hist_leading_phi_HLT_Jet30U_all->Write();
     hist_leading_phi_HLT_Jet50U_all->Write();
     hist_leading_phi_all->Write();

//monitor distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_all->Write();
     hist_sj_leading_pt_HLT_Jet30U_all->Write();
     hist_sj_leading_pt_HLT_Jet50U_all->Write();
     hist_sj_leading_pt_all->Write();

     hist_sj_leading_eta_HLT_Jet15U_all->Write();
     hist_sj_leading_eta_HLT_Jet30U_all->Write();
     hist_sj_leading_eta_HLT_Jet50U_all->Write();
     hist_sj_leading_eta_all->Write();
     hist_sj_leading_eta_HLT_Jet30U_all_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_all_lowpt->Write();
     hist_sj_leading_eta_all_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet30U_all_highpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_all_highpt->Write();
     hist_sj_leading_eta_all_highpt->Write();

     hist_sj_leading_phi_HLT_Jet15U_all->Write();
     hist_sj_leading_phi_HLT_Jet30U_all->Write();
     hist_sj_leading_phi_HLT_Jet50U_all->Write();
     hist_sj_leading_phi_all->Write();
     
//triggered distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_emulated->Write();
     hist_leading_pt_HLT_Jet30U_emulated->Write();
     hist_leading_pt_HLT_Jet50U_emulated->Write();
     hist_leading_pt_emulated->Write();

     hist_leading_central_pt_HLT_Jet15U_emulated->Write();
     hist_leading_central_pt_HLT_Jet30U_emulated->Write();
     hist_leading_central_pt_HLT_Jet50U_emulated->Write();
     hist_leading_central_pt_emulated->Write();

     hist_leading_eta_HLT_Jet15U_emulated->Write();
     hist_leading_eta_HLT_Jet30U_emulated->Write();
     hist_leading_eta_HLT_Jet50U_emulated->Write();
     hist_leading_eta_emulated->Write();
     hist_leading_eta_HLT_Jet30U_emulated_lowpt->Write();
     hist_leading_eta_HLT_Jet50U_emulated_lowpt->Write();
     hist_leading_eta_emulated_lowpt->Write();
     hist_leading_eta_HLT_Jet30U_emulated_highpt->Write();
     hist_leading_eta_HLT_Jet50U_emulated_highpt->Write();
     hist_leading_eta_emulated_highpt->Write();

     hist_leading_phi_HLT_Jet15U_emulated->Write();
     hist_leading_phi_HLT_Jet30U_emulated->Write();
     hist_leading_phi_HLT_Jet50U_emulated->Write();
     hist_leading_phi_emulated->Write();

//triggered distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_emulated->Write();
     hist_sj_leading_pt_HLT_Jet30U_emulated->Write();
     hist_sj_leading_pt_HLT_Jet50U_emulated->Write();
     hist_sj_leading_pt_emulated->Write();

     hist_sj_leading_eta_HLT_Jet15U_emulated->Write();
     hist_sj_leading_eta_HLT_Jet30U_emulated->Write();
     hist_sj_leading_eta_HLT_Jet50U_emulated->Write();
     hist_sj_leading_eta_emulated->Write();
     hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt->Write();
     hist_sj_leading_eta_emulated_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet30U_emulated_highpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_emulated_highpt->Write();
     hist_sj_leading_eta_emulated_highpt->Write();

     hist_sj_leading_phi_HLT_Jet15U_emulated->Write();
     hist_sj_leading_phi_HLT_Jet30U_emulated->Write();
     hist_sj_leading_phi_HLT_Jet50U_emulated->Write();
     hist_sj_leading_phi_emulated->Write();
     
//efficiency distributions with broad binning for dijets
     hist_leading_pt_HLT_Jet15U_eff->Write();
     hist_leading_pt_HLT_Jet30U_eff->Write();
     hist_leading_pt_HLT_Jet50U_eff->Write();
     hist_leading_pt_eff->Write();

     hist_leading_central_pt_HLT_Jet15U_eff->Write();
     hist_leading_central_pt_HLT_Jet30U_eff->Write();
     hist_leading_central_pt_HLT_Jet50U_eff->Write();
     hist_leading_central_pt_eff->Write();

     hist_leading_eta_HLT_Jet15U_eff->Write();
     hist_leading_eta_HLT_Jet30U_eff->Write();
     hist_leading_eta_HLT_Jet50U_eff->Write();
     hist_leading_eta_eff->Write();
     hist_leading_eta_HLT_Jet30U_eff_lowpt->Write();
     hist_leading_eta_HLT_Jet50U_eff_lowpt->Write();
     hist_leading_eta_eff_lowpt->Write();
     hist_leading_eta_HLT_Jet30U_eff_highpt->Write();
     hist_leading_eta_HLT_Jet50U_eff_highpt->Write();
     hist_leading_eta_eff_highpt->Write();

     hist_leading_phi_HLT_Jet15U_eff->Write();
     hist_leading_phi_HLT_Jet30U_eff->Write();
     hist_leading_phi_HLT_Jet50U_eff->Write();
     hist_leading_phi_eff->Write();

//efficiency distributions with broad binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff->Write();
     hist_sj_leading_pt_HLT_Jet30U_eff->Write();
     hist_sj_leading_pt_HLT_Jet50U_eff->Write();
     hist_sj_leading_pt_eff->Write();

     hist_sj_leading_eta_HLT_Jet15U_eff->Write();
     hist_sj_leading_eta_HLT_Jet30U_eff->Write();
     hist_sj_leading_eta_HLT_Jet50U_eff->Write();
     hist_sj_leading_eta_eff->Write();
     hist_sj_leading_eta_HLT_Jet30U_eff_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_eff_lowpt->Write();
     hist_sj_leading_eta_eff_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet30U_eff_highpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_eff_highpt->Write();
     hist_sj_leading_eta_eff_highpt->Write();

     hist_sj_leading_phi_HLT_Jet15U_eff->Write();
     hist_sj_leading_phi_HLT_Jet30U_eff->Write();
     hist_sj_leading_phi_HLT_Jet50U_eff->Write();
     hist_sj_leading_phi_eff->Write();

//monitor distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_all_fine->Write();
     hist_leading_pt_HLT_Jet30U_all_fine->Write();
     hist_leading_pt_HLT_Jet50U_all_fine->Write();
     hist_leading_pt_all_fine->Write();

     hist_leading_central_pt_HLT_Jet15U_all_fine->Write();
     hist_leading_central_pt_HLT_Jet30U_all_fine->Write();
     hist_leading_central_pt_HLT_Jet50U_all_fine->Write();
     hist_leading_central_pt_all_fine->Write();

     hist_leading_eta_HLT_Jet15U_all_fine->Write();
     hist_leading_eta_HLT_Jet30U_all_fine->Write();
     hist_leading_eta_HLT_Jet50U_all_fine->Write();
     hist_leading_eta_all_fine->Write();
     hist_leading_eta_HLT_Jet30U_all_fine_lowpt->Write();
     hist_leading_eta_HLT_Jet50U_all_fine_lowpt->Write();
     hist_leading_eta_all_fine_lowpt->Write();
     hist_leading_eta_HLT_Jet30U_all_fine_highpt->Write();
     hist_leading_eta_HLT_Jet50U_all_fine_highpt->Write();
     hist_leading_eta_all_fine_highpt->Write();

     hist_leading_phi_HLT_Jet15U_all_fine->Write();
     hist_leading_phi_HLT_Jet30U_all_fine->Write();
     hist_leading_phi_HLT_Jet50U_all_fine->Write();
     hist_leading_phi_all_fine->Write();

//monitor distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_all_fine->Write();
     hist_sj_leading_pt_HLT_Jet30U_all_fine->Write();
     hist_sj_leading_pt_HLT_Jet50U_all_fine->Write();
     hist_sj_leading_pt_all_fine->Write();

     hist_sj_leading_eta_HLT_Jet15U_all_fine->Write();
     hist_sj_leading_eta_HLT_Jet30U_all_fine->Write();
     hist_sj_leading_eta_HLT_Jet50U_all_fine->Write();
     hist_sj_leading_eta_all_fine->Write();
     hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt->Write();
     hist_sj_leading_eta_all_fine_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt->Write();
     hist_sj_leading_eta_all_fine_highpt->Write();

     hist_sj_leading_phi_HLT_Jet15U_all_fine->Write();
     hist_sj_leading_phi_HLT_Jet30U_all_fine->Write();
     hist_sj_leading_phi_HLT_Jet50U_all_fine->Write();
     hist_sj_leading_phi_all_fine->Write();
     
//triggered distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_emulated_fine->Write();
     hist_leading_pt_HLT_Jet30U_emulated_fine->Write();
     hist_leading_pt_HLT_Jet50U_emulated_fine->Write();
     hist_leading_pt_emulated_fine->Write();

     hist_leading_central_pt_HLT_Jet15U_emulated_fine->Write();
     hist_leading_central_pt_HLT_Jet30U_emulated_fine->Write();
     hist_leading_central_pt_HLT_Jet50U_emulated_fine->Write();
     hist_leading_central_pt_emulated_fine->Write();

     hist_leading_eta_HLT_Jet15U_emulated_fine->Write();
     hist_leading_eta_HLT_Jet30U_emulated_fine->Write();
     hist_leading_eta_HLT_Jet50U_emulated_fine->Write();
     hist_leading_eta_emulated_fine->Write();
     hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt->Write();
     hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt->Write();
     hist_leading_eta_emulated_fine_lowpt->Write();
     hist_leading_eta_HLT_Jet30U_emulated_fine_highpt->Write();
     hist_leading_eta_HLT_Jet50U_emulated_fine_highpt->Write();
     hist_leading_eta_emulated_fine_highpt->Write();

     hist_leading_phi_HLT_Jet15U_emulated_fine->Write();
     hist_leading_phi_HLT_Jet30U_emulated_fine->Write();
     hist_leading_phi_HLT_Jet50U_emulated_fine->Write();
     hist_leading_phi_emulated_fine->Write();

//triggered distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_emulated_fine->Write();
     hist_sj_leading_pt_HLT_Jet30U_emulated_fine->Write();
     hist_sj_leading_pt_HLT_Jet50U_emulated_fine->Write();
     hist_sj_leading_pt_emulated_fine->Write();

     hist_sj_leading_eta_HLT_Jet15U_emulated_fine->Write();
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine->Write();
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine->Write();
     hist_sj_leading_eta_emulated_fine->Write();
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt->Write();
     hist_sj_leading_eta_emulated_fine_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt->Write();
     hist_sj_leading_eta_emulated_fine_highpt->Write();

     hist_sj_leading_phi_HLT_Jet15U_emulated_fine->Write();
     hist_sj_leading_phi_HLT_Jet30U_emulated_fine->Write();
     hist_sj_leading_phi_HLT_Jet50U_emulated_fine->Write();
     hist_sj_leading_phi_emulated_fine->Write();
     
//efficiency distributions with fine binning for dijets
     hist_leading_pt_HLT_Jet15U_eff_fine->Write();
     hist_leading_pt_HLT_Jet30U_eff_fine->Write();
     hist_leading_pt_HLT_Jet50U_eff_fine->Write();
     hist_leading_pt_eff_fine->Write();

     hist_leading_central_pt_HLT_Jet15U_eff_fine->Write();
     hist_leading_central_pt_HLT_Jet30U_eff_fine->Write();
     hist_leading_central_pt_HLT_Jet50U_eff_fine->Write();
     hist_leading_central_pt_eff_fine->Write();

     hist_leading_eta_HLT_Jet15U_eff_fine->Write();
     hist_leading_eta_HLT_Jet30U_eff_fine->Write();
     hist_leading_eta_HLT_Jet50U_eff_fine->Write();
     hist_leading_eta_eff_fine->Write();
     hist_leading_eta_HLT_Jet30U_eff_fine_lowpt->Write();
     hist_leading_eta_HLT_Jet50U_eff_fine_lowpt->Write();
     hist_leading_eta_eff_fine_lowpt->Write();
     hist_leading_eta_HLT_Jet30U_eff_fine_highpt->Write();
     hist_leading_eta_HLT_Jet50U_eff_fine_highpt->Write();
     hist_leading_eta_eff_fine_highpt->Write();

     hist_leading_phi_HLT_Jet15U_eff_fine->Write();
     hist_leading_phi_HLT_Jet30U_eff_fine->Write();
     hist_leading_phi_HLT_Jet50U_eff_fine->Write();
     hist_leading_phi_eff_fine->Write();

//efficiency distributions with fine binning for leading jet
     hist_sj_leading_pt_HLT_Jet15U_eff_fine->Write();
     hist_sj_leading_pt_HLT_Jet30U_eff_fine->Write();
     hist_sj_leading_pt_HLT_Jet50U_eff_fine->Write();
     hist_sj_leading_pt_eff_fine->Write();

     hist_sj_leading_eta_HLT_Jet15U_eff_fine->Write();
     hist_sj_leading_eta_HLT_Jet30U_eff_fine->Write();
     hist_sj_leading_eta_HLT_Jet50U_eff_fine->Write();
     hist_sj_leading_eta_eff_fine->Write();
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_lowpt->Write();
     hist_sj_leading_eta_eff_fine_lowpt->Write();
     hist_sj_leading_eta_HLT_Jet30U_eff_fine_highpt->Write();
     hist_sj_leading_eta_HLT_Jet50U_eff_fine_highpt->Write();
     hist_sj_leading_eta_eff_fine_highpt->Write();

     hist_sj_leading_phi_HLT_Jet15U_eff_fine->Write();
     hist_sj_leading_phi_HLT_Jet30U_eff_fine->Write();
     hist_sj_leading_phi_HLT_Jet50U_eff_fine->Write();
     hist_sj_leading_phi_eff_fine->Write();

//distributions for the observables
     hist_delta_phi_all->Write();
     hist_delta_eta_all->Write();
     hist_delta_phi_emulated->Write();
     hist_delta_eta_emulated->Write();
     hist_delta_phi_eff->Write();
     hist_delta_eta_eff->Write();

     hist_delta_phi_HLT_Jet15U_all->Write();
     hist_delta_eta_HLT_Jet15U_all->Write();
     hist_delta_phi_HLT_Jet15U_emulated->Write();
     hist_delta_eta_HLT_Jet15U_emulated->Write();
     hist_delta_phi_HLT_Jet15U_eff->Write();
     hist_delta_eta_HLT_Jet15U_eff->Write();

     hist_delta_phi_HLT_Jet30U_all->Write();
     hist_delta_eta_HLT_Jet30U_all->Write();
     hist_delta_phi_HLT_Jet30U_emulated->Write();
     hist_delta_eta_HLT_Jet30U_emulated->Write();
     hist_delta_phi_HLT_Jet30U_eff->Write();
     hist_delta_eta_HLT_Jet30U_eff->Write();

     hist_delta_phi_HLT_Jet50U_all->Write();
     hist_delta_eta_HLT_Jet50U_all->Write();
     hist_delta_phi_HLT_Jet50U_emulated->Write();
     hist_delta_eta_HLT_Jet50U_emulated->Write();
     hist_delta_phi_HLT_Jet50U_eff->Write();
     hist_delta_eta_HLT_Jet50U_eff->Write();

     if (detail) { cout<<"Histograms written sucessfully!"<<endl; }
     
     //close the output file
     data_output->Close();
     
//delete the histograms to avoid memory leak
//
     delete(hist_events);
     delete(hist_events_HLT_Jet15U);
     delete(hist_events_HLT_Jet30U);
     delete(hist_events_HLT_Jet50U);
     delete(hist_events_all);

//monitor distributions with broad binning for dijets
     delete(hist_leading_pt_HLT_Jet15U_all);
     delete(hist_leading_pt_HLT_Jet30U_all);
     delete(hist_leading_pt_HLT_Jet50U_all);
     delete(hist_leading_pt_all);

     delete(hist_leading_central_pt_HLT_Jet15U_all);
     delete(hist_leading_central_pt_HLT_Jet30U_all);
     delete(hist_leading_central_pt_HLT_Jet50U_all);
     delete(hist_leading_central_pt_all);

     delete(hist_leading_eta_HLT_Jet15U_all);
     delete(hist_leading_eta_HLT_Jet30U_all);
     delete(hist_leading_eta_HLT_Jet50U_all);
     delete(hist_leading_eta_all);
     delete(hist_leading_eta_HLT_Jet30U_all_lowpt);
     delete(hist_leading_eta_HLT_Jet50U_all_lowpt);
     delete(hist_leading_eta_all_lowpt);
     delete(hist_leading_eta_HLT_Jet30U_all_highpt);
     delete(hist_leading_eta_HLT_Jet50U_all_highpt);
     delete(hist_leading_eta_all_highpt);

     delete(hist_leading_phi_HLT_Jet15U_all);
     delete(hist_leading_phi_HLT_Jet30U_all);
     delete(hist_leading_phi_HLT_Jet50U_all);
     delete(hist_leading_phi_all);

//monitor distributions with broad binning for leading jet
     delete(hist_sj_leading_pt_HLT_Jet15U_all);
     delete(hist_sj_leading_pt_HLT_Jet30U_all);
     delete(hist_sj_leading_pt_HLT_Jet50U_all);
     delete(hist_sj_leading_pt_all);

     delete(hist_sj_leading_eta_HLT_Jet15U_all);
     delete(hist_sj_leading_eta_HLT_Jet30U_all);
     delete(hist_sj_leading_eta_HLT_Jet50U_all);
     delete(hist_sj_leading_eta_all);
     delete(hist_sj_leading_eta_HLT_Jet30U_all_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_all_lowpt);
     delete(hist_sj_leading_eta_all_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet30U_all_highpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_all_highpt);
     delete(hist_sj_leading_eta_all_highpt);

     delete(hist_sj_leading_phi_HLT_Jet15U_all);
     delete(hist_sj_leading_phi_HLT_Jet30U_all);
     delete(hist_sj_leading_phi_HLT_Jet50U_all);
     delete(hist_sj_leading_phi_all);

//triggered distributions with broad binning for dijets
     delete(hist_leading_pt_HLT_Jet15U_emulated);
     delete(hist_leading_pt_HLT_Jet30U_emulated);
     delete(hist_leading_pt_HLT_Jet50U_emulated);
     delete(hist_leading_pt_emulated);

     delete(hist_leading_central_pt_HLT_Jet15U_emulated);
     delete(hist_leading_central_pt_HLT_Jet30U_emulated);
     delete(hist_leading_central_pt_HLT_Jet50U_emulated);
     delete(hist_leading_central_pt_emulated);

     delete(hist_leading_eta_HLT_Jet15U_emulated);
     delete(hist_leading_eta_HLT_Jet30U_emulated);
     delete(hist_leading_eta_HLT_Jet50U_emulated);
     delete(hist_leading_eta_emulated);
     delete(hist_leading_eta_HLT_Jet30U_emulated_lowpt);
     delete(hist_leading_eta_HLT_Jet50U_emulated_lowpt);
     delete(hist_leading_eta_emulated_lowpt);
     delete(hist_leading_eta_HLT_Jet30U_emulated_highpt);
     delete(hist_leading_eta_HLT_Jet50U_emulated_highpt);
     delete(hist_leading_eta_emulated_highpt);

     delete(hist_leading_phi_HLT_Jet15U_emulated);
     delete(hist_leading_phi_HLT_Jet30U_emulated);
     delete(hist_leading_phi_HLT_Jet50U_emulated);
     delete(hist_leading_phi_emulated);

//triggered distributions with broad binning for leading jet
     delete(hist_sj_leading_pt_HLT_Jet15U_emulated);
     delete(hist_sj_leading_pt_HLT_Jet30U_emulated);
     delete(hist_sj_leading_pt_HLT_Jet50U_emulated);
     delete(hist_sj_leading_pt_emulated);

     delete(hist_sj_leading_eta_HLT_Jet15U_emulated);
     delete(hist_sj_leading_eta_HLT_Jet30U_emulated);
     delete(hist_sj_leading_eta_HLT_Jet50U_emulated);
     delete(hist_sj_leading_eta_emulated);
     delete(hist_sj_leading_eta_HLT_Jet30U_emulated_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_emulated_lowpt);
     delete(hist_sj_leading_eta_emulated_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet30U_emulated_highpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_emulated_highpt);
     delete(hist_sj_leading_eta_emulated_highpt);

     delete(hist_sj_leading_phi_HLT_Jet15U_emulated);
     delete(hist_sj_leading_phi_HLT_Jet30U_emulated);
     delete(hist_sj_leading_phi_HLT_Jet50U_emulated);
     delete(hist_sj_leading_phi_emulated);

//efficiency distributions with broad binning for dijets
     delete(hist_leading_pt_HLT_Jet15U_eff);
     delete(hist_leading_pt_HLT_Jet30U_eff);
     delete(hist_leading_pt_HLT_Jet50U_eff);
     delete(hist_leading_pt_eff);

     delete(hist_leading_central_pt_HLT_Jet15U_eff);
     delete(hist_leading_central_pt_HLT_Jet30U_eff);
     delete(hist_leading_central_pt_HLT_Jet50U_eff);
     delete(hist_leading_central_pt_eff);

     delete(hist_leading_eta_HLT_Jet15U_eff);
     delete(hist_leading_eta_HLT_Jet30U_eff);
     delete(hist_leading_eta_HLT_Jet50U_eff);
     delete(hist_leading_eta_eff);
     delete(hist_leading_eta_HLT_Jet30U_eff_lowpt);
     delete(hist_leading_eta_HLT_Jet50U_eff_lowpt);
     delete(hist_leading_eta_eff_lowpt);
     delete(hist_leading_eta_HLT_Jet30U_eff_highpt);
     delete(hist_leading_eta_HLT_Jet50U_eff_highpt);
     delete(hist_leading_eta_eff_highpt);

     delete(hist_leading_phi_HLT_Jet15U_eff);
     delete(hist_leading_phi_HLT_Jet30U_eff);
     delete(hist_leading_phi_HLT_Jet50U_eff);
     delete(hist_leading_phi_eff);

//efficiency distributions with broad binning for leading jet
     delete(hist_sj_leading_pt_HLT_Jet15U_eff);
     delete(hist_sj_leading_pt_HLT_Jet30U_eff);
     delete(hist_sj_leading_pt_HLT_Jet50U_eff);
     delete(hist_sj_leading_pt_eff);

     delete(hist_sj_leading_eta_HLT_Jet15U_eff);
     delete(hist_sj_leading_eta_HLT_Jet30U_eff);
     delete(hist_sj_leading_eta_HLT_Jet50U_eff);
     delete(hist_sj_leading_eta_eff);
     delete(hist_sj_leading_eta_HLT_Jet30U_eff_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_eff_lowpt);
     delete(hist_sj_leading_eta_eff_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet30U_eff_highpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_eff_highpt);
     delete(hist_sj_leading_eta_eff_highpt);

     delete(hist_sj_leading_phi_HLT_Jet15U_eff);
     delete(hist_sj_leading_phi_HLT_Jet30U_eff);
     delete(hist_sj_leading_phi_HLT_Jet50U_eff);
     delete(hist_sj_leading_phi_eff);

//monitor distributions with fine binning for dijets
     delete(hist_leading_pt_HLT_Jet15U_all_fine);
     delete(hist_leading_pt_HLT_Jet30U_all_fine);
     delete(hist_leading_pt_HLT_Jet50U_all_fine);
     delete(hist_leading_pt_all_fine);

     delete(hist_leading_central_pt_HLT_Jet15U_all_fine);
     delete(hist_leading_central_pt_HLT_Jet30U_all_fine);
     delete(hist_leading_central_pt_HLT_Jet50U_all_fine);
     delete(hist_leading_central_pt_all_fine);

     delete(hist_leading_eta_HLT_Jet15U_all_fine);
     delete(hist_leading_eta_HLT_Jet30U_all_fine);
     delete(hist_leading_eta_HLT_Jet50U_all_fine);
     delete(hist_leading_eta_all_fine);
     delete(hist_leading_eta_HLT_Jet30U_all_fine_lowpt);
     delete(hist_leading_eta_HLT_Jet50U_all_fine_lowpt);
     delete(hist_leading_eta_all_fine_lowpt);
     delete(hist_leading_eta_HLT_Jet30U_all_fine_highpt);
     delete(hist_leading_eta_HLT_Jet50U_all_fine_highpt);
     delete(hist_leading_eta_all_fine_highpt);

     delete(hist_leading_phi_HLT_Jet15U_all_fine);
     delete(hist_leading_phi_HLT_Jet30U_all_fine);
     delete(hist_leading_phi_HLT_Jet50U_all_fine);
     delete(hist_leading_phi_all_fine);

//monitor distributions with fine binning for leading jet
     delete(hist_sj_leading_pt_HLT_Jet15U_all_fine);
     delete(hist_sj_leading_pt_HLT_Jet30U_all_fine);
     delete(hist_sj_leading_pt_HLT_Jet50U_all_fine);
     delete(hist_sj_leading_pt_all_fine);

     delete(hist_sj_leading_eta_HLT_Jet15U_all_fine);
     delete(hist_sj_leading_eta_HLT_Jet30U_all_fine);
     delete(hist_sj_leading_eta_HLT_Jet50U_all_fine);
     delete(hist_sj_leading_eta_all_fine);
     delete(hist_sj_leading_eta_HLT_Jet30U_all_fine_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_all_fine_lowpt);
     delete(hist_sj_leading_eta_all_fine_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet30U_all_fine_highpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_all_fine_highpt);
     delete(hist_sj_leading_eta_all_fine_highpt);

     delete(hist_sj_leading_phi_HLT_Jet15U_all_fine);
     delete(hist_sj_leading_phi_HLT_Jet30U_all_fine);
     delete(hist_sj_leading_phi_HLT_Jet50U_all_fine);
     delete(hist_sj_leading_phi_all_fine);
         
//triggered distributions with fine binning for dijets
     delete(hist_leading_pt_HLT_Jet15U_emulated_fine);
     delete(hist_leading_pt_HLT_Jet30U_emulated_fine);
     delete(hist_leading_pt_HLT_Jet50U_emulated_fine);
     delete(hist_leading_pt_emulated_fine);

     delete(hist_leading_central_pt_HLT_Jet15U_emulated_fine);
     delete(hist_leading_central_pt_HLT_Jet30U_emulated_fine);
     delete(hist_leading_central_pt_HLT_Jet50U_emulated_fine);
     delete(hist_leading_central_pt_emulated_fine);

     delete(hist_leading_eta_HLT_Jet15U_emulated_fine);
     delete(hist_leading_eta_HLT_Jet30U_emulated_fine);
     delete(hist_leading_eta_HLT_Jet50U_emulated_fine);
     delete(hist_leading_eta_emulated_fine);
     delete(hist_leading_eta_HLT_Jet30U_emulated_fine_lowpt);
     delete(hist_leading_eta_HLT_Jet50U_emulated_fine_lowpt);
     delete(hist_leading_eta_emulated_fine_lowpt);
     delete(hist_leading_eta_HLT_Jet30U_emulated_fine_highpt);
     delete(hist_leading_eta_HLT_Jet50U_emulated_fine_highpt);
     delete(hist_leading_eta_emulated_fine_highpt);

     delete(hist_leading_phi_HLT_Jet15U_emulated_fine);
     delete(hist_leading_phi_HLT_Jet30U_emulated_fine);
     delete(hist_leading_phi_HLT_Jet50U_emulated_fine);
     delete(hist_leading_phi_emulated_fine);

//triggered distributions with fine binning for single jet
     delete(hist_sj_leading_pt_HLT_Jet15U_emulated_fine);
     delete(hist_sj_leading_pt_HLT_Jet30U_emulated_fine);
     delete(hist_sj_leading_pt_HLT_Jet50U_emulated_fine);
     delete(hist_sj_leading_pt_emulated_fine);

     delete(hist_sj_leading_eta_HLT_Jet15U_emulated_fine);
     delete(hist_sj_leading_eta_HLT_Jet30U_emulated_fine);
     delete(hist_sj_leading_eta_HLT_Jet50U_emulated_fine);
     delete(hist_sj_leading_eta_emulated_fine);
     delete(hist_sj_leading_eta_HLT_Jet30U_emulated_fine_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_emulated_fine_lowpt);
     delete(hist_sj_leading_eta_emulated_fine_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet30U_emulated_fine_highpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_emulated_fine_highpt);
     delete(hist_sj_leading_eta_emulated_fine_highpt);

     delete(hist_sj_leading_phi_HLT_Jet15U_emulated_fine);
     delete(hist_sj_leading_phi_HLT_Jet30U_emulated_fine);
     delete(hist_sj_leading_phi_HLT_Jet50U_emulated_fine);
     delete(hist_sj_leading_phi_emulated_fine); 

//efficiency distributions with fine binning for dijets
     delete(hist_leading_pt_HLT_Jet15U_eff_fine);
     delete(hist_leading_pt_HLT_Jet30U_eff_fine);
     delete(hist_leading_pt_HLT_Jet50U_eff_fine);
     delete(hist_leading_pt_eff_fine);

     delete(hist_leading_central_pt_HLT_Jet15U_eff_fine);
     delete(hist_leading_central_pt_HLT_Jet30U_eff_fine);
     delete(hist_leading_central_pt_HLT_Jet50U_eff_fine);
     delete(hist_leading_central_pt_eff_fine);

     delete(hist_leading_eta_HLT_Jet15U_eff_fine);
     delete(hist_leading_eta_HLT_Jet30U_eff_fine);
     delete(hist_leading_eta_HLT_Jet50U_eff_fine);
     delete(hist_leading_eta_eff_fine);
     delete(hist_leading_eta_HLT_Jet30U_eff_fine_lowpt);
     delete(hist_leading_eta_HLT_Jet50U_eff_fine_lowpt);
     delete(hist_leading_eta_eff_fine_lowpt);
     delete(hist_leading_eta_HLT_Jet30U_eff_fine_highpt);
     delete(hist_leading_eta_HLT_Jet50U_eff_fine_highpt);
     delete(hist_leading_eta_eff_fine_highpt);

     delete(hist_leading_phi_HLT_Jet15U_eff_fine);
     delete(hist_leading_phi_HLT_Jet30U_eff_fine);
     delete(hist_leading_phi_HLT_Jet50U_eff_fine);
     delete(hist_leading_phi_eff_fine);

//efficiency distributions with fine binning for single jet
     delete(hist_sj_leading_pt_HLT_Jet15U_eff_fine);
     delete(hist_sj_leading_pt_HLT_Jet30U_eff_fine);
     delete(hist_sj_leading_pt_HLT_Jet50U_eff_fine);
     delete(hist_sj_leading_pt_eff_fine);

     delete(hist_sj_leading_eta_HLT_Jet15U_eff_fine);
     delete(hist_sj_leading_eta_HLT_Jet30U_eff_fine);
     delete(hist_sj_leading_eta_HLT_Jet50U_eff_fine);
     delete(hist_sj_leading_eta_eff_fine);
     delete(hist_sj_leading_eta_HLT_Jet30U_eff_fine_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_eff_fine_lowpt);
     delete(hist_sj_leading_eta_eff_fine_lowpt);
     delete(hist_sj_leading_eta_HLT_Jet30U_eff_fine_highpt);
     delete(hist_sj_leading_eta_HLT_Jet50U_eff_fine_highpt);
     delete(hist_sj_leading_eta_eff_fine_highpt);

     delete(hist_sj_leading_phi_HLT_Jet15U_eff_fine);
     delete(hist_sj_leading_phi_HLT_Jet30U_eff_fine);
     delete(hist_sj_leading_phi_HLT_Jet50U_eff_fine);
     delete(hist_sj_leading_phi_eff_fine);

//distributions for the observables
     delete(hist_delta_phi_all);
     delete(hist_delta_eta_all);
     delete(hist_delta_phi_emulated);
     delete(hist_delta_eta_emulated);
     delete(hist_delta_phi_eff);
     delete(hist_delta_eta_eff);

     delete(hist_delta_phi_HLT_Jet15U_all);
     delete(hist_delta_eta_HLT_Jet15U_all);
     delete(hist_delta_phi_HLT_Jet15U_emulated);
     delete(hist_delta_eta_HLT_Jet15U_emulated);
     delete(hist_delta_phi_HLT_Jet15U_eff);
     delete(hist_delta_eta_HLT_Jet15U_eff);

     delete(hist_delta_phi_HLT_Jet30U_all);
     delete(hist_delta_eta_HLT_Jet30U_all);
     delete(hist_delta_phi_HLT_Jet30U_emulated);
     delete(hist_delta_eta_HLT_Jet30U_emulated);
     delete(hist_delta_phi_HLT_Jet30U_eff);
     delete(hist_delta_eta_HLT_Jet30U_eff);

     delete(hist_delta_phi_HLT_Jet50U_all);
     delete(hist_delta_eta_HLT_Jet50U_all);
     delete(hist_delta_phi_HLT_Jet50U_emulated);
     delete(hist_delta_eta_HLT_Jet50U_emulated);
     delete(hist_delta_phi_HLT_Jet50U_eff);
     delete(hist_delta_eta_HLT_Jet50U_eff);

//Sucess confirmation
     if (detail) { cout<<"Done!"<<endl; }
}
