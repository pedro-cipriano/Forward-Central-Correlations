// Pedro Cipriano, Jan 2012
// DESY, CMS
// Last Update: 31 Jan 2013
//
// 

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TPad.h>
#include <TString.h>
#include <TMath.h>
#include <TLegend.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "common_methods.h"

using namespace std;

void compute_trigger_correction(string root_in, string file_out, string output_path = "../output/trigger_correction/", string prefix = "test", string sufix = "_v0", bool issj = false, bool detail = false, bool test = false)
{

//output the configuration
   if (detail) { cout<<"Compute Trigger Correction Configuration"<<endl; }
   if (detail) { cout<<"Root In :              "<<root_in<<endl; }
   if (detail) { cout<<"Root Out :             "<<file_out<<endl; }
   if (detail) { cout<<"Output Path :          "<<output_path<<endl; }
   if (detail) { cout<<"Prefix :               "<<prefix<<endl; }
   if (detail) { cout<<"Sufix :                "<<sufix<<endl; }
   if (detail) { cout<<"Is Single Selection? : "<<issj<<endl; }
   if (detail) { cout<<"Detail level :         "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :            "<<test<<endl; }

//opening the input data files
    if (detail) { cout<<"Opening Root files... "<<endl; }
    TFile *data_in = new TFile( root_in.c_str() );
    if (detail) { cout<<"All files opened sucessfully!"<<endl; }

//setting vars for single jet
    TString hist_mod = "";
    if (issj) { hist_mod = "sj_"; }
    if (issj) { prefix =  hist_mod + prefix; }

//loading histograms
    if (detail) { cout<<"Loading histograms... "<<endl; }
    bool not_loaded = false;
    //HLT_Jet15U pt
    TH1D *pt_all_Jet15 = 0;
    TH1D *pt_emulated_Jet15 = 0;
    TH1D *pt_eff_Jet15 = 0;
    TH1D *pt_all_fine_Jet15 = 0;
    TH1D *pt_emulated_fine_Jet15 = 0;
    TH1D *pt_eff_fine_Jet15 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT15U_all",pt_all_Jet15);
    if (pt_all_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT15U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT15U_emulated",pt_emulated_Jet15);
    if (pt_emulated_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT15U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT15U_eff",pt_eff_Jet15);
    if (pt_eff_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT15U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT15U_all_fine",pt_all_fine_Jet15);
    if (pt_all_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT15U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT15U_emulated_fine",pt_emulated_fine_Jet15);
    if (pt_emulated_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT15U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT15U_eff_fine",pt_eff_fine_Jet15);
    if (pt_eff_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT15U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet15U phi
    TH1D *phi_all_Jet15 = 0;
    TH1D *phi_emulated_Jet15 = 0;
    TH1D *phi_eff_Jet15 = 0;
    TH1D *phi_all_fine_Jet15 = 0;
    TH1D *phi_emulated_fine_Jet15 = 0;
    TH1D *phi_eff_fine_Jet15 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT15U_all",phi_all_Jet15);
    if (phi_all_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT15U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT15U_emulated",phi_emulated_Jet15);
    if (phi_emulated_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT15U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT15U_eff",phi_eff_Jet15);
    if (phi_eff_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT15U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT15U_all_fine",phi_all_fine_Jet15);
    if (phi_all_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT15U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT15U_emulated_fine",phi_emulated_fine_Jet15);
    if (phi_emulated_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT15U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT15U_eff_fine",phi_eff_fine_Jet15);
    if (phi_eff_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT15U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet15U eta
    TH1D *eta_all_Jet15 = 0;
    TH1D *eta_emulated_Jet15 = 0;
    TH1D *eta_eff_Jet15 = 0;
    TH1D *eta_all_fine_Jet15 = 0;
    TH1D *eta_emulated_fine_Jet15 = 0;
    TH1D *eta_eff_fine_Jet15 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT15U_all",eta_all_Jet15);
    if (eta_all_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT15U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT15U_emulated",eta_emulated_Jet15);
    if (eta_emulated_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT15U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT15U_eff",eta_eff_Jet15);
    if (eta_eff_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT15U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT15U_all_fine",eta_all_fine_Jet15);
    if (eta_all_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT15U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT15U_emulated_fine",eta_emulated_fine_Jet15);
    if (eta_emulated_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT15U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT15U_eff_fine",eta_eff_fine_Jet15);
    if (eta_eff_fine_Jet15 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT15U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet15U delta eta and delta phi
    TH1D *delta_phi_all_Jet15 = 0;
    TH1D *delta_phi_emulated_Jet15 = 0;
    TH1D *delta_phi_eff_Jet15 = 0;
    TH1D *delta_eta_all_Jet15 = 0;
    TH1D *delta_eta_emulated_Jet15 = 0;
    TH1D *delta_eta_eff_Jet15 = 0;
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet15U_all",delta_phi_all_Jet15);
    if (delta_phi_all_Jet15 == 0) { cout << "ak5PF_delta_phi_HLT_Jet15U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet15U_emulated",delta_phi_emulated_Jet15);
    if (delta_phi_emulated_Jet15 == 0) { cout << "ak5PF_delta_phi_HLT_Jet15U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet15U_eff",delta_phi_eff_Jet15);
    if (delta_phi_eff_Jet15 == 0) { cout << "ak5PF_delta_phi_HLT_Jet15U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet15U_all",delta_eta_all_Jet15);
    if (delta_eta_all_Jet15 == 0) { cout << "ak5PF_delta_eta_HLT_Jet15U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet15U_emulated",delta_eta_emulated_Jet15);
    if (delta_eta_emulated_Jet15 == 0) { cout << "ak5PF_delta_eta_HLT_Jet15U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet15U_eff",delta_eta_eff_Jet15);
    if (delta_eta_eff_Jet15 == 0) { cout << "ak5PF_delta_eta_HLT_Jet15U_eff not found!" << endl; not_loaded = true; }

    //HLT_Jet30U pt
    TH1D *pt_all_Jet30 = 0;
    TH1D *pt_emulated_Jet30 = 0;
    TH1D *pt_eff_Jet30 = 0;
    TH1D *pt_all_fine_Jet30 = 0;
    TH1D *pt_emulated_fine_Jet30 = 0;
    TH1D *pt_eff_fine_Jet30 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT30U_all",pt_all_Jet30);
    if (pt_all_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT30U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT30U_emulated",pt_emulated_Jet30);
    if (pt_emulated_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT30U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT30U_eff",pt_eff_Jet30);
    if (pt_eff_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT30U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT30U_all_fine",pt_all_fine_Jet30);
    if (pt_all_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT30U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT30U_emulated_fine",pt_emulated_fine_Jet30);
    if (pt_emulated_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT30U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT30U_eff_fine",pt_eff_fine_Jet30);
    if (pt_eff_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT30U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet30U phi
    TH1D *phi_all_Jet30 = 0;
    TH1D *phi_emulated_Jet30 = 0;
    TH1D *phi_eff_Jet30 = 0;
    TH1D *phi_all_fine_Jet30 = 0;
    TH1D *phi_emulated_fine_Jet30 = 0;
    TH1D *phi_eff_fine_Jet30 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT30U_all",phi_all_Jet30);
    if (phi_all_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT30U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT30U_emulated",phi_emulated_Jet30);
    if (phi_emulated_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT30U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT30U_eff",phi_eff_Jet30);
    if (phi_eff_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT30U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT30U_all_fine",phi_all_fine_Jet30);
    if (phi_all_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT30U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT30U_emulated_fine",phi_emulated_fine_Jet30);
    if (phi_emulated_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT30U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT30U_eff_fine",phi_eff_fine_Jet30);
    if (phi_eff_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT30U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet30U eta
    TH1D *eta_all_Jet30 = 0;
    TH1D *eta_emulated_Jet30 = 0;
    TH1D *eta_eff_Jet30 = 0;
    TH1D *eta_all_fine_Jet30 = 0;
    TH1D *eta_emulated_fine_Jet30 = 0;
    TH1D *eta_eff_fine_Jet30 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_all",eta_all_Jet30);
    if (eta_all_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated",eta_emulated_Jet30);
    if (eta_emulated_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_eff",eta_eff_Jet30);
    if (eta_eff_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_all_fine",eta_all_fine_Jet30);
    if (eta_all_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_fine",eta_emulated_fine_Jet30);
    if (eta_emulated_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_fine",eta_eff_fine_Jet30);
    if (eta_eff_fine_Jet30 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet30U eta low pt
    TH1D *eta_all_Jet30_lowpt = 0;
    TH1D *eta_emulated_Jet30_lowpt = 0;
    TH1D *eta_eff_Jet30_lowpt = 0;
    TH1D *eta_all_fine_Jet30_lowpt = 0;
    TH1D *eta_emulated_fine_Jet30_lowpt = 0;
    TH1D *eta_eff_fine_Jet30_lowpt = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_all_lowpt",eta_all_Jet30_lowpt);
    if (eta_all_Jet30_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_all_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_lowpt",eta_emulated_Jet30_lowpt);
    if (eta_emulated_Jet30_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_lowpt",eta_eff_Jet30_lowpt);
    if (eta_eff_Jet30_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_all_fine_lowpt",eta_all_fine_Jet30_lowpt);
    if (eta_all_fine_Jet30_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_all_fine_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_fine_lowpt",eta_emulated_fine_Jet30_lowpt);
    if (eta_emulated_fine_Jet30_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_fine_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_fine_lowpt",eta_eff_fine_Jet30_lowpt);
    if (eta_eff_fine_Jet30_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_fine_lowpt not found!" << endl; not_loaded = true; }

    //HLT_Jet30U eta high pt
    TH1D *eta_all_Jet30_highpt = 0;
    TH1D *eta_emulated_Jet30_highpt = 0;
    TH1D *eta_eff_Jet30_highpt = 0;
    TH1D *eta_all_fine_Jet30_highpt = 0;
    TH1D *eta_emulated_fine_Jet30_highpt = 0;
    TH1D *eta_eff_fine_Jet30_highpt = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_all_highpt",eta_all_Jet30_highpt);
    if (eta_all_Jet30_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_all_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_highpt",eta_emulated_Jet30_highpt);
    if (eta_emulated_Jet30_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_highpt",eta_eff_Jet30_highpt);
    if (eta_eff_Jet30_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_all_fine_highpt",eta_all_fine_Jet30_highpt);
    if (eta_all_fine_Jet30_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_all_fine_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_fine_highpt",eta_emulated_fine_Jet30_highpt);
    if (eta_emulated_fine_Jet30_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_emulated_fine_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_fine_highpt",eta_eff_fine_Jet30_highpt);
    if (eta_eff_fine_Jet30_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT30U_eff_fine_highpt not found!" << endl; not_loaded = true; }

   //HLT_Jet30U delta eta and delta phi
    TH1D *delta_phi_all_Jet30 = 0;
    TH1D *delta_phi_emulated_Jet30 = 0;
    TH1D *delta_phi_eff_Jet30 = 0;
    TH1D *delta_eta_all_Jet30 = 0;
    TH1D *delta_eta_emulated_Jet30 = 0;
    TH1D *delta_eta_eff_Jet30 = 0;
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet30U_all",delta_phi_all_Jet30);
    if (delta_phi_all_Jet30 == 0) { cout << "ak5PF_delta_phi_HLT_Jet30U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet30U_emulated",delta_phi_emulated_Jet30);
    if (delta_phi_emulated_Jet30 == 0) { cout << "ak5PF_delta_phi_HLT_Jet30U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet30U_eff",delta_phi_eff_Jet30);
    if (delta_phi_eff_Jet30 == 0) { cout << "ak5PF_delta_phi_HLT_Jet30U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet30U_all",delta_eta_all_Jet30);
    if (delta_eta_all_Jet30 == 0) { cout << "ak5PF_delta_eta_HLT_Jet30U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet30U_emulated",delta_eta_emulated_Jet30);
    if (delta_eta_emulated_Jet30 == 0) { cout << "ak5PF_delta_eta_HLT_Jet30U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet30U_eff",delta_eta_eff_Jet30);
    if (delta_eta_eff_Jet30 == 0) { cout << "ak5PF_delta_eta_HLT_Jet30U_eff not found!" << endl; not_loaded = true; }


    //HLT_Jet50U pt
    TH1D *pt_all_Jet50 = 0;
    TH1D *pt_emulated_Jet50 = 0;
    TH1D *pt_eff_Jet50 = 0;
    TH1D *pt_all_fine_Jet50 = 0;
    TH1D *pt_emulated_fine_Jet50 = 0;
    TH1D *pt_eff_fine_Jet50 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT50U_all",pt_all_Jet50);
    if (pt_all_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT50U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT50U_emulated",pt_emulated_Jet50);
    if (pt_emulated_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT50U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT50U_eff",pt_eff_Jet50);
    if (pt_eff_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT50U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT50U_all_fine",pt_all_fine_Jet50);
    if (pt_all_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT50U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT50U_emulated_fine",pt_emulated_fine_Jet50);
    if (pt_emulated_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT50U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_HLT50U_eff_fine",pt_eff_fine_Jet50);
    if (pt_eff_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_HLT50U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet50U phi
    TH1D *phi_all_Jet50 = 0;
    TH1D *phi_emulated_Jet50 = 0;
    TH1D *phi_eff_Jet50 = 0;
    TH1D *phi_all_fine_Jet50 = 0;
    TH1D *phi_emulated_fine_Jet50 = 0;
    TH1D *phi_eff_fine_Jet50 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT50U_all",phi_all_Jet50);
    if (phi_all_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT50U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT50U_emulated",phi_emulated_Jet50);
    if (phi_emulated_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT50U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT50U_eff",phi_eff_Jet50);
    if (phi_eff_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT50U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT50U_all_fine",phi_all_fine_Jet50);
    if (phi_all_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT50U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT50U_emulated_fine",phi_emulated_fine_Jet50);
    if (phi_emulated_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT50U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_HLT50U_eff_fine",phi_eff_fine_Jet50);
    if (phi_eff_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_HLT50U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet50U eta
    TH1D *eta_all_Jet50 = 0;
    TH1D *eta_emulated_Jet50 = 0;
    TH1D *eta_eff_Jet50 = 0;
    TH1D *eta_all_fine_Jet50 = 0;
    TH1D *eta_emulated_fine_Jet50 = 0;
    TH1D *eta_eff_fine_Jet50 = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_all",eta_all_Jet50);
    if (eta_all_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated",eta_emulated_Jet50);
    if (eta_emulated_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_eff",eta_eff_Jet50);
    if (eta_eff_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_all_fine",eta_all_fine_Jet50);
    if (eta_all_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_fine",eta_emulated_fine_Jet50);
    if (eta_emulated_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_fine",eta_eff_fine_Jet50);
    if (eta_eff_fine_Jet50 == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_fine not found!" << endl; not_loaded = true; }

    //HLT_Jet50U eta low pt
    TH1D *eta_all_Jet50_lowpt = 0;
    TH1D *eta_emulated_Jet50_lowpt = 0;
    TH1D *eta_eff_Jet50_lowpt = 0;
    TH1D *eta_all_fine_Jet50_lowpt = 0;
    TH1D *eta_emulated_fine_Jet50_lowpt = 0;
    TH1D *eta_eff_fine_Jet50_lowpt = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_all_lowpt",eta_all_Jet50_lowpt);
    if (eta_all_Jet50_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_all_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_lowpt",eta_emulated_Jet50_lowpt);
    if (eta_emulated_Jet50_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_lowpt",eta_eff_Jet50_lowpt);
    if (eta_eff_Jet50_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_all_fine_lowpt",eta_all_fine_Jet50_lowpt);
    if (eta_all_fine_Jet50_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_all_fine_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_fine_lowpt",eta_emulated_fine_Jet50_lowpt);
    if (eta_emulated_fine_Jet50_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_fine_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_fine_lowpt",eta_eff_fine_Jet50_lowpt);
    if (eta_eff_fine_Jet50_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_fine_lowpt not found!" << endl; not_loaded = true; }

    //HLT_Jet50U eta high pt
    TH1D *eta_all_Jet50_highpt = 0;
    TH1D *eta_emulated_Jet50_highpt = 0;
    TH1D *eta_eff_Jet50_highpt = 0;
    TH1D *eta_all_fine_Jet50_highpt = 0;
    TH1D *eta_emulated_fine_Jet50_highpt = 0;
    TH1D *eta_eff_fine_Jet50_highpt = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_all_highpt",eta_all_Jet50_highpt);
    if (eta_all_Jet50_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_all_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_highpt",eta_emulated_Jet50_highpt);
    if (eta_emulated_Jet50_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_highpt",eta_eff_Jet50_highpt);
    if (eta_eff_Jet50_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_all_fine_highpt",eta_all_fine_Jet50_highpt);
    if (eta_all_fine_Jet50_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_all_fine_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_fine_highpt",eta_emulated_fine_Jet50_highpt);
    if (eta_emulated_fine_Jet50_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_emulated_fine_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_fine_highpt",eta_eff_fine_Jet50_highpt);
    if (eta_eff_fine_Jet50_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_HLT50U_eff_fine_highpt not found!" << endl; not_loaded = true; }

    //HLT_Jet50U delta eta and delta phi
    TH1D *delta_phi_all_Jet50 = 0;
    TH1D *delta_phi_emulated_Jet50 = 0;
    TH1D *delta_phi_eff_Jet50 = 0;
    TH1D *delta_eta_all_Jet50 = 0;
    TH1D *delta_eta_emulated_Jet50 = 0;
    TH1D *delta_eta_eff_Jet50 = 0;
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet50U_all",delta_phi_all_Jet50);
    if (delta_phi_all_Jet50 == 0) { cout << "ak5PF_delta_phi_HLT_Jet50U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet50U_emulated",delta_phi_emulated_Jet50);
    if (delta_phi_emulated_Jet50 == 0) { cout << "ak5PF_delta_phi_HLT_Jet50U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_HLT_Jet50U_eff",delta_phi_eff_Jet50);
    if (delta_phi_eff_Jet50 == 0) { cout << "ak5PF_delta_phi_HLT_Jet50U_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet50U_all",delta_eta_all_Jet50);
    if (delta_eta_all_Jet50 == 0) { cout << "ak5PF_delta_eta_HLT_Jet50U_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet50U_emulated",delta_eta_emulated_Jet50);
    if (delta_eta_emulated_Jet50 == 0) { cout << "ak5PF_delta_eta_HLT_Jet50U_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_HLT_Jet50U_eff",delta_eta_eff_Jet50);
    if (delta_eta_eff_Jet50 == 0) { cout << "ak5PF_delta_eta_HLT_Jet50U_eff not found!" << endl; not_loaded = true; }

    //Combination pt
    TH1D *pt_all = 0;
    TH1D *pt_emulated = 0;
    TH1D *pt_eff = 0;
    TH1D *pt_all_fine = 0;
    TH1D *pt_emulated_fine = 0;
    TH1D *pt_eff_fine = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_all",pt_all);
    if (pt_all == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_emulated",pt_emulated);
    if (pt_emulated == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_eff",pt_eff);
    if (pt_eff == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_all_fine",pt_all_fine);
    if (pt_all_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_emulated_fine",pt_emulated_fine);
    if (pt_emulated_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_pt_eff_fine",pt_eff_fine);
    if (pt_eff_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_pt_eff_fine not found!" << endl; not_loaded = true; }

    //Combination phi
    TH1D *phi_all = 0;
    TH1D *phi_emulated = 0;
    TH1D *phi_eff = 0;
    TH1D *phi_all_fine = 0;
    TH1D *phi_emulated_fine = 0;
    TH1D *phi_eff_fine = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_all",phi_all);
    if (phi_all == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_emulated",phi_emulated);
    if (phi_emulated == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_eff",phi_eff);
    if (phi_eff == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_all_fine",phi_all_fine);
    if (phi_all_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_emulated_fine",phi_emulated_fine);
    if (phi_emulated_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_phi_eff_fine",phi_eff_fine);
    if (phi_eff_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_phi_eff_fine not found!" << endl; not_loaded = true; }

    //Combination eta
    TH1D *eta_all = 0;
    TH1D *eta_emulated = 0;
    TH1D *eta_eff = 0;
    TH1D *eta_all_fine = 0;
    TH1D *eta_emulated_fine = 0;
    TH1D *eta_eff_fine = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_all",eta_all);
    if (eta_all == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_emulated",eta_emulated);
    if (eta_emulated == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_eff",eta_eff);
    if (eta_eff == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_all_fine",eta_all_fine);
    if (eta_all_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_all_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_emulated_fine",eta_emulated_fine);
    if (eta_emulated_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_emulated_fine not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_eff_fine",eta_eff_fine);
    if (eta_eff_fine == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_eff_fine not found!" << endl; not_loaded = true; }

    //Combination eta low pt
    TH1D *eta_all_lowpt = 0;
    TH1D *eta_emulated_lowpt = 0;
    TH1D *eta_eff_lowpt = 0;
    TH1D *eta_all_fine_lowpt = 0;
    TH1D *eta_emulated_fine_lowpt = 0;
    TH1D *eta_eff_fine_lowpt = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_all_lowpt",eta_all_lowpt);
    if (eta_all_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_all not_lowpt found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_emulated_lowpt",eta_emulated_lowpt);
    if (eta_emulated_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_emulated_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_eff_lowpt",eta_eff_lowpt);
    if (eta_eff_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_eff_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_all_fine_lowpt",eta_all_fine_lowpt);
    if (eta_all_fine_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_all_fine_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_emulated_fine_lowpt",eta_emulated_fine_lowpt);
    if (eta_emulated_fine_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_emulated_fine_lowpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_eff_fine_lowpt",eta_eff_fine_lowpt);
    if (eta_eff_fine_lowpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_eff_fine_lowpt not found!" << endl; not_loaded = true; }

    //Combination eta high pt
    TH1D *eta_all_highpt = 0;
    TH1D *eta_emulated_highpt = 0;
    TH1D *eta_eff_highpt = 0;
    TH1D *eta_all_fine_highpt = 0;
    TH1D *eta_emulated_fine_highpt = 0;
    TH1D *eta_eff_fine_highpt = 0;
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_all_highpt",eta_all_highpt);
    if (eta_all_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_all not_highpt found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_emulated_highpt",eta_emulated_highpt);
    if (eta_emulated_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_emulated_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_eff_highpt",eta_eff_highpt);
    if (eta_eff_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_eff_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_all_fine_highpt",eta_all_fine_highpt);
    if (eta_all_fine_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_all_fine_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_emulated_fine_highpt",eta_emulated_fine_highpt);
    if (eta_emulated_fine_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_emulated_fine_highpt not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_"+hist_mod+"leading_eta_eff_fine_highpt",eta_eff_fine_highpt);
    if (eta_eff_fine_highpt == 0) { cout << "ak5PF_"+hist_mod+"leading_eta_eff_fine_highpt not found!" << endl; not_loaded = true; }

    //Combination delta eta and delta phi
    TH1D *delta_phi_all = 0;
    TH1D *delta_phi_emulated = 0;
    TH1D *delta_phi_eff = 0;
    TH1D *delta_eta_all = 0;
    TH1D *delta_eta_emulated = 0;
    TH1D *delta_eta_eff = 0;
    data_in->GetObject("ak5PF_delta_phi_all",delta_phi_all);
    if (delta_phi_all == 0) { cout << "ak5PF_delta_phi_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_emulated",delta_phi_emulated);
    if (delta_phi_emulated == 0) { cout << "ak5PF_delta_phi_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_phi_eff",delta_phi_eff);
    if (delta_phi_eff == 0) { cout << "ak5PF_delta_phi_eff not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_all",delta_eta_all);
    if (delta_eta_all == 0) { cout << "ak5PF_delta_eta_all not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_emulated",delta_eta_emulated);
    if (delta_eta_emulated == 0) { cout << "ak5PF_delta_eta_emulated not found!" << endl; not_loaded = true; }
    data_in->GetObject("ak5PF_delta_eta_eff",delta_eta_eff);
    if (delta_eta_eff == 0) { cout << "ak5PF_delta_eta_eff not found!" << endl; not_loaded = true; }

    //Summary regarding the loading of the histograms
    if (detail && !not_loaded) { cout<<"All histograms loaded sucessfully!"<<endl; }
    if (not_loaded)
	{
	cout<<"Some histograms were not loaded! Please check the input rootfile. The routine will terminate now!" << endl;
	return;
	}


//ploting control distributions
//HTL_Jet15U
	plot_2histograms(pt_all_Jet15, "Triggered", pt_emulated_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_leading_pt" + sufix, "top_right", true, detail);
	plot_2histograms(pt_all_fine_Jet15, "Triggered", pt_emulated_fine_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_leading_pt_fine" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_Jet15, "Triggered", phi_emulated_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_leading_phi" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_fine_Jet15, "Triggered", phi_emulated_fine_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_leading_phi_fine" + sufix, "top_right", true, detail);
	plot_2histograms(eta_all_Jet15, "Triggered", eta_emulated_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_leading_eta" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet15, "Triggered", eta_emulated_fine_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_leading_eta_fine" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_phi_all_Jet15, "Triggered", delta_phi_emulated_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_delta_phi" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_eta_all_Jet15, "Triggered", delta_eta_emulated_Jet15, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet15U_delta_eta" + sufix, "bottom_middle", true, detail);

//HTL_Jet30U
	plot_2histograms(pt_all_Jet30, "Triggered", pt_emulated_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_pt" + sufix, "top_right", true, detail);
	plot_2histograms(pt_all_fine_Jet30, "Triggered", pt_emulated_fine_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_pt_fine" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_Jet30, "Triggered", phi_emulated_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_phi" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_fine_Jet30, "Triggered", phi_emulated_fine_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_phi_fine" + sufix, "top_right", true, detail);
	plot_2histograms(eta_all_Jet30, "Triggered", eta_emulated_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_eta" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet30, "Triggered", eta_emulated_fine_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_eta_fine" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_Jet30_lowpt, "Triggered", eta_emulated_Jet30_lowpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_eta_lowpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet30_lowpt, "Triggered", eta_emulated_fine_Jet30_lowpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_eta_fine_lowpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_Jet30_highpt, "Triggered", eta_emulated_Jet30_highpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_eta_highpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet30_highpt, "Triggered", eta_emulated_fine_Jet30_highpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_leading_eta_fine_highpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_phi_all_Jet30, "Triggered", delta_phi_emulated_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_delta_phi" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_eta_all_Jet30, "Triggered", delta_eta_emulated_Jet30, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet30U_delta_eta" + sufix, "bottom_middle", true, detail);

//HTL_Jet50U
	plot_2histograms(pt_all_Jet50, "Triggered", pt_emulated_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_pt" + sufix, "top_right", true, detail);
	plot_2histograms(pt_all_fine_Jet50, "Triggered", pt_emulated_fine_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_pt_fine" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_Jet50, "Triggered", phi_emulated_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_phi" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_fine_Jet50, "Triggered", phi_emulated_fine_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_phi_fine" + sufix, "top_right", true, detail);
	plot_2histograms(eta_all_Jet50, "Triggered", eta_emulated_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_eta" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet50, "Triggered", eta_emulated_fine_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_eta_fine" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_Jet50_lowpt, "Triggered", eta_emulated_Jet50_lowpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_eta_lowpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet50_lowpt, "Triggered", eta_emulated_fine_Jet50_lowpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_eta_fine_lowpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_Jet50_highpt, "Triggered", eta_emulated_Jet50_highpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_eta_highpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine_Jet50_highpt, "Triggered", eta_emulated_fine_Jet50_highpt, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_leading_eta_fine_highpt" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_phi_all_Jet50, "Triggered", delta_phi_emulated_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_delta_phi" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_eta_all_Jet50, "Triggered", delta_eta_emulated_Jet50, "Triggered and Emulated", output_path, "control_" + prefix + "/HLT_Jet50U_delta_eta" + sufix, "bottom_middle", true, detail);

//Combination
	plot_2histograms(pt_all, "Triggered", pt_emulated, "Triggered and Emulated", output_path, "control_" + prefix + "/leading_pt" + sufix, "top_right", true, detail);
	plot_2histograms(pt_all_fine, "Triggered", pt_emulated_fine, "Triggered and Emulated", output_path, "control_" + prefix + "/leading_pt_fine" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all, "Triggered", phi_emulated, "Triggered and Emulated", output_path, "control_" + prefix + "/leading_phi" + sufix, "top_right", true, detail);
	plot_2histograms(phi_all_fine, "Triggered", phi_emulated_fine, "Triggered and Emulated", output_path, "control_" + prefix + "/leading_phi_fine" + sufix, "top_right", true, detail);
	plot_2histograms(eta_all, "Triggered", eta_emulated, "Triggered and Emulated", output_path, "control_" + prefix + "/leading_eta" + sufix, "bottom_middle", true, detail);
	plot_2histograms(eta_all_fine, "Triggered", eta_emulated_fine, "Triggered and Emulated", output_path, "control_" + prefix + "/leading_eta_fine" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_phi_all, "Triggered", delta_phi_emulated, "Triggered and Emulated", output_path, "control_" + prefix + "/delta_phi" + sufix, "bottom_middle", true, detail);
	plot_2histograms(delta_eta_all, "Triggered", delta_eta_emulated, "Triggered and Emulated", output_path, "control_" + prefix + "/delta_eta" + sufix, "bottom_middle", true, detail);

//plot efficiencies
//HLT_Jet15U
	plot_efficiency(pt_eff_Jet15, prefix, "/HLT_Jet15U_leading_pt" + sufix, output_path, "top_left", detail);
	plot_efficiency(pt_eff_fine_Jet15, prefix, "/HLT_Jet15U_leading_pt_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_Jet15, prefix, "/HLT_Jet15U_leading_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_fine_Jet15, prefix, "/HLT_Jet15U_leading_phi_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_Jet15, prefix, "/HLT_Jet15U_leading_eta" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_fine_Jet15, prefix, "/HLT_Jet15U_leading_eta_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_phi_eff_Jet15, prefix, "/HLT_Jet15U_delta_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_eta_eff_Jet15, prefix, "/HLT_Jet15U_delta_eta" + sufix, output_path, "top_left", detail);

//HLT_Jet30U
	plot_efficiency(pt_eff_Jet30, prefix, "/HLT_Jet30U_leading_pt" + sufix, output_path, "top_left", detail);
	plot_efficiency(pt_eff_fine_Jet30, prefix, "/HLT_Jet30U_leading_pt_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_Jet30, prefix, "/HLT_Jet30U_leading_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_fine_Jet30, prefix, "/HLT_Jet30U_leading_phi_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_Jet30, prefix, "/HLT_Jet30U_leading_eta" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_fine_Jet30, prefix, "/HLT_Jet30U_leading_eta_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_phi_eff_Jet30, prefix, "/HLT_Jet30U_delta_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_eta_eff_Jet30, prefix, "/HLT_Jet30U_delta_eta" + sufix, output_path, "top_left", detail);

//HLT_Jet50U
	plot_efficiency(pt_eff_Jet50, prefix, "/HLT_Jet50U_leading_pt" + sufix, output_path, "top_left", detail);
	plot_efficiency(pt_eff_fine_Jet50, prefix, "/HLT_Jet50U_leading_pt_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_Jet50, prefix, "/HLT_Jet50U_leading_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_fine_Jet50, prefix, "/HLT_Jet50U_leading_phi_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_Jet50, prefix, "/HLT_Jet50U_leading_eta" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_fine_Jet50, prefix, "/HLT_Jet50U_leading_eta_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_phi_eff_Jet50, prefix, "/HLT_Jet50U_delta_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_eta_eff_Jet50, prefix, "/HLT_Jet50U_delta_eta" + sufix, output_path, "top_left", detail);

//Combination
	plot_efficiency(pt_eff, prefix, "/leading_pt" + sufix, output_path, "top_left", detail);
	plot_efficiency(pt_eff_fine, prefix, "/leading_pt_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff, prefix, "/leading_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(phi_eff_fine, prefix, "/leading_phi_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff, prefix, "/leading_eta" + sufix, output_path, "top_left", detail);
	plot_efficiency(eta_eff_fine, prefix, "/leading_eta_fine" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_phi_eff, prefix, "/delta_phi" + sufix, output_path, "top_left", detail);
	plot_efficiency(delta_eta_eff, prefix, "/delta_eta" + sufix, output_path, "top_left", detail);


//Fitting
	double fit_eff_Jet30[9*2];
	double fit_eff_fine_Jet30[9*2];
	double fit_eff_Jet50[9*2];
	double fit_eff_fine_Jet50[9*2];

	for (int i=0; i < 9*2; i++)
		{
		fit_eff_Jet30[i] = 0.0;
		fit_eff_fine_Jet30[i] = 0.0;
		fit_eff_Jet50[i] = 0.0;
		fit_eff_fine_Jet50[i] = 0.0;
		}
	
	char *fit_function;
	double range_min30 = 0.0, range_max30 = 0.0, range_min50 = 0.0, range_max50 = 0.0;
	if (sufix == "_v1" && !issj)
		{
		fit_function = "pol8";
		range_min30 = 34.0;
		range_min50 = 36.0;
		range_max30 = 78.0;
		range_max50 = 110.0;
		}
	if (sufix == "_v1" && issj)
		{
		fit_function = "pol8";
		range_min30 = 34.0;
		range_min50 = 34.0;
		range_max30 = 80.0;
		range_max50 = 116.0;
		}
	if (sufix == "_v2" && !issj)
		{
		fit_function = "pol8";
		range_min30 = 34.0;
		range_min50 = 36.0;
		range_max30 = 78.0;
		range_max50 = 110.0;
		}
	if (sufix == "_v2" && issj)
		{
		fit_function = "pol8";
		range_min30 = 34.0;
		range_min50 = 34.0;
		range_max30 = 80.0;
		range_max50 = 116.0;
		}

	//fit_and_plot(pt_eff_Jet30, fit_eff_Jet30, fit_function, 36., 80., prefix, "_Jet30U_leading_pt_fitted" + sufix, output_path, "top_left", detail);
	fit_and_plot(pt_eff_fine_Jet30, fit_eff_fine_Jet30, fit_function, range_min30, range_max30, prefix, "_Jet30U_leading_pt_fine_fitted" + sufix, output_path, "top_left", detail);
	//fit_and_plot(pt_eff_Jet50, fit_eff_Jet50, fit_function, 36., 118., prefix, "_Jet50U_leading_pt_fitted" + sufix, output_path, "top_left", detail);
	fit_and_plot(pt_eff_fine_Jet50, fit_eff_fine_Jet50, fit_function, range_min50, range_max50, prefix, "_Jet50U_leading_pt_fine_fitted" + sufix, output_path, "top_left", detail);

	FILE *file;
   	file = fopen(file_out.c_str(),"w");

	for (int i=0; i <=8; i++)
		{
   		fprintf(file,"%2.17f %2.17f %2.17f %2.17f\n",fit_eff_fine_Jet30[i],fit_eff_fine_Jet30[9+i],fit_eff_fine_Jet50[i],fit_eff_fine_Jet50[9+i]);
		}

   	fclose (file);
   
//close all TFiles
    	if (detail) { cout<<"Closing the files and deleting variables..."<<endl; }
    	data_in->Close();
    	if (detail) { cout<<"Done!"<<endl; }
}
