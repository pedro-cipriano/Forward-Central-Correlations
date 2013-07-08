// Pedro Cipriano, Jun 2013
// DESY, CMS
// Last Update: 17 Jun 2012
//
//


#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "common_methods.h"



void compare_unfolding_results(string data_file, string standard_bbb, string roounfold_bbb, string tunfold, string svd, string bayes, string plots_path, string prefix, bool isMC = false, bool detail = false)
{

//output the configuration
   if (detail) { cout<<"Compare Unfolding Results Configuration"<<endl; }
   if (detail) { cout<<"Uncorrected Data File :             "<<data_file<<endl; }
   if (detail) { cout<<"Standard Bin-By-Bin Results File :  "<<standard_bbb<<endl; }
   if (detail) { cout<<"RooUnfold Bin-By-Bin Results File : "<<roounfold_bbb<<endl; }
   if (detail) { cout<<"TUnfold Results File :              "<<tunfold<<endl; }
   if (detail) { cout<<"SVD Results File :                  "<<svd<<endl; }
   if (detail) { cout<<"Bayes Results File :                "<<bayes<<endl; }
   if (detail) { cout<<"Plots Path :                        "<<plots_path<<endl; }
   if (detail) { cout<<"Prefix :                            "<<prefix<<endl; }
   if (detail) { cout<<"Detail level :                      "<<detail<<endl; }

    if (detail) { cout << "Opening files..." << endl; }
//opens the MC files
    TFile *data = new TFile( data_file.c_str() );
    TFile *sbbb = new TFile( standard_bbb.c_str() );
    TFile *rbbb = new TFile( roounfold_bbb.c_str() );
    TFile *tunf = new TFile( tunfold.c_str() );
    TFile *svdf = new TFile( svd.c_str() );
    TFile *baye = new TFile( bayes.c_str() );

    TString load_prefix = "ak5PF_";
    if (isMC) { load_prefix = "ak5Gen_"; }

//plot delta phi unfolding results
    if (detail) { cout << "Plotting Delta Phi Unfolding Results..." << endl; }

    TH1D *data_delta_phi = 0;
    data->GetObject("ak5PF_delta_phi",data_delta_phi);
    if (data_delta_phi == 0) { cout << "ak5PF_delta_phi not found!" << endl; return; }
    TH1D *sbbb_delta_phi = 0;
    sbbb->GetObject(load_prefix+"delta_phi",sbbb_delta_phi);
    if (sbbb_delta_phi == 0) { cout << load_prefix << "delta_phi not found!" << endl; return; }
    TH1D *rbbb_delta_phi = 0;
    rbbb->GetObject("output_true_delta_phi",rbbb_delta_phi);
    if (rbbb_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }
    TH1D *tunf_delta_phi = 0;
    tunf->GetObject("output_true_delta_phi",tunf_delta_phi);
    if (tunf_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }
    TH1D *svd_delta_phi = 0;
    svdf->GetObject("output_true_delta_phi",svd_delta_phi);
    if (svd_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }
    TH1D *baye_delta_phi = 0;
    baye->GetObject("output_true_delta_phi",baye_delta_phi);
    if (baye_delta_phi == 0) { cout << "resp_delta_phi not found!" << endl; return; }



    plot_six_dist(data_delta_phi, "Uncorrected Data", sbbb_delta_phi, "Standard Bin-By-Bin", rbbb_delta_phi, "RooUnfold Bin-By-Bin", tunf_delta_phi, "TUnfold", svd_delta_phi, "SVD", baye_delta_phi, "Bayes", plots_path, prefix, "delta_phi", "top_left", detail);    
    TH1D *ratio1;
    ratio1 = static_cast<TH1D*>(data_delta_phi->Clone());
    ratio1->Divide(rbbb_delta_phi,sbbb_delta_phi,1.,1.,"B");
    TH1D *ratio2;
    ratio2 = static_cast<TH1D*>(data_delta_phi->Clone());
    ratio2->Divide(tunf_delta_phi,sbbb_delta_phi,1.,1.,"B");
    TH1D *ratio3;
    ratio3 = static_cast<TH1D*>(data_delta_phi->Clone());
    ratio3->Divide(svd_delta_phi,sbbb_delta_phi,1.,1.,"B");
    TH1D *ratio4;
    ratio4 = static_cast<TH1D*>(data_delta_phi->Clone());
    ratio4->Divide(baye_delta_phi,sbbb_delta_phi,1.,1.,"B");
    plot_4histograms(ratio1, "Bin-By-Bin - RooUnfold/Standard", ratio2, "TUnfold/Bin-By-Bin", ratio3, "SVD/Bin-By-Bin", ratio4, "Bayes/Bin-By-Bin", plots_path, prefix+"delta_phi_ratios", "top_right", false, detail);
    //plot_3histograms(ratio1, "Bin-By-Bin - RooUnfold/Standard", ratio3, "SVD/Bin-By-Bin", ratio4, "Bayes/Bin-By-Bin", plots_path, prefix+"delta_phi_ratios", "top_right", false, detail);


}
