// Pedro Cipriano, Dec 2011
// DESY, CMS
// Last Update: 13 Dec 2012
//
// get_pileup_normalization(string *data_in, string data_out, double *data_lumi, int n_files, string vertex_weights = "", TString input_label = "", bool detail = false, bool test = false)
// get the pileup normalization factor

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
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

void get_pileup_normalization(string *data_in, string data_out, double *data_lumi, int n_files, string vertex_weights = "", TString input_label = "", bool detail = false, bool test = false)
{

//output the configuration
   if (detail) { cout<<"Get Pileup Normalization Configuration"<<endl; }
   if (detail) { cout<<"Number of files :           "<<n_files<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :               "<<data_out<<endl; }
   if (detail) { cout<<"Vertex Weights :            "<<vertex_weights<<endl; }
   if (detail) { cout<<"Vertex Weights Input Label: "<<input_label<<endl; }
   if (detail) { cout<<"Detail level :              "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :                 "<<test<<endl; }

//opening the vertex weight file
    bool vertex_weighted = false;
    TH1D *vertex_hist = 0;
    if (vertex_weights != "")
    {
    if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
    TFile *vertex_file = new TFile( vertex_weights.c_str() );
    vertex_file->GetObject("vertex_weights"+input_label,vertex_hist);
    if (vertex_hist != 0) { vertex_weighted = true; }
    }

// declaring and initializing histograms
   TH1D *hist_pileup_normalization_all;
   TH1D *hist_pileup_normalization_selected;

   hist_pileup_normalization_all =  new TH1D("pileup_normalization_all","Pileup normalization in all events;Type;Value", 5, 0, 5);
   hist_pileup_normalization_selected =  new TH1D("pileup_normalization_selected","Pileup normalization in selected events;Type;Value", 8, 0, 8);

   hist_pileup_normalization_all->Sumw2();
   hist_pileup_normalization_all->Fill("Events all",0);
   hist_pileup_normalization_all->Fill("Weights all",0);
   hist_pileup_normalization_all->Fill("Events nopu",0);
   hist_pileup_normalization_all->Fill("Weights nopu",0);
   hist_pileup_normalization_all->Fill("Normalization",0);

   hist_pileup_normalization_selected->Sumw2();
   hist_pileup_normalization_selected->Fill("Events all",0);
   hist_pileup_normalization_selected->Fill("Weights all",0);
   hist_pileup_normalization_selected->Fill("Events nopu",0);
   hist_pileup_normalization_selected->Fill("Weights nopu",0);
   hist_pileup_normalization_selected->Fill("Normalization",0);
   hist_pileup_normalization_selected->Fill("Jets all",0);
   hist_pileup_normalization_selected->Fill("Jets selected",0);
   hist_pileup_normalization_selected->Fill("Jet Ratio",0);

//declaring variables
  double pt_min = 35.0;
  int counter_all = 0, counter_pv = 0, counter_all2 = 0, counter_pv2 = 0;
  double counter_wall = 0.0, counter_wpv = 0.0, counter_wall2 = 0.0, counter_wpv2 = 0.0;
  int counter_jet = 0, jet_tmp, counter_jet2 = 0;
  int nentries = 0;
  double pu_scale = 0.0, pu_scale2 = 0.0, jet_ratio = 0.0;
  double prescale = 0.0;
  double vertex_factor = 1.0;
  double pt = 0.0, eta = 0.0;
  double pt_central, eta_central, pt_forward, eta_forward;

//loop over the files
for (int z = 0; z < n_files; z++)
{
//open the file
  string file = data_in[z];
  cout<<z+1<<"/"<<n_files<<" -> "<<file<<endl;
  TFile *inf = TFile::Open( file.c_str() );

//define the tree branch
  TTree *tr = (TTree*)inf->Get("ak5/ProcessedTree");
  QCDEvent *Event = new QCDEvent();
  TBranch *branch = tr->GetBranch("events");
  branch->SetAddress(&Event);

  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  prescale = data_lumi[z];
  cout<<"Cross-section per event : "<<prescale<<endl;

  if (test) { nentries = 250; }

  int decade = 0;

   for (Int_t i=0;i<nentries;i++) {
     double progress = 10.0*i/(1.0*nentries);
     int k = TMath::FloorNint(progress); 
     if (k > decade) 
      cout<<k*10<<" % "<<endl;
     decade = k; 
     
     tr->GetEntry(i);
     if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1); }

     counter_all++;
     counter_wall = counter_wall + vertex_factor * prescale;
     if (test) { cout<<"Npu : "<<Event->evtHdr().pu()<<endl; }
     if (test) { cout<<"counter all : "<<counter_wall<<endl; }

      if (Event->evtHdr().pu() == 0)
	{
	counter_pv++;
	counter_wpv = counter_wpv + vertex_factor * prescale;
	if (test) { cout<<"counter nopu : "<<counter_wpv<<endl; }
	}

     pt_central = 0.0;
     eta_central = 0.0;
     pt_forward = 0.0;
     eta_forward = 0.0;
     jet_tmp = 0;

    for(unsigned int j=0; j<Event->nGenJets(); j++)
	{
	pt = Event->genjet(j).pt();
	eta = Event->genjet(j).eta();
	if (pt >= pt_min)
		{
		counter_jet++;
		jet_tmp++;
		if (eta <= 2.8 && eta >= -2.8 && pt > pt_central)
			{ pt_central = pt; eta_central = eta; }
		if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward)
			{ pt_forward = pt; eta_forward = eta; }
     		}
	}
     
     if (pt_forward > pt_min && pt_central > pt_min)
     {
     counter_all2++;
     counter_wall2 = counter_wall2 + prescale*vertex_factor;
      if (Event->evtHdr().pu() == 0)
	{
	counter_jet2 = counter_jet2 + jet_tmp;
	counter_pv2++;
	counter_wpv2 = counter_wpv2 + vertex_factor * prescale;
	}
     }

   }

}

  cout<<endl<<endl;
  cout<<"Events read:                               "<<counter_all<<endl;
  cout<<"Events without pileup:                     "<<counter_pv<<endl;
  cout<<"Weights of Events read:                    "<<counter_wall<<endl;
  cout<<"Weights of Events without pileup:          "<<counter_wpv<<endl;
  cout<<"Events selected:                           "<<counter_all2<<endl;
  cout<<"Events selected without pileup:            "<<counter_pv2<<endl;
  cout<<"Weights of Events selected:                "<<counter_wall2<<endl;
  cout<<"Weights of Events selected without pileup: "<<counter_wpv2<<endl;
  cout<<"Total Jets:                                "<<counter_jet<<endl;
  cout<<"Jets in selected events:                   "<<counter_jet2<<endl;

  pu_scale = (double) counter_wall/ (double) counter_wpv;
  pu_scale2 = (double) counter_wall2/ (double) counter_wpv2;
  jet_ratio = (double) counter_jet/ (double) counter_jet2;
  cout<<endl;
  cout<<"Pile up correction (all events) =          "<<pu_scale<<endl;
  cout<<"Pile up correction (selected events) =     "<<pu_scale2<<endl;
  cout<<"Pile up correction (from jets) =           "<<jet_ratio<<endl;
  cout<<endl;

//fill the pileup_normalization_all histogram
     hist_pileup_normalization_all->SetBinContent(1,counter_all);
     hist_pileup_normalization_all->SetBinContent(2,counter_wall);
     hist_pileup_normalization_all->SetBinContent(3,counter_pv);
     hist_pileup_normalization_all->SetBinContent(4,counter_wpv);
     hist_pileup_normalization_all->SetBinContent(5,pu_scale);

//fill the pileup_normalization_selected histogram
     hist_pileup_normalization_selected->SetBinContent(1,counter_all2);
     hist_pileup_normalization_selected->SetBinContent(2,counter_wall2);
     hist_pileup_normalization_selected->SetBinContent(3,counter_pv2);
     hist_pileup_normalization_selected->SetBinContent(4,counter_wpv2);
     hist_pileup_normalization_selected->SetBinContent(5,pu_scale2);
     hist_pileup_normalization_selected->SetBinContent(6,counter_jet);
     hist_pileup_normalization_selected->SetBinContent(7,counter_jet2);
     hist_pileup_normalization_selected->SetBinContent(8,jet_ratio);

     //Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     //write histograms on disc
     hist_pileup_normalization_all->Write();
     hist_pileup_normalization_selected->Write();

     //close the output file
     data_output->Close();

     //delete the histograms to avoid memory leak
     delete(hist_pileup_normalization_all);
     delete(hist_pileup_normalization_selected);

}


void show_pileup_normalization(string p6_z2, string p6_z2_jetmettau, string p6_z2_jetmet, string p6_z2_jet, string p8_4c, string p8_4c_jetmettau, string p8_4c_jetmet, string p8_4c_jet, bool detail = false)
{


//output the configuration
if (detail)
	{
	cout<<"Get Pileup Normalization Configuration"<<endl;
	cout<<"Pythia6 - Tune Z2* :                "<<p6_z2<<endl;
	cout<<"Pythia6 - Tune Z2* JetMETTau_2010A: "<<p6_z2_jetmettau<<endl;
	cout<<"Pythia6 - Tune Z2* JetMET_2010A:    "<<p6_z2_jetmet<<endl;
	cout<<"Pythia6 - Tune Z2* Jet_2010B:       "<<p6_z2_jet<<endl;
	cout<<"Pythia8 - Tune 4C :                 "<<p8_4c<<endl;
	cout<<"Pythia8 - Tune 4C JetMETTau_2010A:  "<<p8_4c_jetmettau<<endl;
	cout<<"Pythia8 - Tune 4C JetMET_2010A:     "<<p8_4c_jetmet<<endl;
	cout<<"Pythia8 - Tune 4C Jet_2010B:        "<<p8_4c_jet<<endl;
	cout<<"Detail level :              "<<detail<<endl;
	}

//opening the files and extracting the histograms
if (detail) { cout<<"Opening the files and extracting the histograms..."<<endl; }

    //p6_z2 weights
    TH1D *weight_p6_z2_all = 0;
    TH1D *weight_p6_z2_selected = 0;
    TFile *p6_z2_file = new TFile( p6_z2.c_str() );
    p6_z2_file->GetObject("pileup_normalization_all",weight_p6_z2_all);
    p6_z2_file->GetObject("pileup_normalization_selected",weight_p6_z2_selected);
    if (weight_p6_z2_all == 0 or weight_p6_z2_selected == 0) { cout << "The file " << p6_z2 << " dont contain valid histograms!" << endl; return; }

    //p6_z2 jetmettau weights
    TH1D *weight_p6_z2_jetmettau_all = 0;
    TH1D *weight_p6_z2_jetmettau_selected = 0;
    TFile *p6_z2_jetmettau_file = new TFile( p6_z2_jetmettau.c_str() );
    p6_z2_jetmettau_file->GetObject("pileup_normalization_all",weight_p6_z2_jetmettau_all);
    p6_z2_jetmettau_file->GetObject("pileup_normalization_selected",weight_p6_z2_jetmettau_selected);
    if (weight_p6_z2_jetmettau_all == 0 or weight_p6_z2_jetmettau_selected == 0) { cout << "The file " << p6_z2_jetmettau << " dont contain valid histograms!" << endl; return; }

    //p6_z2 jetmet weights
    TH1D *weight_p6_z2_jetmet_all = 0;
    TH1D *weight_p6_z2_jetmet_selected = 0;
    TFile *p6_z2_jetmet_file = new TFile( p6_z2_jetmet.c_str() );
    p6_z2_jetmet_file->GetObject("pileup_normalization_all",weight_p6_z2_jetmet_all);
    p6_z2_jetmet_file->GetObject("pileup_normalization_selected",weight_p6_z2_jetmet_selected);
    if (weight_p6_z2_jetmet_all == 0 or weight_p6_z2_jetmet_selected == 0) { cout << "The file " << p6_z2_jetmet << " dont contain valid histograms!" << endl; return; }

    //p6_z2 jet weights
    TH1D *weight_p6_z2_jet_all = 0;
    TH1D *weight_p6_z2_jet_selected = 0;
    TFile *p6_z2_jet_file = new TFile( p6_z2_jet.c_str() );
    p6_z2_jet_file->GetObject("pileup_normalization_all",weight_p6_z2_jet_all);
    p6_z2_jet_file->GetObject("pileup_normalization_selected",weight_p6_z2_jet_selected);
    if (weight_p6_z2_jet_all == 0 or weight_p6_z2_jet_selected == 0) { cout << "The file " << p6_z2_jet << " dont contain valid histograms!" << endl; return; }

    //p8_4c weights
    TH1D *weight_p8_4c_all = 0;
    TH1D *weight_p8_4c_selected = 0;
    TFile *p8_4c_file = new TFile( p8_4c.c_str() );
    p8_4c_file->GetObject("pileup_normalization_all",weight_p8_4c_all);
    p8_4c_file->GetObject("pileup_normalization_selected",weight_p8_4c_selected);
    if (weight_p8_4c_all == 0 or weight_p8_4c_selected == 0) { cout << "The file " << p8_4c << " dont contain valid histograms!" << endl; return; }

    //p8_4c jetmettau weights
    TH1D *weight_p8_4c_jetmettau_all = 0;
    TH1D *weight_p8_4c_jetmettau_selected = 0;
    TFile *p8_4c_jetmettau_file = new TFile( p8_4c_jetmettau.c_str() );
    p8_4c_jetmettau_file->GetObject("pileup_normalization_all",weight_p8_4c_jetmettau_all);
    p8_4c_jetmettau_file->GetObject("pileup_normalization_selected",weight_p8_4c_jetmettau_selected);
    if (weight_p8_4c_jetmettau_all == 0 or weight_p8_4c_jetmettau_selected == 0) { cout << "The file " << p8_4c_jetmettau << " dont contain valid histograms!" << endl; return; }

    //p8_4c jetmet weights
    TH1D *weight_p8_4c_jetmet_all = 0;
    TH1D *weight_p8_4c_jetmet_selected = 0;
    TFile *p8_4c_jetmet_file = new TFile( p8_4c_jetmet.c_str() );
    p8_4c_jetmet_file->GetObject("pileup_normalization_all",weight_p8_4c_jetmet_all);
    p8_4c_jetmet_file->GetObject("pileup_normalization_selected",weight_p8_4c_jetmet_selected);
    if (weight_p8_4c_jetmet_all == 0 or weight_p8_4c_jetmet_selected == 0) { cout << "The file " << p8_4c_jetmet << " dont contain valid histograms!" << endl; return; }

    //p8_4c jet weights
    TH1D *weight_p8_4c_jet_all = 0;
    TH1D *weight_p8_4c_jet_selected = 0;
    TFile *p8_4c_jet_file = new TFile( p8_4c_jet.c_str() );
    p8_4c_jet_file->GetObject("pileup_normalization_all",weight_p8_4c_jet_all);
    p8_4c_jet_file->GetObject("pileup_normalization_selected",weight_p8_4c_jet_selected);
    if (weight_p8_4c_jet_all == 0 or weight_p8_4c_jet_selected == 0) { cout << "The file " << p8_4c_jet << " dont contain valid histograms!" << endl; return; }

//outputing the results
if (detail) { cout<<"Showing the normalization factors..."<<endl; }

cout<<"Pileup Normalization"<<endl;
cout<<"Dataset                            All Events Selected Events"<<endl;
cout<<"Pythia6 - Tune Z2*                 "<<setw(6)<<weight_p6_z2_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p6_z2_selected->GetBinContent(5)<<endl;
cout<<"Pythia6 - Tune Z2* JetMETTau_2010A "<<setw(6)<<weight_p6_z2_jetmettau_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p6_z2_jetmettau_selected->GetBinContent(5)<<endl;
cout<<"Pythia6 - Tune Z2* JetMET_2010A    "<<setw(6)<<weight_p6_z2_jetmet_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p6_z2_jetmet_selected->GetBinContent(5)<<endl;
cout<<"Pythia6 - Tune Z2* Jet_2010B       "<<setw(6)<<weight_p6_z2_jet_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p6_z2_jet_selected->GetBinContent(5)<<endl;
cout<<"Pythia8 - Tune 4C                  "<<setw(6)<<weight_p8_4c_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p8_4c_selected->GetBinContent(5)<<endl;
cout<<"Pythia8 - Tune 4C JetMETTau_2010A  "<<setw(6)<<weight_p8_4c_jetmettau_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p8_4c_jetmettau_selected->GetBinContent(5)<<endl;
cout<<"Pythia8 - Tune 4C JetMET_2010A     "<<setw(6)<<weight_p8_4c_jetmet_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p8_4c_jetmet_selected->GetBinContent(5)<<endl;
cout<<"Pythia8 - Tune 4C Jet_2010B        "<<setw(6)<<weight_p8_4c_jet_all->GetBinContent(5)<<"    "<<setw(6)<<weight_p8_4c_jet_selected->GetBinContent(5)<<endl;

if (detail) { cout<<"Output finished..."<<endl; }
}
