// Pedro Cipriano, Nov 2012
// DESY, CMS
// Last Update: 12 Dec 2012
//
// get_vertex_distribution(string *data_in, string data_out, double *data_lumi, int n_files, string data_type = "DATA", string sel_mode = "1vertex", string vertex_weights = "", TString input_label = "", bool detail = false, bool test = false)
// gets the vertex distribution

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

void get_vertex_distribution(string *data_in, string data_out, double *data_lumi, int n_files, string data_type = "DATA", string sel_mode = "1vertex", string vertex_weights = "", TString input_label = "", bool detail = false, bool test = false)
{

//output the configuration
   if (detail) { cout<<"Get Vertex Distribution Configuration"<<endl; }
   if (detail) { cout<<"Data Type :                 "<<data_type<<endl; }
   if (detail) { cout<<"Selection mode :            "<<sel_mode<<endl; }
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
    if ((vertex_weights != "") and (data_type == "MC_DET" or data_type == "MC_GEN"))
    {
    if (detail) { cout<<"Opening Vertex Weights... "<<endl; }
    TFile *vertex_file = new TFile( vertex_weights.c_str() );
    vertex_file->GetObject("vertex_weights"+input_label,vertex_hist);
    if (vertex_hist != 0) { vertex_weighted = true; }
    }

//setting the prefix
   TString prefix;
   if (data_type == "DATA" or data_type == "MC_DET") { prefix = "ak5PF_"; }
   if (data_type == "MC_GEN") { prefix = "ak5Gen_"; }
   if (detail) { cout<<"Prefix :          "<<prefix<<endl; }

//main cuts
   double pt_min = 35.0;

   TH1D *hist_vertex_selected;
   TH1D *hist_pu_selected;
   TH1D *hist_pvz_selected;
   TH1D *hist_events;

   hist_vertex_selected =  new TH1D(prefix+"vertex_selected","Vertex multiplicity in selected events;Multiplicity;#frac{d#sigma}{dN} [pb]", 15, 0, 15);
   hist_pu_selected =  new TH1D(prefix+"pu_selected","Pileup in selected events;Pileup;#frac{d#sigma}{dPileup} [pb]", 15, 0, 15);
   hist_pvz_selected =  new TH1D(prefix+"pvz_selected","z of the Primary Vertex in selected events;z-position;#frac{d#sigma}{dz} [pb]", 200, -100, 100);
   hist_events =  new TH1D(prefix+"events","Selection Chain;Type;Events or Weight", 6, 0, 6);

   hist_vertex_selected->Sumw2();
   hist_pvz_selected->Sumw2();
   hist_pu_selected->Sumw2();

   hist_events->Sumw2();
   hist_events->Fill("Total Events",0);
   hist_events->Fill("PV Selection",0);
   hist_events->Fill("Trigger Selection",0);
   hist_events->Fill("Selected",0);
   hist_events->Fill("Selected Weights",0);
   hist_events->Fill("Total Jets",0);

//declaring variables
  int counter_hlt = 0, counter_pv = 0, counter_jet = 0, counter_selected = 0, counter_entries = 0;
  double counter_weigth = 0.0;
  int nentries = 0;
  double pu_scale = 0.0;
  double prescale = 0.0, prescale_prov = 0.0;
  bool pv_pass = false, jet_tight = false;
  unsigned int number_of_jets = 1;
  double pt = 0.0, eta = 0.0, phi = 0.0, gen_pt = 0.0, gen_eta = 0.0, gen_phi = 0.0;
  double pt_central = 0.0, eta_central = 0.0, pt_forward = 0.0, eta_forward = 0.0, leading_pt = 0.0;
  double vertex_factor = 1.0;

//setup for the special case with all triggers
// multiple triggers
  int number_triggers = 3;
  string trigger_list[number_triggers];
  int hlt_list[number_triggers];
//  double prescales[number_triggers];
  bool hltPass(false);
//  bool trig_pass[number_triggers];
  trigger_list[0] = "HLT_Jet15U";
  trigger_list[1] = "HLT_Jet30U";
  trigger_list[2] = "HLT_Jet50U";
//  trigger_list[3] = "HLT_DiJetAve15U_8E29";
//  trigger_list[4] = "HLT_DiJetAve30U_8E29";
//  trigger_list[5] = "HLT_DoubleJet15U_ForwardBackward";
//  trigger_list[6] = "HLT_FwdJet20U";

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

//get the trigger information
  TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
  cout<<"Finding trigger mapping: "<<endl;

  if (data_type == "DATA")
  {
  for (int trg = 0; trg < number_triggers; trg++)
	{
	for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
		TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));
		if (ss == trigger_list[trg]) {
			hlt_list[trg] = ibin;
			continue;
			}
		}
  	if (hlt_list[trg] == -1) {
		cout<<"The requested trigger ("<<trigger_list[trg]<<") is not found." << endl;
    		//break;
  		}
	else {
		cout<<trigger_list[trg]<<" --> "<<hlt_list[trg]<<endl;
		}
	}
  }

  prescale = 0.0;
  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  if (data_type == "MC_DET" or data_type == "MC_GEN")
  {
  prescale = data_lumi[z];
  cout<<"Cross-section per event : "<<prescale<<endl;
  }
  if (data_type == "DATA")
  {
  prescale_prov =  1.0 / data_lumi[z];
  cout<<"Luminosity for " << trigger_list[z] << " = " << data_lumi[z] << endl;
  cout<<"Cross-section per event : "<<prescale_prov<<endl;
  }

  //if (test) { nentries = 100; }

  int decade = 0;

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
      if ((sel_mode == "nopileup" ) && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().pu() == 0)    	  { pv_pass = true; }
      if ((sel_mode == "allvertex") && Event->evtHdr().isPVgood() == 1)                                           { pv_pass = true; }
      if ((sel_mode == "1vertex" ) && Event->evtHdr().isPVgood() == 1 && Event->evtHdr().nVtxGood() == 1)    	  { pv_pass = true; }
      
      if (pv_pass == true)
      {
        counter_pv++;
        if (vertex_weighted) { vertex_factor = vertex_hist->GetBinContent(Event->evtHdr().pu() + 1); }

     if (data_type == "DATA")
     {
     // multiple triggers
     hltPass = false;
     prescale = 0.0;
	/*
	for (int trg = 0; trg < 3; trg++)
		{
		prescales[trg] = 1.0;
		trig_pass[trg] = false;
		if (hlt_list[trg] == -1)
			hltPass = true; // no trigger set
        	else
			{
			if (Event->fired(hlt_list[trg]) > 0)
				{
				hltPass = true;
				trig_pass[trg] = true;
				//prescales[trg] = Event->preL1(hlt_list[trg]) * Event->preHLT(hlt_list[trg]);
				}
			}
		}
	*/
     
	if (Event->fired(hlt_list[z]) > 0) { hltPass = true; }
	//cout << "z = " << z << " and trigger = " << hlt_list[z] << " and " << endl;
     }


     if (hltPass and data_type == "DATA")
     	{
     	leading_pt = 0.0;

     	for(unsigned int j=0; j<Event->nPFJets(); j++)
        	{
     		pt = Event->pfjet(j).ptCor();
     		if (leading_pt < pt && Event->pfjet(j).tightID()) { leading_pt = pt; }
     		}

	if (z == 0 and leading_pt > 35.0 and leading_pt < 70.0) { prescale = prescale_prov; }
	if (z == 1 and leading_pt > 70.0 and leading_pt < 110.0) { prescale = prescale_prov; }
	if (z == 2 and leading_pt > 110.0) { prescale = prescale_prov; }

     	//if (leading_pt > 35.0 and leading_pt < 70.0 and trig_pass[0] ) { prescale = prescales[0]; }
     	//if (leading_pt > 70.0 and leading_pt < 110.0 and trig_pass[1] ) { prescale = prescales[1]; }
     	//if (leading_pt > 110.0  and trig_pass[2] ) { prescale = prescales[2];}

     	//cout << "Event : " << i << endl;
     	//cout << "Prescales : " << prescales[0] << " " << prescales[1] << " " << prescales[2] << endl;
     	//cout << "Leading Jet pT : " << leading_pt << " prescale : " << prescale << endl;

// old combination

//     for(int l=0; l<number_triggers; l++)
//	{//
//	probs[l] = 1.0;
//	if (prescales[l] > 0 and trig_pass[l]) { probs[l] = (1 - 1/prescales[l]); }
	//if (prescales[k] > 0) { cout<<"k = "<<k<<" | prob = "<<probs[k]<<" | prescale = "<<prescales[k]<<endl; }
//	}
//     prescale = 1.0 / ( 1.0 - (probs[0] * probs[1] * probs[2]));
     //prescale = 1.0 / ( 1 - (probs[0] * probs[1] * probs[2] * probs[3] * probs[4] * probs[5] * probs[6]));
     //if (prescale > 1) { cout<<"probs = ["<<probs[0]<<" ; "<<probs[1]<<" ; "<<probs[2]<<" ; "<<probs[3]<<" ; "<<probs[4]<<"]"<<endl;   
     //if (prescale > 1) {cout<<"Prescale = "<<prescale<<endl; }

     	}

     if (prescale > 0.0)
	{ //or prescale > 0 hltPass_monitor
        counter_hlt++;
      	//cout<<"Prescale = "<<prescale<<endl;

     pt_central = 0.0;
     eta_central = 0.0;
     pt_forward = 0.0;
     eta_forward = 0.0;

     if (data_type == "DATA" or data_type == "MC_DET") { number_of_jets = Event->nPFJets(); }
     if (data_type == "MC_GEN") { number_of_jets = Event->nGenJets(); }

    for(unsigned int j=0; j<number_of_jets; j++) {
	if (data_type == "DATA" or data_type == "MC_DET")
		{
        	jet_tight = Event->pfjet(j).tightID();
		pt = Event->pfjet(j).ptCor();
		eta = Event->pfjet(j).eta();
		}
	if (data_type == "MC_GEN")
		{
		jet_tight = true;
		pt = Event->genjet(j).pt();
		eta = Event->genjet(j).eta();
		}
	if (data_type == "MC_DET")
		{
		phi = Event->pfjet(j).phi();
		gen_pt = Event->genjet(j).pt();
		gen_eta = Event->genjet(j).eta();
		gen_phi = Event->genjet(j).phi();
		pt = smearpt(phi, eta, pt, gen_phi, gen_eta, gen_pt,test);
		}
     if (pt >= pt_min && jet_tight ) {
     counter_jet++;
     if (eta <= 2.8 && eta >= -2.8 && pt > pt_central)
     { pt_central = pt; eta_central = eta; }
     if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > pt_forward)
     { pt_forward = pt; eta_forward = eta; }
     }
     }
     
     if (pt_forward > pt_min && pt_central > pt_min)
     {
     counter_selected++;
     counter_weigth = counter_weigth + prescale*vertex_factor;
     hist_vertex_selected->Fill(Event->evtHdr().nVtxGood(),prescale*vertex_factor);
     hist_pvz_selected->Fill(Event->evtHdr().PVz(),prescale*vertex_factor);
     if (data_type == "MC_GEN" or data_type == "MC_DET") { hist_pu_selected->Fill(Event->evtHdr().pu(),prescale*vertex_factor); }
     }



      }

   }

}

}


  cout<<endl<<endl;
  cout<<"Events read:                      "<<counter_entries<<endl;
  cout<<"Events after the PV cut:          "<<counter_pv<<endl;
  if (data_type == "DATA") { cout<<"Events after the trigger cut:     "<<counter_hlt<<endl; }
  cout<<"Events Main Selection:            "<<counter_selected<<endl;
  cout<<"Events weight count               "<<counter_weigth<<endl;

     pu_scale = (double) nentries/ (double) counter_pv;
     //double cross_section = (double) counter_weigth / data_lumi[0]; }
     cout<<endl;
     cout<<"Pile up correction =  "<<pu_scale<<endl;
     //cout<<"Dataset luminosity =  "<<data_lumi[0]<<" pb^-1"<<endl;
     //cout<<"Total Cross-section = "<<cross_section<<" pb"<<endl;
     cout<<endl;

//fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_hlt);
     hist_events->SetBinContent(4,counter_selected);
     hist_events->SetBinContent(5,counter_weigth);
     hist_events->SetBinContent(6,counter_jet);

     //normalize histograms
     normalize_histogram(hist_vertex_selected, "Vertex Selected");
     normalize_histogram(hist_pvz_selected, "z-position Selected");
     normalize_histogram(hist_pu_selected, "Pileup Selected");

     //Open the output root file
     TFile *data_output= TFile::Open( data_out.c_str() , "RECREATE");

     //write histograms on disc
     hist_vertex_selected->Write();
     hist_pvz_selected->Write();
     hist_events->Write();

	if (data_type == "MC_GEN" or data_type == "MC_DET")
	{
	hist_pu_selected->Write();
	}

     //close the output file
     data_output->Close();

     //delete the histograms to avoid memory leak
     delete(hist_vertex_selected);
     delete(hist_pvz_selected);
     delete(hist_pu_selected);
     delete(hist_events);

}
