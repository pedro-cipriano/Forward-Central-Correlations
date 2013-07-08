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

void counting_combination_statistics(string *data_in, string data_out, int n_files, string sel_mode = "allvertex", bool detail = false, bool test = false)
{

   double pt_min = 35.0;
   double pt_min_gap = 20.0;
   int HLTJetPtN[3] = {15,30,50};
   int nJetTrig = 3;

//output the configuration
   if (detail) { cout<<"Get Triggered Events Configuration"<<endl; }
   if (detail) { cout<<"Selection mode :  "<<sel_mode<<endl; }
   if (detail) { cout<<"Number of files : "<<n_files<<endl; }
   if (test)   { data_out = data_out + "_test"; }
   if (detail) { cout<<"Output File :     "<<data_out<<endl; }
   if (detail) { cout<<"Pt Min :          "<<pt_min<<endl; }
   if (detail) { cout<<"Detail level :    "<<detail<<endl; }
   if (detail) { cout<<"Test Mode :       "<<test<<endl; }

//declaring the variables
   TString HLTJet[nJetTrig];

   int ihltj[nJetTrig];
   int prescalej[nJetTrig];
   int nentries = 0;
   int counter_entries = 0, counter_pv = 0, counter_hlt = 0, counter_selected = 0, counter_jet = 0;
   int triggered[6], selected[6];
   int events_selected[6][15];
   int index, index1, index2;

   char trigtitle[200];
   
   bool pv_pass = false;
   bool hltPass = false;
   bool hltPassj[nJetTrig];
   bool hard_emission;
   
   double pt, eta, phi;
   double leading_pt, leading_eta, leading_phi;
   double forward_pt, forward_eta, forward_phi;
   double central_pt, central_eta, central_phi;
   double gap_leading_pt, gap_eta, gap_phi;
   double delta_eta;
   double events_loss[2][15];
   double pu_scale, eff, trigger_loss1, trigger_loss2;

//declaring histograms
     TH1D *hist_events;
     TH1D *hist_events_HLT_Jet15U;
     TH1D *hist_events_HLT_Jet30U;
     TH1D *hist_events_HLT_Jet50U;
     TH1D *hist_events_all;
     TH1D *hist_events_inclusive_zones;
     TH1D *hist_events_exclusive_zones;

//declaring histogram binning
     hist_events =  new TH1D("Events","Selection Chain;Type;# Events", 7, 0, 7);
     hist_events_HLT_Jet15U =  new TH1D("Events_HLT_Jet15U","Selection Chain for HLT_Jet15U;Type;# Events", 2, 0, 2);
     hist_events_HLT_Jet30U =  new TH1D("Events_HLT_Jet30U","Selection Chain for HLT_Jet30U;Type;# Events", 2, 0, 2);
     hist_events_HLT_Jet50U =  new TH1D("Events_HLT_Jet50U","Selection Chain for HLT_Jet50U;Type;# Events", 2, 0, 2);
     hist_events_all =  new TH1D("Events_Combined","Selection Chain for combination;Type;# Events", 2, 0, 2);
     hist_events_inclusive_zones =  new TH1D("Events_Inclusive_Zones","Selection Chain for combination Inclusive Zones;Type;# Events", 3, 0, 3);
     hist_events_exclusive_zones =  new TH1D("Events_Exclusive_Zones","Selection Chain for combination Exclusive Zones;Type;# Events", 3, 0, 3);
 
//inicializing the histogram statistics
     hist_events->Sumw2();
     hist_events->Fill("Total Events",0);
     hist_events->Fill("PV Selection",0);
     hist_events->Fill("Trigger Selection",0);
     hist_events->Fill("Selected",0);
     hist_events->Fill("Number of Jets",0);
     hist_events->Fill("Pileup Scale",0);
     hist_events->Fill("Efficiency",0);

     hist_events_HLT_Jet15U->Sumw2();
     hist_events_HLT_Jet15U->Fill("Triggered Events",0);
     hist_events_HLT_Jet15U->Fill("Selected",0);

     hist_events_HLT_Jet30U->Sumw2();
     hist_events_HLT_Jet30U->Fill("Triggered Events",0);
     hist_events_HLT_Jet30U->Fill("Selected",0);

     hist_events_HLT_Jet50U->Sumw2();
     hist_events_HLT_Jet50U->Fill("Triggered Events",0);
     hist_events_HLT_Jet50U->Fill("Selected",0);

     hist_events_all->Sumw2();
     hist_events_all->Fill("Triggered Events",0);
     hist_events_all->Fill("Selected",0);

     hist_events_inclusive_zones->Sumw2();
     hist_events_inclusive_zones->Fill("Triggered Events",0);
     hist_events_inclusive_zones->Fill("Selected",0);
     hist_events_inclusive_zones->Fill("Proportion",0);

     hist_events_exclusive_zones->Sumw2();
     hist_events_exclusive_zones->Fill("Triggered Events",0);
     hist_events_exclusive_zones->Fill("Selected",0);
     hist_events_exclusive_zones->Fill("Proportion",0);

for (int m = 0; m < 15; m++)
	{
	for (int y = 0; y < 6; y++)
		{
		events_selected[y][m] = 0;
		}
	events_loss[0][m] = 0.0;
	events_loss[1][m] = 0.0;
	}

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
      for (int i = 0; i<=nJetTrig+2; i++)
         {
                triggered[i] = 0;
		selected[i] = 0;
         }

      for (int i = 0; i<nJetTrig; i++)
         {
                sprintf(trigtitle, "HLT_Jet%iU",HLTJetPtN[i]);             
		HLTJet[i] = trigtitle;
		//cout << "Testing " << trigtitle << endl;
         }

      for (int i = 0; i<nJetTrig; i++)
         {
                  ihltj[i] = -1 ;
         }

   TH1F *hTrigNames = (TH1F*)inf->Get("ak5/TriggerNames");
   // ------------ Assigning An Integer To Each Trigger Path---------------//
    for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++)
	{
       	TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));

       	for (int ii=0; ii<nJetTrig; ii++)
          {
          if (ss == HLTJet[ii])
		{
                ihltj[ii] = ibin;
                continue;
          	}
          } // for (int ii=0; ii<nJetTrig; ii++)
        } // for(int ibin=0;ibin<55;ibin++)

  // ----------------------- Checking For The Trigger Assignment ------------------- // 
    for (int ij=0; ij<nJetTrig; ij++)
        {
              if (ihltj[ij] == -1)
		{
              	cout<<"The requested trigger ("<<HLTJet[ij]<<") is not found "<<endl;
   //         	break;
               	}
              else
		{
              	cout<<HLTJet[ij]<<" --> "<<ihltj[ij]<<endl;
             	}
      	}
      
  nentries = tr->GetEntries();
  cout<<"Reading TREE: "<<nentries<<" events"<<endl;
  int decade = 0;

  if (test) { nentries = 100000; } //reduced number of read entries, usefull for testing

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
        for (int j=0; j<=nJetTrig+2; j++)
        {
                   hltPassj[j] = false ;
        } // for (int i=0; i<nJetTrig; i++)

  //----------------------- Computing The Prescale Values For a given event ----------------------- //
    for (int j=0; j<nJetTrig; j++)
         {
                 if (ihltj[j] == -1)
                 hltPassj[j] = true; // no trigger set
                 else
                 	{
                 	if (Event->fired(ihltj[j]) > 0)
				{
                 		hltPassj[j] = true;
		 		triggered[j] = triggered[j] + 1;
                 		hltPass = true;
                 		prescalej[j] = Event->preL1(ihltj[j]) * Event->preHLT(ihltj[j]);
        //         		if(i==0)
        //         		cout<<HLTJet[i][j]<<" has prescale = "<<prescalej[i][j]<<endl;
           			}
      		}
      	}
      	

     if (hltPass)
	{
	counter_hlt++;
	triggered[3] = triggered[3] + 1;

     	leading_pt = -10.0;
     	forward_pt = -10.0;
     	central_pt = -10.0;
     	gap_leading_pt = -10.0;
     	forward_eta = -10.0;
     	central_eta = -10.0;
     	leading_eta = -10.0;
     	gap_eta = -10.0;
     	forward_phi = -10.0;
     	central_phi = -10.0;
     	leading_phi = -10.0;
     	gap_phi = -10.0;
	hard_emission = false;
	delta_eta = 0.0;
	index = -1;
	index1 = -1;
	index2 = -1;

    for(unsigned int l=0; l<Event->nPFJets(); l++)
	{
	pt = Event->pfjet(l).ptCor();
	eta = Event->pfjet(l).eta();
	phi = Event->pfjet(l).phi();
     	if (pt >= pt_min && Event->pfjet(l).tightID() )
		{
     		counter_jet++;
     		if (leading_pt < pt) { leading_pt = pt; leading_eta = eta; leading_phi = phi; }
     		if (eta <= 2.8 && eta >= -2.8 && pt > central_pt) { central_pt = pt; central_eta = eta; central_phi = phi; }
     		if (((eta >= 3.2 && eta <= 4.7) || (eta <= -3.2 && eta >= -4.7)) && pt > forward_pt)
     			{
			forward_pt = pt; forward_eta = eta; forward_phi = phi;
			}
     		}
     	}

	if (leading_pt > 35.0 and leading_pt < 70.0 and hltPassj[0])
		{
		hltPassj[4] = true;
		triggered[4] = triggered[4] + 1;
		}
	if (leading_pt > 70.0 and leading_pt < 110.0 and hltPassj[1])
		{
		hltPassj[4] = true;
		triggered[4] = triggered[4] + 1;
		}

	if (leading_pt > 110.0 and hltPassj[2])
		{
		hltPassj[4] = true;
		triggered[4] = triggered[4] + 1;
		}

	if (leading_pt > 35.0 and leading_pt < 70.0 and hltPassj[0] and !hltPassj[1] and !hltPassj[2])
		{
		hltPassj[5] = true;
		triggered[5] = triggered[5] + 1;
		}
	if (leading_pt > 70.0 and leading_pt < 110.0 and !hltPassj[0] and hltPassj[1] and !hltPassj[2])
		{
		hltPassj[5] = true;
		triggered[5] = triggered[5] + 1;
		}

	if (leading_pt > 110.0 and !hltPassj[0] and !hltPassj[1] and hltPassj[2])
		{
		hltPassj[5] = true;
		triggered[5] = triggered[5] + 1;
		}
     
          if (forward_pt > pt_min && central_pt > pt_min)
     		{
     		counter_selected++;
		delta_eta = abs(forward_eta - central_eta);
    		for(unsigned int n=0; n<Event->nPFJets(); n++)
			{
			pt = Event->pfjet(n).ptCor();
			eta = Event->pfjet(n).eta();
			phi = Event->pfjet(n).phi();
     			if (Event->pfjet(n).tightID() )
				{
     				if (central_eta > forward_eta && pt > gap_leading_pt && eta > forward_eta && eta < central_eta )
					{
        				gap_leading_pt = pt; gap_eta = eta; gap_phi = phi;
					}
     				if (central_eta < forward_eta && pt > gap_leading_pt && eta < forward_eta && eta > central_eta )
					{
       					gap_leading_pt = pt; gap_eta = eta; gap_phi = phi;
					}
        			}
     			}

			if (gap_leading_pt > pt_min_gap) { index1 = 1; }
			else { index1 = 2; }
			if (delta_eta >= 0.4 and delta_eta < 2.5) { index2 = 0; }
			if (delta_eta >= 2.5 and delta_eta < 3.5) { index2 = 1; }
			if (delta_eta >= 3.5 and delta_eta < 4.5) { index2 = 2; }
			if (delta_eta >= 4.4 and delta_eta < 7.5) { index2 = 3; }

			index = index1 * 5 + index2;

    		for (int l=0; l<nJetTrig; l++)
         		{
			if (hltPassj[l])
				{
				selected[l] = selected[l] + 1;
				events_selected[l][index] = events_selected[l][index] + 1;
				events_selected[l][index2] = events_selected[l][index2] + 1;
				events_selected[l][index1 * 5 + 4] = events_selected[l][index1 * 5 + 4] + 1;
				events_selected[l][4] = events_selected[l][4] + 1;
				}
			}
		selected[3] = selected[3] + 1;
		events_selected[3][index] = events_selected[3][index] + 1;
		events_selected[3][index2] = events_selected[3][index2] + 1;
		events_selected[3][index1 * 5 + 4] = events_selected[3][index1 * 5 + 4] + 1;
		events_selected[3][4] = events_selected[3][4] + 1;
		if (hltPassj[4])
			{
			selected[4] = selected[4] + 1;
			events_selected[4][index] = events_selected[4][index] + 1;
			events_selected[4][index2] = events_selected[4][index2] + 1;
			events_selected[4][index1 * 5 + 4] = events_selected[4][index1 * 5 + 4] + 1;
			events_selected[4][4] = events_selected[4][4] + 1;
			}
		if (hltPassj[5])
			{
			selected[5] = selected[5] + 1;
			events_selected[5][index] = events_selected[5][index] + 1;
			events_selected[5][index2] = events_selected[5][index2] + 1;
			events_selected[5][index1 * 5 + 4] = events_selected[5][index1 * 5 + 4] + 1;
			events_selected[5][4] = events_selected[5][4] + 1;
			}
		}
	}


}

} 
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
     if (detail) { cout<<"Total Jets :                     "<<counter_jet<<endl; }
     
     //fill the events histogram
     hist_events->SetBinContent(1,counter_entries);
     hist_events->SetBinContent(2,counter_pv);
     hist_events->SetBinContent(3,counter_hlt);
     hist_events->SetBinContent(4,counter_selected);
     hist_events->SetBinContent(7,counter_jet);
     hist_events->SetBinContent(8,pu_scale);
     hist_events->SetBinContent(9,eff);

     //output a summary for HLT_Jet15U
     if (detail) { cout<<"Results for HLT_Jet15U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet15U: "<<triggered[0]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[0]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[0][0]<< "        "<<setw(7)<<events_selected[0][5]<<"       "<<setw(7)<<events_selected[0][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[0][1]<< "        "<<setw(7)<<events_selected[0][6]<<"       "<<setw(7)<<events_selected[0][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[0][2]<< "        "<<setw(7)<<events_selected[0][7]<<"       "<<setw(7)<<events_selected[0][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[0][3]<< "        "<<setw(7)<<events_selected[0][8]<<"       "<<setw(7)<<events_selected[0][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[0][4]<< "        "<<setw(7)<<events_selected[0][9]<<"       "<<setw(7)<<events_selected[0][14]<<endl; }

     //fill the events in HLT_Jet15U histogram
     hist_events_HLT_Jet15U->SetBinContent(1,triggered[0]);
     hist_events_HLT_Jet15U->SetBinContent(2,selected[0]);

     //output a summary for HLT_Jet30U
     if (detail) { cout<<"Results for HLT_Jet30U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet30U:  "<<triggered[1]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[1]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[1][0]<< "        "<<setw(7)<<events_selected[1][5]<<"       "<<setw(7)<<events_selected[1][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[1][1]<< "        "<<setw(7)<<events_selected[1][6]<<"       "<<setw(7)<<events_selected[1][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[1][2]<< "        "<<setw(7)<<events_selected[1][7]<<"       "<<setw(7)<<events_selected[1][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[1][3]<< "        "<<setw(7)<<events_selected[1][8]<<"       "<<setw(7)<<events_selected[1][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[1][4]<< "        "<<setw(7)<<events_selected[1][9]<<"       "<<setw(7)<<events_selected[1][14]<<endl; }

     //fill the events in HLT_Jet30U histogram
     hist_events_HLT_Jet30U->SetBinContent(1,triggered[1]);
     hist_events_HLT_Jet30U->SetBinContent(2,selected[1]);

     //output a summary for HLT_Jet50U
     if (detail) { cout<<"Results for HLT_Jet50U" << endl; }
     if (detail) { cout<<"Total Triggered Events by HLT_Jet50U:  "<<triggered[2]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[2]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[2][0]<< "        "<<setw(7)<<events_selected[2][5]<<"       "<<setw(7)<<events_selected[2][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[2][1]<< "        "<<setw(7)<<events_selected[2][6]<<"       "<<setw(7)<<events_selected[2][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[2][2]<< "        "<<setw(7)<<events_selected[2][7]<<"       "<<setw(7)<<events_selected[2][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[2][3]<< "        "<<setw(7)<<events_selected[2][8]<<"       "<<setw(7)<<events_selected[2][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[2][4]<< "        "<<setw(7)<<events_selected[2][9]<<"       "<<setw(7)<<events_selected[2][14]<<endl; }

     //fill the events in HLT_Jet50U histogram
     hist_events_HLT_Jet50U->SetBinContent(1,triggered[2]);
     hist_events_HLT_Jet50U->SetBinContent(2,selected[2]);

     //output a summary for combination
     if (detail) { cout<<"Results for combination" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[3]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[3]<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[3][0]<< "        "<<setw(7)<<events_selected[3][5]<<"       "<<setw(7)<<events_selected[3][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[3][1]<< "        "<<setw(7)<<events_selected[3][6]<<"       "<<setw(7)<<events_selected[3][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[3][2]<< "        "<<setw(7)<<events_selected[3][7]<<"       "<<setw(7)<<events_selected[3][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[3][3]<< "        "<<setw(7)<<events_selected[3][8]<<"       "<<setw(7)<<events_selected[3][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[3][4]<< "        "<<setw(7)<<events_selected[3][9]<<"       "<<setw(7)<<events_selected[3][14]<<endl; }

     //fill the events in combination histogram
     hist_events_all->SetBinContent(1,triggered[3]);
     hist_events_all->SetBinContent(2,selected[3]);

     trigger_loss1 = (double) selected[4] / (double) selected[3];
     trigger_loss2 = (double) selected[5] / (double) selected[3];

for (int m = 0; m < 15; m++)
	{
		events_loss[0][m] = (double) events_selected[4][m] / (double) events_selected[3][m];
		events_loss[1][m] = (double) events_selected[5][m] / (double) events_selected[3][m];
	}

     //output a summary for combination
     if (detail) { cout<<"Results for combination with inclusive zones" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[4]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[4]<<endl; }
     if (detail) { cout<<"Proportion :                           "<<trigger_loss1<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[4][0]<< "        "<<setw(7)<<events_selected[4][5]<<"       "<<setw(7)<<events_selected[4][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[4][1]<< "        "<<setw(7)<<events_selected[4][6]<<"       "<<setw(7)<<events_selected[4][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[4][2]<< "        "<<setw(7)<<events_selected[4][7]<<"       "<<setw(7)<<events_selected[4][11]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[4][3]<< "        "<<setw(7)<<events_selected[4][8]<<"       "<<setw(7)<<events_selected[4][12]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[4][4]<< "        "<<setw(7)<<events_selected[4][9]<<"       "<<setw(7)<<events_selected[4][14]<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_loss[0][0]<< "        "<<setw(7)<<events_loss[0][5]<<"       "<<setw(7)<<events_loss[0][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_loss[0][1]<< "        "<<setw(7)<<events_loss[0][6]<<"       "<<setw(7)<<events_loss[0][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_loss[0][2]<< "        "<<setw(7)<<events_loss[0][7]<<"       "<<setw(7)<<events_loss[0][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_loss[0][3]<< "        "<<setw(7)<<events_loss[0][8]<<"       "<<setw(7)<<events_loss[0][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_loss[0][4]<< "        "<<setw(7)<<events_loss[0][9]<<"       "<<setw(7)<<events_loss[0][14]<<endl; }

     //fill the events in combination histogram
     hist_events_inclusive_zones->SetBinContent(1,triggered[4]);
     hist_events_inclusive_zones->SetBinContent(2,selected[4]);
     hist_events_inclusive_zones->SetBinContent(3,trigger_loss1);

     //output a summary for combination
     if (detail) { cout<<"Results for combination with exclusive zones" << endl; }
     if (detail) { cout<<"Total Triggered Events :               "<<triggered[5]<<endl; }
     if (detail) { cout<<"Selected :                             "<<selected[5]<<endl; }
     if (detail) { cout<<"Proportion :                           "<<trigger_loss2<<endl; }
     if (detail) { cout<<"Delta eta  Main Scenario  Hard Emmision  Veto Emission"<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_selected[5][0]<< "        "<<setw(7)<<events_selected[5][5]<<"       "<<setw(7)<<events_selected[5][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_selected[5][1]<< "        "<<setw(7)<<events_selected[5][6]<<"       "<<setw(7)<<events_selected[5][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_selected[5][2]<< "        "<<setw(7)<<events_selected[5][7]<<"       "<<setw(7)<<events_selected[5][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_selected[5][3]<< "        "<<setw(7)<<events_selected[5][8]<<"       "<<setw(7)<<events_selected[5][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_selected[5][4]<< "        "<<setw(7)<<events_selected[5][9]<<"       "<<setw(7)<<events_selected[5][14]<<endl; }
     if (detail) { cout<<"0.4 - 2.5  "<<setw(7)<<events_loss[1][0]<< "        "<<setw(7)<<events_loss[1][5]<<"       "<<setw(7)<<events_loss[1][10]<<endl; }
     if (detail) { cout<<"2.5 - 3.5  "<<setw(7)<<events_loss[1][1]<< "        "<<setw(7)<<events_loss[1][6]<<"       "<<setw(7)<<events_loss[1][11]<<endl; }
     if (detail) { cout<<"3.5 - 4.5  "<<setw(7)<<events_loss[1][2]<< "        "<<setw(7)<<events_loss[1][7]<<"       "<<setw(7)<<events_loss[1][12]<<endl; }
     if (detail) { cout<<"4.5 - 7.5  "<<setw(7)<<events_loss[1][3]<< "        "<<setw(7)<<events_loss[1][8]<<"       "<<setw(7)<<events_loss[1][13]<<endl; }
     if (detail) { cout<<"All        "<<setw(7)<<events_loss[1][4]<< "        "<<setw(7)<<events_loss[1][9]<<"       "<<setw(7)<<events_loss[1][14]<<endl; }


     //fill the events in combination histogram
     hist_events_exclusive_zones->SetBinContent(1,triggered[5]);
     hist_events_exclusive_zones->SetBinContent(2,selected[5]);
     hist_events_exclusive_zones->SetBinContent(3,trigger_loss2);

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
     hist_events_inclusive_zones->Write();
     hist_events_exclusive_zones->Write();

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
     delete(hist_events_inclusive_zones);
     delete(hist_events_exclusive_zones);

//Sucess confirmation
     if (detail) { cout<<"Done!"<<endl; }
}
