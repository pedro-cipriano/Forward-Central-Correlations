// Created: Knutsson Aug 2010
// Integrated : Pedro Cipriano Mar 2013 

#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>

#include <iostream>

using namespace std;

TH1F* Gr2Hi(const TGraph* GraphIn){

//init
double xerr[100];
double yerr[100];
double xbinlimits[101];	   
for (Int_t j=0;j<100;j++){
  xerr[j] = 0.0;
  yerr[j] = 0.0;
  xbinlimits[j] = 0.0;
}
xbinlimits[100] = 0.0;


//get points and info
int Npnts= GraphIn->GetN();
cout << "Npoints = " << Npnts << endl;
double *x = GraphIn->GetX();
double *y = GraphIn->GetY();
const char *grName = GraphIn->GetName();


//get errors and create bin limits
for (int j=0;j<Npnts;j++){
  xerr[j] = GraphIn->GetErrorX(j);
  yerr[j] = GraphIn->GetErrorY(j);
  xbinlimits[j]=x[j]-xerr[j];
  //cout << j << " - " << xbinlimits[j] << endl;
}
xbinlimits[Npnts]=x[Npnts-1]+xerr[Npnts-1];  


//Create and fill Histogram
Char_t  hiName[100];
sprintf(hiName,"%s_histo",grName);
TH1F *HistoOut;
HistoOut = new TH1F(hiName, hiName,Npnts, xbinlimits);
for (Int_t j=0;j<Npnts;j++){
  HistoOut->Fill(x[j],y[j]);
  HistoOut->SetBinError(j+1,yerr[j]);
}

return HistoOut;
}

void convert_prediction(string path_in, string path_out, bool detail = false)
{
   if (detail) { cout<<"Converting Prediction Configuration"<<endl; }
   if (detail) { cout<<"Input file:   "<<path_in<<endl; }
   if (detail) { cout<<"Output file : "<<path_out<<endl; }
   if (detail) { cout<<"Detail :      "<<detail<<endl; }

//opening file
if (detail) { cout << "Opening files..." << endl; }
TFile *file_in = 0;
file_in = new TFile( path_in.c_str() );
if (file_in == 0) { cout << "Source file " << path_in << "not found!" << endl; return; }

//converting histograms
    if (detail) { cout << "Converting Delta Phi..." << endl; }
    TH1F *delta_phi;
    TGraph *delta_phi_in = 0;
    delta_phi_in = (TGraph*) file_in->Get("d01-x01-y01");
    if (delta_phi_in == 0) { cout<< "TGraph for delta_phi not found" << endl; return; }
    delta_phi = Gr2Hi(delta_phi_in);
    delta_phi->SetName("ak5Gen_delta_phi");

    if (detail) { cout << "Converting Delta Phi Gap..." << endl; }
    TH1F *delta_phi_gap;
    TGraph *delta_phi_gap_in = 0;
    delta_phi_gap_in = (TGraph*) file_in->Get("d02-x01-y01");
    if (delta_phi_gap_in == 0) { cout<< "TGraph for delta_phi_gap not found" << endl; return; }
    delta_phi_gap = Gr2Hi(delta_phi_gap_in); 
    delta_phi_gap->SetName("ak5Gen_delta_phi_gap");

    if (detail) { cout << "Converting Delta Phi Nogap..." << endl; }
    TH1F *delta_phi_nogap;
    TGraph *delta_phi_nogap_in = 0;
    delta_phi_nogap_in = (TGraph*) file_in->Get("d03-x01-y01");
    if (delta_phi_nogap_in == 0) { cout<< "TGraph for delta_phi_nogap not found" << endl; return; }
    delta_phi_nogap = Gr2Hi(delta_phi_nogap_in);  
    delta_phi_nogap->SetName("ak5Gen_delta_phi_nogap");

    if (detail) { cout << "Creating output file " << path_out << endl; }
    TFile *file_out = new TFile( path_out.c_str() , "RECREATE"); 

    if (detail) { cout << "Writing histograms on file..." << endl; }
    delta_phi->Write();
    delta_phi_gap->Write();
    delta_phi_nogap->Write();

    if (detail) { cout << "Closing files..." << endl; }
    file_in->Close();
    file_out->Close();

    if (detail) { cout << "Deleting variables..." << endl; }
    //delete(delta_phi);
    if (detail) { cout << "Deleting variables..." << endl; }
    //delete(delta_phi_gap);
    //delete(delta_phi_nogap);

    //delete(delta_phi_in);
    //delete(delta_phi_gap_in);
    //delete(delta_phi_nogap_in);
    if (detail) { cout << "Done!" << endl; }
}
