#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TAxis.h"


void nLeptonStack()
{ 

gROOT->SetStyle("Plain");
gStyle->SetOptTitle(0);
gStyle->SetOptStat("nemrou");

TFile *DY = new TFile("DY.root");
TFile *ttbar = new TFile("ttbar.root");
TFile *QCD = new TFile("QCD.root");
TFile *SingleTop = new TFile("SingleTop.root");
TFile *WW = new TFile("WW.root");
TFile *ZW = new TFile("ZW.root");
TFile *ZZ = new TFile("ZZ.root");



TH1D* h2 = (TH1D*)DY->Get("DYnLepton"); 
TH1D* h5 = (TH1D*)ttbar->Get("ttbarnLepton");
TH1D* h3 = (TH1D*)QCD->Get("QCDnLepton");
TH1D* h9 = (TH1D*)ZZ->Get("ZZnLepton");
TH1D* h8 = (TH1D*)ZW->Get("ZWnLepton");
TH1D* h7 = (TH1D*)WW->Get("WWnLepton");
TH1D* h4 = (TH1D*)SingleTop->Get("SingleTopnLepton");

THStack *nLepton = new THStack("nJets","Number of Jets");

	h2->GetXaxis()->SetRange(0.0,20.0);
	h2->SetLineColor(2);
	h2->SetFillColor(2);
	h2->SetMarkerStyle(21);
	h2->SetMarkerColor(2);

	h3->SetLineColor(3);
	h3->SetFillColor(3);
	h3->SetMarkerStyle(21);
	h3->SetMarkerColor(3);
	
	h4->SetLineColor(4);
	h4->SetFillColor(4);
	h4->SetMarkerStyle(21);
	h4->SetMarkerColor(4);

	h5->SetLineColor(5);
	h5->SetFillColor(5);
	h5->SetMarkerStyle(21);
	h5->SetMarkerColor(5);

	h7->SetLineColor(7);
	h7->SetFillColor(7);
	h7->SetMarkerStyle(21);
	h7->SetMarkerColor(7);

	h8->SetLineColor(8);
	h8->SetFillColor(8);
	h8->SetMarkerStyle(21);
	h8->SetMarkerColor(8);

	h9->SetLineColor(9);
	h9->SetFillColor(9);
	h9->SetMarkerStyle(21);
	h9->SetMarkerColor(9);




	nLepton->Add(h3);
	nLepton->Add(h9);
	nLepton->Add(h8);
	nLepton->Add(h7);
	nLepton->Add(h2);
	nLepton->Add(h4);
	nLepton->Add(h5);

	//nLepton->cd();
	nLepton->Draw();
	

      //nLepton->Update();
	nLepton->SaveAs("bJet pT.pdf");	
}


int main()
{
nLeptonStack();
return 0;
}
