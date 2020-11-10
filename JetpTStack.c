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


void JetpTStack()
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



TH1D* h2 = (TH1D*)DY->Get("DYpTJetsS"); 
TH1D* h5 = (TH1D*)ttbar->Get("ttbarpTJetsS");
TH1D* h3 = (TH1D*)QCD->Get("QCDpTJetsS");
TH1D* h9 = (TH1D*)ZZ->Get("ZZpTJetsS");
TH1D* h8 = (TH1D*)ZW->Get("ZWpTJetsS");
TH1D* h7 = (TH1D*)WW->Get("WWpTJetsS");
TH1D* h4 = (TH1D*)SingleTop->Get("SingleToppTJetsS");

THStack *pTJetsS = new THStack("nJets","Number of Jets");
TCanvas *c1= new TCanvas("pT","Transverse Momenta of Jets",800,600);

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




	pTJetsS->Add(h3);
	pTJetsS->Add(h9);
	pTJetsS->Add(h8);
	pTJetsS->Add(h7);
	pTJetsS->Add(h2);
	pTJetsS->Add(h4);
	pTJetsS->Add(h5);

	c1->cd();
	pTJetsS->Draw();
	

        c1->Update();
	c1->SaveAs("bJet pT.pdf");
}


int main()
{
JetpTStack();
return 0;
}
