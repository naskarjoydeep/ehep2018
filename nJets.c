#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TMacro.h"
//#include <iostream.h>
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
//#include "Latex.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TArrow.h"
	



void nJets()
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



TH1D* h2 = (TH1D*)DY->Get("DYnJetsS"); 
TH1D* h5 = (TH1D*)ttbar->Get("ttbarnJetsS");
TH1D* h3 = (TH1D*)QCD->Get("QCDnJetsS");
TH1D* h9 = (TH1D*)ZZ->Get("ZZnJetsS");
TH1D* h8 = (TH1D*)ZW->Get("ZWnJetsS");
TH1D* h7 = (TH1D*)WW->Get("WWnJetsS");
TH1D* h4 = (TH1D*)SingleTop->Get("SingleTopnJetsS");

TCanvas *nJetsS = new TCanvas("nJets","Number of Jets");
//pT->SetLogy();

//h2-> Scale(5.3);
//h3-> Scale(740.);
//h1-> Scale(800000);

//h2-> Scale(0.01);
//h3-> Scale(0.01);
//h1-> Scale(0.01);





h3->SetLineColor(3);
h2->SetLineColor(2);
h9->SetLineColor(9);
h4->SetLineColor(4);
h5->SetLineColor(5);
h8->SetLineColor(8);
h7->SetLineColor(7);


//h1->GetXaxis()->SetRange(0.0,100.0);



h2->Draw();
h3->Draw("same");
h4->Draw("same");
h5->Draw("same");
h7->Draw("same");
h8->Draw("same");
h9->Draw("same");


nJetsS->BuildLegend();



nJetsS -> Print("nJetsOverlay.pdf");
}

int main()
{
nJets();
return 0;
}
