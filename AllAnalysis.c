/*#include "MyAnalysis.h"
#include "Plotter.h"*/
#include <iostream>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include "TH1F.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TAxis.h"

float MET(float MET_px, float MET_py)
	{	float met=sqrt(pow(MET_px,2)+pow(MET_py,2));
		return met;
	}	

float PT(float px, float py)
	{	float pT=sqrt(pow(px,2)+pow(py,2));
		return pT;
	}

float P(float px, float py, float pz)
	{	float p=sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
		return p;
	}
float ETA(float px, float py, float pz)
	{ 	float p=P(px,py,pz);
		float eta=atanh(pz/p);
		return eta;
	}

float PHI(float px, float py)
	{
		float phi=0;
		phi=atan(py/px);
		return phi;			
	}


float HighestpT(int N, float px[N], float py[N], float pz[N])
	{	float pTmax=0;	
		for (int i=0; i<N; i++)
		{	float pt=PT(px[i], py[i]);
			float eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
			if (eta<2.5)	if(pt>pTmax)
				pTmax=pt;
		}
		return pTmax;
	}



float etaHighestpT(int N, float px[N], float py[N], float pz[N])
	{	float pTmax=0;	float eta0=0;
		for (int i=0; i<N; i++)
		{	float pt=PT(px[i], py[i]);
			float eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
			if (eta<2.5)	if(pt>pTmax)
				{pTmax=pt; eta0=eta;}	
		}
		return eta0;
	}

float HT(int N, float px[N], float py[N], float pz[N])
	{	float ht=0;		
		for(int i=0; i<N; i++)
		{
			float p=0;	float eta;
			p=PT(px[i], py[i]);
			eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
			if(p>25 && eta<2.5)
			ht=ht+p;
		}
		return ht;
	}

float ST(int NJet, float pxj[NJet], float pyj[NJet], float pzj[NJet], int NMuon, float pxm[NMuon], float pym[NMuon], float pzm[NMuon])
	{
		float jht, lht, st;	jht=0;	lht=0;	st=0;
		for (int i=0; i<NJet; i++)
		{
			float pt=0;	float eta;
			pt=PT(pxj[i], pyj[i]);
			eta=ETA(pxj[i], pyj[i], pzj[i]);	eta=abs(eta);
			if(pt>25 && eta<2.5)
			jht=jht+pt;
		}
		for (int i=0; i<NMuon; i++)
		{
			float pt=0;	float eta;
			pt=PT(pxm[i], pym[i]);
			eta=ETA(pxm[i], pym[i], pzm[i]);	eta=abs(eta);
			if(pt>25 && eta<2.5)
			lht=lht+pt;
		}
		st=lht+jht;
		return st;
	}



int Factorial(int n)	{	int fact=1;	for (int i=1; i<=n; i++)	fact=fact*i;	return fact;	}

float GetMin(int N, int value[N])
	{	float min=value[0];
		for (int i=1; i<N; i++)
			if (value[i]<min)	min=value[i];
		return min;
	}

float MTab(float Ea, float pxa, float pya, float pza, float Eb, float pxb, float pyb, float pzb)
	{	float etaa=ETA(pxa,pya,pza);
		float etab=ETA(pxb,pyb,pzb);
		float ETa=Ea/cosh(etaa);
		float ETb=Eb/cosh(etab);
		float ET=ETa+ETb;
		float px=pxa+pxb;	float py=pya+pyb;	float pz=pza+pzb;
		float MT=sqrt(abs(pow(ET,2)-pow(px,2)-pow(py,2)));
		return MT;
	}

float M0ab(float Ea, float pxa, float pya, float pza, float Eb, float pxb, float pyb, float pzb)
	{	float E=Ea+Eb;
		float px=pxa+pxb;	float py=pya+pyb;	float pz=pza+pzb;
		float M0=sqrt(abs(pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2)));
		return M0;
	}




float MTcd(float pxa, float pya, float pxb, float pyb)	
	{
		float phil,phin, ptl, ptn;
		phil=PHI(pxa, pya);	
		phin=PHI(pxb, pyb);	
		ptl=PT(pxa, pya);
		ptn=PT(pxb, pyb);
		float mt=0;
		mt=sqrt(2*ptl*ptn*(1-cos(phil-phin)));
		return mt;
	}



float BestMTab(int N, bool JetID[N], float Jet_btag[N], float E[N], float px[N], float py[N], float pz[N])
	{
		int length=N*(N-1)/2;
		float eta, pt;
		float mt[N][N];	for(int i=0; i<N; i++)	for(int j=0; j<N; j++) {mt[i][j]=0;}
		for(int i=0; i<N; i++){
					eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
					pt=PT(px[i], py[i]);					
						if (JetID[i] && Jet_btag[i]<1.76 && eta<2.5 && pt>25)
			for (int j=i+1; j<N; j++){	
					eta=ETA(px[j], py[j], pz[j]);
					pt=PT(px[j], py[j]);	
						if (JetID[j] && Jet_btag[j]<1.76 && eta<2.5 && pt>25)	
						mt[i][j]=MTab(E[i], px[i], py[i], pz[i], E[j], px[j], py[j], pz[j]);
						}	
					}
		float leastdiff=80.385;	//W-Mass
		int pointeri=0; int pointerj=0;		
		for(int i=0; i<N; i++)	
			for(int j=0; j<N; j++) 
			{	float diff;
				diff=mt[i][j]-80.385;	diff=abs(diff);
				if(diff<leastdiff)	{leastdiff=diff;	pointeri=i;	pointerj=j;}
			}
		float bestmt=mt[pointeri][pointerj];
		return bestmt;
	}

float BestM0ab(int N, bool JetID[N], float Jet_btag[N], float E[N], float px[N], float py[N], float pz[N])
	{
		int length=N*(N-1)/2;
		float eta, pt;
		float m0[N][N];	for(int i=0; i<N; i++)	for(int j=0; j<N; j++) {m0[i][j]=0;}
		for(int i=0; i<N; i++){
					eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
					pt=PT(px[i], py[i]);					
						if (JetID[i] && Jet_btag[i]<1.76 && eta<2.5 && pt>25)
			for (int j=i+1; j<N; j++){	
					eta=ETA(px[j], py[j], pz[j]);
					pt=PT(px[j], py[j]);	
						if (JetID[j]&& Jet_btag[i]<1.76 && eta<2.5 && pt>25)	
						m0[i][j]=M0ab(E[i], px[i], py[i], pz[i], E[j], px[j], py[j], pz[j]);
						}	
					}
		float leastdiff=80.385;	//W-Mass
		int pointeri=0; int pointerj=0;		
		for(int i=0; i<N; i++)	
			for(int j=0; j<N; j++) 
			{	float diff;
				diff=m0[i][j]-80.385;	diff=abs(diff);
				if(diff<leastdiff)	{leastdiff=diff;	pointeri=i;	pointerj=j;}
			}
		float bestm0=m0[pointeri][pointerj];
		return bestm0;
	}


float BestMTcd(int N, float Muon_Iso[N], float pxl[N], float pyl[N], float pzl[N], float MEx, float MEy)
	{
		float eta, pt, relIso;
		float mt[N];	for(int i=0; i<N; i++)	{mt[i]=0;}
		for(int i=0; i<N; i++)	{
					eta=ETA(pxl[i], pyl[i], pzl[i]);	eta=abs(eta);
					pt=PT(pxl[i], pyl[i]);	
					relIso=Muon_Iso[i]/pt;				
						if (relIso<0.05 && eta<2.5 && pt>25)
						mt[i]=MTcd(pxl[i], pyl[i], MEx, MEy);							
					}
		float leastdiff=80.385;	//W-Mass
		int pointeri=0;		
		for(int i=0; i<N; i++)
			{	float diff;
				diff=mt[i]-80.385;	diff=abs(diff);
				if(diff<leastdiff)	{leastdiff=diff;	pointeri=i;}
			}
		float bestmt=mt[pointeri];
		return bestmt;
	}


int JetSelect(bool JetID, float px, float py, float pz)
	{
		float pt=0; pt=PT(px,py);
		float eta=ETA(px,py,pz);	eta=abs(eta);
		if (JetID && pt>25 && eta<2.5)	
			return 1;
		else
			return 0;
	}

int bJetSelect(float JetID, float btag ,float px, float py, float pz)
	{	
		int JetSel=JetSelect(JetID, px, py, pz);
		if (JetSel && btag>=1.76)
			 return 1;
		else
			return 0;
	}

int JetCounter(int NJet, bool Jet_ID[NJet], float px[NJet], float py[NJet], float pz[NJet])
	{
		int counter=0;	int JetSel=0;	int i=0;	float pt=0;
		for (int i=0; i<NJet; i++)		
		{	pt=PT(px[i],py[i]);
			JetSel=JetSelect(Jet_ID[i], px[i], py[i], pz[i]); 
			if(JetSel)	counter++;
		}
		return counter;
	}
				
int bJetCounter(int NJet, bool Jet_ID[NJet], float Jet_btag[NJet], float px[NJet], float py[NJet], float pz[NJet])
	{	int counter=0;	int bJetSel=0;
		for(int i=0; i<NJet; i++)	
		{	bJetSel=bJetSelect(Jet_ID[i], Jet_btag[i], px[i], py[i], pz[i]); 
			if(bJetSel)	counter++;
		}
		return counter;
	}

int MuonSelect(float Muon_Iso, float px, float py, float pz)
	{	float relIso;
		float pt=PT(px,py);
		relIso=Muon_Iso/pt;
		float eta=ETA(px,py,pz);	eta=abs(eta);
		if (relIso<0.05 && pt>25 && eta<2.1)	
			return 1;
		else
			return 0;	
	}

int IsoMuonCounter(int NMuon, float Muon_Iso[NMuon], float px[NMuon], float py[NMuon], float pz[NMuon])
	{	int counter=0;	int MuSel=0;
		for(int i=0; i<NMuon; i++)	
		{	MuSel=MuonSelect(Muon_Iso[i], px[i], py[i], pz[i]); 
			if(MuSel)	counter++;
		}
		return counter;
			
	}

float Highest2ndpT(int N, bool JetID[N], float px[N], float py[N], float pz[N])
	{	float pTmax=0;	float pT2nd=0;	int nj;		
			for (int i=0; i<N; i++)
			{	float pt=PT(px[i], py[i]);
				float eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
				if (eta<2.5)	if(pt>pTmax)
				{pT2nd=pTmax;	pTmax=pt;}
			}
		return pT2nd;
	}



float M3abc(float Ea, float pxa, float pya, float pza, float Eb, float pxb, float pyb, float pzb, float Ec, float pxc, float pyc, float pzc)
	{	float E=Ea+Eb+Ec;
		float px=pxa+pxb+pxc;	float py=pya+pyb+pyc;	float pz=pza+pzb+pzc;
		float M3=sqrt(abs(pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2)));
		return M3;
	}

float BestM3abc(int N, bool JetID[N], float Jet_btag[N], float E[N], float px[N], float py[N], float pz[N])
	{
		int length=N*(N-1)/2;
		float eta, x,y,z, pt, transmomenta;
		float trans[N][N][N];	for(int i=0; i<N; i++)	for(int j=0; j<N; j++)	for(int k=0; k<N; k++) {trans[i][j][k]=0;}
		for(int i=0; i<N; i++){
					eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
					pt=PT(px[i], py[i]);
					if (pt>25 && eta<2.5)	{x=px[i];	y=py[i];}					
			for (int j=i+1; j<N; j++){	
					eta=ETA(px[j], py[j], pz[j]);
					pt=PT(px[j], py[j]);	
						if (pt>15 && eta<2.5)	{x=px[i]+px[j];	y=py[i]+py[j];}
						for (int k=j+1; k<N; k++)
						{
							if (pt>5 && eta<2.5)							
							{x=px[i]+px[j]+px[k];	y=py[i]+py[j]+py[k];
							transmomenta=PT(x,y);	
						trans[i][j][k]=transmomenta;}
						}	
						}
					}
		float maxtrans=0;
		int pointeri=0; int pointerj=0; int pointerk=0;		
		for(int i=0; i<N; i++)	
			for(int j=0; j<N; j++)
				for (int k=0; k<N; k++) 
				{	
					if (trans[i][j][k]>maxtrans)	{	pointeri=i; pointerj=j; pointerk=k;	maxtrans=trans[i][j][k];	}
				}
		float bestm3=M3abc(E[pointeri], px[pointeri], py[pointeri], pz[pointeri], E[pointerj], px[pointerj], py[pointerj], pz[pointerj], E[pointerk], px[pointerk], py[pointerk], pz[pointerk]);
		return bestm3;
	}



float TopMass(int N, bool JetID[N], float Jet_btag[N], float E[N], float px[N], float py[N], float pz[N])
	{
		int length=N*(N-1)/2;
		float eta, pt;
		float mt[N][N];	for(int i=0; i<N; i++)	for(int j=0; j<N; j++) {mt[i][j]=0;}
		for(int i=0; i<N; i++){
					eta=ETA(px[i], py[i], pz[i]);	eta=abs(eta);
					pt=PT(px[i], py[i]);					
						if (JetID[i] && Jet_btag[i]<1.76 && eta<2.5 && pt>25)
			for (int j=i+1; j<N; j++){	
					eta=ETA(px[j], py[j], pz[j]);
					pt=PT(px[j], py[j]);	
						if (JetID[j] && Jet_btag[j]<1.76 && eta<2.5 && pt>25)	
						mt[i][j]=MTab(E[i], px[i], py[i], pz[i], E[j], px[j], py[j], pz[j]);
						}	
					}
		float leastdiff=80.385;	//W-Mass
		int pointeri=0; int pointerj=0;		
		for(int i=0; i<N; i++)	
			for(int j=0; j<N; j++) 
			{	float diff;
				diff=mt[i][j]-80.385;	diff=abs(diff);
				if(diff<leastdiff)	{leastdiff=diff;	pointeri=i;	pointerj=j;}
			}
		float top=0;	int pointerk=0;
		for (pointerk=0; pointerk<N;	pointerk++)
		if (Jet_btag[pointerk]>1.76)
			{top=M3abc(E[pointeri], px[pointeri], py[pointeri], pz[pointeri], E[pointerj], px[pointerj], py[pointerj], pz[pointerj], E[pointerk], px[pointerk], py[pointerk], pz[pointerk]);
				break;	}
		return top;

	}



void AllAnalysis()
{

	float Grid[11]={0,25.,50.,75.,100.,125.,160, 210, 250., 350., 500.};
	float Gridht[9]={0, 25, 75, 125, 200, 250, 350, 500, 800};
	float Gridst[9]={0, 25, 100, 200, 250, 300, 350, 450, 800};
	TFile *f1 = new TFile("data.root");
	TTree *t1 = (TTree*)f1->Get("events");
	TH1F *h1a = new TH1F("DataJetpT","Data",30, 0, 450);
	TH1F *h1b = new TH1F("NJet","Data",10, 0., 10.);
	TH1F *h1 = new TH1F("DataTopMass","Data",30, 0, 450);
	TH1F *h1c = new TH1F("DrellYanElectronpT","Drell Yan",30, 0, 450);
	

	TFile *f2 = new TFile("dy.root");
	TTree *t2 = (TTree*)f2->Get("events");
	TH1F *h2a = new TH1F("DrellYanJetpT","Drell Yan",30, 0, 450);
	TH1F *h2b = new TH1F("NJet","Drell Yan",10, 0., 10.);
	TH1F *h2 = new TH1F("DrellYanMuonpT","Drell Yan",30, 0, 450);
	TH1F *h2c = new TH1F("DrellYanElectronpT","Drell Yan",30, 0, 450);

	TFile *f3 = new TFile("qcd.root");
	TTree *t3 = (TTree*)f3->Get("events");
	TH1F *h3a = new TH1F("QCDJetpT","QCD",30, 0, 450);
	TH1F *h3b = new TH1F("QCDNJet","QCD",10, 0., 10.);
	TH1F *h3 = new TH1F("QCDMuonpT","QCD",30, 0, 450);
	TH1F *h3d = new TH1F("QCDElectronpT","Drell Yan",30, 0, 450);

	TFile *f4 = new TFile("single_top.root");
	TTree *t4 = (TTree*)f4->Get("events");
	TH1F *h4a = new TH1F("SingleTopJetpT","Single Top",30, 0, 450);
	TH1F *h4b = new TH1F("SingleTopNJet","SingleTop",10, 0., 10.);
	TH1F *h4 = new TH1F("SingleTopMuonpT","Single Top",30, 0, 450);
	TH1F *h4d = new TH1F("SingleTopElectronpT","Drell Yan",30, 0, 450);

	TFile *f5 = new TFile("ttbar.root");
	TTree *t5 = (TTree*)f5->Get("events");
	TH1F *h5a = new TH1F("ttbarJetpT","ttbar",30, 0, 450);
	TH1F *h5b = new TH1F("ttbarNJet","ttbar",10, 0., 10.);
	TH1F *h5 = new TH1F("ttbarMuonpT","ttbar",30, 0, 450);
	TH1F *h5d = new TH1F("ttbarElectronpT","Drell Yan",30, 0, 450);

	TFile *f6 = new TFile("wjets.root");
	TTree *t6 = (TTree*)f6->Get("events");
	TH1F *h6a = new TH1F("WjetsJetpT","W Jets",30, 0, 450);
	TH1F *h6b = new TH1F("WjetsNJet","WJets",10, 0., 10.);
	TH1F *h6 = new TH1F("WjetsMuonpT","W Jets",30, 0, 450);
	TH1F *h6d = new TH1F("WjetsElectronpT","Drell Yan",30, 0, 450);

	TFile *f7 = new TFile("ww.root");
	TTree *t7 = (TTree*)f7->Get("events");
	TH1F *h7a = new TH1F("WWJetpT","WW",30, 0, 450);
	TH1F *h7b = new TH1F("WWNJet","WW",10, 0., 10.);
	TH1F *h7 = new TH1F("WWMuonpT","WW",30, 0, 450);
	TH1F *h7d = new TH1F("WWElectronpT","Drell Yan",30, 0, 450);

	TFile *f8 = new TFile("wz.root");
	TTree *t8 = (TTree*)f8->Get("events");
	TH1F *h8a = new TH1F("WZJetpT","WZ",30, 0, 450);
	TH1F *h8b = new TH1F("WZNJet","WZ",10, 0., 10.);
	TH1F *h8 = new TH1F("WZMuonpT","WZ",30, 0, 450);
	TH1F *h8d = new TH1F("WZElectronpT","Drell Yan",30, 0, 450);

	TFile *f9 = new TFile("zz.root");
	TTree *t9 = (TTree*)f9->Get("events");
	TH1F *h9a = new TH1F("ZZJetpT","ZZ",30, 0, 450);
	TH1F *h9b = new TH1F("ZZNJet","ZZ",10, 0., 10.);
	TH1F *h9 = new TH1F("ZZMuonpT","ZZ",30, 0, 450);
	TH1F *h9d = new TH1F("ZZElectronpT","Drell Yan",30, 0, 450);	

	THStack *MuonpTStack= new THStack("M3","M3: Top Mass");

	TCanvas *c1= new TCanvas("c1","",800,600);
	float EventWeight, MET_px, MET_py; int NMuon, NJet, n;
	NMuon=1; NJet=1; EventWeight=1.0;	bool triggerIsoMu24, Jet_ID[NJet];	n=1;
	float Muon_Px[NMuon], Muon_Py[NMuon], Muon_Pz[NMuon], Muon_E[NMuon], Muon_ID[NMuon], Muon_Iso[NMuon], Jet_Px[NJet],Jet_Py[NJet],Jet_Pz[NJet],Jet_E[NJet], Jet_btag[NJet] ;
	
	t1-> SetBranchAddress("NJet", &NJet);
	t1-> SetBranchAddress("NMuon", &NMuon);	
	t1-> SetBranchAddress("Muon_Px", &Muon_Px);
	t1-> SetBranchAddress("Muon_Py", &Muon_Py);
	t1-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t1-> SetBranchAddress("Muon_E", &Muon_E);
	t1-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t1-> SetBranchAddress("Jet_Px", &Jet_Px);
	t1-> SetBranchAddress("Jet_Py", &Jet_Py);
	t1-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t1-> SetBranchAddress("Jet_E", &Jet_E);
	t1-> SetBranchAddress("Jet_ID", &Jet_ID);
	t1-> SetBranchAddress("Jet_btag", &Jet_btag);
	t1-> SetBranchAddress("MET_px", &MET_px);
	t1-> SetBranchAddress("MET_py", &MET_py);
	t1-> SetBranchAddress("EventWeight", &EventWeight);
	t1-> SetBranchAddress("triggerIsoMu24",&triggerIsoMu24);
	Int_t n1 = t1->GetEntries();	cout<<n1<<endl;

	for(int i=0; i<n1; i++)
	{	
   		t1->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz);
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz);	
 nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);	

float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);


if (nj>=2 && nb>=1 && met>25 && nm>=1)
	  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h1->Fill(m3, EventWeight);}//cout<<nj<<"\t:"<<JetpT2nd<<endl;}		
	}	
	
	t2-> SetBranchAddress("NJet", &NJet);
	t2-> SetBranchAddress("NMuon", &NMuon);
	t2-> SetBranchAddress("Muon_Px", &Muon_Px);
	t2-> SetBranchAddress("Muon_Py", &Muon_Py);
	t2-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t2-> SetBranchAddress("Muon_E", &Muon_E);
	t2-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t2-> SetBranchAddress("Jet_Px", &Jet_Px);
	t2-> SetBranchAddress("Jet_Py", &Jet_Py);
	t2-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t2-> SetBranchAddress("Jet_E", &Jet_E);
	t2-> SetBranchAddress("Jet_ID", &Jet_ID);
	t2-> SetBranchAddress("Jet_btag", &Jet_btag);
	t2-> SetBranchAddress("MET_px", &MET_px);
	t2-> SetBranchAddress("MET_py", &MET_py);
	t2-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n2 = t2->GetEntries();	cout<<n2<<endl;

	for(int i=0; i<n2; i++)
	{	
   		t2->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, e, pt, p, eta, met;	px=0; py=0; pz=0; float JetpT=0; e=0;  pt=0; eta=0; p=0; met=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);	
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz);
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);

if (nj>=2 && nb>=1 && met>25 && nm>=1)

	{JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h2->Fill(m3, EventWeight);}
	}	

	t3-> SetBranchAddress("NJet", &NJet);
	t3-> SetBranchAddress("NMuon", &NMuon);	
	t3-> SetBranchAddress("Muon_Px", &Muon_Px);
	t3-> SetBranchAddress("Muon_Py", &Muon_Py);
	t3-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t3-> SetBranchAddress("Muon_E", &Muon_E);
	t3-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t3-> SetBranchAddress("Jet_Px", &Jet_Px);
	t3-> SetBranchAddress("Jet_Py", &Jet_Py);
	t3-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t3-> SetBranchAddress("Jet_E", &Jet_E);
	t3-> SetBranchAddress("Jet_ID", &Jet_ID);
	t3-> SetBranchAddress("Jet_btag", &Jet_btag);
	t3-> SetBranchAddress("MET_px", &MET_px);
	t3-> SetBranchAddress("MET_py", &MET_py);
	t3-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n3 = t3->GetEntries();	cout<<n3<<endl;

	for(int i=0; i<n3; i++)
	{	
   		t3->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz); 
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);

	
if (nj>=2 && nb>=1 && met>25 && nm>=1)

	  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h3->Fill(m3, EventWeight);}	
	}
	
	t4-> SetBranchAddress("NJet", &NJet);	
	t4-> SetBranchAddress("NMuon", &NMuon);
	t4-> SetBranchAddress("Muon_Px", &Muon_Px);
	t4-> SetBranchAddress("Muon_Py", &Muon_Py);
	t4-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t4-> SetBranchAddress("Muon_E", &Muon_E);
	t4-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t4-> SetBranchAddress("Jet_Px", &Jet_Px);
	t4-> SetBranchAddress("Jet_Py", &Jet_Py);
	t4-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t4-> SetBranchAddress("Jet_E", &Jet_E);
	t4-> SetBranchAddress("Jet_ID", &Jet_ID);
	t4-> SetBranchAddress("Jet_btag", &Jet_btag);
	t4-> SetBranchAddress("MET_px", &MET_px);
	t4-> SetBranchAddress("MET_py", &MET_py);
	t4-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n4 = t4->GetEntries();	cout<<n4<<endl;

	for(int i=0; i<n4; i++)
	{	
		t4->GetEntry(i);		
	int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz); 
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);

if (nj>=2 && nb>=1 && met>25 && nm>=1)

	
		  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h4->Fill(m3, EventWeight);}	
	}

	t5-> SetBranchAddress("NJet", &NJet);
	t5-> SetBranchAddress("NMuon", &NMuon);
	t5-> SetBranchAddress("Muon_Px", &Muon_Px);
	t5-> SetBranchAddress("Muon_Py", &Muon_Py);
	t5-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t5-> SetBranchAddress("Muon_E", &Muon_E);
	t5-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t5-> SetBranchAddress("Jet_Px", &Jet_Px);
	t5-> SetBranchAddress("Jet_Py", &Jet_Py);
	t5-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t5-> SetBranchAddress("Jet_E", &Jet_E);
	t5-> SetBranchAddress("Jet_ID", &Jet_ID);
	t5-> SetBranchAddress("Jet_btag", &Jet_btag);
	t5-> SetBranchAddress("MET_px", &MET_px);
	t5-> SetBranchAddress("MET_py", &MET_py);
	t5-> SetBranchAddress("EventWeight", &EventWeight);
	t5-> SetBranchAddress("triggerIsoMu24",&triggerIsoMu24);
	Int_t n5 = t5->GetEntries();	cout<<n5<<endl;


	for(int i=0; i<n5; i++)
	{	
   		t5->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); 
	if(triggerIsoMu24)	transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h=0; float s=0; float m3=0; float topm=0;  
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);

float JetpTmax, Jeteta1;	if(triggerIsoMu24) JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta;if(triggerIsoMu24) MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	if(triggerIsoMu24) JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);

if(triggerIsoMu24) {h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz);}
if(triggerIsoMu24)	if (nj>=2 && nb>=1 && met>25 && nm>=1)

	  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       

	h5->Fill(m3, EventWeight);}	
			
	}	

	t6-> SetBranchAddress("NJet", &NJet);	
	t6-> SetBranchAddress("NMuon", &NMuon);
	t6-> SetBranchAddress("Muon_Px", &Muon_Px);
	t6-> SetBranchAddress("Muon_Py", &Muon_Py);
	t6-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t6-> SetBranchAddress("Muon_E", &Muon_E);
	t6-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t6-> SetBranchAddress("Jet_Px", &Jet_Px);
	t6-> SetBranchAddress("Jet_Py", &Jet_Py);
	t6-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t6-> SetBranchAddress("Jet_E", &Jet_E);
	t6-> SetBranchAddress("Jet_ID", &Jet_ID);
	t6-> SetBranchAddress("Jet_btag", &Jet_btag);
	t6-> SetBranchAddress("MET_px", &MET_px);
	t6-> SetBranchAddress("MET_py", &MET_py);
	t6-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n6 = t6->GetEntries();	cout<<n6<<endl;

	for(int i=0; i<n6; i++)
	{	
   		t6->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz); 
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet, Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);
if (nj>=2 && nb>=1 && met>25 && nm>=1)

	  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h6->Fill(m3, EventWeight);}
	}	

	t7-> SetBranchAddress("NJet", &NJet);
	t7-> SetBranchAddress("NMuon", &NMuon);
	t7-> SetBranchAddress("Muon_Px", &Muon_Px);
	t7-> SetBranchAddress("Muon_Py", &Muon_Py);
	t7-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t7-> SetBranchAddress("Muon_E", &Muon_E);
	t7-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t7-> SetBranchAddress("Jet_Px", &Jet_Px);
	t7-> SetBranchAddress("Jet_Py", &Jet_Py);
	t7-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t7-> SetBranchAddress("Jet_E", &Jet_E);
	t7-> SetBranchAddress("Jet_ID", &Jet_ID);
	t7-> SetBranchAddress("Jet_btag", &Jet_btag);
	t7-> SetBranchAddress("MET_px", &MET_px);
	t7-> SetBranchAddress("MET_py", &MET_py);
	t7-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n7 = t7->GetEntries();	cout<<n7<<endl;

	for(int i=0; i<n7; i++)
	{	
   		t7->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz); 
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
if (nj>=2 && nb>=1 && met>25 && nm>=1)

	  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h7->Fill(m3, EventWeight);}	
	}	

	t8-> SetBranchAddress("NJet", &NJet);
	t8-> SetBranchAddress("NMuon", &NMuon);
	t8-> SetBranchAddress("Muon_Px", &Muon_Px);
	t8-> SetBranchAddress("Muon_Py", &Muon_Py);
	t8-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t8-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t8-> SetBranchAddress("Jet_Px", &Jet_Px);
	t8-> SetBranchAddress("Jet_Py", &Jet_Py);
	t8-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t8-> SetBranchAddress("Jet_E", &Jet_E);
	t8-> SetBranchAddress("Jet_ID", &Jet_ID);
	t8-> SetBranchAddress("Jet_btag", &Jet_btag);
	t8-> SetBranchAddress("Muon_E", &Muon_E);
	t8-> SetBranchAddress("MET_px", &MET_px);
	t8-> SetBranchAddress("MET_py", &MET_py);
	t8-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n8 = t8->GetEntries();	cout<<n8<<endl;

	for(int i=0; i<n8; i++)
	{	
   		t8->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz); 
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
if (nj>=2 && nb>=1 && met>25 && nm>=1)	
  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h8->Fill(m3, EventWeight);}	
	}	

	t9-> SetBranchAddress("NJet", &NJet);
	t9-> SetBranchAddress("NMuon", &NMuon);
	t9-> SetBranchAddress("Muon_Px", &Muon_Px);
	t9-> SetBranchAddress("Muon_Py", &Muon_Py);
	t9-> SetBranchAddress("Muon_Pz", &Muon_Pz);
	t9-> SetBranchAddress("Muon_E", &Muon_E);
	t9-> SetBranchAddress("Muon_Iso", &Muon_Iso);	
	t9-> SetBranchAddress("Jet_Px", &Jet_Px);
	t9-> SetBranchAddress("Jet_Py", &Jet_Py);
	t9-> SetBranchAddress("Jet_Pz", &Jet_Pz);
	t9-> SetBranchAddress("Jet_E", &Jet_E);
	t9-> SetBranchAddress("Jet_ID", &Jet_ID);
	t9-> SetBranchAddress("Jet_btag", &Jet_btag);
	t9-> SetBranchAddress("MET_px", &MET_px);
	t9-> SetBranchAddress("MET_py", &MET_py);
	t9-> SetBranchAddress("EventWeight", &EventWeight);
	Int_t n9 = t9->GetEntries();	cout<<n9<<endl;

	for(int i=0; i<n9; i++)
	{	
   		t9->GetEntry(i);
		int j=0; float invmass=0; float transmass=0; float tmass=0;	float px, py, pz, p, eta, pt, met;	px=0; py=0; pz=0; float JetpT=0; eta=0; pt=0; met=0; p=0;
		met=sqrt(pow(MET_px,2)+pow(MET_py,2)); transmass=BestMTab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
 tmass=BestMTcd(NMuon, Muon_Iso, Muon_Px, Muon_Py, Muon_Pz, MET_px, MET_py); invmass=BestM0ab(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz);
		int nb, nj, nm; nb=0; nj=0; nm=0; float h,s, m3, topm; h=HT(NJet, Jet_Px, Jet_Py, Jet_Pz); m3=BestM3abc(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
topm=TopMass(NJet, Jet_ID, Jet_btag, Jet_E, Jet_Px, Jet_Py, Jet_Pz); 
float JetpTmax, Jeteta1;	JetpTmax=HighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz); Jeteta1=etaHighestpT(NJet, Jet_Px, Jet_Py, Jet_Pz);
float MuonpTmax, Muoneta; MuonpTmax=HighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz); Muoneta=etaHighestpT(NMuon, Muon_Px, Muon_Py, Muon_Pz);
float JetpT2nd, Jeteta2;	JetpT2nd=Highest2ndpT(NJet, Jet_ID, Jet_Px, Jet_Py, Jet_Pz);
 s=ST(NJet, Jet_Px, Jet_Py, Jet_Pz, NMuon, Muon_Px, Muon_Py, Muon_Pz); 
 nj=JetCounter(NJet,Jet_ID, Jet_Px, Jet_Py, Jet_Pz); nm=IsoMuonCounter(NMuon,Muon_Iso, Muon_Px, Muon_Py, Muon_Pz);
	nb=bJetCounter(NJet,Jet_ID,Jet_btag, Jet_Px, Jet_Py, Jet_Pz);
if (nj>=2 && nb>=1 && met>25 && nm>=1)

  {JetpTmax=PT(Jet_Px[0], Jet_Py[0]);
 JetpT2nd=PT(Jet_Px[1], Jet_Px[1]);
if(JetpTmax>25 && MuonpTmax>25) if(transmass)       	h9->Fill(m3, EventWeight);}	
	}	

	gStyle->SetOptFit(111111);
	h1->SetLineColor(1);
//	h1->SetFillColor(kRed);
	h1->SetMarkerStyle(20);
	h1->SetMarkerColor(1);

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

	h6->SetLineColor(6);
	h6->SetFillColor(6);
	h6->SetMarkerStyle(21);
	h6->SetMarkerColor(6);

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




	MuonpTStack->Add(h3);
	MuonpTStack->Add(h9);
	MuonpTStack->Add(h8);
	MuonpTStack->Add(h7);
	MuonpTStack->Add(h6);
	MuonpTStack->Add(h2);
	MuonpTStack->Add(h4);
	MuonpTStack->Add(h5);
	
	
	cout<<"\nData: "<<h1->GetEntries();
	cout<<"\nDrell Yan: "<<h2->GetEntries();
	cout<<"\nQCD: "<<h3->GetEntries();
	cout<<"\nSingle Top: "<<h4->GetEntries();
	cout<<"\nttbar: "<<h5->GetEntries();
	cout<<"\nW Jets: "<<h6->GetEntries();
	cout<<"\nWW: "<<h7->GetEntries();
	cout<<"\nWZ: "<<h8->GetEntries();
	cout<<"\nZZ: "<<h9->GetEntries();


	cout<<"\nData: "<<h1->GetEntries()*1;
	cout<<"\nDrell Yan: "<<h2->GetEntries()*0.4389;
	cout<<"\nQCD: "<<h3->GetEntries()*557.5;
	cout<<"\nSingle Top: "<<h4->GetEntries()*0.05482;
	cout<<"\nttbar: "<<h5->GetEntries()*0.2147;
	cout<<"\nW Jets: "<<h6->GetEntries()*1.91;
	cout<<"\nWW: "<<h7->GetEntries()*0.05021;
	cout<<"\nWZ: "<<h8->GetEntries()*0.02077;
	cout<<"\nZZ: "<<h9->GetEntries()*0.003<<endl;
	
	cout<<"\nData: "<<h1->Integral(1,50);
	cout<<"\nDrell Yan: "<<h2->Integral(1,50);
	cout<<"\nQCD: "<<h3->Integral(1,50);
	cout<<"\nSingle Top: "<<h4->Integral(1,50);
	cout<<"\nttbar: "<<h5->Integral(1,50);
	cout<<"\nW Jets: "<<h6->Integral(1,50);
	cout<<"\nWW: "<<h7->Integral(1,50);
	cout<<"\nWZ: "<<h8->Integral(1,50);
	cout<<"\nZZ: "<<h9->Integral(1,50);

	

	//h1->GetYaxis()->SetRangeUser(0,7000);	//gPad->SetLogy();
	h1->GetXaxis()->SetTitle("Mass of Top Quark (GeV/c^{2})");
	h1->GetYaxis()->SetTitle("Events");
	//h1->SetTitle("MuonpTStack");
	c1->cd();
	h1->Draw("pe");
	MuonpTStack->Draw("same HIST");
	//MuonpTStack->Draw("HIST");
	h1->Draw("same pe");
	//gPad->BuildLegend();
	

        c1->Update();
	c1->SaveAs("HT.pdf");	
//	f->Close();
}





