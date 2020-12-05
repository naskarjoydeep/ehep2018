#include "Pythia8/Pythia.h"
// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("ZZjets", &argc, argv);

  // Number of events, generated and listed ones.
  int nEvent    = 50000;			//nEvent=500;
  int nListJets = 5;

  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("WeakDoubleBoson:ffbar2gmZgmZ = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("WeakZ0:gmZmode = 2");
  pythia.init();

  // Common parameters for the two jet finders.
  double etaMax   = 2.5;
  double radius   = 0.7;
  double pTjetMin = 10.;
  // Exclude neutrinos (and other invisible) from stuZZ.
  int    nSel     = 2;
  // Range and granularity of CellJet jet finder.
  int    nEta     = 80;
  int    nPhi     = 64;

  // Set up SlowJet jet finder, with anti-kT clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

  // Set up CellJet jet finder.
  //CellJet cellJet( etaMax, nEta, nPhi, nSel);

 // Create file on which histogram(s) can be saved.
  TFile* outFile1 = new TFile("ZZnJetsS.root", "RECREATE");
  //TFile* outFile2 = new TFile("ZZnJetsC.root", "RECREATE");
  //TFile* outFile3 = new TFile("ZZnJetsD.root", "RECREATE");
  TFile* outFile4 = new TFile("ZZpTJetsS.root", "RECREATE");
  //TFile* outFile5 = new TFile("ZZeTJetsC.root", "RECREATE");
  TFile* outFile6 = new TFile("ZZyJetsS.root", "RECREATE");
  //TFile* outFile7 = new TFile("ZZetaJetsC.root", "RECREATE");
  TFile* outFile8 = new TFile("ZZphiJetsS.root", "RECREATE");
  //TFile* outFile9 = new TFile("ZZphiJetsC.root", "RECREApTTE");
  TFile* outFile10 = new TFile("ZZdistJetsS.root", "RECREATE");
  //TFile* outFile11 = new TFile("ZZdistJetsC.root", "RECREATE");
  TFile* outFile12 = new TFile("ZZpTdiffS.root", "RECREATE");
  //TFile* outFile13 = new TFile("ZZeTdiffC.root", "RECREATE");
  TFile* outFile14 = new TFile("ZZnLeptons.root", "RECREATE");
  TFile* outFile15 = new TFile("ZZLeptonpT.root", "RECREATE");
  TFile* outFile16 = new TFile("ZZLeptonEta.root", "RECREATE");
  TFile* outFile17 = new TFile("ZZMET.root", "RECREATE");


// Book histogram.
  TH1F *nJetsS = new TH1F("ZZnJetsS","ZZ number of jets, SlowJet", 50, -0.5, 49.5);
  //TH1F *nJetsC = new TH1F("ZZnJetsC","ZZ number of jets, CellJet", 50, -0.5, 49.5);
  //TH1F *nJetsD = new TH1F("ZZnJetsD","ZZ number of jets, CellJet - SlowJet", 45, -22.5, 22.5);
  TH1F *pTJetsS = new TH1F("ZZpTJetsS","ZZ pT for jets, SlowJet", 100, 0., 250.);
  //TH1F *eTJetsC = new TH1F("ZZeTJetsC","ZZ eT for jets, CellJet", 100, 0., 250.);
  TH1F *yJetsS = new TH1F("ZZyJetsS","ZZ y for jets, SlowJet", 100, -5., 5.);
  //TH1F *etaJetsC = new TH1F("ZZyJetsC","ZZ eta for jets, CellJet", 100, -5., 5.);
  TH1F *phiJetsS = new TH1F("ZZphiJetsS","ZZ phi for jets, SlowJet", 100, -M_PI, M_PI);
  //TH1F *phiJetsC = new TH1F("ZZphiJetsC","ZZ phi for jets, CellJet", 100, -M_PI, M_PI);
  TH1F *distJetsS = new TH1F("ZZdistJetsS","ZZ R distance between jets, SlowJet", 100, 0., 10.);
  //TH1F *distJetsC = new TH1F("ZZdistJetsC","ZZ R distance between jets, CellJet", 100, 0., 10.);
  TH1F *pTdiffS = new TH1F("ZZpTdiffS","ZZ pT difference, SlowJet", 100, -100., 400.);
  //TH1F *eTdiffC = new TH1F("ZZeTdiffC","ZZ eT difference, CellJet", 100, -100., 400.);
  TH1F *nLepton = new TH1F("ZZnLepton","ZZ number of leptons", 20, -0.5, 19.5);
  //TH1F *m0JetsC = new TH1F("ZZm0","ZZ Invariant Mass of Jets, CellJet", 100, 0., 1000.);
  TH1F *pTLepton = new TH1F("ZZpTLepton","ZZ pT of Leptons", 100, 0., 250.);
  TH1F *etaLepton = new TH1F("ZZetaLepton","ZZ eta of Leptons", 20, -10., 10.);
  TH1F *MET = new TH1F("ZZMET","ZZ Missing Energy", 100, 0., 100.);



 // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {	
    if (!pythia.next()) continue;

    // Fill number of leptons
    int nLep = 0;	float pTmax=0;	int pointer=0;	float Leta=0;	float Ex=0;	float Ey=0;	float Ei=0;	float eta=0;	float phi=0;	float missET=0;
    for (int i = 0; i < pythia.event.size(); ++i)	
      if ( pythia.event[i].isFinal() && (pythia.event[i].isCharged() || pythia.event[i].id()==22) )
								 { 	Ei=pythia.event[i].e();	eta=pythia.event[i].eta();	phi=pythia.event[i].phi();
									Ex+=Ei*cos(phi)/cosh(eta);	Ey+=Ei*sin(phi)/cosh(eta);

								if( pythia.event[i].isLepton() )
        								{	++nLep;
										if(pythia.event[i].pT()>pTmax)
											{	pointer=i;	pTmax=pythia.event[i].pT();	Leta=pythia.event[i].eta();	}
									}
								}
	missET=pow(Ex,2)+ pow(Ey,2);	missET=sqrt(missET);
   	nLepton->Fill( nLep );
	if(pTmax>25)	pTLepton->Fill (pTmax);
   	etaLepton->Fill (Leta);
	MET->Fill(missET);
   						
   // Fill pT and eta of leptons

   	

    // Analyze Slowet jet properties. List first few.
    slowJet. analyze( pythia.event );
    if (iEvent < nListJets) slowJet.list();

    // Fill SlowJet inclusive jet distributions.
    nJetsS->Fill( slowJet.sizeJet() );
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      pTJetsS->Fill( slowJet.pT(i) );
      yJetsS->Fill( slowJet.y(i) );
      phiJetsS->Fill( slowJet.phi(i) );
    }

    // Fill SlowJet distance between jets.
    for (int i = 0; i < slowJet.sizeJet() - 1; ++i)
    for (int j = i +1; j < slowJet.sizeJet(); ++j) {
      double dEta = slowJet.y(i) - slowJet.y(j);
      double dPhi = abs( slowJet.phi(i) - slowJet.phi(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJetsS->Fill( dR );
    }

    // Fill SlowJet pT-difference between jets (to check ordering of list).
    for (int i = 1; i < slowJet.sizeJet(); ++i)
      pTdiffS->Fill( slowJet.pT(i-1) - slowJet.pT(i) );

    // Analyze CellJet jet properties. List first few.
 //   cellJet. analyze( pythia.event, pTjetMin, radius );
 //   if (iEvent < nListJets) cellJet.list();

    // Fill CellJet inclusive jet distributions.
   // nJetsC->Fill( cellJet.size() );
   // for (int i = 0; i < cellJet.size(); ++i) {
   //   m0JetsC->Fill( cellJet.m(i) );
   //   eTJetsC->Fill( cellJet.eT(i) );
   //   etaJetsC->Fill( cellJet.etaWeighted(i) );
   //   phiJetsC->Fill( cellJet.phiWeighted(i) );
    //}

    // Fill CellJet distance between jets.
  //  for (int i = 0; i < cellJet.size() - 1; ++i)
  //  for (int j = i +1; j < cellJet.size(); ++j) {
  //    double dEta = cellJet.etaWeighted(i)
  //      - cellJet.etaWeighted(j);
  //    double dPhi = abs( cellJet.phiWeighted(i)
  //      - cellJet.phiWeighted(j) );
 //     if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
 //     double dR = sqrt( pow2(dEta) + pow2(dPhi) );
 //     distJetsC->Fill( dR );
 //   }

    // Fill CellJet ET-difference between jets (to check ordering of list).
 //   for (int i = 1; i < cellJet.size(); ++i)
 //     eTdiffC->Fill( cellJet.eT(i-1) - cellJet.eT(i) );

    // Compare number of jets for the two finders.
 //   nJetsD->Fill( cellJet.size() - slowJet.sizeJet() );

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();

// Show histogram. Possibility to close it.

  //m0JetsC->Draw(); 
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  nLepton->Draw(); 
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive(); 
  pTLepton->Draw(); 
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive(); 
  etaLepton->Draw(); 
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  MET->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  nJetsS->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  //nJetsC->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  //nJetsD->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  pTJetsS->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  //eTJetsC->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  yJetsS->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  //etaJetsC->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  phiJetsS->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  //phiJetsC->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  distJetsS->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  //distJetsC->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  pTdiffS->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  //eTdiffC->Draw();
  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();
  


  // Save histogram on file and close file.
  //m0JetsC->Write();
  nLepton->Write();
  pTLepton->Write();
  etaLepton->Write();
  nJetsS->Write();
 // nJetsC->Write();
  //nJetsD->Write();
  pTJetsS->Write();
  //eTJetsC->Write();
  yJetsS->Write();
  //etaJetsC->Write();
  phiJetsS->Write();
 // phiJetsC->Write();
  distJetsS->Write();
  //distJetsC->Write();
  pTdiffS->Write();
  //eTdiffC->Write();


  delete outFile1;
  //delete outFile2;
 // delete outFile3;
  delete outFile4;
 // delete outFile5;
  delete outFile6;
  //delete outFile7;
  delete outFile8;
  //delete outFile9;
  delete outFile10;
 // delete outFile11;
  delete outFile12;
 // delete outFile13;
  delete outFile14;
  delete outFile15;
  delete outFile16;
  delete outFile17;
  

// Done.
  return 0;
}

