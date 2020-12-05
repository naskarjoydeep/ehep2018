// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

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
  TApplication theApp("histpTe-100bin", &argc, argv);

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  pythia.readString("Top:all= on");
 //pythia.readString("Top:gg2ttbar = on");
 //pythia.readString("Top:qqbar2ttbar = on");
 //pythia.readString("WeakSingleBoson:ffbar2W"); 
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Beams:eCM = 13000.");
  
  pythia.readString("6:onMode = off");
  pythia.readString("6:onIfAny = 24 5");
  pythia.readString("-6:onMode = off");
  pythia.readString("-6:onIfAny = -24 -5");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = -11 14");
  pythia.readString("-24:onMode = off");
  pythia.readString("-24:onIfAny = 11 -14");
  pythia.init();

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("TopAll.root", "RECREATE");

  // Book histogram.
  TH1F *pTpi = new TH1F("TopAllpTe","TopAll dN/dpT", 100, 0., 200.);
int counter=0;
  // Begin event loop. Generate event; skip if generation aborted.
int iEvent=0;  
for (int iEvent = 0; iEvent < 10000; ++iEvent){
//													while(counter<10000){
    if (!pythia.next()) continue;

 // Loop over particles in event. Find last Z0 copy.
    int ipi = 0;	int pTmax=pythia.event[0].pT();
    for (int i = 0; i < pythia.event.size(); ++i)					
	{						int mother1= pythia.event[i].mother1();	int mother2=pythia.event[i].mother2();

							if ((pythia.event[i].id()==-11 && pythia.event[mother1].id()==24) || (pythia.event[i].id()==11 && pythia.event[mother1].id()==-24))	
								if (pythia.event[i].isFinal()) if(pythia.event[i].pT()>pTmax)
											{pTmax=pythia.event[i].pT();	ipi = i;	counter++;}
	}
    // Fill its pT
	if(ipi)
    pTpi->Fill( pythia.event[ipi].pT() );
  }
  // Statistics on event generation.
  pythia.stat();

  // Show histogram. Possibility to close it.
  pTpi->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();

  // Save histogram on file and close file.
  pTpi->Write();
  delete outFile;
cout<<"\n\n"<<counter<<endl;
  // Done.
  return 0;
}
