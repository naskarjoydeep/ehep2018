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
#include "tgmath.h"

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
  TApplication theApp("dyZhistmTe-", &argc, argv);

   // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
 pythia.readString("Top:gg2ttbar = on");
 pythia.readString("Top:qqbar2ttbar = on");
 //pythia.readString("Top:all = on");
 //pythia.readString("WeakSingleBoson:ffbar2W"); 
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Beams:eCM = 13000.");
  
 
  pythia.readString("6:onMode = off");
  pythia.readString("6:onIfAny = 24 5");
  pythia.readString("-6:onMode = off");
  pythia.readString("-6:onIfAny = -24 -5");
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = -11 -13");
  pythia.readString("-24:onMode = off");
  pythia.readString("-24:onIfAny = 11 13");
  //pythia.readString("-5:onMode = off");
  //pythia.readString("5:onMode = off");



  pythia.init();


  // Create file on which histogram(s) can be saved.
 
    TFile* outFile1 = new TFile("m0TOPALL.root", "RECREATE");
   
  // Book histogram.
  
    TH1F *m0 = new TH1F("Invariant Mass","dN/m0", 100, 0., 150.);	
   
   Hist im("Invariant Mass", 0, 0., 150.);
   
								int mom1, mom2, d2;
								double pTmax;
								double mass;
								double pf[3];	
								double Ei, Ef;
								Vec4 p4f;
								//miss4mass=sqrt(pi.m2Calc()-pf.m2Calc());
								//missmass=sqrt(pow((Ei-Ef),2)-pow((pi[0]-pf[0]),2)- pow((pi[1]-pf[1]),2)-pow((pi[1]-pf[1]),2));
								int pointer;


								int counter=0;
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < 100; ++iEvent) {	mass=0; for (int m=0; m<3; m++) {	pf[m]=0;}	 p4f.reset();	pointer=0;
    if (!pythia.next()) continue;
							
     // Loop over particles in event. Find last Z0 copy.
    pointer = 0;			pTmax=pythia.event[0].pT();
    for (int i = 0; i < pythia.event.size(); ++i)
	 {							
								
						
					int mother1=pythia.event[i].mother1();	int daughter2=pythia.event[mother1].daughter2();
					if (daughter2==i)	daughter2=pythia.event[mother1].daughter1();
															
						if (pythia.event[i].isFinal() && (pythia.event[i].isCharged() || pythia.event[i].id()==22))				
								{
										
							if(pythia.event[i].pT()>pTmax)
											{	pointer=i;	pTmax=pythia.event[i].pT();			
										

									Ef=pythia.event[i].e()+pythia.event[daughter2].e();
									p4f=pythia.event[i].p()+pythia.event[daughter2].p();					
									pf[0]=pythia.event[i].px()+pythia.event[daughter2].px();			
									pf[1]=pythia.event[i].py()+pythia.event[daughter2].py();			
									pf[2]=pythia.event[i].pz()+pythia.event[daughter2].pz();			
									mass=abs(pow(Ef,2)-pow(pf[0],2)- pow(pf[1],2)-pow(pf[2],2));
									mass=sqrt(mass);
								}}	

							

	} if(mass)	counter++;	mom1= pythia.event[pointer].mother1();		mom2= pythia.event[pointer].mother2();	int d2=pythia.event[mom1].daughter2();
													if (d2==pointer)	d2=pythia.event[mom1].daughter1();
					
if (mass<10) 
{	cout<<"\nMass "<<mass<<" with pT "<< pythia.event[pointer].pT()<<" coming from "<<pythia.event[pointer].mother1()<<" : "<<pythia.event[mom1].id()<<" & "<<pythia.event[pointer].mother2()<<" : "<<pythia.event[mom2].id();
				cout<<"\nOther is Daughter "<<d2<<": "<<pythia.event[d2].id()<<" with pT "<<pythia.event[d2].pT()<<endl;
}
     
    // Fill its Invariant Mass.
	if (pointer)
	{	m0->Fill(mass);
       		im.fill(mass);
	}  
}
  // Statistics on event generation.
  pythia.stat();
  // Show histogram. Possibility to close it.
 
  cout << im;
 
  m0->Draw();
  std::cout << "\nDouble click on the histogram window to quit.\n";
  gPad->WaitPrimitive();
  

  // Save histogram on file and close file.
  m0->Write();
  
  delete outFile1;
cout<<"\n\n"<<counter<<endl;
  // Done.
  return 0;
}
