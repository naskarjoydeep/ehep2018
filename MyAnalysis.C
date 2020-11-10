#define MyAnalysis_cxx
// The class definition in MyAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("MyAnalysis.C")
// Root > T->Process("MyAnalysis.C","some options")
// Root > T->Process("MyAnalysis.C+")
//

#include "MyAnalysis.h"
#include <iostream>
#include <TH1F.h>
#include <TLatex.h>

using namespace std;
/////////////////////////////////////////////////////////////////////////////////
//MY Code
/////////////////////////////////////////////////////////////////////////////////


void MyAnalysis::BuildEvent() {
   
   Muons.clear();
   for (int i = 0; i < NMuon; ++i) {
      MyMuon muon(Muon_Px[i], Muon_Py[i], Muon_Pz[i], Muon_E[i]);
      muon.SetIsolation(Muon_Iso[i]);
      muon.SetCharge(Muon_Charge[i]);
      Muons.push_back(muon);
   }
   
   Electrons.clear();
   for (int i = 0; i < NElectron; ++i) {
      MyElectron electron(Electron_Px[i], Electron_Py[i], Electron_Pz[i], Electron_E[i]);
      electron.SetIsolation(Electron_Iso[i]);
      electron.SetCharge(Electron_Charge[i]);
      Electrons.push_back(electron);
   }
   
   Photons.clear();
   for (int i = 0; i < NPhoton; ++i) {
      MyPhoton photon(Photon_Px[i], Photon_Py[i], Photon_Pz[i], Photon_E[i]);
      photon.SetIsolation(Photon_Iso[i]);
      Photons.push_back(photon);
   }
   
   Jets.clear();
   for (int i = 0; i < NJet; ++i) {
      MyJet jet(Jet_Px[i], Jet_Py[i], Jet_Pz[i], Jet_E[i]);
      jet.SetBTagDiscriminator(Jet_btag[i]);
      jet.SetJetID(Jet_ID[i]);
      Jets.push_back(jet);
   }
   
   hadB.SetXYZM(MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, 4.8);
   lepB.SetXYZM(MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz, 4.8);
   hadWq.SetXYZM(MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, 0.0);
   hadWqb.SetXYZM(MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz, 0.0);
   lepWl.SetXYZM(MClepton_px, MClepton_py, MClepton_pz, 0.0);
   lepWn.SetXYZM(MCneutrino_px, MCneutrino_py, MCneutrino_pz, 0.0);
   met.SetXYZM(MET_px, MET_py, 0., 0.);
   
   EventWeight *= weight_factor;
   
}

void MyAnalysis::Begin(TTree * /*tree*/) {
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TString option = GetOption();
   
}

void MyAnalysis::SlaveBegin(TTree * /*tree*/) {
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   
   TString option = GetOption();
   
   h_Mmumu = new TH1F("Mmumu", "Invariant di-muon mass", 60, 60, 120);				//INVARIANT MASS
   h_Mmumu->SetXTitle("M_{#mu#mu}");
   h_Mmumu->Sumw2();
   histograms.push_back(h_Mmumu);
   histograms_MC.push_back(h_Mmumu);

/* 
   h_NMuon = new TH1F("NMuon", "Number of muons", 7, 0, 7);					//NUMBER OF MUONS
   h_NMuon->SetXTitle("No. Muons");
   h_NMuon->Sumw2();
   histograms.push_back(h_NMuon);
   histograms_MC.push_back(h_NMuon);
*/

/*
   h_pTmu = new TH1F("pTmu", "Transverse Momenta", 100, 25, 150);				//TRANSVERSE MOMENTA
   h_pTmu->SetXTitle("pT_{#mu}");
   h_pTmu->Sumw2();
   histograms.push_back(h_pTmu);
   histograms_MC.push_back(h_pTmu);
*/
/*
   h_ymu = new TH1F("ymu", "Rapidity", 20, -10, 10);				//Rapidity
   h_ymu->SetXTitle("y_{#mu}");
   h_ymu->Sumw2();
   histograms.push_back(h_ymu);
   histograms_MC.push_back(h_ymu);
 */

/*
   h_Etmumu = new TH1F("Etmu", "Transverse Energy", 100, 0, 200);				//Transverse Energy
   h_Etmumu->SetXTitle("Et_{#mu}");
   h_Etmumu->Sumw2();
   histograms.push_back(h_Etmumu);
   histograms_MC.push_back(h_Etmumu);
*/
/*
   h_phimu = new TH1F("Phimu", "Angular Distribution", 20, -10, 10);				//Phi
   h_phimu->SetXTitle("#Phi_{#mu}");
   h_phimu->Sumw2();
   histograms.push_back(h_phimu);
   histograms_MC.push_back(h_phimu);
*/
 /*  h_pTmu = new TH1F("pTmu", "Transverse Momenta", 100, 25, 150);				//TRANSVERSE MOMENTA
   h_pTmu->SetXTitle("pT_{#mu}");
   h_pTmu->Sumw2();
   histograms.push_back(h_pTmu);
   histograms_MC.push_back(h_pTmu);
   
   h_pTe = new TH1F("pTmu", "Transverse Momenta", 100, 25, 150);				//TRANSVERSE MOMENTA
   h_pTe->SetXTitle("pT_{#mu}");
   h_pTe->Sumw2();
   histograms.push_back(h_pTe);
   histograms_MC.push_back(h_pTe);
*/
/* 
   h_pTjet = new TH1F("pTjet", "Transverse Momenta", 50, 25, 250);				//TRANSVERSE MOMENTA
   h_pTjet->SetXTitle("p_{T} b-Jet (GeV/c^{2})");
   h_pTjet->Sumw2();
   histograms.push_back(h_pTjet);
   histograms_MC.push_back(h_pTjet);
*/
/*
   h_pTall = new TH1F("pTmu", "Transverse Momenta", 100, 25, 150);				//TRANSVERSE MOMENTA
   h_pTall->SetXTitle("pT_{#mu}");
   h_pTall->Sumw2();
   histograms.push_back(h_pTall);
   histograms_MC.push_back(h_pTall);

   h_NMuon = new TH1F("NJet", "Number of bJets", 7, 0, 7);					//NUMBER OF MUONS
   h_NMuon->SetXTitle("No. b-Jets");
   h_NMuon->Sumw2();
   histograms.push_back(h_NMuon);
   histograms_MC.push_back(h_NMuon);

   h_pTjet = new TH1F("m0", "Transverse Mass of W-boson", 50, 25, 250);				//TRANSVERSE MOMENTA
   h_pTjet->SetXTitle("m_{T} W-boson (GeV/c^{2})");
   h_pTjet->Sumw2();
   histograms.push_back(h_pTjet);
   histograms_MC.push_back(h_pTjet);
*/

   
}
int counter=0;	int runi=0;	int runj=0;		//I ADDED IT
Bool_t MyAnalysis::Process(Long64_t entry) {
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either MyAnalysis::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   
   ++TotalEvents;
   
   GetEntry(entry);
   
   if (TotalEvents % 10000 == 0)
      cout << "Next event -----> " << TotalEvents << endl;
   
   BuildEvent();
   
   double MuonPtCut = 25.;
   double MuonRelIsoCut = 0.05;
   
/* 
    cout << "Jets: " << endl;
      for (vector<MyJet>::iterator it = Jets.begin(); it != Jets.end(); ++it) {
         cout << "pt, eta, phi, btag, id: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", " << it->IsBTagged() << ", " << it->GetJetID()
         << endl;
      }
      cout << "Muons: " << endl;
      for (vector<MyMuon>::iterator it = Muons.begin(); it != Muons.end(); ++it) {
         cout << "pt, eta, phi, iso, charge: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", "
         << it->GetIsolation() << ", " << it->GetCharge() << endl;
      }
      cout << "Electrons: " << endl;
      for (vector<MyElectron>::iterator it = Electrons.begin(); it != Electrons.end(); ++it) {
         cout << "pt, eta, phi, iso, charge: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", "
         << it->GetIsolation() << ", " << it->GetCharge() << endl;
      }
      cout << "Photons: " << endl;
      for (vector<MyPhoton>::iterator it = Photons.begin(); it != Photons.end(); ++it) {
         cout << "pt, eta, phi, iso: " << it->Pt() << ", " << it->Eta() << ", " << it->Phi() << ", " << it->GetIsolation()
         << endl;
      }
   */
   
   //////////////////////////////



// ************************************************************			MY CODE STARTS HERE		*******************************************************************************





   // Exercise 1: Invariant Di-Muon mass
   
   int N_IsoMuon = 0;
   MyMuon *muon1, *muon2;
   float Z0mass;
  runi++; 
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {	runj++;
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
         if (N_IsoMuon == 2) muon2 = &(*jt);
      }
   }
   
   //h_NMuon->Fill(N_IsoMuon, EventWeight);
   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (muon1->Pt()>MuonPtCut && muon2->Pt()>MuonPtCut)
	    if (abs(muon1->Rapidity())<2.1 && abs(muon2->Rapidity())<2.1 )
	 {	Z0mass=(*muon1 + *muon2).M();
         h_Mmumu->Fill(Z0mass, EventWeight);					//M() gives invariant mass
      }
   }	 
   //////////////////////////////
   
   return kTRUE;
}



/*
//////////////////////////////
   // Exercise 2: Kinematic Variables		MUON
	// 2.1 Transverse Momenta pT
   
   int N_IsoMuon = 0;	int N_IsoElectron = 0;
   MyMuon *muon1, *muon2;	
   
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
	 if (N_IsoMuon == 2) muon2 = &(*jt);
      }
   }

   h_NMuon->Fill(N_IsoMuon, EventWeight);

   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (muon1->Pt()>MuonPtCut) {
         h_pTmu->Fill((*muon1).Pt(), EventWeight);					//Perp() OR Pt() gives transverse component of momenta
      }
   }


   //////////////////////////////
   
   return kTRUE;
}
*/
/*

//////////////////////////////
   // Exercise 2: Kinematic Variables	JET
	// 2.1 Transverse Momenta pT
   
   int N_IsoMuon = 0;	int N_IsoElectron = 0;
   MyJet *Jet1, *Jet2;	
   
   for (vector<MyJet>::iterator jt = Jets.begin(); jt != Jets.end(); ++jt) {
      if (jt->GetBTagDiscriminator()>0.8) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) Jet1 = &(*jt);
	 if (N_IsoMuon == 2) Jet2 = &(*jt);
      }
   }

   h_NMuon->Fill(N_IsoMuon, EventWeight);

   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (Jet1->Pt()>25.0) {
         h_pTjet->Fill((*Jet1).Pt(), EventWeight);					//Perp() OR Pt() gives transverse component of momenta
      }
   }


   //////////////////////////////
   
   return kTRUE;
}
*/

/*
//////////////////////////////
   // Exercise X: Mass of Top		incomplete
   
   int N_IsoMuon = 0;
   MyJet *Jet1, *Jet2;
   
   for (vector<MyJet>::iterator jt = Jets.begin(); jt != Jets.end(); ++jt) {
      if (jt->IsBTagged()) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) Jet1 = &(*jt);
	 if (N_IsoMuon == 2) Jet2 = &(*jt);
      }
   }

   h_NMuon->Fill(N_IsoMuon, EventWeight);

   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (Jet1->Pt()>25.0 && Jet1->GetJetID()) {
         h_pTjet->Fill((*Jet1).Pt(), EventWeight);					//Perp() OR Pt() gives transverse component of momenta
      }
   }


   //////////////////////////////
   
   return kTRUE;
}
*/






/*

//////////////////////////////				
   // Exercise 2: Kinematic Variables
	// 2.2 Rapidity
   
   int N_IsoMuon = 0;	
   MyMuon *muon1;
   
	
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
      }
   }
      
   h_NMuon->Fill(N_IsoMuon, EventWeight);
   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (muon1->Pt()>MuonPtCut) {
         h_ymu->Fill((*muon1).Rapidity(), EventWeight);					//Rapidity() gives Rapidity	
			  		 }

   }
   //////////////////////////////
   
   return kTRUE;
}
*/

/*
/////////////////////////////				
   // Exercise 2: Kinematic Variables
	// 2.3 Transverse Energy
   
   int N_IsoMuon = 0;	
   MyMuon *muon1, *muon2;	
   
	
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
	 if (N_IsoMuon == 2) muon2 = &(*jt);
      }
   }
  
   h_NMuon->Fill(N_IsoMuon, EventWeight);
   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (muon1->Pt()>MuonPtCut) {
         h_Etmumu->Fill((*muon1+ *muon2).Et(), EventWeight);					//Et() gives Transverse Energy	
			  		 }

   }
   //////////////////////////////
   
   return kTRUE;
}
*/
/*
//////////////////////////////				
   // Exercise 2: Kinematic Variables
	// 2.2 Phi
   
   int N_IsoMuon = 0;	
   MyMuon *muon1;
   
	
   for (vector<MyMuon>::iterator jt = Muons.begin(); jt != Muons.end(); ++jt) {
      if (jt->IsIsolated(MuonRelIsoCut)) {
         ++N_IsoMuon;
         if (N_IsoMuon == 1) muon1 = &(*jt);
      }
   }
      
   h_NMuon->Fill(N_IsoMuon, EventWeight);
   
   if (N_IsoMuon > 1 && triggerIsoMu24) {
      if (muon1->Pt()>MuonPtCut) {
         h_phimu->Fill((*muon1).Phi(), EventWeight);					//Phi() gives Azimuthal Angle	
			  		 }

   }
   //////////////////////////////
   
   return kTRUE;
}

*/


//////////////////////////////
   // Exercise 2: Kinematic Variables	JET
	// 2.1 Transverse Momenta pT
 /*  
   int N_bJet = 0;	int N_Jet=0;	int N_nonbJet=0;	
   MyJet *Jet1, *Jet2;	
   
   for (vector<MyJet>::iterator jt = Jets.begin(); jt != Jets.end(); ++jt) {
      if(jt->GetJetID()){++N_Jet; if (jt->IsBTagged())++N_bJet; 
	else if (!jt->IsBTagged()) {	++N_nonbJet;					//GetBTagDiscriminator()>1.76 selects b-Jets only	
         if (N_nonbJet == 1) Jet1 = &(*jt);
	 if (N_nonbJet == 2) Jet2 = &(*jt);
      }	}
   }

   
	
   
   if (N_Jet>=2 && N_bJet>= 1 && triggerIsoMu24) {
      if (Jet1->Pt()>25.0 && abs(Jet1->Rapidity())<2.5) {
	 h_NMuon->Fill(N_bJet, EventWeight);
         h_pTjet->Fill((*Jet1+*Jet2).Mt(), EventWeight);					//Perp() OR Pt() gives transverse component of momenta
      }
   }


   //////////////////////////////
   
   return kTRUE;
}


*/


























void MyAnalysis::SlaveTerminate() {
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   
}

void MyAnalysis::Terminate() {
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
}
