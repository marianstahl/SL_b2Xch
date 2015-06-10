/*

    2015-05-28 M. Stahl

    Simple ACLIC compilable script to add new variables to the baryon spectroscopy tuples.

    Use e.g. root -l scripts/NewVariables.cpp++ to run.

*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <exception>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "Riostream.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "RooRealVar.h"
#include "RooDataHist.h"

TString temp;

void NewVariables(){

  const double protonmass = 938.272013; //MeV
  const double pionmass = 139.57018; //MeV
  const double kaonmass = 493.677; //MeV
  //const double muonmass = 105.6583715; //MeV

  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  double p_PT, p_ETA, p_PHI;
  double K_PT, K_ETA, K_PHI;
  double pi_PT, pi_ETA, pi_PHI;
  double Xb_OWNPV_X, Xb_OWNPV_Y, Xb_OWNPV_Z;
  double Xb_ENDVERTEX_X, Xb_ENDVERTEX_Y, Xb_ENDVERTEX_Z;
  double Xb_PT, Xb_ETA, Xb_PHI, Xb_M;
  double Xc_PT, Xc_ETA, Xc_PHI, Xc_M;
  float Added_H_PT[200], Added_H_ETA[200], Added_H_PHI[200];
  int Added_n_Particles;

  gErrorIgnoreLevel = kError;
  TFile *fSLBS = new TFile("/auto/data/mstahl/SLBaryonSpectroscopy/SLBaryonSpectroscopyStrp21.root","read");
  TTree *Xic_tree = (TTree*)gDirectory->Get("Xib02XicMuNu/Xic2pKpi/DecayTree");
  gErrorIgnoreLevel = kPrint;
  Xic_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)  
  Xic_tree->SetBranchStatus("Xib_M",1);
  Xic_tree->SetBranchStatus("Xib_PT",1);
  Xic_tree->SetBranchStatus("Xib_ETA",1);
  Xic_tree->SetBranchStatus("Xib_PHI",1);
  Xic_tree->SetBranchStatus("Xib_OWNPV_X",1);
  Xic_tree->SetBranchStatus("Xib_OWNPV_Y",1);
  Xic_tree->SetBranchStatus("Xib_OWNPV_Z",1);
  Xic_tree->SetBranchStatus("Xib_ENDVERTEX_X",1);
  Xic_tree->SetBranchStatus("Xib_ENDVERTEX_Y",1);
  Xic_tree->SetBranchStatus("Xib_ENDVERTEX_Z",1);

  Xic_tree->SetBranchStatus("Xic_M",1);
  Xic_tree->SetBranchStatus("Xic_PT",1);
  Xic_tree->SetBranchStatus("Xic_ETA",1);
  Xic_tree->SetBranchStatus("Xic_PHI",1);

  Xic_tree->SetBranchStatus("Added_n_Particles",1);
  Xic_tree->SetBranchStatus("Added_H_PT",1);
  Xic_tree->SetBranchStatus("Added_H_ETA",1);
  Xic_tree->SetBranchStatus("Added_H_PHI",1);

  Xic_tree->SetBranchStatus("p_PT",1);
  Xic_tree->SetBranchStatus("p_ETA",1);
  Xic_tree->SetBranchStatus("p_PHI",1);
  Xic_tree->SetBranchStatus("K_PT",1);
  Xic_tree->SetBranchStatus("K_ETA",1);
  Xic_tree->SetBranchStatus("K_PHI",1);
  Xic_tree->SetBranchStatus("pi_PT",1);
  Xic_tree->SetBranchStatus("pi_ETA",1);
  Xic_tree->SetBranchStatus("pi_PHI",1);

  //set the branch addresses
  Xic_tree->SetBranchAddress("Xib_M",&Xb_M);
  Xic_tree->SetBranchAddress("Xib_PT",&Xb_PT);
  Xic_tree->SetBranchAddress("Xib_ETA",&Xb_ETA);
  Xic_tree->SetBranchAddress("Xib_PHI",&Xb_PHI);
  Xic_tree->SetBranchAddress("Xib_OWNPV_X",&Xb_OWNPV_X);
  Xic_tree->SetBranchAddress("Xib_OWNPV_Y",&Xb_OWNPV_Y);
  Xic_tree->SetBranchAddress("Xib_OWNPV_Z",&Xb_OWNPV_Z);
  Xic_tree->SetBranchAddress("Xib_ENDVERTEX_X",&Xb_ENDVERTEX_X);
  Xic_tree->SetBranchAddress("Xib_ENDVERTEX_Y",&Xb_ENDVERTEX_Y);
  Xic_tree->SetBranchAddress("Xib_ENDVERTEX_Z",&Xb_ENDVERTEX_Z);

  Xic_tree->SetBranchAddress("Xic_M",&Xc_M);
  Xic_tree->SetBranchAddress("Xic_PT",&Xc_PT);
  Xic_tree->SetBranchAddress("Xic_ETA",&Xc_ETA);
  Xic_tree->SetBranchAddress("Xic_PHI",&Xc_PHI);

  Xic_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Xic_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Xic_tree->SetBranchAddress("Added_H_ETA",&Added_H_ETA);
  Xic_tree->SetBranchAddress("Added_H_PHI",&Added_H_PHI);

  Xic_tree->SetBranchAddress("p_PT",&p_PT);
  Xic_tree->SetBranchAddress("p_ETA",&p_ETA);
  Xic_tree->SetBranchAddress("p_PHI",&p_PHI);
  Xic_tree->SetBranchAddress("K_PT",&K_PT);
  Xic_tree->SetBranchAddress("K_ETA",&K_ETA);
  Xic_tree->SetBranchAddress("K_PHI",&K_PHI);
  Xic_tree->SetBranchAddress("pi_PT",&pi_PT);
  Xic_tree->SetBranchAddress("pi_ETA",&pi_ETA);
  Xic_tree->SetBranchAddress("pi_PHI",&pi_PHI);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory

  double Xb_CorrM, p_beta, K_beta, pi_beta;
  float Xcpi_CosTheta[200],XcK_CosTheta[200],Xcp_CosTheta[200];
  double p_as_piKpi_M, p_as_KKpi_M, pK_as_pipi_M, pK_as_ppi_M, pKpi_as_K_M, pKpi_as_p_M;

  TFile *f1 = new TFile("/auto/data/mstahl/SLBaryonSpectroscopy/SLBaryonSpectroscopyStrp21_friend.root","RECREATE");
  //f1->mkdir("Xib02XicMuNu/Xic2pKpi");
  //f1->cd("Xib02XicMuNu/Xic2pKpi");
  TTree added_Xic_tree("Xic2pKpi","Xic2pKpi");

  added_Xic_tree.Branch("Xib_CorrM", &Xb_CorrM, "Xib_CorrM/D");
  added_Xic_tree.Branch("p_beta", &p_beta, "p_beta/D");
  added_Xic_tree.Branch("K_beta", &K_beta, "K_beta/D");
  added_Xic_tree.Branch("pi_beta", &pi_beta, "pi_beta/D");
  added_Xic_tree.Branch("Added_n_Particles", &Added_n_Particles, "Added_n_Particles/I");
  added_Xic_tree.Branch("Xcpi_CosTheta", &Xcpi_CosTheta, "Xcpi_CosTheta[Added_n_Particles]/F");
  added_Xic_tree.Branch("XcK_CosTheta", &XcK_CosTheta, "XcK_CosTheta[Added_n_Particles]/F");
  added_Xic_tree.Branch("Xcp_CosTheta", &Xcp_CosTheta, "Xcp_CosTheta[Added_n_Particles]/F");
  added_Xic_tree.Branch("p_as_piKpi_M", &p_as_piKpi_M, "p_as_piKpi_M/D");
  added_Xic_tree.Branch("p_as_KKpi_M", &p_as_KKpi_M, "p_as_KKpi_M/D");
  added_Xic_tree.Branch("pK_as_pipi_M", &pK_as_pipi_M, "pK_as_pipi_M/D");
  added_Xic_tree.Branch("pK_as_ppi_M", &pK_as_ppi_M, "pK_as_ppi_M/D");
  added_Xic_tree.Branch("pKpi_as_K_M", &pKpi_as_K_M, "pKpi_as_K_M/D");
  added_Xic_tree.Branch("pKpi_as_p_M", &pKpi_as_p_M, "pKpi_as_p_M/D");

  UInt_t Xic_nevents = Xic_tree->GetEntries();
  cout << "Entries in Xic tree: " << Xic_nevents << endl;

  for (UInt_t evt = 0; evt < Xic_nevents;evt++) {
    Xic_tree->GetEntry(evt);

    TVector3 dir(Xb_ENDVERTEX_X-Xb_OWNPV_X,Xb_ENDVERTEX_Y-Xb_OWNPV_Y,Xb_ENDVERTEX_Z-Xb_OWNPV_Z);
    TVector3 mom;
    mom.SetPtEtaPhi(Xb_PT,Xb_ETA,Xb_PHI);
    double dmag2 = dir.Mag2();
    double ptprime = 0;
    if ( 0 == dmag2 ) ptprime = mom.Mag();
    else ptprime = (mom - dir * ( mom.Dot( dir ) / dmag2 )).Mag() ;
    Xb_CorrM = sqrt(Xb_M*Xb_M + ptprime*ptprime) + ptprime;

    TLorentzVector Xb;
    Xb.SetPtEtaPhiM(Xb_PT,Xb_ETA,Xb_PHI,Xb_CorrM);
    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);
    for(int i = 0; i < Added_n_Particles; i++){
      TLorentzVector Hpi;
      Hpi.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],pionmass);
      TLorentzVector HK;
      HK.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],kaonmass);
      TLorentzVector Hp;
      Hp.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],protonmass);
      TLorentzVector Xcpi = Hpi + Xc;
      TLorentzVector XcK = HK + Xc;
      TLorentzVector Xcp = Hp + Xc;
      Xcpi.Boost(-Xb.BoostVector());
      Xcpi_CosTheta[i] = cos(Xcpi.Angle(Xb.Vect()));
      XcK.Boost(-Xb.BoostVector());
      XcK_CosTheta[i] = cos(XcK.Angle(Xb.Vect()));
      Xcp.Boost(-Xb.BoostVector());
      Xcp_CosTheta[i] = cos(Xcp.Angle(Xb.Vect()));
    }

    TLorentzVector proton;
    proton.SetPtEtaPhiM(p_PT,p_ETA,p_PHI,protonmass);
    TLorentzVector kaon;
    kaon.SetPtEtaPhiM(K_PT,K_ETA,K_PHI,kaonmass);
    TLorentzVector pion;
    pion.SetPtEtaPhiM(pi_PT,pi_ETA,pi_PHI,pionmass);

    p_beta  = (-proton.P()+kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    K_beta  = ( proton.P()-kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    pi_beta = ( proton.P()+kaon.P()-pion.P())/(proton.P()+kaon.P()+pion.P());

    TLorentzVector p_as_pi;
    p_as_pi.SetVectM(proton.Vect(),pionmass);
    TLorentzVector p_as_K;
    p_as_K.SetVectM(proton.Vect(),kaonmass);

    TLorentzVector K_as_pi;
    K_as_pi.SetVectM(kaon.Vect(),pionmass);
    TLorentzVector K_as_p;
    K_as_p.SetVectM(kaon.Vect(),protonmass);

    TLorentzVector pi_as_K;
    pi_as_K.SetVectM(pion.Vect(),kaonmass);
    TLorentzVector pi_as_p;
    pi_as_p.SetVectM(pion.Vect(),protonmass);

    p_as_piKpi_M = (p_as_pi + kaon + pion).M();
    p_as_KKpi_M = (p_as_K + kaon + pion).M();

    pK_as_pipi_M = (proton + K_as_pi + pion).M();
    pK_as_ppi_M = (proton + K_as_p + pion).M();

    pKpi_as_K_M = (proton + kaon + pi_as_K).M();
    pKpi_as_p_M = (proton + kaon + pi_as_p).M();

    added_Xic_tree.Fill();

  }

  Xic_tree->SetDirectory(0);
  added_Xic_tree.Write();

  fSLBS->cd();
  TTree *Xic0_tree = (TTree*)gDirectory->Get("Xib2Xic0MuNu/Xic02pKKpi/DecayTree");

  double p_P, SSK1_P, SSK2_P, pi_P;
  double SSK1_PT, SSK2_PT, SSK1_ETA, SSK2_ETA, SSK1_PHI, SSK2_PHI;

  Xic0_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)
  Xic0_tree->SetBranchStatus("Xib_M",1);
  Xic0_tree->SetBranchStatus("Xib_PT",1);
  Xic0_tree->SetBranchStatus("Xib_ETA",1);
  Xic0_tree->SetBranchStatus("Xib_PHI",1);
  Xic0_tree->SetBranchStatus("Xib_OWNPV_X",1);
  Xic0_tree->SetBranchStatus("Xib_OWNPV_Y",1);
  Xic0_tree->SetBranchStatus("Xib_OWNPV_Z",1);
  Xic0_tree->SetBranchStatus("Xib_ENDVERTEX_X",1);
  Xic0_tree->SetBranchStatus("Xib_ENDVERTEX_Y",1);
  Xic0_tree->SetBranchStatus("Xib_ENDVERTEX_Z",1);

  Xic0_tree->SetBranchStatus("Xic_M",1);
  Xic0_tree->SetBranchStatus("Xic_PT",1);
  Xic0_tree->SetBranchStatus("Xic_ETA",1);
  Xic0_tree->SetBranchStatus("Xic_PHI",1);

  Xic0_tree->SetBranchStatus("Added_n_Particles",1);
  Xic0_tree->SetBranchStatus("Added_H_PT",1);
  Xic0_tree->SetBranchStatus("Added_H_ETA",1);
  Xic0_tree->SetBranchStatus("Added_H_PHI",1);

  Xic0_tree->SetBranchStatus("p_P",1);
  Xic0_tree->SetBranchStatus("SSK1_P",1);
  Xic0_tree->SetBranchStatus("SSK2_P",1);
  Xic0_tree->SetBranchStatus("pi_P",1);

  Xic0_tree->SetBranchStatus("p_PT",1);
  Xic0_tree->SetBranchStatus("p_ETA",1);
  Xic0_tree->SetBranchStatus("p_PHI",1);
  Xic0_tree->SetBranchStatus("SSK1_PT",1);
  Xic0_tree->SetBranchStatus("SSK1_ETA",1);
  Xic0_tree->SetBranchStatus("SSK1_PHI",1);
  Xic0_tree->SetBranchStatus("SSK2_PT",1);
  Xic0_tree->SetBranchStatus("SSK2_ETA",1);
  Xic0_tree->SetBranchStatus("SSK2_PHI",1);
  Xic0_tree->SetBranchStatus("pi_PT",1);
  Xic0_tree->SetBranchStatus("pi_ETA",1);
  Xic0_tree->SetBranchStatus("pi_PHI",1);

  //set the branch addresses
  Xic0_tree->SetBranchAddress("Xib_M",&Xb_M);
  Xic0_tree->SetBranchAddress("Xib_PT",&Xb_PT);
  Xic0_tree->SetBranchAddress("Xib_ETA",&Xb_ETA);
  Xic0_tree->SetBranchAddress("Xib_PHI",&Xb_PHI);
  Xic0_tree->SetBranchAddress("Xib_OWNPV_X",&Xb_OWNPV_X);
  Xic0_tree->SetBranchAddress("Xib_OWNPV_Y",&Xb_OWNPV_Y);
  Xic0_tree->SetBranchAddress("Xib_OWNPV_Z",&Xb_OWNPV_Z);
  Xic0_tree->SetBranchAddress("Xib_ENDVERTEX_X",&Xb_ENDVERTEX_X);
  Xic0_tree->SetBranchAddress("Xib_ENDVERTEX_Y",&Xb_ENDVERTEX_Y);
  Xic0_tree->SetBranchAddress("Xib_ENDVERTEX_Z",&Xb_ENDVERTEX_Z);

  Xic0_tree->SetBranchAddress("Xic_M",&Xc_M);
  Xic0_tree->SetBranchAddress("Xic_PT",&Xc_PT);
  Xic0_tree->SetBranchAddress("Xic_ETA",&Xc_ETA);
  Xic0_tree->SetBranchAddress("Xic_PHI",&Xc_PHI);

  Xic0_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Xic0_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Xic0_tree->SetBranchAddress("Added_H_ETA",&Added_H_ETA);
  Xic0_tree->SetBranchAddress("Added_H_PHI",&Added_H_PHI);

  Xic0_tree->SetBranchAddress("p_P",&p_P);
  Xic0_tree->SetBranchAddress("SSK1_P",&SSK1_P);
  Xic0_tree->SetBranchAddress("SSK2_P",&SSK2_P);
  Xic0_tree->SetBranchAddress("pi_P",&pi_P);

  Xic0_tree->SetBranchAddress("p_PT",&p_PT);
  Xic0_tree->SetBranchAddress("SSK1_PT",&SSK1_PT);
  Xic0_tree->SetBranchAddress("SSK2_PT",&SSK2_PT);
  Xic0_tree->SetBranchAddress("pi_PT",&pi_PT);
  Xic0_tree->SetBranchAddress("p_ETA",&p_ETA);
  Xic0_tree->SetBranchAddress("SSK1_ETA",&SSK1_ETA);
  Xic0_tree->SetBranchAddress("SSK2_ETA",&SSK2_ETA);
  Xic0_tree->SetBranchAddress("pi_ETA",&pi_ETA);
  Xic0_tree->SetBranchAddress("p_PHI",&p_PHI);
  Xic0_tree->SetBranchAddress("SSK1_PHI",&SSK1_PHI);
  Xic0_tree->SetBranchAddress("SSK2_PHI",&SSK2_PHI);
  Xic0_tree->SetBranchAddress("pi_PHI",&pi_PHI);


  double SSK1_beta, SSK2_beta;

  f1->cd();
  //f1->mkdir("Xib2Xic0MuNu/Xic02pKKpi");
  //f1->cd("Xib2Xic0MuNu/Xic02pKKpi");
  TTree added_Xic0_tree("Xic02pKKpi","Xic02pKKpi");

  added_Xic0_tree.Branch("Xib_CorrM", &Xb_CorrM, "Xib_CorrM/D");
  added_Xic0_tree.Branch("p_beta", &p_beta, "p_beta/D");
  added_Xic0_tree.Branch("SSK1_beta", &SSK1_beta, "SSK1_beta/D");
  added_Xic0_tree.Branch("SSK2_beta", &SSK2_beta, "SSK2_beta/D");
  added_Xic0_tree.Branch("pi_beta", &pi_beta, "pi_beta/D");
  added_Xic0_tree.Branch("Added_n_Particles", &Added_n_Particles, "Added_n_Particles/I");
  added_Xic0_tree.Branch("Xcpi_CosTheta", &Xcpi_CosTheta, "Xcpi_CosTheta[Added_n_Particles]/F");
  added_Xic0_tree.Branch("XcK_CosTheta", &XcK_CosTheta, "XcK_CosTheta[Added_n_Particles]/F");
  added_Xic0_tree.Branch("Xcp_CosTheta", &Xcp_CosTheta, "Xcp_CosTheta[Added_n_Particles]/F");
  added_Xic0_tree.Branch("p_as_piKKpi_M", &p_as_piKpi_M, "p_as_piKKpi_M/D");
  added_Xic0_tree.Branch("p_as_KKKpi_M", &p_as_KKpi_M, "p_as_KKKpi_M/D");

  UInt_t Xic0_nevents = Xic0_tree->GetEntries();
  cout << "Entries in Xic0 tree: " << Xic0_nevents << endl;

  for (UInt_t evt = 0; evt < Xic0_nevents;evt++) {
    Xic0_tree->GetEntry(evt);

    TVector3 dir(Xb_ENDVERTEX_X-Xb_OWNPV_X,Xb_ENDVERTEX_Y-Xb_OWNPV_Y,Xb_ENDVERTEX_Z-Xb_OWNPV_Z);
    TVector3 mom;
    mom.SetPtEtaPhi(Xb_PT,Xb_ETA,Xb_PHI);
    double dmag2 = dir.Mag2();
    double ptprime = 0;
    if ( 0 == dmag2 ) ptprime = mom.Mag();
    else ptprime = (mom - dir * ( mom.Dot( dir ) / dmag2 )).Mag() ;
    Xb_CorrM = sqrt(Xb_M*Xb_M + ptprime*ptprime) + ptprime;

    TLorentzVector Xb;
    Xb.SetPtEtaPhiM(Xb_PT,Xb_ETA,Xb_PHI,Xb_CorrM);
    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);
    for(int i = 0; i < Added_n_Particles; i++){
      TLorentzVector Hpi;
      Hpi.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],pionmass);
      TLorentzVector HK;
      HK.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],kaonmass);
      TLorentzVector Hp;
      Hp.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],protonmass);
      TLorentzVector Xcpi = Hpi + Xc;
      TLorentzVector XcK = HK + Xc;
      TLorentzVector Xcp = Hp + Xc;
      Xcpi.Boost(-Xb.BoostVector());
      Xcpi_CosTheta[i] = cos(Xcpi.Angle(Xb.Vect()));
      XcK.Boost(-Xb.BoostVector());
      XcK_CosTheta[i] = cos(XcK.Angle(Xb.Vect()));
      Xcp.Boost(-Xb.BoostVector());
      Xcp_CosTheta[i] = cos(Xcp.Angle(Xb.Vect()));
    }
    p_beta    = (-p_P+SSK1_P+SSK2_P+pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);
    SSK1_beta = ( p_P-SSK1_P+SSK2_P+pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);
    SSK2_beta = ( p_P+SSK1_P-SSK2_P+pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);
    pi_beta   = ( p_P+SSK1_P+SSK2_P-pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);

    TLorentzVector proton;
    proton.SetPtEtaPhiM(p_PT,p_ETA,p_PHI,protonmass);
    TLorentzVector kaon1;
    kaon1.SetPtEtaPhiM(SSK1_PT,SSK1_ETA,SSK1_PHI,kaonmass);
    TLorentzVector kaon2;
    kaon2.SetPtEtaPhiM(SSK2_PT,SSK2_ETA,SSK2_PHI,kaonmass);
    TLorentzVector pion;
    pion.SetPtEtaPhiM(pi_PT,pi_ETA,pi_PHI,pionmass);

    TLorentzVector p_as_pi;
    p_as_pi.SetVectM(proton.Vect(),pionmass);
    TLorentzVector p_as_K;
    p_as_K.SetVectM(proton.Vect(),kaonmass);

    p_as_piKpi_M = (p_as_pi + kaon1 + kaon2 + pion).M();
    p_as_KKpi_M = (p_as_K + kaon1 + kaon2 + pion).M();

    added_Xic0_tree.Fill();

  }

  added_Xic0_tree.Write();

  fSLBS->cd();
  TTree *Omegac_tree = (TTree*)gDirectory->Get("Omegab2Omegac0MuNu/Omegac2pKKpi/DecayTree");

  Omegac_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)
  Omegac_tree->SetBranchStatus("Omegab_M",1);
  Omegac_tree->SetBranchStatus("Omegab_PT",1);
  Omegac_tree->SetBranchStatus("Omegab_ETA",1);
  Omegac_tree->SetBranchStatus("Omegab_PHI",1);
  Omegac_tree->SetBranchStatus("Omegab_OWNPV_X",1);
  Omegac_tree->SetBranchStatus("Omegab_OWNPV_Y",1);
  Omegac_tree->SetBranchStatus("Omegab_OWNPV_Z",1);
  Omegac_tree->SetBranchStatus("Omegab_ENDVERTEX_X",1);
  Omegac_tree->SetBranchStatus("Omegab_ENDVERTEX_Y",1);
  Omegac_tree->SetBranchStatus("Omegab_ENDVERTEX_Z",1);

  Omegac_tree->SetBranchStatus("Omegac_M",1);
  Omegac_tree->SetBranchStatus("Omegac_PT",1);
  Omegac_tree->SetBranchStatus("Omegac_ETA",1);
  Omegac_tree->SetBranchStatus("Omegac_PHI",1);

  Omegac_tree->SetBranchStatus("Added_n_Particles",1);
  Omegac_tree->SetBranchStatus("Added_H_PT",1);
  Omegac_tree->SetBranchStatus("Added_H_ETA",1);
  Omegac_tree->SetBranchStatus("Added_H_PHI",1);

  Omegac_tree->SetBranchStatus("p_P",1);
  Omegac_tree->SetBranchStatus("SSK1_P",1);
  Omegac_tree->SetBranchStatus("SSK2_P",1);
  Omegac_tree->SetBranchStatus("pi_P",1);

  //set the branch addresses
  Omegac_tree->SetBranchAddress("Omegab_M",&Xb_M);
  Omegac_tree->SetBranchAddress("Omegab_PT",&Xb_PT);
  Omegac_tree->SetBranchAddress("Omegab_ETA",&Xb_ETA);
  Omegac_tree->SetBranchAddress("Omegab_PHI",&Xb_PHI);
  Omegac_tree->SetBranchAddress("Omegab_OWNPV_X",&Xb_OWNPV_X);
  Omegac_tree->SetBranchAddress("Omegab_OWNPV_Y",&Xb_OWNPV_Y);
  Omegac_tree->SetBranchAddress("Omegab_OWNPV_Z",&Xb_OWNPV_Z);
  Omegac_tree->SetBranchAddress("Omegab_ENDVERTEX_X",&Xb_ENDVERTEX_X);
  Omegac_tree->SetBranchAddress("Omegab_ENDVERTEX_Y",&Xb_ENDVERTEX_Y);
  Omegac_tree->SetBranchAddress("Omegab_ENDVERTEX_Z",&Xb_ENDVERTEX_Z);

  Omegac_tree->SetBranchAddress("Omegac_M",&Xc_M);
  Omegac_tree->SetBranchAddress("Omegac_PT",&Xc_PT);
  Omegac_tree->SetBranchAddress("Omegac_ETA",&Xc_ETA);
  Omegac_tree->SetBranchAddress("Omegac_PHI",&Xc_PHI);

  Omegac_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Omegac_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Omegac_tree->SetBranchAddress("Added_H_ETA",&Added_H_ETA);
  Omegac_tree->SetBranchAddress("Added_H_PHI",&Added_H_PHI);

  Omegac_tree->SetBranchAddress("p_P",&p_P);
  Omegac_tree->SetBranchAddress("SSK1_P",&SSK1_P);
  Omegac_tree->SetBranchAddress("SSK2_P",&SSK2_P);
  Omegac_tree->SetBranchAddress("pi_P",&pi_P);

  f1->cd();
  //f1->mkdir("Omegab2Omegac0MuNu/Omegac2pKKpi");
  //f1->cd("Omegab2Omegac0MuNu/Omegac2pKKpi");
  TTree added_Omegac_tree("Omegac2pKKpi","Omegac2pKKpi");

  added_Omegac_tree.Branch("Omegab_CorrM", &Xb_CorrM, "Omegab_CorrM/D");
  added_Omegac_tree.Branch("p_beta", &p_beta, "p_beta/D");
  added_Omegac_tree.Branch("SSK1_beta", &SSK1_beta, "SSK1_beta/D");
  added_Omegac_tree.Branch("SSK2_beta", &SSK2_beta, "SSK2_beta/D");
  added_Omegac_tree.Branch("pi_beta", &pi_beta, "pi_beta/D");
  added_Omegac_tree.Branch("Added_n_Particles", &Added_n_Particles, "Added_n_Particles/I");
  added_Omegac_tree.Branch("Xcpi_CosTheta", &Xcpi_CosTheta, "Xcpi_CosTheta[Added_n_Particles]/F");
  added_Omegac_tree.Branch("XcK_CosTheta", &XcK_CosTheta, "XcK_CosTheta[Added_n_Particles]/F");
  added_Omegac_tree.Branch("Xcp_CosTheta", &Xcp_CosTheta, "Xcp_CosTheta[Added_n_Particles]/F");

  UInt_t Omegac_nevents = Omegac_tree->GetEntries();
  cout << "Entries in Omegac tree: " << Omegac_nevents << endl;

  for (UInt_t evt = 0; evt < Omegac_nevents;evt++) {
    Omegac_tree->GetEntry(evt);

    TVector3 dir(Xb_ENDVERTEX_X-Xb_OWNPV_X,Xb_ENDVERTEX_Y-Xb_OWNPV_Y,Xb_ENDVERTEX_Z-Xb_OWNPV_Z);
    TVector3 mom;
    mom.SetPtEtaPhi(Xb_PT,Xb_ETA,Xb_PHI);
    double dmag2 = dir.Mag2();
    double ptprime = 0;
    if ( 0 == dmag2 ) ptprime = mom.Mag();
    else ptprime = (mom - dir * ( mom.Dot( dir ) / dmag2 )).Mag() ;
    Xb_CorrM = sqrt(Xb_M*Xb_M + ptprime*ptprime) + ptprime;

    TLorentzVector Xb;
    Xb.SetPtEtaPhiM(Xb_PT,Xb_ETA,Xb_PHI,Xb_CorrM);
    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);
    for(int i = 0; i < Added_n_Particles; i++){
      TLorentzVector Hpi;
      Hpi.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],pionmass);
      TLorentzVector HK;
      HK.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],kaonmass);
      TLorentzVector Hp;
      Hp.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],protonmass);
      TLorentzVector Xcpi = Hpi + Xc;
      TLorentzVector XcK = HK + Xc;
      TLorentzVector Xcp = Hp + Xc;
      Xcpi.Boost(-Xb.BoostVector());
      Xcpi_CosTheta[i] = cos(Xcpi.Angle(Xb.Vect()));
      XcK.Boost(-Xb.BoostVector());
      XcK_CosTheta[i] = cos(XcK.Angle(Xb.Vect()));
      Xcp.Boost(-Xb.BoostVector());
      Xcp_CosTheta[i] = cos(Xcp.Angle(Xb.Vect()));
    }
    p_beta    = (-p_P+SSK1_P+SSK2_P+pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);
    SSK1_beta = ( p_P-SSK1_P+SSK2_P+pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);
    SSK2_beta = ( p_P+SSK1_P-SSK2_P+pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);
    pi_beta   = ( p_P+SSK1_P+SSK2_P-pi_P)/(p_P+SSK1_P+SSK2_P+pi_P);

    added_Omegac_tree.Fill();

  }

  added_Omegac_tree.Write();

  clock->Stop();clock->Print();delete clock;
  return;
}
