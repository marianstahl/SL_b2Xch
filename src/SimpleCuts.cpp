/*

    2015-05-25 M. Stahl

    Determine cut and background efficiencies for a given cut string

*/

#include <iostream>
#include <algorithm>
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
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TProfile.h"

#include "TArrow.h"
#include "TLatex.h"

#include "RooRealVar.h"
#include "RooDataHist.h"

#include "../include/Fit_Charm_Baryon.h"
#include "../include/configuration.h"
#include "../include/MyStyle.h"

using namespace std;
using namespace RooFit;

TString temp;

struct cH_mult_cand{
  double mass;
  float ip;
  int index;
};

bool compareByIP(const cH_mult_cand &a, const cH_mult_cand &b)
{
  return a.ip < b.ip;
}

template <class cH> void fill_profile(vector<cH> multiple_cH_candidate, TProfile *Xc_IP);
void SimpleCuts();
template<class T> void make_1D_Plot(T *hist, configuration* myconfig);
void make_2D_Plot(TH2D* hist, configuration* myconfig);
void make_overlay_Plot(TH1D *SS_hist, TH1D *OS_hist, configuration *myconfig);

int main(int argc, char **argv)
{
  SimpleCuts();
  //return 0;
}

void SimpleCuts(){

  //Xic0   --> //p_ProbNNp > 0.1 && SSK1_ProbNNk > 0.25 && SSK2_ProbNNk > 0.25 && pi_ProbNNpi > 0.05 && !(MyFriend.p_as_KKKpi_M > 1855 && MyFriend.p_as_KKKpi_M < 1875)  && Xic_M > 2460 && Xic_M < 2485 && MyFriend.Xib_CorrM > 4500 && MyFriend.Xib_CorrM < 6500 && Added_H_ProbNNpi > 0.4 && Added_H_PT > 200 && ((Added_H_PROBNNPID < 0 && Xic_ID < 0) || (Added_H_PROBNNPID > 0 && Xic_ID > 0)) && Added_CharmH_M < 4000 && TMath::Abs(Added_H_PROBNNPID) == 221 && Added_CharmH_VERTEXCHI2_NEW < 3
  //Omegac --> //p_ProbNNp > 0.3 && pi_ProbNNpi > 0.0 && SSK1_ProbNNk > 0.4 && SSK2_ProbNNk > 0.4 && MyFriend.Omegab_CorrM > 5500 && MyFriend.Omegab_CorrM < 6800 && Omegac_M > 2685 && Omegac_M < 2710 && Added_H_ProbNNpi > 0.4 && Added_H_PT > 200 && ((Added_H_PROBNNPID > 0 && Omegac_ID < 0) || (Added_H_PROBNNPID < 0 && Omegac_ID > 0)) && TMath::Abs(Added_H_PROBNNPID) == 221 && Added_CharmH_VERTEXCHI2_NEW < 3
  //Xic    --> //!(MyFriend.p_as_piKpi_M > 1860 && MyFriend.p_as_piKpi_M < 1880) && !(MyFriend.p_as_piKpi_M > 2005 && MyFriend.p_as_piKpi_M < 2025) && !(MyFriend.p_as_KKpi_M > 1860 && MyFriend.p_as_KKpi_M < 1880) && !(MyFriend.p_as_KKpi_M > 1955 && MyFriend.p_as_KKpi_M < 1985) && p_ProbNNp > 0.28 && K_ProbNNk > 0.15 && pi_ProbNNpi > 0.2 && MyFriend.Xib_CorrM > 4500 && MyFriend.Xib_CorrM < 6500 && p_CosTheta > -0.8 && Added_H_ProbNNpi > 0.4 && Added_H_PT > 200 && ((Added_H_PROBNNPID < 0 && Xic_ID < 0) || (Added_H_PROBNNPID > 0 && Xic_ID > 0)) && Added_CharmH_M < 4000 && TMath::Abs(Added_H_PROBNNPID) == 221 && Added_CharmH_VERTEXCHI2_NEW < 3
  //Xic    -?-> // !(MyFriend.p_as_piKpi_M > 1860 && MyFriend.p_as_piKpi_M < 1880) && !(MyFriend.p_as_piKpi_M > 2005 && MyFriend.p_as_piKpi_M < 2025) && !(MyFriend.p_as_KKpi_M > 1860 && MyFriend.p_as_KKpi_M < 1880) && !(MyFriend.p_as_KKpi_M > 1955 && MyFriend.p_as_KKpi_M < 1985) && p_ProbNNp > 0.28 && K_ProbNNk > 0.15 && pi_ProbNNpi > 0.2 && MyFriend.Xib_CorrM > 4700 && MyFriend.Xib_CorrM < 6300 && ((sqrt(Xic_Dalitz_Kminus_piplus_M2) > 842 && sqrt(Xic_Dalitz_Kminus_piplus_M2) < 942) || (sqrt(Xic_Dalitz_Kminus_pplus_M2) > 1512 && sqrt(Xic_Dalitz_Kminus_pplus_M2) < 1526)) && p_CosTheta > -0.7 && Xic_M > 2460 && Xic_M < 2485 && Added_H_ProbNNk > 0.1 && Added_H_PT > 150 && ((Added_H_PROBNNPID < 0 && Xic_ID > 0) || (Added_H_PROBNNPID > 0 && Xic_ID < 0)) && Added_CharmH_M_kaon < 4000 && TMath::Abs(Added_H_PROBNNPID) == 321 && Added_CharmH_VERTEXCHI2_NEW < 4

  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  bool Xic_sel = true;
  bool Xic0_sel = false;
  bool Omegac_sel = true;
  bool make_IP_profile_plots = false;
  bool clean_vertex = true;
  bool purge_vertex = true;

  const double xicmass = 2469.26; //MeV --> measured by fit
  const double xic0mass = 2472.18; //MeV --> measured by fit

  MyStyle();

  configuration* myconfig = new configuration();
  myconfig->set_version(1);
  myconfig->fill_cs("Basic selection");
  myconfig->set_current_cs("NoXbCut");

  TH1D::AddDirectory(0);
  TProfile::AddDirectory(0);
  TH2D::AddDirectory(0);

  TH1D *Xic_hist = new TH1D("Xic",";;",myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Xic_beta_hist = new TH2D("Xic_beta",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Xic_beta_cut_hist = new TH2D("Xic_beta_cut",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Xic_2D_dump = new TH2D("Xic_p_CosTheta",";M_{inv}(pK#pi) (GeV);#Theta^{CM}_{p}",2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi(),100,-1,1);

  int n_CorrM_bins = 84;
  double CorrM_lo = 3000, CorrM_hi = 7200;
  TH1D *Xib0_hist = new TH1D("Xib0",";M_{corr}(#Xi^{+}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);
  TH1D *Xib_hist = new TH1D("Xib",";M_{corr}(#Xi^{0}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);
  TH1D *Omegab_hist = new TH1D("Omegab",";M_{corr}(#Omega^{0}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);
  TH1D *Xib0_SB_hist = new TH1D("Xib0_SB",";M_{corr}(#Xi^{+}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);
  TH1D *Xib_SB_hist = new TH1D("Xib_SB",";M_{corr}(#Xi^{0}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);
  TH1D *Omegab_SB_hist = new TH1D("Omegab_SB",";M_{corr}(#Omega^{0}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);

  int n_profilebins = 20;
  TProfile *Xic_IP = new TProfile("Xic_IP",";IP (mm);<nTracks>",n_profilebins,0,1," ");
  TProfile *Xic0_IP = new TProfile("Xic0_IP",";IP (mm);<nTracks>",n_profilebins,0,1," ");
  TProfile *Omegac_IP = new TProfile("Omegac_IP",";IP (mm);<nTracks>",n_profilebins,0,1," ");

  double Xicpilo = 2600, Xicpihi = 3800; int nbinsXicpi = 200;
  double XicKlo = 2960, XicKhi = 3560; int nbinsXicK = 150;

  temp.Form(";M_{inv}(#Xi^{#pm}_{c}#pi)-M_{inv}(pK#pi)+M_{PDG}(#Xi^{#pm}_{c}) (GeV); Events/%g MeV",(Xicpihi-Xicpilo)/nbinsXicpi);
  TH1D *XicpiSS_hist = new TH1D("XicpiSS_hist",temp,nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *XicpiOS_hist = new TH1D("Xicpi_hist",temp,nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *XicpiSS_IPfail = new TH1D("XicpiSS_IPfail",temp,nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *XicpiOS_IPfail = new TH1D("Xicpi_IPfail",temp,nbinsXicpi,Xicpilo,Xicpihi);
  temp.Form(";M_{inv}(#Xi^{#pm}_{c}K)-M_{inv}(pK#pi)+M_{PDG}(#Xi^{#pm}_{c}) (GeV); Events/%g MeV",(XicKhi-XicKlo)/nbinsXicK);
  TH1D *XicKSS_hist = new TH1D("XicKSS_hist",temp,nbinsXicK,XicKlo,XicKhi);
  TH1D *XicKOS_hist = new TH1D("XicK_hist",temp,nbinsXicK,XicKlo,XicKhi);
  TH1D *XicKSS_IPfail = new TH1D("XicKSS_IPfail",temp,nbinsXicK,XicKlo,XicKhi);
  TH1D *XicKOS_IPfail = new TH1D("XicK_IPfail",temp,nbinsXicK,XicKlo,XicKhi);
  vector<TH1D*> Xic_hists;Xic_hists.push_back(XicpiSS_hist);Xic_hists.push_back(XicpiOS_hist);Xic_hists.push_back(XicKSS_hist);Xic_hists.push_back(XicKOS_hist);
  vector<TH1D*> Xic_IPfails;Xic_IPfails.push_back(XicpiSS_IPfail);Xic_IPfails.push_back(XicpiOS_IPfail);Xic_IPfails.push_back(XicKSS_IPfail);Xic_IPfails.push_back(XicKOS_IPfail);

  myconfig->set_particle("Xic0");
  TH1D *Xic0_hist = new TH1D("Xic0",";;",myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Xic0_beta_hist = new TH2D("Xic0_beta",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Xic0_beta_cut_hist = new TH2D("Xic0_beta_cut",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Xic0_2D_dump = new TH2D("Xic0_p_CosTheta",";M_{inv}(pK#pi) (GeV);#Theta^{CM}_{p}",2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi(),100,-1,1);

  temp.Form(";M_{inv}(#Xi^{0}_{c}#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Xi^{0}_{c}) (GeV); Events/%g MeV",(Xicpihi-Xicpilo)/nbinsXicpi);
  TH1D *Xic0piSS_hist = new TH1D("Xic0piSS_hist",temp,nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *Xic0piOS_hist = new TH1D("Xic0pi_hist",temp,nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *Xic0piSS_IPfail = new TH1D("Xic0piSS_IPfail",temp,nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *Xic0piOS_IPfail = new TH1D("Xic0pi_IPfail",temp,nbinsXicpi,Xicpilo,Xicpihi);
  temp.Form(";M_{inv}(#Xi^{0}_{c}K)-M_{inv}(pKK#pi)+M_{PDG}(#Xi^{0}_{c}) (GeV); Events/%g MeV",(XicKhi-XicKlo)/nbinsXicK);
  TH1D *Xic0KSS_hist = new TH1D("Xic0KSS_hist",temp,nbinsXicK,XicKlo,XicKhi);
  TH1D *Xic0KOS_hist = new TH1D("Xic0K_hist",temp,nbinsXicK,XicKlo,XicKhi);
  TH1D *Xic0KSS_IPfail = new TH1D("Xic0KSS_IPfail",temp,nbinsXicK,XicKlo,XicKhi);
  TH1D *Xic0KOS_IPfail = new TH1D("Xic0K_IPfail",temp,nbinsXicK,XicKlo,XicKhi);
  vector<TH1D*> Xic0_hists;Xic0_hists.push_back(Xic0piSS_hist);Xic0_hists.push_back(Xic0piOS_hist);Xic0_hists.push_back(Xic0KSS_hist);Xic0_hists.push_back(Xic0KOS_hist);
  vector<TH1D*> Xic0_IPfails;Xic0_IPfails.push_back(Xic0piSS_IPfail);Xic0_IPfails.push_back(Xic0piOS_IPfail);Xic0_IPfails.push_back(Xic0KSS_IPfail);Xic0_IPfails.push_back(Xic0KOS_IPfail);
  TH1D *Xic0pipiSS_hist = new TH1D("Xic0pipiSS_hist",";M_{inv}(#Xi^{0}_{c}#pi#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Xi^{0}_{c}) (GeV); Events/10 MeV",nbinsXicpi,Xicpilo,Xicpihi);
  TH1D *Xic0pipiOS_hist = new TH1D("Xic0pipi_hist",";M_{inv}(#Xi^{0}_{c}#pi#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Xi^{0}_{c}) (GeV); Events/10 MeV",nbinsXicpi,Xicpilo,Xicpihi);
  vector<TH1D*> Xic0pipi_hists;Xic0pipi_hists.push_back(Xic0pipiSS_hist);Xic0pipi_hists.push_back(Xic0pipiOS_hist);

  myconfig->set_particle("Omegac");
  TH1D *Omegac_hist = new TH1D("Omegac",";;",myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Omegac_beta_hist = new TH2D("Omegac_beta",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Omegac_beta_cut_hist = new TH2D("Omegac_beta_cut",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Omegac_2D_dump = new TH2D("Omegac_p_CosTheta",";M_{inv}(pK#pi) (GeV);#Theta^{CM}_{p}",2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi(),100,-1,1);

  double Omegacpilo = 2840, Omegacpihi = 4040; int nbinsOmegacpi = 120;
  double OmegacKlo = 3190, OmegacKhi = 3990; int nbinsOmegacK = 80;
  TH1D *OmegacpiSS_hist = new TH1D("OmegacpiSS_hist",";M_{inv}(#Omega^{0}_{c}#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacpi,Omegacpilo,Omegacpihi);
  TH1D *OmegacpiOS_hist = new TH1D("Omegacpi_hist",";M_{inv}(#Omega^{0}_{c}#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacpi,Omegacpilo,Omegacpihi);
  TH1D *OmegacKSS_hist = new TH1D("OmegacKSS_hist",";M_{inv}(#Omega^{0}_{c}K)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacK,OmegacKlo,OmegacKhi);
  TH1D *OmegacKOS_hist = new TH1D("OmegacK_hist",";M_{inv}(#Omega^{0}_{c}K)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacK,OmegacKlo,OmegacKhi);
  vector<TH1D*> Omegac_hists;Omegac_hists.push_back(OmegacpiSS_hist);Omegac_hists.push_back(OmegacpiOS_hist);Omegac_hists.push_back(OmegacKSS_hist);Omegac_hists.push_back(OmegacKOS_hist);

  TH1D *OmegacpiSS_IPfail = new TH1D("OmegacpiSS_IPfail",";M_{inv}(#Omega^{0}_{c}#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacpi,Omegacpilo,Omegacpihi);
  TH1D *OmegacpiOS_IPfail = new TH1D("Omegacpi_IPfail",";M_{inv}(#Omega^{0}_{c}#pi)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacpi,Omegacpilo,Omegacpihi);
  TH1D *OmegacKSS_IPfail = new TH1D("OmegacKSS_IPfail",";M_{inv}(#Omega^{0}_{c}K)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacK,OmegacKlo,OmegacKhi);
  TH1D *OmegacKOS_IPfail = new TH1D("OmegacK_IPfail",";M_{inv}(#Omega^{0}_{c}K)-M_{inv}(pKK#pi)+M_{PDG}(#Omega^{0}_{c}) (GeV); Events/10 MeV",nbinsOmegacK,OmegacKlo,OmegacKhi);
  vector<TH1D*> Omegac_IPfails;Omegac_IPfails.push_back(OmegacpiSS_IPfail);Omegac_IPfails.push_back(OmegacpiOS_IPfail);Omegac_IPfails.push_back(OmegacKSS_IPfail);Omegac_IPfails.push_back(OmegacKOS_IPfail);

  double p_ProbNNp, pi_ProbNNpi, K_ProbNNk, K2_ProbNNk;
  double Xb_M, Xc_M;
  double Xc_PT, Xc_ETA, Xc_PHI;
  double MisID_D02KKKpi_M, MisID_D2piKpi_M, MisID_Ds2KKpi_M;
  double KStar_M2, Lambda1520_M2;
  double p_CosTheta, p_beta;
  int Added_n_Particles, Added_n_Tracks, Xc_ID;
  float Added_CharmH_M[200], Added_CharmH_M_kaon[200], Added_CharmH_PT[200];
  float Added_H_PT[200], Added_H_ETA[200], Added_H_PHI[200];
  float Added_CharmH_VERTEXCHI2_NEW[200], Added_CharmH_IP_NEW[200];
  float Added_H_PROBNNPID[200], Added_H_ProbNNpi[200], Added_H_ProbNNk[200];

  temp = myconfig->get_tupledir() + "/SLBaryonSpectroscopyStrp21.root";
  gErrorIgnoreLevel = kError;
  TFile *fSLBS = new TFile(temp,"read");
  TTree *Xic_tree = (TTree*)gDirectory->Get("Xib02XicMuNu/Xic2pKpi/DecayTree");
  TTree *Xic0_tree = (TTree*)gDirectory->Get("Xib2Xic0MuNu/Xic02pKKpi/DecayTree");
  TTree *Omegac_tree = (TTree*)gDirectory->Get("Omegab2Omegac0MuNu/Omegac2pKKpi/DecayTree");
  temp = myconfig->get_tupledir() + "SLBaryonSpectroscopyStrp21_friend.root";
  Xic_tree->AddFriend("Xic2pKpi",temp);
  Xic0_tree->AddFriend("Xic02pKKpi",temp);
  Omegac_tree->AddFriend("Omegac2pKKpi",temp);
  gErrorIgnoreLevel = kPrint;

  myconfig->set_particle("Xic");
  Xic_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)
  Xic_tree->SetBranchStatus("p_ProbNNp",1);
  Xic_tree->SetBranchStatus("pi_ProbNNpi",1);
  Xic_tree->SetBranchStatus("K_ProbNNk",1);
  Xic_tree->SetBranchStatus("Xic2pKpi.Xib_CorrM",1);
  Xic_tree->SetBranchStatus("Xic_M",1);
  Xic_tree->SetBranchStatus("Xic_ID",1);
  Xic_tree->SetBranchStatus("Xic2pKpi.p_as_piKpi_M",1);
  Xic_tree->SetBranchStatus("Xic2pKpi.p_as_KKpi_M",1);
  Xic_tree->SetBranchStatus("Xic_Dalitz_Kminus_piplus_M2",1);
  Xic_tree->SetBranchStatus("Xic_Dalitz_Kminus_pplus_M2",1);
  Xic_tree->SetBranchStatus("p_CosTheta",1);
  Xic_tree->SetBranchStatus("Xic2pKpi.p_beta",1);
  Xic_tree->SetBranchStatus("Added_n_Particles",1);
  Xic_tree->SetBranchStatus("Added_n_Tracks",1);
  Xic_tree->SetBranchStatus("Added_CharmH_M",1);
  Xic_tree->SetBranchStatus("Added_CharmH_PT",1);
  Xic_tree->SetBranchStatus("Added_CharmH_M_kaon",1);
  Xic_tree->SetBranchStatus("Added_CharmH_VERTEXCHI2_NEW",1);
  Xic_tree->SetBranchStatus("Added_CharmH_IP_NEW",1);
  Xic_tree->SetBranchStatus("Added_H_PT",1);
  Xic_tree->SetBranchStatus("Added_H_PROBNNPID",1);
  Xic_tree->SetBranchStatus("Added_H_ProbNNpi",1);
  Xic_tree->SetBranchStatus("Added_H_ProbNNk",1);

  Xic_tree->SetBranchAddress("p_ProbNNp",&p_ProbNNp);
  Xic_tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);
  Xic_tree->SetBranchAddress("K_ProbNNk",&K_ProbNNk);
  Xic_tree->SetBranchAddress("Xic2pKpi.Xib_CorrM",&Xb_M);
  Xic_tree->SetBranchAddress("Xic_M",&Xc_M);
  Xic_tree->SetBranchAddress("Xic_ID",&Xc_ID);
  Xic_tree->SetBranchAddress("Xic2pKpi.p_as_piKpi_M",&MisID_D2piKpi_M);
  Xic_tree->SetBranchAddress("Xic2pKpi.p_as_KKpi_M",&MisID_Ds2KKpi_M);
  Xic_tree->SetBranchAddress("Xic_Dalitz_Kminus_piplus_M2",&KStar_M2);
  Xic_tree->SetBranchAddress("Xic_Dalitz_Kminus_pplus_M2",&Lambda1520_M2);
  Xic_tree->SetBranchAddress("p_CosTheta",&p_CosTheta);
  Xic_tree->SetBranchAddress("Xic2pKpi.p_beta",&p_beta);
  Xic_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Xic_tree->SetBranchAddress("Added_n_Tracks",&Added_n_Tracks);
  Xic_tree->SetBranchAddress("Added_CharmH_M",&Added_CharmH_M);
  Xic_tree->SetBranchAddress("Added_CharmH_PT",&Added_CharmH_PT);
  Xic_tree->SetBranchAddress("Added_CharmH_M_kaon",&Added_CharmH_M_kaon);
  Xic_tree->SetBranchAddress("Added_CharmH_VERTEXCHI2_NEW",&Added_CharmH_VERTEXCHI2_NEW);
  Xic_tree->SetBranchAddress("Added_CharmH_IP_NEW",&Added_CharmH_IP_NEW);
  Xic_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Xic_tree->SetBranchAddress("Added_H_PROBNNPID",&Added_H_PROBNNPID);
  Xic_tree->SetBranchAddress("Added_H_ProbNNpi",&Added_H_ProbNNpi);
  Xic_tree->SetBranchAddress("Added_H_ProbNNk",&Added_H_ProbNNk);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory
  RooRealVar Charm_mass("Xc_M","M_{inv}(pK^{-}pi^{+} + c.c.)",myconfig->get_XcMlo(),myconfig->get_XcMhi(),"MeV") ;

  if(Xic_sel){
    UInt_t nXicevents = Xic_tree->GetEntries();
    cout << "Entries in Xic tree: " << nXicevents << endl;

    for (UInt_t evt = 0; evt < nXicevents;evt++) {
      Xic_tree->GetEntry(evt);

      bool vetos = !(1860 < MisID_D2piKpi_M && MisID_D2piKpi_M < 1880) && !(2005 < MisID_D2piKpi_M && MisID_D2piKpi_M < 2025) && !(1860 < MisID_Ds2KKpi_M && MisID_Ds2KKpi_M < 1880) && !(1955 < MisID_Ds2KKpi_M && MisID_Ds2KKpi_M < 1985);
      bool PID_cuts = p_ProbNNp > 0.28 && K_ProbNNk > 0.15 && pi_ProbNNpi > 0.2;
      bool Dalitz_region = (842 < sqrt(KStar_M2) && sqrt(KStar_M2) < 942) || (1505 < sqrt(Lambda1520_M2) && sqrt(Lambda1520_M2) < 1535);
      bool CorrM_cut = true;//4500 < Xb_M && Xb_M < 7200;
      bool basic_selection = vetos && PID_cuts && Dalitz_region && p_CosTheta > -0.9;
      //bool Xic_M_cut = 2460 < Xc_M && Xc_M < 2485;
      bool Xic_M_cut = fabs(Xc_M - xicmass) < 16.5;//6.6 MeV with * 2.5 = 16.5

      Xic_beta_hist->Fill(p_beta,Xc_M);
      if(basic_selection){
        if(CorrM_cut){
          Xic_hist->Fill(Xc_M);
          Xic_beta_cut_hist->Fill(p_beta,Xc_M);
          Xic_2D_dump->Fill(Xc_M,p_CosTheta);
        }// p_CosTheta > -0.7 ?
        if(Xic_M_cut)Xib0_hist->Fill(Xb_M);
        else Xib0_SB_hist->Fill(Xb_M);
      }
      bool gs_selection = basic_selection && CorrM_cut && Xic_M_cut;

      if(!gs_selection)continue;

      vector<cH_mult_cand> multiple_cpi_SS;
      vector<cH_mult_cand> multiple_cpi_OS;
      vector<cH_mult_cand> multiple_cK_SS;
      vector<cH_mult_cand> multiple_cK_OS;
      vector < vector<cH_mult_cand> > multiples;
      for(int ap = 0; ap < Added_n_Particles; ap++){
        bool OS = ((Added_H_PROBNNPID[ap] < 0 && Xc_ID > 0) || (Added_H_PROBNNPID[ap] > 0 && Xc_ID < 0));
        bool added_H_cuts = Added_H_PT[ap] > 150 && Added_CharmH_VERTEXCHI2_NEW[ap] < 4 &&  Added_CharmH_PT[ap] > 2000;
        bool added_K_cuts = Added_H_ProbNNk[ap] > 0.25 && Added_CharmH_M_kaon[ap] < 3960 && TMath::Abs(Added_H_PROBNNPID[ap]) == 321;
        bool added_pi_cuts = Added_H_ProbNNpi[ap] > 0.25 && Added_CharmH_M[ap] < 4000 && TMath::Abs(Added_H_PROBNNPID[ap]) == 221;
        bool piSS_selection = !OS && added_H_cuts && added_pi_cuts;//basic_selection &&
        bool piOS_selection =  OS && added_H_cuts && added_pi_cuts;
        bool KSS_selection  = !OS && added_H_cuts && added_K_cuts;
        bool KOS_selection  =  OS && added_H_cuts && added_K_cuts;
        if(piSS_selection)multiple_cpi_SS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+2467.8,Added_CharmH_IP_NEW[ap]});
        if(piOS_selection)multiple_cpi_OS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+2467.8,Added_CharmH_IP_NEW[ap]});
        if(KSS_selection)multiple_cK_SS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+2467.8,Added_CharmH_IP_NEW[ap]});
        if(KOS_selection)multiple_cK_OS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+2467.8,Added_CharmH_IP_NEW[ap]});
      }
      multiples.push_back(multiple_cpi_SS);
      multiples.push_back(multiple_cpi_OS);
      multiples.push_back(multiple_cK_SS);
      multiples.push_back(multiple_cK_OS);

      for(unsigned int i = 0; i < multiples.size(); i++){
        if(multiples.at(i).size() > 1){
          if(myconfig->get_verbosity() > normal)cout << "we have a multiple candidate" << endl;
          sort(multiples.at(i).begin(), multiples.at(i).end(), compareByIP);
          if(i == 1 && make_IP_profile_plots)fill_profile(multiples.at(i),Xic_IP);
          for(unsigned int j = 0; j < multiples.at(i).size(); j++){
            if(myconfig->get_verbosity() > normal)cout << "IP of sorted candidate # " << j << ": " << multiples.at(i).at(j).ip << endl;
          }
        }
        if(clean_vertex){
          if(multiples.at(i).size() > 1){
            int good_IP_particles = 1;
            int ih1 = 0, ih2 = 0;
            for(unsigned int j = 1; j < multiples.at(i).size(); j++){
              //eliminate Xc*->Xc+nH (where n > 1)
              //Assuming that if n candidates have similar IP, they come from the same vertex
              if(multiples.at(i).at(1).ip/multiples.at(i).at(0).ip < 1.25){
                good_IP_particles++;
                ih1 = multiples.at(i).at(0).index; ih2 = multiples.at(i).at(1).index;
                multiples.at(i).erase(multiples.at(i).begin());//it will always be the first element I hope ?!
                if(purge_vertex)multiples.at(i).clear();
              }
            }
          }
          if(!multiples.at(i).empty())Xic_hists.at(i)->Fill(multiples.at(i).at(0).mass);
        }
        else{
          if(!multiples.at(i).empty())Xic_hists.at(i)->Fill(multiples.at(i).at(0).mass);
          if(multiples.at(i).size() > 1){
            for(unsigned int j = 1; j < multiples.at(i).size(); j++){
              Xic_IPfails.at(i)->Fill(multiples.at(i).at(j).mass);
            }
          }
        }
      }
    }
    make_overlay_Plot(Xib0_SB_hist,Xib0_hist,myconfig);
    if(make_IP_profile_plots)make_1D_Plot(Xic_IP,myconfig);
    make_2D_Plot(Xic_beta_hist,myconfig);
    make_2D_Plot(Xic_beta_cut_hist,myconfig);
    make_2D_Plot(Xic_2D_dump,myconfig);
    make_2D_Plot(Xic_beta_hist,myconfig);
    make_overlay_Plot(XicpiSS_hist,XicpiOS_hist,myconfig);
    make_overlay_Plot(XicKSS_hist,XicKOS_hist,myconfig);
    make_overlay_Plot(XicpiSS_IPfail,XicpiOS_IPfail,myconfig);
    make_overlay_Plot(XicKSS_IPfail,XicKOS_IPfail,myconfig);
    RooDataHist *Xic = new RooDataHist("Xc","Xc",Charm_mass,Import(*Xic_hist)) ;
    Fit_Charm(Xic,myconfig);
  }

  myconfig->set_particle("Xic0");
  Xic0_tree->SetBranchStatus("p_ProbNNp",1);
  Xic0_tree->SetBranchStatus("pi_ProbNNpi",1);
  Xic0_tree->SetBranchStatus("SSK1_ProbNNk",1);
  Xic0_tree->SetBranchStatus("SSK2_ProbNNk",1);
  Xic0_tree->SetBranchStatus("Xic02pKKpi.Xib_CorrM",1);
  Xic0_tree->SetBranchStatus("Xic_M",1);
  Xic0_tree->SetBranchStatus("Xic_PT",1);
  Xic0_tree->SetBranchStatus("Xic_ETA",1);
  Xic0_tree->SetBranchStatus("Xic_PHI",1);
  Xic0_tree->SetBranchStatus("Xic_ID",1);
  Xic0_tree->SetBranchStatus("Xic02pKKpi.p_as_KKKpi_M",1);
  Xic0_tree->SetBranchStatus("p_CosTheta",1);
  Xic0_tree->SetBranchStatus("Xic02pKKpi.p_beta",1);
  Xic0_tree->SetBranchStatus("Added_n_Particles",1);
  Xic0_tree->SetBranchStatus("Added_n_Tracks",1);
  Xic0_tree->SetBranchStatus("Added_CharmH_M",1);
  Xic0_tree->SetBranchStatus("Added_CharmH_PT",1);
  Xic0_tree->SetBranchStatus("Added_CharmH_M_kaon",1);
  Xic0_tree->SetBranchStatus("Added_CharmH_VERTEXCHI2_NEW",1);
  Xic0_tree->SetBranchStatus("Added_CharmH_IP_NEW",1);
  Xic0_tree->SetBranchStatus("Added_H_PT",1);
  Xic0_tree->SetBranchStatus("Added_H_ETA",1);
  Xic0_tree->SetBranchStatus("Added_H_PHI",1);
  Xic0_tree->SetBranchStatus("Added_H_PROBNNPID",1);
  Xic0_tree->SetBranchStatus("Added_H_ProbNNpi",1);
  Xic0_tree->SetBranchStatus("Added_H_ProbNNk",1);

  Xic0_tree->SetBranchAddress("p_ProbNNp",&p_ProbNNp);
  Xic0_tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);
  Xic0_tree->SetBranchAddress("SSK1_ProbNNk",&K_ProbNNk);
  Xic0_tree->SetBranchAddress("SSK2_ProbNNk",&K2_ProbNNk);
  Xic0_tree->SetBranchAddress("Xic02pKKpi.Xib_CorrM",&Xb_M);
  Xic0_tree->SetBranchAddress("Xic_M",&Xc_M);
  Xic0_tree->SetBranchAddress("Xic_PT",&Xc_PT);
  Xic0_tree->SetBranchAddress("Xic_ETA",&Xc_ETA);
  Xic0_tree->SetBranchAddress("Xic_PHI",&Xc_PHI);
  Xic0_tree->SetBranchAddress("Xic_ID",&Xc_ID);
  Xic0_tree->SetBranchAddress("Xic02pKKpi.p_as_KKKpi_M",&MisID_D02KKKpi_M);
  Xic0_tree->SetBranchAddress("p_CosTheta",&p_CosTheta);
  Xic0_tree->SetBranchAddress("Xic02pKKpi.p_beta",&p_beta);
  Xic0_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Xic0_tree->SetBranchAddress("Added_n_Tracks",&Added_n_Tracks);
  Xic0_tree->SetBranchAddress("Added_CharmH_M",&Added_CharmH_M);
  Xic0_tree->SetBranchAddress("Added_CharmH_PT",&Added_CharmH_PT);
  Xic0_tree->SetBranchAddress("Added_CharmH_M_kaon",&Added_CharmH_M_kaon);
  Xic0_tree->SetBranchAddress("Added_CharmH_VERTEXCHI2_NEW",&Added_CharmH_VERTEXCHI2_NEW);
  Xic0_tree->SetBranchAddress("Added_CharmH_IP_NEW",&Added_CharmH_IP_NEW);
  Xic0_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Xic0_tree->SetBranchAddress("Added_H_ETA",&Added_H_ETA);
  Xic0_tree->SetBranchAddress("Added_H_PHI",&Added_H_PHI);
  Xic0_tree->SetBranchAddress("Added_H_PROBNNPID",&Added_H_PROBNNPID);
  Xic0_tree->SetBranchAddress("Added_H_ProbNNpi",&Added_H_ProbNNpi);
  Xic0_tree->SetBranchAddress("Added_H_ProbNNk",&Added_H_ProbNNk);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory

  if(Xic0_sel){

    UInt_t nXic0events = Xic0_tree->GetEntries();
    cout << "Entries in Xic0 tree: " << nXic0events << endl;

    for (UInt_t evt = 0; evt < nXic0events;evt++) {
      Xic0_tree->GetEntry(evt);

      bool vetos = !(1855 < MisID_D02KKKpi_M && MisID_D02KKKpi_M < 1875);
      bool PID_cuts = p_ProbNNp > 0.1 && K_ProbNNk > 0.25 && K2_ProbNNk > 0.25 && pi_ProbNNpi > 0.05;
      bool CorrM_cut = true;//4500 < Xb_M && Xb_M < 7200;
      //bool Xic0_M_cut = 2460 < Xc_M && Xc_M < 2485;
      bool Xic0_M_cut = fabs(Xc_M - xic0mass) < 10.75;//6.6 MeV with * 2.5 = 16.5

      bool basic_selection = vetos && PID_cuts;
      Xic0_beta_hist->Fill(p_beta,Xc_M);
      if(basic_selection){
        if(CorrM_cut){
          Xic0_hist->Fill(Xc_M);
          Xic0_beta_cut_hist->Fill(p_beta,Xc_M);
          Xic0_2D_dump->Fill(Xc_M,p_CosTheta);// p_CosTheta > -0.7 ?
        }
        if(Xic0_M_cut)Xib_hist->Fill(Xb_M);
        else Xib_SB_hist->Fill(Xb_M);
      }
      bool gs_selection = basic_selection && CorrM_cut && Xic0_M_cut;
      if(!gs_selection)continue;
      //Added_H_ProbNNpi > 0.4 && Added_H_PT > 200 && Added_CharmH_M < 4000 && TMath::Abs(Added_H_PROBNNPID) == 221 && Added_CharmH_VERTEXCHI2_NEW < 3

      vector<cH_mult_cand> multiple_cpi_SS;
      vector<cH_mult_cand> multiple_cpi_OS;
      vector<cH_mult_cand> multiple_cK_SS;
      vector<cH_mult_cand> multiple_cK_OS;
      vector < vector<cH_mult_cand> > multiples;
      for(int ap = 0; ap < Added_n_Particles; ap++){
        bool OS = ((Added_H_PROBNNPID[ap] < 0 && Xc_ID > 0) || (Added_H_PROBNNPID[ap] > 0 && Xc_ID < 0));
        bool added_H_cuts = Added_H_PT[ap] > 150 && Added_CharmH_VERTEXCHI2_NEW[ap] < 4 && Added_CharmH_PT[ap] > 2000;
        bool added_K_cuts = Added_H_ProbNNk[ap] > 0.25 && Added_CharmH_M_kaon[ap] < 3960 && TMath::Abs(Added_H_PROBNNPID[ap]) == 321;
        bool added_pi_cuts = Added_H_ProbNNpi[ap] > 0.25 && Added_CharmH_M[ap] < 4000 && TMath::Abs(Added_H_PROBNNPID[ap]) == 221;
        bool piSS_selection = !OS && added_H_cuts && added_pi_cuts;//basic_selection &&
        bool piOS_selection =  OS && added_H_cuts && added_pi_cuts;
        bool KSS_selection  = !OS && added_H_cuts && added_K_cuts;
        bool KOS_selection  =  OS && added_H_cuts && added_K_cuts;
        if(piSS_selection)multiple_cpi_SS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+2470.9,Added_CharmH_IP_NEW[ap],ap});
        if(piOS_selection)multiple_cpi_OS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+2470.9,Added_CharmH_IP_NEW[ap],ap});
        if(KSS_selection)multiple_cK_SS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+2470.9,Added_CharmH_IP_NEW[ap],ap});
        if(KOS_selection)multiple_cK_OS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+2470.9,Added_CharmH_IP_NEW[ap],ap});
      }
      multiples.push_back(multiple_cpi_SS);
      multiples.push_back(multiple_cpi_OS);
      multiples.push_back(multiple_cK_SS);
      multiples.push_back(multiple_cK_OS);

      for(unsigned int i = 0; i < multiples.size(); i++){
        if(multiples.at(i).size() > 1){
          if(myconfig->get_verbosity() > normal)cout << "we have a multiple candidate" << endl;
          sort(multiples.at(i).begin(), multiples.at(i).end(), compareByIP);
          if(i == 0 && make_IP_profile_plots)fill_profile(multiples.at(i),Xic0_IP);
          for(unsigned int j = 0; j < multiples.at(i).size(); j++){
            if(myconfig->get_verbosity() > normal)cout << "IP of sorted candidate # " << j << ": " << multiples.at(i).at(j).ip << endl;
          }
        }
        if(clean_vertex){
          if(multiples.at(i).size() > 1){
            int good_IP_particles = 1;
            int ih1 = 0, ih2 = 0;
            for(unsigned int j = 1; j < multiples.at(i).size(); j++){
              //eliminate Xc*->Xc+nH (where n > 1)
              //Assuming that if n candidates have similar IP, they come from the same vertex
              if(multiples.at(i).at(1).ip/multiples.at(i).at(0).ip < 1.25){
                good_IP_particles++;
                ih1 = multiples.at(i).at(0).index; ih2 = multiples.at(i).at(1).index;
                multiples.at(i).erase(multiples.at(i).begin());//it will always be the first element I hope ?!
                if(purge_vertex)multiples.at(i).clear();
              }
            }
            if(good_IP_particles == 2){
              TLorentzVector Xic_lv, h1_lv, h2_lv;
              Xic_lv.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,2470.9);
              h1_lv.SetPtEtaPhiM(Added_H_PT[ih1],Added_H_ETA[ih1],Added_H_PHI[ih1],139.57);
              h2_lv.SetPtEtaPhiM(Added_H_PT[ih2],Added_H_ETA[ih2],Added_H_PHI[ih2],139.57);
              double Xicpipi_M = (Xic_lv + h1_lv + h2_lv).M();
              if(i < 2)Xic0pipi_hists.at(i)->Fill(Xicpipi_M);
            }
            else{
              for(unsigned int j = 1; j < multiples.at(i).size(); j++){
                Xic0_IPfails.at(i)->Fill(multiples.at(i).at(j).mass);
              }
            }
          }
          if(!multiples.at(i).empty())Xic0_hists.at(i)->Fill(multiples.at(i).at(0).mass);
        }
        else{
          if(!multiples.at(i).empty())Xic0_hists.at(i)->Fill(multiples.at(i).at(0).mass);
          if(multiples.at(i).size() > 1){
            for(unsigned int j = 1; j < multiples.at(i).size(); j++){
              Xic0_IPfails.at(i)->Fill(multiples.at(i).at(j).mass);
            }
          }
        }
      }
    }
    make_overlay_Plot(Xib_SB_hist,Xib_hist,myconfig);
    if(make_IP_profile_plots)make_1D_Plot(Xic0_IP,myconfig);
    make_2D_Plot(Xic0_beta_hist,myconfig);
    make_2D_Plot(Xic0_beta_cut_hist,myconfig);
    make_2D_Plot(Xic0_2D_dump,myconfig);
    make_2D_Plot(Xic0_beta_hist,myconfig);
    make_overlay_Plot(Xic0piSS_hist,Xic0piOS_hist,myconfig);
    make_overlay_Plot(Xic0KSS_hist,Xic0KOS_hist,myconfig);
    make_overlay_Plot(Xic0piSS_IPfail,Xic0piOS_IPfail,myconfig);
    make_overlay_Plot(Xic0KSS_IPfail,Xic0KOS_IPfail,myconfig);
    make_overlay_Plot(Xic0pipiSS_hist,Xic0pipiOS_hist,myconfig);
    Charm_mass.SetTitle("M_{inv}(pK^{-}K^{-}pi^{+}} + c.c.)");
    RooDataHist *Xic0 = new RooDataHist("Xc","Xc",Charm_mass,Import(*Xic0_hist)) ;
    Fit_Charm(Xic0,myconfig);
  }

  myconfig->set_particle("Omegac");
  Omegac_tree->SetBranchStatus("p_ProbNNp",1);
  Omegac_tree->SetBranchStatus("pi_ProbNNpi",1);
  Omegac_tree->SetBranchStatus("SSK1_ProbNNk",1);
  Omegac_tree->SetBranchStatus("SSK2_ProbNNk",1);
  Omegac_tree->SetBranchStatus("Omegac2pKKpi.Omegab_CorrM",1);
  Omegac_tree->SetBranchStatus("Omegac_M",1);
  Omegac_tree->SetBranchStatus("Omegac_ID",1);
  Omegac_tree->SetBranchStatus("p_CosTheta",1);
  Omegac_tree->SetBranchStatus("Omegac2pKKpi.p_beta",1);
  Omegac_tree->SetBranchStatus("Added_n_Particles",1);
  Omegac_tree->SetBranchStatus("Added_n_Tracks",1);
  Omegac_tree->SetBranchStatus("Added_CharmH_M",1);
  Omegac_tree->SetBranchStatus("Added_CharmH_PT",1);
  Omegac_tree->SetBranchStatus("Added_CharmH_M_kaon",1);
  Omegac_tree->SetBranchStatus("Added_CharmH_VERTEXCHI2_NEW",1);
  Omegac_tree->SetBranchStatus("Added_CharmH_IP_NEW",1);
  Omegac_tree->SetBranchStatus("Added_H_PT",1);
  Omegac_tree->SetBranchStatus("Added_H_PROBNNPID",1);
  Omegac_tree->SetBranchStatus("Added_H_ProbNNpi",1);
  Omegac_tree->SetBranchStatus("Added_H_ProbNNk",1);

  Omegac_tree->SetBranchAddress("p_ProbNNp",&p_ProbNNp);
  Omegac_tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);
  Omegac_tree->SetBranchAddress("SSK1_ProbNNk",&K_ProbNNk);
  Omegac_tree->SetBranchAddress("SSK2_ProbNNk",&K2_ProbNNk);
  Omegac_tree->SetBranchAddress("Omegac2pKKpi.Omegab_CorrM",&Xb_M);
  Omegac_tree->SetBranchAddress("Omegac_M",&Xc_M);
  Omegac_tree->SetBranchAddress("Omegac_ID",&Xc_ID);
  Omegac_tree->SetBranchAddress("p_CosTheta",&p_CosTheta);
  Omegac_tree->SetBranchAddress("Omegac2pKKpi.p_beta",&p_beta);
  Omegac_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Omegac_tree->SetBranchAddress("Added_n_Tracks",&Added_n_Tracks);
  Omegac_tree->SetBranchAddress("Added_CharmH_M",&Added_CharmH_M);
  Omegac_tree->SetBranchAddress("Added_CharmH_PT",&Added_CharmH_PT);
  Omegac_tree->SetBranchAddress("Added_CharmH_M_kaon",&Added_CharmH_M_kaon);
  Omegac_tree->SetBranchAddress("Added_CharmH_VERTEXCHI2_NEW",&Added_CharmH_VERTEXCHI2_NEW);
  Omegac_tree->SetBranchAddress("Added_CharmH_IP_NEW",&Added_CharmH_IP_NEW);
  Omegac_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Omegac_tree->SetBranchAddress("Added_H_PROBNNPID",&Added_H_PROBNNPID);
  Omegac_tree->SetBranchAddress("Added_H_ProbNNpi",&Added_H_ProbNNpi);
  Omegac_tree->SetBranchAddress("Added_H_ProbNNk",&Added_H_ProbNNk);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory

  if(Omegac_sel){
    UInt_t nOmegacevents = Omegac_tree->GetEntries();
    cout << "Entries in Omegac tree: " << nOmegacevents << endl;

    for (UInt_t evt = 0; evt < nOmegacevents;evt++) {
      Omegac_tree->GetEntry(evt);

      //p_ProbNNp > 0.3 && pi_ProbNNpi > 0.0 && SSK1_ProbNNk > 0.4 && SSK2_ProbNNk > 0.4 && MyFriend.Omegab_CorrM > 5500 && MyFriend.Omegab_CorrM < 6800 && Omegac_M > 2685 && Omegac_M < 2710 && Added_H_ProbNNpi > 0.4 && Added_H_PT > 200 && ((Added_H_PROBNNPID > 0 && Omegac_ID < 0) || (Added_H_PROBNNPID < 0 && Omegac_ID > 0)) && TMath::Abs(Added_H_PROBNNPID) == 221 && Added_CharmH_VERTEXCHI2_NEW < 3


      bool PID_cuts = p_ProbNNp > 0.3 && K_ProbNNk > 0.4 && K2_ProbNNk > 0.4 && pi_ProbNNpi > 0.00;
      bool CorrM_cut = true;//5500 < Xb_M && Xb_M < 7200;
      bool Omegac_M_cut = 2685 < Xc_M && Xc_M < 2710;
      Omegac_beta_hist->Fill(p_beta,Xc_M);
      if(PID_cuts){
        if(CorrM_cut){
          Omegac_hist->Fill(Xc_M);
          Omegac_beta_cut_hist->Fill(p_beta,Xc_M);
          Omegac_2D_dump->Fill(Xc_M,p_CosTheta);// p_CosTheta > -0.7 ?
        }
        if(Omegac_M_cut)Omegab_hist->Fill(Xb_M);
        else Omegab_SB_hist->Fill(Xb_M);
      }
      bool gs_selection = PID_cuts && CorrM_cut && Omegac_M_cut;
      if(!gs_selection)continue;

      vector<cH_mult_cand> multiple_cpi_SS;
      vector<cH_mult_cand> multiple_cpi_OS;
      vector<cH_mult_cand> multiple_cK_SS;
      vector<cH_mult_cand> multiple_cK_OS;
      vector < vector<cH_mult_cand> > multiples;
      for(int ap = 0; ap < Added_n_Particles; ap++){
        bool OS = ((Added_H_PROBNNPID[ap] < 0 && Xc_ID > 0) || (Added_H_PROBNNPID[ap] > 0 && Xc_ID < 0));
        bool added_H_cuts = Added_H_PT[ap] > 150 && Added_CharmH_VERTEXCHI2_NEW[ap] < 4 && Added_CharmH_PT[ap] > 2000;
        bool added_K_cuts = Added_H_ProbNNk[ap] > 0.25 && Added_CharmH_M_kaon[ap] < 3990 && TMath::Abs(Added_H_PROBNNPID[ap]) == 321;
        bool added_pi_cuts = Added_H_ProbNNpi[ap] > 0.25 && Added_CharmH_M[ap] < 4040 && TMath::Abs(Added_H_PROBNNPID[ap]) == 221;
        bool piSS_selection = !OS && added_H_cuts && added_pi_cuts;//basic_selection &&
        bool piOS_selection =  OS && added_H_cuts && added_pi_cuts;
        bool KSS_selection  = !OS && added_H_cuts && added_K_cuts;
        bool KOS_selection  =  OS && added_H_cuts && added_K_cuts;
        if(piSS_selection)multiple_cpi_SS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+2695.2,Added_CharmH_IP_NEW[ap]});
        if(piOS_selection)multiple_cpi_OS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+2695.2,Added_CharmH_IP_NEW[ap]});
        if(KSS_selection)multiple_cK_SS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+2695.2,Added_CharmH_IP_NEW[ap]});
        if(KOS_selection)multiple_cK_OS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+2695.2,Added_CharmH_IP_NEW[ap]});
      }
      multiples.push_back(multiple_cpi_SS);
      multiples.push_back(multiple_cpi_OS);
      multiples.push_back(multiple_cK_SS);
      multiples.push_back(multiple_cK_OS);

      for(unsigned int i = 0; i < multiples.size(); i++){
        if(multiples.at(i).size() > 1){
          if(myconfig->get_verbosity() > normal)cout << "we have a multiple candidate" << endl;
          sort(multiples.at(i).begin(), multiples.at(i).end(), compareByIP);
          if(i == 0 && make_IP_profile_plots)fill_profile(multiples.at(i),Omegac_IP);
          for(unsigned int j = 0; j < multiples.at(i).size(); j++){
            if(myconfig->get_verbosity() > normal)cout << "IP of sorted candidate # " << j << ": " << multiples.at(i).at(j).ip << endl;
          }
        }
        if(clean_vertex){
          if(multiples.at(i).size() > 1){
            int good_IP_particles = 1;
            int ih1 = 0, ih2 = 0;
            for(unsigned int j = 1; j < multiples.at(i).size(); j++){
              //eliminate Xc*->Xc+nH (where n > 1)
              //Assuming that if n candidates have similar IP, they come from the same vertex
              if(multiples.at(i).at(1).ip/multiples.at(i).at(0).ip < 1.25){
                good_IP_particles++;
                ih1 = multiples.at(i).at(0).index; ih2 = multiples.at(i).at(1).index;
                multiples.at(i).erase(multiples.at(i).begin());//it will always be the first element I hope ?!
                if(purge_vertex)multiples.at(i).clear();
              }
            }
          }
          if(!multiples.at(i).empty())Omegac_hists.at(i)->Fill(multiples.at(i).at(0).mass);
        }
        else{
          if(!multiples.at(i).empty())Omegac_hists.at(i)->Fill(multiples.at(i).at(0).mass);
          if(multiples.at(i).size() > 1){
            for(unsigned int j = 1; j < multiples.at(i).size(); j++){
              Omegac_IPfails.at(i)->Fill(multiples.at(i).at(j).mass);
            }
          }
        }
      }
    }

    make_overlay_Plot(Omegab_SB_hist,Omegab_hist,myconfig);
    if(make_IP_profile_plots)make_1D_Plot(Omegac_IP,myconfig);
    make_2D_Plot(Omegac_beta_hist,myconfig);
    make_2D_Plot(Omegac_beta_cut_hist,myconfig);
    make_2D_Plot(Omegac_2D_dump,myconfig);
    make_2D_Plot(Omegac_beta_hist,myconfig);
    make_overlay_Plot(OmegacpiSS_hist,OmegacpiOS_hist,myconfig);
    make_overlay_Plot(OmegacKSS_hist,OmegacKOS_hist,myconfig);
    make_overlay_Plot(OmegacpiSS_IPfail,OmegacpiOS_IPfail,myconfig);
    make_overlay_Plot(OmegacKSS_IPfail,OmegacKOS_IPfail,myconfig);
    Charm_mass.setRange(myconfig->get_XcMlo(),myconfig->get_XcMhi());
    RooDataHist *Omegac = new RooDataHist("Xc","Xc",Charm_mass,Import(*Omegac_hist)) ;
    Fit_Charm(Omegac,myconfig);
  }

  clock->Stop();clock->Print();
  Xic_tree->SetDirectory(0);
  Xic0_tree->SetDirectory(0);
  Omegac_tree->SetDirectory(0);
  fSLBS->Close();
  cout << "SUCCESS" << endl;
  return;
}

template <class cH>
void fill_profile(vector<cH> multiple_cH_candidate, TProfile *Xc_IP) {
  int n_profilebins = Xc_IP->GetNbinsX();
  double binwidth = Xc_IP->GetXaxis()->GetBinWidth(1);
  vector<int> multiplicity;
  for (int i = 1; i < n_profilebins + 1; i++){
    multiplicity.push_back(0);
    for (int j = 0; j < multiple_cH_candidate.size(); j++){
      if(i < n_profilebins){
        if(multiple_cH_candidate.at(j).ip < (float)(i*binwidth))
          multiplicity[i-1] += 1;
      }
      else multiplicity[i-1] += 1;//get all
    }
  }
  for (int i = 1; i < n_profilebins + 1; i++){
    //cout << "multiplicity below " << i*binwidth << " mm: " << multiplicity.at(i-1) << endl;
    Xc_IP->Fill(Xc_IP->GetXaxis()->GetBinCenter(i),(double)multiplicity.at(i-1));
  }
  return;
}

template<class T>
void make_1D_Plot(T *hist, configuration *myconfig){
  TCanvas *c1 = new TCanvas("c1","C1",10,10,1280,960);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.14);
  gPad->SetTickx();
  gPad->SetTicky();
  //gPad->SetLogy();
  //hist->GetYaxis()->SetTitleOffset(0.8);
  hist->SetLineWidth(2);
  if(static_cast<TString>(hist->GetName()).Contains("_IP")){
    hist->SetLineColor(kBlack);
    hist->SetMarkerColor(kBlack);
  }
  else hist->SetLineColor(kBlue+2);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetRangeUser(0,1.2*hist->GetMaximum());
  if(static_cast<TString>(hist->GetName()).Contains("_IP"))hist->Draw("e1p");
  else hist->Draw("hist");

  TPaveText *blank_box = new TPaveText(1.001-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();
  hist->Draw("axis same");

  if(static_cast<TString>(hist->GetName()).Contains("_IP")){
    TLegend *leg;
    leg = new TLegend(0.5,0.12+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.17+gPad->GetBottomMargin());
    leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
    temp = "#Xi_{c}^{+}+n negative Tracks";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0} + n positive Tracks";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0} + n positive Tracks^{+}";
    leg->AddEntry(hist,temp,"lp");
    leg->Draw();
  }

  temp = myconfig->get_dumpdir()+"/Plots/";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  temp += static_cast<TString>(hist->GetName()) + "_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);
  delete c1;
  return;
}

void make_2D_Plot(TH2D *hist, configuration *myconfig){
  TCanvas *c1 = new TCanvas("c1","C1",10,10,1280,960);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.14);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  gPad->SetTickx();
  gPad->SetTicky();
  //gPad->SetLogy();
  hist->GetYaxis()->SetTitleOffset(0.95);
  hist->Draw("colz");

  TPaveText *blank_box = new TPaveText(1.001-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TPaveText *sbb_box = new TPaveText(0.001+gPad->GetLeftMargin(),1.001-gPad->GetTopMargin(),0.1+gPad->GetLeftMargin(),1.0,"BRNDC");
  sbb_box->SetBorderSize(0);sbb_box->SetFillColor(kWhite);sbb_box->SetTextAlign(12);sbb_box->SetFillStyle(1001);
  sbb_box->AddText(" ");
  sbb_box->Draw();

  hist->Draw("z axis same");

  temp = myconfig->get_dumpdir()+"/Plots/";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  temp += static_cast<TString>(hist->GetName()) + "_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);
  delete c1;
  return;
}

void make_overlay_Plot(TH1D *SS_hist, TH1D *OS_hist, configuration *myconfig){
  TCanvas *c1 = new TCanvas("c1","C1",10,10,1280,960);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.14);
  gPad->SetTickx();
  gPad->SetTicky();
  //gPad->SetLogy();
  //hist->GetYaxis()->SetTitleOffset(0.8);

  double ymax = TMath::Max(SS_hist->GetMaximum(),OS_hist->GetMaximum());

  //SS_hist->SetLineWidth(2);
  SS_hist->SetLineColor(kOrange+1);

  //OS_hist->SetLineWidth(2);
  OS_hist->SetLineColor(kBlue+2);

  OS_hist->GetXaxis()->SetLabelSize(0.05);
  OS_hist->GetYaxis()->SetLabelSize(0.05);
  OS_hist->GetXaxis()->SetTitleSize(0.05);
  OS_hist->GetYaxis()->SetTitleSize(0.05);
  OS_hist->GetYaxis()->SetRangeUser(0,1.2*ymax);

  //SS_hist->SetFillColorAlpha(kOrange+1, 0.35);
  //SS_hist->SetFillColor(kOrange+1);
  OS_hist->Draw("hist");
  SS_hist->Draw("histsame");

  TPaveText *blank_box = new TPaveText(1.001-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();
  OS_hist->Draw("axis same");

  if(myconfig->get_particle().Contains("Xic") && static_cast<TString>(OS_hist->GetName()).Contains("pi")){

    c1->Update();
    double Xic_2645_mass = 2646;
    double Xic_2790_mass = 2790 - 107.5;//107.5 : Xic' - Xic mass diffrrence
    double LQCD_Xicc_mass = 3610;

    TArrow Xic_arrows;
    Xic_arrows.SetAngle(40);
    Xic_arrows.SetLineWidth(2);
    Xic_arrows.DrawArrow(Xic_2645_mass,0.85*gPad->GetUymax(),Xic_2645_mass,0.88*gPad->GetUymax(),0.03,"<|");
    Xic_arrows.DrawArrow(Xic_2790_mass,0.6*gPad->GetUymax(),Xic_2790_mass,0.7*gPad->GetUymax(),0.03,"<|");
    Xic_arrows.DrawArrow(LQCD_Xicc_mass,0.10*gPad->GetUymax(),LQCD_Xicc_mass,0.2*gPad->GetUymax(),0.03,"<|");

    TLatex Xic_labels;Xic_labels.SetTextAlign(12);Xic_labels.SetTextFont(42);Xic_labels.SetTextSize(0.05);
    Xic_labels.DrawLatex(Xic_2645_mass-10,0.93*gPad->GetUymax(),"#Xi_{c}(2645)");
    Xic_labels.DrawLatex(Xic_2790_mass-10,0.75*gPad->GetUymax(),"#Xi_{c}(2790)#rightarrow#Xi_{c}'#pi");
    Xic_labels.SetTextAlign(22);
    Xic_labels.DrawLatex(LQCD_Xicc_mass,0.25*gPad->GetUymax(),"LQCD #Xi_{cc}");

  }

  if(myconfig->get_particle().CompareTo("Xic") == 0 && static_cast<TString>(OS_hist->GetName()).Contains("K")){

    c1->Update();

    double Marco_1_QVal = 40.0;
    double Marco_2_QVal = 90.0;
    double Marco_3_QVal = 105.0;
    double Marco_4_QVal = 128.0;
    double Marco_5_QVal = 160.0;

    double peak1_mass = Marco_1_QVal + 2469 + 494;
    double peak2_mass = Marco_2_QVal + 2469 + 494;
    double peak3_mass = Marco_3_QVal + 2469 + 494;
    double peak4_mass = Marco_4_QVal + 2469 + 494;
    double peak5_mass = Marco_5_QVal + 2469 + 494;

    TArrow OmegacStar_arrows;
    OmegacStar_arrows.SetAngle(40);
    OmegacStar_arrows.SetLineWidth(2);
    OmegacStar_arrows.DrawArrow(peak1_mass,0.85*gPad->GetUymax(),peak1_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak2_mass,0.85*gPad->GetUymax(),peak2_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak3_mass,0.85*gPad->GetUymax(),peak3_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak4_mass,0.85*gPad->GetUymax(),peak4_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak5_mass,0.85*gPad->GetUymax(),peak5_mass,0.88*gPad->GetUymax(),0.03,"<|");

    TLatex OmegacStar_labels;OmegacStar_labels.SetTextAlign(21);OmegacStar_labels.SetTextFont(42);OmegacStar_labels.SetTextSize(0.05);
    OmegacStar_labels.DrawLatex(peak1_mass,0.90*gPad->GetUymax(),"1");
    OmegacStar_labels.DrawLatex(peak2_mass,0.90*gPad->GetUymax(),"2");
    OmegacStar_labels.DrawLatex(peak3_mass,0.90*gPad->GetUymax(),"3");
    OmegacStar_labels.DrawLatex(peak4_mass,0.90*gPad->GetUymax(),"4");
    OmegacStar_labels.DrawLatex(peak5_mass,0.90*gPad->GetUymax(),"5");

  }

  TLegend *leg;
  if(static_cast<TString>(SS_hist->GetName()).Contains("SB"))leg = new TLegend(gPad->GetLeftMargin()+0.06,0.82-gPad->GetTopMargin(),0.40-gPad->GetLeftMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.82-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  if(static_cast<TString>(OS_hist->GetName()).Contains("IPfail"))leg->SetHeader("Multiple candidates");
  if(static_cast<TString>(OS_hist->GetName()).Contains("pi")){
    temp = "#Xi_{c}^{#pm}#pi^{#mp} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}#pi^{#mp} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}#pi^{#mp} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("K")){
    temp = "#Xi_{c}^{#pm}K^{#mp} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}K^{#mp} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}K^{#mp} + c.c.";
  }
  else{
    temp = "from #Xi_{c}^{#pm} signal";
    if(myconfig->get_particle().Contains("Xic0"))temp = "from #Xi_{c}^{0} signal";
    if(myconfig->get_particle().Contains("Omega"))temp = "from #Omega_{c}^{0} signal";
  }
  leg->AddEntry(OS_hist,temp,"l");
  if(static_cast<TString>(OS_hist->GetName()).Contains("pi")){
    temp = "#Xi_{c}^{#pm}#pi^{#pm} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}#pi^{#pm} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}#pi^{#pm} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("K")){
    temp = "#Xi_{c}^{#pm}K^{#pm} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}K^{#pm} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}K^{#pm} + c.c.";
  }
  else temp = "from sidebands";
  leg->AddEntry(SS_hist,temp,"l");
  leg->Draw();

  temp = myconfig->get_dumpdir()+"/Plots/";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  temp += static_cast<TString>(OS_hist->GetName()) + "_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);
  delete c1;
  return;
}
