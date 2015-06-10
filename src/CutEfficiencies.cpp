/*

    2015-05-25 M. Stahl

    Determine cut and background efficiencies for a given cut string

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

#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataHist.h"

#include "../include/Fit_Charm_Baryon.h"
#include "../include/configuration.h"
#include "../include/MyStyle.h"

TString temp;

void get_efficiencies(configuration* myconfig);

//http://root.cern.ch/download/doc/primer/ROOTPrimer.html#root-macros
//#ifndef __CINT__
int main(int argc, char **argv)
{
  TString jobname;
  for(int i = 1; i < argc; i++){
    //if statement so that a jobname without whitespace stays without whitespace
    if(i > 1)jobname += " " + TString(argv[i]);
    else jobname += TString(argv[i]);
  }
  configuration* myconfig = new configuration(jobname);

  get_efficiencies(myconfig);

  delete myconfig;
  return 0;
}
//#endif

void get_efficiencies(configuration* myconfig){

  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  MyStyle();

  TString cuts = myconfig->get_cutstring();
  //get the program to work with empty relation operators and values (e.g. for TISTOS efficiencies)
  //make the program work with single cut values
  TString variable, relop, value;
  unsigned int ncuts = 0;
  double lo_cut = 0, hi_cut = 0;
  //don't append histos to files
  TH1D::AddDirectory(false);
  TH1D *pass_hist = new TH1D("pass_hist",";;",myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH1D *fail_hist = (TH1D*)pass_hist->Clone("fail_hist");
  vector<TH1D*> pass_hists;
  vector<TH1D*> fail_hists;
  vector<double> cutvalues;

  if(cuts.Tokenize("#")->GetEntries() == 3){
    variable = ((TObjString*)(cuts.Tokenize("#")->At(0)))->String();
    relop = ((TObjString*)(cuts.Tokenize("#")->At(1)))->String();
    value = ((TObjString*)(cuts.Tokenize("#")->At(2)))->String();
    if(value.Tokenize(":")->GetEntries() == 2){
      ncuts = static_cast<unsigned int>(((TObjString*)(value.Tokenize(":")->At(1)))->String().Atoi() - 1);
      lo_cut = ((TObjString*)((((TObjString*)(value.Tokenize(":")->At(0)))->String()).Tokenize(",")->At(0)))->String().Atof();
      hi_cut = ((TObjString*)((((TObjString*)(value.Tokenize(":")->At(0)))->String()).Tokenize(",")->At(1)))->String().Atof();
      double stepsize = (hi_cut-lo_cut)/(ncuts+1);
      for(unsigned int i = 0; i < ncuts; i++){
        temp.Form("pass_hist_%d",i);
        pass_hists.emplace_back((TH1D*)pass_hist->Clone(temp));
        temp.Form("fail_hist_%d",i);
        fail_hists.emplace_back((TH1D*)fail_hist->Clone(temp));
        double the_value = lo_cut+(i+1)*stepsize;
        cutvalues.push_back(the_value);
        temp.Form(myconfig->get_particle() + "_" + variable + "_" + relop + "_%g",the_value);
        myconfig->fill_cs(temp);
      }
    }
    else myconfig->fill_cs(myconfig->get_particle() + "_" + myconfig->get_cutstring());//need to refine this else a bit
  }
  else myconfig->fill_cs(myconfig->get_particle() + "_" + myconfig->get_cutstring());

  double cutvar = 0, mass;
  bool boolcutvar = true;

  temp = myconfig->get_tupledir() + "/SLBaryonSpectroscopyStrp21.root";
  gErrorIgnoreLevel = kError;
  TFile *fSLBS = new TFile(temp,"read");
  if(myconfig->get_particle().Contains("Xic0"))temp = "Xib2Xic0MuNu/Xic02pKKpi/DecayTree";
  else if(myconfig->get_particle().Contains("Omegac"))temp = "Omegab2Omegac0MuNu/Omegac2pKKpi/DecayTree";
  else temp = "Xib02XicMuNu/Xic2pKpi/DecayTree";
  TTree *SLBS_tree = (TTree*)gDirectory->Get(temp);
  gErrorIgnoreLevel = kPrint;
  SLBS_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)
  SLBS_tree->SetBranchStatus(variable,1);
  if(myconfig->get_particle().Contains("Omegac"))temp = "Omegac_M";
  else temp = "Xic_M";
  SLBS_tree->SetBranchStatus(temp,1);
  if(ncuts > 0)SLBS_tree->SetBranchAddress(variable,&cutvar);
  else SLBS_tree->SetBranchAddress(variable,&boolcutvar);
  SLBS_tree->SetBranchAddress(temp,&mass);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory

  UInt_t nevents = SLBS_tree->GetEntries();
  cout << "Entries in tree: " << nevents << endl;

  //Lambda to parse the relation operator
  auto pass_cut = [=] (double cutvalue, double cutvariable) {
    if(relop.CompareTo("GT")==0)return cutvariable > cutvalue;//}  {cout << cutvariable << " " << cutvalue << endl;
    if(relop.CompareTo("LT")==0)return cutvariable < cutvalue;
    if(relop.CompareTo("EQ")==0)return cutvariable == cutvalue;
    if(relop.CompareTo("NE")==0)return cutvariable != cutvalue;
    if(relop.CompareTo("GE")==0)return cutvariable >= cutvalue;
    if(relop.CompareTo("LE")==0)return cutvariable <= cutvalue;
  };

  //for (UInt_t evt = 0; evt < nevents;evt++) {
  for (UInt_t evt = 0; evt < 100000;evt++) {
    SLBS_tree->GetEntry(evt);
    if(cutvalues.empty())
      boolcutvar ? pass_hist->Fill(mass) : fail_hist->Fill(mass);
    else{
      for(unsigned int i = 0; i < ncuts; i++){
        //TODO:try to implement that fail histos are filled automatically after having been filled once without evaluating the lambda over and over again
        pass_cut(cutvalues.at(i),cutvar) ? static_cast<TH1D*>(pass_hists.at(i))->Fill(mass) : static_cast<TH1D*>(fail_hists.at(i))->Fill(mass);
        //pass_cut(cutvalues.at(i),cutvar) ? cout << "pass " << mass << endl : cout << "fail " << mass << endl;
      }
    }
  }

  SLBS_tree->SetDirectory(0);
  fSLBS->Close();
  delete fSLBS;

  auto fill_ASCII = [&] (double *efficiencies, unsigned int iter) {
    ofstream *ofile = new ofstream;
    temp = myconfig->get_dumpdir() + "/Raw_numbers";
    if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);
    temp = myconfig->get_dumpdir() + "/Raw_numbers/eff_" + myconfig->get_cs_at(iter);
    ofile->open(temp);
    *ofile << efficiencies[0] << "\t" << efficiencies[1] << "\t" << efficiencies[2] << "\t" << efficiencies[3] << endl;
    ofile->close();
    delete ofile;
  };

  RooRealVar Xc_M("Xc_M","M_{inv}(pK^{-}pi^{+}} + c.c.)",myconfig->get_XcMlo(),myconfig->get_XcMhi(),"MeV") ;
  if(myconfig->get_particle().Contains("Xic0") || myconfig->get_particle().Contains("Omegac"))Xc_M.SetTitle("M_{inv}(pK^{-}K^{-}pi^{+}} + c.c.)");

  if(cutvalues.empty()){
    RooDataHist *Xc_cut_pass = new RooDataHist("Xc_cut_pass","Xc_cut_pass",Xc_M,Import(*pass_hist)) ;
    RooDataHist *Xc_cut_fail = new RooDataHist("Xc_cut_fail","Xc_cut_fail",Xc_M,Import(*fail_hist)) ;
    myconfig->set_current_cs(myconfig->get_cs_at(0));
    double *efficiencies = Fit_Charm(Xc_cut_pass,Xc_cut_fail,myconfig);
    if(!myconfig->is_debug())fill_ASCII(efficiencies,0);
    delete Xc_cut_pass;
    delete Xc_cut_fail;
  }
  else{
    for(unsigned int i = 0; i < ncuts; i++){
      RooDataHist *Xc_cut_pass = new RooDataHist("Xc_cut_pass","Xc_cut_pass",Xc_M,Import(*pass_hists.at(i))) ;
      RooDataHist *Xc_cut_fail = new RooDataHist("Xc_cut_fail","Xc_cut_fail",Xc_M,Import(*fail_hists.at(i))) ;
      myconfig->set_current_cs(myconfig->get_cs_at(i));
      double *efficiencies = Fit_Charm(Xc_cut_pass,Xc_cut_fail,myconfig);
      if(!myconfig->is_debug())fill_ASCII(efficiencies,i);
      delete Xc_cut_pass;
      delete Xc_cut_fail;
    }
  }
  clock->Stop();clock->Print();delete clock;
  return;
}
