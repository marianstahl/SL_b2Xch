/*

    2015-05-24 M. Stahl

    Fitter to calculate cut efficiencies for ground state charm baryons from semileptonic Xb -> Xc mu nu transitions
    TODO :
          - Clean up convergence check. E.g. code it in a while loop with function calls and an external counter provided by the configuration class
          - Make fitparameters acessible from python options, so that no recompiling is needed to change them

*/

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "TObject.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLine.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooArgList.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGlobalFunc.h"

#include "configuration.h"
//#include "AsymCalc.h"

using namespace RooFit;
using namespace std;

//global variables
extern TString temp;
bool simple_fits = false;
int cpus_for_fit = 1;

//functions to check if the fit converged and to re-fit the samples
bool convcheck(RooFitResult *fitres, vector<RooRealVar*> vars,vector<RooRealVar*> posis, RooAbsPdf *model, RooDataHist *data, configuration* myconfig);
bool errcheck(vector<RooRealVar*> posis);
bool paratlimit(vector<RooRealVar*> vars);
void setnewparranges(vector<RooRealVar*> vars);

double* Fit_Charm(RooDataHist *datapass, RooDataHist *datafail, configuration* myconfig){

  double binsize = (myconfig->get_XcMhi() - myconfig->get_XcMlo())/myconfig->get_NBinsForFit();

  double pull_threshold = 4.5;

  // (common) Mean Jpsi
  double MXcstart = 2470; double MXclo = 2460; double MXchi = 2478;
  if(myconfig->get_particle().Contains("Omegac")){MXcstart = 2695; MXclo = 2685; MXchi = 2705;}
  // sigma's
  double sigma1start = 8; double sigma1lo = 3; double sigma1hi = 20;
  //double sigma2start = 55; double sigma2lo = 10; double sigma2hi = 100.0;
  //double wfstart = 1.6; double wflo = 1.3; double wfhi = 2.1;

  //Background
  double chebastart = -0.3; double chebalo = -0.9; double chebahi = 0.9;
  double chebbstart = -0.03; double chebblo = -0.3; double chebbhi = 0.3;
  double chebcstart = 0.002; double chebclo = -0.1; double chebchi = 0.1;
  //double taustart = -0.001; double taulo = -0.1; double tauhi = 0.1;

  //normalizations
  //double sfracstart = 0.4; double sfraclo = 0.3; double sfrachi = 0.75;
  double Xcsig1start = 15000; double Xcsiglo = 500; double Xcsighi = 25000;
  double bkgstart = 100000; double bkglo = 100; double bkghi = 1e+7;

  //cout << "building model.." << endl;

  // --- Observable ---
  RooRealVar Xc_M("Xc_M","M_{inv}(pK^{-}pi^{+}} + c.c.)",myconfig->get_XcMlo(),myconfig->get_XcMhi(),"MeV") ;
  if(myconfig->get_particle().Contains("Xic0") || myconfig->get_particle().Contains("Omegac"))Xc_M.SetTitle("M_{inv}(pK^{-}K^{-}pi^{+}} + c.c.)");

  // --- Signal ---
  RooRealVar *Xc_sigmean = new RooRealVar("#mu_{Xc}","#mu_{Xc}",MXcstart,MXclo,MXchi) ;
  RooRealVar *Xc_sigma1 = new RooRealVar("#sigma_{Xc,1}","#sigma_{Xc,1}",sigma1start,sigma1lo,sigma1hi) ;//,sigma1lo,sigma1hi
  RooGaussian Xc_gauss1("Xc_gauss1","Gaussian PDF 1",Xc_M,*Xc_sigmean,*Xc_sigma1) ;

  //RooRealVar *widthfactor = new RooRealVar("widthfactor","widthfactor",wfstart,wflo,wfhi);
  //RooFormulaVar *Xc_sigma2 = new RooFormulaVar("gen","@0*@1",RooArgList(*Xc_sigma1,*widthfactor));
  //RooGaussian Xc_gauss2("Xc_gauss2","Gaussian PDF 2",Xc_M,*Xc_sigmean,*Xc_sigma2) ;
  //RooRealVar *signalfraction = new RooRealVar("signalfraction", "signalfraction",sfracstart,sfraclo,sfrachi);
  //RooAddPdf *signal = new RooAddPdf("signal","s",RooArgList(Xc_gauss1,Xc_gauss2),*signalfraction) ;

  // --- Background ---
  RooRealVar *cheba = new RooRealVar("cheba","chebBackground parameter a",chebastart,chebalo,chebahi ) ;
  RooRealVar *chebb = new RooRealVar("chebb","chebBackground parameter b",chebbstart,chebblo,chebbhi) ;
  RooRealVar *chebc = new RooRealVar("chebc","chebBackground parameter c",chebcstart,chebclo,chebchi) ;
  RooArgList bkg_par_list(*cheba,*chebb,*chebc);
  //RooRealVar *tau = new RooRealVar("tau","#tau",taustart,taulo,tauhi);

  RooChebychev chebbkg("chebbkg","cheb Background PDF",Xc_M,bkg_par_list);
  //RooExponential expbkg("expbkg","exp Background PDF",Xc_M,*tau);

  // --- Signal + Background ---
  RooRealVar *signal_yield = new RooRealVar("N_{Xc}","Xc signal yield",Xcsig1start,Xcsiglo,Xcsighi);
  RooRealVar *background_yield = new RooRealVar("N_{bkg}","background yield",bkgstart,bkglo,bkghi);
  //RooExtendPdf* sig = new RooExtendPdf("sig","Extendend signal PDF",*signal,*signal_yield);
  RooExtendPdf* sig = new RooExtendPdf("sig","Extendend signal PDF",Xc_gauss1,*signal_yield);
  RooExtendPdf* bkg = new RooExtendPdf("bkg","Extendend background PDF",chebbkg,*background_yield);
  //RooExtendPdf* bkg = new RooExtendPdf("bkg","Extendend background PDF",expbkg,*background_yield);
  RooAddPdf* sum = new RooAddPdf("sum","Full Model pdf", RooArgList(*sig, *bkg));

  RooRealVar *effz = new RooRealVar("effz","effz",0.95,0.0,1.00);
  RooFormulaVar *signal_yield_fail = new RooFormulaVar("signal_yield_fail","@0*((1/@1)-1)",RooArgList(*signal_yield,*effz));
  RooExtendPdf* sig_fail = new RooExtendPdf("sig_fail","Extendend signal PDF",Xc_gauss1,*signal_yield_fail);
  //RooExtendPdf* sig_fail = new RooExtendPdf("sig_fail","Extendend signal PDF",*signal,*signal_yield_fail);
  RooAddPdf* sum_fail = new RooAddPdf("sum_fail","Full Model pdf", RooArgList(*sig_fail,*bkg));

  //Some things for convergence checks
  //the vector vars_for_convcheck contains the parameters which are checked for convergence. Took signal and background yield away, to not spoil the efficiencies
  vector<RooRealVar*> vars_for_convcheck;
  vars_for_convcheck.push_back(Xc_sigmean);vars_for_convcheck.push_back(Xc_sigma1);
  //vars_for_convcheck.push_back(widthfactor);vars_for_convcheck.push_back(signalfraction);
  vars_for_convcheck.push_back(cheba);vars_for_convcheck.push_back(chebb);vars_for_convcheck.push_back(chebc);

  vector<RooRealVar*> posis;
  posis.push_back(signal_yield);posis.push_back(background_yield);
  vector<RooRealVar*> fail_posis;
  fail_posis.push_back(background_yield);fail_posis.push_back(effz);

  float pavetextsize = 0.04;
  float labeltextsize = 0.05;

  //if(myconfig->get_verbosity() > 0) fitres = sum->fitTo(*plusdatapass,Strategy(2),Save(true),NumCPU(cpus_for_fit)) ;//Minos(true) ,NumCPU(2)
  //else fitres = sum->fitTo(*plusdatapass,Strategy(2),Save(true),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
  RooFitResult *FitResPass = sum->fitTo(*datapass,Strategy(2),Save(true),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;

  cout << "Fit done!" << endl;

  if(simple_fits == false){
    if(convcheck(FitResPass,vars_for_convcheck,posis,sum,datapass,myconfig)==false){
      cerr << "\033[0;31mError while fitting plusdatapass \033[0m" << endl << endl << endl;
      //return NULL;
    }
  }

  TCanvas* c1 = new TCanvas("c1","Fit pass",10,10,1280,2048) ;
  float padsplitmargin = 1e-6;
  float Left_margin  = 0.15;
  float Right_margin = 0.005;
  float Bigy_margin = 0.13;
  float padsplit = 0.25;//height of the pull histo
  c1->Divide(1,2,padsplitmargin,padsplitmargin);
  c1->cd(1);

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad1 = new TPad("pad1","pad1",padsplitmargin,padsplit,1-padsplitmargin,1-padsplitmargin);
  pad1->Draw();
  pad1->cd();
  pad1->SetBorderSize(0);
  pad1->SetLeftMargin(Left_margin);
  pad1->SetRightMargin(Right_margin);
  pad1->SetTopMargin(Bigy_margin/(1-padsplit));
  pad1->SetBottomMargin(padsplitmargin);
  pad1->SetTickx();
  pad1->SetTicky();
  //pad1->SetLogy();
  c1->Update();

  // --- Plot ---
  RooPlot* XcFramePass = Xc_M.frame(Title(" ")) ;

  datapass->plotOn(XcFramePass,Binning(myconfig->get_NBinsForFit()),MarkerColor(1),LineColor(1)) ;
  sum->plotOn(XcFramePass,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines()) ;
  sum->plotOn(XcFramePass,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4)) ;
  sum->plotOn(XcFramePass,LineColor(4),LineWidth(2)) ;
  RooHist* passpullhist = XcFramePass->pullHist() ;
  int bins_gt_threshold = 0;
  if(simple_fits == false){
    double *pass_pulls = passpullhist->GetY();
    for(int ij = 0; ij < myconfig->get_NBinsForFit(); ij++){
      if(pass_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in datapass \033[0m" << endl << endl;
      //return NULL;
    }
  }
  datapass->plotOn(XcFramePass,Binning(myconfig->get_NBinsForFit()),MarkerColor(1),LineColor(1)) ;

  double NsigPass = signal_yield->getVal();
  double NbkgPass = background_yield->getVal();

  TPaveText *passbox = new TPaveText(0.64,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  passbox->SetName("paramBox");
  passbox->SetBorderSize(0);passbox->SetFillColor(kWhite);passbox->SetTextAlign(12);passbox->SetTextSize(pavetextsize/(1-padsplit));passbox->SetTextFont(42);passbox->SetFillStyle(0);
  //passbox->AddText(methlabel);
  //passbox->AddText("matched #mu^{+} probe legs");
  temp.Form("N_{#Xi^{#pm}_{c}} = %.0f #pm %.0f", NsigPass,signal_yield->getError());
  if(myconfig->get_particle().Contains("Xic0"))temp.Form("N_{#Xi^{0}_{c}} = %.0f #pm %.0f", NsigPass,signal_yield->getError());
  if(myconfig->get_particle().Contains("Omegac"))temp.Form("N_{#Omega^{0}_{c}} = %.0f #pm %.0f", NsigPass,signal_yield->getError());
  passbox->AddText(temp);
  temp.Form("N_{Bkg} = %.0f #pm %.0f", NbkgPass,background_yield->getError());
  passbox->AddText(temp);
  XcFramePass->addObject(passbox) ;

  TPaveText *passcutstringbox = new TPaveText(pad1->GetLeftMargin()+0.03,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),0.4,1-pad1->GetTopMargin()-0.03,"BRNDC");
  passcutstringbox->SetBorderSize(0);passcutstringbox->SetFillColor(kWhite);passcutstringbox->SetTextAlign(12);passcutstringbox->SetTextSize(pavetextsize/(1-padsplit));passcutstringbox->SetTextFont(42);passcutstringbox->SetFillStyle(0);
  temp.Form("Passed "+ myconfig->get_current_cs());
  passcutstringbox->AddText(temp);
  XcFramePass->addObject(passcutstringbox);

  TPaveText *pt3 = new TPaveText(pad1->GetLeftMargin(),1.002-pad1->GetTopMargin(),pad1->GetLeftMargin()+0.1,1,"BRNDC");
  pt3->SetBorderSize(0);pt3->SetFillColor(kWhite);pt3->SetFillStyle(1001);pt3->SetTextColor(kWhite);
  pt3->AddText(" ");
  XcFramePass->addObject(pt3);

  TPaveText *pt4 = new TPaveText(1.002-pad1->GetRightMargin(),pad1->GetBottomMargin(),1,pad1->GetBottomMargin()+0.05,"BRNDC");
  pt4->SetBorderSize(0);pt4->SetFillColor(kWhite);pt4->SetFillStyle(1001);pt4->SetTextColor(kWhite);
  pt4->AddText(" ");
  XcFramePass->addObject(pt4);

  XcFramePass->GetXaxis()->SetNdivisions(305) ;
  char yaxislabel[100];
  if(XcFramePass->GetMaximum() > 1e+3 && XcFramePass->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(XcFramePass->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(XcFramePass->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  XcFramePass->GetYaxis()->SetTitle(yaxislabel) ;
  XcFramePass->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  XcFramePass->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  XcFramePass->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  XcFramePass->GetYaxis()->SetRangeUser(0.1,XcFramePass->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  XcFramePass->Draw() ;

  c1->cd(1);

  TPad *pad2 = new TPad("pad2","pad2",padsplitmargin,padsplitmargin,1-padsplitmargin,padsplit);
  pad2->SetBorderSize(0);
  pad2->SetLeftMargin(Left_margin);
  pad2->SetRightMargin(Right_margin);
  pad2->SetTopMargin(padsplitmargin);
  pad2->SetBottomMargin(padsplitmargin);
  pad2->SetTickx();
  pad2->SetTicky();

  RooPlot* PullFramePass = Xc_M.frame(Title(" ")) ;

  Color_t pullcolor = kAzure+3;
  passpullhist->SetLineColor(pullcolor);
  passpullhist->SetMarkerColor(pullcolor);

  PullFramePass->GetXaxis()->SetTickLength(XcFramePass->GetXaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  PullFramePass->GetXaxis()->SetNdivisions(305);
  //h18->GetYaxis()->SetTickLength(Jpsi_Mframeplus->GetYaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  PullFramePass->GetYaxis()->SetLabelSize(XcFramePass->GetYaxis()->GetLabelSize()/(padsplit/(1-padsplit)));
  PullFramePass->GetYaxis()->SetTitleSize(XcFramePass->GetYaxis()->GetTitleSize()/(padsplit/(1-padsplit)));
  PullFramePass->GetYaxis()->SetNdivisions(205);
  PullFramePass->GetYaxis()->SetTitleOffset(1.2*padsplit);
  PullFramePass->GetYaxis()->SetTitle("Pull ");

  TLine* zeroline = new TLine(myconfig->get_XcMlo(),0,myconfig->get_XcMhi(),0);
  zeroline->SetLineStyle(2);
  zeroline->SetLineColor(1);
  TLine* m3sigma = new TLine(myconfig->get_XcMlo(),-3,myconfig->get_XcMhi(),-3);
  m3sigma->SetLineStyle(2);
  m3sigma->SetLineColor(1);
  TLine* p3sigma = new TLine(myconfig->get_XcMlo(),3,myconfig->get_XcMhi(),3);
  p3sigma->SetLineStyle(2);
  p3sigma->SetLineColor(1);
  PullFramePass->addObject(zeroline,"L") ;
  PullFramePass->addObject(m3sigma,"L") ;
  PullFramePass->addObject(p3sigma,"L") ;
  TPaveText *pt5 = new TPaveText(1.002-pad2->GetRightMargin(),pad2->GetBottomMargin(),1,pad2->GetBottomMargin()+0.05,"BRNDC");
  pt5->SetBorderSize(0);pt5->SetFillColor(kWhite);pt5->SetFillStyle(1001);pt5->SetTextColor(kWhite);
  pt5->AddText(" ");
  PullFramePass->addObject(pt5);
  PullFramePass->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);

  pad2->Draw();
  pad2->cd();
  PullFramePass->addPlotable(passpullhist,"P") ;
  PullFramePass->Draw();
  c1->Update();


  RooFitResult *fitResFail = sum_fail->fitTo(*datafail,Strategy(2),Optimize(0),Save(true),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false)) ;
  /*if(myconfig->get_verbosity() > 0) fitResFail = sum_fail->fitTo(*datafail,Strategy(2),Save(true),NumCPU(cpus_for_fit)) ;
  else fitResFail = sum_fail->fitTo(*datafail,Strategy(2),Save(true),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false)) ;//,EvalErrorWall(false)*/

  if(simple_fits == false){
    if(!convcheck(fitResFail,vars_for_convcheck,fail_posis,sum_fail,datafail,myconfig)){
      cerr << "\033[0;31mError while fitting plusdatafail\033[0m" << endl << endl << endl;
      //return NULL;
    }
    if(effz->getMax() > 1)effz->setMax(1);
    if(effz->getMax() < 0)effz->setMin(0);
  }

  /*RooAbsReal *neg_log_likelihood = sum_fail->createNLL(*plusdatafail,Strategy(2),Save(true),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false)) ;
    RooMinuit(*neg_log_likelihood).migrad() ;*/

  /*double reconsigplus = signal_yield->getVal();
    double reconbkgplus = background_yield->getVal();*/    
  double cutefficiency = effz->getVal();
  /*double lo_delta_epsilonplus = fabs(effz->getErrorLo());
    double hi_delta_epsilonplus = effz->getErrorHi();*/
  double delta_ceff = effz->getError();
  double NbkgFail = background_yield->getVal();
  double bkgefficiency = NbkgPass/(NbkgPass+NbkgFail);
  double delta_beff = sqrt((bkgefficiency*(1-bkgefficiency))/(NbkgPass+NbkgFail));

  c1->cd(2);

  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  float lowerpadsplit = padsplit+(Bigy_margin/(1-padsplit));

  TPad *pad3 = new TPad("pad3","pad3",padsplitmargin,lowerpadsplit,1-padsplitmargin,1-padsplitmargin);
  pad3->Draw();
  pad3->cd();
  pad3->SetBorderSize(0);
  pad3->SetLeftMargin(Left_margin);
  pad3->SetRightMargin(Right_margin);
  pad3->SetTopMargin(padsplitmargin);
  pad3->SetBottomMargin(padsplitmargin);
  pad3->SetTickx();
  pad3->SetTicky();
  //pad3->SetLogy();
  c1->Update();

  // --- Plot ---
  RooPlot* XcFrameFail = Xc_M.frame(Title(" ")) ;

  datafail->plotOn(XcFrameFail,Binning(myconfig->get_NBinsForFit()),MarkerColor(1),LineColor(1));
  sum_fail->plotOn(XcFrameFail,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines());
  sum_fail->plotOn(XcFrameFail,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4));
  sum_fail->plotOn(XcFrameFail,LineColor(4),LineWidth(2));
  RooHist* failpullhist = XcFrameFail->pullHist();
  if(simple_fits == false){
    double *fail_pulls = failpullhist->GetY();
    bins_gt_threshold = 0;
    for(int ij = 0; ij < myconfig->get_NBinsForFit(); ij++){
      if(fail_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high in datafail \033[0m" << endl << endl;
      //return NULL;
    }
  }
  datafail->plotOn(XcFrameFail,Binning(myconfig->get_NBinsForFit()),MarkerColor(1),LineColor(1));

  TPaveText *effbox = new TPaveText(0.60,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  effbox->SetName("paramBox");
  effbox->SetBorderSize(0);effbox->SetFillColor(kWhite);effbox->SetTextAlign(12);effbox->SetTextSize(pavetextsize/(1-padsplit));effbox->SetTextFont(42);effbox->SetFillStyle(0);
  //effbox->AddText(methlabel);
  //effbox->AddText("Fail cut");
  temp.Form( "#varepsilon_{sig} = %.2f #pm %.2f %%", 100*cutefficiency,100*delta_ceff);
  effbox->AddText(temp);
  temp.Form( "#varepsilon_{bkg} = %.2f #pm %.2f %%", 100*bkgefficiency,100*delta_beff);
  effbox->AddText(temp);

  TPaveText *failcutstringbox = new TPaveText(pad1->GetLeftMargin()+0.03,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),0.4,1-pad1->GetTopMargin()-0.03,"BRNDC");
  failcutstringbox->SetBorderSize(0);failcutstringbox->SetFillColor(kWhite);failcutstringbox->SetTextAlign(12);failcutstringbox->SetTextSize(pavetextsize/(1-padsplit));failcutstringbox->SetTextFont(42);failcutstringbox->SetFillStyle(0);
  temp.Form("Fail "+ myconfig->get_current_cs());
  failcutstringbox->AddText(temp);
  XcFramePass->addObject(failcutstringbox);

  XcFrameFail->addObject(effbox) ;
  //add the other PaveTexts
  XcFrameFail->addObject(failcutstringbox) ;
  XcFrameFail->addObject(pt3);
  XcFrameFail->addObject(pt4);

  if(XcFrameFail->GetMaximum() > 1e+3 && XcFramePass->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(XcFrameFail->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(XcFrameFail->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);

  XcFrameFail->GetXaxis()->SetNdivisions(305) ;
  XcFrameFail->GetYaxis()->SetTitle(yaxislabel) ;
  XcFrameFail->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  XcFrameFail->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  XcFrameFail->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  //if(meth == 1)Jpsi_Mframeplus_fail->GetYaxis()->SetRangeUser(0.1,Jpsi_Mframeplus_fail->GetMaximum());
  //else
  XcFrameFail->GetYaxis()->SetRangeUser(0.1,1.4*XcFrameFail->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  XcFrameFail->Draw() ;

  /*TCanvas *c35 = new TCanvas("c35","Effz likelihood",1280,960);
    RooNLLVar *nll = new RooNLLVar("nll","nll",*sum_fail,*plusdatafail);
    RooPlot *effzframe = effz->frame(Title(" "),Range(0.98,1));
    RooAbsReal* pll = neg_log_likelihood->createProfile(effz) ;
    nll->plotOn(effzframe,ShiftToZero(),Precision(1e-4));
    pll->plotOn(effzframe,ShiftToZero(),Precision(1e-4),LineColor(2),LineStyle(kDashed));
    effzframe->Draw();
    c35->SaveAs("NLL.pdf");*/

  c1->cd(2);

  failpullhist->SetLineColor(pullcolor);
  failpullhist->SetMarkerColor(pullcolor);

  TPad *pad6 = (TPad*)pad2->Clone("pad6");
  pad6->Draw();
  pad6->cd();
  RooPlot* PullFrameFail = (RooPlot*)PullFramePass->Clone("PullFrameFail");
  PullFrameFail->addPlotable(failpullhist,"P") ;
  PullFrameFail->GetYaxis()->SetRangeUser(-pull_threshold,pull_threshold);
  PullFrameFail->Draw();
  c1->Update();

  temp = myconfig->get_dumpdir()+"/Fits";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);
  temp += "/" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);

  delete c1;  

  //this is the important stuff
  double* efficiencies;
  efficiencies = new double[4];
  efficiencies[0] = cutefficiency;
  if(cutefficiency + delta_ceff < 1 && cutefficiency - delta_ceff > 0)efficiencies[1] = delta_ceff;
  else {
    if(cutefficiency + delta_ceff > 1){
      efficiencies[1] = delta_ceff;
      if(myconfig->get_verbosity() > 0)cout << "efficiency + upper uncertainty exceeds 1, be careful" << endl;
    }
    if(cutefficiency - delta_ceff < 0){
      cout << "efficiency - lower uncertainty smaller than 0, this shouldn't happen - return NULL!" << endl;
      return NULL;
    }
  }
  efficiencies[2] = bkgefficiency;
  efficiencies[3] = delta_beff;

  cout << "\033[0;33mCut efficiency = ( " << 100*cutefficiency << " +- "<< 100*efficiencies[1] << " ) %" << endl;
  cout << "Background efficiency = ( " << 100*bkgefficiency << " +- "<< 100*delta_beff << " ) %\033[0m" << endl;

  delete passbox;
  delete passcutstringbox;  
  delete failcutstringbox;
  delete passpullhist;
  delete failpullhist;
  delete effbox;  
  delete zeroline;
  delete m3sigma;
  delete p3sigma;
  delete pt3;
  delete pt4;
  delete pt5;  
  return efficiencies;

}


void Fit_Charm(RooDataHist *data, configuration* myconfig){

  double binsize = (myconfig->get_XcMhi() - myconfig->get_XcMlo())/myconfig->get_NBinsForFit();

  double pull_threshold = 4.5;

  // (common) Mean Jpsi
  double MXcstart = 2470; double MXclo = 2460; double MXchi = 2478;
  if(myconfig->get_particle().Contains("Omegac")){MXcstart = 2695; MXclo = 2685; MXchi = 2705;}
  // sigma's
  double sigma1start = 8; double sigma1lo = 2; double sigma1hi = 20;
  //double sigma2start = 55; double sigma2lo = 10; double sigma2hi = 100.0;
  double wfstart = 1.6; double wflo = 1.3; double wfhi = 2.1;

  //Background
  double chebastart = -0.3; double chebalo = -0.9; double chebahi = 0.9;
  double chebbstart = -0.03; double chebblo = -0.3; double chebbhi = 0.3;
  double chebcstart = 0.002; double chebclo = -0.1; double chebchi = 0.1;
  //double taustart = -0.001; double taulo = -0.1; double tauhi = 0.1;

  //normalizations
  double sfracstart = 0.4; double sfraclo = 0.3; double sfrachi = 0.75;
  double Xcsig1start = 15000; double Xcsiglo = 0.0; double Xcsighi = 25000;
  double bkgstart = 100000; double bkglo = 100; double bkghi = 1e+7;

  //cout << "building model.." << endl;

  // --- Observable ---
  RooRealVar Xc_M("Xc_M","M_{inv}(pK^{-}pi^{+} + c.c.)",myconfig->get_XcMlo(),myconfig->get_XcMhi(),"GeV") ;
  if(myconfig->get_particle().Contains("Xic0") || myconfig->get_particle().Contains("Omegac"))Xc_M.SetTitle("M_{inv}(pK^{-}K^{-}pi^{+} + c.c.)");

  // --- Signal ---
  RooRealVar *Xc_sigmean = new RooRealVar("#mu_{Xc}","#mu_{Xc}",MXcstart,MXclo,MXchi) ;
  RooRealVar *Xc_sigma1 = new RooRealVar("#sigma_{Xc,1}","#sigma_{Xc,1}",sigma1start,sigma1lo,sigma1hi) ;//,sigma1lo,sigma1hi
  RooGaussian Xc_gauss1("Xc_gauss1","Gaussian PDF 1",Xc_M,*Xc_sigmean,*Xc_sigma1) ;

  RooRealVar *widthfactor = new RooRealVar("widthfactor","widthfactor",wfstart,wflo,wfhi);
  RooFormulaVar *Xc_sigma2 = new RooFormulaVar("gen","@0*@1",RooArgList(*Xc_sigma1,*widthfactor));
  RooGaussian Xc_gauss2("Xc_gauss2","Gaussian PDF 2",Xc_M,*Xc_sigmean,*Xc_sigma2) ;
  RooRealVar *signalfraction = new RooRealVar("signalfraction", "signalfraction",sfracstart,sfraclo,sfrachi);
  RooAddPdf *signal = new RooAddPdf("signal","s",RooArgList(Xc_gauss1,Xc_gauss2),*signalfraction) ;

  // --- Background ---
  RooRealVar *cheba = new RooRealVar("cheba","chebBackground parameter a",chebastart,chebalo,chebahi ) ;
  RooRealVar *chebb = new RooRealVar("chebb","chebBackground parameter b",chebbstart,chebblo,chebbhi) ;
  RooRealVar *chebc = new RooRealVar("chebc","chebBackground parameter c",chebcstart,chebclo,chebchi) ;
  RooArgList bkg_par_list(*cheba,*chebb,*chebc);
  //RooRealVar *tau = new RooRealVar("tau","#tau",taustart,taulo,tauhi);

  RooChebychev chebbkg("chebbkg","cheb Background PDF",Xc_M,bkg_par_list);
  //RooExponential expbkg("expbkg","exp Background PDF",Xc_M,*tau);

  // --- Signal + Background ---
  RooRealVar *signal_yield = new RooRealVar("N_{Xc}","Xc signal yield",Xcsig1start,Xcsiglo,Xcsighi);
  RooRealVar *background_yield = new RooRealVar("N_{bkg}","background yield",bkgstart,bkglo,bkghi);
  RooExtendPdf* sig;
  if(myconfig->get_particle().Contains("Omegac"))sig = new RooExtendPdf("sig","Extendend signal PDF",Xc_gauss1,*signal_yield);
  else sig = new RooExtendPdf("sig","Extendend signal PDF",*signal,*signal_yield);
  //RooExtendPdf* sig = new RooExtendPdf("sig","Extendend signal PDF",Xc_gauss1,*signal_yield);
  RooExtendPdf* bkg = new RooExtendPdf("bkg","Extendend background PDF",chebbkg,*background_yield);
  //RooExtendPdf* bkg = new RooExtendPdf("bkg","Extendend background PDF",expbkg,*background_yield);
  RooAddPdf* sum = new RooAddPdf("sum","Full Model pdf", RooArgList(*sig, *bkg));

  //Some things for convergence checks
  //the vector vars_for_convcheck contains the parameters which are checked for convergence. Took signal and background yield away, to not spoil the efficiencies
  vector<RooRealVar*> vars_for_convcheck;
  vars_for_convcheck.push_back(Xc_sigmean);vars_for_convcheck.push_back(Xc_sigma1);
  if(!myconfig->get_particle().Contains("Omegac"))vars_for_convcheck.push_back(widthfactor);vars_for_convcheck.push_back(signalfraction);
  vars_for_convcheck.push_back(cheba);vars_for_convcheck.push_back(chebb);vars_for_convcheck.push_back(chebc);

  vector<RooRealVar*> posis;
  posis.push_back(signal_yield);posis.push_back(background_yield);

  float pavetextsize = 0.04;
  float labeltextsize = 0.05;

  RooFitResult *FitRes = sum->fitTo(*data,Strategy(2),Save(true),NumCPU(cpus_for_fit),Optimize(0)) ;//PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)

  if(simple_fits == false){
    if(convcheck(FitRes,vars_for_convcheck,posis,sum,data,myconfig)==false){
      cerr << "\033[0;31mError while fitting plusdatapass \033[0m" << endl << endl << endl;
      //return NULL;
    }
  }

  TCanvas* c1 = new TCanvas("c1","Fit pass",10,10,1280,1280) ;
  float padsplitmargin = 1e-6;
  float Left_margin  = 0.15;
  float Top_margin  = 0.02;
  float Right_margin = 0.01;
  float Bigy_margin = 0.15;
  float padsplit = 0.32;//height of the pull histo

  //Must control the margins in the individual pads and not the subpads of the canvas, because of the axis labels
  gPad->SetTopMargin(padsplitmargin);
  gPad->SetBottomMargin(padsplitmargin);
  gPad->SetRightMargin(padsplitmargin);
  gPad->SetLeftMargin(padsplitmargin);

  TPad *pad1 = new TPad("pad1","pad1",padsplitmargin,padsplit,1-padsplitmargin,1-padsplitmargin);
  pad1->Draw();
  pad1->cd();
  pad1->SetBorderSize(0);
  pad1->SetLeftMargin(Left_margin);
  pad1->SetRightMargin(Right_margin);
  pad1->SetTopMargin(Top_margin);
  pad1->SetBottomMargin(padsplitmargin);
  pad1->SetTickx();
  pad1->SetTicky();
  //pad1->SetLogy();
  c1->Update();

  // --- Plot ---
  RooPlot* XcFrame = Xc_M.frame(Title(" ")) ;

  data->plotOn(XcFrame,Binning(myconfig->get_NBinsForFit()),MarkerColor(1),LineColor(1)) ;
  sum->plotOn(XcFrame,Components("bkg"),FillColor(5),LineColor(5),LineWidth(0),DrawOption("F"),VLines()) ;
  sum->plotOn(XcFrame,Components("bkg"),LineStyle(kDashed),LineWidth(2),LineColor(4)) ;
  sum->plotOn(XcFrame,LineColor(4),LineWidth(2)) ;
  RooHist* pullhist = XcFrame->pullHist() ;
  int bins_gt_threshold = 0;
  if(simple_fits == false){
    double *pass_pulls = pullhist->GetY();
    for(int ij = 0; ij < myconfig->get_NBinsForFit(); ij++){
      if(pass_pulls[ij] > pull_threshold)bins_gt_threshold++;
    }
    if(bins_gt_threshold >= 2){
      cerr << "\033[0;31mPulls too high \033[0m" << endl << endl;
      //return NULL;
    }
  }
  data->plotOn(XcFrame,Binning(myconfig->get_NBinsForFit()),MarkerColor(1),LineColor(1)) ;

  double Nsig = signal_yield->getVal();
  double Nbkg = background_yield->getVal();

  TPaveText *statbox = new TPaveText(0.64,1-pad1->GetTopMargin()-0.03-(4*pavetextsize/(1-padsplit)),1-pad1->GetRightMargin()-0.03,1-pad1->GetTopMargin()-0.03,"BRNDC");
  statbox->SetName("paramBox");
  statbox->SetBorderSize(0);statbox->SetFillColor(kWhite);statbox->SetTextAlign(12);statbox->SetTextSize(pavetextsize/(1-padsplit));statbox->SetTextFont(42);statbox->SetFillStyle(0);
  //passbox->AddText(methlabel);
  //passbox->AddText("matched #mu^{+} probe legs");
  temp.Form("N_{#Xi^{#pm}_{c}} = %.0f #pm %.0f", Nsig,signal_yield->getError());
  if(myconfig->get_particle().Contains("Xic0"))temp.Form("N_{#Xi^{0}_{c}} = %.0f #pm %.0f", Nsig,signal_yield->getError());
  if(myconfig->get_particle().Contains("Omegac"))temp.Form("N_{#Omega^{0}_{c}} = %.0f #pm %.0f", Nsig,signal_yield->getError());
  statbox->AddText(temp);
  temp.Form("N_{Bkg} = %.0f #pm %.0f", Nbkg,background_yield->getError());
  statbox->AddText(temp);
  if(myconfig->get_particle().Contains("Omegac"))temp.Form("#sigma = %.1f #pm %.1f MeV", Xc_sigma1->getVal(),Xc_sigma1->getError());
  else{
    double a_tmp = pow(signalfraction->getVal()*Xc_sigma1->getError(),2);
    double sig2_val = Xc_sigma1->getVal()*widthfactor->getVal();
    double b_tmp = pow(Xc_sigma1->getVal()*signalfraction->getError()-sig2_val*signalfraction->getError(),2);
    double sig2_err = sqrt(pow(widthfactor->getVal()*Xc_sigma1->getError(),2)+pow(widthfactor->getError()*Xc_sigma1->getVal(),2));
    double c_tmp = pow(sig2_err-sig2_err*signalfraction->getVal(),2);
    temp.Form("#bar{#sigma} = %.1f #pm %.1f MeV", signalfraction->getVal()*Xc_sigma1->getVal()+(1-signalfraction->getVal())*sig2_val,sqrt(a_tmp+b_tmp+c_tmp));
  }
  statbox->AddText(temp);
  XcFrame->addObject(statbox) ;

  TPaveText *Selectionbox = new TPaveText(pad1->GetLeftMargin()+0.05,1-pad1->GetTopMargin()-0.06-(1*pavetextsize/(1-padsplit)),0.4,1-pad1->GetTopMargin()-0.06,"BRNDC");
  Selectionbox->SetBorderSize(0);Selectionbox->SetFillColor(kWhite);Selectionbox->SetTextAlign(12);Selectionbox->SetTextSize(pavetextsize/(1-padsplit));Selectionbox->SetTextFont(42);Selectionbox->SetFillStyle(0);
  temp.Form(myconfig->get_cs_at(0));
  Selectionbox->AddText(temp);
  XcFrame->addObject(Selectionbox);

  TPaveText *pt3 = new TPaveText(pad1->GetLeftMargin(),1.002-pad1->GetTopMargin(),pad1->GetLeftMargin()+0.1,1,"BRNDC");
  pt3->SetBorderSize(0);pt3->SetFillColor(kWhite);pt3->SetFillStyle(1001);pt3->SetTextColor(kWhite);
  pt3->AddText(" ");
  XcFrame->addObject(pt3);

  TPaveText *pt4 = new TPaveText(1.002-pad1->GetRightMargin(),pad1->GetBottomMargin(),1,pad1->GetBottomMargin()+0.05,"BRNDC");
  pt4->SetBorderSize(0);pt4->SetFillColor(kWhite);pt4->SetFillStyle(1001);pt4->SetTextColor(kWhite);
  pt4->AddText(" ");
  XcFrame->addObject(pt4);

  //XcFrame->GetXaxis()->SetNdivisions(305) ;
  char yaxislabel[100];
  if(XcFrame->GetMaximum() > 1e+3 && XcFrame->GetMaximum() < 1e+6)sprintf(yaxislabel,"10^{3} Events/%.0f MeV",binsize);
  if(XcFrame->GetMaximum() > 1e+6)sprintf(yaxislabel,"10^{6} Events/%.0f MeV",binsize);
  if(XcFrame->GetMaximum() < 1e+3) sprintf(yaxislabel,"Events/%.0f MeV",binsize);
  //sprintf(yaxislabel,"weighted Events/%.0f MeV",binsize);
  XcFrame->GetYaxis()->SetTitle(yaxislabel) ;
  XcFrame->GetYaxis()->SetTitleSize(labeltextsize/(1-padsplit)) ;
  XcFrame->GetYaxis()->SetLabelSize(labeltextsize/(1-padsplit)) ;
  XcFrame->GetYaxis()->SetTitleOffset(1.2*(1-padsplit)) ;
  XcFrame->GetYaxis()->SetRangeUser(0.1,XcFrame->GetMaximum());
  //Jpsi_Mframeplus->GetYaxis()->SetRangeUser(TMath::Exp(0.9*TMath::Log(Jpsi_Mframeplus->GetMinimum())),TMath::Exp(1.1*TMath::Log(Jpsi_Mframeplus->GetMaximum())));
  XcFrame->Draw() ;

  c1->cd();

  TPad *pad2 = new TPad("pad4","pad4",padsplitmargin,padsplitmargin,1-padsplitmargin,padsplit);
  pad2->SetBorderSize(0);
  pad2->SetLeftMargin(Left_margin);
  pad2->SetRightMargin(Right_margin);
  pad2->SetTopMargin(padsplitmargin);
  pad2->SetBottomMargin(Bigy_margin/padsplit);
  pad2->SetTickx();
  pad2->SetTicky();

  RooPlot* PullFrame = Xc_M.frame(Title(" ")) ;

  Color_t pullcolor = kAzure+3;
  pullhist->SetLineColor(pullcolor);
  pullhist->SetMarkerColor(pullcolor);

  PullFrame->GetXaxis()->SetTitleOffset(1.05);
  PullFrame->GetXaxis()->SetTitleSize(labeltextsize/padsplit);
  PullFrame->GetXaxis()->SetLabelSize(labeltextsize/padsplit);
  PullFrame->GetXaxis()->SetTitle(Xc_M.getTitle(true));
  PullFrame->GetXaxis()->SetTickLength(XcFrame->GetXaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  //h18->GetYaxis()->SetTickLength(Jpsi_Mframeplus->GetYaxis()->GetTickLength()/(padsplit/(1-padsplit)));
  PullFrame->GetYaxis()->SetLabelSize(labeltextsize/padsplit);
  PullFrame->GetYaxis()->SetTitleSize(labeltextsize/padsplit);
  //PullFrame->GetXaxis()->SetNdivisions(305);
  PullFrame->GetYaxis()->SetNdivisions(5);
  PullFrame->GetYaxis()->SetTitleOffset(1.2*padsplit);
  PullFrame->GetYaxis()->SetTitle("Pull ");

  TLine* zeroline = new TLine(myconfig->get_XcMlo(),0,myconfig->get_XcMhi(),0);
  zeroline->SetLineStyle(2);
  zeroline->SetLineColor(1);
  TLine* m3sigma = new TLine(myconfig->get_XcMlo(),-3,myconfig->get_XcMhi(),-3);
  m3sigma->SetLineStyle(2);
  m3sigma->SetLineColor(1);
  TLine* p3sigma = new TLine(myconfig->get_XcMlo(),3,myconfig->get_XcMhi(),3);
  p3sigma->SetLineStyle(2);
  p3sigma->SetLineColor(1);
  PullFrame->addObject(zeroline,"L") ;
  PullFrame->addObject(m3sigma,"L") ;
  PullFrame->addObject(p3sigma,"L") ;
  TPaveText *pt5 = new TPaveText(1.002-pad2->GetRightMargin(),pad2->GetBottomMargin(),1,pad2->GetBottomMargin()+0.05,"BRNDC");
  pt5->SetBorderSize(0);pt5->SetFillColor(kWhite);pt5->SetFillStyle(1001);pt5->SetTextColor(kWhite);
  pt5->AddText(" ");
  PullFrame->addObject(pt5);
  PullFrame->GetYaxis()->SetRangeUser(-4.8,4.8);

  pad2->Draw();
  pad2->cd();
  PullFrame->addPlotable(pullhist,"P") ;
  PullFrame->Draw();
  c1->Update();

  temp = myconfig->get_dumpdir()+"/Fits";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);
  temp += "/" + myconfig->get_particle() + "_simple_cuts_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);

  delete Xc_sigmean;
  delete Xc_sigma1;
  delete widthfactor;
  delete Xc_sigma2;
  delete signalfraction;
  delete signal;
  delete cheba;
  delete chebb;
  delete chebc;
  delete signal_yield;
  delete background_yield;
  delete sig;
  delete bkg;
  delete sum;
  delete c1;
  delete statbox;
  delete Selectionbox;
  delete zeroline;
  delete m3sigma;
  delete p3sigma;
  delete pt3;
  delete pt4;  
  delete pt5;  
  return;
}



bool convcheck(RooFitResult *fitres, vector<RooRealVar*> vars_for_convcheck, vector<RooRealVar*> posis, RooAbsPdf *model, RooDataHist *data, configuration* myconfig){

  /*
Checks if the fitresult of model to data given the parameters vars_for_convcheck and the parameters of special interest posis converged and refits automatically in case not.
The non-convergence considered here has 3 main symptomps: either someting went wrong in the calculation of error matrices/covariances, the fit din't converge at all, or the parameters reached their limits.

The first 2 problems (case A) can be tackled by simply refitting the data with different, more careful options (minimizers) of the fit strategy.
To be more accurate: if something went wrong in the calculation of error matrices/covariances, it is likely to be fixed by a more careful treatment of the errors by the MINOS routine,
which takes rather long and sometimes has problems with fits where the start parameters are not chosen carefully. In our case, the starting parameters are the result from the previous fit, so this should be fine.
If the fit didn't converge at all, new start paramers should be chosen, which is automatically the case here (it is assumed that the non-converged fit is at least a bit better than the model with only the start paraeters).
The third problem (parameters reached limits, case B) can be tackled by defining a new range for these parameters. Currently all other parameters are also constrained to make the fit run faster.

It turned out that an alternating query of the cases gives best performance. I.e. check A first, if it failed, check then B first, if that failed, check again A first.

Currenty 3 levels of checks and refits are implemented (maybe this could be solved more general, i.e. give the number of refits as an option and be able to refit forever and a day ;) )
    */

  double maxedm = 1e-3;
  if ((fitres -> status() != 0  && fitres->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
    cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: this fit did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;
    //Refit I case A: status or uncertainties failed. This can be severe. Try to refit with different options
    if ((fitres -> status() != 0 && fitres->edm() > maxedm)|| errcheck(posis) == false){
      cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options\033[0m" << endl << endl;
      RooFitResult *fitres_refitstat;
      if(myconfig->get_verbosity() > 0)fitres_refitstat = model->fitTo(*data,Save(true),Strategy(2),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit));
      else fitres_refitstat = model->fitTo(*data,Save(true),Strategy(2),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit),PrintLevel(-1),Verbose(false));
      if ((fitres_refitstat -> status() != 0  && fitres_refitstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
        cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 1st attempt to re-fit (statuscheck or bad uncertainties) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;
        //Refit II case AB: parameter reached limit. Will be solved by manipulation of the parameterranges
        if (paratlimit(vars_for_convcheck) == false){
          cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
          setnewparranges(vars_for_convcheck);
          RooFitResult *fitres_refitstatlim;
          if(myconfig->get_verbosity() > 0)fitres_refitstatlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit)) ;
          else fitres_refitstatlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
          if ((fitres_refitstatlim -> status() != 0 && fitres_refitstatlim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (parameter reached limits) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

            //Refit III case ABA: status or uncertainties failed. This can be severe. Try to refit with different options
            if ((fitres_refitstatlim -> status() != 0 && fitres_refitstatlim->edm() > maxedm) || errcheck(posis) == false){
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options\033[0m" << endl << endl;
              RooFitResult *fitres_refitstatlimstat;
              if(myconfig->get_verbosity() > 0)fitres_refitstatlimstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit));
              else fitres_refitstatlimstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refitstatlimstat -> status() != 0 && fitres_refitstatlimstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case ABB: parameter reached limit. Will be solved by manipulation of the parameterranges
            else{
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refitstat2lim;
              if(myconfig->get_verbosity() > 0) fitres_refitstat2lim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit)) ;
              else fitres_refitstat2lim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
              if ((fitres_refitstat2lim -> status() != 0 && fitres_refitstat2lim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
        //Refit II case AA: status or uncertainties failed. This can be severe. Try to refit with different options
        else{
          cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. 2nd try...\033[0m" << endl << endl;
          RooFitResult *fitres_refit2stat;
          if(myconfig->get_verbosity() > 0) fitres_refit2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit));
          else fitres_refit2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
          if ((fitres_refit2stat -> status() != 0  && fitres_refit2stat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (statuscheck or bad uncertainties) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;
            //Refit III case AAB: parameter reached limit. Will be solved by manipulation of the parameterranges
            if (paratlimit(vars_for_convcheck) == false){
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refit2statlim;
              if(myconfig->get_verbosity() > 0) fitres_refit2statlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit)) ;
              else fitres_refit2statlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
              if ((fitres_refit2statlim -> status() != 0 && fitres_refit2statlim->edm() > maxedm)|| errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case AAA: status or uncertainties failed. This can be severe. Try to refit with different options
            else{
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. Last try!\033[0m" << endl << endl;
              RooFitResult *fitres_refit3stat;
              if(myconfig->get_verbosity() > 0) fitres_refit3stat = model->fitTo(*data,Save(true),Minimizer("Minuit","migradimproved"),Strategy(2),NumCPU(cpus_for_fit));
              else fitres_refit3stat = model->fitTo(*data,Save(true),Minimizer("Minuit","migradimproved"),Strategy(2),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refit3stat -> status() != 0 && fitres_refit3stat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
      }
      return true;
    }
    //Refit I case B: parameter reached limit. Will be solved by manipulation of the parameterranges
    else{
      cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
      setnewparranges(vars_for_convcheck);
      RooFitResult *fitres_refitlim;
      if(myconfig->get_verbosity() > 0) fitres_refitlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit)) ;
      else fitres_refitlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
      if ((fitres_refitlim -> status() != 0 && fitres_refitlim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
        cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 1st attempt to re-fit (parameter reached limits) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

        //Refit II case BA: status or uncertainties failed. This can be severe. Try to refit with different options
        if ((fitres_refitlim -> status() != 0 && fitres_refitlim->edm() > maxedm) || errcheck(posis) == false){
          cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. 2nd try...\033[0m" << endl << endl;
          RooFitResult *fitres_refitlimstat;
          if(myconfig->get_verbosity() > 0) fitres_refitlimstat = model->fitTo(*data,Save(true),Strategy(2),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit));
          else fitres_refitlimstat = model->fitTo(*data,Save(true),Strategy(2),InitialHesse(true),Extended(true),NumCPU(cpus_for_fit),PrintLevel(-1),Verbose(false));
          if ((fitres_refitlimstat -> status() != 0 && fitres_refitlimstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (statuscheck or bad uncertainties) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

            //Refit III case BAB: parameter reached limit. Will be solved by manipulation of the parameterranges
            if (paratlimit(vars_for_convcheck) == false){
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refitlimstatlim;
              if(myconfig->get_verbosity() > 0) fitres_refitlimstatlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit)) ;
              else fitres_refitlimstatlim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
              if ((fitres_refitlimstatlim -> status() != 0 && fitres_refitlimstatlim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case BAA: status or uncertainties failed. This can be severe. Try to refit with different options
            else{
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options. Last try!\033[0m" << endl << endl;
              RooFitResult *fitres_refitlim2stat;
              if(myconfig->get_verbosity() > 0) fitres_refitlim2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit));
              else fitres_refitlim2stat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refitlim2stat -> status() != 0 && fitres_refitlim2stat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
        //Refit II case BB: parameter reached limit. Will be solved by manipulation of the parameterranges
        else{
          cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
          setnewparranges(vars_for_convcheck);
          RooFitResult *fitres_refit2lim;
          if(myconfig->get_verbosity() > 0) fitres_refit2lim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit)) ;
          else fitres_refit2lim = model->fitTo(*data,Save(true),Strategy(2),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1)) ;
          if ((fitres_refit2lim -> status() != 0 && fitres_refit2lim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
            cout << endl << endl << "\033[0;31m>>>>>>>>>> ERROR: 2nd attempt to re-fit (parameter reached limits) did not converge! <<<<<<<<<<\033[0m" << endl << endl << endl;

            //Refit III case BBA: status or uncertainties failed. This can be severe. Try to refit with different options
            if ((fitres_refit2lim -> status() != 0 && fitres_refit2lim->edm() > maxedm) || errcheck(posis) == false){
              cout << "\033[0;33mFit result returned bad status or bad uncertainties. Retry the same fit with different options\033[0m" << endl << endl;
              RooFitResult *fitres_refit2limstat;
              if(myconfig->get_verbosity() > 0) fitres_refit2limstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit));
              else fitres_refit2limstat = model->fitTo(*data,Save(true),Minimizer("Minuit2","migrad"),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refit2limstat -> status() != 0 && fitres_refit2limstat->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (statuscheck or bad uncertainties) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
            //Refit III case BBB: parameter reached limit. Will be solved by manipulation of the parameterranges
            else{
              cout << "\033[0;33mAt least one parameter reached its limits. Retry the same fit with different parameter ranges\033[0m" << endl << endl;
              setnewparranges(vars_for_convcheck);
              RooFitResult *fitres_refit3lim;
              if(myconfig->get_verbosity() > 0) fitres_refit3lim = model->fitTo(*data,Save(true),Strategy(1),NumCPU(cpus_for_fit));
              else fitres_refit3lim = model->fitTo(*data,Save(true),Strategy(1),NumCPU(cpus_for_fit),PrintLevel(-1),Warnings(false),Verbose(false),PrintEvalErrors(-1));
              if ((fitres_refit3lim -> status() != 0 && fitres_refit3lim->edm() > maxedm) || errcheck(posis) == false || paratlimit(vars_for_convcheck) == false) {
                cerr << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                cout << endl << endl << "\033[1;31m>>>>>>>>>> ERROR: Last attempt to fit  (parameter reached limits) did not converge! There is no hope for this one. Something has to be done by hand!\033[0m" << endl;
                return false;
              }
              return true;
            }
          }
          return true;
        }
      }
      return true;
    }
  }
  return true;
}

bool errcheck(vector<RooRealVar*> posis){
  cout << endl << "\033[0;34mChecking size of the uncertainties for the parameters of special interest...";

  // Rule of thumb could be Number of free parameters ...
  // Dealing with weighted spectra, so the error from SumW2Errors can get large
  double insaneerrors = 32.0;

  // Figure of merit is the Poissonian standard deviation sqrt(N)
  for (vector<RooRealVar*>::const_iterator siter = posis.begin() ;  posis.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    if( (iter->getError() > insaneerrors*sqrt(iter->getVal()) ) || (iter->getError() < (1/insaneerrors)*sqrt(iter->getVal()) ) ){
      cout << endl << "\033[0;31mParameter: "<< iter->GetName() << " has a unexpectedly high or low uncertainty ( Delta N/sqrt(N) = " << iter->getError()/sqrt(iter->getVal()) << " )\033[0m" << endl;
      return false;
    }
  }
  cout << " OK\033[0m" << endl;
  return true;
}

bool paratlimit(vector<RooRealVar*> vars){

  double convergencecriterion = 0.001;

  cout << endl << "Checking for convergence..." << endl<<endl;

  for (vector<RooRealVar*>::const_iterator siter = vars.begin() ;  vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    if(iter->isConstant())continue;
    cout << "Checking parameter "<< iter->GetName() << " in range [" << iter->getMin() + convergencecriterion*(iter->getMax() - iter->getMin()) << "," << iter->getMax() - convergencecriterion*(iter->getMax() - iter->getMin()) << "]";

    // if the value with its errors reaches the lower edge or the higher edge (+epsilon = convergencecriterion*total range) of the range, the fit has not converged
    if( ( (iter->getVal() - iter->getError()) < iter->getMin() + convergencecriterion*(iter->getMax() - iter->getMin()) ) ||
        ((iter->getVal() + iter->getError()) > iter->getMax() - convergencecriterion*(iter->getMax() - iter->getMin()) ) ){
      cout << endl << "\033[0;31mConvergence check failed at: "<< iter->GetName() << " = " << iter->getVal() << " +- " << iter->getError() << "\033[0m" << endl;
      return false;
    }
    else cout << " ... OK " << endl;
  }
  cout << endl << "\033[0;32mConvergence check passed! " << "\033[0m" << endl << endl;
  return true;
}

void setnewparranges(vector<RooRealVar*> vars){
  //double errfornewrange = 0.02;
  //double errangecheck = 2; // size of error compared to the initial range
  double newrange = 6; // sigma
  double convergencecriterion = 0.001;
  cout << "Assigning new ranges for some parameter..." << endl;
  for (vector<RooRealVar*>::const_iterator siter = vars.begin() ;  vars.end() != siter; ++siter) {
    RooRealVar *iter = *siter;
    //This checks if a parameter reached it's lower limit
    if( (iter->getVal() - iter->getError()) < iter->getMin() + convergencecriterion*(iter->getMax() - iter->getMin())){
      //If this guy has reached it's limit, expand its limits towards lower values
      iter->setRange(iter->getMin()-newrange*iter->getError(),iter->getMax());//be generous ;)
      cout << "New range for parameter " << iter->GetName() << " = [" << iter->getMin() <<"," << iter->getMax() << "]" << endl;
      //continue;
    }
    if ((iter->getVal() + iter->getError()) > iter->getMax() - convergencecriterion*(iter->getMax() - iter->getMin()) ){
      //If this guy has reached it's limit, expand its limits towards higher values
      iter->setRange(iter->getMin(),iter->getMax()+newrange*iter->getError());//be generous ;)
      cout << "New range for parameter " << iter->GetName() << " = [" << iter->getMin() <<"," << iter->getMax() << "]" << endl;
      //continue;
    }
    //This is to constrain the other parameters
    // if the error of the parameter is > 0 and reasonably small compared to the initial parameter range, set a new range
    /*if (iter->getError() > 0 && iter->getError()/iter->getVal() > errfornewrange && iter->getError() < errangecheck*(iter->getMax() - iter->getMin())){
            iter->setRange(iter->getVal()-newrange*iter->getError(),iter->getVal()+newrange*iter->getError());
            cout << "New range for parameter " << iter->GetName() << " = [" << iter->getMin() <<"," << iter->getMax() << "]" << endl;
        }*/
  }
  return;
}
