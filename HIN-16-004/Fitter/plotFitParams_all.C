#include <vector>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TString.h"
#include "TLatex.h"
#include "TMathText.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/bin.h"

using namespace std;
void plotBkgOrder(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC);
void getUnfoldingInput(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC, bool statErr);
double readSyst(const char* systfile, double zedmin, double zedmax, double rapmin, double rapmax);
void plotMCPars_easyButQuickVersionForNow_temporary(const char* workDirName,
						    const char* rapRegion,
						    const char* DSTag, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
						    const char* fitType, // "mass", "ctau"...
						    bool wantPureSMC, // =false,
						    const char* applyCorr, // = "",
						    bool applyJEC // =false
						    )
{
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  //double zbins1624 [] = {0.16, 0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};

  int nzbins = 0;
  double zedmin, zedmax;
  if (strcmp(rapRegion,"1624")) {
    nzbins = sizeof(zbins016)/sizeof(double)-1;
    zedmin = zbins016[0];
    zedmax = zbins016[nzbins];
  }
  else {
    nzbins = sizeof(zbins1624)/sizeof(double)-1;
    zedmin = zbins1624[0];
    zedmax = zbins1624[nzbins];
  }
  TCanvas* c = new TCanvas ("c","",1000,800);

  gSystem->mkdir(Form("Output/%s/DataFits_%s/%s/%s/MCPars",workDirName, rapRegion, fitType, DSTag));
  ofstream fileOut (Form("Output/%s/DataFits_%s/%s/%s/MCPars/ParsEvol.txt",workDirName, rapRegion, fitType, DSTag));
  fileOut << "zmin  zmax  n  nerr  alpha  alphaerr"<<endl;

  TH1F* aevol = new TH1F ("aevol",";z(J/#psi);alpha",nzbins, (strcmp(rapRegion,"1624")?zbins016:zbins1624));
  TH1F* nevol = new TH1F ("nevol",";z(J/#psi);n", nzbins, (strcmp(rapRegion,"1624")?zbins016:zbins1624));
  //TH1F* aevolnw = new TH1F ("aevolnw",";z(J/#psi);alpha",10, 0, 1);
  //TH1F* nevolnw = new TH1F ("nevolnw",";z(J/#psi);n", 10, 0, 1);

  TH1F* atot = new TH1F ("atot",";z(J/#psi);alpha",1000, zedmin, zedmin+0.03);
  TH1F* ntot = new TH1F ("ntot",";z(J/#psi);n", 1000, zedmin, zedmin+0.03);
  //TH1F* atotnw = new TH1F ("atotnw",";z(J/#psi);alpha",10, 0, 0.00001);
  //TH1F* ntotnw = new TH1F ("ntotnw",";z(J/#psi);n", 10, 0, 0.00001);


  TString treeFileName = Form ("Output/%s/DataFits_%s/%s/%s/result/tree_allvars.root",workDirName, rapRegion, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl; 
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
    f = new TFile(treeFileName);
    if (!f) return;
  }
  TString ShapeTag = "";
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float n, n_errL,n_errH;
  float alpha, alpha_errL,alpha_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("N_Jpsi_val",&val);
  tr->SetBranchAddress("N_Jpsi_errL",&errL);
  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("n_Jpsi_val",&n);
  tr->SetBranchAddress("n_Jpsi_errL",&n_errL);
  tr->SetBranchAddress("n_Jpsi_errH",&n_errH);
  tr->SetBranchAddress("alpha_Jpsi_val",&alpha);
  tr->SetBranchAddress("alpha_Jpsi_errL",&alpha_errL);
  tr->SetBranchAddress("alpha_Jpsi_errH",&alpha_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);

  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    if (zmin < zedmin+0.02 && zmax == 1){
      atot->SetBinContent(atot->FindBin(zmin+0.0001),alpha);
      atot->SetBinError(atot->FindBin(zmin+0.0001),alpha_errL);
      ntot->SetBinContent(ntot->FindBin(zmin+0.0001),n);
      ntot->SetBinError(ntot->FindBin(zmin+0.0001),n_errL);
      fileOut<< zmin << "  " << zmax << "  " <<n<<"  "<<n_errL<<"  "<<alpha<<"  "<<alpha_errL<<endl;
      ShapeTag = jpsiName;
    }
    else
      { 
	avr=avr+alpha;
	avrn=avrn+n;
	tot++;
	aevol->SetBinContent(aevol->FindBin(zmin+0.02),alpha);
	aevol->SetBinError(aevol->FindBin(zmin+0.02),alpha_errL);
	nevol->SetBinContent(nevol->FindBin(zmin+0.02),n);
	nevol->SetBinError(nevol->FindBin(zmin+0.02),n_errL);

	fileOut<< zmin << "  " << zmax << "  " <<n<<"  "<<n_errL<<"  "<<alpha<<"  "<<alpha_errL<<endl;
      }
  }
  fileOut.close();
  TLegend* leg = NULL;
  TPaveText* tbox = NULL;
  TLine * ave = NULL;
  //TLine * avenw = NULL;

  aevol->SetMinimum(0);
  aevol->SetMaximum(4);
  nevol->SetMinimum(0);
  nevol->SetMaximum(4);
  atot->SetMinimum(0);
  atot->SetMaximum(4);
  ntot->SetMinimum(0);
  ntot->SetMaximum(4);

  avr=avr/tot;
  ave= new TLine(zedmin,avr,1,avr);
  aevol->SetMarkerColor(kRed);
  aevol->SetMarkerStyle(33);
  aevol->SetLineColor(kRed+2);
  aevol->SetOption("E1");
  //aevol->SetFillColor(kRed+2);

  atot->SetMarkerColor(kRed);
  atot->SetMarkerStyle(29);
  atot->SetMarkerSize(2);
  atot->SetLineColor(kRed+2);
  atot->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);

  leg = new TLegend(0.2, 0.6, 0.4, 0.8);
  leg->AddEntry(aevol, "diff.", "lep");
  leg->AddEntry(atot, "int.", "lep");
  leg->AddEntry(ave, "average", "l");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  tbox = new TPaveText(0.2,0.8,0.4,0.9, "BRNDC");
  tbox->AddText(strcmp(DSTag,"MCJPSINOPR")?(strcmp(DSTag,"MCJPSIPR")?"Data":"Prompt J/#psi MC"):"Nonprompt J/#psi MC");
  tbox->AddText(Form("%.1f < |y| < %.1f", ymin, ymax));
  tbox->SetBorderSize(0);
  tbox->SetFillColor(0);
  tbox->SetFillStyle(0);


  c->cd();
  aevol->Draw();
  ave->Draw("same");
  atot->Draw("same E1");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/alphaEvol_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/alphaEvol_%s_%s.pdf",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/alphaEvol_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));

  avrn=avrn/tot;
  ave = new TLine(zedmin,avrn,1, avrn);

  nevol->SetMarkerColor(kRed);
  nevol->SetMarkerStyle(33);
  nevol->SetLineColor(kRed+2);
  nevol->SetOption("E1");
  //nevol->SetFillColor(kRed+2);                                                                                                                                                                       

  ntot->SetMarkerColor(kRed);
  ntot->SetMarkerStyle(29);
  ntot->SetMarkerSize(2);
  ntot->SetLineColor(kRed+2);
  ntot->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);

  c->cd();
  nevol->Draw();
  ntot->Draw("same E1");
  ave->Draw("same");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/nEvol_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/nEvol_%s_%s.pdf",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/MCPars/nEvol_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
  f->Close();
  delete f;delete atot; delete ntot; delete aevol; delete nevol; delete c; delete leg; delete tbox; delete ave;
}


void plotBkgOrder(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC) {
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};

  int nzbins = 0;
  double zedmin, zedmax;
  if (strcmp(rapRegion,"1624")) {
    nzbins = sizeof(zbins016)/sizeof(double)-1;
    zedmin = zbins016[0];
    zedmax = zbins016[nzbins];
  }
  else {
    nzbins = sizeof(zbins1624)/sizeof(double)-1;
    zedmin = zbins1624[0];
    zedmax = zbins1624[nzbins];
  }

  TCanvas* c = new TCanvas ("c","",1000,800);
  TH1F* bkgOrd = NULL;
  if (strcmp(rapRegion,"1624"))
    bkgOrd = new TH1F ("bkgOrd",";z(J/#psi);background order", 7, 0.02, 1);
  else 
    bkgOrd = new TH1F ("bkgOrd",";z(J/#psi);background order", 5, 0, 1);

  string bkgPol [] = {"Uniform", "Chebychev1", "Chebychev2", "Chebychev3", "Chebychev4", "Chebychev5", "Chebychev6"};
  string bkgExp [] = {"Uniform", "ExpChebychev1", "ExpChebychev2", "ExpChebychev3", "ExpChebychev4", "ExpChebychev5", "ExpChebychev6"};
  gSystem->mkdir(Form("Output/%s/DataFits_%s/%s/%s/fitsPars",workDirName, rapRegion, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/DataFits_%s/%s/%s/result/tree_allvars.root",workDirName, rapRegion, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl; 
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
    f = new TFile(treeFileName);
    if (!f) return;
  }
  TString ShapeTag = "";
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float n, n_errL,n_errH;
  float alpha, alpha_errL,alpha_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  tr->SetBranchAddress("N_Jpsi_val",&val);
  tr->SetBranchAddress("N_Jpsi_errL",&errL);
  tr->SetBranchAddress("N_Jpsi_errH",&errH);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);
  //tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("n_Jpsi_val",&n);
  tr->SetBranchAddress("n_Jpsi_errL",&n_errL);
  tr->SetBranchAddress("n_Jpsi_errH",&n_errH);
  tr->SetBranchAddress("alpha_Jpsi_val",&alpha);
  tr->SetBranchAddress("alpha_Jpsi_errL",&alpha_errL);
  tr->SetBranchAddress("alpha_Jpsi_errH",&alpha_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);
  int ord = 0;
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    for (int j = 0; j<7; j++){
      if (bkgName == bkgPol[j]) { ord = j; ShapeTag = "PolChebychev"; break;}
      else if (bkgName == bkgExp [j]) {ord = j; ShapeTag = "ExpChebychev"; break;}
    }
    cout<<"[INFO] z: "<<zmin<<"-"<<zmax<< "bkg: " << bkgName <<endl;
    //if (!(zmin < zedmin+0.02 && zmax == 1)){
      bkgOrd->SetBinContent(bkgOrd->FindBin(zmin+0.001),ord);
      bkgOrd->SetBinError(bkgOrd->FindBin(zmin+0.001), 0.0001);
      //}
  }

    TLatex *  text2 = new TLatex(0.175 ,0.8,strcmp(ShapeTag,"PolChebychev")?"Exp. Chebychev bkg.":"Pol. Chebychev bkg.");
    text2->SetNDC();
    text2->SetTextFont(42);
    text2->SetTextSize(0.05);
    text2->SetLineWidth(2);

    TLatex *  text3 = new TLatex(0.21 ,0.75, Form("%.1f < |y| < %.1f", ymin, ymax));
    text3->SetNDC();
    text3->SetTextFont(42);
    text3->SetTextSize(0.05);
    text3->SetLineWidth(2);

    TLatex *  text = new TLatex(0.75 ,0.8,"CMS");
    text->SetNDC();
    text->SetTextFont(42);
    text->SetTextSize(0.06708595);
    text->SetLineWidth(5);

    TLatex *  text1 = new TLatex(0.7 ,0.72,"Preliminary");
    text1->SetNDC();
    text1->SetTextFont(42);
    text1->SetTextSize(0.05);
    text1->SetLineWidth(2);

    TLatex *  text4 = new TLatex(0.5 ,0.91,"pp 27.39 pb^{-1} (5.02 TeV)");
    text4->SetNDC();
    text4->SetTextFont(42);
    text4->SetTextSize(0.05);
    text4->SetLineWidth(2);


    bkgOrd->GetYaxis()->SetRangeUser(0, 3);
    bkgOrd->SetMarkerColor(kMagenta+3);
    bkgOrd->SetMarkerStyle(33);
    bkgOrd->SetMarkerSize(3);
    bkgOrd->SetLineColor(kMagenta+2);
    //bkgOrd->SetOption("E1");

    c->cd();
    bkgOrd->Draw("EP");
    text->Draw("same");
    text1->Draw("same");
    text2->Draw("same");
    text3->Draw("same");
    text4->Draw("same");
    c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.png",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.pdf",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    c->SaveAs(Form("Output/%s/DataFits_%s/%s/%s/fitsPars/bkgOrder_%s_%s.root",workDirName, rapRegion, fitType, DSTag, ShapeTag.Data(), rapRegion));
    f->Close();
    delete f; delete bkgOrd;
}


void getUnfoldingInput(const char* workDirName, const char* rapRegion, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool applyJEC, bool statErr) {
  gStyle->SetOptStat(0);

  double zbins016 [] = {0.3, 0.44, 0.58, 0.72, 0.86, 1.0};
  //double zbins1624 [] = {0.16, 0.3, 0.44, 0.58, 0.72, 0.86, 1.0};

  //double zbins016 [] = {0.2, 0.4, 0.6, 0.8, 1.0};
  double zbins1624 [] = {0.2, 0.4, 0.6, 0.8, 1.0};


  string binTag = workDirName;
  if (binTag.find("midJtPt")!=std::string::npos) binTag = "midJtPt";
  else if (binTag.find("lowJtPt")!=std::string::npos) binTag = "lowJtPt";
  else if (binTag.find("highJtPt")!=std::string::npos) binTag = "highJtPt";

  double prSyst;
  double nprSyst;
  string systName [] = {"ctauBkg", "ctauErr", "ctauRes", "ctauTrue", "massBkg", "massSig", "AccEffMisMod", "tnpmuidSyst", "AccEffStat", "tnpstaStat", "tnpstaSyst", "tnptrgStat", "tnptrgSyst", "tnpbinned", "tnptrkSyst", "tnpmuidStat"};

  int nSyst = sizeof(systName)/sizeof(systName[0]);

  int nzbins = 0;
  double zedmin, zedmax;
  if (strcmp(rapRegion,"1624")) {
    nzbins = sizeof(zbins016)/sizeof(double)-1;
    zedmin = zbins016[0];
    zedmax = zbins016[nzbins];
  }
  else {
    nzbins = sizeof(zbins1624)/sizeof(double)-1;
    zedmin = zbins1624[0];
    zedmax = zbins1624[nzbins];
  }

  TH1F* prNhist = NULL;//new TH1F ("prNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);
  TH1F* nprNhist = NULL;//new TH1F ("nprNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);

  if (!strcmp(rapRegion,"016")){
    prNhist = new TH1F ("prNhist",";z(J/#psi);N(J/#psi)", 7, 0.02, 1);
    nprNhist = new TH1F ("nprNhist",";z(J/#psi);N(J/#psi)", 7, 0.02, 1);
  }
  else {
    prNhist = new TH1F ("prNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);
    nprNhist = new TH1F ("nprNhist",";z(J/#psi);N(J/#psi)", 5, 0, 1);
  }
  gSystem->mkdir(Form("Output/%s/DataFits_%s/%s/%s/fitsPars",workDirName, rapRegion, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/DataFits_%s/%s/%s/result/tree_allvars.root",workDirName, rapRegion, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl;
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(Form("%s/DataFits_%s", workDirName, rapRegion), DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
    f = new TFile(treeFileName);
    if (!f) return;
  }

  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL=0, bfrac_errH=0;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);
  //tr->SetBranchAddress("N_Jpsi_val",&val);
  //tr->SetBranchAddress("N_Jpsi_errL",&errL);
  //tr->SetBranchAddress("N_Jpsi_errH",&errH);
  tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);                                                                                                                                                 
  tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);                                                                                                                                             
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);
  int ord = 0;
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    //if (zmax < zmin+0.22){
    prSyst = 0;
    nprSyst = 0;
    for (int j=0; j<nSyst ; j++) {
      double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_%s_NJpsi_prompt_PP_%s.csv", binTag.c_str(), rapRegion, systName[i].c_str()), zmin, zmax, ymin, ymax);
      double v2 = readSyst(Form("../Fitter/Systematics/csv/syst_%s_%s_NJpsi_nonprompt_PP_%s.csv", binTag.c_str(), rapRegion, systName[i].c_str()), zmin, zmax, ymin, ymax);
      prSyst=sqrt(pow(prSyst,2)+pow(v1,2));
      nprSyst=sqrt(pow(nprSyst,2)+pow(v2,2));
    }
      prNhist->SetBinContent(prNhist->FindBin(zmin+0.001),val*(1-bfrac));
      if (statErr)
	prNhist->SetBinError(prNhist->FindBin(zmin+0.001), val*(1-bfrac)*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
      else
	prNhist->SetBinError(prNhist->FindBin(zmin+0.001), val*(1-bfrac)*prSyst);

      nprNhist->SetBinContent(nprNhist->FindBin(zmin+0.001),val*bfrac);
      if (statErr)
	nprNhist->SetBinError(nprNhist->FindBin(zmin+0.001), val*bfrac*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
      else 
	nprNhist->SetBinError(nprNhist->FindBin(zmin+0.001),val*bfrac*nprSyst);
      //}
  }
  TFile* fsave = new TFile (Form("Output/%s/DataFits_%s/%s/%s/fitsPars/unfoldingInput_%s_rap%s_%s.root", workDirName, rapRegion, fitType, DSTag, binTag.c_str(), rapRegion, statErr?"statErr":"systErr"),"RECREATE");
  fsave->ls();
  prNhist->Write(Form("prHist_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  nprNhist->Write(Form("nprHist_%s_rap%s_%s", binTag.c_str(), rapRegion, statErr?"statErr":"systErr"));
  fsave->Close();
  delete prNhist; delete nprNhist; delete fsave; delete f;
}


void plotXS(const char* workDirName, const char* DSTag, const char* fitType, bool wantPureSMC, const char* applyCorr, bool statErr) {
  gStyle->SetOptStat(0);
  double ptbins [] = {0, 1.6, 2.4};//{6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35, 50}; //{0, 1.6, 2.4};
  double nbins = sizeof(ptbins)/sizeof(double)-1;

  string systName [] = {"ctauBkg", "ctauErr", "ctauRes", "ctauTrue", "massBkg", "massSig", "AccEffMisMod", "tnpmuidSyst", "AccEffStat", "tnpstaStat", "tnpstaSyst", "tnptrgStat", "tnptrgSyst", "tnpbinned", "tnptrkSyst", "tnpmuidStat"};
  int nSyst = sizeof(systName)/sizeof(systName[0]);

  string lumiTag = workDirName;

  gSystem->mkdir(Form("Output/%s/%s/%s/fitsPars",workDirName, fitType, DSTag));
  TString treeFileName = Form ("Output/%s/%s/%s/result/tree_allvars.root",workDirName, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl;
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(workDirName, DSTag,"", fitType, wantPureSMC, applyCorr, 1);
    f = new TFile(treeFileName);
    if (!f) return;
  }
  ofstream prOut (Form("Output/%s/%s/%s/fitsPars/prXC_vsPt.csv",workDirName, fitType, DSTag));
  ofstream nprOut (Form("Output/%s/%s/%s/fitsPars/nprXC_vsPt.csv",workDirName, fitType, DSTag));

  TH1F* prXS = new TH1F ("prXS", ";pt;XS", nbins, ptbins); prXS->Sumw2();
  TH1F* nprXS = new TH1F ("nprXS", ";pt;XS", nbins, ptbins); nprXS->Sumw2();

  TH1F* prTot = new TH1F ("prTot", "tot XS for pt[6.5-50];whatever;XS", 1, 0, 10); prTot->Sumw2();
  TH1F* nprTot = new TH1F ("nprTot", "tot XS for pt[6.5-50];whatever;XS", 1, 0, 10); nprTot->Sumw2();

  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL=0, bfrac_errH=0;
  float correl=0;
  int ival=-999;
  char collSystem[5];
  char jpsiName[50];
  char bkgName[50];
  float avr=0;
  float avrn=0;
  int tot=0;
  tr->SetBranchAddress("zmin",&zmin);
  tr->SetBranchAddress("zmax",&zmax);
  tr->SetBranchAddress("ptmin",&ptmin);
  tr->SetBranchAddress("ptmax",&ptmax);
  tr->SetBranchAddress("ymin",&ymin);
  tr->SetBranchAddress("ymax",&ymax);
  tr->SetBranchAddress("centmin",&centmin);
  tr->SetBranchAddress("centmax",&centmax);

  tr->SetBranchAddress("N_Jpsi_parLoad_mass",&val);                                                                                                                                                 
  tr->SetBranchAddress("N_Jpsi_parLoad_mass_err",&errL);                                                                                                                                             
  tr->SetBranchAddress("collSystem",collSystem);
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  tr->SetBranchAddress("jpsiName",&jpsiName);
  tr->SetBranchAddress("bkgName",&bkgName);
  int ord = 0;
  int ntr = tr->GetEntries();
  double deltaPt = 1.0;
  double deltaRap = 1.0;
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);

    if (lumiTag.find("muonJSON")!=std::string::npos || lumiTag.find("18025")!=std::string::npos)
      lumi = 28.0e6;
    else if (lumiTag.find("oldLumi")!=std::string::npos) 
      lumi = 25.8e6;
    cout<<"[INFO] lumi = "<<lumi<<endl;
    deltaPt = 1.0;//ptmax - ptmin;
    deltaRap = ymax - ymin;
    double normfact = 1.0;//(1.0/(lumi*deltaPt*deltaRap*1e-3));
    int xBin = prXS->FindBin((ymin+ymax)/2);
    //if ((ptmax-ptmin)<20)
    //{

    double prSyst = 0;
    double nprSyst = 0;
    for (int j=0; j<nSyst ; j++) {
      double v1 = readSyst(Form("../Fitter/Systematics/csv/syst_NoJets_total_NJpsi_prompt_PP_%s.csv", systName[i].c_str()), zmin, zmax, ymin, ymax);
      double v2 = readSyst(Form("../Fitter/Systematics/csv/syst_NoJets_total_NJpsi_nonprompt_PP_%s.csv", systName[i].c_str()), zmin, zmax, ymin, ymax);
      prSyst=sqrt(pow(prSyst,2)+pow(v1,2));
      nprSyst=sqrt(pow(nprSyst,2)+pow(v2,2));
    }
    
    prXS->SetBinContent(xBin, val*(1-bfrac)*normfact);
    if (statErr)
      prXS->SetBinError(xBin, val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    else 
      prXS->SetBinError(xBin, val*(1-bfrac)*prSyst);

    nprXS->SetBinContent(xBin, val*bfrac*normfact);
    if (statErr)
      nprXS->SetBinError(xBin, val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    else nprXS->SetBinError(xBin, val*bfrac*prSyst);

    normfact = 1.0;
    cout<<"pr XS = " << prXS->GetBinContent(xBin)/normfact<<" +- "<<prXS->GetBinError(xBin)/normfact<< " for "<<ptmin<<" < pt < "<<ptmax<<endl;
    cout<<"npr XS = " << nprXS->GetBinContent(xBin)/normfact<<" +- "<<nprXS->GetBinError(xBin)/normfact << " for "<<ptmin<<" < pt < "<<ptmax <<endl;
    prOut<<ptmin<<"  "<<ptmax<<"   "<<prXS->GetBinContent(xBin)/normfact<<"  "<<prXS->GetBinError(xBin)/normfact<<endl;
    nprOut<<ptmin<<"  "<<ptmax<<"   "<<nprXS->GetBinContent(xBin)/normfact<<"  "<<nprXS->GetBinError(xBin)/normfact<<endl;
    //}
    //else if (ptmax>35)
    //{
    //prTot->SetBinContent(prTot->FindBin(5), val*(1-bfrac)*normfact);
    //prTot->SetBinError(prTot->FindBin(5), val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    //nprTot->SetBinContent(nprTot->FindBin(5), val*bfrac*normfact);
    //nprTot->SetBinError(nprTot->FindBin(5), val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)));
    
    //cout<<"pr XS = " << val*(1-bfrac)*normfact <<" +- "<< val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)) << " for "<<ptmin<<" < pt < "<<ptmax<<endl;
    //cout<<"npr XS = " << val*bfrac*normfact <<" +- "<< val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)) << " for "<<ptmin<<" < pt < "<<ptmax<<endl;
    //prOut<<ptmin<<"  "<<ptmax<<"  "<< val*(1-bfrac)*normfact <<"  "<< val*(1-bfrac)*normfact*sqrt(pow(errL/val,2)-2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2)) <<endl;
    //nprOut<<ptmin<<"  "<<ptmax<<"  "<<val*bfrac*normfact<<"  "<< val*bfrac*normfact*sqrt(pow(errL/val,2)+2*correl*errL*bfrac_errL/(val*bfrac)+pow(bfrac_errL/bfrac,2))<<endl;
    //}
  }
  prOut.close();
  nprOut.close();
  prXS->Draw();
  TFile* fsave = new TFile (Form("Output/%s/%s/%s/fitsPars/XSPlot_%s.root", workDirName, fitType, DSTag, statErr?"statErr":"systErr"), "RECREATE");
  fsave->ls();
  prXS->Write("prtotXS");
  nprXS->Write("nprtotXS");
  prTot->Write("prtotIntXS");
  nprTot->Write("nprtotIntXS");
  fsave->Close();
  //delete prNhist; delete nprNhist; delete fsave; delete f;
}//end of plotXS function



void compXSPt(bool isPr, bool isDiff) {
  gStyle->SetOptStat(0);

  double ptbins [] ={6.5, 7.5, 8.5, 9.5, 11, 13, 15, 17.5, 20, 25, 30, 35, 50};
  int nptbins = sizeof(ptbins)/sizeof(double);
  double centbins [] = {0,10};

  double bins [15];
  int nbins;
  if (isDiff){
    nbins = sizeof(ptbins)/sizeof(double);
    for (int i =0; i<nbins; i++)
      bins[i] = ptbins [i];
  }
  else {
    nbins = sizeof(centbins)/sizeof(double);
    for (int i =0; i<nbins; i++)
      bins[i] = centbins [i];
  }
  nbins = nbins-1;

  TH1F* corrHist = new TH1F("corrHist","", nbins, bins); corrHist->Sumw2();
  TH1F* refHist = new TH1F("refHist","", nbins, bins); refHist->Sumw2();
  TH1F* goldHist = NULL;

  ifstream accEff;
  accEff.open(Form("../MyEfficiency/RatioStudy/AccEff16025_%s.csv", isDiff?"pt_rap0024":"cent_rap0024"));

  Float_t x,y,z,w;
  Int_t nlines = 0;

  while (1)
    {
      accEff >> x >> y >> z >> w;
      if (!accEff.good()) break;
      if (isPr) corrHist->SetBinContent(corrHist->FindBin((x+y)*0.5), z);
      else corrHist->SetBinContent(corrHist->FindBin((x+y)*0.5), w);
      corrHist->SetBinError(corrHist->FindBin((x+y)*0.5), 0.001);
      nlines++;
    }
  accEff.close();

  double refVal = 0;
  double refErr = 0;

  ifstream ref;
  ref.open(Form("Input/%sXS16025_pt.csv", isPr?"pr":"npr"));
  nlines = 0;
  while (1)
    {
      ref >> x >> y >> z >> w;
      if (!ref.good()) break;
      if (isDiff)
	{
	  refHist->SetBinContent(refHist->FindBin((x+y)*0.5), z);
	  refHist->SetBinError(refHist->FindBin((x+y)*0.5), w);
	}
      else { 
	refVal = refVal+z; 
	refErr = refErr+(w/z)*(w/z);
      }
      nlines++;
    }
  ref.close();
  //cout << "[INFO] intVal = "<<refVal<<" +- "<<refVal*sqrt(refErr)<<" intVal/(dpt) = "<<refVal/(50-6.5)<<" +- "<<refVal*sqrt(refErr)<<endl;
  refVal = refVal/(50-6.5);
  refErr = refVal*sqrt(refErr);
  if (!isDiff)
    refHist->SetBinContent(refHist->FindBin(5), refVal);

  TFile* goldFile = TFile::Open("Output/DataFits_totalN_y024/ctauMass/DATA/fitsPars/XSPlot.root");

  goldHist = (TH1F*) goldFile->Get(Form("%stot%sXS",isPr?"pr":"npr", isDiff?"":"Int"));

  goldHist->Divide(corrHist);

  refHist->SetLineColor(kRed+2);
  goldHist->SetLineColor(kBlue+2);

  refHist->SetMarkerColor(kRed);
  goldHist->SetMarkerColor(kBlue);

  refHist->SetMarkerStyle(4);
  goldHist->SetMarkerStyle(4);

  refHist->SetMarkerSize(1);
  goldHist->SetMarkerSize(1);

  TLegend* leg = new TLegend (0.5, 0.7, 0.9, 0.8);
  leg->AddEntry(refHist, "HIN-16-025 results", "lep");
  leg->AddEntry(goldHist, "HIN-18-012 (evt-by-evt)", "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TH1F* axisHist = new TH1F("axisHist","", nbins, bins);
  axisHist->SetTitle(Form("%s XS comparison", isPr?"prompt":"nonprompt"));
  axisHist->GetYaxis()->SetLimits(0,(goldHist->GetMaximum()<1)?1:(1.2*goldHist->GetMaximum()));
  goldHist->SetTitle(Form("%s XS comparison", isPr?"prompt":"nonprompt"));
  refHist->SetTitle(Form("%s XS comparison", isPr?"prompt":"nonprompt"));

  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  if (!isDiff) axisHist->Draw();
  else refHist->Draw();
  refHist->Draw("same");
  goldHist->Draw("same");
  leg->Draw("same");

  gSystem->mkdir("Output/XSComparison");
  c->SaveAs(Form("Output/XSComparison/%s%scomparison.pdf", isPr?"pr":"npr", isDiff?"Diff":"Int"));

  if (!isDiff) cout<<"[INFO] xs(tot18-012) = "<<goldHist->GetBinContent(1)<<endl;//" xs(18012) = "<<muonHist->GetBinContent(muonHist->FindBin(6.5))<<endl;
  goldFile->Close();

}//end of compXS function

void plotXStot(bool isPr) {
  gStyle->SetOptStat(0);

  double ybins [] = {0, 1.6, 2.4};

  TFile* totFile = TFile::Open("Output/DataFits/DataFits_NoJets/DataFits_total/ctauMass/DATA/fitsPars/XSPlot_statErr.root");
  TFile* totSystFile = TFile::Open("Output/DataFits/DataFits_NoJets/DataFits_total/ctauMass/DATA/fitsPars/XSPlot_systErr.root");
  TFile* mcFile = TFile::Open(Form("Output/MCResults/mcResult_%s_underflowOff.root", isPr?"prompt":"nonprompt"));

  TH1F* totHist = (TH1F*) totFile->Get(Form("%stotXS",isPr?"pr":"npr"));
  TH1F* totSystHist = (TH1F*) totSystFile->Get(Form("%stotXS",isPr?"pr":"npr"));

  TH1F* midMCHist = (TH1F*) mcFile->Get("zDist_mid");
  TH1F* midMCTotHist =(TH1F*) mcFile->Get("Ntot_mid");
  TH1F* fwdMCHist = (TH1F*) mcFile->Get("zDist_fwd");
  TH1F* fwdMCTotHist =(TH1F*) mcFile->Get("Ntot_fwd");
  TH1F* mcHist = new TH1F("mcHist", "", 2, ybins);

  TFile* jetMidFile = TFile::Open(Form("Output/unfData/results/unfResult_%s_mid_wStats.root",isPr?"prompt":"nonprompt"));
  TFile* jetMidSystFile = TFile::Open(Form("Output/unfData/results/systErrs_%s_Mid_Total.root",isPr?"Prompt":"NonPrompt"));
  TFile* jetFwdFile = TFile::Open(Form("Output/unfData/results/unfResult_%s_fwd_wStats.root",isPr?"prompt":"nonprompt"));
  TFile* jetFwdSystFile =TFile::Open(Form("Output/unfData/results/systErrs_%s_Fwd_Total.root",isPr?"Prompt":"NonPrompt"));

  TH1F* midHist = (TH1F*) jetMidFile->Get("zUnf"); 
  TH1F* fwdHist = (TH1F*)jetFwdFile->Get("zUnf");
  TH1F* midSystHist = (TH1F*) jetMidSystFile->Get("totalSyst");
  TH1F* fwdSystHist = (TH1F*)jetFwdSystFile->Get("totalSyst");

  double valMid =0;
  double errMid =0;
  double systMid =0;
  double valFwd =0;
  double errFwd =0;
  double systFwd =0;
  double valMCMid =0;
  double valMCFwd =0;
  double errMCMid =0;
  double errMCFwd =0;

  for (int j =0; j<4 ; j++){
    valFwd =valFwd + fwdHist->GetBinContent(fwdHist->FindBin(0.2+j*(0.2)));
    errFwd = sqrt(errFwd*errFwd + fwdHist->GetBinError(fwdHist->FindBin(0.2+j*0.2))*fwdHist->GetBinError(fwdHist->FindBin(0.2+j*0.2)));
    systFwd = sqrt(systFwd*systFwd + fwdSystHist->GetBinError(fwdSystHist->FindBin(0.2+j*0.2))*fwdSystHist->GetBinError(fwdSystHist->FindBin(0.2+j*0.2)));
    valMCFwd= valMCFwd + fwdMCHist->GetBinContent(fwdMCHist->FindBin(0.2+j*(0.2)));
    errMCFwd= sqrt(errMCFwd*errMCFwd + fwdMCHist->GetBinError(fwdMCHist->FindBin(0.2+j*0.2))*fwdMCHist->GetBinError(fwdMCHist->FindBin(0.2+j*0.2)));
  }
  for (int i =0; i<5 ; i++){
    valMid = valMid+midHist->GetBinContent(midHist->FindBin(0.44+(i*0.14)));
    errMid = sqrt(errMid*errMid + midHist->GetBinError(midHist->FindBin(0.44+(i*0.16)))*midHist->GetBinError(midHist->FindBin(0.46+(i*0.16))));
    systMid = sqrt(systMid*systMid + midSystHist->GetBinError(midSystHist->FindBin(0.44+(i*0.16)))*midSystHist->GetBinError(midSystHist->FindBin(0.46+(i*0.16))));
    valMCMid= valMCMid + midMCHist->GetBinContent(midMCHist->FindBin(0.44+(i*0.14)));
    errMCMid= sqrt(errMCMid*errMCMid + midMCHist->GetBinError(midMCHist->FindBin(0.44+(i*0.16)))*midMCHist->GetBinError(midMCHist->FindBin(0.44+(i*0.14))));
  }

  int bin =0;
  bin = totHist->FindBin(1);
  cout<<"totXS in mid rapidity = "<<totHist->GetBinContent(bin)/(27.39e3*1.6)<<", statErr =  "<<totHist->GetBinError(bin)/(27.39e3*1.6)<<", systErr "<<totSystHist->GetBinError(bin)/(27.39e3*1.6)<<endl;
  errMid = errMid/valMid;
  systMid = systMid/valMid;
  valMid = valMid/totHist->GetBinContent(bin);
  errMid = sqrt(errMid*errMid + totHist->GetBinError(bin)*totHist->GetBinError(bin)/(totHist->GetBinContent(bin)*totHist->GetBinContent(bin)));
  errMid = errMid*valMid;
  systMid = sqrt(systMid*systMid + totSystHist->GetBinError(bin)*totSystHist->GetBinError(bin)/(totSystHist->GetBinContent(bin)*totSystHist->GetBinContent(bin)));
  systMid = systMid*valMid;
  totHist->SetBinContent(bin, valMid);
  totHist->SetBinError(bin, errMid);
  totSystHist->SetBinContent(bin, valMid);
  totSystHist->SetBinError(bin, systMid);


  errMCMid = errMCMid/valMCMid;
  valMCMid = valMCMid/midMCTotHist->GetBinContent(midMCTotHist->FindBin(0.0001));
  errMCMid = sqrt(errMCMid*errMCMid + midMCTotHist->GetBinError(midMCTotHist->FindBin(0.0001))*midMCTotHist->GetBinError(midMCTotHist->FindBin(0.0001))/(midMCTotHist->GetBinContent(midMCTotHist->FindBin(0.0001))*midMCTotHist->GetBinContent(midMCTotHist->FindBin(0.0001))));
  errMCMid=errMCMid*valMCMid;
  mcHist->SetBinContent(bin, valMCMid);
  mcHist->SetBinError(bin, errMCMid);

  bin = totHist->FindBin(2);
  cout<<"totXS in fwd rapidity = "<<totHist->GetBinContent(bin)/(27.39e3*1.6)<<", statErr =  "<<totHist->GetBinError(bin)/(27.39e3*1.6)<<", systErr "<<totSystHist->GetBinError(bin)/(27.39e3*1.6)<<endl;
  errFwd = errFwd/valFwd;
  systFwd = systFwd/valFwd;
  valFwd = valFwd/totHist->GetBinContent(bin);
  errFwd = sqrt(errFwd*errFwd + totHist->GetBinError(bin)*totHist->GetBinError(bin)/(totHist->GetBinContent(bin)*totHist->GetBinContent(bin)));
  errFwd = errFwd*valFwd;
  systFwd = sqrt(systFwd*systFwd + totSystHist->GetBinError(bin)*totSystHist->GetBinError(bin)/(totSystHist->GetBinContent(bin)*totSystHist->GetBinContent(bin)));
  systFwd = systFwd*valFwd;
  totHist->SetBinContent(bin, valFwd);
  totHist->SetBinError(bin, errFwd);
  totSystHist->SetBinContent(bin, valFwd);
  totSystHist->SetBinError(bin, systFwd);

  errMCFwd = errMCFwd/valMCFwd;
  valMCFwd = valMCFwd/fwdMCTotHist->GetBinContent(fwdMCTotHist->FindBin(0.0001));
  errMCFwd = sqrt(errMCFwd*errMCFwd + fwdMCTotHist->GetBinError(fwdMCTotHist->FindBin(0.0001))*fwdMCTotHist->GetBinError(fwdMCTotHist->FindBin(0.0001))/(fwdMCTotHist->GetBinContent(fwdMCTotHist->FindBin(0.0001))*fwdMCTotHist->GetBinContent(fwdMCTotHist->FindBin(0.0001))));
  errMCFwd=errMCFwd*valMCFwd;
  mcHist->SetBinContent(bin, valMCFwd);
  mcHist->SetBinError(bin, errMCFwd);

  gSystem->mkdir("Output/XSComparison");
  TFile *fsave = new TFile(Form("Output/XSComparison/%sXSRatio.root", isPr?"pr":"npr"), "RECREATE");
  totHist->Write("XS");
  totSystHist->Write("XS_syst");
  mcHist->Write("XS_mc");
  fsave->Close();
}//end of compXStot

void drawXStot() {
  gStyle->SetOptStat(0);

  TFile *prf = TFile::Open("Output/XSComparison/prXSRatio.root"); 
  TFile *nprf = TFile::Open("Output/XSComparison/nprXSRatio.root");
  TH1F *prHist = (TH1F*) prf->Get("XS");
  TH1F *nprHist = (TH1F*) nprf->Get("XS");
  TH1F *prSystHist = (TH1F*) prf->Get("XS_syst");
  TH1F *nprSystHist = (TH1F*) nprf->Get("XS_syst");
  TH1F *prMCHist = (TH1F*) prf->Get("XS_mc");
  TH1F *nprMCHist = (TH1F*) nprf->Get("XS_mc");

  prHist->Scale(100);
  nprHist->Scale(100);
  prSystHist->Scale(100);
  nprSystHist->Scale(100);
  prMCHist->Scale(100);
  nprMCHist->Scale(100);

  prHist->SetMarkerColor(kMagenta+2);
  prHist->SetMarkerStyle(kFullDiamond);
  prHist->SetMarkerSize(2);
  prHist->SetLineColor(kMagenta+2);
  prHist->SetLineWidth(2);

  nprHist->SetMarkerColor(kCyan+2);
  nprHist->SetMarkerStyle(kFullDiamond);
  nprHist->SetMarkerSize(2);
  nprHist->SetLineColor(kCyan+2);
  nprHist->SetLineWidth(2);

  prSystHist->SetMarkerColor(kMagenta+2);
  prSystHist->SetMarkerStyle(kFullDiamond);
  prSystHist->SetMarkerSize(0);
  prSystHist->SetLineColor(kMagenta+2);
  prSystHist->SetFillColorAlpha(kMagenta-5, 0.5);
  prSystHist->SetLineWidth(2);

  nprSystHist->SetMarkerColor(kCyan+2);
  nprSystHist->SetMarkerStyle(kFullDiamond);
  nprSystHist->SetMarkerSize(0);
  nprSystHist->SetLineColor(kCyan+2);
  nprSystHist->SetFillColorAlpha(kCyan-5, 0.5);
  nprSystHist->SetLineWidth(2);

  prMCHist->SetMarkerColor(kMagenta+2);
  prMCHist->SetMarkerStyle(kOpenDiamond);
  prMCHist->SetMarkerSize(2);
  prMCHist->SetLineColor(kMagenta+2);
  prMCHist->SetLineWidth(2);

  nprMCHist->SetMarkerColor(kCyan+2);
  nprMCHist->SetMarkerStyle(kOpenDiamond);
  nprMCHist->SetMarkerSize(2);
  nprMCHist->SetLineColor(kCyan+2);
  nprMCHist->SetLineWidth(2);

  double ybins [] = {0, 1.6, 2.4};
  TH1F *axisHist = new TH1F("axisHist","", 2, ybins);
  axisHist->GetYaxis()->SetRangeUser(0.01, 500);
  axisHist->GetYaxis()->SetTitle("#sigma(J/#psi-Jet)/#sigma_{tot}(J/#psi) (%)");
  axisHist->GetXaxis()->SetTitle("|y(J/#psi)|");
  axisHist->GetXaxis()->SetNdivisions(505);
  axisHist->GetXaxis()->CenterTitle(true);
  axisHist->GetYaxis()->CenterTitle(true);

  TLegend* leg = new TLegend(0.62,0.62,0.88,0.80);

  leg->AddEntry(prHist, "prompt data", "lp");
  leg->AddEntry(nprHist, "nonprompt data", "lp");
  leg->AddEntry(prMCHist, "prompt MC", "lp");
  leg->AddEntry(nprMCHist, "nonprompt MC", "lp");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  TLatex *  text2 = new TLatex(0.160 ,0.72, "6.5 < p_{T} < 35 GeV for |y| < 1.6");
  text2->SetNDC();
  text2->SetTextFont(42);
  text2->SetTextSize(0.03);
  text2->SetLineWidth(2);
  
  TLatex *  text3 = new TLatex(0.160 ,0.67, "3 < p_{T} < 35 GeV for 1.6 < |y| < 2.4");
  text3->SetNDC();
  text3->SetTextFont(42);
  text3->SetTextSize(0.03);
  text3->SetLineWidth(2);
  
  TLatex *  text = new TLatex(0.18 ,0.82,"CMS");
  text->SetNDC();
  text->SetTextFont(42);
  text->SetTextSize(0.06);
  text->SetLineWidth(8);
  
  TLatex *  text1 = new TLatex(0.16 ,0.77,"Preliminary");
  text1->SetNDC();
  text1->SetTextFont(42);
  text1->SetTextSize(0.04);
  text1->SetLineWidth(2);
  
  TLatex *  text4 = new TLatex(0.5 ,0.91,"pp 27.39 pb^{-1} (5.02 TeV)");
  text4->SetNDC();
  text4->SetTextFont(42);
  text4->SetTextSize(0.04);
  text4->SetLineWidth(2);
  

  TCanvas* c = new TCanvas("c", "", 1000, 1000);
  axisHist->Draw();
  prSystHist->Draw("E2 same");
  prHist->Draw("E1 same");
  prMCHist->Draw("E1 same");
  nprSystHist->Draw("E2 same");
  nprHist->Draw("E1 same");
  nprMCHist->Draw("E1 same");
  leg->Draw("same");
  text2->Draw("same");
  text3->Draw("same");
  text->Draw("same");
  text1->Draw("same");
  text4->Draw("same");
  //c->SaveAs("Output/XSComparison/totXSPlot.pdf");
  c->SetLogy();
  c->SaveAs("Output/XSComparison/totXSPlot_logScale.pdf");
}

double readSyst(const char* systfile, double zedmin, double zedmax, double rapmin, double rapmax) {
  double ans;
  ans = 0;
  ifstream file(systfile);
  if (!(file.good())) return ans;

  string systname; getline(file,systname);

  string line;
  double zmin=0, zmax=0, ymin=0, ymax=0, ptmin=0, ptmax=0, centmin=0, centmax=0, value=0;

  while (file.good()) {
    getline(file,line);
    if (line.size()==0) break;
    TString tline(line.c_str());
    TString t; Int_t from = 0, cnt=0;
    while (tline.Tokenize(t, from , ",")) {
      t.Strip(TString::kBoth,' ');
      value = atof(t.Data());
      if (cnt==0) ymin = atof(t.Data());
      else if (cnt==1) ymax = value;
      else if (cnt==2) ptmin = value;
      else if (cnt==3) ptmax = value;
      else if (cnt==4) zmin = value;
      else if (cnt==5) zmax = value;
      else if (cnt==6) centmin = value;
      else if (cnt==7) centmax = value;
      else if (cnt>8) {
	cout << "Warning, too many fields, I'll take the last one." << endl;
	continue;
      }
      cnt++;
    }
    //if (!(zmin == 0.4 && zmax == 1.0) && !(zmin == 0.2 && zmax == 1.0)) ///// to not take the integrated results
    if (zmin<zedmin+0.001 && zmin>zedmin-0.001 && zmax<zedmax+0.001 && zmax>zedmax-0.001 && ymin<rapmin+0.001 && ymin>rapmin-0.001 && ymax<rapmax+0.001 && ymax>rapmax-0.001)
      //ans.push_back(value);
      ans = value;
  }
  file.close();
  if (ans == 0) cout <<"[WARNING] systematic = 0;"<<endl;
  if (ans >10) cout <<"[WARNING] huge systematic in "<<systfile<<endl;
  return ans;
}

