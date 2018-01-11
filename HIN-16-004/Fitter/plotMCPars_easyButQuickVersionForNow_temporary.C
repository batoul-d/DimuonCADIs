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
#include "Macros/Utilities/resultUtils.h"

#include "Macros/Utilities/bin.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>

using namespace std;

void plotMCPars_easyButQuickVersionForNow_temporary(const char* workDirName,
						    const char* DSTag, //="DATA", // Data Set tag can be: "DATA","MCPSI2SP", "MCJPSIP" ...
						    const char* fitType, // "mass", "ctau"...
						    bool wantPureSMC, // =false,
						    const char* applyCorr, // = "",
						    bool applyJEC // =false
						    )
{
  gStyle->SetOptStat(0);


  TCanvas* c = new TCanvas ("c","",1000,800);

  gSystem->mkdir(Form("Output/%s/%s/%s/MCPars",workDirName, fitType, DSTag));
  ofstream fileOut (Form("Output/%s/%s/%s/MCPars/ParsEvol.dat",workDirName, fitType, DSTag));
  fileOut << "zmin  zmax  n  nerr  alpha  alphaerr"<<endl;

  TH1F* aevol = new TH1F ("aevol",";z(J/#psi);alpha",10, 0, 1);
  TH1F* nevol = new TH1F ("nevol",";z(J/#psi);n", 10, 0, 1);
  //TH1F* aevolnw = new TH1F ("aevolnw",";z(J/#psi);alpha",10, 0, 1);
  //TH1F* nevolnw = new TH1F ("nevolnw",";z(J/#psi);n", 10, 0, 1);

  TH1F* atot = new TH1F ("atot",";z(J/#psi);alpha",10, 0, 0.00001);
  TH1F* ntot = new TH1F ("ntot",";z(J/#psi);n", 10, 0, 0.00001);
  //TH1F* atotnw = new TH1F ("atotnw",";z(J/#psi);alpha",10, 0, 0.00001);
  //TH1F* ntotnw = new TH1F ("ntotnw",";z(J/#psi);n", 10, 0, 0.00001);


  TString treeFileName = Form ("Output/%s/%s/%s/result/tree_allvars.root",workDirName, fitType, DSTag);
  cout << "[INFO] extracting MC parameters from "<<treeFileName<<endl; 
  TFile *f = new TFile(treeFileName);
  if (!f || !f->IsOpen()) {
    cout << "[INFO] tree file not found! creating the result trees."<<endl;
    results2tree(workDirName, DSTag,"", fitType, wantPureSMC, applyCorr, applyJEC);
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

    if (zmin==0 && zmax==1) {
      atot->SetBinContent(atot->FindBin(zmin),alpha);
      atot->SetBinError(atot->FindBin(zmin),alpha_errL);
      ntot->SetBinContent(ntot->FindBin(zmin),n);
      ntot->SetBinError(ntot->FindBin(zmin),n_errL);
      fileOut<< zmin << "  " << zmax << "  " <<n<<"  "<<n_errL<<"  "<<alpha<<"  "<<alpha_errL<<endl;
      ShapeTag = Form("rap%.0f%.0f_%s", ymin*10, ymax*10, jpsiName);
    }    
    else if (zmin>0.2)
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
  ave= new TLine(0,avr,1,avr);
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
  c->SaveAs(Form("Output/%s/%s/%s/MCPars/alphaEvol_%s.png",workDirName, fitType, DSTag, ShapeTag.Data()));
  c->SaveAs(Form("Output/%s/%s/%s/MCPars/alphaEvol_%s.pdf",workDirName, fitType, DSTag, ShapeTag.Data()));
  c->SaveAs(Form("Output/%s/%s/%s/MCPars/alphaEvol_%s.root",workDirName, fitType, DSTag, ShapeTag.Data()));

  avrn=avrn/tot;
  ave = new TLine(0,avrn,1, avrn);

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
  c->SaveAs(Form("Output/%s/%s/%s/MCPars/nEvol_%s.png",workDirName, fitType, DSTag, ShapeTag.Data()));
  c->SaveAs(Form("Output/%s/%s/%s/MCPars/nEvol_%s.pdf",workDirName, fitType, DSTag, ShapeTag.Data()));
  c->SaveAs(Form("Output/%s/%s/%s/MCPars/nEvol_%s.root",workDirName, fitType, DSTag, ShapeTag.Data()));
  f->Close();
  delete f;delete atot; delete ntot; delete aevol; delete nevol; delete c; delete leg; delete tbox; delete ave;
}
