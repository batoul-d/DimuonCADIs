#include <vector>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TArrow.h"
#include "TString.h"
#include "TLatex.h"
#include "TMathText.h"

using namespace std;

void plotMCPars_easyButQuickVersionForNow_temporary()
{
  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas ("c","",1000,800);
  TFile* fsave = new TFile ("Output/MCPars/MCPars.root", "RECREATE");

  ///////// with JEC ///// alpha and n free
  float  nVal_npr_w_016 [] = {1.176, 0.716, 1.142, 0.93, 1.191, 1.239, 1.4, 1.334, 1.533, 2.765};
  float  nErr_npr_w_016 [] = {0.325, 1.8, 59, 0.201, 0.089, 0.091, 0.117, 0.15, 0.241, 1.577};
  float  aVal_npr_w_016 [] = {2.041, 1.748, 1.924, 2.125, 2.012, 2.023, 1.978, 2.072, 2.075, 1.931};
  float  aErr_npr_w_016 [] = {0.214, 2.14, 42.99, 0.156, 0.058, 0.057, 0.066, 0.086, 0.118, 0.4};

  float  nVal_npr_nw_016 [] = {1.23, 1.81, 1.004, 0.964, 1.199, 1.24, 1.418, 1.355, 1.554, 2.702};
  float  nErr_npr_nw_016 [] = {0.041, 2.39, 0.195, 0.094, 0.076, 0.081, 0.107, 0.135, 0.222, 1.307};
  float  aVal_npr_nw_016 [] = {2.035, 1.228, 1.948, 2.104, 2.008, 2.026, 1.971, 2.062, 2.07, 1.956};
  float  aErr_npr_nw_016 [] = {0.025, 0.842, 0.142, 0.07, 0.049, 0.05, 0.059, 0.076, 0.107, 0.261};

  float  nVal_pr_w_016 [] = {};
  float  nErr_pr_w_016 [] = {};
  float  aVal_pr_w_016 [] = {};
  float  aErr_pr_w_016 [] = {};

  float  nVal_pr_nw_016 [] = {};
  float  nErr_pr_nw_016 [] = {};
  float  aVal_pr_nw_016 [] = {};
  float  aErr_pr_nw_016 [] = {};

  // forward rap                                                                                                                                                                                            
  float  nVal_npr_w_1624 [] = {};
  float  nErr_npr_w_1624 [] = {};
  float  aVal_npr_w_1624 [] = {};
  float  aErr_npr_w_1624 [] = {};
  
  float  nVal_npr_nw_1624 [] = {};
  float  nErr_npr_nw_1624 [] = {};
  float  aVal_npr_nw_1624 [] = {};
  float  aErr_npr_nw_1624 [] = {};
  
  float  nVal_pr_w_1624 [] = {};
  float  nErr_pr_w_1624 [] = {};
  float  aVal_pr_w_1624 [] = {};
  float  aErr_pr_w_1624 [] = {};
  
  float  nVal_pr_nw_1624 [] = {};
  float  nErr_pr_nw_1624 [] = {};
  float  aVal_pr_nw_1624 [] = {};
  float  aErr_pr_nw_1624 [] = {};
  
  /////// without JEC ////// alpha and n free
  // mid rap 
  //float  nVal_npr_w_016 [] = {1.23, 1.13, 1.5, 1.03, 1.183, 1.365, 1.35, 1.8, 2.56, 4.99};
  //float  nErr_npr_w_016 [] = {0.26, 8.5, 0.0083, 0.15, 0.079, 0.088, 0.1, 0.18, 0.57, 5.62};
  //float  aVal_npr_w_016 [] = {2.022, 1.14, 1.7, 2.06, 2.0421, 1.97, 2.01, 1.92, 1.92, 1.85};
  //float  aErr_npr_w_016 [] = {0.166, 5.3, 1.3, 0.11, 0.052, 0.05, 0.06, 0.076, 0.16, 0.3};

  //float  nVal_npr_nw_016 [] = {1.27,3,1.04,1.06,1.16,1.37,1.34,1.83,2.57,5};
  //float  nErr_npr_nw_016 [] = {0.03,2.93,0.19,0.084,0.07,0.077,0.09,0.16,0.47,4.81};
  //float  aVal_npr_nw_016 [] = {2.022,0.57,1.94,2.06,2.05,1.96,2.02,1.91,1.93,1.83};
  //float  aErr_npr_nw_016 [] = {0.022,0.34,0.13,0.059,0.045,0.044,0.05,0.07,0.13,0.28};

  //float  nVal_pr_w_016 [] = {1.07, 2.78, 1, 0.63, 0.41, 1.25, 1.1, 1.4, 1.86, 5};
  //float  nErr_pr_w_016 [] = {0.37, 149.6, 0.37, 1.94, 0.3, 0.46, 0.16, 0.24, 0.46, 15.6};
  //float  aVal_pr_w_016 [] = {2.12, 2.99, 2.23, 2.4, 2.46, 1.86, 1.91, 1.89, 2, 1.64};
  //float  aErr_pr_w_016 [] = {0.26, 122.79, 1.66, 1.73, 0.344, 0.31, 0.1, 0.13, 0.2, 0.15};

  //float  nVal_pr_nw_016 [] = {1.19, 2.5, 0.14, 0.72, 0.43, 1.23, 1.1, 1.37, 1.82, 5};
  //float  nErr_pr_nw_016 [] = {0.11, 1.55, 0.15, 0.35, 0.23, 0.3, 0.15, 0.24, 0.32, 3.67};
  //float  aVal_pr_nw_016 [] = {2.06, 2.94, 2.99, 2.37, 2.45, 1.89, 1.91, 1.91, 2.01, 1.62};
  //float  aErr_pr_nw_016 [] = {0.07, 1.96, 2.13, 0.28, 0.25, 0.19, 0.1, 0.13, 0.13, 0.15};
  
  // forward rap
  //float  nVal_npr_w_1624 [] = {0.979, 0.0001, 1.153, 0.922, 0.807, 0.559, 1.712, 2.007, 3.00};
  //float  nErr_npr_w_1624 [] = {0.104, 0.185, 0.747, 0.217, 0.1911, 0.222, 0.42, 0.913, 10.454};
  //float  aVal_npr_w_1624 [] = {2.266, 2.693, 1.993, 2.206, 2.334, 2.648, 1.998, 2.060, 1.905};
  //float  aErr_npr_w_1624 [] = {0.062, 0.394, 1.51, 0.134, 0.123, 0.162, 0.163, 0.294, 1.703};

  //float  nVal_npr_nw_1624 [] = {0.999, 0.001, 1.288, 0.905, 0.774, 0.543, 1.679, 2.166, 3};
  //float  nErr_npr_nw_1624 [] = {0.124, 2.11, 0.709, 0.217, 0.169, 0.230, 0.335, 0.533, 0.286};
  //float  aVal_npr_nw_1624 [] = {2.280, 2.69, 1.977, 2.231, 2.372, 2.663, 2.014, 1.998, 1.885};
  //float  aErr_npr_nw_1624 [] = {0.073, 0.39, 0.384, 0.139, 0.111, 0.176, 0.132, 0.149, 0.06};

  //float  nVal_pr_w_1624 [] = {1.606, 2.99, 3, 0.476, 2.902, 0.076, 0.22, 2.61, 0.56};
  //float  nErr_pr_w_1624 [] = {0.633, 1.60, 96, 0.08, 2.206, 53.6, 0.419, 1.66, 2.18};
  //float  aVal_pr_w_1624 [] = {2.016, 3, 1.25, 2.417, 1.418, 3, 2.853, 1.699, 3};
  //float  aErr_pr_w_1624 [] = {0.273, 0.2, 4.86, 0.0647, 0.519, 68.6, 0.507, 0.518, 3.5};

  //float  nVal_pr_nw_1624 [] = {1.545, 2.994, 2.999, 0.445, 2.983, 0.073, 0.191, 2.789, 0.58};
  //float  nErr_pr_nw_1624 [] = {0.486, 2.81, 2.81, 0.025, 0.495, 0.668, 2.711, 0.393, 2.653, 0.59};
  //float  aVal_pr_nw_1624 [] = {2.047, 3, 1.235, 2.36, 1.418, 2.999, 2.88, 1.674, 2.985};
  //float  aErr_pr_nw_1624 [] = {0.215, 1.67, 0.024, 0.41, 0.048, 2.5, 0.431, 0.286, 2.74};


  TH1F* aevol = new TH1F ("aevol",";z(J/#psi);alpha",10, 0, 1);
  TH1F* nevol = new TH1F ("nevol",";z(J/#psi);n", 10, 0, 1);

  TH1F* aevolnw = new TH1F ("aevolnw",";z(J/#psi);alpha",10, 0, 1);
  TH1F* nevolnw = new TH1F ("nevolnw",";z(J/#psi);n", 10, 0, 1);

  TH1F* atot = new TH1F ("atot",";z(J/#psi);alpha",10, 0, 0.00001);
  TH1F* ntot = new TH1F ("ntot",";z(J/#psi);n", 10, 0, 0.00001);
  TH1F* atotnw = new TH1F ("atotnw",";z(J/#psi);alpha",10, 0, 0.00001);
  TH1F* ntotnw = new TH1F ("ntotnw",";z(J/#psi);n", 10, 0, 0.00001);

  TLegend* leg = NULL;
  TPaveText* tbox = NULL;
  float avr = nVal_npr_w_016[2]+nVal_npr_w_016[3]+nVal_npr_w_016[4]+nVal_npr_w_016[5]+nVal_npr_w_016[6]+nVal_npr_w_016[7]+nVal_npr_w_016[8];
  avr = avr/7;
  TLine * ave = NULL;
  TLine * avenw = NULL;

  aevol->SetMinimum(0);
  aevol->SetMaximum(4);
  nevol->SetMinimum(0);
  nevol->SetMaximum(4);
  aevolnw->SetMinimum(0);
  aevolnw->SetMaximum(4);
  nevolnw->SetMinimum(0);
  nevolnw->SetMaximum(4);
  atot->SetMinimum(0);
  atot->SetMaximum(4);
  ntot->SetMinimum(0);
  ntot->SetMaximum(4);
  atotnw->SetMinimum(0);
  atotnw->SetMaximum(4);
  ntotnw->SetMinimum(0);
  ntotnw->SetMaximum(4);

  atot->SetBinContent(atot->FindBin(0), aVal_npr_w_016[0]);
  atot->SetBinError(atot->FindBin(0),aErr_npr_w_016[0]);

  ntot->SetBinContent(ntot->FindBin(0), nVal_npr_w_016[0]);
  ntot->SetBinError(ntot->FindBin(0),nErr_npr_w_016[0]);

  atotnw->SetBinContent(atotnw->FindBin(0), aVal_npr_nw_016[0]);
  atotnw->SetBinError(atotnw->FindBin(0),aErr_npr_nw_016[0]);

  ntotnw->SetBinContent(ntotnw->FindBin(0), nVal_npr_nw_016[0]);
  ntotnw->SetBinError(ntotnw->FindBin(0),nErr_npr_nw_016[0]);

  for (int i=2; i<10; i++)
    {
      aevol->SetBinContent(aevol->FindBin(i*0.1), aVal_npr_w_016[i]);
      aevol->SetBinError(aevol->FindBin(i*0.1),aErr_npr_w_016[i]);
      
      nevol->SetBinContent(nevol->FindBin(i*0.1), nVal_npr_w_016[i]);
      nevol->SetBinError(nevol->FindBin(i*0.1),nErr_npr_w_016[i]);

      aevolnw->SetBinContent(aevolnw->FindBin(i*0.1), aVal_npr_nw_016[i]);
      aevolnw->SetBinError(aevolnw->FindBin(i*0.1),aErr_npr_nw_016[i]);

      nevolnw->SetBinContent(nevolnw->FindBin(i*0.1), nVal_npr_nw_016[i]);
      nevolnw->SetBinError(nevolnw->FindBin(i*0.1),nErr_npr_nw_016[i]);
    }

  avr = aVal_npr_w_016[2]+aVal_npr_w_016[3]+aVal_npr_w_016[4]+aVal_npr_w_016[5]+aVal_npr_w_016[6]+aVal_npr_w_016[7]+aVal_npr_w_016[8];
  avr=avr/7;
  ave= new TLine(0,avr,1,avr);
  avr = aVal_npr_nw_016[2]+aVal_npr_nw_016[3]+aVal_npr_nw_016[4]+aVal_npr_nw_016[5]+aVal_npr_nw_016[6]+aVal_npr_nw_016[7]+aVal_npr_nw_016[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);
  aevol->SetMarkerColor(kRed);
  aevol->SetMarkerStyle(33);
  aevol->SetLineColor(kRed+2);
  aevol->SetOption("E1");
  //aevol->SetFillColor(kRed+2);

  aevolnw->SetMarkerColor(kBlue);
  aevolnw->SetMarkerStyle(33);
  aevolnw->SetMarkerStyle(kBlue+2);
  aevolnw->SetOption("E1");
  //aevolnw->SetFillColor(kBlue+2);
  atot->SetMarkerColor(kRed);
  atot->SetMarkerStyle(29);
  atot->SetMarkerSize(2);
  atot->SetLineColor(kRed+2);
  atot->SetOption("E1");

  atotnw->SetMarkerColor(kBlue);
  atotnw->SetMarkerStyle(29);
  atotnw->SetMarkerSize(2);
  atotnw->SetLineColor(kBlue+2);
  atotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  leg = new TLegend(0.15, 0.6, 0.35, 0.8);
  leg->AddEntry(aevol, "diff. with corr.", "lep");
  leg->AddEntry(aevolnw, "diff. without corr.", "lep");
  leg->AddEntry(atot, "int. with corr", "lep");
  leg->AddEntry(atotnw, "int. without corr", "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  tbox = new TPaveText(0.15,0.8,0.35,0.9, "BRNDC");
  tbox->AddText("npr MC");
  tbox->AddText("|y| < 1.6");
  tbox->SetBorderSize(0);
  tbox->SetFillColor(0);
  tbox->SetFillStyle(0);

  fsave->cd();
  aevol->Write("alpha_npr_w_016");
  aevolnw->Write("alpha_npr_nw_016");
  c->cd();
  aevol->Draw();
  aevolnw->Draw("same E1");
  ave->Draw("same");
  avenw->Draw("same");
  atot->Draw("same E1");
  atotnw->Draw("same E1");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/alpha_npr_016.png");

  avr = nVal_npr_w_016[2]+nVal_npr_w_016[3]+nVal_npr_w_016[4]+nVal_npr_w_016[5]+nVal_npr_w_016[6]+nVal_npr_w_016[7]+nVal_npr_w_016[8];
  avr=avr/7;
  ave = new TLine(0,avr,1, avr);
  avr = nVal_npr_nw_016[2]+nVal_npr_nw_016[3]+nVal_npr_nw_016[4]+nVal_npr_nw_016[5]+nVal_npr_nw_016[6]+nVal_npr_nw_016[7]+nVal_npr_nw_016[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);

  nevol->SetMarkerColor(kRed);
  nevol->SetMarkerStyle(33);
  nevol->SetLineColor(kRed+2);
  nevol->SetOption("E1");
  //nevol->SetFillColor(kRed+2);                                                                                                                                                                            
  nevolnw->SetMarkerColor(kBlue);
  nevolnw->SetMarkerStyle(33);
  nevolnw->SetMarkerStyle(kBlue+2);
  nevolnw->SetOption("E1");
  //nevolnw->SetFillColor(kBlue+2);

  ntot->SetMarkerColor(kRed);
  ntot->SetMarkerStyle(29);
  ntot->SetMarkerSize(2);
  ntot->SetLineColor(kRed+2);
  ntot->SetOption("E1");
  ntotnw->SetMarkerColor(kBlue);
  ntotnw->SetMarkerStyle(29);
  ntotnw->SetMarkerSize(2);
  ntotnw->SetLineColor(kBlue+2);
  ntotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  fsave->cd();
  nevol->Write("n_npr_w_016");
  nevolnw->Write("n_npr_nw_016"); 
  c->cd();
  nevol->Draw();
  nevolnw->Draw("same E1");
  ntot->Draw("same E1");
  ntotnw->Draw("same");
  ave->Draw("same");
  avenw->Draw("same");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/n_npr_016.png");



  atot->SetBinContent(atot->FindBin(0), aVal_pr_w_016[0]);
  atot->SetBinError(atot->FindBin(0),aErr_pr_w_016[0]);

  ntot->SetBinContent(ntot->FindBin(0), nVal_pr_w_016[0]);
  ntot->SetBinError(ntot->FindBin(0),nErr_pr_w_016[0]);

  atotnw->SetBinContent(atotnw->FindBin(0), aVal_pr_nw_016[0]);
  atotnw->SetBinError(atotnw->FindBin(0),aErr_pr_nw_016[0]);

  ntotnw->SetBinContent(ntotnw->FindBin(0), nVal_pr_nw_016[0]);
  ntotnw->SetBinError(ntotnw->FindBin(0),nErr_pr_nw_016[0]);

  for (int i=2; i<10; i++)
    {
      aevol->SetBinContent(aevol->FindBin(i*0.1), aVal_pr_w_016[i]);
      aevol->SetBinError(aevol->FindBin(i*0.1),aErr_pr_w_016[i]);
      
      nevol->SetBinContent(nevol->FindBin(i*0.1), nVal_pr_w_016[i]);
      nevol->SetBinError(nevol->FindBin(i*0.1),nErr_pr_w_016[i]);

      aevolnw->SetBinContent(aevolnw->FindBin(i*0.1), aVal_pr_nw_016[i]);
      aevolnw->SetBinError(aevolnw->FindBin(i*0.1),aErr_pr_nw_016[i]);

      nevolnw->SetBinContent(nevolnw->FindBin(i*0.1), nVal_pr_nw_016[i]);
      nevolnw->SetBinError(nevolnw->FindBin(i*0.1),nErr_pr_nw_016[i]);
    }

  avr = aVal_pr_w_016[2]+aVal_pr_w_016[3]+aVal_pr_w_016[4]+aVal_pr_w_016[5]+aVal_pr_w_016[6]+aVal_pr_w_016[7]+aVal_pr_w_016[8];
  avr=avr/7;
  ave= new TLine(0,avr,1,avr);
  avr = aVal_pr_nw_016[2]+aVal_pr_nw_016[3]+aVal_pr_nw_016[4]+aVal_pr_nw_016[5]+aVal_pr_nw_016[6]+aVal_pr_nw_016[7]+aVal_pr_nw_016[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);
  aevol->SetMarkerColor(kRed);
  aevol->SetMarkerStyle(33);
  aevol->SetLineColor(kRed+2);
  aevol->SetOption("E1");
  //aevol->SetFillColor(kRed+2);

  aevolnw->SetMarkerColor(kBlue);
  aevolnw->SetMarkerStyle(33);
  aevolnw->SetMarkerStyle(kBlue+2);
  aevolnw->SetOption("E1");
  //aevolnw->SetFillColor(kBlue+2);
  atot->SetMarkerColor(kRed);
  atot->SetMarkerStyle(29);
  atot->SetMarkerSize(2);
  atot->SetLineColor(kRed+2);
  atot->SetOption("E1");

  atotnw->SetMarkerColor(kBlue);
  atotnw->SetMarkerStyle(29);
  atotnw->SetMarkerSize(2);
  atotnw->SetLineColor(kBlue+2);
  atotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  leg = new TLegend(0.15, 0.6, 0.35, 0.8);
  leg->AddEntry(aevol, "diff. with corr.", "lep");
  leg->AddEntry(aevolnw, "diff. without corr.", "lep");
  leg->AddEntry(atot, "int. with corr", "lep");
  leg->AddEntry(atotnw, "int. without corr", "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  tbox = new TPaveText(0.15,0.8,0.35,0.9, "BRNDC");
  tbox->AddText("pr MC");
  tbox->AddText(" |y| < 1.6");
  tbox->SetBorderSize(0);
  tbox->SetFillColor(0);
  tbox->SetFillStyle(0);

  fsave->cd();
  aevol->Write("alpha_pr_w_016");
  aevolnw->Write("alpha_pr_nw_016");
  c->cd();
  aevol->Draw();
  aevolnw->Draw("same E1");
  ave->Draw("same");
  avenw->Draw("same");
  atot->Draw("same E1");
  atotnw->Draw("same E1");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/alpha_pr_016.png");

  avr = nVal_pr_w_016[2]+nVal_pr_w_016[3]+nVal_pr_w_016[4]+nVal_pr_w_016[5]+nVal_pr_w_016[6]+nVal_pr_w_016[7]+nVal_pr_w_016[8];
  avr=avr/7;
  ave = new TLine(0,avr,1, avr);
  avr = nVal_pr_nw_016[2]+nVal_pr_nw_016[3]+nVal_pr_nw_016[4]+nVal_pr_nw_016[5]+nVal_pr_nw_016[6]+nVal_pr_nw_016[7]+nVal_pr_nw_016[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);

  nevol->SetMarkerColor(kRed);
  nevol->SetMarkerStyle(33);
  nevol->SetLineColor(kRed+2);
  nevol->SetOption("E1");
  //nevol->SetFillColor(kRed+2);                                                                                                                                                                            
  nevolnw->SetMarkerColor(kBlue);
  nevolnw->SetMarkerStyle(33);
  nevolnw->SetMarkerStyle(kBlue+2);
  nevolnw->SetOption("E1");
  //nevolnw->SetFillColor(kBlue+2);

  ntot->SetMarkerColor(kRed);
  ntot->SetMarkerStyle(29);
  ntot->SetMarkerSize(2);
  ntot->SetLineColor(kRed+2);
  ntot->SetOption("E1");
  ntotnw->SetMarkerColor(kBlue);
  ntotnw->SetMarkerStyle(29);
  ntotnw->SetMarkerSize(2);
  ntotnw->SetLineColor(kBlue+2);
  ntotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  fsave->cd();
  nevol->Write("n_pr_w_016");
  nevolnw->Write("n_pr_nw_016"); 
  c->cd();
  nevol->Draw();
  nevolnw->Draw("same E1");
  ntot->Draw("same E1");
  ntotnw->Draw("same");
  ave->Draw("same");
  avenw->Draw("same");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/n_pr_016.png");


  atot->SetBinContent(atot->FindBin(0), aVal_npr_w_1624[0]);
  atot->SetBinError(atot->FindBin(0),aErr_npr_w_1624[0]);

  ntot->SetBinContent(ntot->FindBin(0), nVal_npr_w_1624[0]);
  ntot->SetBinError(ntot->FindBin(0),nErr_npr_w_1624[0]);

  atotnw->SetBinContent(atotnw->FindBin(0), aVal_npr_nw_1624[0]);
  atotnw->SetBinError(atotnw->FindBin(0),aErr_npr_nw_1624[0]);

  ntotnw->SetBinContent(ntotnw->FindBin(0), nVal_npr_nw_1624[0]);
  ntotnw->SetBinError(ntotnw->FindBin(0),nErr_npr_nw_1624[0]);

  for (int i=2; i<10; i++)
    {
      aevol->SetBinContent(aevol->FindBin(i*0.1), aVal_npr_w_1624[i]);
      aevol->SetBinError(aevol->FindBin(i*0.1),aErr_npr_w_1624[i]);
      
      nevol->SetBinContent(nevol->FindBin(i*0.1), nVal_npr_w_1624[i]);
      nevol->SetBinError(nevol->FindBin(i*0.1),nErr_npr_w_1624[i]);

      aevolnw->SetBinContent(aevolnw->FindBin(i*0.1), aVal_npr_nw_1624[i]);
      aevolnw->SetBinError(aevolnw->FindBin(i*0.1),aErr_npr_nw_1624[i]);

      nevolnw->SetBinContent(nevolnw->FindBin(i*0.1), nVal_npr_nw_1624[i]);
      nevolnw->SetBinError(nevolnw->FindBin(i*0.1),nErr_npr_nw_1624[i]);
    }

  avr = aVal_npr_w_1624[2]+aVal_npr_w_1624[3]+aVal_npr_w_1624[4]+aVal_npr_w_1624[5]+aVal_npr_w_1624[6]+aVal_npr_w_1624[7]+aVal_npr_w_1624[8];
  avr=avr/7;
  ave= new TLine(0,avr,1,avr);
  avr = aVal_npr_nw_1624[2]+aVal_npr_nw_1624[3]+aVal_npr_nw_1624[4]+aVal_npr_nw_1624[5]+aVal_npr_nw_1624[6]+aVal_npr_nw_1624[7]+aVal_npr_nw_1624[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);
  aevol->SetMarkerColor(kRed);
  aevol->SetMarkerStyle(33);
  aevol->SetLineColor(kRed+2);
  aevol->SetOption("E1");
  //aevol->SetFillColor(kRed+2);

  aevolnw->SetMarkerColor(kBlue);
  aevolnw->SetMarkerStyle(33);
  aevolnw->SetMarkerStyle(kBlue+2);
  aevolnw->SetOption("E1");
  //aevolnw->SetFillColor(kBlue+2);
  atot->SetMarkerColor(kRed);
  atot->SetMarkerStyle(29);
  atot->SetMarkerSize(2);
  atot->SetLineColor(kRed+2);
  atot->SetOption("E1");

  atotnw->SetMarkerColor(kBlue);
  atotnw->SetMarkerStyle(29);
  atotnw->SetMarkerSize(2);
  atotnw->SetLineColor(kBlue+2);
  atotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  leg = new TLegend(0.15, 0.6, 0.35, 0.8);
  leg->AddEntry(aevol, "diff. with corr.", "lep");
  leg->AddEntry(aevolnw, "diff. without corr.", "lep");
  leg->AddEntry(atot, "int. with corr", "lep");
  leg->AddEntry(atotnw, "int. without corr", "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  tbox = new TPaveText(0.15,0.8,0.35,0.9, "BRNDC");
  tbox->AddText("npr MC");
  tbox->AddText("1.6 < |y| < 2.4");
  tbox->SetBorderSize(0);
  tbox->SetFillColor(0);
  tbox->SetFillStyle(0);

  fsave->cd();
  aevol->Write("alpha_npr_w_1624");
  aevolnw->Write("alpha_npr_nw_1624");
  c->cd();
  aevol->Draw();
  aevolnw->Draw("same E1");
  ave->Draw("same");
  avenw->Draw("same");
  atot->Draw("same E1");
  atotnw->Draw("same E1");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/alpha_npr_1624.png");

  avr = nVal_npr_w_1624[2]+nVal_npr_w_1624[3]+nVal_npr_w_1624[4]+nVal_npr_w_1624[5]+nVal_npr_w_1624[6]+nVal_npr_w_1624[7]+nVal_npr_w_1624[8];
  avr=avr/7;
  ave = new TLine(0,avr,1, avr);
  avr = nVal_npr_nw_1624[2]+nVal_npr_nw_1624[3]+nVal_npr_nw_1624[4]+nVal_npr_nw_1624[5]+nVal_npr_nw_1624[6]+nVal_npr_nw_1624[7]+nVal_npr_nw_1624[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);

  nevol->SetMarkerColor(kRed);
  nevol->SetMarkerStyle(33);
  nevol->SetLineColor(kRed+2);
  nevol->SetOption("E1");
  //nevol->SetFillColor(kRed+2);                                                                                                                                                                            
  nevolnw->SetMarkerColor(kBlue);
  nevolnw->SetMarkerStyle(33);
  nevolnw->SetMarkerStyle(kBlue+2);
  nevolnw->SetOption("E1");
  //nevolnw->SetFillColor(kBlue+2);

  ntot->SetMarkerColor(kRed);
  ntot->SetMarkerStyle(29);
  ntot->SetMarkerSize(2);
  ntot->SetLineColor(kRed+2);
  ntot->SetOption("E1");
  ntotnw->SetMarkerColor(kBlue);
  ntotnw->SetMarkerStyle(29);
  ntotnw->SetMarkerSize(2);
  ntotnw->SetLineColor(kBlue+2);
  ntotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  fsave->cd();
  nevol->Write("n_npr_w_1624");
  nevolnw->Write("n_npr_nw_1624"); 
  c->cd();
  nevol->Draw();
  nevolnw->Draw("same E1");
  ntot->Draw("same E1");
  ntotnw->Draw("same");
  ave->Draw("same");
  avenw->Draw("same");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/n_npr_1624.png");



  atot->SetBinContent(atot->FindBin(0), aVal_pr_w_1624[0]);
  atot->SetBinError(atot->FindBin(0),aErr_pr_w_1624[0]);

  ntot->SetBinContent(ntot->FindBin(0), nVal_pr_w_1624[0]);
  ntot->SetBinError(ntot->FindBin(0),nErr_pr_w_1624[0]);

  atotnw->SetBinContent(atotnw->FindBin(0), aVal_pr_nw_1624[0]);
  atotnw->SetBinError(atotnw->FindBin(0),aErr_pr_nw_1624[0]);

  ntotnw->SetBinContent(ntotnw->FindBin(0), nVal_pr_nw_1624[0]);
  ntotnw->SetBinError(ntotnw->FindBin(0),nErr_pr_nw_1624[0]);

  for (int i=2; i<10; i++)
    {
      aevol->SetBinContent(aevol->FindBin(i*0.1), aVal_pr_w_1624[i]);
      aevol->SetBinError(aevol->FindBin(i*0.1),aErr_pr_w_1624[i]);
      
      nevol->SetBinContent(nevol->FindBin(i*0.1), nVal_pr_w_1624[i]);
      nevol->SetBinError(nevol->FindBin(i*0.1),nErr_pr_w_1624[i]);

      aevolnw->SetBinContent(aevolnw->FindBin(i*0.1), aVal_pr_nw_1624[i]);
      aevolnw->SetBinError(aevolnw->FindBin(i*0.1),aErr_pr_nw_1624[i]);

      nevolnw->SetBinContent(nevolnw->FindBin(i*0.1), nVal_pr_nw_1624[i]);
      nevolnw->SetBinError(nevolnw->FindBin(i*0.1),nErr_pr_nw_1624[i]);
    }

  avr = aVal_pr_w_1624[2]+aVal_pr_w_1624[3]+aVal_pr_w_1624[4]+aVal_pr_w_1624[5]+aVal_pr_w_1624[6]+aVal_pr_w_1624[7]+aVal_pr_w_1624[8];
  avr=avr/7;
  ave= new TLine(0,avr,1,avr);
  avr = aVal_pr_nw_1624[2]+aVal_pr_nw_1624[3]+aVal_pr_nw_1624[4]+aVal_pr_nw_1624[5]+aVal_pr_nw_1624[6]+aVal_pr_nw_1624[7]+aVal_pr_nw_1624[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);
  aevol->SetMarkerColor(kRed);
  aevol->SetMarkerStyle(33);
  aevol->SetLineColor(kRed+2);
  aevol->SetOption("E1");
  //aevol->SetFillColor(kRed+2);

  aevolnw->SetMarkerColor(kBlue);
  aevolnw->SetMarkerStyle(33);
  aevolnw->SetMarkerStyle(kBlue+2);
  aevolnw->SetOption("E1");
  //aevolnw->SetFillColor(kBlue+2);
  atot->SetMarkerColor(kRed);
  atot->SetMarkerStyle(29);
  atot->SetMarkerSize(2);
  atot->SetLineColor(kRed+2);
  atot->SetOption("E1");

  atotnw->SetMarkerColor(kBlue);
  atotnw->SetMarkerStyle(29);
  atotnw->SetMarkerSize(2);
  atotnw->SetLineColor(kBlue+2);
  atotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  leg = new TLegend(0.15, 0.6, 0.35, 0.8);
  leg->AddEntry(aevol, "diff. with corr.", "lep");
  leg->AddEntry(aevolnw, "diff. without corr.", "lep");
  leg->AddEntry(atot, "int. with corr", "lep");
  leg->AddEntry(atotnw, "int. without corr", "lep");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  tbox = new TPaveText(0.15,0.8,0.35,0.9, "BRNDC");
  tbox->AddText("pr MC");
  tbox->AddText("1.6 < |y| < 2.4");
  tbox->SetBorderSize(0);
  tbox->SetFillColor(0);
  tbox->SetFillStyle(0);

  fsave->cd();
  aevol->Write("alpha_pr_w_1624");
  aevolnw->Write("alpha_pr_nw_1624");
  c->cd();
  aevol->Draw();
  aevolnw->Draw("same E1");
  ave->Draw("same");
  avenw->Draw("same");
  atot->Draw("same E1");
  atotnw->Draw("same E1");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/alpha_pr_1624.png");

  avr = nVal_pr_w_1624[2]+nVal_pr_w_1624[3]+nVal_pr_w_1624[4]+nVal_pr_w_1624[5]+nVal_pr_w_1624[6]+nVal_pr_w_1624[7]+nVal_pr_w_1624[8];
  avr=avr/7;
  ave = new TLine(0,avr,1, avr);
  avr = nVal_pr_nw_1624[2]+nVal_pr_nw_1624[3]+nVal_pr_nw_1624[4]+nVal_pr_nw_1624[5]+nVal_pr_nw_1624[6]+nVal_pr_nw_1624[7]+nVal_pr_nw_1624[8];
  avr=avr/7;
  avenw = new TLine(0,avr, 1, avr);

  nevol->SetMarkerColor(kRed);
  nevol->SetMarkerStyle(33);
  nevol->SetLineColor(kRed+2);
  nevol->SetOption("E1");
  //nevol->SetFillColor(kRed+2);                                                                                                                                                                            
  nevolnw->SetMarkerColor(kBlue);
  nevolnw->SetMarkerStyle(33);
  nevolnw->SetMarkerStyle(kBlue+2);
  nevolnw->SetOption("E1");
  //nevolnw->SetFillColor(kBlue+2);

  ntot->SetMarkerColor(kRed);
  ntot->SetMarkerStyle(29);
  ntot->SetMarkerSize(2);
  ntot->SetLineColor(kRed+2);
  ntot->SetOption("E1");
  ntotnw->SetMarkerColor(kBlue);
  ntotnw->SetMarkerStyle(29);
  ntotnw->SetMarkerSize(2);
  ntotnw->SetLineColor(kBlue+2);
  ntotnw->SetOption("E1");

  ave->SetLineColor(kRed);
  ave->SetLineStyle(2);
  avenw->SetLineColor(kBlue);
  avenw->SetLineStyle(2);

  fsave->cd();
  nevol->Write("n_pr_w_1624");
  nevolnw->Write("n_pr_nw_1624"); 
  c->cd();
  nevol->Draw();
  nevolnw->Draw("same E1");
  ntot->Draw("same E1");
  ntotnw->Draw("same");
  ave->Draw("same");
  avenw->Draw("same");
  leg->Draw("same");
  tbox->Draw("same");
  c->SaveAs("Output/MCPars/n_pr_1624.png");

  fsave->Close();
}
