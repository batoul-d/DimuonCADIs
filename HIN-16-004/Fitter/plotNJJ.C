#include "Macros/CMS/CMS_lumi.C"
#include "Macros/CMS/tdrstyle.C"
#include "Macros/Utilities/bin.h"
#include "Macros/Utilities/EVENTUTILS.h"
#include "Macros/Utilities/initClasses.h"
#include "Macros/Utilities/resultUtils.h"
#include "Macros/Utilities/texUtils.h"
#include "Macros/Utilities/monster.h"
#include "Systematics/syst.h"

#include <vector>
#include <map>
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

bool doprompt       = true;
bool dononprompt    = false;
bool is18XXX        = true;
bool doLogPt        = false;
bool  isPreliminary = false;
string nameTag_base = "_prompt";    // can put here e.g. "_prompt", "_nonprompt", ...
const bool useNcoll = false; // false -> use TAA / NMB, true -> use Ncoll / lumiPbPb
//////////////////
// DECLARATIONS //
//////////////////

void printOptions();
void setOptions(bool adoprompt, bool adononprompt, bool ais18XXX, bool adoLogPt, string anameTag_base="");
void plotNJJ(vector<anabin> thecats, string xaxis, string outputDir);
void plotGraphNJJ(map<anabin, TGraphAsymmErrors*> theGraphs, string xaxis, string outputDir);
int color(int i);
int markerstyle(int i);
string nameTag;

class njj_input {
public:
  double npp;
  double dnpp;
  double systpp;
  double naa;
  double dnaa_stat;
  double systaa;
  double taa;
  double lumipp;
  double lumiaa;
  double ncoll;
  double bfracpp;
  double dbfracpp;
  double systbfracpp;
  double bfracaa;
  double dbfracaa;
  double systbfracaa;
  double correlaa;
  double correlpp;
  syst   statpp;
};


/////////////////////////////////////////////
// MAIN FUNCTIONS TO BE CALLED BY THE USER //
/////////////////////////////////////////////


void plotZed(string workDirName, string poiname, int iplot) {
  
  string xaxis = "zed";
  vector<anabin> theCats;  

  // 10 z bins in 3 rap intervals 
  //if (iplot==0) {
  //theCats.push_back(anabin(0.0,1.0,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.0,0.1,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.1,0.2,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.2,0.3,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.3,0.4,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.4,0.5,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.5,0.6,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.6,0.7,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.7,0.8,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.8,0.9,0,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.9,1.0,0,2.4,6.5,50,0,200));
  //}
  if (iplot==0)
    theCats.push_back(anabin(0,1.0,0,2.4,6.5,50,0,200));
  if (iplot==1)
    theCats.push_back(anabin(0,1.0,0,1.6,6.5,50,0,200));
  if (iplot==2)
    theCats.push_back(anabin(0,1.0,1.6,2.4,3,50,0,200));
  
  //if (iplot==1) {
  //theCats.push_back(anabin(0.0,1.0,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.0,0.1,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.1,0.2,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.2,0.3,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.3,0.4,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.4,0.5,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.5,0.6,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.6,0.7,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.7,0.8,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.8,0.9,0,1.6,6.5,50,0,200));
  //theCats.push_back(anabin(0.9,1.0,0,1.6,6.5,50,0,200));
  //}

  //if (iplot==2) {
  //theCats.push_back(anabin(0.0,1.0,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.0,0.1,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.1,0.2,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.2,0.3,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.3,0.4,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.4,0.5,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.5,0.6,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.6,0.7,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.7,0.8,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.8,0.9,1.6,2.4,6.5,50,0,200));
  //theCats.push_back(anabin(0.9,1.0,1.6,2.4,6.5,50,0,200));
  //}
  
  
  nameTag = nameTag_base + Form("_%i",iplot);
  
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
  
};

void plotPt(string workDirName, string poiname, int iplot) {
  
  string xaxis = "pt";
  vector<anabin> theCats;  

  // 3 rap intervals integrated in z
  if (iplot==0) {
    theCats.push_back(anabin(0,1,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0,1,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0,1,1.6,2.4,3,50,0,200));
  }
  
  
  nameTag = nameTag_base + Form("_%i",iplot);
  
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
  
};

void plotCent(string workDirName, string poiname, int iplot) {
  
  string xaxis = "cent";
  vector<anabin> theCats;
  
  // 3 rapidity intervals integrated in z
  if (iplot==0) {
    theCats.push_back(anabin(0,1,0,2.4,6.5,50,0,200));
    theCats.push_back(anabin(0,1,0,1.6,6.5,50,0,200));
    theCats.push_back(anabin(0,1,1.6,2.4,3,50,0,200));
  }
   
  nameTag = nameTag_base + Form("_%i",iplot);
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
};


void plotRap(string workDirName, string poiname) {
  
  string xaxis = "rap";
  vector<anabin> theCats;
  
  theCats.push_back(anabin(0,1,0,2.4,6.5,50,0,200));
  
  nameTag = nameTag_base;
  if (poiname.find("NJJ")!=std::string::npos) plotNJJ(theCats,xaxis,workDirName);
  else {
    cout << "[ERROR] : The observable you want to plot is not supported" << endl;
    return;
  }
};


void plotAll(string workDirName, string poiname) {
//  if (dononprompt) nameTag_base = "_nonprompt";
  if (!doprompt && !dononprompt) nameTag_base = "";
  //workDirName= workDirName+"/mass/DATA";
  
  if (is18XXX)
  {
    plotZed(workDirName,poiname,2);
    plotPt(workDirName,poiname,2);
  }
};

void doAllplots(bool is18XXX=false) {
  
  if (is18XXX)
  {
    setOptions(true, false, true, false);
    printOptions();
    plotAll("test/Datafits_1624/mass/DATA","NJJ");
  }

};


/////////////////////
// OTHER FUNCTIONS //
/////////////////////

void plotNJJ(vector<anabin> thecats, string xaxis, string outputDir){
  if (doprompt && dononprompt) {
    cout << "ERROR you can't set both doprompt and dononprompt to true." << endl;
    return;
  }
  float jtptmin=0, jtptmax=0;
  if (outputDir.find("mid")!=std::string::npos) {jtptmin = 25; jtptmax=35;}
  if (outputDir.find("low")!=std::string::npos) {jtptmin = 15; jtptmax=25;}
  if (outputDir.find("high")!=std::string::npos) {jtptmin = 35; jtptmax=45;}
  TString poi("NJpsi");
  if (doprompt) poi = "NJpsi_prompt";
  if (dononprompt) poi = "NJpsi_nonprompt";
  
  string plotLabel = "";
  TFile *f = new TFile(treeFileName(outputDir.c_str(),""));
  if (!f || !f->IsOpen()) {
    results2tree(outputDir.c_str(),"");
    f = new TFile(treeFileName(outputDir.c_str(),""));
    if (!f) return;
  }
  TTree *tr = (TTree*) f->Get("fitresults");
  if (!tr) return;

  TH1F * nprRes = new TH1F ("npRres", ";z(J/#psi);dN(J/#psi-jet)/N", 5, 0.0, 1);
  TH1F * prRes = new TH1F ("prRes", ";z(J/#psi);dN(J/#psi-jet)/N", 5, 0.0, 1);

  //nprRes->GetYaxis()->SetLimits(0, 1);
  //prRes->GetYaxis()->SetLimits(0, 1);

  TString sTag("16025");
  if (is18XXX) sTag = "18XXX";
  map<anabin, njj_input> theVars_inputs;

  vector<double> x, ex, y, ey;
  float zmin, zmax, ptmin, ptmax, ymin, ymax, centmin, centmax;
  float /*eff, acc,*/ lumi, taa, ncoll;
  float val, errL=0, errH=0;
  float bfrac, bfrac_errL,bfrac_errH;
  float correl=0;
  int ival=-999;
  char collSystem[5];
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
  //if (!dononprompt)
  //{
  //tr->SetBranchAddress("eff_val",&eff);
  //tr->SetBranchAddress("acc_val",&acc);
  //}
  //else
  //{
  //tr->SetBranchAddress("effnp_val",&eff);
  //tr->SetBranchAddress("accnp_val",&acc);
  //}
  tr->SetBranchAddress("lumi_val",&lumi);
  tr->SetBranchAddress("taa_val",&taa);
  tr->SetBranchAddress("ncoll_val",&ncoll);
  tr->SetBranchAddress("b_Jpsi_val",&bfrac);
  tr->SetBranchAddress("b_Jpsi_errL",&bfrac_errL);
  tr->SetBranchAddress("b_Jpsi_errH",&bfrac_errH);
  tr->SetBranchAddress("correl_N_Jpsi_vs_b_Jpsi_val",&correl);
  
  int ntr = tr->GetEntries();
  for (int i=0; i<ntr; i++) {
    tr->GetEntry(i);
    
    anabin thebin(zmin, zmax, ymin, ymax, ptmin, ptmax, centmin, centmax);
    
    bool ispp = (TString(collSystem)=="PP");
    ispp = true;

    if (ispp){    
      theVars_inputs[thebin].npp = val;
      theVars_inputs[thebin].dnpp = errL;
      theVars_inputs[thebin].bfracpp = bfrac;
      theVars_inputs[thebin].dbfracpp = bfrac_errL;
      //theVars_inputs[thebin].systpp = syst_PP[thebin].value;
      //syst thestat_PP; thestat_PP.name = "stat_PP"; thestat_PP.value = errL/val;
      //stat_PP[thebin] = thestat_PP;
      //theVars_inputs[thebin].statpp = thestat_PP;
      theVars_inputs[thebin].lumipp = lumi;
      theVars_inputs[thebin].correlpp = correl;
      cout<<"Njj = "<<val<<"  dNjj = "<<errL<<endl;
      if (i==0)
	{
	  plotLabel = Form("_rap%.0f%.0f_pt%.0f%.f",(ymin*10.0), (ymax*10.0), (ptmin*10.0), (ptmax*10.0));
	}
      if (i!=0) {
	prRes->SetBinContent(prRes->FindBin(zmin+0.03),val*(1-bfrac));
	prRes->SetBinError(prRes->FindBin(zmin+0.03),(val*(1-bfrac)*sqrt(pow((errL/val),2) - 2*correl*errL*bfrac_errL/(val*bfrac) + pow(bfrac_errL/bfrac,2))));
        nprRes->SetBinContent(nprRes->FindBin(zmin+0.03),val*bfrac);
        nprRes->SetBinError(nprRes->FindBin(zmin+0.03),(val*(1-bfrac)*sqrt(pow((errL/val),2) + 2*correl*errL*bfrac_errL/(val*bfrac) + pow(bfrac_errL/bfrac,2))));
      }
    }
  }

  //nprRes->Scale(1.0/(nprRes->Integral()));
  nprRes->SetMarkerColor(2);
  nprRes->SetMarkerStyle(33);
  //prRes->Scale(1.0/(prRes->Integral()));
  prRes->SetMarkerColor(2);
  prRes->SetMarkerStyle(33);
  TPaveText* tb = new TPaveText(0.15,0.6,0.4,0.8,"BRNDC");
  tb->AddText("Nonprompt J/#psi");
  if (ymin ==0)
    tb->AddText(Form("|y| < %.1f",ymax));
  else 
    tb->AddText(Form("%.1f < |y| < %.1f", ymin, ymax));
  tb->AddText(Form("%.1f < p_{T}(J/#psi) < %.1f GeV/c", ptmin, ptmax));
  tb->AddText(Form("%.1f < p_{T}(jet) < %.1f GeV/c", jtptmin, jtptmax));
  tb->SetBorderSize(0);
  tb->SetFillColor(0);

  TPaveText* tb1 = new TPaveText(0.15,0.6,0.4,0.8,"BRNDC");
  tb1->AddText("Prompt J/#psi");
  if (ymin ==0)
    tb1->AddText(Form("|y| < %.1f",ymax));
  else 
    tb1->AddText(Form("%.1f < |y| < %.1f", ymin, ymax));
  tb1->AddText(Form("%.1f < p_{T}(J/#psi) < %.1f GeV/c", ptmin, ptmax));
  tb1->AddText(Form("%.1f < p_{T}(jet) < %.1f GeV/c", jtptmin, jtptmax));
  tb1->SetBorderSize(0);
  tb1->SetFillColor(0);

  gSystem->mkdir("Results");
  TCanvas *cres = new TCanvas ("cres","",1000,800);
  cres->cd();
  nprRes->Draw();
  tb->Draw("same");
  cres->SaveAs(Form("Results/plot_DATA_NPr%s_jtpt%.0f%.0f.root", plotLabel.c_str(), jtptmin, jtptmax));

  prRes->Draw();
  tb1->Draw("same");
  cres->SaveAs(Form("Results/plot_DATA_Pr%s_jtpt%.0f%.0f.root", plotLabel.c_str(), jtptmin, jtptmax));

  TFile *fsave = new TFile (Form("Results/histo_DATA%s_jtpt%.0f%.0f.root", plotLabel.c_str(), jtptmin, jtptmax),"RECREATE");
  prRes->Write("PRHist");
  nprRes->Write("NPRHist");
  fsave->Close();
  map<anabin, vector<anabin> > theBins;
  map<anabin, vector<double> > theVarsBinned_pp;
  map<anabin, vector<double> > theVarsBinned_stat_pp;;
  map<anabin, TGraphAsymmErrors* > theGraphs_pp;

  //map<anabin, njj_input > theResults18XXX; //In case they are needed to compute the Psi2S NJJ

  // initialize the maps
  for (vector<anabin>::const_iterator it=thecats.begin(); it!=thecats.end(); it++) {
    theBins[*it] = vector<anabin>();
    theVarsBinned_pp[*it] = vector<double>();
  }
  
  for (map<anabin, njj_input>::const_iterator it=theVars_inputs.begin(); it!=theVars_inputs.end(); it++) {
    anabin thebin = it->first;
    anabin thebinOrig = it->first; // Original bin to retrieve results later if needed (cause binok() will overwrite thebin)
    njj_input s = it->second;
    if (!binok(thecats,xaxis,thebin)) continue;
    anabin thebinPP = it->first; thebinPP.setcentbin(binI(0,200));
    njj_input spp = theVars_inputs[thebinPP];
    if (s.npp<=0) continue;
    theBins[thebin].push_back(it->first);
    
    double npp  = spp.npp;
    double dnpp = spp.dnpp;

    if (doprompt) {
    npp = spp.npp*(1.-spp.bfracpp);
    dnpp = npp*sqrt(pow(spp.dnpp/spp.npp,2)
                    - 2.*spp.correlpp*spp.dnpp*spp.dbfracpp/(spp.npp*spp.bfracpp)
                    + pow(spp.dbfracpp/spp.bfracpp,2));
    }
    if (dononprompt) {
    npp = spp.npp*spp.bfracpp;
    dnpp = npp*sqrt(pow(spp.dnpp/spp.npp,2)
                    + 2.*spp.correlpp*spp.dnpp*spp.dbfracpp/(spp.npp*spp.bfracpp)
                    + pow(spp.dbfracpp/spp.bfracpp,2));
    }

    double njj = npp;
    double dnjj = njj>0 ? dnpp:0;

    theVarsBinned_pp[thebin].push_back(njj);
    theVarsBinned_stat_pp[thebin].push_back(dnjj);
  }

    int cnt=0;
    for (vector<anabin>::const_iterator it=thecats.begin();it!=thecats.end(); it++){
      int n = theBins[*it].size();
      if (n==0) {
	cout << "Error, nothing found for this category" << endl;
	theGraphs_pp[*it] = NULL;
	continue;
      }
      theGraphs_pp[*it] = new TGraphAsymmErrors(n);
      theGraphs_pp[*it]->SetName(Form("bin_%i_pp",cnt));
      for (int i=0; i<n; i++) {
	double x=0, exl=0, exh=0,  ypp=0, eylpp=0, eyhpp=0;
	double low=0, high=0;
	anabin thebin = theBins[*it][i];
	ypp = theVarsBinned_pp[*it][i];
	if (xaxis=="pt" || xaxis=="rap" || xaxis=="zed") {
	  if (xaxis=="pt") {
	    low= thebin.ptbin().low();
	    high = thebin.ptbin().high();
	  } else if (xaxis=="rap"){
	    low= thebin.rapbin().low();
	    high = thebin.rapbin().high();
	  }
	  else {
	    low= thebin.zbin().low();
	    high = thebin.zbin().high();
	  }
	  x = (low+high)/2.;
	  exh = (high-low)/2.;
	  exl = (high-low)/2.;
	}
	eylpp = fabs(theVarsBinned_stat_pp[*it][i]);
	eyhpp = eylpp;
	theGraphs_pp[*it]->SetPoint(i,x,ypp);
	theGraphs_pp[*it]->SetPointError(i,exl,exh,eylpp,eyhpp);
      }
      cnt++;
    }
  plotGraphNJJ(theGraphs_pp, xaxis, outputDir);
}

void plotGraphNJJ(map<anabin, TGraphAsymmErrors*> theGraphs, string xaxis, string outputDir) {
  //setTDRStyle();
  const char* ylabel = "N_{JJ}";
  int intervals2Plot = theGraphs.size();

  vector<anabin> theCats;
  
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",600,600);
  gStyle->SetOptStat(0);
  // the axes
  TH1F *haxes=NULL; TLine line;
  if (xaxis=="pt") {
    haxes = new TH1F("haxes","", 1, 0, 50);
    //line = TLine(0,1,50,1);
  }
  if (xaxis=="zed"){
    haxes = new TH1F("haxes","", 1, 0, 1);
    //line = TLine(0,1,1,1);
  }
  if (xaxis=="rap") {
    haxes = new TH1F("haxes","",1,0,2.4);
    haxes->GetXaxis()->SetNdivisions(306,false);
    //line = TLine(0,1,2.4,1);
  }
  if (xaxis=="cent") {
    haxes = new TH1F("haxes","",1,0,420);
    haxes->GetXaxis()->SetTickLength(gStyle->GetTickLength("X"));
    line = TLine(0,1,420,1);
  }
  haxes->GetYaxis()->SetRangeUser(0,1.5);
  haxes->GetYaxis()->SetTitle(ylabel);
  const char* xlabel = (xaxis=="pt") ? "p_{T} (GeV/c)" : ((xaxis=="rap") ? "|y|" : "z");
  haxes->GetXaxis()->SetTitle(xlabel);
  haxes->GetXaxis()->CenterTitle(true);
  haxes->Draw();

  double xshift=0.025;
  TLegend *tleg(0x0);
  if (xaxis!="cent" && intervals2Plot == 2) tleg = new TLegend(0.44,0.50,0.76,0.62);
  else if (xaxis=="cent" && intervals2Plot == 2) tleg = new TLegend(0.19,0.16,0.51,0.28);
  else if (dononprompt && intervals2Plot == 3) tleg = new TLegend(0.56,0.47,0.88,0.62);
  else if (doprompt && intervals2Plot == 3) tleg = new TLegend(0.19,0.49,0.51,0.64);
  else tleg = new TLegend(0.56,0.42,0.88,0.62);
  tleg->SetBorderSize(0);
  tleg->SetFillStyle(0);
  tleg->SetTextFont(42);
  tleg->SetTextSize(0.04);
  bool drawLegend(true);
  
  const char* xname = (xaxis=="cent") ? "Centrality" : (xaxis=="pt" ? "\\pt" : (xaxis=="rap"?"$|y|$":"z"));
  gSystem->mkdir(Form("Output/%s/tex/", outputDir.c_str()), kTRUE);
  char texname[2048]; sprintf(texname, "Output/%s/tex/result_%s_NJJ_%s%s%s_%s.tex",outputDir.c_str(), "JPsi",xaxis.c_str(),nameTag.c_str(), (xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : "", "noCorr");
  string yname("$\\N (\\Jpsi-Jet)$");
  inittex(texname, xname, yname);
  
  int cnt=0;
  map<anabin, TGraphAsymmErrors*>::const_iterator it=theGraphs.begin();
  for (; it!=theGraphs.end(); it++) {
    anabin thebin = it->first;
    TGraphAsymmErrors* tg = it->second;
    if (!tg) continue;
    
    theCats.push_back(thebin);
    int style = cnt;
    int colorI = cnt;
    int colorF = color(colorI)-11;

    if (intervals2Plot==2)
      {
	style = 4;
	cnt==0 ? colorI = 9 : colorI =4;
	colorF = color(colorI)-11;
      }
    if (intervals2Plot==3)
      {
	style = 0;
	colorI = cnt+6;
	colorF = color(colorI)-10;
      }
    if (intervals2Plot>3)
      {
	style = cnt+1;
	colorI = cnt+1;
	colorF = color(colorI)-11;
      }
    
    tg->SetMarkerStyle(markerstyle(style));
    tg->SetMarkerColor(color(colorI));
    tg->SetLineColor(color(colorI));
    //tg_syst->SetLineColor(color(colorI));
    //tg_syst->SetFillColorAlpha(colorF, 0.5);
    if (markerstyle(style) == kFullStar) tg->SetMarkerSize(2.3);
    else if (markerstyle(style) == kFullDiamond) tg->SetMarkerSize(2.2);
    else if (markerstyle(style) == kFullCross) tg->SetMarkerSize(2.0);
    else tg->SetMarkerSize(1.5);
    tg->SetLineWidth(tg->GetLineWidth()*2);
    
    //if (xaxis=="cent") {
      // do not plot wide centrality bins
      //prune(tg, tg_syst);
      //}
    bool plot = true;

    if (plot)
      {
	double x(0.), dx(0.), y(0.), dy_low(0.), dy_high(0.);
	double rightA = 0.;
	int centminGlob(0),centmaxGlob(200);
	if (xaxis=="cent") {
	  dx = 10;
	  rightA = 420.;
	} else if (xaxis=="pt") {
	  dx = 0.625;
	  if (intervals2Plot == 1)
	    {
	      rightA = 50.;
	      dx = 1.25;
	    }
	  else
	    {
	      rightA = 30.;
	      dx = 0.65;
	    }
	  if (intervals2Plot == 3)
	    {
	      centminGlob = it->first.centbin().low();
	      centmaxGlob = it->first.centbin().high();
	    }
	} else if (xaxis=="rap") {
	  dx = 0.06;
	  rightA = 2.4;
	}
	else if (xaxis=="zed"){
	  dx=0.05;
	  rightA = 1;
	}
	x = rightA - (2*dx*cnt + dx);
	//if ( cnt>0) x = rightA - (2*dx + dx);
	y = 1;
      
	anabin thebinglb(it->first.zbin().low(),
			 it->first.zbin().high(),
			 it->first.rapbin().low(),
			 it->first.rapbin().high(),
			 it->first.ptbin().low(),
			 it->first.ptbin().high(),
			 centminGlob,centmaxGlob);
	thebinglb.print();
	TBox *tbox = new TBox(x-dx,y-dy_low,x+dx,y+dy_high);
	tbox->SetFillColorAlpha(colorF, 1);
	tbox->SetLineColor(color(colorI));
	tbox->Draw("l");
	TBox *tboxl = (TBox*) tbox->Clone("tboxl");
	tboxl->SetFillStyle(0);
	tboxl->Draw("l");

	gStyle->SetEndErrorSize(5);
	tg->Draw("P");
      }    

    TLatex tl;
    double tlx = 0.25; //0.92;
    double tly = 0.80; //0.69;
    tl.SetTextFont(42); // consistent font for symbol and plain text
    if (plot)
      {
        TString zlabel = Form("%g < z < %g GeV/c",it->first.zbin().low(), it->first.zbin().high());
        TString raplabel = Form("%.1f < |y| < %.1f",it->first.rapbin().low(),it->first.rapbin().high());
        if (it->first.rapbin().low()<0.1) raplabel = Form("|y| < %.1f",it->first.rapbin().high());
        TString ptlabel = Form("%g < p_{T} < %g GeV/c",it->first.ptbin().low(), it->first.ptbin().high());
        TString centlabel = Form("%i-%i%s",(int) (it->first.centbin().low()/2.), (int) (it->first.centbin().high()/2.), "%");
        
        if (xaxis == "pt")
	  {
	    if (is18XXX) tleg->AddEntry(tg, raplabel, "p");
	  }
        if (xaxis == "cent")
	  {
	    if (is18XXX) tleg->AddEntry(tg, Form("#splitline{%s}{%s}",raplabel.Data(),ptlabel.Data()), "p");
	  }
        if (xaxis == "rap")
	  {
	    if (intervals2Plot > 3) tleg->AddEntry(tg, ptlabel, "p");
	    else if (intervals2Plot == 1) drawLegend = false;
	    else tl.DrawLatexNDC(tlx,tly,Form("#splitline{%s}{Cent. %s}",ptlabel.Data(),centlabel.Data()));
	  }
	if (xaxis =="zed")
	  if(is18XXX) tleg->AddEntry(tg, raplabel, "p");
      }
    ostringstream oss;
    oss.precision(1); oss.setf(ios::fixed);
    oss << "$" << it->first.rapbin().low() << "<|y|<" << it->first.rapbin().high() << "$, ";
    if (xaxis == "pt") oss << (int) (it->first.zbin().low()) << "\\% - " << (int) (it->first.zbin().high()) << "\\%";
    if (xaxis == "cent" || xaxis == "rap" || xaxis == "zed") oss << "$" << it->first.ptbin().low() << "<\\pt<" << it->first.ptbin().high() << "\\GeVc $";
    
    addline(texname,oss.str());
    //printGraph(tg, texname);
   
    cnt++;
  }

  if (drawLegend) tleg->Draw();
  line.Draw();
  
  if(xaxis!="cent" && intervals2Plot > 3)
  {
    //      TLatex *tex = new TLatex(0.21,0.86,"Cent. 0-100%");
    TLatex *tex = new TLatex(0.2,0.78,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 2 && !is18XXX)
  {
    TLatex *tex = new TLatex(0.2,0.78,"1.8 < |y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot > 3)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="cent" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"6.5 < p_{T} < 50 GeV/c");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  if(xaxis=="pt" && intervals2Plot == 3)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="pt" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"|y| < 2.4");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"Cent. 0-100%");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  if(xaxis=="pt" && is18XXX)
  {
    TLatex *tex = new TLatex(0.2,0.78,"Cent. 0-100%");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  if(xaxis=="rap" && intervals2Plot == 1)
  {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    
    TLatex *tex2 = new TLatex(0.2,0.73,"Cent. 0-100%");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
  }
  if(xaxis=="zed" && is18XXX)
    {
    TLatex *tex = new TLatex(0.2,0.78,"6.5 < p_{T} < 50 GeV/c");
    tex->SetNDC();
    tex->SetTextSize(0.044);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

    TLatex *tex2 = new TLatex(0.2,0.73,"|y| < 2.4");
    tex2->SetNDC();
    tex2->SetTextSize(0.044);
    tex2->SetTextFont(42);
    tex2->SetLineWidth(2);
    tex2->Draw();
    }
  
  TLatex tl;

  double tlx = 0.20; //0.92;
  double tly = 0.85; //0.69;

  tl.SetTextFont(42); // consistent font for symbol and plain text
  tl.SetTextSize(0.057); 
  if (doprompt) tl.DrawLatexNDC(tlx,tly, "Prompt J/#psi");
  if (dononprompt) tl.DrawLatexNDC(tlx,tly,"J/#psi from b hadrons");
  tl.SetTextSize(0.046);
  
  int iPos = 33;
  if (xaxis=="cent") CMS_lumi( (TPad*) gPad, 1061, iPos, "", isPreliminary );
  else CMS_lumi( (TPad*) gPad, 106, iPos, "", isPreliminary );
  
  c1->cd();
  c1->Update();
  c1->RedrawAxis();
  gSystem->mkdir(Form("Output/%s/plot/RESULT/root/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/root/result_%s_NJJ_%s%s%s.root",outputDir.c_str(), "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : ""));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/png/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/png/result_%s_NJJ_%s%s%s.png",outputDir.c_str(), "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : ""));
  gSystem->mkdir(Form("Output/%s/plot/RESULT/pdf/", outputDir.c_str()), kTRUE);
  c1->SaveAs(Form("Output/%s/plot/RESULT/pdf/result_%s_NJJ_%s%s%s.pdf",outputDir.c_str(),  "JPsi", xaxis.c_str(), nameTag.c_str(),(xaxis=="pt") ? (doLogPt ? "_logX" :"_linearX") : ""));
  
  delete tleg;
  delete haxes;
  delete c1;
  
  // close tex
  closetex(texname);
  cout << "Closed " << texname << endl;
}


int color(int i) {
  if (i==0) return kRed+2;
  else if (i==1) return kBlue+2;
  else if (i==2) return kMagenta+2;
  else if (i==3) return kCyan+2;
  else if (i==4) return kGreen+2;
  else if (i==5) return kOrange+2;
  else if (i==6) return kRed+1;
  else if (i==7) return kYellow+1;
  else if (i==8) return kAzure+1;
  else if (i==9) return kBlack;
  else return kBlack;
}

int markerstyle(int i) {
  if (i==0) return kFullSquare;
  else if (i==1) return kFullCircle;
  else if (i==2) return kFullStar;
  else if (i==3) return kFullCross;
  else if (i==4) return kFullDiamond;
  else if (i==5) return kOpenSquare;
  else if (i==6) return kOpenCircle;
  else if (i==7) return kOpenStar;
  else if (i==8) return kFullTriangleDown;
  else return kOpenCross;
}

void setOptions(bool adoprompt, bool adononprompt, bool ais18XXX, bool adoLogPt, string anameTag_base){
  doprompt = adoprompt;
  dononprompt = adononprompt;
  is18XXX = ais18XXX;
  doLogPt = adoLogPt;
  nameTag_base = anameTag_base;
  if (doprompt) nameTag_base += "_prompt";
  if (dononprompt) nameTag_base += "_nonprompt";
  if (is18XXX)  nameTag_base += "_18XXX";
}

void printOptions() {
  cout <<
  "doprompt = " << doprompt << ", " <<
  "dononprompt = " << dononprompt << ", " <<
  "is18XXX = " << is18XXX << ", " <<
  "doLogPt = " << doLogPt << ", " <<
  "nameTag_base = \"" << nameTag_base << "\"" <<
  endl;
}



