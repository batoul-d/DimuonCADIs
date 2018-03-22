#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TPaveText.h>
#include "RooCategory.h"
#include "RooNumIntConfig.h"
#include "RooPlotable.h"
#include <TUnfold.h>
#include <TLorentzVector.h>
#include <vector>
#include <TRandom.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TEfficiency.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TStyle.h"
// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

using namespace std;

void GetMCResults (bool isPr=true, bool isMid=true, bool underflowOff=true);

void GetMCResults_all ()
{
  GetMCResults(true, true, true);
  GetMCResults(true, true, false);
  GetMCResults(true, false, true);
  GetMCResults(true, false, false);

  GetMCResults(false, true, true);
  GetMCResults(false, true, false);
  GetMCResults(false, false, true);
  GetMCResults(false, false, false);
}

void GetMCResults (bool isPr, bool isMid, bool underflowOff)
{
  TFile *treeFile = TFile::Open(Form("/data_CMS/cms/mnguyen/jPsiJet/mc/%s/v9_ext/merged_HiForestAOD.root", isPr?"prompt":"nonprompt"));
  TTree* oniaTree = (TTree*) treeFile->Get("hionia/myTree");
  TTree* jetTree = (TTree*) treeFile->Get("ak4PFJetAnalyzer/t");
  oniaTree->AddFriend(jetTree);

  TH1F* zHist = new TH1F ("zHist","", 5, 0, 1); zHist->Sumw2();
  double zed;
  double drmin;

  TTree* fChain;

  Int_t           Gen_QQ_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_QQ_mupl_4mom;
  TClonesArray    *Gen_QQ_mumi_4mom;
  Int_t           ngen;
  Float_t         genpt[28];   //[ngen]
  Float_t         geneta[28];   //[ngen]
  Float_t         geny[28];   //[ngen]
  Float_t         genphi[28];   //[ngen]
  Float_t         genm[28];   //[ngen]

  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_QQ_mupl_4mom;   //!
  TBranch        *b_Gen_QQ_mumi_4mom;   //!
  TBranch        *b_ngen;   //!
  TBranch        *b_genpt;   //!
  TBranch        *b_geneta;   //!
  TBranch        *b_geny;   //!
  TBranch        *b_genphi;   //!
  TBranch        *b_genm;   //!

  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;

  if (!oniaTree) { cout<<"[ERROR] no tree found"<<endl; return;}

  fChain = oniaTree;
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
  if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
  if (fChain->GetBranch("ngen")) fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
  if (fChain->GetBranch("genpt")) fChain->SetBranchAddress("genpt", genpt, &b_genpt);
  if (fChain->GetBranch("geneta")) fChain->SetBranchAddress("geneta", geneta, &b_geneta);
  if (fChain->GetBranch("geny")) fChain->SetBranchAddress("geny", geny, &b_geny);
  if (fChain->GetBranch("genphi")) fChain->SetBranchAddress("genphi", genphi, &b_genphi);
  if (fChain->GetBranch("genm")) fChain->SetBranchAddress("genm", genm, &b_genm);

  Long64_t nentries =fChain->GetEntries();
  //nentries = 2000000;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
      //Long64_t ientry = LoadTree(jentry);
      //if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      fChain->GetEntry(jentry);
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  zed = -1;
	  drmin = 0.5;

	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  if (GenQQ4mom->M()<2.6 || GenQQ4mom->M()>3.5) continue;
	  if (GenQQ4mom->Pt() < 3 || GenQQ4mom->Pt() > 35) continue;
	  if (isMid && GenQQ4mom->Pt() < 6.5) continue;
	  if (abs(GenQQ4mom->Rapidity())>2.4) continue;
	  if (isMid && abs(GenQQ4mom->Rapidity())>1.6) continue;
	  if (!isMid && abs(GenQQ4mom->Rapidity())<1.6) continue;
	  
	  for (int iJet=0; iJet<ngen; iJet++)
	    {
	      if (genpt[iJet]<25 || genpt[iJet]>35) continue;
	      if (abs(geny[iJet])>2.4) continue;

	      TLorentzVector v_jet;
	      v_jet.SetPtEtaPhiM(genpt[iJet], geneta[iJet], genphi[iJet], genm[iJet]);
	      if (GenQQ4mom->DeltaR (v_jet)<=drmin)
		{
		  drmin = GenQQ4mom->DeltaR (v_jet);
		  zed = GenQQ4mom->Pt()*1.0/genpt[iJet];
		}
	      
	    }// end of gen jet loop
	  if (zed == -1) continue;
	  if (zed > 1 && zed <= 1.000001) zed = 0.9999999;
	  if (underflowOff && isMid && zed<0.4) continue;
	  if (underflowOff && !isMid && zed<0.2) continue;
	  zHist->Fill(zed);
	} //end of genQQ loop 
    }//end of events loop

  zHist->Scale(1.0/zHist->Integral("width"));
  //zHist->Scale(1.0/0.2);
  gSystem->mkdir("Output/MCResults");
  TFile *fsave = new TFile(Form("Output/MCResults/mcResult_%s_%s_%s.root",isPr?"prompt":"nonprompt", isMid?"mid":"fwd", underflowOff?"underflowOff":"all"),"RECREATE");
  zHist->Write("zDistMC");
  fsave->Close();
}
