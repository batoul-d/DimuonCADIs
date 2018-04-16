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

#include <vector> 
//#ifdef __MAKECINT__ 
//#pragma link C++ class vector<vector<float> >+; 
//#endif

using namespace std;

Double_t ptbins2D []= {3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50};
Double_t ybins2D []= {-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

void GetMCResults (bool isPr=true, bool underflowOff=true, bool plotBpt=false);
Bool_t isGlobalMuonInAccept2015 (TLorentzVector* Muon);

void GetMCResults_all ()
{
  GetMCResults(true, true, false);
  //GetMCResults(true, false, false);

  GetMCResults(false,  true, false);
  //GetMCResults(false,  false, false);

  GetMCResults(false, true, true);
  //GetMCResults(false, false, true);
}

void GetMCResults (bool isPr, bool underflowOff, bool plotBpt)
{
  int ny2D = sizeof(ybins2D)/sizeof(double)-1;
  int npt2D = sizeof(ptbins2D)/sizeof(double)-1;

  double sefer=0;

  TFile *treeFile = TFile::Open(Form("/data_CMS/cms/mnguyen/jPsiJet/mc/%s/merged_HiForestAOD.root", isPr?"prompt/v9_ext":"nonprompt/v9_ext_bHadronGen"));
  cout << Form("[INFO] Reading tree from file /data_CMS/cms/mnguyen/jPsiJet/mc/%s/merged_HiForestAOD.root", isPr?"prompt/v9_ext":"nonprompt/v9_ext_bHadronGen")<<endl;

  TTree* oniaTree = (TTree*) treeFile->Get("hionia/myTree");
  TTree* jetTree = (TTree*) treeFile->Get("ak4PFJetAnalyzer/t");
  TTree* bTree(0x0);
  bTree = (TTree*) treeFile->Get("bHadronAna/hi");
  if (!bTree) cout<<"b tree not found"<<endl;
  oniaTree->AddFriend(jetTree);
  if (plotBpt)  oniaTree->AddFriend(bTree);

  TH1F* zMidHist = new TH1F ("zMidHist","", 7, 0.02, 1); zMidHist->Sumw2();
  TH1F* zFwdHist = new TH1F ("zFwdHist","", 5, 0, 1); zFwdHist->Sumw2();
  TH1F* nMidHist = new TH1F ("nMidHist","", 10, 0, 1); nMidHist->Sumw2();
  TH1F* nFwdHist = new TH1F ("nFwdHist","", 10, 0, 1); nFwdHist->Sumw2();

  cout <<"[INFO] Importing AccFiles to correct"<<endl;
  TFile *corrFile = TFile::Open(Form("../MyEfficiency/FilesAccxEff/Acc/%sAccHists.root",isPr?"pr":"npr"));
  TH2F* corrNum = (TH2F*) corrFile->Get("hnum_2d_nominal");
  TH2F* corrDeno =(TH2F*) corrFile->Get("hdeno_2d");

  TEfficiency* accCorr = new TEfficiency("accCorr", "Acc(y,pt); y; pt; eff", ny2D, ybins2D, npt2D, ptbins2D);
  accCorr->SetStatisticOption(TEfficiency::kBBayesian);
  accCorr->SetPassedHistogram(*corrNum,"f");
  accCorr->SetTotalHistogram(*corrDeno,"f");

  double accWeight;
  double zed;
  double drmin;

  TTree* fChain;

  Int_t           Gen_QQ_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_QQ_mupl_4mom;
  TClonesArray    *Gen_QQ_mumi_4mom;
  Int_t           ngen;
  Float_t         genpt[99];   //[ngen]
  Float_t         geneta[99];   //[ngen]
  Float_t         geny[99];   //[ngen]
  Float_t         genphi[99];   //[ngen]
  Float_t         genm[99];   //[ngen]
  //vector<vector<float> >    *jtbHadronPt; //
  std::vector<float>   *pt =0;
  std::vector<float>   *eta =0;
  std::vector<float>   *phi = 0;
  Int_t           mult;


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
  //TBranch        *b_jtbHadronPt;
  TBranch        *b_mult;
  TBranch        *b_pt = 0;
  TBranch        *b_eta = 0;
  TBranch        *b_phi = 0;

  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;
  //jtbHadronPt = 0;

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
  if (fChain->GetBranch("mult")) fChain->SetBranchAddress("mult", &mult, &b_mult);
  if (fChain->GetBranch("pt")) fChain->SetBranchAddress("pt", &pt, &b_pt);
  if (fChain->GetBranch("eta")) fChain->SetBranchAddress("eta", &eta, &b_eta);
  if (fChain->GetBranch("phi")) fChain->SetBranchAddress("phi", &phi, &b_phi);
  
  
  Long64_t nentries =fChain->GetEntries();
  //nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<Form(" for (%s, %s, %s)",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"plot z(B)":"plot z(Jpsi)")<<endl;
      nb = fChain->GetEntry(jentry);   
      nbytes += fChain->GetEntry(jentry);
      fChain->GetEntry(jentry);
      
      for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	{
	  zed = -1;
	  drmin = 0.5;
	  TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	  TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);

	  if (GenQQ4mom->M()<2.6 || GenQQ4mom->M()>3.5) continue;
	  if (GenQQ4mom->Pt() < 3 || GenQQ4mom->Pt() > 35) continue;
	  if (abs(GenQQ4mom->Rapidity())>2.4) continue;
	  if (!isGlobalMuonInAccept2015(GenQQmupl) || !isGlobalMuonInAccept2015(GenQQmumi)) continue;

	  accWeight = 1.0/(accCorr->GetEfficiency(accCorr->FindFixBin(GenQQ4mom->Rapidity(), GenQQ4mom->Pt())));

	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    nFwdHist->Fill(sefer, accWeight);
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt() > 6.5)
	    nMidHist->Fill(sefer, accWeight);

	  for (int iJet=0; iJet<ngen; iJet++)
	    {
	      if (genpt[iJet]<25 || genpt[iJet]>35) continue;
	      if (abs(geny[iJet])>2.4) continue;

	      TLorentzVector v_jet;
	      v_jet.SetPtEtaPhiM(genpt[iJet], geneta[iJet], genphi[iJet], genm[iJet]);
	      if (GenQQ4mom->DeltaR (v_jet)<=drmin)
		{
		  drmin = GenQQ4mom->DeltaR (v_jet);
		  if (!plotBpt)
		    zed = GenQQ4mom->Pt()*1.0/genpt[iJet];
		  if (plotBpt) 
		    {
		      int ibestB = -1;
		      float drminb = 999;

		      for (int ib = 0; ib<mult; ib++) {
			float deta = eta->at(ib);
			deta = eta->at(ib) - geneta[iJet];
			float dphi = phi->at(ib) - genphi[iJet];
			dphi = acos(cos(dphi));
			float dR = sqrt(deta*deta + dphi*dphi);
			if (dR < drminb) 
			  {
			    drminb = dR;
			    ibestB = ib;
			  }
		      }// end of b loop 
		      if (ibestB > 0)
			zed = pt->at(ibestB)*1.0/genpt[iJet];
		      if (ibestB > 0 && drminb>0.4) cout <<"[WARNING] Bbig dR(b-jet) = "<<drminb<<"for pt(b) = "<<pt->at(ibestB)<<endl;
		    }// end of b cond
		}
	      
	    }// end of gen jet loop
	  if (zed == -1) continue;
	  if (zed > 1 && zed <= 1.000001) zed = 0.9999999;
	  if (underflowOff  && zed<0.2) continue;
	  if (abs(GenQQ4mom->Rapidity())>1.6)
	    zFwdHist->Fill(zed, accWeight);
	  if (underflowOff  && zed<0.44) continue;
	  if (abs(GenQQ4mom->Rapidity())<1.6 && GenQQ4mom->Pt()>6.5)
	    zMidHist->Fill(zed, accWeight);
	} //end of genQQ loop 
    }//end of events loop
  //if (underflowOff)
  //{
  //zMidHist->Scale(1.0/zMidHist->Integral("width"));
  //zFwdHist->Scale(1.0/zFwdHist->Integral("width"));
  //}
  //zHist->Scale(1.0/0.2);
  gSystem->mkdir("Output/MCResults");
  TFile *fsave = new TFile(Form("Output/MCResults/mcResult_%s_%s%s.root",isPr?"prompt":"nonprompt", underflowOff?"underflowOff":"all", plotBpt?"_bHadronPt":""),"RECREATE");
  zMidHist->Write("zDist_mid");
  nMidHist->Write("Ntot_mid");
  zFwdHist->Write("zDist_fwd");
  nFwdHist->Write("Ntot_fwd");
  fsave->Close();
}

Bool_t isGlobalMuonInAccept2015 (TLorentzVector* Muon)
{
  return (fabs(Muon->Eta()) < 2.4 &&
	  ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
	   (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
	   (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
}
