#define myTree_cxx
#define _USE_MATH_DEFINES
#include "myTree.h"
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

#include "tnp_weight.h"
using namespace std;
using namespace  RooFit;
//using namespace RooPlot;

TLorentzVector* matchReco;
TLorentzVector* matchGen;
int nj=0;
int nt=0;
int np=0;
int nn=0;
int ntp=0;
int ntn=0;
Float_t jpsi_m;
Float_t jpsi_pt;
Float_t jpsi_eta;
Float_t jpsi_phi;
Float_t jpsi_rap;
Float_t dr;
Float_t dphi;
Float_t dphimin;
Float_t deta;
Float_t drmin; 
Float_t z=100;
int triggerIndex_PP =0;
int k;
Float_t weight;
Float_t tnp_weight;
Double_t ptbins []={3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 24, 27, 30, 35, 40, 45, 50};
int nptbins = ((sizeof(ptbins)/sizeof(double))-1);
Double_t etabins []={-2.4, -2, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
int netabins = ((sizeof(etabins)/sizeof(double))-1);

Double_t ybins []={0, 0.4, 0.8, 1.2, 1.6, 2, 2.4};
int nybins = ((sizeof(ybins)/sizeof(double))-1);

RooRealVar* mass = new RooRealVar("Mass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
RooRealVar* genmass = new RooRealVar("genMass","#mu#mu mass", 2.6, 3.4, "GeV/c^{2}");
RooRealVar* ctau = new RooRealVar("ctau","c_{#tau}", -2.0, 5.0, "mm");
RooRealVar* genctau = new RooRealVar("genctau","gen c_{#tau}", -2.0, 5.0, "mm");
RooRealVar* zed = new RooRealVar("z", "z_{J/#psi}", 0, 1);
RooRealVar* r = new RooRealVar("r", "#DeltaR", 0, 0.5);
RooRealVar* genzed = new RooRealVar("genz", "genz_{J/#psi}", 0, 1);
RooRealVar* w = new RooRealVar("w","MC weight", 0.0, 10000000.0);
RooRealVar* genr = new RooRealVar("genr", "#genDeltaR", 0, 0.5);
RooRealVar* pt = new RooRealVar("pt", "pt", 6.5, 25.5);
RooRealVar* rap = new RooRealVar("rap", "rapidity", -2.4, 2.4);
RooRealVar* gen_pt = new RooRealVar("gen_pt", "pt for gen j/psi", 6.5, 25.5);
RooRealVar* gen_rap = new RooRealVar("gen_rap", "rapidity for gen j/psi", -2.4, 2.4);

RooArgSet* recoset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);
RooArgSet*  unwset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap);
RooArgSet*  genset = new RooArgSet (*genmass, *genctau, *genzed, *genr, *gen_pt, *gen_rap);
RooArgSet* prset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);
RooArgSet* nprset = new RooArgSet(*mass, *ctau, *zed, *r, *pt, *rap, *w);

RooDataSet* data = new RooDataSet ("data", "data for reconstructed J/#psi", *recoset, WeightVar(*w));
RooDataSet* unwdata = new RooDataSet ("unwdata", "unweighted data", *unwset);
RooDataSet* gendata = new RooDataSet ("gendata", "data for generated J/#psi", *genset);
RooDataSet* prdata = new RooDataSet ("prdata", "data for prompt J/#psi", *prset, WeightVar(*w));
RooDataSet* nprdata = new RooDataSet ("nprdata", "data for non prompt J/#psi", *nprset, WeightVar(*w));

void myTree::EffCalc()
{
  if (isMc)
    {
      TH1F* ptg = new TH1F ("ptg", "N_{gen} vs p_{T}; p_{T}; N_{total}", nptbins, ptbins);
      TH1F* ptr = new TH1F ("ptr", "N_{reco} vs p_{T}; p_{T}; N_{reco}", nptbins, ptbins);
      TH1F* rapr = new TH1F ("rapr","N_{reco} vs #eta; #eta; N_{reco}", netabins, etabins);
      TH1F* rapg = new TH1F ("rapg","N_{gen} vs #eta; #eta; N_{total}", netabins, etabins);
      TH2F* ptrapg = new TH2F ("ptrapg", "N_{gen} vs p_{T} and y; y; p_{T}; N_{total}", netabins, etabins, nptbins, ptbins);
      TH2F* ptrapr = new TH2F ("ptrapr", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);

      TH2F* num_binned = new TH2F ("num_", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_plus1sig = new TH2F ("num_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_minus1sig = new TH2F ("num_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid_sta = new TH2F ("num_muid_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid = new TH2F ("num_muid", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid_plus1sig = new TH2F ("num_muid_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_muid_minus1sig = new TH2F ("num_muid_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_sta = new TH2F ("num_sta", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_sta_plus1sig = new TH2F ("num_sta_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_sta_minus1sig = new TH2F ("num_sta_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_trk_plus1sig = new TH2F ("num_trk_plus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);
      TH2F* num_trk_minus1sig = new TH2F ("num_trk_minus1sig", "N_{reco} vs p_{T} and y; y; p_{T}; N_{reco}", netabins, etabins, nptbins, ptbins);

      ptg->Sumw2(); ptr->Sumw2(); rapr->Sumw2(); rapg->Sumw2(); ptrapg->Sumw2(); ptrapr->Sumw2(); num_binned->Sumw2(); num_plus1sig->Sumw2(); num_minus1sig->Sumw2(); num_muid_sta->Sumw2(); num_muid->Sumw2(); num_muid_plus1sig->Sumw2(); num_muid_minus1sig->Sumw2(); num_sta->Sumw2(); num_sta_plus1sig->Sumw2(); num_sta_minus1sig->Sumw2(); num_trk_plus1sig->Sumw2(); num_trk_minus1sig->Sumw2();

      Long64_t nentries =fChain->GetEntries();
      //nentries = 2000000;
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m=GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();
	      if (jpsi_pt>3 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
		  ptg->Fill(jpsi_pt);
		  rapg->Fill(jpsi_rap);
		  ptrapg->Fill(jpsi_rap, jpsi_pt);
		}
	    }

	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
		  TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m=RecoQQ4mom->M();
		  if (
		      jpsi_pt > 3  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) && // if it matches the trigger 
		      (isMatchedRecoDiMuon(iQQ))
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  /////nominal
			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  ptr->Fill(jpsi_pt,tnp_weight);
			  rapr->Fill(jpsi_rap,tnp_weight);
			  ptrapr->Fill(jpsi_rap, jpsi_pt, tnp_weight);

			  /////systematics
			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-10) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-10) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_binned->Fill(jpsi_rap, jpsi_pt, tnp_weight);

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * 
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid_sta->Fill(jpsi_rap, jpsi_pt, tnp_weight);
 
			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_muid_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_muid_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_muid_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_sta->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-1) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-1) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_sta_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_sta_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),-2) * tnp_weight_sta_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),-2) *
			    tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
			  num_sta_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(-1) * tnp_weight_trk_pp(-1);
			  num_trk_plus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight); 

			  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) *
			    tnp_weight_trk_pp(-2) * tnp_weight_trk_pp(-2);
			  num_trk_minus1sig->Fill(jpsi_rap, jpsi_pt, tnp_weight);
			}
		    }
		}
	    }
	}
      TEfficiency* gptef = new TEfficiency("gptef", "AccxEff(pt)", nptbins, ptbins);
      gptef->SetStatisticOption(TEfficiency::kBBayesian);
      gptef->SetPassedHistogram(*ptr,"f");
      gptef->SetTotalHistogram(*ptg,"f");

      TEfficiency* grapef = new TEfficiency("grapef", "AccxEff(y)", netabins, etabins);
      grapef->SetStatisticOption(TEfficiency::kBBayesian);
      grapef->SetPassedHistogram(*rapr,"f");
      grapef->SetTotalHistogram(*rapg,"f");

      /////nominal 2D
      TEfficiency* eff_nom = new TEfficiency("eff_nom", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_nom->SetStatisticOption(TEfficiency::kBBayesian);
      eff_nom->SetPassedHistogram(*ptrapr,"f");
      eff_nom->SetTotalHistogram(*ptrapg,"f");
      eff_nom->SetName("hcorr_Jpsi_PP");

      ////systematics
      TEfficiency* eff_binned = new TEfficiency("eff_binned", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_binned->SetStatisticOption(TEfficiency::kBBayesian);
      eff_binned->SetPassedHistogram(*num_binned,"f");
      eff_binned->SetTotalHistogram(*ptrapg,"f");
      eff_binned->SetName("hcorr_binned");

      TEfficiency* eff_plus1sig = new TEfficiency("eff_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_plus1sig->SetPassedHistogram(*num_plus1sig,"f");
      eff_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_plus1sig->SetName("hcorr_plus1sig");

      TEfficiency* eff_minus1sig = new TEfficiency("eff_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_minus1sig->SetPassedHistogram(*num_minus1sig,"f");
      eff_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_minus1sig->SetName("hcorr_minus1sig");

      TEfficiency* eff_muid_sta = new TEfficiency("eff_muid_sta", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid_sta->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid_sta->SetPassedHistogram(*num_muid_sta,"f");
      eff_muid_sta->SetTotalHistogram(*ptrapg,"f");
      eff_muid_sta->SetName("hcorr_muid_sta");

      TEfficiency* eff_muid = new TEfficiency("eff_muid", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid->SetPassedHistogram(*num_muid,"f");
      eff_muid->SetTotalHistogram(*ptrapg,"f");
      eff_muid->SetName("hcorr_muid");

      TEfficiency* eff_muid_plus1sig = new TEfficiency("eff_muid_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid_plus1sig->SetPassedHistogram(*num_muid_plus1sig,"f");
      eff_muid_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_muid_plus1sig->SetName("hcorr_muid_plus1sig");

      TEfficiency* eff_muid_minus1sig = new TEfficiency("eff_muid_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_muid_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_muid_minus1sig->SetPassedHistogram(*num_muid_minus1sig,"f");
      eff_muid_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_muid_minus1sig->SetName("hcorr_muid_minus1sig");

      TEfficiency* eff_sta = new TEfficiency("eff_sta", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_sta->SetStatisticOption(TEfficiency::kBBayesian);
      eff_sta->SetPassedHistogram(*num_sta,"f");
      eff_sta->SetTotalHistogram(*ptrapg,"f");
      eff_sta->SetName("hcorr_sta");

      TEfficiency* eff_sta_plus1sig = new TEfficiency("eff_sta_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_sta_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_sta_plus1sig->SetPassedHistogram(*num_sta_plus1sig,"f");
      eff_sta_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_sta_plus1sig->SetName("hcorr_sta_plus1sig");

      TEfficiency* eff_sta_minus1sig = new TEfficiency("eff_sta_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_sta_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_sta_minus1sig->SetPassedHistogram(*num_sta_minus1sig,"f");
      eff_sta_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_sta_minus1sig->SetName("hcorr_sta_minus1sig");

      TEfficiency* eff_trk_plus1sig = new TEfficiency("eff_trk_plus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_trk_plus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_plus1sig->SetPassedHistogram(*num_trk_plus1sig,"f");
      eff_trk_plus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_trk_plus1sig->SetName("hcorr_trk_plus1sig");

      TEfficiency* eff_trk_minus1sig = new TEfficiency("eff_trk_minus1sig", "AccxEff(y,pt); y; pt; eff", netabins, etabins, nptbins, ptbins);
      eff_trk_minus1sig->SetStatisticOption(TEfficiency::kBBayesian);
      eff_trk_minus1sig->SetPassedHistogram(*num_trk_minus1sig,"f");
      eff_trk_minus1sig->SetTotalHistogram(*ptrapg,"f");
      eff_trk_minus1sig->SetName("hcorr_trk_minus1sig");


      TFile* ef (0x0);
      if (isPr)
	ef = new TFile ("../Fitter/Input/pr_correction_AccEff.root", "RECREATE");
      else
	ef = new TFile ("../Fitter/Input/npr_correction_AccEff.root", "RECREATE");


      gptef->Write("effVsPt");
      grapef->Write("effVsRap");
      ptrapg->Write("hcorr_his_deno");
      ptrapr->Write("hcorr_his_num");
      eff_nom->Write("hcorr_Jpsi_PP");

      eff_binned->Write("hcorr_binned"); 
      eff_plus1sig->Write("hcorr_plus1sig"); 
      eff_minus1sig->Write("hcorr_minus1sig"); 
      eff_muid_sta->Write("hcorr_muid_sta"); 
      eff_muid->Write("hcorr_muid"); 
      eff_muid_plus1sig->Write("hcorr_muid_plus1sig"); 
      eff_muid_minus1sig->Write("hcorr_muid_minus1sig"); 
      eff_sta->Write("hcorr_sta"); 
      eff_sta_plus1sig->Write("hcorr_sta_plus1sig"); 
      eff_sta_minus1sig->Write("hcorr_sta_minus1sig"); 
      eff_trk_plus1sig->Write("hcorr_trk_plus1sig"); 
      eff_trk_minus1sig->Write("hcorr_trk_minus1sig");

      ef->Close();
    }
  else 
    cout<< "[ERROR] this is data and not MC"<<endl;
}

void myTree::ANEffPlots()
{
  if (isMc)
    {
      TH1F* Accpthist = new TH1F ("Accpthist", ";p_{T}^{#mu#mu} (GeV/c); Acceptance",nptbins,ptbins);
      TH1F* Effpthist = new TH1F("Effpthist", ";p_{T}^{#mu#mu} (GeV/c); Efficiency",nptbins,ptbins);
      TH1F* Accyhist = new TH1F("Accyhist", ";|y|^{#mu#mu}; Acceptance",nybins,ybins);
      TH1F* Effyhist = new TH1F("Effyhist", ";|y|^{#mu#mu}; Efficiency",nybins,ybins);
      TH1F* AccEffpthist = new TH1F("AccEffpthist", ";p_{T}^{#mu#mu} (GeV/c); Acc x Eff",nptbins,ptbins);
      TH1F* AccEffyhist = new TH1F("AccEffyhist", ";|y|^{#mu#mu}; Acc x Eff",nybins,ybins);

      TCanvas* c = new TCanvas("c","",1000,800);
      TH1F* Acc_pt_016_num = new TH1F ("Acc_pt_016_num", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TH1F* Acc_pt_016_deno = new TH1F ("Acc_pt_016_deno", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TH1F* Acc_pt_1624_num = new TH1F ("Acc_pt_1624_num", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TH1F* Acc_pt_1624_deno = new TH1F ("Acc_pt_1624_deno", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TH1F* Acc_y_num = new TH1F ("Acc_y_num", ";|y|^{#mu#mu}; Acceptance", nybins, ybins);
      TH1F* Acc_y_deno = new TH1F ("Acc_y_deno", ";|y|^{#mu#mu}; Acceptance", nybins, ybins);
      TH2F* Acc_ypt_num = new TH2F ("Acc_ypt_num", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      TH2F* Acc_ypt_deno = new TH2F ("Acc_ypt_deno", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      
      TH1F* Eff_pt_016_num = new TH1F ("Eff_pt_016_num", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TH1F* Eff_pt_016_deno = new TH1F ("Eff_pt_016_deno", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TH1F* Eff_pt_1624_num = new TH1F ("Eff_pt_1624_num", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TH1F* Eff_pt_1624_deno = new TH1F ("Eff_pt_1624_deno", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TH1F* Eff_y_num = new TH1F ("Eff_y_num", ";|y|^{#mu#mu}; Efficiency", nybins, ybins);
      TH1F* Eff_y_deno = new TH1F ("Eff_y_deno", ";|y|^{#mu#mu}; Efficiency", nybins, ybins);
      TH2F* Eff_ypt_num = new TH2F ("Eff_ypt_num", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      TH2F* Eff_ypt_deno = new TH2F ("Eff_ypt_deno", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      
      TEfficiency* Acc_pt_016 = new TEfficiency ("Acc_pt_016", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TEfficiency* Acc_pt_1624 = new TEfficiency ("Acc_pt_1624", ";p_{T}^{#mu#mu} (GeV/c); Acceptance", nptbins, ptbins);
      TEfficiency* Acc_y = new TEfficiency ("Acc_y", ";y^{#mu#mu}; Acceptance", nybins, ybins);
      TEfficiency* Acc_ypt = new TEfficiency ("Acc_ypt", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      
      TEfficiency* Eff_pt_016 = new TEfficiency ("Eff_pt_016", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TEfficiency* Eff_pt_1624 = new TEfficiency ("Eff_pt_1624", ";p_{T}^{#mu#mu} (GeV/c); Efficiency", nptbins, ptbins);
      TEfficiency* Eff_y = new TEfficiency ("Eff_y", ";y^{#mu#mu}; Efficiency", nybins, ybins);
      TEfficiency* Eff_ypt = new TEfficiency ("Eff_ypt", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      
      TEfficiency* AccEff_pt_016 = new TEfficiency ("AccEff_pt_016", ";p_{T}^{#mu#mu} (GeV/c); Acc x Eff", nptbins, ptbins);
      TEfficiency* AccEff_pt_1624 = new TEfficiency ("AccEff_pt_1624", ";p_{T}^{#mu#mu} (GeV/c); Acc x Eff", nptbins, ptbins);
      TEfficiency* AccEff_y = new TEfficiency ("AccEff_y", ";y^{#mu#mu}; Acc x Eff", nybins, ybins);
      TEfficiency* AccEff_ypt = new TEfficiency ("AccEff_ypt", ";|y|^{#mu#mu};p_{T}^{#mu#mu} (GeV/c)", nybins, ybins, nptbins, ptbins);
      
      Long64_t nentries = fChain->GetEntries();
      //nentries = 2000000;
      
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++)
	{
	  if (jentry%1000000==0) cout<<"[INFO] "<<jentry<<"/"<<nentries<<endl;
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      if ( abs(GenQQ4mom->Rapidity()) < 2.4 && abs (GenQQ4mom->Rapidity()) > 1.6 && GenQQ4mom->Pt() > 3)
		{
		  Acc_pt_1624_deno->Fill(GenQQ4mom->Pt());
		  if(areGenMuonsInAcceptance2015(iQQ))
		    {
		      Acc_pt_1624_num->Fill(GenQQ4mom->Pt());
		      Eff_pt_1624_deno->Fill(GenQQ4mom->Pt());
		    }
		}
	      if (abs(GenQQ4mom->Rapidity()) < 1.6 && GenQQ4mom->Pt() > 6.5)
		{
		  Acc_pt_016_deno->Fill(GenQQ4mom->Pt());
		  if(areGenMuonsInAcceptance2015(iQQ))
		    {
		      Acc_pt_016_num->Fill(GenQQ4mom->Pt());
		      Eff_pt_016_deno->Fill(GenQQ4mom->Pt());
		    }
		}
	    if (abs(GenQQ4mom->Rapidity()) < 2.4 && GenQQ4mom->Pt() > 6.5)
	      {
		Acc_y_deno->Fill(abs(GenQQ4mom->Rapidity()));
		if(areGenMuonsInAcceptance2015(iQQ))
		  {
		    Acc_y_num->Fill(abs(GenQQ4mom->Rapidity()));
                    Eff_y_deno->Fill(abs(GenQQ4mom->Rapidity()));
		  }
	      }
	    if (abs(GenQQ4mom->Rapidity()) < 2.4 && GenQQ4mom->Pt() > 3)
	      {
	       	Acc_ypt_deno->Fill(abs(GenQQ4mom->Rapidity()), GenQQ4mom->Pt());
		if (areGenMuonsInAcceptance2015(iQQ))
		  {
		    Acc_ypt_num->Fill(abs(GenQQ4mom->Rapidity()), GenQQ4mom->Pt());
		    Eff_ypt_deno->Fill(abs(GenQQ4mom->Rapidity()), GenQQ4mom->Pt());
		  }
	      }
	  }
	  for(int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
	    {
	      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	      TLorentzVector *RecoQQmupl4mom = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iQQ);
	      TLorentzVector *RecoQQmumi4mom = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iQQ);
	      if(
		 (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts                                                                                                                
		 (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts                                                                                                                  
		 (isTriggerMatch(iQQ, triggerIndex_PP)) && // if it matches the trigger                                                                                                                
		 (isMatchedRecoDiMuon(iQQ))
		 )
		{
		  tnp_weight = tnp_weight_trg_pp(RecoQQmupl4mom->Pt(),RecoQQmupl4mom->Eta(),0) * tnp_weight_trg_pp(RecoQQmumi4mom->Pt(),RecoQQmumi4mom->Eta(),0) * tnp_weight_trk_pp(0) * tnp_weight_trk_pp(0);
		  if (abs(RecoQQ4mom->Rapidity()) < 2.4 && abs (RecoQQ4mom->Rapidity()) > 1.6 && RecoQQ4mom->Pt() > 3)
		    Eff_pt_1624_num->Fill(RecoQQ4mom->Pt(), tnp_weight);
		  if(abs(RecoQQ4mom->Rapidity()) < 1.6 && RecoQQ4mom->Pt() > 6.5)
		    Eff_pt_016_num->Fill(RecoQQ4mom->Pt(), tnp_weight);
		  if(abs(RecoQQ4mom->Rapidity()) < 2.4 && RecoQQ4mom->Pt() > 6.5)
		    Eff_y_num->Fill(abs(RecoQQ4mom->Rapidity()),tnp_weight);
		  if(abs(RecoQQ4mom->Rapidity()) < 2.4 && RecoQQ4mom->Pt() > 3)
		    Eff_ypt_num->Fill(abs(RecoQQ4mom->Rapidity()), RecoQQ4mom->Pt(),tnp_weight);
		}
	    }
	}
      
      if(TEfficiency::CheckConsistency(*Acc_pt_1624_num,*Acc_pt_1624_deno))
	Acc_pt_1624 = new TEfficiency (*Acc_pt_1624_num,*Acc_pt_1624_deno);
      
      if(TEfficiency::CheckConsistency(*Acc_pt_016_num,*Acc_pt_016_deno))
	Acc_pt_016 = new TEfficiency (*Acc_pt_016_num,*Acc_pt_016_deno);
      
      if(TEfficiency::CheckConsistency(*Acc_y_num,*Acc_y_deno))
	Acc_y = new TEfficiency (*Acc_y_num,*Acc_y_deno);
      
      if(TEfficiency::CheckConsistency(*Acc_ypt_num,*Acc_ypt_deno))
	Acc_ypt = new TEfficiency (*Acc_ypt_num,*Acc_ypt_deno);

      if(TEfficiency::CheckConsistency(*Eff_pt_1624_num,*Eff_pt_1624_deno))
	Eff_pt_1624 = new TEfficiency (*Eff_pt_1624_num,*Eff_pt_1624_deno);
      
      if(TEfficiency::CheckConsistency(*Eff_pt_016_num,*Eff_pt_016_deno))
	Eff_pt_016 = new TEfficiency (*Eff_pt_016_num,*Eff_pt_016_deno);

      if(TEfficiency::CheckConsistency(*Eff_y_num,*Eff_y_deno))
	Eff_y = new TEfficiency (*Eff_y_num,*Eff_y_deno);

      if(TEfficiency::CheckConsistency(*Eff_ypt_num,*Eff_ypt_deno))
	Eff_ypt = new TEfficiency (*Eff_ypt_num,*Eff_ypt_deno);

      if(TEfficiency::CheckConsistency(*Eff_pt_1624_num,*Acc_pt_1624_deno))
        AccEff_pt_1624 = new TEfficiency (*Eff_pt_1624_num,*Acc_pt_1624_deno);

      if(TEfficiency::CheckConsistency(*Eff_pt_016_num,*Acc_pt_016_deno))
        AccEff_pt_016 = new TEfficiency (*Eff_pt_016_num,*Acc_pt_016_deno);

      if(TEfficiency::CheckConsistency(*Eff_y_num,*Acc_y_deno))
        AccEff_y = new TEfficiency (*Eff_y_num,*Acc_y_deno);

      if(TEfficiency::CheckConsistency(*Eff_ypt_num,*Acc_ypt_deno))
        AccEff_ypt = new TEfficiency (*Eff_ypt_num,*Acc_ypt_deno);

      TFile* Corr = new TFile(Form("ANPlots/%sAccEffCorr.root", isPr?"pr":"npr"),"RECREATE");
      Corr->cd();
      Acc_pt_1624->Write("Acc_pt_1624");
      Acc_pt_016->Write("Acc_pt_016");
      Acc_y->Write("Acc_y");
      Acc_ypt->Write("Acc_ypt");
      Eff_pt_1624->Write("Eff_pt_1624");
      Eff_pt_016->Write("Eff_pt_016");
      Eff_y->Write("Eff_y");
      Eff_ypt->Write("Eff_ypt");
      AccEff_pt_1624->Write("AccEff_pt_1624");
      AccEff_pt_016->Write("AccEff_pt_016");
      AccEff_y->Write("AccEff_y");
      AccEff_ypt->Write("AccEff_ypt");
      Corr->Close();
      
      c->cd();
      gStyle->SetOptStat(0);
      //Acc_pt_1624->SetMinimum(0);
      //Acc_pt_1624->SetMaximum(1.2);
      //Acc_pt_016->SetMinimum(0);
      //Acc_pt_016->SetMaximum(1.2);
      //Acc_y->SetMinimum(0);
      //Acc_y->SetMaximum(1.2);
      //Eff_pt_1624->SetMinimum(0);
      //Eff_pt_1624->SetMaximum(1.2);
      //Eff_pt_016->SetMinimum(0);
      //Eff_pt_016->SetMaximum(1.2);
      //Eff_y->SetMinimum(0);
      //Eff_y->SetMaximum(1.2);

      Accpthist->SetMinimum(0);
      Accpthist->SetMaximum(1.2);
      Effpthist->SetMinimum(0);
      Effpthist->SetMaximum(1.2);
      Accyhist->SetMinimum(0);
      Accyhist->SetMaximum(1.2);
      Effyhist->SetMinimum(0);
      Effyhist->SetMaximum(1.2);
      AccEffpthist->SetMinimum(0);
      AccEffpthist->SetMaximum(1.2);
      AccEffyhist->SetMinimum(0);
      AccEffyhist->SetMaximum(1.2);

      TLine* pt1 = new TLine(3,1,50,1);
      pt1->SetLineColor(kRed);
      pt1->SetLineStyle(2);
      TLine* y1 = new TLine(0,1,2.4,1);
      y1->SetLineColor(kRed);
      y1->SetLineStyle(2);
      
      Acc_pt_1624->SetMarkerStyle(27);
      Acc_pt_1624->SetMarkerColor(40);
      Acc_pt_1624->SetMarkerSize(1);
      Acc_pt_1624->SetLineColor(39);
      Acc_pt_016->SetMarkerStyle(30);
      Acc_pt_016->SetMarkerColor(46);
      Acc_pt_016->SetMarkerSize(1);
      Acc_pt_016->SetLineColor(48);
      Acc_y->SetMarkerStyle(26);
      Acc_y->SetMarkerColor(8);
      Acc_y->SetMarkerSize(1);
      Acc_y->SetLineColor(29);

      Eff_pt_1624->SetMarkerStyle(27);
      Eff_pt_1624->SetMarkerColor(40);
      Eff_pt_1624->SetMarkerSize(1);
      Eff_pt_1624->SetLineColor(39);
      Eff_pt_016->SetMarkerStyle(30);
      Eff_pt_016->SetMarkerColor(46);
      Eff_pt_016->SetMarkerSize(1);
      Eff_pt_016->SetLineColor(48);
      Eff_y->SetMarkerStyle(26);
      Eff_y->SetMarkerColor(8);
      Eff_y->SetMarkerSize(1);
      Eff_y->SetLineColor(29);
      
      AccEff_pt_1624->SetMarkerStyle(27);
      AccEff_pt_1624->SetMarkerColor(40);
      AccEff_pt_1624->SetMarkerSize(1);
      AccEff_pt_1624->SetLineColor(39);
      AccEff_pt_016->SetMarkerStyle(30);
      AccEff_pt_016->SetMarkerColor(46);
      AccEff_pt_016->SetMarkerSize(1);
      AccEff_pt_016->SetLineColor(48);
      AccEff_y->SetMarkerStyle(26);
      AccEff_y->SetMarkerColor(8);
      AccEff_y->SetMarkerSize(1);
      AccEff_y->SetLineColor(29);

      TLegend *leg = new TLegend(0.6, 0.775, 0.8, 0.875);
      leg->AddEntry(Acc_pt_016, "|y| < 1.6", "lep");
      leg->AddEntry(Acc_pt_1624, "1.6 < |y| < 2.4", "lep");
      leg->SetBorderSize(1);
      leg->SetFillStyle(0);
      
      TPaveText *tbox = new TPaveText(0.15,0.8,0.35,0.9, "BRNDC");
      tbox->AddText(Form("%s J/#psi",isPr?"Prompt":"Nonprompt"));
      tbox->SetBorderSize(0);
      tbox->SetFillColor(0);
      tbox->SetFillStyle(0);

      Accpthist->Draw();    
      Acc_pt_1624->Draw("same");
      Acc_pt_016->Draw("same");
      leg->Draw("same");
      tbox->Draw("same");
      pt1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccVsPt.pdf", isPr?"pr":"npr"));

      Effpthist->Draw();
      Eff_pt_1624->Draw("same");
      Eff_pt_016->Draw("same");
      leg->Draw("same");
      tbox->Draw("same");
      pt1->Draw("same");
      c->SaveAs(Form("ANPlots/%sEffVsPt.pdf", isPr?"pr":"npr"));
      
      AccEffpthist->Draw();
      AccEff_pt_1624->Draw("same");
      AccEff_pt_016->Draw("same");
      leg->Draw("same");
      tbox->Draw("same");
      pt1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccEffVsPt.pdf", isPr?"pr":"npr"));

      tbox->AddText("6.5 < p_{T} < 50 GeV/c");
      
      Accyhist->Draw();
      Acc_y->Draw("same");
      tbox->Draw("same");
      y1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccVsY.pdf", isPr?"pr":"npr"));

      Effyhist->Draw();
      Eff_y->Draw("same");
      tbox->Draw("same");
      y1->Draw("same");
      c->SaveAs(Form("ANPlots/%sEffVsY.pdf", isPr?"pr":"npr"));
      
      AccEffyhist->Draw();
      AccEff_y->Draw("same");
      tbox->Draw("same");
      y1->Draw("same");
      c->SaveAs(Form("ANPlots/%sAccEffVsY.pdf", isPr?"pr":"npr"));

      Acc_ypt->Draw("COLZ");
      c->SaveAs(Form("ANPlots/%sAcc2D.pdf", isPr?"pr":"npr"));
      Eff_ypt->Draw("COLZ");
      c->SaveAs(Form("ANPlots/%sEff2D.pdf", isPr?"pr":"npr"));
      AccEff_ypt->Draw("COLZ");
      c->SaveAs(Form("ANPlots/%sAccEff2D.pdf", isPr?"pr":"npr"));
    }
  else {
    cout << "[ERROR] this is data not MC" << endl;
  }
}
void myTree::ClosureTest()
{
  TH1F* rpt = new TH1F ("rpt", "pt distribution at reco level", 10, 6.5, 26.5);
  TH1F* gpt = new TH1F ("gpt", "pt distribution at gen level", 10, 6.5, 26.5);
  TFile* f (0x0);
  if (isPr)
    f= TFile::Open("prEff.root");
  else
    f= TFile::Open ("nprEff.root");
  TEfficiency* eff = (TEfficiency*) f->Get("ptrap");

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  if (isMc)
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m = GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
		  gpt->Fill(jpsi_pt);
		}
	    }
	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  weight=1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			  rpt->Fill(jpsi_pt, weight);
			}
		    }
		}
	    }
	}
      TFile* testfile (0x0);
      if (isPr)
	testfile = new TFile ("prClosureTest.root","RECREATE");
      else
	testfile = new TFile ("nprClosureTest.root","RECREATE");
      rpt->Write("recopt");
      gpt->Write("genpt");
      testfile->Close();
    }
  else
    cout<< "this is data and not MC"<<endl;
}

void myTree::Loop()
{
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas ("c1", "", 1000, 800);
  TH1F* tot = new TH1F ("tot", ";p_{t}(#mu#mu)(GeV/c);Events", 25, ptbins);
  TH1F* jpjt = new TH1F ("jpjt", ";p_{t}(#mu#mu/jet)(GeV/c);Events", 25, ptbins);
  TH1F* ptot = new TH1F ("ptot", ";p_{t}(pr #mu#mu)(GeV/c);Events", 25, ptbins);
  TH1F* ntot = new TH1F ("ntot", ";p_{t}(npr #mu#mu)(GeV/c);Events", 25, ptbins);
  TH1F* pjpjt = new TH1F ("pjpjt", ";p_{t}(pr #mu#mu/jet)(GeV/c);Events", 25, ptbins);
  TH1F* njpjt = new TH1F ("njpjt", ";p_{t}(npr #mu#mu/jet)(GeV/c);Events", 25, ptbins);
  TF1* lcut = new TF1 ("lcut", "0.013+0.25/x", 6.5, 50);
  if (fChain == 0) return;
  //genzed->setBins(10);
  //mass->setBins(40);
  //zed->setBins(10);
  //rap->setBins(12);
  //gen_rap->setBins(12);
  Bool_t ismatch;
  TFile* f (0x0);
  TFile* prf = TFile::Open("prEff.root");
  TFile* nprf = TFile::Open ("nprEff.root");

  TEfficiency* preff = (TEfficiency*) prf->Get("ptrap");
  TEfficiency* npreff = (TEfficiency*) prf->Get("ptrap");
  TEfficiency* eff(0x0);
  if (isPr)
    eff=preff;
  else
    eff=npreff;

  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  if (isMc)
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  ismatch=false;
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for (int iQQ=0; iQQ<Gen_QQ_size; iQQ++)
	    {
	      TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	      jpsi_m=GenQQ4mom->M();
	      jpsi_pt = GenQQ4mom->Pt();
	      jpsi_rap = GenQQ4mom->Rapidity();
	      drmin = 1000;

	      if (jpsi_pt>6.5 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5)
		{
  
		  for (Long64_t ijet=0; ijet<ngen; ijet++)
		    {
		      if (abs(geny[ijet])<2.4 && genpt[ijet]>20)
			{
			  TLorentzVector v_jet;
			  v_jet.SetPtEtaPhiM(genpt[ijet], geneta[ijet], genphi[ijet], genm[ijet]);
			  dphi= GenQQ4mom->DeltaPhi(v_jet);
			  dr = GenQQ4mom->DeltaR (v_jet);
			  if (dr<=drmin)
			    {
			      drmin=dr;
			      dphimin=dphi;
			      deta=(jpsi_eta-geneta[ijet]);
			      z= jpsi_pt/genpt[ijet];
			    }
			}
		    }
		  if (drmin<0.5)
		    {
		      genmass->setVal(jpsi_m);
		      gen_pt->setVal(jpsi_pt);
		      gen_rap->setVal(jpsi_rap);  
		      genzed->setVal(z);
		      genr->setVal(drmin);
		      genctau->setVal(Gen_QQ_ctau[iQQ]);
		      gendata->add(*genset);
		      ismatch=true;
		    }
		}
	    }
	  if (!ismatch) continue;
	  if ( HLT_HIL1DoubleMu0ForPPRef_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) 
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();
		  drmin = 1000;

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  nt++;
			  tot->Fill(jpsi_pt, (1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)))));
			  for (Long64_t ijet=0; ijet<nref; ijet++)
			    {
			      if ((isPr && Reco_QQ_ctau[iQQ] < lcut->Eval(jpsi_pt))
			      || (!isPr && Reco_QQ_ctau[iQQ] > lcut->Eval(jpsi_pt)))
				{
				  if (abs(jty[ijet])<2.4 && jtpt[ijet]>50)
				    {
				      TLorentzVector v_jet;
				      v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
				      dphi= RecoQQ4mom->DeltaPhi(v_jet);
				      dr = RecoQQ4mom->DeltaR (v_jet);
				      if (dr<=drmin)
					{
					  drmin=dr;
					  dphimin=dphi;
					  deta=(jpsi_eta-jteta[ijet]);
					  z= jpsi_pt/jtpt[ijet];
					}
				    }
				}
			    }
			  if (drmin<0.5)
			    {
			      nj++;
			      jpjt->Fill(jpsi_pt, (1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)))));
			      mass->setVal(jpsi_m);
			      pt->setVal(jpsi_pt);
			      rap->setVal(jpsi_rap);
			      zed->setVal(z);
			      r->setVal(drmin);
			      ctau->setVal(Reco_QQ_ctau[iQQ]);
			      weight=(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			      if (weight==0)
				cout <<"rap= "<<jpsi_rap<<"  pt= "<<jpsi_pt<< "  e= "<<weight<<endl;
			      w->setVal(1.0/weight);
			      data->add(*recoset, w->getVal());
			      unwdata->add(*unwset);
			    }
			}
		    }
		}
	    }
	}
      cout<<"[INFO} N_{#mu#mu/Jets} = "<<nj<<endl<<"[INFO] N_{total} = "<<nt<<endl;
      for (int ibin=0; ibin <25; ibin++)
	{
	  tot->SetBinContent(ibin, ((tot->GetBinContent(ibin)*1.0)/tot->GetBinWidth(ibin)));
	  jpjt->SetBinContent(ibin, ((jpjt->GetBinContent(ibin)*1.0)/jpjt->GetBinWidth(ibin)));
	}
      gSystem->cd("ptdist");
      c1->cd();
      tot->SetMarkerColor(6);
      tot->SetMarkerStyle(33);
      jpjt->SetMarkerColor(kBlue);
      jpjt->SetMarkerStyle(33);

      TLegend* l1 = new TLegend (0.6, 0.6, 0.8, 0.7);
      l1->AddEntry(tot, "N_{#mu#mu tot}", "lep");
      l1->AddEntry(jpjt, "N_{#mu#mu in jets}", "lep");
      l1->SetBorderSize(0);

      TPaveText* tbox = new TPaveText(0.6,0.4,0.8,0.5, "BRNDC");
      tbox->AddText(Form("N_{#mu#mu tot} = %f", tot->GetEntries()));
      tbox->AddText(Form("N_{#mu#mu in jets} = %f", jpjt->GetEntries()));
      tbox->AddText(Form("N_{in jets}/N_{tot} = %f ", (nj*1.0)/nt));
      tbox->SetBorderSize(0);
      tbox->SetFillColor(0);

      tot->Draw();
      if (isPr)
	{
	  c1->SaveAs("prMcNtot.png");
	  c1->SaveAs("prMcNtot.C");
	}
      else
	{
	  c1->SaveAs("nprMcNtot.png");
	  c1->SaveAs("nprMcNtot.C");
	}
      jpjt->Draw();
      if (isPr)
	{
	  c1->SaveAs("prMcNjet.png");
	  c1->SaveAs("prMcNjet.C");
	}
      else
	{
	  c1->SaveAs("nprMcNjet.png");
	  c1->SaveAs("nprMcNjet.C");
	}
      tot->Draw();
      jpjt->Draw("same");
      l1->Draw("same");
      tbox->Draw("same");
      if (isPr)
	{
	  c1->SaveAs("prMc.png");
	  c1->SaveAs("prMc.C");
	}
      else
	{
	  c1->SaveAs("nprMc.png");
	  c1->SaveAs("nprMc.C");
	}
      gSystem->cd("..");
    }
  else
    {
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;

	  if ( HLT_HIL1DoubleMu0_v1 && pPAprimaryVertexFilter)
	    {
	      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
		{
		  TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
		  jpsi_pt = RecoQQ4mom->Pt();
		  jpsi_eta = RecoQQ4mom->Eta();
		  jpsi_rap = RecoQQ4mom->Rapidity();
		  jpsi_m = RecoQQ4mom->M();
		  drmin = 1000;

		  if (
		      jpsi_pt > 6.5  &&
		      (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		      (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		      (isTriggerMatch(iQQ, triggerIndex_PP)) // if it matches the trigger 
		      )
		    {
		      if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
			{
			  nt++;
			  tot->Fill(jpsi_pt, (1.0/eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt))));
			  if(Reco_QQ_ctau[iQQ] < lcut->Eval(jpsi_pt))
			    {
			      ntp++;
			      ptot->Fill(jpsi_pt, (1.0/preff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt))));
			    }
			  else
			    {
			      ntn++;
			      ntot->Fill(jpsi_pt, (1.0/npreff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt))));
			    }
			  for (Long64_t ijet=0; ijet<nref; ijet++)
			    {
	  
			      if (abs(jty[ijet])<2.4 && jtpt[ijet]>50)
				{
				  TLorentzVector v_jet;
				  v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
				  dphi= RecoQQ4mom->DeltaPhi(v_jet);
				  dr = RecoQQ4mom->DeltaR (v_jet);
				  if (dr<=drmin)
				    {
				      drmin=dr;
				      dphimin=dphi;
				      deta=(jpsi_eta-jteta[ijet]);
				      z= jpsi_pt/jtpt[ijet];
				    }
				}
			    }
			  if (drmin<0.5)
			    {
			      nj++;
			      jpjt->Fill(jpsi_pt, (1.0/eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt))));
			      mass->setVal(jpsi_m);
			      pt->setVal(jpsi_pt);
			      rap->setVal(jpsi_rap);
			      zed->setVal(z);
			      r->setVal(drmin);
			      weight=(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			      if (weight==0)
				{
				  cout <<"rap= "<<jpsi_rap<<"  pt= "<<jpsi_pt<< "  e= "<<weight<<endl;
				  weight=1;
				}
			      w->setVal(1.0/weight);
			      data->add(*recoset, w->getVal());
			      unwdata->add(*unwset);
			      if(Reco_QQ_ctau[iQQ] < lcut->Eval(jpsi_pt))
				{
				  np++;
				  pjpjt->Fill(jpsi_pt, (1.0/preff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt))));
				  weight=preff->GetEfficiency(preff->FindFixBin(jpsi_rap, jpsi_pt));
				  if (weight==0)
				    weight=1;
				  w->setVal(1.0/weight);
				  prdata->add(*prset, w->getVal());
				}
			      else
				{
				  nn++;
				  njpjt->Fill(jpsi_pt, (1.0/npreff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt))));
				  weight=npreff->GetEfficiency(npreff->FindFixBin(jpsi_rap, jpsi_pt));
				  if (weight==0)
				    weight=1;
				  w->setVal(1.0/weight);
				  nprdata->add(*nprset, w->getVal());
				}
			    }
			}
		    }
		}
	    }
	}
      cout<<"[INFO} N_{#mu#mu/Jets} = "<<nj<<"   N_{pr/Jets} = "<<np<<"   N_{npr/Jets} = "<<nn<<endl<<"[INFO] N_{total} = "<<nt<<"   N_{pr} = "<<ntp<<"   N_{npr} = "<<ntn<<endl;

      for (int ibin=0; ibin <25; ibin++)
	{
	  tot->SetBinContent(ibin, ((tot->GetBinContent(ibin)*1.0)/tot->GetBinWidth(ibin)));
	  jpjt->SetBinContent(ibin, ((jpjt->GetBinContent(ibin)*1.0)/jpjt->GetBinWidth(ibin)));
	  ptot->SetBinContent(ibin, ((ptot->GetBinContent(ibin)*1.0)/ptot->GetBinWidth(ibin)));
	  pjpjt->SetBinContent(ibin, ((pjpjt->GetBinContent(ibin)*1.0)/pjpjt->GetBinWidth(ibin)));
	  ntot->SetBinContent(ibin, ((ntot->GetBinContent(ibin)*1.0)/ntot->GetBinWidth(ibin)));
	  njpjt->SetBinContent(ibin, ((njpjt->GetBinContent(ibin)*1.0)/njpjt->GetBinWidth(ibin)));
	}

      gSystem->cd("ptdist");
      c1->cd();
      tot->SetMarkerColor(6);
      tot->SetMarkerStyle(33);
      jpjt->SetMarkerColor(kBlue);
      jpjt->SetMarkerStyle(33);

      ptot->SetMarkerColor(6);
      ptot->SetMarkerStyle(33);
      pjpjt->SetMarkerColor(kBlue);
      pjpjt->SetMarkerStyle(33);

      ntot->SetMarkerColor(6);
      ntot->SetMarkerStyle(33);
      njpjt->SetMarkerColor(kBlue);
      njpjt->SetMarkerStyle(33);

      TLegend* l1 = new TLegend (0.6, 0.6, 0.8, 0.7);
      l1->AddEntry(tot, "N_{#mu#mu tot}", "lep");
      l1->AddEntry(jpjt, "N_{#mu#mu in jets}", "lep");
      l1->SetBorderSize(0);

      TLegend* l2 = new TLegend (0.6, 0.6, 0.8, 0.7);
      l2->AddEntry(ptot, "N_{pr#mu#mu tot}", "lep");
      l2->AddEntry(pjpjt, "N_{pr#mu#mu in jets}", "lep");
      l2->SetBorderSize(0);

      TLegend* l3 = new TLegend (0.6, 0.6, 0.8, 0.7);
      l3->AddEntry(ntot, "N_{npr#mu#mu tot}", "lep");
      l3->AddEntry(njpjt, "N_{npr#mu#mu in jets}", "lep");
      l3->SetBorderSize(0);

      TPaveText* tbox = new TPaveText(0.6,0.4,0.8,0.5, "BRNDC");
      tbox->AddText(Form("N_{#mu#mu tot} = %f", tot->GetEntries()));
      tbox->AddText(Form("N_{#mu#mu in jets} = %f", jpjt->GetEntries()));
      tbox->AddText(Form("N_{in jets}/N_{tot} = %f ", jpjt->GetEntries()*1.0/tot->GetEntries()));
      tbox->SetBorderSize(0);
      tbox->SetFillColor(0);

      TPaveText* tbox1 = new TPaveText(0.6,0.4,0.8,0.5, "BRNDC");
      tbox1->AddText(Form("N_{pr#mu#mu tot} = %f", ptot->GetEntries()));
      tbox1->AddText(Form("N_{pr#mu#mu in jets} = %f", pjpjt->GetEntries()));
      tbox1->AddText(Form("N_{in jets}/N_{tot} = %f ", pjpjt->GetEntries()*1.0/ptot->GetEntries()));
      tbox1->SetBorderSize(0);
      tbox1->SetFillColor(0);

      TPaveText* tbox2 = new TPaveText(0.6,0.4,0.8,0.5, "BRNDC");
      tbox2->AddText(Form("N_{npr#mu#mu tot} = %f", ntot->GetEntries()));
      tbox2->AddText(Form("N_{npr#mu#mu in jets} = %f", njpjt->GetEntries()));
      tbox2->AddText(Form("N_{in jets}/N_{tot} = %f ", (njpjt->GetEntries()*1.0)/ntot->GetEntries()));
      tbox2->SetBorderSize(0);
      tbox2->SetFillColor(0);

      tot->Draw();
      c1->SaveAs("dataTot.png");
      c1->SaveAs("dataTot.C");
      jpjt->Draw();
      c1->SaveAs("dataJet.png");
      c1->SaveAs("dataJet.C");
      tot->Draw();
      jpjt->Draw("same");
      l1->Draw("same");
      tbox->Draw("same");
      c1->SaveAs("data.png");
      c1->SaveAs("data.C");

      ptot->Draw();
      c1->SaveAs("prdataTot.png");
      c1->SaveAs("prdataTot.C");
      pjpjt->Draw();
      c1->SaveAs("prdataJet.png");
      c1->SaveAs("prdataJet.C");
      ptot->Draw();
      pjpjt->Draw("same");
      l2->Draw("same");
      tbox1->Draw("same");
      c1->SaveAs("prdata.png");
      c1->SaveAs("prdata.C");

      ntot->Draw();
      c1->SaveAs("nprdataTot.png");
      c1->SaveAs("nprdataTot.C");
      njpjt->Draw();
      c1->SaveAs("nprdataJet.png");
      c1->SaveAs("nprdataJet.C");
      ntot->Draw();
      njpjt->Draw("same");
      l3->Draw("same");
      tbox2->Draw("same");
      c1->SaveAs("nprdata.png");
      c1->SaveAs("nprdata.C");

      gSystem->cd("..");
    }	    
  TFile *datafile (0x0);
  if (isMc)
    {
      if (isPr)
	datafile = new TFile ("prMcDatasets","RECREATE");
      else 
	datafile = new TFile ("nprMcDatasets","RECREATE");
      datafile->cd();
      gendata->Write("genData");
      data->Write("recoData");
      unwdata->Write("unwData");
    }
  else
    {
      if (isPr)
	datafile = new TFile ("datasets_pr","RECREATE");
      else
	datafile = new TFile ("datasets_npr","RECREATE");
      datafile->cd();
      data->Write("data");
      unwdata->Write("unwData");
      prdata->Write("prData");
      nprdata->Write("nprData");
    }
      datafile->Write();
      datafile->Close();
      delete datafile;
}

void myTree::Plot()
{
  Float_t nv[12];
  Float_t av[12];
  Double_t zbins [] = {0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1};
  TFile* mcfile (0x0);
  TFile* datafile (0x0);
  TH1F* gn = new TH1F ("gn", "values of n from mc", 12, 0, 1.1);
  TH1F* ga = new TH1F ("ga", "values of alpha from mc", 12, 0, 1.1);
  TH1F* nz = new TH1F ("nz", "N_{J/#psi} in z bins", 7, zbins);
  //TH1F* nz = new TH1F ("nz", "N_{J/#psi} in z bins", 10, 0, 1);
  TH1F* np = new TH1F ("np", "N_{J/#psi} in z bins", 7, zbins);
  TAxis* xaxis(0x0);
  if (isPr)
    {
      mcfile = TFile::Open("prMcDatasets");
      datafile = TFile::Open("datasets_pr");
    }
  else
    {
      mcfile = TFile::Open("nprMcDatasets");
      datafile = TFile::Open("datasets_npr");
    }
  gendata= (RooDataSet*) mcfile->Get("genData");
  prdata= (RooDataSet*) datafile->Get("prData");
  nprdata= (RooDataSet*) datafile->Get("nprData");
  if (isMc)
    {
      data= (RooDataSet*) mcfile->Get("recoData");
      unwdata= (RooDataSet*) mcfile->Get("unwData");
    }
  else
    {  
      data= (RooDataSet*) datafile->Get("data");
      unwdata= (RooDataSet*) datafile->Get("unwData");
    }
  genzed->setBins(10);
  mass->setBins(40);
  zed->setBins(10);
  rap->setBins(12);
  gen_rap->setBins(12);
  pt->setBins(19);
  gen_pt->setBins(19);
  r->setMax(2);	
  //r->setBins(10);
  RooDataSet* d = (RooDataSet*) data->reduce("Mass>2.6 && Mass<3.5 && r<0.5");
  RooDataSet* unwd = (RooDataSet*) unwdata->reduce("Mass>2.6 && Mass<3.5 && r<0.5");
  RooDataSet* gend = (RooDataSet*) gendata->reduce("genMass>2.6 && genMass<3.5 && genr<0.5");

  RooPlot* mframe = mass->frame();
  RooPlot* comp = mass->frame();
  RooPlot* zframe = zed->frame();
  RooPlot* rframe = r->frame();
  RooPlot* ptframe = pt->frame();
  RooPlot* rapframe = rap->frame();

  RooRealVar* mean = new RooRealVar ("mean", "mean", 3.096, 3, 3.2);
  RooRealVar* sigma = new RooRealVar ("sigma", "sigma", 0.04, 0, 0.1);
  RooRealVar* sigma1 = new RooRealVar("sigma1", "sigma1", 0.02, 0, 0.1);
  RooRealVar* alpha (0x0);
  RooRealVar* n (0x0);


  if (isMc) 
    {
      alpha = new RooRealVar ("alpha", "alpha", 1.5, 0, 4);
      n = new RooRealVar("n", "n", 1.5, 0.2, 4);
    }
  else
    {
      alpha = new RooRealVar ("alpha", "alpha", 1.94);
      n = new RooRealVar("n", "n", 1.39);
    }

  RooCBShape cb ("cb", "cb", *mass, *mean, *sigma, *alpha, *n);
  RooCBShape cb1 ("cb1", "cb1", *mass, *mean, *sigma1, *alpha, *n);

  RooRealVar c0("c0","coefficient #0", -0.1, -5, 5);
  RooRealVar c1("c1","coefficient #1", 0.1, -5, 5);
  //RooRealVar c2("c2","coefficient #1", 0.1, -5, 5);
  RooChebychev bkg("bkg","background p.d.f.", *mass, RooArgList(c0, c1)) ;

  RooRealVar fb ("fb", "background nb", 1000, 0, data->numEntries());
  RooRealVar fs ("fs", "signal nb", 10000, 0, data->numEntries());
  RooRealVar fcb ("fcb", "cb", 0.1, 0, 1);

  mass->setRange("window", 2.9, 3.2) ;

  RooExtendPdf ecb ("ecb", "ecb", cb, fcb);

  RooAddPdf masspeak ("masspeak", "fit of mass peak", RooArgList(cb, cb1), fcb);
  RooExtendPdf es ("es", "es", masspeak, fs, "window");
  RooExtendPdf eb ("eb", "eb", bkg, fb);
  RooAddPdf* model;
  if (isMc)
    model = new RooAddPdf ("model", "model", RooArgList(cb, cb1), fcb);
  else 
    model = new RooAddPdf ("model", "model", RooArgList(es, eb));

  d->plotOn(mframe,DataError(RooAbsData::SumW2));
  d->plotOn(rframe);
  //rframe->Draw();
  d->plotOn(ptframe, DataError(RooAbsData::SumW2));
  d->plotOn(rapframe, DataError(RooAbsData::SumW2));

  model->fitTo(*d, Extended(!isMc), SumW2Error(kTRUE));
  if (isMc)
    {
      xaxis=gn->GetXaxis();
      gn->SetBinContent(xaxis->FindBin(1), n->getValV());
      xaxis=ga->GetXaxis();
      ga->SetBinContent(xaxis->FindBin(1), alpha->getValV());
    }
  model->plotOn(mframe, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
  model->plotOn(mframe, Name("signal"), Components(RooArgSet(masspeak)), LineStyle(kDashed));
  model->plotOn(mframe, Name("total"), LineColor(kRed));
  //model->paramOn(mframe);


  TFile* fsave (0x0);
  if (isMc)
    {
      if (isPr)
	fsave = new TFile ("prMcplots.root","RECREATE");
      else
	fsave = new TFile ("nprMcplots.root", "RECREATE");
    }
  else
    {
      if (isPr)
	fsave = new TFile("prdataplots.root", "RECREATE");
      else 
	fsave = new TFile("nprdataplots.root", "RECREATE");
    }

  RooWorkspace* w (0x0);
  if (isMc)
    {
      for (int i =2; i<9; i++)
	{
	  RooPlot* mzframe = mass->frame();
	  w = new RooWorkspace(Form("w%d",i));
	  RooDataSet* dz(0x0);
	  if (i==2)
	  dz = (RooDataSet*) d->reduce("z>0 && z<0.3");
	  else if (i==8)
	    dz = (RooDataSet*) d->reduce("z>0.8 && z<=1");
	  else 
	    dz = (RooDataSet*) d->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
	  w->import(*dz);
	  w->import(*model);
	  w->data("data")->plotOn(mzframe, DataError(RooAbsData::SumW2));
	  w->pdf("model")->fitTo(*w->data("data"),  Extended(!isMc), SumW2Error(kTRUE));
	  w->pdf("model")->plotOn(mzframe);
	  if (isMc)
	    {
	      xaxis=gn->GetXaxis();
	      gn->SetBinContent(xaxis->FindBin(i*0.1+0.05), w->var("n")->getValV());
	      xaxis=ga->GetXaxis();
	      ga->SetBinContent(xaxis->FindBin(i*0.1+0.05), w->var("alpha")->getValV());
	      xaxis=nz->GetXaxis();
	      nz->SetBinContent(xaxis->FindBin(i*0.1+0.05), w->data("data")->sumEntries());
	      nz->SetBinError(xaxis->FindBin(i*0.1+0.05), sqrt(w->data("data")->sumEntries())); 
	    }
	  TPaveText* tbox = new TPaveText(0.15,0.6,0.4,0.8, "BRNDC");
	  tbox->AddText("Monte Carlo");
	  if (isPr)
	    tbox->AddText("prompt J/#psi");
	  else
	    tbox->AddText("b-decay J/#psi");
	  if (i==2)
	  tbox-> AddText("0 #leq z < 0.3");
	  else if (i==8)
	  tbox->AddText("0.8 #leq z #leq 1");
	  else 
	    tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	  tbox->AddText("p_{t}(J/#psi)>6.5GeV");
	  tbox->AddText("p_{t}(jets)>50GeV");
	  tbox->AddText("|y|<2.4");
	  mzframe->addObject(tbox);
	  tbox->SetBorderSize(0);
	  tbox->SetFillColor(0);
	  //w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("massframe%d%d",i,(i+1)));
	}
    }

  else
    {
      for (int i =2; i<9; i++)
	{
	  w = new RooWorkspace(Form("prw%d",i));
	  RooPlot* mzframe = mass->frame();
	  RooDataSet* prdz (0x0);
	  if (i==2)
	    prdz=(RooDataSet*) prdata->reduce("z>=0 && z<0.3");
	  else if (i==8)
	    prdz=(RooDataSet*) prdata->reduce("z>=0.8 && z<=1");
	  else 
	    prdz= (RooDataSet*) prdata->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
	  w->import(*prdz);
	  w->import(*model);
	      if (i==3)
		{
		  w->var("c0")->setVal(0.5);
		  w->var("c0")->setConstant(false);
		  w->var("c0")->setRange(-3, 3);
		  w->var("c1")->setVal(0.2);
		  w->var("c1")->setConstant(false);
		  w->var("c1")->setRange(-3, 3);
		  w->var("fs")->setVal(15000);
		  w->var("fs")->setConstant(false);
		  w->var("fb")->setVal(1000);
		  w->var("fb")->setConstant(false);
		  w->var("fcb")->setVal(0.8);
		  w->var("fcb")->setConstant(false);
		  w->var("fcb")->setMin(0.3);
		  //w->var("sigma")->setRange(0.0, 0.05);
		  //w->var("sigma1")->setRange(0.01, 0.2);
		  //w->var("c0")->setRange(-1, 1);
		  //w->var("c1")->setRange(-1, 3);
		  //w->var("fs")->setMax(50000);
		  //w->var("fcb")->setRange(0.05, 0.6);
		}

	      //if  (i ==6)
	      //{
	      //}
	  w->data("prdata")->plotOn(mzframe, DataError(RooAbsData::SumW2));
	  w->pdf("model")->fitTo(*w->data("prdata"),  Extended(!isMc), SumW2Error(kTRUE));
	  //w->pdf("model")->plotOn(mzframe);

	  w->pdf("model")->plotOn(mzframe, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
	  w->pdf("model")->plotOn(mzframe, Name("signal"), Components(RooArgSet(masspeak)), LineStyle(kDashed));
	  w->pdf("model")->plotOn(mzframe, Name("total"), LineColor(kRed));
	  w->data("prdata")->plotOn(mzframe, DataError(RooAbsData::SumW2));

	  xaxis=nz->GetXaxis();
	  nz->SetBinContent(xaxis->FindBin(i*0.1+0.05), w->var("fs")->getValV());
	  nz->SetBinError(xaxis->FindBin(i*0.1+0.05), w->var("fs")->getError());

	  TPaveText* tbox = new TPaveText(0.15,0.6,0.4,0.8, "BRNDC");
	  tbox->AddText("prompt J/#psi");
	  if(i==2)
	    tbox->AddText("0 #leq z < 0.3");
	  else if (i==8)
	    tbox->AddText("0.8 #leq z #leq 1");
	  else 
	    tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	  tbox->AddText("p_{t}(J/#psi)>6.5GeV");
	  tbox->AddText("p_{t}(jets)>50GeV");
	  tbox->AddText("|y|<2.4");
	  mzframe->addObject(tbox);
	  tbox->SetBorderSize(0);
	  tbox->SetFillColor(0);
	  TLegend* l = new TLegend (0.7,0.6,1,0.8);
	  l->AddEntry(mzframe->RooPlot::findObject("signal"), "J/#psi", "L");
	  l->AddEntry(mzframe->RooPlot::findObject("background"), "background", "f");
	  l->AddEntry(mzframe->RooPlot::findObject("total"), "Total fit", "L");
	  l->SetBorderSize(0);
	  mzframe->addObject(l);
	  w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("prmfr%d%d",i,(i+1)));
	}

      for (int i =2; i<9; i++)
	{
	  w = new RooWorkspace(Form("nprw%d",i));
	  RooPlot* mzframe = mass->frame();
	  RooDataSet* nprdz (0x0);
	  if (i==2)
	    nprdz=(RooDataSet*) nprdata->reduce("z>=0 && z<0.3");
	  else if (i==8)
	    nprdz=(RooDataSet*) nprdata->reduce("z>=0.8 && z<=1");
	  else 
	    nprdz=(RooDataSet*) nprdata->reduce(Form("z >= %f && z<%f", i*0.1, (i+1)*0.1));
	  w->import(*nprdz);
	      w->import(*model);
	      if (i==2)
		{
		  //w->var("sigma")->setVal(0.08);
		  //w->var("sigma1")->setVal(0.1);
		  w->var("sigma")->setRange(0.0, 0.05);
		  w->var("sigma1")->setRange(0.02, 0.1);
		  //w->var("c0")->setVal(-0.003);
		  w->var("c0")->setRange(-1, 1);
		  //w->var("c1")->setVal(1);
		  w->var("c1")->setRange(-1, 3);
		  //w->var("fs")->setVal(20000);
		  w->var("fs")->setMax(100000);
		  //w->var("fb")->setVal(1500);
		  //w->var("fcb")->setVal(0.4);
		  w->var("fcb")->setMin(0.3);
		 }
	      if (i==8)
		{
		  w->var("sigma")->setRange(0.0, 0.05);
		  w->var("sigma1")->setRange(0.02, 0.1);
		  w->var("c0")->setRange(-1, 1);
		  w->var("c1")->setRange(-1, 3);
		  w->var("fs")->setMax(10000);
		  w->var("fcb")->setMin(0.3);
		}
	      w->data("nprdata")->plotOn(mzframe, DataError(RooAbsData::SumW2));
	      w->pdf("model")->fitTo(*w->data("nprdata"),  Extended(!isMc), SumW2Error(kTRUE));
	      //w->pdf("model")->plotOn(mzframe);
	      w->pdf("model")->plotOn(mzframe, Name("background"), Components(RooArgSet(bkg)),DrawOption("F"), FillColor(kGray));
	      w->pdf("model")->plotOn(mzframe, Name("signal"), Components(RooArgSet(masspeak)), LineStyle(kDashed));
	      w->pdf("model")->plotOn(mzframe, Name("total"), LineColor(kRed));
	      w->data("nprdata")->plotOn(mzframe, DataError(RooAbsData::SumW2));

	      xaxis=np->GetXaxis();
	      np->SetBinContent(xaxis->FindBin(i*0.1+0.05), w->var("fs")->getValV());
	      np->SetBinError(xaxis->FindBin(i*0.1+0.05), w->var("fs")->getError());
	  TPaveText* tbox = new TPaveText(0.15,0.6,0.4,0.8, "BRNDC");
	  tbox->AddText("b-decay J/#psi");
	  if(i==2)
	    tbox->AddText("0 #leq z < 0.3");
	  else if (i==8)
	    tbox->AddText("0.8 #leq z #leq 1");
	  else 
	    tbox-> AddText (Form("%.1f #leq z < %.1f", i*0.1, (i+1)*0.1));
	  tbox->AddText("p_{t}(J/#psi)>6.5GeV");
	  tbox->AddText("p_{t}(jets)>50GeV");
	  tbox->AddText("|y|<2.4");
	  mzframe->addObject(tbox);
	  tbox->SetBorderSize(0);
	  tbox->SetFillColor(0);
	  TLegend* l = new TLegend (0.7,0.6,1,0.8);
	  l->AddEntry(mzframe->RooPlot::findObject("signal"), "J/#psi", "L");
	  l->AddEntry(mzframe->RooPlot::findObject("background"), "background", "f");
	  l->AddEntry(mzframe->RooPlot::findObject("total"), "Total fit", "L");
	  l->SetBorderSize(0);
	  mzframe->addObject(l);
	  w->pdf("model")->paramOn(mzframe);
	  mzframe->Write(Form("nprmfr%d%d",i,(i+1)));
	}
    }



  mframe->Write("massframe");
  rframe->Write("rframe");
  ptframe->Write("recopt");
  rapframe->Write("recorap");
 if (isMc)
   {
     RooPlot* genzframe= genzed->frame();
     gend->plotOn(genzframe, MarkerColor(kBlue), LineColor(kOrange+2), LineStyle(1));
     //TLegend* l = new TLegend (0.2,0.6,0.4,0.8);
     //l->AddEntry(zframe->RooPlot::findObject("recod"), "N_{reco J/#psi} with corrections", "L");
     //l->AddEntry(genzframe->RooPlot::findObject("gend"), "N_{reco J/#psi} with corrections", "L");
     //l->SetBorderSize(0);
     //genzframe->addObject(l);
     //genzframe->Write("genz");
     gn->Write("nv");
     ga->Write("av");
   }
 else
   {
     RooPlot* prmass=mass->frame();
     RooPlot* nprmass=mass->frame();
     prdata->plotOn(prmass, DataError(RooAbsData::SumW2));
     nprdata->plotOn(nprmass, DataError(RooAbsData::SumW2));
     prmass->Write("prmassframe");
     nprmass->Write("nprmassframe");
   }

 nz->Scale(1.0/(nz->Integral()));
 //RooDataHist* zd = new RooDataHist ("zd", "N_{J/#psi} in z bins", *zed, nz);
 if (isMc)
   {
     //d->plotOn(zframe, Name("zmc"), MarkerStyle(kFullSquare), MarkerColor(kBlue));
     //zframe->Draw();
     //zframe->Write("zframe");
     //zd= new RooDataHist ("nzd", "N_{J/#psi} in z bins", *zed, nz);
     //zd->plotOn(zframe, Name("ze"), MarkerColor(kBlue));
     //zframe->Write("nz");
     nz->Write("nz");
   }
 else
   {
     //RooPlot* nprzfr=zed->frame();
     np->Scale(1.0/(np->Integral()));
     //zd->plotOn(zframe, Name("prz"), MarkerColor(kRed));
     //zd= new RooDataHist ("nzd", "N_{J/#psi} in z bins", *zed, np);
     //zframe->Write("prz");
     //zd->plotOn(nprzfr, Name("nprz"), MarkerColor(kRed));
     //nprzfr->Write("nprz");
     nz->Write("prz");
     np->Write("nprz");
   }
 fsave->Close();
}



Bool_t myTree::isTriggerMatch (Int_t iRecoQQ, Int_t TriggerBit) 
  {
    Bool_t cond = true;
    cond = cond && ( (HLTriggers&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) ); 
    cond = cond && ( (Reco_QQ_trig[iRecoQQ]&((ULong64_t)pow(2, TriggerBit))) == ((ULong64_t)pow(2, TriggerBit)) );
    return cond;
  };

Bool_t myTree::isGlobalMuonInAccept2015 (TLorentzVector* Muon) 
  {
  return (fabs(Muon->Eta()) < 2.4 &&
          ((fabs(Muon->Eta()) < 1.2 && Muon->Pt() >= 3.5) ||
           (1.2 <= fabs(Muon->Eta()) && fabs(Muon->Eta()) < 2.1 && Muon->Pt() >= 5.77-1.89*fabs(Muon->Eta())) ||
           (2.1 <= fabs(Muon->Eta()) && Muon->Pt() >= 1.8)));
  };

Bool_t myTree::areMuonsInAcceptance2015 (Int_t iRecoQQ)
  {
    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoQQ);
    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoQQ);
    return ( isGlobalMuonInAccept2015(RecoQQmupl) && isGlobalMuonInAccept2015(RecoQQmumi) );
  };

Bool_t myTree::areGenMuonsInAcceptance2015 (Int_t iGenQQ)
  {
    TLorentzVector *GenQQmupl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenQQ);
    TLorentzVector *GenQQmumi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenQQ);
    return (isGlobalMuonInAccept2015(GenQQmupl) && isGlobalMuonInAccept2015(GenQQmumi));
  };  
  
Bool_t myTree::passQualityCuts2015 (Int_t iRecoQQ) 
  {
    Bool_t cond = true;
    cond = cond && (Reco_QQ_mumi_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 1)));
    cond = cond && (Reco_QQ_mumi_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 3)));    
    // cond = cond && (Reco_QQ_mumi_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mumi_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mumi_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mumi_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mumi_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mumi_dz[iRecoQQ]) < 20.);
    
    cond = cond && (Reco_QQ_mupl_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 1)));
    cond = cond && (Reco_QQ_mupl_SelectionType[iRecoQQ]&((ULong64_t)pow(2, 3)));
    // cond = cond && (Reco_QQ_mupl_highPurity[iRecoQQ]);
    cond = cond && (Reco_QQ_mupl_isGoodMuon[iRecoQQ]==1);
    cond = cond && (Reco_QQ_mupl_nTrkWMea[iRecoQQ] > 5);
    cond = cond && (Reco_QQ_mupl_nPixWMea[iRecoQQ] > 0);
    cond = cond && (fabs(Reco_QQ_mupl_dxy[iRecoQQ]) < 0.3);
    cond = cond && (fabs(Reco_QQ_mupl_dz[iRecoQQ]) < 20.);    
    cond = cond && (Reco_QQ_VtxProb[iRecoQQ] > 0.01);
    
    return cond;
  };


Double_t myTree::deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

Bool_t myTree::isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR)
{
  TLorentzVector* RecoMuonpl = (TLorentzVector*) Reco_QQ_mupl_4mom->At(iRecoDiMuon);
  TLorentzVector* RecoMuonmi = (TLorentzVector*) Reco_QQ_mumi_4mom->At(iRecoDiMuon);
  
  bool isMatched(false);
  int iGenMuon(0);
  while ( !isMatched && (iGenMuon < Gen_QQ_size) )
  {
    TLorentzVector *GenMuonpl = (TLorentzVector*)Gen_QQ_mupl_4mom->At(iGenMuon);
    TLorentzVector *GenMuonmi = (TLorentzVector*)Gen_QQ_mumi_4mom->At(iGenMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  ) 
      {
	isMatched = true;
	matchGen = (TLorentzVector*) Gen_QQ_4mom->At(iGenMuon);
      }
    iGenMuon++;
  }
  return isMatched;
};



Bool_t myTree::isMatchedGenDiMuon(int iGenDiMuon, double maxDeltaR)
{
  TLorentzVector* GenMuonpl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuonmi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iGenDiMuon);
  TLorentzVector* GenMuon = (TLorentzVector*) Gen_QQ_4mom->At(iGenDiMuon);

  bool isMatched(false);
  int iRecoMuon(0);
  while ( !isMatched && (iRecoMuon < Reco_QQ_size) )
  {
    TLorentzVector *RecoMuonpl = (TLorentzVector*)Reco_QQ_mupl_4mom->At(iRecoMuon);
    TLorentzVector *RecoMuonmi = (TLorentzVector*)Reco_QQ_mumi_4mom->At(iRecoMuon);
    double dRpl = deltaR(GenMuonpl,RecoMuonpl);
    double dRmi = deltaR(GenMuonmi,RecoMuonmi);
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR) ) 
      {
	isMatched = true;
	matchReco = (TLorentzVector*) Reco_QQ_4mom->At(iRecoMuon);
      }
    iRecoMuon++;
  }
  
  return isMatched;
};
void myTree::Unfolding()
{
  TH2D* mcpt = new TH2D ("mcpt","pt distribution in MC; gen pt; reco pt", 44, 6.5, 50.5, 44, 6.5, 50.5);
  TH1D* datapt = new TH1D ("datapt","pt distribution in data ; pt",44, 6.5, 50.5);
  TH1D* recopt = new TH1D ("recopt","reco pt distribution ; pt",44, 6.5, 50.5);
  TH1D* genpt = new TH1D ("genpt","gen pt distribution; pt",44, 6.5, 50.5);
  Float_t jet_pt;
  TFile* f (0x0);
  if (isPr)
    f= TFile::Open("prEff.root");
  else 
    f= TFile::Open("nprEff.root");

  TEfficiency* eff = (TEfficiency*) f->Get("ptrap");
  f->Close();
  if (isMc)
    {
      TFile* fsave(0x0);
      if (isPr)
	fsave = new TFile ("prMcunfold1.root","RECREATE");
      else
	fsave = new TFile ("nprMcunfold1.root", "Recreate");

      TString filename;
      if (isPr)
	filename ="prTruthMeasured";
      else
	filename ="nprTruthMeasured";

      ofstream file_out(filename);
      file_out<< "gen      reco"<<endl;
      Long64_t nentries = fChain->GetEntries();
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<(nentries/2);jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for(int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
	    {
	      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	      jpsi_pt = RecoQQ4mom->Pt();
	      jpsi_eta = RecoQQ4mom->Eta();
	      jpsi_rap = RecoQQ4mom->Rapidity();
	      jpsi_m = RecoQQ4mom->M();
	      drmin = 1000;

	      if (
		  jpsi_pt > 6.5  &&
		  (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		  (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		  (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		  )
		{
		  if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
		    {
		      for (Long64_t ijet=0; ijet<nref; ijet++)
			{
	  
			  if (abs(jty[ijet])<2.4 && jtpt[ijet]>15)
			    {
			      TLorentzVector v_jet;
			      v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
			      dphi= RecoQQ4mom->DeltaPhi(v_jet);
			      dr = RecoQQ4mom->DeltaR (v_jet);
			      if (dr<=drmin)
				{
				  drmin=dr;
				  //z= jpsi_pt/jtpt[ijet];
				  //jet_pt=jtpt[ijet];
				}
			    }
			}
		      if (drmin<0.5 && isMatchedRecoDiMuon(iQQ))
			{
			  weight=1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			  mcpt->Fill(matchGen->Pt(), jpsi_pt, weight);
			  recopt->Fill(jpsi_pt, weight);
			  genpt->Fill(matchGen->Pt(), weight);
			  file_out << matchGen->Pt()<<"  "<<jpsi_pt<< endl;
			}
		    }
		}
	    }
	}
      file_out.close();
      mcpt->Write("ptdist2D");
      recopt->Write("recopt1D");
      genpt->Write("genpt1D");
      fsave->Close();
    }
  else
    {
      Long64_t nentries =fChain->GetEntries();
      Long64_t nbytes = 0, nb = 0;
      for (Long64_t jentry=0; jentry<nentries;jentry++) 
	{
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;


	  for(int iQQ=0; iQQ<Reco_QQ_size; iQQ++)
	    {
	      TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
	      jpsi_pt = RecoQQ4mom->Pt();
	      jpsi_eta = RecoQQ4mom->Eta();
	      jpsi_rap = RecoQQ4mom->Rapidity();
	      jpsi_m = RecoQQ4mom->M();
	      drmin = 1000;

	      if (
		  jpsi_pt > 6.5  &&
		  (areMuonsInAcceptance2015(iQQ))&&  // 2015 Global Muon Acceptance Cuts
		  (passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
		  (isTriggerMatch(iQQ, triggerIndex_PP)) //&&// if it matches the trigger 
		  )
		{
		  if (Reco_QQ_sign[iQQ]==0 && abs(jpsi_rap)<2.4 && jpsi_m>2.6 && jpsi_m<3.5) 
		    {
		      for (Long64_t ijet=0; ijet<nref; ijet++)
			{
	  
			  if (abs(jty[ijet])<2.4 && jtpt[ijet]>15)
			    {
			      TLorentzVector v_jet;
			      v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
			      dphi= RecoQQ4mom->DeltaPhi(v_jet);
			      dr = RecoQQ4mom->DeltaR (v_jet);
			      if (dr<=drmin)
				{
				  drmin=dr;
				  //z= jpsi_pt/jtpt[ijet];
				  //jet_pt=jtpt[ijet];
				}
			    }
			}
		      if (drmin<0.5)
			{
			  weight=1.0/(eff->GetEfficiency(eff->FindFixBin(jpsi_rap, jpsi_pt)));
			  datapt->Fill(jpsi_pt, weight);
			}
		    }
		}
	    }
	}
      TFile* fsave(0x0);
      if (isPr)
	fsave = new TFile ("prdataunfold.root","RECREATE");
      else
	fsave = new TFile ("nprdataunfold.root", "RECREATE");
      datapt->Write("ptdist");
      fsave->Close();
    }
}

void myTree::JetPtRange()
{ 
  gSystem->mkdir("JtPtRange");
  gSystem->cd("JtPtRange");
  int imatch=0;
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "", 1000,800);
  TH1F* genptdist = new TH1F ("genptdist", ";p_{T}(jet) (GeV/c)", 100, 10, 50);
  TH1F* recoptdist = new TH1F ("recoptdist", ";p_{T}(jet) (GeV/c)", 100, 10, 50);
  TH1F* ptres = new TH1F ("ptres", ";p_{T}^{gen}-p_{T}^{reco}(jet) (GeV/c)", 40, -20, 20);
  Long64_t nentries = 1000000; //fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;  

      for (Long64_t ijet=0; ijet<nref; ijet++)
	{
	  
	  if (abs(jty[ijet])<1.6 /*&& abs(jty[ijet])<2.1*/ && jtpt[ijet]>25 && jtpt[ijet]<40)
	    {
	      recoptdist->Fill(jtpt[imatch]);
	      genptdist->Fill(refpt[imatch]);
	      ptres->Fill(refpt[imatch]-jtpt[imatch]); 
	    }
	}
    }

  TString filename;
  if (isPr)
    filename ="prJtPt2540_mid.dat";
  else
    filename ="nprJtPt2540_mid.dat";

  ofstream file_out(filename);

  for (int i=0; i<genptdist->GetSize(); i++)
    {
      file_out<< "for min pt(jet) = "<<genptdist->GetXaxis()->GetBinLowEdge(i)<< " Ngen/Nreco = "<<genptdist->Integral(i, genptdist->GetSize()-1)*1.0/recoptdist->Integral()<<endl;
      cout<< "for min pt(jet) = "<<genptdist->GetXaxis()->GetBinLowEdge(i)<< " Ngen/Nreco = "<<genptdist->Integral(i, genptdist->GetSize()-1)*1.0/recoptdist->Integral()<<endl;
    }
  file_out.close();
  c->cd();
  ptres->Draw();
  if (isPr)
    c->SaveAs("prJtPtRes2540_mid.pdf");
  else
    c->SaveAs("nprJtPtRes2540_mid.pdf");
  recoptdist->SetLineColor(4);
  recoptdist->Draw();
  genptdist->SetLineColor(2);
  genptdist->Draw("same");
  TLegend* l = new TLegend (0.65,0.7,0.85,0.8);
  l->AddEntry(genptdist, "gen", "lp");
  l->AddEntry(recoptdist, "reco", "lp");
  l->SetBorderSize(0);
  l->Draw("same");

  TPaveText* tbox0 = new TPaveText(0.6,0.5,0.8,0.7, "BRNDC");
  if (isPr)
    tbox0->AddText("prompt MC");
  else
    tbox0->AddText("nonprompt MC");
  tbox0->AddText("25 < p_{T}^{reco}(jet) < 40 GeV/c");
  tbox0->AddText("|y(jet)| < 1.6");
  tbox0->SetBorderSize(0);
  tbox0->SetFillColor(0);
  tbox0->Draw("same");
  if (isPr)
    c->SaveAs("prJtPtDist2540_mid.pdf");
  else
    c->SaveAs("nprJtPtDist2540_mid.pdf");
  gSystem->cd("..");
}
