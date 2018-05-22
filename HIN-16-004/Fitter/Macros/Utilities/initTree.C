//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jan 16 22:46:58 2016 by ROOT version 5.34/32
// from TTree myTree/My TTree of dimuons
// found on file: OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262620_263511_noCUT.root
//////////////////////////////////////////////////////////

#ifndef initTree_C
#define initTree_C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

Int_t           fCurrent; //!current Tree number in a TChain

// Declaration of leaf types
UInt_t          eventNb;
UInt_t          runNb;
UInt_t          LS;
Float_t         zVtx;
Int_t           Centrality;
Int_t           nTrig;
Int_t           trigPrescale[99];   //[nTrig]
ULong64_t       HLTriggers;
Int_t           Npix;
Int_t           NpixelTracks;
Int_t           Ntracks;
Float_t         SumET_HF;
Float_t         SumET_HFplus;
Float_t         SumET_HFminus;
Float_t         SumET_HFplusEta4;
Float_t         SumET_HFminusEta4;
Float_t         SumET_ET;
Float_t         SumET_EE;
Float_t         SumET_EB;
Float_t         SumET_EEplus;
Float_t         SumET_EEminus;
Float_t         SumET_ZDC;
Float_t         SumET_ZDCplus;
Float_t         SumET_ZDCminus;
Int_t           nEP;
Float_t         rpAng[99];   //[nEP]
Float_t         rpSin[99];   //[nEP]
Float_t         rpCos[99];   //[nEP]
Int_t           Reco_QQ_size;
Int_t           Reco_QQ_type[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_sign[99];   //[Reco_QQ_size]
TClonesArray    *Reco_QQ_4mom;
TClonesArray    *Reco_QQ_mupl_4mom;
TClonesArray    *Reco_QQ_mumi_4mom;
ULong64_t       Reco_QQ_trig[99];   //[Reco_QQ_size]
ULong64_t       Reco_QQ_mupl_trig[99];   //[Reco_QQ_size]
ULong64_t       Reco_QQ_mumi_trig[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_isCowboy[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctau[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctau3D[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauErr3D[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauTrue[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_ctauTrue3D[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_VtxProb[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_dca[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_MassErr[99];   //[Reco_QQ_size]
TClonesArray    *Reco_QQ_vtx;
Int_t           Reco_QQ_Ntrk[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_SelectionType[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_SelectionType[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_isGoodMuon[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_isGoodMuon[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_highPurity[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_highPurity[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_TrkMuArb[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_TrkMuArb[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mupl_TMOneStaTight[99];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_TMOneStaTight[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nPixValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nPixValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nMuValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nMuValHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nTrkHits[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nTrkHits[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_normChi2_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_normChi2_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_normChi2_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_normChi2_global[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nPixWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nPixWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_nTrkWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nTrkWMea[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mupl_StationsMatched[99];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_StationsMatched[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dxy[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dxy[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dxyErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dxyErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dz[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dz[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dzErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dzErr[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_pt_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_pt_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_pt_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_pt_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_ptErr_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_ptErr_inner[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_ptErr_global[99];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_ptErr_global[99];   //[Reco_QQ_size]
Int_t           Reco_mu_size;
Int_t           Reco_mu_type[99];   //[Reco_mu_size]
Int_t           Reco_mu_charge[99];   //[Reco_mu_size]
TClonesArray    *Reco_mu_4mom;
ULong64_t       Reco_mu_trig[99];   //[Reco_mu_size]
Bool_t          Reco_mu_isGoodMuon[99];   //[Reco_mu_size]
Bool_t          Reco_mu_highPurity[99];   //[Reco_mu_size]
Bool_t          Reco_mu_TrkMuArb[99];   //[Reco_mu_size]
Bool_t          Reco_mu_TMOneStaTight[99];   //[Reco_mu_size]
Int_t           Reco_mu_nPixValHits[99];   //[Reco_mu_size]
Int_t           Reco_mu_nMuValHits[99];   //[Reco_mu_size]
Int_t           Reco_mu_nTrkHits[99];   //[Reco_mu_size]
Float_t         Reco_mu_normChi2_inner[99];   //[Reco_mu_size]
Float_t         Reco_mu_normChi2_global[99];   //[Reco_mu_size]
Int_t           Reco_mu_nPixWMea[99];   //[Reco_mu_size]
Int_t           Reco_mu_nTrkWMea[99];   //[Reco_mu_size]
Int_t           Reco_mu_StationsMatched[99];   //[Reco_mu_size]
Float_t         Reco_mu_dxy[99];   //[Reco_mu_size]
Float_t         Reco_mu_dxyErr[99];   //[Reco_mu_size]
Float_t         Reco_mu_dz[99];   //[Reco_mu_size]
Float_t         Reco_mu_dzErr[99];   //[Reco_mu_size]
Float_t         Reco_mu_pt_inner[99];   //[Reco_mu_size]
Float_t         Reco_mu_pt_global[99];   //[Reco_mu_size]
Float_t         Reco_mu_ptErr_inner[99];   //[Reco_mu_size]
Float_t         Reco_mu_ptErr_global[99];   //[Reco_mu_size]
Int_t           Gen_QQ_size;
Int_t           Gen_QQ_type[99];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_4mom;
Int_t           Gen_QQ_momId[99];   //[Gen_QQ_size]
Float_t         Gen_QQ_ctau[99];   //[Gen_QQ_size]
Float_t         Gen_QQ_ctau3D[99];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_mupl_4mom;
TClonesArray    *Gen_QQ_mumi_4mom;
Int_t           Gen_mu_size;
Int_t           Gen_mu_type[99];   //[Gen_mu_size]
Int_t           Gen_mu_charge[99];   //[Gen_mu_size]
TClonesArray    *Gen_mu_4mom;

Int_t           evt;
Float_t         b;
Int_t           nref;
Float_t         rawpt[56];   //[nref]
Float_t         jtpt[56];   //[nref]
Float_t         jteta[56];   //[nref]
Float_t         jty[56];   //[nref]
Float_t         jtphi[56];   //[nref]
Float_t         jtpu[56];   //[nref]
Float_t         jtm[56];   //[nref]
Float_t         jtarea[56];   //[nref]
Float_t         jtPfCHF[56];   //[nref]
Float_t         jtPfNHF[56];   //[nref]
Float_t         jtPfCEF[56];   //[nref]
Float_t         jtPfNEF[56];   //[nref]
Float_t         jtPfMUF[56];   //[nref]
Int_t           jtPfCHM[56];   //[nref]
Int_t           jtPfNHM[56];   //[nref]
Int_t           jtPfCEM[56];   //[nref]
Int_t           jtPfNEM[56];   //[nref]
Int_t           jtPfMUM[56];   //[nref]
Float_t         jttau1[56];   //[nref]
Float_t         jttau2[56];   //[nref]
Float_t         jttau3[56];   //[nref]
Float_t         discr_jetID_cuts[56];   //[nref]
Float_t         discr_jetID_bdt[56];   //[nref]
Float_t         discr_fr01[56];   //[nref]
Float_t         trackMax[56];   //[nref]
Float_t         trackSum[56];   //[nref]
Int_t           trackN[56];   //[nref]
Float_t         trackHardSum[56];   //[nref]
Int_t           trackHardN[56];   //[nref]
Float_t         chargedMax[56];   //[nref]
Float_t         chargedSum[56];   //[nref]
Int_t           chargedN[56];   //[nref]
Float_t         chargedHardSum[56];   //[nref]
Int_t           chargedHardN[56];   //[nref]
Float_t         photonMax[56];   //[nref]
Float_t         photonSum[56];   //[nref]
Int_t           photonN[56];   //[nref]
Float_t         photonHardSum[56];   //[nref]
Int_t           photonHardN[56];   //[nref]
Float_t         neutralMax[56];   //[nref]
Float_t         neutralSum[56];   //[nref]
Int_t           neutralN[56];   //[nref]
Float_t         hcalSum[56];   //[nref]
Float_t         ecalSum[56];   //[nref]
Float_t         eMax[56];   //[nref]
Float_t         eSum[56];   //[nref]
Int_t           eN[56];   //[nref]
Float_t         muMax[56];   //[nref]
Float_t         muSum[56];   //[nref]
Int_t           muN[56];   //[nref]
Float_t         discr_ssvHighEff[56];   //[nref]
Float_t         discr_ssvHighPur[56];   //[nref]
Float_t         discr_csvV1[56];   //[nref]
Float_t         discr_csvV2[56];   //[nref]
Float_t         discr_muByIp3[56];   //[nref]
Float_t         discr_muByPt[56];   //[nref]
Float_t         discr_prob[56];   //[nref]
Float_t         discr_probb[56];   //[nref]
Float_t         discr_tcHighEff[56];   //[nref]
Float_t         discr_tcHighPur[56];   //[nref]
Float_t         ndiscr_ssvHighEff[56];   //[nref]
Float_t         ndiscr_ssvHighPur[56];   //[nref]
Float_t         ndiscr_csvV1[56];   //[nref]
Float_t         ndiscr_csvV2[56];   //[nref]
Float_t         ndiscr_muByPt[56];   //[nref]
Float_t         pdiscr_csvV1[56];   //[nref]
Float_t         pdiscr_csvV2[56];   //[nref]
Int_t           nsvtx[56];   //[nref]
Int_t           svtxntrk[56];   //[nref]
Float_t         svtxdl[56];   //[nref]
Float_t         svtxdls[56];   //[nref]
Float_t         svtxdl2d[56];   //[nref]
Float_t         svtxdls2d[56];   //[nref]
Float_t         svtxm[56];   //[nref]
Float_t         svtxpt[56];   //[nref]
Float_t         svtxmcorr[56];   //[nref]
Int_t           nIPtrk[56];   //[nref]
Int_t           nselIPtrk[56];   //[nref]
Float_t         mue[56];   //[nref]
Float_t         mupt[56];   //[nref]
Float_t         mueta[56];   //[nref]
Float_t         muphi[56];   //[nref]
Float_t         mudr[56];   //[nref]
Float_t         muptrel[56];   //[nref]
Int_t           muchg[56];   //[nref]
Int_t           beamId1;
Int_t           beamId2;
Float_t         pthat;
Float_t         refpt[56];   //[nref]
Float_t         refeta[56];   //[nref]
Float_t         refy[56];   //[nref]
Float_t         refphi[56];   //[nref]
Float_t         refm[56];   //[nref]
Float_t         refarea[56];   //[nref]
Float_t         reftau1[56];   //[nref]
Float_t         reftau2[56];   //[nref]
Float_t         reftau3[56];   //[nref]
Float_t         refdphijt[56];   //[nref]
Float_t         refdrjt[56];   //[nref]
Float_t         refparton_pt[56];   //[nref]
Int_t           refparton_flavor[56];   //[nref]
Int_t           refparton_flavorForB[56];   //[nref]
Float_t         genChargedSum[56];   //[nref]
Float_t         genHardSum[56];   //[nref]
Float_t         signalChargedSum[56];   //[nref]
Float_t         signalHardSum[56];   //[nref]
Int_t           subid[56];   //[nref]
Int_t           ngen;
Int_t           genmatchindex[28];   //[ngen]
Float_t         genpt[28];   //[ngen]
Float_t         geneta[28];   //[ngen]
Float_t         geny[28];   //[ngen]
Float_t         gentau1[28];   //[ngen]
Float_t         gentau2[28];   //[ngen]
Float_t         gentau3[28];   //[ngen]
Float_t         genphi[28];   //[ngen]
Float_t         genm[28];   //[ngen]
Float_t         gendphijt[28];   //[ngen]
Float_t         gendrjt[28];   //[ngen]
Int_t           gensubid[28];   //[ngen]

Int_t           Onia2MuMuPAT;
Int_t           ana_step;
Int_t           pHBHENoiseFilterResultProducer;
Int_t           HBHENoiseFilterResult;
Int_t           HBHENoiseFilterResultRun1;
Int_t           HBHENoiseFilterResultRun2Loose;
Int_t           HBHENoiseFilterResultRun2Tight;
Int_t           HBHEIsoNoiseFilterResult;
Int_t           pPAprimaryVertexFilter;
Int_t           pBeamScrapingFilter;
Int_t           pVertexFilterCutG;
Int_t           pVertexFilterCutGloose;
Int_t           pVertexFilterCutGtight;
Int_t           pVertexFilterCutGplus;
Int_t           pVertexFilterCutE;
Int_t           pVertexFilterCutEandG;


// List of branches
TBranch        *b_eventNb;   //!
TBranch        *b_runNb;   //!
TBranch        *b_LS;   //!
TBranch        *b_zVtx;   //!
TBranch        *b_Centrality;   //!
TBranch        *b_nTrig;   //!
TBranch        *b_trigPrescale;   //!
TBranch        *b_HLTriggers;   //!
TBranch        *b_Npix;   //!
TBranch        *b_NpixelTracks;   //!
TBranch        *b_Ntracks;   //!
TBranch        *b_SumET_HF;   //!
TBranch        *b_SumET_HFplus;   //!
TBranch        *b_SumET_HFminus;   //!
TBranch        *b_SumET_HFplusEta4;   //!
TBranch        *b_SumET_HFminusEta4;   //!
TBranch        *b_SumET_ET;   //!
TBranch        *b_SumET_EE;   //!
TBranch        *b_SumET_EB;   //!
TBranch        *b_SumET_EEplus;   //!
TBranch        *b_SumET_EEminus;   //!
TBranch        *b_SumET_ZDC;   //!
TBranch        *b_SumET_ZDCplus;   //!
TBranch        *b_SumET_ZDCminus;   //!
TBranch        *b_nEP;   //!
TBranch        *b_rpAng;   //!
TBranch        *b_rpSin;   //!
TBranch        *b_rpCos;   //!
TBranch        *b_Reco_QQ_size;   //!
TBranch        *b_Reco_QQ_type;   //!
TBranch        *b_Reco_QQ_sign;   //!
TBranch        *b_Reco_QQ_4mom;   //!
TBranch        *b_Reco_QQ_mupl_4mom;   //!
TBranch        *b_Reco_QQ_mumi_4mom;   //!
TBranch        *b_Reco_QQ_trig;   //!
TBranch        *b_Reco_QQ_mupl_trig;   //!
TBranch        *b_Reco_QQ_mumi_trig;   //!
TBranch        *b_Reco_QQ_isCowboy;   //!
TBranch        *b_Reco_QQ_ctau;   //!
TBranch        *b_Reco_QQ_ctauErr;   //!
TBranch        *b_Reco_QQ_ctau3D;   //!
TBranch        *b_Reco_QQ_ctauErr3D;   //!
TBranch        *b_Reco_QQ_ctauTrue;   //!
TBranch        *b_Reco_QQ_ctauTrue3D;   //!
TBranch        *b_Reco_QQ_VtxProb;   //!
TBranch        *b_Reco_QQ_dca;   //!
TBranch        *b_Reco_QQ_MassErr;   //!
TBranch        *b_Reco_QQ_vtx;   //!
TBranch        *b_Reco_QQ_Ntrk;   //!
TBranch        *b_Reco_QQ_mupl_SelectionType; //!
TBranch        *b_Reco_QQ_mumi_SelectionType; //!
TBranch        *b_Reco_QQ_mupl_isGoodMuon;   //!
TBranch        *b_Reco_QQ_mumi_isGoodMuon;   //!
TBranch        *b_Reco_QQ_mupl_highPurity;   //!
TBranch        *b_Reco_QQ_mumi_highPurity;   //!
TBranch        *b_Reco_QQ_mupl_TrkMuArb;   //!
TBranch        *b_Reco_QQ_mumi_TrkMuArb;   //!
TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
TBranch        *b_Reco_QQ_mupl_normChi2_inner;   //!
TBranch        *b_Reco_QQ_mumi_normChi2_inner;   //!
TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
TBranch        *b_Reco_QQ_mupl_dxy;   //!
TBranch        *b_Reco_QQ_mumi_dxy;   //!
TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
TBranch        *b_Reco_QQ_mupl_dz;   //!
TBranch        *b_Reco_QQ_mumi_dz;   //!
TBranch        *b_Reco_QQ_mupl_dzErr;   //!
TBranch        *b_Reco_QQ_mumi_dzErr;   //!
TBranch        *b_Reco_QQ_mupl_pt_inner;   //!
TBranch        *b_Reco_QQ_mumi_pt_inner;   //!
TBranch        *b_Reco_QQ_mupl_pt_global;   //!
TBranch        *b_Reco_QQ_mumi_pt_global;   //!
TBranch        *b_Reco_QQ_mupl_ptErr_inner;   //!
TBranch        *b_Reco_QQ_mumi_ptErr_inner;   //!
TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
TBranch        *b_Reco_mu_size;   //!
TBranch        *b_Reco_mu_type;   //!
TBranch        *b_Reco_mu_charge;   //!
TBranch        *b_Reco_mu_4mom;   //!
TBranch        *b_Reco_mu_trig;   //!
TBranch        *b_Reco_mu_isGoodMuon;   //!
TBranch        *b_Reco_mu_highPurity;   //!
TBranch        *b_Reco_mu_TrkMuArb;   //!
TBranch        *b_Reco_mu_TMOneStaTight;   //!
TBranch        *b_Reco_mu_nPixValHits;   //!
TBranch        *b_Reco_mu_nMuValHits;   //!
TBranch        *b_Reco_mu_nTrkHits;   //!
TBranch        *b_Reco_mu_normChi2_inner;   //!
TBranch        *b_Reco_mu_normChi2_global;   //!
TBranch        *b_Reco_mu_nPixWMea;   //!
TBranch        *b_Reco_mu_nTrkWMea;   //!
TBranch        *b_Reco_mu_StationsMatched;   //!
TBranch        *b_Reco_mu_dxy;   //!
TBranch        *b_Reco_mu_dxyErr;   //!
TBranch        *b_Reco_mu_dz;   //!
TBranch        *b_Reco_mu_dzErr;   //!
TBranch        *b_Reco_mu_pt_inner;   //!
TBranch        *b_Reco_mu_pt_global;   //!
TBranch        *b_Reco_mu_ptErr_inner;   //!
TBranch        *b_Reco_mu_ptErr_global;   //!
TBranch        *b_Gen_QQ_size;   //!
TBranch        *b_Gen_QQ_type;   //!
TBranch        *b_Gen_QQ_4mom;   //!
TBranch        *b_Gen_QQ_momId;   //!
TBranch        *b_Gen_QQ_ctau;   //!
TBranch        *b_Gen_QQ_ctau3D;   //!
TBranch        *b_Gen_QQ_mupl_4mom;   //!
TBranch        *b_Gen_QQ_mumi_4mom;   //!
TBranch        *b_Gen_mu_size;   //!
TBranch        *b_Gen_mu_type;   //!
TBranch        *b_Gen_mu_charge;   //!
TBranch        *b_Gen_mu_4mom;   //!

TBranch        *b_evt;   //!
TBranch        *b_b;   //!
TBranch        *b_nref;   //!
TBranch        *b_rawpt;   //!
TBranch        *b_jtpt;   //!
TBranch        *b_jteta;   //!
TBranch        *b_jty;   //!
TBranch        *b_jtphi;   //!
TBranch        *b_jtpu;   //!
TBranch        *b_jtm;   //!
TBranch        *b_jtarea;   //!
TBranch        *b_jtPfCHF;   //!
TBranch        *b_jtPfNHF;   //!
TBranch        *b_jtPfCEF;   //!
TBranch        *b_jtPfNEF;   //!
TBranch        *b_jtPfMUF;   //!
TBranch        *b_jtPfCHM;   //!
TBranch        *b_jtPfNHM;   //!
TBranch        *b_jtPfCEM;   //!
TBranch        *b_jtPfNEM;   //!
TBranch        *b_jtPfMUM;   //!
TBranch        *b_jttau1;   //!
TBranch        *b_jttau2;   //!
TBranch        *b_jttau3;   //!
TBranch        *b_discr_jetID_cuts;   //!
TBranch        *b_discr_jetID_bdt;   //!
TBranch        *b_discr_fr01;   //!
TBranch        *b_trackMax;   //!
TBranch        *b_trackSum;   //!
TBranch        *b_trackN;   //!
TBranch        *b_trackHardSum;   //!
TBranch        *b_trackHardN;   //!
TBranch        *b_chargedMax;   //!
TBranch        *b_chargedSum;   //!
TBranch        *b_chargedN;   //!
TBranch        *b_chargedHardSum;   //!
TBranch        *b_chargedHardN;   //!
TBranch        *b_photonMax;   //!
TBranch        *b_photonSum;   //!
TBranch        *b_photonN;   //!
TBranch        *b_photonHardSum;   //!
TBranch        *b_photonHardN;   //!
TBranch        *b_neutralMax;   //!
TBranch        *b_neutralSum;   //!
TBranch        *b_neutralN;   //!
TBranch        *b_hcalSum;   //!
TBranch        *b_ecalSum;   //!
TBranch        *b_eMax;   //!
TBranch        *b_eSum;   //!
TBranch        *b_eN;   //!
TBranch        *b_muMax;   //!
TBranch        *b_muSum;   //!
TBranch        *b_muN;   //!
TBranch        *b_discr_ssvHighEff;   //!
TBranch        *b_discr_ssvHighPur;   //!
TBranch        *b_discr_csvV1;   //!
TBranch        *b_discr_csvV2;   //!
TBranch        *b_discr_muByIp3;   //!
TBranch        *b_discr_muByPt;   //!
TBranch        *b_discr_prob;   //!
TBranch        *b_discr_probb;   //!
TBranch        *b_discr_tcHighEff;   //!
TBranch        *b_discr_tcHighPur;   //!
TBranch        *b_ndiscr_ssvHighEff;   //!
TBranch        *b_ndiscr_ssvHighPur;   //!
TBranch        *b_ndiscr_csvV1;   //!
TBranch        *b_ndiscr_csvV2;   //!
TBranch        *b_ndiscr_muByPt;   //!
TBranch        *b_pdiscr_csvV1;   //!
TBranch        *b_pdiscr_csvV2;   //!
TBranch        *b_nsvtx;   //!
TBranch        *b_svtxntrk;   //!
TBranch        *b_svtxdl;   //!
TBranch        *b_svtxdls;   //!
TBranch        *b_svtxdl2d;   //!
TBranch        *b_svtxdls2d;   //!
TBranch        *b_svtxm;   //!
TBranch        *b_svtxpt;   //!
TBranch        *b_svtxmcorr;   //!
TBranch        *b_nIPtrk;   //!
TBranch        *b_nselIPtrk;   //!
TBranch        *b_mue;   //!
TBranch        *b_mupt;   //!
TBranch        *b_mueta;   //!
TBranch        *b_muphi;   //!
TBranch        *b_mudr;   //!
TBranch        *b_muptrel;   //!
TBranch        *b_muchg;   //!

TBranch        *b_beamId1;   //!
TBranch        *b_beamId2;   //!
TBranch        *b_pthat;   //!
TBranch        *b_refpt;   //!
TBranch        *b_refeta;   //!
TBranch        *b_refy;   //!
TBranch        *b_refphi;   //!
TBranch        *b_refm;   //!
TBranch        *b_refarea;   //!
TBranch        *b_reftau1;   //!
TBranch        *b_reftau2;   //!
TBranch        *b_reftau3;   //!
TBranch        *b_refdphijt;   //!
TBranch        *b_refdrjt;   //!
TBranch        *b_refparton_pt;   //!
TBranch        *b_refparton_flavor;   //!
TBranch        *b_refparton_flavorForB;   //!
TBranch        *b_genChargedSum;   //!
TBranch        *b_genHardSum;   //!
TBranch        *b_signalChargedSum;   //!
TBranch        *b_signalHardSum;   //!
TBranch        *b_subid;   //!
TBranch        *b_ngen;   //!
TBranch        *b_genmatchindex;   //!
TBranch        *b_genpt;   //!
TBranch        *b_geneta;   //!
TBranch        *b_geny;   //!
TBranch        *b_gentau1;   //!
TBranch        *b_gentau2;   //!
TBranch        *b_gentau3;   //!
TBranch        *b_genphi;   //!
TBranch        *b_genm;   //!
TBranch        *b_gendphijt;   //!
TBranch        *b_gendrjt;   //!
TBranch        *b_gensubid;   //!

TBranch        *b_Onia2MuMuPAT;   //!                                                                                                                                                                
TBranch        *b_ana_step;   //!                                                                                                                                                                    
TBranch        *b_pHBHENoiseFilterResultProducer;   //!                                                                                                                                              
TBranch        *b_HBHENoiseFilterResult;   //!                                                                                                                                                       
TBranch        *b_HBHENoiseFilterResultRun1;   //!                                                                                                                                                   
TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!                                                                                                                                              
TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!                                                                                                                                              
TBranch        *b_HBHEIsoNoiseFilterResult;   //!                                                                                                                                                    
TBranch        *b_pPAprimaryVertexFilter;   //!                                                                                                                                                      
TBranch        *b_pBeamScrapingFilter;   //!                                                                                                                                                         
TBranch        *b_pVertexFilterCutG;   //!                                                                                                                                                           
TBranch        *b_pVertexFilterCutGloose;   //!                                                                                                                                                      
TBranch        *b_pVertexFilterCutGtight;   //!                                                                                                                                                      
TBranch        *b_pVertexFilterCutGplus;   //!                                                                                                                                                       
TBranch        *b_pVertexFilterCutE;   //!                                                                                                                                                           
TBranch        *b_pVertexFilterCutEandG;   //!


string TreeName("hionia/myTree");
string jetTreeName("ak4PFJetAnalyzer/t");
string skimTreeName("skimanalysis/HltTree");

TTree *htr;
TTree *jtr;
TTree *str;

void initTree(TChain *tree)
{
  std::cout << "[INFO] Initializing TTree " << TreeName.c_str() << std::endl;
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed). 

  TChain   *fChain;   //!pointer to the analyzed TTree or TChain

  // Set object pointer
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  Reco_QQ_vtx = 0;
  Reco_mu_4mom = 0;
  Gen_QQ_4mom = 0;
  Gen_QQ_mupl_4mom = 0;
  Gen_QQ_mumi_4mom = 0;
  Gen_mu_4mom = 0;
  // Set branch addresses and branch pointers

  if (!tree) return;
  fChain = tree;
  fCurrent = -1;

  if (fChain->GetBranch("eventNb")) fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  if (fChain->GetBranch("runNb")) fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
  if (fChain->GetBranch("LS")) fChain->SetBranchAddress("LS", &LS, &b_LS);
  if (fChain->GetBranch("zVtx")) fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  if (fChain->GetBranch("Centrality")) fChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  if (fChain->GetBranch("nTrig")) fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
  if (fChain->GetBranch("trigPrescale")) fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
  if (fChain->GetBranch("HLTriggers")) fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  if (fChain->GetBranch("Npix")) fChain->SetBranchAddress("Npix", &Npix, &b_Npix);
  if (fChain->GetBranch("NpixelTracks")) fChain->SetBranchAddress("NpixelTracks", &NpixelTracks, &b_NpixelTracks);
  if (fChain->GetBranch("Ntracks")) fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
  if (fChain->GetBranch("SumET_HF")) fChain->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
  if (fChain->GetBranch("SumET_HFplus")) fChain->SetBranchAddress("SumET_HFplus", &SumET_HFplus, &b_SumET_HFplus);
  if (fChain->GetBranch("SumET_HFminus")) fChain->SetBranchAddress("SumET_HFminus", &SumET_HFminus, &b_SumET_HFminus);
  if (fChain->GetBranch("SumET_HFplusEta4")) fChain->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4, &b_SumET_HFplusEta4);
  if (fChain->GetBranch("SumET_HFminusEta4")) fChain->SetBranchAddress("SumET_HFminusEta4", &SumET_HFminusEta4, &b_SumET_HFminusEta4);
  if (fChain->GetBranch("SumET_ET")) fChain->SetBranchAddress("SumET_ET", &SumET_ET, &b_SumET_ET);
  if (fChain->GetBranch("SumET_EE")) fChain->SetBranchAddress("SumET_EE", &SumET_EE, &b_SumET_EE);
  if (fChain->GetBranch("SumET_EB")) fChain->SetBranchAddress("SumET_EB", &SumET_EB, &b_SumET_EB);
  if (fChain->GetBranch("SumET_EEplus")) fChain->SetBranchAddress("SumET_EEplus", &SumET_EEplus, &b_SumET_EEplus);
  if (fChain->GetBranch("SumET_EEminus")) fChain->SetBranchAddress("SumET_EEminus", &SumET_EEminus, &b_SumET_EEminus);
  if (fChain->GetBranch("SumET_ZDC")) fChain->SetBranchAddress("SumET_ZDC", &SumET_ZDC, &b_SumET_ZDC);
  if (fChain->GetBranch("SumET_ZDCplus")) fChain->SetBranchAddress("SumET_ZDCplus", &SumET_ZDCplus, &b_SumET_ZDCplus);
  if (fChain->GetBranch("SumET_ZDCminus")) fChain->SetBranchAddress("SumET_ZDCminus", &SumET_ZDCminus, &b_SumET_ZDCminus);
  if (fChain->GetBranch("nEP")) fChain->SetBranchAddress("nEP", &nEP, &b_nEP);
  if (fChain->GetBranch("rpAng")) fChain->SetBranchAddress("rpAng", rpAng, &b_rpAng);
  if (fChain->GetBranch("rpSin")) fChain->SetBranchAddress("rpSin", rpSin, &b_rpSin);
  if (fChain->GetBranch("rpCos")) fChain->SetBranchAddress("rpCos", rpCos, &b_rpCos);
  if (fChain->GetBranch("Reco_QQ_size")) fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  if (fChain->GetBranch("Reco_QQ_type")) fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
  if (fChain->GetBranch("Reco_QQ_sign")) fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  if (fChain->GetBranch("Reco_QQ_4mom")) fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) fChain->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) fChain->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  if (fChain->GetBranch("Reco_QQ_trig")) fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  if (fChain->GetBranch("Reco_QQ_mupl_trig")) fChain->SetBranchAddress("Reco_QQ_mupl_trig", Reco_QQ_mupl_trig, &b_Reco_QQ_mupl_trig);
  if (fChain->GetBranch("Reco_QQ_mumi_trig")) fChain->SetBranchAddress("Reco_QQ_mumi_trig", Reco_QQ_mumi_trig, &b_Reco_QQ_mumi_trig);
  if (fChain->GetBranch("Reco_QQ_isCowboy")) fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
  if (fChain->GetBranch("Reco_QQ_ctau")) fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
  if (fChain->GetBranch("Reco_QQ_ctauErr")) fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
  if (fChain->GetBranch("Reco_QQ_ctau3D")) fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
  if (fChain->GetBranch("Reco_QQ_ctauTrue")) fChain->SetBranchAddress("Reco_QQ_ctauTrue", Reco_QQ_ctauTrue, &b_Reco_QQ_ctauTrue);
  if (fChain->GetBranch("Reco_QQ_ctauTrue3D")) fChain->SetBranchAddress("Reco_QQ_ctauTrue3D", Reco_QQ_ctauTrue3D, &b_Reco_QQ_ctauTrue3D);
  if (fChain->GetBranch("Reco_QQ_VtxProb")) fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  if (fChain->GetBranch("Reco_QQ_dca")) fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
  if (fChain->GetBranch("Reco_QQ_MassErr")) fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
  if (fChain->GetBranch("Reco_QQ_vtx")) fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
  if (fChain->GetBranch("Reco_QQ_Ntrk")) fChain->SetBranchAddress("Reco_QQ_Ntrk", Reco_QQ_Ntrk, &b_Reco_QQ_Ntrk);
  if (fChain->GetBranch("Reco_QQ_mupl_SelectionType")) fChain->SetBranchAddress("Reco_QQ_mupl_SelectionType", Reco_QQ_mupl_SelectionType, &b_Reco_QQ_mupl_SelectionType);
  if (fChain->GetBranch("Reco_QQ_mumi_SelectionType")) fChain->SetBranchAddress("Reco_QQ_mumi_SelectionType", Reco_QQ_mumi_SelectionType, &b_Reco_QQ_mumi_SelectionType);
  if (fChain->GetBranch("Reco_QQ_mupl_isGoodMuon")) fChain->SetBranchAddress("Reco_QQ_mupl_isGoodMuon", Reco_QQ_mupl_isGoodMuon, &b_Reco_QQ_mupl_isGoodMuon);
  if (fChain->GetBranch("Reco_QQ_mumi_isGoodMuon")) fChain->SetBranchAddress("Reco_QQ_mumi_isGoodMuon", Reco_QQ_mumi_isGoodMuon, &b_Reco_QQ_mumi_isGoodMuon);
  if (fChain->GetBranch("Reco_QQ_mupl_highPurity")) fChain->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
  if (fChain->GetBranch("Reco_QQ_mumi_highPurity")) fChain->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);
  if (fChain->GetBranch("Reco_QQ_mupl_TrkMuArb")) fChain->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
  if (fChain->GetBranch("Reco_QQ_mumi_TrkMuArb")) fChain->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);
  if (fChain->GetBranch("Reco_QQ_mupl_TMOneStaTight")) fChain->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  if (fChain->GetBranch("Reco_QQ_mumi_TMOneStaTight")) fChain->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
  if (fChain->GetBranch("Reco_QQ_mupl_nPixValHits")) fChain->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
  if (fChain->GetBranch("Reco_QQ_mumi_nPixValHits")) fChain->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
  if (fChain->GetBranch("Reco_QQ_mupl_nMuValHits")) fChain->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
  if (fChain->GetBranch("Reco_QQ_mumi_nMuValHits")) fChain->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
  if (fChain->GetBranch("Reco_QQ_mupl_nTrkHits")) fChain->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  if (fChain->GetBranch("Reco_QQ_mumi_nTrkHits")) fChain->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  if (fChain->GetBranch("Reco_QQ_mupl_normChi2_inner")) fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", Reco_QQ_mupl_normChi2_inner, &b_Reco_QQ_mupl_normChi2_inner);
  if (fChain->GetBranch("Reco_QQ_mumi_normChi2_inner")) fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", Reco_QQ_mumi_normChi2_inner, &b_Reco_QQ_mumi_normChi2_inner);
  if (fChain->GetBranch("Reco_QQ_mupl_normChi2_global")) fChain->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  if (fChain->GetBranch("Reco_QQ_mumi_normChi2_global")) fChain->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  if (fChain->GetBranch("Reco_QQ_mupl_nPixWMea")) fChain->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  if (fChain->GetBranch("Reco_QQ_mumi_nPixWMea")) fChain->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  if (fChain->GetBranch("Reco_QQ_mupl_nTrkWMea")) fChain->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  if (fChain->GetBranch("Reco_QQ_mumi_nTrkWMea")) fChain->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  if (fChain->GetBranch("Reco_QQ_mupl_StationsMatched")) fChain->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
  if (fChain->GetBranch("Reco_QQ_mumi_StationsMatched")) fChain->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
  if (fChain->GetBranch("Reco_QQ_mupl_dxy")) fChain->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  if (fChain->GetBranch("Reco_QQ_mumi_dxy")) fChain->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  if (fChain->GetBranch("Reco_QQ_mupl_dxyErr")) fChain->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  if (fChain->GetBranch("Reco_QQ_mumi_dxyErr")) fChain->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  if (fChain->GetBranch("Reco_QQ_mupl_dz")) fChain->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  if (fChain->GetBranch("Reco_QQ_mumi_dz")) fChain->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  if (fChain->GetBranch("Reco_QQ_mupl_dzErr")) fChain->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  if (fChain->GetBranch("Reco_QQ_mumi_dzErr")) fChain->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  if (fChain->GetBranch("Reco_QQ_mupl_pt_inner")) fChain->SetBranchAddress("Reco_QQ_mupl_pt_inner", Reco_QQ_mupl_pt_inner, &b_Reco_QQ_mupl_pt_inner);
  if (fChain->GetBranch("Reco_QQ_mumi_pt_inner")) fChain->SetBranchAddress("Reco_QQ_mumi_pt_inner", Reco_QQ_mumi_pt_inner, &b_Reco_QQ_mumi_pt_inner);
  if (fChain->GetBranch("Reco_QQ_mupl_pt_global")) fChain->SetBranchAddress("Reco_QQ_mupl_pt_global", Reco_QQ_mupl_pt_global, &b_Reco_QQ_mupl_pt_global);
  if (fChain->GetBranch("Reco_QQ_mumi_pt_global")) fChain->SetBranchAddress("Reco_QQ_mumi_pt_global", Reco_QQ_mumi_pt_global, &b_Reco_QQ_mumi_pt_global);
  if (fChain->GetBranch("Reco_QQ_mupl_ptErr_inner")) fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_inner", Reco_QQ_mupl_ptErr_inner, &b_Reco_QQ_mupl_ptErr_inner);
  if (fChain->GetBranch("Reco_QQ_mumi_ptErr_inner")) fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_inner", Reco_QQ_mumi_ptErr_inner, &b_Reco_QQ_mumi_ptErr_inner);
  if (fChain->GetBranch("Reco_QQ_mupl_ptErr_global")) fChain->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
  if (fChain->GetBranch("Reco_QQ_mumi_ptErr_global")) fChain->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
  if (fChain->GetBranch("Reco_mu_size")) fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  if (fChain->GetBranch("Reco_mu_type")) fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
  if (fChain->GetBranch("Reco_mu_charge")) fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
  if (fChain->GetBranch("Reco_mu_4mom")) fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  if (fChain->GetBranch("Reco_mu_trig")) fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  if (fChain->GetBranch("Reco_mu_isGoodMuon")) fChain->SetBranchAddress("Reco_mu_isGoodMuon", Reco_mu_isGoodMuon, &b_Reco_mu_isGoodMuon);
  if (fChain->GetBranch("Reco_mu_highPurity")) fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  if (fChain->GetBranch("Reco_mu_TrkMuArb")) fChain->SetBranchAddress("Reco_mu_TrkMuArb", Reco_mu_TrkMuArb, &b_Reco_mu_TrkMuArb);
  if (fChain->GetBranch("Reco_mu_TMOneStaTight")) fChain->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  if (fChain->GetBranch("Reco_mu_nPixValHits")) fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  if (fChain->GetBranch("Reco_mu_nMuValHits")) fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  if (fChain->GetBranch("Reco_mu_nTrkHits")) fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  if (fChain->GetBranch("Reco_mu_normChi2_inner")) fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
  if (fChain->GetBranch("Reco_mu_normChi2_global")) fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  if (fChain->GetBranch("Reco_mu_nPixWMea")) fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  if (fChain->GetBranch("Reco_mu_nTrkWMea")) fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  if (fChain->GetBranch("Reco_mu_StationsMatched")) fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  if (fChain->GetBranch("Reco_mu_dxy")) fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  if (fChain->GetBranch("Reco_mu_dxyErr")) fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  if (fChain->GetBranch("Reco_mu_dz")) fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  if (fChain->GetBranch("Reco_mu_dzErr")) fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  if (fChain->GetBranch("Reco_mu_pt_inner")) fChain->SetBranchAddress("Reco_mu_pt_inner", Reco_mu_pt_inner, &b_Reco_mu_pt_inner);
  if (fChain->GetBranch("Reco_mu_pt_global")) fChain->SetBranchAddress("Reco_mu_pt_global", Reco_mu_pt_global, &b_Reco_mu_pt_global);
  if (fChain->GetBranch("Reco_mu_ptErr_inner")) fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
  if (fChain->GetBranch("Reco_mu_ptErr_global")) fChain->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
  if (fChain->GetBranch("Gen_QQ_size")) fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  if (fChain->GetBranch("Gen_QQ_type")) fChain->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
  if (fChain->GetBranch("Gen_QQ_4mom")) fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  if (fChain->GetBranch("Gen_QQ_momId")) fChain->SetBranchAddress("Gen_QQ_momId", Gen_QQ_momId, &b_Gen_QQ_momId);
  if (fChain->GetBranch("Gen_QQ_ctau")) fChain->SetBranchAddress("Gen_QQ_ctau", Gen_QQ_ctau, &b_Gen_QQ_ctau);
  if (fChain->GetBranch("Gen_QQ_ctau3D")) fChain->SetBranchAddress("Gen_QQ_ctau3D", Gen_QQ_ctau3D, &b_Gen_QQ_ctau3D);
  if (fChain->GetBranch("Gen_QQ_mupl_4mom")) fChain->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
  if (fChain->GetBranch("Gen_QQ_mumi_4mom")) fChain->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
  if (fChain->GetBranch("Gen_mu_size")) fChain->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  if (fChain->GetBranch("Gen_mu_type")) fChain->SetBranchAddress("Gen_mu_type", Gen_mu_type, &b_Gen_mu_type);
  if (fChain->GetBranch("Gen_mu_charge")) fChain->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);
  if (fChain->GetBranch("Gen_mu_4mom")) fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  if (fChain->GetBranch("evt")) fChain->SetBranchAddress("evt", &evt, &b_evt);
  if (fChain->GetBranch("b")) fChain->SetBranchAddress("b", &b, &b_b);
  if (fChain->GetBranch("nref")) fChain->SetBranchAddress("nref", &nref, &b_nref);
  if (fChain->GetBranch("rawpt")) fChain->SetBranchAddress("rawpt", rawpt, &b_rawpt);
  if (fChain->GetBranch("jtpt")) fChain->SetBranchAddress("jtpt", jtpt, &b_jtpt);
  if (fChain->GetBranch("jteta")) fChain->SetBranchAddress("jteta", jteta, &b_jteta);
  if (fChain->GetBranch("jty")) fChain->SetBranchAddress("jty", jty, &b_jty);
  if (fChain->GetBranch("jtphi")) fChain->SetBranchAddress("jtphi", jtphi, &b_jtphi);
  if (fChain->GetBranch("jtpu")) fChain->SetBranchAddress("jtpu", jtpu, &b_jtpu);
  if (fChain->GetBranch("jtm")) fChain->SetBranchAddress("jtm", jtm, &b_jtm);
  if (fChain->GetBranch("jtarea")) fChain->SetBranchAddress("jtarea", jtarea, &b_jtarea);
  if (fChain->GetBranch("jtPfCHF")) fChain->SetBranchAddress("jtPfCHF", jtPfCHF, &b_jtPfCHF);
  if (fChain->GetBranch("jtPfNHF")) fChain->SetBranchAddress("jtPfNHF", jtPfNHF, &b_jtPfNHF);
  if (fChain->GetBranch("jtPfCEF")) fChain->SetBranchAddress("jtPfCEF", jtPfCEF, &b_jtPfCEF);
  if (fChain->GetBranch("jtPfNEF")) fChain->SetBranchAddress("jtPfNEF", jtPfNEF, &b_jtPfNEF);
  if (fChain->GetBranch("jtPfMUF")) fChain->SetBranchAddress("jtPfMUF", jtPfMUF, &b_jtPfMUF);
  if (fChain->GetBranch("jtPfCHM")) fChain->SetBranchAddress("jtPfCHM", jtPfCHM, &b_jtPfCHM);
  if (fChain->GetBranch("jtPfNHM")) fChain->SetBranchAddress("jtPfNHM", jtPfNHM, &b_jtPfNHM);
  if (fChain->GetBranch("jtPfCEM")) fChain->SetBranchAddress("jtPfCEM", jtPfCEM, &b_jtPfCEM);
  if (fChain->GetBranch("jtPfNEM")) fChain->SetBranchAddress("jtPfNEM", jtPfNEM, &b_jtPfNEM);
  if (fChain->GetBranch("jtPfMUM")) fChain->SetBranchAddress("jtPfMUM", jtPfMUM, &b_jtPfMUM);
  if (fChain->GetBranch("jttau1")) fChain->SetBranchAddress("jttau1", jttau1, &b_jttau1);
  if (fChain->GetBranch("jttau2")) fChain->SetBranchAddress("jttau2", jttau2, &b_jttau2);
  if (fChain->GetBranch("jttau3")) fChain->SetBranchAddress("jttau3", jttau3, &b_jttau3);
  if (fChain->GetBranch("discr_jetID_cuts")) fChain->SetBranchAddress("discr_jetID_cuts", discr_jetID_cuts, &b_discr_jetID_cuts);
  if (fChain->GetBranch("discr_jetID_bdt")) fChain->SetBranchAddress("discr_jetID_bdt", discr_jetID_bdt, &b_discr_jetID_bdt);
  if (fChain->GetBranch("discr_fr01")) fChain->SetBranchAddress("discr_fr01", discr_fr01, &b_discr_fr01);
  if (fChain->GetBranch("trackMax")) fChain->SetBranchAddress("trackMax", trackMax, &b_trackMax);
  if (fChain->GetBranch("trackSum")) fChain->SetBranchAddress("trackSum", trackSum, &b_trackSum);
  if (fChain->GetBranch("trackN")) fChain->SetBranchAddress("trackN", trackN, &b_trackN);
  if (fChain->GetBranch("trackHardSum")) fChain->SetBranchAddress("trackHardSum", trackHardSum, &b_trackHardSum);
  if (fChain->GetBranch("trackHardN")) fChain->SetBranchAddress("trackHardN", trackHardN, &b_trackHardN);
  if (fChain->GetBranch("chargedMax")) fChain->SetBranchAddress("chargedMax", chargedMax, &b_chargedMax);
  if (fChain->GetBranch("chargedSum")) fChain->SetBranchAddress("chargedSum", chargedSum, &b_chargedSum);
  if (fChain->GetBranch("chargedN")) fChain->SetBranchAddress("chargedN", chargedN, &b_chargedN);
  if (fChain->GetBranch("chargedHardSum")) fChain->SetBranchAddress("chargedHardSum", chargedHardSum, &b_chargedHardSum);
  if (fChain->GetBranch("chargedHardN")) fChain->SetBranchAddress("chargedHardN", chargedHardN, &b_chargedHardN);
  if (fChain->GetBranch("photonMax")) fChain->SetBranchAddress("photonMax", photonMax, &b_photonMax);
  if (fChain->GetBranch("photonSum")) fChain->SetBranchAddress("photonSum", photonSum, &b_photonSum);
  if (fChain->GetBranch("photonN")) fChain->SetBranchAddress("photonN", photonN, &b_photonN);
  if (fChain->GetBranch("photonHardSum")) fChain->SetBranchAddress("photonHardSum", photonHardSum, &b_photonHardSum);
  if (fChain->GetBranch("photonHardN")) fChain->SetBranchAddress("photonHardN", photonHardN, &b_photonHardN);
  if (fChain->GetBranch("neutralMax")) fChain->SetBranchAddress("neutralMax", neutralMax, &b_neutralMax);
  if (fChain->GetBranch("neutralSum")) fChain->SetBranchAddress("neutralSum", neutralSum, &b_neutralSum);
  if (fChain->GetBranch("neutralN")) fChain->SetBranchAddress("neutralN", neutralN, &b_neutralN);
  if (fChain->GetBranch("hcalSum")) fChain->SetBranchAddress("hcalSum", hcalSum, &b_hcalSum);
  if (fChain->GetBranch("ecalSum")) fChain->SetBranchAddress("ecalSum", ecalSum, &b_ecalSum);
  if (fChain->GetBranch("eMax")) fChain->SetBranchAddress("eMax", eMax, &b_eMax);
  if (fChain->GetBranch("eSum")) fChain->SetBranchAddress("eSum", eSum, &b_eSum);
  if (fChain->GetBranch("eN")) fChain->SetBranchAddress("eN", eN, &b_eN);
  if (fChain->GetBranch("muMax")) fChain->SetBranchAddress("muMax", muMax, &b_muMax);
  if (fChain->GetBranch("muSum")) fChain->SetBranchAddress("muSum", muSum, &b_muSum);
  if (fChain->GetBranch("muN")) fChain->SetBranchAddress("muN", muN, &b_muN);
  //if (fChain->GetBranch("discr_ssvHighEff")) fChain->SetBranchAddress("discr_ssvHighEff", discr_ssvHighEff, &b_discr_ssvHighEff);
  //if (fChain->GetBranch("discr_ssvHighPur")) fChain->SetBranchAddress("discr_ssvHighPur", discr_ssvHighPur, &b_discr_ssvHighPur);
  //if (fChain->GetBranch("discr_csvV1")) fChain->SetBranchAddress("discr_csvV1", discr_csvV1, &b_discr_csvV1);
  //if (fChain->GetBranch("discr_csvV2")) fChain->SetBranchAddress("discr_csvV2", discr_csvV2, &b_discr_csvV2);
  //if (fChain->GetBranch("discr_muByIp3")) fChain->SetBranchAddress("discr_muByIp3", discr_muByIp3, &b_discr_muByIp3);
  //if (fChain->GetBranch("discr_muByPt")) fChain->SetBranchAddress("discr_muByPt", discr_muByPt, &b_discr_muByPt);
  //if (fChain->GetBranch("discr_prob")) fChain->SetBranchAddress("discr_prob", discr_prob, &b_discr_prob);
  //if (fChain->GetBranch("discr_probb")) fChain->SetBranchAddress("discr_probb", discr_probb, &b_discr_probb);
  //if (fChain->GetBranch("discr_tcHighEff")) fChain->SetBranchAddress("discr_tcHighEff", discr_tcHighEff, &b_discr_tcHighEff);
  //if (fChain->GetBranch("discr_tcHighPur")) fChain->SetBranchAddress("discr_tcHighPur", discr_tcHighPur, &b_discr_tcHighPur);
  //if (fChain->GetBranch("ndiscr_ssvHighEff")) fChain->SetBranchAddress("ndiscr_ssvHighEff", ndiscr_ssvHighEff, &b_ndiscr_ssvHighEff);
  //if (fChain->GetBranch("ndiscr_ssvHighPur")) fChain->SetBranchAddress("ndiscr_ssvHighPur", ndiscr_ssvHighPur, &b_ndiscr_ssvHighPur);
  //if (fChain->GetBranch("ndiscr_csvV1")) fChain->SetBranchAddress("ndiscr_csvV1", ndiscr_csvV1, &b_ndiscr_csvV1);
  //if (fChain->GetBranch("ndiscr_csvV2")) fChain->SetBranchAddress("ndiscr_csvV2", ndiscr_csvV2, &b_ndiscr_csvV2);
  //if (fChain->GetBranch("ndiscr_muByPt")) fChain->SetBranchAddress("ndiscr_muByPt", ndiscr_muByPt, &b_ndiscr_muByPt);
  //if (fChain->GetBranch("pdiscr_csvV1")) fChain->SetBranchAddress("pdiscr_csvV1", pdiscr_csvV1, &b_pdiscr_csvV1);
  //if (fChain->GetBranch("pdiscr_csvV2")) fChain->SetBranchAddress("pdiscr_csvV2", pdiscr_csvV2, &b_pdiscr_csvV2);
  //if (fChain->GetBranch("nsvtx")) fChain->SetBranchAddress("nsvtx", nsvtx, &b_nsvtx);
  //if (fChain->GetBranch("svtxntrk")) fChain->SetBranchAddress("svtxntrk", svtxntrk, &b_svtxntrk);
  //if (fChain->GetBranch("svtxdl")) fChain->SetBranchAddress("svtxdl", svtxdl, &b_svtxdl);
  //if (fChain->GetBranch("svtxdls")) fChain->SetBranchAddress("svtxdls", svtxdls, &b_svtxdls);
  //if (fChain->GetBranch("svtxdl2d")) fChain->SetBranchAddress("svtxdl2d", svtxdl2d, &b_svtxdl2d);
  //if (fChain->GetBranch("svtxdls2d")) fChain->SetBranchAddress("svtxdls2d", svtxdls2d, &b_svtxdls2d);
  //if (fChain->GetBranch("svtxm")) fChain->SetBranchAddress("svtxm", svtxm, &b_svtxm);
  //if (fChain->GetBranch("svtxpt")) fChain->SetBranchAddress("svtxpt", svtxpt, &b_svtxpt);
  //if (fChain->GetBranch("svtxmcorr")) fChain->SetBranchAddress("svtxmcorr", svtxmcorr, &b_svtxmcorr);
  //if (fChain->GetBranch("nIPtrk")) fChain->SetBranchAddress("nIPtrk", nIPtrk, &b_nIPtrk);
  //if (fChain->GetBranch("nselIPtrk")) fChain->SetBranchAddress("nselIPtrk", nselIPtrk, &b_nselIPtrk);
  if (fChain->GetBranch("mue")) fChain->SetBranchAddress("mue", mue, &b_mue);
  if (fChain->GetBranch("mupt")) fChain->SetBranchAddress("mupt", mupt, &b_mupt);
  if (fChain->GetBranch("mueta")) fChain->SetBranchAddress("mueta", mueta, &b_mueta);
  if (fChain->GetBranch("muphi")) fChain->SetBranchAddress("muphi", muphi, &b_muphi);
  if (fChain->GetBranch("mudr")) fChain->SetBranchAddress("mudr", mudr, &b_mudr);
  if (fChain->GetBranch("muptrel")) fChain->SetBranchAddress("muptrel", muptrel, &b_muptrel);
  if (fChain->GetBranch("muchg")) fChain->SetBranchAddress("muchg", muchg, &b_muchg);

  if (fChain->GetBranch("beamId1")) fChain->SetBranchAddress("beamId1", &beamId1, &b_beamId1);
  if (fChain->GetBranch("beamId2")) fChain->SetBranchAddress("beamId2", &beamId2, &b_beamId2);
  if (fChain->GetBranch("pthat")) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
  if (fChain->GetBranch("refpt")) fChain->SetBranchAddress("refpt", refpt, &b_refpt);
  if (fChain->GetBranch("refeta")) fChain->SetBranchAddress("refeta", refeta, &b_refeta);
  if (fChain->GetBranch("refy")) fChain->SetBranchAddress("refy", refy, &b_refy);
  if (fChain->GetBranch("refphi")) fChain->SetBranchAddress("refphi", refphi, &b_refphi);
  if (fChain->GetBranch("refm")) fChain->SetBranchAddress("refm", refm, &b_refm);
  if (fChain->GetBranch("refarea")) fChain->SetBranchAddress("refarea", refarea, &b_refarea);
  if (fChain->GetBranch("reftau1")) fChain->SetBranchAddress("reftau1", reftau1, &b_reftau1);
  if (fChain->GetBranch("reftau2")) fChain->SetBranchAddress("reftau2", reftau2, &b_reftau2);
  if (fChain->GetBranch("reftau3")) fChain->SetBranchAddress("reftau3", reftau3, &b_reftau3);
  if (fChain->GetBranch("refdphijt")) fChain->SetBranchAddress("refdphijt", refdphijt, &b_refdphijt);
  if (fChain->GetBranch("refdrjt")) fChain->SetBranchAddress("refdrjt", refdrjt, &b_refdrjt);
  if (fChain->GetBranch("refparton_pt")) fChain->SetBranchAddress("refparton_pt", refparton_pt, &b_refparton_pt);
  if (fChain->GetBranch("refparton_flavor")) fChain->SetBranchAddress("refparton_flavor", refparton_flavor, &b_refparton_flavor);
  if (fChain->GetBranch("refparton_flavorForB")) fChain->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
  if (fChain->GetBranch("genChargedSum")) fChain->SetBranchAddress("genChargedSum", genChargedSum, &b_genChargedSum);
  if (fChain->GetBranch("genHardSum")) fChain->SetBranchAddress("genHardSum", genHardSum, &b_genHardSum);
  if (fChain->GetBranch("signalChargedSum")) fChain->SetBranchAddress("signalChargedSum", signalChargedSum, &b_signalChargedSum);
  if (fChain->GetBranch("signalHardSum")) fChain->SetBranchAddress("signalHardSum", signalHardSum, &b_signalHardSum);
  if (fChain->GetBranch("subid")) fChain->SetBranchAddress("subid", subid, &b_subid);
  if (fChain->GetBranch("ngen")) fChain->SetBranchAddress("ngen", &ngen, &b_ngen);
  if (fChain->GetBranch("genmatchindex")) fChain->SetBranchAddress("genmatchindex", genmatchindex, &b_genmatchindex);
  if (fChain->GetBranch("genpt")) fChain->SetBranchAddress("genpt", genpt, &b_genpt);
  if (fChain->GetBranch("geneta")) fChain->SetBranchAddress("geneta", geneta, &b_geneta);
  if (fChain->GetBranch("geny")) fChain->SetBranchAddress("geny", geny, &b_geny);
  if (fChain->GetBranch("gentau1")) fChain->SetBranchAddress("gentau1", gentau1, &b_gentau1);
  if (fChain->GetBranch("gentau2")) fChain->SetBranchAddress("gentau2", gentau2, &b_gentau2);
  if (fChain->GetBranch("gentau3")) fChain->SetBranchAddress("gentau3", gentau3, &b_gentau3);
  if (fChain->GetBranch("genphi")) fChain->SetBranchAddress("genphi", genphi, &b_genphi);
  if (fChain->GetBranch("genm")) fChain->SetBranchAddress("genm", genm, &b_genm);
  if (fChain->GetBranch("gendphijt")) fChain->SetBranchAddress("gendphijt", gendphijt, &b_gendphijt);
  if (fChain->GetBranch("gendrjt")) fChain->SetBranchAddress("gendrjt", gendrjt, &b_gendrjt);
  if (fChain->GetBranch("gensubid")) fChain->SetBranchAddress("gensubid", gensubid, &b_gensubid);

  //skimanalysis tree
  if (fChain->GetBranch("Onia2MuMuPAT")) fChain->SetBranchAddress("Onia2MuMuPAT", &Onia2MuMuPAT, &b_Onia2MuMuPAT);
  if (fChain->GetBranch("ana_step")) fChain->SetBranchAddress("ana_step", &ana_step, &b_ana_step);
  if (fChain->GetBranch("pHBHENoiseFilterResultProducer")) fChain->SetBranchAddress("pHBHENoiseFilterResultProducer", &pHBHENoiseFilterResultProducer, &b_pHBHENoiseFilterResultProducer);
  if (fChain->GetBranch("HBHENoiseFilterResult")) fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
  if (fChain->GetBranch("HBHENoiseFilterResultRun1")) fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
  if (fChain->GetBranch("HBHENoiseFilterResultRun2Loose")) fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
  if (fChain->GetBranch("HBHENoiseFilterResultRun2Tight")) fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
  if (fChain->GetBranch("HBHEIsoNoiseFilterResult")) fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);
  if (fChain->GetBranch("pPAprimaryVertexFilter")) fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
  if (fChain->GetBranch("pBeamScrapingFilter")) fChain->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter, &b_pBeamScrapingFilter);
  if (fChain->GetBranch("pVertexFilterCutG")) fChain->SetBranchAddress("pVertexFilterCutG", &pVertexFilterCutG, &b_pVertexFilterCutG);
  if (fChain->GetBranch("pVertexFilterCutGloose")) fChain->SetBranchAddress("pVertexFilterCutGloose", &pVertexFilterCutGloose, &b_pVertexFilterCutGloose);
  if (fChain->GetBranch("pVertexFilterCutGplus")) fChain->SetBranchAddress("pVertexFilterCutGtight", &pVertexFilterCutGtight, &b_pVertexFilterCutGtight);
  if (fChain->GetBranch("pVertexFilterCutGplus")) fChain->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus, &b_pVertexFilterCutGplus);
  if (fChain->GetBranch("pVertexFilterCutE")) fChain->SetBranchAddress("pVertexFilterCutE", &pVertexFilterCutE, &b_pVertexFilterCutE);
  if (fChain->GetBranch("pVertexFilterCutEandG")) fChain->SetBranchAddress("pVertexFilterCutEandG", &pVertexFilterCutEandG, &b_pVertexFilterCutEandG);
}
#endif // #ifndef initTree_C
