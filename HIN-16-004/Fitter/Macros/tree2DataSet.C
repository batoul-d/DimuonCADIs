// -*- C++ -*-
//
// Package:    Fitter
//
/*
 Description: TTree to RooDataSet converter.
 Implementation:
 This program creates two RooDataSets (opposite-sign and same-sign dimuons) from an Onia Tree.
 */
// Original Author:  Andre Stahl,
//         Created:  Feb 26 19:08 CET 2016
//
//

#include "Utilities/initTree.C"
#include "Utilities/EVENTUTILS.h"
#include "Utilities/initClasses.h"
#include "TEfficiency.h"
#include "TObjArray.h"

map<int, double>   fCentMap; // map for centrality-Ncoll mapping
double             fCentBinning[200];
int                fCentBins;
TObjArray*         fcorrArray = NULL; // Array with the 2D correction for weighting
TH2F*              corrHist = NULL;
double drmin = 0.5;
int matchGR= 0;
string  findMyTree(string FileName);
string  findJetTree(string FileName);
string  findSkimTree(string FileName);
bool    getTChain(TChain* fChain, TChain* jChain, TChain* sChain, vector<string> FileNames);
void    iniBranch(TChain* fChain,bool isMC=false);
bool    checkDS(RooDataSet* DS, string DSName);
double  deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon);
bool    isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR=0.03);
double  getNColl(int centr, bool isPP);
double  getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP);
bool    readCorrection(const char* file);
void    setCentralityMap(const char* file);
float   jecCorr(double jtPt, double rawPt, double jpsiPt);

bool tree2DataSet(RooWorkspace& Workspace, vector<string> InputFileNames, string DSName, string OutputFileName, bool UpdateDS=false)
{
  RooDataSet* dataOS = NULL; RooDataSet* dataSS = NULL; RooDataSet* dataOSNoBkg = NULL;
  bool isMC = false;
  if (DSName.find("MC")!=std::string::npos) isMC =true;
  
  bool isPbPb = false;
  if (DSName.find("PbPb")!=std::string::npos) isPbPb =true;
  int triggerIndex_PP   = PP::HLT_HIL1DoubleMu0_v1;
  int triggerIndex_PbPb = HI::HLT_HIL1DoubleMu0_v1;
  int CentFactor = 1;
  
  bool usePeriPD = false;
  if (InputFileNames[0].find("HIOniaPeripheral30100")!=std::string::npos) {
    cout << "[INFO] Working with Peripheral PbPb PD" << endl;
    usePeriPD = true;
    triggerIndex_PbPb = HI::HLT_HIL1DoubleMu0_2HF_Cent30100_v1;
  }
  
  bool applyWeight = false;
  //if (isMC && isPbPb) applyWeight = true;
  if (isMC && DSName.find("JPSIPR")!=std::string::npos) applyWeight = true;  

  bool isPureSDataset = false;
  if (OutputFileName.find("_PureS")!=std::string::npos) isPureSDataset = true;

  bool applyWeight_Corr = false;
  if ( (OutputFileName.find("_AccEff")!=std::string::npos) || (OutputFileName.find("_lJpsiEff")!=std::string::npos) ) applyWeight_Corr = true;
  //if(applyWeight == true) applyWeight_Corr = false;
  TString corrName = "";
  TString corrFileName = "";
  if (OutputFileName.find("_AccEff")!=std::string::npos)
  {
    corrFileName = "correction_AccEff.root";
    corrName = "AccEff";
  }
  else if (OutputFileName.find("_lJpsiEff")!=std::string::npos)
  {
    corrFileName = "correction_lJpsiEff.root";
    corrName = "lJpsiEff";
  }

  bool applyJEC = false;
  if (OutputFileName.find("_JEC")!=std::string::npos) applyJEC = true;

  bool createDS = ( gSystem->AccessPathName(OutputFileName.c_str()) || UpdateDS );

  if ( !gSystem->AccessPathName(OutputFileName.c_str()) ) {
    cout << "[INFO] Loading RooDataSet from " << OutputFileName << endl;
    
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"READ");
    if (isMC && isPureSDataset) {
      if (applyWeight_Corr) {
	dataOSNoBkg = (RooDataSet*)DBFile->Get(Form("dOS_%s_NoBkg_%s%s", DSName.c_str(),corrName.Data(), (applyJEC?"_JEC":"")));
	if (checkDS(dataOSNoBkg, DSName)==false) { createDS = true; }
      }
      else {
	dataOSNoBkg = (RooDataSet*)DBFile->Get(Form("dOS_%s_NoBkg%s", DSName.c_str(), (applyJEC?"_JEC":"")));
	if (checkDS(dataOSNoBkg, DSName)==false) { createDS = true; }
      }  
    } 
    else if (applyWeight_Corr) {
      dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s_%s%s", DSName.c_str(),corrName.Data(), (applyJEC?"_JEC":"")));
      if (checkDS(dataOS, DSName)==false) { createDS = true; }
    }
    else {
      dataOS = (RooDataSet*)DBFile->Get(Form("dOS_%s%s", DSName.c_str(),(applyJEC?"_JEC":"")));
      if (checkDS(dataOS, DSName)==false) { createDS = true; }
      dataSS = (RooDataSet*)DBFile->Get(Form("dSS_%s%s", DSName.c_str(),(applyJEC?"_JEC":"")));
      if (checkDS(dataSS, DSName)==false) { createDS = true; }
    }
    DBFile->Close(); delete DBFile;
  }
  //cout<< "FileName: "<< OutputFileName << " DSName: " << DSName.c_str() << " CorrName: " << corrName.Data() << " or " <<corrName << endl;
  if (createDS) {
    cout << "[INFO] Creating " << (isPureSDataset ? "pure signal " : "") << "RooDataSet for " << DSName << endl;
    TreeName = findMyTree(InputFileNames[0]); if(TreeName==""){return false;}
    jetTreeName = findJetTree(InputFileNames[0]); //if(jetTreeName==""){return false;}
    skimTreeName = findSkimTree(InputFileNames[0]);
    TChain* theTree = new TChain(TreeName.c_str(),"");
    TChain* jetTree = new TChain(jetTreeName.c_str(),"");
    TChain* skimTree = new TChain(skimTreeName.c_str(),"");
    if(!getTChain(theTree, jetTree, skimTree, InputFileNames)){ return false; }     // Import files to TChain
    initTree(theTree);                         // Initialize the Tree
    iniBranch(theTree,isMC);                   // Initialize the Branches

    RooRealVar* mass    = new RooRealVar("invMass","#mu#mu mass", 1.0, 6.0, "GeV/c^{2}");
    RooRealVar* zed     = new RooRealVar("zed", "z_{J/#psi}", -0.1, 2);
    RooRealVar* ctau    = new RooRealVar("ctau","c_{#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ctauN    = new RooRealVar("ctauN","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauTrue = new RooRealVar("ctauTrue","c_{#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ctauNRes = new RooRealVar("ctauNRes","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauRes = new RooRealVar("ctauRes","c_{#tau}", -100000.0, 100000.0, "");
    RooRealVar* ctauErr = new RooRealVar("ctauErr","#sigma_{c#tau}", -100000.0, 100000.0, "mm");
    RooRealVar* ptQQ    = new RooRealVar("pt","#mu#mu p_{T}", -1.0, 10000.0, "GeV/c");
    RooRealVar* rapQQ   = new RooRealVar("rap","#mu#mu y", -2.5, 2.5, "");
    RooRealVar* ptJet    = new RooRealVar("jetpt","Jet p_{T}", -1.0, 10000.0, "GeV/c");
    RooRealVar* rapJet   = new RooRealVar("jetrap","Jet y", -2.5, 2.5, "");
    RooRealVar* cent    = new RooRealVar("cent","centrality", -1.0, 1000.0, "");
    RooRealVar* weight  = new RooRealVar("weight","MC weight", 0.0, 10000000.0, "");
    RooRealVar* weightCorr   = new RooRealVar("weightCorr","Data correction weight", 0.0, 10000000.0, "");
    RooArgSet*  cols    = NULL;
    
    if (applyWeight && !applyWeight_Corr)
    {
      setCentralityMap(Form("%s/Input/CentralityMap_PbPb2015.txt",gSystem->ExpandPathName(gSystem->pwd())));
      if (isMC) {
        cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weight);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
	cols->add(*ptJet);
	cols->add(*rapJet);
      } else {
        cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent, *weight);
        cols->add(*ctauN);
	cols->add(*ptJet);
	cols->add(*rapJet);
      }
      dataOS = new RooDataSet(Form("dOS_%s%s", DSName.c_str(), (applyJEC?"_JEC":"")), "dOS", *cols, WeightVar(*weight), StoreAsymError(*mass));
      dataSS = new RooDataSet(Form("dSS_%s%s", DSName.c_str(), (applyJEC?"_JEC":"")), "dSS", *cols, WeightVar(*weight), StoreAsymError(*mass));
      if (isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg%s", DSName.c_str(), (applyJEC?"_JEC":"")), "dOSNoBkg", *cols, WeightVar(*weight), StoreAsymError(*mass));
    }
    else if (applyWeight_Corr)
    {
      if (isMC) {
        cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent, *weightCorr);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
	cols->add(*ptJet);
	cols->add(*rapJet);
      } else {
        cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ptQQ, *rapQQ,*cent, *weightCorr);
        cols->add(*ctauN);
	cols->add(*ptJet);
	cols->add(*rapJet);
      }
      if (!readCorrection(Form("%s/Input/%s",gSystem->ExpandPathName(gSystem->pwd()),corrFileName.Data()))){ return false; }
      dataOS = new RooDataSet(Form("dOS_%s_%s%s", DSName.c_str(),corrName.Data(), (applyJEC?"_JEC":"")), "dOS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      if (isMC && isPureSDataset)
	dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg_%s%s", DSName.c_str(),corrName.Data(),(applyJEC?"_JEC":"")), "dOSNoBkg", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      //      dataSS = new RooDataSet(Form("dSS_%s", DSName.c_str()), "dSS", *cols, WeightVar(*weightCorr), StoreAsymError(*mass));
      cout<<"[INFO] "<<corrName<<" applied!"<<endl;
    }
    else
    {
      if (isMC) {
        cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ctauTrue, *ptQQ, *rapQQ, *cent);
        cols->add(*ctauNRes);
        cols->add(*ctauRes);
	cols->add(*ptJet);
	cols->add(*rapJet);
      } else {
        cols = new RooArgSet(*mass, *zed, *ctau, *ctauErr, *ptQQ, *rapQQ, *cent);
        cols->add(*ctauN);
	cols->add(*ptJet);
	cols->add(*rapJet);
      }  
      dataOS = new RooDataSet(Form("dOS_%s%s", DSName.c_str(),(applyJEC?"_JEC":"")), "dOS", *cols, StoreAsymError(*mass));
      dataSS = new RooDataSet(Form("dSS_%s%s", DSName.c_str(), (applyJEC?"_JEC":"")), "dSS", *cols, StoreAsymError(*mass));
      if (isMC && isPureSDataset) dataOSNoBkg = new RooDataSet(Form("dOS_%s_NoBkg%s", DSName.c_str(),(applyJEC?"_JEC":"")), "dOSNoBkg", *cols, StoreAsymError(*mass));
    }
    if (applyWeight) cout<<"[INFO] pt weights applied!"<<endl;


    ////////////////////////////////////
    Long64_t nentries = theTree->GetEntries();
    //nentries = 2000000;
    float normF = 0.;
    if (isMC && isPbPb)
    {
      cout << "[INFO] Computing sum of weights for " << nentries << " nentries" << endl;
      
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        if (jentry%1000000==0) cout << "[INFO] " << jentry << "/" << nentries << endl;
        
        if (theTree->LoadTree(jentry)<0) break;
        if (theTree->GetTreeNumber()!=fCurrent) {
          fCurrent = theTree->GetTreeNumber();
          cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
        }
        
        theTree->GetEntry(jentry);
       // normF += theTree->GetWeight()*getNColl(Centrality,!isPbPb);
        normF += theTree->GetWeight();  
    }
      
      normF = nentries/normF;
    }
    // creating the tree to use in the unfolding
    string fl = InputFileNames[0];
    //cout<<"[INFO] nb of files "<<InputFileNames.size()<<endl;
    if (fl.find("ext")!=std::string::npos && InputFileNames.size()>1) fl = "";
    else if (fl.find("ext")!=std::string::npos && InputFileNames.size()==1) fl = "_ext";
    else if (fl.find("pthat15")!=std::string::npos) fl = "_pthat15";
    else if (fl.find("pthat25")!=std::string::npos) fl = "_pthat25";
    else if (fl.find("pthat35")!=std::string::npos) fl = "_pthat35";
    else if (fl.find("pthat45")!=std::string::npos) fl = "_pthat45";

    gSystem->mkdir("TreesForUnfolding");
    string trUnfFileName = Form("TreesForUnfolding/tree_%s%s%s%s%s.root", DSName.c_str(), (isPureSDataset?"_NoBkg":""), (applyWeight_Corr?Form("_%s",corrName.Data()):""), (applyJEC?"_JEC":""), (applyWeight? fl.c_str():""));

    TFile * trUnfFile = new TFile (trUnfFileName.c_str(),"RECREATE");
    //trUnfFile->cd();
    TTree* trUnf = new TTree ("trUnf","tree used for the unfolding");
    Int_t evtNb; Float_t jp_pt; Float_t jp_rap; Float_t jp_eta; Float_t jp_mass; Float_t jp_phi; Float_t jp_l; 
    Float_t jp_gen_pt; Float_t jp_gen_rap; Float_t jp_gen_eta; Float_t jp_gen_phi; 
    Float_t jt_pt; Float_t jt_rap; Float_t jt_eta; Float_t jt_phi; Float_t jt_CHF; Float_t jt_NHF; Float_t jt_CEF; Float_t jt_NEF; Float_t jt_MUF; Float_t jt_CHM; Float_t jt_NHM; Float_t jt_CEM; Float_t jt_NEM; Float_t jt_MUM; 
    Float_t jt_ref_pt; Float_t jt_ref_rap; Float_t jt_ref_eta; Float_t jt_ref_phi; Float_t jt_gen_chargedSum; Float_t jt_gen_hardSum; Float_t jt_signal_chargedSum; Float_t jt_signal_hardSum;
    Float_t z; Float_t gen_z; 
    Float_t corr_AccEff; Float_t pt_hat; Float_t corr_ptw;

    cout<< "[INFO] Creating the tree to use in the unfolding" << endl;
    trUnf->Branch("evtNb", &evtNb, "evtNb/I");
    trUnf->Branch("jp_pt", &jp_pt, "jp_pt/F");
    trUnf->Branch("jp_rap", &jp_rap, "jp_rap/F");
    trUnf->Branch("jp_eta", &jp_eta, "jp_eta/F");
    trUnf->Branch("jp_phi", &jp_phi, "jp_phi/F");
    trUnf->Branch("jp_mass", &jp_mass, "jp_mass/F");
    trUnf->Branch("jp_l", &jp_l, "jp_l/F");
    trUnf->Branch("jt_pt", &jt_pt, "jt_pt/F");
    trUnf->Branch("jt_rap", &jt_rap, "jt_rap/F");
    trUnf->Branch("jt_eta", &jt_eta, "jt_eta/F");
    trUnf->Branch("jt_phi", &jt_phi, "jt_phi/F");
    trUnf->Branch("jt_CHF", &jt_CHF, "jt_CHF/F");
    trUnf->Branch("jt_NHF", &jt_NHF, "jt_NHF/F");
    trUnf->Branch("jt_CEF", &jt_CEF, "jt_CEF/F");
    trUnf->Branch("jt_NEF", &jt_NEF, "jt_NEF/F");
    trUnf->Branch("jt_MUF", &jt_MUF, "jt_MUF/F");
    trUnf->Branch("jt_CHM", &jt_CHM, "jt_CHM/F");
    trUnf->Branch("jt_NHM", &jt_NHM, "jt_NHM/F");
    trUnf->Branch("jt_CEM", &jt_CEM, "jt_CEM/F");
    trUnf->Branch("jt_NEM", &jt_NEM, "jt_NEM/F");
    trUnf->Branch("jt_MUM", &jt_MUM, "jt_MUM/F");
    trUnf->Branch("z", &z, "z/F");
    trUnf->Branch("corr_AccEff", &corr_AccEff, "corr_AccEff/F");
    if (isMC) {
      trUnf->Branch("corr_ptw", &corr_ptw, "corr_ptw/F");
      trUnf->Branch("pt_hat", &pt_hat, "pt_hat/F");
      trUnf->Branch("jp_gen_pt", &jp_gen_pt, "jp_gen_pt/F");
      trUnf->Branch("jp_gen_rap", &jp_gen_rap, "jp_gen_rap/F");
      trUnf->Branch("jp_gen_eta", &jp_gen_eta, "jp_gen_eta/F");
      trUnf->Branch("jp_gen_phi", &jp_gen_phi, "jp_gen_phi/F");
      trUnf->Branch("jt_ref_pt", &jt_ref_pt, "jt_ref_pt/F");
      trUnf->Branch("jt_ref_rap", &jt_ref_rap, "jt_ref_rap/F");
      trUnf->Branch("jt_ref_eta", &jt_ref_eta, "jt_ref_eta/F");
      trUnf->Branch("jt_ref_phi", &jt_ref_phi, "jt_ref_phi/F");
      trUnf->Branch("jt_gen_chargedSum", &jt_gen_chargedSum, "jt_gen_chargedSum/F");
      trUnf->Branch("jt_gen_hardSum", &jt_gen_hardSum, "jt_gen_hardSum/F");
      trUnf->Branch("jt_signal_chargedSum", &jt_signal_chargedSum, "jt_signal_chargedSum/F");
      trUnf->Branch("jt_signal_hardSum", &jt_gen_hardSum, "jt_signal_hardSum/F");
      trUnf->Branch("gen_z", &gen_z, "gen_z/F");
    }

    cout << "[INFO] Starting to process " << nentries << " nentries" << endl;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if (jentry%1000000==0) cout << "[INFO] " << jentry << "/" << nentries <<endl;//", runNb = " <<runNb<< endl;
      
      pPAprimaryVertexFilter=1;
      pBeamScrapingFilter=1;

      if (theTree->LoadTree(jentry)<0) break;
      if (theTree->GetTreeNumber()!=fCurrent) {
        fCurrent = theTree->GetTreeNumber();
        cout << "[INFO] Processing Root File: " << InputFileNames[fCurrent] << endl;
      }
      Reco_QQ_4mom->Clear();
      Reco_QQ_mumi_4mom->Clear();
      Reco_QQ_mupl_4mom->Clear();
      if (isMC) {
        Gen_QQ_mumi_4mom->Clear();
        Gen_QQ_mupl_4mom->Clear();
      }
      theTree->GetEntry(jentry);

      for (int iQQ=0; iQQ<Reco_QQ_size; iQQ++) {
	drmin= 0.5;
	zed->setVal(100);
	z=100;
	ptJet->setVal(1);
	rapJet->setVal(100);

        TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(iQQ);
        mass->setVal(RecoQQ4mom->M());
        if (theTree->GetBranch("Reco_QQ_ctau3D")) { ctau->setVal(Reco_QQ_ctau3D[iQQ]); }
        else if (theTree->GetBranch("Reco_QQ_ctau")) { ctau->setVal(Reco_QQ_ctau[iQQ]); }
        else { cout << "[ERROR] No ctau information found in the Onia Tree" << endl; }
        if (theTree->GetBranch("Reco_QQ_ctauErr3D")) { ctauErr->setVal(Reco_QQ_ctauErr3D[iQQ]); }
        else if (theTree->GetBranch("Reco_QQ_ctauErr")) { ctauErr->setVal(Reco_QQ_ctauErr[iQQ]); }
        else { cout << "[ERROR] No ctauErr information found in the Onia Tree" << endl; }
        
        ctauN->setVal(ctau->getVal()/ctauErr->getVal());
        ptQQ->setVal(RecoQQ4mom->Pt());
        rapQQ->setVal(RecoQQ4mom->Rapidity());
        cent->setVal(Centrality*CentFactor);
	jp_pt = RecoQQ4mom->Pt();
	jp_rap = RecoQQ4mom->Rapidity();
	jp_eta = RecoQQ4mom->Eta();
	jp_phi = RecoQQ4mom->Phi();
	jp_mass = RecoQQ4mom->M();
	jp_l = Reco_QQ_ctau[iQQ];
	pt_hat = pthat;
	corr_ptw = 1;

	for (Long64_t ijet=0; ijet<nref; ijet++)
	{
	    TLorentzVector v_jet;
	    v_jet.SetPtEtaPhiM(jtpt[ijet], jteta[ijet], jtphi[ijet], jtm[ijet]);
	    if (RecoQQ4mom->DeltaR (v_jet)<=drmin)
	      {
		drmin = RecoQQ4mom->DeltaR (v_jet);
		if (applyJEC){
		  zed->setVal(RecoQQ4mom->Pt()/jecCorr(jtpt[ijet], rawpt[ijet], RecoQQ4mom->Pt()));
		  if (zed->getVal() > 1 && zed->getVal() <= 1.000001) zed->setVal (0.9999999);
		  ptJet->setVal(jecCorr(jtpt[ijet], rawpt[ijet], RecoQQ4mom->Pt()));
		  jt_pt = jecCorr(jtpt[ijet], rawpt[ijet], RecoQQ4mom->Pt());
		  z = jp_pt/jt_pt;
		  if (z > 1 && z <= 1.000001) z = 0.9999999;
		}
		else {
		  zed->setVal(RecoQQ4mom->Pt()/jtpt[ijet]);
		  if (zed->getVal() > 1 && zed->getVal() <= 1.000001) zed->setVal(0.9999999);
		  ptJet->setVal(jtpt[ijet]);
		  jt_pt = jtpt[ijet];
		  z = RecoQQ4mom->Pt()/jt_pt;
		  if (z > 1 && z <= 1.000001) z = 0.9999999;
		}
		rapJet->setVal(jty[ijet]);
		jt_rap = jty[ijet];
		jt_eta = jteta[ijet];
		jt_phi = jtphi[ijet];
		jt_CHF = jtPfCHF[ijet];
		jt_NHF = jtPfNHF[ijet];
		jt_CEF = jtPfCEF[ijet];
		jt_NEF = jtPfNEF[ijet];
		jt_MUF = jtPfMUF[ijet];
		jt_CHM = jtPfCHM[ijet];
		jt_NHM = jtPfNHM[ijet];
		jt_CEM = jtPfCEM[ijet];
		jt_NEM = jtPfNEM[ijet];
		jt_MUM = jtPfMUM[ijet];

		if (isMC) {
		  jt_ref_pt = refpt[ijet];
		  jt_ref_rap = refy[ijet];
		  jt_ref_eta = refeta[ijet];
		  jt_ref_phi = refphi[ijet];
		  jt_gen_chargedSum = genChargedSum[ijet];
		  jt_gen_hardSum = genHardSum[ijet];
		  jt_signal_chargedSum = signalChargedSum[ijet];
		  jt_signal_hardSum = signalHardSum[ijet];
		}
	      }
	  }

        if (isMC) {
          if (theTree->GetBranch("Reco_QQ_ctauTrue3D")) { ctauTrue->setVal(Reco_QQ_ctauTrue3D[iQQ]); }
          else if (theTree->GetBranch("Reco_QQ_ctauTrue")) { ctauTrue->setVal(Reco_QQ_ctauTrue[iQQ]); }
          else { cout << "[ERROR] No ctauTrue information found in the Onia Tree" << endl; }
          ctauNRes->setVal( (ctau->getValV() - ctauTrue->getValV())/(ctauErr->getValV()) );
          ctauRes->setVal( (ctau->getValV() - ctauTrue->getValV()) );
        }

        if (applyWeight && !applyWeight_Corr){
	  double w = theTree->GetWeight();
	  if (isMC && isPbPb) w = w*normF;//*getNColl(Centrality,!isPbPb)*normF;
	  weight->setVal(w);
        }
        else if (applyWeight_Corr)
	  { double Corr = 1;
	    Corr = getCorr(RecoQQ4mom->Rapidity(),RecoQQ4mom->Pt(),RecoQQ4mom->M(),!isPbPb);
	    double wCorr = 1/Corr;
	    //// add pt weights for prompt MC
	    if (applyWeight)
	      {
		if (pthat >= 15 && pthat < 25)  corr_ptw = 0.168329;
		else if (pthat >= 25 && pthat < 35) corr_ptw = 0.0238285;
		else if (pthat >= 35 && pthat < 45) corr_ptw = 0.00517143;
		else if (pthat >= 45) corr_ptw = 0.00150836;
	      }
	    corr_AccEff = wCorr;
	    weightCorr->setVal(wCorr*corr_ptw);
	  }
        if (
            ( RecoQQ::areMuonsInAcceptance2015(iQQ) ) &&  // 2015 Global Muon Acceptance Cuts
            ( RecoQQ::passQualityCuts2015(iQQ)) &&  // 2015 Soft Global Muon Quality Cuts
            ( isPbPb ? (RecoQQ::isTriggerMatch(iQQ,triggerIndex_PbPb) || (usePeriPD ? RecoQQ::isTriggerMatch(iQQ,HI::HLT_HIL1DoubleMu0_2HF0_Cent30100_v1) : (RecoQQ::isTriggerMatch(iQQ,HI::HLT_HIL1DoubleMu0_2HF_v1) || RecoQQ::isTriggerMatch(iQQ,HI::HLT_HIL1DoubleMu0_2HF0_v1)))) :
              RecoQQ::isTriggerMatch(iQQ, triggerIndex_PP) )     // if PbPb && !periPD then (HLT_HIL1DoubleMu0_v1 || HLT_HIL1DoubleMu0_2HF_v1)
            )
	  {
	    if (pPAprimaryVertexFilter)
	    {
	      if (pBeamScrapingFilter)
		  { 
		    if (Reco_QQ_sign[iQQ]==0) { // Opposite-Sign dimuons
		      if (isMC && isPureSDataset && isMatchedRecoDiMuon(iQQ)) {
			if (applyWeight_Corr)
			  dataOSNoBkg->add(*cols, weightCorr->getVal()); //Signal-only dimuons
			else 
			  dataOSNoBkg->add(*cols, (applyWeight ? weight->getVal() : 1.0)); // Signal-only dimuons
			TLorentzVector *GenQQ4mom = (TLorentzVector*) Gen_QQ_4mom->At(matchGR);
			jp_gen_pt = GenQQ4mom->Pt();
			jp_gen_rap = GenQQ4mom->Rapidity();
			jp_gen_eta = GenQQ4mom->Eta();
			jp_gen_phi = GenQQ4mom->Phi();
			gen_z = jp_gen_pt/jt_ref_pt;
			if (gen_z > 1 && gen_z <= 1.000001) gen_z = 0.9999999;
		      }
		      else if (isMC && isPureSDataset && !isMatchedRecoDiMuon(iQQ))
			gen_z = -1;
		      
		      else if (applyWeight_Corr) dataOS->add(*cols,weightCorr->getVal()); //Signal and background dimuons
		      else dataOS->add(*cols, ( applyWeight ? weight->getVal() : 1.0)); //Signal and background dimuons
		      evtNb = jentry;
		      if (isMC && isPureSDataset && gen_z >= 0 && z < 100)
			trUnf->Fill();
		      
		      else if (!isPureSDataset && z < 100)
			trUnf->Fill();
		    }
		    else { // Like-Sign dimuons
		      if (!isPureSDataset && !applyWeight_Corr ) dataSS->add(*cols, ( applyWeight  ? weight->getVal() : 1.0));
		    }
		  }
	      }
	  }
      }
    }

    // Close the TChain and all its pointers
    delete Reco_QQ_4mom; delete Reco_QQ_mumi_4mom; delete Reco_QQ_mupl_4mom; delete Gen_QQ_mumi_4mom; delete Gen_QQ_mupl_4mom;
    theTree->Reset(); delete theTree;
    // Save the tree for the unfolding
    cout<<"[INFO] Saving the tree for the unfolding" << endl;
    trUnfFile->cd();
    trUnf->Write("treeForUnfolding");
    trUnfFile->Close(); delete trUnfFile;
    // Save all the datasets
    TFile *DBFile = TFile::Open(OutputFileName.c_str(),"RECREATE");
    DBFile->cd();
    if (isMC && isPureSDataset &&  applyWeight_Corr) {
      dataOSNoBkg->Write(Form("dOS_%s_NoBkg_%s%s", DSName.c_str(),corrName.Data(),(applyJEC?"_JEC":"")));
    }
    else if (isMC && isPureSDataset && !applyWeight_Corr) {
      dataOSNoBkg->Write(Form("dOS_%s_NoBkg%s", DSName.c_str(),(applyJEC?"_JEC":"")));
    }
    else if (!isPureSDataset && applyWeight_Corr) {
      dataOS->Write(Form("dOS_%s_%s%s", DSName.c_str(),corrName.Data(),(applyJEC?"_JEC":"")));
    }
    else 
    {
      dataOS->Write(Form("dOS_%s%s", DSName.c_str(),(applyJEC?"_JEC":"")));
      dataSS->Write(Form("dSS_%s%s", DSName.c_str(),(applyJEC?"_JEC":"")));
    }
    DBFile->Write(); DBFile->Close(); delete DBFile;
  }
  
  // Import datasets to workspace
  if (isMC && isPureSDataset)
  {
    if (!dataOSNoBkg) { cout << "[ERROR] " << DSName << "_NoBkg was not found" << endl; return false; }
    Workspace.import(*dataOSNoBkg);
  }
  else if (applyWeight_Corr)
  {
    if(!dataOS) { cout << "[ERROR] " << DSName << "_" << corrName.Data() << " was not found" << endl; return false; }
    Workspace.import(*dataOS);
  }
  else
  {
    if(!dataOS || !dataSS) { cout << "[ERROR] " << DSName << " was not found" << endl; return false; }
    Workspace.import(*dataOS);
    Workspace.import(*dataSS);
  }
  
  // delete the local datasets
  delete dataSS; delete dataOS; delete dataOSNoBkg;
  
  // delete the correction array
  if (fcorrArray) delete fcorrArray; delete corrHist;
  return true;
};

string findMyTree(string FileName)
{
  TFile *f = TFile::Open(FileName.c_str(), "READ");
  string name = "";
  if(f->GetListOfKeys()->Contains("hionia")) name = "hionia/myTree";
  else if(f->GetListOfKeys()->Contains("myTree")) name = "myTree";
  else { cout << "[ERROR] myTree was not found in: " << FileName << endl;}
  f->Close(); delete f;
  return name;
};

string  findJetTree(string FileName)
{
  TFile *f = TFile::Open(FileName.c_str(), "READ");
  string name = "";
  if(f->GetListOfKeys()->Contains("ak4PFJetAnalyzer")) name = "ak4PFJetAnalyzer/t";
  else if(f->GetListOfKeys()->Contains("t")) name = "t";
  else { cout << "[ERROR] t was not found in: " << FileName << endl; }
  f->Close(); delete f;
  return name;
};

string  findSkimTree(string FileName)
{
  TFile *f = TFile::Open(FileName.c_str(), "READ");
  string name = "";
  if(f->GetListOfKeys()->Contains("skimanalysis")) name = "skimanalysis/HltTree";
  else if(f->GetListOfKeys()->Contains("HltTree")) name = "HltTree";
  else { cout << "[ERROR] HltTree was not found in: " << FileName << endl; }
  f->Close(); delete f;
  return name;
};

bool getTChain(TChain *fChain, TChain *jChain, TChain *sChain, vector<string> FileNames)
{
  cout << "[INFO] Extrating TTree " << TreeName.c_str() << endl;
  for (vector<string>::iterator FileName = FileNames.begin() ; FileName != FileNames.end(); ++FileName){
    cout << "[INFO] Adding TFile " << FileName->c_str() << endl;
    fChain->Add(Form("%s/%s", FileName->c_str(),  TreeName.c_str()));
    jChain->Add(Form("%s/%s", FileName->c_str(),  jetTreeName.c_str()));
    sChain->Add(Form("%s/%s", FileName->c_str(),  skimTreeName.c_str()));
  }
  if (jChain)
    fChain->AddFriend(jChain);
  if (sChain)
    fChain->AddFriend(sChain);

  if (!fChain /*|| !jChain*/) { cout << "[ERROR] fChain was not created, some input files are missing" << endl; return false; }
  return true;
};

void iniBranch(TChain* fChain, bool isMC)
{
  cout << "[INFO] Initializing Branches of " << TreeName.c_str() << endl;
  if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->GetBranch("Reco_QQ_4mom")->SetAutoDelete(false);      }
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(false); }
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(false); }
  if (isMC) 
    {
      if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->GetBranch("Gen_QQ_mupl_4mom")->SetAutoDelete(false); }
      if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->GetBranch("Gen_QQ_mumi_4mom")->SetAutoDelete(false); }
    }
  fChain->SetBranchStatus("*",0);
  RecoQQ::iniBranches(fChain);
  if (fChain->GetBranch("runNb"))             { fChain->SetBranchStatus("runNb",1);             }
  if (fChain->GetBranch("Centrality"))        { fChain->SetBranchStatus("Centrality",1);        }
  if (fChain->GetBranch("Reco_QQ_size"))      { fChain->SetBranchStatus("Reco_QQ_size",1);      }
  if (fChain->GetBranch("Reco_QQ_sign"))      { fChain->SetBranchStatus("Reco_QQ_sign",1);      }
  if (fChain->GetBranch("Reco_QQ_4mom"))      { fChain->SetBranchStatus("Reco_QQ_4mom",1);      }
  if (fChain->GetBranch("Reco_QQ_mupl_4mom")) { fChain->SetBranchStatus("Reco_QQ_mupl_4mom",1); }
  if (fChain->GetBranch("Reco_QQ_mumi_4mom")) { fChain->SetBranchStatus("Reco_QQ_mumi_4mom",1); }
  if (fChain->GetBranch("Reco_QQ_ctau3D"))    { fChain->SetBranchStatus("Reco_QQ_ctau3D",1);    }
  if (fChain->GetBranch("Reco_QQ_ctauErr3D")) { fChain->SetBranchStatus("Reco_QQ_ctauErr3D",1); }
  if (fChain->GetBranch("Reco_QQ_ctau"))      { fChain->SetBranchStatus("Reco_QQ_ctau",1);      }
  if (fChain->GetBranch("Reco_QQ_ctauErr"))   { fChain->SetBranchStatus("Reco_QQ_ctauErr",1);   }
  if (fChain->GetBranch("nref"))              { fChain->SetBranchStatus("nref",1);              }
  if (fChain->GetBranch("jtpt"))              { fChain->SetBranchStatus("jtpt",1);              }
  if (fChain->GetBranch("rawpt"))             { fChain->SetBranchStatus("rawpt",1);             }
  if (fChain->GetBranch("jteta"))             { fChain->SetBranchStatus("jteta",1);             }
  if (fChain->GetBranch("jty"))               { fChain->SetBranchStatus("jty",1);               }
  if (fChain->GetBranch("jtphi"))             { fChain->SetBranchStatus("jtphi",1);             }
  if (fChain->GetBranch("jtm"))               { fChain->SetBranchStatus("jtm",1);               }
  if (fChain->GetBranch("jtPfCHF"))           { fChain->SetBranchStatus("jtPfCHF",1);           }
  if (fChain->GetBranch("jtPfNHF"))           { fChain->SetBranchStatus("jtPfNHF",1);           }
  if (fChain->GetBranch("jtPfCEF"))           { fChain->SetBranchStatus("jtPfCEF",1);           }
  if (fChain->GetBranch("jtPfNEF"))           { fChain->SetBranchStatus("jtPfNEF",1);           }
  if (fChain->GetBranch("jtPfMUF"))           { fChain->SetBranchStatus("jtPfMUF",1);           }
  if (fChain->GetBranch("jtPfCHM"))           { fChain->SetBranchStatus("jtPfCHM",1);           }
  if (fChain->GetBranch("jtPfNHM"))           { fChain->SetBranchStatus("jtPfNHM",1);           }
  if (fChain->GetBranch("jtPfCEM"))           { fChain->SetBranchStatus("jtPfCEM",1);           }
  if (fChain->GetBranch("jtPfNEM"))           { fChain->SetBranchStatus("jtPfNEM",1);           }
  if (fChain->GetBranch("jtPfMUM"))           { fChain->SetBranchStatus("jtPfMUM",1);           }
  if (fChain->GetBranch("pPAprimaryVertexFilter")){ fChain->SetBranchStatus("pPAprimaryVertexFilter",1);}
  if (fChain->GetBranch("pBeamScrapingFilter")){ fChain->SetBranchStatus("pBeamScrapingFilter",1);}

  if (isMC)
  {
    if (fChain->GetBranch("pthat"))            { fChain->SetBranchStatus("pthat",1);            }
    if (fChain->GetBranch("Gen_QQ_4mom"))      { fChain->SetBranchStatus("Gen_QQ_4mom",1);      }
    if (fChain->GetBranch("Gen_QQ_size"))      { fChain->SetBranchStatus("Gen_QQ_size",1);      }
    if (fChain->GetBranch("Gen_QQ_mupl_4mom")) { fChain->SetBranchStatus("Gen_QQ_mupl_4mom",1); }
    if (fChain->GetBranch("Gen_QQ_mumi_4mom")) { fChain->SetBranchStatus("Gen_QQ_mumi_4mom",1); }
    if (fChain->GetBranch("Reco_QQ_ctauTrue3D")) { fChain->SetBranchStatus("Reco_QQ_ctauTrue3D",1); }
    if (fChain->GetBranch("Reco_QQ_ctauTrue")) { fChain->SetBranchStatus("Reco_QQ_ctauTrue",1); }
    if (fChain->GetBranch("refpt"))            { fChain->SetBranchStatus("refpt",1);            }
    if (fChain->GetBranch("refeta"))           { fChain->SetBranchStatus("refeta",1);           }
    if (fChain->GetBranch("refy"))             { fChain->SetBranchStatus("refy",1);             }
    if (fChain->GetBranch("refphi"))           { fChain->SetBranchStatus("refphi",1);           }
    if (fChain->GetBranch("refm"))             { fChain->SetBranchStatus("refm",1);             }
    if (fChain->GetBranch("genChargedSum"))    { fChain->SetBranchStatus("genChargedSum",1);    }
    if (fChain->GetBranch("genHardSum"))       { fChain->SetBranchStatus("genHardSum",1);       }
    if (fChain->GetBranch("signalChargedSum")) { fChain->SetBranchStatus("signalChargedSum",1); }
    if (fChain->GetBranch("signalHardSum"))       { fChain->SetBranchStatus("signalHardSum",1);    }
  }
};

bool checkDS(RooDataSet* DS, string DSName)
{
  bool incCent     = (DSName.find("PbPb")!=std::string::npos);
  bool incCtauTrue = (DSName.find("MC")!=std::string::npos);
  const RooArgSet* row = DS->get();
  if (
      (row->find("invMass")!=0) &&
      (row->find("pt")!=0)      &&
      (row->find("ctau")!=0)    &&
      (row->find("ctauErr")!=0) &&
      (row->find("zed")!=0)     &&
      (row->find("jetpt")!=0)   &&
      (incCent     ? row->find("cent")!=0     : true) &&
      (incCtauTrue ? row->find("ctauTrue")!=0 : true) &&
      (incCtauTrue ? row->find("ctauRes")!=0 : true) &&
      (incCtauTrue ? row->find("ctauNRes")!=0 : true) &&
      (incCtauTrue ? true : row->find("ctauN")!=0)
      ) 
    { return true; }
  else 
    { cout << "[WARNING] Original dataset: " << DS->GetName() << " is corrupted, will remake it!" << endl; }

  return false;
};

double deltaR(TLorentzVector* GenMuon, TLorentzVector* RecoMuon)
{
  double dEta = RecoMuon->Eta() - GenMuon->Eta();
  double dPhi = TVector2::Phi_mpi_pi(RecoMuon->Phi() - GenMuon->Phi());
  return ((double) TMath::Sqrt( (dEta*dEta) + (dPhi*dPhi) ) );
};

bool isMatchedRecoDiMuon(int iRecoDiMuon, double maxDeltaR)
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
    if ( (dRpl < maxDeltaR) && (dRmi < maxDeltaR)  ) isMatched = true;
    matchGR = iGenMuon;
    iGenMuon++;
  }
  return isMatched;
};

double getNColl(int centr, bool isPP)
{
  // Returns the corresponding Ncoll value to the "centr" centrality bin
  
  if ( isPP ) return 1.;
  
  int normCent = TMath::Nint(centr/2.);
  
  int lcent = 0;
  int ucent = 0;
  for ( int i = 0 ; i < fCentBins ; i++ )
  {
    ucent = fCentBinning[i];
    if ( (normCent >= lcent) && (normCent < ucent) ) return fCentMap[ucent];
    else lcent = ucent;
  }
  return 1.;
};

void setCentralityMap(const char* file)
{
  // Creates a mapping between centrality and Ncoll, based on a text file (taken from: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHeavyIonCentrality)
  
  if ( strlen(file) > 0 )
  {
    char line[1024];
    ifstream in(file);
    float lcent;
    float ucent;
    float Ncoll;
    
    fCentBins = 0;
    while ( in.getline(line,1024,'\n'))
    {
      sscanf(line,"%f %f %f",&lcent,&ucent,&Ncoll);
      
      fCentMap[ucent] = Ncoll;
      fCentBinning[fCentBins++] = ucent;
    }
    if ( fCentBins == 0 ) std::cout << "[INFO] No centrality map could be defined: The file provided is empty" << std::endl;
    else std::cout << "[INFO] Defining centrality map" << std::endl;
  }
  else
  {
    fCentBins = 0;
    std::cout << "[INFO] No centrality map could be defined: No file provided" << std::endl;
  }
};

bool readCorrection(const char* file)
{
  TFile *froot = new TFile(file,"READ");
  if (!froot)
  {
    cout << "[ERROR] File "<< file << " for correction of events not found" << endl;
    return false;
  }
  TList* lcorr = froot->GetListOfKeys();
  TIter nextCorr(lcorr);
  
  fcorrArray = new TObjArray();
  fcorrArray->SetOwner(kTRUE);
  
  TObjString* fname(0x0);
  while ( (fname = static_cast<TObjString*>(nextCorr.Next())) )
  {
  TEfficiency* h = static_cast<TEfficiency*>(froot->FindObjectAny(fname->GetString().Data()));

  TString sName(h->GetName());
  if ( sName.Contains("hcorr") ){ 
  fcorrArray->Add(h->Clone());
  cout<<"[INFO] Adding "<< sName << " to the correction array"<<endl;
  }
  else cout << "[WARNING] Correction histo " << sName.Data() << " not according to naming convention. Not included in correction array" << endl;
  }
  
  if (!(fcorrArray->GetSize()>0))
  {
  cout << "[ERROR] Correction array empty: No corrections found." << endl;
  return false;
  }
  delete lcorr;
  //froot->Close(); delete froot;
  return true;
};

double getCorr(Double_t rapidity, Double_t pt, Double_t mass, bool isPP)
{
  const char* collName = "PbPb";
  const char* massName = "Jpsi";
  if (isPP) collName = "PP";
  
    if (!fcorrArray)
    {
    cout << "[ERROR] No correction array exist" << endl;
    return 0;
    }

  Double_t corr = 1.;
  if (!strcmp(massName,"Interp"))
  {
    TEfficiency* corrHistoJpsi = static_cast<TEfficiency*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s",collName)));
    //TH2* corrHistoJpsi = static_cast<TH2*>(fcorrArray->FindObject("ptrapg_clone"));
    //TH2* corrHistoPsi2S = static_cast<TH2*>(fcorrArray->FindObject(Form("hcorr_Psi2S_%s",collName)));
    if (!corrHistoJpsi) // || !corrHistoPsi2S)
      {
	std::cout << "[Error] No histogram provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
	return 1.;
      }
    
    Int_t binJpsi = corrHistoJpsi->FindFixBin(rapidity, pt); // removed fabs()
    Double_t corrJpsi = corrHistoJpsi->GetEfficiency(binJpsi);
    
    //Int_t binPsi2S = corrHistoPsi2S->FindBin(fabs(rapidity), pt); //removed fabs()
    //Double_t corrPsi2S = corrHistoPsi2S->GetBinContent(binPsi2S);
    
    //corr = ((corrPsi2S - corrJpsi)/(3.5-3.3))*(mass-3.3) + corrJpsi;
    corr = corrJpsi;
  }
  else
    {
      //TH2* corrHisto = static_cast<TH2*>(fcorrArray->FindObject("ptrapg_clone"));
      //TF1  *bfrac = new TF1("bfrac","0.360511-0.348448*pow(x,0.5)+0.145442*x+-0.0134963*pow(x,1.5)", 3, 50);// Francois' polynomial
      
      TF1  *bfrac = new TF1("bfrac","exp(-2.74079+0.211476*pow(x,1)-0.007024*pow(x,2)+(7.90067e-05)*pow(x,3))", 3, 50);
      TEfficiency* prcorrHisto = static_cast<TEfficiency*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s_pr",collName)));
      TEfficiency* nprcorrHisto = static_cast<TEfficiency*>(fcorrArray->FindObject(Form("hcorr_Jpsi_%s_npr",collName)));
      double prEff = 1.0; double nprEff = 1.0; double bf = 1.0;
      if (!prcorrHisto || !nprcorrHisto)
	{
	  std::cout << "[Error] pr or npr histogram not provided for correction of " << collName << " " << massName << ". Weight set to 1." << std::endl;
	  return 1.;
	}
      if (pt > 3 && pt < 35 && abs(rapidity) < 2.4){
	prEff = prcorrHisto->GetEfficiency(prcorrHisto->FindFixBin(rapidity, pt));
	nprEff = nprcorrHisto->GetEfficiency(nprcorrHisto->FindFixBin(rapidity, pt));
	bf = bfrac->Eval(pt);
	corr = bf*nprEff + (1-bf)*prEff;
	//if (bf > 0.5 && bf <0.6) std::cout<<"[INFO] dimuon with pt = "<<pt<<", bfrac = "<<bf<<", prEff = "<<prEff<<", nprEff = "<<nprEff<<", corr = "<<corr<<endl;
      }
    }
  if(corr<0.00001) corr=1.0;
  return corr;
};

float jecCorr(double jtPt, double rawPt, double jpsiPt)
{
  return ( (1-(jpsiPt/rawPt))*jtPt + ((jpsiPt/rawPt)/(jtPt/rawPt))*jtPt );
}
