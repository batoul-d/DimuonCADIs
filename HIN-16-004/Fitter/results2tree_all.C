#ifndef results2tree_all_C
#define results2tree_all_C

#include <iostream>

#include "TString.h"

#include "results2tree.C"

using namespace std;

void results2tree_all()
{ 
  TString workdir("");
  ///////////////////////////////
  //        Nominal            //
  ///////////////////////////////
  workdir = "DataFits/DataFits_midJtPt/DataFits_016";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = Form("DataFits/DataFits_midJtPt/DataFits_1624");
  // results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  //low and high Jt pt
  workdir = "DataFits/DataFits_lowJtPt/DataFits_016";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = Form("DataFits/DataFits_lowJtPt/DataFits_1624");
  // results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits/DataFits_highJtPt/DataFits_016";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1); 
  workdir = Form("DataFits/DataFits_highJtPt/DataFits_1624");
  // results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);

  ///////////////////////////////
  //   Background variations   //
  ///////////////////////////////
  workdir = "DataFits_bkgExp/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_bkgExp/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_pval25/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_pval25/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_pval100/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_pval100/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);

  ///////////////////////////////
  //    Signal variations      //
  ///////////////////////////////
  workdir = "DataFits_constrained/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_constrained/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_CBGauss/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_CBGauss/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);

  ///////////////////////////////
  //    CtauErr variations     //
  ///////////////////////////////
  workdir = "DataFits_ctauErrTot/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_ctauErrTot/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);

  ///////////////////////////////
  //   CtauTrue variations     //
  ///////////////////////////////
  workdir = "DataFits_ctauReco/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_ctauReco/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  
  ///////////////////////////////
  //   CtauRes variations      //
  ///////////////////////////////
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  //workdir = "";
  //results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);

  ///////////////////////////////
  //   CtauBkg variations      //
  ///////////////////////////////
  workdir = "DataFits_ctauBkgTemplate/DataFits_midJtPt/DataFits_016";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);
  workdir = "DataFits_ctauBkgTemplate/DataFits_midJtPt/DataFits_1624";
  results2tree(workdir.Data(), "DATA", "", "ctauMass", 0, "AccEff", 1);

  return;
}

#endif
