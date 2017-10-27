#include <string>
#include <algorithm>
#include "/afs/cern.ch/user/h/hinzmann/unfoldscript/src/RooUnfold.h"
#include "/afs/cern.ch/user/h/hinzmann/unfoldscript/src/RooUnfoldInvert.h"
#include "/afs/cern.ch/user/h/hinzmann/unfoldscript/src/RooUnfoldBayes.h"
#include "/afs/cern.ch/user/h/hinzmann/unfoldscript/src/RooUnfoldSvd.h"
#include "/afs/cern.ch/user/h/hinzmann/unfoldscript/src/RooUnfoldResponse.h"
#include <vector>
#include <cmath>
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TH1.h"
#include "TAttText.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TColor.h"
#include "TStyle.h"

//*********************************************************************************************
//**********************************************************************************************
std::string rootFiles[2] = {"Pasversion4Fig5and6Syst12FinalAgainData.root", "Pasversion4Fig5and6SystJESWResp.root"};

std::string outFileName = "test.root";

//std::string rootFiles[2] = {"Pasversion4Fig5and6Syst12FinalAgainData.root", "Pasversion4Fig5and6SystJESWResp.root"};
//std::string rootFiles[2] = {"Pasversion4Fig5and6Syst12FinalAgainData.root", "Pasversion4Fig5and6SystJERWResp.root"};
//std::string rootFiles[2] = {"Pasversion4Fig5and6Syst12FinalAgainData.root", "Pasversion4Fig5and6SystTrackptWResp.root"};
//std::string rootFiles[2] = {"Pasversion4Fig5and6Syst12FinalAgainData.root", "Pasversion4Fig5and6SystTrackreconeffWResp.root"};
//wtptHerwigcorraxis_Other2varModbinallpTmass.root
//std::string rootFiles[2] = {"Data12Finalcorraxis_Other2varModbinrebinallpTmass.root", "wtptFlatPy6JERunfoldcorraxis_Other2varModbindijetmassrebin.root"};
//wtptFlatPy6JESunfoldcorraxis_Other2varModbinrebin.root
//std::string rootFiles[2] = {"Data12Finalcorraxis_Other2varModbinrebinallpTmass.root", "wtptFlatPy6reconeffunfoldcorraxis_Other2varModbindijetmassrebin.root"};
//std::string rootFiles[2] = {"Data12Finalcorraxis_Other2varModbinrebinallpTmass.root", "wtptFlatPy6reconeffunfoldcorraxis_Other2varModbinrebin.root"};
//std::string rootFiles[2] = {"Data12finalcorraxis_Other2varModbinrebin.root", "wtptHerwigcorraxis_Other2varModbin.root"};
//std::string rootFiles[2] = {"Data12finalcorraxis_Other2varModbindijetmassrebin.root", "wtptFlatPy6corraxis_Other2varModbindijetmassrebin.root"}; 
//std::string rootFiles[2] = {"Data12finalcorraxis_Other2varModbinrebin.root", "wtptFlatPy6JESunfoldcorraxis_Other2varModbinrebin.root"};
//std::string rootFiles[2] = {"12Datacorraxis_Other2varModbinrebin.root", "wtptFlatPy6JESunfoldcorraxis_Other2varModbinrebin.root"};
//std::string rootFiles[2] = {"12Datacorraxis_Other2varModbinrebin.root", "wtptFlatPy6reconeffunfoldcorraxis_Other2varModbinrebin.root"};
//std::string varname[12] = {"jetchargeL1k6", "jetchargeL2k6", "GenchargeL1k6", "GenchargeL2k6", "recochargeL1k6vsgenchargeL1k6", "recochargeL2k6vsgenchargeL2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetchargeT1k6pt1000", "jetcharge2k6", "GenchargeT1k6pt1000", "Gencharge2k6", "recochargeT1k6vsgenchargeT1k6pt1000", "recocharge2k6vsgencharge2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetchargeT1k6", "jetchargeT2k6", "GenchargeT1k6", "GenchargeT2k6", "recochargeT1k6vsgenchargeT1k6", "recochargeT2k6vsgenchargeT2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetchargeT1", "jetcharge2", "GenchargeT1", "Gencharge2", "recochargeT1vsgenchargeT1", "recocharge2vsgencharge2", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetcharge1k6pt400", "jetcharge2", "Gencharge1k6pt400", "Gencharge2", "recocharge1k6vsgencharge1k6pt400", "recocharge2vsgencharge2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetchargeT1", "jetchargeT2", "GenchargeT1", "GenchargeT2", "recochargeT1vsgenchargeT1", "recochargeT2vsgenchargeT2", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetcharge1effplus", "jetcharge2", "Gencharge1", "Gencharge2", "recocharge1vsgencharge1effplus", "recocharge2vsgencharge2", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};

//std::string varname[12] = {"jetcharge1k6pt400", "jetcharge2k6pt400", "Gencharge1k6pt400", "Gencharge2k6pt400", "recocharge1k6vsgencharge1k6pt400", "recocharge2k6vsgencharge2k6pt400", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};

std::string varname[12] = {"jetchargeT1k6pt1000", "jetchargeL2k6pt400", "GenchargeT1k6pt1000", "GenchargeL2k6pt400", "recochargeT1k6vsgenchargeT1k6pt1000", "recochargeL2k6vsgenchargeL2k6pt400", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};


//name convention for eff and ptres
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//std::string varname[12] = {"jetchargeT1k6pt1000ptresminus", "jetcharge2k6", "GenchargeT1k6pt1000", "Gencharge2k6", "recochargeT1k6vsgenchargeT1k6pt1000ptresminus", "recocharge2k6vsgencharge2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//name convention for jer and jes and ptres
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//std::string varname[12] = {"jetcharge1k6jesuppt400", "jetcharge2k6", "Gencharge1k6pt400", "Gencharge2k6", "recocharge1k6vsgencharge1k6jesuppt400", "recocharge2k6vsgencharge2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//std::string varname[12] = {"jetcharge1k6", "jetcharge2k6", "Gencharge1k6", "Gencharge2k6", "recocharge1k6vsgencharge1k6", "recocharge2k6vsgencharge2k6", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetcharge1k3", "jetcharge2k3", "Gencharge1k3", "Gencharge2k3", "recocharge1k3vsgencharge1k3", "recocharge2k3vsgencharge2k3", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};
//std::string varname[12] = {"jetcharge1k3", "jetcharge2k3", "Gencharge1k3", "Gencharge2k3", "recocharge1k3vsgencharge1k3", "recocharge2k3vsgencharge2k3", "pT1", "pT2", "genpT1", "genpT2", "genpT1vsrecopT1", "genpT2vsrecopT2"};

struct histos {
  histos() {h_meas = h_data = h_true = 0; h_res = 0;}
  TH1D     *h_meas, *h_true, *h_data, *h_trdata, *h_gen;
  TH2D     *h_res, *h_res2;
};

//*******************************************************************************************

void jetcharge_Rebin(TH2D*, TH1D*, TH1D*, bool debug=false);
void ResetBinContent(TH1D*, bool debug=false);
void ChangeBin(char* name, TH2D* h_res, TH1D* h_meas, TH1D* h_true,
	       TH1D* h_data, TH1D* h_trdata, TH1D* h_gen, histos & out,
	       bool debug=false);
void compareResponse(TH2D*, TH2D*);

//********************************************************************************************

void compareResponse(TH2D *hR1, TH2D *hR2) {
  int nm = hR1->GetNbinsX();
  int nt = hR1->GetNbinsY();
  
  if (nm != hR2->GetNbinsX() || nt != hR2->GetNbinsY()) {
    std::cout << "Incompatibe # of bins X:" << nm << "|" << hR2->GetNbinsX() << " Y:" << nt << ":" << hR2->GetNbinsY() << std::endl;
    return;
  }
  
  bool ok(true);
  
  for (int i1=1; i1<=nm; ++i1) {
    for (int i2=1; i2<=nt; ++i2) {
      double diff = hR1->GetBinContent(i1,i2) - hR2->GetBinContent(i1,i2);
      if (std::abs(diff) > 0.001) ok = false;
    }
  }
  std::cout << hR1->GetTitle() << " compared to " << hR2->GetTitle() << " identical test " << ok << std::endl;
}

void jetcharge_Rebin(TH2D *h_res, TH1D *h_meas, TH1D *h_tr, bool debug) {
  
  if (debug)  std::cout << h_meas << " " << h_tr << " " << h_res << " " << std::endl;
  
  int nm = h_meas->GetNbinsX();
  int nt = (h_tr > 0) ? (h_tr->GetNbinsX()) : 0;
  TH1D *h_meas2 = (TH1D*)h_meas->Clone();
  h_meas->Reset();
  
  for (int i1=1; i1<=nm; ++i1) {
    double nmes=0, wmes=0;
    for (int i2=1; i2<=nt; ++i2) {
      
      nmes += h_res->GetBinContent(i1,i2);
      wmes += (h_res->GetBinError(i1,i2))*(h_res->GetBinError(i1,i2));
    }
    
    h_meas->SetBinContent(i1,nmes);
    h_meas->SetBinError(i1,std::sqrt(wmes));
    
    if (debug) std::cout << "Bin[" << i1 << "] = " << h_meas2->GetBinContent(i1) << "/" << h_meas->GetBinContent(i1) << " " << h_meas2->GetBinError(i1) << "/" << h_meas->GetBinError(i1) << std::endl;
  }
  
  if (h_tr != 0) {
    
    TH1D *h_tr2 = (TH1D*)h_tr->Clone();
    h_tr->Reset();
    
    for (int i2=1; i2<=nt; ++i2) {
      double nmes=0, wmes=0;
      for (int i1=1; i1<=nm; ++i1) {
	nmes += h_res->GetBinContent(i1,i2);
	wmes += (h_res->GetBinError(i1,i2))*(h_res->GetBinError(i1,i2));
      }
      
      h_tr->SetBinContent(i2,nmes);
      h_tr->SetBinError(i2,std::sqrt(wmes));
      h_tr->Draw();
      
      if (debug) std::cout << "Bin[" << i2 << "] = " << h_tr2->GetBinContent(i2) << "/" << h_tr->GetBinContent(i2) << " " << h_tr2->GetBinError(i2) << "/" << h_tr->GetBinError(i2) << std::endl;
    }
  }

  std::cout << " hmeas bins " << h_meas->GetNbinsX() << std::endl;
}

//*********************************************************************************

void ChangeBin(char* nameIn, TH2D* h_res, TH1D* h_meas, TH1D* h_true,
	       TH1D* h_data, TH1D* h_trdata, TH1D* h_gen, histos &out, bool debug) {
  
  if (debug)  std::cout << "ChangeBin for " << nameIn << " res " << h_res
                        << " meas " << h_meas << " true " << h_true << " data "
                        << h_data << std::endl;
  
  int nm = h_meas->GetNbinsX();

  std::cout << " Again h_meas bin " << nm << std::endl;
  
  int nt = h_true->GetNbinsX();
  
  int    nmStart(nm), nmEnd(0), ntStart(nt), ntEnd(0);
  double xmRange[200], xtRange[200];
  
  if (debug) std::cout << "Measured Histogram : " << h_meas->GetName() << std::endl;
  

  for (int i=1; i<=nm; ++i) {
    
    xmRange[i-1] = h_meas->GetXaxis()->GetBinLowEdge(i);
    xmRange[i]   = h_meas->GetXaxis()->GetBinUpEdge(i);
    
    if (h_meas->GetBinContent(i) > 0) {
      if (i < nmStart) nmStart = i;
      if (i > nmEnd)   nmEnd   = i;
    }
    
    if (debug) std::cout << "Bin[" << i << "] " << xmRange[i-1] << ":" << xmRange[i] << " Content " << h_meas->GetBinContent(i) << std::endl;

    std::cout << " nmstart " << nmStart << std::endl;
    std::cout << " nmEnd " << nmEnd << std::endl;
  }

  if (debug) std::cout << "True Histogram : " << h_true->GetName() << std::endl;
  
  for (int i=1; i<=nt; ++i) {
    
    xtRange[i-1] = h_true->GetXaxis()->GetBinLowEdge(i);
    xtRange[i]   = h_true->GetXaxis()->GetBinUpEdge(i);
    
    if (h_true->GetBinContent(i) > 0) {
      
      if (i < ntStart) ntStart = i;
      if (i > ntEnd)   ntEnd   = i;
    }
    
    if (debug) std::cout << "Bin[" << i << "] " << xtRange[i-1] << ":" << xtRange[i] << " Content " << h_true->GetBinContent(i) << std::endl;
  }
  
  if (debug) std::cout << "Measured: Total " << nm << " start " << nmStart << " end " << nmEnd << " True:Total " << nt << " start " << ntStart << " end " << ntEnd << std::endl;

    
  if (nmStart > ntStart) ntStart = nmStart;
  else                   nmStart = ntStart;
  if (nmEnd < ntEnd)     ntEnd   = nmEnd;
  else                   nmEnd   = ntEnd;


  int nmTot = (nmEnd - nmStart + 1);
  std::cout << " nmTot " << nmTot << std::endl;
  int ntTot = (ntEnd - ntStart + 1);
  
  if (debug) std::cout << "Measured: Total " << nmTot << " start " << nmStart << " end " << nmEnd << " True:Total " << ntTot << " start " << ntStart << " end " << ntEnd << std::endl;
  
  double binsM[200], binsT[200];
  
  for (int i=0; i<= nmTot; ++i) binsM[i] = xmRange[nmStart+i-1];
  for (int i=0; i<= ntTot; ++i) binsT[i] = xtRange[ntStart+i-1];
  
  char name[200], htit[400];
  
  sprintf (name, "h_meas_%s", nameIn);
  out.h_meas = new TH1D(name, h_meas->GetXaxis()->GetTitle(), nmTot, binsM);

  std::cout << " bins of out.h_meas " << out.h_meas->GetNbinsX() << std::endl;

  
  for (int i1=1; i1<=nmTot; ++i1) {
    
    out.h_meas->SetBinContent(i1, h_meas->GetBinContent(i1+nmStart-1));
    out.h_meas->SetBinError(i1, h_meas->GetBinError(i1+nmStart-1));
    
    if (debug) std::cout << "Meas Bin[" << i1 << "] = " << out.h_meas->GetBinContent(i1) << " +- " << out.h_meas->GetBinError(i1) << std::endl;
  }
  
  int hmeasbinx = out.h_meas->GetNbinsX();
  std::cout << " Again bins of out.h_meas " << out.h_meas->GetNbinsX() << std::endl;
  
  for (int i1=1; i1<=hmeasbinx; ++i1) {
    if (debug) std::cout << "Again Meas Bin[" << i1 << "] = " << out.h_meas->GetBinContent(i1) << " +- " << out.h_meas->GetBinError(i1) << std::endl;
  }
  
  sprintf (name, "h_true_%s", nameIn);
  out.h_true = new TH1D(name, h_true->GetXaxis()->GetTitle(), ntTot, binsT);
  if (h_gen != 0) {
    sprintf (name, "h_gen_%s", nameIn);
    out.h_gen = new TH1D(name, h_gen->GetXaxis()->GetTitle(), ntTot, binsT);
  }
  
  for (int i1=1; i1<=ntTot; ++i1) {
    
    out.h_true->SetBinContent(i1,h_true->GetBinContent(i1+ntStart-1));
    out.h_true->SetBinError(i1,h_true->GetBinError(i1+ntStart-1));
    
    if (h_gen != 0) {
      out.h_gen->SetBinContent(i1,h_gen->GetBinContent(i1+ntStart-1));
      out.h_gen->SetBinError(i1,h_gen->GetBinError(i1+ntStart-1));
    }
    
    if (debug) std::cout << "True Bin[" << i1 << "] = " << out.h_true->GetBinContent(i1) << " +- " << out.h_true->GetBinError(i1) << std::endl;
  }
  
  sprintf (name, "h_res_%s", nameIn);
  out.h_res = new TH2D(name, h_res->GetXaxis()->GetTitle(), nmTot, binsM, ntTot, binsT);
  
  for (int i2=1; i2<=ntTot; ++i2) {                                                                                                                 
    for (int i1=1; i1<=nmTot; ++i1) {
      
      out.h_res->SetBinContent(i1,i2,h_res->GetBinContent(nmStart+i1-1,ntStart+i2-1));
      out.h_res->SetBinError(i1,i2,h_res->GetBinError(nmStart+i1-1,ntStart+i2-1));
      
      if (debug) std::cout << "Reso Bin[" << i1 << "," << i2 << "] = " << out.h_res->GetBinContent(i1,i2) << " +- " << out.h_res->GetBinError(i1,i2) << std::endl;
    }
  }
  
  int hresbinx = out.h_res->GetXaxis()->GetNbins();
  int hresbiny = out.h_res->GetYaxis()->GetNbins();
  
  for (int i2=1; i2<=hresbinx; ++i2) {
    for (int i1=1; i1<=hresbiny; ++i1) {
      
      out.h_res->GetBinContent(i1,i2);
      if (debug) std::cout << "Again Reso Bin[" << i1 << "," << i2 << "] = " << out.h_res->GetBinContent(i1,i2) << " +- " << out.h_res->GetBinError(i1,i2) << std::endl;
    }
  }
  
  if (h_data != 0) {
    
    sprintf (name, "h_data_%s", nameIn);
    out.h_data = new TH1D(name, h_data->GetXaxis()->GetTitle(), nmTot, binsM);
    for (int i1=1; i1<=nmTot; ++i1) {
      out.h_data->SetBinContent(i1, h_data->GetBinContent(i1+nmStart-1));
      out.h_data->SetBinError(i1, h_data->GetBinError(i1+nmStart-1));
      if (debug) std::cout << "Data Bin[" << i1 << "] = " << out.h_data->GetBinContent(i1) << " +- " << out.h_data->GetBinError(i1) << std::endl;
    }
  }

  std::cout << " h_data binning " << out.h_data->GetNbinsX() << std::endl;
  std::cout << " Again again bins of out.h_meas " << out.h_meas->GetNbinsX() << std::endl;
}

//**************************************************************************************

void ResetBinContent(TH1D *histo, bool debug) {
  
  double value = 0.0, error = 0.0, errfac, width;
  histo->Scale(1/histo->Integral());
  for (int k=1; k<=histo->GetNbinsX(); k++) {
    value = histo->GetBinContent(k);
    error = histo->GetBinError(k);
    width = histo->GetBinWidth(k);
    if (debug) std::cout << "bin " << k << " content " << value << " width " << width  << std::endl;
    histo->SetBinContent(k, value/width);
    histo->SetBinError(k, error/width);
  }
}

//*****************************************************************************************
//**********************************************************************************************
//*****************************************************************

void StartUnfolding(int var=0, int algo=0) {
  
  char hname[100],hname1[100],hname2[100],title[100], title1[100];
  double chi2Cut(0.01);
  
  TH1D *hMeas(0), *hData(0);
  TH1D *hTr(0), *hGen(0);
  TH2D *hResponse(0), *hR2(0);
  
  int var1 = var + 4;
  
  for (int i=0; i<2; i++)
    {
      
      TFile *File = TFile::Open(rootFiles[i].c_str());
      std::cout << "[" << i << "] " << rootFiles[i] << std:: endl;
      
      if (!File) {
	
	std::cout << "Cannot find file " << rootFiles<< std::endl;
      } else {
	sprintf (hname, "%s", varname[var].c_str());
	sprintf (hname1, "%s", varname[var+2].c_str());
	sprintf (hname2, "%s", varname[var1].c_str());
	
	sprintf (title, "Unfolding of jetcharge1");
	sprintf (title1, "Unfolding of jetcharge1");
	
	std::cout << "Loading " << hname << std::endl;
	if (i==0)
	  {
	    hData = (TH1D*)(File->FindObjectAny(hname));
	    hData->SetName(((std::string(hData->GetName()))+"0").c_str());
	  }
	
	if (i==1)
	  {
	    hMeas = (TH1D*)(File->FindObjectAny(hname));
	    hMeas->SetName(((std::string(hMeas->GetName()))+"1").c_str());
	    hTr   = (TH1D*)(File->FindObjectAny(hname1)); 
	    hTr->SetName(((std::string(hTr->GetName()))+"1").c_str());
            hResponse  = (TH2D*)(File->FindObjectAny(hname2));
	    hResponse->SetName(((std::string(hResponse->GetName()))+"1").c_str());
	  }  
      }
    }

  std::cout << " hMeas bins " << hMeas->GetNbinsX() << std::endl;
  //jetcharge_Rebin(hResponse, hMeas, hTr, false);
  
  histos out;
  
  ChangeBin(title, hResponse, hMeas, hTr, hData, 0, hGen, out, true);
  RooUnfoldResponse response(out.h_meas, out.h_true, out.h_res, title, title);
  
  //TVectorD cov = RooUnfold::ErecoV(0);
  //cov->Draw();
  //RooUnfoldResponse response(0, out.h_true, out.h_res, title, title);
  //Unfolding code taking into account different algorithms
  
  char histname[100];
  
  //*****************************Different Algorithms****************************************************************
  
  if (algo == 0) //Bayes Method
    {
      TH1D *h_unfolded[10];
      
      for (int j=0; j<10; j++)
        {
          h_unfolded[j] = 0;
        }
      
      for (int it=1; it<=10; ++it) {
	    RooUnfoldBayes unfold_bayes (&response, out.h_data, it);
	    unfold_bayes.IncludeSystematics(2);
	    h_unfolded[it-1] = (TH1D*)unfold_bayes.Hreco(RooUnfold::kCovariance);
	    sprintf (histname, "jetcharge1 with bayesian method iteration %d", it);
	    h_unfolded[it-1]->SetName(histname);
	    if (it == 5)
	      {
		TMatrixD covmatr = unfold_bayes.Ereco(RooUnfold::kCovariance);
                covmatr.Print();
	      }
	    //matr = unfold_bayes.Ereco(2);
      }
      
      TFile* File1 = new TFile(outFileName.c_str(), "RECREATE");
      
      for (int k=0; k<10; k++)
      {
	if (k == 4)
	  h_unfolded[k]->Write("unfolded");
      }
      covmatr.Write("cov");
      out.h_res->Write();
      out.h_meas->Write();
      out.h_true->Write();
      out.h_data->Write();
      
      File1->Close();
      
    }
  
  if (algo == 1) //SVD Method
    {
      
      TH1D *h_unfolded[10];
      for (int j=0; j<10; j++)
        {
          h_unfolded[j] = 0;
        }
      
      for (int kreg=1; kreg<=10; ++kreg) {
      	RooUnfoldSvd unfold_svd (&response, out.h_data, kreg);
	sprintf (histname, "jetcharge1 with svd method kreg %d", kreg);
	h_unfolded[kreg-1] = (TH1D*)unfold_svd.Hreco(RooUnfold::kCovToy);
        h_unfolded[kreg-1]->SetName(histname);
      }

      TFile* File1 = new TFile("dataunfoldPy6_Q1k3FNAL_SVD.root", "RECREATE");
      
      for (int k=0; k<10; k++)
	{
	  h_unfolded[k]->Write();
	}
      
      out.h_res->Write();
      out.h_meas->Write();
      out.h_true->Write();
      out.h_data->Write();
      File1->Close();
    }

  if (algo == 2) // Matrix Inversion Method 
    {
      
      TH1D *h_unfolded(0);
      RooUnfoldInvert unfold_invert (&response, out.h_data);
      sprintf (histname, "jetcharge1 with matrix inversion method kreg");
       h_unfolded = (TH1D*)unfold_invert.Hreco(RooUnfold::kCovToy);
      h_unfolded->SetName(histname);
      TFile* File1 = new TFile("dataunfoldPy6_Q1FNAL_INV.root", "RECREATE");
      h_unfolded->Write();
      out.h_res->Write();
      out.h_meas->Write();
      out.h_true->Write();
      out.h_data->Write();
      File1->Close();
    }
}

void chargeunfold()
{
  varname[0] = "jetcharge1";
  varname[2] = "Gencharge1";
  varname[4] = "recocharge1vsgencharge1";
  outFileName="jetcharge1.root";
  StartUnfolding();
  varname[0] = "jetcharge1k6";
  varname[2] = "Gencharge1k6";
  varname[4] = "recocharge1k6vsgencharge1k6";
  outFileName="jetcharge1k6.root";
  StartUnfolding();
  varname[0] = "jetcharge1k3";
  varname[2] = "Gencharge1k3";
  varname[4] = "recocharge1k3vsgencharge1k3";
  outFileName="jetcharge1k3.root";
  StartUnfolding();
  varname[0] = "jetchargeL1";
  varname[2] = "GenchargeL1";
  varname[4] = "recochargeL1vsgenchargeL1";
  outFileName="jetchargeL1.root";
  StartUnfolding();
  varname[0] = "jetchargeL1k6";
  varname[2] = "GenchargeL1k6";
  varname[4] = "recochargeL1k6vsgenchargeL1k6";
  outFileName="jetchargeL1k6.root";
  StartUnfolding();
  varname[0] = "jetchargeL1k3";
  varname[2] = "GenchargeL1k3";
  varname[4] = "recochargeL1k3vsgenchargeL1k3";
  outFileName="jetchargeL1k3.root";
  StartUnfolding();
  varname[0] = "jetchargeT1";
  varname[2] = "GenchargeT1";
  varname[4] = "recochargeT1vsgenchargeT1";
  outFileName="jetchargeT1.root";
  StartUnfolding();
  varname[0] = "jetchargeT1k6";
  varname[2] = "GenchargeT1k6";
  varname[4] = "recochargeT1k6vsgenchargeT1k6";
  outFileName="jetchargeT1k6.root";
  StartUnfolding();
  varname[0] = "jetchargeT1k3";
  varname[2] = "GenchargeT1k3";
  varname[4] = "recochargeT1k3vsgenchargeT1k3";
  outFileName="jetchargeT1k3.root";
  StartUnfolding();
  varname[0] = "jetcharge1k6pt400";
  varname[2] = "Gencharge1k6pt400";
  varname[4] = "recocharge1k6vsgencharge1k6pt400";
  outFileName="jetcharge1k6pt400.root";
  StartUnfolding();
  varname[0] = "jetcharge1k6pt700";
  varname[2] = "Gencharge1k6pt700";
  varname[4] = "recocharge1k6vsgencharge1k6pt700";
  outFileName="jetcharge1k6pt700.root";
  StartUnfolding();
  varname[0] = "jetcharge1k6pt1000";
  varname[2] = "Gencharge1k6pt1000";
  varname[4] = "recocharge1k6vsgencharge1k6pt1000";
  outFileName="jetcharge1k6pt1000.root";
  StartUnfolding();
  varname[0] = "jetchargeL1k6pt400";
  varname[2] = "GenchargeL1k6pt400";
  varname[4] = "recochargeL1k6vsgenchargeL1k6pt400";
  outFileName="jetchargeL1k6pt400.root";
  StartUnfolding();
  varname[0] = "jetchargeL1k6pt700";
  varname[2] = "GenchargeL1k6pt700";
  varname[4] = "recochargeL1k6vsgenchargeL1k6pt700";
  outFileName="jetchargeL1k6pt700.root";
  StartUnfolding();
  varname[0] = "jetchargeL1k6pt1000";
  varname[2] = "GenchargeL1k6pt1000";
  varname[4] = "recochargeL1k6vsgenchargeL1k6pt1000";
  outFileName="jetchargeL1k6pt1000.root";
  StartUnfolding();
  varname[0] = "jetchargeT1k6pt400";
  varname[2] = "GenchargeT1k6pt400";
  varname[4] = "recochargeT1k6vsgenchargeT1k6pt400";
  outFileName="jetchargeT1k6pt400.root";
  StartUnfolding();
  varname[0] = "jetchargeT1k6pt700";
  varname[2] = "GenchargeT1k6pt700";
  varname[4] = "recochargeT1k6vsgenchargeT1k6pt700";
  outFileName="jetchargeT1k6pt700.root";
  StartUnfolding();
  varname[0] = "jetchargeT1k6pt1000";
  varname[2] = "GenchargeT1k6pt1000";
  varname[4] = "recochargeT1k6vsgenchargeT1k6pt1000";
  outFileName="jetchargeT1k6pt1000.root";
  StartUnfolding();
}
