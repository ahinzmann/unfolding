from ROOT import *
import ROOT
import array, math
import os
from math import *

if __name__=="__main__":

    print "start ROOT"
    gROOT.Reset()
    gROOT.SetStyle("Plain")
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetTitleOffset(1.2,"Y")
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetPadBottomMargin(0.20)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.23)
    gStyle.SetMarkerSize(2.5)
    gStyle.SetHistLineWidth(1)
    gStyle.SetStatFontSize(0.020)
    gStyle.SetTitleSize(0.05, "XYZ")
    gStyle.SetLabelSize(0.04, "XYZ")
    gStyle.SetLegendBorderSize(0)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetEndErrorSize(5)

    stops = [0.0, 0.5, 1.0]
    red   = [0.0, 1.0, 1.0]
    green = [0.0, 1.0, 0.0]
    blue  = [1.0, 1.0, 0.0]
    s = array.array('d', stops)
    r = array.array('d', red)
    g = array.array('d', green)
    b = array.array('d', blue)
    TColor.CreateGradientColorTable(len(s), s, r, g, b, 999)
    gStyle.SetNumberContours(999)
    #gStyle.SetPalette(1)

    print "start CMS_lumi"

    gROOT.LoadMacro("CMS_lumi.C");
    iPeriod = 2;	#// 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV 
    iPos = 1;
    
    plots=[
      ("jetcharge1", "Leading-jet Q^{#kappa=1.0} [|e|]",1,""),
      ("jetcharge1k6", "Leading-jet Q^{#kappa=0.6} [|e|]",1.5,""),
      ("jetcharge1k3", "Leading-jet Q^{#kappa=0.3} [|e|]",2.5,""),
      ("jetchargeL1", "Leading-jet Q_{L}^{#kappa=1.0} [|e|]",1,""),
      ("jetchargeL1k6", "Leading-jet Q_{L}^{#kappa=0.6} [|e|]",0.8,""),
      ("jetchargeL1k3", "Leading-jet Q_{L}^{#kappa=0.3} [|e|]",0.7,""),
      ("jetchargeT1", "Leading-jet Q_{T}^{#kappa=1.0} [|e|]",1.0,""),
      ("jetchargeT1k6", "Leading-jet Q_{T}^{#kappa=0.6} [|e|]",0.7,""),
      ("jetchargeT1k3", "Leading-jet Q_{T}^{#kappa=0.3} [|e|]",0.7,""),
      ("jetcharge1k6pt400", "Leading-jet Q^{#kappa=0.6} [|e|]",1.5,"400 < p_{T} < 700 GeV"),
      ("jetcharge1k6pt700", "Leading-jet Q^{#kappa=0.6} [|e|]",1.5,"700 < p_{T} < 1000 GeV"),
      ("jetcharge1k6pt1000", "Leading-jet Q^{#kappa=0.6} [|e|]",1.5,"p_{T} > 1000 GeV"),
      ("jetchargeL1k6pt400", "Leading-jet Q_{L}^{#kappa=0.6} [|e|]",0.8,"400 < p_{T} < 700 GeV"),
      ("jetchargeL1k6pt700", "Leading-jet Q_{T}^{#kappa=0.6} [|e|]",0.8,"700 < p_{T} < 1000 GeV"),
      ("jetchargeL1k6pt1000", "Leading-jet Q_{T}^{#kappa=0.6} [|e|]",0.8,"p_{T} > 1000 GeV"),
      ("jetchargeT1k6pt400", "Leading-jet Q_{T}^{#kappa=0.6} [|e|]",0.7,"400 < p_{T} < 700 GeV"),
      ("jetchargeT1k6pt700", "Leading-jet Q_{T}^{#kappa=0.6} [|e|]",0.7,"700 < p_{T} < 1000 GeV"),
      ("jetchargeT1k6pt1000", "Leading-jet Q_{T}^{#kappa=0.6} [|e|]",0.7,"p_{T} > 1000 GeV"),
    ]
    for plot,label,axis,caption in plots:
    
      f=TFile.Open(plot+".root")

      t1=f.Get("unfolded")

      c = TCanvas("combined", "combined", 0, 0, 300, 300)
    
      t1.GetXaxis().SetTitle(label)
      t1.GetXaxis().SetRangeUser(-axis,axis)
      t1.SetTitle("")
      t1.Draw("hist")
      
      l=TLegend(0.2,0.75,0.5,0.85,caption)
      l.SetFillStyle(0)
      l.SetTextSize(0.04)
      l.Draw("same")
        
      #// writing the lumi information and the CMS "logo"
      CMS_lumi( c, iPeriod, iPos );
    
      c.SaveAs(plot+"_unfolded.pdf")
      c.SaveAs(plot+"_unfolded.eps")
    
      m=f.Get("cov")

      c = TCanvas("combined", "combined", 0, 0, 300, 300)
      
      t=TH2D("covariance","",t1.GetXaxis().GetNbins(),t1.GetXaxis().GetXbins().GetArray(),t1.GetXaxis().GetNbins(),t1.GetXaxis().GetXbins().GetArray())
      print plot
      print "bin &",
      for x in range(t1.GetXaxis().GetNbins()):
       if t.GetXaxis().GetBinCenter(x+1)>-axis and t.GetXaxis().GetBinCenter(x+1)<axis:
        print t.GetXaxis().GetBinCenter(x+1),"&",
      print "\\\\ \hline"
      for x in range(t1.GetXaxis().GetNbins()):
       if t.GetXaxis().GetBinCenter(x+1)>-axis and t.GetXaxis().GetBinCenter(x+1)<axis:
        print t.GetXaxis().GetBinCenter(x+1),"&",
       for y in range(t1.GetXaxis().GetNbins()):
        t.SetBinContent(x+1,y+1,m[y][x]/t1.Integral()/t1.Integral()/t1.GetXaxis().GetBinWidth(x+1)/t1.GetXaxis().GetBinWidth(y+1))
	if t.GetXaxis().GetBinCenter(x+1)>-axis and t.GetXaxis().GetBinCenter(x+1)<axis and t.GetXaxis().GetBinCenter(y+1)>-axis and t.GetXaxis().GetBinCenter(y+1)<axis:
	  print("%.3g &" % (t.GetBinContent(x+1,y+1))),
       print "\\\\"
      for x in range(t1.GetXaxis().GetNbins()):
       for y in range(t1.GetXaxis().GetNbins()):
	upx=t.GetXaxis().GetBinUpEdge(x+1)
        upy=t.GetYaxis().GetBinUpEdge(y+1)
        if plot=="jetcharge1k3":
         if str(t.GetXaxis().GetBinLowEdge(x+1)) in ["-2.88","-2.76","-2.64","-2.4","-2.28","-2.16","2.16","2.28","2.4","2.52","2.64","2.76","2.88"]: continue
         if str(t.GetYaxis().GetBinLowEdge(y+1)) in ["-2.88","-2.76","-2.64","-2.4","-2.28","-2.16","2.16","2.28","2.4","2.52","2.64","2.76","2.88"]: continue
         if t.GetXaxis().GetBinLowEdge(x+1) in [-2.52]:
	   upx=-2.04
         if t.GetXaxis().GetBinLowEdge(x+1) in [2.04]:
           upx=2.52
         if t.GetYaxis().GetBinLowEdge(y+1) in [-2.52]:
           upy=-2.04
         if t.GetYaxis().GetBinLowEdge(y+1) in [2.04]:
           upy=2.52
        if plot=="jetchargeL1k6" or plot=="jetchargeL1k6pt400" or plot=="jetchargeL1k6pt700" or plot=="jetchargeL1k6pt1000":
         if str(t.GetXaxis().GetBinLowEdge(x+1)) in ["-0.95","0.8","0.95"]: continue
         if str(t.GetYaxis().GetBinLowEdge(y+1)) in ["-0.95","0.8","0.95"]: continue
        if plot=="jetchargeL1k3" or plot=="jetchargeT1k6" or plot=="jetchargeT1k3" or plot=="jetchargeT1k6pt400":
         if str(t.GetXaxis().GetBinLowEdge(x+1)) in ["-0.95","-0.9","-0.8","0.7","0.8","0.9","-0.95"]: continue
         if str(t.GetYaxis().GetBinLowEdge(y+1)) in ["-0.95","-0.9","-0.8","0.7","0.8","0.9","-0.95"]: continue
	print("%.3g TO %.3g; %.3g TO %.3g; %.3g;" % (t.GetXaxis().GetBinLowEdge(x+1),upx,t.GetYaxis().GetBinLowEdge(y+1),upy,t.GetBinContent(x+1,y+1)))
      t.GetXaxis().SetTitle(label)
      t.GetYaxis().SetTitle(label)
      t.GetXaxis().SetRangeUser(-axis,axis)
      t.GetYaxis().SetRangeUser(-axis,axis)
      t.GetZaxis().SetRangeUser(-max(-t.GetMinimum(),t.GetMaximum()),max(-t.GetMinimum(),t.GetMaximum()))
      t.Draw("colz")
      
      l.Draw("same")
        
      #// writing the lumi information and the CMS "logo"
      CMS_lumi( c, iPeriod, iPos );
    
      c.SaveAs(plot+"_covariance.pdf")
      c.SaveAs(plot+"_covariance.eps")
    