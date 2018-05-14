from ROOT import *
from os import mkdir, chdir
import copy, math
import tdrstyle
from array import array

gROOT.SetBatch()
tdrstyle.setTDRStyle()

myColor = [1,kRed-4,kBlue-4,kGreen+3,kMagenta-4, kOrange+7,kCyan-7, kGray, kRed-7, kBlue-7,kMagenta-9]
mySysColor = [kYellow-7,kGreen-7]
mySysColor = [kOrange+1,kGreen+1]
myMarker = [20,21,23,26,32,1,1,1]
myLine = [1,1,1,2,2,7,7]
myLine = [1]*10
myLine = [1,2,3,4,8,7,7]
myLine = [1]*10

def reBin(h1,h2):
  tmp = h1.Clone(h2.GetName())
  tmp.Reset()
  for x in xrange(h2.GetNbinsX()):
    v = h2.GetBinContent(x+1)
    b = h2.GetBinCenter(x+1)
    tmp.Fill(b,v)
  return tmp

def setOnBin(onList,hist):
  for ibin, xbin in enumerate(onList):
    if xbin == 0: hist.SetBinContent(ibin+1,0)

def drawCut(cut):
  if cut == "dr_pt3": t_cut = ""
  if cut == "sdr_pt3": t_cut = "#Delta R_{23} < 1.0"
  if cut == "ldr_pt3": t_cut = "1.0 < #Delta R_{23}"
  if cut == "dr_lpt3": t_cut = "p_{T3}/p_{T2} < 0.3"
  if cut == "dr_hpt3": t_cut = "0.6 < p_{T3}/p_{T2}"
  if cut == "dr_mpt3": t_cut = "0.3 < p_{T3}/p_{T2} < 0.6"
  return [t_cut]

def drawLumi(text_l, x, y):
  cmsText = TLatex()
  cmsText.SetNDC()
  cmsText.SetTextFont(61)
  cmsText.SetTextSize(0.035)
  cmsText.DrawLatex(x-0.03,y+0.05,"CMS")
  
  extraText = TLatex()
  extraText.SetNDC()
  extraText.SetTextFont(52)
  extraText.SetTextSize(0.03)
  extraText.DrawLatex(x+0.05,y+0.05,text_l[0]) 

  lumiText = TLatex()
  lumiText.SetNDC()
  lumiText.SetTextFont(42)
  lumiText.SetTextSize(0.03)
  lumiText.DrawLatex(x+0.55,y+0.05,text_l[1])
  #lumiText.DrawLatex(x+0.426,y+0.05,text_l[1])

  cutText =TLatex()
  cutText.SetNDC()
  cutText.SetTextFont(42)
  cutText.SetTextSize(0.03)
  cutText.DrawLatex(x,y-0.45,text_l[2])
  #cutText.DrawLatex(x,y-0.03,text_l[2])

  cutText2 =TLatex()
  cutText2.SetNDC()
  cutText2.SetTextFont(42)
  cutText2.SetTextSize(0.03)
  cutText2.DrawLatex(x,y-0.5,"510 GeV < p_{T1} < 2500 GeV")
  #cutText2.DrawLatex(x,y-0.08,"510 GeV < p_{T1} < 2500 GeV")

def drawLumi2(text_l, x, y):
  cmsText = TLatex()
  cmsText.SetNDC()
  cmsText.SetTextFont(61)
  cmsText.SetTextSize(0.035)
  cmsText.DrawLatex(x-0.03,y+0.05,"CMS")

  extraText = TLatex()
  extraText.SetNDC()
  extraText.SetTextFont(52)
  extraText.SetTextSize(0.03)
  extraText.DrawLatex(x+0.05,y+0.05,text_l[0])

  lumiText = TLatex()
  lumiText.SetNDC()
  lumiText.SetTextFont(42)
  lumiText.SetTextSize(0.03)
  lumiText.DrawLatex(x+0.426,y+0.05,text_l[1])

  cutText =TLatex()
  cutText.SetNDC()
  cutText.SetTextFont(42)
  cutText.SetTextSize(0.03)
  cutText.DrawLatex(x,y-0.03,text_l[2])

  cutText2 =TLatex()
  cutText2.SetNDC()
  cutText2.SetTextFont(42)
  cutText2.SetTextSize(0.03)
  cutText2.DrawLatex(x,y-0.08,"510 GeV < p_{T1} < 2500 GeV")

def calChi2(histList):
  chi2L =[] 
  mV = []
  mE = []
  ent = 0
  for x in xrange(len(histList)-1):
    chi2L.append(0.0)
    mV.append(0.0)
    mE.append(0.0)
  for x in xrange(histList[0].GetNbinsX()):
    dV = histList[0].GetBinContent(x+1)
    dE = histList[0].GetBinError(x+1)
    for mc in xrange(len(histList)-1):
      mV[mc] = histList[mc+1].GetBinContent(x+1)
      mE[mc] = histList[mc+1].GetBinError(x+1)
    if dV == 0 or dE == 0 : continue
    for mc in xrange(len(histList)-1):
      chi2L[mc] += (dV - mV[mc])**2/(dE**2 + mE[mc]**2)
    ent += 1
  return chi2L, ent-1, ent

def setError(x):
  for bin in xrange(x.GetNbinsX()):
    er = math.sqrt(x.GetBinContent(bin+1))
    x.SetBinError(bin+1,er)
    #x.Sumw2()
  return x

def normHist(hist):
  intg = hist.Integral()
  hist.Scale(1.0/intg)
  for x in xrange(hist.GetNbinsX()):
    binW = hist.GetBinWidth(x+1)
    binC = hist.GetBinContent(x+1)
    binE = hist.GetBinError(x+1)
    hist.SetBinContent(x+1, binC/binW)
    hist.SetBinError(x+1, binE/binW)

class plotMaker:
  def __init__(self,cx_size = 600, cy_size = 700) :
    gROOT.Reset()
    self.pad = TCanvas("c1", "c1", cx_size, cy_size)
    self.drawLog = False
    self.sysDraw = False
    self.sysName = []
    self.rMax = 1.51
    self.rMin = 0.49
    self.norDraw = True
    self.histList = []
    self.divHistList = []
    self.mPlotBMargin = 0.01
    self.mPlotTMargin = 0.1
    self.divPlotSet = True
    self.myColor = myColor
    self.myMarker = myMarker
    self.mySysColor = mySysColor    
    self.myLine = myLine
    self.allSysList = []
    self.padList = []
    self.padH = []
    self.padS = []
    self.sysGRList = []
    self.sysGRDivList = []
    self.divTitle = "MC / Data"
    self.mP1 = 1.1
    self.mP2 = 0.05
    self.mP3 = 0.035
    self.mP4 = 1.8
    self.mP5 = 1.8

  def setDivTitle(self, title):
    self.divTitle = title

  def setRange(self, max, min):
    self.rMax = max+0.01
    self.rMin = min-0.01

  def setNorDraw(self, norSet):
    self.norDraw = norSet

  def setLogByName(self, name):
    tmp = name.split("_")    
    logList = ["pt", "met", "mass"]
    try:
      if logList.index(tmp): self.drawLog = True
    except: self.drawLog = False
    if name.endswith("jet3_pt_jet2_pt"):
      self.setLog_ = False

  def setHistList(self, histList):
    self.histList = histList
    self.histStSet() 
    tmp0 = self.histList[0]
    for i,x in enumerate(self.histList):
      tmp1 = x.Clone(x.GetName()+"_div{}".format(i)) 
      tmp1.Divide(tmp0)
      self.divHistList.append(tmp1)

  def setColor(self, colorList):
    self.myColor = colorList
  def setMarker(self, markerList):
    self.myMarker = markerList
  def setLine(self, lineList):
    self.myLine = lineList

  def histStSet(self):
    for i, x in enumerate(self.histList):
      x.SetLineColor(self.myColor[i])
      x.SetLineWidth(2)
      x.SetMarkerStyle(self.myMarker[i])
      x.SetMarkerColor(self.myColor[i])
      x.SetLineStyle(self.myLine[i])
      x.SetStats(0)
      if self.norDraw:
        normHist(x)

  def legMaker(self, nHist):
    self.le = TLegend(0.45, 0.9-nHist*0.06, 0.99, 0.9)
    self.le.SetTextSize(0.035)
    self.le.SetFillStyle(0)
    self.le.SetFillColor(kWhite)
    self.le.SetBorderSize(0)
    return self.le

  def findMaxMin(self):
    tmpMax = 0 
    tmpMin = 10
    for x in self.histList:
      hMax = x.GetMaximum()
      hMin = x.GetMinimum()
      if tmpMax < hMax: tmpMax = hMax
      if tmpMin > hMin: tmpMin = hMin
    self.maxV = tmpMax*1.8
    self.minV = 0.0001
    if tmpMin == 0: self.minV = 0.05

  def sysGR(self):
    self.sysGRList = []
    nBin = self.histList[0].GetNbinsX()
    xBin = []
    xEr = []
    yBin = []
    for x in xrange(nBin):
      xBin.append(self.histList[0].GetBinCenter(x+1)) 
      xEr.append(self.histList[0].GetBinWidth(x+1)/2.0)
      yBin.append(self.histList[0].GetBinContent(x+1))
    for i, sysList in enumerate(self.sysSet[1]):
      yEr = []
      for ibin, sysV in enumerate(sysList):
        yEr.append(float(sysV)*yBin[ibin])
      tmpGR = TGraphAsymmErrors(nBin, array("d", xBin), array("d", yBin), array("d", xEr), array("d", xEr), array("d", yEr), array("d", yEr))
      tmpGR.SetFillColor(self.mySysColor[i])
      tmpGR.SetMarkerStyle(1)
      tmpGR.SetMarkerColor(self.mySysColor[i])
      tmpGR.SetLineColor(self.mySysColor[i])
      self.sysGRList.append(copy.deepcopy(tmpGR))
      del(tmpGR)
    return self.sysGRList

  def sysGRDiv(self):
    self.sysGRDivList = []
    nBin = self.histList[0].GetNbinsX()
    xBin = []
    xEr = []
    yBin = []
    for x in xrange(nBin):
      xBin.append(self.histList[0].GetBinCenter(x+1)) 
      xEr.append(self.histList[0].GetBinWidth(x+1)/2.0)
      yBin.append(1.0)
    for i, sysList in enumerate(self.sysSet[1]):
      yEr = []
      for ibin, sysV in enumerate(sysList):
        yEr.append(float(sysV))
      tmpGR = TGraphAsymmErrors(nBin, array("d", xBin), array("d", yBin), array("d", xEr), array("d", xEr), array("d", yEr), array("d", yEr))
      tmpGR.SetFillColor(self.mySysColor[i])
      tmpGR.SetMarkerStyle(1)
      tmpGR.SetMarkerColor(self.mySysColor[i])
      tmpGR.SetLineColor(self.mySysColor[i])
      self.sysGRDivList.append(copy.deepcopy(tmpGR))
      del(tmpGR)
    return self.sysGRDivList

  def setSysDraw(self, sysSet):
    self.sysSet = sysSet
    self.sysDraw = sysSet[0]
    self.sysName = sysSet[2]

  def histTxt(self, tag):
    outBin = open(tag+"_"+self.histList[0].GetName()+"_bin.txt","w")
    for x in xrange(self.histList[0].GetNbinsX()):
      outBin.write("%02d bin : %.5f\n"%(x+1, self.histList[0].GetBinContent(x+1)))
    outBin.close()    
    if len(self.histList) > 2:
      for i,x in enumerate(self.divHistList):
        outTxt = open(tag+"_"+self.histList[i].GetName()+".txt","w")
        for y in xrange(x.GetNbinsX()):
          if self.histList[0].GetBinContent(y+1) == 0:
            tmpEnt = 1
          else: tmpEnt = self.histList[0].GetBinContent(y+1)
          outTxt.write("%02d bin value, error, statistics : %3.8f, %3.8f, %3.8f\n"%(y+1, x.GetBinContent(y+1), x.GetBinError(y+1), self.histList[0].GetBinError(y+1)/tmpEnt))
        outTxt.close()
    chiR = calChi2(self.histList)
    outChi = open(tag+"_"+self.histList[0].GetName()+"_chi2.txt", "w")
    for x in chiR[0]:
      outChi.write("chi2 : %f, %f\n"%(x/chiR[1], x))
  def mPlotSet(self, setList):
    self.mPlotTMargin = setList[0]
    self.mPlotBMargin = setList[1]
    self.mP1 = setList[2]
    self.mP2 = setList[3]
    self.mP3 = setList[4]
    self.mP4 = setList[5]
    self.mP5 = setList[6]

  def mPlot(self, s, plotHist):
    self.pad.cd()
    padt = TPad("t", "t", s[0], s[1], s[2], s[3])
    padt.Draw()
    padt.cd()
    self.setLogByName(plotHist[0].GetName())
    padt.SetLogy(self.drawLog)
    #padt.SetLogy(1)
    padt.SetBottomMargin(self.mPlotBMargin)
    padt.SetTopMargin(self.mPlotTMargin)
    
    plotHist[0].GetYaxis().SetTitleOffset(self.mP1)
    plotHist[0].GetYaxis().SetTitleSize(self.mP2) 
    plotHist[0].GetYaxis().SetLabelSize(self.mP3)
    plotHist[0].GetXaxis().SetLabelOffset(self.mP4)
    plotHist[0].GetXaxis().SetTitleOffset(self.mP5)
    plotHist[0].SetMaximum(self.maxV/1.5)
    #if self.maxV/1.5 < 5.0:
    plotHist[0].SetMaximum(3.55)
    #plotHist[0].SetMaximum(5.55)
    plotHist[0].SetMinimum(self.minV) 
    plotHist[0].Draw("HIST][")
    if self.sysDraw:
      self.sysGR()
      for x in reversed(self.sysGRList):
        x.Draw("P2")
    plotHist[0].Draw("same")
    plotHist[0].Draw("HIST][ same")
    #plotHist[0].Draw("HIST same")
    
    for x in plotHist[1:]:
      x.Draw("][ same") 
      x.Draw("Hist][ same")
      # x.Draw("Hist same")
    plotHist[0].Draw("AXIS same") 

    self.pad.cd()

  def dPlot(self, s, plotHist):
    self.pad.cd()
    padb = TPad("b", "b", s[0], s[1], s[2], s[3])
    padb.Draw()
    padb.cd()
    #padb.SetGrid()
    padb.SetTopMargin(0.018)
    padb.SetBottomMargin(0.3)

    plotHist[0].Reset()
    plotHist[0].SetTitle("")
    plotHist[0].SetYTitle(self.divTitle)
    plotHist[0].SetMaximum(self.rMax)
    plotHist[0].SetMinimum(self.rMin)
    plotHist[0].GetYaxis().SetTitleSize(0.08)
    plotHist[0].GetYaxis().SetTitleOffset(0.6)
    plotHist[0].GetYaxis().SetLabelSize(0.07)
    plotHist[0].GetXaxis().SetTitleSize(0.1)
    plotHist[0].GetXaxis().SetLabelSize(0.07)
    #plotHist[0].Draw("][")
    #plotHist[0].Draw("HIST")
    plotHist[0].Draw()

    if self.sysDraw:
      self.sysGRDiv()
      for x in reversed(self.sysGRDivList):
        x.Draw("P2") 
    for x in plotHist[1:]:
      x.Draw("same")
      x.Draw("HIST same")
    plotHist[0].Draw("AXIS same")
    self.pad.cd()

  def sysPlot(self, s,m,frame, sysN,sGR):
    self.pad.cd()
    self.padList.append(TPad(sysN, sysN, s[0], s[1], s[2], s[3]))
    self.padList[-1].SetTopMargin(0.018)
    self.padList[-1].SetBottomMargin(m)
    self.padList[-1].SetGrid()
    self.padList[-1].Draw()
    self.padList[-1].cd()
    frame.Reset()
    frame.SetTitle("")
    frame.SetYTitle(sysN)
    frame.SetMaximum(self.rMax)
    frame.SetMinimum(self.rMin)
    if m != 0.018: coe = 1.0-0.17-(0.17*(0.2-0.018))
    else: coe = 1.0
    frame.GetYaxis().SetTitleSize(0.11*coe*1.5)
    frame.GetYaxis().SetTitleOffset(0.37+(1-coe)*coe*0.5)
    frame.GetYaxis().SetLabelSize(0.12*coe)
    frame.GetXaxis().SetTitleSize(0.2)
    frame.GetXaxis().SetTitleOffset(0.8)
    frame.GetXaxis().SetLabelSize(0.1)
    frame.GetXaxis().SetLabelOffset(0.01)
    frame.Draw()
    sGR.Draw("P2")
    frame.Draw("AXIS same")
    self.pad.cd()
    self.padList[-1].Draw()
   
  def allSys(self, allSysList):
    self.allSysList = allSysList   

  def listSysGR(self, list):
    nBin = self.histList[0].GetNbinsX()
    xBin = []
    xEr = []
    yBin = []
    yEr = []
    for x in xrange(nBin):
      xBin.append(self.histList[0].GetBinCenter(x+1))
      xEr.append(self.histList[0].GetBinWidth(x+1)/2.0)
      yBin.append(1.0)
      yEr.append(float(list[1][x]))
      tmpGR = TGraphAsymmErrors(nBin, array("d", xBin), array("d", yBin), array("d", xEr), array("d", xEr), array("d", yEr), array("d", yEr))
      tmpGR.SetFillColor(self.mySysColor[0])
      tmpGR.SetMarkerStyle(1)
      tmpGR.SetMarkerColor(self.mySysColor[0])
      tmpGR.SetLineColor(self.mySysColor[0])
    self.sGR = tmpGR

  def compPlot(self, histList, histName, tag, sysSet, text):
    self.setHistList(histList)
    self.histName = histName
    self.setSysDraw(sysSet)
    self.findMaxMin()
    self.mPlot([0.0, 1.0, 1.0, 0.35], self.histList)
    self.dPlot([0.0, 0.35, 1.0, 0.0], self.divHistList)
    le = TLegend(0.5, 0.94-10.0*0.026, 0.8, 0.95)
    #le = TLegend(0.2, 0.9-6.0*0.026, 0.7, 0.9)
    le.SetHeader(histName[0])
    le.SetTextSize(0.023)
    le.SetTextSize(0.033)
    le.SetFillStyle(0)
    le.SetFillColor(kWhite)
    le.SetBorderSize(0)
    for i, x in enumerate(self.histList):
      le.AddEntry(x, histName[i+1])
    for i, x in enumerate(self.sysGRList):
      le.AddEntry(x, sysSet[2][i])
    le.Draw()
    if text:
      drawLumi(text[0], text[1], text[2])
    storeN = tag+"_"+histList[0].GetName()
    self.pad.SaveAs(storeN+".png")
    self.pad.SaveAs(storeN+".eps")
    self.histTxt(tag) 

  def compPlotWOD(self, histList, histName, tag, sysSet, text):
    self.setHistList(histList)
    self.histName = histName
    self.setSysDraw(sysSet)
    self.findMaxMin()
    self.mPlotSet([0.1, 0.1, 1.8, 0.03, 0.035, 0.0, 1.1])
    self.mPlot([0.0, 1.0, 1.0, 0.0], self.histList)
    #le = TLegend(0.2, 0.9-6.0*0.026, 0.7, 0.9)
    le = TLegend(0.55, 0.94-6.0*0.026, 0.7, 0.9)
    le.SetHeader(histName[0])
    le.SetTextSize(0.021)
    le.SetFillStyle(0)
    le.SetFillColor(kWhite)
    le.SetBorderSize(0)
    for i, x in enumerate(self.histList):
      le.AddEntry(x, histName[i+1])
    for i, x in enumerate(self.sysGRList):
      le.AddEntry(x, sysSet[2][i])
    le.Draw()
    if text:
      drawLumi(text[0], text[1], text[2])
    storeN = tag+"_"+histList[0].GetName()
    self.pad.SaveAs(storeN+".png")
    self.pad.SaveAs(storeN+".eps")

  def sysSPlot(self, histList, tag,text):
    self.setHistList(histList)
    for x in xrange(len(self.allSysList)):
      self.listSysGR(self.allSysList[x])
      self.padH.append(self.divHistList[0].Clone(self.allSysList[x][0]))
      self.padS.append(self.sGR)
      if x == len(self.allSysList)-1:
        self.sysPlot([0.0, 0.92-x*0.14, 1.0, 0.92-((x+1)*0.14+0.14*(0.2-0.018))],0.2, self.padH[-1], self.allSysList[x][0] ,self.padS[-1] )
      else:
        self.sysPlot([0.0, 0.92-x*0.14, 1.0, 0.92-(x+1)*0.14],0.018, self.padH[-1], self.allSysList[x][0] ,self.padS[-1] )
    cutText =TLatex()
    cutText.SetNDC()
    cutText.SetTextFont(42)
    cutText.SetTextSize(0.05)
    cutText.DrawLatex(0.2,0.95,text[0][2])
    store_n = tag+"_"+histList[0].GetName()
    self.pad.SaveAs(store_n+".png")
    self.pad.SaveAs(store_n+".eps")

  def sysSPlot2(self, histList, tag,text):
    self.setHistList(histList)
    for x in xrange(len(self.allSysList)):
      self.listSysGR(self.allSysList[x])
      self.padH.append(self.divHistList[0].Clone(self.allSysList[x][0]))
      self.padS.append(self.sGR)
      if x == len(self.allSysList)-1:
        self.sysPlot([0.0, 0.92-x*0.14, 1.0, 0.92-((x+1)*0.14+0.14*(0.2-0.018))],0.2, self.padH[-1], self.allSysList[x][0] ,self.padS[-1] )
      else:
        self.sysPlot([0.0, 0.92-x*0.14, 1.0, 0.92-(x+1)*0.14],0.018, self.padH[-1], self.allSysList[x][0] ,self.padS[-1] )
    cutText =TLatex()
    cutText.SetNDC()
    cutText.SetTextFont(42)
    cutText.SetTextSize(0.05)
    cutText.DrawLatex(0.2,0.95,text[0][2])
    store_n = tag+"_"+histList[0].GetName()
    self.pad.SaveAs(store_n+".png")
    self.pad.SaveAs(store_n+".eps")


