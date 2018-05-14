from ROOT import *
from os import mkdir, chdir, system
import copy
from array import array
from plotMaker import *
gROOT.SetBatch()

if __name__ == '__main__':
  root_f_list = ["Unfold_hist_MADGRAPH_iter7.root", "QCD_Pt_TuneCUETP8M1_13TeV_pythia8_hist.root","QCD_powheg_ct10_pythia8_cuetp8m1_13TeV_hist_noPS.root","QCD_powheg_ct10_pythia8_cuetp8m1_13TeV_hist.root", "QCD_HT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_hist.root"]
  tr_list = []
  for x in root_f_list:
    tr_list.append(TFile(x))
  chdir("./sysDraw")
  plotList = [["del_r23", "dr_lpt3", "#DeltaR_{23}"],["del_r23", "dr_hpt3", "#DeltaR_{23}"], ["jet3_pt_jet2_pt", "sdr_pt3", "p_{T3}/p_{T2}"], ["jet3_pt_jet2_pt", "ldr_pt3", "p_{T3}/p_{T2}"], ["del_r23", "dr_lpt3", "#DeltaR_{23}"],["del_r23", "dr_hpt3", "#DeltaR_{23}"], ["jet3_pt_jet2_pt", "sdr_pt3", "p_{T3}/p_{T2}"], ["jet3_pt_jet2_pt", "ldr_pt3", "p_{T3}/p_{T2}"]]
  for v in plotList:
    h_list = []
    nom_key = "non_eta_high_pt_"+v[1]+"_nom_"+v[0]
    rivet_key = "non_eta_high_pt_"+v[1]+"_rivet_"+v[0]
    data = copy.deepcopy(tr_list[0].Get(nom_key.replace("_nom_", "_roounfold_res_")+"_rd"))
    gen_pyt =copy.deepcopy(tr_list[1].Get(nom_key.replace("_nom_","_genJet_")))
    gen_mad = copy.deepcopy(tr_list[4].Get(nom_key.replace("_nom_","_genJet_")))
 
    mc1 = copy.deepcopy(setError(tr_list[2].Get(rivet_key)))
    mc2 = copy.deepcopy(setError(tr_list[3].Get(rivet_key)))
    h_list.append(data)
    h_list.append(gen_pyt)
    #h_list.append(mc1)
    h_list.append(mc2)
    h_list.append(gen_mad)
    #h_list.append(mad_noPS)
    h_list[0].SetXTitle(v[2])
    h_list[0].SetYTitle("#frac{1}{#sigma} #frac{d#sigma}{d%s}"%v[2])
    #normHist(h_list[0])
    #normHist(h_list[2])

    #print [h_list[2].GetBinContent(x+1)/h_list[0].GetBinContent(x+1) for x in xrange(h_list[0].GetNbinsX())]

    sys_r = open("./sysTxt/"+nom_key.replace("nom", "sys_tot")+".txt","r")
    sys_sta = sys_r.readline().split(":")[1].split(",")
    sys_sys = sys_r.readline().split(":")[1].split(",")
    sys_tot = sys_r.readline().split(":")[1].split(",")
    sysSet = [True, [sys_sys, sys_tot], ["sys.", "sys. + sta."]]
    #sysSet = [False, [sys_sys, sys_tot], ["sys.", "sys. + sta."]]
    text = [["Preliminary", "2.29 fb^{-1} (13 TeV)"]+drawCut(v[1]), 0.2, 0.9]
    myPlot = plotMaker()
    myPlot.setRange(1.3,0.7)
    #myPlot.compPlot(h_list, ["","Unfolded Data","PYTHIA8 (CUETP8M1)","POWHEG + PYTHIA8 (CUETP8M1)","POWHEG + PYTHIA8 without PS", "MADGRAPH + PYTHIA8 (CUETP8M1)", "MADGRAPH + PYTHIA8 without PS"], "SYS_DRAW", sysSet, text)
    #myPlot.compPlot(h_list, ["","Unfolded Data","PYTHIA8 (CUETP8M1)","POWHEG + PYTHIA8 (CUETP8M1)","POWHEG + PYTHIA8 without PS", "MADGRAPH + PYTHIA8 (CUETP8M1)"], "SYS_DRAW", sysSet, text)
    #myPlot.compPlot(h_list, ["","Unfolded Data","PYTHIA (LO2jets+PS)","POWHEG (NLO2jets)","POWHEG (NLO2jets+PS)", "MADGRAPH (LO4jets+PS)"], "SYS_DRAW_13TeV", sysSet, text)
    myPlot.compPlot(h_list, ["","Data","PYTHIA (LO2jets+PS)","POWHEG (NLO2jets+PS)", "MADGRAPH (LO4jets+PS)"], "SYS_DRAW_13TeV", sysSet, text)
    #myPlot.compPlot(h_list, ["","Unfolded Data","PYTHIA (LO2jets+PS)", "MADGRAPH (LO4jets+PS)"], "SYS_DRAW_13TeV", sysSet, text)
