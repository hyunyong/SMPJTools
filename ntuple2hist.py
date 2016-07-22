#!/usr/bin/env python
import ROOT, os, sys, copy
from array import array

ROOT.gROOT.SetBatch()
pi = ROOT.TMath.Pi()

in_f = sys.argv[1]

out_f = in_f.replace(".root", "_hist.root").split("/")[-1]
in_rf = ROOT.TFile(in_f)
mc = True

if not in_f.startswith("QCD"):
  mc = False

def getMCW13TeV(in_f):
  in_f = in_f.split("/")[-1]
  luminosity = 2294.7
  mad_cross = [2.110e+08, 1717000, 351300, 31630, 6802, 1206, 120.4, 25.25]
  mad_ent = [82110945, 18822308, 17040053, 19701790, 15575940, 5085104, 3952170, 1981228]
  mad_name = ["100to200","200to300","300to500", "500to700", "700to1000", "1000to1500","1500to2000", "2000toInf"]
  pyt_cross = [61018300000, 5887580000, 1837410000, 140932000, 19204300, 2762530, 471100, 117276, 7823, 648.2, 186.9, 32.293, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445] 
  pyt_ent = [6991536, 6842800, 38387715, 9808025, 9775360, 6953590, 6848223, 6918748, 5968960, 3977770, 3979884, 3973224, 2953982, 395725, 393760, 398452, 391108]
  pyt_name = ["5to10","10to15","15to30","30to50","50to80","80to120","120to170","170to300","300to470","470to600","600to800","800to1000","1000to1400","1400to1800","1800to2400","2400to3200","3200toInf"]

  if in_f.find("QCD") < 0:
    print "It is Not MC File!"
    return 1.0
  if in_f.find("madgraphMLM-pythia8") > 0:
    print "It is MADGRAPH"
    index = -1
    for i,x in enumerate(mad_name):
      if in_f.startswith("QCD_HT"+x): index = i
    if index < 0:
      print "Can Not Finde File Index in MADGRAPH!"
      return 0
    return luminosity*mad_cross[index]/mad_ent[index]
  else:
    print "It is PYTHIA8"
    index = -1
    for i,x in enumerate(pyt_name):
      if in_f.startswith("QCD_Pt_"+x): index = i
    if index < 0:
      print "Can Not Finde File Index in PYTHIA8!"
      return 0
    return luminosity*pyt_cross[index]/pyt_ent[index]

def hist_maker(name, title, bin_set, x_name, y_name, tr, br, w,s = 1.0):
  if bin_set[2] == 2500:
    bin = []
    for x in xrange(15):
      bin.append(2500.0/2.0*float(x)/15.0)
    bin.append(1500)
    bin.append(2000)
    bin.append(2500)
    pt_bin = array('d', bin)
    hist = ROOT.TH1F(name, title, len(pt_bin)-1, pt_bin)
  elif len(bin_set)>3:
    bin_a = array('d', bin_set)
    hist = ROOT.TH1F(name, title, len(bin_set)-1, bin_a)
  else:
    hist = ROOT.TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
  hist.GetXaxis().SetTitle(x_name)
  hist.GetYaxis().SetTitle(y_name)
  hist.Sumw2()
  print br, w
  #aaaaa =str(br)+">>"++str(name)
  #bbbbb =str(w)
  tr.Project(name, br, w)
  #tr.Draw("%s>>%s" % (br, name), "1.0*%s" % w, "goff")
  #tr.Draw(aaaaa, bbbbb, "goff")
  hist.Scale(s)
  return hist

def prof_maker(name, title, bin_set, x_name, y_name, tr, br, w):
  prof = ROOT.TProfile(name, title, bin_set[0], bin_set[1], bin_set[3])
  prof.GetXaxis().SetTitle(x_name)
  prof.GetYaxis().SetTitle(y_name)
  prof.Sumw2()
  tr.Project(name, br, w)
  return prof

def hist2_maker(name, title, bin_set, x_name, y_name, tr, br_x, br_y, w, s = 1.0):
  if len(bin_set[0]) > 3:
    bin_x =array('d', bin_set[0])
    bin_y =array('d', bin_set[1])
    hist2 = ROOT.TH2D(name, title, len(bin_x)-1, bin_x, len(bin_y)-1, bin_y)
  else:
    hist2 = ROOT.TH2D(name, title, bin_set[0][0], bin_set[0][1], bin_set[0][2], bin_set[1][0], bin_set[1][1], bin_set[1][2])
  hist2.GetXaxis().SetTitle(x_name)
  hist2.GetYaxis().SetTitle(y_name)
  hist2.Sumw2()
  tr.Project(name, br_y+":"+br_x, w)
  hist2.Scale(s)
  return hist2

## event selection
r_mass_cut = "(raw_mass > 220)*"
dp12_cut = "(abs(abs(del_phi12)-{})<1.0)*".format(pi)
eta_cut = "(abs(jet1_eta)<2.5)*"
met_cut = "(metSig<0.3)*"
pt3 = "(jet3_pt/jet2_pt < 0.9)*"
d_cut = r_mass_cut+dp12_cut+eta_cut+met_cut+pt3

## event classification
n_eta = "(abs(jet2_eta)<2.5)*(jet3_pt>30)*"#(jet3_pt>30)*(jet3_p<100)*"#(abs(jet3_eta)>0.08)*"
l_eta = "(abs(jet2_eta)>0.0)*(abs(jet2_eta)<0.8)*(jet3_pt>30)*"#(jet3_pt>30)*(jet3_p<100)*"#(abs(jet3_eta)>0.08)*"
m_eta = "(abs(jet2_eta)>0.8)*(abs(jet2_eta)<1.5)*(jet3_pt>30)*"#(jet3_pt>30)*(jet3_p<200)*"#(abs(jet3_eta)>0.9)*"
h_eta = "(abs(jet2_eta)>1.5)*(abs(jet2_eta)<2.5)*(jet3_pt>30)*"#(jet3_pt>30)*(jet3_p<300)*"#(abs(jet3_eta)>1.92)*"

h_pt = "(jet1_pt>510)*(jet1_pt<2500)*(hlt_450_pass == 1)*"
m_pt = "(jet1_pt>400)*(jet1_pt<500)*(hlt_320_pass == 1)*"
l_pt = "(jet1_pt>170)*(jet1_pt<350)*(hlt_140_pass == 1)*"

dr_cut = "(del_r23>0.4)*(del_r23<1.5)*"
dr_small = "(del_r23>0.4)*(del_r23<1.0)*"
dr_large = "(del_r23>1.0)*(del_r23<1.5)*"

pt3_non = "(jet3_pt/jet2_pt < 0.9)*"
pt3_low = "(jet3_pt/jet2_pt < 0.3)*"
pt3_high = "(jet3_pt/jet2_pt > 0.6)*"

eta_bin = ["non_eta", "low_eta", "medium_eta", "high_eta"]
eta_bin_cut = [n_eta, l_eta, m_eta, h_eta]
pt_bin = ["low_pt", "medium_pt", "high_pt"]
pt_bin_cut = [l_pt, m_pt, h_pt]
dr_bin = ["dr", "sdr", "ldr"]
dr_bin_cut = [dr_cut, dr_small, dr_large]
pt3_bin = ["non", "lpt3", "hpt3"]
pt3_bin_cut = [pt3_non, pt3_low, pt3_high]

beta_l = ["beta23", "del_eta23", "del_phi23", "del_r23", "jet3_pt/jet2_pt"]
gbeta_l = ["gbeta", "gdel_eta", "gdel_phi", "gdel_r", "gjet3_pt/gjet2_pt"]
gen_beta_l = ["gen_beta23", "gen_del_eta23", "gen_del_phi23", "gen_del_r23", "gen_jet3_pt/gen_jet2_pt"]
beta_bin = [[18, 0, pi], [30, -3, 3], [30, -3, 3], [20, 0.0, 2], [0.0,0.1,0.2,0.3,0.4,0.6,0.9,1.0]]

beta_l = ["del_r23", "jet3_pt/jet2_pt"]
gbeta_l = ["gdel_r", "gjet3_pt/gjet2_pt"]
gen_beta_l = ["gen_del_r23", "gen_jet3_pt/gen_jet2_pt"]
beta_bin = [[20, 0.0, 2], [0.0,0.1,0.2,0.3,0.4,0.6,0.9,1.0]]

jet_l = ["pt", "eta", "phi"]
jet_bin = [[30,0,2500],[30,-3,3],[30,-pi,pi]]
ev_l = ["njet", "metSig", "nvtx", "raw_mass", "del_phi12", "pileupWeight"]  
ev_bin = [[30, 0, 30],[10, 0, 1], [50, 0, 50],[50,0,5000],[30,-pi,pi], [20,0,2]]

cut_l = []
cut_name = []

mc_ss = getMCW13TeV(in_f)
mc_s = 1.0
for i, eta in enumerate(eta_bin_cut):
  for j, pt in enumerate(pt_bin_cut):
    for k, dr in enumerate(dr_bin_cut):
      for l, pt3 in enumerate(pt3_bin_cut):
        cut_l.append(d_cut+eta+pt+dr+pt3+"(1.0)")
        cut_name.append(eta_bin[i]+"_"+pt_bin[j]+"_"+dr_bin[k]+"_"+pt3_bin[l])
if mc:
    sys_e = ["nom", "jer","jer_u", "jer_d", "jar", "pu_u", "pu_d"]
    sys_e = ["nom"]
else:
    sys_e = ["nom", "jes_u", "jes_d"]    
hist_l = []

for sys in sys_e:
  if mc:
    e_w = "({}*pileupWeight)".format(mc_ss)
    if sys.startswith("pu"):
      if sys == "pu_u":
        e_w = "(1.0*pileupWeight_up)"
      else:
        e_w = "(1.0*pileupWeight_dn)"
  else:
    e_w = "(1.0)"
  for i, cut in enumerate(cut_l):
    if sys.startswith("pu"):
      tr = in_rf.Get("cc/nom").Clone(cut_name[i])
    else:
      tr = in_rf.Get("cc/%s"%sys).Clone(cut_name[i])
    #tr.Draw(">>eventList", cut)
    #el = ROOT.gDirectory.Get("eventList")
    #tr.SetEventList(el)
 
    for bi, beta_loop in enumerate(beta_l):
      name = cut_name[i]+"_"+sys+"_"+beta_loop
      name = name.replace("/","_")
      title = name
      bin_set = beta_bin[bi]
      x_name = beta_loop
      y_name = "count"
      br = beta_loop
      hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w, mc_s)))
    for ji in xrange(3):
      for jii, jet_loop in enumerate(jet_l):
        name = cut_name[i]+"_"+sys+"_jet%d_"%(ji+1)+jet_loop
        title = name
        bin_set = jet_bin[jii]
        x_name = jet_loop
        y_name = "count"
        br = "jet%d_"%(ji+1)+jet_loop
        hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w,mc_s)))
    for ei, ev_loop in enumerate(ev_l):
      name = cut_name[i]+"_"+sys+"_"+ev_loop
      title = name
      bin_set = ev_bin[ei]
      x_name = ev_loop
      y_name = "count"
      br = ev_loop
      hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w, mc_s)))
    # unfolding part
    if (sys == "nom" and mc):
      for bi, beta_loop in enumerate(beta_l):

        print beta_loop
        name = cut_name[i]+"_gen_"+beta_loop
        name = name.replace("/","_")
        title = name
        bin_set = beta_bin[bi]
        x_name = gbeta_l[bi]
        y_name = "count"
        br = gbeta_l[bi]
        hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w, mc_s)))

        name = cut_name[i]+"_fake_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist_maker(name, title, beta_bin[bi], beta_l[bi], "count", tr, beta_l[bi], e_w+"*((gdel_r1>0.4)+(gdel_r2>0.4)+(gdel_r3>0.4)+(%s==-99))"%beta_l[bi], mc_s)))
        name = cut_name[i]+"_miss_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist_maker(name, title, beta_bin[bi], beta_l[bi], "count", tr, gbeta_l[bi], e_w+"*((gdel_r1>0.4)+(gdel_r2>0.4)+(gdel_r3>0.4)+(%s==-99))"%gbeta_l[bi], mc_s)))

        name = cut_name[i]+"_resM_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist2_maker(name, title, [bin_set,bin_set],"RECO", "GEN", tr, beta_loop, gbeta_l[bi], e_w+"*(gdel_r1<0.4)*(gdel_r2<0.4)*(gdel_r3<0.4)", mc_s)))

        name = cut_name[i]+"_ct_preco_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist_maker(name, title, beta_bin[bi], beta_l[bi], "count", tr, beta_l[bi], e_w+"*(Entry$%2==0)", mc_s)))
        name = cut_name[i]+"_ct_treco_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist_maker(name, title, beta_bin[bi], beta_l[bi], "count", tr, beta_l[bi], e_w+"*(Entry$%2==1)", mc_s)))

        name = cut_name[i]+"_ct_pgen_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist_maker(name, title, beta_bin[bi], beta_l[bi], "count", tr, gbeta_l[bi], e_w+"*(Entry$%2==0)", mc_s)))
        name = cut_name[i]+"_ct_tgen_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist_maker(name, title, beta_bin[bi], beta_l[bi], "count", tr, gbeta_l[bi], e_w+"*(Entry$%2==1)", mc_s)))
  
        name = cut_name[i]+"_ct_resM_"+beta_loop
        name = name.replace("/","_")
        title = name
        hist_l.append(copy.deepcopy(hist2_maker(name, title, [bin_set,bin_set],"RECO", "GEN", tr, beta_loop, gbeta_l[bi], e_w+"*(gdel_r1<0.4)*(gdel_r2<0.4)*(gdel_r3<0.4)*(Entry$%2==0)", mc_s)))
        if beta_loop == "jet3_pt/jet2_pt":
          name = cut_name[i]+"_pt_swp_"+beta_loop
          name = name.replace("/","_")
          title = name
          bin_set = beta_bin[bi]
          x_name = beta_loop
          y_name = "count"
          br = beta_loop
          hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, "(gjet3_pt > gjet2_pt)", mc_s)))
          hist_l.append(copy.deepcopy(hist2_maker(name.replace("_pt_swp_","_pt3_bin_"), title, [bin_set, [10,30,500]], x_name, "Jet_{3} p_{T}", tr, beta_loop, "jet3_pt", e_w, mc_s)))
          hist_l.append(copy.deepcopy(hist2_maker(name.replace("_pt_swp_","_pt2_bin_"), title, [bin_set, [10,30,500]], x_name, "Jet_{2} p_{T}", tr, beta_loop, "jet2_pt", e_w, mc_s)))
          hist_l.append(copy.deepcopy(hist2_maker(name.replace("_pt_swp_","_pt3_bin_swp_"), title, [bin_set, [10,30,500]], x_name, "Jet_{3} p_{T}", tr, beta_loop, "jet3_pt", "(gjet3_pt > gjet2_pt)", mc_s)))
          hist_l.append(copy.deepcopy(hist2_maker(name.replace("_pt_swp_","_pt2_bin_swp_"), title, [bin_set, [10,30,500]], x_name, "Jet_{2} p_{T}", tr, beta_loop, "jet2_pt", "(gjet3_pt > gjet2_pt)", mc_s)))
    del tr
    del el

# gen jet 
if mc:
  ## event selection
  r_mass_cut = "(gen_raw_mass > 220)*"
  dp12_cut = "(abs(abs(gen_del_phi12)-%f)<1.0)*"%pi
  eta_cut = "(abs(gen_jet1_eta)<2.5)*"
  pt3 = "(gen_jet3_pt/gen_jet2_pt < 0.9)*"
  d_cut = r_mass_cut+dp12_cut+eta_cut+pt3
  
  ## event classification
  l_eta = "(abs(gen_jet2_eta)>0.0)*(abs(gen_jet2_eta)<0.8)*(gen_jet3_pt>30)*"
  m_eta = "(abs(gen_jet2_eta)>0.8)*(abs(gen_jet2_eta)<1.5)*(gen_jet3_pt>30)*"
  h_eta = "(abs(gen_jet2_eta)>1.5)*(abs(gen_jet2_eta)<2.5)*(gen_jet3_pt>30)*"
  eta_bin_cut = [l_eta, m_eta, h_eta]

  h_pt = "(gen_jet1_pt>510)*(gen_jet1_pt<2500)*"
  m_pt = "(gen_jet1_pt>400)*(gen_jet1_pt<500)*"
  l_pt = "(gen_jet1_pt>170)*(gen_jet1_pt<350)*"
  pt_bin_cut = [l_pt, m_pt, h_pt]

  dr_cut = "(gen_del_r23>0.4)*(gen_del_r23<1.5)*"
  dr_small = "(gen_del_r23>0.4)*(gen_del_r23<1.0)*"
  dr_large = "(gen_del_r23>1.0)*(gen_del_r23<1.5)*"
  dr_bin_cut = [dr_cut, dr_small, dr_large]

  pt3_non = "(gen_jet3_pt/gen_jet2_pt < 0.9)*"
  pt3_low = "(gen_jet3_pt/gen_jet2_pt < 0.3)*"
  pt3_high = "(gen_jet3_pt/gen_jet2_pt > 0.6)*"
  pt3_bin_cut = [pt3_non, pt3_low, pt3_high]

  
  cut_l = []
  cut_name = []
  
  for i, eta in enumerate(eta_bin_cut):
    for j, pt in enumerate(pt_bin_cut):
      for k, dr in enumerate(dr_bin_cut):
        for l, pt3 in enumerate(pt3_bin_cut):
          cut_l.append(d_cut+eta+pt+dr+pt3+"(1.0)")
          cut_name.append(eta_bin[i]+"_"+pt_bin[j]+"_"+dr_bin[k]+"_"+pt3_bin[l])
  print cut_l  
  for i, cut in enumerate(cut_l):
    tr = in_rf.Get("cc/nom").Clone(cut_name[i])
    tr.Draw(">>eventList", cut)
    el = ROOT.gDirectory.Get("eventList")
    tr.SetEventList(el)
    e_w = "(1.0)"
    for bi, beta_loop in enumerate(beta_l):
      name = cut_name[i]+"_genJet_"+beta_loop
      name = name.replace("/","_")
      title = "Gen Jet "+name
      bin_set = beta_bin[bi]
      x_name = beta_loop
      y_name = "count"
      br = gen_beta_l[bi]
      hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w, mc_s)))
    for ji in xrange(3):
      for jii, jet_loop in enumerate(jet_l):
        name = cut_name[i]+"_genJet_jet%d_"%(ji+1)+jet_loop
        title = "Gen Jet "+name
        bin_set = jet_bin[jii]
        x_name = jet_loop
        y_name = "count"
        br = "gen_jet%d_"%(ji+1)+jet_loop
        hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w, mc_s)))
    del tr
    #del el

 


# out file write
out_rf = ROOT.TFile(out_f, "RECREATE")
for x in hist_l:
  print x.GetName(), x.GetEntries(), x.GetEffectiveEntries()
  x.Write()
out_rf.Write()
out_rf.Close()
in_rf.Close()
