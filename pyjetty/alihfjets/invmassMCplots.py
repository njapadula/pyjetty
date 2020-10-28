#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as mp
import ROOT
import yaml
import math
ROOT.gROOT.SetBatch(True)

from pyjetty.mputils import treereader

parser = argparse.ArgumentParser(description='D0 inv mass')
parser.add_argument('-f', '--ifile', help='input rootfile', type=str, default=None, required=True)
parser.add_argument('-e', '--energy', help='system energy', type=str, default='5TeV', required=False)
parser.add_argument('-s', '--system', help='collision system', type=str, default='pp', required=False) 
parser.add_argument('-t', '--trailer', help='special trailer', type=str, default='', required=False)
args = parser.parse_args()

def make_var_hist(name, nbins, xarray):
	h = ROOT.TH1F(name, name, nbins, xarray)
	return h

main_dir = "/home/software/users/napadula/"
plot_dir = main_dir + "plots/"
root_dir = main_dir + "rootfiles/"
ofile = "InvMassMC_ptbins_" + args.energy + "_" + args.system + "_" + args.trailer + ".root"
png = args.system + args.energy + args.trailer + ".png"

nbins = 8
pt_low = [3, 4, 5, 6, 7, 8, 10, 12]
pt_high = [4, 5, 6, 7, 8, 10, 12, 15]
pt_mid = [3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13.5]
pt_edge = array('d',[3, 4, 5, 6, 7, 8, 10, 12, 15])

root_filename = root_dir + ofile
rootfile = ROOT.TFile.Open(root_filename, "READ")

branch_names = ['inv_mass', 'pt_cand', 'pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p']
cut_names = ['pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p']
cut_type = ['greater', 'greater', 'less', 'abs', 'less', 'greater']
cut_value = [[0.7, 0.7, 0.03, 0.8, -0.0002, 0.9],
[0.7, 0.7, 0.03, 0.8, -0.0002, 0.9],
[0.7, 0.7, 0.03, 0.8, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.8, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.8, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.9, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.9, -0.00005, 0.85],
[0.6, 0.6, 0.03, 1, 10, 0.8]]

signal = []
bg = []
sob = []
signif = []

hists = []
hnames = ["mean_vs_pt", "sigma_vs_pt", "raw_yield", "background", "subtracted", "SB", "significance"]

h_invmass = []

for p in range(len(pt_low)):
	name = "InvMass_{}_to_{}".format(pt_low[p],pt_high[p])
	h_invmass.insert(p, rootfile.Get(name))
	h_invmass[p].Sumw2()
	h_invmass[p].SetDirectory(0)

	print(name)

h_genD0_pt_yield = rootfile.Get("h_genD0_pt_yield")
h_genD0_pt_yield.Sumw2()
h_genD0_pt_yield.SetDirectory(0)

rootfile.Close()

for n in range(len(hnames)):
	hists.insert(n, make_var_hist(hnames[n], nbins, pt_edge))
	hists[n].Sumw2()

for p in range(nbins):
	#bgfunc = ROOT.TF1("bgfunc", "pol1", 1, 3)
	bgfunc = ROOT.TF1("bgfunc", "expo", 1, 3)
	#p1 = ROOT.TF1("p1", "pol1", 1.65, 1.8)
	p1 = ROOT.TF1("p1", "expo", 1.65, 1.8)
	g1 = ROOT.TF1("g1", "gaus", 1.82, 1.89)
	#f1 = ROOT.TF1("f1", "pol1(0)+gaus(2)", 1.65, 2.1) 
	f1 = ROOT.TF1("f1", "expo(0)+gaus(2)", 1.65, 2.1) 
	g1.SetParLimits(1, 1.864, 1.87)
	g1.SetParLimits(2, 0.005, 0.03)
	f1.SetParLimits(3, 1.864, 1.87)
	f1.SetParLimits(4, 0.005, 0.03)
	par = [0, 0, 0, 0, 0]
	#f1.FixParameter(4, 1.869)

	h_invmass[p].Fit("p1","R")
	par[0:2] = [p1.GetParameter(0), p1.GetParameter(1)]
	h_invmass[p].Fit("g1","R")
	par[2:] = [g1.GetParameter(0), g1.GetParameter(1), g1.GetParameter(2)]
	print(par)
	#f1.FixParameter(3, 1.868)
	f1.SetParameters(par[0], par[1], par[2], 1.868, par[4])
	h_invmass[p].Fit("f1","R")

	polconst = f1.GetParameter(0)
	polslope = f1.GetParameter(1)
	#polquad = f1.GetParameter(2)
	mean = f1.GetParameter(3)
	sigma = f1.GetParameter(4)
	mean_err = f1.GetParError(3)
	sig_err = f1.GetParError(4)

	bgfunc.SetParameters(polconst, polslope)

	bin = hists[0].GetXaxis().FindBin(pt_mid[p])
	hists[0].SetBinContent(bin, mean)
	hists[0].SetBinError(bin, mean_err)

	hists[1].SetBinContent(bin, sigma)
	hists[1].SetBinError(bin, sig_err)

	bin_low = h_invmass[p].GetXaxis().FindBin(mean - 2*sigma)
	bin_high = h_invmass[p].GetXaxis().FindBin(mean + 2*sigma)
	counts = h_invmass[p].Integral(bin_low, bin_high)
	count_err = math.sqrt(counts)
	hists[2].SetBinContent(bin, counts)
	hists[2].SetBinError(bin, count_err)

	lowbin9 = h_invmass[p].GetXaxis().GetBinLowEdge(h_invmass[p].GetXaxis().FindBin(mean - 9*sigma))
	lowbin4 = h_invmass[p].GetXaxis().GetBinUpEdge(h_invmass[p].GetXaxis().FindBin(mean - 4*sigma))
	highbin9 = h_invmass[p].GetXaxis().GetBinUpEdge(h_invmass[p].GetXaxis().FindBin(mean + 9*sigma))
	highbin4 = h_invmass[p].GetXaxis().GetBinLowEdge(h_invmass[p].GetXaxis().FindBin(mean + 4*sigma))
	lowpeak = h_invmass[p].GetXaxis().GetBinLowEdge(h_invmass[p].GetXaxis().FindBin(mean - 2*sigma))
	highpeak = h_invmass[p].GetXaxis().GetBinUpEdge(h_invmass[p].GetXaxis().FindBin(mean + 2*sigma))

	bgcount = h_invmass[p].Integral(h_invmass[p].GetXaxis().FindBin(mean - 9*sigma), h_invmass[p].GetXaxis().FindBin(mean - 4*sigma)) + h_invmass[p].Integral(h_invmass[p].GetXaxis().FindBin(mean + 4*sigma), h_invmass[p].GetXaxis().FindBin(mean + 9*sigma))
	bgcount_err = math.sqrt(bgcount)
	
	bglow = bgfunc.Integral(lowbin9, lowbin4)
	bghigh = bgfunc.Integral(highbin4, highbin9)
	bgpeak = bgfunc.Integral(lowpeak, highpeak)
	alpha = bgpeak/(bglow + bghigh)

	bgcount = bgcount*alpha
	bgcount_err = bgcount_err*alpha

	hists[3].SetBinContent(bin, bgcount)
	hists[3].SetBinError(bin, bgcount_err)

	signal.insert(p, (counts-bgcount)/0.9545)
	bg.insert(p, bgcount)
	sob.insert(p, signal[p]/bg[p])
	signif.insert(p, signal[p]/math.sqrt(signal[p]+bg[p]))

	print("signal: ", signal[p])
	print("bg: ", bgcount)
	print("s/b: ", sob[p])
	print("significance: ", signif[p])



hists[4].Add(hists[2])
hists[4].Add(hists[3], -1)
hists[4].Scale(1/0.9545)
hists[5].Add(hists[4])
hists[5].Divide(hists[3])
rootfile = ROOT.TFile(root_filename, 'UPDATE')
for n in range(len(hnames)):
	hists[n].Write()
rootfile.Close()


#Draw the histograms and save the canvases
t = ROOT.TLatex()
t.SetNDC()
t.SetTextSize(18)
t.SetTextFont(43)
t.SetTextAlign(13)
c = []
for i in range(nbins):
	c.append(ROOT.TCanvas("c"+str(i), "c"+str(i), 700, 700))
	ROOT.gStyle.SetOptStat(0)
	h_invmass[i].SetXTitle("inv mass")
	h_invmass[i].GetXaxis().SetRangeUser(1.7, 2.05)
	h_invmass[i].Draw("E0")
	t.DrawLatex(0.65, 0.86, "D^{0} counts (2#sigma): %.0f" %signal[i])
	t.DrawLatex(0.65, 0.81, "S/B: %.02f" %sob[i])
	t.DrawLatex(0.65, 0.76, "Significance: %.01f" %signif[i])
	t.Draw("same")
	c[i].SaveAs(plot_dir + "InvMass_pt_" + str(pt_low[i]) + "_to_" + str(pt_high[i]) + png)	

cc = ROOT.TCanvas("c2","c2", 700, 700)
h = ROOT.TH1F("h1", "", 16, -0.5, 15.5)
h.GetYaxis().SetRangeUser(1, 100000)
cc.SetLogy()
ROOT.gStyle.SetOptStat(0)
h.SetXTitle("pT")
h.SetYTitle("raw counts")
h.Draw()
hists[2].SetMarkerStyle(20)
hists[2].SetMarkerColor(4)
hists[2].Draw("P same")

colors = [1, 2, 4]
#for i in range(3):
#	hists[i+2].SetMarkerStyle(20)
#	hists[i+2].SetMarkerColor(colors[i])
#	hists[i+2].Draw("P same")

cc.SaveAs(plot_dir + "yields_" + png)

can = ROOT.TCanvas("can","can", 700, 700)
hh = ROOT.TH1F("hh", "Signal over Background", 16, -0.5, 15.5)
hh.GetYaxis().SetRangeUser(0, 4)
ROOT.gStyle.SetOptStat(0)
hh.SetXTitle("pT")
hh.SetYTitle("Signal/Background")
hh.Draw()

hists[5].SetMarkerStyle(21)
hists[5].SetMarkerColor(1)
hists[5].Draw("P same")
can.SaveAs(plot_dir + "SB_" + png) 

yax = ["mean", "sigma"]
ccc = ROOT.TCanvas("c3", "c3", 1000, 500)
ccc.Divide(2,1)
for i in range(2):
	ccc.cd(i+1)
	hists[i].SetXTitle("pT")
	hists[i].SetYTitle(yax[i])
	hists[i].Draw("E0")

ccc.SaveAs(plot_dir + "mean_sigma_" + png)

