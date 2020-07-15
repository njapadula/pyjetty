#! /usr/bin/env python

import sys
import os
import argparse
from array import *
import numpy as mp
import ROOT
import yaml
import math

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
ofile = "InvMass" + args.system + args.energy + args.trailer + ".root"
png = args.system + args.energy + args.trailer + ".png"

nbins = 4
pt_low = [2, 4, 6, 8]
pt_high = [4, 6, 8, 50]
pt_mid = [3, 5, 7, 29]
pt_edge = array('d',[2, 4, 6, 8, 50])

root_filename = root_dir + ofile
rootfile = ROOT.TFile(root_filename, 'RECREATE')
rootfile.Close()

bgfunc = ROOT.TF1("bgfunc", "pol1", 1, 3)
p1 = ROOT.TF1("p1", "pol1", 1.65, 1.8)
g1 = ROOT.TF1("g1", "gaus", 1.84, 1.8)
f1 = ROOT.TF1("f1", "pol1(0)+gaus(2)", 1.65, 2.1) 
f1.SetParLimits(2, 0, 100000000)
f1.SetParLimits(3, 1.82, 1.9)
f1.SetParLimits(4, 0.01, 0.2)
par = [0, 0, 0, 0, 0]

hists = []
hnames = ["mean_vs_pt", "sigma_vs_pt", "raw_yield", "background", "subtracted", "SB"]

for n in range(len(hnames)):
	hists.insert(n, make_var_hist(hnames[n], nbins, pt_edge))
	hists[n].Sumw2()

h = []

for p in range(len(pt_low)):
	print(p)
	name = "InvMass_{}_to_{}".format(pt_low[p],pt_high[p])
	h.insert(p, ROOT.TH1F(name, name, 45, 1.65, 2.1))
	h[p].Sumw2()

	print(name)

	tr = treereader.RTreeReader(tree_name='d0', branches = ['inv_mass','pt_cand','pt_prong0','pt_prong1','dca','cos_t_star','imp_par_prong0','imp_par_prong1', 'imp_par_prod', 'cos_p'], file_name = args.ifile)


	for i in range(tr.tree.GetEntries()):
		tr.tree.GetEntry(i)
		if (tr.pt_cand[0] >= pt_low[p]) and (tr.pt_cand[0] < pt_high[p]) and (tr.inv_mass[0] > 0):
			if(args.system!="PbPb"):
				h[p].Fill(tr.inv_mass[0])
			elif(tr.pt_prong0[0] > 0.6 and tr.pt_prong1[0] > 0.6 and tr.dca[0] < 0.03 and abs(tr.cos_t_star[0]) < 0.8 and abs(tr.imp_par_prong0[0]) < 0.1 and abs(tr.imp_par_prong1[0]) < 0.1 and tr.imp_par_prod[0]<-0.0001 and tr.cos_p[0]>0.9):
				h[p].Fill(tr.inv_mass[0])
			

	h[p].Fit("p1","R")
	par[0:2] = [p1.GetParameter(0), p1.GetParameter(1)]
	h[p].Fit("g1","R")
	par[2:5] = [g1.GetParameter(0), g1.GetParameter(1), g1.GetParameter(2)]
	print(par)
	f1.SetParameters(par[0], par[1], par[2], par[3], par[4])
	h[p].Fit("f1","R")
	# Write to a root file
	rootfile = ROOT.TFile(root_filename, 'UPDATE')
	h[p].Write()
	rootfile.Close()

	polconst = f1.GetParameter(0)
	polslope = f1.GetParameter(1)
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

	bin_low = h[p].GetXaxis().FindBin(mean - 2*sigma)
	bin_high = h[p].GetXaxis().FindBin(mean + 2*sigma)
	counts = h[p].Integral(bin_low, bin_high)
	count_err = math.sqrt(counts)
	hists[2].SetBinContent(bin, counts)
	hists[2].SetBinError(bin, count_err)

	lowbin9 = h[p].GetXaxis().GetBinLowEdge(h[p].GetXaxis().FindBin(mean - 9*sigma))
	lowbin4 = h[p].GetXaxis().GetBinUpEdge(h[p].GetXaxis().FindBin(mean - 4*sigma))
	highbin9 = h[p].GetXaxis().GetBinUpEdge(h[p].GetXaxis().FindBin(mean + 9*sigma))
	highbin4 = h[p].GetXaxis().GetBinLowEdge(h[p].GetXaxis().FindBin(mean + 4*sigma))
	lowpeak = h[p].GetXaxis().GetBinLowEdge(h[p].GetXaxis().FindBin(mean - 2*sigma))
	highpeak = h[p].GetXaxis().GetBinUpEdge(h[p].GetXaxis().FindBin(mean + 2*sigma))

	bgcount = h[p].Integral(h[p].GetXaxis().FindBin(mean - 9*sigma), h[p].GetXaxis().FindBin(mean - 4*sigma)) + h[p].Integral(h[p].GetXaxis().FindBin(mean + 4*sigma), h[p].GetXaxis().FindBin(mean + 9*sigma))
	bgcount_err = math.sqrt(bgcount)
	
	bglow = bgfunc.Integral(lowbin9, lowbin4)
	bghigh = bgfunc.Integral(highbin4, highbin9)
	bgpeak = bgfunc.Integral(lowpeak, highpeak)
	alpha = bgpeak/(bglow + bghigh)

	bgcount = bgcount*alpha
	bgcount_err = bgcount_err*alpha

	hists[3].SetBinContent(bin, bgcount)
	hists[3].SetBinError(bin, bgcount_err)


#hists[2].Scale(1/0.9545)
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
c = ROOT.TCanvas("c1", "c1", 700, 700)
c.Divide(2,2)
for i in range(4):
	c.cd(i+1)
	h[i].SetXTitle("inv mass")
	h[i].Draw("E0")

c.SaveAs(plot_dir + "InvMass" + png)	

cc = ROOT.TCanvas("c2","c2", 700, 700)
h = ROOT.TH1F("h1", "h1", 51, -0.5, 50.5)
h.GetYaxis().SetRangeUser(1, 100000)
cc.SetLogy()
ROOT.gStyle.SetOptStat(0)
h.SetXTitle("pT")
h.SetYTitle("raw counts")
h.Draw()

colors = [1, 2, 4]
for i in range(3):
	hists[i+2].SetMarkerStyle(20)
	hists[i+2].SetMarkerColor(colors[i])
	hists[i+2].Draw("P same")

cc.SaveAs(plot_dir + "yields" + png)

can = ROOT.TCanvas("can","can", 700, 700)
hh = ROOT.TH1F("hh", "Signal over Background", 51, -0.5, 50.5)
hh.GetYaxis().SetRangeUser(0, 0.5)
ROOT.gStyle.SetOptStat(0)
hh.SetXTitle("pT")
hh.SetYTitle("Signal/Background")
hh.Draw()

hists[5].SetMarkerStyle(21)
hists[5].SetMarkerColor(1)
hists[5].Draw("P same")
can.SaveAs(plot_dir + "SB" + png) 

yax = ["mean", "sigma"]
ccc = ROOT.TCanvas("c3", "c3", 1000, 500)
ccc.Divide(2,1)
for i in range(2):
	ccc.cd(i+1)
	hists[i].SetXTitle("pT")
	hists[i].SetYTitle(yax[i])
	hists[i].Draw("E0")

ccc.SaveAs(plot_dir + "mean_sigma" + png)

