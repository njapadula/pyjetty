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
parser.add_argument('-t', '--trailer', help='information trailer', type=str, default='', required=False)
args = parser.parse_args()

pt_low = 3
pt_high = 100

main_dir = "/home/software/users/napadula/"
root_dir = main_dir + "rootfiles/"
plot_dir = main_dir + "plots/"
ofile = "InvMass_pt_" + str(pt_low) + "_" + args.system + args.energy + "_" + args.trailer + ".root"
png = ".png"

root_filename = root_dir + ofile
rootfile = ROOT.TFile(root_filename, 'RECREATE')
rootfile.Close()

cut_h = []
#branch_names = ['inv_mass', 'pt_cand', 'pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prong0', 'imp_par_prong1', 'imp_par_prod', 'cos_p']
#cut_names = ['pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prong0', 'imp_par_prong1', 'imp_par_prod', 'cos_p', 'allcuts']
#cut_type = ['greater', 'greater', 'less', 'abs', 'abs', 'abs', 'less', 'greater', 'all']
#cut_value = [0.6, 0.6, 0.03, 0.8, 0.1, 0.1, -0.0001, 0.9, 0]	
#how_many_cuts = 8
branch_names = ['inv_mass', 'pt_cand', 'pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p']
cut_names = ['pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p', 'allcuts']
cut_type = ['greater', 'greater', 'less', 'abs', 'less', 'greater', 'all']
cut_value = [0.6, 0.6, 0.03, 0.8, -0.0001, 0.9, 0]	
how_many_cuts = 6

cut_string = '(pt_cand)>' + str(pt_low) + ', '

for n in range(how_many_cuts):
	abs_string = ''
	less_string = '<'
	if cut_type[n]=='abs':
		abs_string = 'abs' 
	if cut_type[n]=='greater':
		less_string = '>'
	cut_string += abs_string + '(' + cut_names[n] + ')' + less_string + str(cut_value[n]) + ', '
print(cut_string)

f = ROOT.TFile(args.ifile, 'READ')
d0 = f.Get("d0")

for n in range(len(cut_names)-1):
	this_cut = cut_names[n] + ' >> hnew'
	d0.Draw(this_cut, "pt_cand>3", "goff")
	hnew = ROOT.gROOT.FindObject('hnew')
	hnew.SetDirectory(0)
	cut_h.insert(n, ROOT.TH1F())
	cut_h[n] = hnew.Clone()
	cut_h[n].SetDirectory(0)
	cut_h[n].SetName(cut_names[n]+"_dist")

f.Close()

tr = treereader.RTreeReader(tree_name='d0', branches = branch_names, file_name = args.ifile)

def make_cut(d0tree, cut, type, value):
	if type == 'greater':
		if getattr(d0tree, cut)[0] > value:
			return True
	if type == 'less':
		if getattr(d0tree, cut)[0] < value:
			return True
	if type == 'abs':
		if abs(getattr(d0tree, cut)[0]) < value:
			return True

	return False

#Fill histograms
invmass_h = ROOT.TH1F("InvMassAbovepT"+str(pt_low), cut_string, 45, 1.65, 2.10)
invmass_h.Sumw2()	
for i in range(tr.tree.GetEntries()):
	tr.tree.GetEntry(i)
	if(i % 100000 == 0):
		print(i)
	if((tr.pt_cand[0] >= pt_low) and (tr.pt_cand[0] < pt_high) and (tr.inv_mass[0] > 0)):
		flag = 0
		for n in range(how_many_cuts):
			if make_cut(tr, cut_names[n], cut_type[n], cut_value[n]):
				flag = flag + 1		
		if flag == how_many_cuts: 
			invmass_h.Fill(tr.inv_mass[0])

bgfunc = ROOT.TF1("bgfunc", "expo", 1, 3)
p1 = ROOT.TF1("p1", "expo", 1.65, 1.8)
g1 = ROOT.TF1("g1", "gaus", 1.82, 1.89)
f1 = ROOT.TF1("f1", "expo(0)+gaus(2)", 1.65, 2.1)
f1.SetParLimits(3, 1.864, 1.872)
f1.SetParLimits(4, 0.01, 0.03)
par = [0, 0, 0, 0, 0]

print("entries: ", invmass_h.GetEntries())

if (invmass_h.GetEntries()<15000):
	print("Fitting with pol1")
	p1.FixParameter(2, 0)
	f1.FixParameter(2, 0)	

invmass_h.Fit("p1","R")
par[0:2] = [p1.GetParameter(0), p1.GetParameter(1)]
invmass_h.Fit("g1","R")
par[2:] = [g1.GetParameter(0), g1.GetParameter(1), g1.GetParameter(2)]
f1.SetParameters(par[0], par[1], par[2], par[3], par[4])
invmass_h.Fit("f1","R")

polconst = f1.GetParameter(0)
polslope = f1.GetParameter(1)
mean = f1.GetParameter(3)
sigma = f1.GetParameter(4)
mean_err = f1.GetParError(3)
sig_err = f1.GetParError(4)

#bgfunc.SetParameters(polconst, polslope, polquad)
bgfunc.SetParameters(polconst, polslope)

bin_low = invmass_h.GetXaxis().FindBin(mean - 2*sigma)
bin_high = invmass_h.GetXaxis().FindBin(mean + 2*sigma)
counts = invmass_h.Integral(bin_low, bin_high)
count_err = math.sqrt(counts)
print("counts: ", counts)

lowbin9 = invmass_h.GetXaxis().GetBinLowEdge(invmass_h.GetXaxis().FindBin(mean - 9*sigma))
lowbin4 = invmass_h.GetXaxis().GetBinUpEdge(invmass_h.GetXaxis().FindBin(mean - 4*sigma))
highbin9 = invmass_h.GetXaxis().GetBinUpEdge(invmass_h.GetXaxis().FindBin(mean + 9*sigma))
highbin4 = invmass_h.GetXaxis().GetBinLowEdge(invmass_h.GetXaxis().FindBin(mean + 4*sigma))
lowpeak = invmass_h.GetXaxis().GetBinLowEdge(invmass_h.GetXaxis().FindBin(mean - 2*sigma))
highpeak = invmass_h.GetXaxis().GetBinUpEdge(invmass_h.GetXaxis().FindBin(mean + 2*sigma))

bgcount = invmass_h.Integral(invmass_h.GetXaxis().FindBin(mean - 9*sigma), invmass_h.GetXaxis().FindBin(mean - 4*sigma)) + invmass_h.Integral(invmass_h.GetXaxis().FindBin(mean + 4*sigma), invmass_h.GetXaxis().FindBin(mean + 9*sigma))
bgcount_err = math.sqrt(bgcount)

bglow = bgfunc.Integral(lowbin9, lowbin4)
bghigh = bgfunc.Integral(highbin4, highbin9)
bgpeak = bgfunc.Integral(lowpeak, highpeak)
alpha = bgpeak/(bglow + bghigh)

print("bglow: ", bglow)
print("bghigh: ", bghigh)
print("bgpeak: ", bgpeak)
print("alpha: ", alpha)

bgcount = bgcount * alpha
bgcount_err = bgcount_err * alpha

signal = counts - bgcount
signal = signal/0.9545
s_over_b = signal/bgcount
sig = signal/(math.sqrt(signal+bgcount))

print("signal: ", signal)
print("bg: ", bgcount)
print("s/b: ", s_over_b)
print("significance: ", sig)

rootfile = ROOT.TFile(root_filename, 'UPDATE')
invmass_h.Write()
for n in range(len(cut_names)-1):
	cut_h[n].Write()
rootfile.Close()

c = ROOT.TCanvas("c1", "c1", 700, 700)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
invmass_h.SetXTitle("inv mass")
invmass_h.SetTitle(args.energy + " " + args.system + " " + args.trailer)
invmass_h.Draw("E0")

c.SaveAs(plot_dir + "InvMass_" + args.energy + "_" + args.system + "_" + args.trailer + png)
