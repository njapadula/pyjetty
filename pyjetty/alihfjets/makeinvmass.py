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
import utils

from pyjetty.mputils import treereader

parser = argparse.ArgumentParser(description='D0 inv mass')
parser.add_argument('-f', '--ifile', help='input rootfile', type=str, default=None, required=True)
parser.add_argument('-e', '--energy', help='system energy', type=str, default='5TeV', required=False)
parser.add_argument('-s', '--system', help='collision system', type=str, default='pp', required=False)
parser.add_argument('-t', '--trailer', help='special trailer', type=str, default='', required=False)
parser.add_argument('-r', '--runevents', help='run over tree', type=int,  default=1, required=False)
args = parser.parse_args()

global main_dir, plot_dir, root_dir, ofile, png

main_dir = "/home/software/users/napadula/"
plot_dir = main_dir + "plots/"
root_dir = main_dir + "rootfiles/"
ofile = "HistsFromTree_" + args.energy + "_" + args.system + "_" + args.trailer + ".root"
png = args.system + args.energy + args.trailer + ".png"

def run_over_tree(hnamelist=[]):
	root_filename = root_dir + ofile
	rootfile = ROOT.TFile(root_filename, 'RECREATE')
	rootfile.Close()

	h_invmass_2D = []

	for h in range(len(hnamelist)):
		name = "InvMass_vs_pt_" + hnamelist[h]

		h_invmass_2D.insert(h, utils.make_hist(name, 55, 1.6, 2.15, 40, 0, 40))

	utils.fill_hist_from_tree('d0', utils.branch_names, args.ifile, h_invmass_2D, utils.cut_names, utils.cut_type, utils.cut_value_13TeV_Note, hnamelist)

	rootfile = ROOT.TFile(root_filename, 'UPDATE')
	for h in range(len(hnamelist)):
		h_invmass_2D[h].Write()
	rootfile.Close()

def open_hist(histname):
	root_filename = root_dir + ofile
	rootfile = ROOT.TFile.Open(root_filename, "READ")	

	name = "InvMass_vs_pt_" + histname
	#histname = "InvMass_vs_pt_reflection"
	hist = utils.get_hist_from_file(rootfile, name)

	rootfile.Close()
	return hist

def fit_and_plot(histname, fittype, fitrange=[], fitparam=[], fitparlim=[]):
	h_invmass = []
	c = []

	hist = open_hist(histname)

	print(utils.nbins)

	for p in range(utils.nbins):
		name = "InvMass" + histname + str(utils.pt_low[p]) + "_to_" + str(utils.pt_high[p])
	
		thefit = utils.set_fit_function(fittype, "thefit", fitrange, fitparam, fitparlim)

		ptlow = utils.get_bin(hist, utils.pt_low[p], False)
		pthigh = utils.get_bin(hist, utils.pt_high[p], False)
		
		h_invmass.insert(p, hist.ProjectionX(name, ptlow, pthigh-1))
		h_invmass[p].Fit("thefit","R")
	
		c.insert(p, ROOT.TCanvas("c"+str(p), "c"+str(p), 700, 700))
		ROOT.gStyle.SetOptStat(0)
		title_name = "Inv Mass " + histname + " {} < pt < {}".format(utils.pt_low[p], utils.pt_high[p])
		utils.setup_hist_to_draw(h_invmass[p], title_name, "inv mass", "", [1.7, 2.05])
		h_invmass[p].Draw("E0")
		title_name = "Inv Mass {} {} < pt < {}".format("Reflections",utils.pt_low[p],utils.pt_high[p])
		c[p].SaveAs(plot_dir + histname + "_InvMass_pt_" + str(utils.pt_low[p]) + "_to_" + str(utils.pt_high[p]) + png)


def integrate_hist(hist, xvalues=[], yvalues=[]):
	if xvalues:
		if yvalues:
			integral_val = hist.Integral(utils.get_bin(hist, xvalues[0]), utils.get_bin(hist, xvalues[1]), utils.get_bin(hist, yvalues[0], False), utils.get_bin(hist, yvalues[1], False))
		else:
			integral_val = hist.Integral(utils.get_bin(hist, xvalues[0]), utils.get_bin(hist, xvalues[1]))
	else:
		integral_val = hist.Integral()

	return integral_val

def main():
	#print("WTF")
	#run_over_tree(["signal","reflection"])
	hists = ["signal","reflection"]
	hinvmass=[]
	for i, hist in enumerate(hists):
		hinvmass.insert(i, open_hist(hist))
	#h_refoversig = make_var_hist("ReflectionPercentage", utils.nbins, utils.pt_edge)
	#xrange = [utils.get_bin(hinvmass[0], 1.7), utils.get_bin(hinvmass[0], 2.05)]
	for p in range(utils.nbins):
		sig_int = integrate_hist(hinvmass[0], [1.7, 2.05], [utils.pt_low[p], utils.pt_high[p]])
		ref_int = integrate_hist(hinvmass[1], [1.7, 2.05], [utils.pt_low[p], utils.pt_high[p]])
		ref_per = ref_int/sig_int
		ref_err = (ref_int/sig_int)*((math.sqrt(ref_int)/ref_int) + (math.sqrt(sig_int)/sig_int))
		print("{} +/- {}".format(ref_per, ref_err))
			

	#fit_and_plot("reflection", "dgaus", [1.7, 2.05], [0, 1.87, 0.05, 0, 1.87, 0.1])
	fit_and_plot("reflection", "dgaus", [1.7, 2.05], [20, 1.8, 0.5, 20, 1.8, 0.1])	

if __name__ == '__main__':
	main()
