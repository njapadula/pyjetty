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
rootfile = ROOT.TFile(root_filename, 'RECREATE')
rootfile.Close()

branch_names = ['inv_mass', 'pt_cand', 'pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p']
gen_branch_names = ['pt_cand']
cut_names = ['pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p']
cut_type = ['greater', 'greater', 'less', 'abs', 'less', 'greater']
cut_value_var = [[0.7, 0.7, 0.03, 0.8, -0.0002, 0.9],
[0.7, 0.7, 0.03, 0.8, -0.0002, 0.9],
[0.7, 0.7, 0.03, 0.8, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.8, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.8, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.9, -0.00005, 0.85],
[0.7, 0.7, 0.03, 0.9, -0.00005, 0.85],
[0.6, 0.6, 0.03, 1, 10, 0.8]]

cut_value = [[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9],
[0.6, 0.6, 0.03, 0.8, -0.0001, 0.9]]

h_invmass = []

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

h_genD0_pt_yield = ROOT.TH1F("h_genD0_pt_yield", "h_genD0_pt_yield", nbins, pt_edge)
h_genD0_pt_yield.Sumw2()

for p in range(len(pt_low)):
	name = "InvMass_{}_to_{}".format(pt_low[p],pt_high[p])
	h_invmass.insert(p, ROOT.TH1F(name, name, 45, 1.65, 2.1))
	h_invmass[p].Sumw2()

	print(name)

tr = treereader.RTreeReader(tree_name='d0', branches = branch_names, file_name = args.ifile)
	
for i in range(tr.tree.GetEntries()):
	tr.tree.GetEntry(i)
	if(i % 100000 == 0):
		print(i)
	for n in range(nbins):
		if (tr.pt_cand[0] >= pt_low[n]) and (tr.pt_cand[0] < pt_high[n]):
			for c in range(len(cut_names)):
				flag = make_cut(tr, cut_names[c], cut_type[c], cut_value[n][c])
				if not flag:
					break
			if flag:
				h_invmass[n].Fill(tr.inv_mass[0])


trg = treereader.RTreeReader(tree_name='d0_gen', branches = gen_branch_names, file_name = args.ifile)
for i in range(trg.tree.GetEntries()):
	trg.tree.GetEntry(i)
	h_genD0_pt_yield.Fill(trg.pt_cand[0])

			
rootfile = ROOT.TFile(root_filename, 'UPDATE')
for n in range(nbins):
	h_invmass[n].Write()
h_genD0_pt_yield.Write()
rootfile.Close()

