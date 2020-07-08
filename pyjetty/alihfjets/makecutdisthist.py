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

main_dir = "/home/software/users/napadula/"
root_dir = main_dir + "rootfiles/"
ofile = "InvMassCutSelectionHists" + args.system + args.energy + args.trailer + ".root"

root_filename = root_dir + ofile
rootfile = ROOT.TFile(root_filename, 'RECREATE')
rootfile.Close()

invmass_h = []
cut_h = []
branch_names = ['inv_mass', 'pt_cand', 'pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prong0', 'imp_par_prong1', 'imp_par_prod', 'cos_p']
cut_names = ['pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prong0', 'imp_par_prong1', 'imp_par_prod', 'cos_p', 'allcuts']
cut_type = ['greater', 'greater', 'less', 'abs', 'abs', 'abs', 'less', 'greater', 'all']
cut_value = [0.6, 0.6, 0.03, 0.8, 0.1, 0.1, -0.0001, 0.9, 0]	

pt_low = 5
pt_high = 100

for n in range(len(cut_names)):
	invmass_h.insert(n, ROOT.TH1F("InvMass_"+cut_names[n], "InvMass_"+cut_names[n], 45, 1.65, 2.10))
	invmass_h[n].Sumw2()	

f = ROOT.TFile(args.ifile, 'READ')
d0 = f.Get("d0")

for n in range(len(cut_names)-1):
	this_cut = cut_names[n] + ' >> hnew'
	d0.Draw(this_cut, "pt_cand>5", "goff")
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
for i in range(tr.tree.GetEntries()):
	tr.tree.GetEntry(i)
	if(i % 100000 == 0):
		print(i)
	if((tr.pt_cand[0] >= pt_low) and (tr.pt_cand[0] < pt_high) and (tr.inv_mass[0] > 0)):
		flag = 0
		for n in range(len(cut_names)-1):
			if make_cut(tr, cut_names[n], cut_type[n], cut_value[n]):
				invmass_h[n].Fill(tr.inv_mass[0])
				flag = flag + 1		
		if flag == (len(cut_names)-1): 
			invmass_h[-1].Fill(tr.inv_mass[0])


rootfile = ROOT.TFile(root_filename, 'UPDATE')
for n in range(len(cut_names)-1):
	invmass_h[n].Write()
	cut_h[n].Write()
invmass_h[len(cut_names)-1].Write()
#hnew.Write()
#cut_h.Write()
rootfile.Close()
