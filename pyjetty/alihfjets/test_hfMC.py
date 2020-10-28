#!/usr/bin/env python3

import argparse
import os
import pyjetty.alihfjets.hf_MC_io as hfdio
from pyjetty.mputils import perror, pinfo, pwarning, treewriter
import ROOT
ROOT.gROOT.SetBatch(True)

class HFAnalysisInvMass(hfdio.HFAnalysis):
	def __init__(self, **kwargs):
		self.fout = None
		super(HFAnalysisInvMass, self).__init__(**kwargs)
		self.fout = ROOT.TFile(self.name+'.root', 'recreate')
		self.fout.cd()
		self.tw = treewriter.RTreeWriter(tree_name='d0', fout=self.fout)
		self.tw_gen = treewriter.RTreeWriter(tree_name='d0_gen', fout=self.fout)
		

	def analysis(self, df):
		for index, row in df.iterrows():
			for c in df.columns:
				self.tw.fill_branch(c, row[c])
			self.tw.fill_tree()

	def analysis_gen(self, df):
		for index, row in df.iterrows():
			for c in df.columns:
				self.tw_gen.fill_branch(c, row[c])
			self.tw_gen.fill_tree()
				
	def finalize(self):
		self.fout.Write()
		self.fout.Close()
		pinfo(self.fout.GetName(), 'written.')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='D0 analysis on alice data', prog=os.path.basename(__file__))
	parser.add_argument('-f', '--flist', help='file list to process', type=str, default=None, required=True)
	parser.add_argument('-n', '--nfiles', help='max n files to process', type=int, default=0, required=False)
	parser.add_argument('-o', '--output', help="output name / file name in the end", type=str, default='test_hfana')
	args = parser.parse_args()

	hfaio = hfdio.HFAnalysisIO()

	hfa = HFAnalysisInvMass(name = args.output)
	hfa.add_selection_range('pt_cand', 2, 1e3)
	hfa.add_selection_range_abs('z_vtx_reco', 10)
	hfa.add_selection_range('pt_prong0', 0.5, 1e3)
	hfa.add_selection_range('pt_prong1', 0.5, 1e3)
	hfa.add_selection_range_abs('eta_cand', 0.8)
	hfa.add_selection_nsig('nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0', 3, -900)

	hfaio.add_analysis(hfa)

	# hfaio.load_file("./AnalysisResults.root")
	# hfaio.execute_analyses()
	hfaio.execute_analyses_on_file_list(args.flist, args.nfiles)

	hfa.finalize()
