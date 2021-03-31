from pyjetty.mputils import MPBase, pwarning, pinfo, perror
import random
import uproot
import pandas as pd
import fastjet as fj
import fjext
import os
import tqdm
# from numba import jit
# import numexpr

class HFAnalysis(MPBase):
	def __init__(self, **kwargs):
		super(HFAnalysis, self).__init__(**kwargs)
		self.selection = []
		self.df_selection = None
		self.query_strings = []
		self.query_string = ''
		self.callback = None

	def add_selection_equal(self, what, val):
		self.selection.append([what, val, None, 0])
		self.query_strings.append('({} == {})'.format(what, val))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_range(self, what, minv, maxv):
		self.selection.append([what, minv, maxv, 1])
		self.query_strings.append('({} > {}) & ({} < {})'.format(what, minv, what, maxv))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_range_abs(self, what, val):
		self.selection.append([what, val, None, 2])
		self.query_strings.append('({} > {}) & ({} < {})'.format(what, -val, what, +val))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_cond(self, what, val1, val2):
		self.selection.append([what, val1, val2, 3])
		self.query_strings.append('((({} > {}) & ({} < {})) | ({} < {}))'.format(what, -val2, what, +val2, what, val1))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def add_selection_nsig(self, tpcp0, tofp0, tpck1, tofk1, tpcp1, tofp1, tpck0, tofk0, val1, val2):
	#	self.selection.append([tpcp0, tofp0, None, 4])
		self.query_strings.append('((abs({}) < {} & (abs({}) < {} | {} < {}) & abs({}) < {} & (abs({}) < {} | {} < {})) | (abs({}) < {} & (abs({}) < {} | {} < {}) & abs({}) < {} & (abs({}) < {} | {} < {})))'.format(tpcp0, val1, tofp0, val1, tofp0, val2, tpck1, val1, tofk1, val1, tofk1, val2, tpcp1, val1, tofp1, val1, tofp1, val2, tpck0, val1, tofk0, val1, tofk0, val2))
		self.query_string = '(' + ' & '.join(self.query_strings) + ')'

	def compile_selection(self, df):
		self.df_selection = True
		for c in self.selection:
			if c[3] == 0:
				self.df_selection = (self.df_selection) & (df[c[0]] == c[1])
			if c[3] == 1:
				self.df_selection = (self.df_selection) & (df[c[0]] > c[1]) & (df[c[0]] < c[2])
			if c[3] == 2:
				self.df_selection = (self.df_selection) & (df[c[0]] > -c[1]) & (df[c[0]] < c[1])
			#if c[3] == 3:
			#	self.df_selection = (self.df_selection) & (((df[c[0]] > -c[2]) & (df[c[0]] < c[2])) | (df[c[0]] < c[1]))

	def analyze(self, df):
		self.compile_selection(df)
		_df = df[self.df_selection]
		self.analysis(_df)
		if self.callback is not None:
			self.callback(df['ev_id'].values[0])

	def analyze_slower(self, df):
		_df = df.query(self.query_string)
		#print(self.query_string)
		self.analysis(_df)

	def analyze_gen(self, df):
		_df = df
		self.analysis_gen(_df)

	# analysis on the single data frame
	# this is something specific to user - overload this one
	def analysis(self, df):
		if len(df) > 0:
			print (df)


class HFAnalysisIO(MPBase):
	def __init__(self, **kwargs):
		self.configure_from_args(tree_name='PWGHF_TreeCreator/tree_D0')
		self.event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
		self.configure_from_args(gen_tree_name='PWGHF_TreeCreator/tree_D0_gen')
		#self.gen_tree_name = 'PWGHF_TreeCreator/tree_D0_gen'
		super(HFAnalysisIO, self).__init__(**kwargs)
		self.analyses = []
		self.track_df_grouped = None
		self.d0_gen = None

	def reset_analyses_list(self):
		self.analyses = []

	def add_analysis(self, a):
		self.analyses.append(a)

	def load_file(self, path):
		if not os.path.exists(path):
			pwarning('[w] file', path, 'does not exists.')
			return
		try:
			event_tree = uproot.open(path)[self.event_tree_name]
		except:
			pwarning('error getting', self.event_tree_name, 'from file:', path)
			return False
		if not event_tree:
			perror('Tree {} not found in file {}'.format(self.event_tree_name, path))
			return False
		event_df_orig = event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
		event_df_orig.reset_index(drop=True)
		event_df = event_df_orig.query('is_ev_rej == 0')
		event_df.reset_index(drop=True)
		# Load gen tree into df
		try:
			gen_tree = uproot.open(path)[self.gen_tree_name]
		except:
			pwarning('error getting', self.gen_tree_name, 'from file:', path)
			return False
		if not gen_tree:
			perror('Tree {} not found in file {}'.format(gen_tree_name, path))
			return False
		gen_df_orig = gen_tree.pandas.df(['run_number', 'ev_id', 'pt_cand', 'eta_cand', 'cand_type'])
		gen_df_orig.sort_values(by=['run_number','ev_id'], inplace=True)
		df_genruns = gen_df_orig[['run_number','ev_id']].copy()
		# Load track tree into dataframe
		try:
			track_tree = uproot.open(path)[self.tree_name]
		except:
			pwarning('error getting', self.tree_name, 'from file:', path)
			return False
		if not track_tree:
			perror('Tree {} not found in file {}'.format(tree_name, path))
			return False
		track_df_orig = track_tree.pandas.df(['run_number', 'ev_id', 'inv_mass', 'pt_cand', 'pt_prong0', 'pt_prong1', 'dca', 'cos_t_star', 'imp_par_prod', 'cos_p', 'cand_type', 'imp_par_prong0', 'imp_par_prong1', 'norm_dl_xy', 'eta_cand', 'nsigTPC_Pi_0', 'nsigTOF_Pi_0', 'nsigTPC_K_1', 'nsigTOF_K_1', 'nsigTPC_Pi_1', 'nsigTOF_Pi_1', 'nsigTPC_K_0', 'nsigTOF_K_0'])
		# Merge event info into track tree
		track_df = pd.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
		track_df.sort_values(by=['run_number','ev_id'], inplace=True)
		df_d0runs = track_df[['run_number','ev_id']].copy()
		df_runs = pd.merge(df_d0runs, df_genruns, on=['run_number','ev_id'])
		df_runs.drop_duplicates(keep='first', inplace=True)
		gen_df_orig = pd.merge(gen_df_orig, df_runs, on=['run_number','ev_id'])
		track_df = pd.merge(track_df, df_runs, on=['run_number','ev_id'])
		self.d0_gen = gen_df_orig.groupby(['run_number','ev_id'])
		self.track_df_grouped = track_df.groupby(['run_number','ev_id'])
		return True

	def execute_analyses(self):
		[self.track_df_grouped.apply(a.analyze_slower) for a in self.analyses]
		[self.d0_gen.apply(a.analyze_gen) for a in self.analyses]

	def update_status(self, mark):
		if mark != self.pbar2_mark:
			self.pbar2_mark = mark
			self.pbar2.update(1)

	def execute_analyses_on_file_list(self, file_list, nfiles=0):
		self.pbar2 = tqdm.tqdm(mininterval=20, maxinterval=60)
		self.pbar2_mark = None
		for a in self.analyses:
			a.callback = self.update_status
		print()
		if os.path.exists(file_list):
			with open(file_list) as f:
				files = f.readlines()
			if int(nfiles) > 0:
				files = files[:nfiles]
			for f in files:
				fn = f.strip('\n')
				pinfo('+file:', fn)
			for f in tqdm.tqdm(files):
				fn = f.strip('\n')
				if self.load_file(fn):
					self.execute_analyses()
			self.pbar2.close()
		else:
			perror('file list does not exist', file_list)
		pinfo('done.')

	def __def__(self):
		self.track_df_grouped = None
		self.d0_gen = None
