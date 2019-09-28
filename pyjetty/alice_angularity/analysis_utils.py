#!/usr/bin/env python3

"""
  Analysis utilities for jet analysis with track dataframe.
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import math

# Data analysis and plotting
from collections import defaultdict
import uproot
import pandas
import numpy as np
import ROOT

# Fastjet via python (from external library fjpydev)
import fastjet as fj
import fjext

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

#---------------------------------------------------------------
def analysis_utils():

  print('These are some analysis utilities for jet analysis with track dataframe')

#---------------------------------------------------------------
# Convert ROOT TTree to pandas dataframe
# Return merged track+event dataframe from a given input file
# Returned dataframe has one row per jet constituent:
#     run_number, ev_id, ParticlePt, ParticleEta, ParticlePhi
#---------------------------------------------------------------
def load_dataframe(inputFile, tree_name):
  
  # Load event tree into dataframe, and apply event selection
  event_tree_name = 'PWGHF_TreeCreator/tree_event_char'
  event_tree = uproot.open(inputFile)[event_tree_name]
  if not event_tree:
    print('Tree {} not found in file {}'.format('tree_event_char', inputFile))
  event_df_orig = event_tree.pandas.df(['run_number', 'ev_id', 'z_vtx_reco','is_ev_rej'])
  event_df_orig.reset_index(drop=True)
  event_df = event_df_orig.query('is_ev_rej == 0')
  event_df.reset_index(drop=True)

  # Load track tree into dataframe
  track_tree_name = 'PWGHF_TreeCreator/{}'.format(tree_name)
  track_tree = uproot.open(inputFile)[track_tree_name]
  if not track_tree:
    print('Tree {} not found in file {}'.format(tree_name, inputFile))
  track_df_orig = track_tree.pandas.df()

  # Merge event info into track tree
  track_df = pandas.merge(track_df_orig, event_df, on=['run_number', 'ev_id'])
  return track_df

#---------------------------------------------------------------
# Transform the track dataframe into a SeriesGroupBy object
# of fastjet particles per event.
#---------------------------------------------------------------
def group_fjparticles(track_df):

  print('Transform the track dataframe into a series object of fastjet particles per event...')

  # (i) Group the track dataframe by event
  #     track_df_grouped is a DataFrameGroupBy object with one track dataframe per event
  track_df_grouped = track_df.groupby(['run_number','ev_id'])
  
  # (ii) Transform the DataFrameGroupBy object to a SeriesGroupBy of fastjet particles
  df_fjparticles = track_df_grouped.apply(get_fjparticles)
  
  return df_fjparticles

#---------------------------------------------------------------
# Return fastjet:PseudoJets from a given track dataframe
#---------------------------------------------------------------
def get_fjparticles(df_tracks):
  
  # Use swig'd function to create a vector of fastjet::PseudoJets from numpy arrays of pt,eta,phi
  fj_particles = fjext.vectorize_pt_eta_phi(df_tracks['ParticlePt'].values, df_tracks['ParticleEta'].values, df_tracks['ParticlePhi'].values)
  
  return fj_particles

#---------------------------------------------------------------
# Check if det-jet passes acceptance criteria
#---------------------------------------------------------------
def is_det_jet_accepted(jet_det):

  for track in jet_det.constituents():

    if track.pt() > 100.:
      return False

  return True

#---------------------------------------------------------------
# Normalize a histogram by its integral
#---------------------------------------------------------------
def scale_by_integral(h):
  
  if h.GetSumw2() is 0:
    h.Sumw2()
  
  integral = h.Integral()
  if integral > 0:
    h.Scale(1./integral)
  else:
    print('Integral is 0, check for problem')

#---------------------------------------------------------------
# Plot and save a 1D histogram
#---------------------------------------------------------------
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False):
  
  h.SetLineColor(1)
  h.SetLineWidth(1)
  h.SetLineStyle(1)
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  ROOT.gPad.SetLeftMargin(0.15)
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  ROOT.gPad.SetLeftMargin(0.15)

  h.Draw(drawOptions)
  c.SaveAs(outputFilename)
  c.Close()

#---------------------------------------------------------------
# Initialize 3D dictionary with keys being lists of binned lambda values.
# { (pTbin1): { (jetR1): { (beta1): [lambdabin1, lambdabin2, ... ],
#                          (beta2): [ ... ], 
#                          ... },
#               (jetR2): { ... },
#               ... },
#   (pTbin2:) { ... },
#   ... }
#---------------------------------------------------------------
def init_jet_lambda(pTbins, jetRs, betas, nbins):
  jet_lambda = defaultdict(lambda : defaultdict(dict))
  for pTbin in pTbins:
    for jetR in jetRs:
      for beta in betas:
        jet_lambda[pTbin][jetR][beta] = [0] * self.n_lambda_bins
  return jet_lambda

#---------------------------------------------------------------
# Calculate the jet substructure variable lambda for given jet, beta, kappa, and jet R
#---------------------------------------------------------------
def calc_lambda_single_jet(jet, b, k, jet_R):

  # Calculate & sum lambda for all jet constituents
  lambda_bk = 0
  jet_pT = np.sqrt(jet.px()**2 + jet.py()**2)
  for constituent in jet.consituents():
    deltaR = np.sqrt( (jet.rap() - constituent.rap())**2 + (jet.phi() - constituent.phi())**2 )
    constit_pT = np.sqrt(constituent.px()**2 + constituent.py()**2)
    lambda_bk += ( constit_pT / jet_pT )**k * (deltaR / jet_R)**b

  return lambda_bk

#----------------------------------------------------------------------
if __name__ == '__main__':
  
  analysis_utils()