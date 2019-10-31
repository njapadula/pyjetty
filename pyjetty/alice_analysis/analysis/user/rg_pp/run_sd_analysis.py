#! /usr/bin/env python

import sys
import os
import argparse
import itertools
from array import *
import numpy
import math
import ROOT
import yaml

from pyjetty.alice_analysis.analysis.base import common_base
from pyjetty.alice_analysis.analysis.base import analysis_utils
import roounfold_sd

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class run_sd_analysis(common_base.common_base):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, config_file='', **kwargs):
    super(run_sd_analysis, self).__init__(**kwargs)
    self.config_file = config_file
    
    # Initialize utils class
    self.utils = analysis_utils.analysis_utils()
    
    # Initialize yaml config
    self.initialize_config()

    print(self)
  
  #---------------------------------------------------------------
  # Initialize config file into class members
  #---------------------------------------------------------------
  def initialize_config(self):
    
    # Read config file
    with open(self.config_file, 'r') as stream:
      config = yaml.safe_load(stream)
    
    # Set list of observables
    self.observables = config['observables']
    
    # Set which analysis steps to perform
    self.do_unfolding = config['do_unfolding']
    self.do_systematics = config['do_systematics']
    self.do_plot_final_result = config['do_plot_final_result']
    self.force_rebin_response=config['force_rebin']
    
    # Retrieve list of SD grooming settings
    self.jetR_list = config['jetR']
    self.sd_config_dict = config['SoftDrop']
    self.sd_config_list = list(self.sd_config_dict.keys())
    self.sd_settings = [[self.sd_config_dict[name]['zcut'], self.sd_config_dict[name]['beta']] for name in self.sd_config_list]
    
    for observable in self.observables:
    
      if observable == 'theta_g':
        xtitle = '#it{#theta}_{g}'
        ytitle = '#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#it{#theta}_{g}}'
      if observable == 'zg':
        xtitle = '#it{z}_{g}'
        ytitle = '#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#it{z}_{g}}'
      setattr(self, 'xtitle_{}'.format(observable), xtitle)
      setattr(self, 'ytitle_{}'.format(observable), ytitle)

    # Retrieve histogram binnings for each SD setting
    for i, sd_setting in enumerate(self.sd_settings):
        
      zcut = sd_setting[0]
      beta = sd_setting[1]
      sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
      config_name = self.sd_config_list[i]
      
      pt_bins_truth = (self.sd_config_dict[config_name]['pt_bins_truth'])
      
      n_pt_bins_truth = len(pt_bins_truth) - 1
      setattr(self, 'n_pt_bins_truth_{}'.format(sd_label), n_pt_bins_truth)

      truth_pt_bin_array = array('d',pt_bins_truth)
      setattr(self, 'truth_pt_bin_array_{}'.format(sd_label), truth_pt_bin_array)

      if 'theta_g' in self.observables:
      
        rg_bins_truth = (self.sd_config_dict[config_name]['rg_bins_truth'])
      
        n_rg_bins_truth = len(rg_bins_truth) - 1
        setattr(self, 'n_rg_bins_truth_{}'.format(sd_label), n_rg_bins_truth)
        
        truth_rg_bin_array = array('d',rg_bins_truth)
        setattr(self, 'truth_rg_bin_array_{}'.format(sd_label), truth_rg_bin_array)
      
      if 'zg' in self.observables:
        
        zg_bins_truth = (self.sd_config_dict[config_name]['zg_bins_truth'])

        n_zg_bins_truth = len(zg_bins_truth) - 1
        setattr(self, 'n_zg_bins_truth_{}'.format(sd_label), n_zg_bins_truth)
        
        truth_zg_bin_array = array('d',zg_bins_truth)
        setattr(self, 'truth_zg_bin_array_{}'.format(sd_label), truth_zg_bin_array)

    # List of systematic variations to perform
    self.systematics_list = config['systematics_list']
    
    # Load paths to processing output, to be unfolded
    self.main_data = config['main_data']
    self.main_response = config['main_response']
    
    if 'kTrackEff' in self.systematics_list:
      self.trkeff_data = config['trkeff_data']
      self.trkeff_response = config['trkeff_response']
          
    # Create output dirs
    self.file_format = config['file_format']
    self.output_dir = config['output_dir']
    for observable in self.observables:
      self.create_output_dirs(observable)

  #---------------------------------------------------------------
  # Create a set of output directories for a given observable
  #---------------------------------------------------------------
  def create_output_dirs(self, observable):
    
    output_dir = os.path.join(self.output_dir, observable)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    self.create_output_subdir(observable, output_dir, 'main')
      
    if 'kTrackEff' in self.systematics_list:
      self.create_output_subdir(observable, output_dir, 'trkeff')
    if 'kPrior1' in self.systematics_list:
      self.create_output_subdir(observable, output_dir, 'prior1')
    if 'kPrior2' in self.systematics_list:
      self.create_output_subdir(observable, output_dir, 'prior2')
    if 'kTruncation' in self.systematics_list:
      self.create_output_subdir(observable, output_dir, 'truncation')
    if 'kBinning' in self.systematics_list:
      self.create_output_subdir(observable, output_dir, 'binning')

    if self.do_systematics:
      output_dir_systematics = self.create_output_subdir(observable, output_dir, 'systematics')
      sys_root_filename = os.path.join(output_dir_systematics, 'fSystematics.root')
      fSystematics = ROOT.TFile(sys_root_filename, 'RECREATE')
  
    if self.do_plot_final_result:
      output_dir_final = self.create_output_subdir(observable, output_dir, 'final_results')
      final_result_root_file = os.path.join(output_dir_final, 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_file, 'RECREATE')

  #---------------------------------------------------------------
  # Create a single output subdirectory
  #---------------------------------------------------------------
  def create_output_subdir(self, observable, output_dir, name):
    
    output_subdir = os.path.join(output_dir, name)
    setattr(self, 'output_dir_{}_{}'.format(name, observable), output_subdir)
    if not os.path.isdir(output_subdir):
      os.makedirs(output_subdir)

    return output_subdir

  #---------------------------------------------------------------
  # Main processing function
  #---------------------------------------------------------------
  def run_sd_analysis(self):

    for observable in self.observables:
    
      if self.do_unfolding:
        self.perform_unfolding(observable)
      
      for jetR in self.jetR_list:
        for sd_setting in self.sd_settings:
          
          zcut = sd_setting[0]
          beta = sd_setting[1]
          sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)

          if self.do_systematics:
            self.compute_systematics(observable, jetR, sd_label, zcut, beta)

          if self.do_plot_final_result:
            self.plot_final_result(observable, jetR, sd_label, zcut, beta)

        self.plot_final_result_overlay(observable, jetR)
      
  #----------------------------------------------------------------------
  def perform_unfolding(self, observable):
    print('Perform unfolding for all systematic variations: {} ...'.format(observable))
    
    # Main result
    output_dir = getattr(self, 'output_dir_main_{}'.format(observable))
    analysis_main = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir))
    analysis_main.roounfold_sd()

    # Tracking efficiency variation
    if 'kTrackEff' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_trkeff_{}'.format(observable))
      analysis_trkeff = roounfold_sd.roounfold_sd(observable, self.trkeff_data, self.trkeff_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir))
      analysis_trkeff.roounfold_sd()

    # Prior variation 1
    if 'kPrior1' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_prior1_{}'.format(observable))
      analysis_prior1 = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir), power_law_offset=0.5)
      analysis_prior1.roounfold_sd()
    
    # Prior variation 2
    if 'kPrior2' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_prior2_{}'.format(observable))
      analysis_prior2 = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir), power_law_offset=-0.5)
      analysis_prior2.roounfold_sd()

    # Truncation variation
    if 'kTruncation' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_truncation_{}'.format(observable))
      analysis_truncation = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir), truncation=True)
      analysis_truncation.roounfold_sd()
        
    # Binning variation
    if 'kBinning' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_binning_{}'.format(observable))
      analysis_binning = roounfold_sd.roounfold_sd(observable, self.main_data, self.main_response, self.config_file, output_dir, self.file_format, rebin_response=self.check_rebin_response(output_dir), binning=True)
      analysis_binning.roounfold_sd()

  #----------------------------------------------------------------------
  def check_rebin_response(self, output_dir):
    
    rebin_response = True
    response_path = os.path.join(output_dir, 'response.root')
    if os.path.exists(response_path):
      rebin_response = False
      if self.force_rebin_response:
        print('Response {} exists -- force re-create...'.format(response_path))
        rebin_response = True
      else:
        print('Response {} exists -- don\'t re-create.'.format(response_path))
    else:
      print('Response {} doesn\'t exist -- create it...'.format(response_path))

    return rebin_response

  #----------------------------------------------------------------------
  def compute_systematics(self, observable, jetR, sd_label, zcut, beta):
    print('Compute systematics for {}: R = {}, {} ...'.format(observable, jetR, sd_label))

    # Get main result
    output_dir = getattr(self, 'output_dir_main_{}'.format(observable))
    path_main = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
    fMain = ROOT.TFile(path_main, 'READ')
    
    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, observable, jetR)

    name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
    hMain = fMain.Get(name)
    hMain.SetDirectory(0)
    setattr(self, name, hMain)
    
    # Tagging rate histogram
    name = 'hTaggingFractions_R{}_{}'.format(jetR, sd_label)
    hTaggingFractions = fMain.Get(name)
    hTaggingFractions.SetDirectory(0)
    setattr(self, name, hTaggingFractions)

    # Regularization parameter +2
    name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final+2)
    hRegParam1 = fMain.Get(name)
    hRegParam1.SetDirectory(0)
    setattr(self, name, hRegParam1)
    
    # Regularization parameter -2
    name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final-2)
    hRegParam2 = fMain.Get(name)
    hRegParam2.SetDirectory(0)
    setattr(self, name, hRegParam2)
    
    # Get trkeff result
    if 'kTrackEff' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_trkeff_{}'.format(observable))
      path_trkeff = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fTrkEff = ROOT.TFile(path_trkeff, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hTrkEff = fTrkEff.Get(name)
      hTrkEff.SetDirectory(0)
      setattr(self, '{}_trkeff'.format(name), hTrkEff)
    
    # Get prior result
    if 'kPrior1' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_prior1_{}'.format(observable))
      path_prior1 = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fPrior1 = ROOT.TFile(path_prior1, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hPrior1 = fPrior1.Get(name)
      hPrior1.SetDirectory(0)
      setattr(self, '{}_prior1'.format(name), hPrior1)
    
    # Get prior result
    if 'kPrior2' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_prior2_{}'.format(observable))
      path_prior2 = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fPrior2 = ROOT.TFile(path_prior2, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hPrior2 = fPrior2.Get(name)
      hPrior2.SetDirectory(0)
      setattr(self, '{}_prior2'.format(name), hPrior2)
    
    # Get truncation result
    if 'kTruncation' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_truncation_{}'.format(observable))
      path_truncation = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fTruncation = ROOT.TFile(path_truncation, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hTruncation = fTruncation.Get(name)
      hTruncation.SetDirectory(0)
      setattr(self, '{}_truncation'.format(name), hTruncation)
    
    # Get binning result
    if 'kBinning' in self.systematics_list:
      output_dir = getattr(self, 'output_dir_binning_{}'.format(observable))
      path_binning = os.path.join(output_dir, 'fResult_R{}_{}.root'.format(jetR, sd_label))
      fBinning = ROOT.TFile(path_binning, 'READ')
      name = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
      hBinning = fBinning.Get(name)
      hBinning.SetDirectory(0)
      setattr(self, '{}_binning'.format(name), hBinning)
    
    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
    
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.compute_sd_observable_systematic(observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth)

  #----------------------------------------------------------------------
  def get_sd_observable_distribution(self, jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth):
  
    h2D = getattr(self, name2D)
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = h2D.ProjectionY() # Better to use ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) ?
    h.SetName(name1D)
    h.SetDirectory(0)
    
    name = 'hTaggingFractions_R{}_{}'.format(jetR, sd_label)
    hTaggingFrac = getattr(self, name)
    x = (min_pt_truth + max_pt_truth)/2.
    fraction_tagged =  hTaggingFrac.GetBinContent(hTaggingFrac.FindBin(x))
    setattr(self, '{}_fraction_tagged'.format(name1D), fraction_tagged)

    n_jets_tagged = h.Integral(1, h.GetNbinsX())
    n_jets_inclusive = n_jets_tagged/fraction_tagged
    h.Scale(1./n_jets_inclusive, 'width')
  
    setattr(self, name1D, h)

    return h

  #----------------------------------------------------------------------
  def compute_sd_observable_systematic(self, observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth):
    
    #------------------------------------
    # Get 1D histograms
    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin (cross-check this)
    
    reg_param_final = self.utils.get_reg_param(self.sd_settings, self.sd_config_list, self.sd_config_dict, sd_label, observable, jetR)
    
    # Get main histogram
    name2D = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final)
    name1D = 'hMain_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hMain = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)

    # Get reg param +2
    name2D = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final+2)
    name1D = 'hRegParam1_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hRegParam1 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get reg param -2
    name2D = 'hUnfolded_{}_R{}_{}_{}'.format(observable, jetR, sd_label, reg_param_final-2)
    name1D = 'hRegParam2_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hRegParam2 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)

    # Get trk eff
    if 'kTrackEff' in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_trkeff'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hTrkEff_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hTrkEff = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get prior1
    if 'kPrior1' in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_prior1'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hPrior1_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hPrior1 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get prior2
    if 'kPrior2' in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_prior2'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hPrior2_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hPrior2 = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get truncation
    if 'kTruncation' in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_truncation'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hTruncation_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hTruncation = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    # Get binning
    if 'kBinning' in self.systematics_list:
      name2D = 'hUnfolded_{}_R{}_{}_{}_binning'.format(observable, jetR, sd_label, reg_param_final)
      name1D = 'hBinning_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hBinning = self.get_sd_observable_distribution(jetR, sd_label, name2D, name1D, min_pt_truth, max_pt_truth)
    
    #------------------------------------
    # Compute systematics

    # Reg param +2
    name = 'hSystematic_{}_RegParam1_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hSystematic_RegParam1 = hMain.Clone()
    hSystematic_RegParam1.SetName(name)
    hSystematic_RegParam1.Divide(hRegParam1)
    self.change_to_per(hSystematic_RegParam1)
    setattr(self, name, hSystematic_RegParam1)

    # Reg param -2
    name = 'hSystematic_{}_RegParam2_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hSystematic_RegParam2 = hMain.Clone()
    hSystematic_RegParam2.SetName(name)
    hSystematic_RegParam2.Divide(hRegParam2)
    self.change_to_per(hSystematic_RegParam2)
    setattr(self, name, hSystematic_RegParam2)
    
    name = 'hSystematic_{}_RegParam_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hSystematic_RegParam = self.build_average(hSystematic_RegParam1, hSystematic_RegParam2)
    setattr(self, name, hSystematic_RegParam)
    
    output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
    outputFilename = os.path.join(output_dir, 'hSystematic_RegParam_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
    #self.utils.plot_hist(hSystematic_RegParam, outputFilename, 'P E')
    
    # Prior 1
    hSystematic_Prior1 = None
    if 'kPrior1' in self.systematics_list:
      name = 'hSystematic_{}_Prior1_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_Prior1 = hMain.Clone()
      hSystematic_Prior1.SetName(name)
      hSystematic_Prior1.Divide(hPrior1)
      self.change_to_per(hSystematic_Prior1)
      setattr(self, name, hSystematic_Prior1)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_Prior1_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
      #self.utils.plot_hist(hSystematic_Prior1, outputFilename, 'P E')

    # Prior 2
    hSystematic_Prior2 = None
    if 'kPrior2' in self.systematics_list:
      name = 'hSystematic_{}_Prior2_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_Prior2 = hMain.Clone()
      hSystematic_Prior2.SetName(name)
      hSystematic_Prior2.Divide(hPrior2)
      self.change_to_per(hSystematic_Prior2)
      setattr(self, name, hSystematic_Prior2)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_Prior2_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
      #self.utils.plot_hist(hSystematic_Prior2, outputFilename, 'P E')

    # Truncation
    hSystematic_Truncation = None
    if 'kTruncation' in self.systematics_list:
      name = 'hSystematic_{}_Truncation_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_Truncation = hMain.Clone()
      hSystematic_Truncation.SetName(name)
      hSystematic_Truncation.Divide(hTruncation)
      self.change_to_per(hSystematic_Truncation)
      setattr(self, name, hSystematic_Truncation)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_Truncation_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
      #self.utils.plot_hist(hSystematic_Truncation, outputFilename, 'P E')
    
    # Binning
    hSystematic_Binning = None
    if 'kBinning' in self.systematics_list:
      name = 'hSystematic_{}_Binning_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_Binning = hMain.Clone()
      hSystematic_Binning.SetName(name)
      hSystematic_Binning.Divide(hBinning)
      self.change_to_per(hSystematic_Binning)
      setattr(self, name, hSystematic_Binning)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_Binning_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
      #self.utils.plot_hist(hSystematic_Binning, outputFilename, 'P E')
    
    # Trk eff
    if 'kTrackEff' in self.systematics_list:
      name = 'hSystematic_{}_TrkEff_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      hSystematic_TrkEff = hMain.Clone()
      hSystematic_TrkEff.SetName(name)
      hSystematic_TrkEff.Divide(hTrkEff)
      self.change_to_per(hSystematic_TrkEff)
      setattr(self, name, hSystematic_TrkEff)
      
      output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
      outputFilename = os.path.join(output_dir, 'hSystematic_TrkEff_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
      #self.utils.plot_hist(hSystematic_TrkEff, outputFilename, 'P E')

    # Add uncertainties in quadrature
    hSystematic_Total = self.add_in_quadrature(hSystematic_RegParam, hSystematic_Prior1, hSystematic_Prior2, hSystematic_Truncation, hSystematic_Binning, hSystematic_TrkEff)
    output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
    outputFilename = os.path.join(output_dir, 'hSystematic_Total_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
    #self.utils.plot_hist(hSystematic_Total, outputFilename, 'P E')

    name = 'hResult_{}_systotal_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    hResult_sys = hMain.Clone()
    hResult_sys.SetName(name)
    hResult_sys.SetDirectory(0)
    self.AttachErrToHist(hResult_sys, hSystematic_Total)
    setattr(self, name, hResult_sys)
      
    # Plot systematic uncertainties
    h_list = [hSystematic_RegParam, hSystematic_Prior1, hSystematic_Prior2, hSystematic_Truncation, hSystematic_Binning, hSystematic_TrkEff]
    self.plot_systematic_uncertainties(observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, h_list, hSystematic_Total)
      
  #----------------------------------------------------------------------
  def plot_systematic_uncertainties(self, observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, h_list, h_total):
  
    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()
    
    name = 'cSys_R{}_{}_{}-{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    if observable == 'theta_g':
      n_bins_truth = getattr(self, 'n_rg_bins_truth_{}'.format(sd_label))
      truth_bin_array = getattr(self, 'truth_rg_bin_array_{}'.format(sd_label))
    if observable == 'zg':
      n_bins_truth = getattr(self, 'n_zg_bins_truth_{}'.format(sd_label))
      truth_bin_array = getattr(self, 'truth_zg_bin_array_{}'.format(sd_label))

    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_bins_truth, truth_bin_array)
    myBlankHisto.SetNdivisions(505)
    xtitle = getattr(self, 'xtitle_{}'.format(observable))
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle('Systematic uncertainty (%)')
    myBlankHisto.SetMaximum(1.5*h_total.GetMaximum())
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")
    
    leg = ROOT.TLegend(0.67,0.65,0.8,0.92)
    self.utils.setup_legend(leg,0.04)

    for i, h in enumerate(h_list):
      if i == 0:
        h.SetMarkerStyle(20)
        label = 'Reg. param'
      if i == 1:
        h.SetMarkerStyle(21)
        label = 'Prior1'
      if i == 2:
        h.SetMarkerStyle(22)
        label = 'Prior2'
      if i == 3:
        h.SetMarkerStyle(23)
        label = 'Truncation'
      if i == 4:
        h.SetMarkerStyle(33)
        label = 'Binning'
      if i == 5:
        h.SetMarkerStyle(34)
        label = 'Tracking efficiency'
      
      h.SetMarkerSize(1.5)
      h.SetMarkerColor(600-5+i)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(600-5+i)

      h.DrawCopy('P X0 same')

      leg.AddEntry(h, label, 'Pe')

    h_total.SetLineStyle(1)
    h_total.SetLineColor(1)
    h_total.SetLineWidth(2)
    h_total.DrawCopy('same hist')
    leg.AddEntry(h_total, 'Total', 'l')

    leg.Draw()

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR)
    text_latex.DrawLatex(0.25, 0.85, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.78, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'z_{cut} = ' + str(zcut) + '   #beta = ' + str(beta)
    text_latex.DrawLatex(0.25, 0.71, text)

    output_dir = getattr(self, 'output_dir_systematics_{}'.format(observable))
    outputFilename = os.path.join(output_dir, 'hSystematics_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format))
    c.SaveAs(outputFilename)
    c.Close()

    sys_root_filename = os.path.join(output_dir, 'fSystematics.root')
    fSystematics = ROOT.TFile(sys_root_filename, 'UPDATE')
    h_total.Write()
    fSystematics.Close()
      
  #----------------------------------------------------------------------
  def add_in_quadrature(self, h1, h2, h3, h4, h5, h6):
  
    h_new = h1.Clone()
    h_new.SetName('{}_new'.format(h1.GetName()))
    
    for i in range(1, h_new.GetNbinsX()+1):
      value1 = h1.GetBinContent(i)
      value2 = h2.GetBinContent(i)
      value3 = h3.GetBinContent(i)
      value4 = h4.GetBinContent(i)
      value5 = h5.GetBinContent(i)
      value6 = h6.GetBinContent(i)

      new_value = math.sqrt(value1*value1 + value2*value2  + value3*value3  + value4*value4  + value5*value5  + value6*value6)
    
      h_new.SetBinContent(i, new_value)
    
    return h_new
  
  #----------------------------------------------------------------------
  def plot_final_result(self, observable, jetR, sd_label, zcut, beta):
    print('Plot final results for {}: R = {}, {} ...'.format(observable, jetR, sd_label))

    self.utils.set_plotting_options()
    ROOT.gROOT.ForceStyle()

    # Loop through pt slices, and compute systematics for each 1D theta_g distribution
    n_pt_bins_truth = getattr(self, 'n_pt_bins_truth_{}'.format(sd_label))
    truth_pt_bin_array = getattr(self, 'truth_pt_bin_array_{}'.format(sd_label))
  
    for bin in range(1, n_pt_bins_truth-3):
      min_pt_truth = truth_pt_bin_array[bin]
      max_pt_truth = truth_pt_bin_array[bin+1]
      
      self.plot_observable(observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, plot_pythia=False)
      self.plot_observable(observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, plot_pythia=True)

  #----------------------------------------------------------------------
  def plot_final_result_overlay(self, observable, jetR):
  
    # Plot overlay of different beta, for fixed pt bin
    if observable == 'theta_g':
      ymin_ratio = 0.
      ymax_ratio = 2.8
    if observable == 'zg':
      ymin_ratio = 0.3
      ymax_ratio = 1.6
    
    self.plot_observable_overlay_beta(observable, jetR, 20., 40., plot_pythia=True, plot_ratio = True, ymin_ratio=ymin_ratio, ymax_ratio=ymax_ratio)
    self.plot_observable_overlay_beta(observable, jetR, 40., 60., plot_pythia=True, plot_ratio = True, ymin_ratio=ymin_ratio, ymax_ratio=ymax_ratio)
    self.plot_observable_overlay_beta(observable, jetR, 60., 80., plot_pythia=True, plot_ratio = True, ymin_ratio=ymin_ratio, ymax_ratio=ymax_ratio)

  #----------------------------------------------------------------------
  def plot_observable_overlay_beta(self, observable, jetR, min_pt_truth, max_pt_truth, plot_pythia=False, plot_ratio=False, ymin_ratio=0., ymax_ratio=2.):
    
    name = 'cResult_overlay_R{}_allpt_{}-{}'.format(jetR, min_pt_truth, max_pt_truth)
    if plot_ratio:
      c = ROOT.TCanvas(name, name, 600, 650)
    else:
      c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    if plot_ratio:
      pad1 = ROOT.TPad('myPad', 'The pad',0,0.3,1,1)
    else:
      pad1 = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    pad1.SetLeftMargin(0.2)
    pad1.SetTopMargin(0.07)
    pad1.SetRightMargin(0.04)
    pad1.SetBottomMargin(0.13)
    if plot_ratio:
      pad1.SetBottomMargin(0.)
    pad1.Draw()
    pad1.cd()
    
    if observable == 'theta_g':
      xmin = 0.
      xmax = 1.
      ymax = 4.
    if observable == 'zg':
      xmin = 0.
      xmax = 0.5
      ymax = 11.
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', 1, xmin, xmax)
    myBlankHisto.SetNdivisions(505)
    xtitle = getattr(self, 'xtitle_{}'.format(observable))
    ytitle = getattr(self, 'ytitle_{}'.format(observable))
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    if plot_ratio:
      myBlankHisto.SetMinimum(2e-4) # Don't draw 0 on top panel
      myBlankHisto.GetYaxis().SetTitleSize(0.065)
      myBlankHisto.GetYaxis().SetTitleOffset(1.4)
      myBlankHisto.GetYaxis().SetLabelSize(0.06)
    myBlankHisto.Draw("E")

    # Plot ratio
    if plot_ratio:
      
      c.cd()
      pad2 = ROOT.TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
      pad2.SetTopMargin(0)
      pad2.SetBottomMargin(0.4)
      pad2.SetLeftMargin(0.2)
      pad2.SetRightMargin(0.04)
      pad2.Draw()
      pad2.cd()
      
      myBlankHisto2 = myBlankHisto.Clone("myBlankHisto_C")
      myBlankHisto2.SetYTitle("#frac{Data}{PYTHIA}")
      myBlankHisto2.SetXTitle(xtitle)
      myBlankHisto2.GetXaxis().SetTitleSize(30)
      myBlankHisto2.GetXaxis().SetTitleFont(43)
      myBlankHisto2.GetXaxis().SetTitleOffset(4.)
      myBlankHisto2.GetXaxis().SetLabelFont(43)
      myBlankHisto2.GetXaxis().SetLabelSize(25)
      myBlankHisto2.GetYaxis().SetTitleSize(20)
      myBlankHisto2.GetYaxis().SetTitleFont(43)
      myBlankHisto2.GetYaxis().SetTitleOffset(2.2)
      myBlankHisto2.GetYaxis().SetLabelFont(43)
      myBlankHisto2.GetYaxis().SetLabelSize(25)
      myBlankHisto2.GetYaxis().SetNdivisions(505)
      myBlankHisto2.GetYaxis().SetRangeUser(ymin_ratio, ymax_ratio)
      myBlankHisto2.Draw()
    
      line = ROOT.TLine(0,1,xmax,1)
      line.SetLineColor(920+2)
      line.SetLineStyle(2)
      line.Draw()

    pad1.cd()
    myLegend = ROOT.TLegend(0.66,0.65,0.8,0.85)
    self.utils.setup_legend(myLegend,0.035)
      
    # Retrieve histogram binnings for each SD setting
    for i, sd_setting in enumerate(self.sd_settings):
      
      zcut = sd_setting[0]
      beta = sd_setting[1]
      sd_label = 'zcut{}_B{}'.format(self.utils.remove_periods(zcut), beta)
      
      if observable == 'theta_g':
        n_rg_bins_truth = getattr(self, 'n_rg_bins_truth_{}'.format(sd_label))
        truth_rg_bin_array = getattr(self, 'truth_rg_bin_array_{}'.format(sd_label))
      if observable == 'zg':
        n_rg_bins_truth = getattr(self, 'n_zg_bins_truth_{}'.format(sd_label))
        truth_rg_bin_array = getattr(self, 'truth_zg_bin_array_{}'.format(sd_label))
      
      if i == 0:
        marker = 20
        marker_pythia = marker+4
        color = 600-6
      if i == 1:
        marker = 21
        marker_pythia = marker+4
        color = 632-4
      if i == 2:
        marker = 33
        marker_pythia = 27
        color = 416-2

      pad1.cd()
      if plot_pythia:
        
        fPythia_name = '/Users/jamesmulligan/alidock/theta_g/Pythia_new/pythia.root'
        fPythia = ROOT.TFile(fPythia_name, 'READ')
        hname = 'histogram_h_{}_B{}_{}-{}'.format(observable, beta, int(min_pt_truth), int(max_pt_truth))
        hPythia = fPythia.Get(hname)
        hPythia.SetDirectory(0)
        
        n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
        n_jets_tagged = hPythia.Integral(hPythia.FindBin(truth_rg_bin_array[0]), hPythia.GetNbinsX())
        fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
        hPythia.Scale(1./n_jets_inclusive, 'width')
        
        plot_errors = False
        if plot_errors:
          hPythia.SetMarkerSize(0)
          hPythia.SetMarkerStyle(0)
          hPythia.SetMarkerColor(color)
          hPythia.SetFillColor(color)
          hPythia.DrawCopy('E3 same')
        else:
          hPythia.SetLineColor(color)
          hPythia.SetLineColorAlpha(color, 0.5)
          hPythia.SetLineWidth(4)
          hPythia.DrawCopy('L hist same')

      h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth))
      h_sys.SetLineColor(0)
      h_sys.SetFillColor(color)
      h_sys.SetFillColorAlpha(color, 0.3)
      h_sys.SetFillStyle(1001)
      h_sys.SetLineWidth(0)
      h_sys.DrawCopy('E2 same')
    
      name = 'hMain_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
      fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
      h = getattr(self, name)
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(marker)
      h.SetMarkerColor(color)
      h.SetLineStyle(1)
      h.SetLineWidth(2)
      h.SetLineColor(color)
      h.DrawCopy('PE X0 same')
        
      myLegend.AddEntry(h, 'ALICE pp #beta={}'.format(beta), 'pe')

      if plot_ratio:
        pad2.cd()
        
        hRatioSys = h_sys.Clone()
        hRatioSys.SetName('{}_Ratio'.format(h_sys.GetName()))
        hRatioSys.Divide(hPythia)
        hRatioSys.SetLineColor(0)
        hRatioSys.SetFillColor(color)
        hRatioSys.SetFillColorAlpha(color, 0.3)
        hRatioSys.SetFillStyle(1001)
        hRatioSys.SetLineWidth(0)
        hRatioSys.DrawCopy('E2 same')

        hRatioStat = h.Clone()
        hRatioStat.SetName('{}_Ratio'.format(h.GetName()))
        hRatioStat.Divide(hPythia)
        hRatioStat.SetMarkerSize(1.5)
        hRatioStat.SetMarkerStyle(marker)
        hRatioStat.SetMarkerColor(color)
        hRatioStat.SetLineStyle(1)
        hRatioStat.SetLineWidth(2)
        hRatioStat.SetLineColor(color)
        hRatioStat.DrawCopy('PE X0 same')
        
    pad1.cd()
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      myLegend.AddEntry(hPythia, 'PYTHIA8 Monash2013', 'l')
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE Preliminary'
    text_latex.DrawLatex(0.25, 0.87, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.81, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'Charged jets   anti-#it{k}_{T}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.75, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{R} = ' + str(jetR) + '  |#it{#eta}_{jet}| < 0.5' + '  #it{z}_{cut} = ' + str(zcut)
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.69, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.25, 0.63, text)
    
    myLegend.Draw()

    name = 'h_{}_R{}_{}-{}{}'.format(observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'h_{}_R{}_{}-{}_Pythia{}'.format(observable, self.utils.remove_periods(jetR), int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results_{}'.format(observable))
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

  #----------------------------------------------------------------------
  def plot_observable(self, observable, jetR, sd_label, zcut, beta, min_pt_truth, max_pt_truth, plot_pythia=False):
    
    name = 'cResult_R{}_{}_{}-{}'.format(jetR, sd_label, min_pt_truth, max_pt_truth)
    c = ROOT.TCanvas(name, name, 600, 450)
    c.Draw()
    
    c.cd()
    myPad = ROOT.TPad('myPad', 'The pad',0,0,1,1)
    myPad.SetLeftMargin(0.2)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.04)
    myPad.SetBottomMargin(0.13)
    myPad.Draw()
    myPad.cd()
    
    if observable == 'theta_g':
      n_rg_bins_truth = getattr(self, 'n_rg_bins_truth_{}'.format(sd_label))
      truth_rg_bin_array = getattr(self, 'truth_rg_bin_array_{}'.format(sd_label))
      ymax = 5.
    if observable == 'zg':
      n_rg_bins_truth = getattr(self, 'n_zg_bins_truth_{}'.format(sd_label))
      truth_rg_bin_array = getattr(self, 'truth_zg_bin_array_{}'.format(sd_label))
      ymax = 12.
    myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_rg_bins_truth, truth_rg_bin_array)
    myBlankHisto.SetNdivisions(505)
    xtitle = getattr(self, 'xtitle_{}'.format(observable))
    ytitle = getattr(self, 'ytitle_{}'.format(observable))
    myBlankHisto.SetXTitle(xtitle)
    myBlankHisto.GetYaxis().SetTitleOffset(1.5)
    myBlankHisto.SetYTitle(ytitle)
    myBlankHisto.SetMaximum(ymax)
    myBlankHisto.SetMinimum(0.)
    myBlankHisto.Draw("E")

    if plot_pythia:
    
      plot_pythia_from_response = False
      plot_pythia_from_mateusz = True
      if plot_pythia_from_response:
        hPythia = self.get_pythia_from_response(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
        n_jets_inclusive = hPythia.Integral(0, hPythia.GetNbinsX()+1)
        n_jets_tagged = hPythia.Integral(hPythia.FindBin(truth_rg_bin_array[0]), hPythia.GetNbinsX())
        fraction_tagged_pythia =  n_jets_tagged/n_jets_inclusive
        hPythia.Scale(1./n_jets_inclusive, 'width')
        hPythia.SetFillStyle(0)
        hPythia.SetMarkerSize(1.5)
        hPythia.SetMarkerStyle(21)
        hPythia.SetMarkerColor(2)
        hPythia.SetLineColor(2)
        hPythia.SetLineWidth(1)
        hPythia.Draw('E2 same')
    
      if plot_pythia_from_mateusz:
        fPythia_name = '/Users/jamesmulligan/alidock/theta_g/Pythia_new/pythia.root'
        fPythia = ROOT.TFile(fPythia_name, 'READ')
        hname = 'histogram_h_{}_B{}_{}-{}'.format(observable, beta, int(min_pt_truth), int(max_pt_truth))
        hPythia2 = fPythia.Get(hname)

        n_jets_inclusive2 = hPythia2.Integral(0, hPythia2.GetNbinsX()+1)
        n_jets_tagged2 = hPythia2.Integral(hPythia2.FindBin(truth_rg_bin_array[0]), hPythia2.GetNbinsX())
        fraction_tagged_pythia =  n_jets_tagged2/n_jets_inclusive2
        hPythia2.Scale(1./n_jets_inclusive2, 'width')
        hPythia2.SetFillStyle(0)
        hPythia2.SetMarkerSize(1.5)
        hPythia2.SetMarkerStyle(21)
        hPythia2.SetMarkerColor(1)
        hPythia2.SetLineColor(1)
        hPythia2.SetLineWidth(1)
        hPythia2.Draw('E2 same')
    
    color = 600-6
    h_sys = getattr(self, 'hResult_{}_systotal_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth))
    h_sys.SetLineColor(0)
    h_sys.SetFillColor(color)
    h_sys.SetFillColorAlpha(color, 0.3)
    h_sys.SetFillStyle(1001)
    h_sys.SetLineWidth(0)
    h_sys.DrawCopy('E2 same')
    
    name = 'hMain_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth)
    fraction_tagged = getattr(self, '{}_fraction_tagged'.format(name))
    h = getattr(self, name)
    h.SetMarkerSize(1.5)
    h.SetMarkerStyle(20)
    h.SetMarkerColor(color)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(color)
    h.DrawCopy('PE X0 same')
  
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'ALICE Preliminary'
    text_latex.DrawLatex(0.57, 0.87, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'pp #sqrt{#it{s}} = 5.02 TeV'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.8, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = 'R = ' + str(jetR) + '   | #eta_{jet}| < 0.5'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.73, text)

    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.66, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text = '#it{z}_{cut} = ' + str(zcut) + '  #beta = ' + str(beta)
    text_latex.SetTextSize(0.045)
    text_latex.DrawLatex(0.57, 0.59, text)
    
    text_latex = ROOT.TLatex()
    text_latex.SetNDC()
    text_latex.SetTextSize(0.04)
    text = '#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged
    text_latex.DrawLatex(0.57, 0.52, text)
    
    if plot_pythia:
      text_latex = ROOT.TLatex()
      text_latex.SetNDC()
      text_latex.SetTextSize(0.04)
      text = ('#it{f}_{tagged}^{data} = %3.3f' % fraction_tagged) + (', #it{f}_{tagged}^{pythia} = %3.3f' % fraction_tagged_pythia)
      text_latex.DrawLatex(0.57, 0.52, text)

    myLegend = ROOT.TLegend(0.25,0.7,0.5,0.85)
    self.utils.setup_legend(myLegend,0.035)
    myLegend.AddEntry(h, 'ALICE pp', 'pe')
    myLegend.AddEntry(h_sys, 'Sys. uncertainty', 'f')
    if plot_pythia:
      if plot_pythia_from_response:
        myLegend.AddEntry(hPythia, 'PYTHIA Monash2013', 'pe')
      if plot_pythia_from_mateusz:
        myLegend.AddEntry(hPythia2, 'PYTHIA Monash2013', 'pe')
    myLegend.Draw()

    name = 'hUnfolded_R{}_{}_{}-{}{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    if plot_pythia:
      name = 'hUnfolded_R{}_{}_{}-{}_Pythia{}'.format(self.utils.remove_periods(jetR), sd_label, int(min_pt_truth), int(max_pt_truth), self.file_format)
    output_dir = getattr(self, 'output_dir_final_results_{}'.format(observable))
    outputFilename = os.path.join(output_dir, name)
    c.SaveAs(outputFilename)
    c.Close()

    # Write result to ROOT file
    if not plot_pythia:
      final_result_root_filename = os.path.join(output_dir, 'fFinalResults.root')
      fFinalResults = ROOT.TFile(final_result_root_filename, 'UPDATE')
      h.Write()
      fFinalResults.Close()


  #----------------------------------------------------------------------
  def get_pythia_from_response(self, observable, jetR, sd_label, min_pt_truth, max_pt_truth):
  
    output_dir = getattr(self, 'output_dir_main_{}'.format(observable))
    file = os.path.join(output_dir, 'response.root')
    f = ROOT.TFile(file, 'READ')

    if observable == 'theta_g':
      thn_name = 'hResponse_JetPt_ThetaG_R{}_{}_rebinned'.format(jetR, sd_label)
    if observable == 'zg':
      thn_name = 'hResponse_JetPt_zg_R{}_{}_rebinned'.format(jetR, sd_label)
    
    thn = f.Get(thn_name)
    thn.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)

    h = thn.Projection(3)
    h.SetName('hPythia_{}_R{}_{}_{}-{}'.format(observable, jetR, sd_label, min_pt_truth, max_pt_truth))
    h.SetDirectory(0)

    return h

  #----------------------------------------------------------------------
  def change_to_per(self, h):
    
    for bin in range(0, h.GetNbinsX()):
      content = h.GetBinContent(bin)
      content_new = math.fabs(1-content)
      h.SetBinContent(bin, content_new*100)

  #----------------------------------------------------------------------
  def build_average(self, h1, h2, takeMaxDev=False):
  
    h_avg = h1.Clone()
    h_avg.SetName('{}_avg'.format(h1.GetName()))
  
    for i in range(1, h_avg.GetNbinsX()+1):
      value1 = h1.GetBinContent(i)
      value2 = h2.GetBinContent(i)
      avg =  0.5*(value1 + value2)
      
      if takeMaxDev:
        if value1>value2:
          avg = value1
        else:
          avg = value2
    
      h_avg.SetBinContent(i, avg)
  
    return h_avg

  #----------------------------------------------------------------------
  def AttachErrToHist(self, h, hPercError):
  
    #Fill array with lower bin edges of data histogram
    for bin in range(1, h.GetNbinsX()+1):
      content = h.GetBinContent(bin)
      perErr = hPercError.GetBinContent(bin)
      h.SetBinError(bin, content*perErr*0.01)

#----------------------------------------------------------------------
if __name__ == '__main__':

  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold theta_g distribution')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')

  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('configFile: \'{0}\''.format(args.configFile))
  
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  analysis = run_sd_analysis(config_file = args.configFile)
  analysis.run_sd_analysis()
