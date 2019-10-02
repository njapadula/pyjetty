#! /usr/bin/env python

import sys
import os
import argparse
import itertools
from array import *
import numpy
import ROOT
import yaml

# Load the RooUnfold library
ROOT.gSystem.Load("$ALIBUILD_WORK_DIR/slc7_x86-64/RooUnfold/latest/lib/libRooUnfold.so")

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

# Suppress a lot of standard output
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Set plotting options
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

analysis_utils = analysis_utils.analysis_utils()

###########################################################################################
###########################################################################################
def roounfold_rg(input_file_data, input_file_response, config_file, output_dir, file_format):

  fData = ROOT.TFile(input_file_data)
  fResponse = ROOT.TFile(input_file_response)
  
  # Create output dir for unfolding histograms and result
  if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
  
  # Read config file
  with open(config_file, 'r') as stream:
    config = yaml.safe_load(stream)
  
  jetR_list = config['jetR']
  beta_dict = config['beta']
  beta_list = list(beta_dict.keys())
  reg_param_final = config['reg_param']

  min_pt_reported = 20
  max_pt_reported = 80
  
  #--------------------------------------------------------------

  # Define an empty dictionary to store the final unfolded histograms from each label
  unfolded_dict = {}
  
  for jetR in jetR_list:
    for beta in beta_list:
      
      unfoldSingleOutputList(fData, fResponse, unfolded_dict, jetR, beta, beta_dict, reg_param_final, min_pt_reported, max_pt_reported, output_dir, file_format)

###################################################################################################
# Unfold jet spectrum from a single output list
###################################################################################################
def unfoldSingleOutputList(fData, fResponse, unfolded_dict, jetR, beta, beta_dict, reg_param_final, min_pt_reported, max_pt_reported, output_dir, file_format):

  print('jetR = {},  beta = {}'.format(jetR, beta))
  
  # Get data jet spectrum
  name = 'hThetaG_JetPt_R{}_B{}_Rebinned'.format(jetR, beta)
  hData_PerBin = fData.Get(name)
  hData_PerBin.Sumw2()
  hData_PerBin.GetXaxis().SetRangeUser(0., 100.)
  outputFilename = os.path.join(output_dir, 'hData_R{}_B{}{}'.format(analysis_utils.remove_periods(jetR), beta, file_format))
  plotHist(hData_PerBin, outputFilename, 'colz', False, True)
  
  # Get RooUnfoldResponse object
  name = 'roounfold_response_R{}_B{}'.format(jetR, beta)
  response = fResponse.Get(name)
  response.UseOverflow(False)
  
  # Plot various slices of the response matrix (from the THn)
  plot_RM_slices(fResponse, jetR, beta, beta_dict, output_dir, file_format)
  
  # Plot the kinematic efficiency from the response THn
  plot_kinematic_efficiency(fResponse, jetR, beta, beta_dict, output_dir, file_format)
  # Can either pass this to the RooUnfoldResponse object, or else can apply it afterward
  # (Need to think...)
  
  # Get MC-det and MC-truth 2D projections for unfolding closure test
  hMC_Det = get_MCdet2D(fResponse, jetR, beta)
  hMC_Truth = get_MCtruth2D(fResponse, jetR, beta)

  # Set prior
  # Get the (matched) truth-level jet spectrum from the THn (for prior)
  #hPrior = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  # Normalize the response matrix by setting each pT-truth projection to the intended prior distribution
  # Keep response matrix as per-bin probabilities (i.e. don't scale by bin width)
  # Scale also hJetSpectrumTrueUncutPerBin accordingly, soas to preserve the kinematic efficiency
  #setPrior(hResponseMatrix, hJetSpectrumTrueUncutPerBin, outputDir, fileFormat)
  #
  # (Need to think how to modify the prior...)

  # Unfold spectrum
  if hData_PerBin and response and hMC_Det and hMC_Truth:
    
    unfoldJetSpectrum(hData_PerBin, response, unfolded_dict, jetR, beta, beta_dict, reg_param_final, min_pt_reported, max_pt_reported, hMC_Det, hMC_Truth, output_dir, file_format)

#################################################################################################
# Unfold jet spectrum
#################################################################################################
def unfoldJetSpectrum(hData_PerBin, response, unfolded_dict, jetR, beta, beta_dict, reg_param_final, min_pt_reported, max_pt_reported, hMC_Det, hMC_Truth, output_dir, file_format):
  
  regularizationParamName = "n_iter"
  
  # Create canvas to superpose iterations on one plot, to examine convergence
  c1 = ROOT.TCanvas('c1_{}_{}'.format(jetR, beta),'c1_{}_{}: histos'.format(jetR, beta),600,450)
  c1.cd()
  c1.SetLogy()
  ROOT.gPad.SetLeftMargin(0.15)
  leg = ROOT.TLegend(0.6,0.5,0.88,0.83,"{} Unfolding".format(type))
  leg.SetFillColor(10)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  
  # Loop over values of regularization parameter
  for i in range(1, reg_param_final + 3):
    
    # Set up the Bayesian unfolding object
    unfoldBayes = ROOT.RooUnfoldBayes(response, hData_PerBin, i)
    #unfoldBayes.SetNToys(1000)
    
    # Perform the unfolding
    print('Unfolding with {} = {}'.format(regularizationParamName, i))
    errorType = ROOT.RooUnfold.kCovToy
    hUnfolded = unfoldBayes.Hreco(errorType) # Produces the truth distribution, with errors, PerBin (will scale by bin width below, after refolding checks)
    hUnfolded.SetName('hUnfolded_R{}_B{}_{}'.format(jetR, beta, i))
    
    # Plot Pearson correlation coeffs for each k, to get a measure of the correlation between the bins
    covarianceMatrix = unfoldBayes.Ereco(errorType) # Get the covariance matrix
    #plotCorrelationCoefficients(covarianceMatrix, i, output_dir, file_format)
    
    unfolded_dict[i] = hUnfolded

  plot_unfolded_rg(unfolded_dict, jetR, beta, beta_dict, reg_param_final, regularizationParamName, min_pt_reported, max_pt_reported, output_dir, file_format)

  plot_unfolded_pt(unfolded_dict, jetR, beta, beta_dict, reg_param_final, regularizationParamName, min_pt_reported, max_pt_reported, output_dir, file_format)

  #--------------------------------------------------------------
  
  # Smear MC truth spectrum according to the error bars on the measured spectrum, for closure test
  measuredErrors = getMeasuredErrors(hData_PerBin)
  smearSpectrum(hMC_Det, measuredErrors)

  # Loop over values of regularization parameter to do unfolding checks
  for i in range(1, reg_param_final + 3):
    
    # Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
    plot_refolding_test(response, unfolded_dict[i], hData_PerBin, i, regularizationParamName, jetR, beta, beta_dict, output_dir, file_format)

    # Unfold the smeared det-level result with response, and compare to truth-level MC.
    plot_closure_test(response, hMC_Det, hMC_Truth, i, regularizationParamName, jetR, beta, beta_dict, output_dir, file_format)

#################################################################################################
# Plot various slices of the response matrix (from the THn)
#################################################################################################
def plotRegParamSystematic():
  
  #Plot the spectra comparing only k=+1 and k-1 to the main result
  xRangeMin = min_pt_reported
  xRangeMax = max_pt_reported
  yAxisTitle = "#frac{d#sigma}{dp_{T}} [mb (GeV/c)^{-1}]"
  ratioYAxisTitle = "Ratio to k={}".format(reg_param)
  outputFilename = os.path.join(output_dir, "hJetSpectraUnfoldedRatio" + file_format)
  legendTitle = "{} Unfolding".format(type)
  hLegendLabel = "k = {}".format(reg_param-1)
  h2LegendLabel = "k = {}".format(reg_param)
  h3LegendLabel = "k = {}".format(reg_param+1)
  if hLowerkResult and hMainResult and hHigherkResult:
    # To get sensible error bars, assume main result has no errors
    for bin in range(1, hMainResult.GetNbinsX() + 1):
      hMainResult.SetBinError(bin, 0)
    plotSpectra(hLowerkResult, hMainResult, hHigherkResult, 1., xRangeMin, xRangeMax, yAxisTitle, ratioYAxisTitle, outputFilename, "", legendTitle, hLegendLabel, h2LegendLabel, h3LegendLabel)

##################################################################################################
# Set prior in repsonse matrix
# RooUnfold takes the prior as the truth-axis projection of the response matrix.
# After normalizing the response matrix to our desired prior, RooUnfold will then
# automatically normalize the response matrix to conserve the number of jets (i.e. each bin of
# truth-axis projection normalized to 1).
##################################################################################################
def setPrior(hResponseMatrix, hJetSpectrumTrueUncutPerBin, outputDir, fileFormat):
  
  priorFolder = os.path.join(outputDir, "PriorProperties/")
  if not os.path.exists(priorFolder):
    os.makedirs(priorFolder)
  
  # Plot response matrix before normalization
  outputFilename = os.path.join(priorFolder, "hResponseMatrixBeforeNormalization" + fileFormat)
  plotHist(hResponseMatrix, outputFilename, "colz")
  
  # Make projection onto pT-true axis (y-axis), and scale appropriately
  hTruthProjectionBefore = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hTruthProjectionBefore.SetName("hTruthProjectionBefore")

  #hResponseMatrix has a detector min pT restriction
  outputFilename = os.path.join(priorFolder, "hTruthProjectionBeforeNormalization" + fileFormat)
  plotHist(hTruthProjectionBefore, outputFilename, "hist", True)
  
  #hJetSpectrumTrueUncutPerBin has no detector pT restriction
  outputFilename = os.path.join(priorFolder, "hJetSpectrumTrueUncutPerBin_asPrior{}".format(fileFormat))
  plotHist(hJetSpectrumTrueUncutPerBin, outputFilename,"hist", True)
  
  # Check kinematic efficiency before
  hKinematicEfficiencyBefore = hTruthProjectionBefore.Clone()
  hKinematicEfficiencyBefore.SetName("hKinematicEfficiencyBefore")
  hKinematicEfficiencyBefore.Divide(hTruthProjectionBefore, hJetSpectrumTrueUncutPerBin, 1., 1., "B")
  hKinematicEfficiencyBefore.SetMarkerStyle(21)
  outputFilename = os.path.join(priorFolder, "hKinematicEfficiencyBefore{}".format(fileFormat))
  plotHist(hKinematicEfficiencyBefore, outputFilename, "P E")
  
  # Loop through truth-level bins (of RM and uncut truth spectrum), and normalize each bin (of RM and uncut truth spectrum) to the prior
  nBinsY = hResponseMatrix.GetNbinsY() # pT-gen
  nBinsX = hResponseMatrix.GetNbinsX() # pT-det

  # Plot truth projection after normalization (i.e. prior distribution)
  hTruthProjectionAfter = hResponseMatrix.ProjectionY("_py",1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hTruthProjectionAfter.SetName("hTruthProjectionAfterNormalization")
  outputFilename = os.path.join(priorFolder, "hTruthProjectionAfterNormalization" + fileFormat)
  plotHist(hTruthProjectionAfter, outputFilename, "hist", True)
  
  # Plot response matrix after normalization
  outputFilename = os.path.join(priorFolder, "hResponseMatrixAfterNormalization" + fileFormat)
  plotHist(hResponseMatrix, outputFilename, "colz")
  
  #Prior used for the unfolding
  outputFilename = os.path.join(priorFolder, "hPriorForUnfolding{}".format(fileFormat))
  plotHist(hJetSpectrumTrueUncutPerBin, outputFilename,"hist", True)
  
  # Check kinematic efficiency after normalization
  hKinematicEfficiencyAfter = hTruthProjectionAfter.Clone()
  hKinematicEfficiencyAfter.SetName("hKinematicEfficiencyAfter")
  hKinematicEfficiencyAfter.Divide(hTruthProjectionAfter, hJetSpectrumTrueUncutPerBin, 1., 1., "B")
  hKinematicEfficiencyAfter.SetMarkerStyle(21)
  outputFilename = os.path.join(priorFolder, "hKinematicEfficiencyAfter{}".format(fileFormat))
  plotHist(hKinematicEfficiencyAfter, outputFilename, "P E")

#################################################################################################
# Plot various slices of the response matrix (from the THn)
#################################################################################################
def plot_unfolded_rg(unfolded_dict, jetR, beta, beta_dict, reg_param_final, regularizationParamName, min_pt_reported, max_pt_reported, output_dir, file_format):
  
  binning_dict = beta_dict[beta]
  pt_bins_truth = (binning_dict['pt_bins_truth'])
  rg_bins_truth = (binning_dict['rg_bins_truth'])
  n_pt_bins_truth = len(pt_bins_truth) - 1
  n_rg_bins_truth = len(rg_bins_truth) - 1
  truth_pt_bin_array = array('d',pt_bins_truth)
  truth_rg_bin_array = array('d',rg_bins_truth)
  
  for bin in range(1, n_pt_bins_truth-1):
    min_pt_truth = pt_bins_truth[bin]
    max_pt_truth = pt_bins_truth[bin+1]
    
    plot_rg(unfolded_dict, jetR, beta, reg_param_final, regularizationParamName, min_pt_truth, max_pt_truth, n_rg_bins_truth, truth_rg_bin_array, output_dir, file_format)

#################################################################################################
# Plot various slices of the response matrix (from the THn)
#################################################################################################
def plot_rg(unfolded_dict, jetR, beta, reg_param_final, regularizationParamName, min_pt_truth, max_pt_truth, n_rg_bins_truth, truth_rg_bin_array, output_dir, file_format):
  
  setOptions()
  ROOT.gROOT.ForceStyle()
  
  name = 'cResult_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth)
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
  
  myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_rg_bins_truth, truth_rg_bin_array)
  myBlankHisto.SetNdivisions(505)
  myBlankHisto.SetXTitle('#theta_{g}')
  myBlankHisto.GetYaxis().SetTitleOffset(1.5)
  myBlankHisto.SetYTitle('#frac{1}{#it{N}_{jets, inc}} #frac{d#it{N}}{d#theta_{g}}')
  myBlankHisto.SetMaximum(3)
  myBlankHisto.SetMinimum(0.)
  myBlankHisto.Draw("E")
  
  leg = ROOT.TLegend(0.75,0.65,0.88,0.92)
  setupLegend(leg,0.04)
  
  for i in range(1, reg_param_final + 3):
    
    h2D = unfolded_dict[i]
    h2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    #h = h2D.ProjectionY('{}_py'.format(h2D.GetName()), 1, h2D.GetNbinsX()) # exclude under- and over-flow bins
    h = h2D.ProjectionY()
    
    # Normalize by integral, i.e. N_jets,inclusive in this pt-bin (cross-check this)
    integral = h.Integral()
    h.Scale(1./integral, 'width')
    
    if i == 1:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(20)
    elif i == 2:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(21)
    elif i == 3:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(22)
    elif i == 4:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(23)
    elif i == 5:
      h.SetMarkerSize(2)
      h.SetMarkerStyle(33)
    elif i == 6:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(34)
    elif i == 7:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(35)
    elif i == 8:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(36)
    else:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(19)
    h.SetMarkerColor(600-6+i)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(600-6+i)
    
    h.DrawCopy('PE X0 same')
    
    label = '{} = {}'.format(regularizationParamName, i)
    leg.AddEntry(h, label, 'Pe')

  leg.Draw()

  text_latex = ROOT.TLatex()
  text_latex.SetNDC()
  text = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth)
  text_latex.DrawLatex(0.45, 0.85, text)

  text_latex = ROOT.TLatex()
  text_latex.SetNDC()
  text = 'R = ' + str(jetR) + '   #beta = ' + str(beta)
  text_latex.DrawLatex(0.45, 0.75, text)
  
  outputFilename = os.path.join(output_dir, 'hUnfolded_R{}_B{}_{}-{}{}'.format(analysis_utils.remove_periods(jetR), beta, min_pt_truth, max_pt_truth, file_format))
  c.SaveAs(outputFilename)
  c.Close()

#################################################################################################
# Plot various slices of the response matrix (from the THn)
#################################################################################################
def plot_unfolded_pt(unfolded_dict, jetR, beta, beta_dict, reg_param_final, regularizationParamName, min_pt_reported, max_pt_reported, output_dir, file_format):
  
  binning_dict = beta_dict[beta]
  pt_bins_truth = (binning_dict['pt_bins_truth'])
  n_pt_bins_truth = len(pt_bins_truth) - 1
  truth_pt_bin_array = array('d',pt_bins_truth)
  
  setOptions()
  ROOT.gROOT.ForceStyle()
  
  name = 'cResultPt_R{}_B{}'.format(jetR, beta)
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
  
  myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_pt_bins_truth, truth_pt_bin_array)
  myBlankHisto.SetNdivisions(505)
  myBlankHisto.SetXTitle('#it{p}_{T, ch jet}')
  myBlankHisto.GetYaxis().SetTitleOffset(2.2)
  myBlankHisto.SetYTitle('#frac{dN}{d#it{p}_{T, ch jet}}')
  myBlankHisto.SetMaximum(5000)
  myBlankHisto.SetMinimum(0.)
  myBlankHisto.Draw("E")
  
  leg = ROOT.TLegend(0.75,0.65,0.88,0.92)
  setupLegend(leg,0.04)
  
  for i in range(1, reg_param_final + 3):
    
    h2D = unfolded_dict[i]
    h2D.GetXaxis().SetRangeUser(5., 120.)
    h = h2D.ProjectionX()
    
    h.Scale(1., 'width')
    
    if i == 1:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(20)
    elif i == 2:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(21)
    elif i == 3:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(22)
    elif i == 4:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(23)
    elif i == 5:
      h.SetMarkerSize(2)
      h.SetMarkerStyle(33)
    elif i == 6:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(34)
    elif i == 7:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(35)
    elif i == 8:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(36)
    else:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(19)
    h.SetMarkerColor(600-6+i)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(600-6+i)
    
    h.DrawCopy('PE X0 same')
    
    label = '{} = {}'.format(regularizationParamName, i)
    leg.AddEntry(h, label, 'Pe')
  
  leg.Draw()

  text_latex = ROOT.TLatex()
  text_latex.SetNDC()
  text = 'R = ' + str(jetR) + '   #beta = ' + str(beta)
  text_latex.DrawLatex(0.45, 0.75, text)

  outputFilename = os.path.join(output_dir, 'hUnfoldedPt_R{}_B{}{}'.format(analysis_utils.remove_periods(jetR), beta, file_format))
  c.SaveAs(outputFilename)
  c.Close()

###################################################################################################
# Plot kinematic efficiency
# The kinematic efficiency is the ratio:
#   Numerator: 2D truth-level projection [pt-true, rg-true] using no cut on det-level
#   Denominator: 2D truth-level projection [pt-true, rg-true] using [pt-det, rg-det] cut on det-level
###################################################################################################
def plot_kinematic_efficiency(fResponse, jetR, beta, beta_dict, output_dir, file_format):
  
  # (pt-det, pt-true, theta_g-det, theta_g-true)
  name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
  hResponse = fResponse.Get(name).Clone()
  hResponse.SetName('hResponse_KinEff_JetPt_ThetaG_R{}_B{}'.format(jetR, beta))
  
  # Denominator -- by default, under/over-flow bins are included in projection
  hDenominator = hResponse.Projection(3, 1)
  hDenominator.SetName('{}_Denominator'.format(hDenominator.GetName()))
  
  # Numerator -- cut on det-level input binning
  binning_dict = beta_dict[beta]
  pt_bins_det = (binning_dict['pt_bins_det'])
  rg_bins_det = (binning_dict['rg_bins_det'])
  min_pt_det = pt_bins_det[0]
  max_pt_det = pt_bins_det[-1]
  min_rg_det = rg_bins_det[0]
  max_rg_det = rg_bins_det[-1]
  hResponse.GetAxis(0).SetRangeUser(min_pt_det, max_pt_det)
  hResponse.GetAxis(2).SetRangeUser(min_rg_det, max_rg_det)
  hNumerator = hResponse.Projection(3, 1)
  hNumerator.SetName('{}_Numerator'.format(hNumerator.GetName()))
  
  hKinematicEfficiency = hNumerator.Clone()
  hKinematicEfficiency.SetName('hKinematicEfficiency_R{}_B{}'.format(jetR, beta))
  hKinematicEfficiency.Divide(hDenominator)
  outputFilename = os.path.join(output_dir, 'hKinematicEfficiency2D_R{}_B{}{}'.format(analysis_utils.remove_periods(jetR), beta, file_format))
  plotHist(hKinematicEfficiency, outputFilename, "colz")
  
  plot_kinematic_efficiency_projections(hKinematicEfficiency, jetR, beta, beta_dict, output_dir, file_format)

#################################################################################################
# Unfold jet spectrum
#################################################################################################
def plot_kinematic_efficiency_projections(hKinematicEfficiency2D, jetR, beta, beta_dict, output_dir, file_format):
  
  binning_dict = beta_dict[beta]
  pt_bins_truth = (binning_dict['pt_bins_truth'])
  rg_bins_truth = (binning_dict['rg_bins_truth'])
  n_pt_bins_truth = len(pt_bins_truth) - 1
  n_rg_bins_truth = len(rg_bins_truth) - 1
  truth_pt_bin_array = array('d',pt_bins_truth)
  truth_rg_bin_array = array('d',rg_bins_truth)
  
  setOptions()
  ROOT.gROOT.ForceStyle()
  
  name = 'cKinEff_R{}_B{}'.format(jetR, beta)
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
  
  myBlankHisto = ROOT.TH1F('myBlankHisto','Blank Histogram', n_rg_bins_truth, truth_rg_bin_array)
  myBlankHisto.SetNdivisions(505)
  myBlankHisto.SetXTitle('#theta_{g}')
  myBlankHisto.GetYaxis().SetTitleOffset(1.5)
  myBlankHisto.SetYTitle('#varepsilon_{kin}')
  myBlankHisto.SetMaximum(2)
  myBlankHisto.SetMinimum(0.)
  myBlankHisto.Draw("E")
  
  leg = ROOT.TLegend(0.6,0.65,0.72,0.92)
  setupLegend(leg,0.04)
  
  for bin in range(1, n_pt_bins_truth-1):
    min_pt_truth = pt_bins_truth[bin]
    max_pt_truth = pt_bins_truth[bin+1]
    
    hKinematicEfficiency2D.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
    h = hKinematicEfficiency2D.ProjectionY()
    h.SetName('hKinematicEfficiency_R{}_B{}_{}-{}'.format(jetR, beta, min_pt_truth, max_pt_truth))
    
    if bin == 1:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(20)
    elif bin == 2:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(21)
    elif bin == 3:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(22)
    elif bin == 4:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(23)
    elif bin == 5:
      h.SetMarkerSize(2)
      h.SetMarkerStyle(33)
    else:
      h.SetMarkerSize(1.5)
      h.SetMarkerStyle(19)
    h.SetMarkerColor(600-6+bin)
    h.SetLineStyle(1)
    h.SetLineWidth(2)
    h.SetLineColor(600-6+bin)
    
    h.DrawCopy('P X0 same')
    
    label = str(min_pt_truth) + ' < #it{p}_{T, ch jet} < ' + str(max_pt_truth) + ' GeV/#it{c}'
    leg.AddEntry(h, label, 'Pe')
  
  leg.Draw()

  text_latex = ROOT.TLatex()
  text_latex.SetNDC()
  text = 'R = ' + str(jetR) + '   #beta = ' + str(beta)
  text_latex.DrawLatex(0.3, 0.85, text)

  outputFilename = os.path.join(output_dir, 'hKinematicEfficiency_R{}_B{}{}'.format(analysis_utils.remove_periods(jetR), beta, file_format))
  c.SaveAs(outputFilename)
  c.Close()

#################################################################################################
# Plot various slices of the response matrix (from the THn)
#################################################################################################
def plot_RM_slices(fResponse, jetR, beta, beta_dict, output_dir, file_format):
  
  output_dir_RM = os.path.join(output_dir, 'RM')
  if not os.path.isdir(output_dir_RM):
    os.makedirs(output_dir_RM)
  
  # (pt-det, pt-true, theta_g-det, theta_g-true)
  name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
  hResponse = fResponse.Get(name)
  
  # Fix pt-true, and plot the 2D theta_g response
  binning_dict = beta_dict[beta]
  pt_bins_truth = binning_dict['pt_bins_truth']
  n_pt_bins_truth = len(pt_bins_truth) - 1
  
  for bin in range(1, n_pt_bins_truth-1):
    min_pt_truth = pt_bins_truth[bin]
    max_pt_truth = pt_bins_truth[bin+1]
    
    plot_ThetaG_Response(jetR, beta, min_pt_truth, max_pt_truth, hResponse, output_dir_RM, file_format)

#################################################################################################
# Plot 2D theta_g response for a fixed range of pt-truth
#################################################################################################
def plot_ThetaG_Response(jetR, beta, min_pt_truth, max_pt_truth, hResponse, output_dir, file_format):
  
  hResponse4D = hResponse.Clone()
  hResponse4D.SetName('{}_{}_{}'.format(hResponse4D.GetName(), min_pt_truth, max_pt_truth))
  
  hResponse4D.GetAxis(1).SetRangeUser(min_pt_truth, max_pt_truth)
  hResponse_ThetaG = hResponse4D.Projection(3,2)
  hResponse_ThetaG.SetName('hResponse_ThetaG_R{}_B{}_{}_{}'.format(jetR, beta, min_pt_truth, max_pt_truth))
  
  hResponse_ThetaG_Normalized = normalizeResponseMatrix(hResponse_ThetaG)
  
  text = str(min_pt_truth) + ' < #it{p}_{T, ch jet}^{true} < ' + str(max_pt_truth)
  
  outputFilename = os.path.join(output_dir, '{}{}'.format(hResponse_ThetaG.GetName(), file_format))
  plotHist(hResponse_ThetaG_Normalized, outputFilename, 'colz', False, True, text)

################################################################################################
# Normalize response matrix
# Normalize the truth projection to 1
################################################################################################
def normalizeResponseMatrix(hResponseMatrix):
  
  # Make projection onto true axis (y-axis), and scale appropriately
  hTruthProjectionBefore = hResponseMatrix.ProjectionY('{}_py'.format(hResponseMatrix.GetName()),1,hResponseMatrix.GetNbinsX()) # Do exclude under and overflow bins
  hTruthProjectionBefore.SetName("hTruthProjectionBefore")
  
  # Loop through truth-level bins, and apply normalization factor to all bins.
  nBinsY = hResponseMatrix.GetNbinsY() # truth
  nBinsX = hResponseMatrix.GetNbinsX() # det
  for truthBin in range(1,nBinsY+1):
    normalizationFactor = hTruthProjectionBefore.GetBinContent(truthBin)
    if normalizationFactor > 0:
      truthBinCenter = hTruthProjectionBefore.GetXaxis().GetBinCenter(truthBin)
      
      for detBin in range(1,nBinsX+1):
        binContent = hResponseMatrix.GetBinContent(detBin, truthBin)
        hResponseMatrix.SetBinContent(detBin, truthBin, binContent/normalizationFactor)

  return hResponseMatrix

#################################################################################################
# Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
#################################################################################################
def plot_refolding_test(response, hUnfolded, hData_PerBin, i, regularizationParamName, jetR, beta, beta_dict, output_dir, file_format):
  
  # In principle should use two statistically independent response matrices
  
  output_dir_refolding = os.path.join(output_dir, 'Refolding')
  if not os.path.isdir(output_dir_refolding):
    os.makedirs(output_dir_refolding)
  
  hFoldedTruth = response.ApplyToTruth(hUnfolded) # Produces folded distribution PerBin (unfolded spectrum is also PerBin at the moment)
  
  binning_dict = beta_dict[beta]
  pt_bins_det = (binning_dict['pt_bins_det'])
  n_pt_bins_det = len(pt_bins_det) - 1
  det_pt_bin_array = array('d',pt_bins_det)
  
  for bin in range(1, n_pt_bins_det):
    min_pt_det = pt_bins_det[bin]
    max_pt_det = pt_bins_det[bin+1]
    
    plot_refolded_slice(hFoldedTruth, hData_PerBin, i, jetR, beta, regularizationParamName, min_pt_det, max_pt_det, output_dir_refolding, file_format)

#################################################################################################
# Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
#################################################################################################
def plot_refolded_slice(hFoldedTruth, hData_PerBin, i, jetR, beta, regularizationParamName, min_pt_det, max_pt_det, output_dir, file_format):
  
  hFoldedTruth.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
  hFolded_rg = hFoldedTruth.ProjectionY()
  hFolded_rg.SetName('hFolded_rg_R{}_B{}_{}_{}-{}'.format(jetR, beta, i, min_pt_det, max_pt_det))
  
  hData_PerBin.GetXaxis().SetRangeUser(min_pt_det, max_pt_det)
  hData_rg = hData_PerBin.ProjectionY()
  hData_rg.SetName('hData_rg_R{}_B{}_{}_{}-{}'.format(jetR, beta, i, min_pt_det, max_pt_det))
  
  yAxisTitle = '#frac{1}{N_{jet,inc}} #frac{dN}{d#theta_{g}}'
  legendTitle = ''
  h1LegendLabel = 'Folded truth, {} = {}'.format(regularizationParamName,i)
  h2LegendLabel = 'Measured pp'
  ratioYAxisTitle = 'Folded truth / Measured'
  outputFilename = os.path.join(output_dir, 'hFoldedTruth_R{}_B{}_{}-{}_{}{}'.format(analysis_utils.remove_periods(jetR), beta, min_pt_det, max_pt_det, i, file_format))
  plot_rg_ratio(hFolded_rg, hData_rg, None, yAxisTitle, ratioYAxisTitle, min_pt_det, max_pt_det, jetR, beta, outputFilename, 'width', legendTitle, h1LegendLabel, h2LegendLabel)

#################################################################################################
# Closure test
#################################################################################################
def plot_closure_test(response, hMC_Det, hMC_Truth, i, regularizationParamName, jetR, beta, beta_dict, output_dir, file_format):
  
  output_dir_closure = os.path.join(output_dir, 'ClosureTest')
  if not os.path.isdir(output_dir_closure):
    os.makedirs(output_dir_closure)

  # Unfold smeared det-level spectrum with RM
  unfold2 = ROOT.RooUnfoldBayes(response, hMC_Det, i)
  hUnfolded2 = unfold2.Hreco() # Produces the truth distribution, with errors, PerBin d

  binning_dict = beta_dict[beta]
  pt_bins_truth = (binning_dict['pt_bins_truth'])
  n_pt_bins_truth = len(pt_bins_truth) - 1
  truth_pt_bin_array = array('d',pt_bins_truth)

  for bin in range(1, n_pt_bins_truth-1):
    min_pt_truth = pt_bins_truth[bin]
    max_pt_truth = pt_bins_truth[bin+1]
    
    plot_closure_slice(hUnfolded2, hMC_Truth, i, jetR, beta, regularizationParamName, min_pt_truth, max_pt_truth, output_dir_closure, file_format)

#################################################################################################
# Apply RM to unfolded result, and check that I obtain measured spectrum (simple technical check)
#################################################################################################
def plot_closure_slice(hUnfolded, hMC_Truth, i, jetR, beta, regularizationParamName, min_pt_truth, max_pt_truth, output_dir, file_format):
  
  hUnfolded.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
  hUnfolded_rg = hUnfolded.ProjectionY()
  hUnfolded_rg.SetName('hUnfolded_rg_R{}_B{}_{}_{}-{}'.format(jetR, beta, i, min_pt_truth, max_pt_truth))
  
  hMC_Truth.GetXaxis().SetRangeUser(min_pt_truth, max_pt_truth)
  hMCTruth_rg = hMC_Truth.ProjectionY()
  hMCTruth_rg.SetName('hMCTruth_rg_R{}_B{}_{}_{}-{}'.format(jetR, beta, i, min_pt_truth, max_pt_truth))
  
  yAxisTitle = '#frac{1}{N_{jet,inc}} #frac{dN}{d#theta_{g}}'
  legendTitle = ''
  h1LegendLabel = 'Unfolded MC-det, {} = {}'.format(regularizationParamName,i)
  h2LegendLabel = 'MC-truth'
  ratioYAxisTitle = 'Unfolded MC det / Truth'
  outputFilename = os.path.join(output_dir, 'hClosure_R{}_B{}_{}-{}_{}{}'.format(analysis_utils.remove_periods(jetR), beta, min_pt_truth, max_pt_truth, i, file_format))
  plot_rg_ratio(hUnfolded_rg, hMCTruth_rg, None, yAxisTitle, ratioYAxisTitle, min_pt_truth, max_pt_truth, jetR, beta, outputFilename, 'width', legendTitle, h1LegendLabel, h2LegendLabel)

#################################################################################################
# Get errors from measured spectrum, stored as dictionary {bin:error}
#################################################################################################
def getMeasuredErrors(h):
  
  dictErrors = {}
  for binx in range(0, h.GetNbinsX()+1):
    for biny in range(0, h.GetNbinsY()+1):
      global_bin = h.GetBin(binx, biny)
      content = h.GetBinContent(global_bin)
      error = h.GetBinError(global_bin)
      if content!=0:
        dictErrors[global_bin] = error/content
      #print('Bin ({},{}) with content {} has error {}'.format(binx, biny, content, error/content))
      else:
        print('    Error: check your historgram ranges - something might be wrong')
        print('    Bin {}, content {}, error {}'.format(global_bin, content, error))
        dictErrors[global_bin] = 0
  return dictErrors

#################################################################################################
# Smear spectrum according to the error bars on the measured spectrum
#################################################################################################
def smearSpectrum(h, measuredErrors):
  
  # Zero errors in all bins
  for binx in range(0, h.GetNbinsX() + 1):
    for biny in range(0, h.GetNbinsY() + 1):
      global_bin = h.GetBin(binx, biny)
      h.SetBinError(global_bin, 0)
  
  # Loop through relevant bins and smear content according to errors in measured data, and set new errors
  r = ROOT.TRandom3(0)
  for bin, error in measuredErrors.items():
    content = h.GetBinContent(bin)
    errorNew = error * content
    contentNew = content + r.Gaus(0, errorNew)
    h.SetBinContent(bin, contentNew)
    h.SetBinError(bin, errorNew)
    #print('Setting bin {} to have relative error {}'.format(bin, errorNew/contentNew))

#################################################################################################
# Get MC-det 2D projection
#################################################################################################
def get_MCdet2D(fResponse, jetR, beta):
  
  # (pt-det, pt-true, theta_g-det, theta_g-true)
  name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
  hResponse = fResponse.Get(name)
  
  hResponse4D = hResponse.Clone()
  hResponse4D.SetName('{}_clone'.format(hResponse4D.GetName()))
  
  hMC_Det = hResponse4D.Projection(2,0)
  hMC_Det.SetName('hMC_Det_R{}_B{}'.format(jetR, beta))
  return hMC_Det

#################################################################################################
# Get MC-det 2D projection
#################################################################################################
def get_MCtruth2D(fResponse, jetR, beta):
  
  # (pt-det, pt-true, theta_g-det, theta_g-true)
  name = 'hResponse_JetPt_ThetaG_R{}_B{}'.format(jetR, beta)
  hResponse = fResponse.Get(name)
  
  hResponse4D = hResponse.Clone()
  hResponse4D.SetName('{}_clone'.format(hResponse4D.GetName()))
  
  hMC_Truth = hResponse4D.Projection(3,1)
  hMC_Truth.SetName('hMC_Truth_R{}_B{}'.format(jetR, beta))
  return hMC_Truth

#################################################################################################
# Plot spectra and ratio of h (and h3, if supplied) to h2
#################################################################################################
def plot_rg_ratio(h, h2, h3, yAxisTitle, ratioYAxisTitle, min_pt_det, max_pt_det, jetR, beta, outputFilename, scalingOptions = "", legendTitle = "",hLegendLabel = "", h2LegendLabel = "", h3LegendLabel = "", yRatioMax = 2.2):
  
  c = ROOT.TCanvas("c","c: pT",800,850)
  c.cd()
  pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
  pad1.SetBottomMargin(0)
  pad1.SetLeftMargin(0.15)
  pad1.SetRightMargin(0.05)
  pad1.SetTopMargin(0.05)
  pad1.Draw()
  pad1.cd()
  
  h.SetLineColor(1)
  h.SetLineWidth(2)
  h.SetLineStyle(1)
  
  h.Scale(1./h.Integral(), scalingOptions)
  h.GetYaxis().SetTitle(yAxisTitle)
  h.GetYaxis().SetTitleSize(0.06)
  h.GetYaxis().SetRangeUser(0., 3.)
  h.GetYaxis().SetLabelFont(43)
  h.GetYaxis().SetLabelSize(20)
  xAxisTitle = h.GetXaxis().GetTitle()
  h.GetXaxis().SetTitle("")
  
  h2.SetLineColor(4)
  h2.SetLineWidth(2)
  h2.SetLineStyle(1)
  h2.Scale(1./h2.Integral(), scalingOptions)
  
  h.Draw("hist same E")
  h2.Draw("hist same E")
  
  if h3:
    h3.SetLineColor(2)
    h3.SetLineWidth(2)
    h3.SetLineStyle(1)
    h3.Scale(1./h3.Integral(), scalingOptions)
    h3.Draw("hist same")
  
  c.cd()
  pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
  pad2.SetTopMargin(0)
  pad2.SetBottomMargin(0.35)
  pad2.SetLeftMargin(0.15)
  pad2.SetRightMargin(0.05)
  pad2.Draw()
  pad2.cd()

  # plot ratio h/h2
  hRatio = h.Clone()
  hRatio.Divide(h2)
  hRatio.SetMarkerStyle(21)
  hRatio.SetMarkerColor(1)
  
  hRatio.GetXaxis().SetTitleSize(30)
  hRatio.GetXaxis().SetTitleFont(43)
  hRatio.GetXaxis().SetTitleOffset(4.)
  hRatio.GetXaxis().SetLabelFont(43)
  hRatio.GetXaxis().SetLabelSize(20)
  hRatio.GetXaxis().SetTitle(xAxisTitle)
  
  hRatio.GetYaxis().SetTitle(ratioYAxisTitle)
  hRatio.GetYaxis().SetTitleSize(20)
  hRatio.GetYaxis().SetTitleFont(43)
  hRatio.GetYaxis().SetTitleOffset(2.2)
  hRatio.GetYaxis().SetLabelFont(43)
  hRatio.GetYaxis().SetLabelSize(20)
  hRatio.GetYaxis().SetNdivisions(505)
  min= hRatio.GetBinContent(hRatio.GetMinimumBin())
  max= hRatio.GetBinContent(hRatio.GetMaximumBin())
  #automatic zoom-in for a very small scatter of the points
  if min>0.5 and max<1.5:
    hRatio.GetYaxis().SetRangeUser(0.5,1.5)
  elif yRatioMax>2:
    hRatio.GetYaxis().SetRangeUser(0,yRatioMax)
  else:
    hRatio.GetYaxis().SetRangeUser(2-yRatioMax,yRatioMax)

  hRatio.Draw("P E")

  # plot ratio h3/h2
  if h3:
    hRatio3 = h3.Clone()
    hRatio3.Divide(h2)
    hRatio3.SetMarkerStyle(21)
    hRatio3.SetMarkerColor(2)
    hRatio3.Draw("P E same")
  
  pad1.cd()
  
  leg2 = ROOT.TLegend(0.55,0.7,0.88,0.93,legendTitle)
  leg2.SetFillColor(10)
  leg2.SetBorderSize(0)
  leg2.SetFillStyle(0)
  leg2.SetTextSize(0.04)
  leg2.AddEntry(h, hLegendLabel, "l")
  if h3:
    leg2.AddEntry(h3, h3LegendLabel, "l")
  if h2:
    leg2.AddEntry(h2, h2LegendLabel, "l")
  leg2.Draw("same")

  text_latex = ROOT.TLatex()
  text_latex.SetNDC()
  text = str(min_pt_det) + ' < #it{p}_{T,ch jet} < ' + str(max_pt_det)
  text_latex.DrawLatex(0.25, 0.85, text)

  text_latex = ROOT.TLatex()
  text_latex.SetNDC()
  text = 'R = ' + str(jetR) + '   #beta = ' + str(beta)
  text_latex.DrawLatex(0.25, 0.75, text)

  c.SaveAs(outputFilename)
  c.Close()

###################################################################################################
# Plot basic histogram
###################################################################################################
def plotHist(h, outputFilename, drawOptions = "", setLogy = False, setLogz = False, text = ""):
  
  c = ROOT.TCanvas("c","c: hist",600,450)
  c.cd()
  c.cd().SetLeftMargin(0.15)
  
  if setLogy:
    c.SetLogy()
  if setLogz:
    c.SetLogz()
  h.DrawCopy(drawOptions)

  if text:
    textFit = ROOT.TLatex()
    textFit.SetTextSize(0.04)
    textFit.SetNDC()
    textFit.DrawLatex(0.6,0.8,text)

  c.SaveAs(outputFilename)
  c.Close()

###################################################################################################
# Set legend parameters
###################################################################################################
def setupLegend(leg, textSize):
  
  leg.SetTextFont(42);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetMargin(0.25);
  leg.SetTextSize(textSize);
  leg.SetEntrySeparation(0.5);

###################################################################################################
# Configure style options
###################################################################################################
def setOptions():
  
  font = 42
  
  ROOT.gStyle.SetFrameBorderMode(0)
  ROOT.gStyle.SetFrameFillColor(0)
  ROOT.gStyle.SetCanvasBorderMode(0)
  ROOT.gStyle.SetPadBorderMode(0)
  ROOT.gStyle.SetPadColor(10)
  ROOT.gStyle.SetCanvasColor(10)
  ROOT.gStyle.SetTitleFillColor(10)
  ROOT.gStyle.SetTitleBorderSize(1)
  ROOT.gStyle.SetStatColor(10)
  ROOT.gStyle.SetStatBorderSize(1)
  ROOT.gStyle.SetLegendBorderSize(1)
  
  ROOT.gStyle.SetDrawBorder(0)
  ROOT.gStyle.SetTextFont(font)
  ROOT.gStyle.SetStatFont(font)
  ROOT.gStyle.SetStatFontSize(0.05)
  ROOT.gStyle.SetStatX(0.97)
  ROOT.gStyle.SetStatY(0.98)
  ROOT.gStyle.SetStatH(0.03)
  ROOT.gStyle.SetStatW(0.3)
  ROOT.gStyle.SetTickLength(0.02,"y")
  ROOT.gStyle.SetEndErrorSize(3)
  ROOT.gStyle.SetLabelSize(0.05,"xyz")
  ROOT.gStyle.SetLabelFont(font,"xyz")
  ROOT.gStyle.SetLabelOffset(0.01,"xyz")
  ROOT.gStyle.SetTitleFont(font,"xyz")
  ROOT.gStyle.SetTitleOffset(1.2,"xyz")
  ROOT.gStyle.SetTitleSize(0.045,"xyz")
  ROOT.gStyle.SetMarkerSize(1)
  ROOT.gStyle.SetPalette(1)
  
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptFit(0)

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  parser = argparse.ArgumentParser(description='Unfold theta_g distribution')
  parser.add_argument('-d', '--inputFileData', action='store',
                      type=str, metavar='inputFileData',
                      default='AnalysisResults.root',
                      help='Path of AnalysisResults.root file containing spectrum to be unfolded')
  parser.add_argument('-r', '--inputFileResponse', action='store',
                      type=str, metavar='inputFileResponse',
                      default='AnalysisResults.root',
                      help='Path of AnalysisResults.root file containing response matrix')
  parser.add_argument('-c', '--configFile', action='store',
                      type=str, metavar='configFile',
                      default='analysis_config.yaml',
                      help='Path of config file for analysis')
  parser.add_argument('-o', '--outputDir', action='store',
                      type=str, metavar='outputDir',
                      default='./unfolding_output/',
                      help='Output directory for plots to be written to')
  parser.add_argument('-i', '--imageFormat', action='store',
                      type=str, metavar='imageFormat',
                      default='.pdf',
                      help='Image format to save plots in, e.g. \".pdf\" or \".png\"')
                      
  # Parse the arguments
  args = parser.parse_args()
  
  print('Configuring...')
  print('inputFileData: \'{0}\''.format(args.inputFileData))
  print('inputFileResponse: \'{0}\''.format(args.inputFileResponse))
  print('configFile: \'{0}\''.format(args.configFile))
  print('ouputDir: \'{0}\''.format(args.outputDir))
  print('imageFormat: \'{0}\''.format(args.imageFormat))
  
  # If invalid inputFile is given, exit
  if not os.path.exists(args.inputFileData):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileData))
    sys.exit(0)
  if not os.path.exists(args.inputFileResponse):
    print('File \"{0}\" does not exist! Exiting!'.format(args.inputFileResponse))
    sys.exit(0)
  # If invalid configFile is given, exit
  if not os.path.exists(args.configFile):
    print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
    sys.exit(0)

  roounfold_rg(input_file_data = args.inputFileData, input_file_response = args.inputFileResponse, config_file = args.configFile, output_dir = args.outputDir, file_format = args.imageFormat)