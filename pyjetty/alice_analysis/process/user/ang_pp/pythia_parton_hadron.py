#!/usr/bin/env python

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import copy
import argparse
import os
import numpy as np

from pyjetty.mputils import *

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

from pyjetty.alice_analysis.process.user.ang_pp.helpers import lambda_beta_kappa

def main():
    parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly',
                                     prog=os.path.basename(__file__))
    pyconf.add_standard_pythia_args(parser)
    # Could use --py-seed
    parser.add_argument('--user-seed', help='PYTHIA starting seed', default=1111, type=int)
    parser.add_argument('--output', default="output.root", type=str)
    parser.add_argument('--sd_beta', help='SoftDrop beta', default=None, type=float)
    parser.add_argument('--jetR', help='Jet radius/resolution parameter', default=0.4, type=float)
    parser.add_argument('--no-match-level', help="Save simulation for only one level with " + \
                        "no matching. Options: 'p', 'h', 'ch'", default=None, type=str)
    parser.add_argument('--p-ch-MPI', help="Match between parton level (no MPI) and charged " + \
                        "hadron level (with MPI) for ALICE response characterization",
                        action='store_true', default=False)
    args = parser.parse_args()

    level = args.no_match_level
    if args.p_ch_MPI:
        if level:
            print("ERROR: --no-match-level and --p-ch-MPI cannot be set simultaneously.")
            exit(1)
    if level not in [None, 'p', 'h', 'ch']:
        print("ERROR: Unrecognized type %s. Please use 'p', 'h', or 'ch'" % args.type_only)
        exit(1)

    # Angularity beta values
    betas = [1, 1.5, 2, 3]

    if args.user_seed < 0:
        args.user_seed = 1111
    pinfo('user seed for pythia', args.user_seed)
    # mycfg = ['PhaseSpace:pThatMin = 100']
    mycfg = ['Random:setSeed=on', 'Random:seed={}'.format(args.user_seed)]
    mycfg.append('HadronLevel:all=off')

    # Have at least 1 event
    if args.nev < 1:
        args.nev = 1

    if args.p_ch_MPI:
        d = vars(args)
        d['py_noMPI'] = True
        pythia_noMPI = pyconf.create_and_init_pythia_from_args(args, mycfg)
        args_MPI = copy.deepcopy(args)
        d = vars(args_MPI)
        d['py_noMPI'] = False
        pythia = pyconf.create_and_init_pythia_from_args(args_MPI, mycfg)
    else:
        pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

    # print the banner first
    fj.ClusterSequence.print_banner()
    print()
    # set up our jet definition and a jet selector
    jet_R0 = args.jetR
    jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
    print(jet_def)

    # hadron level - ALICE
    max_eta_hadron = 0.9
    pwarning('max eta for particles after hadronization set to', max_eta_hadron)
    parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
    jet_selector = fj.SelectorPtMin(5.0) #& fj.SelectorAbsEtaMax(max_eta_hadron - jet_R0)

    max_eta_parton = max_eta_hadron + 2. * jet_R0
    pwarning('max eta for partons set to', max_eta_parton)
    parts_selector_p = fj.SelectorAbsEtaMax(max_eta_parton)

    if not args.sd_beta is None:
        args.output = args.output.replace('.root', '_sdbeta{}.root'.format(args.sd_beta))
    outf = ROOT.TFile(args.output, 'recreate')
    outf.cd()
    t = ROOT.TTree('t', 't')
    tw = RTreeWriter(tree=t)

    count1 = 0
    count2 = 0

    # event loop
    for iev in range(args.nev):  #tqdm.tqdm(range(args.nev)):
        if not pythia.next():
            if args.p_ch_MPI:
                pythia_noMPI.next()
            continue
        if args.p_ch_MPI and not pythia_noMPI.next():
            continue

        parts_pythia_p = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
        parts_pythia_p_selected = parts_selector_p(parts_pythia_p)

        if args.p_ch_MPI:
            parts_pythia_p = pythiafjext.vectorize_select(pythia_noMPI, [pythiafjext.kFinal], 0, True)
            parts_pythia_p_selected = parts_selector_p(parts_pythia_p)

        hstatus = pythia.forceHadronLevel()
        if not hstatus:
            #pwarning('forceHadronLevel false event', iev)
            continue
        # parts_pythia_h = pythiafjext.vectorize_select(
        #     pythia, [pythiafjext.kHadron, pythiafjext.kCharged])
        parts_pythia_h = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], 0, True)
        parts_pythia_h_selected = parts_selector_h(parts_pythia_h)

        parts_pythia_hch = pythiafjext.vectorize_select(
            pythia, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
        parts_pythia_hch_selected = parts_selector_h(parts_pythia_hch)

        # pinfo('debug partons...')
        # for p in parts_pythia_p_selected:
        #   pyp = pythiafjext.getPythia8Particle(p)
        #   print(pyp.name())
        # pinfo('debug hadrons...')
        # for p in parts_pythia_h_selected:
        #   pyp = pythiafjext.getPythia8Particle(p)
        #   print(pyp.name())
        # pinfo('debug ch. hadrons...')
        # for p in parts_pythia_hch_selected:
        #   pyp = pythiafjext.getPythia8Particle(p)
        #   print(pyp.name())

        # parts = pythiafjext.vectorize(pythia, True, -1, 1, False)
        jets_p = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_p)))
        jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h)))
        jets_ch = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_hch)))

        if level:  # Only save info at one level w/o matching
            jets = locals()["jets_%s" % level]
            for jet in jets:
                tw.fill_branch('iev', iev)
                tw.fill_branch(level, jet)
                kappa = 1
                for beta in betas:
                    label = str(beta).replace('.', '')
                    tw.fill_branch('l_%s_%s' % (level, label),
                                   lambda_beta_kappa(jet, jet_R0, beta, kappa))
            continue

        if not args.sd_beta is None:
            sd = fjcontrib.SoftDrop(args.sd_beta, 0.1, jet_R0)

        for i,jchh in enumerate(jets_ch):

            # match hadron (full) jet
            drhh_list = []
            for j, jh in enumerate(jets_h):
                drhh = jchh.delta_R(jh)
                if drhh < jet_R0 / 2.:
                    drhh_list.append((j,jh))
            if len(drhh_list) != 1:
                count1 += 1
            else:  # Require unique match
                j, jh = drhh_list[0]

                # match parton level jet
                dr_list = []
                for k, jp in enumerate(jets_p):
                    dr = jh.delta_R(jp)
                    if dr < jet_R0 / 2.:
                        dr_list.append((k, jp))
                if len(dr_list) != 1:
                    count2 += 1
                else:
                    k, jp = dr_list[0]

                    # pwarning('event', iev)
                    # pinfo('matched jets: ch.h:', jchh.pt(), 'h:', jh.pt(),
                    #       'p:', jp.pt(), 'dr:', dr)

                    tw.fill_branch('iev', iev)
                    tw.fill_branch('ch', jchh)
                    tw.fill_branch('h', jh)
                    tw.fill_branch('p', jp)

                    kappa = 1
                    for beta in betas:
                        label = str(beta).replace('.', '')
                        tw.fill_branch("l_ch_%s" % label,
                                       lambda_beta_kappa(jchh, jet_R0, beta, kappa))
                        tw.fill_branch("l_h_%s" % label,
                                       lambda_beta_kappa(jh, jet_R0, beta, kappa))
                        tw.fill_branch("l_p_%s" % label,
                                       lambda_beta_kappa(jp, jet_R0, beta, kappa))

                    if not args.sd_beta is None:
                        jchh_sd = sd.result(jchh)
                        jchh_sd_info = fjcontrib.get_SD_jet_info(jchh_sd)
                        jh_sd = sd.result(jh)
                        jh_sd_info = fjcontrib.get_SD_jet_info(jh_sd)
                        jp_sd = sd.result(jp)
                        jp_sd_info = fjcontrib.get_SD_jet_info(jp_sd)

                        tw.fill_branch('p_zg', jp_sd_info.z)
                        tw.fill_branch('p_Rg', jp_sd_info.dR)
                        tw.fill_branch('p_thg', jp_sd_info.dR/jet_R0)
                        tw.fill_branch('p_mug', jp_sd_info.mu)

                        tw.fill_branch('h_zg', jh_sd_info.z)
                        tw.fill_branch('h_Rg', jh_sd_info.dR)
                        tw.fill_branch('h_thg', jh_sd_info.dR/jet_R0)
                        tw.fill_branch('h_mug', jh_sd_info.mu)

                        tw.fill_branch('ch_zg', jchh_sd_info.z)
                        tw.fill_branch('ch_Rg', jchh_sd_info.dR)
                        tw.fill_branch('ch_thg', jchh_sd_info.dR/jet_R0)
                        tw.fill_branch('ch_mug', jchh_sd_info.mu)

            #print("  |-> SD jet params z={0:10.3f} dR={1:10.3f} mu={2:10.3f}".format(
            #    sd_info.z, sd_info.dR, sd_info.mu))

    tw.fill_tree()
    if args.p_ch_MPI:
        pythia_noMPI.stat()
    pythia.stat()

    print("%i jets cut at first match criteria; %i jets cut at second match criteria." % 
          (count1, count2))

    outf.Write()
    outf.Close()

if __name__ == '__main__':
    main()
