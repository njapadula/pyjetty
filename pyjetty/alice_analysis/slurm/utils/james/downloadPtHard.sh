#!/bin/bash
#
# This script downloads Pt-hard bin data, and merges and scales appropriately using Pt-hard weights computed from histograms filled
# in AliAnalysisTaskPWGJEQA.
#
# The script has two modes, depending whether run-by-run QA is desired -- set flag doRunByRun.
#
# The script:
#   (1) Downloads output analysis files per Pt-hard bin and run (at the moment this assumes run-by-run output from train)
#   (2) Merges output files for each Pt-hard bin (per run, if run-by-run mode)
#   (3) Calls python script scalePtHardHistos.py to re-weight histograms in each Pt-hard bin file (per run, if run-by-run mode)
#   (4) Sums all weighted Pt-hard bins into a single final output file (per run, if run-by-run mode)
# To use:
#   Fill in appropriate fields in config. Then, execute script from anywhere,
#   and it will write out to folder TrainOutput in the specified outputDir.
#
# Author: James Mulligan <james.mulligan@yale.edu>

# LHC15o
# RUNLIST="246810 246751 246994 245954 246991 246012 245952 245833 245949 246766 246765 246989 246003 246053 245705 246089 245831 246763 246276 245702 246185 246495 246052 246676 246225 246809 246153 246275 245700 246760 246493 246182 246049 246001 246675 245829 246759 246808 246222 245692 246152 246984 246758 246048 246807 246181 246087 246272 246757 246434 246151 246042 246488 246804 246948 246982 246846 246928 246115 246945 246180 246805 246217 246851 245683 246847 246750 246271 246487 246431 246844 246037 246178 246980 246845 246428 246424 245923 246113 246036 245064 246392 244982 244918 244975 244980 244983 245061 245066 245068 246390 246391 245152 245151 245146 245145 245232 245231 246648 246583 246575 246553 245793 245775 245738 245148 245411 245454 245259 245409 245507 245554 245452 245505 245407 245353 245450 245504 245545 245544 245446 245349 245501 245401 245543 245441 245347 245397 245542 245497 245346 245540 245439 245496 245396 245345 245535 245343 245785 245766 245759 245752 245731 245729 246871 246870 246867 246865 245410"

# LHC16j5 (EMCal good)
# RUNLIST="246945 246846 246845 246844 246810 246809 246808 246807 246805 246804 246766 246765 246760 246759 246758 246757 246751 246750 246495 246493 246488 246487 246434 246424 246272 246271 246225 246222 246217 246115 246113 246089 246087 246053 246052 246037 246003 246001 245954 245952 245949 245833 245831 245829 245705 245702 245700 245683"

# LHC18b8
#RUNLIST="282343 282342 282341 282340 282314 282313 282312 282309 282307 282306 282305 282304 282303 282302 282247 282230 282229 282227 282224 282206 282189 282147 282146 282127 282126 282125 282123 282122 282120 282119 282118 282099 282098 282078 282051 282050 282031 282025 282021 282016 282008 282367 282366 282365"

# LHC19f4
RUNLIST="297588 297540 297481 297413 297333 297221 297117 296850 296694 296690 296512 296414 296379 296377 296269 296240 296060 295816 295712 295586"
SUBDIRS="0001 0002 0003 0004 0005 0006 0007 0008 0009 0010 0011 0012 0013 0014 0015 0016 0017 0018 0019 0020"

outputDir="/mnt/rstorage/alice/data/LHC19f4/403/child_2"
PtHardBins=20

data="sim"
year="2019"
period="LHC19f4_2"
trainName="HF_TreeCreator"
trainPWG="PWGHF"
trainNumber="403_20200404-2156_child_2"
filename="AnalysisResults.root"

PREFIX="/alice/${data}/${year}/${period}/"
SUFFIX="/${trainPWG}/${trainName}/${trainNumber}"

# Create outputDir and cd into it
if [ ! -d $outputDir/TrainOutput ]; then
  mkdir -p $outputDir/TrainOutput
fi
cd $outputDir/TrainOutput
echo "output dir: $outputDir/TrainOutput"

# Copy the files
for bin in $(seq 19 $PtHardBins);
do
  for RUN in $RUNLIST
  do
    for SUBDIR in $SUBDIRS
    do
      if [ ! -d "$bin/$RUN/$SUBDIR" ]; then
        mkdir -p "$bin/$RUN/$SUBDIR"
      fi        
      alien_cp alien://$PREFIX$bin/$RUN$SUFFIX/$SUBDIR/AnalysisResults.root $bin/$RUN/$SUBDIR
    done
  done
done
