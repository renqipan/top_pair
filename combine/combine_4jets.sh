#!/bin/bash
source env.sh
cd $YourCombine_CMSSW_src_Dir #/directory of higgs combine analysis
eval `scramv1 runtime -sh`
cd $YourDatacard13TeV_Dir #/ directory of datacard
cd datacard2/

text2workspace.py ttbar_4jets.txt -o workspace_ttbar_4jets.root --PO doStage0 --PO doacttbar -P HiggsAnalysis.CombinedLimit.stagex_ttwc3:stagex_ttwc3 -m 125 --X-allow-no-background  -v 7 
text2workspace.py ttbar_4jets.txt -o workspace_ttbar_4jets_fcp.root --PO doStage0 --PO doacttbar --PO dofcp -P HiggsAnalysis.CombinedLimit.stagex_ttwc3:stagex_ttwc3 -m 125 --X-allow-no-background  -v 7 

# -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=40 --alignEdges=1 -P x --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_x -v 3 -m 125 --setParameterRanges x=-1.0,1.0 --saveInactivePOI=1
combine -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=60 --alignEdges=1 -P y --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_y -v 3 -m 125 --setParameterRanges y=-2.0,2.0 --saveInactivePOI=1
combine -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=60 --alignEdges=1 -P z --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_z -v 3 -m 125 --setParameterRanges z=-1.0,3.0 --saveInactivePOI=1
combine -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=60 --alignEdges=1 -P k --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_k -v 3 -m 125 --setParameterRanges k=-2.0,2.0 --saveInactivePOI=1

combine -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=3600 --alignEdges=1  -P k -P y --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_2Dyk -v 3 -m 125 --setParameterRanges y=-2.0,2.0:k=-2.0,2.0 --saveInactivePOI=1
combine -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=3600 --alignEdges=1  -P z -P k --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_2Dzk -v 3 -m 125 --setParameterRanges z=-1.0,3.0:k=-2.0,2.0 --saveInactivePOI=1
combine -M MultiDimFit workspace_ttbar_4jets.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=3600 --alignEdges=1  -P y -P z --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_2Dyz -v 3 -m 125 --setParameterRanges y=-2.0,2.0:z=-1.0,3.0 --saveInactivePOI=1
combine -M MultiDimFit workspace_ttbar_4jets_fcp.root -S 1 -t -1 --expectSignal=1 --algo=grid --points=40 --alignEdges=1 -P fcp --floatOtherPOIs=1 --X-rtd TMCSO_AdaptivePseudoAsimov=10 -n ttbar_4jets_fcp -v 3 -m 125 --setParameterRanges fcp=-1.0,1.0:y=0.0,0.0 --saveInactivePOI=1

rename higgsCombine limit_ *.root
rename MultiDimFit.mH125. '' *.root
cd ../
#  -S [ --systematics ] arg (=1) Include constrained systematic uncertainties, 
#-S 0 will ignore systematics constraint terms in the datacard.
#   -t [ --toys ] arg (=0)  Number of Toy MC extractions
# --algo arg (=none)   Algorithm to compute uncertainties
#--grid arg  Use the specified file containing a grid of SamplingDistributions for the limit
# --X-rtd arg Define some constants to be used at runtime (for debugging purposes). 
# -v [ --verbose ] arg (=0) Verbosity level (-1 = very quiet; 0 =quiet, 1 = verbose, 2+ = debug) 
# --alignEdges arg (=0) Align the grid points such that the endpoints of the ranges are included
