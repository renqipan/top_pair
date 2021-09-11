#!/bin/bash
#voms-proxy-init --voms cms -valid 192:00 -out ~/temp/x509up
# sumbit asssinment: condor_submit condor.sub
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh
cd /afs/cern.ch/user/r/repan/work/top_pair/condor_dilepton/Chunk$2
var=$1
temp=${var%%-pythia8*}
process=${temp:1}
#dasgoclient --query "file dataset=$1" > ${process}.txt
#dasgoclient -query="file dataset=$1" > ${process}.txt
# root -l -q -b ../get_info.cpp"(\"$1\")"
root -l -q -b ../get_dilepton_condor.cpp"(\"${process}.txt\")"

