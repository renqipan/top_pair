#!/bin/bash
#voms-proxy-init --voms cms -valid 192:00 -out ~/temp/x509up
# sumbit asssinment: condor_submit condor.sub
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh
cd /afs/cern.ch/user/r/repan/work/top_pair/condor_dilepton/condor_out/$1
file=$(ls *.txt)
root -l -q -b ../../get_dilepton_condor.cpp"(\"$file\")"

if rootfile=$(ls *.root) 2> /dev/null
then
	if [ $rootfile =~ "TTTo" ]
	then
		echo $rootfile
		root -l -q -b ../../add_weight_branch.c"(\"$rootfile\")"
	fi
else
	echo "the root file doesn't exist. task unfinished"
fi


#var=$1
#if [[ $var =~ "_pythia" ]]
#then
#	temp=${var%%_pythia8*}
#else
#    temp=${var%%-pythia8*}
#fi
#process=${temp:1}
#dasgoclient --query "file dataset=$1" > ${process}.txt
#dasgoclient -query="file dataset=$1" > ${process}.txt
# root -l -q -b ../get_info.cpp"(\"$1\")"