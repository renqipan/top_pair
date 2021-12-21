#!/bin/bash
#voms-proxy-init --voms cms -valid 192:00 -out ~/temp/x509up
# sumbit asssinment: condor_submit condor.sub
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/thisroot.sh
mkdir -p myout
output=$PWD/myout
echo "output: $output"
cd /afs/cern.ch/user/r/repan/work/top_pair/condor_semi/condor_out/$1
file=$(ls *.txt)
dir=$(cat $file)
dir="root://cms-xrd-global.cern.ch/"$dir
echo "input file: $dir"
cd /afs/cern.ch/user/r/repan/work/top_pair/CMSSW_10_6_19_patch2/src/PhysicsTools/NanoAODTools/crab
eval `scramv1 runtime -sh`   
python crab_script.py $dir $output
if sys_root=$(ls $output| grep .root);then
	cd $output
	inputFile=${file%.txt*}
	inputFile=${inputFile}.root
	mv $sys_root $inputFile
	cd /afs/cern.ch/user/r/repan/work/top_pair/condor_semi/condor_out/$1
	root -l -q -b ../../get_info_3jet_condor.cpp"(\"$output\",\"$inputFile\")"
    eos="/afs/cern.ch/user/r/repan/work/top_pair/condor_semi/output"
	if rootfile=$(ls $output/*.root) 2> /dev/null
	then
		if [[ $rootfile =~ "TTTo" ]]
		then
			echo $rootfile
			root -l -q -b ../../add_weight_branch.c"(\"$rootfile\")"
            mv $rootfile $eos
        else 
            mv $rootfile $eos
		fi
	else
		echo "Events reconstruction failed. task unfinished"
	fi
else
	echo "Running crab_script.py failed!!!"
fi
rm -rf $output

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
