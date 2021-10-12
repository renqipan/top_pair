#!/bin/bash
# bash script to run MCFM-8.3
output=EW_kappa
Cpq3=(-1.0 -0.5 0.5 1.0)
Cpu=(-1.0 -0.5 0.5 1.0)
#ReCup=(1.0 1.5 2.0 3.0)
#ImCup=(0.0 0.5 1.0 2.0)
ReCup=(1.0 2.0 3.0 4.0)
ImCup=(0.0 1.0 2.0)
for a3 in ${ReCup[*]} 
do
for a4 in ${ImCup[*]}
do
	echo "Current parameters are: $a3 $a4."
	dir=ci00$a3$a4
	rm -rf $dir
	mkdir -p $dir
	cp input.DAT $dir
	cd $dir
	sed -i "130s/+0.0/$a3/" input.DAT
	sed -i "131s/+0.0/$a4/" input.DAT
	sed -i "18s/ci0000/$dir/" input.DAT
	cd ../
	nohup ./mcfm_omp $dir/ input.DAT >/dev/null  2>log &
done
done

flag="true"
while [ $flag = "true" ]
do
	ps | grep "mcfm" 
	if [ $? -eq 0 ]; then
		echo "the program is runing."
		sleep 1m
	else 
		echo "the program is finished."
		flag="false"
	fi
done

ps | grep "mcfm" 
if [ $? -ne 0 ]; then
	echo "the program is finished."
	echo "run the .C files via CERN ROOT." 
	echo "move the .root files to a new directory."
	rm -rf $output
	mkdir -p $output
	files=$(ls |grep ci00)
	for filename in $files
	do
		cd $filename
		sed -i '172s/ci.*_ci/ci/' *.C
		vim -e -c 'g/SetBinError.*NaN/.-2,.d' -c 'wq' *.C 
		root -l -q -b *.C
		mv *.root ../$output/
		cd ../ 
	done
	echo "all files are copied to $output."
#	cd EW_files/
#	rename 's/.*ci/ci/' *.root
#	cd ../
	rm -rf ci00*
	echo "EW corrections files are storied in $output."
fi
