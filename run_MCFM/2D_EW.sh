#!/bin/bash
# bash script to run MCFM-8.3
Cpq3=(0 1)
Cpu=(0 1)
ReCup=(1 2)
ImCup=(0 1)
m34min=(345 360 380 400 450 500 550 600 650 700 750 800 900 1200)
m34max=(360 380 400 450 500 550 600 650 700 750 800 900 1200 2000)
#m34min=(345 360 )
#m34max=(360 380 )
output=EW_2Dcorrection
i=0
for a1 in ${Cpq3[*]}
do 
for a2 in ${Cpu[*]}
do	
for a3 in ${ReCup[*]} 
do
for a4 in ${ImCup[*]}
do
for((k=0;k<14;k++))
do
	echo "Current parameters are: $a1 $a2 $a3 $a4."
	echo "mtt range: ${m34min[k]}, ${m34max[k]}"
	coupl=ci$a1$a2$a3$a4
	mtt=${m34min[k]}_${m34max[k]}
	file=${coupl}_${mtt}

	rm -rf $file
	mkdir -p $file
	cp input.DAT $file
	cd $file
	sed -i "50s/0.0/${m34min[k]}/" input.DAT
	sed -i "51s/14000.0/${m34max[k]}/" input.DAT 
	sed -i "128s/+0.0/$a1/" input.DAT
	sed -i "129s/+0.0/$a2/" input.DAT
	sed -i "130s/+0.0/$a3/" input.DAT
	sed -i "131s/+0.0/$a4/" input.DAT
	sed -i "18s/ci0000/$file/" input.DAT
	echo "$file is prepared."
	cd ../
	nohup ./mcfm_omp $file input.DAT >/dev/null  2>log &
	let i=i+1
	if [ $i -eq 20 ]; then
		sleep 60m
		let i=0
	fi
done
done
done
done
done

flag="true"
while [ $flag = "true" ]
do
	ps | grep "mcfm" 
	if [ $? -eq 0 ]; then
		echo "the program is runing."
		sleep 3m
	else 
		echo "the MCFM is finished."
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
	cp draw.cpp $output
	files=$(ls |grep ci)
	for filename in $files
	do
		cd $filename
		sed -i '172s/ci.*_ci/ci/' *.C
		# delete lines containing "NaN" and its bin
		vim -e -c 'g/SetBinError.*NaN/.-2,.d' -c 'wq' *.C 
		root -l -q -b *.C
		cp *.root ../$output/
		cd ../ 
	done
	echo "all root files are copied to $output."
	cd $output
	for couplings in `ls ci* | sed "s/_.*//"| sort -u`
	do
		root -l -q -b draw.cpp"(\"$couplings\")"
	done
 	rm -rf ci*.root
	cd ../
	echo "EW corrections files are storied in $output."
    	rm -rf ci*
fi


