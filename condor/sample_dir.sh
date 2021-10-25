#writen by Renqi Pan in August 30,2021
input="dataset.txt"
i=0
for dataset in `cat $input`
do 
	echo $dataset
	var="Chunk$i"
	rm -rf $var
    mkdir $var
    cd $var
	mkdir condor_out
	if [[ $dataset =~ "_pythia" ]]
	then
		temp=${dataset%%_pythia8*}
	else
    	temp=${dataset%%-pythia8*}
    fi
    process=${temp:1}
    #dasgoclient --query "file dataset=$dataset" > ${process}.txt
    dasgoclient -query="file dataset=$dataset" > ${process}.txt
	cd ../
	let i=i+1
done

for dataset in `cat $input`
do 
	echo $dataset
	var="Chunk$i"
    cd $var
    if [[ $dataset =~ "_pythia" ]]
	then
		temp=${dataset%%_pythia8*}
	else
    	temp=${dataset%%-pythia8*}
    fi
    process=${temp:1}
	total=$(cat *.txt | wc -l)
	lines=$(($total/7))
	remind=$(($total%7))
	num_txt=1
	num_line=1
	touch ${process}_{1..8}.txt
	for line in $(cat *.txt)
	do
		if [ $num_line -le $lines ]; then
			txt_name=${process}_$(num_txt).txt
			echo line > $txt_name
			let num_line=num_line+1

		fi
		if [ $num_line -eq $(($lines+1)) ]; then
			let num_txt=num_txt+1
			let num_line=1
			echo $txt_name
		fi

	done

	cd ../
	let i=i+1
done