#writen by Renqi Pan in Oct 22,2021
# to write dataset files into .txt and divide
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

#source adjust_entry.sh

echo "divide each dataset into several pieces "
for var in `ls | grep Chunk`
do 
    cd $var
    process=$( ls *.txt )
    process=${process%%.txt*}
    if [[ $process =~ "TTToSemiLeptonic" ]]
    then 
	    total_ttsemi=$(cat *.txt | wc -l)
    fi
	total=$(cat *.txt | wc -l)
    divide=$total  #edit this line  
	lines=$(($total / $divide ))
	remind=$(($total % $divide ))
	echo "total=$total line=$lines divide=$divide remind=$remind"
    num_txt=1
	num_line=1
	#touch ${process}_{1..8}.txt
	if [[ $total -le $divide ]];then
		for line in $(cat *.txt)
		do
          mkdir ../${process}_${num_txt}
          cd ../${process}_${num_txt}
		  txt_name=${process}_${num_txt}.txt
		  echo $line >> $txt_name
		  #echo $txt_name
		  let num_txt=num_txt+1
          cd ../$var

		done

	fi
	cd ../
done
out=condor_out
rm -rf $out
mkdir $out
exp="mv *_{1.."$total_ttsemi"} "$out
eval $exp 2> /dev/null
rm -rf Chunk*
ls $out  > condor_list.txt
echo "directories are written into condor_list.txt"
