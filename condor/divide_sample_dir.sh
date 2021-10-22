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

source adjust_entry.sh

echo "divide each dataset into several pieces"
i=0
divide=8
for dataset in `cat $input`
do 
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
	echo "total=$total lines=$lines divide=$divide remind=$remind"
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
		  echo $txt_name
		  let num_txt=num_txt+1
          cd ../$var

		done

	fi
	if [[ $total -gt $divide ]]; then
		for line in $(cat *.txt)
		do  
			if [[ ( $num_line -le $lines && $num_txt -le $(($divide-1)) ) ||$num_txt -eq $divide ]]; then
				txt_name=${process}_${num_txt}.txt
                mkdir -p ../${process}_${num_txt}
                cd ../${process}_${num_txt}
				echo $line >> $txt_name
				let num_line=num_line+1

			fi
			if [[ $num_line -eq $(($lines+1)) || ( $num_txt -eq $divide && $num_line -eq $(($remind+1)) ) ]]; then
				let num_txt=num_txt+1
				let num_line=1
				echo $txt_name
			fi
            cd ../$var

		done
	fi

	cd ../
	let i=i+1
done
out=condor_out
rm -rf $out
mkdir $out
exp="mv *_{1.."$divide"} "$out
eval $exp
rm -rf Chunk*
ls $out  > condor_list.txt
echo "directories are written into condor_list.txt"
