#writen by Renqi Pan in Oct 22,2021
#to adjust total entries of dataset
echo "adjust total entries of dataset."
for dir in $( ls | grep Chunk)
do   
	cd $dir
	dataset=$(ls | grep .txt)
    echo $dir
    echo $dataset
    count=$(cat $dataset | wc -l)
    echo "number of lines: $count"
	if [[ $dataset =~ "TTToSemiLeptonic" ]]
	then
		entries=50267258
		flag=true
	elif [[ $dataset =~ "TTTo2L2Nu" ]]
	then
		entries=12199607
		flag=true
	elif [[ $dataset =~ "TTToHadronic" ]]
	then
		entries=51780528
		flag=true
	elif [[ $dataset =~ "WJetsToLNu_HT-100To200" ]]
    then
		entries=184360900
		flag=true
	elif [[ $dataset =~ "WJetsToLNu_HT-200To400" ]]
	then
		entries=49278900
		flag=true
	else
		flag=false
	fi
	if $flag
	then
		echo "number of entries is adjusted to $entries."
		echo "delete extra root files"
		sum=0
		num_line=0
		for line in $(cat $dataset)
		do
			if [ $sum -lt $entries ]
			then 
				run='root -l -q -b ../read_entry.cpp"(\"$line\")"'
				txt=$(eval $run 2> .out | grep has)
                num=${txt#* }
				num=${num% *}
				let sum=sum+num
				let num_line=num_line+1
              #  echo "sum: $sum, num_line: $num_line, num: $num"
			else
				break
			fi
		done
		head -n $num_line $dataset > tmp.txt && mv tmp.txt $dataset
		echo $dir
		echo "line $(($num_line+1)) to the last line are deleted"
	fi
    echo ""
	cd ../	

done
