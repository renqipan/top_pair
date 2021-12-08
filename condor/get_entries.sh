dir=output2
echo "get entries from samples"
cd $dir
for tt in $(ls *.root | grep "TTTo")
do
	echo $tt
	run='root -l -q -b ../get_entries.cpp"(\"$tt\")"'
	if txt=$(eval $run 2> .out | grep has)
	then	
	    num=${txt#* }
		num=${num% *}
		echo "has $num entries."
		echo ""
	else
		echo "$tt read failed"
	fi
done
cd ../
