dir=output
txt=condor_sec.txt
rm -rf $dir
mkdir $dir
cd condor_out
for var in `ls`
do 
    cp $var/run.log ../$dir/${var}_run.log 2>/dev/null
    cp $var/run.err ../$dir/${var}_run.err 2>/dev/null
	if  cp $var/run.out ../$dir/${var}_run.out 2>/dev/null
	then 
		if  mv $var/*.root ../$dir 2>/dev/null
		then 
			if cat $var/run.err | grep "error"
			then
				echo "$dir have errors !!!!"
				echo $var >> $txt
			fi
		else
			echo "$var didn't succed"
			echo $var >> $txt
				
		fi	

	else
		echo "$var didn't finish"
		echo $var >> $txt
	fi
    
#   rm -rf $var
done
cd ../$dir
cat *_run.out > all2.out
cat *_run.err > all2.err
cat *_run.log > all2.log
rm *_run.*
echo "if any see any anomalous tips, please have a cheack,"
echo "The failed tasks(if have) are summaried in $txt "
