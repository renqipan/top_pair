dir=output2
txt=condor_sec.txt
#rm -rf $dir
mkdir -p $dir
cd condor_out
for var in `ls`
do 
    cp $var/run.log ../$dir/${var}_run.log 2>/dev/null
    cp $var/run.err ../$dir/${var}_run.err 2>/dev/null
	if  cp $var/run.out ../$dir/${var}_run.out 2>/dev/null
	then 
		if  mv $var/*.root ../$dir 2>/dev/null
		then 
			flag=true	
	    else
		    echo "$var didn't succeed"
		    echo $var >> $txt
		fi	

	else
		echo "$var didn't finish"
		echo $var >> $txt
        rm $var/*.root 
	fi
    
#   rm -rf $var
done
mv $txt ..
cd ../$dir
cat *_run.out > all.out
cat *_run.err > all.err
cat *_run.log > all.log
rm *_run.*
echo "if any see any anomalous tips, please have a cheack,"
echo "The failed tasks(if have) are summaried in $txt "
cd ..
