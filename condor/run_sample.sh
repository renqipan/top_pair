input="dataset.txt"
i=0
for process in `cat $input`
do 
	echo $process
	var="Chunk$i"
	mkdir $var
        cd $var
	mkdir condor_out
	let i=i+1
	cd ../
done
