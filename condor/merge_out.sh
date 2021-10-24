dir=output
rm -rf $dir
mkdir $dir
input="dataset.txt"
i=0
for process in `cat $input`
do 
    var=Chunk$i
    cp $var/condor_out/run.out $dir/${var}_run.out 2>/dev/null
    cp $var/condor_out/run.err $dir/${var}_run.err 2>/dev/null
    cp $var/condor_out/run.log $dir/${var}_run.log 2>/dev/null
    mv $var/new*.root $dir 2>/dev/null || :
#    rm -rf $var
    let i=i+1
done
cd $dir
cat *_run.out > all2.out
cat *_run.err > all2.err
cat *_run.log > all2.log
rm *run*

