dir=output
rm -rf $dir
mkdir $dir
cd condor_out
for var in `ls`
do 
    cp $var/run.out ../$dir/${var}_run.out 2>/dev/null
    cp $var/run.err ../$dir/${var}_run.err 2>/dev/null
    cp $var/run.log ../$dir/${var}_run.log 2>/dev/null
    mv $var/new*.root ../$dir 2>/dev/null || :
#    rm -rf $var
done
cd $dir
cat *_run.out > all2.out
cat *_run.err > all2.err
cat *_run.log > all2.log
rm *run*

