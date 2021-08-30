#written by Renqi Pan in August 30,2021
rm -rf out_put
mkdir  out_put
input="dataset.txt"
i=0
for process in `cat $input`
do 
    var=Chunk$i
    cp $var/condor_out/run.out out_put/${var}_run.out 2>/dev/null
    cp $var/condor_out/run.err out_put/${var}_run.err 2>/dev/null
    cp $var/condor_out/run.log out_put/${var}_run.log 2>/dev/null
#    mv $var/new*.root out_put 2>/dev/null || :
#    rm -rf $var
    let i=i+1
done
cd out_put
cat *_run.out > all.out
cat *_run.err > all.err
cat *_run.log > all.log
rm *run*
