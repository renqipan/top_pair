dir=output
echo "add weights to tt samples"
cd $dir
for tt in $(ls *.root | grep "TTTo")
do
	root -l -q -b ../add_weight_branch.c"(\"$tt\")"
done
cd ../
