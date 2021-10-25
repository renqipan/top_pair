if rootfile=$(ls semi_*.root) 2> /dev/null
then
	if [ $rootfile =~ "TTTo" ]
	then
		echo $rootfile
		root -l -q -b ./add_weight_branch.c"(\"$rootfile\")"
	fi
else
	echo "the root file doesn't exist. task unfinished"
fi
