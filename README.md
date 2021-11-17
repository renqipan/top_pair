ttbar production at LHC in the SMEFT framewrok  
1. get data from nanoaod
2. reconstruct top qurak pairs through maximum likelihood
2. estimate ttbar seliletopnic reconstruction efficiency
3. run all signal and background dataset through lxplus condor
4. scan the likelihood to get the sensitivity  

Steps to run the code: 
1. visit the directory in lxplus:  
**cd /afs/cern.ch/user/r/repan/work/top_pair/condor**  
2. get voms-proxy  
**voms-proxy-init --voms cms -valid 192:00 -out ~/temp/x509up**    
3. run get_info_3jet_condor.cpp and add_weight_branch.c through condor  
**source divide_sample_dir.sh**  
**cd condor_sample/**  
**condor_submit condor.sub**  
4. check condor processing and merge results into a directory  
**condor_q**  
**source check.sh** # check processing  
**source merge_out.sh** #copy all root files to one directory 
5. add EW weights to root files(if didn't add)  
**cd /afs/cern.ch/user/r/repan/work/top_pair/condor**  
**source add_weight.sh**
6. wirte histogram distribution of siganl and background into  
a root file and preapre a datacard(.txt)  
**root -l -q -b prepare.cpp**  
7. write a python script to paramerize the model(signal),  
put it in HiggsAnalysis/CombinedLimit/python and then rebuild    
**cd ~/tth_cms/CMSSW_8_1_0/src**  
**cmsenv**  
**cd HiggsAnalysis/CombinedLimit**  
**scramv1 b clean; scramv1 b**  
8. get likelihood scan through higgs combie tool  
**source combine_3jets.sh**  
**source combine_semi.sh**  


