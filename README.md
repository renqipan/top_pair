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
3. run get_info_3jet_condor.cpp through condor  
**source sample_dir.sh**  
**cd condor_sample/**  
**condor_submit condor.sub**  
4. check condor processing and merge results into a directory  
**condor_q**  
**source merge_out.sh**  
5. add EW weights to root files  
**cd /afs/cern.ch/user/r/repan/work/top_pair/condor**  
**root -l -q -b add_weight_branch.c**
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


