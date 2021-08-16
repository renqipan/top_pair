#authenticate your grid certificate
voms-proxy-init --rfc --voms cms -valid 192:00
dasgoclient -query="dataset=/TTToSemiLeptonic_TuneCP*_13TeV-powheg-pythia8/*/NANOAODSIM"

# copy a specific root file to local storage
xrdcp root://xrootd-cms.infn.it//store/mc/RunIISummer19UL16NanoAOD/TTToSemiLeptonic_TuneCP5CR1_13TeV-powheg-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v13-v1/240000/5F888D3E-628B-C54C-A2D9-BFCA76580FB4.root .

# query dataset name of a root file
dasgoclient --query "dataset file=/store/mc/RunIISummer16MiniAODv3/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/120000/9AD694A9-E9EA-E811-B81E-6CC2173DAD00.root"
# query summary of a dataset
dasgoclient --query "summary dataset=/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"
# query availible root files of a dataset
dasgoclient --query "file dataset=/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/MINIAODSIM"

# copy a specific root file to local storage
xrdcp root://xrootd-cms.infn.it//store/mc/RunIISummer19UL18MiniAODv2/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/270000/005B1938-1AA2-8D48-A2B3-59800B0812CA.root ./
dasgoclient -query="dataset=/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/swertz-TopNanoAODv6p1_2018-0d1d4920f08f56d048ece029b873a2cc/USER"

#ttbar 2018 files
dasgoclient --query "file dataset=/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19_ext3-v1/NANOAODSIM"
#address
root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv5/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano1June2019_102X_upgrade2018_realistic_v19_ext3-v1/20000/B25E3DCB-2DD8-FD4E-9FBD-8AA9015657CD.root
#copy root file
xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv5/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano1June2019_102X_upgrade2018_realistic_v19_ext3-v1/20000/B25E3DCB-2DD8-FD4E-9FBD-8AA9015657CD.root ~/work/top_pair/
#cern box address: /eos/user/r/repan