from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = '17_WJetsToLNu_HT-800To1200'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'

tempfilelist = ['Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFPuppi.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK4PFPuppi.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFPuppi.txt']
filepath ='/afs/cern.ch/work/q/qiguo/B2G/work/03/03_10/CMSSW_10_2_18/src/ExoDiBosonResonances/EDBRTreeMaker/test/'
filelist = []
for i in tempfilelist :
   filelist.append(filepath+i)

config.JobType.inputFiles = filelist
#config.JobType.inputFiles = ['Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFchs.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK8PFPuppi.txt','Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFPuppi.txt','Fall17_17Nov2017_V32_MC_L2Relative_AK4PFPuppi.txt','Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFPuppi.txt','L1PrefiringMaps_new.root']

config.JobType.psetName    = 'analysis.py'

config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =2
config.Data.totalUnits = -1
config.Data.publication = False

config.Data.outLFNDirBase='/store/user/qiguo/crab/VVV/17sample'

config.Data.outputDatasetTag = '17_WJetsToLNu_HT-800To1200'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERNBOX'
