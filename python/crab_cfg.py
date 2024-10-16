import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Test'
config.General.workArea = 'v1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'gensim.py'
config.JobType.maxMemoryMB = 2000

config.section_("Data")
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 2000 #number of files per jobs
config.Data.totalUnits =  -1 #number of event
config.Data.outLFNDirBase = '/store/user/user_name/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Test'

config.section_("Site")
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = []
