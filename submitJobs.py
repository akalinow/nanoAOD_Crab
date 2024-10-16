#!/usr/bin/env python3

import os, re
import math
import urllib
import glob

import CRABClient
from CRABAPI.RawCommand import crabCommand

import sys
sys.path.insert(0,'./python')

from crab_cfg import *

import FWCore.ParameterSet.Config as cms

#########################################
#########################################
def prepareCMSSWCfg(inputFileList):
    
    process = cms.Process('PyROOTProxy')
    process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(inputFileList),
                            inputCommands = cms.untracked.vstring('drop *')
    )   

    out = open('tmpConfig.py','w')
    out.write(process.dumpPython())
    out.close()
#########################################
#########################################
def prepareCrabCfg(dataset,
                   eventsPerJob,
                   outLFNDirBase,
                   storage_element,
                   outputDatasetTag,
                   outputFileList,
                   runLocal,
                   localFiles):
    
    shortName = dataset.split("/")[1]
    if dataset.split("/")[2].find("Run202")!=-1:
        shortName += "_"+dataset.split("/")[2]
    requestName = shortName+"_PyROOT_"+outputDatasetTag
    requestName = requestName.replace("-","_")

    inputFileList = []
    
    if runLocal:
            inputFileList = localFiles

    prepareCMSSWCfg(inputFileList)       
                                                              
    ##Modify CRAB3 configuration
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'tmpConfig.py'
    config.JobType.disableAutomaticOutputCollection = True
    config.JobType.scriptExe = 'python/runPyROOT.py'
    config.JobType.inputFiles = ['python/analysis.py', 'xml/FrameworkJobReport.xml']
    config.JobType.outputFiles = outputFileList
    
    config.General.requestName = requestName
    config.General.workArea = "Crab_tasks"

    config.Data.inputDBS = 'global'
    config.Data.inputDataset = dataset
    config.Data.outLFNDirBase = outLFNDirBase+outputDatasetTag
    config.Data.publication = True
    config.Data.outputDatasetTag = outputDatasetTag
    
    config.Site.storageSite = storage_element
    
    config.Data.unitsPerJob = eventsPerJob
    config.Data.totalUnits = -1

    out = open('crabTmp.py','w')
    out.write(config.pythonise_())
    out.close()

    if runLocal:
        os.system("cp tmpConfig.py PSet.py")
        os.system("python3 python/runPyROOT.py")

    else:
        os.system("crab submit -c crabTmp.py")

    os.system("rm -f tmpConfig.py crabTmp.py PSet.py")
#########################################
#########################################
##Those are the steering parameters
eventsPerJob = int(1E7)
outLFNDirBase = "/store/user/akalinow/Data/WW/"
storage_element="T3_CH_CERNBOX"

outputDatasetTag = "test3"
outputFileList = ['vbstt.root', 'vbsttee.root',  'vbsttmue.root',  'vbsttmumu.root']

#List od DBS datasets to be analyzed
datasets = {'/EGamma0/Run2023C-22Sep2023_v1-v1/NANOAOD',
            #'/EGamma1/Run2023C-22Sep2023_v1-v1/NANOAOD'
            }

runLocal = False
##List of locally accesible files to be analysed by running on local machine.
path = "/home/akalinow/scratch/CMS/Data/EGamma0/Run2023C-22Sep2023_v1-v1/NANOAOD/"
localFiles = glob.glob(path+"*.root")
localFiles = ["file:"+aFile for aFile in localFiles]
########################################################
########################################################
for aDataset in datasets:
    prepareCrabCfg(dataset=aDataset,                   
                   eventsPerJob=eventsPerJob,
                   outLFNDirBase = outLFNDirBase, 
                   storage_element=storage_element,
                   outputDatasetTag = outputDatasetTag,
                    outputFileList = outputFileList,
                   runLocal=runLocal,
                   localFiles=localFiles)
########################################################

