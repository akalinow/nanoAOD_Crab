## Introduction

This is a script for submitting python code for nanoAOD analysis on CMS samples available in [DAS](https://cmsweb.cern.ch/das/). The sctipst subn it a Crab job that will execute a function

```
analyse(fileList)
```

from [analysis.py](analysis.py) file.

## Installation instructions:

* fetch this repository:

``` 
git clone https://github.com/akalinow/nanoAOD_Crab
```
## Run instructions

* setup any CMSSW area to have ROOT available
```
cd CMSSW_X_Y_Z/src
cmsenv
cd -
```

The jobs are submitted with a single command:

```
cd nanoAOD_Crab
./submitJobs.py
```

The [submitJobs.py](submitJobs.py) script contains following control [parameters](https://github.com/akalinow/nanoAOD/blob/main/submitJobs.py#L91-L105)

* **eventsPerJob** - number of events to be generated per job. 
* **outLFNDirBase** - Logical File Name base directory for storing the output on SE, for example */store/user/akalinow/Data/*
  The results will be stored in subdirectories created automatically.
* **storage_element** - name of the SE, for example *T3_CH_CERNBOX*
* **outputDatasetTag** - tag added to the name of the dataset in the name the job, for example *akalinow_crab_EGamma0_Run2023C_22Sep2023_v1_v1_PyROOT_test2*,
                         where the tag is *test2*
* **outputFileList** - list of files produced by the analysis. Those files will be saved to SE.

* **runLocal** - a flag for running a local test. If set to *True* the Crab job will not be submitted. The python script will be run on a list of local files
* **localFiles** - a list of local files. Used only when `runLocal=True`