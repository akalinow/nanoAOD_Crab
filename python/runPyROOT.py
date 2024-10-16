#!/usr/bin/env python3
import os.path
import ROOT
from analysis import *

import sys
sys.path.insert(0,'./')

#import PSet and extract the input files list
from PSet import process
fileList = process.source.fileNames

#write list to file
with open("LFN_fileList.txt", "w") as file:
    for item in fileList:
        file.write("%s\n" % item)

# run edmFileUtil
os.system("edmFileUtil -d -F LFN_fileList.txt > PFN_fileList.txt")

#read back the list
fileList = []
with open("PFN_fileList.txt", "r") as file:
    for line in file:
        fileList.append(line.strip())

print(fileList)

analyse(fileList)



