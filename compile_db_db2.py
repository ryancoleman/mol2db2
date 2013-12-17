#!/usr/bin/env python

#another script (compile_db_db2.py) to make actual db2.gz files.

import string
import sys
import os
import shutil

#db2file = "db.db2.gz"
#db2file = "db.more.db2.gz"
#db2file = "db.more.0.6.1200.db2.gz"
db2file = "db.new.db2.gz"

def readOtherDict(dictFileName):
  '''reads a dict of the format X->[1, 2, etc]'''
  bigDict = {}
  dictFile = open(dictFileName, 'r')
  try:
    for line in dictFile:
      tokens = string.split(line)
      if tokens[0] not in bigDict.keys():
        bigDict[tokens[0]] = tokens[1:]
      else: #need to append
        bigDict[tokens[0]].extend(tokens[1:]) #adds other ligand/decoy sets
  except StopIteration:
    pass #EOF okay
  dictFile.close()
  return bigDict

def div10k(temp):
  '''used to break into 10000 size chunks'''
  return temp / 10000

target = os.getcwd()
sizeLimit = 100000000
tempName = "newTEMP." #what a great name
newTempDir = os.path.join(os.getcwd(), "db2_inputs") #where newTEMP dirs end up
outTempDir = os.path.join(os.getcwd(), "db2_outputs") #where db2.gz files end up
if not os.path.exists(outTempDir):
  os.mkdir(outTempDir)
id2name = {} #E01234567 newTEMP.1 newTEMP.2
id2tar = {}  #E01234567 pde4a_lig pde4b_lig
dictFileName1 = os.path.join(os.getcwd(), "db2_name") #master file
dictFileName2 = os.path.join(os.getcwd(), "db2_tar") #master file
if os.path.exists(dictFileName1): #file exists, read stuff from it
  id2name = readOtherDict(dictFileName1)
if os.path.exists(dictFileName2): #file exists, read stuff from it
  id2tar = readOtherDict(dictFileName2)
targetList = []
for values in id2tar.itervalues():
  for aValue in values:
    if aValue not in targetList:
      targetList.append(aValue)
for aTarget in targetList:
  print aTarget
  targetCount = 1
  targetOuts = []
  targetOut = os.path.join(outTempDir, 
                       aTarget + "-" + string.zfill(targetCount, 5) + ".db2")
  targetOuts.append(targetOut)
  for name, targets in id2tar.iteritems():
    if aTarget in targets: #the right target
      try:
        curSize = os.path.getsize(targetOut)
      except OSError:
        curSize = 0 #doesn't exist yet, can't have a size
      if sizeLimit < curSize:  
        targetCount += 1 #increment for next
        targetOut = os.path.join(outTempDir, 
                       aTarget + "-" + string.zfill(targetCount, 5) + ".db2")
        targetOuts.append(targetOut)
      try:
        for aTempDir in id2name[name]:
          subDir = str(div10k(long(string.split(aTempDir, '.')[1])))
          db2gzFile = os.path.join(newTempDir, subDir, aTempDir, db2file)
          os.popen("zcat " + db2gzFile + " >> " +  targetOut)
      except KeyError:
        pass #some ligands don't have save.mol2 files from the build process
  for targetCount, aDb2File in enumerate(targetOuts):
    if 0 == targetCount % 2: #hack to run 2 gzips at once
      os.popen("gzip -9 " + aDb2File + " & ")
    else:
      os.popen("gzip -9 " + aDb2File)
