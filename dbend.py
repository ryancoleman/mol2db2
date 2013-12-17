#!/usr/bin/env python2.7

#this script is a replacement for dbend.csh for db2.gz files, much more sensible

import string
import sys
import os
import shutil
import glob

db2file = "db.db2.bz2"

sizeLimit = 100000000
outTempDir = os.getcwd()
aTarget = "dbout"
targetCount = 1
targetOuts = []
targetOut = os.path.join(outTempDir, 
                       aTarget + "-" + string.zfill(targetCount, 5) + ".db2")
print targetOut
targetOuts.append(targetOut)
for oneDb2file in glob.iglob(os.path.join("*", "TEMP*.mol2", db2file)):
  try:
    curSize = os.path.getsize(targetOut)
  except OSError:
    curSize = 0 #doesn't exist yet, can't have a size
  if sizeLimit < curSize:  
    targetCount += 1 #increment for next
    targetOut = os.path.join(outTempDir, 
                       aTarget + "-" + string.zfill(targetCount, 5) + ".db2")
    print targetOut
    targetOuts.append(targetOut)
  os.popen("bzcat " + oneDb2file + " >> " +  targetOut)
for targetCount, aDb2File in enumerate(targetOuts):
  if 0 == targetCount % 2: #hack to run 2 gzips at once
    print "gzipping " + aDb2File
    os.popen("gzip -9 " + aDb2File + " & ")
  else:
    print "gzipping " + aDb2File
    os.popen("gzip -9 " + aDb2File)
