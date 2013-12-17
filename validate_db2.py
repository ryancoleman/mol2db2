#!/usr/bin/env python2.7

#looks through directories, finds db.db2.gz files, validates that they are 
#good db2 files. currently this means length > 0 and ending in a E line.
#if they are bad move them to bad.db.db2.gz name

#run from the directory with subdirs on sgehead

import string, sys, os
import gzip
import shutil

#outName = "db.db2.gz"
#outName = "db.more.db2.gz"
#outName = "db.more.0.6.1200.db2.gz"
#outName = "db.newtest.db2.gz"
#outName = "db.new.db2.gz"
#outName = "db.new.0.6.1200.db2.gz"
#outName = "db.more.0.6.1500.db2.gz"
#outName = "db.more.0.5.1500.db2.gz"
#outName = "db.more.0.4.2000.db2.gz"
#outName = "db.more.0.4.2000.E.db2.gz"
outName = "db.more.0.4.2000.cl.db2.gz"

def runThisDir(dirName, saveDir):
  if os.path.exists(os.path.join(dirName, outName)):
    bad = False
    #print "running ", dirName
    #os.chdir(dirName)
    if 0 == os.path.getsize(os.path.join(dirName, outName)):
      bad = True
    if not bad: #yet
      tempRead = gzip.GzipFile(os.path.join(dirName, outName), 'r')
      try:
        for line in tempRead:
          if line[0] == 'S':
            if int(string.split(line)[1]) > 2000: #more than 2k sets not allowed
              bad = True
              break #no reason to do anymore
      except StopIteration: #EOF
        pass
      except IOError: #CRC check failing because the gz file format is bad
        bad = True
      tempRead.close()
      if (not bad) and string.strip(line) != "E": #also catches 0-length entries
        bad = True
    if bad:
      print dirName, " was bad"
      shutil.move(os.path.join(dirName, outName), \
                  os.path.join(dirName, "bad." + outName))
    else:
      pass
      #print dirName, " was good"

if -1 != string.find(sys.argv[0], "validate_db2.py"):
  saveDir = os.getcwd()
  workedOnce = False
  for dirName in os.listdir(saveDir):
    if dirName.find("TEMP") > 0:
      workedOnce = True
      dirNameFull = saveDir + "/" + dirName
      runThisDir(dirNameFull, saveDir)  
  if not workedOnce:  #in a subdir for dud
    for dirName in os.listdir(saveDir):
      try:
        for subDirName in os.listdir(dirName):
          if subDirName.find("TEMP") > 0:
            dirNameFull = saveDir + "/" + dirName + "/" + subDirName
            runThisDir(dirNameFull, saveDir)
      except OSError:
        pass #not a big deal, wasn't a directory

