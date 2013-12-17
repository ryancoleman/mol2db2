#!/usr/bin/env python

#another script (compile_db_db2_sgearray.py) to make actual db2.gz files.
#makes a bunch of jobs for the cluster

import string
import sys
import os
import shutil

#db2file = "db.db2.gz"
#db2file = "db.more.db2.gz"
#db2file = "db.more.0.6.1200.db2.gz"
#db2file = "db.new.db2.gz"
#db2file = "db.new.0.6.1200.db2.gz"
#db2file = "db.more.0.6.1500.db2.gz"
#db2file = "db.more.0.5.1500.db2.gz"
#db2file = "db.more.0.4.2000.db2.gz"
#db2file = "db.more.0.4.2000.E.db2.gz"
db2file = "db.more.0.4.2000.cl.db2.gz"

sgeFile = '''#$ -S /bin/csh
#$ -cwd
#$ -q all.q
#$ -j yes
'''

scriptFileStart = '''#!/usr/bin/env python

#another script (compile_db_db2_sge.py) to make actual db2.gz files.
#makes a bunch of jobs for the cluster

import string
import sys
import os
import shutil

'''

scriptFileEnd = '''

sizeLimit = 100000000
tempName = "newTEMP." #what a great name
newTempDir = os.path.join(rootDir, "db2_inputs") #where newTEMP dirs end up
outTempDir = os.path.join(rootDir, db2outName) #where db2.gz files end up
targetCount = 1
targetOuts = []
targetOut = os.path.join(outTempDir, 
                       aTarget + "-" + string.zfill(targetCount, 5) + ".db2")
targetOuts.append(targetOut)
for db2gzFile in db2gzFiles:
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
    os.popen("zcat " + db2gzFile + " >> " +  targetOut)
  except KeyError:
    pass #some ligands don't have save.mol2 files from the build process
for targetCount, aDb2File in enumerate(targetOuts):
  os.popen("gzip -9 " + aDb2File)
'''


def readOtherDict(dictFileName):
  '''reads a dict of the format X->[1, 2, etc]'''
  bigDict = {}
  dictFile = open(dictFileName, 'r')
  try:
    for count, line in enumerate(dictFile):
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
rootDir = os.getcwd()
sizeLimit = 100000000
db2outName = "db2_outputs_0.4_cl"
tempName = "newTEMP." #what a great name
newTempDir = os.path.join(rootDir, "db2_inputs") #where newTEMP dirs end up
outTempDir = os.path.join(rootDir, db2outName) #where db2.gz files end up
scriptDir = os.path.join(rootDir, "compile_scripts") #where scripts go to build
if not os.path.exists(outTempDir):
  os.mkdir(outTempDir)
if not os.path.exists(scriptDir):
  os.mkdir(scriptDir)
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
  db2gzFiles = []
  for name, targets in id2tar.iteritems():
    if aTarget in targets: #the right target
      try:
        for aTempDir in id2name[name]:
          subDir = str(div10k(long(string.split(aTempDir, '.')[1])))
          db2gzFile = os.path.join(newTempDir, subDir, aTempDir, db2file)
          db2gzFiles.append(db2gzFile)
      except KeyError:
        pass #some ligands don't have save.mol2 files from the build process
  #make script now. 
  scriptName = "compile_" + aTarget + db2file  + "_db2.py"
  sgeName = "sge_" + aTarget + db2file + "_db2.csh"
  scriptFile = open(os.path.join(scriptDir, scriptName), 'w')
  scriptFile.write(scriptFileStart)
  #needs aTarget, db2gzFiles to be set between these
  scriptFile.write("db2outName = '" + db2outName + "'\n")
  scriptFile.write("aTarget = '" + aTarget + "'\n")
  scriptFile.write("rootDir = '" + rootDir + "'\n")
  scriptFile.write("db2gzFiles = ['")
  for aDb2gzFile in db2gzFiles[:-1]:
    scriptFile.write(aDb2gzFile + "', '")
  scriptFile.write(db2gzFiles[-1] + "' ]\n")
  scriptFile.write(scriptFileEnd)
  scriptFile.close()
  sgeFileOut = open(os.path.join(scriptDir, sgeName), 'w')
  sgeFileOut.write(sgeFile + "\n")
  sgeFileOut.write('#$ -o ' + 
            os.path.join(scriptDir, aTarget + db2file + ".test.log.txt") + "\n")
  sgeFileOut.write("cd " + scriptDir + "\n")
  sgeFileOut.write('chmod 755 ' + scriptName + '\n')
  sgeFileOut.write(scriptName + '\n')
  sgeFileOut.close()
  os.popen('qsub ' + os.path.join(scriptDir, sgeName))
