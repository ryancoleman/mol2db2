#!/usr/bin/env python2.7

#stupid script to run mol2db2 on a directory of files that are setup
# compiles whole thing into one test_this.db2.gz file at the end
# makes scripts and runs everything on the cluster

#run from the directory with subdirs on sgehead

mol2db2 = "/raid1/people/rgc/Source/bks_src/mol2db2_3.0_save/mol2db2.py"
scriptName = "mol2db2array.csh"
prefix = "newTEMP."
#inputMol2 = "db.mol2.gz"
#outputDb2 = "db.db2.gz"
#inputMol2 = "db.more.mol2.gz"
#outputDb2 = "db.more.db2.gz"
#inputMol2 = "db.more.0.6.1200.mol2.gz"
#outputDb2 = "db.more.0.6.1200.db2.gz"
#inputMol2 = "db.mol2.gz"
#outputDb2 = "db.new.db2.gz"
#inputMol2 = "db.more.0.6.1200.mol2.gz"
#outputDb2 = "db.new.0.6.1200.db2.gz"
#inputMol2 = "db.more.0.6.1500.mol2.gz"
#outputDb2 = "db.more.0.6.1500.db2.gz"
#inputMol2 = "db.more.0.5.1500.mol2.gz"
#outputDb2 = "db.more.0.5.1500.db2.gz"
#inputMol2 = "db.more.0.4.2000.mol2.gz"
#outputDb2 = "db.more.0.4.2000.db2.gz"
#inputMol2 = "db.more.0.4.2000.E.mol2.gz"
#outputDb2 = "db.more.0.4.2000.E.db2.gz"
inputMol2 = "db.more.0.4.2000.mol2.gz"
outputDb2 = "db.more.0.4.2000.cl.db2.gz"

cpulimit = str(7200) #seconds mol2db2 is allowed to run

numRunTogether = 24

import string, sys, os

def runTheseDirs(listDirNums, saveDir, parentDir):
  listToRun = []
  for dirNum in listDirNums:
    if not os.path.exists(os.path.join(parentDir, prefix + dirNum, outputDb2)):
      listToRun.append(dirNum)
  print parentDir, len(listToRun)
  if len(listToRun) > 0: #actually have jobs to run
    totalJobs = (len(listToRun) / numRunTogether) #no plus one?
    os.chdir(parentDir)
    scriptFile = open(scriptName, 'w') #open for writing
    scriptFile.write("#$ -S /bin/csh\n")
    scriptFile.write("#$ -cwd\n")
    scriptFile.write("#$ -q all.q\n")
    scriptFile.write("#$ -t 1-" + str(totalJobs) + "\n")
    scriptFile.write('#$ -l hostname="node-[56]-*"\n') #have to only use highmem
    scriptFile.write("#$ -j yes\n")
    scriptFile.write("#$ -o /dev/null\n") #these are idiotically static
    scriptFile.write("#$ -e /dev/null\n") #these are idiotocally static
    #scriptFile.write("#$ -o test.log.txt.$SGE_TASK_ID\n")
    for count in range(totalJobs):
      oneCount = count + 1
      scriptFile.write("if ($SGE_TASK_ID == " + str(oneCount) + ") then\n")
      lastOne = min((count * numRunTogether) + numRunTogether, len(listToRun))
      scriptFile.write("  set TASK_DIR = ( ")
      for dirCount in xrange(count * numRunTogether, lastOne):
        scriptFile.write(prefix + listToRun[dirCount] + " ")
      scriptFile.write(" ) \n") #end of taskdir line
      scriptFile.write("endif\n")
    scriptFile.write("foreach onetask ( $TASK_DIR ) \n")
    scriptFile.write("  cd $onetask\n")
    scriptFile.write("  rm test.log.txt\n")
    scriptFile.write("  limit cputime "+ cpulimit + "\n") #this is dumb, but use for now
    scriptFile.write("  " + mol2db2 + " -m " + inputMol2 + " -o " + outputDb2+ \
                     " -v -a >& test.log.txt \n") #reset is default here
    scriptFile.write("  unlimit cputime\n")
    scriptFile.write("  cd ..\n")
    scriptFile.write("end\n")
    scriptFile.close()
    os.popen("qsub " + scriptName)
  os.chdir(saveDir)
  
saveDir = os.getcwd()
workedOnce = False
listDirs = []
for dirName in os.listdir(saveDir):
  if dirName.find("TEMP") >= 0:
    workedOnce = True
    dirNumber = string.split(string.strip(dirName), ".")[1]
    listDirs.append(dirNumber)
if len(listDirs) > 0:
  runTheseDirs(listDirs, saveDir, saveDir + "/")
if not workedOnce:  #in a subdir for dud
  for dirName in os.listdir(saveDir):
    try:
      listDirs = []
      for subDirName in os.listdir(dirName):
        if subDirName.find("TEMP") >= 0:
          dirNumber = string.split(string.strip(subDirName), ".")[1]
          listDirs.append(dirNumber)
      if len(listDirs) > 0:
        runTheseDirs(listDirs, saveDir, saveDir + "/" + dirName + "/")
    except OSError:
      pass #not a big deal, wasn't a directory

