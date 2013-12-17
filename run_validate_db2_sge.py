#!/usr/bin/env python2.7

#runs as an sge job all subdirs to check for valid db2 files
#run from the directory with subdirs on sgehead
#Ryan G. Coleman, BKS lab. 2010.

#IMPORTANT NOTE: doing this pretty much saturates disk reads/writes. 
#running 42 of these jobs does anyway. don't run this without monitoring the 
#disk usage. talk to your friendly sysadmin before starting fires.

import string, sys, os

scriptName = "validate_sge.csh"

validate = "/raid1/people/rgc/Source/bks_src/mol2db2_3.0_save/validate_db2.py"

prefix = "newTEMP."

def runTheseDirs(listDirs, saveDir, parentDir):
  for subDir in listDirs:
    print os.path.join(parentDir, subDir)
    os.chdir(os.path.join(parentDir, subDir))
    scriptFile = open(scriptName, 'w') #open for writing
    scriptFile.write("#$ -S /bin/csh\n")
    scriptFile.write("#$ -cwd\n")
    scriptFile.write("#$ -q all.q\n")
    scriptFile.write('#$ -l hostname="node-[56]-*"\n') #have to only use highmem
    scriptFile.write("#$ -j yes\n")
    scriptFile.write("#$ -o validate_log.txt\n") 
    scriptFile.write(validate + " \n") #reset is default here
    scriptFile.close()
    os.popen("qsub " + scriptName)
  os.chdir(saveDir)

if -1 != string.find(sys.argv[0], "run_validate_db2_sge.py"):
  saveDir = os.getcwd()
  workedOnce = False
  listDirs = []
  for dirName in os.listdir(saveDir):
    if dirName.find(prefix) >= 0:
      workedOnce = True
      listDirs.append(saveDirName)
      break
  if len(listDirs) > 0:
    runTheseDirs(listDirs, saveDir, os.path.dirname(saveDir))
  if not workedOnce:  #in a subdir for dud
    listDirs = []
    for dirName in os.listdir(saveDir):
      try:
        for subDirName in os.listdir(dirName):
          if subDirName.find(prefix) >= 0:
            listDirs.append(dirName)
            break
      except OSError:
        pass #not a big deal, wasn't a directory
    if len(listDirs) > 0:
      runTheseDirs(listDirs, saveDir, saveDir)

