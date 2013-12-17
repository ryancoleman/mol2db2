#!/usr/bin/env python2.7

#Ryan G. Coleman, Brian K. Shoichet Lab

import string
import sys
import optparse
import sybyl2dock
import mol2
import hierarchy
import solv
import clash
import hydrogens
import time
import math


def mol2db2(options):
  '''function that does all the actual work you may want to do to convert a
  mol2 file and solv file into a db2 file.'''
  if options.timeit:
    timeStart = time.time()
  else:
    timeStart = None
  mol2data = mol2.Mol2(options.mol2file, nameFileName=options.namefile)
  mol2data.convertDockTypes(options.atomtypefile)
  mol2data.addColors(options.colortablefile)
  solvdata = solv.Solv(options.solvfile)
  if options.verbose:
    print "names:", mol2data.name, mol2data.protName, mol2data.smiles, \
                   mol2data.longname
    print "dock atom types:", mol2data.atomType
    print "dock atom type numbers:", mol2data.dockNum
    print "dock color type numbers:", mol2data.colorNum
  clashDecider = clash.Clash(options.clashfile)
  hydrogenRotater = hydrogens.Hydrogens(options.hydrogenfile)
  if options.timeit:
    timeReadIn = time.time()
    print "time to read in files:", timeReadIn-timeStart
  #0th step is to do the hydrogen operation directly on the mol2data.
  #done here so the hierarchy knows exactly how many conformations it must deal
  #with, not some estimate. and so each hierarchy has the max # of allowed sets
  #again, without guessing.
  if options.rotateh or options.reseth:
    hydrogenRotater.findTerminalHydrogens(mol2data)
    if options.verbose:
      print mol2data.hydrogensToRotate, " hydrogens need rotated"
    if options.reseth and 0 < mol2data.hydrogensToRotate:
      hydrogenRotater.resetHydrogens(mol2data)
    if options.rotateh and 0 < mol2data.hydrogensToRotate:
      mol2data = hydrogenRotater.rotateHydrogens(mol2data)
  if options.timeit:
    hydTime = time.time()
    print "time to move hydrogens:", hydTime-timeReadIn
  if options.verbose:
    print len(mol2data.atomXyz), " conformations in input"

  #Generators to pipeline hierarchy generation
  def hierarchyDataGenerator():
    try:
      yield hierarchy.Hierarchy(mol2data, clashDecider, \
                        tolerance=options.tolerance, \
                        verbose=options.verbose, timeit=options.timeit, \
                        limitset=options.limitset, \
                        limitconf=options.limitconf, \
                        limitcoord=options.limitcoord, solvdata=solvdata)
    except hierarchy.TooBigError as limitError:
      breaks = hierarchy.computeBreaks(limitError, options)
      origConfsPer = max(len(mol2data.atomXyz)/(breaks + 1), 1) #at least 1
      if options.verbose:
        print "splitting original input conformations into", breaks + 1, "parts"
      for snap in xrange(breaks + 1): #extra one to do leftovers
        newMol2data = mol2data.copy() #copy orig before hydrogen rotations
        first = origConfsPer * snap
        last = origConfsPer * (snap + 1)
        newMol2data.keepConfsOnly(first, last)
        #print "atomXyz breaks", first, last, len(newMol2data.atomXyz), len(mol2data.atomXyz)
        if len(newMol2data.atomXyz) > 0:
          hierarchyData = hierarchy.Hierarchy(newMol2data, clashDecider, \
                          tolerance=options.tolerance, \
                          verbose=options.verbose, timeit=options.timeit, \
                          solvdata=solvdata)
          yield hierarchyData

  if options.timeit:
    timeHier = time.time()
    print "time to (start) construction hierarchy (subtotal):", timeHier-hydTime
  return timeStart, hierarchyDataGenerator()

def mol2db2writeDb2(options, timeStart, hierarchyDatas):
  '''does the writing of the output files. separate so you can make but
  not write. requires timeStart (can be None if no times wanted) and
  all the hierarchyDatas from the mol2db2 function'''
  if options.timeit:
    timeHier = time.time()
  # Clear db2gzfile
  open(options.db2gzfile, 'wb').close()
  for hierarchyData in hierarchyDatas:
    hierarchyData.write(db2gzFileName=options.db2gzfile, \
                        verbose=options.verbose,\
                        timeit=options.timeit, limitset=options.limitset, \
                        writeMode='a') #append so we don't overwrite
    yield hierarchyData
  if options.timeit:
    timeHierOut = time.time()
    print "time to save db2:", timeHierOut-timeHier
  if options.timeit:
    timeFinal = time.time()
    print "time to do everything (total):", timeFinal-timeStart

if -1 != string.find(sys.argv[0], "mol2db2.py"):
  description = "Convert multimol2 and solv files into db2 files"
  usage = "%prog [options]"
  parser = optparse.OptionParser(usage=usage, description=description)
  parser.add_option("-v", "--verbose",
                    action="store_true", dest="verbose", default=False,
                    help="lots of debugging output")
  parser.add_option("-a", "--time",
                    action="store_true", dest="timeit", default=False,
                    help="timing output")
  parser.add_option("-m", "--mol2", type="string", action="store", \
                    dest="mol2file", default="db.mol2.gz", \
                    help="multimol2 input file, (default: %default)")
  parser.add_option("-n", "--namemol2", type="string", action="store", \
                    dest="namefile", default="name.txt", \
                    help="mol2 name file, (default: %default)")
  parser.add_option("-s", "--solv", type="string", action="store", \
                    dest="solvfile", default="db.solv", \
                    help="solv input file, (default: %default)")
  parser.add_option("-o", "--db", type="string", action="store", \
                    dest="db2gzfile", default="db.db2.gz", \
                    help="db2.gz output file, (default: %default)")
  parser.add_option("-t", "--atomtype", type="string", action="store", \
                    dest="atomtypefile", default=None, \
                    help="atom type conversion file, (default: %default)")
  parser.add_option("-c", "--colortable", type="string", action="store", \
                    dest="colortablefile", default=None, \
                    help="color table conversion file, (default: %default)")
  parser.add_option("-d", "--clash", type="string", action="store", \
                    dest="clashfile", default=None, \
                    help="clash checking parameter file, (default: %default)")
  parser.add_option("-y", "--hydrogens", type="string", action="store", \
                    dest="hydrogenfile", default=None, \
                   help="terminal hydrogen parameter file, (default: %default)")
  parser.add_option("-r", "--noreseth",
                    action="store_false", dest="reseth", default=True,
                    help="reset the planar terminal hydrogen positions," + \
                         " (default: %default)")
  parser.add_option("-z", "--norotateh",
                    action="store_false", dest="rotateh", default=True,
                    help="rotate the terminal hydrogen positions," + \
                         " (default: %default)")
  #yes you can do reset and rotate in any combination!
  parser.add_option("-x", "--disttol", type="float", action="store", \
                    dest="tolerance", default=0.001, \
                    help="distance tolerance in angstroms, (default: %default)")
  parser.add_option("-l", "--limitset", type="long", action="store", \
                    dest="limitset", default=500, \
                    help="limit on the number of sets, (default: %default)")
  parser.add_option("--limitconf", type="long", action="store", \
                    dest="limitconf", default=2000, \
                    help="limit on the number of confs, (default: %default)")
  parser.add_option("--limitcoord", type="long", action="store", \
                    dest="limitcoord", default=40000, \
                    help="limit on the number of coords, (default: %default)")
  options, args = parser.parse_args() #default reads from argv[1:]
  if 0 != len(args):
    parser.error("mol2db2.py takes no positional arguments" +
                 "  Use --help for more information")
  else:
    #print options
    if options.verbose:
      print "verbose debugging requested."
      print options
    timeStart, hierarchyDatas = mol2db2(options) #main program call
    finishedHierarchyDatas = mol2db2writeDb2(options, timeStart, hierarchyDatas) #write output
    for hierarchyData in finishedHierarchyDatas:
      if options.verbose:
        print "Cleaning up"
      del hierarchyData
