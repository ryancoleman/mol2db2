#!/usr/bin/env python2.7

#use numeric utility to compute principal components of points
#written in Sharp lab by Ryan Coleman 2008-9
#edited, extended in Shoichet Lab by Ryan Coleman 2011-2015
#let penn & ucsf spend millions on lawyers, i've got work to do.

try: #to use numeric (old)
  from Matrix import Matrix
  from LinearAlgebra import eigenvectors
except ImportError: #use numpy (new)
  import numpy.oldnumeric as Numeric
  from numpy import matrix as Matrix
  from numpy.linalg import eig as eigenvectors
except ImportError:
  print "you do not have numpy or numeric installed under this python version"
  exit(1)
import geometry_basic #for getAverage of points, and dot product
import operator #for sorting tricks

def pca2d(pointList):
  '''sets up the pca for a list of points in 2d. solves eig problem'''
  matrixList = [[0,0],[0,0]] #init to 0
  avgPt = geometry_basic.getAverageArbitraryDimension(pointList, 2)
  diffs = [0,0]
  for point in pointList:
    for index in range(2):
      diffs[index] = point[index] - avgPt[index]
    #matrix is old plus a2 ab
    #                   ba b2 
    #          only compute upper diagonal now
    matrixList[0][0] += diffs[0]*diffs[0] #just hardcode it
    matrixList[0][1] += diffs[0]*diffs[1]
    matrixList[1][1] += diffs[1]*diffs[1]
  #make symmetric
  matrixList[1][0] = matrixList[0][1]
  actualMatrix = Matrix(matrixList)
  val,vec = eigenvectors(actualMatrix)
  return val, vec

def pca3d(pointList):
  '''sets up the pca for a list of points in 3d. solves eig problem'''
  matrixList = [[0,0,0],[0,0,0],[0,0,0]] #init to 0
  avgPt = geometry_basic.getAverage(pointList)
  diffs = [0,0,0]
  for point in pointList:
    for index in range(3):
      diffs[index] = point[index] - avgPt[index]
    #matrix is old plus a2 ab ac
    #                   ba b2 bc
    #                   ca cb c2 only compute upper diagonal now
    matrixList[0][0] += diffs[0]*diffs[0] #just hardcode it
    matrixList[0][1] += diffs[0]*diffs[1]
    matrixList[0][2] += diffs[0]*diffs[2]
    matrixList[1][1] += diffs[1]*diffs[1]
    matrixList[1][2] += diffs[1]*diffs[2]
    matrixList[2][2] += diffs[2]*diffs[2]
  #make symmetric
  matrixList[1][0] = matrixList[0][1]
  matrixList[2][0] = matrixList[0][2]
  matrixList[2][1] = matrixList[1][2]
  actualMatrix = Matrix(matrixList)
  val,vec = eigenvectors(actualMatrix)
  return val, vec

def flatten(pointListList):
  '''flattens an input list of lists of lists into a list of lists, removing the
  least significant list (lowest order)'''
  #let's go ahead and flatten the input list!
  flattenedList = []
  for pointList in pointListList:
    flatTemp = [pts for sublist in pointList for pts in sublist]
    flattenedList.append(flatTemp)
  return flattenedList

def pcaN3d(pointListList):
  '''sets up the pca for a list of list of points in 3d. solves eig problem.
  the pointListList is a list of sets of points of the same length.'''
  length = len(pointListList[0]) * 3 #know 3d points times length of list of pts
  oneList = [0. for count in range(length)] #bunch of 0.0s
  matrixList = [oneList[:] for count in range(length)] #copy list of 0.0s
  #let's go ahead and flatten the input list!
  flattenedList = flatten(pointListList)
  avgPt = geometry_basic.getAverageArbitraryDimension(flattenedList, dimension=length)
  diffs = oneList[:]
  for point in flattenedList: #this point has length dimensions
    for index in xrange(length):
      diffs[index] = point[index] - avgPt[index]
    #matrix is old plus a2 ab ac ad etc
    #                   ba b2 bc bd etc
    #                   ca cb c2 cd etc
    #                   da db dc d2 etc
    #                   etc etc etc etc etc (etc)
    #only compute upper diagonal now
    #can't get away with hardcoding anymore
    for row in xrange(length):
      for col in xrange(row, length):
        matrixList[row][col] += diffs[row] * diffs[col]
  #make symmetric
  for row in xrange(length):
    for col in xrange(0, row):
      matrixList[row][col] = matrixList[col][row]
  #print "mat" #debugging madness
  #for row in xrange(length):
  #  for col in xrange(length):
  #    print matrixList[row][col],
  #  print " "
  actualMatrix = Matrix(matrixList)
  val,vec = eigenvectors(actualMatrix)
  return val, vec

def findLongestProjectedDirection(pointListList):
  '''does pca, gets the eigenresult, find the longest direction, returns it.
  done for a list of 3d points, projects them onto 1d using the eigenvector 
  and dot product.'''
  eigenvalues, eigenvectors = pcaN3d(pointListList) 
  maxVal, maxIndex = 0,0
  try:
    for index in range(len(eigenvalues)):
      if maxVal < eigenvalues[index]:
        maxVal = eigenvalues[index]
        maxIndex = index
    maxVec = eigenvectors[maxIndex]
  except TypeError: #caused by complex numbers
    for index in range(len(eigenvalues)):
      if maxVal < eigenvalues[index].real: #real prevents imaginary problems
        maxVal = eigenvalues[index].real
        maxIndex = index
    maxVec = eigenvectors[maxIndex].real
  try:
    maxVecRet = maxVec.tolist() #stupid numpy new2010
  except AttributeError:
    maxVecRet = maxVec
  if 1 == len(maxVecRet): 
    maxVecRet = maxVecRet[0]
  return maxVecRet

def findBiggestGapSplit(projectedPts):
  '''finds the biggest gap (diff between adjacent pts), splits on either side'''
  sortedPts = projectedPts[:] #copy
  sortedPts.sort()
  biggestGap, biggestIndex = 0., -1
  for index in xrange(len(sortedPts) - 1):
    diff = sortedPts[index + 1] - sortedPts[index]
    if diff > biggestGap:
      biggestIndex = index
      biggestGap = diff
  print biggestIndex, biggestGap
  breakPoint = sortedPts[biggestIndex]
  #print "index, length of split", biggestIndex, len(sortedPts) #debugging
  splits = ([], []) 
  for index, point in enumerate(projectedPts):
    if point <= breakPoint:
      splits[0].append(index)
    else:
      splits[1].append(index)
  return splits

def findBisectiveSplit(projectedPts, avgPoint):
  '''for a given list of points, find the best split, based on average.'''
  splits = ([],[]) 
  for index, point in enumerate(projectedPts):
    if point <= avgPoint:
      splits[0].append(index)
    else:
      splits[1].append(index)
  if 0 == len(splits[1]): #nothing in second half. try splitting w/o =
    splits = ([],[]) 
    for index, point in enumerate(projectedPts):
      if point < avgPoint: #just < not <=
        splits[0].append(index)
      else:
        splits[1].append(index)
  return splits

def findProjectAndSplit(pointListList, altSplit=False):
  '''finds the eigs, projects the points, splits into 2 sub lists. returns
  indices from orig pointListList split into 2 groups. 
  altsplit true means biggest gap splitting
  altsplit false means bisective splitting (based on average)
  maybe want a 3rd option to split on median (equal subgroups)'''
  maxVecRet = findLongestProjectedDirection(pointListList)
  flattenedList = flatten(pointListList)
  projectedPts = []
  for pointList in flattenedList:
    projectedPt = geometry_basic.dot(maxVecRet, pointList)
    try:
      projectedPt = projectedPt.real #in case it is complex
    except AttributeError:
      pass #this is okay, not a complex number
    projectedPts.append(projectedPt)
  if not altSplit: #use bisective splitting, i.e. split based on mean
    avgPoint = geometry_basic.getAverage1(projectedPts)
    indicesSplit = findBisectiveSplit(projectedPts, avgPoint)
  else:
    indicesSplit = findBiggestGapSplit(projectedPts)
  return indicesSplit

def findLongestDirection(pointList):
  '''does pca, gets the eigenresult, find the longest direction, returns it'''
  eigenvalues, eigenvectors = pca3d(pointList)
  maxVal, maxIndex = 0,0
  try:
    for index in range(len(eigenvalues)):
      if maxVal < eigenvalues[index]:
        maxVal = eigenvalues[index]
        maxIndex = index
    maxVec = eigenvectors[maxIndex]
  except TypeError: #caused by complex numbers
    for index in range(len(eigenvalues)):
      if maxVal < eigenvalues[index].real: #real prevents imaginary problems
        maxVal = eigenvalues[index].real
        maxIndex = index
    maxVec = eigenvectors[maxIndex].real
  try:
    maxVecRet = maxVec.tolist() #stupid numpy new2010
  except AttributeError:
    maxVecRet = maxVec
  return maxVecRet
    
def sortDirections(pointList):
  '''finds the eigenvalues and eigenvectors, sorts based on eigenvalue and 
  returns list of direction vectorsin descending order'''
  eigenvalues, eigenvectors = pca3d(pointList)
  maxVal, maxIndex = 0, 0
  #print eigenvalues, eigenvectors
  try:
    newEigList = []
    for index in xrange(len(eigenvalues)):
      newEigList.append((eigenvalues[index], eigenvectors[index]))
    newEigList.sort(key=operator.itemgetter(0))
    newEigList.reverse()
  except TypeError: #imaginary problem
    newEigList = []
    for index in xrange(len(eigenvalues)):
      newEigList.append((eigenvalues[index].real, eigenvectors[index].real))
    newEigList.sort(key=operator.itemgetter(0))
    newEigList.reverse()
  retEigVecList = []
  for eigenvalue, eigenvector in newEigList:
    try:
      newEig = eigenvector.tolist() #stupid numpy new 2010 [0] before?
      retEigVecList.append(newEig)
    except AttributeError: #numpy is actually not being run, not a problem
      retEigVecList.append(eigenvector)
  return retEigVecList

def findLongestDimension(pointList):
  '''calls direction, returns length in that direction'''
  direction = findLongestDirection(pointList)
  #direction is a unit vector, so can do scalar projection, i.e.
  #dot product of each point with direction gives the length in that direction
  #and is negative if in the opposite direction, so just find max-min and return
  firstPoint = pointList[0]
  try:
    if 1 == len(direction):
      direction = direction[0] #numpy/numeric return differently
  except TypeError:
    pass
  min = geometry_basic.dot(firstPoint, direction)
  max = min #same for now
  for point in pointList[1:]: #already done first point
    newScalar = geometry_basic.dot(point, direction)
    try:
      if newScalar < min:
        min = newScalar
      if newScalar > max:
        max = newScalar
    except TypeError:
      newScalar = newScalar.real
      min = min.real
      max = max.real
      if newScalar < min:
        min = newScalar
      if newScalar > max:
        max = newScalar
  return max-min   

def findDimensions(pointList):
  '''finds the dimension in the 3 principal directions, longest first'''
  dimensions = []
  directions = sortDirections(pointList)
  for direction in directions:
    firstPoint = pointList[0]
    try:
      if 1 == len(direction):
        direction = direction[0] #stupid numpy
    except TypeError:
      pass
    min = geometry_basic.dot(firstPoint, direction)
    max = min #same for now
    for point in pointList[1:]: #already done first point
      newScalar = geometry_basic.dot(point, direction)
      try:
        if newScalar < min:
          min = newScalar
        if newScalar > max:
          max = newScalar
      except TypeError:
        newScalar = newScalar.real
        min = min.real
        max = max.real
        if newScalar < min:
          min = newScalar
        if newScalar > max:
          max = newScalar
    dimensions.append(max - min)
  return dimensions 

#disable testing now
'''
import sys,string
if -1 != string.find(sys.argv[0], "pca.py"):
  print findLongestDimension([[1,1,1],[0,0,0],[-1,-1,-1]]) #stupid test
  print findLongestDimension([[10,10,10],[9,9,9],[8,8,8]]) #stupid test
  print findLongestDimension([[10,9,1],[9,9,9],[8,7,8]]) #stupid test
  print findDimensions([[10,9,1],[9,9,9],[8,7,8],[1,1,1]]) #stupid test
  #print pcaN3d([[[10,9,2],[9,9,7],[8,6,8],[1,1,4]],[[10,9,1],[9,9,9],[8,7,8],[1,1,1]]])
  #print findLongestProjectedDirection([[[10,4,2],[9,8,7],[8,6,8],[1,2,4]],[[10,9,1],[9,9,9],[8,7,8],[1,1,1]],[[-10,-9,-1],[-9,-9,-9],[-8,-7,-8],[-1,-1,-1]],[[90,90,10],[90,80,80],[70,70,80],[8,9,10]]])
  print findProjectAndSplit([[[10,4,2],[9,8,7],[8,6,8],[1,2,4]],[[10,9,1],[9,9,9],[8,7,8],[1,1,1]],[[-10,-9,-1],[-9,-9,-9],[-8,-7,-8],[-1,-1,-1]],[[90,90,10],[90,80,80],[70,70,80],[8,9,10]]])
  print findProjectAndSplit([[[10,4,2],[9,8,7],[8,6,8],[1,2,4]],[[10,9,1],[9,9,9],[8,7,8],[1,1,1]],[[-10,-9,-1],[-9,-9,-9],[-8,-7,-8],[-1,-1,-1]]])
'''
