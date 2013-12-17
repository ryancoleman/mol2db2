#ryan g coleman, ryangc@mail.med.upenn.edu
#copyright 2006-7 ryan g coleman, kim sharp crystal.med.upenn.edu
#geometric primitives like distance functions and such

import math

def distL2(a,b):
  '''no error checking, very fast, should use everywhere'''
  sum = 0.
  for count in xrange(len(a)):
    sum += (b[count]-a[count])**2.
  return math.sqrt(sum) # is this faster than **0.5?

def distL2Squared3(a, b):
  '''no error checking, unrolled loop'''
  return (b[0] - a[0])**2. + (b[1] - a[1])**2. + (b[2] - a[2])**2.

def distL2Squared(a,b):
  '''no error checking, very fast, should use everywhere, doesn't square root'''
  sum = 0.
  for count in xrange(len(a)):
    sum += (b[count]-a[count])**2.
  return sum

def dist(a,b,metric='L2'):
  '''a and b should be lists of equal length (any dimension)
  calculates distance needed and returns it (L1,L2,LINF,L2SQUARED).
  these new versions are twice the speed of using list comprehensions.'''
  if metric == 'L2':
    sum = 0.
    for count in xrange(len(a)):
      sum += (b[count]-a[count])**2.
    return sum**0.5
  elif metric == 'LINF':
    max = 0.
    for count in xrange(len(a)):
      new = abs(b[count]-a[count])
      if new > max:
        max = new
    return max
  elif metric == 'L2SQUARED':
    sum = 0.
    for count in xrange(len(a)):
      sum += (b[count]-a[count])**2.
    return sum
  elif metric == 'L1':
    sum = 0.
    for count in xrange(len(a)):
      sum += abs(b[count]-a[count])
    return sum

def longestAndMeanDist(pts):
  '''given a list of points, finds the largest distance between any 2. also 
  finds mean distance between all pairs. returns both, in that order.'''
  longestDist = 0.
  sumDists, countDists = 0., 0
  for indexOne, ptOne in enumerate(pts):
    for ptTwo in pts[indexOne+1:]: #no duplicates, minimal looping
      thisDist = distL2(ptOne, ptTwo)
      longestDist = max(thisDist, longestDist)
      sumDists += thisDist
      countDists += 1
  return longestDist, sumDists/float(countDists)

def getAngle(a,b):
  '''helper function for triangle interior, returns angle between two vectors'''
  #print "a =", a
  #print "b =", b
  #return math.acos(max(-1.,min(1.,(sum( a[i]*b[i] for i in range(len(a)) ))/(((sum( x**2. for x in a ))**0.5)*((sum( x**2. for x in b ))**0.5)))))
  ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2] #all inlined for speed
  aSquared = a[0]**2. + a[1]**2. + a[2]**2.
  bSquared = b[0]**2. + b[1]**2. + b[2]**2.
  #ab = 0. #tons of debugging here
  #aSquared = 0.
  #bSquared = 0.
  #for index in xrange(len(a)):
  #  ab += a[index] * b[index]
  #  aSquared += a[index]**2.
  #  bSquared += b[index]**2.
  return math.acos(max(-1.,min(1.,(ab)/(((aSquared)**0.5)*((bSquared)**0.5)))))

def calcTriAreaList(abc):
  '''uses heron's formula'''
  a,b,c=abc #unpack
  dists = [distL2(a,b),distL2(b,c),distL2(a,c)]
  s = (dists[0] + dists[1] + dists[2])*0.5
  triArea = (s*(s-dists[0])*(s-dists[1])*(s-dists[2]))**(0.5)
  return triArea

def calcTriArea(a,b,c): #3 points in 3d
  '''uses heron's formula'''
  dists = [distL2(a,b),distL2(b,c),distL2(a,c)]
  s = (dists[0] + dists[1] + dists[2])*0.5
  triArea = (s*(s-dists[0])*(s-dists[1])*(s-dists[2]))**(0.5)
  return triArea

def getVector(a, b):
  '''does a-b, returns'''
  return [a[i]-b[i] for i in range(len(a))]

def getNormalVector(a, b):
  '''normal(a-b)'''
  return normalizeVector(getVector(a, b))

def length(vector):
  '''vector length'''
  total = 0.
  for coord in vector:
    total += coord**2.
  total = total**0.5
  return total

def normalizeVector(vector):
  '''divides each by the total components squared'''
  total = 0.
  for coord in vector:
    total += coord**2.
  total = total**0.5
  newVect = []
  for coord in vector:
    newVect.append(coord/total)
  return newVect

def dot(x, y):
  '''gives dot product of two vectors of any dimension, assumes same length'''
  dot = 0.
  for index in range(len(x)):
    dot += x[index] * y[index]
  return dot

def cross(x, y):
  '''gives cross product of two vectors'''
  return [ x[1]*y[2] - x[2]*y[1],
           x[2]*y[0] - x[0]*y[2],
           x[0]*y[1] - x[1]*y[0] ]

def getDihedralUnited(all):
  '''list of 4 xyzs, gets the dihedral'''
  return getDihedral(all[0], all[1], all[2], all[3])

def getDihedral(a, b, c, d):
  '''4 xyzs, gets the dihedral'''
  cross1 = normalizeVector(cross( \
                 getNormalVector(a, b), \
                 getNormalVector(b, c)))
  cross2 = normalizeVector(cross( \
                 getNormalVector(b, c), \
                 getNormalVector(c, d)))
  try:
    dihedral1 = math.acos(dot(cross1, cross2))
  except ValueError:
    dihedral1 = 0.0 #sometimes the dot ends up a tiny bit above 1.0
  #have to figure out +- direction
  planeD = calculatePlaneD(cross1, b)
  planeFull = (cross1[0], cross1[1], cross1[2], planeD)
  if not checkPlaneSide(planeFull, d):
    dihedral1 = -dihedral1
  return dihedral1

def calculatePlaneD(normal, pointOnP):
  '''calculates the d of a plane where d = -ax -by -cz where normal = a,b,c
  and point on plane = x,y,z'''
  return -normal[0]*pointOnP[0]-normal[1]*pointOnP[1]-normal[2]*pointOnP[2]

def checkPlaneSide(plane, point):
  '''plane is normal + D (from function calculatePlaneD). sees if point is 
  in the direction of normal or not, return boolean'''
  sign = plane[0]*point[0]+plane[1]*point[1]+plane[2]*point[2]+plane[3]
  if sign >= 0:
    return True
  else:
    return False

def rotateAboutLine(aIn, dIn, xyz, theta):
  '''rotates the point xyz about the line d-a to an angle of theta radians'''
  #based on http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
  #first we have to constrain theta to be within -pi to +pi
  while theta < math.pi:
    theta += 2 * math.pi
  while theta > math.pi:
    theta -= 2 * math.pi
  da = getVector(dIn, aIn) #line through a and d
  #break down and just use the worst notation ever. someone punch me in the face
  a, b, c = aIn #unpack many things
  d, e, f = dIn
  u, v, w = da
  x, y, z = xyz
  #shortcuts
  uvw = length(da)
  uvw2 = uvw*uvw
  #long stupid equations
  newX =(a * (v**2. + w**2.) + u * (- b * v - c * w + u * x + v * y + w * z) + \
         (- a * (v**2. + w**2.) + u * (b * v + c * w - v * y - w * z) + \
            x * (v**2. + w**2.)) * math.cos(theta) + \
         (- c * v + b * w - w * y + v * z) * math.sin(theta) * uvw) / uvw2
  newY =(b * (u**2. + w**2.) + v * (- a * u - c * w + u * x + v * y + w * z) + \
         (- b * (u**2. + w**2.) + v * (a * u + c * w - u * x - w * z) + \
            y * (u**2. + w**2.)) * math.cos(theta) + \
         (c * u - a * w + w * x - u * z) * math.sin(theta) * uvw) / uvw2
  newZ =(c * (v**2. + u**2.) + w * (- a * u - b * v + u * x + v * y + w * z) + \
         (- c * (v**2. + u**2.) + w * (a * u + b * v - u * x - v * y) + \
            z * (v**2. + u**2.)) * math.cos(theta) + \
         (- b * u + a * v - v * x + u * y) * math.sin(theta) * uvw) / uvw2
  return newX, newY, newZ

def getTriNormalList(united):
  return getTriNormal(united[0], united[1], united[2])

def getTriNormal(a,b,c,firstTime=True):
  '''a,b and c are triange points in clockwise order, returns normal vector 
  that points out. returns NORMALIZED vector now. or 0s.'''
  #find a-b and c-b
  #vecAB = normalizeVector(getVector(a,b))
  #vecCB = normalizeVector(getVector(c,b))
  vecAB = getVector(a,b)
  vecCB = getVector(c,b)
  #does the cross product, that's all there is to it
  normal = cross(vecAB,vecCB)
  #only enter this part if all 0 and if first time being called
  if not firstTime: #has been called recursively. don't check 0s.don't normalize
    return normal
  elif firstTime and normal[0] == 0. and normal[1] == 0. and normal[2] == 0.:
    '''this is a big problem. attempt to call after permuting values'''
    newNor = getTriNormal(b,c,a,firstTime=False) #still maintains clockwise
    if newNor[0] == 0. and newNor[1] == 0. and newNor[2] == 0.:
      lastNo = getTriNormal(c,a,b,firstTime=False) #again
      #if this is zero we still have to return it
      if lastNo[0] == 0. and lastNo[1] == 0. and lastNo[2] == 0.:
        return lastNo #0s knowingly returned
      else:
        return normalizeVector(lastNo)
    else:
      return normalizeVector(newNor)
  else:
    return normalizeVector(normal)
  
def getAverage(listPoints):
  '''averages any number of 3d points passed in as list'''  
  average = [0.,0.,0.]
  for point in listPoints:
    for index in xrange(len(average)):
      average[index] += point[index]
  for index in xrange(len(average)):
    average[index] /= len(listPoints)
  return average

def getAverage1(listPoints):
  '''averages any number of 1d points passed in as list'''  
  average = 0.
  for point in listPoints:
    average += point
  average /= len(listPoints)
  return average

def getAverageArbitraryDimension(listPoints, dimension=2):
  '''averages any number of nD points passed in as list'''  
  average = [0. for count in xrange(dimension)]
  for point in listPoints:
    for index in xrange(len(average)):
      average[index] += point[index]
  for index in xrange(len(average)):
    average[index] /= len(listPoints)
  return average

def findMinsMaxsSpheres(spheres):
  '''goes through all spheres, finds the min and max in each dimension.
  spheres are expected in [x,y,z,r] format'''
  if 0 == len(spheres):
    return False,False #indicates failure
  mins, maxs = [], []
  for xyz in range(3):
    mins.append(spheres[0][xyz] - spheres[0][3]) #x - radius then y-rad, z-rad
    maxs.append(spheres[0][xyz] + spheres[0][3]) #x + radius then y+rad, z+rad
  for sphere in spheres[1:]: #already did the first
    for xyz in range(3):
      mins[xyz] = min(mins[xyz], sphere[xyz]-sphere[3]) 
      maxs[xyz] = max(maxs[xyz], sphere[xyz]+sphere[3]) 
  return mins, maxs

def planeDistToOrigin(normal):
  '''uses formula from http://mathworld.wolfram.com/Plane.html
  normal is a,b,c,d of plane
  dist = d / ((a^2 + b^2 + c^2) ^ (1/2))'''
  a,b,c,d = normal #unpack tuple for laziness
  return d / ((a**2. + b**2. + c**2.) ** 0.5)

def fixNormalZeros(vector):
  '''if all 0s, return unchanged, that's fine.
  if 1 or 2 0s, permute a tiny bit so there are no 0s. normalize and return'''
  alpha = 0.0000000000000000001
  if vector[0] == 0. and vector[1] == 0. and vector[2] == 0.:
    return vector #all zeros
  elif vector[0] == 0. or vector[1] == 0. or vector[2] == 0.:
    newVec = vector[:] #copy
    if vector[0] == 0.:
      newVec[0] += alpha
    if vector[1] == 0.:
      newVec[1] += alpha
    if vector[2] == 0.:
      newVec[2] += alpha
    return normalizeVector(newVec)
  else:
    return vector #no zeros

def withinTolerance(pointA, pointB, tolerance):
  '''trying to make something fast to check if pointA and pointB are within 
  the tolerance of each other. 
  exact distance function (l2, l1, linf) not a big deal'''
  if abs(pointA[0] - pointB[0]) < tolerance:
    if abs(pointA[1] - pointB[1]) < tolerance:
      if abs(pointA[2] - pointB[2]) < tolerance:
        return True
  return False

