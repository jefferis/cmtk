#!/usr/bin/python
import io
import os
import struct
import math
import scipy
import scipy.linalg
from scipy import matrix
import numpy
import resample


sliverRoot = '/Users/mphasak/Development/sliver07'
testDir = os.path.join(sliverRoot, 'resamp')
rawDir = os.path.join(sliverRoot, 'resamp')
comDir = os.path.join(sliverRoot, 'com')
dimsDir = os.path.join(sliverRoot, 'dims')
cropDir = os.path.join(sliverRoot, 'cropfiles')

volumeIDs = ['01','02','03','04','05','06','07','08','09','10',\
             '11','12','13','14','15','16','17','18','19','20']
goodIDs = ['01','02','03','04','05','06','07','08','09','10',\
           '11','12','13','15','16','17','18','19']

# Description of the image data format.
# This will need to be determined on the fly
# (or ideally rendered moot by the use of a
# real NRRD-reader library), but for now
# I'm hard-coding little-endian short (which
# is two bytes per value)
fmt = '<h'      # < indicates little-endian, h indicates short
ctypelength = 2 # length of the ctype encoding the image data

def getCentersOfMass():
  centersOfMass = []
  fileList = os.listdir(comDir)
  for filename in fileList:
    if filename.endswith(".loc"):
      comFile = os.path.join(comDir, filename)
      with open(comFile) as f:
        for line in f:
          line = line.rstrip('\n')
          items = line.split(" ")
          items = map(float, items)
          centersOfMass.append(items)
#      print filename
  return centersOfMass

def getVolumeDimsInVoxels():
  dims = []
  fileList = os.listdir(dimsDir)
  for filename in fileList:
    if filename.endswith(".voxels"):
      dimFile = os.path.join(dimsDir, filename)
      with open(dimFile) as f:
        for line in f:
          line = line.rstrip('\n')
          items = line.split(" ")
          items = map(float, items)
          dims.append(items)
#      print filename
  return dims

def getVolumeDimsInMM():
  dims = []
  fileList = os.listdir(dimsDir)
  for filename in fileList:
    if filename.endswith(".volsize"):
      dimFile = os.path.join(dimsDir, filename)
      with open(dimFile) as f:
        for line in f:
          line = line.rstrip('\n')
          items = line.split(" ")
          items = map(float, items)
          dims.append(items)
#      print filename
  return dims

def getVoxelSizesInMM():
  dims = []
  fileList = os.listdir(dimsDir)
  for filename in fileList:
    if filename.endswith(".voxelsize"):
      dimFile = os.path.join(dimsDir, filename)
      with open(dimFile) as f:
        for line in f:
          line = line.rstrip('\n')
          items = line.split(" ")
          items = map(float, items)
          dims.append(items)
#      print filename
  return dims

def normalizeCOMs(coms, dims):
  assert len(coms) == len(dims)
  normalizedCOMs = []
  for i in range(0, len(coms)):
    assert len(coms[i]) == len(dims[i])
    normalizedCOM = []
    for dim in range(0, len(coms[i])):
      normalizedCOM.append(coms[i][dim]/dims[i][dim])
    normalizedCOMs.append(normalizedCOM)
  return normalizedCOMs  

def scaleVoxelDims(sourceDimsMM, targetCubeDim):
  scaledVoxelDims = []
  for i in range(0, len(sourceDimsMM)):
    scaledVoxelDim = []
    for dim in range(0, len(sourceDimsMM[i])):
      scaledVoxelDim.append(sourceDimsMM[i][dim] / targetCubeDim)
    scaledVoxelDims.append(scaledVoxelDim)
  return scaledVoxelDims  

def getTrainingIntensityVectors(resampleTo):
  """
  This function may be hard to read. It's parsing the binary
  data of a NRRD volume, and relies on some unusual string
  concatenation syntax and the 'struct' class, which has
  an entirely different meaning than 'struct' in C: here
  it's a class providing the 'unpack' function for turning
  a string representation of a ctype (short, float, etc.) 
  into a python float.
   (This could be cleaned us by finding and using an actual
  proper NRRD library.)
  """
  intensityVectors = []
  fileList = os.listdir(rawDir)
  for filename in fileList:
    if filename.endswith("-%d.raw" % resampleTo):
      rawFile = os.path.join(rawDir, filename)
      rawFileHandle = io.FileIO(rawFile, 'r') # FileIO opens in binary mode
      bytes = rawFileHandle.read() # a list of uninterpreted bytes
      # Unpack the bytes:
      intensities = []
      for i in range(0, len(bytes)):
        if ((i % ctypelength)==0):
          # Current raw item is the concatenation of bytes from
          # list index i to index i+ctypelength.
          currentDatum = ''.join(byte for byte in bytes[i:i+ctypelength])
          # "Unpack" raw bytes to a numeric value:
          intensities.append(struct.unpack(fmt, currentDatum)[0])
      intensityVectors.append(intensities)
#      print filename
  return intensityVectors

def getTestingIntensityVectors():
  intensityVectors = []
  fileList = os.listdir(testDir)
  for filename in fileList:
    if filename.endswith(".raw"):
      rawFile = os.path.join(testDir, filename)
      rawFileHandle = io.FileIO(rawFile, 'r') # FileIO opens in binary mode
      bytes = rawFileHandle.read() # a list of uninterpreted bytes
      # Unpack the bytes:
      intensities = []
      for i in range(0, len(bytes)):
        if ((i % ctypelength)==0):
          # Current raw item is the concatenation of bytes from
          # list index i to index i+ctypelength.
          currentDatum = ''.join(byte for byte in bytes[i:i+ctypelength])
          # "Unpack" raw bytes to a numeric value:
          intensities.append(struct.unpack(fmt, currentDatum)[0])
      intensityVectors.append(intensities)
#      print filename
  return intensityVectors


def pseudoInverse(A):
  return numpy.linalg.pinv(A)

def computeProposedCOM(A, B, x, row): 
  selectedVector = A[row]
  proposedCOM = numpy.dot(selectedVector, x)
  return proposedCOM
 
def computeResidualZ(A, B, x, ommittedRow):
  testVector = A[ommittedRow]
  actualCOM = B[ommittedRow]
  proposedCOM = numpy.dot(testVector, x)
  r = actualCOM[2] - proposedCOM[2]
  #pz = proposedCOM[2]
  #az = actualCOM[2]
  #zd = pz - az
  #print pz, az, zd
  return r

def computeResidual(A, B, x, ommittedRow):
  testVector = A[ommittedRow]
  actualCOM = B[ommittedRow]
  proposedCOM = numpy.dot(testVector, x)
  residualVector = actualCOM - proposedCOM
  squaredResiduals = map(lambda a: a*a, residualVector)
  r = math.sqrt(sum(squaredResiduals))
  return r

def leaveOneOut(A, B):
  assert len(A) == len(B)
  rowRange = range(0,len(A))
  errs = []
  for leaveOut in rowRange:
    currentSubA = []
    currentSubB = []
    for i in rowRange:
      if i != leaveOut:
        currentSubA.append(A[i])
        currentSubB.append(B[i])
    x = numpy.dot(pseudoInverse(currentSubA), currentSubB)
#    computeProposedCOM(A, B, x, leaveOut)
#    r = computeResidual(A, B, x, leaveOut)
    r = computeResidualZ(A, B, x, leaveOut)
    errs.append(r)
  return errs  

rf = 5
coms = getCentersOfMass()
dimsMM = getVolumeDimsInMM()
dimsVox = getVolumeDimsInVoxels()
voxelSizes = getVoxelSizesInMM()

def computeCOMs(rf):
  trainA = getTrainingIntensityVectors(rf)
  testA = getTrainingIntensityVectors(rf)
  b = normalizeCOMs(coms, dimsMM)
  x = numpy.dot(pseudoInverse(trainA), b)
  newComs = []
  for i in range(0,len(testA)):
    newComs.append(computeProposedCOM(testA,b,x,i))
  return newComs

def mmPtToVoxel(mmPt, imageIdx):
  vox = []
  for i in range(0,3):
    vox.append(int(mmPt[i] * dimsVox[imageIdx][i]))
  return vox

def computeBoundingBoxes(rf):
  newComs = computeCOMs(rf)
  boxheightPx = 210
  boxes = []
  for i in range(0,len(goodIDs)):
    boxCenterVox = mmPtToVoxel(newComs[i],i)
    centerZ = boxCenterVox[2]
    topZ = min(512, centerZ + int(boxheightPx/2))
    botZ = max(0, centerZ - int(boxheightPx/2))
    boxes.append([0,0,botZ,512,512,topZ])
  return boxes  

def saveBoundingBoxes(rf):
  boxes = computeBoundingBoxes(rf)
  
  if not os.path.exists(cropDir):
    os.makedirs(cropDir)

  for i in range(0, len(goodIDs)):
    filename = "liver-orig0%s.crop" % goodIDs[i]
    cropFile = os.path.join(cropDir, filename)
    with open(cropFile, 'w') as f:
      for p in range(0,len(boxes[i])):
        if (p > 0): f.write(",") 
        f.write(str(boxes[i][p]))
      f.close()


#scaledVoxelSizes_coll = {}
#A_coll = {}
#b_coll = {}
#x_coll = {}
#errs_coll = {}
#avgerrs = []
#stddevs = []
#resampleFactors = [30, 25, 22, 21, 20, 19, 18, 16, 15, 14, 12, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
#for rf in resampleFactors:
#  scaledVoxelSizes_coll[rf] = scaleVoxelDims(dimsMM, rf)
#  resample.reformat(rf, scaledVoxelSizes) 
#  A_coll[rf] = getTrainingIntensityVectors(rf)
#  b_coll[rf] = normalizeCOMs(coms, dimsMM)
#  x_coll[rf] = numpy.dot(pseudoInverse(A_coll[rf]), b_coll[rf])
#  errs_coll[rf] = leaveOneOut(A_coll[rf], b_coll[rf])
#  avgerr = numpy.average(errs_coll[rf])
#  stddev = numpy.std(errs_coll[rf])
#  avgerrs.append(avgerr)
#  stddevs.append(stddev)
#  print "%d,%f,%f" % (rf, avgerr, stddev)
