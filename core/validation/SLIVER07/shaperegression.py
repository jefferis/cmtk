
import io
import os
import struct
import scipy
import scipy.linalg
from scipy import matrix
import numpy

pathRoot = '/Users/mphasak/Development/sliver07'
rawDir = os.path.join(pathRoot, 'resamp30')
comDir = os.path.join(pathRoot, 'com')

volumeIDs = ['01','02','03','04','05','06','07','08','09','10',\
             '11','12','13','14','15','16','17','18','19','20']
resampleFactor = "8"

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
      print filename
  return centersOfMass

def getIntensityVectors():
  intensityVectors = []
  fileList = os.listdir(rawDir)
  for filename in fileList:
    if filename.endswith("-30.raw"):
      rawFile = os.path.join(rawDir, filename)
      rawFileHandle = io.FileIO(rawFile, 'r') # FileIO opens in binary mode
      bytes = rawFileHandle.read() # a list of uninterpreted bytes
      # Unpack the bytes:
      intensities = []
      for i in range(0, len(bytes)):
        if ((i%ctypelength)==0):
          # Current raw item is the concatenation of bytes from
          # list index i to index i+ctypelength:
          currentDatum = reduce((lambda a,b: a+b),bytes[i:i+ctypelength])
          # "Unpack" raw bytes to a numeric value:
          intensities.append(struct.unpack(fmt, currentDatum)[0])
      intensityVectors.append(intensities)
      print filename
  return intensityVectors


def pseudoInverse(A):
  return numpy.linalg.pinv(A)

def computeResidual(A, B, x, ommittedRow):
  testVector = A[ommittedRow]
  actualCOM = B[ommittedRow]
  proposedCOM = testVector * x
  residualVector = actualCOM - proposedCOM


def leaveOneOut(A, B):
  assert len(A) == len(B)
  rowRange = range(0,len(A))
  for leaveOut in rowRange:
    currentSubA = []
    currentSubB = []
    for i in rowRange:
      if i != leaveOut:
        currentSubA.append(A[i])
        currentSubB.append(B[i])
    print currentSubA
    currentSubA = matrix(currentSubA)
    currentSubB = matrix(currentSubB)
    print currentSubA.shape
    print pseudoInverse(currentSubA).shape
    print currentSubB.shape
    #x = pseudoInverse(currentSubA) * currentSubB
    print x
  

