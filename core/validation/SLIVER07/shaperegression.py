#!/usr/bin/python
import io
import os
import struct
import math
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
    if filename.endswith("-30.raw"):
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
      print filename
  return intensityVectors


def pseudoInverse(A):
  return numpy.linalg.pinv(A)

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
  for leaveOut in rowRange:
    currentSubA = []
    currentSubB = []
    for i in rowRange:
      if i != leaveOut:
        currentSubA.append(A[i])
        currentSubB.append(B[i])
    x = numpy.dot(pseudoInverse(currentSubA), currentSubB)
    r = computeResidual(A, B, x, leaveOut)
    print r
  

