#!/usr/bin/python

import os

appRoot = '/Users/mphasak/Development/cmtk/cur/build/bin'
sliverRoot = '/Users/mphasak/Development/sliver07'

convert = os.path.join(appRoot,"convert")
reformatx = os.path.join(appRoot,"reformatx")
mk_analyze_hdr = os.path.join(appRoot,"mk_analyze_hdr")

# The liver volume IDs
ids = ['01','02','03','04','05','06','07','08','09','10',\
       '11','12','13','14','15','16','17','18','19','20']

goodIDs = ['01','02','03','04','05','06','07','08','09','10',\
       '11','12','13','15','16','17','18','19']

# These two volumes didn't resample well so leave them out for now
badResamples = ['14','20']

# number of Z voxels in each volume
zdims = [183,64,79,212,319,111,251,228,210,191,\
         388,220,145,129,394,151,121,245,335,183]

def reformat(toGridSize, toVoxelSize):
  inDir = os.path.join(sliverRoot, "training")
  gridDir = os.path.join(sliverRoot, "grid" + str(toGridSize))
  #outDir = os.path.join(sliverRoot, "resamp" + str(toGridSize))
  outDir = os.path.join(sliverRoot, "resamp")
 
  if not os.path.exists(gridDir):
    os.makedirs(gridDir)
  if not os.path.exists(gridDir):
    os.makedirs(outDir)

  for i in range(0,len(goodIDs)):
    id = goodIDs[i]
    if id not in badResamples:
      
      inFile = os.path.join(inDir, "liver-orig0%s.nhdr" % id)
      gridFile = os.path.join(gridDir,"liver-orig0%s-%dgrid.nhdr" % (id,toGridSize))
      outFile = os.path.join(outDir,"liver-orig0%s-%d.nhdr" % (id,toGridSize))
     
      #divXBy = divYBy = 512 / float(toGridSize)
      #divZBy = zdims[i] / float(toGridSize) 
     
      # % here is like sprintf in C
      newGridSize = "%d,%d,%d" % (toGridSize, toGridSize, toGridSize)  # voxels
      newVolumeSizeMM = "%f,%f,%f" % (toVoxelSize[i][0],toVoxelSize[i][1],toVoxelSize[i][2]) # mm

      cmd_mk_grid = ( mk_analyze_hdr + 
                      " -D " + newGridSize +
                      " -V " + newVolumeSizeMM + 
                      " " +  gridFile )

      cmd_reformat = ( reformatx + " -o " + outFile +
                       " --floating " + inFile +
                       " " + gridFile )

      cmd_unzip = "gunzip -f " + outDir + "/*gz"

      
      if not os.path.exists(gridFile):
        print cmd_mk_grid
        os.system(cmd_mk_grid)
        print cmd_reformat
        os.system(cmd_reformat)
        print cmd_unzip
        os.system(cmd_unzip)

#reformat(30)
