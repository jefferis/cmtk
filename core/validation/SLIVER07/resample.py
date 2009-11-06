#!/usr/bin/python

import os

appRoot = '/Users/mphasak/Development/cmtk/cur/build/bin'
sliverRoot = '/Users/mphasak/Development/sliver07'

convert = os.path.join(appRoot,"convert")
reformatx = os.path.join(appRoot,"reformatx")
mk_analyze_hdr = os.path.join(appRoot,"mk_analyze_hdr")

# The liver volume IDs
ids = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20']

# These two volumes didn't resample well so leave them out for now
badResamples = ['14','20']

# number of Z voxels in each volume
zdims = [183,64,79,212,319,111,251,228,210,191,388,220,145,129,394,151,121,245,335,183]

def reformat(to):
  inDir = "training"
  gridDir = "grid" + str(to)
  outDir = "resamp" + str(to)
  
  for i in range(0,len(ids)):
    id = ids[i]
    if id not in badResamples:
      
      infile = os.path.join(sliverRoot, inDir, "liver-orig0" + id + ".nhdr")
      gridfile = os.path.join(sliverRoot, gridDir,"liver-orig0" + id + "-30grid.nhdr")
      outfile = os.path.join(sliverRoot, outDir,"liver-orig0" + id + "-30.nhdr")
     
      divXBy = divYBy = 512 / float(to)
      divZBy = zdims[i] / float(to) 
     
      # % here is like sprintf in C
      newGridSize = "%f,%f,%f" % (to, to, to)  # voxels
      newVolumeSizeMM = "%f,%f,%f" % (divXBy, divYBy, divZBy) # mm

      cmd_mk_grid = ( mk_analyze_hdr + 
                      " -D " + newGridSize +
                      " -V " + newVolumeSizeMM + 
                      " " +  gridfile )

      cmd_reformat = ( reformatx + " -o " + outfile +
                       " --floating " + infile +
                       " " + gridfile )

      print cmd_mk_grid
      print cmd_reformat
      
      #os.system(cmd_mk_grid)
      #os.system(cmd_reformat)

#reformat(30)
