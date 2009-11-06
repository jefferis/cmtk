#!/usr/bin/python

import os

appRoot = '/Users/mphasak/Development/cmtk/cur/build/bin/'

convert = os.path.join(appRoot,"convert")
reformatx = os.path.join(appRoot,"reformatx")
mk_analyze_hdr = os.path.join(appRoot,"mk_analyze_hdr")

# The liver volume IDs
ids = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20']

# These two volumes didn't resample well so leave them out for now
badResamples = ['14','20']

# Factors for converting the Z thickness to 30
reformat30factors = [6.0999999999999996, 2.1333333333333333, 2.6333333333333333, 7.0666666666666664, 10.633333333333333, 3.7000000000000002, 8.3666666666666671, 7.5999999999999996, 7.0, 6.3666666666666663, 12.933333333333334, 7.333333333333333, 4.833333333333333, 4.2999999999999998, 13.133333333333333, 5.0333333333333332, 4.0333333333333332, 8.1666666666666661, 11.166666666666666, 6.0999999999999996]


def reformat30():
  inDir = "training"
  gridDir = "grid30"
  outDir = "resamp30"
  for i in range(0,len(ids)):
    id = ids[i]
    if id not in badResamples:
      divBy = str(reformat30factors[i])
      
      infile = os.path.join(inDir, "liver-orig0" + id + ".nhdr")
      gridfile = os.path.join(gridDir,"liver-orig0" + id + "-30grid.nhdr")
      outfile = os.path.join(outDir,"liver-orig0" + id + "-30.nhdr")
      
      cmd_mk_grid = mk_analyze_hdr + " -D 30,30,30 -V 17.06,17.06," + divBy +\
                            " " +  gridfile

      cmd_reformat = reformatx + " -o " + outfile +\
                        " --floating " + infile +\
                        " " + gridfile

      print cmd_mk_grid
      print cmd_reformat
      
      os.system(cmd_mk_grid)
      os.system(cmd_reformat)

reformat30()
