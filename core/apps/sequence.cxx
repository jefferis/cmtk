/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkValueSequence.h>

#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#include <math.h>
#include <stdio.h>

#include <cfloat>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <vector>

double MaxThreshold = 0;
bool UseMaxThreshold = false;
bool AbsoluteValues = false;

const char *OutputFormat = "%.6f";

const char *WriteHistogram = NULL;

double HistogramMin = 0;
bool HistogramMinSet = false;

double HistogramMax = 100;
bool HistogramMaxSet = false;

int HistogramBins = 1000;

int doMain(const int argc, const char *argv[]) {
  try {
    cmtk::CommandLine cl;
    cl.SetProgramInfo(cmtk::CommandLine::PRG_TITLE, "Value sequence");
    cl.SetProgramInfo(cmtk::CommandLine::PRG_DESCR,
                      "Analyze sequence of numerical values, which is read "
                      "from standard input");

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption(Key('t', "thresh"), &MaxThreshold,
                 "Maximum value threshold. All values above are ignored.",
                 &UseMaxThreshold);
    cl.AddSwitch(Key('a', "abs"), &AbsoluteValues, true,
                 "Use absolute values.");

    cl.BeginGroup("Histogram", "Histogram Options");
    cl.AddOption(Key("histogram-min"), &HistogramMin,
                 "Minimum of the histogram value range. All values below this "
                 "will be counted in the first histogram bin.",
                 &HistogramMinSet);
    cl.AddOption(Key("histogram-max"), &HistogramMax,
                 "Maximum of the histogram value range. All values above this "
                 "will be counted in the last histogram bin.",
                 &HistogramMaxSet);
    cl.AddOption(Key("histogram-bins"), &HistogramBins,
                 "Number of histogram bins.");
    cl.EndGroup();

    cl.BeginGroup("Output", "Output Options");
    cl.AddOption(Key('f', "format"), &OutputFormat,
                 "Output number format in printf() style.");
    cl.AddOption(
        Key("write-histogram"), &WriteHistogram,
        "Path for optional histogram output in comma-separated (CSV) format.");
    cl.EndGroup();

    cl.Parse(argc, argv);
  } catch (const cmtk::CommandLine::Exception &e) {
    cmtk::StdErr << "Command line parse error: " << e << "\n";
    return 1;
  }

  cmtk::ValueSequence<double> seq;
  std::list<double> list;

  unsigned int countOverThreshold = 0;
  double f;
  while (!std::cin.eof()) {
    std::cin >> f;

    if (!finite(f)) break;

    if (AbsoluteValues) f = fabs(f);

    if (UseMaxThreshold && (f > MaxThreshold))
      ++countOverThreshold;
    else {
      seq.Proceed(f);
      list.push_back(f);
    }

    f = std::numeric_limits<double>::signaling_NaN();
  }

  // Check for empty input and simply exit in that case
  if (list.empty()) {
    return 0;
  }

  // Print value counts
  const size_t totalNumberOfValues = seq.GetNValues() + countOverThreshold;
  printf("Number of Values:\t%d\n", (int)totalNumberOfValues);
  printf("Values over Threshold:\t%u (%.2f%%)\n", countOverThreshold,
         100.0 * countOverThreshold / totalNumberOfValues);

  // Format and print simple statistics
  char format[120];
  snprintf(format, sizeof(format),
           "\nSTAT\tMin\tMax\tMean\tStdDev\nSTATval\t%s\t%s\t%s\t%s\n",
           OutputFormat, OutputFormat, OutputFormat, OutputFormat);
  printf(format, seq.GetMinimum(), seq.GetMaximum(), seq.GetAverage(),
         sqrt(seq.GetVariance()));

  // Create sorted list for percentile computation
  std::vector<double> sorted(list.begin(), list.end());
  std::sort(sorted.begin(), sorted.end());

  printf("\nPERC");
  const int percentiles[] = {5, 10, 25, 50, 75, 90, 95, -1};
  for (size_t idx = 0; percentiles[idx] > 0; ++idx) {
    printf("\t%d", percentiles[idx]);
  }

  snprintf(format, sizeof(format), "\t%s", OutputFormat);
  printf("\nPERCval");
  for (size_t idx = 0; percentiles[idx] > 0; ++idx) {
    if (percentiles[idx] == 50) {
      const size_t medianIdx = sorted.size() / 2;
      if (sorted.size() & 1) {
        printf(format, sorted[medianIdx]);
      } else {
        printf(format, 0.5 * (sorted[medianIdx] + sorted[1 + medianIdx]));
      }
    } else {
      printf(format, sorted[(size_t)(sorted.size() * 0.01 * percentiles[idx])]);
    }
  }
  printf("\n");

  // Create histogram is desired
  if (WriteHistogram) {
    // If no histogram minimum given, use minimum data value (first in sorted
    // vector)
    if (!HistogramMinSet) HistogramMin = sorted[0];

    // If no histogram maximum given, use maximum data value (last in sorted
    // vector)
    if (!HistogramMaxSet) HistogramMax = sorted[sorted.size() - 1];

    std::ofstream hstream(WriteHistogram);
    if (!hstream.good()) {
      cmtk::StdErr << "ERROR: could not open file " << WriteHistogram
                   << " for writing the histogram.\n";
      throw cmtk::ExitException(1);
    }

    hstream << "bin_min,bin_max,bin_count\n";

    const double binWidth = (HistogramMax - HistogramMin) / HistogramBins;
    std::vector<double>::const_iterator sample = sorted.begin();
    for (int bin = 0; bin < HistogramBins; ++bin) {
      const double binFrom = HistogramMin + bin * binWidth;
      const double binTo = HistogramMin + (1 + bin) * binWidth;

      // go through sorted vector and add to this bin whatever belongs
      int countBin = 0;
      for (; (sample != sorted.end()) && (*sample < binTo);
           ++countBin, ++sample) {
      }

      // last bin -- add whatever is left
      if (bin == HistogramBins - 1)
        for (; (sample != sorted.end()); ++countBin, ++sample) {
        }

      hstream << binFrom << "," << binTo << "," << countBin << "\n";
    }
  }

  return 0;
}

#include "cmtkSafeMain"
