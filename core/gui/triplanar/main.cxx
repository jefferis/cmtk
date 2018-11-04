/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#include <System/cmtkStackBacktrace.h>
namespace cmtk {
static StackBacktrace StackBacktraceInstance;
}

#include <Qt/cmtkQtProgress.h>
#include <Qt/cmtkQtTriplanarViewer.h>

#include <IO/cmtkStudy.h>
#include <System/cmtkCommandLine.h>

#include <qapplication.h>

#include <string.h>
#include <list>
#include <string>

int doMain(int argc, char *argv[]) {
  QApplication app(argc, argv);

  cmtk::QtTriplanarViewer *viewer = new cmtk::QtTriplanarViewer;
  if (viewer) {
    cmtk::QtProgress progressInstance(viewer);
    progressInstance.SetProgressWidgetMode(cmtk::QtProgress::PROGRESS_DIALOG);
    cmtk::Progress::SetProgressInstance(&progressInstance);

    if ((argc > 1)) {
      if (!strcmp(argv[1], "--exec")) {
        viewer->hide();
        return viewer->ExecuteBatchMode(argc - 2, argv + 2);
      } else if (!strcmp(argv[1], "--xml")) {
        exit(1);
      }
    }

    for (int i = 1; i < argc; ++i) {
      viewer->slotAddStudy(argv[i]);
    }
    viewer->show();

    return app.exec();
  }

  return 0;
}

int main(const int argc, char *argv[]) {
#ifdef _MSC_VER
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  cmtk::Threads::CheckEnvironment();  // need this to check for
                                      // "CMTK_NUM_THREADS" and constrain OpenMP
                                      // accordingly

#ifdef CMTK_BUILD_STACKTRACE
  cmtk::StackBacktrace::Static();
#endif

  int exitCode = 0;
  try {
    exitCode = doMain(argc, argv);
  } catch (const cmtk::ExitException &ex) {
    exitCode = ex.ExitCode();
  }
  return exitCode;
}
