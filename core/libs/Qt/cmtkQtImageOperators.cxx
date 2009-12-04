/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkQtImageOperators.h>

#include <math.h>

#include <qinputdialog.h>

#include <cmtkQtProgress.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QMenu*
QtImageOperators::CreatePopupMenu()
{
  QMenu* operatorsMenu = new QMenu;
  operatorsMenu->addAction( "&Median Filter...", this, SLOT( slotOperatorMedian() ) );
  operatorsMenu->addAction( "&Histogram Equalization...", this, SLOT( slotOperatorHistEq() ) );
  operatorsMenu->addAction( "&Sobel Edge Filter", this, SLOT( slotOperatorSobel() ) );
  operatorsMenu->addSeparator();

  QMenu* algOperatorsMenu = operatorsMenu->addMenu( "&Algebraic" );
  algOperatorsMenu->addAction( "&abs()", this, SLOT( slotOperatorAbs() ) );
  algOperatorsMenu->addAction( "&log()", this, SLOT( slotOperatorLog() ) );
  algOperatorsMenu->addAction( "&exp()", this, SLOT( slotOperatorExp() ) );

  return operatorsMenu;
}

void
QtImageOperators::slotOperatorsMedian()
{
  if ( this->StudyDataValid() ) 
    {
    bool ok;
    int radius = QInputDialog::getInt( this->Parent, "Median Filter", "Neighborhood radius:",  1, 1, 5, 1, &ok );
    if ( ok )
      {
      // user entered something and pressed OK
      if ( this->ProgressInstance )
	this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_DIALOG );
      (*(this->CurrentStudy))->GetVolume()->ApplyMedianFilter( radius );
      emit dataChanged( *(this->CurrentStudy) );
      } 
    else
      {
      // user pressed Cancel
      }
    }
}

void
QtImageOperators::slotOperatorsSobel()
{
  if ( this->StudyDataValid() ) 
    {
    if ( this->ProgressInstance )
      this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );
    (*(this->CurrentStudy))->GetVolume()->ApplySobelFilter();
    emit dataChanged( *(this->CurrentStudy) );
    }
}

void
QtImageOperators::slotOperatorsHistEq()
{
  if ( this->StudyDataValid() ) 
    {
    if ( this->ProgressInstance )
      this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );
    bool ok;
    int bins = QInputDialog::getInt( this->Parent, "Histogram Equalization", "Number of Histogram Bins:", 256, 2, 256, 1, &ok );
    if ( ok ) 
      {
      // user entered something and pressed OK
      if ( this->ProgressInstance )
	this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_DIALOG );
      (*(this->CurrentStudy))->GetVolume()->GetData()->HistogramEqualization( bins );
      emit dataChanged( *(this->CurrentStudy) );
      } 
    else
      {
      // user pressed Cancel
      }
    }
}

void
QtImageOperators::slotOperatorsAbs()
{
  if ( this->StudyDataValid() ) 
    {
    (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunction( fabs );
    emit dataChanged( *(this->CurrentStudy) );
    }
}

void
QtImageOperators::slotOperatorsLog()
{
  if ( this->StudyDataValid() ) 
    {
    (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunction( log );
    emit dataChanged( *(this->CurrentStudy) );
    }
}

void
QtImageOperators::slotOperatorsExp()
{
  if ( this->StudyDataValid() ) 
    {
    (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunction( exp );
    emit dataChanged( *(this->CurrentStudy) );
    }
}

} // namespace cmtk
