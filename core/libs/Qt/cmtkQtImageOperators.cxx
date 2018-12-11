/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkQtImageOperators.h"

#include <math.h>

#include <qinputdialog.h>

#include <Qt/cmtkQtProgress.h>

#include <Base/cmtkTypedArrayFunctionHistogramEqualization.h>
#include <Base/cmtkDataGridFilter.h>
#include <Base/cmtkMathFunctionWrappers.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QMenu*
QtImageOperators::CreatePopupMenu()
{
  QMenu* operatorsMenu = new QMenu;
  operatorsMenu->setTitle( "&Operators" );
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
QtImageOperators::slotOperatorMedian()
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

      (*(this->CurrentStudy))->GetVolume()->SetData( DataGridFilter( (*(this->CurrentStudy))->GetVolume() ).GetDataMedianFiltered( radius ) );
      
      emit dataChanged( *(this->CurrentStudy) );
      } 
    else
      {
      // user pressed Cancel
      }
    }
}

void
QtImageOperators::slotOperatorSobel()
{
  if ( this->StudyDataValid() ) 
    {
    if ( this->ProgressInstance )
      this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );

    (*(this->CurrentStudy))->GetVolume()->SetData( DataGridFilter( (*(this->CurrentStudy))->GetVolume() ).GetDataSobelFiltered() );

    emit dataChanged( *(this->CurrentStudy) );
    }
}

void
QtImageOperators::slotOperatorHistEq()
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
      (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunctionObject( TypedArrayFunctionHistogramEqualization( (*(*(this->CurrentStudy))->GetVolume()->GetData()), bins ) );
      emit dataChanged( *(this->CurrentStudy) );
      } 
    else
      {
      // user pressed Cancel
      }
    }
}

void
QtImageOperators::slotOperatorAbs()
{
  if ( this->StudyDataValid() ) 
    {
    (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Abs );
    emit dataChanged( *(this->CurrentStudy) );
    }
}

void
QtImageOperators::slotOperatorLog()
{
  if ( this->StudyDataValid() ) 
    {
    (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Log );
    emit dataChanged( *(this->CurrentStudy) );
    }
}

void
QtImageOperators::slotOperatorExp()
{
  if ( this->StudyDataValid() ) 
    {
    (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunctionDouble( cmtk::Wrappers::Exp );
    emit dataChanged( *(this->CurrentStudy) );
    }
}

} // namespace cmtk
