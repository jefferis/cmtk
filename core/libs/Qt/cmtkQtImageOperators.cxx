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
#include <Q3PopupMenu>

#include <cmtkQtProgress.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

Q3PopupMenu*
QtImageOperators::CreatePopupMenu()
{
  Q3PopupMenu* algOperatorsMenu = new Q3PopupMenu;
  algOperatorsMenu->insertItem( "&abs()", OPERATORS_MENU_ABS );
  algOperatorsMenu->insertItem( "&log()", OPERATORS_MENU_LOG );
  algOperatorsMenu->insertItem( "&exp()", OPERATORS_MENU_EXP );
  QObject::connect( algOperatorsMenu, SIGNAL( activated( int ) ),
		    this, SLOT( slotOperatorsMenu( int ) ) );

  Q3PopupMenu* operatorsMenu = new Q3PopupMenu;
  operatorsMenu->insertItem( "&Median Filter...", OPERATORS_MENU_MEDIAN );
  operatorsMenu->insertItem( "&Histogram Equalization...", OPERATORS_MENU_HISTOGRAM );
  operatorsMenu->insertItem( "&Sobel Edge Filter", OPERATORS_MENU_SOBEL );
  operatorsMenu->insertSeparator();
  operatorsMenu->insertItem( "&Algebraic", algOperatorsMenu );
  QObject::connect( operatorsMenu, SIGNAL( activated( int ) ), this, SLOT( slotOperatorsMenu( int ) ) );

  return operatorsMenu;
}

void
QtImageOperators::slotOperatorsMenu( int command )
{
  switch ( command )
    {
    case OPERATORS_MENU_MEDIAN: 
    {
    if ( this->StudyDataValid() ) 
      {
      bool ok;
      int radius = QInputDialog::getInteger( "Median Filter", "Neighborhood radius:",  1, 1, 5, 1, &ok, this->Parent );
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
    break;
    }
    case OPERATORS_MENU_SOBEL: 
    {
    if ( this->StudyDataValid() ) 
      {
      if ( this->ProgressInstance )
	this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );
      (*(this->CurrentStudy))->GetVolume()->ApplySobelFilter();
      emit dataChanged( *(this->CurrentStudy) );
      }
    break;
    }
    case OPERATORS_MENU_HISTOGRAM: 
    {
    if ( this->StudyDataValid() ) 
      {
      if ( this->ProgressInstance )
	this->ProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );
      bool ok;
      int bins = QInputDialog::getInteger( "Histogram Equalization", "Number of Histogram Bins:", 256, 2, 256, 1, &ok, this->Parent );
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
    break;
    case OPERATORS_MENU_ABS: 
    case OPERATORS_MENU_LOG: 
    case OPERATORS_MENU_EXP: 
    {
    if ( this->StudyDataValid() ) 
      {
      switch ( command ) 
	{
	case OPERATORS_MENU_ABS: 
	  (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunction( fabs );
	  emit dataChanged( *(this->CurrentStudy) );
	  break;
	case OPERATORS_MENU_LOG: 
	  (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunction( log );
	  emit dataChanged( *(this->CurrentStudy) );
	  break;
	case OPERATORS_MENU_EXP: 
	  (*(this->CurrentStudy))->GetVolume()->GetData()->ApplyFunction( exp );
	  emit dataChanged( *(this->CurrentStudy) );
	  break;
	}
      }
    break;
    }
    }
    }
}

} // namespace cmtk
