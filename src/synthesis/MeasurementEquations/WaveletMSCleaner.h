//# Copyright (C) 1997-2010
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id:  $

#ifndef SYNTHESIS_WAVELETMSCLEANER_H
#define SYNTHESIS_WAVELETMSCLEANER_H

//# Includes
#include <synthesis/MeasurementEquations/MatrixCleaner.h>
#include <deque>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/casa/Utilities/CountedPtr.h>

namespace casa { //# NAMESPACE CASA - BEGIN

class WaveletMSCleaner : public MatrixCleaner
{
public:
  // Create a cleaner : default constructor
  WaveletMSCleaner();

  // The destructor does nothing special.
  ~WaveletMSCleaner();
  
  void setWaveletControl(const casacore::Float waveletdummyparam) {itsWaveletDummyParam=waveletdummyparam;}
  
protected:

  casacore::Float itsWaveletDummyParam;
  
  void makeScale(casacore::Matrix<casacore::Float>& scale, const casacore::Float& scaleSize) override;
  float wavelet(float rad);

};

} //# NAMESPACE CASA - END

#endif
