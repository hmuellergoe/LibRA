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

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/IO/ArrayIO.h>
#include <casacore/casa/BasicMath/Math.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Logging/LogIO.h>

#include <synthesis/MeasurementEquations/MatrixCleaner.h>
#include <synthesis/TransformMachines2/Utils.h>

#include<synthesis/MeasurementEquations/WaveletMSCleaner.h>

using namespace casacore;
using namespace std::chrono;
namespace casa { //# NAMESPACE CASA - BEGIN
WaveletMSCleaner::WaveletMSCleaner():
  MatrixCleaner()
{}

WaveletMSCleaner::~WaveletMSCleaner()
{
  destroyScales();
  if(!itsMask.null()) itsMask=0;
}

// Make a single scale size image
void WaveletMSCleaner::makeScale(Matrix<Float>& iscale, const Float& scaleSize) 
{
  
  Int nx=iscale.shape()(0);
  Int ny=iscale.shape()(1);
  //Matrix<Float> iscale(nx, ny);
  iscale=0.0;
  
  Double refi=nx/2;
  Double refj=ny/2;
  
  if(scaleSize==0.0) {
    iscale(Int(refi), Int(refj)) = 1.0;
  }
  else {
    AlwaysAssert(scaleSize>0.0,AipsError);

    Int mini = max( 0, (Int)(refi-scaleSize));
    Int maxi = min(nx-1, (Int)(refi+scaleSize));
    Int minj = max( 0, (Int)(refj-scaleSize));
    Int maxj = min(ny-1, (Int)(refj+scaleSize));

    Float ypart=0.0;
    Float volume=0.0;
    Float rad2=0.0;
    Float rad=0.0;

    for (Int j=minj;j<=maxj;j++) {
      ypart = square( (refj - (Double)(j)) / scaleSize );
      for (Int i=mini;i<=maxi;i++) {
		rad2 =  ypart + square( (refi - (Double)(i)) / scaleSize );
		if (rad2 < 1.0) {
			if (rad2 <= 0.0) {
				rad = 0.0;
			} else {
				rad = sqrt(rad2);
			}
			iscale(i,j) = dogwavelet(rad, scaleSize);
			volume += iscale(i,j);
		} else {
			iscale(i,j) = 0.0;
		}
      }
    }
    iscale/=volume;
  }
}

// Calculate the spheroidal function
Float WaveletMSCleaner::wavelet(Float rad) {
  return 0.0;
}

Float WaveletMSCleaner::dogwavelet(Float rad, Float scaleSize) {
	
  if (rad <= 0) {
    return 1.0;
  } else if (rad >= 5.0) {
    return 0.0;
  }
  
    Float sigma_1 = scaleSize;
    Float sigma_2 = 1.6 * scaleSize; //Research showed that 1.6 was the ideal difference between STD values
    
    Float gaussian_1 = (1/(2*sigma_1*sigma_1*3.141)) * (exp(-0.5 * pow(rad/sigma_1 , 2.0)));
    Float gaussian_2 = (1/(2*sigma_2*sigma_2*3.141)) * (exp(-0.5 * pow(rad/sigma_2 , 2.0)));
    
    Float DoG_wavelet = gaussian_1 - gaussian_2;
  
    return DoG_wavelet;
}

} //# NAMESPACE CASA - END
