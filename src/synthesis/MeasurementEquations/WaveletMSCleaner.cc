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
void WaveletMSCleaner::makeScaleDiff(Matrix<Float>& iscale, const Float& scaleSize1, const Float& scaleSize2) 
{
  
  Int nx=iscale.shape()(0);
  Int ny=iscale.shape()(1);
  //Matrix<Float> iscale(nx, ny);
  iscale=0.0;
  
  Double refi=nx/2;
  Double refj=ny/2;
   
  if(scaleSize1==0.0) {
	iscale(Int(refi), Int(refj)) = 1.0;
  }
  else {
	AlwaysAssert(scaleSize1>0.0,AipsError);
	AlwaysAssert(scaleSize2>0.0,AipsError);
	AlwaysAssert(scaleSize2>=scaleSize1,AipsError);

	// Note that the larger scale always needs to be ScaleSize2
	Int mini = max( 0, (Int)(refi-5*scaleSize2));
	Int maxi = min(nx-1, (Int)(refi+5*scaleSize2));
	Int minj = max( 0, (Int)(refj-5*scaleSize2));
	Int maxj = min(ny-1, (Int)(refj+5*scaleSize2));

	Float ypart=0.0;
	Float volume=0.0;
	Float rad2=0.0;
	Float rad=0.0;

	for (Int j=minj;j<=maxj;j++) {
	  ypart = square( (refj - (Double)(j)) );
	  for (Int i=mini;i<=maxi;i++) {
		rad2 =  ypart + square( (refi - (Double)(i)) );
		if (rad2 < 1.0) {
			if (rad2 <= 0.0) {
				rad = 0.0;
			} else {
				rad = sqrt(rad2);
			}
			if (scaleSize1 == scaleSize2)
				iscale(i,j) = gaussian(rad, scaleSize1); //last Gaussian component
			else
				iscale(i,j) = dogwavelet(rad, scaleSize1, scaleSize2);
			//volume += iscale(i,j);
		} else {
			iscale(i,j) = 0.0;
		}
	  }
	}
	//iscale/=volume;
  }
}

// Calculate the spheroidal function
Float WaveletMSCleaner::wavelet(Float rad) {
  return 0.0;
}

Float WaveletMSCleaner::dogwavelet(Float rad, Float scaleSize1, Float scaleSize2) {
	
  if (rad <= 0) {
    return 1.0;
  } else if (rad >= 5.0) {
    return 0.0;
  }
  
    Float sigma_1 = scaleSize1;
    Float sigma_2 = scaleSize2;
    
    Float gaussian_1 = (1/(2*sigma_1*sigma_1*3.141)) * (exp(-0.5 * pow(rad/sigma_1 , 2.0)));
    Float gaussian_2 = (1/(2*sigma_2*sigma_2*3.141)) * (exp(-0.5 * pow(rad/sigma_2 , 2.0)));
    
    Float DoG_wavelet = gaussian_2 - gaussian_1;
  
    return DoG_wavelet;
}

Float WaveletMSCleaner::gaussian(Float rad, Float scaleSize) {
	
  if (rad <= 0) {
    return 1.0;
  } else if (rad >= 5.0) {
    return 0.0;
  }
  
    Float sigma = scaleSize;
    
    Float gaussian = (1/(2*sigma*sigma*3.141)) * (exp(-0.5 * pow(rad/sigma , 2.0)));
    
    return gaussian;
}

// We calculate all the scales and the corresponding convolutions
// and cross convolutions.

Bool WaveletMSCleaner::setscales(const Vector<Float>& scaleSizes)
{
  LogIO os(LogOrigin("MatrixCleaner", "setscales()", WHERE));

  Int scale;

  defineScales(scaleSizes);

  // Residual, psf, and mask, plus cross terms
  // e.g. for 5 scales this is 45. for 6 it is 60.
  Int nImages=3*itsNscales+itsNscales*(itsNscales+1);
  os << "Expect to use "  << nImages << " scratch images" << LogIO::POST;

  // Now we can update the size of memory allocated
  itsMemoryMB=0.5*Double(HostInfo::memoryTotal()/1024)/Double(nImages);
  os << "Maximum memory allocated per image "  << itsMemoryMB << "MB" << LogIO::POST;

  itsDirtyConvScales.resize(itsNscales);
  itsScaleMasks.resize(itsNscales);
  itsScaleXfrs.resize(itsNscales);
  itsPsfConvScales.resize((itsNscales+1)*(itsNscales+1));
  for(scale=0; scale<itsNscales;scale++) {
    itsDirtyConvScales[scale].resize();
    itsScaleMasks[scale].resize();
    itsScaleXfrs[scale].resize();
  }
  for(scale=0; scale<((itsNscales+1)*(itsNscales+1));scale++) {
    itsPsfConvScales[scale].resize();
  }

  AlwaysAssert(!itsDirty.null(), AipsError);

  FFTServer<Float,Complex> fft(itsDirty->shape());

  Matrix<Complex> dirtyFT;
  fft.fft0(dirtyFT, *itsDirty);

  ///Having problem with fftw with openmp
  //#pragma parallel default(shared) private(scale) firstprivate(fft)
  {
    //#pragma omp for 
    for (scale=0; scale<itsNscales;scale++) {
      os << "Calculating scale image and Fourier transform for scale " << scale << LogIO::POST;
      //cout << "Calculating scale image and Fourier transform for scale " << scale << endl;
      itsScales[scale] = Matrix<Float>(itsDirty->shape());
      //AlwaysAssert(itsScales[scale], AipsError);
      // First make the scale
      if(scale == itsNscales - 1)
		makeScaleDiff(itsScales[scale], itsScaleSizes(scale), itsScaleSizes(scale));
	  else
		makeScaleDiff(itsScales[scale], itsScaleSizes(scale), itsScaleSizes(scale+1));
      itsScaleXfrs[scale] = Matrix<Complex> ();
       fft.fft0(itsScaleXfrs[scale], itsScales[scale]);
    }
  }
  
  // Now we can do all the convolutions
  Matrix<Complex> cWork;
  for (scale=0; scale<itsNscales;scale++) {
    os << "Calculating convolutions for scale " << scale << LogIO::POST;
    
    // PSF * scale
     itsPsfConvScales[scale] = Matrix<Float>(itsDirty->shape());

    cWork=((*itsXfr)*(itsScaleXfrs[scale]));

    fft.fft0((itsPsfConvScales[scale]), cWork, false);
    fft.flip(itsPsfConvScales[scale], false, false);
    
    itsDirtyConvScales[scale] = Matrix<Float>(itsDirty->shape());
    cWork=((dirtyFT)*(itsScaleXfrs[scale]));
    fft.fft0(itsDirtyConvScales[scale], cWork, false);
    fft.flip(itsDirtyConvScales[scale], false, false);
    ///////////
    /*
    {
      String axisName = "TabularDoggies1";
      String axisUnit = "km";
      Double crval = 10.12;
      Double crpix = -128.32;
      Double cdelt = 3.145;
      TabularCoordinate tab1(crval, cdelt, crpix, axisUnit, axisName);
      TabularCoordinate tab2(crval, cdelt, crpix, axisUnit, "Dogma");
      CoordinateSystem csys;
      csys.addCoordinate(tab1);
      csys.addCoordinate(tab2);
      PagedImage<Float> limage(itsPsfConvScales[scale].shape(), csys, "psfconvscale_"+String::toString(scale));
      limage.put(itsPsfConvScales[scale]);
    }
    */
      ///////////

    for (Int otherscale=scale;otherscale<itsNscales;otherscale++) {
      
      AlwaysAssert(index(scale, otherscale)<Int(itsPsfConvScales.nelements()),
		   AipsError);
      
      // PSF *  scale * otherscale
       itsPsfConvScales[index(scale,otherscale)] =Matrix<Float>(itsDirty->shape());
      cWork=((*itsXfr)*conj(itsScaleXfrs[scale])*(itsScaleXfrs[otherscale]));
      fft.fft0(itsPsfConvScales[index(scale,otherscale)], cWork, false);
      //fft.flip(*itsPsfConvScales[index(scale,otherscale)], false, false);
    }
  }

  itsScalesValid=true;

  if (!itsMask.null()) {
    makeScaleMasks();
  }

  return true;
}

void WaveletMSCleaner::makePsfScales(){
  LogIO os(LogOrigin("MatrixCleaner", "makePsfScales()", WHERE));
  if(itsNscales < 1)
    throw(AipsError("Scales have to be set"));
  if(itsXfr.null())
    throw(AipsError("Psf is not defined"));
  destroyScales();
  itsScales.resize(itsNscales, true);
  itsScaleXfrs.resize(itsNscales, true);
  itsPsfConvScales.resize((itsNscales+1)*(itsNscales+1), true);
  FFTServer<Float,Complex> fft(psfShape_p);
  Int scale=0;
  for(scale=0; scale<itsNscales;scale++) {
    itsScales[scale] = Matrix<Float>(psfShape_p);
    if(scale == itsNscales - 1)
		makeScaleDiff(itsScales[scale], itsScaleSizes(scale), itsScaleSizes(scale));
	else
		makeScaleDiff(itsScales[scale], itsScaleSizes(scale), itsScaleSizes(scale+1));
    itsScaleXfrs[scale] = Matrix<Complex> ();
    fft.fft0(itsScaleXfrs[scale], itsScales[scale]);
  }
  Matrix<Complex> cWork;
  
  for (scale=0; scale<itsNscales;scale++) {
    os << "Calculating convolutions for scale " << scale << LogIO::POST;
    //PSF * scale
    itsPsfConvScales[scale] = Matrix<Float>(psfShape_p);
    cWork=((*itsXfr)*(itsScaleXfrs[scale])*(itsScaleXfrs[scale]));
    //cout << "shape "  << cWork.shape() << "   " << itsPsfConvScales[scale].shape() << endl;

    fft.fft0((itsPsfConvScales[scale]), cWork, false);
    fft.flip(itsPsfConvScales[scale], false, false);
    
    //cout << "psf scale " << scale << " " << max(itsPsfConvScales[scale]) << " " << min(itsPsfConvScales[scale]) << endl;

    for (Int otherscale=scale;otherscale<itsNscales;otherscale++) {
      
      AlwaysAssert(index(scale, otherscale)<Int(itsPsfConvScales.nelements()),
		   AipsError);
      
      // PSF *  scale * otherscale
      itsPsfConvScales[index(scale,otherscale)] =Matrix<Float>(psfShape_p);
      cWork=((*itsXfr)*(itsScaleXfrs[scale])*(itsScaleXfrs[otherscale]));
      fft.fft0(itsPsfConvScales[index(scale,otherscale)], cWork, false);
      //For some reason this complex->real fft  does not need a flip ...may be because conj(a)*a is real
      //fft.flip(*itsPsfConvScales[index(scale,otherscale)], false, false);
    }
  }
  
  itsScalesValid=true;

}

} //# NAMESPACE CASA - END
