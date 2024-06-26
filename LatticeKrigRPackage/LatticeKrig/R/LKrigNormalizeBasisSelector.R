# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKrigNormalizeBasisSelector <- function(LKinfo, Level, x1, verbose){
  
  # Extracting number of basis functions from LKinfo 
  basisNum <- max(LKinfo$latticeInfo$mxDomain[Level,1], LKinfo$latticeInfo$mxDomain[Level,2])
  
  # Dimensions of original data
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  maxDimension <- max(nr, nc)
  miniGridSize <- 4 * basisNum
  
  # Method selection
  # If coarse grid size is less than the size of the data, use FFT
  if (Level < 4 && miniGridSize < maxDimension){
    if (verbose){
      cat("Using FFT Interpolation method for level", Level, fill = TRUE)
    }
    # Old (for variable sizes)
    #wght <- LKrigNormalizeBasisFFTInterpolate(LKinfo, Level, x1)
  
    # For scale factor 
    wght <- LKrigNormalizeBasisFFTInterpolate_scale(LKinfo, Level, x1)
  }
  
  # When coarse grid size gets too big (too many basis functions)
  else {
    if (verbose){
      cat("Using Kronecker method for level", Level, fill = TRUE)
    }
    wght <- LKrigNormalizeBasisFast(LKinfo,  Level,  x1)
  }
  
  return(wght)
}
