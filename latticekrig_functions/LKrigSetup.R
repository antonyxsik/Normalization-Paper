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

LKrigSetup <- function(x = NULL,
                       nlevel=NULL,
# the following attributes of the covariance models can 
# vary in their class depending on the geometry and if they are
# stationary or not. 
                       alpha=NA, alphaObject=NULL, nu = NULL,
                       a.wght = NA, a.wghtObject=NULL,
# arguments that set the initial lattice size
# for a square rectangle domain there are NCxNC lattice points
# with NC.buffer added on the 4 edges.
# NC must be specified 
# defaults for NC.buffer are set in setDefaultsLKinfo
                       NC = NULL, NC.buffer=NULL,
# default normalization method is exact as FFT, Kronecker, 
# and both have more limited use cases and depend on the 
# geometry of the problem
                       normalize=TRUE, normalizeMethod = "exact",
                       lambda = NA, sigma = NA, rho = NA, rho.object = NULL,
                       latticeInfo=NULL, basisInfo=NULL, 
# default geometry is a rectangular domain with Euclidean distance 
# (see NC above)
                       LKGeometry="LKRectangle",
# even if the geometry is different than rectangle many of the following defaults
# still are appropriate                       
                       distance.type="Euclidean",
# some natural defaults for the basis functions (including the default to be radial)                      
                       BasisFunction = "WendlandFunction", overlap=2.5,                            
                       V=NULL,
                       BasisType= "Radial",
# some natural defaults for the fixed part of the model
# NOTE: these can be modified according to geometry in 
# setDefaultsLKinfo below.
                       fixedFunction="LKrigDefaultFixedFunction",    
                       fixedFunctionArgs = list(m=2),
                       collapseFixedEffect = FALSE,
# defaults for sparse matrix size.                        
                       max.points=NULL, mean.neighbor=50, choleskyMemory=NULL,
# useful for debugging                       
                       verbose = FALSE, noCheck=FALSE,
                       returnCall = FALSE,
                       dense=FALSE, 
# these additional arguments will just be added as a list to the LKinfo object as setupArgs
                          ... ) { 
#
# create initial shell of LKinfo object with classes and common parameters
   setupArgs<- list( ...)
# the initial LKinfo object list that will be modified and added to 
# below
   LKinfo<-  list(      x = x,
                   nlevel = nlevel,
                    alpha = alpha,
              alphaObject = alphaObject,
                   a.wght = a.wght,
             a.wghtObject = a.wghtObject,
                       NC = NC,
                NC.buffer = NC.buffer,
                       nu = nu,
                normalize = normalize,
          normalizeMethod = normalizeMethod,
                   lambda = lambda,
                    sigma = sigma,
                      rho = rho,
               rho.object = rho.object,
              LKGeometry  = LKGeometry,
            distance.type = distance.type, 
            BasisFunction = BasisFunction,
                  overlap = overlap,                            
                        V = V,
                BasisType = BasisType,
            fixedFunction = fixedFunction,
        fixedFunctionArgs = fixedFunctionArgs,
      collapseFixedEffect = collapseFixedEffect,
               max.points = max.points,
            mean.neighbor = mean.neighbor,
           choleskyMemory = choleskyMemory,
                setupArgs =  setupArgs,
                    dense = dense
        
                 ) 
# 
    LKinfo$basisInfo<- list(             BasisType = BasisType,
                                     BasisFunction = BasisFunction,
                                           overlap = overlap,
                                        max.points = max.points,
                                     mean.neighbor = mean.neighbor,
                                                 V = V
                              )
   if( verbose){
     temp<- LKinfo
     class( temp) <- NULL
     print( as.list(temp))
   }
# set the geometry class -- this will determine what
# functions are applied below in the "Setup" calls
   class(LKinfo) <- c( "LKinfo", LKGeometry)
# the next function is used to reset obvious defaults for
# a particular geometry. The default method here is to 
# do nothing but for most geometries it makes sense to 
# set some options to default values.
# or to add some additional information to the LKinfo object.   
   LKinfo<- setDefaultsLKinfo( LKinfo )
   if( verbose){
    	cat("----- After call to setDefaultsLKinfo -----", fill=TRUE)
     temp<- LKinfo
     class( temp) <- NULL
     print( as.list(temp))
   }
#    
# Create information to construct the lattice at each level based on the
# geometry of spatial domain, number of levels, and the data locations
# this function is overloaded and is determine by LKGeometry
# Note that extra arguments given to LKrigSetup are passed through to
# this function setting up lattice.
# there is no default method
# conceptually the following call is jsut LKrigSetupLattice( LKinfo, x, verbose, ...)
# but this does not work because the additional arguments in list( ...) are 
# not handled correctly
   LKinfo$latticeInfo<- do.call( "LKrigSetupLattice",
                 c( list( object = LKinfo, verbose = verbose ),
                    setupArgs)
                  ) 
#
# reformat, modify, and check the parameters for the Markox random field/ GP model
#  fix up the alpha parameters
#   the default method is probably adequate for most geometries and SARs
    LKinfo$alpha<- LKrigSetupAlpha(LKinfo)
    if( verbose){
        print(alpha)}
# fix up the a.wght parameters specfic geometries might need 
# a specific function here.
    LKinfo$a.wght<-LKrigSetupAwght(LKinfo)
# set lambda if sigma and rho are passed.
    if (is.na(lambda[1])) {
        lambda <- sigma^2/rho
        LKinfo$lambda<- lambda
    }  
# Note saving call argument in return allows the function 
#      to be re evaluated
    if ( returnCall){
        LKinfo$call<- match.call()
          }
    else{ 
        LKinfo$call<- NULL
        }
# Finally make some generic checks that all basic components of
# LKinfo class are included
   if( !noCheck){
    LKinfoCheck(LKinfo)
  }
   return(LKinfo)
}