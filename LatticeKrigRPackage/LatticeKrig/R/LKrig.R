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

LKrig <- function(x, y,
                 weights = NULL,
                       Z = NULL,
                  LKinfo = NULL, 
                   iseed = NA, 
	                NtrA = 20, 
            use.cholesky = NULL,
         return.cholesky = TRUE,
                	   X = NULL, 
                	   U = NULL,
                	  wX = NULL,
                	  wU = NULL,
          return.wXandwU = TRUE,
               	            ...,
             getVarNames = TRUE,
               	 verbose = FALSE) {
# if LKinfo is missing create it from passed arguments 
# if it is passed update this object with the ... arguments
# for example a new lambda value can be passed in this part.
# LKinfo is a list that describes the LatticeKrig model.
#   the list(...) only pertains to LKinfo arguments  
	if ( is.null(LKinfo) ) {
		LKinfo <- do.call("LKrigSetup", c(list(x = x ), list(...),
		 list(verbose = verbose)))
	}
	else{
		LKinfo<- do.call("LKinfoUpdate", c(list(LKinfo=LKinfo), list(...)) )	
	}
	
	if( verbose){
		cat(" ", fill=TRUE)
		cat("LKrig: updated LKinfo object", fill=TRUE)
		print(LKinfo)
		
	}
#  
# the variable names are used to generate column labels if 
# those are missing. 
# the getVarNames switch is important -- if
# LKrig is called via the do.call function then strange things
# happen with the substitute function and I understand is an unresolved 
# aspect of Rli
#
  if( getVarNames){
    xName<- as.character( substitute( x) )
    ZName<- as.character( substitute( Z) )
    UName<- as.character( substitute( U) )
    xName<- tail( xName, 1)
    ZName<- tail( ZName, 1)
    UName<- tail( UName, 1)
    
# just take last component
    
  }
  else{
    xName<- "xVar"
    ZName<- "ZVar"
    UName<- "UVar"
  }  
#
# create the initial parts of LKrig object
# this list is added to as the computation proceeds 
# using the device  object<- list( object, newStuff)
# and the full object is only obtained at the end 
# NOTE default for weights are just 1's and filled in by 
# the next call    

    object<- createLKrigObject( x, y, 
                              weights = weights,
                                    Z = Z,
                                    X = X,
                                    U = U, 
                               LKinfo = LKinfo,
                                xName = xName, 
                                ZName = ZName,
                                UName = UName, 
                              verbose = verbose)
  
    nObs <-  nrow( object$y )
    nReps <- ncol( object$y )
# for readablity make a local copy of LKinfo
# but don't change it in this function! 
    LKinfo<- object$LKinfo                                              	 	
	# Begin computations ....
	# weighted observation vector
    wy <- sqrt(object$weights) * object$y
   
# create matrix for fixed part of model    
# Spatial drift matrix -- default is assumed to be linear in coordinates (m=2)
# and includes possible covariate(s) -- the Z matrix.
# the choice of fixed part of the model is controlled in LKinfo
# (see also LKrigSetup)
   if (is.null(wU)) {
   	wU<- LKrigMakewU( object,  verbose=verbose)
   }
   
# some column indices to keep track of fixed part of the model	
# NOTE nZ <= nt because Z is a subset of U
    object$nt <- ifelse( is.null(ncol(wU)), 0, ncol(wU))
# create matrix for random part of model (basis functions)
#  wX is the matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
# and multiplied by square root of diagonal weight matrix
# this can be a large matrix if not encoded in sparse format.
if (is.null(wX)) {
      timewX<- system.time(
          wX<- LKrigMakewX( object, verbose=verbose)
          )	
	}
else{
	timewX<- rep(0,5)
}	
 
 #   Precision matrix of the lattice process
#   inverse of Q is proportional to the covariance matrix of the Markov Random Field
timeQ<-system.time(
        Q <- LKrig.precision(LKinfo, verbose=verbose)
        )
if( LKinfo$dense){
  if( !is.null(  use.cholesky)){
    stop("Can not update (use.cholesky) with dense matrix option")
  }
  Q<- spam2full(Q)
  wX<- spam2full(wX)
}

if( verbose & (!LKinfo$dense) ) {
  cat("LKrig: Nonzero entries in Q:", length(Q@entries), fill=TRUE)		
}

# G is the regularized (ridge) regression matrix that is 
# the key to the entire algorithm:
timeM<- system.time(	
	G <- t(wX) %*% wX + LKinfo$lambda * (Q)
	)
	if( verbose){
	  if( !LKinfo$dense){
	    cat("LKrig: Nonzero entries in M:", length(G@entries), fill=TRUE)	
	  }
	  else{
	    cat( "Dense matrix methods used", fill=TRUE)
	  }
	}	
#  Find Cholesky square root of M
#  This is where the heavy lifting happens!  M is in sparse, spam format so
#  by the overloading this is a sparse cholesky decomposition.
#  if this function has been coded efficiently this step should dominate
#  all other computations ...
#  If a previous sparse cholesky decoposition is passed then the
#  pattern of sparseness is used for the decoposition.
#  This can speed the computation as the symbolic decomposition part of the
#  sparse Cholesky is a nontrivial step. The condition is that
#  the current 'M' matrix  has the same sparse pattern as that
#  which  will result in the same sparse pattern of factorization 
#  as when chol was applied to the case passed in 'use.cholesky'
if (is.null(use.cholesky)) {
	timeChol<- system.time(
		GCholesky <- chol(G, memory = LKinfo$choleskyMemory)
		)
	} else {	
			timeChol<- system.time(
		GCholesky <- update.spam.chol.NgPeyton(use.cholesky, G)
		)
	}

if( !LKinfo$dense){
   nonzero.entries<- length(GCholesky@entries)
}
   else{
     nonzero.entries<-NA
   }


  if( verbose){
     	cat("LKrig: nonzero entries of GCholesky:",nonzero.entries, fill=TRUE)
  }
# use GCholesky to find coefficients of estimate
# Note that this functions also finds an important piece of the likelihood (quad.form)
  timeCoef<- system.time(
	out1 <- LKrig.coef(GCholesky, wX, wU, wy,
	               LKinfo$lambda, 
	   collapseFixedEffect = LKinfo$collapseFixedEffect, 
	               verbose=verbose)
	)
	
# Note collapseFixedEffect added as component here in the return
# finding coefficients
# fill in names of the fixed coefficients 
# fill in names of fixed model coefficients
	
	rownames(out1$d.coef)<- colnames( wU )
#
	object <- c(object, out1)
	if( verbose){
	  cat("fixed model coefficients", fill=TRUE)
	  cat( object$d.coef, fill=TRUE)
	}
	
# compute predicted values  and residuals
	wfitted.values <- (wX %*% out1$c.coef)
	if ( !is.null(wU) ) {
		wfitted.values.fixed <- (wU %*% out1$d.coef)
		wfitted.values <- wfitted.values.fixed + wfitted.values
	}
# X and U actully include the weights so need to divide these
# out to get fitted values	
	object$fitted.values<- wfitted.values/sqrt(object$weights)		
# For reference: fitted.values <- predict.LKrig(object, x, Znew = object$Z)
# but at this point it is less efficient because X will be recomputed.
	object$residuals <- object$y - object$fitted.values	
# find likelihood
timeLike<- system.time(	
#	out2 <- LKrig.lnPlikeOLD(GCholesky, Q, wy,
#	             object$residuals, object$weights,
#	             LKinfo)
  out2 <- LKrig.lnPlike(GCholesky, Q, object$quad.form,
                           nObs, nReps,
                           object$weights,LKinfo)
  )
	if( verbose){
		cat("Likelihood/MLE list:",  fill=TRUE)
		print( out2)
	}   
	object <- c(object, out2)
	
# estimate trace of hat matrix (effective degrees of freedom)
# by Monte Carlo if NtrA greater than zero
timeTrA<- system.time(
	if (NtrA > 0) {
		out3 <- LKrig.traceA(GCholesky, wX, wU, LKinfo$lambda, object$weights, NtrA, iseed = iseed)
		# find GCV using this trace estimate
		n<- length( object$weights)
		out3$GCV = (sum(object$weights * (object$residuals)^2)/n)/(1 - out3$trA.est/n)^2
	} else {
		out3 <- list(trA.est = NA, trA.SE = NA, GCV = NA)
	}
	)
	object <- c(object, out3)
	
# create the table of times for individual function calls
timingTable<- rbind(timewX, timeQ, timeM, timeChol, timeCoef, timeLike,  timeTrA)
timingTable<- timingTable[,1:3]
timingTable <- rbind( timingTable, colSums(timingTable))

# last of required arguments to LKrig object
object <- c(object,
 list( 
        lambda.fixed = LKinfo$lambda, 
     nonzero.entries = nonzero.entries,
                call = match.call(), 
         timingLKrig = timingTable )
    )
# finally add in some large matrices that can be reused if only 
# lambda is varied  
# NOTE to keep the code with the same components the name Mc is kept although
# this is actually the cholesky decomposition of the G matrix as in the LKrig
# article. 
if (return.cholesky) {
		object$Mc <- GCholesky
	}
if (return.wXandwU) {
		object$wX <- wX
		object$wU <- wU
	}
# refresh the class 
    class(object) <- "LKrig"
# all done!  
	return(object)
}

