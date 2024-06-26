

LKrigNormalizeBasisFFTInterpolate_scale <- function(LKinfo, Level, x1){
  # This functions evaluates the variance of the basis functions on a coarser grid, 
  # then uses 2D interpolation via FFT in order to smooth/interpolate the variance 
  # up to the size of the original grid that we were working with. Should provide a 
  # significant computational speedup. 
  # big grid N must be a multiple of small grid n, N= Mn - M + 1
  
  cat("Using FFT normalization method for level", Level, fill = TRUE)
  
  
  
  # Extracting important information from LKinfo 
  bounds <- cbind(c(min(LKinfo$x[,1]), max(LKinfo$x[,1])), 
                  c(min(LKinfo$x[,2]), max(LKinfo$x[,2])))
  basisNum_big <- max(LKinfo$latticeInfo$mxDomain[Level,1], 
                      LKinfo$latticeInfo$mxDomain[Level,2])
  basisNum_small <- min(LKinfo$latticeInfo$mxDomain[Level,1], 
                        LKinfo$latticeInfo$mxDomain[Level,2])
  gridOrientation <- which.max(c(LKinfo$latticeInfo$mxDomain[Level,1],
                                 LKinfo$latticeInfo$mxDomain[Level,2]))
  
  buffer <- LKinfo$NC.buffer
  alphaNum <- LKinfo$alpha[Level]
  awght <- LKinfo$a.wght[Level]
  
  #dimensions of original data, also helpful for shift parameter
  nr <- length(unique(x1[,1]))
  nc <- length(unique(x1[,2]))
  maxDimension <- max(nr, nc)
  minDimension <- min(nr, nc)
  
  # Setting up a new, single level LKrig object using the extracted info 
  LKinfoNew <- LKrigSetup(bounds, nlevel = 1, NC = basisNum_big, NC.buffer = buffer, 
                          alpha = alphaNum, a.wght = awght, normalize = FALSE) 
  
  # Setting a default coarse grid size based on the number of basis functions 
  # MINIMUM VALUE is 2 * basisNum - 1
  # NOTE: can play with this for accuracy
  scalefactor <- 4
  miniGridSize_big <- scalefactor * basisNum_big - 3
  miniGridSize_small <- scalefactor * basisNum_small - 3
  
  if (miniGridSize_big >= maxDimension || miniGridSize_small >= minDimension) {
    stop("Warning: Minimum coarse grid based on the number of basis functions is 
         greater than the size of the data. This method is not appropriate here. 
         Either choose less basis functions, or choose a different method, 
         such as exactKronecker or fast. See help file on 
         LKrigNormalizeBasisFFTInterpolate.")
  }
  
  else{
    
    # Creating the actual grid
    # if the first row of basis funcs is larger, the y is the bigger side
    if (gridOrientation == 1){
      gridList<- list( x= seq( bounds[1,1],bounds[2,1],length.out = miniGridSize_big),
                       y= seq( bounds[1,2],bounds[2,2],length.out = miniGridSize_small) )
    }
    # if the second row of basis funcs is larger, the x is the bigger side
    if (gridOrientation == 2){
      gridList<- list( x= seq( bounds[1,1],bounds[2,1],length.out = miniGridSize_small),
                       y= seq( bounds[1,2],bounds[2,2],length.out = miniGridSize_big) )
    }
    
    sGrid<- make.surface.grid(gridList)
    
    
    # Calling LKrig.cov to evaluate the variance on the coarse grid
    # this is the initial, small variance calculation that we will upsample
    
    #exact 
    #roughMat <- as.surface(sGrid, LKrig.cov(sGrid, LKinfo = LKinfoNew, marginal=TRUE ))[["z"]]
    
    #exact but kronecker (should add if statement for kroneck conditions)
    roughMat <- as.surface(sGrid, LKrigNormalizeBasisFast(LKinfo = LKinfoNew, Level = 1, x = sGrid))[["z"]]
    

    # FFT step: taking the fft of the small variance
    # reliant on fftwtools package
    fftStep <- fftw2d((roughMat)/length(roughMat))

    # Helpful dimensions
    snr <- nrow(fftStep) 
    snc <- ncol(fftStep) 
    
    M <- ceiling(maxDimension/miniGridSize_big)
    
    NR <- M * snr
    NC <- M * snc
    
    # Instantiate empty matrix
    temp <- matrix(0, nrow = NR, ncol = NC)
    
    # Helpful for indexing later 
    bigindY <- 1:ceiling(snr/2) 
    bigindX <- 1:ceiling(snc/2) 
    indY <- 1:(snr/2) 
    indX <- 1:(snc/2) 
    # helpful offsets
    bigOffsetY <- (NR - floor(snr/2)) 
    bigOffsetX <- (NC - floor(snc/2))
    smallOffsetY <- (snr - floor(snr/2)) 
    smallOffsetX <- (snc - floor(snc/2))
    
    # Stuffing the small FFT result into the large matrix of zeroes 
    temp[bigindY, bigindX] <- fftStep[bigindY, bigindX] #top left corner
    temp[bigindY, (indX + bigOffsetX)] <- fftStep[bigindY, (indX + smallOffsetX)] #top right corner
    temp[(indY + bigOffsetY), bigindX] <- fftStep[(indY + smallOffsetY), bigindX] #bottom left corner 
    temp[(indY + bigOffsetY), (indX + bigOffsetX)] <- fftStep[(indY + smallOffsetY), (indX + smallOffsetX)] #bottom right corner
    
    # takes the IFFT of the modified big matrix to return our interpolated/upsampled variance
    # again reliant on fftwtools package
    wght <- Re(fftw2d(temp, inverse = 1))
    
    # trim interpolated array to be within bounds of source image
    wght<- wght[1: (NR - M + 1), 1: (NC - M + 1)]

    # Shifting due to the periodicity assumed by the FFT
    #wght <- LKrig.shift.matrix(wght, shift.row = yShift, shift.col = xShift, periodic = c(TRUE, TRUE))
    
    #the next steps are done to account for oddly sized/missing data
    # the fft calculation will produce a full grid
    # but only certain observations need a variance associated with them
    
    lat_min <- min(x1[,1]) # minimum latitude
    lat_max <- max(x1[,1]) # maximum latitude
    lon_min <- min(x1[,2]) # minimum longitude
    lon_max <- max(x1[,2]) # maximum longitude
    
    # the grid resolution is known, but if one needs to calculate steps based on nr and nc:
    lat_step <- (lat_max - lat_min) / (nr - 1)
    lon_step <- (lon_max - lon_min) / (nc - 1)
    
    # Map lat-lon to grid indices
    row_indices <- round((x1[,1] - lat_min) / lat_step) + 1
    col_indices <- round((x1[,2] - lon_min) / lon_step) + 1
    
    # Ensure indices are within the bounds and adjust if needed
    row_indices <- pmin(pmax(row_indices, 1), nr)
    col_indices <- pmin(pmax(col_indices, 1), nc)
    
    # Extract the actual values associated with existing data using the indices
    values <- wght[cbind(row_indices, col_indices)]
    
    # vectorize the output matrix to be compatible with LKrig.basis
    wght <- c(values)
  }
  
  return (wght)
  
}
