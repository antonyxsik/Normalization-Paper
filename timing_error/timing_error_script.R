#SETUP ####################################################################
###########################################################################


# this script takes long to run on a powerful laptop (roughly 5 nights)
# contains all experiments that create dataframes for figure_generation.Rmd


#Make sure to set correct working directory
setwd("C:/Users/anton/Desktop/Normalization_Paper_Github/Normalization-Paper")


#UNCOMMENT IF YOU HAVE NOT INSTALLED LATTICEKRIG 
# install.packages("LatticeKrig")


#libraries
library(LatticeKrig)
library(spam64)
library(fftwtools)


#LOAD THE RDA FILE AND ADD TO IT
load("timing_error/dataframes/normalization_times_error.rda")


###########################################################################
# THIS IS PURE TIMING FOR LARGE MATRICES, PROVIDES BOTH ERROR AND TIME
###########################################################################

normalization_timer <- function(sideLengths, basisFuncList, numBuffer, l, numIters, 
                                df, checkpointPath){
  
  # sideLengths is a list of side lengths of a large square (for simplicity) matrix
  # that we want to time. Ex: c(200, 400, 1000, 2000)
  # 
  # basisfunclist is a list of basis function numbers 
  # 
  # numBuffer is an integer number of additional outside basis functions that reduce 
  # edge effects. Ex: 5
  # 
  # l we should always set to 1. it is the integer number of levels in the 
  # LatticeKrig model.
  #
  # numIters is an integer number of times that we will go through the list
  # we go through it multiple times to assess variability and calculate 
  # statistics like averages. Ex: 10
  #
  # df is the dataframe object that we will be adding rows to
  # Ex: df_timing
  # 
  # checkPointPath is a way of temporarily saving the df 
  
  if(file.exists(checkpointPath)){
    load(checkpointPath)
    df <- df
  }
  
  xr <- c(-1,1)
  yr <- c(-1,1)
  
  # Iteration is the outermost loop
  for (iteration in 1:numIters) {
    
    # Loop through different numbers of basis functions
    for (numBasis in basisFuncList) {
      
      for (size in sideLengths){ 
        
        # Setting up the grid/x1
        gridList <- list(x = seq(-1, 1, length.out = size),
                         y = seq(-1, 1, length.out = size))
        x1 <- make.surface.grid(gridList)
        
        
        #no normalization 
        LKinfo_none <- LKrigSetup(cbind(xr, yr), nlevel=l, NC= numBasis, NC.buffer = numBuffer,
                                  a.wght= 4.05, normalize= FALSE) 
        
        #exact normalization (slow)
        LKinfo_exact <- LKrigSetup(cbind(xr, yr), nlevel=l, NC= numBasis, NC.buffer = numBuffer,
                                   a.wght= 4.05, normalize= TRUE, normalizeMethod = "exact")
        
        #exact Kronecker normalization (faster)
        LKinfo_kroneck <- LKrigSetup(cbind(xr, yr), nlevel=l, NC= numBasis, NC.buffer = numBuffer,
                                     a.wght= 4.05, normalize= TRUE, normalizeMethod = "exactKronecker")
        
        #approximate fft normalization (fastest)
        LKinfo_fft <- LKrigSetup(cbind(xr, yr), nlevel=l, NC= numBasis, NC.buffer = numBuffer,
                                 a.wght= 4.05, normalize= TRUE, normalizeMethod = "fftInterpolation")
        
        
        #storing the times
        none_time <- system.time(LKrig.basis(x1, LKinfo_none, verbose = FALSE))[[3]]
        exact_time <- system.time(LKrig.basis(x1, LKinfo_exact, verbose = FALSE))[[3]]
        kroneck_time <- system.time(LKrig.basis(x1, LKinfo_kroneck, verbose = FALSE))[[3]]
        fft_time <- system.time(LKrig.basis(x1, LKinfo_fft, verbose = FALSE))[[3]]
        
        
        #doing these to calculate exact to fft error
        kroneck_var <- LKrigNormalizeBasisFast(LKinfo_none,  Level=l,  x=x1)
        fft_var <- LKrigNormalizeBasisFFTInterpolate(LKinfo_none, Level=l, x1=x1)
        
        
        #percent error calculations (average and max absolute)
        AvgPercent <- mean((kroneck_var - fft_var)/kroneck_var) * 100
        MaxPercent <- max(abs(kroneck_var - fft_var)/kroneck_var) * 100
        
        
        # Add a row to the data frame
        df <- rbind(df, c(iteration, size, numBasis, none_time, exact_time, kroneck_time, fft_time, 
                          AvgPercent, MaxPercent))
        
        print(paste(numBasis, size))
        
        #checkpoint save
        save(df, file = checkpointPath)
        
        # Update the column names to include iteration, size, and number of basis functions
        colnames(df) <- c("Iteration", "Size", "NumBasis", "None", "Exact", "Kronecker", "FFT", 
                          "Mean % Error", "Max % Error")
      }
    }
  }
  
  return(df)
}



# Actual params for timing experiment
varying_basis_funcs <- c(15, 25, 35, 50, 100)
varying_sidelength <- c(500, 750, 1000, 1250, 1500, 2000)
numBuffer <- 10
l <- 1
numIters <- 1
tempPath <- "timing_error/dataframes/normalization_times_error.rda"
df_timing_error <- data.frame()
system.time(result_df <- normalization_timer(sideLengths = varying_sidelength, 
                                 basisFuncList = varying_basis_funcs, 
                                 numBuffer = numBuffer, 
                                 l = l, 
                                 numIters = numIters, 
                                 df = df_timing_error, 
                                 checkpointPath = tempPath))

print(result_df)
result_df <- result_df[order(result_df$NumBasis, result_df$Size),]

save(df, file = "timing_error/dataframes/normalization_times_error.rda")
###########################################################################
###########################################################################





###########################################################################
# THIS EXPERIMENT FOCUSES ON SEEING HOW MUCH THE INITIAL COARSE GRID SIZE IMPACTS THE ERROR
###########################################################################

#function for testing
fftSmoother <- function(LKinfo, sGrid, nr, nc, plots = FALSE, col = turbo(256)){
  
  #fft step and dimensions
  roughMat <- as.surface(sGrid, LKrig.cov(sGrid, LKinfo = LKinfo, marginal= TRUE ))[["z"]]
  fftStep <- fftw2d((roughMat)/length(roughMat))
  snr <- nrow(fftStep)
  snc <- ncol(fftStep)
  
  #instantiate empty matrix
  temp <- matrix(0, nrow = nr, ncol = nc)
  
  indY <- 1:(snr/2)
  indX <- 1:(snc/2)
  bigOffsetY <- (nr - floor(snr/2))
  bigOffsetX <- (nc - floor(snc/2))
  smallOffsetY <- (snr - floor(snr/2))
  smallOffsetX <- (snc - floor(snc/2))
  
  
  temp[indY, indX] <- fftStep[indY, indX] #top left corner
  temp[indY, (indX + bigOffsetX)] <- fftStep[indY, (indX + smallOffsetX)] #top right corner
  temp[(indY + bigOffsetY), indX] <- fftStep[(indY + smallOffsetY), indX] #bottom left corner 
  temp[(indY + bigOffsetY), (indX + bigOffsetX)] <- fftStep[(indY + smallOffsetY), (indX + smallOffsetX)] #bottom right corner
  
  #takes the inverse of the modified big matrix to return a smoothed version of the input matrix
  smoothMat <- Re(fftw2d(temp, inverse = 1))
  
  #shifting based on how large the size-up is (not sure which one to use, ceiling is safer imo)
  #yShift <- ceiling(((nr/snr) - 1)/2)
  #xShift <- ceiling(((nc/snc) - 1)/2)
  
  #uncomment the next 3 lines
  yShift <- floor((nr/snr)/2)
  xShift <- floor((nc/snc)/2)
  
  smoothMat <- LKrig.shift.matrix(smoothMat, shift.row = yShift, shift.col = xShift, periodic = c(TRUE, TRUE))
  
  if (plots == TRUE){
    
    #plotting params
    par(mfrow = c(2,2))
    
    #capture the axis of the original grid (the -1 to 1 part) (take from the object)
    xboundLeft <- as.numeric(LKinfo$x[1,1])
    xboundRight <- as.numeric(LKinfo$x[2,1])
    yBoundLow <- as.numeric(LKinfo$x[1,2])
    yBoundHigh <- as.numeric(LKinfo$x[2,2])
    
    smallXaxis <- seq(xboundLeft,xboundRight, length.out = snc)
    smallYaxis <- seq(yBoundLow,yBoundHigh, length.out = snr)
    bigXaxis <- seq(xboundLeft,xboundRight, length.out = nc)
    bisYaxis <- seq(yBoundLow,yBoundHigh, length.out = nr)
    
    imagePlot(smallXaxis, smallYaxis, roughMat, col = col)
    imagePlot(smallXaxis, smallYaxis, log(Mod(fftStep) + 0.000000001), col = col)
    #necessary because log(0) will not work (for plotting)
    imagePlot(bigXaxis, bisYaxis, log(Mod(temp) + 0.000000001), col = col)
    imagePlot(bigXaxis, bisYaxis, smoothMat, col = col)
    
    #plotting params
    par(mfrow = c(1,1))
    
    #delete this
    imagePlot(bigXaxis, bisYaxis, smoothMat, col = col)
  }
  
  return(smoothMat)
}

df_CoarseGrid <- data.frame()
  
proveItWorks <- function(ORIGINALGRID, BIGGRID, NC = 25, NC.buffer = 10, AWGHT = 4.05){
    # create a simple 1 level model
        xr<- c( -1,1)
        yr<- c( -1,1)
    LKinfo<- LKrigSetup(cbind(xr, yr), nlevel=1, NC=NC, NC.buffer = NC.buffer, 
                        a.wght= AWGHT, normalize= FALSE)
    
    #CHANGE HERE 
    gridList<- list( x= seq( -1,1,length.out= ORIGINALGRID),
                     y= seq( -1,1,length.out= ORIGINALGRID) )
    sGrid<- make.surface.grid(gridList)
    look<- LKrig.cov(sGrid, LKinfo = LKinfo, marginal=TRUE )
    #imagePlot( as.surface( sGrid, look), main = "Small Original", col = turbo(256))
    
    
    Time <- system.time(tester <- fftSmoother(LKinfo, sGrid, nr = BIGGRID, nc = BIGGRID, plots = FALSE))[[3]]
    
    par(mfrow = c(1,3), cex.axis = 0.6)
    #imagePlot(tester, col = turbo(256))
    #title("Smoothed")
    
    
    gridList1 <- list( x= seq( -1,1,length.out=BIGGRID),
                       y= seq( -1,1,length.out=BIGGRID) )
    
    sGrid1 <- make.surface.grid(gridList1)
    RawTime <- system.time(look1 <- LKrig.cov(sGrid1, LKinfo = LKinfo, marginal=TRUE ))[[3]]
    #imagePlot( as.surface( sGrid1, look1) , main = "Big Original", col = turbo(256))
    
    #imagePlot((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]], col = turbo(256), main = "Error (%)")
    
    MeanErrorPercent <- mean(((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]]))* 100
    MaxErrorPercent <- max(abs((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]]))* 100
    
    df_CoarseGrid <- rbind(df_CoarseGrid, c(ORIGINALGRID, BIGGRID, NC, NC.buffer, AWGHT, MeanErrorPercent, 
                                            MaxErrorPercent, Time, RawTime))
    
    colnames(df_CoarseGrid) <- c("Coarse Grid", "Fine Grid", "NumBasis", "NumBuffer", "A.wght", "MeanErrorPercent", "MaxErrorPercent", 
                                 "FFT Time (s)", "Raw Time (s)")
    
    print(MaxErrorPercent)
    print(MeanErrorPercent)
    par(mfrow = c(1,1))
    
  return(df_CoarseGrid)
}

for (i in c(250, 500, 740, 750)){
  print(i)
  df_CoarseGrid <- proveItWorks(ORIGINALGRID = i, BIGGRID = 750, NC = 10, NC.buffer = 10, AWGHT = 4.05)
}
  
print(df_CoarseGrid)
save(df_CoarseGrid, file = "timing_error/dataframes/coarsegrid_error.rda")
###########################################################################
###########################################################################



###########################################################################
# THIS EXPERIMENT FOCUSES ON SEEING HOW MUCH THE BUFFER BASIS FUNCTIONS IMPACT THE ERROR
###########################################################################

df_Buffer <- data.frame()

proveItWorksBuf <- function(ORIGINALGRID, BIGGRID, NC = 25, NC.buffer = 10, AWGHT = 4.05){
  # create a simple 1 level model
  xr<- c( -1,1)
  yr<- c( -1,1)
  LKinfo<- LKrigSetup(cbind(xr, yr), nlevel=1, NC=NC, NC.buffer = NC.buffer, 
                      a.wght= AWGHT, normalize= FALSE)
  
  #CHANGE HERE 
  gridList<- list( x= seq( -1,1,length.out= ORIGINALGRID),
                   y= seq( -1,1,length.out= ORIGINALGRID) )
  sGrid<- make.surface.grid(gridList)
  look<- LKrig.cov(sGrid, LKinfo = LKinfo, marginal=TRUE )
  #imagePlot( as.surface( sGrid, look), main = "Small Original", col = turbo(256))
  
  
  Time <- system.time(tester <- fftSmoother(LKinfo, sGrid, nr = BIGGRID, nc = BIGGRID, plots = FALSE))[[3]]
  
  par(mfrow = c(1,3), cex.axis = 0.6)
  #imagePlot(tester, col = turbo(256))
  #title("Smoothed")
  
  
  gridList1 <- list( x= seq( -1,1,length.out=BIGGRID),
                     y= seq( -1,1,length.out=BIGGRID) )
  
  sGrid1 <- make.surface.grid(gridList1)
  RawTime <- system.time(look1 <- LKrig.cov(sGrid1, LKinfo = LKinfo, marginal=TRUE ))[[3]]
  #imagePlot( as.surface( sGrid1, look1) , main = "Big Original", col = turbo(256))
  
  #imagePlot((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]], col = turbo(256), main = "Error (%)")
  
  MeanErrorPercent <- mean(((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]]))* 100
  MaxErrorPercent <- max(abs((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]]))* 100
  
  df_Buffer <- rbind(df_Buffer, c(ORIGINALGRID, BIGGRID, NC, NC.buffer, AWGHT, MeanErrorPercent, 
                                          MaxErrorPercent, Time, RawTime))
  
  colnames(df_Buffer) <- c("Coarse Grid", "Fine Grid", "NumBasis", "NumBuffer", "A.wght", "MeanErrorPercent", "MaxErrorPercent", 
                               "FFT Time (s)", "Raw Time (s)")
  
  print(MaxErrorPercent)
  print(MeanErrorPercent)
  par(mfrow = c(1,1))
  
  return(df_Buffer)
}

for (i in c(30, 40)){
  print(i)
  df_Buffer <- proveItWorksBuf(ORIGINALGRID = 21, BIGGRID = 750, NC = 10, NC.buffer = i, AWGHT = 4.05)
}

print(df_Buffer)
save(df_Buffer, file = "timing_error/dataframes/buffer_error.rda")
###########################################################################
###########################################################################





###########################################################################
# THIS EXPERIMENT FOCUSES ON SEEING HOW MUCH THE A.WGHT IMPACTS THE ERROR 
###########################################################################

df_awght <- data.frame()

proveItWorksBuf <- function(ORIGINALGRID, BIGGRID, NC = 25, NC.buffer = 10, AWGHT = 4.05){
  # create a simple 1 level model
  xr<- c( -1,1)
  yr<- c( -1,1)
  LKinfo<- LKrigSetup(cbind(xr, yr), nlevel=1, NC=NC, NC.buffer = NC.buffer, 
                      a.wght= AWGHT, normalize= FALSE)
  
  #CHANGE HERE 
  gridList<- list( x= seq( -1,1,length.out= ORIGINALGRID),
                   y= seq( -1,1,length.out= ORIGINALGRID) )
  sGrid<- make.surface.grid(gridList)
  look<- LKrig.cov(sGrid, LKinfo = LKinfo, marginal=TRUE )
  #imagePlot( as.surface( sGrid, look), main = "Small Original", col = turbo(256))
  
  
  Time <- system.time(tester <- fftSmoother(LKinfo, sGrid, nr = BIGGRID, nc = BIGGRID, plots = FALSE))[[3]]
  
  par(mfrow = c(1,3), cex.axis = 0.6)
  #imagePlot(tester, col = turbo(256))
  #title("Smoothed")
  
  
  gridList1 <- list( x= seq( -1,1,length.out=BIGGRID),
                     y= seq( -1,1,length.out=BIGGRID) )
  
  sGrid1 <- make.surface.grid(gridList1)
  RawTime <- system.time(look1 <- LKrig.cov(sGrid1, LKinfo = LKinfo, marginal=TRUE ))[[3]]
  #imagePlot( as.surface( sGrid1, look1) , main = "Big Original", col = turbo(256))
  
  #imagePlot((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]], col = turbo(256), main = "Error (%)")
  
  MeanErrorPercent <- mean(((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]]))* 100
  MaxErrorPercent <- max(abs((as.surface(sGrid1, look1)[["z"]] - tester)/as.surface(sGrid1, look1)[["z"]]))* 100
  
  df_awght <- rbind(df_awght, c(ORIGINALGRID, BIGGRID, NC, NC.buffer, AWGHT, MeanErrorPercent, 
                                  MaxErrorPercent, Time, RawTime))
  
  colnames(df_awght) <- c("Coarse Grid", "Fine Grid", "NumBasis", "NumBuffer", "A.wght", "MeanErrorPercent", "MaxErrorPercent", 
                           "FFT Time (s)", "Raw Time (s)")
  
  print(MaxErrorPercent)
  print(MeanErrorPercent)
  par(mfrow = c(1,1))
  
  return(df_awght)
}

for (i in c(4.005, 4.05, 4.1, 4.2, 4.5, 5)){
  print(i)
  df_awght <- proveItWorksBuf(ORIGINALGRID = 21, BIGGRID = 750, NC = 10, NC.buffer = 10, AWGHT = i)
}

print(df_awght)
save(df_awght, file = "timing_error/dataframes/awght_error.rda")
###########################################################################
###########################################################################