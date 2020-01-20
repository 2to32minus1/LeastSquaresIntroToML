##############################################################
### Copyright (c) Richard Creamer 2019 - All Rights Reserved
### Inquiries: email 2to32minus1@gmail.com
##############################################################

### Globals
fontScale = 1.2    # Global font size variable
lf        = 1.2    # Enlarge font of plot labels
verbose   = FALSE  # Verbose printing
lw        = 2      # Line width
xMargin   = 1      # X-axis margin for plots
yMargin   = 1      # Y-axis margin for plots
de        = 0.01   # Small constant for estimating derivatives numerically
smlPts    = 1      # Small  data point/circle
medPts    = 1.5    # Medium data point/circle
lrgPts    = 2.0    # Large  data point/circle

### Run this function to generate this plot
runFunc = function() {
	dw = 1000
	dh = 750	
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot48/"
	savePngVerbose( paste(path, "Plot48.png", sep=""), plotGaussianCostFuncExt, w=1400, h=1400, un="px" )
}

### Plot extended Gaussian cost function without convergence trails
plotGaussianCostFuncExt = function( nLevels=250, drawCrossHairs=TRUE ) {
	set.seed( 3 )
	gd = makeGaussianTestData()
	trueS = gd$trueSigma
	trueA = gd$trueAmpl	
	maxA = trueA * 1.5
	minA = -maxA
	maxS = trueS * 3.5
	minS = -maxS
	gridAmpls  = seq( minA * 0.9, maxA * 0.9, length=400 )
	gridSigmas = seq( minS * 0.9, maxS * 0.9, length=400 )
	plotTitle = paste( "1D Gaussian Cost Function Contours - Extended Range\n", sprintf("Model: Amplitude = %.1f  Sigma = %.1f", trueA, trueS), sep="" )
	
	z = getGaussErrorGrid( meanGaussSsdErr, gridAmpls, gridSigmas, gd$x, gd$y )
	contour( gridAmpls, gridSigmas, t(z), nlevels=nLevels, xlab="", ylab="", main=plotTitle, cex.main=1.3*fontScale )
	addAxisLabels( xLabel="Amplitude", yLabel="Sigma", cexVal=1.7 )
	
	if ( drawCrossHairs ) {
		x1 = gpl()$xLim[1]
		x2 = gpl()$xLim[2]
		y1 = gpl()$yLim[1]
		y2 = gpl()$yLim[2]
		lines( c(x1, x2), c( trueS, trueS ), col="blue", lwd=lw )
		lines( c(trueA, trueA), c(y1, y2), col="blue", lwd=lw )
	}	
}	

### Shared use case Gaussian data set generator used by several cooperating functions
makeGaussianTestData = function( trueAmpl=5, trueSigma=2, aGuess=sqrt(trueAmpl), sGuess=sqrt(trueSigma), nPts=40 ) {
	set.seed( 3 )
	x = seq( -5*trueSigma, 5*trueSigma, length=nPts )
	y = gauss1dFunc( x, trueAmpl, trueSigma ) + rnorm( length(x), sd=0.2 )
	guesses = list( a=aGuess, s=sGuess )
	tuple = list( trueAmpl=trueAmpl, trueSigma=trueSigma, x=x, y=y, nPts=nPts, guesses=guesses )
	return( tuple )
}

### Compute mean SSD Gaussian cost function over a grid of Amplitude and Sigma values
### Note: returned matrix may need to be transposed for some functions such as contour()
### TODO: Try to use outer()
getGaussErrorGrid = function( errorFunc, ampls, sigmas, x, y ) {
	nr = length( sigmas )
	nc = length( ampls )
	z = matrix( nrow=nr, ncol=nc )	
	for ( row in 1:nr ) {     # loop over y-coordinates (sigma)
		sigma = sigmas[row]
		for ( col in 1:nc ) { # loop over x-coordinates (amplitude)
			ampl = ampls[col]
			z[row, col] = errorFunc( x, y, a=ampl, sigma=sigma )
		}
	}
	return( z )
}

### 1D Gaussian
gauss1dFunc = function( x, a, sigma ) {
	return( a * exp( -0.5 * (x/sigma)^2 ) )
}

### Gaussian error/cost function - SSD
meanGaussSsdErr = function( x, y, a, sigma ) {
	return( (1/length(x)) * sum( (gauss1dFunc( x, a, sigma ) - y)^2 ) )
}

### Get Plot Limits : gpl()$xLim[1] --> left x-coord
gpl = function() {
	u = par( "usr" )
	return( list( xLim=u[1:2], yLim=u[3:4] ) )
}

### Convience method to add x/y axis labels
addAxisLabels = function( xLabel="x", yLabel="y", cexVal=1.3 ) {
	mtext( text=xLabel, side=1, line=2.5, cex=cexVal )
	mtext( text=yLabel, side=2, line=2.5, cex=cexVal )
}

### Save to PNG file, specify width and height
savePngVerbose = function( path, plotFunc, w=512, h=512, un="px", doCopyright=TRUE, ... ) {
	png( filename = path, type="cairo", units=un, width=w, height=h, pointsize=12, res=96 )
	plotFunc( ... )
	if ( doCopyright )
		addCopyright()
	dev.off()
}

### Add copyright notice to plot via text()
addCopyright = function() {
	mtext( "Copyright \uA9 2019 Richard Creamer - All Rights Reserved", side=4, line=0, adj=0, cex=1.1 )
	mtext( "Email: 2to32minus1@gmail.com", side=4, line=1, adj=0, cex=1.1 )
}

runFunc()
