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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot46/"
	savePngVerbose( paste(path, "Plot46.png", sep=""), plotGaussErrorSurface, w=dw, h=dw, un="px" )
}

### 3D plot error surface for range of Amplitude and Sigma values
plotGaussErrorSurface = function( xyAngle=0, vertAngle=15 ) {
	set.seed( 3 )
	a = 5
	s = 2
	x = seq( -10, 10, length=100 )
	y = 1 + gauss1dFunc( x, a, s ) + rnorm( length(x), sd=0.2 )
	minAmpl = 0
	maxAmpl = 10
	minSigma = 0.1
	maxSigma = 6
	nr = 11
	nc = 11
	ampls = seq( minAmpl, maxAmpl, length.out=nc )
	sigmas = seq( minSigma, maxSigma, length.out=nr )
	z = getGaussErrorGrid( meanGaussSsdErr, ampls, sigmas, x, y )
	xLim1 = floor( minAmpl )
	xLim2 = ceiling( maxAmpl )
	yLim1 = floor( minSigma )
	yLim2 = ceiling( maxSigma )
	xform = persp( ampls, sigmas, t(z), xlim=c(xLim1,xLim2), ylim=c(yLim1,yLim2), cex.lab=lf,
				   main="1D Gaussian Cost Function\nMean Fit Error (z) vs. Amplitude & Sigma (x,y)",
				   xlab="Amplitude", ylab="Sigma", zlab="\nMean Error",
				   cex.lab=1.2,
				   theta=xyAngle, phi=vertAngle, expand = 0.8, shade=0.5,
				   col = "lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=1.3 )
}

### 1D Gaussian
gauss1dFunc = function( x, a, sigma ) {
	return( a * exp( -0.5 * (x/sigma)^2 ) )
}

### Gaussian error/cost function - SSD
meanGaussSsdErr = function( x, y, a, sigma ) {
	return( (1/length(x)) * sum( (gauss1dFunc( x, a, sigma ) - y)^2 ) )
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
