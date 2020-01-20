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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot2/"
	savePngVerbose( paste(path, "Plot2.png", sep=""), plotGaussErrorSurfaceStereoPairGraphic, w=1200, h=400, un="px", doCopyright=FALSE )	
}

### Cross-eyes pair of Gaussian cost function surface
plotGaussErrorSurfaceStereoPairGraphic = function() {
	set.seed( 3 )
	gd = makeGaussianTestData( trueAmpl=5, trueSigma=2 )
	a = gd$trueAmpl
	s = gd$trueSigma
	x = gd$x
	y = gd$y
	xSpan = 44
	ySpan = 24
	ctrAmpl  = 0
	ctrSigma = 0
	minAmpl  = ctrAmpl  - xSpan/2
	maxAmpl  = ctrAmpl  + xSpan/2
	minSigma = ctrSigma - ySpan/2
	maxSigma = ctrSigma + ySpan/2
	aspect = xSpan/ySpan
	nr = 20
	nc = as.integer( nr *xSpan / ySpan )
	ampls = seq( minAmpl, maxAmpl, length.out=nc )
	sigmas = seq( minSigma, maxSigma, length.out=nr )
	z = getGaussErrorGrid( meanGaussSsdErr, ampls, sigmas, x, y )
	
	m = matrix( c( 1, 2 ), byrow=TRUE, nrow=1, ncol=2 )
	layout( m, widths = c( 1, 1 ), heights = c(1 ) )
	
	xyAngle=-1.5
	zAngle=20
	
	parSave = par( mar=c(0,1,0,1), fg="white", col.main="white", bg = rgb(31/255, 78/255, 121/255) )
	
	persp( ampls, sigmas, t(z), xlab="", ylab="", zlab="", main="", cex.lab=1.3*fontScale, axes=FALSE,
		   theta=xyAngle+3, phi=zAngle, expand=0.5, shade=0.3, col="lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=fontScale )
	plotTitle = sprintf( "1D Gaussian Least-Squares Cost Function\nSSD Fit Error vs. Amplitude & Sigma\nLeft Eye", s, a )
	persp( ampls, sigmas, t(z), xlab="", ylab="", zlab="", main="",  cex.lab=1.3*fontScale, axes=FALSE,
		   theta=xyAngle, phi=zAngle, expand=0.5, shade=0.3, col="lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=fontScale )
	
	par( parSave )
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

### 1D Gaussian
gauss1dFunc = function( x, a, sigma ) {
	return( a * exp( -0.5 * (x/sigma)^2 ) )
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
