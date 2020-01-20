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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot70/"
	savePngVerbose( paste(path, "Plot70.png", sep=""), plotLineCostFuncSurface, w=dw, h=dh, un="px" )		
}

### Plot line cost function 3D surface
plotLineCostFuncSurface = function() {
	set.seed( 3 )
	trueB = 2
	trueM = 3
	coeff = c( trueB, trueM )
	x = seq( 0, 10, length=100 )
	y = evalPoly( x, coeff ) + rnorm( length(x), sd=2 )
	nc=20
	nr=20
	slopes = seq( trueM - 0.5, trueM + 0.5, length=nr )
	intercepts = seq( trueB - 2, trueB + 2, length=nc )
	z = getLineErrorGrid( meanLineSsdErr, slopes, intercepts, x, y )
	
	parSave = par( mar=c( 3, 6, 5, 4 ) )
	
	plotTitle = sprintf( "Line Cost Function\nTrue Slope: %.1f  True y-Intercept: %.1f", trueM, trueB )
	persp( intercepts, slopes, t(z), theta=40, phi=10, zlab="Mean SSD Error", xlab="y-Intercept", ylab="Slope", main=plotTitle,
		   shade=0.2, col="lightblue", expand=0.5, cex.lab=1.3, cex.main=1.3*fontScale, cex.axis=1.3, ticktype="detailed" )
	
	par( parSave )
}

### Compute line fit cost function value over grid of slopes and intercepts
### errorFunc = the specific cost function
### Note: returned matrix may need to be transposed for some functions such as contour()
getLineErrorGrid = function( errorFunc, slopes, intercepts, x, y ) {
	nr = length( slopes )
	nc = length( intercepts )
	z = matrix( nrow=nr, ncol=nc )	
	for ( row in 1:nr ) {     # loop over y-coordinates (slopes)
		slope = slopes[row]
		for ( col in 1:nc ) { # loop over x-coordinates (y-intercepts)
			yInt = intercepts[col]
			z[row, col] = errorFunc( x, y, slope=slope, yInt=yInt )
		}
	}
	return( z )
}

### Line fit error function
meanLineSsdErr = function( x, y, slope, yInt ) {
	return( (1/length(x)) * sum( ( yInt + slope*x - y )^2 ) )
}

### Compute f(x) for a polynomial and a vector of x-coordinates
evalPoly = function( x, coeff ) {
	if ( length( coeff ) < 1 ) return( c(0) )  	
	termSum = 0
	for ( i in 1:length(coeff) ) {
		termSum = termSum + coeff[i] * x^(i-1)
	}
	return( termSum )
}

### Add copyright notice to plot via text()
addCopyright = function() {
	mtext( "Copyright \uA9 2019 Richard Creamer - All Rights Reserved", side=4, line=0, adj=0, cex=1.1 )
	mtext( "Email: 2to32minus1@gmail.com", side=4, line=1, adj=0, cex=1.1 )
}

### Save to PNG file, specify width and height
savePngVerbose = function( path, plotFunc, w=512, h=512, un="px", doCopyright=TRUE, ... ) {
	png( filename = path, type="cairo", units=un, width=w, height=h, pointsize=12, res=96 )
	plotFunc( ... )
	if ( doCopyright )
		addCopyright()
	dev.off()
}

runFunc() 
