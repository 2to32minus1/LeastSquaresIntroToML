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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot36/"
	savePngVerbose( paste(path, "Plot36.png", sep=""), plot15PointLineBruteForceFit, w=dw, h=dw, un="px" )
}

### Main purpose of this function is to write a CSV file with grid of SSD fit error vs. m and b
plot15PointLineBruteForceFit = function() {
	csvPath = "D:/Code/R/LeastSquaresPlay/BruteForceLineErrorGrid.csv"
	set.seed( 3 )	
	trueM = 1.3 # slope
	trueB = 2.7 # y-intercept
	coeff = c( trueB, trueM )
	nPts = 15
	sigma = 1.0
	x = 1:nPts
	y = evalPoly( x, coeff ) + rnorm( length(x), sd=sigma )	
	xLimits = axisLimits( x )
	yLimits = axisLimits( y )	
	lmFit = lm( y ~ x )
	fitB = lmFit$coefficients[1]
	fitM = lmFit$coefficients[2]
	fitMeanSsdErr = (1/nPts) * sum( (fitM*x + fitB - y)^2 )
	plot( x, y, axes = FALSE, xlim = xLimits, ylim = yLimits,
		  cex.main = 1.3*fontScale, cex=lrgPts, lwd=lw+1,
		  xlab="", ylab="",
		  main=sprintf( "Fit Line to Linear Model Data (noise added)\nTrueM: %.1f  TrueB: %.1f", trueM, trueB ) )
	addAxisLabels()	
	grid()
	box()
	abline( lmFit, lwd=lw+1, col="blue" )
	axis( side = 1, at = xLimits[1]:xLimits[2], cex.axis=1.3 )
	axis( side = 2, at = yLimits[1]:yLimits[2], las=2, cex.axis=1.3 )	
	
	# Write data to csv file w/grid of SSD error values
	nc = 9
	nr = 9
	mVals = seq( trueM - 2, trueM + 2, length=nr )
	bVals = seq( trueB - 2, trueB + 2, length=nc )
	errVals = matrix( nrow=nr+1, ncol=nc+1, data=0 )
	for ( row in 2:(nr+1) ) {
		m = mVals[row-1]
		for ( col in 2:(nc+1) ) {
			b = bVals[col-1]
			errVals[row, col] = (1/nPts)*sum( (m*x + b - y)^2 )
		}
	}
	errVals[1,2:(nc+1)] = bVals
	errVals[2:(nr+1),1] = mVals
	write.table( errVals, file=csvPath, col.names=FALSE, row.names=FALSE, sep="," )
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

### Convience method to add x/y axis labels
addAxisLabels = function( xLabel="x", yLabel="y", cexVal=1.3 ) {
	mtext( text=xLabel, side=1, line=2.5, cex=cexVal )
	mtext( text=yLabel, side=2, line=2.5, cex=cexVal )
}

### Compute plot axis limits to fit range of data
axisLimits = function( v, margin=0 ) {
	lowerLimit = floor( min(v) ) - margin
	upperLimit = ceiling( max(v) ) + margin
	c( lowerLimit, upperLimit )
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
