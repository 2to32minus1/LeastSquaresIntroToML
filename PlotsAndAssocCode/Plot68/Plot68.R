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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot68/"
	savePngVerbose( paste(path, "Plot68.png", sep=""), plotInvParabolaWithTangents, w=dh, h=dh, un="px" )
}

### Draw inverted parabola with tangent lines
plotInvParabolaWithTangents = function( lineLength=1 ) {
	coeff = c( 3, 0, -1 )
	x = c( -2, 0, 2 )
	y = evalPoly( x, coeff )
	xLimits = axisLimits( x, xMargin )
	yLimits = axisLimits( y, yMargin )
	polyDerivCoeff = getPolyDerivCoeff( coeff )
	plotTitle = expression( paste( "y = ", -x^2 + 3, ",  y' = -2x" ) )
	plot( x, y, main=plotTitle, axes=FALSE, xlim = xLimits, ylim = yLimits, asp=1, cex.lab=lf,
		  xlab="", ylab="", cex=lrgPts, lwd=lw+1, cex.main=1.4*fontScale )
	addAxisLabels()	
	axis( side = 1, at = xLimits[1]:xLimits[2], cex.axis=1.3 )
	axis( side = 2, at = yLimits[1]:yLimits[2], cex.axis=1.3 )
	box()
	grid()
	xc = ( xLimits[1] + xLimits[2] )/2
	text( -3.2, 4, "Tangent Slopes: 4, 0, -4", cex=1.4*fontScale, adj=0 )
	xPred = seq( x[1] - 1, x[length(x)] + 1, length.out=200 )
	yPred = evalPoly( xPred, coeff )
	lines( xPred, yPred, col="blue", lwd=lw )
	plotTangents( x, coeff, lineLength=lineLength, tanOffset=0.2 )
	text( -2.5, -1, "m = 4", cex=1.4*fontScale, adj=1 )
	text( -0.4, 3.4, "m = 0", cex=1.4*fontScale, adj=0 )
	text( 2.5, -1, "m = -4", cex=1.4*fontScale,adj=0 )
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

### Compute polynomial coefficients for derivative of argument polynomial
getPolyDerivCoeff = function( coeff ) {
	if ( length(coeff) < 2 ) return( c(0) )
	derivCoeff = c()
	for ( i in 2:length( coeff ) ) {
		termPower = i - 1
		derivCoeff[ length( derivCoeff ) + 1 ] = termPower * coeff[i]
	}
	return( derivCoeff )
}

### Plot tangent(s) to polynomial at x[i], y[i]
### Will not look correct unless both x and y axes are same scale (1:1 aspect ratio)
plotTangents = function( x, coeff, lineLength=1, tanOffset=0.2, col="red", lwd=2 ) {
	if ( length(coeff) < 2 ) return()  	
	y = evalPoly( x, coeff )
	polyDerivCoeff = getPolyDerivCoeff( coeff )	
	tanLineSlopes = evalPoly( x, polyDerivCoeff )
	tanLineYInts = y - tanLineSlopes * x
	for ( i in 1:length(x) ) {
		lineCoeff = c( tanLineYInts[i], tanLineSlopes[i] )
		plotTangentLine( x[i], y[i], lineCoeff, lineLength=lineLength, tanOffset=tanOffset, col=col, lwd=lwd )
	}
}

### Get tangentl line ends point about center point (x,y)
getTanLineEndPts = function( x, y, lineCoeff, lineLength=1.0, tanOffset=0.2 ) {
	x1 = x - 1
	x2 = x + 1
	y1 = evalPoly( x1, lineCoeff )
	y2 = evalPoly( x2, lineCoeff )
	xComp = x2 - x1
	yComp = y2 - y1
	vMag = sqrt( xComp^2 + yComp^2 )
	vUnit = c( xComp/vMag, yComp/vMag )
	vOffset = rotVecCcw( vUnit[1], vUnit[2], 90 )
	ox = vOffset[1] * tanOffset
	oy = vOffset[2] * tanOffset
	xComp = xComp/vMag * lineLength
	yComp = yComp/vMag * lineLength
	fromX = x - xComp/2 - ox
	fromY = y - yComp/2 - oy
	toX = x + xComp/2 - ox
	toY = y + yComp/2 - oy	
	return( list(fromX=fromX, fromY=fromY, toX=toX, toY=toY) )
}

### Plot one tangent line segment
plotTangentLine = function( x, y, lineCoeff, lineLength=1, tanOffset=0.2, col="red", lwd=1 ) {
	endPts = getTanLineEndPts( x, y, lineCoeff, lineLength=lineLength, tanOffset=tanOffset )
	lines( c(endPts$fromX, endPts$toX), c(endPts$fromY, endPts$toY), col=col, lwd=lwd )
}

### Rotate vector CCW
rotVecCcw = function( xComp, yComp, deg ) {
	theta = deg/180 * pi
	rotXComp = cos(theta) * xComp + sin(theta) * yComp;
	rotYComp = -sin(theta) * xComp + cos(theta) * yComp;
	return( c(rotXComp, rotYComp) )
}

### Compute plot axis limits to fit range of data
axisLimits = function( v, margin=0 ) {
	lowerLimit = floor( min(v) ) - margin
	upperLimit = ceiling( max(v) ) + margin
	c( lowerLimit, upperLimit )
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
