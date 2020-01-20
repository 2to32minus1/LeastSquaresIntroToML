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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot38/"
	savePngVerbose( paste(path, "Plot38.png", sep=""), plotCubicWithDeriv, w=1.2*600, h=1.2*900, un="px" )	
}

### Plot 5th order polynomial - hard coded example/function
plotCubicWithDeriv = function() {
	x = c( 1, 5, 6, 10 )
	y = c( 8, 3, 5, 1 )
	clrs = list( cubicClr="blue", derivClr="red", vLineClr="chartreuse4", tanClr="tan" )
	legText = c( "Cubic", "Cubic's Derivative", "Derivative Zeros", "Cubic Tangents at Derivative's Zeros" )
	legClrs = c( clrs$cubicClr, clrs$derivClr, clrs$vLineClr, clrs$tanClr )
	plotTitle = "Cubic Polynomial and its Derivative"
	plot( x, y, main=plotTitle, axes = FALSE, asp=1, cex.lab=lf, xlab="", ylab="", cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale, type="n" )	
	addAxisLabels()
	box(); grid()
	axis( side = 1, at = 1:10, cex.axis=1.3 )
	axis( side = 2, at = -3:11, las=2, cex.axis=1.3 )
	lmFit = lm( y ~ x + I(x^2) + I(x^3) )		
	predX = seq( x[1] - 1, x[length(x)] + 1, length.out=400 )
	predY = evalPoly( predX, lmFit$coefficients )	
	lines( predX, predY, col=clrs$cubicClr, lwd=lw+1 )
	cubicCoeff = lmFit$coefficients
	derivCoeff = getPolyDerivCoeff( cubicCoeff )
	derivPredX = seq( 1, 10, length.out = 400 )
	derivPredY = evalPoly( derivPredX, derivCoeff )
	lines( derivPredX, derivPredY, col=clrs$derivClr, lwd=lw+1 )
	derivRoots = qe( derivCoeff )
	cubicZeroSlopeX = c( derivRoots )
	cubicZeroSlopeY = evalPoly( cubicZeroSlopeX, cubicCoeff )
	plotTangents( cubicZeroSlopeX, cubicCoeff, lineLength=0.7, tanOffset=0.3, lwd=lw+1, col=clrs$tanClr )
	lines( c( gpl()$xLim[1], gpl()$xLim[2] ), c(0,0), col="black", lwd=lw+1 ) # horiz x-axis line
	x1 = derivRoots[1]
	x2 = derivRoots[2]
	y1 = evalPoly( x1, cubicCoeff )
	y2 = evalPoly( x2, cubicCoeff )
	lines( c( x1, x1 ), c( y1, -3 ), type="l", lty=2, col=clrs$vLineClr, lwd=lw+1 )
	lines( c( x2, x2 ), c( y2, -3 ), type="l", lty=2, col=clrs$vLineClr, lwd=lw+1 )
	legend( x="topright", y="top", legend=legText, col=legClrs, lty=c(1,1,2,1), lwd=lw+1, cex=1.3*fontScale, inset = c(0.06,0.06) )
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

### Rotate vector CCW
rotVecCcw = function( xComp, yComp, deg ) {
	theta = deg/180 * pi
	rotXComp = cos(theta) * xComp + sin(theta) * yComp;
	rotYComp = -sin(theta) * xComp + cos(theta) * yComp;
	return( c(rotXComp, rotYComp) )
}

### Plot one tangent line segment
plotTangentLine = function( x, y, lineCoeff, lineLength=1, tanOffset=0.2, col="red", lwd=1 ) {
	endPts = getTanLineEndPts( x, y, lineCoeff, lineLength=lineLength, tanOffset=tanOffset )
	lines( c(endPts$fromX, endPts$toX), c(endPts$fromY, endPts$toY), col=col, lwd=lwd )
}

### Naive implementation of quadratic equation solver: ax^2 + bx + c = 0
### Only special cases handled: 1) radicand < 0, and 2) only 1 root
qe = function( coeff ) {
	C = coeff[1]
	B = coeff[2]
	A = coeff[3]
	D = B^2 - 4*A*C
	if ( D < 0 ) {
		return( c() )
	} else {
		r1 = ( -B + sqrt(D) )/(2*A) # (A == 0) ? boom!
		r2 = ( -B - sqrt(D) )/(2*A)
		if ( r1 == r2 )
			return( c(r1) )
		return( c(r1, r2) )
	}
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

### Get Plot Limits : gpl()$xLim[1] --> left x-coord
gpl = function() {
	u = par( "usr" )
	return( list( xLim=u[1:2], yLim=u[3:4] ) )
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
