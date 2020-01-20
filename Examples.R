##############################################################
### Copyright (c) Richard Creamer 2019 - All Rights Reserved
### Inquiries: email 2to32minus1@gmail.com
##############################################################

source( "d:/code/R/LeastSquaresPlay/Util.R" )

#######################
# Top-level / exported functions below
#######################

### Draw cubic polynomial with tangents
plotCubicWithTangents = function( lineLength=1 ) {
	x = c( 1, 5, 6, 4 )
	y = c( 4, 2, 3, 1 )
	xLimits = axisLimits( x, xMargin )
	yLimits = axisLimits( y, yMargin )
	plot( x, y, axes=FALSE, asp=1, cex.lab=lf, xlab="", ylab="", cex=lrgPts, lwd=lw+1 )
	addAxisLabels()
	box(); grid()
	axis( side = 1, at = xLimits[1]:xLimits[2], cex.axis=1.3 )
	axis( side = 2, at = yLimits[1]:yLimits[2], las=2, cex.axis=1.3 )	
	lmFit = lm( y ~ x + I(x^2) + I(x^3) )
	xPred = seq( xLimits[1], xLimits[2], length.out=200 )
	yPred = evalPoly( xPred, lmFit$coefficients )
	lines( xPred, yPred, col="blue", lwd=lw )
	plotTangents( x, lmFit$coefficients, lineLength=0.3, tanOffset=0.1 )
	if ( verbose )
		print( lmFit$coefficients )
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

### Plot 5th order polynomial - hard coded example
plot5thOrderPoly = function() {
	x=c( 1, 2.8, 4.9, 6.4, 8, 10 )
	y=c( 1, 3.4, 2.3, 6.3, 5, 10 )	
	plot( x, y, axes = FALSE, asp=1, cex.lab=lf, xlab="", ylab="", cex=lrgPts, lwd=lw+1 )
	addAxisLabels()	
	box(); grid()
	axis( side = 1, at = 0:11, cex.axis=1.3 )
	axis( side = 2, at = 0:11, las=2, cex.axis=1.3 )	
	lmFit = lm( y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) )	
	predX = seq( x[1] - 1, x[length(x)] + 1, length.out=400 )
	predY = evalPoly( predX, lmFit$coefficients )	
	lines( predX, predY, col=4, lwd=lw )
	plotTangents( x, lmFit$coefficients, lineLength=0.5, tanOffset=0.2 )	
	return( c( x=x, y=y, lmFit=lmFit ) )
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

### Draw 2-point line -> exact fit
plot2PointLine = function() {
	set.seed( 3 )	
	trueM = 1.3
	trueB = 2.7
	coeff = c( trueB, trueM )
	nPts = 2	
	x = c( 1, nPts )
	y = evalPoly( x, coeff )
	xLimits = axisLimits( x )
	yLimits = axisLimits( y )	
	plot( x, y, axes = FALSE, xlim = xLimits, ylim = yLimits, cex.lab=lf, xlab="", ylab="",
		  main = "Fit Line to Linear Model Data: 2 points --> no/zero fit error", cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale )
	addAxisLabels()	
	axis( side = 1, at = xLimits[1]:xLimits[2], cex.axis=1.3 )
	axis( side = 2, at = yLimits[1]:yLimits[2], las=2, cex.axis=1.3 )
	box()
	grid()
	lmFit = lm( y ~ x )
	abline( lmFit, col = "blue", lwd=lw )
	e = compPolyFitError( x, y, lmFit$coefficients )
	addLegend( coeff, lmFit$coefficients, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError )
}

### Main purpose of this function is to write a CSV file with grid of SSD fit error vs. m and b
plot15PointLineBruteForceFit = function() {
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
#	print( sprintf( "fitB = %.3f  fitM = %.3f  fitError = %.2f", fitB, fitM, fitMeanSsdErr ) )

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
	write.table( errVals, file="D:/Code/R/LeastSquaresPlay/BruteForceLineErrorGrid.csv", col.names=FALSE, row.names=FALSE, sep="," )
}

### Draw 15-point quadratic + noise, fit quadratic
plot15PointPoly2DataFitPoly2 = function() {
	set.seed( 3 )	
	coeff = c( 1, 1.3, 2.2 )
	nPts = 15
	sigma = 20
	x = 1:nPts
	y = evalPoly( x, coeff ) + rnorm( length(x), sd=sigma )	
	xLimits = axisLimits( x )
	yLimits = axisLimits( y )	
	plot( x, y, axes = FALSE, xlim = xLimits, ylim = yLimits, cex.lab=lf, xlab="", ylab="",
		  main = "Fit 2nd Order Polynomial to Quadratic Model Data  (noise added)", cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale )	
	addAxisLabels()	
	grid()
	box()
	axis( side = 1, at = xLimits[1]:xLimits[2], cex.axis=1.3 )
	axis( side = 2, at = seq( yLimits[1], yLimits[2], 50 ), las=2, cex.axis=1.3 )
	
	lmFit = lm( y ~ x + I(x^2) )
	
	predX = seq( x[1],x[length(x)],length.out=500 )
	predY = evalPoly( predX, lmFit$coefficients )
	lines( predX, predY, col=4, lwd=lw )
	segments( x, y, x, predict( lmFit ), col = "red", lwd=lw )
	
	e = compPolyFitError( x, y, lmFit$coefficients )
	addLegend( coeff, lmFit$coefficients, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError )
}

### Fit Nth order polynomial to linear data
plotNthOrderPolyToLinearData = function( nPts=15, N=1, modelM=1.3, modelB=2.7, noiseSigma=1.5, doPlot=FALSE ) {
	set.seed( 3 )
	
	if ( N < 1 || N > 11 ) { print( sprintf( "Error: N (%d) outside range 1:11", N) ); return() }
	if ( nPts < (N + 1) )  { print( sprintf( "Error: nPts (%d) must be > N (%d)", nPts, N ) ); return() }
	
	coeff = c( modelB, modelM )
	x = 1:nPts
	y = evalPoly( x, coeff ) + rnorm( length(x), sd=noiseSigma )

	# Build appropriate formula, example: "y ~ x + I(x^2) + I(x^3) + ..."
	polyFormula = makePolyFormula( polyOrder=N )
	
	lmFit = lm( as.formula(polyFormula) )
	e = compPolyFitError( x, y, lmFit$coefficients )
	if ( doPlot ) {
		plotTitle = sprintf( "Fit Polynomial of Order %d to Linear Data\nNoise Added", N );
		plot( x, y, cex.lab=lf, xlab="", ylab="", main=plotTitle, cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale )
		addAxisLabels()
		grid()
		box()
		leftX = gpl()$xLim[1]
		rightX = gpl()$xLim[2]
		predX = seq( leftX, rightX, length.out=200 )
		predY = evalPoly( predX, lmFit$coefficients )
		lines( predX, predY, col=4, lwd=lw )
		segments( x, y, x, predict(lmFit), col = "red", lwd=lw )  # Note: predict(lmFit) same as evalPoly(x, lmFit$coefficients)
		if ( N > 9 ) {
			addLegend( coeff, lmFit$coefficients, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError, wrapLimit=4 )
		}
		else {
			addLegend( coeff, lmFit$coefficients, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError )
		}
	}
	return( list(nPts=nPts, x=x, y=y, lmFit=lmFit, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError) )
}

### Simulate temperature vs. time
tempVsTime = function() {
	set.seed( 3 )
	# Make some fake data
	trueCoeff = c( 44, 6, -0.23 )
	x = 1:10
	y = evalPoly( x, trueCoeff ) + rnorm( length(x), sd=1.5 )
	# Plot it
	plot( x, y, main="Temperature vs. Minute", cex.lab=lf, xlab="", ylab="", cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale, cex.axis=1.3 )
	addAxisLabels( xLabel="Minute", yLabel="Temperature" )	
	# Fit a quadratic curve to simulated real-world data
	lmFit = lm( y ~ x + I(x^2) )
	
	predX = seq( x[1],x[length(x)],length.out=200 )
	predY = evalPoly( predX, lmFit$coefficients )
	lines( predX, predY, col=4, lwd=lw )
	
	e = compPolyFitError( x, y, lmFit$coefficients )
	addLegend( trueCoeff, lmFit$coefficients, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError )
}

### Play with non-linear curve fitting using nlm() built-in optimizer
fitGaussNls = function() {
	set.seed( 3 )
	trueMu = 57.5
	trueSigma = 2
	truePeakHt = 180
	tuple = getSimGammaData( mu=trueMu, sigma=trueSigma, peakHt=truePeakHt )
	x = tuple$x
	y = tuple$y
	data = data.frame( x, y )
	guesses = list( c0=150, c1=-1, c2=0.3, peakHt=200, sigma=3, mu=55 )
	f = y ~ c0 + c1*x + c2*x^2 + gauss(x, mu, sigma, peakHt)
	model= nls( f, data, start=guesses )
	coeff = coef( model )
	c0 = coeff["c0"]
	c1 = coeff["c1"]
	c2 = coeff["c2"]
	peakHt = coeff["peakHt"]
	sigma = coeff["sigma"]
	mu = coeff["mu"]
	xPred = seq( x[1], x[length(x)], length.out=500 )
	yPred = evalGammaFunc( xPred, c0, c1, c2, mu, sigma, peakHt )
	xLab = "Photon Energy"
	yLab = "Photon Count"
	plotTitle = paste( "Simulated Gamma Spectroscopy: Quadratic Background + Gaussian Peak" )
	plot( x, y, ylim=c(0,300), main=plotTitle, xlab="", ylab="", cex.lab=lf, cex=medPts, lwd=lw+1, cex.axis=1.3, cex.main=1.3*fontScale )
	addAxisLabels( xLabel=xLab, yLabel=yLab )	
	lines( xPred, yPred, col="blue",lwd=lw )
	yGauss = gauss( xPred, mu, sigma, peakHt )
	lines( xPred, yGauss, col="red", lwd=lw )
	text( 0, 10, "Extracted Gaussian", pos=4, cex=fontScale, pch=15 )
	legLine1 = paste( "y ==", quadPlusGaussianMathText(c0,c1,c2,peakHt,sigma, mu), sep="" )
	legLine2 = paste( "y ==", getGaussianMathText( peakHt, sigma, mu ), sep="" )
	legLines = c( parse(text=legLine1), parse(text=legLine2) )
	legCols = c( "blue", "red"  )
	legend( "topleft", legend=legLines, col=legCols, lwd=lw, cex=1.4*fontScale, inset=c(.0,.0),box.col="transparent" )
	textLine0 = "Gradient Descent Summary (nls):"
	textLine1 = sprintf( "True Gaussian Model:      Amplitude=%.2f Sigma=%.2f Mu=%.2f", truePeakHt, trueSigma, trueMu )
	textLine2 = sprintf( "Init. Gaussian Guesses:   Amplitude=%.2f Sigma=%.2f Mu=%.2f", guesses$peakHt, guesses$sigma, guesses$mu )
	textLine3 = sprintf( "Fitted Gaussian Model:    Amplitude=%.2f Sigma=%.2f Mu=%.2f", peakHt, sigma, mu )
	textLines = paste( textLine0, "\n", textLine1, "\n", textLine2, "\n", textLine3, sep="" )
	text( 0, 55, textLines, cex=fontScale,adj=0 )
	text( 75, 20, "Scientists need this", adj=0, cex=fontScale )
	mathText = "Model:  y == c[0] + c[1]*x + c[2]*x^2 + a*e^frac(-(x-mu)^2,2*sigma^2)"
	text( 61, 190, parse(text=mathText), cex=1.6*fontScale, adj=0 )
	arrows( 74, 20, 63, 20,col="blue",lwd=lw, angle=15 )
}

### Plot inverted parabola with a tangent line at vertex
plotInvParabolaWithVertexTangent = function() {
	f = function( x ) { return( -(x-3)*(x-3) + 5 ) }
	x = -1:7
	y = f(x)
	xLimits = c( min(x), max(x) )
	yLimits = c( min(y), max(y) + 6 )
	plot( x, y, axes=FALSE, ylim=yLimits, xlim=xLimits, main="Parabola + Tangent at Vertex", cex.lab=lf,
		  xlab="", ylab="", cex=lrgPts, lwd=lw+1, cex.main=1.4*fontScale, cex.axis=1.3 )
	addAxisLabels()	

	axis( side=1, at=xLimits[1]:xLimits[2], cex.axis=1.3 )
	axis( side=2, at=yLimits[1]:yLimits[2], las=2, cex.axis=1.3 )
	box()
	xVertex = 3
	yVertex = f( xVertex )
	xs = seq( x[1], x[length(x)], length.out=500 )
	ys = f(xs)
	lines( xs, ys, col=4, lwd=lw )
	lines( c( xVertex - 3, xVertex + 3), c( yVertex, yVertex ), col="red", lwd=lw+1 )
	lines( c(3,3), c(gpl()$yLim[1],f(3)), col="chartreuse4", lwd=lw+1, lty=2 )
	eqText = expression( paste( "y = ", -x^2 + 6*x - 4, ",  y' = -2x + 6,  y'(3) = 0" ) )
	text( 3, 11, eqText, cex=1.4*fontScale )
	text( 3, 9, "x-coord of maxima is where deriv and tangent slope = 0", cex=1.4*fontScale )
}

### Draw inverted gaussian surface with a few gradient vectors with proportional magnitudes
gradDemo = function( zFunc=elevFunc, xyAngle=30, vertAngle=45, sigma=3, peakHt=10, col="white", lwd=lw, expand=0.2 ) {
	peakHt = 0.5
	x = seq( -10, 10, length=21 )
	y = x
	z = outer( x, y, elevFunc, sigma, peakHt )
	xform = persp( x, y, t(z), xlim=c(-10,10), ylim=c(-10,10), cex.lab=1.3, zlab="Elevation",
				  main="Gaussian Gradients - Direction + Magnitude",
				  theta = xyAngle, phi = vertAngle, expand = 0.3, shade=0.4,
				  col = "lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=1.3 )
	arrowLen = 80
	x1 = c( -4, 0, 4, -4,  0,  4, 7, -7, -5.8, 5.8, -5.8,  5.8 )
	y1 = c( 4,  7, 4, -4, -7, -4, 0,  0,  5.8, 5.8, -5.8, -5.8 )
	z1 = elevFunc( x1, y1, sigma, peakHt )
	gradX = -gaussDzDx( x1, y1, sigma, peakHt )
	gradY = -gaussDzDy( x1, y1, sigma, peakHt )
	mags = sqrt( gradX^2 + gradY^2 )
	gradX = arrowLen * gradX
	gradY = arrowLen * gradY
	x2 = x1 + gradX
	y2 = y1 + gradY
	z2 = z1	
	pt1 = trans3d( x1, y1, z1, xform )
	pt2 = trans3d( x2, y2, z2, xform )
	arrows( pt1$x, pt1$y, pt2$x, pt2$y, col=col, lwd=lwd, angle=15, length = 0.125 )
}

### Draw inverted gaussian surface with a few gradient vectors with constant magnitudes (dir only)
gradDemo2 = function( zFunc=elevFunc, xyAngle=30, vertAngle=50, sigma=3, peakHt=10, col="white", lwd=lw, expand=0.2 ) {
	peakHt = 0.5
	x = seq( -10, 10, length=21 )
	y = x
	z = outer( x, y, elevFunc, sigma, peakHt )
	xform = persp( x, y, t(z), xlim=c(-10,10), ylim=c(-10,10), cex.lab=1.3, zlab="Elevation",
				  main="Gaussian Gradients - Direction Only",
				  theta = xyAngle, phi = vertAngle, expand = 0.3, shade=0.4,
				  col = "lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=1.3 )
	arrowLen = 1.4
	x1 = c()
	y1 = c()
	byVal=3
	for ( i in seq(-9,9,by=byVal) ) {
		for ( j in seq(-9,9,by=byVal) ) {
			x1[length(x1) + 1] = i
			y1[length(y1) + 1] = j
		}
	}
	z1 = elevFunc( x1, y1, sigma, peakHt )
	gradX = -gaussDzDx( x1, y1, sigma, peakHt )
	gradY = -gaussDzDy( x1, y1, sigma, peakHt )
	mags = sqrt( gradX^2 + gradY^2 )
	gradX = arrowLen * gradX/mags
	gradY = arrowLen * gradY/mags
	x2 = x1 + gradX
	y2 = y1 + gradY
	z2 = z1	
	pt1 = trans3d( x1, y1, z1, xform )
	pt2 = trans3d( x2, y2, z2, xform )
	arrows( pt1$x, pt1$y, pt2$x, pt2$y, col=col, lwd=lwd, angle=15, length = 0.125 )
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

### 2D plot error contours for range of Amplitude and Sigma values
### Amplitude maps to x-axis, Sigma maps to y-axis
plotGaussErrorContours = function( nLevels=30, customTitle="", drawGradVecs=FALSE, normGradVecs=TRUE, gd=NULL ) {
	set.seed( 3 )
	if ( is.null( gd ) )
		gd = makeGaussianTestData()
	a = gd$trueAmpl
	s = gd$trueSigma
	x = gd$x
	y = gd$y
	minAmpl  = a * 0.1;	maxAmpl  = a * 1.7
	minSigma = s * 0.1;	maxSigma = s * 1.7
	# Adjust ranges so plot is square
	deltaA = maxAmpl - minAmpl
	deltaS = maxSigma - minSigma
	maxDelta = if ( deltaA > deltaS ) deltaA else deltaS
	if ( deltaA > deltaS ) {
		maxSigma = minSigma + maxDelta
	} else {
		maxAmpl = minAmpl + maxDelta
	}
	nr = 200
	nc = 200
	ampls = seq( minAmpl, maxAmpl, length.out=nc )
	sigmas = seq( minSigma, maxSigma, length.out=nr )
	z = getGaussErrorGrid( meanGaussSsdErr, ampls, sigmas, x, y )
	xLim1 = minAmpl
	xLim2 = maxAmpl
	yLim1 = minSigma
	yLim2 = maxSigma
	plotTitle = ""
	if ( !drawGradVecs ) {
		plotTitle = sprintf( "1D Gaussian Cost Function Contours\nTrue Model: Amplitude = %.1f  Sigma = %.1f (noise added)", a, s )
	} else if ( normGradVecs ) {
		plotTitle = sprintf( "1D Gaussian Cost Function Contours w/Gradient (Unit) Vectors\nTrue Model: Amplitude = %.1f  Sigma = %.1f (noise added)", a, s )
	} else {
		plotTitle = sprintf( "1D Gaussian Cost Function Contours w/Gradient Vectors\nTrue Model: Amplitude = %.1f  Sigma = %.1f (noise added)", a, s )
	}
	if ( customTitle != "" )
		plotTitle = customTitle
	xform = contour( ampls, sigmas, t(z), axes=FALSE, cex.lab=lf, main=plotTitle, asp=1.0,
					xlab="", ylab="", labcex=1.2, nlevels=nLevels, cex.main=1.3*fontScale, cex.axis=1.3 )
	xAxisMin = ceiling( gpl()$xLim[1] )
	xAxisMax = floor( gpl()$xLim[2] )
	yAxisMin = ceiling( gpl()$yLim[1] )
	yAxisMax = floor( gpl()$yLim[2] )
	axis( side=1, at=xAxisMin:xAxisMax, cex.axis=1.3 )
	axis( side=2, at=yAxisMin:yAxisMax, cex.axis=1.3 )
	addAxisLabels( xLabel="Amplitude", yLabel="Sigma" )
	box()
	xLimits = axisLimits( x )
	yLimits = axisLimits( y )
	lines( c(-1,12), c( s, s ), col="blue", lwd=lw )
	lines( c( a, a), c(-1,10), col="blue", lwd=lw )
	if ( drawGradVecs ) {
		gradPtsX1 = c()
		gradPtsY1 = c()
		gradPtsX2 = c()
		gradPtsY2 = c()
		gridSize=13
		for ( sigma in seq( from=sigmas[2], to=sigmas[length(sigmas)-1], length=gridSize ) ) {
			for ( ampl in seq( from=ampls[2], to=ampls[length(ampls)-1], length=gridSize ) ) {
				gradPtsX1[length(gradPtsX1)+1] = ampl
				gradPtsY1[length(gradPtsY1)+1] = sigma
				da = gaussDeDa( x, y, ampl, sigma )
				ds = gaussDeDs( x, y, ampl, sigma )
				mag = sqrt( da*da + ds*ds )
				if ( normGradVecs ) {
					da = da/mag/4 # hard-coded case-specific scale factor, unit vecs a bit too long for plot
					ds = ds/mag/4
				} else {
					da = da * 0.22 # hard-coded case-specific scale factor
					ds = ds * 0.22
				}
				gradPtsX2[length(gradPtsX2)+1] = ampl  - da
				gradPtsY2[length(gradPtsY2)+1] = sigma - ds
			}
		}
		points( gradPtsX1, gradPtsY1, col="blue", lwd=lw )
		arrows( gradPtsX1, gradPtsY1, gradPtsX2, gradPtsY2, col="blue", lwd=lw, cex=medPts, angle=15, length=0.125 )
	}
}

### Cross-eyes pair of Gaussian cost function surface
plotGaussErrorSurfaceStereoPair = function() {
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

	parSave = par( mar=c( 1, 4, 5, 4 ) )
	xyAngle=-1.5
	zAngle=20
	plotTitle = sprintf( "1D Gaussian Least-Squares Cost Function\nSSD Fit Error vs. Amplitude & Sigma\nRight Eye", s, a )
	persp( ampls, sigmas, t(z), xlab="Amplitude", ylab="", zlab="", main=plotTitle, cex.lab=1.3*fontScale,
		   theta=xyAngle+3, phi=zAngle, expand=0.5, shade=0.3, col="lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=fontScale )
	plotTitle = sprintf( "1D Gaussian Least-Squares Cost Function\nSSD Fit Error vs. Amplitude & Sigma\nLeft Eye", s, a )
	persp( ampls, sigmas, t(z), xlab="Amplitude", ylab="", zlab="", main=plotTitle,  cex.lab=1.3*fontScale,
		   theta=xyAngle, phi=zAngle, expand=0.5, shade=0.3, col="lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=fontScale )
	par( parSave )
}

### 2D plot error contours for WIDE range of Amplitude and Sigma values + stereo pair of cost function surface
### Amplitude maps to x-axis, Sigma maps to y-axis
plotGaussErrorContoursStereo = function( nLevels=170, customTitle="", drawCrossHairs=TRUE ) {
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
	nr = 200
	nc = as.integer( nr *xSpan / ySpan )
	ampls = seq( minAmpl, maxAmpl, length.out=nc )
	sigmas = seq( minSigma, maxSigma, length.out=nr )
	z = getGaussErrorGrid( meanGaussSsdErr, ampls, sigmas, x, y )
	plotTitle = sprintf( "1D Gaussian Cost Function Mean SSD Error Contours\nTrue Model: Amplitude = %.1f  Sigma = %.1f (noise added)", a, s )
	if ( customTitle != "" )
		plotTitle = customTitle
	
	m = matrix( c( 1, 1, 2, 3 ), nrow=2, byrow=TRUE )
	layout( m, widths = c( 1, 1), heights = c(2, 1 ) )
	xform = contour( ampls, sigmas, t(z), axes=FALSE, cex.lab=lf, main=plotTitle, # asp=1.0,
					 xlab="", ylab="", labcex=1.3, nlevels=nLevels, cex.main=1.3*fontScale, cex.axis=1.3 )

	xTickVals = seq( minAmpl, maxAmpl, 2 )
	yTickVals = seq( minSigma, maxSigma, 2 )
	axis( side=1, at=xTickVals, cex.axis=1.3 )
	axis( side=2, at=yTickVals, cex.axis=1.3 )
	addAxisLabels( xLabel="Amplitude", yLabel="Sigma" )
	box()
	
	if ( drawCrossHairs ) {
		x1 = gpl()$xLim[1]
		x2 = gpl()$xLim[2]
		y1 = gpl()$yLim[1]
		y2 = gpl()$yLim[2]
		lines( c(x1, x2), c( s, s ), col="blue", lwd=lw )
		lines( c(a, a), c(y1, y2), col="blue", lwd=lw )
	}

	sampleFreq = 13
	perspX = ampls[ seq( 1, length(ampls), by=sampleFreq ) ]
	perspY = sigmas[ seq( 1, length(sigmas), by=sampleFreq ) ]
	perspZ = z[ seq( 1, length(sigmas), by=sampleFreq ), seq( 1, length(ampls), by=sampleFreq ) ]
	
	parSave = par( mar=c( 3.5, 2, 1, 2 ) )
	xyAngle=-1.5
	zAngle=20
	persp( perspX, perspY, t(perspZ), xlab="Amplitude", ylab="Sigma", zlab="", main="Right Eye: Cost Function Value = SSD Error", cex.lab=1.3*fontScale,
		   theta=xyAngle+3, phi=zAngle, expand=0.5, shade=0.3, col="lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=fontScale )
	persp( perspX, perspY, t(perspZ), xlab="Amplitude", ylab="Sigma", zlab="", main="Left Eye: Cost Function Value = SSD Error",  cex.lab=1.3*fontScale,
		   theta=xyAngle, phi=zAngle, expand=0.5, shade=0.3, col="lightblue", ticktype="detailed", cex.main=1.3*fontScale, cex.axis=fontScale )
	par( parSave )
}

### 2D plot error contours for WIDE range of Amplitude and Sigma values
### Amplitude maps to x-axis, Sigma maps to y-axis
plotGaussErrorContoursMono = function( nLevels=160, customTitle="", drawCrossHairs=TRUE ) {
	set.seed( 3 )
	gd = makeGaussianTestData( trueAmpl=5, trueSigma=2 )
	a = gd$trueAmpl
	s = gd$trueSigma
	x = gd$x
	y = gd$y
	xSpan = 36
	ySpan = 20
	ctrAmpl  = 0
	ctrSigma = 0
	minAmpl  = ctrAmpl  - xSpan/2
	maxAmpl  = ctrAmpl  + xSpan/2
	minSigma = ctrSigma - ySpan/2
	maxSigma = ctrSigma + ySpan/2
	aspect = xSpan/ySpan
	nr = 200
	nc = as.integer( nr *xSpan / ySpan )
	ampls = seq( minAmpl, maxAmpl, length.out=nc )
	sigmas = seq( minSigma, maxSigma, length.out=nr )
	z = getGaussErrorGrid( meanGaussSsdErr, ampls, sigmas, x, y )
	plotTitle = sprintf( "1D Gaussian Cost Function Mean SSD Error Contours\nTrue Model: Amplitude = %.1f  Sigma = %.1f (noise added)", a, s )
	if ( customTitle != "" )
		plotTitle = customTitle

	m = matrix( c( 1, 1, 1, 1, 2, 2 ), ncol=2, byrow=TRUE )
	layout( m, widths = c( 1, 1), heights = c(1, 1, 1.2 ) )
	
	xform = contour( ampls, sigmas, t(z), main=plotTitle, asp=1.0, axes=FALSE,
					 xlab="", ylab="", cex.lab=lf, labcex=1.3, cex.main=1.5*fontScale, cex.axis=1.3, nlevels=nLevels )

	tickFreq = 2
	xTickVals = seq( minAmpl, maxAmpl, by=tickFreq )
	yTickVals = seq( minSigma, maxSigma, by=tickFreq )
	axis( side=1, at=xTickVals, cex.axis=1.3*fontScale )
	axis( side=2, at=yTickVals, cex.axis=1.3*fontScale )
	addAxisLabels( xLabel="Amplitude", yLabel="Sigma" )
	
	box()
	
	if ( drawCrossHairs ) {
		x1 = gpl()$xLim[1]
		x2 = gpl()$xLim[2]
		y1 = gpl()$yLim[1]
		y2 = gpl()$yLim[2]
		lines( c(x1, x2), c( s, s ), col="blue", lwd=lw )
		lines( c( a, a), c(y1, y2), col="blue", lwd=lw )
	}
	
	sampleFreq = 7
	perspX = ampls[ seq( 1, length(ampls), by=sampleFreq ) ]
	perspY = sigmas[ seq( 1, length(sigmas), by=sampleFreq ) ]
	perspZ = z[ seq( 1, length(sigmas), by=sampleFreq ), seq( 1, length(ampls), by=sampleFreq ) ]
	
	parSave = par( mar=c( 5, 2, 1, 2 ) )
	
	xyAngle = 0
	zAngle = 20
	persp( perspX, perspY, t(perspZ), xlab="Amplitude", ylab="Sigma", zlab="", main="Cost Function Value = SSD Error",
		   cex.lab=1.5*fontScale, cex.main=1.5*fontScale, cex.axis=1.3*fontScale, 
		   theta=xyAngle, phi=zAngle, expand=0.5, shade=0.4, col="lightblue", ticktype="detailed" )

	par( parSave )
}

### Nonlinear least-squares gradient-descent curve fitting for Gaussian
fitGauss1d = function( learnRate = 0.5, nIter=100, nPts=40, drawPointLines=FALSE, doPlot=TRUE, numDerivs=FALSE, plotContours=FALSE, gd=NULL ) {
	set.seed( 3 )
	a = s = x = y = guesses = NULL
	if ( is.null( gd ) )
		gd = makeGaussianTestData()
	a = gd$trueAmpl
	s = gd$trueSigma
	trueAmpl = gd$trueAmpl
	trueSigma = gd$trueSigma
	x = gd$x
	y = gd$y
	nPts=gd$nPts
	guesses = gd$guesses
	aEst = guesses$a
	sEst = guesses$s
	da = 0
	ds = 0
	amplLegs = c( aEst )
	sigmaLegs = c( sEst )
	if ( nIter > 0 ) {
		# Basic implementation of Gradient Descent
		for ( i in 1:nIter ) {
			da = if ( numDerivs ) gaussDeDaNum( x, y, aEst, sEst ) else gaussDeDa( x, y, aEst, sEst )
			ds = if ( numDerivs ) gaussDeDsNum( x, y, aEst, sEst ) else gaussDeDs( x, y, aEst, sEst )
			aEst = aEst - learnRate * da
			sEst = sEst - learnRate * ds
			amplLegs[ length(amplLegs) + 1 ] = aEst
			sigmaLegs[ length(sigmaLegs) + 1 ] = sEst
			if ( verbose )
				print( sprintf( "Mean ABS fit error: %f\n", meanGaussAbsErr( x, y, aEst, sEst ) ) )
		}
	}
	finalMeanAbsFitError = meanGaussAbsErr( x, y, aEst, sEst )
	if ( doPlot) {
		plotTitle = "Gradient Descent Fit 1D Gaussian using Analytic Partial Derivatives"
		if ( numDerivs )
			plotTitle = "Gradient Descent Fit 1D Gaussian using Numerical Est. Partial Derivatives"
		plot( x, y, cex.lab=lf, xlab="", ylab="", main=plotTitle, cex=medPts, lwd=lw+1, cex.main=1.3*fontScale, cex.axis=1.3 )
		addAxisLabels()	
		predX = seq( x[1], x[length(x)], length.out=200 )
		predY = gauss1dFunc( predX, aEst, sEst )
		lines( predX, predY, col="blue", lwd=lw+1 )
		if ( drawPointLines )
			segments( x, y, x, gauss1dFunc( x, aEst, sEst ), col = "red" )
		legText = c()		
		legText[length(legText)+1] = sprintf( "Num Iterations   %d", nIter )
		legText[length(legText)+1] = sprintf( "Learn Rate:      %.3f", learnRate )
	       legText[length(legText)+1] = sprintf( "True Amplitude:  %.2f", a )
		legText[length(legText)+1] = sprintf( "True Sigma:      %.2f", s )
		legText[length(legText)+1] = sprintf( "Start Amplitude: %.2f", guesses$a )
		legText[length(legText)+1] = sprintf( "Start Sigma:     %.2f", guesses$s )
		legText[length(legText)+1] = sprintf( "Final Amplitude: %.2f", aEst )
		legText[length(legText)+1] = sprintf( "Final Sigma:     %.2f", sEst )
		legText[length(legText)+1] = sprintf( "Mean ABS Error:  %.2f", meanGaussAbsErr( x, y, aEst, sEst ) )
		legLtys = seq( from=0, to=0, length=length(legText) )
		parSave = par( family="mono", font=2 )		
		legend( "topleft", inset=c(0,0), legend=legText, lty=legLtys, box.col="transparent", cex=1.3*fontScale, bg="transparent" )
		par( parSave )
	}
	if ( plotContours ) {
		plotGaussErrorContours( gd=gd )
		lines( amplLegs, sigmaLegs, col="red", lwd=lw )
	}
	return( list( x=x, y=y, a=aEst, s=sEst, meanAbsFitError=finalMeanAbsFitError, amplLegs=amplLegs, sigmaLegs=sigmaLegs, lr=learnRate ) )
}

### Plot Gaussian fit step vectors for multiple starting points
plotGaussianParamStepVectors = function( nLevels=70 ) {
	set.seed( 3 )
	gd = makeGaussianTestData()
	trueS = gd$trueSigma
	trueA = gd$trueAmpl	
	a1 = trueA * 0.2
	a2 = trueA * 1.5
	s1 = trueS * 0.3
	s2 = trueS * 3.5
	startAmpls =  c( a1, a1, a2, a2 )
	startSigmas = c( s1, s2, s1, s2 )
	endAmpls = c()
	endSigmas = c()
	trailColors = c( "red", "blue", "chartreuse4", "purple1" )
	custTitle = paste( "1D Gaussian Cost Function Contours + (4) Gradient Descent Convergence Parameter Trails\n",
					   sprintf("True Model: Amplitude = %.1f  Sigma = %.1f  (noise added)", trueA, trueS ) )
	plotGaussErrorContours( nLevels=nLevels, customTitle = custTitle, drawGradVecs=FALSE, gd=gd )
	for ( i in 1:length( startAmpls ) ) {
		gd$guesses = list(a=startAmpls[i], s=startSigmas[i] )
		tuple = fitGauss1d( doPlot=FALSE, gd=gd )
		aLegs = tuple$amplLegs
		sLegs = tuple$sigmaLegs
		nLegs = length( aLegs )
		for ( j in 2:nLegs ) {
			x1 = aLegs[j-1]
			y1 = sLegs[j-1]
			x2 = aLegs[j]
			y2 = sLegs[j]
			arrows( x1, y1, x2, y2, col=trailColors[i], lwd=lw+1, angle=19, length=0.125 )
		}
		endAmpls[i] = aLegs[nLegs]
		endSigmas[i] = sLegs[nLegs]
	}
	# Draw beg/end circles last (so on top)
	for ( i in 1:length( startAmpls ) ) {
		points( startAmpls[i], startSigmas[i],  col=trailColors[i], cex=1.5, lwd=lw+1 )
		points( endAmpls[i],     endSigmas[i],  col=trailColors[i], cex=(i+0.5), lwd=lw+1 )
	}
}

### Plot Gaussian fit step vectors for multiple starting points
plotGaussianParamStepVectorsEx = function( nLevels=250, gridSize=3, learnRate=0.5, nIter=100 ) {
	set.seed( 3 )
	gd = makeGaussianTestData()
	trueS = gd$trueSigma
	trueA = gd$trueAmpl	
	maxA = trueA * 1.5
	minA = -maxA
	maxS = trueS * 3.5
	minS = -maxS
	ampls  = seq( minA * 0.9, maxA * 0.9, length=gridSize )
	sigmas = seq( minS * 0.9, maxS * 0.9, length=gridSize )
	gridAmpls  = seq( minA * 0.9, maxA * 0.9, length=400 )
	gridSigmas = seq( minS * 0.9, maxS * 0.9, length=400 )
	startAmpls = c()
	startSigmas = c()
	endAmpls = c()
	endSigmas = c()
	plotTitle = paste( "1D Gaussian Cost Function Contours + Gradient Descent Convergence Parameter Trails\n",
					   sprintf("Model: Amplitude = %.1f  Sigma = %.1f  gridSize = %d  learnRate = %.2f nIter = %d", 
					   		trueA, trueS, gridSize, learnRate, nIter ), sep="" )

	z = getGaussErrorGrid( meanGaussSsdErr, gridAmpls, gridSigmas, gd$x, gd$y )
	contour( gridAmpls, gridSigmas, t(z), nlevels=nLevels, xlab="", ylab="", main=plotTitle )
	addAxisLabels( xLabel="Amplitude", yLabel="Sigma", cexVal=1.7 )
	
	for ( s in sigmas ) {
		for ( a in ampls ) {
			startAmpls[length(startAmpls)+1] = a
			startSigmas[length(startSigmas)+1] = s
			gd$guesses = list( a=a, s=s )
			
			tuple = fitGauss1d( doPlot=FALSE, gd=gd, learnRate=learnRate, nIter=nIter )
			
			aLegs = tuple$amplLegs
			sLegs = tuple$sigmaLegs
			nLegs = length( aLegs )
			endAmpls[length(endAmpls)+1] = aLegs[nLegs]
			endSigmas[length(endSigmas)+1] = sLegs[nLegs]
			for ( j in 2:nLegs ) {
				x1 = aLegs[j-1]
				y1 = sLegs[j-1]
				x2 = aLegs[j]
				y2 = sLegs[j]
				arrows( x1, y1, x2, y2, col="blue", angle=19, lwd=lw, length=0.125 )
			}
		}
	}
	# Draw beg/end circles last (so on top)
	for ( i in 1:length( startAmpls ) ) {
		points( startAmpls[i], startSigmas[i],  col="chartreuse4", cex=1.5, lwd=lw+1 )
		points( endAmpls[i], endSigmas[i],  col="red2", cex=1.5, lwd=lw )
	}
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

### Call fitGauss1d() with increasing nIter values to get fitError as f( nIter ); plot result
plotGauss1dErrorVsnIter = function( gd=NULL ) {
	if ( is.null( gd ) )
		gd = makeGaussianTestData()
	meanAbsFitError = c()
	iter = seq( from=1, to=50, by=2 )
	plotGaussErrorContours( gd=gd )
	learnRate = 0.5
	for ( i in iter ) {
		meanAbsError = fitGauss1d( nIter=i, doPlot=FALSE, learnRate=learnRate )$meanAbsFitError
		meanAbsFitError[length(meanAbsFitError) + 1] = meanAbsError
	}
	plot( iter, meanAbsFitError, main="Gradient Descent Fit 1D Gaussian\nFit Error vs. Number Iterations",
		  xlab="", ylab="", cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale, cex.axis=1.3 )
	addAxisLabels( xLabel="Number Iterations", yLabel="Mean ABS Fit Error" )
	lines( iter, meanAbsFitError, col="blue", lwd=lw )
	text( 25, 0.63, sprintf("Learn Rate = %.02f", learnRate ), cex=1.3*fontScale )
}

### Use Gradient Descent and other methods to fit a line to data
fitLineGradientDescentVsOthers = function( nIter=100, learnRate=.01, nPts=10, drawContourPlot=FALSE, subTitle="" ) {
	set.seed( 3 )
	lmColor = "blue"
	exactDerivGdColor = "red"
	numDerivGdColor = "green3"
	trueM = 1
	trueB = 1
	coeff = c( trueB, trueM )
	x = seq( from=0, to=1, length=nPts )
	y = trueM * x + trueB + rnorm( length(x), sd=0.05 )
	xLimits = axisLimits( x )
	yLimits = axisLimits( y )	
	guesses = list( yInt=1.5, slope=2 )
	bEst = guesses$yInt
	mEst = guesses$slope
	bEstNum = guesses$yInt
	mEstNum = guesses$slope
	dm = 0
	db = 0
	histM = c( mEst )
	histB = c( bEst )
	numHistM = c( mEst )
	numHistB = c( bEst )
	# Gradient descent iteration using analytic and numeric derivatives
	for ( i in 1:nIter ) {
		# Analytic derivatives
		dm      = meanLineDeDm( x, y, mEst, bEst )
		db      = meanLineDeDb( x, y, mEst, bEst )
		mEst    = mEst - learnRate * dm
		bEst    = bEst - learnRate * db
		histM[length(histM)+1] = mEst
		histB[length(histB)+1] = bEst
		# Numerica estimates of derivatives
		dmNum   = lineDeDmNum( x, y, mEstNum, bEstNum ) # numerically estimated derivatives
		dbNum   = lineDeDbNum( x, y, mEst, bEstNum )
		mEstNum = mEstNum - learnRate * dmNum
		bEstNum = bEstNum - learnRate * dbNum
		numHistM[length(numHistM)+1] = mEstNum
		numHistB[length(numHistB)+1] = bEstNum
	}

	defTitle = "Fit a Line: Gradient Descent vs. Other Methods"
	plotTitle = if ( subTitle != "" ) sprintf( "%s\n%s", defTitle, subTitle ) else defTitle
	plot( x, y, main=plotTitle, cex=medPts,
		  xlab="", ylab="", lwd=lw+1, cex.main=1.3*fontScale, cex.axis=1.3 )	
	addAxisLabels( "x", "y" )
	
	# Now employ other methods to fit line
	ne = normalEq( x, y )               # 1) Normal Equation
	kr = fitLineKramersRule( x, y )     # 2) Solve linear system via Kramer's Rule
	lmFit = lm( y ~ x )                 # 3) Built-in lm() function
	lmB = lmFit$coefficients[1]
	lmM = lmFit$coefficients[2]
	lmMeanAbsErr =  (1/length(x))*sum( abs( lmFit$residuals ) )
	abline( bEst, mEst, col=exactDerivGdColor, lwd=lw+1, lty=4 )
	abline( bEstNum, mEstNum, col=numDerivGdColor, lwd=lw+1, lty=2 )
	abline( lmB, lmM, col=lmColor, lwd=lw+1, lty=3 )
	adMeanAbsErr = meanLineAbsErr( x, y, mEst, bEst)
	ndMeanAbsErr = meanLineAbsErr( x, y, mEstNum, bEstNum)
	krMeanAbsErr = meanLineAbsErr( x, y, kr$m, kr$b )
	neMeanAbsErr = meanLineAbsErr( x, y, ne[2], ne[1] )
	if ( verbose ) {
		print( sprintf( "Run params:" ) )
		print( sprintf( "    trueM     = %.2f", trueM ) )
		print( sprintf( "    trueB     = %.2f", trueB ) )
		print( sprintf( "    guessB    = %.2f", guesses$yInt ) )
		print( sprintf( "    guessM    = %.2f", guesses$slope ) )
		print( sprintf( "    nPts      = %d",   nPts ) )
		print( sprintf( "    learnRate = %.4f", learnRate ) )
		print( sprintf( "    nIter     = %d",   nIter ) )
		print( "Results:" )
		print( sprintf( "    True Model (noise added): slope = %.3f  yInt = %.3f", trueM, trueB ) )
		print( sprintf( "    Analytic Derivs:          slope = %.3f  yInt = %.3f  Mean ABS Error = %.3f", mEst,    bEst,    adMeanAbsErr ) )
		print( sprintf( "    Numeric Derivs:           slope = %.3f  yInt = %.3f  Mean ABS Error = %.3f", mEstNum, bEstNum, ndMeanAbsErr ) )
		print( sprintf( "    Kramer's Rule:            slope = %.3f  yInt = %.3f  Mean ABS Error = %.3f", kr$b,    kr$m,    krMeanAbsErr ) )
		print( sprintf( "    R lm():                   slope = %.3f  yInt = %.3f  Mean ABS Error = %.3f", lmB,     lmM,     lmMeanAbsErr ) )
		print( sprintf( "    Normal Equation:          slope = %.3f  yInt = %.3f  Mean ABS Error = %.3f", ne[1],   ne[2],   neMeanAbsErr ) )
	}	
	parSave = par( family="mono", font=2 )
	legText = c( sprintf( "Grad Desc Analytic Derivs: Mean ABS Err = %.4f", adMeanAbsErr ),
				 sprintf( "Grad Desc Numeric Derivs:  Mean ABS Err = %.4f", ndMeanAbsErr ),
				 sprintf( "R lm():                    Mean ABS Err = %.4f", lmMeanAbsErr ),
				 sprintf( "Normal Equation:           Mean ABS Err = %.4f", neMeanAbsErr ),
				 sprintf( "Kramer's Rule:             Mean ABS Err = %.4f", krMeanAbsErr ),
				 sprintf( "True/Model Slope: = %.2f", trueM ),
				 sprintf( "True/Model y-Int: = %.2f", trueB ),
				 sprintf( "lm() Slope:       = %.2f", lmM ),
				 sprintf( "lm() y-Int:       = %.2f", lmB ),
				 sprintf( "Grad Desc Slope   = %.2f", mEst ),
				 sprintf( "Grad Desc y-Int   = %.2f", bEst ),
				 sprintf( "# Iterations      = %d",   nIter ),
				 sprintf( "# Points          = %d",   nPts ),
				 sprintf( "Learn Rate        = %.3f", learnRate ) )
	legClrs = c( exactDerivGdColor, numDerivGdColor, lmColor, "transparent", "transparent", "transparent",
				 "transparent", lmColor, lmColor, exactDerivGdColor, exactDerivGdColor, "transparent", "transparent", "transparent" )
	legLtys = c( 4, 2, 3, 0, 0, 0, 0, 3, 3, 4, 4, 0, 0, 0 )
	legend( "topleft", legend=legText, col=legClrs, lty=legLtys, lwd=lw, inset = c(0.0,0.0),
		   bg="transparent", box.col="transparent", cex=1.2*fontScale, text.font=2 )
	par( parSave ) # Restore prior font
	
	if ( drawContourPlot )
		plotLineGradDescParamTrails( x, y, meanLineSsdErr, trueM, trueB, numHistM, numHistB,
									 takeSqrtError=FALSE, nIter=nIter, meanAbsFitError=adMeanAbsErr )	
}

### Slightly more reusable cost function plot + parameter evolution trails
### y-intercept maps to x-axis, slope maps to y-axis
plotLineGradDescParamTrails = function( x, y, errorFunc, modelM, modelB, histM, histB, takeSqrtError=TRUE, nIter, meanAbsFitError ) {
	nc = 200
	nr = 200	
	# Make contour plot region slightly larger than M and B parameter trails
	minM = min( histM ) * 0.8
	maxM = max( histM ) * 1.2
	minB = min( histB ) * 0.8
	maxB = max( histB ) * 1.2
	# Adjust these limits so that ***if plot is square*** the trails are perpendicular to contour lines
	deltaM = maxM - minM
	deltaB = maxB - minB
	gap = abs( deltaM - deltaB )
	if ( deltaM > deltaB ) {
		minB = minB - gap/2
		maxB = maxB + gap/2
	} else {
		minM = minM - gap/2
		maxM = maxM + gap/2
	}
	# Generate cost function sampling x and y coordinates
	intercepts = seq( minB, maxB, length.out=nc )
	slopes = seq( minM, maxM, length.out=nr )
	# Fill z 2D matrix with const function values for (y-intercept, slope) x,y pairs
	z = getLineErrorGrid( errorFunc, slopes, intercepts, x, y )
	if ( takeSqrtError )
		z = sqrt( z )
	# Create contour plot
	xLim = c( min(intercepts), max(intercepts) )
	yLim = c( min(slopes), max(slopes) )
	plotTitle = sprintf( "Fit Line using Gradient Descent\nCost Function Contours with Parameter Convergence Trail" )
	tz = t( z ) # contour() wants z first dim -> x-axis and second dim -> y-axis	
	contour( intercepts, slopes, tz, main=plotTitle, xlim=xLim, ylim=yLim, cex.lab=lf,
			 cex.main=1.3*fontScale, cex=1.3*fontScale, lwd=lw,
			 nlevels=50, labcex=1.3, xlab="", ylab="", cex.axis=1.3 )	
	addAxisLabels( "y-Intercept", "Slope", cexVal=1.5 )
	# Draw 'true' model (y-intercept, slope) crosshairs
	# Final converged gradient desccent params will != 'true' m and be because rand error added to model
	pl = gpl()
	lines( c( pl$xLim[1], pl$xLim[2] ), c( modelM, modelM ), col="blue", lwd=lw ) # Draw crosshair @ (b,m)
	lines( c( modelB, modelB ), c( pl$yLim[1], pl$yLim[2] ), col="blue", lwd=lw )
	for ( i in 2:length(histB) )
		arrows( histB[i-1], histM[i-1], histB[i], histM[i], col="red", lwd=lw, length=0.16, angle=17 )
	# Draw start/end point (initial guess + final converged params)
	x1 = histB[1];               y1 = histM[1]
	x2 = histB[length(histB)];   y2 = histM[length(histM)]	
	points( c(x1), c(y1), col="red", lwd=lw+1, cex=2 )
	points( c(x2), c(y2), col="red", lwd=lw+1, cex=2 )
	shift = 0.03
	text( x1 + shift, y1 + shift, sprintf( "(%.1f, %.1f)", x1, y1 ), cex=fontScale, adj=0, col="red" )
	text( x2 + shift, y2 + 0, sprintf( "(%.3f, %.3f)", x2, y2 ), cex=fontScale, adj=0, col="red" )
	lt1 = sprintf( "True model: slope = %.1f  y-int = %.1f", modelM, modelB )
	lt2 = "Convergence trail from init guesses"
	lt3 = sprintf( "# Iterations = %d", nIter )
	lt4 = sprintf( "Mean ABS Error = %.3f", meanAbsFitError )
	legText = c( lt1, lt2, lt3, lt4 )
	legClrs = c( "blue", "red", "transparent", "transparent" )
	# Below some hard-coded values for example
	legend( "topleft", legend=legText, col=legClrs, lwd=lw, cex=fontScale*1.2, inset=c(0.01,0.01) )
}

### Least-squares line fit cost function contour plot
plotLineFitErrorContours = function( plotCrossHairs=TRUE, plotLineFit=TRUE, plotErrorContour=TRUE, takeSquareRoot=FALSE, nLevels=50 ) {
	# Create fake data
	set.seed( 3 )
	trueB = 2
	trueM = 3
	coeff = c( trueB, trueM )
	x = seq( 0, 10, length=100 )
	y = evalPoly( x, coeff ) + rnorm( length(x), sd=2 )
	# Fit line so we know where center of bowl is located (best y-intercept & slope)
	lmFit = lm( y ~ x )
	fitB = lmFit$coefficients[1]
	fitM = lmFit$coefficients[2]
	# Draw line fit plot
	if ( plotLineFit ) {
		plot( x, y, main=sprintf( "Fit Line Using lm()\nTrue Slope = %.2f True y-Intercept = %.2f", trueM, trueB ),
			  xlab="", ylab="", cex=lrgPts, lwd=lw+1, cex.main=1.3*fontScale )
		abline( fitB, fitM, lwd=lw, col="blue" )
		
		e = compPolyFitError( x, y, lmFit$coefficients )
		addLegend( coeff, lmFit$coefficients, meanAbsErr=e$meanAbsError, meanSsdErr=e$meanSsdError )
		addAxisLabels()	
		grid()		
	}
	# Compute grid of cost function values for 'z' surface
	nc = 200
	nr = 200
	intercepts = seq( trueB - 2, trueB + 2, length.out=nc )
	slopes = seq( trueM - 0.5, trueM + 0.5, length.out=nr )
	z = getLineErrorGrid( meanLineSsdErr, slopes, intercepts, x, y )
	if ( plotErrorContour ) {
		# Create contour plot
		xLim = c( min(intercepts), max(intercepts) )
		yLim = c( min(slopes), max(slopes) )
		plotTitle = sprintf( "Least-Squares SSD Line Cost Function Contour Plot\nTrue Slope = %.2f True y-Intercept = %.2f (noise added)",
							 trueM, trueB)
		if ( takeSquareRoot )
			z = sqrt( z )
		contour( intercepts, slopes, t(z), main=plotTitle,
				 cex.axis=1.3, # Axis numeric label font scale factor
				 cex.main=1.3*fontScale, # Plot title font size
				 nlevels=nLevels,
				 labcex=1.4, # Contour line numeric labels,
				 xlab="", ylab="" )
		addAxisLabels( "y-Intercept", "Slope", cexVal=1.6 )		
		if ( plotCrossHairs ) {
			pl = gpl()
			lines( c(pl$xLim[1], pl$xLim[2]), c(trueM, trueM), col="blue", lwd=lw )
			lines( c(trueB, trueB), c(pl$yLim[1], pl$yLim[2]), col="blue", lwd=lw )
			lines( c(pl$xLim[1], pl$xLim[2]), c(fitM, fitM), col="red", lwd=lw )
			lines( c(fitB, fitB), c(pl$yLim[1], pl$yLim[2]), col="red", lwd=lw )
			legLine1 = "True Params (noise added)"
			legLine2 = "Fitted Params"
			legend( 0, 2.95, legend=c(legLine1,legLine2),col=c("blue","red"), lwd=lw, cex=fontScale )
		}
	}
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

# Basic R lm() and plot() example
basicPlotLmExample = function() {
	set.seed( 3 )  # For reproducibility
	modelM = 2
	modelB = 1
	x = 1:10
	y = modelM * x + modelB + rnorm( length(x), 1.5 ) # Add random noise
	lmFit = lm( y ~ x )
	fitB = lmFit$coefficients[1]
	fitM = lmFit$coefficients[2]
	plotTitle = paste( "R lm() and plot() Example\n",
					   sprintf( "modelM = %.2f  modelB = %.2f  fitM = %.2f fitB = %.2f",
					   		 modelM, modelB, fitM, fitB ), sep="" )
	plot( x, y, main=plotTitle, xlab="x-axis", ylab="y-axis",
		  cex.main=1.3,  # 'main' title font scale/size
		  cex.axis=1.3,  # Axis tick label font scale/size
		  cex.lab=1.3,   # Axis label font
		  cex=2.0,       # Data point size/diameter
		  lwd=2 )        # Data point line thickness scale
	abline( lmFit, lwd=2, col="blue" )
}

# Mockup lens parabolic auto-focus parameter estimation
lensFocusMockup = function() {
	x = c( 2, 4, 6, 8, 10, 12 )
	y = c( 10, 6, 5, 4.8, 6.25, 10.3 )
	plotTitle = "Lens Auto-Focus Parabolic Estimation Mockup"
	plot( x, y, ylim=c(0, max(y)+2), main=plotTitle, ylab="Edge Width (px)", xlab="Focus Motor Position", cex.main=1.3, cex.axis=1.3, cex.lab=1.3, cex=2.0, lwd=2 )
	lmFit = lm( y ~ x + I(x^2) )
	coeff = lmFit$coefficients
	xPred = seq( x[1] - 1, x[length(x)] + 1, length=100 )
	yPred = evalPoly( xPred, coeff )
	lines( xPred, yPred, col="blue", lwd=lw )
	a = coeff[3]; b = coeff[2]; c = coeff[1]
	focusPointX = -b/(2*a)
	focusPointY = a * focusPointX^2 + b * focusPointX + c
	minY = gpl()$yLim[1]
	arrows( focusPointX, minY, focusPointX, focusPointY, col="magenta", lwd=lw, angle=15, length = 0.3 )
}

### Run fitLineGradientDescentVsOthers() with different call parameters and save -> png files
runVsOthersAndSave = function( path ) { # nIter, learnRate, nPts
	f = fitLineGradientDescentVsOthers
	fName = paste( path, "fitLineGradientDescentVsOthers_", sep="" )
	
	numPts = 50; width=1200; height = 1200
	
	# Vary learnRate for 100 points
	subTitle = "Various Learn Rates"
	savePngVerbose( paste( fName, "_a.png", sep="" ), f, w=width, h=height, un="px", nIter=20, learnRate=0.05, nPts=numPts, subTitle=subTitle )
	savePngVerbose( paste( fName, "_b.png", sep="" ), f, w=width, h=height, un="px", nIter=20, learnRate=0.1,  nPts=numPts, subTitle=subTitle )
	savePngVerbose( paste( fName, "_c.png", sep="" ), f, w=width, h=height, un="px", nIter=20, learnRate=0.5,  nPts=numPts, subTitle=subTitle )

	# Vary nIter for 100 points
	subTitle = "Various # Iterations"
	savePngVerbose( paste( fName, "_d.png", sep="" ), f, w=width, h=height, un="px", nIter=20,  learnRate=0.05, nPts=numPts, subTitle=subTitle )
	savePngVerbose( paste( fName, "_e.png", sep="" ), f, w=width, h=height, un="px", nIter=100, learnRate=0.05, nPts=numPts, subTitle=subTitle )
	savePngVerbose( paste( fName, "_f.png", sep="" ), f, w=width, h=height, un="px", nIter=500, learnRate=0.05, nPts=numPts, subTitle=subTitle )	
	
	# This time trigger a contour param trail plot for save PNG target
	savePngVerbose( paste( fName, "_j.png", sep="" ), f, w=height, h=height, un="px", nIter=500, learnRate=0.2, nPts=numPts, drawContourPlot=TRUE )	
}

### Save plots as PNG files
runAllAndSave = function() {	
	dw = 1000
	dh = 750
	path = "d:/code/R/LeastSquaresPlay/PngPlots/"
	
	runVsOthersAndSave( path )
	
	savePngVerbose( paste(path, "plot15PointLineBruteForceFit.png", sep=""), plot15PointLineBruteForceFit, w=dw, h=dw, un="px" )	
	savePngVerbose( paste(path, "tempVsTime.png", sep=""), tempVsTime, w=dw, h=dw, un="px" )	
	savePngVerbose( paste(path, "plotLineCostFuncSurface.png", sep=""), plotLineCostFuncSurface, w=dw, h=dh, un="px" )		
	savePngVerbose( paste(path, "plotGaussianCostFuncExt.png", sep=""), plotGaussianCostFuncExt, w=1400, h=1400, un="px" )
	
	doStepVecTests = TRUE
	svSize = 1600
	# Learn rate default (0.5) nIter default (100)	
	if ( doStepVecTests ) {
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_001_LearnRate_0.5.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=1 )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_003_LearnRate_0.5.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=3 )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_005_LearnRate_0.5.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=5 )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_007_LearnRate_0.5.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=7 )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_011_LearnRate_0.5.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=11 )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_021_LearnRate_0.5.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=21 )	
		learnRate = 0.1 ; nIter = 1000
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_001_LearnRate_0.1.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=1, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_003_LearnRate_0.1.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=3, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_005_LearnRate_0.1.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=5, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_007_LearnRate_0.1.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=7, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_011_LearnRate_0.1.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=11, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_021_LearnRate_0.1.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=21, learnRate=learnRate, nIter=nIter )	
		learnRate = 0.25 ; nIter = 1000
	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_001_LearnRate_0.25.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=1, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_003_LearnRate_0.25.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=3, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_005_LearnRate_0.25.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=5, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_007_LearnRate_0.25.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=7, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_011_LearnRate_0.25.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=11, learnRate=learnRate, nIter=nIter )	
		savePngVerbose( paste(path, "plotGaussianParamStepVectorsEx_Grid_021_LearnRate_0.25.png", sep=""), 
						plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=21, learnRate=learnRate, nIter=nIter )	
	}
	
	gd = makeGaussianTestData()
	
	savePngVerbose( paste(path, "plotGaussErrorContours_NoVecs.png", sep=""), 
					plotGaussErrorContours, w=1200, h=1280, un="px", nLevels=110, drawGradVecs=FALSE, normGradVecs=TRUE, gd=gd )
	savePngVerbose( paste(path, "plotGaussErrorContours_NormVecs.png", sep=""), 
					plotGaussErrorContours, w=1200, h=1280, un="px", nLevels=110, drawGradVecs=TRUE, normGradVecs=TRUE, gd=gd )
	savePngVerbose( paste(path, "plotGaussErrorContours_UnNormVecs.png", sep=""), 
					plotGaussErrorContours,  w=1200, h=1280, un="px", nLevels=110, drawGradVecs=TRUE, normGradVecs=FALSE, gd=gd )

	savePngVerbose( paste(path, "plotGaussErrorContoursStereo_ContourAndPerspStereo.png", sep=""), 
					plotGaussErrorContoursStereo,  w=1400, h=1400, un="px" )
	savePngVerbose( paste(path, "plotGaussErrorContoursMono_ContourAndPerspMono.png", sep=""), 
					plotGaussErrorContoursMono, w=1000, h=1280, un="px" )
	savePngVerbose( paste(path, "plotGaussErrorSurfaceStereoPair.png", sep=""), 
					plotGaussErrorSurfaceStereoPair, w=1280, h=768, un="px" )	
	
	savePngVerbose( paste(path, "plotGaussianParamStepVectors.png", sep=""), 
					plotGaussianParamStepVectors, w=1200, h=1280, un="px", nLevels=120 )
	
	savePngVerbose( paste(path, "plotGaussErrorSurface.png", sep=""), 
					plotGaussErrorSurface, w=dw, h=dw, un="px" )
	savePngVerbose( paste(path, "plotGauss1dErrorVsnIter.png", sep=""), 
					plotGauss1dErrorVsnIter, w=dw, h=dh, un="px" )	
	
	savePngVerbose( paste(path, "fitGaussNls.png", sep=""), fitGaussNls, w=1280, h=800, un="px" )
	
	savePngVerbose( paste(path, "fitGauss1d_Iter_000.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=0, gd=gd )
	savePngVerbose( paste(path, "fitGauss1d_Iter_001.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=1, gd=gd )
	savePngVerbose( paste(path, "fitGauss1d_Iter_005.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=5, gd=gd )
	savePngVerbose( paste(path, "fitGauss1d_Iter_010.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=10, gd=gd )
	savePngVerbose( paste(path, "fitGauss1d_Iter_020.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=20, gd=gd )
	savePngVerbose( paste(path, "fitGauss1d_Iter_100.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=100, gd=gd )

	savePngVerbose( paste(path, "fitGauss1dNumDerivs_Iter_000.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=0,    numDerivs=TRUE, gd=gd )
	savePngVerbose( paste(path, "fitGauss1dNumDerivs_Iter_001.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=1,    numDerivs=TRUE, gd=gd )
	savePngVerbose( paste(path, "fitGauss1dNumDerivs_Iter_005.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=5,    numDerivs=TRUE, gd=gd )
	savePngVerbose( paste(path, "fitGauss1dNumDerivs_Iter_010.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=10,   numDerivs=TRUE, gd=gd )
	savePngVerbose( paste(path, "fitGauss1dNumDerivs_Iter_020.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=20,   numDerivs=TRUE, gd=gd )
	savePngVerbose( paste(path, "fitGauss1dNumDerivs_Iter_100.png", sep=""), 
					fitGauss1d, w=1280, h=dh, un="px", nIter=100,  numDerivs=TRUE, gd=gd )
	
	savePngVerbose( paste(path, "gradDemo.png", sep=""), gradDemo, w=dw, h=dw, un="px" )
	savePngVerbose( paste(path, "gradDemo2.png", sep=""), gradDemo2, w=dw, h=dw, un="px" )
	
	savePngVerbose( paste(path, "plotCubicWithTangents.png", sep=""), plotCubicWithTangents, w=1.2*600, h=1.2*900, un="px" )
	savePngVerbose( paste(path, "plotCubicWithDeriv.png", sep=""), plotCubicWithDeriv, w=1.2*600, h=1.2*900, un="px" )	
	savePngVerbose( paste(path, "plot5thOrderPoly.png", sep=""), plot5thOrderPoly, w=dw, h=dh, un="px" )
	savePngVerbose( paste(path, "plotInvParabolaWithTangents.png", sep=""), plotInvParabolaWithTangents, w=dh, h=dh, un="px" )
	savePngVerbose( paste(path, "plotInvParabolaWithVertexTangent.png", sep=""), plotInvParabolaWithVertexTangent, w=dw, h=dw, un="px" )
	savePngVerbose( paste(path, "plot2PointLine.png", sep=""), plot2PointLine, w=dw, h=dw, un="px" )
	
	savePngVerbose( paste(path, "FitLineToLinearData.png", sep=""), 
					plotNthOrderPolyToLinearData,   w=dw, h=dw, un="px", N=1, doPlot=TRUE )
	savePngVerbose( paste(path, "Fit2ndOrderPolyToLinearData.png", sep=""), 
					plotNthOrderPolyToLinearData,   w=dw, h=dw, un="px", N=2, doPlot=TRUE )
	savePngVerbose( paste(path, "Fit5thOrderPolyToLinearData.png", sep=""), 
					plotNthOrderPolyToLinearData,   w=dw, h=dw, un="px", N=5, doPlot=TRUE )
	savePngVerbose( paste(path, "Fit9thOrderPolyToLinearData.png", sep=""), 
					plotNthOrderPolyToLinearData,   w=dw, h=dw, un="px", N=9, doPlot=TRUE )
	savePngVerbose( paste(path, "Fit11thOrderPolyToLinearData.png", sep=""), 
					plotNthOrderPolyToLinearData, w=dw, h=dw, un="px", N=11, doPlot=TRUE, nPts=15 )
	savePngVerbose( paste(path, "Fit11thOrderPolyToLinearData13Pts.png", sep=""), 
					plotNthOrderPolyToLinearData, w=dw, h=dw, un="px", N=11, doPlot=TRUE, nPts=13 )
	savePngVerbose( paste(path, "Fit11thOrderPolyToLinearData30Pts.png", sep=""), 
					plotNthOrderPolyToLinearData, w=dw, h=dw, un="px", N=11, doPlot=TRUE, nPts=30 )
	
	savePngVerbose( paste(path, "plotLineFitErrorContours_LineFit.png", sep=""), 
					plotLineFitErrorContours, w=dw, h=dw, un="px", plotErrorContour=FALSE )
	savePngVerbose( paste(path, "plotLineFitErrorContours_Contour.png", sep=""), 
					plotLineFitErrorContours, w=dw, h=dw, un="px", plotLineFit=FALSE )
	savePngVerbose( paste(path, "plotLineFitErrorContours_ContourNoCrossHairs.png", sep=""), 
					plotLineFitErrorContours, w=dw, h=dw, un="px", plotLineFit=FALSE, plotCrossHairs=FALSE )
	savePngVerbose( paste(path, "plot15PointPoly2DataFitPoly2.png", sep=""), 
					plot15PointPoly2DataFitPoly2, w=800, h=800, un="px" )	
	
	savePngVerbose( paste(path, "genGaussianPlot.png", sep=""), genGaussianPlot, w=600, h=800, un="px" )	
	savePngVerbose( paste(path, "basicPlotLmExample.png", sep=""), basicPlotLmExample, w=900, h=900, un="px" )		
	savePngVerbose( paste(path, "lensFocusMockup.png", sep=""), lensFocusMockup, w=900, h=600, un="px" )		
}

