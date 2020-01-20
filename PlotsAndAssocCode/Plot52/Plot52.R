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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot52/"
	svSize = 1600
	# Learn rate default (0.5) nIter default (100)	
	learnRate = 0.25
	nIter = 1000
	savePngVerbose( paste(path, "Plot52.png", sep=""), 
					plotGaussianParamStepVectorsEx, w=svSize, h=svSize, un="px", gridSize=1, learnRate=learnRate, nIter=nIter )	
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

### Shared use case Gaussian data set generator used by several cooperating functions
makeGaussianTestData = function( trueAmpl=5, trueSigma=2, aGuess=sqrt(trueAmpl), sGuess=sqrt(trueSigma), nPts=40 ) {
	set.seed( 3 )
	x = seq( -5*trueSigma, 5*trueSigma, length=nPts )
	y = gauss1dFunc( x, trueAmpl, trueSigma ) + rnorm( length(x), sd=0.2 )
	guesses = list( a=aGuess, s=sGuess )
	tuple = list( trueAmpl=trueAmpl, trueSigma=trueSigma, x=x, y=y, nPts=nPts, guesses=guesses )
	return( tuple )
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

### Gaussian error/cost function - abs() of line lengths (not squared)
meanGaussAbsErr = function( x, y, a, sigma ) {
	return( (1/length(x)) * sum( abs(gauss1dFunc( x, a, sigma ) - y) ) )
}

### Gaussian error/cost function - SSD
meanGaussSsdErr = function( x, y, a, sigma ) {
	return( (1/length(x)) * sum( (gauss1dFunc( x, a, sigma ) - y)^2 ) )
}

### Analytical/exact partial derivative of error/cost function wrt Amplitude
gaussDeDa = function( x, y, a, sigma ) {
	term1 = sum( gauss1dFunc(x, a, sigma) * exp( -0.5 * (x/sigma)^2) )
	term2 = sum( y * exp( -0.5 * (x/sigma)^2) )
	return( 2 * ( term1 - term2 )/length(x) )
}

### Analytical/exact partial derivative of error/cost function wrt Sigma
gaussDeDs = function( x, y, a, sigma ) {
	v = 2*a/sigma^3
	term1 = sum( gauss1dFunc(x, a, sigma) * x^2 * exp( -0.5 * (x/sigma)^2) )
	term2 = sum( x^2 * y * exp( -0.5 * (x/sigma)^2) )
	return( v * ( term1 - term2 )/length(x) )
}

### Numerical estimation of partial derivative of error/cost function wrt Amplitude
gaussDeDaNum = function( x, y, a, sigma ) {
	e1 = meanGaussSsdErr( x, y, a, sigma )
	e2 = meanGaussSsdErr( x, y, a + de, sigma )
	return( (e2 - e1)/de )
}

### Numerical estimation of partial derivative of error/cost function wrt Sigma
gaussDeDsNum = function( x, y, a, sigma ) {
	e1 = meanGaussSsdErr( x, y, a, sigma )
	e2 = meanGaussSsdErr( x, y, a, sigma + de )
	return( (e2 - e1)/de )
}

### Get Plot Limits : gpl()$xLim[1] --> left x-coord
gpl = function() {
	u = par( "usr" )
	return( list( xLim=u[1:2], yLim=u[3:4] ) )
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
