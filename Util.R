##############################################################
### Copyright (c) Richard Creamer 2019 - All Rights Reserved
### Inquiries: email 2to32minus1@gmail.com
##############################################################

# TODO: use image( t( z^0.2 ), useRaster=TRUE, axes=FALSE, col=gray(seq(0,1,length=255) ) )

# this works to put quadratic + gaussian into a legend text with math rendering
# legend( "bottomright", legend=parse(text=c( "1+2*x^2+3.24*x^3+e^(-x^2/2*sigma^2)" ) ), inset=c(0.01,0.01) )

### Globals - Constant(s) and Helper functions
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

#######################
# Helper functions below
#######################

### Return string suitable for plot titles with actual exponents (not '^')
### Examples: [ "1.23", "1.23-2.21x", "1.23-2.21x+4.98x^2" ]
prettyPoly = function( coeff, nSigFigs=3, wrapLimit=3 ) {
	rc = signif( coeff, nSigFigs )
	numTerms = length( coeff )
	terms = c()
	for ( termNum in 1:numTerms ) {
		if ( termNum == 1 ) {
			s = paste( rc[termNum], sep="" )
		} else {
			s = paste( abs(rc[termNum]), sep="" )
		}
		if ( termNum > 1 )
			s = paste( s, "*x", sep="" )
		if ( termNum > 2 )
			s = paste( s, "^", (termNum-1), sep="" )
		terms[ length(terms) + 1 ] = s
	}
	wrapLines = list()
	catStr = ""
	for ( termNum in 1:numTerms ) {
		if ( termNum > 1 && ( termNum - 1 ) %% wrapLimit == 0 ) {
			wrapLines[length(wrapLines) + 1] = catStr
			catStr = ""
		}
		if ( termNum == 1 ) {
			catStr = paste( catStr, terms[termNum], sep="" )
		} else {
			if ( rc[termNum] < 0 ) {
				catStr = paste( catStr, " - ", terms[termNum], sep="" )
			} else {
				catStr = paste( catStr, " + ", terms[termNum], sep="" )
			}
		}
	}
	if ( catStr != "" )
		wrapLines[length(wrapLines) + 1] = catStr
	return( wrapLines )
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

### Add a legend to a plot suitable for this deck's polynomials/purpose
addLegend = function( coeff, fitCoeff, meanAbsErr=NULL, meanSsdErr=NULL, fontSz=1.35*fontScale, pos="topleft", wrapLimit=5 ) {
	legList = list()
	
	# Generate and append legend text lines for 'true model'
	truePolyStrList = prettyPoly( coeff, nSigFigs=3, wrapLimit=wrapLimit )
	legList[1] = paste( "'True model': y ==", truePolyStrList[1], sep="" )
	if ( length( truePolyStrList ) > 1 ) {
		for ( i in 2:length( truePolyStrList ) ) {
			legList[length(legList)+1] = paste( "'     '", truePolyStrList[i], sep="" )
		}
	}
	
	# Generate and append legend text lines for 'fitted model'
	fitPolyStrList = prettyPoly( fitCoeff, nSigFigs=3, wrapLimit=wrapLimit )
	legList[length(legList) + 1] = paste( "'Fitted model': y ==", fitPolyStrList[1], sep="" )
	if ( length( fitPolyStrList ) > 1 ) {
		for ( i in 2:length( fitPolyStrList ) ) {
			legList[length(legList)+1] = paste( "'     '", fitPolyStrList[i], sep="" )
		}
	}
	
	if ( !is.null(meanAbsErr) ) {
		legError = paste( "'Mean ABS Error': ", sprintf( "%.2f", meanAbsErr ), sep="" )
		legList[ length( legList ) + 1 ] = legError
	}
	if ( !is.null(meanSsdErr) ) {
		legError = paste( "'Mean SSD Error': ", sprintf( "%.2f", meanSsdErr ), sep="" )
		legList[ length( legList ) + 1 ] = legError
	}
	
	legend( pos, bty="n",inset=c(0.005,0.01),
			legend=parse(text=legList), col=c("transparent","transparent"),
			cex=fontSz, pch=15, y.intersp=1.1 )
}

### Save to PNG file
savePng = function( path, plotFunc ) {
	png( filename = path, type="cairo", units="in", width=5, height=5, pointsize=12, res=96 )
	plotFunc()
	addCopyright()
	dev.off()
}

### Save to PNG file, specify width and height
savePngVerbose = function( path, plotFunc, w=512, h=512, un="px", doCopyright=TRUE, ... ) {
	png( filename = path, type="cairo", units=un, width=w, height=h, pointsize=12, res=96 )
	plotFunc( ... )
	if ( doCopyright )
		addCopyright()
	dev.off()
}

### Convience method to add x/y axis labels
addAxisLabels = function( xLabel="x", yLabel="y", cexVal=1.3 ) {
	mtext( text=xLabel, side=1, line=2.5, cex=cexVal )
	mtext( text=yLabel, side=2, line=2.5, cex=cexVal )
}

### Line fit cost function derivative wrt y-intercept
meanLineDeDb = function( x, y, m, b ) {
	return( (1/length(x))*2*sum( (m*x + b) - y ) )
}

### Line fit cost function derivative wrt slope
meanLineDeDm = function( x, y, m, b ) {
	return( (1/length(x))*2*sum( ( (m*x + b) - y ) * x ) )
}

### Add copyright notice to plot via text()
addCopyright = function() {
	mtext( "Copyright \uA9 2019 Richard Creamer - All Rights Reserved", side=4, line=0, adj=0, cex=1.1 )
	mtext( "Email: 2to32minus1@gmail.com", side=4, line=1, adj=0, cex=1.1 )
}

### Compute partial derivative of 1D line fit cost function wrt slope
lineDeDmNum = function( x, y, slope, yInt ) {
	delta = 0.0000001
	e1 = meanLineSsdErr( x, y, slope, yInt )
	e2 = meanLineSsdErr( x, y, slope + delta, yInt )
	return ( (e2 - e1)/delta )
}

### Compute partial derivative of 1D line fit cost function wrt slope
lineDeDbNum = function( x, y, slope, yInt ) {
	delta = 0.0000001
	e1 = meanLineSsdErr( x, y, slope, yInt )
	e2 = meanLineSsdErr( x, y, slope, yInt + delta )
	return ( (e2 - e1)/delta )
}

### Fit line via linear system of equations + Kramer's rule
fitLineKramersRule = function( x, y ) {
	n = length(x)
	d = n * sum(x^2) - (sum(x))^2
	b = (1/d) * ( sum(y) * sum(x^2) - sum(x) * sum(x*y) )
	m = (1/d) * ( n * sum(x*y) - sum(x) * sum(y) )
	return( list( b=b, m=m ) )
}

### Linear regression using Normal Equation; works for 1+ sample features
normalEq = function( x, y ) {
	X = cbind( 1, x ) # cbind() = 'column bind': create matrix: prepend col of 1's to x
	return( solve( t(X) %*% X ) %*% t(X) %*% y )
}

### Get Plot Limits : gpl()$xLim[1] --> left x-coord
gpl = function() {
	u = par( "usr" )
	return( list( xLim=u[1:2], yLim=u[3:4] ) )
}

### Gaussian related code

### Cubic background + gaussian peak function
bkgPlusGaussian = function( x, bkgCoeff, mu, sigma, peakHt ) {
	set.seed( 3 )
	return( evalPoly( x, bkgCoeff ) + 		        # background
				gauss( x, mu, sigma, peakHt ) +     # Gaussian
				rnorm( length(x), sd=7 ) )          # random noise
}

### Evaluate gamma sim given x and curve params
evalGammaFunc = function( x, quadC0, quadC1, quadC2, mu, sigma, peakHt ) {
	return( evalPoly( x, c(quadC0, quadC1, quadC2) ) + gauss( x, mu, sigma, peakHt ) )
}

### Gaussian function
gauss = function( x, mu, sigma, peakHt ) {
	return( peakHt * exp( -0.5 * ((x-mu)/sigma)^2 ) )
}

### Get simulated Gamma plot data: x[], y[]
getSimGammaData = function( mu=mu, sigma=sigma, peakHt=peakHt) {
	bkgCoeff = getQuadraticBkgCoeff()
	x = 1:100
	y = bkgPlusGaussian( x, bkgCoeff, mu, sigma, peakHt )
	return( list(x=x, y=y ) )
}

### Return coeff for a cubic polynomial simulating gamma background curve
getQuadraticBkgCoeff = function(drawPlot = FALSE) {
	x = c( 0, 20, 40, 80, 100 )
	y = c( 170, 135, 110, 65, 50 )
	lmFit = lm( y ~ x + I(x^2) )
	quadCoeff = lmFit$coefficients
	if ( drawPlot ) {
		plot( x, y, xlim=c(0, 100), ylim=c(0,350), cex.lab=lf, xlab="", ylab="" )
		addAxisLabels()	
		xPred = seq( x[1], x[length(x)], length.out=200 )
		yPred = evalPoly( xPred, quadCoeff )
		lines( xPred, yPred, col=4, lwd=lw )
	}
	return( quadCoeff ) # coefficients of fitted 2nd order polynomial
}

### Returns an elevation value as f(x,y) for down-oriented gaussian w/base at z=bias
elevFunc = function( x, y, sigma=3, peakHt=0.5, bias=1 ) {
	r = sqrt( x^2 + y^2 )
	elev = bias - peakHt * exp( -0.5 * (r/sigma)^2 )
	return( elev )
}

### Exact analytical value for dz/dx for inverted gaussian
gaussDzDx = function( x, y, sigma, peakHt ) {
	return( ( (peakHt*x)/sigma^2 ) * exp( -(x^2 + y^2)/(2*sigma^2) ) )
}

### Exact analytical value for dz/dy for inverted gaussian
gaussDzDy = function( x, y, sigma, peakHt ) {
	return( ( (peakHt*y)/sigma^2 ) * exp( -(x^2 + y^2)/(2*sigma^2) ) )
}	

### Analytical/exact partial derivative of error/cost function wrt Amplitude
gaussDeDa = function( x, y, a, sigma ) {
	term1 = sum( gauss1dFunc(x, a, sigma) * exp( -0.5 * (x/sigma)^2) )
	term2 = sum( y * exp( -0.5 * (x/sigma)^2) )
	return( 2 * ( term1 - term2 )/length(x) )
}

### Numerical estimation of partial derivative of error/cost function wrt Amplitude
gaussDeDaNum = function( x, y, a, sigma ) {
	e1 = meanGaussSsdErr( x, y, a, sigma )
	e2 = meanGaussSsdErr( x, y, a + de, sigma )
	return( (e2 - e1)/de )
}

### Analytical/exact partial derivative of error/cost function wrt Sigma
gaussDeDs = function( x, y, a, sigma ) {
	v = 2*a/sigma^3
	term1 = sum( gauss1dFunc(x, a, sigma) * x^2 * exp( -0.5 * (x/sigma)^2) )
	term2 = sum( x^2 * y * exp( -0.5 * (x/sigma)^2) )
	return( v * ( term1 - term2 )/length(x) )
}

### Numerical estimation of partial derivative of error/cost function wrt Sigma
gaussDeDsNum = function( x, y, a, sigma ) {
	e1 = meanGaussSsdErr( x, y, a, sigma )
	e2 = meanGaussSsdErr( x, y, a, sigma + de )
	return( (e2 - e1)/de )
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

### Get quadratic + gaussian math-renderable text string; must call parse() on ret val
quadPlusGaussianMathText = function( c0, c1, c2, ampl, sigma, mu ) {
	polyStr = prettyPoly( c(c0, c1, c2) )
	expStr = getGaussianMathText( ampl, sigma, mu )
	retStr = paste( polyStr, "+", expStr, sep="" )
	return( retStr )
}

### Get gaussian math-renderable text string; must call parse() on ret val
getGaussianMathText = function( ampl, sigma, mu ) {
	if ( mu == 0 ) {
		expStr = sprintf( "%.2f*e^frac(-x^2, 2*(%.2f)^2)", ampl, sigma )		
	}
	else {
		expStr = sprintf( "%.2f*e^frac(-(x-%.2f)^2,2*(%.2f)^2)", ampl, mu, sigma )
	}
	return( expStr )
}

### Generate small Gaussian plot
genGaussianPlot = function( x1=-5, x2=5, n=21, ampl=100, sigma=1, xLab="x", yLab="y" ) {
	x = seq( x1, x2, length=n )
	y = ampl * exp( -(x-mean(x))^2/(2*sigma^2) )
	plot( x,y,xlab="", ylab="", cex=medPts, lwd=lw+1, cex.axis=1.3 )
	addAxisLabels( xLab, yLab )
	xPred = seq( x[1], x[length(x)], length=n*10 )
	yPred = ampl * exp( -(xPred-mean(xPred))^2/(2*sigma^2) )
	lines( xPred, yPred, col="blue", lwd=lw )
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

### Print a Gaussian data set
printGaussData = function( gd ) {
	print( "Gaussian Data Set:" )
	print( sprintf( "trueAmpl:  %.2f", gd$trueAmpl ) )
	print( sprintf( "trueSigma: %.2f", gd$trueSigma ) )
	print( sprintf( "nPts:      %d",   gd$nPts ) )
	print( sprintf( "Guesses$a: %.2f", gd$guesses$a ) )
	print( sprintf( "Guesses$s: %.2f", gd$guesses$s ) )
	print( "x:" ); print( gd$x )
	print( "y:" ); print( gd$y )
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

### Line fit error function: abs() not diff^2
meanLineAbsErr = function( x, y, slope, yInt ) {
	return( (1/length(x)) * sum( abs( yInt + slope*x - y ) ) )
}

### Line fit error function
meanLineSsdErr = function( x, y, slope, yInt ) {
	return( (1/length(x)) * sum( ( yInt + slope*x - y )^2 ) )
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

### Temp test function, delete before sharing source file on web
imgTest = function( fName="lenna.grn" ) {
	istImgDir = "D:/AutogenyData/Images/IstImages/"
	path = paste( istImgDir, fName, sep="" )
	f = file( path, "rb" )
	pixels = readBin( f, integer(), n=262144, size=1, signed=FALSE )
	close( f )
	
	im1 = matrix( data=pixels,nrow=512, ncol=512 )
	image( im1, useRaster=TRUE, axes=FALSE, col=gray(seq(0,1,length=256)) )
	im2 = im1[ , 512:1]
	image( im2, useRaster=TRUE, axes=FALSE, col=gray(seq(0,1,length=256)) )
}

### Temp test function, delete before sharing source file on web
demTest = function( demDir="d:/AutogenyData/Gis/Dems/UsgsExampleDems/", fName="susanville.dem", nr=601, nc=601 ) {
	path = paste( demDir, fName, sep="" )
	f = file( path, "rb" )
	pixels = data=readBin( f, integer(), n=nr*nc, size=2, signed=FALSE )
	im16Bits = matrix( data=pixels,nrow=nr, ncol=nc )
	im16BitsFlipped = im16Bits[ , nc:1]
	im8BitsFlipped = im16BitsFlipped/256
	image( im8BitsFlipped, useRaster=TRUE, axes=FALSE, col=gray(seq(0,1,length=256)) )
	contour( 0:(nc-1), 0:(nr-1), t(im16BitsFlipped), nlevels=30 )
	persp( 0:(nc-1), 0:(nr-1), t(im8BitsFlipped), theta=0, phi=30, ltheta=-120, expand=0.1, shade=0.9, border=NA, box=FALSE, col="lightblue", scale=FALSE )
	close( f )
}

# Build Nth order polynomial formula string
makePolyFormula = function( indepVar="x", depVar="y", polyOrder=1 ) {
	formulaStr = paste( depVar, " ~ ", indepVar, sep="" )
	if ( polyOrder > 1 ) {
		for ( i in 2:polyOrder )
			formulaStr = paste( formulaStr, " + I(", indepVar, "^", i, ")", sep="" )
	}
	return( formulaStr )
}

# Compute errors, ABS and SSD, return tuple
compPolyFitError = function( x, y, coeff ) {
	predY = evalPoly( x, coeff )
	meanAbsError = (1/length(x)) * sum( abs(predY - y) )
	meanSsdError = (1/length(x)) * sum( (predY - y)^2 )
	return( list( meanAbsError=meanAbsError, meanSsdError=meanSsdError ) )
}

# Call plot function w/blue background
blueBkgPlot = function( pf, ... ) {
	parSave = par( fg="white", col.main="white", bg = rgb(31/255, 78/255, 121/255) )
	pf(...)
	par( parSave )
}

### Gen stereo graphic
genStereoGraphic = function() {
	path = "d:/code/R/LeastSquaresPlay/PngPlots/"
	savePngVerbose( paste(path, "stereoGraphic.png", sep=""), plotGaussErrorSurfaceStereoPairGraphic, w=1200, h=400, un="px", doCopyright=FALSE )	
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

