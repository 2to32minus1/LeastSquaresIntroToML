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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot21/"
	savePngVerbose( paste(path, "Plot21.png", sep=""), fitGaussNls, w=1280, h=800, un="px" )
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

### Get simulated Gamma plot data: x[], y[]
getSimGammaData = function( mu=mu, sigma=sigma, peakHt=peakHt) {
	bkgCoeff = getQuadraticBkgCoeff()
	x = 1:100
	y = bkgPlusGaussian( x, bkgCoeff, mu, sigma, peakHt )
	return( list(x=x, y=y ) )
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

### Cubic background + gaussian peak function
bkgPlusGaussian = function( x, bkgCoeff, mu, sigma, peakHt ) {
	set.seed( 3 )
	return( evalPoly( x, bkgCoeff ) + 		        # background
				gauss( x, mu, sigma, peakHt ) +     # Gaussian
				rnorm( length(x), sd=7 ) )          # random noise
}

### Get simulated Gamma plot data: x[], y[]
getSimGammaData = function( mu=mu, sigma=sigma, peakHt=peakHt) {
	bkgCoeff = getQuadraticBkgCoeff()
	x = 1:100
	y = bkgPlusGaussian( x, bkgCoeff, mu, sigma, peakHt )
	return( list(x=x, y=y ) )
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

### Evaluate gamma sim given x and curve params
evalGammaFunc = function( x, quadC0, quadC1, quadC2, mu, sigma, peakHt ) {
	return( evalPoly( x, c(quadC0, quadC1, quadC2) ) + gauss( x, mu, sigma, peakHt ) )
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
