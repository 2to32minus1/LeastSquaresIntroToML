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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot7/"
	savePngVerbose( paste(path, "Plot7.png", sep=""), 
					plotNthOrderPolyToLinearData, w=dw, h=dw, un="px", N=11, doPlot=TRUE, nPts=13 )	
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

### Get Plot Limits : gpl()$xLim[1] --> left x-coord
gpl = function() {
	u = par( "usr" )
	return( list( xLim=u[1:2], yLim=u[3:4] ) )
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

# Compute errors, ABS and SSD, return tuple
compPolyFitError = function( x, y, coeff ) {
	predY = evalPoly( x, coeff )
	meanAbsError = (1/length(x)) * sum( abs(predY - y) )
	meanSsdError = (1/length(x)) * sum( (predY - y)^2 )
	return( list( meanAbsError=meanAbsError, meanSsdError=meanSsdError ) )
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

### Compute f(x) for a polynomial and a vector of x-coordinates
evalPoly = function( x, coeff ) {
	if ( length( coeff ) < 1 ) return( c(0) )  	
	termSum = 0
	for ( i in 1:length(coeff) ) {
		termSum = termSum + coeff[i] * x^(i-1)
	}
	return( termSum )
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
