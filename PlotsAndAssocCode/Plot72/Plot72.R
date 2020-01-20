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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot72/"
	savePngVerbose( paste(path, "Plot72.png", sep=""), 
					plotLineFitErrorContours, w=dw, h=dw, un="px", plotLineFit=FALSE, plotCrossHairs=FALSE )
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

### Compute f(x) for a polynomial and a vector of x-coordinates
evalPoly = function( x, coeff ) {
	if ( length( coeff ) < 1 ) return( c(0) )  	
	termSum = 0
	for ( i in 1:length(coeff) ) {
		termSum = termSum + coeff[i] * x^(i-1)
	}
	return( termSum )
}

# Compute errors, ABS and SSD, return tuple
compPolyFitError = function( x, y, coeff ) {
	predY = evalPoly( x, coeff )
	meanAbsError = (1/length(x)) * sum( abs(predY - y) )
	meanSsdError = (1/length(x)) * sum( (predY - y)^2 )
	return( list( meanAbsError=meanAbsError, meanSsdError=meanSsdError ) )
}

### Get Plot Limits : gpl()$xLim[1] --> left x-coord
gpl = function() {
	u = par( "usr" )
	return( list( xLim=u[1:2], yLim=u[3:4] ) )
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
