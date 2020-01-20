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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot22/"
	f = fitLineGradientDescentVsOthers	
	numPts = 50; width=1200; height = 1200	
	subTitle = "Various Learn Rates"
	savePngVerbose( paste( path, "Plot22.png", sep="" ), f, w=width, h=height, un="px", nIter=20, learnRate=0.05, nPts=numPts, subTitle=subTitle )
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

### Line fit error function: abs() not diff^2
meanLineAbsErr = function( x, y, slope, yInt ) {
	return( (1/length(x)) * sum( abs( yInt + slope*x - y ) ) )
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
