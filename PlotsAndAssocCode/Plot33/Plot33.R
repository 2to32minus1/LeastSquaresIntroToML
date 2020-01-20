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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot33/"
	savePngVerbose( paste(path, "Plot33.png", sep=""), lensFocusMockup, w=900, h=600, un="px" )
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

### Compute f(x) for a polynomial and a vector of x-coordinates
evalPoly = function( x, coeff ) {
	if ( length( coeff ) < 1 ) return( c(0) )  	
	termSum = 0
	for ( i in 1:length(coeff) ) {
		termSum = termSum + coeff[i] * x^(i-1)
	}
	return( termSum )
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
