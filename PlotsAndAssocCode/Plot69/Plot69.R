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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot69/"
	savePngVerbose( paste(path, "Plot69.png", sep=""), plotInvParabolaWithVertexTangent, w=dw, h=dw, un="px" )
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
