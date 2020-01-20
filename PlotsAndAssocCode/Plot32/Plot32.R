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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot32/"
	savePngVerbose( paste(path, "Plot32.png", sep=""), gradDemo2, w=dw, h=dw, un="px" )
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
