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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot30/"
	savePngVerbose( paste(path, "Plot30.png", sep=""), genGaussianPlot, w=600, h=800, un="px" )	
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
