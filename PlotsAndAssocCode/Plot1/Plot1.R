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
	path = "D:/Code/R/LeastSquaresPlay/PlotsAndAssocCode/Plot1/"
	savePngVerbose( paste(path, "Plot1.png", sep=""), basicPlotLmExample, w=900, h=900, un="px" )	
}

### Basic R lm() and plot() example
basicPlotLmExample = function() {
	set.seed( 3 )  # For reproducibility
	modelM = 2
	modelB = 1
	x = 1:10
	y = modelM * x + modelB + rnorm( length(x), 1.5 ) # Add random noise
	lmFit = lm( y ~ x )
	fitB = lmFit$coefficients[1]
	fitM = lmFit$coefficients[2]
	plotTitle = paste( "R lm() and plot() Example\n",
					   sprintf( "modelM = %.2f  modelB = %.2f  fitM = %.2f fitB = %.2f",
					   		 modelM, modelB, fitM, fitB ), sep="" )
	plot( x, y, main=plotTitle, xlab="x-axis", ylab="y-axis",
		  cex.main=1.3,  # 'main' title font scale/size
		  cex.axis=1.3,  # Axis tick label font scale/size
		  cex.lab=1.3,   # Axis label font
		  cex=2.0,       # Data point size/diameter
		  lwd=2 )        # Data point line thickness scale
	abline( lmFit, lwd=2, col="blue" )
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
