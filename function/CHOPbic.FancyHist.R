# Draw a fancy histogram as in http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=151
# Use the original source codes with modification
CHOPbic.FancyHist<-function(x, main='Histogram', xlab='', ylab='', freq=TRUE, breaks='Sturges', ...) {
# x       Numeric, data to be drawn
# main, xlab, ylab    Character, same as in hist()
# breaks, freq        Same as in hist()

# axis labels always horizontal
# margin smaller than default
par( las = 1, mar = c(5,5,3,3))  

# calculate histogram (but do not display)
h <- hist( x, plot = FALSE, breaks = breaks )

# calculate density estimation
d <- density( x ) 

if (freq) ylab='Count' else ylab='Frequency'

# just set the scene and axes
plot( h , border = NA, freq = freq, xlab = xlab, ylab = ylab, main= main, cex.lab=1.25, cex.main=1.5) 

# grab the limits of the region
usr <- par( "usr" )
ncolors <- 100
dy <- ( usr[4] - usr[3] ) / ncolors

# create the colors
colors <- colorRampPalette( c("yellow","orange","red") )(ncolors)

# horizontal line
abline( h = axTicks(2) , col = "gray", lwd = .5 )

# for each color, clip into a region and redraw the histogram with that color
for( i in 1:ncolors){
  clip( usr[1], usr[2], usr[3] + (i-1) * dy, usr[3] + i*dy )
  plot( h, add = TRUE, axes = FALSE, ylab = "", xlab = "", col = colors[i], border = NA, freq = freq)
}
# reset the clipping area. See ?clip
do.call( clip, as.list( usr) )

# just to get the boxes right
plot( h, add = TRUE, lwd = .5 , freq = freq, xlab = "", ylab = "", axes = FALSE )

if (freq) d[[2]]<-d[[2]]*length(x)/length(h[[1]])

# adds a tick and translucent line to display the density estimate
lines( d, lwd = 3, col = "#22222288" )

# add a rug (also using translucency)
rug(x, col = "#00000088" ) 

# box around the plot
box()

}


