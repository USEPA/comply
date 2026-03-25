# 1. Load necessary libraries
library(akima)



a <- "/work/ROMO/users/kpc/comply/dev/diag.level4.txt"
comply <- read.csv(file=a,skip=0,header=T)
comply$SECTOR <- trimws(comply$SECTOR)
comply$RADNAME <- trimws(comply$RADNAME)
comply$DOSE <- comply$IMMERSION + comply$INHALATION + comply$SURFACE + comply$MEAT + comply$MILK + comply$VEG
head(comply)

one <- subset(comply,STACK==1 & RADNAME=="I-131")
lookup <- c(N=0,NNE=22.5,NE=45,ENE=67.5,E=90,ESE=112.5,SE=135,SSE=157.5,S=180,SSW=202.5,SW=225,WSW=247.5,W=270,WNW=292.5,NW=315,NNW=337.5)
one$angle <- lookup[one$SECTOR]
one$x <- one$DISTANCE * sin(one$angle * 3.1416 / 180) # *pi/180 to convert degrees to radiens
one$y <- one$DISTANCE * cos(one$angle * 3.1416 / 180)
head(one)

# Create example irregular data
x_data <- one$x
y_data <- one$y
z_data <- one$DOSE
dat <- data.frame(x = one$x, y = one$y, z = z_data)

# Define the new grid points where you want to interpolate
#x_out <- seq(min(x_data), max(x_data), length.out = 120)
#y_out <- seq(min(y_data), max(y_data), length.out = 20)
x_out <- seq(-200, 200, by=2)
y_out <- x_out

# Perform bilinear interpolation
# The 'linear = TRUE' argument specifies bilinear interpolation
interpolated_data <- with(dat, akima::interp(x, y, z, xo = x_out, yo = y_out, linear = TRUE, extrap = TRUE))

# The result is a list with components x, y, and a matrix z
# You can use the results for plotting or further analysis

### PLOT
require(fields)

my.colors <- colorRampPalette(c("white","purple","deepskyblue2","green","yellow","orange","red","brown"))
maxval <- max(one$DOSE)*1.1
ncols <- 50

fname <- paste("/home/kbaker/comply/spatialplot1.png",sep="")
#pdf(file=fname,width=8,height=8)
png(file=fname,width=800,height=600)
par(mfrow=c(1,1))
par(mar=c(2.5,2.5,1.5,4))
### MAR (bottom, left, top, right)

#interpolated_data[interpolated_data > maxval] <- maxval
image.plot(interpolated_data,zlim=c(0,maxval),xlim=c(-85,85),ylim=c(-85,85),col=my.colors(ncols),legend.args=list(text="mrem/year",col="black",cex=2,side=2,line=0.5))

breaks <- seq(0,maxval,by=(maxval-0)/ncols)
nbreak <- length(breaks)+1
my.rsq.col <- my.colors(ncols)
#one$DOSE[one$DOSE > maxval] <- maxval
cut.rsq <- cut(one$DOSE,breaks)

points(one$x,one$y,col="black",bg=my.rsq.col[cut.rsq],pch=21,cex=2)
#contour(interpolated_data, add = TRUE)


box(col="black")
dev.off()

maxval