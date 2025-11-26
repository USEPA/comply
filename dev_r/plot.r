#Read in input namelist file
a <- "/work/ROMO/users/kpc/comply/dev_r/input.level3.txt"
lines <- readLines(a)
for (line in lines) {
        print(line)
}

#Get descriptive information
projectname <- substr(lines[1],21,500)
#populate arrays for each stack
latitude <- as.numeric(as.character(substr(lines[8],21,500)))
longitude <- as.numeric(as.character(substr(lines[9],21,500)))
distance <- as.numeric(as.character(substr(lines[13],21,500)))
distanceveg <- as.numeric(as.character(substr(lines[14],21,500)))
distancecow <- as.numeric(as.character(substr(lines[15],21,500)))


a <- "output.bydistance.csv"
one <- read.csv(file=a,skip=0,header=T)
one$RADNAME <- one$Nuclide
one$DOSE <- one$Total.Dose..mrem.y.
one$DISTANCE <- one$Distance.m.
head(one)

temp <- one[c("RADNAME")]
temp2 <- unique(temp)
rads <- temp2$RADNAME
rads
nrads <- length(rads)
nrads

colors <- c("red","blue","orange")

#####################################################
#####################################################


fname <- paste("qaplot.png",sep="")
#pdf(file=fname,width=8,height=8)
png(file=fname,width=600,height=400)
par(mfrow=c(1,1))
par(mar=c(4,4,2,1))
### MAR (bottom, left, top, right)

labtext <- "mrem/yr"
lim <- max(one$DOSE)

temp <- one
plot(temp$DISTANCE,temp$DOSE,xlab="Distance from source (m)",ylab="Dose (mrem/yr)",col="white",cex=0.1,ylim=c(0,lim),xlim=c(0,300))

legend("topright",c(rads),pch=c(1,1,1),col=c(colors),cex=1.25,bg="white",ncol=1)
legend("top",c("COMPLY v1.7"),ncol=1,cex=1.1,bty='n',col="black")

for (i in 1:length(rads)) {
temp <- subset(one,RADNAME==rads[i])
#points(temp$DISTANCE,temp$ALTDOSE,pch=19,cex=0.66,lwd=1,col=colors[i])
points(temp$DISTANCE,temp$DOSE,pch=1,cex=2,lwd=1,col=colors[i])

}

abline(h=10,col="lightgray",lwd=2)
abline(v=distance,col="darkgray",lwd=2,lty=2)
abline(v=distanceveg,col="darkgray",lwd=2,lty=2)
abline(v=distancecow,col="darkgray",lwd=2,lty=2)


#lines(temp$week,temp$AOD_ob,pch=19,lwd=4,col="black")



dev.off()

