#====================================#
#### COMPLYr v2.0                 ####
#### February 26, 2026 version    ####
#====================================#


#Read in input namelist file
a <- "/work/ROMO/users/kpc/comply/dev_r/input.level4r.txt"
lines <- readLines(a)
for (line in lines) {
	print(line)
}

#Get descriptive information
projectname <- substr(lines[1],21,500)
companyname <- substr(lines[2],21,500)
facilityname <- substr(lines[3],21,500)
address <- substr(lines[4],21,500) 
contact <- substr(lines[5],21,500)
email <- substr(lines[6],21,500)

#Get number of stacks
nstacks <- as.numeric(as.character(substr(lines[7],21,500)))
print(paste("Number of stacks=",nstacks))


allout <- NULL
alldiagout <- NULL
totaldose <- 0.0

#start loop over number of stacks
for (nstack in 1:nstacks) {

#first line for each stack is a spacer that should be ignored
latitude <- as.numeric(as.character(substr(lines[7+2+29*(nstack-1)],21,500)))
longitude <- as.numeric(as.character(substr(lines[7+3+29*(nstack-1)],21,500)))
releaseht <- as.numeric(as.character(substr(lines[7+4+29*(nstack-1)],21,500)))
exitvelocity <- as.numeric(as.character(substr(lines[7+5+29*(nstack-1)],21,500)))
    exitvelocity[is.na(exitvelocity) | exitvelocity <= 0] <- 0.3
diameter <- as.numeric(as.character(substr(lines[7+6+29*(nstack-1)],21,500)))
stacktemp <- as.numeric(as.character(substr(lines[7+7+29*(nstack-1)],21,500)))
airtemp <- as.numeric(as.character(substr(lines[7+8+29*(nstack-1)],21,500)))
verticalmomentum <- substr(lines[7+9+29*(nstack-1)],21,500)
rademis <- as.numeric(as.character(substr(lines[7+10+29*(nstack-1)],21,500)))
metinfo <- substr(lines[7+11+29*(nstack-1)],21,500)
pctcalm <- as.numeric(as.character(substr(lines[7+12+29*(nstack-1)],21,500)))
nwindbins <- as.numeric(as.character(substr(lines[7+13+29*(nstack-1)],21,500)))

#set up arrays for sectors
windbinname <- rep(nwindbins,0.0)
windbinfreq <- rep(nwindbins,0.0)
windbinspeed <- rep(nwindbins,0.0)
windbinsrdist <- rep(nwindbins,0.0)
sectorname <- rep(nwindbins,0.0)
conc <- rep(nwindbins,0.0)
ihdose <- rep(nwindbins,0.0)

maxdose <- 0.0
maxsector <- 0
maxdist <- 0.0

#start loop over number of wind sector bins
 for (nnuc in 1:nwindbins) {
 tempx <- substr(lines[7+13+29*(nstack-1)+nnuc],21,500)
 temp <- unlist(strsplit(tempx,","))
 windbinname[nnuc] <- as.character(temp[1])
 windbinfreq[nnuc] <- as.numeric(as.character(temp[2]))
 windbinspeed[nnuc] <- as.numeric(as.character(temp[3]))
 windbinsrdist[nnuc] <- as.numeric(as.character(temp[4]))
 sectorname[nnuc] <- as.character(temp[5]) 

 if(windbinsrdist[nnuc]>maxdist) { maxdist <- windbinsrdist[nnuc] }

#calculate air concentration for each sector bin
 HEFF <- releaseht

 if (releaseht > 0.0) {

 RHOAIR <- 0.0735 * ( 80.0 + 460.0 ) / ( airtemp + 460.0 )
 RHOAIR <- RHOAIR * 454.0 / ( 0.3048^3 )
 RHOSTK <- RHOAIR * ( airtemp + 460.0 ) / ( stacktemp + 460.0 )

 AREA <- 3.14 * ( diameter^2.0 ) / 4.0
 STKVEL <- exitvelocity / AREA
 R <- STKVEL / windbinspeed[nnuc]
 FM <- ( RHOSTK / RHOAIR ) * STKVEL * exitvelocity / 3.14
 TERM1 <- 75.0 * FM / ( ( ( 1.0 + 3.0 / R ) * windbinspeed[nnuc] )^2.0 )
 TERM2 <- releaseht + 7.79 * sqrt(FM) / ( ( 1.0 + 3.0 / R ) * windbinspeed[nnuc] )

  if (verticalmomentum == "YES") {
  DELHMO <- 0.93 * (TERM1^(3/7)) * (TERM2^(1/7))
  } else {
  DELHMO <- 0.0
  }

 FB <- 9.81 * ( RHOAIR - RHOSTK ) * exitvelocity / ( 3.14 * RHOAIR )
 DELHBU <- 23.67 * ( ( FB / ( windbinspeed[nnuc]^3.0 ) )^0.6 ) * ( releaseht + 195.1 * FB / ( windbinspeed[nnuc]^3.0))^0.4

 HEFF <- releaseht + max(DELHMO, DELHBU, na.rm=TRUE)
 
 } #end condition where release height > 0

 SIGMAZ <- (0.06 * windbinsrdist[nnuc]) / sqrt(1.0 + 0.0015 * windbinsrdist[nnuc])

 CAPSIG <- sqrt(SIGMAZ^2 + 0.0 / 3.14159265)

 zz <- 0.5 * ( (HEFF / CAPSIG)^2 )

   if (zz > 0.001) {
   pp <- 2.032 * exp(-zz) / (windbinsrdist[nnuc] * CAPSIG)
   } else {
   pp <- 2.032 / (windbinsrdist[nnuc] * CAPSIG)
   }

 conc[nnuc] <- ((windbinfreq[nnuc]/(1.0-pctcalm)) *rademis*pp) / windbinspeed[nnuc]
 if (!is.na(conc[nnuc]) && conc[nnuc] < 1e-32) {
  conc[nnuc] <- 1e-32
 }

 #calculate inhalation dose

 INHFAC <- .0092 *   100.   * 1000.
 TINIT <- 11.8
 T <- TINIT + windbinsrdist[nnuc] / ( windbinspeed[nnuc] * 60. )

 EQILOUT <- 1 - 0.0479 * exp(-T / 4.39) - 2.1963 * exp(-T / 38.6) + 1.2442 * exp(-T / 28.4)
 EQILOUT <- min(EQILOUT, 1.0)  # caps at 1; preserves NA
 EQILIN <-  0.35 * ( 1. + EQILOUT )
 EQILIN <- min(EQILIN, 0.7)  # caps at 0.7; preserves NA
 INHDOS <- conc[nnuc] *  1.E12  *  (1./1000.)  * (1./100.) * (8760./170.) * INHFAC
 FINDOOR <- .75

 ihdose[nnuc] <- INHDOS * ( FINDOOR * EQILIN + (1.0 - FINDOOR) * EQILOUT )

 if (ihdose[nnuc] > maxdose) { 
	 maxsector <- nnuc
	 maxdose <- ihdose[nnuc]
 }

 } #end loop over wind sector bins

 for (nnuc in 1:nwindbins) {

temp <- c(nstack,rademis,nnuc,maxsector,sectorname[nnuc],windbinsrdist[nnuc],ihdose[nnuc],conc[nnuc],windbinfreq[nnuc],windbinspeed[nnuc])
allout <- rbind(allout,temp)
 }

 totaldose <- totaldose + ihdose[maxsector]





 for (nnuc in 1:nwindbins) {

## calc conc and dose for diagnostic output file
 receptors <- seq(50,maxdist,50)

 for (dist in 1:length(receptors)) {

 SIGMAZ <- (0.06 * receptors[dist]) / sqrt(1.0 + 0.0015 * receptors[dist])
 CAPSIG <- sqrt(SIGMAZ^2 + 0.0 / 3.14159265)

 zz <- 0.5 * ( (HEFF / CAPSIG)^2 )

   if (zz > 0.001) {
   pp <- 2.032 * exp(-zz) / (receptors[dist] * CAPSIG)
   } else {
   pp <- 2.032 / ( receptors[dist] * CAPSIG)
   }

 conc[nnuc] <- ((windbinfreq[nnuc]/(1.0-pctcalm)) *rademis*pp) / windbinspeed[nnuc]
 if (!is.na(conc[nnuc]) && conc[nnuc] < 1e-32) {
  conc[nnuc] <- 1e-32
 }

 INHFAC <- .0092 *   100.   * 1000.
 TINIT <- 11.8
 T <- TINIT + receptors[dist] / ( windbinspeed[nnuc] * 60. )

 EQILOUT <- 1 - 0.0479 * exp(-T / 4.39) - 2.1963 * exp(-T / 38.6) + 1.2442 * exp(-T / 28.4)
 EQILOUT <- min(EQILOUT, 1.0)  # caps at 1; preserves NA
 EQILIN <-  0.35 * ( 1. + EQILOUT )
 EQILIN <- min(EQILIN, 0.7)  # caps at 0.7; preserves NA
 INHDOS <- conc[nnuc] *  1.E12  *  (1./1000.)  * (1./100.) * (8760./170.) * INHFAC
 FINDOOR <- .75

 ihdose[nnuc] <- INHDOS * ( FINDOOR * EQILIN + (1.0 - FINDOOR) * EQILOUT )

 temp <- c(nstack,rademis,nnuc,maxsector,sectorname[nnuc],receptors[dist],ihdose[nnuc],conc[nnuc])
 alldiagout <- rbind(alldiagout,temp)
 } #end loop over diagnostic distances
 } #end loop over sectors

} #end loop over number of stacks


colnames(allout) <- c("Stack","Emissions (Ci/y)","Sector","Max Sector","Sector Name","Distance (m)","Total Dose (mrem/y)","Conc (Ci/m3)","Wind Frequency","Wind Speed (m/s)")
allout 

write.table(allout, file="output.level4r.csv", sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

colnames(alldiagout) <- c("Stack","Emissions (Ci/y)","Sector","Max Sector","Sector Name","Distance (m)","Total Dose (mrem/y)","Conc (Ci/m3)")

write.table(alldiagout, file="output.level4r.bydist.csv", sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

print(totaldose) #COMPLIANCE VALUE!


