#====================================#
#### COMPLY v2.0 (levels 2 and 3) ####
#### November 26, 2025 version    ####
#====================================#

#Read in dose coefficient input files

a <- "/work/ROMO/users/kpc/comply/dev_r/comply171_original_input_data.csv"
one <- read.csv(file=a,skip=0,header=T)
one$xradname <- one$Nuclide
one$wbd <- one$WBdose.y 
one$vegdose <- one$VegDose
one$cowdose <- one$CowMilkDose
head(one)

four <- one[c("xradname","wbd","vegdose","cowdose")]
rm(one)
head(four) 


#Read in input namelist file
a <- "/work/ROMO/users/kpc/comply/dev_r/input.level3.txt"
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
phone <- substr(lines[6],21,500)

#Get number of stacks
nstacks <- as.numeric(as.character(substr(lines[7],21,500)))
print(paste("Number of stacks=",nstacks))

#start loop over number of stacks
for (nstack in 1:nstacks) {

#populate arrays for each stack
latitude <- as.numeric(as.character(substr(lines[8],21,500)))
longitude <- as.numeric(as.character(substr(lines[9],21,500)))
releaseht <- as.numeric(as.character(substr(lines[10],21,500)))
bldheight <- as.numeric(as.character(substr(lines[11],21,500)))
bldwidth <- as.numeric(as.character(substr(lines[12],21,500)))
distance <- as.numeric(as.character(substr(lines[13],21,500)))
distanceveg <- as.numeric(as.character(substr(lines[14],21,500)))
distancecow <- as.numeric(as.character(substr(lines[15],21,500)))
windspeed <- as.numeric(as.character(substr(lines[16],21,500)))

#Get number of radionuclides
nnucs <- as.numeric(as.character(substr(lines[17],21,500)))
print(paste("Number of radionuclides=",nnucs))

#set up arrays for number of radionuclides
radname <- rep(nnucs,0.0)
rademis <- rep(nnucs,0.0)
dose <- rep(nnucs,0.0)
vdose <- rep(nnucs,0.0)
cdose <- rep(nnucs,0.0)
totaldose <- rep(nnucs,0.0)

 #start loop over number of radionuclides
 for (nnuc in 1:nnucs) {
 tempx <- substr(lines[17+nnuc],21,500)
 temp <- unlist(strsplit(tempx,","))
 radname[nnuc] <- as.character(temp[1])
 rademis[nnuc] <- as.numeric(as.character(temp[2]))

 } #end loop over number of radionuclides in a given stack
} #end loop over number of stacks


#calculate whole body dose at nearest human receptor

	area <- bldwidth * bldheight 
        xcrit <- 2.5 * sqrt(area)
	xarea <- min(bldwidth,bldheight)

       	for (nnuc in 1:nnucs) {

        if ( distance < xcrit ) { 

	conc <- (0.25*rademis[nnuc])/(3.14*windspeed*1.0*xarea)

        } else {

        sigmaz <- (0.06*distance)/sqrt(1.0+(0.0015*distance))
	capsig <- sqrt((sigmaz^2) + (area/3.14))
        pp <- 2.032 / (distance*capsig)
	conc <- (0.25*rademis[nnuc]*pp) / windspeed

        }

	test <- subset(four,trimws(xradname)==radname[nnuc])
        dose[nnuc] = as.numeric(as.character(test$wbd)) * conc

	}

#calculate vegetable dose at garden receptor location

        area <- bldwidth * bldheight
        xcrit <- 2.5 * sqrt(area)
        xarea <- min(bldwidth,bldheight)

        for (nnuc in 1:nnucs) {

        if ( distanceveg < xcrit ) {

        conc <- (0.25*rademis[nnuc])/(3.14*windspeed*1.0*xarea)

        } else {

        sigmaz <- (0.06*distanceveg)/sqrt(1.0+(0.0015*distanceveg))
        capsig <- sqrt((sigmaz^2) + (area/3.14))
        pp <- 2.032 / (distanceveg*capsig)
        conc <- (0.25*rademis[nnuc]*pp) / windspeed

        }

        test <- subset(four,trimws(xradname)==radname[nnuc])
        vdose[nnuc] = as.numeric(as.character(test$vegdose)) * conc

        }


#calculate ingested milk/meat at farm receptor location

        area <- bldwidth * bldheight
        xcrit <- 2.5 * sqrt(area)
        xarea <- min(bldwidth,bldheight)

        for (nnuc in 1:nnucs) {

        if ( distancecow < xcrit ) {

        conc <- (0.25*rademis[nnuc])/(3.14*windspeed*1.0*xarea)

        } else {

        sigmaz <- (0.06*distancecow)/sqrt(1.0+(0.0015*distancecow))
        capsig <- sqrt((sigmaz^2) + (area/3.14))
        pp <- 2.032 / (distancecow*capsig)
        conc <- (0.25*rademis[nnuc]*pp) / windspeed

        }

        test <- subset(four,trimws(xradname)==radname[nnuc])
        cdose[nnuc] = as.numeric(as.character(test$cowdose)) * conc

        }

#calculate cumulative dose over all pathways

	allnucsdose <- 0.0
        for (nnuc in 1:nnucs) {
        print(radname[nnuc])
        totaldose[nnuc] = dose[nnuc] + vdose[nnuc] + cdose[nnuc]
	allnucsdose <- allnucsdose + totaldose[nnuc]
        print(totaldose[nnuc])
	}

        print(allnucsdose)

#Output data

allout <- NULL
 for (nnuc in 1:nnucs) {

temp <- c(radname[nnuc],rademis[nnuc],totaldose[nnuc],dose[nnuc],vdose[nnuc],cdose[nnuc])
allout <- rbind(allout,temp)
 }

colnames(allout) <- c("Nuclide","Emissions (Ci/y)","Total Dose (mrem/y)","Inhalation Dose (mrem/y)","Vegetable Dose (mrem/y)","Cow+Milk Dose (mrem/y)")
allout 

write.table(allout, file="output.level3.csv", sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)


######QA data calculate total dose by distances out to 2 km 

        alldata <- NULL

        receptors <- seq(5,2000,5)
        receptors

        area <- bldwidth * bldheight
        xcrit <- 2.5 * sqrt(area)
        xarea <- min(bldwidth,bldheight)

        for (nnuc in 1:nnucs) {

        for (dist in 1:length(receptors)) { 

        distance <- receptors[dist]

        if ( distance < xcrit ) {

        conc <- (0.25*rademis[nnuc])/(3.14*windspeed*1.0*xarea)

        } else {

        sigmaz <- (0.06*distance)/sqrt(1.0+(0.0015*distance))
        capsig <- sqrt((sigmaz^2) + (area/3.14))
        pp <- 2.032 / (distance*capsig)
        conc <- (0.25*rademis[nnuc]*pp) / windspeed

        } #end condition on critical distance

        test <- subset(four,trimws(xradname)==radname[nnuc])

	dose[nnuc] = as.numeric(as.character(test$wbd)) * conc
        vdose[nnuc] = as.numeric(as.character(test$vegdose)) * conc
        cdose[nnuc] = as.numeric(as.character(test$cowdose)) * conc
        total = dose[nnuc]+vdose[nnuc]+cdose[nnuc]
        alldata <- rbind(alldata,c(radname[nnuc],distance,total,dose[nnuc],vdose[nnuc],cdose[nnuc]))

        } #end loop over receptors

        } #end loop over nuclides

colnames(alldata) <- c("Nuclide","Distance(m)","Total Dose (mrem/y)","Inhalation Dose (mrem/y)","Vegetable Dose (mrem/y)","Cow+Milk Dose (mrem/y)")

write.table(alldata, file="output.level3.bydistance.csv", sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)

