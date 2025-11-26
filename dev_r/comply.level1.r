#====================================#
#### COMPLY v2.0 (levels 2 and 3) ####
#### November 26, 2025 version    ####
#====================================#

#Read in dose coefficient input files

a <- "/work/ROMO/users/kpc/comply/dev_r/comply171_original_possession_tables.csv"
one <- read.csv(file=a,skip=0,header=T)
one$xradname <- one$Nuclide
head(one)

four <- one[c("xradname","Gas","Liquid","Solid")]
rm(one)
head(four) 


#Read in input namelist file
a <- "/work/ROMO/users/kpc/comply/dev_r/input.level1.txt"
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

#Get number of radionuclides
nnucs <- as.numeric(as.character(substr(lines[7],21,500)))
print(paste("Number of radionuclides=",nnucs))

#set up arrays for number of radionuclides
radname <- rep(nnucs,0.0)
radtype <- rep(nnucs,0.0)
rademis <- rep(nnucs,0.0)
radpossession <- rep(nnucs,0.0)

 #start loop over number of radionuclides
 for (nnuc in 1:nnucs) {
 tempx <- substr(lines[7+nnuc],21,500)
 temp <- unlist(strsplit(tempx,","))
 radname[nnuc] <- as.character(temp[1])
 radtype[nnuc] <- as.character(temp[2])
 rademis[nnuc] <- as.numeric(as.character(temp[3]))

 } #end loop over number of radionuclides in a given stack

#compare possession limits

       	for (nnuc in 1:nnucs) {

	test <- subset(four,trimws(xradname)==radname[nnuc])

	if (radtype[nnuc]=='G') {
        radpossession[nnuc] = as.numeric(as.character(test$Gas))
        }

        if (radtype[nnuc]=='L') {
        radpossession[nnuc] = as.numeric(as.character(test$Liquid))
        }

	if (radtype[nnuc]=='S') {
        radpossession[nnuc] = as.numeric(as.character(test$Solid))
        }


	}

#calculate total possession

        total <- 0.0
        for (nnuc in 1:nnucs) {
        print(radname[nnuc])
        total <- total + rademis[nnuc]/radpossession[nnuc]
        }

	a <- paste('You possess ',round(total,1),' times the allowable amount given in the possession table.')
        print(a)


#Output data

allout <- NULL
 for (nnuc in 1:nnucs) {

temp <- c(radname[nnuc],radtype[nnuc],rademis[nnuc],radpossession[nnuc],rademis[nnuc]/radpossession[nnuc])
allout <- rbind(allout,temp)
 }

colnames(allout) <- c("Nuclide","Phase","Possession Amount (Ci/y)","Legal Limit (Ci/y)","Percent of Limit")
allout 

write.table(allout, file="output.level1.csv", sep=",",quote=FALSE,row.names=FALSE,col.names=TRUE)



