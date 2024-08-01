##Photosynthesis and Respiration code October 2019
##created by Danielle Becker 03/02/20
##edited by Danielle Becker 04/02/20

#clear lists
rm(list=ls())

#Read in required libraries

##### Include Versions of libraries
remotes::install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('tidyverse')


# get the file path
setwd("Respirometry/")
path.p<-"../Respirometry/Data/Oct_TT" #the location of all your respirometry files

getwd()

# bring in the all files for on big respirometry file
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#generate a 3 column dataframe with specific column names
Photo.R<- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=4))
colnames(Photo.R) <- c("fragment.ID.full","Intercept", "umol.L.sec","Temp.C")
View(Photo.R)

#Load your respiration data file, with all the times, etc.
Sample.Info <- read.csv(file="Data/resp_data_TT.csv", header=T) #read in sample.info data
view(Sample.Info)

# load surface area data
SA <- read.csv(file=paste0(path.p,"/../sample_info_TT.csv"), header=T) #read sample.info data

#change the volume of water to numeric
as.numeric(SA$volume)

View(SA)

# joint the sample info and surface area and volume measurements
Sample.Info<-left_join(Sample.Info, SA)

View(Sample.Info)


# make start and stop times real times
Sample.Info$start.time <- as.POSIXct(Sample.Info$start.time,format="%H:%M:%S", tz = "") #convert time from character to time
Sample.Info$stop.time <- as.POSIXct(Sample.Info$stop.time,format="%H:%M:%S", tz = "") #convert time from character to time


# for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
for(i in 1:length(file.names.full)) { # for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
  
  #find the lines in sample info that have the same file name that is being brought in
  FRow<-which(Sample.Info$fragment.ID==strsplit(file.names[i],'.csv'))
  
  # read in the O2 data one by one
  Photo.Data1 <-read.csv(file.path(path.p,file.names.full[i]), skip = 1, header=T) # skips the first line
  Photo.Data1  <- Photo.Data1[,c("Time","Value","Temp")] #subset columns of interest
  Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  Photo.Data1 <- na.omit(Photo.Data1) #omit NA from data frame
  
  
  # clean up some of the data
  n<-dim(Photo.Data1)[1] # length of full data
  Photo.Data1 <-Photo.Data1[240:(n-3),] #start at data point ~4 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Photo.Data1)[1] #list length of trimmed data
  Photo.Data1$sec <- (1:n) #set seconds by one from start to finish of run in a new column
  
  
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  rename <- sub(".csv","", file.names[i]) # remove all the extra stuff in the file name
  
  pdf(paste0("Output/",rename,"thinning.pdf")) # open the graphics device
  
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot (empty plot to fill) data as a function of time
  usr  <-  par('usr') # extract the size of the figure margins
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) # put a grey background on the plot
  whiteGrid() # make a grid
  box() # add a box around the plot
  points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1) # add the x axis
  axis(2, las=1) # add the y-axis
  
  # Thin the data to make the code run faster
  Photo.Data.orig <-Photo.Data1 #save original unthinned data
  Photo.Data1 <-  thinData(Photo.Data1 ,by=20)$newData1 #thin data by every 20 points for all the O2 values
  Photo.Data1$sec <- as.numeric(rownames(Photo.Data1 )) #maintain numeric values for time
  Photo.Data1$Temp<-NA # add a new column to fill with the thinned data
  Photo.Data1$Temp <-  thinData(Photo.Data.orig,xy = c(1,3),by=20)$newData1[,2] #thin data by every 20 points for the temp values
  
  # plot the thinned data
  plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "eq", "pc")
  Regs  <-  rankLocReg(xall=Photo.Data1$sec, yall=Photo.Data1$Value, alpha=0.5, method="pc", verbose=TRUE)  
  
  # add the regression data
  plot(Regs)
  dev.off()
  
  
  # fill in all the O2 consumption and rate data
  Photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[i,1] <- rename #stores the file name in the Date column
  Photo.R[i,4] <- mean(Photo.Data1$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
  #Photo.R[i,5] <- PR[j] #stores whether it is photosynthesis or respiration
  
  
  # rewrite the file everytime... I know this is slow, but it will save the data that is already run
}
write.csv(Photo.R, 'Output/Photo.R.csv')  

#read in Photo.R file so dont need to run entire for loop again
Photo.R <- read.csv('Output/Photo.R.csv')

# Calculate P and R rate
Photo.R$fragment.ID.full<-Photo.R$fragment.ID
Photo.R$fragment.ID<-NULL
View(Photo.R)

Photo.R<-left_join(Photo.R, Sample.Info)
View(Photo.R)

#Convert sample volume to mL
Photo.R$volume <- Photo.R$volume/1000 #calculate volume

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Photo.R$umol.sec <- Photo.R$umol.L.sec*Photo.R$volume

#Account for blank rate by temperature
#convert character columns to factors
Photo.R <- Photo.R %>%
  mutate_if(sapply(., is.character), as.factor)
View(Photo.R)

#make the blank column a factor
Photo.R$BLANK<-ifelse(Photo.R$treatment=='BLANK', 1,0)
Photo.R$BLANK<-as.factor(Photo.R$BLANK)
View(Photo.R)

#aggregate certain columns together for easier viewing
photo.blnk <- aggregate(umol.sec ~ species*temp.Cat*light_dark*BLANK, data=Photo.R, mean)

#pull out only the blanks
#photo.blnk<-photo.blnk[photo.blnk$Species=='BK',]
photo.blnk<-photo.blnk[photo.blnk$BLANK==1,]

# remove the blank column
photo.blnk$BLANK<-NULL

# rename the blank rate
colnames(photo.blnk)[4]<-'blank.rate'  

# join the blank data with the rest of the data
Photo.R<-left_join(Photo.R, photo.blnk)
View(Photo.R)

# subtract the blanks
Photo.R$umol.sec.corr<-Photo.R$umol.sec-Photo.R$blank.rate

View(Photo.R)

#### Normalize to SA (surface area)#####

#Calculate net P and R
Photo.R$umol.cm2.hr <- (Photo.R$umol.sec.corr*3600)/Photo.R$surf.area.cm2 #mmol cm-2 hr-1

# remove NAs and blanks
Photo.R<-Photo.R[Photo.R$BLANK==0,]

#make respiration positive
Photo.R$umol.cm2.hr<- Photo.R$umol.cm2.hr

# log the rates
Photo.R$Rate.ln<-log(Photo.R$umol.cm2.hr+0.1)

#remove empty rows
Photo.R<-Photo.R[-which(is.na(Photo.R$fragment.ID.full)),]

#making the respiration values positive (pull out data for dark treaments)
Photo.R$umol.cm2.hr[Photo.R$light_dark=="dark"]<-Photo.R$umol.cm2.hr[Photo.R$light_dark=="dark"]*-1 
lessthan <- which(Photo.R$light_dark=="dark" & Photo.R$umol.cm2.hr < 0)
Photo.R$umol.cm2.hr[lessthan] <- 0
View(Photo.R)

ggplot(Photo.R, aes(x=Temp.C, y=umol.cm2.hr,group = fragment.ID, col = fragment.ID))+
  geom_line(size=2)+
  geom_point()+  
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*light_dark, labeller = labeller(.multi_line = FALSE))+
  ggsave(filename = "Output/lowvshigh_curves.pdf", device = "pdf", width = 10, height = 10)
 
write.csv(Photo.R, 'Output/TT_Rates.csv') # export all the uptake rates

View(Photo.R)


###calculating gross photosynthesis from data frame###

#make ifelse statements to assign light treatments as NP and dark treatments as resp
#light will be assigned NP for net photosynthesis 
Photo.R$rate.type <-ifelse(Photo.R$light_dark=='light', "NP", "R")
Photo.R$rate.type<-as.factor(Photo.R$rate.type)

View(Photo.R)

#rename fragment ID 
Photo.R$individual.ID <- str_split(Photo.R$fragment.ID, "_", n = Inf, simplify = TRUE)[,1]
Photo.R$individual.ID <- as.factor(Photo.R$individual.ID)

#rename dataframe to work from Photo.R2 to preserve Photo.R
#make column for GP and group by fragment ID and temp to keep R and NP together
Photo.R2 <- Photo.R %>% 
  filter(temp.Cat <= 37) %>%
  group_by(individual.ID, temp.Cat, treatment) %>% 
  summarize(rates = sum(umol.cm2.hr), Temp.C=mean(Temp.C)) %>%
  mutate(rate.type="GP", light_dark="L") %>%
  rename(umol.cm2.hr = rates) %>%
  mutate(fragment.ID=paste0(individual.ID, "_", light_dark))

#rename light darks
Photo.R2$light_dark <- ifelse(Photo.R2$light_dark=="L", "light", "dark")

#make the negative gphoto values zeros
Photo.R2$umol.cm2.hr[Photo.R2$umol.cm2.hr<0]<-0

Photo.R <- Photo.R[,c("individual.ID", "temp.Cat","treatment", "Temp.C", "umol.cm2.hr", "rate.type", "light_dark", "fragment.ID")]
view(Photo.R)

#left join Photo.R and Photo.R2
Photo.T <- rbind(Photo.R, as.data.frame(Photo.R2))
View(Photo.T)


write.csv(Photo.T, '../Thermal_Performance/Data/Photo.T.csv') # export all the uptake rates

#check all rate types (GP, NP and R) curves
ggplot(Photo.T, aes(x=Temp.C, y=umol.cm2.hr, group = individual.ID, col = individual.ID))+
  geom_line(size=2)+
  geom_point()+  
  theme_bw ()+
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*rate.type, labeller = labeller(.multi_line = FALSE))+ 
  ggsave(filename = "../Respirometry/Output/initialcurves.pdf", device = "pdf", width = 10, height = 10)

#check three various rate types (GP, NP and R) individaul curves
Photo.T.GP <- Photo.T %>%
  filter(rate.type =="GP")

ggplot(Photo.T.GP, aes(x=Temp.C, y=umol.cm2.hr, group = individual.ID, col = individual.ID))+
  geom_line(size=2)+
  geom_point()+  
  theme_bw ()+
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*rate.type*fragment.ID, labeller = labeller(.multi_line = FALSE))+ 
  ggsave(filename = "../Respirometry/Output/GPcurves.pdf", device = "pdf", width = 10, height = 10)

Photo.T.NP <- Photo.T %>%
  filter(rate.type =="NP")

ggplot(Photo.T.NP, aes(x=Temp.C, y=umol.cm2.hr, group = individual.ID, col = individual.ID))+
  geom_line(size=2)+
  geom_point()+  
  theme_bw ()+
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*rate.type*fragment.ID, labeller = labeller(.multi_line = FALSE))+ 
  ggsave(filename = "../Respirometry/Output/NPcurves.pdf", device = "pdf", width = 10, height = 10)

Photo.T.R <- Photo.T %>%
  filter(rate.type =="R")

ggplot(Photo.T.R, aes(x=Temp.C, y=umol.cm2.hr, group = individual.ID, col = individual.ID))+
  geom_line(size=2)+
  geom_point()+  
  theme_bw ()+
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*rate.type*fragment.ID, labeller = labeller(.multi_line = FALSE))+ 
  ggsave(filename = "../Respirometry/Output/Rcurves.pdf", device = "pdf", width = 10, height = 10)


