setwd("~C:/Users/flo_f/OneDrive - University of Rhode Island/GitHub/Coral_Priming_Experiments_Summer_2024/data/environmental/Tris") #set working directory

Calib.Data <-read.table("20240724.csv", header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
model <-lm(mVTris ~ Ttris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
coe <- coef(model) #extracts the coeffecients
R2<-summary(model)$r.squared
pdf("C:/Users/flo_f/OneDrive - University of Rhode Island/GitHub/Coral_Priming_Experiments_Summer_2024/output/environmental/TrisCal/20240724_triscal.pdf")
plot(mVTris ~ Ttris, data=Calib.Data)
abline(lm(mVTris ~ Ttris, data=Calib.Data))
legend('topleft', legend = bquote(R^2 == .(format(R2, digits = 3))), bty='n')

dev.off()

