library(rgbif)
library(terra)
library(geodata)
library(predicts)
setwd("data")

#occ = sp_occurrence("Lacerta","schreiberi") #download occurences from GBIF
#occ = sp_occurrence("Lacerta","schreiberi")
#write.csv(occ,"occurences.csv")
occ = read.csv("occurences.csv")


# download daily average temperatures from worldclim averaged by month
tavg <- worldclim_global(var = 'tmax', res = 10, path = ".")

# plots map of daily average temperature for July (month 7) with European extent
plot(tavg[[7]], ext=ext_eur)

# download daily solar radiation from worldclim averaged by month
srad <- worldclim_global(var = 'srad', res = 10, path = ".")

# plots map of daily solar radiation for July with European extent
plot(srad[[7]], ext=ext_eur)

# sample climate data for occurences
tavg_sample = extract(tavg[[7]],vect(occ))
srad_sample = extract(srad[[7]],vect(occ))
plot(srad_sample[,2],tavg_sample[,2])

world_data <- world(path = ".") #download world countries limits
ext_eur <- c(-15,45,35,72) # coordinates of European extent
ext_eur2 <- c(-15,20,35,50) # more reduced extent for pseudoabsences

# plot the map of Europe
plot(world_data,
     xlab = "Longitude", ylab = "Latitude", axes = TRUE, ext=ext_eur)

# add the species occurrences
points(occ$lon,occ$lat,cex=.1,col="blue")


# create pseudo-absences
eumask<-crop(tavg[[1]], ext(ext_eur))
absences<-backgroundSample(eumask,nrow(occ),vect(occ),tryf=10)
points(absences,cex=.1,col="red")


# plot absences in climate space
tavg_abs = extract(tavg[[7]],vect(absences))
srad_abs = extract(srad[[7]],vect(absences))

plot(srad_sample[,2],tavg_sample[,2])
points(srad_abs[,2],tavg_abs[,2],col="red")

allpts <- rbind(cbind(occ$lon,occ$lat,rep(1,nrow(occ))),
                cbind(absences,rep(0,nrow(absences))))
colnames(allpts)[3]<-"pres"

bioclim <- worldclim_global(var = 'bioc', res = 10, path = ".")
bioclim_sample<-extract(bioclim,allpts[,1:2])

allpts_bioclim <- cbind(allpts,bioclim_sample)

predictor_vars <- names(allpts_bioclim)[4:ncol(allpts_bioclim)]
regres<-lapply(predictor_vars, function(predictor) {
  # Create formula for single regression
  formula <- as.formula(paste("pres", "~", predictor))

  # Fit the linear model
  lm_result <- glm(formula, data = allpts_bioclim,  family="binomial")

  # Return summary of the model
  summary(lm_result)
})

# Name each element in the list with the corresponding predictor variable
names(regres) <- predictor_vars

out<-sapply(regres, function(x)  1 - (x$deviance / x$null.deviance))
print(as.data.frame(out[rev(order(out))]))
tail(allpts_bioclim_plot$pres)

allpts_bioclim_plot<-allpts_bioclim
allpts_bioclim_plot$pres<-as.factor(allpts_bioclim_plot$pres)
ggplot(allpts_bioclim_plot, aes(x=wc2.1_10m_bio_19,y=wc2.1_10m_bio_4,shape=pres, col=pres)) +
  geom_point(size = 3, alpha=0.5) + scale_color_manual(values = c("0" = "red", "1" = "blue"))
ggplot(allpts_bioclim_plot, aes(x=wc2.1_10m_bio_19,y=wc2.1_10m_bio_3,shape=pres, col=pres)) +
  geom_point(size = 3, alpha=0.5) + scale_color_manual(values = c("0" = "red", "1" = "blue"))

cor_mat=cor(allpts_bioclim[,-c(1:3)],use="complete.obs",method = "spearman")
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

#pairs(bioclim_sample[1:100,-1])
#env=envelope(bioclim_sample[-1])
#plot(env)



