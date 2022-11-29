rm(list = ls())

#load packages
library(spatialreg)
library(spdep)
library(igraph)
library(classInt)
library(RColorBrewer)

#load Fridenship dataset
idfriend <- "1dwX4kKlx-ctkU0JyH74p3r1w-jJdTAqi"
friendshiplazega <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", idfriend))
friendshiplazega <- as.matrix(friendshiplazega)
friendshiplazegac <- friendshiplazega
str(friendshiplazega)
dim(friendshiplazega)

#Reading attributes also provided by Lazega
idattributes <- "1e0GtrRS5PFFNdnd1e4fJcjeuBZ6deF7g"
datattrout <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", idattributes))

#Normalize the dataset. Transforming the matrices to spatial form (row normalized) as shown in following equation 
friendshiplazega <-friendshiplazega /rowSums(friendshiplazega)
summary(rowSums(friendshiplazega))
#6 Na's means there are 6 layers didn't have friends with anybody.

# Replacing potential NAN to zeros as shown in equation (22)
friendshiplazega[is.na(friendshiplazega)]<-0
listwAd<-mat2listw(friendshiplazega)

#we are allowing the isolated individuals in the weighted list,zero.policy=T
summary(listwAd, zero.policy=T)
str(listwAd) #data structure, two attributes, one is connections, another is weights
 
#Q1 Testing influence. Is there evidence of stronger dependence in the per hour rate compared to the fees brought in 1990 amounts?
moran.test(datattrout$HrRATE90,listwAd, zero.policy=TRUE)
moran.test(datattrout$FeesCollec90,listwAd, zero.policy=TRUE)
#morean test for the individual hour rate, p<0.001, Moran I is positive 0.521, indicating that we have enough evidence to conclude that the our outcomes consistently change with our peers, either our outcomes and peers outcomes all over the mean or all under the mean.
#morean test for the individual fees brought in 1990, p<0.001, Moran I is positive 0.455, indicating that we have enough evidence to conclude that the fees we brought are similar with our peers, we either have consistently fees brought in 1990 over the mean, or we and our peers fees brought in 1990 all under the mean value.
#And compared with peers influence on individual's hour rate, our peers have less impact on fees brought in 1990, since 0.521>0.455. Therefore, there is evidence of stronger dependence in the per hour rate compared to the fees brought in 1990 amounts.

#########
#Visualizing Local Moran’s I clusters and outliers 
#this is for the empty model
#########
mtHrRATE90 <- moran.test(datattrout$HrRATE90, listwAd, zero.policy=TRUE)
label_x = "Individual Per Hour Rate"
label_y = "Lagged Peer's Per Hour Rate"

mpHrRATE90 <- moran.plot(datattrout$HrRATE90, listwAd, zero.policy=T,
                 labels=datattrout$id, xlab = label_x, ylab = label_y)
title(main="Moran’s Plot for Per Hour Rate", cex.main=2, col.main="grey11", 
      font.main=2, sub=paste("Plot includes 71 lawyers (Moran’s I = ", 
                             round(mtHrRATE90$ estimate[1], 3), ", p < .0001)", sep=""), cex.sub=1.15, col.sub="grey11", font.sub=2,)
#x axis means individual's per hour rate, y axis means our peers' mean per hour rate
#the dots with label means they may not meaningful under the HH/LH/LL/HL region, they may indicate some special mechanism, you could explore more.

mtFeesCollec90 <- moran.test(datattrout$FeesCollec90, listwAd, zero.policy=TRUE)
label_x = "Individual Fees Brought in 1990"
label_y = "Lagged Peer's Fees Brought in 1990"

mpFeesCollec90 <- moran.plot(datattrout$FeesCollec90, listwAd, zero.policy=T,
                         labels=datattrout$id, xlab = label_x, ylab = label_y)
title(main="Moran’s Plot for Fees Brought in 1990", cex.main=2, col.main="grey11", 
      font.main=2, sub=paste("Plot includes 71 lawyers (Moran’s I = ", 
                             round(mtFeesCollec90$ estimate[1], 3), ", p < .0001)", sep=""), cex.sub=1.15, col.sub="grey11", font.sub=2,)


#Q2 Did you have to modify the empty model to correct for spatial dependence? That is,  was the dependence issue addressed with the empty regression models? Explain the likely mechanisms behind this process.
#SAR (Simultaneous Autoregressive), uses a regression on the values from the other neighboring areas, as captured by the weight matrix (wij), also used in Moran's I procedures, to account for spatial/outcome dependence.Since this outcome dependece is the source of RSODA, the inclusion of λwijy in the model should yeild model residuals that are i.i.d.
#Procedually, we may estimate the null or empty model relying on SAR,then we may extract the model residuals and test for whether their tend to co-vary in the same direction as their neighbors' residuals do, with respect to the residual mean.

##SAR procedures for social dependence for outcome 1 - Per Hour Rate
Hrrate90Ad <- spautolm(formula = HrRATE90 ~ 1, data = data.frame(datattrout), listw = listwAd)
summary(Hrrate90Ad)
#lambda λ = 0.904, means this is the average value of per hour rate of our peers. p-value<0.001, the λ is significant, indicating that the 1 unit change in our peers' hour rate, our rate will consistently change 0.903 unit under the same direction.

#Saving the residuals
hNULL<-residuals(Hrrate90Ad)
summary(hNULL)
#Testing residuals for sp dependence
moran.test(hNULL,listwAd, zero.policy=TRUE)
#since the p-value for Moran I test is 0.430 > 0.05, indicating that the residual is i.i.d. Now we can say that the dependence issue addressed with the empty regression friendship models in terms of per hour rate.

##SAR procedures for social dependence for outcome 2 - Fees Brought in 1990
fees90Ad <- spautolm(formula = FeesCollec90 ~ 1, data = data.frame(datattrout), listw = listwAd)
summary(fees90Ad)
#lambda λ = 0.871, means this is the average value of fees brought in 1990 of our peers. p-value<0.001, the λ is significant, indicating that the 1 unit change in our peers' fees brought in 1990, our fees will consistently change 0.871 unit under the same direction.

#Saving the residuals
hNULL2<-residuals(fees90Ad)
#Testing residuals for sp dependence
moran.test(hNULL2, listwAd, zero.policy=TRUE)
#since the p-value for Moran I test is 0.446 > 0.05, indicating that the residual is i.i.d. Now we can say that the dependence issue addressed with the empty regression friendship models in terms of fees brought in 1990.

#Q3 Are there any lagged indicators that you may be interested in testing?

#First run the feature importance analysis
#install.packages("Boruta")
library(Boruta)

set.seed(123)
boruta.train <- Boruta(HrRATE90 ~. , data.frame(datattrout), doTrace = 2)
print(boruta.train)

plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])
names(lz) <- colnames(boruta.train$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.7)
#the plots show that yearsfirm and partner status are importance features.


#Getting socially lagged indicators- whether my friends years at firm will impact my salary
datattrout$lag.yearsfirmAd <- lag.listw(listwAd, datattrout$yearsfirm, zero.policy=T, na.action=na.omit)
#Network visualization
Hrrate90Adwithlag1 <- spautolm(formula = HrRATE90 ~ lag.yearsfirmAd , data = data.frame(datattrout), listw = listwAd)
summary(Hrrate90Adwithlag1)
#interpretation: Whether my peers years at firm is significantly associated with my salary with p-value<0.001. For each additional year my friend has been with the company, my salary will increase by $5.46$.

#Getting socially lagged indicators- whether my friends are partners 
datattrout$lag.partnersAd <- lag.listw(listwAd, datattrout$partner, zero.policy=T, na.action=na.omit)
#Network visualization
Hrrate90Adwithlag2 <- spautolm(formula = HrRATE90 ~ lag.partnersAd , data = data.frame(datattrout), listw = listwAd)
summary(Hrrate90Adwithlag2)
#interpretation: Whether my peers are partners status is significantly associated with my salary with p-value<0.001. If 100% of my friends are all partners, my salary goes up by 94.78$. If 60% of my friends are partners, my salary goes up by 94.78*0.6=56.87$.

#Q4 After removing isolates:a. How many higher order neighbors did you find? b.Does this change by outcome?

#How are disconnected unit influencing these results and how to remove them?
#########
sub_listNAd <- subset(listwAd[[2]], subset=card(listwAd[[2]])> 0)
sub_listNAd
summary(sub_listNAd)
str(sub_listNAd) #no weights, only have connections
#add weights,using style="W" (the normalized version of weights)
sub_listWAd <- nb2listw(sub_listNAd, glist=NULL, style="W", zero.policy=NULL) 
#making this a spatial points dataframe
sub_datattrout <- subset(datattrout, subset=card(listwAd[[2]]) > 0)
dim(sub_datattrout)

#b. Does this change by outcome? Test with new sub_datatrout
#Q5 How do the Moran's I estimates change when removing isolate?
moran.test(sub_datattrout$HrRATE90,sub_listWAd, zero.policy=TRUE)
moran.test(sub_datattrout$FeesCollec90,sub_listWAd, zero.policy=TRUE)
#The results changed after we removed the isolated samples. Specifically, the Moran I for per hour rate is 0.684 with p-value<0.001, the Moran I for the fees brought in 1990 is 0.528 with p-value<0.001. The Moran I values all increased after removing isolated individuals.

# How many neighbors of neighbors (higher order neighboring structures) do we have to account for to establish a data driven selection of neighboring structures?
#a.How many higher order neighbors did you find?
plot.spcor(sp.correlogram(sub_listNAd, sub_datattrout$HrRATE90, order = 6, method = "I", zero.policy=T), xlab = "Social lags", main = "Social correlogram: Autocorrelation with CIs")

plot.spcor(sp.correlogram(sub_listNAd, sub_datattrout$FeesCollec90, order = 6, method = "I", zero.policy=T), xlab = "Social lags", main = "Social correlogram: Autocorrelation with CIs")
#The first 2 order is enough for us to establish a data driven selection of neighboring structures in terms of per hour rate, and first 2 order for fees brought in 1990 as well.



