rm(list = ls())

#load packages
library(spatialreg)
library(spdep)
library(igraph)
library(classInt)
library(RColorBrewer)
library(splitstackshape)
library(magrittr)

#read table as my datattrout
setwd("~/Data/SSNA")
a<-read.csv("scopus1.csv")
a

dim(a)
authors<-as.data.frame(a[,1])
head(authors)
colnames(authors)<-"AU"

#split authors name
a1<-cSplit(authors, splitCols = "AU", sep = ",", direction = "wide", drop = FALSE)

#Here we just drop the original first column
a1<-a1[,-1]
#y <- data.frame(table(unlist(a1)))

#read it as a matrix
mat <- as.matrix(a1)
mat<-tolower(mat)
dim(mat)

#create author-publication two-mode edgelist
a1 <- cbind(a$EID, a1)

mat <- as.matrix(a1)

edgelist_two_mode <- cbind(mat[, 1], c(mat[, -1]))
edgelist_two_mode
edgelist_two_mode <- edgelist_two_mode[!is.na(edgelist_two_mode[,2]), ]
g2<- graph.edgelist(edgelist_two_mode[, 2:1], directed = TRUE)  #this is from author to publication

V(g2)$type <- V(g2)$name %in% edgelist_two_mode[ ,2]
g2

mat2 <- t(get.incidence(g2)) #rectangular matix
dim(mat2) #158 human, 56 publications

mat2_to_1<-mat2%*%t(mat2)
mat2_to_1
dim(mat2_to_1)
summary(diag(mat2_to_1)) #the number of authors and the number of publications.

y<-diag(mat2_to_1)
diag(mat2_to_1) <- 0

#Normalize the dataset. Transforming the matrices to spatial form (row normalized) as shown in following equation 
coauthorship <- mat2_to_1 /rowSums(mat2_to_1)
summary(rowSums(coauthorship))
#9 Na's means there are 9 authors are single authors.

# Replacing potential NA to zeros as shown 
coauthorship[is.na(coauthorship)]<-0
listwCoau<-mat2listw(coauthorship)
summary(rowSums(coauthorship))

#we are allowing the isolated individuals in the weighted list,zero.policy=T
summary(listwCoau, zero.policy=T)
str(listwCoau) #data structure, two attributes, one is connections, another is weights

#Q1 Test for outcome dependence using author's number of publication as a function of her/his coauthors' number of publications
moran.test(y,listwCoau, zero.policy=TRUE)
#Moran I test for the authors' number of publication:Moran I is 0.309 with p value < 0.001, indicating that we have enough evidence to conclude that author's number of publication consistently change with their coauthors' outcomes, either their publications all over the mean or all under the mean.

y1 <- data.frame(id=names(y), y=y)
head(y1)

#Q1 After Removing isolates
sub_listNCoau <- subset(listwCoau[[2]], subset=card(listwCoau[[2]])> 0)
sub_listNCoau
summary(sub_listNCoau)
str(sub_listNCoau) #no weights, only have connections
#add weights,using style="W" (the normalized version of weights)
sub_listWCoau <- nb2listw(sub_listNCoau, glist=NULL, style="W", zero.policy=NULL) 
#making this a spatial points dataframe
sub_sortAu <- subset(y1, subset=card(listwCoau[[2]]) > 0)
dim(sub_sortAu)

moran.test(sub_sortAu$y,sub_listWCoau, zero.policy=TRUE)
#After removing the isolated authors, the Moran I is 0.245 with p <0.001, indicating that coauthor's number of publications have impact on author's number of publication. Compared with the above analysis, the Moran I value decreased a little bit.

#Q2 Also visualize this influence level, what may be the mechanism to explain these clusters and outliers?
mtpubnum <- moran.test(y, listwCoau, zero.policy=TRUE)
label_x = "Author's Number of Publications"
label_y = "Lagged Coauthors' Number of Publications"

mppubnum <- moran.plot(y, listwCoau, zero.policy=T,
                         labels=names(y), xlab = label_x, ylab = label_y)
title(main="Moran’s Plot for the Number of Publication", cex.main=2, col.main="grey11", 
      font.main=2, sub=paste("Plot includes 158 Authors (Moran’s I = ",
                             round(mtpubnum$ estimate[1], 3), ", p < .0001)", sep=""), cex.sub=1.15, col.sub="grey11", font.sub=2,)
#x axis means author's number of publications, y axis means coauthors' mean number of publication
#The dots with label means they may not meaningful under the HH/LH/LL/HL region, they may indicate some special mechanism, ....you could explore more.

#Q3 Test how many higher order neighbors should we account for in this framework?
plot.spcor(sp.correlogram(sub_listNCoau, sub_sortAu$y, order = 3, method = "I", zero.policy=T), xlab = "Social lags", main = "Social correlogram: Autocorrelation with CIs")
#As it shown in the plot, we should only account for the first 2 order to in terms of number of publications.

#Q4, Q5 If needed create a higher order weight matrix and test for autocorrelation. Finally, test whether an empty model may be enough to address social dependence.
PubnumCoau <- spautolm(formula = y ~ 1, data = data.frame(y), listw = listwCoau)
summary(PubnumCoau)

# lambda = 0.209, means this is the average number of publications of coauthors. Since p value < 0.05, the lambda is significant, indicating that 1 more paper coauthor published, the author's publication number will change 0.209 unit under the same direction.

#Saving the residuals
hNULL<-residuals(PubnumCoau)
summary(hNULL)
#Testing residuals for sp dependence
moran.test(hNULL,listwCoau, zero.policy=TRUE)
#Since the p-value for Moran I test is 0.517 > 0.05, indicating that the residual is i.i.d. Now we can say that the dependence issue addressed with the empty regression coauthorship models in terms of number of publications.


