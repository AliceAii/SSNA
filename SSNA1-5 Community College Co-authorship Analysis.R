rm(list=ls())

#load packages
library(splitstackshape)
library(igraph)
library(dplyr)
#install.packages('tinytex') 
#tinytex::install_tinytex()

#read table
setwd("~/Data/SSNA")
a<-read.csv("scopus1.csv")
a
#Keeping only column with author relationships
authors<-as.data.frame(a[,1])
head(authors)
colnames(authors)<-"AU"

#As can be seen in the file the separator of interest is :
a1<-cSplit(authors, splitCols = "AU", sep = ",", direction = "wide", drop = FALSE)
#Here we just drop the original first column
a1<-a1[,-1]
class(a1)

#read it as a matrix
mat <- as.matrix(a1)
mat

mat<-tolower(mat)

#dimision result 56, 7
dim(mat)

# create edgelist _ all combinations except first author with itself
a1<-mat
edgelist1<-matrix(NA, 1, 2)#empty matrix two columns
for (i in 1:(ncol(a1)-1)) {
  edgelist11 <- cbind(a1[, i], c(a1[, -c(1:i)]))
  edgelist1 <- rbind(edgelist1,edgelist11)
  edgelist1<-edgelist1[!is.na(edgelist1[,2]),]
  edgelist1<-edgelist1[edgelist1[,2]!="",]
}
dim(edgelist1)  #231, 2

#plot author relationship
plot(graph.edgelist(edgelist1))

g<- graph.edgelist(edgelist1, directed = FALSE)


#create weight add weight to E(g.c)
E(g)$weight	<- 1#must step
g.c <- simplify(g)
E(g.c)$weight 

links<-as.data.frame(cbind(get.edgelist(g.c), E(g.c)$weight))

head(links,20)
dim(links)
links$V3<-as.numeric(links$V3)

links<- links[order(links$V3, decreasing=T),]
links
              
###########################################
###		Code for second deliverable     ###
###########################################

#Add publication as first column
a1 <- cbind(a$Source.title, a1)
              
mat <- as.matrix(a1)

edgelist_two_mode <- cbind(mat[, 1], c(mat[, -1]))
edgelist_two_mode
edgelist_two_mode <- edgelist_two_mode[!is.na(edgelist_two_mode[,2]), ]
g2<- graph.edgelist(edgelist_two_mode[, 2:1], directed = FALSE)  #this is from author to publication

V(g2)$type <- V(g2)$name %in% edgelist_two_mode[ ,1]
g2

degree(g2)[V(g2)$type==F]

#sort(degree(g2)[V(g2)$type==F], decreasing =T)
#Authors
sort1 <-sort(degree(g2)[V(g2)$type==F], decreasing =T)
head_print(sort1, 10)

#does this mean 146 authors, 36 journals
table(V(g2)$type)

#Transformation 
links2<-as.data.frame(cbind(get.edgelist(g2)))
head(links2)

#Publications
sort(degree(g2)[V(g2)$type==T], decreasing = T)
head(sort(degree(g2)[V(g2)$type==T], decreasing = T))


###########################################
###		Code for third deliverable     ###
###########################################
#for deliverable 3
# Transformations
# mat2 <- get.adjacency(g2) #squre matrix

mat2 <- get.incidence(g2) #rectangular matix
g2
dim(mat2) #146 human, 36 publications

mat2_to_1<-mat2%*%t(mat2)
mat2_to_1
dim(mat2_to_1)
diag(mat2_to_1) #the number of authors and the number of publications. can be used in the first deliverable

g

g.c

g2 #(182=146+36)

###########################################
###		Code for forth deliverable     ###
###########################################

#identify the top six most central human actors
#one-mode 
evcent(g)$vector
plot(g)

dta1<-data.frame(ev=evcent(g)$vector, deg=degree(g), bet=(betweenness(g, normalized=F)/max(betweenness(g, normalized = F))), clo=closeness(g)) 
dta1

#two mode
evcent(g2)$vector
plot(g2)

dta2<-data.frame(ev=evcent(g2)$vector, deg=degree(g2), bet=(betweenness(g2, normalized=F)/max(betweenness(g2, normalized = F))), clo=closeness(g2)) 
dta2

#centrality table
one_mode_df <-read.table( header=T, sep='|', text='
Author|Total Degree|Eigenvector|Betweenness|Closeness
Hu S.|14|1.000|0.875|0.111
Park T.|14|1.000|0.875|0.111
Bertrand Jones T.|7|0.672|0.000|0.071
Woods C. S.|7|0.672|0.000|0.071
Tandberg D.|4|0.392|0.000|0.071
Hu X.|4|0.306|0.000|0.071' )
print(one_mode_df)

two_mode_df <-read.table( header=T, sep='|', text='
Author|Total Degree|Eigenvector|Betweenness|Closeness
Hu S.|4|0.608|0.673|0.008
Park T.|4|0.608|0.673|0.008
Bertrand Jones T.|2|0.395|0.055|0.006
Woods C. S.|2|0.395|0.055|0.006
Yonah Meiselman A.|1|0.250|0.000|0.005
Uretsky M. C.|1|0.250|0.000|0.005' )
print(two_mode_df)

###########################################
###		Code for fifth deliverable     ###
###########################################


cent<-data.frame(bet=betweenness(g, normalized=F)/max(betweenness(g, normalized=F)),eig=evcent(g)$vector)
cent

rownames(cent)<-rownames(dta1)
cent

#add residual, add coef
residuals<-lm(eig~bet,data=cent)$residuals
cent<-transform(cent,res=residuals)
coeffs<-as.data.frame(coef(lm(eig~bet,data=cent)))

summary(lm(cent$bet~cent$eig))

library(ggplot2) 

#code for Figure 1 eig~bet
pdf("Figure1_Key Actor Analysis for Co-authorship.pdf ", width = 15.0, height = 15.0)

figure1<-ggplot(cent,aes(x=bet,y=eig, label=rownames(cent),colour=eig, size=(1/abs(res))))+xlab("Betweenness Centrality")+ylab("Eigenvector Centrality")

figure1 + 
  geom_point() + 
  geom_text(hjust=1.5, vjust=1)+
  labs(title="Key Actor Analysis for Co-authorship") + 
  geom_abline(intercept = coeffs[1,], slope = coeffs[2,],colour = "red", size = 2,alpha=.25) + 
  theme(legend.position = "none")
dev.off()

#code for figure 2
set.seed(47)
l<-layout.fruchterman.reingold(g, niter=5000)
V(g)$name<-rownames(dta1)
V(g)$size<-abs((cent$bet)/max(cent$bet))*15#The divisor is the highest betweenness
nodes<-V(g)$name 

#set 5% threshold
x<-quantile(cent$eig, .95)
nodes[which(abs(cent$eig)<(x))]<-NA

table(is.na(nodes))

pdf("Figure2_Key Actor Analysis for Co-authorship.pdf ", width = 15.0, height = 15.0)

plot(g,layout=l,vertex.label=nodes, vertex.label.dist=0.25,
     vertex.label.color="red",edge.width=1)
title(main="Key Actor Analysis of Co-authorship", sub="Key actors weighted by BC, names are top 5% EC", col.main="black", col.sub="black", cex.sub=1.2,cex.main=2,font.sub=2)

dev.off()
