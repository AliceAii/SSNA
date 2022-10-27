rm(list=ls())

#######################################################################
#######Code for Sixth Deliverable, Two-mode Interactive Network########
#######################################################################

#load packages
library(splitstackshape)
library(igraph)
library(networkD3)
library(magrittr)
library(htmlwidgets)
library(htmltools)

#read table
setwd("~/Data/SSNA")
a<-read.csv("scopus1.csv")
a

#clean authors name, split authors name
a$Authors<- gsub(" Jr.,", "", 
                 gsub(" II.,", "", 
                      gsub(" Jr.", "",
                           gsub(" M.S.", "",
                                gsub(" M.S.,", "",
                                     gsub(" II.", "", a[,1]))))))
a$Authors<-gsub("\\.,", ";", a$Authors)

authorname<- as.data.frame(a[,1])
colnames(authorname)<- "authorname"

AUNA<-cSplit(authorname, splitCols="authorname", 
             sep=";", direction="wide", drop=FALSE)
AUNA<-AUNA[,-1]
mat_AUNA<- as.matrix(AUNA)
mat_AUNA<- tolower(mat_AUNA)
dim(mat_AUNA) #the same number with AUID

#link publication with author name
EID_AUNA <-cbind(a$EID, mat_AUNA)
mat_EID_AUNA <- as.matrix(EID_AUNA)

#create edgelist_two_mode 
edgelist_two_mode <- cbind(mat_EID_AUNA[,1], c(mat_EID_AUNA[,-1]))
edgelist_two_mode <- edgelist_two_mode[!is.na(edgelist_two_mode[,2]), ]
df_authorname<- as.data.frame(edgelist_two_mode)
colnames(df_authorname)<- c('EID', 'AUNA')

#prepare whole dataset
a1<-merge(a, df_authorname,by="EID")
names(a1)

#Retrieve the author-publication connections to be saved under a graph "g2"
g2<- graph.data.frame(a1[, c("AUNA", "EID", "Cited.by", "Title")])

#Add the two-mode structure to the graph "g"
V(g2)$type <- V(g2)$name %in% a1[,c("AUNA")] 
table(V(g2)$type)[2]
i<-table(V(g2)$type)[2]

sort<- sort(degree(g2)[V(g2)$type==T], decreasing = T)

#Gets centrality measures
cent<-data.frame(bet=betweenness(g2, normalized=T, directed = FALSE)/max(betweenness(g2, normalized=T, directed = FALSE)),eig=evcent(g2)$vector, degree=degree(g2, mode="total")) 
cent$name<-rownames(cent) #author name in this case
head(cent);tail(cent)

cent$bet[1:i]<-cent$bet[1:i]/max(cent$bet[1:i])
cent$eig[1:i]<-cent$eig[1:i]/max(cent$eig[1:i])
summary(cent[1:i,])
cent$bet[(i+1):nrow(cent)]<-cent$bet[(i+1):nrow(cent)]/max(cent$bet[(i+1):nrow(cent)])
cent$eig[(i+1):nrow(cent)]<-cent$eig[(i+1):nrow(cent)]/max(cent$eig[(i+1):nrow(cent)])
summary(cent[(i+1):nrow(cent),])

#Gets time invariant actor attributes
actors<- a1[!duplicated(a1$AUNA), c("AUNA", "Affiliations")]
actors$pubnum<-sort
actors$pubnum<- as.numeric(actors$pubnum)

#Interactive Visualization
V(g2)$label<-V(g2)$name
V(g2)$name<-1:length(V(g2)) 

#Prepare two datasets (nodes,links) for HTML visualizations (like NetworkD3) 
#Gets edgelist from graph, also any other attribute at the edge level to be included in the mapping
links<-as.data.frame(cbind(get.edgelist(g2), E(g2)$Cited.by, E(g2)$Title))
#Needs to be numeric
links$V1<-as.numeric(as.character(links$V1))
links$V2<-as.numeric(as.character(links$V2))
links$V3<-as.numeric(as.character(links$V3))
str(links)

colnames(links)<-c("source","target", "value", "Title")

#Counts begin at zero in computer programming
links[,1:2]<-(links[,1:2]-1)
dim(links)

#Adding attributes at the actor level
names(actors)
V(g2)$Affiliations<- actors$Affiliations[match(V(g2)$label, actors$AUNA)]
V(g2)$pubnum<- actors$pubnum[match(V(g2)$label, actors$AUNA)]
V(g2)$eig<- cent$eig[match(V(g2)$label, cent$name)]
V(g2)$degree<- cent$degree[match(V(g2)$label, cent$name)]
V(g2)$size<- cent$bet[match(V(g2)$label, cent$name)]

V(g2)$Affiliations[1:i] <- ifelse(is.na(V(g2)$Affiliations[1:i]), "Affiliationsn_NA", V(g2)$Affiliations[1:i])
V(g2)$pubnum[1:i] <- ifelse(is.na(V(g2)$pubnum[1:i]), "pubnum_NA", V(g2)$pubnum[1:i])

summary(V(g2)$pubnum)
summary(V(g2)$eig)
summary(V(g2)$degree)
summary(V(g2)$size)

#create nodes dataset
nodes <- data.frame(name= c(paste(" Author Name: ", V(g2)$label[1:i], 
                                "; Publication number: ", V(g2)$pubnum[1:i], 
                                "; Affiliations: ", V(g2)$Affiliations[1:i], sep=""), 
                            V(g2)$label[(i+1):length(V(g2)$label)]), 
                    pubnum = V(g2)$pubnum,
                    eigencentrality = abs(V(g2)$eig), 
                    degree = abs(V(g2)$degree),
                    size=abs(V(g2)$size))
nodes$name<-as.character(nodes$name)

head(nodes); tail(nodes)

#create nodes groups based on the number of publication 
nodes$group<-NA
nodes$group[1:i]<-cut(nodes$pubnum[1:i], c(0,1.1,2.1,3.1,4.1), 
                     c('1 Paper','2 Paper', '3 Paper', '4 Paper'))
table(is.na(nodes$group))
table(nodes$group)

counts<-data.frame(table(nodes$group))
counts$labels <- paste(counts$Var1, ", N= ", counts$Freq, sep="")
nodes$groups <- counts$labels[match(nodes$group, counts$Var1)] 
head(nodes)

#get network visulization
netviz<-forceNetwork(Links = links, Nodes = nodes,
                     Source = 'source', Target = 'target',
                     NodeID = 'name',
                     Group = "groups", # color nodes by group calculated earlier
                     charge = -30, # node repulsion
                     linkDistance = JS("function(d) { return d.linkDistance; }"),#JS("function(d){return d.value}"),
                     linkWidth = JS("function(d) { return Math.sqrt(d.value)*2; }"),
                     opacity = 0.8,
                     Value = "value",
                     Nodesize = 'size', 
                     radiusCalculation = JS("Math.sqrt(d.nodesize*30)+4"),
                     zoom = T, 
                     fontSize=14,
                     bounded= F,
                     legend= TRUE,
                     colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10)"))

HTMLaddons <- 
  "function(el, x) { 
d3.select('body').style('background-color', ' #212f3d ')
d3.selectAll('.legend text').style('fill', 'white') 
 d3.selectAll('.link').append('svg:title')
      .text(function(d) { return 'Grade course: ' + d.value + ', Campus: ' + d.campus ; })
  var options = x.options;
  var svg = d3.select(el).select('svg')
  var node = svg.selectAll('.node');
  var link = svg.selectAll('link');
  var mouseout = d3.selectAll('.node').on('mouseout');
  function nodeSize(d) {
    if (options.nodesize) {
      return eval(options.radiusCalculation);
    } else {
      return 6;
    }
  }

  
d3.selectAll('.node').on('click', onclick)

  function onclick(d) {
    if (d3.select(this).on('mouseout') == mouseout) {
      d3.select(this).on('mouseout', mouseout_clicked);
    } else {
      d3.select(this).on('mouseout', mouseout);
    }
  }

  function mouseout_clicked(d) {
    node.style('opacity', +options.opacity);
    link.style('opacity', +options.opacity);

    d3.select(this).select('circle').transition()
      .duration(750)
      .attr('r', function(d){return nodeSize(d);});
    d3.select(this).select('text').transition()
	
      .duration(1250)
      .attr('x', 0)
      .style('font', options.fontSize + 'px ');
  }

}
"
netviz$x$links$linkDistance <- (1/links$value)*1500
netviz$x$links$Title <- links$Title
onRender(netviz, HTMLaddons)
