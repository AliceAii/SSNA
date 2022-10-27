rm(list=ls())

#######################################################################
#######Code for 7 Deliverable, One-mode Interactive Network########
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

#clean authors name, splite authors name
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

mat <- t(get.incidence(g2))

mat <- t(get.incidence(g2))
mat <- mat%*%t(mat)
dim(mat)

g<- graph.adjacency(mat, mode = "undirected", weighted = TRUE, diag = F)

#Gets centrality measures
cent<-data.frame(bet=betweenness(g, normalized=T, directed = FALSE)/max(betweenness(g, normalized=T, directed = FALSE)),eig=evcent(g)$vector, degree=degree(g, mode="total")) 
cent$name<-rownames(cent) #authors in this case
cent$pubnum <- diag(mat)[match(cent$name, rownames(mat))]
head(cent);tail(cent)

summary(cent)

actors<- a1[!duplicated(a1$AUNA), c("AUNA", "Affiliations")]
head(actors)
sort<- sort(degree(g2)[V(g2)$type==T], decreasing = T)
actors$pubnum<-sort
actors$pubnum<- as.numeric(actors$pubnum)

#Interactive Visualization
V(g)$label<-V(g)$name
V(g)$name<-1:length(V(g)) 

#create links
links<-as.data.frame(cbind(get.edgelist(g), E(g)$weight))

#Needs to be numeric
links$V1<-as.numeric(as.character(links$V1))

links$V2<-as.numeric(as.character(links$V2))
str(links)

links$V3<-round(as.numeric(as.character(links$V3)),3)
head(links)
colnames(links)<-c("source","target", "value")

#Counts begin at zero in computer programming
links[,1:2]<-(links[,1:2]-1)

dim(links)

colnames(links)<-c("source","target", "value")

names(actors)
V(g)$Affiliations<- actors$Affiliations[match(V(g)$label, actors$AUNA)]
V(g)$pubnum<- actors$pubnum[match(V(g)$label, actors$AUNA)]
V(g)$size<- cent$bet[match(V(g)$label, cent$name)]

summary(V(g)$size)
summary(V(g)$pubnum)
summary(V(g)$Affiliations)


nodes <- data.frame(name= c(paste("Author Name: ", V(g)$label,
                                  "; Publications: ", V(g)$pubnum, 
                                  "; Affiliations: ", V(g)$Affiliations)), 
                    pubnum = V(g)$pubnum,
                    size=abs(V(g)$size)) 
nodes$name<-as.character(nodes$name)

head(nodes);tail(nodes)

#create groups
nodes$group<-NA
nodes$group<-cut(nodes$pubnum, c(0,1.1,2.1,3.1,4.1), 
                      c('1 Paper','2 Paper', '3 Paper', '4 Paper'))
table(is.na(nodes$group))
table(nodes$group)

counts<-data.frame(table(nodes$group))
counts$labels <- paste(counts$Var1, ", N= ", counts$Freq, sep="")
nodes$groups <- counts$labels[match(nodes$group, counts$Var1)] 
head(nodes)

netviz2<-forceNetwork(Links = links, Nodes = nodes,
                     Source = 'source', Target = 'target',
                     NodeID = 'name',
                     Group = "groups",
                     charge = -30, # node repulsion
                     linkDistance = JS("function(d) { return d.linkDistance; }"),#JS("function(d){return d.value}"),
                     linkWidth = JS("function(d) { return Math.sqrt(d.value)*2; }"),
                     opacity = 0.8,
                     Value = "value",
                     Nodesize = 'size', 
                     radiusCalculation = JS("Math.sqrt(d.nodesize*2)+4"),
                     zoom = T, 
                     fontSize=14,
                     bounded= F,
                     legend= TRUE,
                     colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10)"))

HTMLaddons2 <- 
  "function(el, x) { 
d3.select('body').style('background-color', ' #f2c9d3 ')
d3.selectAll('.legend text').style('fill', 'white') 
 d3.selectAll('.link').append('svg:title')
      .text(function(d) { return 'Number of Co-publish time: ' + d.value ; })
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
netviz2$x$links$linkDistance <- (1/links$value)*400
netviz2$x$links$value <- links$value
onRender(netviz2, HTMLaddons2)

