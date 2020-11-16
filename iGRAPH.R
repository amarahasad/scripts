                       #########################################      Set CYTOSCAPE  ##########################################
### PRE_Process
# 1
 # # open terminal -> go to cytoscape directory -> and type this command (java -Xmx512M -jar cytoscape.jar -p plugins)
 # Double-click on the icon created by the installer or by running cytoscape.sh from the command line (Linux or Mac OS X) or double-clicking cytoscape.bat (Windows). Alternatively, you can pass the .jar file to Java directly using the command java -Xmx512M -jar cytoscape.jar -p plugins. The -Xmx512M flag tells java to allocate more memory for Cytoscape and the -p plugins option tells cytoscape to load all of the plugins in the plugins directory. Loading the plugins is important because many key features like layouts, filters and the attribute browser are included with Cytoscape as plugins in the plugins directory.

# 2
# CytoscapeRPC installation consists of copying the plugin's jar file to the 'Plugins' directory of your installed Cytoscape application, and restarting Cytoscape. You should see CytoscapeRPC in Cytoscape's Plugins menu; click to activate the plugin, which starts the XMLRPC server, and which listens for commands from R (or elsewhere). The default and usual choice is to communicate over port 9000 on localhost. You can choose a differenet port from the Plugins->CytoscapeRPC activate menu item. If you choose a different port, be sure to use that port number when you call the RCytoscape constructor - the default value on in the constructor argument list is also 9000.

                       #########################################      Set R CYTOSCAPE  ##########################################

# source("https://bioconductor.org/biocLite.R")  -----------------------------if not installed already
# biocLite("RCytoscape")
# browseVignettes("RCytoscape")
#install.packages("igraph")
                       
library ('RCytoscape')
cy = CytoscapeConnection ()
pluginVersion (cy)
rm(list = ls())

library("igraph")
library("plyr")

datafile <- read.csv("pair_up_over.csv", header = TRUE, sep = ",")
dat<-datafile[,c(3,4,6)]
#sam<-dat[sample(nrow(dat), 201), ]
dataSet<-dat
colnames(dataSet)=c('V1','V2','V3')
gD <- simplify(graph.data.frame(dataSet, directed=FALSE))
vcount(gD)
ecount(gD)
degAll <- degree(gD, v = V(gD), mode = "all")
betAll <- betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
rm(betAll)
dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
gD <- set.vertex.attribute(gD, "degree", index = V(gD), value = degAll)
gD <- set.vertex.attribute(gD, "betweenness", index = V(gD), value = betAll.norm)
summary(gD)
F1 <- function(x) {data.frame(V4 = dsAll[which(V(gD)$name == as.character(x$V1)), which(V(gD)$name == as.character(x$V2))])}
dataSet.ext <- ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))
colnames(dataSet)=c('V1','V2','V3')
F1 <- function(x) {data.frame(V4 = dsAll[which(V(gD)$name == as.character(x$V1)), which(V(gD)$name == as.character(x$V2))])}
dataSet.ext <- ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))
gD <- set.edge.attribute(gD, "weight", index = E(gD), value = 0)
gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
summary(gD)
gD.cyt <- igraph.to.graphNEL(gD)
gD.cyt <- initNodeAttribute(gD.cyt, 'degree', 'numeric', 0)
gD.cyt <- initNodeAttribute(gD.cyt, 'betweenness', 'numeric', 0)
gD.cyt <- initEdgeAttribute (gD.cyt, "weight", 'integer', 0)
gD.cyt <- initEdgeAttribute (gD.cyt, "similarity", 'numeric', 0)
gDCW <- new.CytoscapeWindow("SAMPLE NETWORK", graph = gD.cyt, overwriteWindow = TRUE)
displayGraph(gDCW)
cy <- CytoscapeConnection()
hlp <-getLayoutNames(cy)
setLayoutProperties (gDCW, hlp[18], list (edge_attribute = 'similarity', iterations = 1000))
layoutNetwork(gDCW, hlp[18])
setDefaultBackgroundColor(gDCW, '#FFFFFF')
setDefaultEdgeColor(gDCW, '#CDC9C9')
setDefaultEdgeLineWidth(gDCW, 4)
setDefaultNodeBorderColor(gDCW, '#000000')
setDefaultNodeBorderWidth(gDCW, 3)
setDefaultNodeShape(gDCW, 'ellipse')
setDefaultNodeColor(gDCW, '#87CEFA')
setDefaultNodeSize(gDCW, 60)
setDefaultNodeFontSize(gDCW, 20)
setDefaultNodeLabelColor(gDCW, '#000000')
redraw(gDCW)
setNodeColorRule(gDCW, 'degree', c(min(degAll), mean(degAll), max(degAll)), c('#F5DEB3', '#FFA500', '#FF7F50', '#FF4500', '#FF0000'), mode = 'interpolate')
setNodeSizeRule(gDCW, 'betweenness', c(min(betAll.norm), mean(betAll.norm), max(betAll.norm)), c(30, 45, 60, 80, 100), mode = 'interpolate')
setEdgeColorRule(gDCW, 'weight', c(min(as.numeric(dataSet.ext$V3)), mean(as.numeric(dataSet.ext$V3)), max(as.numeric(dataSet.ext$V3))), c('#FFFF00', '#00FFFF', '#00FF7F', '#228B22', '#006400'), mode='interpolate')
redraw (gDCW)