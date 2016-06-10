##########
#plotTheNelGraph:plot Using graphNEL
#Params
#gR: the graph
#factorsE: the factors 
#' @import graph
#' @import org.Hs.eg.db 
#' @import Rgraphviz
#' @import AnnotationDbi
#private
plotPathway2Colors <- function(pathway.i, subclass, name = "A pathway",  
                               twocolors = c("lightgreen", "#FEB24C"), twoshapes = c("box", "box"),
                               twofontcolor = c("black", "black") , fontsize = 20, numberofchar = 15, showEdgeData = TRUE){
  
  entrezOrganism = "org.Hs.eg"
  
  entrez2Sym <- org.Hs.egSYMBOL2EG
  mapped_genes <- AnnotationDbi::mappedkeys(entrez2Sym)
  entrez2sym_l <- AnnotationDbi::as.list(entrez2Sym[mapped_genes])
  
  find_sym <- function(x){
    posi <- which( entrez2sym_l == x)
    if(length(posi)==0){
      posi <- x
    }else{
      posi <- names(posi)[1]
    }
    posi     
  }
  
  entrezID <- graph::nodes(pathway.i)
  #entrezID <- sub("hsa:","", entrezID)
  
  symbols <-lapply(entrezID, find_sym )
  symbols <- substr(symbols, 1, numberofchar)
  symbols <- unlist(symbols)
  graph::nodes(pathway.i) <- symbols
  edgeLabels <- unlist(edgeWeights(pathway.i))
  edgesFromTo <- names(edgeLabels)
  
  edgeLabels <- unlist(edgeData(pathway.i))
  
  nodeType <- 1 + (graph::nodes(pathway.i) %in% subclass)
  nA = Rgraphviz::makeNodeAttrs(pathway.i, 
                                height = "1", #width = "1",
                                fillcolor = twocolors[nodeType], shape = twoshapes[nodeType], 
                                fontcolor = twofontcolor[nodeType], fontsize = fontsize )
  eAttrs <- list()
  edgesFromTo <- gsub("\\.", "~", edgesFromTo) 
  edgesFromTo -> names(edgeLabels)
  if(showEdgeData){
    eAttrs$label <- edgeLabels
  }
  edgeColors <- edgeLabels
  for(i in seq_along(edgeColors)){
    if(edgeColors[i] == "repression"){
      edgeColors[i]<- "red"
    } else 
      edgeColors[i]<- "black"
  }
  eAttrs$color <- edgeColors
  att = list(graph = list(rankdir = "LR", rank = ""), node = list(fixedsize = FALSE))
  plot(pathway.i, main = name, attrs = att, nodeAttrs = nA, edgeAttrs = eAttrs)
  legend("topright", legend = c("genes", "microRNAs"), 
         #col = twocolors,
         pch = c(19,19),
         cex = 0.9)
  legend("right", legend = c("activation", "repression"), 
         col = c("black", "red"),
         pch = c('-','-'),
         lwd = 3,
         cex = 0.9)
}

