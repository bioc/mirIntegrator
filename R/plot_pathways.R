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
                               twocolors = c("lightgreen", "black"), twoshapes = c("box", "box"),
                               twofontcolor = c("black", "white") , fontsize = 10, numberofchar = 15, showEdgeData = TRUE){
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
  if(showEdgeData){
    edgeLabels <- unlist(edgeData(pathway.i))
  }
  nodeType <- 1 + (graph::nodes(pathway.i) %in% subclass)
  nA = Rgraphviz::makeNodeAttrs(pathway.i,fixedSize=FALSE, 
                     height = "1", width = "1",
                     fillcolor = twocolors[nodeType], shape = twoshapes[nodeType], 
                     fontcolor = twofontcolor[nodeType], fontsize = fontsize )
  eAttrs <- list()
  edgesFromTo <- gsub("\\.", "~", edgesFromTo) 
  edgesFromTo -> names(edgeLabels)
  eAttrs$label <- edgeLabels
  att = list(graph = list(rankdir = "LR", rank = ""))
  plot(pathway.i, main = name, attrs = att, nodeAttrs = nA, edgeAttrs = eAttrs)
  legend("bottomright", legend = c("original", "new"), 
         col = twocolors,
         pch = c(19,19), title = " ",
         lwd = 2, 
         cex = 0.5)
}

