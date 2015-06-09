##########
#' Get the smallest pathway
#' 
#' @description
#' Find the pathway with the fewer number of nodes among a list of pathways.
#' This simple function is an example of how to navigate the genes on 
#' a list of pathways. 
#' 
#' @param pathways A list of graph::graphNEL objects.
#' 
#' @return
#' The index of the pathway with fewer number of nodes. 
#' 
#' @author
#' Diana Diaz <dmd at wayne dot edu>
#' 
#' @examples
#' data(augmented_pathways)
#' smallest_pathway(augmented_pathways)
#' smallest_pathway
#' 
#' @import graph
#' @export
smallest_pathway <- function(pathways){
  min <- Inf
  j <- 0
  for(i in seq_along(pathways)){ 
    n_nodes <- length(graph::nodes(pathways[[i]]))
    if(min > n_nodes){
      min <- n_nodes
      j <- i
    }
  }
  j
}


##########
#' Plotting of augmented pathway
#' 
#' @description
#' Functions for plotting a particular augmented pathway. In the plot,
#' miRNAs that were
#' added to the original pathway are differentiated from proteins that were 
#' originally in the
#' pathway. Blue boxes represent the proteins that were part of the
#' original
#' pathway, and black boxes represent the miRNAs that were added during 
#' augmentation.
#'  
#' @param original_pathway A
#' graph::graphNEL object where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent the genes involved in the 
#' pathway and edges represent the biological interactions (activation or
#' repression) among those genes. 
#' @param augmented_pathway A
#' graph::graphNEL object where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent genes 
#' and miRNAs involved in the 
#' pathway and edges represent the biological interactions (activation or
#' repression) among them.
#' @param pathway_name The name of the pathway.
#' 
#' @return
#' A plot of one augmented pathway with the new nodes
#' highlighted in black.
#' 
#' @author
#' Diana Diaz <dmd at wayne dot edu>
#' 
#' @examples
#' data(augmented_pathways)
#' data(kegg_pathways)
#' data(names_pathways) 
#' 
#' plot_augmented_pathway(kegg_pathways[[18]], augmented_pathways[[18]],
#'                       pathway_name = names_pathways[[18]])
#' 
#' @import graph
#' @export
plot_augmented_pathway <- function(original_pathway, augmented_pathway, pathway_name = " "){
  oldnodes = graph::nodes(original_pathway)
  newnodes = graph::nodes(augmented_pathway) 
  augmentednodes = setdiff(newnodes, oldnodes)
  plotPathway2Colors(pathway.i = augmented_pathway, subclass = augmentednodes, 
                     name = paste(pathway_name, "augmented pathway."))
}


##########
#' Export augmented pathways to pdf
#' 
#' @description
#' This function creates a pdf file with plottings of a list of augmented 
#' pathways.
#' 
#' @param original_pathways A list of 
#' graph::graphNEL objects where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent the genes involved in the 
#' pathway and edges represent the biological interactions (activation or
#' repression) among those genes (activation or repression). 
#' @param augmented_pathways A list of 
#' graph::graphNEL objects where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent genes 
#' and miRNAs involved in the 
#' pathway and edges represent the biological interactions (activation or 
#' repression) among them.
#' @param pathway_names A list of names of the pathways named by '<pathway_ID>'.
#' @param file The name of the file where the plots will be saved.
#' 
#' @return 
#' A pdf file with the plottings of the augmented pathways.
#' 
#' @author
#' Diana Diaz <dmd at wayne dot edu>
#' 
#' @examples
#' data(augmented_pathways)
#' data(kegg_pathways)
#' data(names_pathways)
#' #The following instruction writes a pfd with three pathways
#' pathways2pdf(kegg_pathways[18:20],augmented_pathways[18:20],
#'              names_pathways[18:20], "three_pathways.pdf")
#' #The following instruction writes a pfd with all the pathways:
#' #NOTE: It may take time.
#' # pathways2pdf(kegg_pathways,augmented_pathways, 
#' #              names_pathways, "all_pathways.pdf")
#' 
#' @import ROntoTools 
#' @import graph
#' @export
pathways2pdf <- function(original_pathways,
                         augmented_pathways,pathway_names,file){
  pdf( file=file, onefile=TRUE)  
  opar <- par()
  par(mfrow=c(1,1))
  for(i in seq_along(original_pathways)){ 
    pathway.i = original_pathways[[i]]
    newpathw.i = augmented_pathways[[i]]
    name.i = pathway_names[i]
    plot_augmented_pathway(pathway.i, newpathw.i, name.i)
  }
  par(opar)
  graphics.off()
}




#' @import ggplot2
#private
plotLines <- function (line1, line2, line3,lab1, lab2, lab3,xlab, ylab ){
  #1:
  df <- data.frame(x = rep(seq_along(line1), 3), y = c(line1, line2, line3), 
                   colors = rep(c(lab1, lab2, lab3), each=length(line1)) )
  plines <- ggplot2::ggplot(data = df, ggplot2::aes(x=x, y=y, col=colors)) + ggplot2::geom_line() + ggplot2::xlab(xlab)  + ggplot2::ylab(ylab)
  plines + ggplot2::guides(fill=ggplot2::guide_legend(title=NULL))
  plines + ggplot2::theme(legend.title=ggplot2::element_blank())
  plines
}

##########
#' Plotting the change in pathways order
#' 
#' @description
#' Function for plotting a lines plot of the difference in pathways' order.
#' The resultant plot shows the comparison between the order of the original 
#' pathways 
#' and the order of the augmented pathways. It also contains a line with
#' the order difference (order of the augmented pathways minus 
#' order of the original pathways). The order of a biological pathway is 
#' the number of genes that are involved in it.
#'  
#' @param original_pathways A list of 
#' graph::graphNEL objects where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent the genes involved in the 
#' pathway and edges represent the biological interactions (activation or
#' repression) among those genes (activation or repression). 
#' @param augmented_pathways A list of 
#' graph::graphNEL objects where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent genes 
#' and miRNAs involved in the 
#' pathway and edges represent the biological interactions (activation or 
#' repression) among them.
#' @param pathway_names A list of names of the pathways named by '<pathway_ID>'.
#' 
#' @return
#' A lines plot of the comparison of pathways order.
#' 
#' @author
#' Diana Diaz <dmd at wayne dot edu>
#' 
#' @examples
#' data(augmented_pathways)
#' data(kegg_pathways)
#' data(names_pathways)
#' plot_change(kegg_pathways,augmented_pathways, names_pathways)
#' 
#' @import ROntoTools 
#' @import graph
#' @export
plot_change <- function(original_pathways,augmented_pathways,pathway_names){
  microRNAadded <- data.frame(name = pathway_names, genesOriginal = NA, genesAugmented = NA, microRNAadded  = NA)
  for(i in seq_along(original_pathways)){ 
    pathway.i = original_pathways[[i]]
    newpathw.i = augmented_pathways[[i]]
    genesOriginal = length(graph::nodes(pathway.i))
    genesAugmentes = length(graph::nodes(newpathw.i))
    microRNAaddedn = genesAugmentes - genesOriginal
    microRNAadded[i,2:4] = c(genesOriginal,genesAugmentes, microRNAaddedn)
  }
  microRNAadded <- microRNAadded[with(microRNAadded, order(genesOriginal)),]
  
  plines <- plotLines(line1 = microRNAadded$genesOriginal,   line2 = microRNAadded$genesAugmented,
                      line3 = microRNAadded$microRNAadded,lab1 = "Original", lab2 = "Augmented",
                      lab3 = "Difference",   xlab = "Pathways",  ylab = "Number of genes")
  plines
}