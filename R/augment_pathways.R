##########
# @description
# Find the MIRs that targets each gene on each pathway 
# 
# @param ent.gen.i The entrez id of the gene i
#
# @return 
# Get the miRNA that regulate these genes
# 
# @author
# Diana Diaz <dmd at wayne dot edu>
# 
# @seealso \code{\link{peRes}} \code{\link{summary}}
#' @import ROntoTools 
#' @import graph
# private
mirTarget <- function(ent.gen.i, targets_db)
{
  mirTarget.i <- targets_db[targets_db$Target.ID %in% ent.gen.i,"miRNA"]
  mirTarget.i <- as.character(unique(mirTarget.i))
  return(mirTarget.i)
}



# Find the targets of the mir.i and merge them to the pathway.i
#' @import graph
# private
addMirEffect <- function(mir.i, ent.gen.i, 
                         newpathway.i, targets_db , 
                         typeOfRelationMirnas = "repression" )
{
  mir.i <- unique(mir.i)
  alltargets <- targets_db[targets_db$miRNA %in% mir.i,"Target.ID"]
  alltargets <- unique(as.character(alltargets)) 
  targets.i <- intersect(ent.gen.i, alltargets)
  newpathway.i <- graph::addNode(node=mir.i, object=newpathway.i)
  newpathway.i <- graph::addEdge(from=mir.i, to=as.character(targets.i),
                                 graph=newpathway.i)
  edgeData(newpathway.i, from=mir.i, to=as.character(targets.i),
           attr="subtype") <- typeOfRelationMirnas
  newpathway.i
}
# Private
generatePathway.i <- function(pathway.i, targets)
{
  ent.gen.i <- graph::nodes(pathway.i)
  mir.i <- mirTarget(ent.gen.i, targets)
  for(i in seq_along(mir.i)){  
    pathway.i <- addMirEffect(mir.i= mir.i[i], 
                                 ent.gen.i, pathway.i, targets)
  }
  pathway.i
}


##########
#' Produce augmented pathways
#' 
#' @description
#' This function takes each pathway of the input list of signaling pathways and
#' adds the miRNAs that are related to it.
#'  
#' @param original_pathways A list of 
#' graph::graphNEL objects where each of the nodes is named with '<gene_ID>'.
#' Gene IDs used to identify the nodes must be the same gene IDs used to
#' identify the genes on the miRNA-target 
#' interactions data.frame, \code{targets_db}.
#' i.e. If the genes are identified by Entrez ID on the \code{original_pathways}
#' graph::graphNEL list, then the \code{targets_db} 
#' data.frame must identify the genes by Entrez ID as well.
#' Nodes of each graph::graphNEL represent the genes involved in the 
#' pathway and edges represent the biological interactions (activation or
#' repression) among those genes (activation or repression). 
#' @param targets_db A data.frame with
#' columns: 'miRNA' which names the miRNAs
#' and 'Target.ID' which gives the gene ID of the target gene.
#' The Gene IDs used to identify the "Target.ID" column 
#' must be the same gene IDs used on the nodes of the \code{original_pathways}.
#' i.e. If the genes are identified by Entrez ID on the \code{original_pathways}
#' graph::graphNEL list, then the \code{targets_db} 
#' data.frame must identify the genes by Entrez ID as well.
#' 
#' @return
#' Gene signaling pathways augmented with miRNA interactions.
#' This is a list of 
#' graph::graphNEL objects where each of the nodes is named with '<gene_ID>'.
#' Nodes of each graph::graphNEL represent genes 
#' and miRNAs involved in the 
#' pathway and edges represent the biological interactions (activation or 
#' repression) among them.
#' 
#' @author
#' Diana Diaz <dmd at wayne dot edu>
#' 
#' @examples
#' data(kegg_pathways)
#' data(mirTarBase)
#' kegg_pathways <- kegg_pathways[1:5] #delete this for augmenting all pathways.
#' augmented_pathways <- integrate_mir(kegg_pathways, mirTarBase)
#' 
#' @import ROntoTools
#' @export
integrate_mir <- function(original_pathways, targets_db)
{
  isGraphNEL <- function(x) is(x, "graphNEL")
  if( ("miRNA" %in% names(targets_db)) && ("Target.ID" %in% names(targets_db))
     && is.data.frame(targets_db) 
     && all(vapply(original_pathways, isGraphNEL, logical(1), USE.NAMES=FALSE))  )
  {
    targets_db <- targets_db[,c("miRNA","Target.ID")]
    genPath <- function(X){
      generatePathway.i(X,targets_db)
    }  
    newpathways <- lapply(original_pathways, genPath)
    return(newpathways)
  } else {
    stop("'targets_db' data.frame must contain the columns:
          'Target.ID' and 'miRNA'.")
  }
}
