# The following script downloads the kegg_pathways dataset
# included in the mirIntegrator package. 
# This dataset is included only for demostration purposes.

library("ROntoTools")
kegg_pathways_O <- keggPathwayGraphs("hsa")

set.gene.id <- function(pathway.i)
{
  genes.i <- nodes(pathway.i)
  ent.gen.i <- gsub(".*\\:","",genes.i)
  nodes(pathway.i) <- ent.gen.i
  pathway.i
}

kegg_pathways <- lapply(kegg_pathways_O, set.gene.id)

save(kegg_pathways, file = "kegg_pathways.rda")
library(tools)
resaveRdaFiles("kegg_pathways.rda")
