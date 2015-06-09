# The following script downloads the names_pathways list
# included in the mirIntegrator package. 
# This list is included only for demostration purposes.

library("ROntoTools")
names_pathways <- keggPathwayNames("hsa")
save(names_pathways, file = "names_pathways.rda")
library(tools)
resaveRdaFiles("names_pathways.rda")

