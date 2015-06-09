# The following script constructs the augmented_pathways dataset
# included in the mirIntegrator package. 
# This dataset is included only for demostration purposes.

require("mirIntegrator")
data(mirTarBase)
data(kegg_pathways)
augmented_pathways <- integrate_mir(kegg_pathways, mirTarBase)  

save(augmented_pathways, file = "augmented_pathways.rda")
library(tools)
resaveRdaFiles("augmented_pathways.rda")

