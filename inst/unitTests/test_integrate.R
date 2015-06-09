test_integrate <- function(){
  set.seed(42L)
  require("graph")
  data(kegg_pathways)
  data(mirTarBase)
  kegg_pathways <- kegg_pathways[18:19]
  augmented_pathways2 <- integrate_mir(kegg_pathways, mirTarBase)  
  checkEquals(length(nodes(augmented_pathways2[[1]])), 20)
  checkEquals(length(nodes(augmented_pathways2$"path:hsa04130")), 77)
  data(augmented_pathways)
  checkEquals(all.equal.list(augmented_pathways[18:19],
                             augmented_pathways2), TRUE)
}


