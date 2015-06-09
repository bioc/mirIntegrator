test_datasets <- function(){
  data("mirTarBase", "augmented_pathways", "GSE43592_miRNA", "GSE43592_mRNA",
       "kegg_pathways")
  checkEquals(dim(mirTarBase)[1], 39091)
  checkEquals(dim(mirTarBase)[2], 9)
  checkEquals(length(augmented_pathways), 149)
  checkEquals(dim(GSE43592_miRNA)[1], 881)
  checkEquals(dim(GSE43592_mRNA)[1], 19611)
  checkEquals(length(kegg_pathways), 149)
}