\name{augmented_pathways}
\alias{augmented_pathways}
\docType{data}
\title{
Signaling pathways augmented with miRNA.
}
\description{
Human signaling KEGG pathways augmented with validated miRNA-target 
interactions from mirTarBase using the mirIntegrator package. These 
interactions represent the biological miRNA 
repression over its target genes and are included in the model as negative
links.
}
\usage{data("augmented_pathways")}
\value{
 A list of graphNEL objects where each graph 
 is a pathway that were augmented with miRNA-target interactions.  The name of 
 each pathway 
 is its KEGG pathway identifier.
}
\source{
Generated using the mirIntegrator package. 
A script that constructs the augmented_pathways object may be found in 
'inst/scripts/get_augmented_pathways.R', see the example.
}
\seealso{
  \code{\link{mirTarBase}} and
  \code{\link{kegg_pathways}}
}
\references{
M. Kanehisa and S. Goto, \cite{KEGG: Kyoto Encyclopedia of Genes and Genomes},
Nucleic Acids Research, vol. 28, pp. 27-30, January 2000.

S.-D. Hsu, Y.-T. Tseng, S. Shrestha, Y.-L. Lin, A. Khaleel, C.-H. Chou, C.-F. 
Chu, H.-Y. Huang, C.-M. Lin, S.-Y. Ho, T.-Y. Jian, F.-M. Lin, T.-H. Chang, S.-L. 
Weng, K.-W. Liao, I.-E. Liao, C.-C. Liu, and H.-D. Huang, \cite{miRTarBase 
update 2014: an information resource for experimentally validated miRNA-target 
interactions}, Nucleic Acids Research, vol. 42,
pp. D78 - D85, Jan. 2014.
}
\examples{
data(augmented_pathways)
head(augmented_pathways)

script <- system.file("scripts", "get_augmented_pathways.R", 
                      package = "mirIntegrator")
script
readLines(script)
}
\keyword{datasets, pathway, microRNA, miRNA}
