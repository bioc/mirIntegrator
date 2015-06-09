\name{kegg_pathways}
\alias{kegg_pathways}
\docType{data}
\title{
List of KEGG signaling pathways of human.
}
\description{
This dataset contains 149 KEGG signaling pathways of human. The original pathways were parsed to a list of graphNEL objects using the ROntoTools package. The original KEGG pathways were published by Kanehisa Laboratories, release 73.0+/01-03, Jan 2015.
}
\usage{data("kegg_pathways")}
\value{
 A list of graphNEL objects where each graph
 represents one KEGG signaling pathway. The name of each pathway 
 is its KEGG pathway identifier.
}
\source{
Obtained using the ROntoTools package Version 1.2.0 with KEGG database
release 73.0+/01-03, Jan 2015.
A script that constructs the kegg_pathways object may be found in 
'inst/scripts/get_kegg_pathways.R', see the example.
}
\references{
M. Kanehisa and S. Goto, \cite{KEGG: Kyoto Encyclopedia of Genes and Genomes},
Nucleic Acids Research, vol. 28, pp. 27-30, January 2000.

C. Voichita, M. Donato, and S. Draghici, \cite{Incorporating gene significance 
in the impact analysis of signaling pathways}, in 2012 11th International 
Conference on Machine Learning and Applications (ICMLA), vol. 1, pp. 126-131, 
Dec. 2012.
}
\examples{
data(kegg_pathways)
head(kegg_pathways)

script <- system.file("scripts", "get_kegg_pathways.R", 
                      package = "mirIntegrator")
script
readLines(script)
}
\keyword{datasets}
