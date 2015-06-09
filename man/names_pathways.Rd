\name{names_pathways}
\alias{names_pathways}
\docType{data}
\title{
List of KEGG signaling pathways' names.
}
\description{
Names of the KEGG signaling pathways in human obtained with the 
ROntoTools package. The original KEGG pathways were published by 
Kanehisa Laboratories, release 73.0+/01-03, Jan 2015.
}
\usage{data("names_pathways")}
\value{
A list of KEGG signaling pathways' names. 
}
\source{
Obtained using the ROntoTools package Version 1.2.0 with KEGG database
release 73.0+/01-03, Jan 2015. 
A script that constructs the names_pathways object may be found in 
'inst/scripts/get_names_pathways.R', see the example.
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
 data(names_pathways)
 head(names_pathways)
 
 script <- system.file("scripts", "get_names_pathways.R", 
                      package = "mirIntegrator")
 script
 readLines(script)
}
\keyword{datasets}
