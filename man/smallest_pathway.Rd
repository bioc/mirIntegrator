% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/view_pathways.R
\name{smallest_pathway}
\alias{smallest_pathway}
\title{Get the smallest pathway}
\usage{
smallest_pathway(pathways)
}
\arguments{
\item{pathways}{A list of graph::graphNEL objects.}
}
\value{
The index of the pathway with fewer number of nodes.
}
\description{
Find the pathway with the fewer number of nodes among a list of pathways.
This simple function is an example of how to navigate the genes on 
a list of pathways.
}
\examples{
data(augmented_pathways)
smallest_pathway(augmented_pathways)
smallest_pathway

}
\author{
Diana Diaz <dmd at wayne dot edu>
}

