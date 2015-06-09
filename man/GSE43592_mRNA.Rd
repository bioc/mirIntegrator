\name{GSE43592_mRNA}
\alias{GSE43592_mRNA}
\docType{data}
\title{
Top table of preprocessed mRNA of GSE43592 dataset. 
}
\description{
A data.frame with the Log fold change and p-value of preprocessed mRNA 
expression of GSE43592 dataset.
}
\usage{data(GSE43592_mRNA)}
\value{
  A data frame with 19611 mRNAs with the following 8 variables:
  entre, ID, logFC, AveExpr, t, P.Value, adj.P.Val, B.
}
\source{
Raw data obtained from
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43592} and 
preprocessed with the limma package version 3.24.0. 
}
\references{
M. Jernas, C. Malmestrom, M. Axelsson, I. Nookaew, H. Wadenvik, J. Lycke, and 
B. Olsson,\cite{MicroRNA regulate immune pathways in t-cells in multiple 
sclerosis (MS),} BMC immunology, vol. 14, p. 32, 2013.

Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). 
\cite{limma powers differential expression analyses for RNA-sequencing and microarray studies.} Nucleic Acids Research, 43(7), pp. e47.
}
\examples{
data(GSE43592_mRNA)
head(GSE43592_mRNA)
}
\keyword{datasets}
