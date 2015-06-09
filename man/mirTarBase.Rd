\name{mirTarBase}
\alias{mirTarBase}
\docType{data}
\title{
MicroRNA-target interactions in human.
}
\description{
Dataset of miRNA-target interactions in human obtained from mirTarBase
release 4.5: Nov. 1, 2013.
}
\usage{data(mirTarBase)}
\format{
  A data.frame with 39083 interactions and 9 variables.
  The columns needed for this package are:
  \itemize{
    \item \code{miRNA}: which contains the miRNA ID,
   \item \code{Target.ID}: contains the entrez ID of the gene targeted by the 
  \code{miRNA}
  }
}
\value{
  A data.frame with human miRNA-targets interactions
}
\details{
 This dataset is licensed by its authors (Hsu et al.), see
 \url{http://mirtarbase.mbc.nctu.edu.tw/cache/download/LICENSE}.
}
\source{
Downloaded from \url{http://mirtarbase.mbc.nctu.edu.tw/} on 4/1/2015. A script 
which downloads the file and constructs the
mirTarBase object may be found in 'inst/scripts/get_mirTarBase.R', see the 
example.
}
\references{
S.-D. Hsu, Y.-T. Tseng, S. Shrestha, Y.-L. Lin, A. Khaleel, C.-H. Chou, C.-F. 
Chu, H.-Y. Huang, C.-M. Lin, S.-Y. Ho, T.-Y. Jian, F.-M. Lin, T.-H. Chang, S.-L. 
Weng, K.-W. Liao, I.-E. Liao, C.-C. Liu, and H.-D. Huang, \cite{miRTarBase 
update 2014: an information resource for experimentally validated miRNA-target 
interactions}, Nucleic Acids Research, vol. 42,
pp. D78 - D85, Jan. 2014.
}
\examples{
data(mirTarBase)
head(mirTarBase)

script <- system.file("scripts", "get_mirTarBase.R", 
                      package = "mirIntegrator")
script
readLines(script)
}
\keyword{datasets}
