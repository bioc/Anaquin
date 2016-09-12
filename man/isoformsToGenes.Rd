\name{isoformsToGenes}
\alias{isoformsToGenes}
\title{Convert sequin isoforms to genes}
\description{
Convert RnaQuin sequin isoforms to genes. Due to alternative splicing,
there will be smaller number of genes than isoforms.
}
\usage{
isoformsToGenes(names)
}
\arguments{
  \item{names}{Isoform names}
}
\details{
For example:

  R1_11_1 to R1_11
  R1_11_2 to R1_11
  R1_12_2 to R1_12
}
\value{
A list of RnaQuin sequin genes.
}
\examples{
data(UserGuideData_5.4.5.1)
isoformsToGenes(row.names(UserGuideData_5.4.5.1))
}

\author{Ted Wong \email{t.wong@garvan.org.au}}