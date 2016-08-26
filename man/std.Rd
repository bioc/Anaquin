\name{std}
\alias{std}
\alias{std,AnaquinData-method}
\title{Return standard deviation}
\description{
Return standard deviation from \code{AnaquinData} data set.
}
\usage{
std(object, ...)
\S4method{std}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} object and return the standard
deviation in the data set.
}
\value{
List of standard deviation for each sequin. Return \code{NULL} if unavailable.
}
\examples{
data('seqDESeq2')

# Sequin names
names <- row.names(seqDESeq2)

# Expected LFC
input <- log2(seqDESeq2$ExpLFC)

# Measured LFC
measured <- log2(seqDESeq2$ObsLFC)

# Standard deviation
sd <- seqDESeq2$SD

anaquin <- AnaquinData(analysis='PlotLinear',
                           seqs=names,
                          input=input,
                            std=sd,
                       measured=measured)
std(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}