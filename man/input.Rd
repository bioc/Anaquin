\name{input}
\alias{input}
\alias{input,AnaquinData-method}
\title{Return input concentration}
\description{
Return input concentarion from data set created by \code{AnaquinData}.
}
\usage{
input(object, ...)
\S4method{input}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} object and return the input
concentration in the data set.
}
\value{
List of input concentration for each sequin. Return \code{NULL} if unavailable.
}
\examples{
data('seqCufflinks')

# Sequin names
names <- row.names(seqCufflinks)

# Input concentration
input <- log2(seqCufflinks$InputConcent)

# Measured FPKM
measured <- log2(seqCufflinks[,2:4])

anaquin <- AnaquinData(analysis='PlotLinear',
                           seqs=names,
                          input=input,
                       measured=measured)
input(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}