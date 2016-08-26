\name{seqs}
\alias{seqs}
\alias{seqs,AnaquinData-method}
\title{Return sequin names}
\description{
Return sequin names from \code{AnaquinData} data set.
}
\usage{
seqs(object, ...)
\S4method{seqs}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} object and return the sequin names
in the data set.
}
\value{
List of sequin names.
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
seqs(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}