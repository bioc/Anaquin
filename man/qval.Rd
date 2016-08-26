\name{qval}
\alias{qval}
\alias{qval,AnaquinData-method}
\title{Return q-value probability}
\description{
Return q-value probability from data set created by \code{AnaquinData}.
}
\usage{
qval(object, ...)
\S4method{qval}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}  
}
\details{
This function takes an \code{AnaquinData} object and return it's q-value
probability.
}
\value{
List of q-value proability for each sequin. Return \code{NULL} if unavailable.
}
\examples{
data('seqDESeq2')

# Sequin names
seqs <- row.names(seqDESeq2)

# Expected log-fold
ratio <- seqDESeq2$ExpLFC

# Measured average abundance
measured <- seqDESeq2$Mean

# P-value
pv <- seqDESeq2$Pval

# Q-value
qv <- seqDESeq2$Qval

anaquin <- AnaquinData(analysis='PlotLODR',
                           seqs=seqs,
                      measured=measured,
                         ratio=ratio,
                          pval=pv,
                          qval=qv)
qval(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}