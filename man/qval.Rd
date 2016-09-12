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
data(UserGuideData_5.6.3)

# Sequin names
seqs <- row.names(UserGuideData_5.6.3)

# Expected log-fold
ratio <- UserGuideData_5.6.3$ExpLFC

# Measured average abundance
measured <- UserGuideData_5.6.3$Mean

# P-value
pv <- UserGuideData_5.6.3$Pval

# Q-value
qv <- UserGuideData_5.6.3$Qval

anaquin <- AnaquinData(analysis='PlotLODR',
                           seqs=seqs,
                      measured=measured,
                         ratio=ratio,
                          pval=pv,
                          qval=qv)
qval(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}