\name{measured}
\alias{measured}
\alias{measured,AnaquinData-method}
\title{Return measured abundance}
\description{
Return measured abundance from data set created by \code{AnaquinData}.
}
\usage{
measured(object, ...)
\S4method{measured}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} object and return the measured
abundance.
}
\value{
List of measured abundance for each sequin. Return \code{NULL} if unavailable.
}
\examples{
data(UserGuideData_5.4.6.3)

# Sequin names
seqs <- row.names(UserGuideData_5.4.6.3)

# Input concentration
input <- log2(UserGuideData_5.4.6.3$InputConcent)

# Measured FPKM
FPKM <- log2(UserGuideData_5.4.6.3[,2:4])

anaquin <- AnaquinData(analysis='PlotLinear',
                           seqs=seqs,
                          input=input,
                       measured=FPKM)
measured(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}