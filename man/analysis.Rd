\name{analysis}
\alias{analysis}
\alias{analysis,AnaquinData-method}
\title{Return type of analysis}
\usage{
analysis(object, ...)

\S4method{analysis}{AnaquinData}(object)
}
\arguments{
\item{...}{Not used}
\item{object}{S4 object created by \code{AnaquinData}}
}
\value{
Type of analysis, eg: \code{PlotLinear}
}
\description{
Return type of analysis from data set created by \code{AnaquinData}
}
\details{
This function takes an \code{AnaquinData} object and return the type of
analysis.
}
\examples{
# Sequin names
names <- row.names(seqCufflinks)

#Input concentration
input <- log2(seqCufflinks$InputConcent)

# Measured FPKM
measured <- log2(seqCufflinks[,2:4])

anaquin <- AnaquinData(analysis='PlotLinear',
                           seqs=names,
                          input=input,
                       measured=measured)
analysis(anaquin)

}

