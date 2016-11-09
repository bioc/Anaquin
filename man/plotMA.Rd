\name{plotMA}
\alias{plotMA}
\title{Create MA plot}
\description{
Create an MA plot between the M (log ratios) and A (mean average) for
endogenous features and sequins.
}
\usage{
plotMA(data, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{AnaquinData}. It needs to
              define information in \bold{Details}.}
  \item{...}{Reserved for internal testing}
}
\details{

\code{plotMA} requires the following data inputs from \code{AnaquinData}.

\tabular{rlll}{
    \tab \code{seqs}    \tab List of sequin identifiers (eg. R2_11_2).\cr
    \tab \code{mean}    \tab Log mean average\cr
    \tab \code{lfc}     \tab Log-fold ratio\cr
}
Create an MA plot between the M (log ratios) and A (mean average) for
endogenous features and sequins.

This function can used for differential isoform/gene analysis. The
distribution of the sequins are compared with the samples.
}

\value{
The functions returns a `\code{ggplot2} object for the plot.
}

\examples{
library(Anaquin)

data('UserGuideData_5.6.3')
data <- AnaquinData(analysis='PlotMA',
                        seqs=row.names(UserGuideData_5.6.3),
                        mean=log2(UserGuideData_5.6.3$Mean),
                        lfc=UserGuideData_5.6.3$ObsLFC)

plotMA(data)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}