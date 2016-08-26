\name{plotLogistic}
\alias{plotLogistic}
\title{Plot logistic model for sequins}
\description{
Plot a GLM logistic model between sequins against measured variable
(eg: sensitivity).
}
\usage{
plotLogistic(data, title, xlab, ylab, showLOA, threshold, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{AnaquinData}. It needs to define
              information in \bold{Details}.}
  \item{title}{Label of the plot.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{showLOA}{Show limit of assembly?}
  \item{threshold}{Threshold required for limit-of-assembly (LOA)}
  \item{...}{Reserved for internal testing}
}
\value{
The functions does not return anything but it prints a scatter plot with a fitted logistic function.
}
\details{

\code{plotLogistic} requires the following inputs from
\code{AnaquinData}.

\tabular{rlll}{
    \tab \code{seqs}     \tab Sequin names\cr
    \tab \code{input}    \tab Input concentration in attomol/ul\cr
    \tab \code{measured} \tab Measurement (eg: sensitivity)\cr
}

The analysis specifies expected abundance (eg: input concentration) on the
x-axis, and measured abundance (eg: sensitivity) on the y-axis. Curve
fitting is performed by non-linear least square fitting. LOA is defined as
the least abundant sequin reaching sensitivity of \code{threshold}.
}

\examples{
data('seqCuffcompare')

title <- 'Assembly Plot'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'Sensitivity'

# Sequin names
seqs <- row.names(seqCuffcompare)

# Input concentration
input <- log2(seqCuffcompare$InputConcent)

# Measured sensitivity
measured <- seqCuffcompare$Sn

anaquin <- AnaquinData(analysis='PlotLogistic',
                           seqs=seqs,
                          input=input,
                       measured=measured)

plotLogistic(anaquin, title=title, xlab=xlab, ylab=ylab, showLOA=TRUE)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}
