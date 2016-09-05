\name{plotLinear}
\alias{plotLinear}
\title{Plot linear model for sequins}
\description{
Plot a linear model between sequins against measured abundance (eg: FPKM).
}
\usage{
plotLinear(data, title, xlab, ylab, showSD, showLOQ, xBreaks, yBreaks, showAxis, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{AnaquinData}. It needs to define
              information in \bold{Details}.}
  \item{title}{Label of the plot.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{xBreaks}{Breaks for the x-axis}
  \item{yBreaks}{Breaks for the y-axis}
  \item{showSD}{Show standard deviation bars vertically? Default to
                \code{TRUE}.}
  \item{showLOQ}{Show limit-of-quantification? Default to \code{TRUE}.}
  \item{showAxis}{Show x-axis and y-axis? Default to \code{TRUE}.}
  \item{...}{Reserved for internal testing}  
}
\value{
The functions does not return anything but it prints a scatter plot.
}
\details{

\code{plotLinear} requires the following inputs from \code{AnaquinData}.

\tabular{rlll}{
    \tab \code{seqs}   \tab Sequin names (eg: R1_1_1)\cr
    \tab \code{input}   \tab Input concentration in attomol/ul\cr
    \tab \code{measured} \tab Measured abundance (eg: FPKM)\cr
}

The function plots expected abundance on the x-axis, and measured abundance on
the y-axis. Generally, the expected abundance is simply the input
concentration, although other measures (eg: expected allele frequency) is also
possible. The plot builds a linear regression between the two variables, and
reports the standard statistics; R2, correlation and regression parameters on
the graph.

}

\examples{
data('seqCufflinks')

title <- 'Gene Expression'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'FPKM (log2)'

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

plotLinear(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}