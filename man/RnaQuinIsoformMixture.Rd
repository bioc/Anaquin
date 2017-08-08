\docType{data}
\name{RnaQuinIsoformMixture}
\alias{RnaQuinIsoformMixture}
\title{RnaQuin mixture (isoform level)}
\format{Data frame:
\itemize{
  \item Name: Sequin name
  \item Length: Sequin length
  \item MixA: Input concentration for mixture A
  \item MixB: Input concentration for mixture B
}}
\usage{
data(RnaQuinIsoformMixture)
}
\description{
Individual sequins are combined across a range of precise concentrations to
formulate mixtures. By modulating the concentration at which each sequin is
present in the mixture, we can emulate quantitative features of genome biology.

This is the mixture A and B in \code{RnaQuin}. File name is \code{A.R.5.csv} on
\code{http://www.sequins.xyz}.
}
\value{
Data frame with columns defined in \code{Format}.
}