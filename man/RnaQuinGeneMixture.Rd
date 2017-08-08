\docType{data}
\name{RnaQuinGeneMixture}
\alias{RnaQuinGeneMixture}
\title{RnaQuin mixture (gene level)}
\format{Data frame:
\itemize{
  \item Name: Sequin name
  \item Length: Gene length
  \item MixA: Input concentration for mixture A
  \item MixB: Input concentration for mixture B
}}
\usage{
data(RnaQuinGeneMixture)
}
\description{
Individual sequins are combined across a range of precise concentrations to
formulate mixtures. By modulating the concentration at which each sequin is
present in the mixture, we can emulate quantitative features of genome biology.

This is the mixture A and B in \code{RnaQuin}. File name is \code{A.R.6.csv} on
\code{http://www.sequins.xyz}.
}
\value{
Data frame with columns defined in \code{Format}.
}