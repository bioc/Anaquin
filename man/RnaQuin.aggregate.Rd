\name{RnaQuin.aggregate}
\alias{RnaQuin.aggregate}
\alias{RnaQuin.aggregate}
\title{Aggregate isoform concentration}
\description{
Aggregate isoform concentration for genes
}
\usage{
RnaQuin.aggregate(isos, concent)
}
\arguments{
  \item{isos}{List of RnaQuin isoform names}
  \item{concent}{List of input concentration for each RnaQuin isoform}
}
\details{
This function provides a convenient mechanism to add up the RnaQuin input
concentration for isoforms to the gene level. For example, if we have the
following isoforms: \code{R1_101_1}: 10 attomol/ul and \code{R1_101_2}:
15 attomol/ul.

\code{R1_101_1} and \code{R1_101_2} forms the splice variants of a gene.
Their input concentration would add up to 25 attomol/ul for
\code{R1_101}.
}
\value{
Data frame with two columns. The first column gives the RnaQuin gene names.
The second column gives the aggregated gene level input concentation
attomol/ul. 
}
\examples{
data('RnaQuinMixtureA')

# Aggregate input concentation for isoforms into genes
RnaQuin.aggregate(RnaQuinMixtureA$ID, RnaQuinMixtureA$MixA)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}