\name{RnaQuin.iso2gen}
\alias{RnaQuin.iso2gen}
\alias{RnaQuin.iso2gen}
\title{Convert RnaQuin isoforms to genes}
\description{
Convert RnaQuin isoform names to gene names.
}
\usage{
RnaQuin.iso2gen(isos)
}
\arguments{
  \item{isos}{List of character strings for the isoform names}
}
\details{
This function provideds a convenient mechanism to convert RnaQuin isoform names
to gene names.
}
\value{
List of RnaQuin gene names for each input isoform.
}
\examples{
# This could also come from a count matrix table
x <- c('R1_1_1', 'R1_1_2', 'R1_1_3', 'R1_2_1')

# Convert RnaQuin isoforms to genes
RnaQuin.iso2gen(x)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}