\name{RnaQuin.isoforms}
\alias{RnaQuin.isoforms}
\alias{RnaQuin.isoforms}
\title{Filter RnaQuin isoforms}
\description{
Filter RnaQuin isoforms from a list of isoforms names.
}
\usage{
RnaQuin.isoforms(isos)
}
\arguments{
  \item{isos}{List of character strings for the isoform names}
}
\details{
This function provideds a convenient mechanism to filter RnaQuin isoforms from
a list of isoform names. A most common usage is for filtering a isoform expression
table.
}
\value{
List of RnaQuin isoforms from the input. If no sequin isoforms is found, the function
returns an empty list.
}
\examples{
# This could also come from a count matrix table
x <- c('R1_11_1', 'This_Is_Not_Sequin', 'R1_12_2')

RnaQuin.isoforms(x)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}