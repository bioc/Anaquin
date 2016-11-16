\name{RnaQuin.genes}
\alias{RnaQuin.genes}
\alias{RnaQuin.genes}
\title{Filter RnaQuin genes}
\description{
Filter RnaQuin genes from a list of genes names.
}
\usage{
RnaQuin.genes(genes)
}
\arguments{
  \item{genes}{List of character strings for the gene names}
}
\details{
This function provideds a convenient mechanism to filter RnaQuin genes from
a list of gene names. A most common usage is for filtering a gene expression
table.
}
\value{
List of RnaQuin genes from the input. If no sequin gene is found, the function
returns an empty list.
}
\examples{
# This could also come from a count matrix table
x <- c('R1_11', 'R1_12')

RnaQuin.genes(x)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}