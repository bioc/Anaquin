\name{ratio}
\alias{ratio}
\alias{ratio,AnaquinData-method}
\title{Return ratio}
\description{
Return ratio from data set created by \code{AnaquinData}.
}
\usage{
ratio(object, ...)
\S4method{ratio}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} data set object and return its
sequin ratio.
}
\value{
List of sequin ratio for each sequin. Return \code{NULL} if unavailable.
}
\examples{
data('seqDESeq2')

# Sequin names
seqs <- row.names(seqDESeq2)

# Expected log-fold
ratio <- seqDESeq2$ExpLFC

# How the ROC curves are ranked
score <- 1-seqDESeq2$Pval

# Classified labels (TP/FP)
label <- seqDESeq2$Label

anaquin <- AnaquinData(analysis='PlotROC',
                           seqs=seqs,
                          ratio=ratio,
                          score=score,
                          label=label)
ratio(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}