\name{label}
\alias{label}
\alias{label,AnaquinData-method}
\title{Return classified labels}
\description{
Return classified labels from data set created by \code{AnaquinData}.
}
\usage{
label(object, ...)
\S4method{label}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} object and return the classified
labels.
}
\value{
List of classified labels for each sequin. Return \code{NULL} if unavailable.
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
label(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}