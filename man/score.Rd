\name{score}
\alias{score}
\alias{score,AnaquinData-method}
\title{Return scores}
\description{
Return scores from data set created by \code{AnaquinData}.
}
\usage{
score(object, ...)
\S4method{score}{AnaquinData}(object)
}
\arguments{
  \item{object}{S4 object created by \code{AnaquinData}}
  \item{...}{Not used}
}
\details{
This function takes an \code{AnaquinData} object and return the ROC ranking
scores.
}
\value{
List of scores for each sequin. Return \code{NULL} if unavailable.
}
\examples{
data(UserGuideData_5.6.3)

# Sequin names
seqs <- row.names(UserGuideData_5.6.3)

# Expected log-fold
ratio <- UserGuideData_5.6.3$ExpLFC

# How the ROC curves are ranked
score <- 1-UserGuideData_5.6.3$Pval

# Classified labels (TP/FP)
label <- UserGuideData_5.6.3$Label

anaquin <- AnaquinData(analysis='PlotROC',
                           seqs=seqs,
                          ratio=ratio,
                          score=score,
                          label=label)
score(anaquin)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}