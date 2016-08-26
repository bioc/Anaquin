\name{plotROC}
\alias{plotROC}
\title{Create ROC plot }
\description{
Create a receiver operating characteristic (ROC) plot. The curve models true
positive rate (TPR) against false positive rate (FPR) at various threshold
settings.
}
\usage{
plotROC(data, refRats, title, legTitle, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{AnaquinData}. It needs to
              define information in \bold{Details}.}
  \item{refRats}{Reference ratio groups}
  \item{title}{Label of the plot. Default to \code{NULL}.}
  \item{legTitle}{Title of the legend. Default to \code{Ratio}.}
  \item{...}{Reserved for internal testing}
}
\details{

\code{plotROC} requires the following inputs from \code{AnaquinData}.

\tabular{rlll}{
    \tab \code{seqs}  \tab Sequin names\cr
    \tab \code{label} \tab Classified labels ('TP' or 'FP')\cr
    \tab \code{score} \tab How the ROC points should be ranked\cr    
    \tab \code{ratio} \tab Expected ratio; eg: expected log-fold ratio\cr
}

The ROC plot illustrates the performance of an experiment as the scoring
threshold is varied. The curve is created by plotting the true positive
rate (TPR) against the false positive rate (FPR) at various threshold settings.

The plot provides diagnose performance to select possibly optimal models and
to discard suboptimal ones. In particularly, the AUC statistics indicate the
performance of the model relatively to a random experiment (AUC 0.5). 
}

\value{
The functions does not return anything but it prints a ROC plot and it's AUC
statistics.
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

plotROC(anaquin, title='ROC Plot', refRats=0)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}