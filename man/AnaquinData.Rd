\name{AnaquinData}
\alias{AnaquinData}
\title{Create Anaquin dataset}
\description{
Create an Anaquin dataset. The resulting S4 object is required for every
Anaquin analysis.
}
\usage{
AnaquinData(analysis, ...)
}
\arguments{
  \item{analysis}{Type of analysis}
  \item{...}{Context specific inputs}
}
\details{
Every Anaquin analysis requires an input dataset. The function enables that
is `\code{AnaquinData}. AnaquinData is a S4 constructor, inputs vary depending
on the analysis.

The constructor requires the following compulsory inputs:

\tabular{ll}{
   \code{analysis} \tab Type of analysis \cr
   \code{seqs}    \tab Sequin names  \cr
}

\code{analysis} gives the type of the analysis; how a statistical model should
build for sequins. The current release supports the following analysis:

\tabular{ll}{
   \code{PlotLinear}   \tab Linear model with sequins \cr
   \code{PlotLogistic} \tab GLM logistic model with sequins  \cr
   \code{PlotROC}      \tab ROC analysis with sequins  \cr
   \code{PlotLODR}     \tab LODR (LOESS) analysis with sequins  \cr
}

\code{seqs} gives the sequin names; critical information otherwise no
statistical model for sequins is possible. \code{AnaquinData} will give
an error message unless both \code{analysis} and \code{seqs} are given.

The constructor also accepts the following optional inputs:

\tabular{ll}{
   \code{sd}       \tab Standard deviation \cr
   \code{pval}     \tab P-value probability  \cr
   \code{qval}     \tab Q-value probability  \cr
   \code{ratio}    \tab Expected sequin ratio (eg: allele frequency)  \cr
   \code{input}    \tab Input concentration (attomol/ul)  \cr
   \code{measured} \tab Measured variable (eg: FPKM)  \cr
   \code{label}    \tab Classified labels (eg: TP, FP)  \cr
   \code{score}    \tab How to rank ROC points? \cr
}

It is important to check the documentation for individual analysis on what the
required inputs. We will illustrate an example for plotting a scatter plot
between input concentration against measured FPKM values.
}
\value{
An S4 object of class \code{AnaquinData}.
}
\examples{
#
# Example 1: Create Anaquin dataset for PlotLogistic
#

data('seqCuffcompare')

title <- 'Assembly Plot'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'Sensitivity'

# Sequin names
seqs <- row.names(seqCuffcompare)

# Input concentration
input <- log2(seqCuffcompare$InputConcent)

# Measured sensitivity
measured <- seqCuffcompare$Sn

anaquin <- AnaquinData(analysis='PlotLogistic',
                           seqs=seqs,
                          input=input,
                       measured=measured)

plotLogistic(anaquin, title=title, xlab=xlab, ylab=ylab, showLOA=TRUE)

#
# Example 2: Create Anaquin dataset for PlotLinear
#

data(seqCufflinks)

title <- 'Gene Expression'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'FPKM (log2)'

# Sequin names
seqs <- row.names(seqCufflinks)

# Input concentration
input <- log2(seqCufflinks$InputConcent)

# Measured FPKKM
measured <- log2(seqCufflinks[,2:4])

anaquin <- AnaquinData(analysis='PlotLinear',
                           seqs=seqs,
                          input=input,
                       measured=measured)

plotLinear(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)

#
# Example 3: Create Anaquin dataset for plotROC
#

data('seqDESeq2')

# Sequin names
seqs <- row.names(seqDESeq2)

# Expected log-fold
ratio <- seqDESeq2$ExpLFC

# How ROC points are ranked (scoring function)
score <- 1-seqDESeq2$Pval

# Classified labels (TP/FP)
label <- seqDESeq2$Label

anaquin <- AnaquinData(analysis='PlotROC',
                           seqs=seqs,
                          ratio=ratio,
                          score=score,
                          label=label)

plotROC(anaquin, title='ROC Plot', refRats=0)

#
# Example 4: Create Anaquin dataset for plotLODR
#

data('seqDESeq2')

xlab  <- 'Average Counts'
ylab  <- 'P-value'
title <- 'LODR Curves'

# Sequin names
seqs <- row.names(seqDESeq2)

# Measured average mean
measured <- seqDESeq2$Mean

# Expected log-fold
ratio <- seqDESeq2$ExpLFC

# P-value
pval <- seqDESeq2$Pval

# Q-value
qval <- seqDESeq2$Qval

anaquin <- AnaquinData(analysis='PlotLODR',
                           seqs=seqs,
                       measured=measured,
                          ratio=ratio,
                           pval=pval,
                           qval=qval)

plotLODR(anaquin, xlab=xlab, ylab=ylab, title=title, FDR=0.1)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}